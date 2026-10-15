# =============================================================================
# phenoIMPACT | Environmental moderators of phenological plasticity:
#               WITHIN- vs BETWEEN-species decomposition
# Version 1.1 | 2026-10-14
#
# PURPOSE
#   Refit the THREE final 15-km phenology-plasticity models using a Mundlak
#   decomposition of the four environmental moderators:
#
#     moderator_population =
#       moderator_within_species + moderator_species_mean
#
#   This separates:
#     WITHIN species:
#       do populations of the SAME species living in different environmental
#       contexts differ in phenological plasticity?
#
#     BETWEEN species:
#       do species occupying different mean environmental contexts differ in
#       phenological plasticity?
#
# ANALYSES
#   1) onset                  | 60-d temperature window
#   2) offset_univoltine      | 90-d window + additive annual ONSET_mean_z
#   3) offset_multivoltine    | 90-d window + additive annual ONSET_mean_z
#
# MODERATORS
#   photoperiod
#   background temperature
#   temperature predictability
#   warming trend
#
# MODEL
#   phenology ~
#     anomaly * (
#       photo_within + photo_between +
#       background_within + background_between +
#       predictability_within + predictability_between +
#       trend_within + trend_between
#     )
#     [+ ONSET_mean_z for offsets]
#     + SAME random effects as source model
#
# Spatial structure is unchanged:
#   spatial = "off"
#   spatiotemporal = "iid"
#   time = YEAR
#   SAME 15-km mesh
#
# IMPORTANT
#   - Same fitted rows as the FINAL-QC source models.
#   - Same anomaly windows and original moderator scaling.
#   - Species means use EQUAL POPULATION weight.
#   - Within + between reconstruct the original moderator EXACTLY.
#   - Before fitting, the script verifies that setting within = between
#     coefficients reproduces the source fixed-effect linear predictor.
#   - No latitude/longitude terms are included in this PRIMARY model.
#   - No causal mediation is claimed.
#
# SOURCE MODELS
#   Onset:
#     current accepted primary 15-km spatial model.
#
#   Offsets:
#     FINAL-QC onset-adjusted 15-km fits after excluding known
#     OFFSET_mean < LAST_PEAK failures, from:
#       output/phenology_plasticity/offset_with_onset_qc_lastpeak_15km/
#       run_fbf8b3101079
#
# RUN
#   source(file.choose(), encoding = "UTF-8")
#   pheno_env_mundlak$audit_sources()
#   env_mundlak_job <- pheno_env_mundlak$run(monitor = FALSE)
#
# MONITOR
#   pheno_env_mundlak$status()
#
# The three fits run sequentially to protect RAM.
# Completed analyses are reused if the same deterministic run is restarted.
# =============================================================================


pheno_env_mundlak <- local({

  VERSION <- "environment_mundlak_qc_final_spatial_15km_1.2"

  ROOT <- "E:/phenoIMPACT project/code/phenoIMPACT"

  OFFSET_RUN <- file.path(
    ROOT,
    "output", "phenology_plasticity",
    "offset_with_onset_qc_lastpeak_15km",
    "run_fbf8b3101079"
  )

  ANALYSES <- c(
    "onset",
    "offset_univoltine",
    "offset_multivoltine"
  )

  MOD_LABELS <- c(
    photoperiod = "Photoperiod",
    background = "Background temperature",
    predictability = "Temperature predictability",
    trend = "Warming trend"
  )

  state <- new.env(parent = emptyenv())


  # ---------------------------------------------------------------------------
  # Utilities
  # ---------------------------------------------------------------------------

  need <- function(ok, msg) {
    if (!isTRUE(ok)) stop(msg, call. = FALSE)
  }

  csv <- function(x, p) {
    utils::write.csv(x, p, row.names = FALSE, na = "")
  }

  atomic_save <- function(x, p) {

    dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)

    tmp <- tempfile(
      ".writing_",
      tmpdir = dirname(p),
      fileext = ".rds"
    )

    on.exit(unlink(tmp), add = TRUE)

    saveRDS(x, tmp, compress = "gzip")

    bak <- paste0(p, ".previous")

    if (file.exists(bak)) unlink(bak)

    had <- file.exists(p)

    if (had) {
      need(
        file.rename(p, bak),
        paste("Cannot retain previous file:", p)
      )
    }

    if (!file.rename(tmp, p)) {
      if (had && file.exists(bak)) file.rename(bak, p)
      stop("Cannot finalize: ", p, call. = FALSE)
    }

    if (file.exists(bak)) unlink(bak)

    invisible(p)
  }


  phase <- function(
    cfg,
    stage,
    analysis = "",
    detail = "",
    terminal = FALSE
  ) {

    rec <- list(
      stage = stage,
      analysis = analysis,
      detail = detail,
      time = Sys.time(),
      terminal = terminal
    )

    atomic_save(
      rec,
      file.path(cfg$out, "progress.rds")
    )

    message(
      "[", format(rec$time, "%Y-%m-%d %H:%M:%S"), "] ",
      stage,
      if (nzchar(analysis)) paste0(" | ", analysis) else "",
      if (nzchar(detail)) paste0(" | ", detail) else ""
    )
  }


  check_packages <- function() {

    pkgs <- c(
      "sdmTMB", "TMB", "Matrix", "fmesher",
      "lme4", "digest", "callr", "ps"
    )

    missing <- pkgs[
      !vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
    ]

    need(
      !length(missing),
      paste(
        "Missing package(s):",
        paste(missing, collapse = ", ")
      )
    )

    need(
      as.character(utils::packageVersion("sdmTMB")) == "1.1.0",
      "Use the original sdmTMB 1.1.0 library."
    )
  }


  disk_check <- function(path, minimum_GiB = 50) {

    d <- ps::ps_disk_usage(path)

    need(
      nrow(d) == 1L && is.finite(d$available),
      "Cannot check free disk space."
    )

    free <- as.numeric(d$available) / 1024^3

    message(
      "Disk preflight: ",
      round(free, 1),
      " GiB available."
    )

    need(
      free >= minimum_GiB,
      paste(
        "Only", round(free, 1),
        "GiB free; require at least",
        minimum_GiB, "GiB."
      )
    )

    invisible(free)
  }


  text_formula <- function(x) {
    paste(
      deparse(x, width.cutoff = 500L),
      collapse = " "
    )
  }


  # ---------------------------------------------------------------------------
  # Analysis specification
  # ---------------------------------------------------------------------------

  spec <- function(analysis) {

    need(
      analysis %in% ANALYSES,
      "Unknown analysis."
    )

    if (analysis == "onset") {

      w <- 60L

      list(
        response = "ONSET_mean",
        window = w,
        anomaly = paste0("clim_anomaly_tw", w),
        control = character(),
        orientation = -1,
        orientation_label = "advancement"
      )

    } else if (analysis == "offset_univoltine") {

      w <- 90L

      list(
        response = "OFFSET_mean",
        window = w,
        anomaly = paste0("clim_anomaly_tw", w),
        control = "ONSET_mean_z",
        orientation = -1,
        orientation_label = "advancement"
      )

    } else {

      w <- 90L

      list(
        response = "OFFSET_mean",
        window = w,
        anomaly = paste0("clim_anomaly_tw", w),
        control = "ONSET_mean_z",
        orientation = 1,
        orientation_label = "delay / extension"
      )
    }
  }


  moderators <- function(w) {

    c(
      photoperiod = paste0("photo_tw", w),
      background = paste0("clim_background_tw", w),
      predictability = paste0("clim_predictability_tw", w),
      trend = paste0("clim_trend_tw", w)
    )
  }


  # ---------------------------------------------------------------------------
  # Locate source onset model
  # ---------------------------------------------------------------------------

  locate_onset <- function(root) {

    spatial_root <- file.path(
      root,
      "output", "phenology_plasticity",
      "spatial_models"
    )

    dirs <- list.dirs(
      spatial_root,
      recursive = FALSE,
      full.names = TRUE
    )

    dirs <- dirs[
      startsWith(
        basename(dirs),
        "primary__onset_"
      )
    ]

    required <- c(
      "prepared_input.rds",
      "fit_spatial_field.rds",
      "metadata.rds"
    )

    dirs <- dirs[
      vapply(
        dirs,
        function(d) {
          all(
            file.exists(
              file.path(d, required)
            )
          )
        },
        logical(1)
      )
    ]

    need(
      length(dirs) == 1L,
      paste(
        "Expected exactly one primary onset spatial model; found",
        length(dirs)
      )
    )

    normalizePath(
      dirs,
      winslash = "/",
      mustWork = TRUE
    )
  }


  # ---------------------------------------------------------------------------
  # Read exact FINAL source model/input for one analysis
  # ---------------------------------------------------------------------------

  read_source <- function(cfg, analysis) {

    s <- spec(analysis)

    if (analysis == "onset") {

      d <- locate_onset(cfg$root)

      p <- readRDS(
        file.path(d, "prepared_input.rds")
      )

      saved <- readRDS(
        file.path(d, "fit_spatial_field.rds")
      )

      fit <- saved$fit

      need(
        inherits(fit, "sdmTMB"),
        "Onset source fit is not an sdmTMB object."
      )

      need(
        identical(p$response, s$response) &&
          identical(p$anomaly, s$anomaly),
        "Onset source response/window mismatch."
      )

      data <- p$data
      formula <- p$formula
      mesh <- fit$spde

      source_path <- file.path(
        d,
        "fit_spatial_field.rds"
      )

    } else {

      folder <- file.path(
        cfg$offset_run,
        analysis
      )

      input_file <- file.path(
        folder,
        "input.rds"
      )

      fit_file <- file.path(
        folder,
        "fit_with_onset.rds"
      )

      completed_file <- file.path(
        folder,
        "completed.rds"
      )

      need(
        all(
          file.exists(
            c(
              input_file,
              fit_file,
              completed_file
            )
          )
        ),
        paste(
          "Final onset-adjusted source missing for",
          analysis
        )
      )

      done <- readRDS(completed_file)

      need(
        identical(done$status, "FIT_OK"),
        paste(
          analysis,
          "onset-adjusted source is not FIT_OK."
        )
      )

      input <- readRDS(input_file)
      saved <- readRDS(fit_file)

      fit <- saved$fit

      need(
        inherits(fit, "sdmTMB"),
        paste(
          analysis,
          "source fit is not an sdmTMB object."
        )
      )

      data <- input$data
      formula <- input$formula
      mesh <- input$mesh

      source_path <- fit_file
    }

    need(
      identical(fit$family$family, "gaussian") &&
        identical(fit$family$link, "identity"),
      "Expected Gaussian identity-link source."
    )

    need(
      s$response %in% names(data) &&
        s$anomaly %in% names(data),
      "Source data lack response or anomaly."
    )

    need(
      all(
        moderators(s$window) %in% names(data)
      ),
      "Source data lack one or more moderators."
    )

    if (length(s$control)) {
      need(
        all(s$control %in% names(data)),
        "Offset source lacks ONSET_mean_z."
      )
    }

    list(
      fit = fit,
      data = data,
      formula = formula,
      mesh = mesh,
      source_path = normalizePath(
        source_path,
        winslash = "/",
        mustWork = TRUE
      )
    )
  }


  # ---------------------------------------------------------------------------
  # Recover source fixed effects + covariance
  # ---------------------------------------------------------------------------

  source_fixed <- function(source) {

    Xnames <- colnames(
      stats::model.matrix(
        lme4::nobars(source$formula),
        source$data
      )
    )

    ix <- which(
      names(source$fit$sd_report$par.fixed) == "b_j"
    )

    need(
      length(ix) == length(Xnames),
      "Cannot align source fixed coefficients."
    )

    b <- stats::setNames(
      as.numeric(
        source$fit$sd_report$par.fixed[ix]
      ),
      Xnames
    )

    V <- as.matrix(
      source$fit$sd_report$cov.fixed[
        ix, ix,
        drop = FALSE
      ]
    )

    dimnames(V) <- list(
      Xnames,
      Xnames
    )

    need(
      all(is.finite(b)) &&
        all(is.finite(V)) &&
        all(diag(V) > 0),
      "Invalid source fixed-effect covariance."
    )

    list(
      beta = b,
      V = V
    )
  }


  # ---------------------------------------------------------------------------
  # Exact random-effects string from source formula
  # ---------------------------------------------------------------------------

  random_terms_text <- function(formula) {

    bars <- lme4::findbars(formula)

    need(
      length(bars) > 0L,
      "Source formula has no random terms."
    )

    paste(
      vapply(
        bars,
        function(x) {
          paste0(
            "(",
            paste(
              deparse(
                x,
                width.cutoff = 500L
              ),
              collapse = " "
            ),
            ")"
          )
        },
        character(1)
      ),
      collapse = " + "
    )
  }


  # ---------------------------------------------------------------------------
  # Decompose moderators within/between species
  # ---------------------------------------------------------------------------

  decompose_moderators <- function(
    data,
    s
  ) {

    mods <- moderators(s$window)

    required <- c(
      "SPECIES",
      "SITE_ID",
      unname(mods)
    )

    need(
      all(required %in% names(data)),
      "Source data lack SPECIES/SITE_ID/moderators."
    )

    d <- data

    d$SPECIES <- as.character(d$SPECIES)
    d$SITE_ID <- as.character(d$SITE_ID)

    pop_key <- interaction(
      d$SPECIES,
      d$SITE_ID,
      drop = TRUE,
      lex.order = TRUE
    )

    # Verify each moderator is a population-context variable,
    # i.e. invariant across years for the same species x site.
    for (v in unname(mods)) {

      ranges <- tapply(
        d[[v]],
        pop_key,
        function(x) {
          diff(range(x))
        }
      )

      need(
        max(
          abs(ranges),
          na.rm = TRUE
        ) < 1e-10,
        paste(
          "Moderator varies among years within SPECIES x SITE_ID:",
          v,
          "\nThis requires a different decomposition."
        )
      )
    }

    pop <- unique(
      d[c(
        "SPECIES",
        "SITE_ID",
        unname(mods)
      )]
    )

    need(
      !anyDuplicated(
        pop[c("SPECIES", "SITE_ID")]
      ),
      "Population moderator table is duplicated."
    )

    scaling <- list()

    for (nm in names(mods)) {

      v <- mods[[nm]]

      sp_mean <- tapply(
        pop[[v]],
        pop$SPECIES,
        mean
      )

      between_name <- paste0(
        nm,
        "_species_mean"
      )

      within_name <- paste0(
        nm,
        "_within_species"
      )

      pop[[between_name]] <-
        unname(
          sp_mean[pop$SPECIES]
        )

      pop[[within_name]] <-
        pop[[v]] -
        pop[[between_name]]

      reconstruction <-
        pop[[within_name]] +
        pop[[between_name]] -
        pop[[v]]

      need(
        max(abs(reconstruction)) < 1e-10,
        paste(
          "Within/between reconstruction failed:",
          nm
        )
      )

      within_sum <- rowsum(
        pop[[within_name]],
        pop$SPECIES,
        reorder = FALSE
      )

      need(
        max(abs(within_sum)) < 1e-8,
        paste(
          "Species-centering failed:",
          nm
        )
      )

      scaling[[nm]] <- data.frame(
        moderator = nm,
        source_variable = v,
        source_SD = stats::sd(pop[[v]]),
        within_SD = stats::sd(
          pop[[within_name]]
        ),
        between_population_weighted_SD =
          stats::sd(
            pop[[between_name]]
          ),
        stringsAsFactors = FALSE
      )
    }

    join_key <- paste(
      d$SPECIES,
      d$SITE_ID,
      sep = "\r"
    )

    pop_join_key <- paste(
      pop$SPECIES,
      pop$SITE_ID,
      sep = "\r"
    )

    m <- match(
      join_key,
      pop_join_key
    )

    need(
      !anyNA(m),
      "Population decomposition join failed."
    )

    new_cols <- grep(
      "_(within_species|species_mean)$",
      names(pop),
      value = TRUE
    )

    for (v in new_cols) {
      d[[v]] <- pop[[v]][m]
    }

    list(
      data = d,
      populations = pop,
      scaling = do.call(
        rbind,
        scaling
      )
    )
  }


  # ---------------------------------------------------------------------------
  # Build Mundlak formula
  # ---------------------------------------------------------------------------

  make_formula <- function(
    source_formula,
    s
  ) {

    fixed_moderators <- c(
      "photoperiod_within_species",
      "photoperiod_species_mean",
      "background_within_species",
      "background_species_mean",
      "predictability_within_species",
      "predictability_species_mean",
      "trend_within_species",
      "trend_species_mean"
    )

    fixed <- paste0(
      s$anomaly,
      " * (",
      paste(
        fixed_moderators,
        collapse = " + "
      ),
      ")"
    )

    if (length(s$control)) {
      fixed <- paste(
        fixed,
        "+",
        paste(
          s$control,
          collapse = " + "
        )
      )
    }

    rhs_random <- random_terms_text(
      source_formula
    )

    stats::as.formula(
      paste0(
        s$response,
        " ~ ",
        fixed,
        " + ",
        rhs_random
      ),
      env = baseenv()
    )
  }


  # ---------------------------------------------------------------------------
  # Map source coefficients to exact equivalent start
  # ---------------------------------------------------------------------------

  map_start_beta <- function(
    source_beta,
    new_formula,
    data,
    s
  ) {

    new_names <- colnames(
      stats::model.matrix(
        lme4::nobars(new_formula),
        data
      )
    )

    beta <- stats::setNames(
      rep(0, length(new_names)),
      new_names
    )

    # Shared terms: intercept, anomaly, ONSET_mean_z.
    shared <- intersect(
      names(source_beta),
      names(beta)
    )

    beta[shared] <- source_beta[shared]

    old_mods <- moderators(s$window)

    for (nm in names(old_mods)) {

      old <- old_mods[[nm]]

      within <- paste0(
        nm,
        "_within_species"
      )

      between <- paste0(
        nm,
        "_species_mean"
      )

      need(
        old %in% names(source_beta),
        paste(
          "Missing source moderator coefficient:",
          old
        )
      )

      beta[within] <- source_beta[old]
      beta[between] <- source_beta[old]

      old_int <- intersect(
        c(
          paste0(s$anomaly, ":", old),
          paste0(old, ":", s$anomaly)
        ),
        names(source_beta)
      )

      need(
        length(old_int) == 1L,
        paste(
          "Missing source anomaly interaction:",
          old
        )
      )

      for (newvar in c(within, between)) {

        new_int <- intersect(
          c(
            paste0(s$anomaly, ":", newvar),
            paste0(newvar, ":", s$anomaly)
          ),
          names(beta)
        )

        need(
          length(new_int) == 1L,
          paste(
            "Missing new anomaly interaction:",
            newvar
          )
        )

        beta[new_int] <-
          source_beta[old_int]
      }
    }

    beta
  }


  # ---------------------------------------------------------------------------
  # Numerical checks
  # ---------------------------------------------------------------------------

  fit_checks <- function(fit) {

    grad <- tryCatch(
      as.numeric(
        fit$tmb_obj$gr(
          fit$model$par
        )
      ),
      error = function(e) NA_real_
    )

    max_grad <- if (
      length(grad) &&
      all(is.finite(grad))
    ) {
      max(abs(grad))
    } else {
      Inf
    }

    convergence <- if (
      !is.null(fit$model$convergence)
    ) {
      fit$model$convergence
    } else {
      NA_integer_
    }

    pd <- isTRUE(
      fit$sd_report$pdHess
    )

    ok <-
      identical(convergence, 0L) &&
      pd &&
      is.finite(max_grad) &&
      max_grad < 0.001

    data.frame(
      convergence_code = convergence,
      pdHess = pd,
      max_abs_gradient = max_grad,
      numerical_checks_ok = ok,
      logLik = as.numeric(
        stats::logLik(fit)
      ),
      AIC = stats::AIC(fit),
      BIC = stats::BIC(fit),
      nobs = nrow(fit$data),
      stringsAsFactors = FALSE
    )
  }


  # ---------------------------------------------------------------------------
  # Fixed-effect table
  # ---------------------------------------------------------------------------

  fixed_table <- function(
    fit,
    formula,
    data
  ) {

    terms <- colnames(
      stats::model.matrix(
        lme4::nobars(formula),
        data
      )
    )

    ix <- which(
      names(fit$sd_report$par.fixed) == "b_j"
    )

    need(
      length(ix) == length(terms),
      "Cannot align fitted fixed effects."
    )

    b <- stats::setNames(
      as.numeric(
        fit$sd_report$par.fixed[ix]
      ),
      terms
    )

    V <- as.matrix(
      fit$sd_report$cov.fixed[
        ix, ix,
        drop = FALSE
      ]
    )

    dimnames(V) <- list(
      terms,
      terms
    )

    se <- sqrt(diag(V))

    tab <- data.frame(
      term = terms,
      estimate = as.numeric(b),
      std.error = as.numeric(se),
      conf.low = as.numeric(
        b - 1.96 * se
      ),
      conf.high = as.numeric(
        b + 1.96 * se
      ),
      p_Wald = as.numeric(
        2 * stats::pnorm(
          -abs(b / se)
        )
      ),
      stringsAsFactors = FALSE
    )

    list(
      beta = b,
      V = V,
      table = tab
    )
  }


  # ---------------------------------------------------------------------------
  # Linear contrast helper
  # ---------------------------------------------------------------------------

  contrast <- function(
    beta,
    V,
    weights,
    label
  ) {

    cvec <- rep(
      0,
      length(beta)
    )

    names(cvec) <- names(beta)

    need(
      all(names(weights) %in% names(beta)),
      paste(
        "Missing contrast terms:",
        label
      )
    )

    cvec[names(weights)] <-
      weights

    est <- sum(
      cvec * beta
    )

    variance <- as.numeric(
      crossprod(
        cvec,
        V %*% cvec
      )
    )

    need(
      is.finite(variance) &&
        variance >= -1e-10,
      paste(
        "Invalid contrast variance:",
        label
      )
    )

    se <- sqrt(
      max(variance, 0)
    )

    data.frame(
      contrast = label,
      estimate = est,
      std.error = se,
      conf.low = est - 1.96 * se,
      conf.high = est + 1.96 * se,
      p_Wald = 2 * stats::pnorm(
        -abs(est / se)
      ),
      stringsAsFactors = FALSE
    )
  }


  # ---------------------------------------------------------------------------
  # Prepare one analysis
  # ---------------------------------------------------------------------------

  prepare_one <- function(
    cfg,
    analysis
  ) {

    phase(
      cfg,
      "PREPARING",
      analysis,
      "Reading FINAL-QC source model"
    )

    s <- spec(analysis)

    source <- read_source(
      cfg,
      analysis
    )

    fx <- source_fixed(
      source
    )

    dec <- decompose_moderators(
      source$data,
      s
    )

    formula <- make_formula(
      source$formula,
      s
    )

    X <- stats::model.matrix(
      lme4::nobars(formula),
      dec$data
    )

    need(
      all(is.finite(X)),
      paste(
        analysis,
        "new design contains non-finite values."
      )
    )

    need(
      qr(X)$rank == ncol(X),
      paste(
        analysis,
        "new fixed-effect design is rank deficient."
      )
    )

    beta_start <- map_start_beta(
      fx$beta,
      formula,
      dec$data,
      s
    )

    X_old <- stats::model.matrix(
      lme4::nobars(source$formula),
      source$data
    )

    eta_old <- as.numeric(
      X_old %*% fx$beta
    )

    eta_new_start <- as.numeric(
      X %*% beta_start
    )

    reconstruction_error <- max(
      abs(
        eta_old -
        eta_new_start
      )
    )

    need(
      reconstruction_error < 1e-8,
      paste(
        analysis,
        "within/between start does not exactly reproduce source fixed predictor."
      )
    )

    start <- source$fit$parlist

    need(
      is.list(start) &&
        "b_j" %in% names(start),
      paste(
        analysis,
        "source parlist unavailable."
      )
    )

    start$b_j <- unname(
      beta_start
    )

    prepared <- list(
      version = VERSION,
      analysis = analysis,
      spec = s,
      data = dec$data,
      populations = dec$populations,
      scaling = dec$scaling,
      formula = formula,
      source_formula = source$formula,
      source_path = source$source_path,
      mesh = source$mesh,
      start = start,
      source_beta = fx$beta,
      source_V = fx$V,
      beta_start = beta_start,
      reconstruction_error = reconstruction_error
    )

    out <- file.path(
      cfg$out,
      analysis
    )

    dir.create(
      out,
      recursive = TRUE,
      showWarnings = FALSE
    )

    atomic_save(
      prepared,
      file.path(
        out,
        "prepared_input.rds"
      )
    )

    csv(
      dec$scaling,
      file.path(
        out,
        "moderator_component_scaling.csv"
      )
    )

    csv(
      data.frame(
        n_rows = nrow(dec$data),
        n_populations =
          nrow(dec$populations),
        n_species =
          length(
            unique(
              dec$populations$SPECIES
            )
          ),
        reconstruction_max_abs_error =
          reconstruction_error,
        fixed_design_columns =
          ncol(X),
        fixed_design_rank =
          qr(X)$rank,
        stringsAsFactors = FALSE
      ),
      file.path(
        out,
        "preparation_checks.csv"
      )
    )

    csv(
      dec$populations,
      file.path(
        out,
        "population_environment_components.csv"
      )
    )

    # Predictor collinearity at population level.
    component_vars <- grep(
      "_(within_species|species_mean)$",
      names(dec$populations),
      value = TRUE
    )

    R <- stats::cor(
      dec$populations[
        component_vars
      ],
      use = "complete.obs"
    )

    csv(
      data.frame(
        predictor = rownames(R),
        R,
        check.names = FALSE
      ),
      file.path(
        out,
        "environment_component_correlations.csv"
      )
    )

    vif <- lapply(
      component_vars,
      function(v) {

        others <- setdiff(
          component_vars,
          v
        )

        mod <- stats::lm(
          stats::reformulate(
            others,
            response = v
          ),
          data = dec$populations
        )

        r2 <- summary(mod)$r.squared

        data.frame(
          predictor = v,
          R2 = r2,
          VIF = 1 / (1 - r2),
          stringsAsFactors = FALSE
        )
      }
    )

    csv(
      do.call(rbind, vif),
      file.path(
        out,
        "environment_component_VIF.csv"
      )
    )

    writeLines(
      c(
        paste(
          "SOURCE:",
          text_formula(source$formula)
        ),
        "",
        paste(
          "MUNDLAK:",
          text_formula(formula)
        ),
        "",
        paste(
          "Exact source fixed predictor reconstruction error:",
          format(
            reconstruction_error,
            scientific = TRUE
          )
        ),
        "",
        "Species means use equal population weight.",
        "Original moderator scaling is retained.",
        "Within + between = original moderator exactly."
      ),
      file.path(
        out,
        "model_specification.txt"
      )
    )

    invisible(prepared)
  }


  # ---------------------------------------------------------------------------
  # Export inferential results
  # ---------------------------------------------------------------------------

  export_results <- function(
    fit,
    p,
    out
  ) {

    s <- p$spec

    fx <- fixed_table(
      fit,
      p$formula,
      p$data
    )

    checks <- fit_checks(
      fit
    )

    csv(
      checks,
      file.path(
        out,
        "fit_checks.csv"
      )
    )

    csv(
      fx$table,
      file.path(
        out,
        "fixed_effects.csv"
      )
    )

    csv(
      data.frame(
        term = rownames(fx$V),
        fx$V,
        check.names = FALSE
      ),
      file.path(
        out,
        "fixed_effect_covariance.csv"
      )
    )

    # Extract the 8 anomaly x moderator-component interactions.
    rows <- list()
    diffs <- list()

    for (nm in names(MOD_LABELS)) {

      within <- paste0(
        nm,
        "_within_species"
      )

      between <- paste0(
        nm,
        "_species_mean"
      )

      within_term <- intersect(
        c(
          paste0(s$anomaly, ":", within),
          paste0(within, ":", s$anomaly)
        ),
        fx$table$term
      )

      between_term <- intersect(
        c(
          paste0(s$anomaly, ":", between),
          paste0(between, ":", s$anomaly)
        ),
        fx$table$term
      )

      need(
        length(within_term) == 1L &&
          length(between_term) == 1L,
        paste(
          "Cannot identify Mundlak interaction terms:",
          nm
        )
      )

      for (scale in c(
        within_species = within_term,
        between_species = between_term
      )) {

        z <- fx$table[
          match(
            scale,
            fx$table$term
          ),
          ,
          drop = FALSE
        ]

        z$analysis <- p$analysis
        z$moderator <- MOD_LABELS[[nm]]
        z$effect_scale <-
          if (
            identical(
              names(scale),
              "within_species"
            )
          ) {
            "within species"
          } else {
            "between species"
          }

        # Biological orientation:
        # positive = stronger onset/uni-offset advancement,
        # or stronger multi-offset delay/extension.
        z$estimate_biological <-
          s$orientation *
          z$estimate

        z$conf.low_biological <-
          min(
            s$orientation *
              c(
                z$conf.low,
                z$conf.high
              )
          )

        z$conf.high_biological <-
          max(
            s$orientation *
              c(
                z$conf.low,
                z$conf.high
              )
          )

        z$biological_response <-
          s$orientation_label

        rows[[length(rows) + 1L]] <- z
      }

      d <- contrast(
        fx$beta,
        fx$V,
        weights = stats::setNames(
          c(1, -1),
          c(
            within_term,
            between_term
          )
        ),
        label = paste0(
          nm,
          "_within_minus_between"
        )
      )

      d$analysis <- p$analysis
      d$moderator <- MOD_LABELS[[nm]]
      d$estimate_biological <-
        s$orientation *
        d$estimate

      diffs[[length(diffs) + 1L]] <- d
    }

    key <- do.call(
      rbind,
      rows
    )

    key <- key[
      c(
        "analysis",
        "moderator",
        "effect_scale",
        "term",
        "estimate",
        "std.error",
        "conf.low",
        "conf.high",
        "p_Wald",
        "estimate_biological",
        "conf.low_biological",
        "conf.high_biological",
        "biological_response"
      )
    ]

    csv(
      key,
      file.path(
        out,
        "environment_within_between_effects.csv"
      )
    )

    diff_tab <- do.call(
      rbind,
      diffs
    )

    csv(
      diff_tab,
      file.path(
        out,
        "within_minus_between_contrasts.csv"
      )
    )

    # Rescale interaction effects by observed SD of within/between components.
    scaled_rows <- list()

    for (i in seq_len(nrow(key))) {

      mod_key <- names(MOD_LABELS)[
        MOD_LABELS ==
          key$moderator[i]
      ]

      sd_row <- p$scaling[
        p$scaling$moderator ==
          mod_key,
        ,
        drop = FALSE
      ]

      component_sd <- if (
        key$effect_scale[i] ==
          "within species"
      ) {
        sd_row$within_SD
      } else {
        sd_row$between_population_weighted_SD
      }

      scaled_rows[[i]] <- data.frame(
        analysis = key$analysis[i],
        moderator = key$moderator[i],
        effect_scale = key$effect_scale[i],
        component_SD = component_sd,
        raw_interaction_estimate =
          key$estimate[i],
        biological_effect_per_observed_component_SD =
          key$estimate_biological[i] *
          component_sd,
        biological_low_per_observed_component_SD =
          key$conf.low_biological[i] *
          component_sd,
        biological_high_per_observed_component_SD =
          key$conf.high_biological[i] *
          component_sd,
        stringsAsFactors = FALSE
      )
    }

    csv(
      do.call(
        rbind,
        scaled_rows
      ),
      file.path(
        out,
        "effects_per_observed_component_SD.csv"
      )
    )

    # Source-vs-Mundlak nesting via exact equality constraints is documented
    # by preparation reconstruction; do not use a naive LRT because equality
    # constraints require a dedicated constrained refit for exact inference.

    sanity <- tryCatch(
      capture.output(
        sdmTMB::sanity(
          fit,
          gradient_thresh = .001
        )
      ),
      error = function(e) {
        paste(
          "Sanity export failed:",
          conditionMessage(e)
        )
      }
    )

    writeLines(
      sanity,
      file.path(
        out,
        "sanity.txt"
      )
    )

    tryCatch(
      csv(
        as.data.frame(
          sdmTMB::tidy(
            fit,
            effects = "ran_pars",
            conf.int = TRUE
          )
        ),
        file.path(
          out,
          "random_parameters.csv"
        )
      ),
      error = function(e) {
        writeLines(
          conditionMessage(e),
          file.path(
            out,
            "random_parameter_export_error.txt"
          )
        )
      }
    )

    list(
      checks = checks,
      key = key,
      differences = diff_tab
    )
  }


  # ---------------------------------------------------------------------------
  # Fit one analysis
  # ---------------------------------------------------------------------------

  fit_one <- function(
    cfg,
    analysis
  ) {

    out <- file.path(
      cfg$out,
      analysis
    )

    p <- readRDS(
      file.path(
        out,
        "prepared_input.rds"
      )
    )

    final_file <- file.path(
      out,
      "fit_mundlak.rds"
    )

    done_file <- file.path(
      out,
      "completed.rds"
    )

    if (
      file.exists(done_file) &&
      file.exists(final_file)
    ) {

      done <- readRDS(done_file)

      need(
        identical(
          done$signature,
          cfg$signature
        ),
        "Completed cache signature mismatch."
      )

      phase(
        cfg,
        "REUSING",
        analysis,
        done$status
      )

      return(done$status)
    }

    if (file.exists(final_file)) {

      saved <- readRDS(final_file)

      need(
        identical(
          saved$signature,
          cfg$signature
        ),
        "Saved fit signature mismatch."
      )

      phase(
        cfg,
        "EXPORTING_CACHED",
        analysis,
        "Optimization already completed"
      )

      ex <- export_results(
        saved$fit,
        p,
        out
      )

      status <- if (
        isTRUE(
          ex$checks$numerical_checks_ok
        )
      ) {
        "FIT_OK"
      } else {
        "NUMERICAL_REVIEW_REQUIRED"
      }

      atomic_save(
        list(
          signature = cfg$signature,
          status = status,
          time = Sys.time()
        ),
        done_file
      )

      return(status)
    }

    phase(
      cfg,
      "FITTING",
      analysis,
      paste(
        nrow(p$data),
        "rows | Mundlak environmental decomposition"
      )
    )

    options(
      sdmTMB.cores = 1L
    )

    TMB::openmp(
      n = 1L,
      DLL = "sdmTMB"
    )

    warnings <- character()

    catcher <- function(w) {

      warnings <<- unique(
        c(
          warnings,
          conditionMessage(w)
        )
      )

      message(
        "WARNING: ",
        conditionMessage(w)
      )

      invokeRestart(
        "muffleWarning"
      )
    }

    set.seed(20261014L)

    fit <- withCallingHandlers(
      sdmTMB::sdmTMB(
        formula = p$formula,
        data = p$data,
        mesh = p$mesh,
        time = "YEAR",
        family = stats::gaussian(
          link = "identity"
        ),
        spatial = "off",
        spatiotemporal = "iid",
        reml = FALSE,
        silent = FALSE,
        control = sdmTMB::sdmTMBcontrol(
          start = p$start,
          multiphase = FALSE,
          nlminb_loops = 1L,
          newton_loops = 1L,
          parallel = 1L,
          get_joint_precision = FALSE,
          collapse_spatial_variance = FALSE
        )
      ),
      warning = catcher
    )

    # Synchronize all point-estimate caches before saving.
    fit$tmb_obj$fn(
      fit$model$par
    )

    fit$parlist <-
      fit$tmb_obj$env$parList(
        par =
          fit$tmb_obj$env$last.par.best
      )

    fit$last.par.best <-
      fit$tmb_obj$env$last.par.best

    atomic_save(
      list(
        version = VERSION,
        signature = cfg$signature,
        analysis = analysis,
        fit = fit,
        warnings = warnings
      ),
      final_file
    )

    writeLines(
      unique(warnings),
      file.path(
        out,
        "warnings.txt"
      )
    )

    phase(
      cfg,
      "EXPORTING",
      analysis,
      "Numerical checks and within/between effects"
    )

    ex <- export_results(
      fit,
      p,
      out
    )

    status <- if (
      isTRUE(
        ex$checks$numerical_checks_ok
      )
    ) {
      "FIT_OK"
    } else {
      "NUMERICAL_REVIEW_REQUIRED"
    }

    atomic_save(
      list(
        signature = cfg$signature,
        status = status,
        time = Sys.time()
      ),
      done_file
    )

    atomic_save(
      fit$parlist,
      file.path(
        out,
        "parameters.rds"
      )
    )

    writeLines(
      capture.output(
        sessionInfo()
      ),
      file.path(
        out,
        "sessionInfo.txt"
      )
    )

    rm(fit)
    invisible(gc())

    phase(
      cfg,
      status,
      analysis,
      "Completed"
    )

    status
  }


  # ---------------------------------------------------------------------------
  # Combine outputs
  # ---------------------------------------------------------------------------

  combine_outputs <- function(
    cfg,
    statuses
  ) {

    effects <- list()
    differences <- list()
    checks <- list()
    scaling <- list()

    for (analysis in ANALYSES) {

      out <- file.path(
        cfg$out,
        analysis
      )

      f1 <- file.path(
        out,
        "environment_within_between_effects.csv"
      )

      f2 <- file.path(
        out,
        "within_minus_between_contrasts.csv"
      )

      f3 <- file.path(
        out,
        "fit_checks.csv"
      )

      f4 <- file.path(
        out,
        "moderator_component_scaling.csv"
      )

      if (file.exists(f1)) {
        effects[[analysis]] <-
          utils::read.csv(
            f1,
            stringsAsFactors = FALSE
          )
      }

      if (file.exists(f2)) {
        differences[[analysis]] <-
          utils::read.csv(
            f2,
            stringsAsFactors = FALSE
          )
      }

      if (file.exists(f3)) {
        x <- utils::read.csv(
          f3,
          stringsAsFactors = FALSE
        )
        x$analysis <- analysis
        checks[[analysis]] <- x
      }

      if (file.exists(f4)) {
        x <- utils::read.csv(
          f4,
          stringsAsFactors = FALSE
        )
        x$analysis <- analysis
        scaling[[analysis]] <- x
      }
    }

    if (length(effects)) {
      csv(
        do.call(rbind, effects),
        file.path(
          cfg$out,
          "ALL_environment_within_between_effects.csv"
        )
      )
    }

    if (length(differences)) {
      csv(
        do.call(rbind, differences),
        file.path(
          cfg$out,
          "ALL_within_minus_between_contrasts.csv"
        )
      )
    }

    if (length(checks)) {
      csv(
        do.call(rbind, checks),
        file.path(
          cfg$out,
          "ALL_fit_checks.csv"
        )
      )
    }

    if (length(scaling)) {
      csv(
        do.call(rbind, scaling),
        file.path(
          cfg$out,
          "ALL_moderator_component_scaling.csv"
        )
      )
    }

    csv(
      statuses,
      file.path(
        cfg$out,
        "workflow_status.csv"
      )
    )
  }


  # ---------------------------------------------------------------------------
  # Worker
  # ---------------------------------------------------------------------------

  worker <- function(cfg) {

    tryCatch({

      check_packages()
      disk_check(
        cfg$out,
        cfg$minimum_free_GiB
      )

      phase(
        cfg,
        "PREPARING",
        "",
        "Preparing all three analyses"
      )

      for (analysis in ANALYSES) {

        prep_file <- file.path(
          cfg$out,
          analysis,
          "prepared_input.rds"
        )

        if (!file.exists(prep_file)) {
          prepare_one(
            cfg,
            analysis
          )
        }
      }

      statuses <- data.frame(
        analysis = character(),
        status = character(),
        stringsAsFactors = FALSE
      )

      # Sequential order: onset -> uni offset -> multi offset
      for (analysis in ANALYSES) {

        disk_check(
          cfg$out,
          20
        )

        status <- fit_one(
          cfg,
          analysis
        )

        statuses <- rbind(
          statuses,
          data.frame(
            analysis = analysis,
            status = status,
            stringsAsFactors = FALSE
          )
        )

        csv(
          statuses,
          file.path(
            cfg$out,
            "workflow_status_partial.csv"
          )
        )

        gc()
      }

      combine_outputs(
        cfg,
        statuses
      )

      writeLines(
        c(
          "ENVIRONMENTAL MODERATORS OF PHENOLOGICAL PLASTICITY — WITHIN/BETWEEN SPECIES",
          "",
          "Three FINAL-QC 15-km Gaussian spatial phenology models:",
          "Offsets use the LAST_PEAK QC source (OFFSET_mean >= LAST_PEAK).",
          "  onset",
          "  offset_univoltine + annual ONSET_mean_z",
          "  offset_multivoltine + annual ONSET_mean_z",
          "",
          "Each of the four environmental moderators is decomposed into:",
          "  population deviation from its species mean (within species)",
          "  species mean environmental context (between species)",
          "",
          "Species means use equal population weight.",
          "The original standardized moderator scale is retained.",
          "Within + between exactly reconstructs the original moderator.",
          "",
          "Primary inference:",
          "  anomaly x moderator_within_species",
          "asks whether populations of the SAME species differ in plasticity",
          "according to environmental context.",
          "",
          "No latitude/longitude is included in these primary fits.",
          "No causal mediation is claimed.",
          "",
          "ALL_environment_within_between_effects.csv is the main result table.",
          "",
          "Positive biological estimates mean:",
          "  onset: stronger advancement",
          "  univoltine offset: stronger advancement",
          "  multivoltine offset: stronger delay / late-season extension"
        ),
        file.path(
          cfg$out,
          "README_environment_mundlak.txt"
        )
      )

      review_files <- list.files(
        cfg$out,
        recursive = TRUE,
        full.names = FALSE
      )

      review_files <- review_files[
        grepl(
          "\\.(csv|txt)$",
          review_files,
          ignore.case = TRUE
        )
      ]

      zip_file <- file.path(
        cfg$out,
        "environment_mundlak_qc_final_results_to_review.zip"
      )

      if (
        requireNamespace(
          "zip",
          quietly = TRUE
        )
      ) {

        zip::zipr(
          zip_file,
          files = review_files,
          root = cfg$out,
          mode = "mirror",
          include_directories = FALSE
        )

      } else {

        old <- getwd()
        on.exit(
          setwd(old),
          add = TRUE
        )

        setwd(cfg$out)

        utils::zip(
          zipfile = zip_file,
          files = review_files
        )
      }

      final_status <- if (
        all(
          statuses$status == "FIT_OK"
        )
      ) {
        "COMPLETED"
      } else {
        "COMPLETED_WITH_NUMERICAL_REVIEW"
      }

      phase(
        cfg,
        final_status,
        "",
        zip_file,
        terminal = TRUE
      )

      list(
        status = final_status,
        output = cfg$out,
        zip = zip_file
      )

    }, error = function(e) {

      original_message <- conditionMessage(e)
      original_call <- tryCatch(
        paste(deparse(conditionCall(e), width.cutoff = 500L), collapse = " "),
        error = function(e2) "<condition call unavailable>"
      )

      message(
        "\nORIGINAL_WORKER_ERROR: ", original_message,
        "\nORIGINAL_WORKER_CALL: ", original_call
      )

      diagnostic_text <- c(
        paste0("ORIGINAL_WORKER_ERROR: ", original_message),
        paste0("ORIGINAL_WORKER_CALL: ", original_call),
        paste0("TIME: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
      )

      try(
        writeLines(
          diagnostic_text,
          file.path(cfg$out, "fatal_error.txt")
        ),
        silent = TRUE
      )

      try(
        phase(
          cfg,
          "STOPPED_WITH_ERROR",
          "",
          original_message,
          terminal = TRUE
        ),
        silent = TRUE
      )

      list(
        status = "STOPPED_WITH_ERROR",
        error = original_message,
        error_call = original_call,
        output = cfg$out
      )
    })
  }


  # ---------------------------------------------------------------------------
  # Dump worker engine
  # ---------------------------------------------------------------------------

  dump_engine <- function(out) {

    e <- environment(dump_engine)

    nm <- ls(
      e,
      all.names = TRUE
    )

    fn <- nm[
      vapply(
        nm,
        function(x) {
          is.function(
            get(
              x,
              envir = e
            )
          )
        },
        logical(1)
      )
    ]

    path <- file.path(
      out,
      "environment_mundlak_engine.R"
    )

    dump(
      c(
        "VERSION",
        "ROOT",
        "OFFSET_RUN",
        "ANALYSES",
        "MOD_LABELS",
        fn
      ),
      file = path,
      envir = e,
      control = "all"
    )

    normalizePath(
      path,
      winslash = "/",
      mustWork = TRUE
    )
  }


  # ---------------------------------------------------------------------------
  # Public source/decomposition audit -- NO fitting
  # ---------------------------------------------------------------------------

  audit_sources <- function(
    project_root = ROOT,
    offset_run = OFFSET_RUN
  ) {

    check_packages()

    project_root <- normalizePath(
      project_root,
      winslash = "/",
      mustWork = TRUE
    )

    offset_run <- normalizePath(
      offset_run,
      winslash = "/",
      mustWork = TRUE
    )

    need(
      grepl("offset_with_onset_qc_lastpeak_15km", offset_run, fixed = TRUE) &&
        basename(offset_run) == "run_fbf8b3101079",
      paste0(
        "Audit is locked to the FINAL-QC LAST_PEAK source: ",
        "offset_with_onset_qc_lastpeak_15km/run_fbf8b3101079"
      )
    )

    cfg <- list(
      root = project_root,
      offset_run = offset_run
    )

    rows <- list()

    for (analysis in ANALYSES) {

      s <- spec(analysis)
      source <- read_source(cfg, analysis)
      fx <- source_fixed(source)
      dec <- decompose_moderators(source$data, s)
      new_formula <- make_formula(source$formula, s)
      beta_start <- map_start_beta(
        fx$beta,
        new_formula,
        dec$data,
        s
      )

      X0 <- stats::model.matrix(
        lme4::nobars(source$formula),
        source$data
      )

      X1 <- stats::model.matrix(
        lme4::nobars(new_formula),
        dec$data
      )

      eta_error <- max(
        abs(
          as.numeric(X0 %*% fx$beta) -
            as.numeric(X1 %*% beta_start)
        )
      )

      need(
        is.finite(eta_error) && eta_error < 1e-8,
        paste(
          analysis,
          "Mundlak constrained start does not reproduce the source fixed predictor."
        )
      )

      pop <- unique(
        dec$data[c("SPECIES", "SITE_ID")]
      )

      conv <- if (
        !is.null(source$fit$model$convergence) &&
          length(source$fit$model$convergence) == 1L
      ) {
        as.integer(source$fit$model$convergence)
      } else {
        NA_integer_
      }

      pd <- isTRUE(source$fit$sd_report$pdHess)

      grad <- tryCatch(
        {
          g <- source$fit$gradients
          if (length(g) && all(is.finite(g))) max(abs(g)) else NA_real_
        },
        error = function(e) NA_real_
      )

      rows[[analysis]] <- data.frame(
        analysis = analysis,
        response = s$response,
        window_days = s$window,
        anomaly = s$anomaly,
        n_rows = nrow(source$data),
        n_populations = nrow(pop),
        n_species = length(unique(as.character(source$data$SPECIES))),
        n_sites = length(unique(as.character(source$data$SITE_ID))),
        n_years = length(unique(source$data$YEAR)),
        source_convergence = conv,
        source_pdHess = pd,
        source_max_abs_gradient = grad,
        source_objective = as.numeric(source$fit$model$objective),
        constrained_fixed_eta_max_abs_error = eta_error,
        source_file = source$source_path,
        stringsAsFactors = FALSE
      )
    }

    audit <- do.call(rbind, rows)
    rownames(audit) <- NULL

    print(
      audit[
        c(
          "analysis",
          "n_rows",
          "n_populations",
          "n_species",
          "n_sites",
          "n_years",
          "window_days",
          "source_convergence",
          "source_pdHess",
          "source_max_abs_gradient",
          "constrained_fixed_eta_max_abs_error"
        )
      ],
      row.names = FALSE
    )

    cat("\nSource files:\n")
    for (i in seq_len(nrow(audit))) {
      cat(
        audit$analysis[i], ":\n  ",
        audit$source_file[i],
        "\n",
        sep = ""
      )
    }

    invisible(audit)
  }


  # ---------------------------------------------------------------------------
  # Public runner
  # ---------------------------------------------------------------------------

  run <- function(
    project_root = ROOT,
    offset_run = OFFSET_RUN,
    minimum_free_GiB = 50,
    monitor = FALSE
  ) {

    check_packages()

    current <- state$job

    if (
      !is.null(current) &&
      current$process$is_alive()
    ) {

      message(
        "Environmental Mundlak workflow is already running."
      )

      return(
        invisible(current)
      )
    }

    project_root <- normalizePath(
      project_root,
      winslash = "/",
      mustWork = TRUE
    )

    offset_run <- normalizePath(
      offset_run,
      winslash = "/",
      mustWork = TRUE
    )

    need(
      grepl("offset_with_onset_qc_lastpeak_15km", offset_run, fixed = TRUE) &&
        basename(offset_run) == "run_fbf8b3101079",
      paste0(
        "This workflow is locked to the FINAL-QC LAST_PEAK source: ",
        "offset_with_onset_qc_lastpeak_15km/run_fbf8b3101079"
      )
    )

    # Run all source/decomposition checks synchronously before launching a
    # multi-hour background fit.
    invisible(
      audit_sources(
        project_root = project_root,
        offset_run = offset_run
      )
    )

    # Deterministic fingerprint from source files + scientific settings.
    onset_dir <- locate_onset(
      project_root
    )

    source_files <- c(
      file.path(
        onset_dir,
        "prepared_input.rds"
      ),
      file.path(
        onset_dir,
        "fit_spatial_field.rds"
      ),
      file.path(
        offset_run,
        "offset_univoltine",
        "input.rds"
      ),
      file.path(
        offset_run,
        "offset_univoltine",
        "fit_with_onset.rds"
      ),
      file.path(
        offset_run,
        "offset_multivoltine",
        "input.rds"
      ),
      file.path(
        offset_run,
        "offset_multivoltine",
        "fit_with_onset.rds"
      )
    )

    need(
      all(file.exists(source_files)),
      "One or more required source files are missing."
    )

    hashes <- unname(
      tools::md5sum(
        source_files
      )
    )

    signature <- digest::digest(
      list(
        version = VERSION,
        source_files = normalizePath(
          source_files,
          winslash = "/"
        ),
        md5 = hashes,
        analyses = ANALYSES,
        moderator_labels = MOD_LABELS,
        species_mean_weighting =
          "equal population",
        geography_in_primary = FALSE
      ),
      algo = "sha256"
    )

    out_root <- file.path(
      project_root,
      "output",
      "phenology_plasticity",
      "environment_mundlak_qc_final_15km"
    )

    dir.create(
      out_root,
      recursive = TRUE,
      showWarnings = FALSE
    )

    out <- file.path(
      out_root,
      paste0(
        "run_",
        substr(signature, 1, 12)
      )
    )

    dir.create(
      out,
      recursive = TRUE,
      showWarnings = FALSE
    )

    out <- normalizePath(
      out,
      winslash = "/",
      mustWork = TRUE
    )

    cfg <- list(
      version = VERSION,
      signature = signature,
      root = project_root,
      offset_run = offset_run,
      out = out,
      minimum_free_GiB =
        minimum_free_GiB,
      source_files =
        normalizePath(
          source_files,
          winslash = "/"
        ),
      source_md5 = hashes,
      libpath = .libPaths()
    )

    cfg$engine <- dump_engine(
      out
    )

    atomic_save(
      cfg,
      file.path(
        out,
        "run_config.rds"
      )
    )

    csv(
      data.frame(
        file = cfg$source_files,
        md5 = cfg$source_md5
      ),
      file.path(
        out,
        "source_manifest.csv"
      )
    )

    phase(
      cfg,
      "QUEUED",
      "",
      paste(
        "Three sequential fits:",
        paste(
          ANALYSES,
          collapse = " -> "
        )
      )
    )

    log_file <- file.path(
      out,
      paste0(
        "controller_",
        format(
          Sys.time(),
          "%Y%m%d_%H%M%S"
        ),
        ".log"
      )
    )

    child <- callr::r_bg(
      function(engine, cfg) {

        e <- new.env(
          parent = globalenv()
        )

        sys.source(
          engine,
          envir = e
        )

        need_same <- identical(
          unname(
            tools::md5sum(
              cfg$source_files
            )
          ),
          cfg$source_md5
        )

        if (!need_same) {
          stop(
            "Source files changed after configuration.",
            call. = FALSE
          )
        }

        e$worker(cfg)
      },
      args = list(
        engine = cfg$engine,
        cfg = cfg
      ),
      libpath = cfg$libpath,
      stdout = log_file,
      stderr = "2>&1",
      supervise = TRUE,
      user_profile = FALSE,
      system_profile = FALSE,
      wd = project_root,
      env = c(
        callr::rcmd_safe_env(),
        OMP_NUM_THREADS = "1",
        OPENBLAS_NUM_THREADS = "1",
        MKL_NUM_THREADS = "1"
      )
    )

    job <- list(
      process = child,
      out = out,
      log = log_file,
      started = Sys.time()
    )

    state$job <- job

    options(
      phenoimpact.env.mundlak.job =
        job
    )

    message(
      "Environmental Mundlak workflow launched."
    )

    message(
      "Output: ",
      out
    )

    message(
      "Monitor with: pheno_env_mundlak$status()"
    )

    if (monitor) {
      status(job)
    }

    invisible(job)
  }


  # ---------------------------------------------------------------------------
  # Status
  # ---------------------------------------------------------------------------

  status <- function(job = NULL) {

    if (is.null(job)) {
      job <- state$job
    }

    if (is.null(job)) {
      job <- getOption(
        "phenoimpact.env.mundlak.job"
      )
    }

    need(
      !is.null(job),
      "No environmental Mundlak job registered."
    )

    progress <- file.path(
      job$out,
      "progress.rds"
    )

    if (file.exists(progress)) {

      p <- tryCatch(
        readRDS(progress),
        error = function(e) NULL
      )

      if (!is.null(p)) {
        print(p)
      }
    }

    alive <- job$process$is_alive()

    cat(
      "\nAlive:",
      alive,
      "\n"
    )

    partial <- file.path(
      job$out,
      "workflow_status_partial.csv"
    )

    if (file.exists(partial)) {

      cat(
        "\nCompleted fits:\n"
      )

      print(
        utils::read.csv(
          partial,
          stringsAsFactors = FALSE
        ),
        row.names = FALSE
      )
    }

    if (!alive) {

      cat(
        "Exit status:",
        job$process$get_exit_status(),
        "\n"
      )

      result_file <- file.path(
        job$out,
        "ALL_environment_within_between_effects.csv"
      )

      if (file.exists(result_file)) {

        cat(
          "\nWithin-between environmental effects:\n"
        )

        print(
          utils::read.csv(
            result_file,
            stringsAsFactors = FALSE
          ),
          row.names = FALSE
        )
      }

      zip_file <- file.path(
        job$out,
        "environment_mundlak_qc_final_results_to_review.zip"
      )

      if (file.exists(zip_file)) {

        cat(
          "\nFILE TO SEND:\n",
          zip_file,
          "\n",
          sep = ""
        )
      }

      error_file <- file.path(
        job$out,
        "fatal_error.txt"
      )

      if (file.exists(error_file)) {

        cat(
          "\nError:\n",
          paste(
            readLines(
              error_file,
              warn = FALSE
            ),
            collapse = "\n"
          ),
          "\n",
          sep = ""
        )
      }
    }

    invisible(job)
  }


  list(
    audit_sources = audit_sources,
    run = run,
    status = status
  )
})


run_environment_mundlak <- function(...) {
  pheno_env_mundlak$run(...)
}


status_environment_mundlak <- function(...) {
  pheno_env_mundlak$status(...)
}


# Explicit FINAL-QC aliases
pheno_env_mundlak_qc <- pheno_env_mundlak
run_environment_mundlak_qc_final <- function(...) pheno_env_mundlak_qc$run(...)
status_environment_mundlak_qc_final <- function(...) pheno_env_mundlak_qc$status(...)
