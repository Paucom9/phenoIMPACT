# =============================================================================
# phenoIMPACT | WITHIN-BETWEEN (MUNDLAK) plasticity-demography bridge
# Version 1.0 | 2026-10-13
#
# SCIENTIFIC QUESTION
#   Is the association between phenological plasticity and population trend
#   expressed WITHIN species (more-plastic populations versus less-plastic
#   populations of the same species), BETWEEN species, or both?
#
# TWO MODELS, fitted sequentially:
#   1) multivoltine
#   2) univoltine
#
# FIXED-EFFECT MODEL
#   ABUND_INDEX ~ year_decade * (
#       latitude_within_species_10deg +
#       longitude_within_species_10deg +
#       onset_plasticity_within_species_z +
#       onset_plasticity_species_mean_z +
#       offset_plasticity_within_species_z +
#       offset_plasticity_species_mean_z
#   )
#
# RANDOM / RESIDUAL STRUCTURE
#   - independent species intercept and temporal slope
#   - independent site intercept and temporal slope
#   - population AR(1), preserving calendar gaps
#   - annual IID Matern spatial field on the existing 15-km mesh
#   - Gamma/log likelihood
#
# PARAMETERIZATION
#   Plasticity variables retain the ORIGINAL group-wide z-score units used in
#   the final bridge models. For each plasticity variable P:
#
#       P_population = P_within_species + mean(P)_species
#
#   No second standardization is used in fitting. This makes the existing bridge
#   model exactly nested in this model under the constraints:
#
#       beta_within = beta_between
#
#   for the onset and offset main effects and their year interactions.
#
# INPUTS
#   Final geography + plasticity bridge run:
#     bridge_geography_plasticity_spatial_15km/run_20261008_100948_5724
#
#   Multivoltine warm start:
#     accepted FIT_OK bridge fit
#
#   Univoltine warm start:
#     accepted FIT_OK numerical adjudication/refinement of the bridge fit
#
# IMPORTANT
#   - Same annual rows, abundance values, plasticity estimates, coordinates,
#     likelihood, AR(1), random effects, mesh and annual field as the final bridge.
#   - Current bridge fits are warm starts only. All parameters are re-estimated.
#   - This is a within-between association model, NOT a causal mediation model.
#   - Checkpoints are saved during optimization. The output directory is
#     deterministic, so rerunning the same command resumes/reuses the same run.
#
# RUN
#   source(file.choose(), encoding = "UTF-8")
#   mundlak_job <- pheno_mundlak$run(monitor = FALSE)
#
# MONITOR
#   pheno_mundlak$status()
#
# Keep R open. Only one large fit runs at a time.
# =============================================================================


pheno_mundlak <- local({

  VERSION <- "mundlak_geography_plasticity_spatial_15km_1.0"
  SOURCE_BRIDGE_VERSION <- "bridge_geography_plasticity_spatial_15km_1.0"
  REFINED_UNIVOLTINE_VERSION <- "univoltine_bridge_refinement_1.0"

  ROOT <- "E:/phenoIMPACT project/code/phenoIMPACT"

  BRIDGE_RUN <- paste0(
    ROOT,
    "/output/population_trends/bridge_geography_plasticity_spatial_15km/",
    "run_20261008_100948_5724"
  )

  GROUPS <- c("multivoltine", "univoltine")
  GRAD_TOL <- 0.05

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
        paste("Cannot retain previous:", p)
      )
    }

    if (!file.rename(tmp, p)) {

      if (had && file.exists(bak)) {
        file.rename(bak, p)
      }

      stop("Cannot finalize: ", p, call. = FALSE)
    }

    if (file.exists(bak)) unlink(bak)

    invisible(p)
  }


  phase <- function(
    cfg,
    stage,
    detail = "",
    group = "",
    terminal = FALSE
  ) {

    z <- list(
      stage = stage,
      group = group,
      detail = detail,
      time = Sys.time(),
      terminal = terminal
    )

    atomic_save(z, cfg$progress)

    message(
      "[", format(z$time, "%Y-%m-%d %H:%M:%S"), "] ",
      stage,
      if (nzchar(group)) paste0(" | ", group) else "",
      if (nzchar(detail)) paste0(" | ", detail) else ""
    )
  }


  check_packages <- function() {

    pkgs <- c(
      "TMB", "Matrix", "fmesher", "sf",
      "callr", "ps", "digest"
    )

    missing <- pkgs[
      !vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
    ]

    need(
      !length(missing),
      paste("Missing package(s):", paste(missing, collapse = ", "))
    )
  }


  disk_check <- function(path, min_gib = 25, verbose = TRUE) {

    d <- ps::ps_disk_usage(path)

    need(
      nrow(d) == 1L && is.finite(d$available),
      "Cannot check free disk space."
    )

    free <- as.numeric(d$available) / 1024^3

    if (verbose) {
      message("Disk preflight: ", round(free, 1), " GiB free.")
    }

    need(
      free >= min_gib,
      paste(
        "Only", round(free, 1),
        "GiB free; require at least",
        min_gib, "GiB."
      )
    )

    invisible(free)
  }


  fixed_table <- function(beta, V) {

    se <- sqrt(diag(V))
    z <- beta / se

    data.frame(
      term = names(beta),
      estimate = as.numeric(beta),
      std.error = as.numeric(se),
      conf.low = as.numeric(beta - 1.96 * se),
      conf.high = as.numeric(beta + 1.96 * se),
      z = as.numeric(z),
      p_Wald = as.numeric(2 * stats::pnorm(-abs(z))),
      stringsAsFactors = FALSE
    )
  }


  # ---------------------------------------------------------------------------
  # Current accepted bridge fits
  # ---------------------------------------------------------------------------

  current_fits <- function(bridge_run) {

    bridge_run <- normalizePath(
      bridge_run,
      winslash = "/",
      mustWork = TRUE
    )

    cfg_file <- file.path(
      bridge_run,
      "run_config.rds"
    )

    need(
      file.exists(cfg_file),
      "Bridge run_config.rds is missing."
    )

    bridge_cfg <- readRDS(cfg_file)

    need(
      identical(
        bridge_cfg$version,
        SOURCE_BRIDGE_VERSION
      ),
      "Unexpected source bridge workflow version."
    )

    prepared_multi <- file.path(
      bridge_run,
      "prepared_bridge_multivoltine.rds"
    )

    prepared_uni <- file.path(
      bridge_run,
      "prepared_bridge_univoltine.rds"
    )

    need(
      all(file.exists(c(prepared_multi, prepared_uni))),
      "One or both prepared bridge inputs are missing."
    )

    p_multi <- readRDS(prepared_multi)
    p_uni <- readRDS(prepared_uni)

    multi_fit <- file.path(
      bridge_run,
      "multivoltine__bridge",
      "fit.rds"
    )

    need(
      file.exists(multi_fit),
      "Accepted multivoltine bridge fit is missing."
    )

    multi <- readRDS(multi_fit)

    need(
      identical(multi$version, SOURCE_BRIDGE_VERSION) &&
        identical(multi$signature, p_multi$signature) &&
        isTRUE(multi$checks$numerical_checks_ok),
      "Multivoltine source bridge fit is not the accepted FIT_OK model."
    )

    refinement_root <- file.path(
      bridge_run,
      "refinement_univoltine_bridge"
    )

    need(
      dir.exists(refinement_root),
      "Refined/adjudicated univoltine bridge directory is missing."
    )

    dirs <- list.dirs(
      refinement_root,
      recursive = FALSE,
      full.names = TRUE
    )

    good <- dirs[vapply(dirs, function(d) {

      result_file <- file.path(d, "result.rds")
      fit_file <- file.path(
        d,
        "univoltine__bridge",
        "fit.rds"
      )

      if (!all(file.exists(c(result_file, fit_file)))) {
        return(FALSE)
      }

      result <- tryCatch(
        readRDS(result_file),
        error = function(e) NULL
      )

      fit <- tryCatch(
        readRDS(fit_file),
        error = function(e) NULL
      )

      !is.null(result) &&
        !is.null(fit) &&
        identical(result$status, "FIT_OK") &&
        identical(fit$version, REFINED_UNIVOLTINE_VERSION) &&
        identical(
          fit$source_bridge_signature,
          p_uni$signature
        ) &&
        isTRUE(fit$checks$numerical_checks_ok)

    }, logical(1))]

    need(
      length(good) == 1L,
      paste(
        "Expected exactly one accepted refined univoltine bridge fit; found",
        length(good)
      )
    )

    uni_fit <- file.path(
      good,
      "univoltine__bridge",
      "fit.rds"
    )

    list(
      bridge_run = bridge_run,
      bridge_config = bridge_cfg,
      source_run = normalizePath(
        bridge_cfg$source_run,
        winslash = "/",
        mustWork = TRUE
      ),
      prepared = list(
        multivoltine = normalizePath(
          prepared_multi,
          winslash = "/",
          mustWork = TRUE
        ),
        univoltine = normalizePath(
          prepared_uni,
          winslash = "/",
          mustWork = TRUE
        )
      ),
      fits = list(
        multivoltine = normalizePath(
          multi_fit,
          winslash = "/",
          mustWork = TRUE
        ),
        univoltine = normalizePath(
          uni_fit,
          winslash = "/",
          mustWork = TRUE
        )
      )
    )
  }


  # ---------------------------------------------------------------------------
  # Existing compiled likelihood engine
  # ---------------------------------------------------------------------------

  load_dll <- function(source_run) {

    dll <- TMB::dynlib(
      file.path(
        source_run,
        "compiled",
        "trend_spde"
      )
    )

    need(
      file.exists(dll),
      paste("Missing source DLL:", dll)
    )

    if (!"trend_spde" %in% names(getLoadedDLLs())) {
      dyn.load(dll)
    }

    TMB::openmp(
      n = 1L,
      DLL = "trend_spde"
    )

    invisible(dll)
  }


  sparse_TMB <- function(x) {

    methods::as(
      methods::as(
        methods::as(x, "dMatrix"),
        "generalMatrix"
      ),
      "TsparseMatrix"
    )
  }


  make_objective <- function(data, parameters) {

    for (nm in c("A", "M0", "M1", "M2")) {
      data[[nm]] <- sparse_TMB(data[[nm]])
    }

    TMB::MakeADFun(
      data = data,
      parameters = parameters,
      random = c(
        "b_species",
        "b_site",
        "ar_obs",
        "field"
      ),
      DLL = "trend_spde",
      silent = TRUE,
      inner.control = list(
        maxit = 1000,
        trace = FALSE
      )
    )
  }


  # ---------------------------------------------------------------------------
  # Descriptive VIF: population-level predictor collinearity
  # ---------------------------------------------------------------------------

  predictor_vif <- function(pop, vars) {

    out <- lapply(vars, function(v) {

      others <- setdiff(vars, v)

      f <- stats::as.formula(
        paste(
          v,
          "~",
          paste(others, collapse = " + ")
        )
      )

      mod <- stats::lm(f, data = pop)

      r2 <- summary(mod)$r.squared

      data.frame(
        predictor = v,
        R2_from_other_predictors = r2,
        VIF = 1 / (1 - r2),
        stringsAsFactors = FALSE
      )
    })

    do.call(rbind, out)
  }


  # ---------------------------------------------------------------------------
  # Linear contrasts and joint Wald tests
  # ---------------------------------------------------------------------------

  linear_contrast <- function(
    beta,
    V,
    weights,
    label,
    scale,
    interpretation
  ) {

    need(
      all(names(weights) %in% names(beta)),
      paste("Missing coefficient for contrast:", label)
    )

    cvec <- rep(0, length(beta))
    names(cvec) <- names(beta)
    cvec[names(weights)] <- weights

    estimate <- sum(cvec * beta)

    variance <- as.numeric(
      crossprod(
        cvec,
        V %*% cvec
      )
    )

    need(
      is.finite(variance) &&
        variance >= -1e-10,
      paste("Invalid contrast variance:", label)
    )

    se <- sqrt(max(0, variance))

    data.frame(
      contrast = label,
      scale = scale,
      estimate = estimate,
      std.error = se,
      conf.low = estimate - 1.96 * se,
      conf.high = estimate + 1.96 * se,
      p_Wald = 2 * stats::pnorm(
        -abs(estimate / se)
      ),
      ratio_of_decadal_factors = exp(estimate),
      ratio_conf.low = exp(estimate - 1.96 * se),
      ratio_conf.high = exp(estimate + 1.96 * se),
      percent_difference = 100 * expm1(estimate),
      percent_low = 100 * expm1(estimate - 1.96 * se),
      percent_high = 100 * expm1(estimate + 1.96 * se),
      interpretation = interpretation,
      stringsAsFactors = FALSE
    )
  }


  joint_wald <- function(
    beta,
    V,
    contrast_matrix,
    label,
    interpretation
  ) {

    C <- as.matrix(contrast_matrix)

    need(
      ncol(C) == length(beta),
      paste("Contrast matrix has wrong width:", label)
    )

    need(
      all(is.finite(C)),
      paste("Contrast matrix is non-finite:", label)
    )

    q <- as.numeric(C %*% beta)
    VC <- C %*% V %*% t(C)

    need(
      all(is.finite(VC)),
      paste("Joint covariance is non-finite:", label)
    )

    rank <- qr(VC)$rank

    need(
      rank == nrow(C),
      paste("Joint contrast is rank deficient:", label)
    )

    statistic <- as.numeric(
      crossprod(
        q,
        solve(VC, q)
      )
    )

    data.frame(
      test = label,
      Wald_chisq = statistic,
      df = nrow(C),
      p_Wald_joint = stats::pchisq(
        statistic,
        df = nrow(C),
        lower.tail = FALSE
      ),
      interpretation = interpretation,
      stringsAsFactors = FALSE
    )
  }


  contrast_row <- function(beta_names, weights) {

    x <- rep(0, length(beta_names))
    names(x) <- beta_names

    need(
      all(names(weights) %in% beta_names),
      paste(
        "Contrast references missing coefficient(s):",
        paste(
          setdiff(names(weights), beta_names),
          collapse = ", "
        )
      )
    )

    x[names(weights)] <- weights
    x
  }


  # ---------------------------------------------------------------------------
  # Prepare ONE within-between model
  # ---------------------------------------------------------------------------

  prepare_group <- function(cfg, group) {

    phase(
      cfg,
      "PREPARING",
      "Constructing exact within-between plasticity decomposition",
      group
    )

    source_prepared_file <- cfg$prepared[[group]]
    source_fit_file <- cfg$current_fits[[group]]

    need(
      file.exists(source_prepared_file) &&
        file.exists(source_fit_file),
      paste(group, "source bridge inputs are missing.")
    )

    source_prepared <- readRDS(source_prepared_file)
    source_fit <- readRDS(source_fit_file)

    need(
      identical(
        source_prepared$version,
        SOURCE_BRIDGE_VERSION
      ),
      paste(group, "unexpected prepared bridge version.")
    )

    source_bridge_signature <- if (
      !is.null(source_fit$source_bridge_signature)
    ) {
      source_fit$source_bridge_signature
    } else {
      source_fit$signature
    }

    need(
      identical(
        source_bridge_signature,
        source_prepared$signature
      ),
      paste(group, "source fit and prepared bridge signatures disagree.")
    )

    need(
      isTRUE(source_fit$checks$numerical_checks_ok),
      paste(group, "source bridge fit is not FIT_OK.")
    )

    d <- as.data.frame(source_prepared$frame)
    pop <- as.data.frame(source_prepared$populations)
    sites <- as.data.frame(source_prepared$sites)

    required_frame <- c(
      "SPECIES",
      "SITE_ID",
      "pop_id",
      "YEAR",
      "year_decade",
      "ABUND_INDEX",
      "latitude_within_species_10deg",
      "longitude_within_species_10deg",
      "onset_plasticity_z",
      "offset_plasticity_bio_z",
      "onset_plasticity_within_species_z",
      "onset_plasticity_species_mean_z",
      "offset_plasticity_within_species_z",
      "offset_plasticity_species_mean_z"
    )

    required_pop <- setdiff(
      required_frame,
      c("YEAR", "year_decade", "ABUND_INDEX")
    )

    need(
      all(required_frame %in% names(d)),
      paste(group, "prepared annual frame lacks within-between variables.")
    )

    need(
      all(required_pop %in% names(pop)),
      paste(group, "population table lacks within-between variables.")
    )

    d$SPECIES <- as.character(d$SPECIES)
    d$SITE_ID <- as.character(d$SITE_ID)
    d$pop_id <- as.character(d$pop_id)

    pop$SPECIES <- as.character(pop$SPECIES)
    pop$SITE_ID <- as.character(pop$SITE_ID)
    pop$pop_id <- as.character(pop$pop_id)

    need(
      !anyDuplicated(pop$pop_id),
      paste(group, "population table has duplicate pop_id.")
    )

    # Exact decomposition identity in the original group-wide z-score units.
    onset_identity <- max(abs(
      pop$onset_plasticity_z -
        (
          pop$onset_plasticity_within_species_z +
            pop$onset_plasticity_species_mean_z
        )
    ))

    offset_identity <- max(abs(
      pop$offset_plasticity_bio_z -
        (
          pop$offset_plasticity_within_species_z +
            pop$offset_plasticity_species_mean_z
        )
    ))

    need(
      is.finite(onset_identity) &&
        onset_identity < 1e-10,
      paste(group, "onset within-between identity failed.")
    )

    need(
      is.finite(offset_identity) &&
        offset_identity < 1e-10,
      paste(group, "offset within-between identity failed.")
    )

    # Species-centred components must sum to zero with equal population weight.
    need(
      max(abs(rowsum(
        pop$onset_plasticity_within_species_z,
        pop$SPECIES,
        reorder = FALSE
      ))) < 1e-8,
      paste(group, "onset within-species centring failed.")
    )

    need(
      max(abs(rowsum(
        pop$offset_plasticity_within_species_z,
        pop$SPECIES,
        reorder = FALSE
      ))) < 1e-8,
      paste(group, "offset within-species centring failed.")
    )

    # -----------------------------------------------------------------------
    # Within-between fixed-effect design
    # -----------------------------------------------------------------------

    X <- stats::model.matrix(
      ~ year_decade * (
        latitude_within_species_10deg +
        longitude_within_species_10deg +
        onset_plasticity_within_species_z +
        onset_plasticity_species_mean_z +
        offset_plasticity_within_species_z +
        offset_plasticity_species_mean_z
      ),
      data = d
    )

    expected <- c(
      "(Intercept)",
      "year_decade",
      "latitude_within_species_10deg",
      "longitude_within_species_10deg",
      "onset_plasticity_within_species_z",
      "onset_plasticity_species_mean_z",
      "offset_plasticity_within_species_z",
      "offset_plasticity_species_mean_z",
      "year_decade:latitude_within_species_10deg",
      "year_decade:longitude_within_species_10deg",
      "year_decade:onset_plasticity_within_species_z",
      "year_decade:onset_plasticity_species_mean_z",
      "year_decade:offset_plasticity_within_species_z",
      "year_decade:offset_plasticity_species_mean_z"
    )

    need(
      identical(colnames(X), expected),
      paste(
        group,
        "unexpected within-between fixed-effect design:\n",
        paste(colnames(X), collapse = "\n")
      )
    )

    need(
      qr(X)$rank == ncol(X),
      paste(group, "within-between design is rank deficient.")
    )

    data <- source_prepared$data
    data$X <- unname(X)

    need(
      length(data$y) == nrow(d),
      paste(group, "response/frame length mismatch.")
    )

    # -----------------------------------------------------------------------
    # Exact nested warm start from the accepted bridge
    #
    # At beta_within = beta_between = beta_original, the new fixed linear
    # predictor equals the accepted bridge predictor exactly because
    # P = P_within + P_between.
    # -----------------------------------------------------------------------

    source_beta <- source_fit$beta

    required_source_terms <- c(
      "(Intercept)",
      "year_decade",
      "latitude_within_species_10deg",
      "longitude_within_species_10deg",
      "onset_plasticity_z",
      "offset_plasticity_bio_z",
      "year_decade:latitude_within_species_10deg",
      "year_decade:longitude_within_species_10deg",
      "year_decade:onset_plasticity_z",
      "year_decade:offset_plasticity_bio_z"
    )

    need(
      all(required_source_terms %in% names(source_beta)),
      paste(group, "accepted bridge fixed effects have unexpected names.")
    )

    beta0 <- setNames(
      rep(0, ncol(X)),
      colnames(X)
    )

    directly_shared <- c(
      "(Intercept)",
      "year_decade",
      "latitude_within_species_10deg",
      "longitude_within_species_10deg",
      "year_decade:latitude_within_species_10deg",
      "year_decade:longitude_within_species_10deg"
    )

    beta0[directly_shared] <- source_beta[directly_shared]

    beta0[c(
      "onset_plasticity_within_species_z",
      "onset_plasticity_species_mean_z"
    )] <- source_beta["onset_plasticity_z"]

    beta0[c(
      "offset_plasticity_within_species_z",
      "offset_plasticity_species_mean_z"
    )] <- source_beta["offset_plasticity_bio_z"]

    beta0[c(
      "year_decade:onset_plasticity_within_species_z",
      "year_decade:onset_plasticity_species_mean_z"
    )] <- source_beta["year_decade:onset_plasticity_z"]

    beta0[c(
      "year_decade:offset_plasticity_within_species_z",
      "year_decade:offset_plasticity_species_mean_z"
    )] <- source_beta["year_decade:offset_plasticity_bio_z"]

    parameters <- source_fit$parameters
    parameters$beta <- unname(beta0)

    need(
      nrow(parameters$b_species) ==
        length(unique(d$SPECIES)),
      paste(group, "species random-effect dimension mismatch.")
    )

    need(
      nrow(parameters$b_site) ==
        nrow(sites),
      paste(group, "site random-effect dimension mismatch.")
    )

    need(
      length(parameters$ar_obs) ==
        nrow(d),
      paste(group, "AR1 latent dimension mismatch.")
    )

    # -----------------------------------------------------------------------
    # Predictor support, component scales and collinearity
    # -----------------------------------------------------------------------

    species_one <- pop[
      !duplicated(pop$SPECIES),
      ,
      drop = FALSE
    ]

    component_scales <- data.frame(
      group = group,
      component = c(
        "onset_within_species",
        "onset_between_species",
        "offset_within_species",
        "offset_between_species"
      ),
      fitted_variable = c(
        "onset_plasticity_within_species_z",
        "onset_plasticity_species_mean_z",
        "offset_plasticity_within_species_z",
        "offset_plasticity_species_mean_z"
      ),
      SD_for_descriptive_rescaling = c(
        stats::sd(pop$onset_plasticity_within_species_z),
        stats::sd(species_one$onset_plasticity_species_mean_z),
        stats::sd(pop$offset_plasticity_within_species_z),
        stats::sd(species_one$offset_plasticity_species_mean_z)
      ),
      SD_weighting = c(
        "equal population",
        "equal species",
        "equal population",
        "equal species"
      ),
      fitted_unit = paste0(
        "one original group-wide plasticity SD; ",
        "no component re-standardization in the model"
      ),
      stringsAsFactors = FALSE
    )

    need(
      all(is.finite(
        component_scales$SD_for_descriptive_rescaling
      )) &&
        all(
          component_scales$SD_for_descriptive_rescaling > 0
        ),
      paste(group, "one or more within-between components are constant.")
    )

    predictors <- c(
      "latitude_within_species_10deg",
      "longitude_within_species_10deg",
      "onset_plasticity_within_species_z",
      "onset_plasticity_species_mean_z",
      "offset_plasticity_within_species_z",
      "offset_plasticity_species_mean_z"
    )

    cor_pop <- stats::cor(
      pop[predictors],
      use = "complete.obs"
    )

    cor_out <- data.frame(
      predictor = rownames(cor_pop),
      cor_pop,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )

    vif_out <- predictor_vif(
      pop,
      predictors
    )

    # Variance decomposition of each original plasticity variable.
    variance_decomposition <- data.frame(
      group = group,
      plasticity = c("onset", "offset"),
      total_variance = c(
        stats::var(pop$onset_plasticity_z),
        stats::var(pop$offset_plasticity_bio_z)
      ),
      within_population_variance = c(
        stats::var(pop$onset_plasticity_within_species_z),
        stats::var(pop$offset_plasticity_within_species_z)
      ),
      between_population_weighted_variance = c(
        stats::var(pop$onset_plasticity_species_mean_z),
        stats::var(pop$offset_plasticity_species_mean_z)
      ),
      stringsAsFactors = FALSE
    )

    variance_decomposition$proportion_variance_within <-
      variance_decomposition$within_population_variance /
      variance_decomposition$total_variance

    variance_decomposition$proportion_variance_between <-
      variance_decomposition$between_population_weighted_variance /
      variance_decomposition$total_variance

    # -----------------------------------------------------------------------
    # Signature and exports
    # -----------------------------------------------------------------------

    signature <- digest::digest(
      list(
        version = VERSION,
        group = group,
        source_bridge_signature = source_prepared$signature,
        source_fit_path = normalizePath(
          source_fit_file,
          winslash = "/",
          mustWork = TRUE
        ),
        source_fit_objective =
          source_fit$optimizer$objective,
        source_fit_size =
          file.info(source_fit_file)$size,
        source_fit_mtime =
          as.character(file.info(source_fit_file)$mtime),
        X = X,
        response = data$y
      ),
      algo = "sha256"
    )

    prepared <- list(
      version = VERSION,
      group = group,
      signature = signature,
      source_bridge_signature =
        source_prepared$signature,
      source_data_signature =
        source_prepared$source_data_signature,
      frame = d,
      populations = pop,
      sites = sites,
      years = source_prepared$years,
      X = X,
      data = data,
      parameters = parameters,
      source_prepared_file =
        normalizePath(
          source_prepared_file,
          winslash = "/",
          mustWork = TRUE
        ),
      source_fit_file =
        normalizePath(
          source_fit_file,
          winslash = "/",
          mustWork = TRUE
        ),
      source_fit_version = source_fit$version,
      source_fit_beta = source_fit$beta,
      source_fit_V = source_fit$V,
      source_fit_checks = source_fit$checks,
      source_fit_objective =
        source_fit$optimizer$objective,
      source_fit_eta = source_fit$eta,
      component_scales = component_scales
    )

    prepared_file <- file.path(
      cfg$out,
      paste0(
        "prepared_mundlak_",
        group,
        ".rds"
      )
    )

    atomic_save(
      prepared,
      prepared_file
    )

    csv(
      component_scales,
      file.path(
        cfg$out,
        paste0(
          "component_scales_",
          group,
          ".csv"
        )
      )
    )

    csv(
      variance_decomposition,
      file.path(
        cfg$out,
        paste0(
          "plasticity_variance_decomposition_",
          group,
          ".csv"
        )
      )
    )

    csv(
      cor_out,
      file.path(
        cfg$out,
        paste0(
          "predictor_correlations_",
          group,
          ".csv"
        )
      )
    )

    csv(
      vif_out,
      file.path(
        cfg$out,
        paste0(
          "predictor_VIF_",
          group,
          ".csv"
        )
      )
    )

    data.frame(
      group = group,
      n_observations = nrow(d),
      n_populations = nrow(pop),
      n_species = length(unique(d$SPECIES)),
      n_sites = nrow(sites),
      min_year = min(d$YEAR),
      max_year = max(d$YEAR),
      onset_identity_max_abs_error = onset_identity,
      offset_identity_max_abs_error = offset_identity,
      maximum_population_level_VIF =
        max(vif_out$VIF),
      source_fit = normalizePath(
        source_fit_file,
        winslash = "/",
        mustWork = TRUE
      ),
      mundlak_signature = signature,
      stringsAsFactors = FALSE
    )
  }


  # ---------------------------------------------------------------------------
  # Fit ONE within-between model
  # ---------------------------------------------------------------------------

  # ---------------------------------------------------------------------------
  # Fit ONE bridge model
  # ---------------------------------------------------------------------------

  fit_group <- function(cfg, group) {

    load_dll(cfg$source_run)

    prepared_file <- file.path(
      cfg$out,
      paste0(
        "prepared_mundlak_",
        group,
        ".rds"
      )
    )

    need(
      file.exists(prepared_file),
      paste(
        "Prepared within-between input missing for",
        group
      )
    )

    p <- readRDS(prepared_file)

    out_group <- file.path(
      cfg$out,
      paste0(
        group,
        "__mundlak"
      )
    )

    dir.create(
      out_group,
      recursive = TRUE,
      showWarnings = FALSE
    )

    fit_file <- file.path(
      out_group,
      "fit.rds"
    )

    if (file.exists(fit_file)) {

      fit <- readRDS(fit_file)

      need(
        identical(fit$version, VERSION) &&
          identical(fit$signature, p$signature),
        paste(group, "cached within-between fit mismatch.")
      )

      phase(
        cfg,
        "REUSING",
        "Completed within-between fit found",
        group
      )

      return(
        if (isTRUE(fit$checks$numerical_checks_ok)) {
          "FIT_OK_CACHED"
        } else {
          "NUMERICAL_REVIEW_REQUIRED"
        }
      )
    }

    data <- p$data

    # -----------------------------------------------------------------------
    # Verify that the constrained starting point exactly reproduces the
    # accepted bridge model before fitting the expanded within-between model.
    # -----------------------------------------------------------------------

    replay_file <- file.path(
      out_group,
      "nested_start_replay.csv"
    )

    if (!file.exists(replay_file)) {

      phase(
        cfg,
        "VERIFYING_NESTED_START",
        paste0(
          "Checking beta_within = beta_between reproduces ",
          "the accepted bridge likelihood"
        ),
        group
      )

      replay_obj <- make_objective(
        data,
        p$parameters
      )

      replay_value <- replay_obj$fn(
        replay_obj$par
      )

      replay_report <- replay_obj$report(
        replay_obj$env$last.par
      )

      eta_error <- max(abs(
        as.numeric(replay_report$eta) -
          as.numeric(p$source_fit_eta)
      ))

      objective_error <- abs(
        replay_value -
          p$source_fit_objective
      )

      replay <- data.frame(
        group = group,
        source_bridge_objective =
          p$source_fit_objective,
        mapped_mundlak_objective =
          replay_value,
        absolute_objective_error =
          objective_error,
        eta_max_abs_error =
          eta_error,
        equivalent =
          is.finite(replay_value) &&
          objective_error < 0.001 &&
          is.finite(eta_error) &&
          eta_error < 1e-5,
        stringsAsFactors = FALSE
      )

      csv(
        replay,
        replay_file
      )

      try(
        TMB::FreeADFun(replay_obj),
        silent = TRUE
      )

      rm(
        replay_obj,
        replay_report
      )

      invisible(gc())

      need(
        isTRUE(replay$equivalent),
        paste(
          group,
          "within-between constrained start did not reproduce the accepted bridge."
        )
      )

    } else {

      replay <- utils::read.csv(
        replay_file,
        stringsAsFactors = FALSE
      )

      need(
        nrow(replay) == 1L &&
          isTRUE(replay$equivalent),
        paste(group, "saved nested-start replay is invalid.")
      )
    }

    pars <- p$parameters

    warm_file <- file.path(
      out_group,
      "last_evaluation.rds"
    )

    opt_file <- file.path(
      out_group,
      "optimized_parameters.rds"
    )

    warm <- if (file.exists(warm_file)) {

      z <- readRDS(warm_file)

      if (
        identical(z$version, VERSION) &&
        identical(z$signature, p$signature)
      ) {
        z
      } else {
        NULL
      }

    } else {
      NULL
    }

    optimized <- if (file.exists(opt_file)) {

      z <- readRDS(opt_file)

      if (
        identical(z$version, VERSION) &&
        identical(z$signature, p$signature)
      ) {
        z
      } else {
        NULL
      }

    } else {
      NULL
    }

    if (!is.null(warm)) {
      pars <- warm$parameters
    }

    if (!is.null(optimized)) {
      pars <- optimized$parameters
    }

    obj <- make_objective(
      data,
      pars
    )

    on.exit(
      try(
        TMB::FreeADFun(obj),
        silent = TRUE
      ),
      add = TRUE
    )

    if (
      !is.null(warm) &&
      is.null(optimized)
    ) {
      obj$par <- warm$par
    }

    if (!is.null(optimized)) {
      obj$par <- optimized$optimizer$par
    }

    warnings <- character()

    capture_warning <- function(w) {

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

      invokeRestart("muffleWarning")
    }

    current_parameters <- function(par) {

      obj$fn(par)

      x <- obj$env$parList(
        par = obj$env$last.par
      )

      need(
        all(vapply(
          x,
          function(v) {
            is.numeric(v) &&
              all(is.finite(v))
          },
          logical(1)
        )),
        paste(group, "non-finite parameter arrays.")
      )

      x
    }

    if (is.null(optimized)) {

      n_fn <- 0L
      n_gr <- 0L
      last_save <- Sys.time() - 3600

      objective <- function(par) {

        val <- obj$fn(par)

        n_fn <<- n_fn + 1L

        if (
          is.finite(val) &&
          (
            n_fn == 1L ||
            as.numeric(
              difftime(
                Sys.time(),
                last_save,
                units = "secs"
              )
            ) >= 90
          )
        ) {

          pars_now <- current_parameters(par)

          atomic_save(
            list(
              version = VERSION,
              signature = p$signature,
              par = par,
              parameters = pars_now,
              objective = val,
              time = Sys.time()
            ),
            warm_file
          )

          phase(
            cfg,
            "FITTING",
            paste(
              "fn", n_fn,
              "| gr", n_gr,
              "| objective", signif(val, 10),
              "| checkpoint saved"
            ),
            group
          )

          last_save <<- Sys.time()
        }

        val
      }

      gradient <- function(par) {
        n_gr <<- n_gr + 1L
        obj$gr(par)
      }

      phase(
        cfg,
        "FITTING",
        paste(
          nrow(p$frame),
          "annual rows;",
          nrow(p$populations),
          "populations;",
          length(unique(p$frame$SPECIES)),
          "species;",
          ncol(p$X),
          "fixed coefficients"
        ),
        group
      )

      set.seed(
        if (group == "multivoltine") {
          202610131L
        } else {
          202610132L
        }
      )

      opt <- withCallingHandlers(
        stats::nlminb(
          obj$par,
          objective,
          gradient,
          control = list(
            iter.max = cfg$iter,
            eval.max = cfg$eval,
            rel.tol = 1e-10
          )
        ),
        warning = capture_warning
      )

      if (opt$convergence != 0L) {

        phase(
          cfg,
          "REFINING",
          paste(
            "nlminb code",
            opt$convergence,
            "-> one BFGS refinement"
          ),
          group
        )

        second <- withCallingHandlers(
          stats::optim(
            opt$par,
            objective,
            gradient,
            method = "BFGS",
            control = list(
              maxit = 300L,
              reltol = 1e-10
            )
          ),
          warning = capture_warning
        )

        if (
          is.finite(second$value) &&
          second$value <=
            opt$objective + 1e-6
        ) {

          opt <- list(
            par = second$par,
            objective = second$value,
            convergence = second$convergence,
            message = "BFGS refinement",
            first_optimizer = opt,
            counts = second$counts
          )
        }
      }

      final_value <- obj$fn(
        opt$par
      )

      need(
        is.finite(final_value) &&
          abs(
            final_value -
              opt$objective
          ) < 0.01,
        paste(group, "final objective mismatch.")
      )

      atomic_save(
        list(
          version = VERSION,
          signature = p$signature,
          parameters = current_parameters(
            opt$par
          ),
          optimizer = opt
        ),
        opt_file
      )

    } else {

      phase(
        cfg,
        "RECOVERING_UNCERTAINTY",
        "Optimizer already completed; no repeated optimization",
        group
      )

      opt <- optimized$optimizer

      value <- obj$fn(
        opt$par
      )

      need(
        is.finite(value) &&
          abs(
            value -
              opt$objective
          ) < 0.01,
        paste(group, "saved optimum replay failed.")
      )
    }

    # -----------------------------------------------------------------------
    # Hessian / covariance
    # -----------------------------------------------------------------------

    phase(
      cfg,
      "UNCERTAINTY",
      "Computing Hessian/covariance",
      group
    )

    sd <- withCallingHandlers(
      TMB::sdreport(
        obj,
        par.fixed = opt$par,
        getJointPrecision = FALSE,
        getReportCovariance = FALSE
      ),
      warning = capture_warning
    )

    gradient <- as.numeric(
      obj$gr(opt$par)
    )

    V_full <- sd$cov.fixed

    beta_ids <- which(
      names(opt$par) == "beta"
    )

    need(
      length(beta_ids) == ncol(p$X),
      paste(group, "beta mapping changed.")
    )

    beta <- setNames(
      as.numeric(
        opt$par[beta_ids]
      ),
      colnames(p$X)
    )

    V <- V_full[
      beta_ids,
      beta_ids,
      drop = FALSE
    ]

    dimnames(V) <- list(
      names(beta),
      names(beta)
    )

    scaled_gradient <- if (
      all(is.finite(gradient)) &&
      is.matrix(V_full) &&
      all(is.finite(V_full))
    ) {

      sqrt(max(
        0,
        as.numeric(
          crossprod(
            gradient,
            V_full %*% gradient
          )
        )
      ))

    } else {
      Inf
    }

    obj$fn(opt$par)

    parameters <- current_parameters(
      opt$par
    )

    report <- obj$report(
      obj$env$last.par
    )

    extent <- max(
      diff(range(p$sites$x_km)),
      diff(range(p$sites$y_km))
    )

    good <-
      opt$convergence == 0L &&
      isTRUE(sd$pdHess) &&
      is.matrix(V) &&
      all(is.finite(V)) &&
      all(diag(V) > 0) &&
      is.finite(scaled_gradient) &&
      scaled_gradient < GRAD_TOL &&
      is.finite(opt$objective) &&
      opt$objective <= p$source_fit_objective + 0.001

    checks <- data.frame(
      group = group,
      convergence_code =
        opt$convergence,
      pdHess = isTRUE(sd$pdHess),
      max_abs_gradient =
        max(abs(gradient)),
      gradient_covariance_norm =
        scaled_gradient,
      numerical_checks_ok = good,
      logLik = -opt$objective,
      parameter_df = length(opt$par),
      AIC =
        2 * opt$objective +
        2 * length(opt$par),
      objective_change_from_constrained_bridge =
        opt$objective - p$source_fit_objective,
      nobs = nrow(p$frame),
      spatial_range_km =
        as.numeric(report$range),
      spatial_sd =
        as.numeric(report$spatial_sd),
      range_greater_than_1_5_extent =
        as.numeric(report$range) >
        1.5 * extent,
      spatial_sd_below_0_001 =
        as.numeric(report$spatial_sd) <
        0.001,
      AR1_rho =
        as.numeric(report$rho),
      Gamma_shape =
        as.numeric(report$shape),
      stringsAsFactors = FALSE
    )

    fit <- list(
      version = VERSION,
      signature = p$signature,
      source_bridge_signature =
        p$source_bridge_signature,
      source_data_signature =
        p$source_data_signature,
      parameters = parameters,
      optimizer = opt,
      beta = beta,
      V = V,
      cov_fixed = V_full,
      checks = checks,
      eta = as.numeric(report$eta),
      field_at_sites = report$spatial,
      group = group,
      warnings = unique(warnings)
    )

    atomic_save(
      fit,
      fit_file
    )

    csv(
      checks,
      file.path(
        out_group,
        "fit_checks.csv"
      )
    )

    writeLines(
      unique(warnings),
      file.path(
        out_group,
        "warnings.txt"
      )
    )

    fx <- fixed_table(
      beta,
      V
    )

    csv(
      fx,
      file.path(
        out_group,
        "fixed_effects.csv"
      )
    )

    csv(
      data.frame(
        term = rownames(V),
        V,
        check.names = FALSE
      ),
      file.path(
        out_group,
        "fixed_effect_covariance.csv"
      )
    )

    # -----------------------------------------------------------------------
    # Focal temporal effects
    # -----------------------------------------------------------------------

    lat_term <-
      "year_decade:latitude_within_species_10deg"

    lon_term <-
      "year_decade:longitude_within_species_10deg"

    onset_within_term <-
      "year_decade:onset_plasticity_within_species_z"

    onset_between_term <-
      "year_decade:onset_plasticity_species_mean_z"

    offset_within_term <-
      "year_decade:offset_plasticity_within_species_z"

    offset_between_term <-
      "year_decade:offset_plasticity_species_mean_z"

    key_terms <- c(
      lat_term,
      lon_term,
      onset_within_term,
      onset_between_term,
      offset_within_term,
      offset_between_term
    )

    key <- fx[
      match(key_terms, fx$term),
      ,
      drop = FALSE
    ]

    need(
      !anyNA(key$term),
      paste(group, "one or more focal temporal terms are missing.")
    )

    key$group <- group

    key$effect_scale <- c(
      "within species",
      "within species",
      "within species",
      "between species",
      "within species",
      "between species"
    )

    key$unit <- c(
      "+10 degrees latitude",
      "+10 degrees longitude",
      "+1 original group-wide SD in population onset deviation",
      "+1 original group-wide SD in species mean onset plasticity",
      if (group == "univoltine") {
        "+1 original group-wide SD in population offset-advancement deviation"
      } else {
        "+1 original group-wide SD in population offset-delay deviation"
      },
      if (group == "univoltine") {
        "+1 original group-wide SD in species mean offset advancement"
      } else {
        "+1 original group-wide SD in species mean offset delay"
      }
    )

    key$ratio_of_decadal_growth_factors <-
      exp(key$estimate)

    key$ratio_conf.low <-
      exp(key$conf.low)

    key$ratio_conf.high <-
      exp(key$conf.high)

    key$percent_difference_in_decadal_growth_factor <-
      100 * expm1(key$estimate)

    key$percent_low <-
      100 * expm1(key$conf.low)

    key$percent_high <-
      100 * expm1(key$conf.high)

    key <- key[c(
      "group",
      "effect_scale",
      "term",
      "unit",
      "estimate",
      "std.error",
      "conf.low",
      "conf.high",
      "p_Wald",
      "ratio_of_decadal_growth_factors",
      "ratio_conf.low",
      "ratio_conf.high",
      "percent_difference_in_decadal_growth_factor",
      "percent_low",
      "percent_high"
    )]

    csv(
      key,
      file.path(
        out_group,
        "key_within_between_effects.csv"
      )
    )

    # -----------------------------------------------------------------------
    # Effects rescaled to one observed component SD.
    # The fitted model itself remains in the original group-wide z units.
    # -----------------------------------------------------------------------

    component_term_map <- data.frame(
      fitted_variable = c(
        "onset_plasticity_within_species_z",
        "onset_plasticity_species_mean_z",
        "offset_plasticity_within_species_z",
        "offset_plasticity_species_mean_z"
      ),
      term = c(
        onset_within_term,
        onset_between_term,
        offset_within_term,
        offset_between_term
      ),
      stringsAsFactors = FALSE
    )

    scaled <- merge(
      p$component_scales,
      component_term_map,
      by = "fitted_variable",
      sort = FALSE
    )

    scaled <- merge(
      scaled,
      fx,
      by = "term",
      sort = FALSE
    )

    scaled$estimate_per_component_SD <-
      scaled$estimate *
      scaled$SD_for_descriptive_rescaling

    scaled$SE_per_component_SD <-
      scaled$std.error *
      scaled$SD_for_descriptive_rescaling

    scaled$low_per_component_SD <-
      scaled$conf.low *
      scaled$SD_for_descriptive_rescaling

    scaled$high_per_component_SD <-
      scaled$conf.high *
      scaled$SD_for_descriptive_rescaling

    scaled$ratio_per_component_SD <-
      exp(scaled$estimate_per_component_SD)

    scaled$percent_per_component_SD <-
      100 * expm1(
        scaled$estimate_per_component_SD
      )

    csv(
      scaled,
      file.path(
        out_group,
        "plasticity_effects_per_component_SD.csv"
      )
    )

    # -----------------------------------------------------------------------
    # Within-versus-between contrasts
    # -----------------------------------------------------------------------

    differences <- rbind(
      linear_contrast(
        beta,
        V,
        weights = c(
          setNames(1, onset_within_term),
          setNames(-1, onset_between_term)
        ),
        label = "onset_within_minus_between",
        scale = "contextual difference",
        interpretation = paste0(
          "Difference between the within-species onset-plasticity association ",
          "and the between-species onset-plasticity association."
        )
      ),
      linear_contrast(
        beta,
        V,
        weights = c(
          setNames(1, offset_within_term),
          setNames(-1, offset_between_term)
        ),
        label = "offset_within_minus_between",
        scale = "contextual difference",
        interpretation = paste0(
          "Difference between the within-species offset-plasticity association ",
          "and the between-species offset-plasticity association."
        )
      )
    )

    differences$group <- group

    differences <- differences[c(
      "group",
      setdiff(
        names(differences),
        "group"
      )
    )]

    csv(
      differences,
      file.path(
        out_group,
        "within_between_differences.csv"
      )
    )

    # -----------------------------------------------------------------------
    # Joint Wald tests
    # -----------------------------------------------------------------------

    beta_names <- names(beta)

    C_within <- rbind(
      contrast_row(
        beta_names,
        setNames(1, onset_within_term)
      ),
      contrast_row(
        beta_names,
        setNames(1, offset_within_term)
      )
    )

    C_between <- rbind(
      contrast_row(
        beta_names,
        setNames(1, onset_between_term)
      ),
      contrast_row(
        beta_names,
        setNames(1, offset_between_term)
      )
    )

    C_contextual <- rbind(
      contrast_row(
        beta_names,
        c(
          setNames(1, onset_within_term),
          setNames(-1, onset_between_term)
        )
      ),
      contrast_row(
        beta_names,
        c(
          setNames(1, offset_within_term),
          setNames(-1, offset_between_term)
        )
      )
    )

    C_geography <- rbind(
      contrast_row(
        beta_names,
        setNames(1, lat_term)
      ),
      contrast_row(
        beta_names,
        setNames(1, lon_term)
      )
    )

    joint_tests <- rbind(
      joint_wald(
        beta,
        V,
        C_within,
        label = "within_species_plasticity_temporal_effects",
        interpretation = paste0(
          "Joint 2-df test that the within-species onset and offset ",
          "plasticity associations with abundance trend are both zero."
        )
      ),
      joint_wald(
        beta,
        V,
        C_between,
        label = "between_species_plasticity_temporal_effects",
        interpretation = paste0(
          "Joint 2-df test that the between-species onset and offset ",
          "plasticity associations with abundance trend are both zero."
        )
      ),
      joint_wald(
        beta,
        V,
        C_contextual,
        label = "within_equals_between_for_onset_and_offset",
        interpretation = paste0(
          "Joint 2-df test of equality between within- and between-species ",
          "temporal plasticity associations."
        )
      ),
      joint_wald(
        beta,
        V,
        C_geography,
        label = "temporal_geography",
        interpretation = paste0(
          "Joint 2-df test of year:latitude and year:longitude after ",
          "within-between plasticity adjustment."
        )
      )
    )

    joint_tests$group <- group

    joint_tests <- joint_tests[c(
      "group",
      setdiff(
        names(joint_tests),
        "group"
      )
    )]

    csv(
      joint_tests,
      file.path(
        out_group,
        "joint_Wald_tests.csv"
      )
    )

    # -----------------------------------------------------------------------
    # Onset-offset balance separately within and between species
    # -----------------------------------------------------------------------

    if (group == "univoltine") {

      within_coordinated_weights <- c(
        setNames(1, onset_within_term),
        setNames(1, offset_within_term)
      )

      within_extension_weights <- c(
        setNames(1, onset_within_term),
        setNames(-1, offset_within_term)
      )

      between_coordinated_weights <- c(
        setNames(1, onset_between_term),
        setNames(1, offset_between_term)
      )

      between_extension_weights <- c(
        setNames(1, onset_between_term),
        setNames(-1, offset_between_term)
      )

    } else {

      # Multivoltine offset_plasticity_bio_z is positive towards delay.
      # Therefore offset advancement = -offset_bio.
      within_coordinated_weights <- c(
        setNames(1, onset_within_term),
        setNames(-1, offset_within_term)
      )

      within_extension_weights <- c(
        setNames(1, onset_within_term),
        setNames(1, offset_within_term)
      )

      between_coordinated_weights <- c(
        setNames(1, onset_between_term),
        setNames(-1, offset_between_term)
      )

      between_extension_weights <- c(
        setNames(1, onset_between_term),
        setNames(1, offset_between_term)
      )
    }

    balance <- rbind(
      linear_contrast(
        beta,
        V,
        within_coordinated_weights,
        label = "coordinated_advancement",
        scale = "within species",
        interpretation = paste0(
          "Within species: onset and offset shift towards coordinated ",
          "whole-season advancement, controlling the extension axis and geography."
        )
      ),
      linear_contrast(
        beta,
        V,
        within_extension_weights,
        label = "season_extension",
        scale = "within species",
        interpretation = paste0(
          "Within species: stronger onset advancement combined with relatively ",
          "later termination, controlling coordinated advancement and geography."
        )
      ),
      linear_contrast(
        beta,
        V,
        between_coordinated_weights,
        label = "coordinated_advancement",
        scale = "between species",
        interpretation = paste0(
          "Between species: species-mean coordinated whole-season advancement, ",
          "controlling the extension axis and geography."
        )
      ),
      linear_contrast(
        beta,
        V,
        between_extension_weights,
        label = "season_extension",
        scale = "between species",
        interpretation = paste0(
          "Between species: species-mean onset advancement combined with ",
          "relatively later termination, controlling coordinated advancement ",
          "and geography."
        )
      )
    )

    balance$group <- group

    balance <- balance[c(
      "group",
      setdiff(
        names(balance),
        "group"
      )
    )]

    csv(
      balance,
      file.path(
        out_group,
        "onset_offset_balance_within_between.csv"
      )
    )

    # -----------------------------------------------------------------------
    # Compare the original bridge coefficient with within and between effects
    # -----------------------------------------------------------------------

    source_fx <- fixed_table(
      p$source_fit_beta,
      p$source_fit_V
    )

    source_terms <- c(
      onset =
        "year_decade:onset_plasticity_z",
      offset =
        "year_decade:offset_plasticity_bio_z"
    )

    comparison <- rbind(
      data.frame(
        group = group,
        plasticity = "onset",
        source_bridge_estimate =
          source_fx$estimate[
            match(
              source_terms["onset"],
              source_fx$term
            )
          ],
        source_bridge_low =
          source_fx$conf.low[
            match(
              source_terms["onset"],
              source_fx$term
            )
          ],
        source_bridge_high =
          source_fx$conf.high[
            match(
              source_terms["onset"],
              source_fx$term
            )
          ],
        within_estimate =
          beta[onset_within_term],
        within_low =
          fx$conf.low[
            match(onset_within_term, fx$term)
          ],
        within_high =
          fx$conf.high[
            match(onset_within_term, fx$term)
          ],
        between_estimate =
          beta[onset_between_term],
        between_low =
          fx$conf.low[
            match(onset_between_term, fx$term)
          ],
        between_high =
          fx$conf.high[
            match(onset_between_term, fx$term)
          ],
        stringsAsFactors = FALSE
      ),
      data.frame(
        group = group,
        plasticity = "offset",
        source_bridge_estimate =
          source_fx$estimate[
            match(
              source_terms["offset"],
              source_fx$term
            )
          ],
        source_bridge_low =
          source_fx$conf.low[
            match(
              source_terms["offset"],
              source_fx$term
            )
          ],
        source_bridge_high =
          source_fx$conf.high[
            match(
              source_terms["offset"],
              source_fx$term
            )
          ],
        within_estimate =
          beta[offset_within_term],
        within_low =
          fx$conf.low[
            match(offset_within_term, fx$term)
          ],
        within_high =
          fx$conf.high[
            match(offset_within_term, fx$term)
          ],
        between_estimate =
          beta[offset_between_term],
        between_low =
          fx$conf.low[
            match(offset_between_term, fx$term)
          ],
        between_high =
          fx$conf.high[
            match(offset_between_term, fx$term)
          ],
        stringsAsFactors = FALSE
      )
    )

    csv(
      comparison,
      file.path(
        out_group,
        "source_bridge_vs_within_between.csv"
      )
    )

    # -----------------------------------------------------------------------
    # Nested likelihood comparison against the constrained bridge
    # -----------------------------------------------------------------------

    source_checks <- p$source_fit_checks

    need(
      nrow(source_checks) == 1L,
      paste(group, "unexpected source bridge check table.")
    )

    LR <- 2 * (
      checks$logLik -
      source_checks$logLik
    )

    df_added <-
      checks$parameter_df -
      source_checks$parameter_df

    nested <- data.frame(
      group = group,
      constrained_model =
        "geography_plus_overall_plasticity_bridge",
      expanded_model =
        "geography_plus_within_between_plasticity",
      logLik_constrained =
        source_checks$logLik,
      logLik_expanded =
        checks$logLik,
      AIC_constrained =
        source_checks$AIC,
      AIC_expanded =
        checks$AIC,
      delta_AIC_expanded_minus_constrained =
        checks$AIC -
        source_checks$AIC,
      LR_statistic =
        max(0, LR),
      added_parameters =
        df_added,
      p_LRT = if (
        is.finite(LR) &&
        LR >= -1e-5 &&
        df_added > 0
      ) {
        stats::pchisq(
          max(0, LR),
          df = df_added,
          lower.tail = FALSE
        )
      } else {
        NA_real_
      },
      interpretation = paste0(
        "The constrained bridge imposes equal within- and between-species ",
        "coefficients for onset and offset main effects and year interactions."
      ),
      stringsAsFactors = FALSE
    )

    csv(
      nested,
      file.path(
        out_group,
        "nested_within_between_LRT.csv"
      )
    )

    phase(
      cfg,
      if (good) {
        "FIT_OK"
      } else {
        "NUMERICAL_REVIEW_REQUIRED"
      },
      paste(
        "onset within",
        signif(beta[onset_within_term], 5),
        "| onset between",
        signif(beta[onset_between_term], 5),
        "| offset within",
        signif(beta[offset_within_term], 5),
        "| offset between",
        signif(beta[offset_between_term], 5),
        "| scaled gradient",
        signif(scaled_gradient, 5)
      ),
      group
    )

    if (good) {
      "FIT_OK"
    } else {
      "NUMERICAL_REVIEW_REQUIRED"
    }
  }


  # ---------------------------------------------------------------------------
  # Combine outputs
  # ---------------------------------------------------------------------------

  # ---------------------------------------------------------------------------
  # Combine outputs
  # ---------------------------------------------------------------------------

  combine_outputs <- function(cfg, statuses) {

    specs <- list(
      key = c(
        file = "key_within_between_effects.csv",
        combined = "MUNDLAK_key_effects.csv"
      ),
      scaled = c(
        file = "plasticity_effects_per_component_SD.csv",
        combined = "MUNDLAK_effects_per_component_SD.csv"
      ),
      differences = c(
        file = "within_between_differences.csv",
        combined = "MUNDLAK_within_between_differences.csv"
      ),
      joint = c(
        file = "joint_Wald_tests.csv",
        combined = "MUNDLAK_joint_Wald_tests.csv"
      ),
      balance = c(
        file = "onset_offset_balance_within_between.csv",
        combined = "MUNDLAK_onset_offset_balance.csv"
      ),
      source_comparison = c(
        file = "source_bridge_vs_within_between.csv",
        combined = "MUNDLAK_source_bridge_comparison.csv"
      ),
      nested = c(
        file = "nested_within_between_LRT.csv",
        combined = "MUNDLAK_nested_LRT.csv"
      ),
      checks = c(
        file = "fit_checks.csv",
        combined = "MUNDLAK_fit_checks.csv"
      )
    )

    for (nm in names(specs)) {

      rows <- list()

      for (group in GROUPS) {

        path <- file.path(
          cfg$out,
          paste0(group, "__mundlak"),
          specs[[nm]][["file"]]
        )

        if (file.exists(path)) {
          rows[[group]] <- utils::read.csv(
            path,
            stringsAsFactors = FALSE,
            check.names = FALSE
          )
        }
      }

      if (length(rows)) {
        csv(
          do.call(rbind, rows),
          file.path(
            cfg$out,
            specs[[nm]][["combined"]]
          )
        )
      }
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
  # Lock helpers
  # ---------------------------------------------------------------------------

  lock_owner_alive <- function(lock) {

    owner_file <- file.path(
      lock,
      "owner.rds"
    )

    owner <- tryCatch(
      readRDS(owner_file),
      error = function(e) NULL
    )

    if (
      is.null(owner) ||
      is.null(owner$pid) ||
      is.null(owner$create_time)
    ) {
      return(FALSE)
    }

    tryCatch(
      ps::ps_is_running(
        ps::ps_handle(
          owner$pid,
          time = owner$create_time
        )
      ),
      error = function(e) FALSE
    )
  }


  release_lock <- function(cfg) {

    if (
      !is.null(cfg$lock) &&
      dir.exists(cfg$lock)
    ) {

      owner <- tryCatch(
        readRDS(
          file.path(
            cfg$lock,
            "owner.rds"
          )
        ),
        error = function(e) NULL
      )

      if (
        !is.null(owner) &&
        identical(
          owner$token,
          cfg$lock_token
        )
      ) {
        unlink(
          cfg$lock,
          recursive = TRUE
        )
      }
    }

    invisible(NULL)
  }


  # ---------------------------------------------------------------------------
  # Worker
  # ---------------------------------------------------------------------------

  worker <- function(cfg) {

    check_packages()
    disk_check(cfg$out, cfg$min_disk)
    load_dll(cfg$source_run)

    tryCatch({

      phase(
        cfg,
        "PREPARING",
        "Preparing multivoltine and univoltine within-between models"
      )

      summaries <- lapply(
        GROUPS,
        function(group) {
          prepare_group(
            cfg,
            group
          )
        }
      )

      csv(
        do.call(rbind, summaries),
        file.path(
          cfg$out,
          "mundlak_data_summary.csv"
        )
      )

      statuses <- data.frame(
        group = character(),
        status = character(),
        stringsAsFactors = FALSE
      )

      for (group in GROUPS) {

        disk_check(
          cfg$out,
          min_gib = 10,
          verbose = FALSE
        )

        status <- fit_group(
          cfg,
          group
        )

        statuses <- rbind(
          statuses,
          data.frame(
            group = group,
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

        invisible(gc())
      }

      combine_outputs(
        cfg,
        statuses
      )

      writeLines(
        c(
          "phenoIMPACT | WITHIN-BETWEEN PLASTICITY-DEMOGRAPHY BRIDGE",
          "",
          paste("Source bridge run:", cfg$bridge_run),
          paste("Likelihood/mesh source:", cfg$source_run),
          "",
          "Question:",
          "Does phenological plasticity predict abundance trends within species,",
          "between species, or both, after controlling within-species latitude",
          "and longitude?",
          "",
          "Model:",
          "ABUND_INDEX ~ year_decade * (latitude_within_species_10deg +",
          "longitude_within_species_10deg +",
          "onset_plasticity_within_species_z +",
          "onset_plasticity_species_mean_z +",
          "offset_plasticity_within_species_z +",
          "offset_plasticity_species_mean_z),",
          "with the same Gamma/log likelihood, species/site temporal slopes,",
          "population AR1 and annual IID Matern field as the accepted bridge.",
          "",
          "Plasticity units:",
          "All plasticity components retain the original group-wide z-score",
          "scale. No component is re-standardized in the fitted model.",
          "Therefore P = P_within + P_species_mean exactly, and the accepted",
          "bridge model is nested under beta_within = beta_between.",
          "",
          "Primary interpretation:",
          "- year:onset_within and year:offset_within test population differences",
          "  within the same species.",
          "- year:onset_species_mean and year:offset_species_mean test differences",
          "  between species in their mean plasticity.",
          "- within_between_differences tests whether the two scales differ.",
          "",
          "This is an association model, not a causal mediation analysis.",
          "",
          "Main combined outputs:",
          "MUNDLAK_key_effects.csv",
          "MUNDLAK_effects_per_component_SD.csv",
          "MUNDLAK_within_between_differences.csv",
          "MUNDLAK_joint_Wald_tests.csv",
          "MUNDLAK_onset_offset_balance.csv",
          "MUNDLAK_source_bridge_comparison.csv",
          "MUNDLAK_nested_LRT.csv",
          "MUNDLAK_fit_checks.csv"
        ),
        file.path(
          cfg$out,
          "README_mundlak.txt"
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

      zip_path <- file.path(
        cfg$out,
        "mundlak_results_to_review.zip"
      )

      if (requireNamespace("zip", quietly = TRUE)) {

        zip::zipr(
          zip_path,
          files = review_files,
          root = cfg$out,
          mode = "mirror",
          include_directories = FALSE
        )

      } else {

        old <- getwd()
        on.exit(setwd(old), add = TRUE)

        setwd(cfg$out)

        utils::zip(
          zipfile = zip_path,
          files = review_files
        )
      }

      terminal_status <- if (
        all(grepl("^FIT_OK", statuses$status))
      ) {
        "COMPLETED"
      } else {
        "COMPLETED_WITH_NUMERICAL_REVIEW"
      }

      result <- list(
        version = VERSION,
        signature = cfg$run_signature,
        status = terminal_status,
        output = cfg$out,
        zip = zip_path,
        finished = Sys.time()
      )

      atomic_save(
        result,
        file.path(
          cfg$out,
          "result.rds"
        )
      )

      phase(
        cfg,
        terminal_status,
        zip_path,
        terminal = TRUE
      )

      result

    }, error = function(e) {

      writeLines(
        conditionMessage(e),
        file.path(
          cfg$out,
          "fatal_error.txt"
        )
      )

      try(
        phase(
          cfg,
          "STOPPED_WITH_ERROR",
          conditionMessage(e),
          terminal = TRUE
        ),
        silent = TRUE
      )

      list(
        status = "STOPPED_WITH_ERROR",
        error = conditionMessage(e),
        output = cfg$out
      )
    })
  }


  worker_entry <- function(cfg) {

    on.exit(
      release_lock(cfg),
      add = TRUE
    )

    worker(cfg)
  }


  # ---------------------------------------------------------------------------
  # Dump self-contained worker engine
  # ---------------------------------------------------------------------------

  dump_engine <- function(out) {

    e <- environment(dump_engine)

    n <- ls(
      e,
      all.names = TRUE
    )

    fun <- n[
      vapply(
        n,
        function(k) {
          is.function(
            get(k, envir = e)
          )
        },
        logical(1)
      )
    ]

    path <- file.path(
      out,
      "mundlak_workflow_engine.R"
    )

    dump(
      c(
        "VERSION",
        "SOURCE_BRIDGE_VERSION",
        "REFINED_UNIVOLTINE_VERSION",
        "ROOT",
        "BRIDGE_RUN",
        "GROUPS",
        "GRAD_TOL",
        fun
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
  # Public runner
  # ---------------------------------------------------------------------------

  run <- function(
    project_root = ROOT,
    bridge_run = BRIDGE_RUN,
    monitor = FALSE,
    minimum_free_GiB = 25,
    optimizer_iterations = 1200L,
    optimizer_evaluations = 1800L
  ) {

    check_packages()

    current <- state$job

    if (
      !is.null(current) &&
      !is.null(current$process) &&
      current$process$is_alive()
    ) {
      message("Within-between workflow already running.")
      return(invisible(current))
    }

    project_root <- normalizePath(
      project_root,
      winslash = "/",
      mustWork = TRUE
    )

    sources <- current_fits(
      bridge_run
    )

    source_paths <- c(
      file.path(
        sources$bridge_run,
        "run_config.rds"
      ),
      unlist(
        sources$prepared,
        use.names = FALSE
      ),
      unlist(
        sources$fits,
        use.names = FALSE
      ),
      TMB::dynlib(
        file.path(
          sources$source_run,
          "compiled",
          "trend_spde"
        )
      ),
      file.path(
        sources$source_run,
        "mesh.rds"
      )
    )

    need(
      all(file.exists(source_paths)),
      "One or more deterministic-run source files are missing."
    )

    source_metadata <- data.frame(
      path = normalizePath(
        source_paths,
        winslash = "/",
        mustWork = TRUE
      ),
      bytes = file.info(source_paths)$size,
      modified = as.character(
        file.info(source_paths)$mtime
      ),
      stringsAsFactors = FALSE
    )

    run_signature <- digest::digest(
      list(
        version = VERSION,
        source_metadata = source_metadata,
        optimizer_iterations =
          as.integer(optimizer_iterations),
        optimizer_evaluations =
          as.integer(optimizer_evaluations),
        gradient_threshold = GRAD_TOL
      ),
      algo = "sha256"
    )

    out_root <- file.path(
      project_root,
      "output",
      "population_trends",
      "mundlak_geography_plasticity_spatial_15km"
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
        substr(run_signature, 1, 12)
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

    result_file <- file.path(
      out,
      "result.rds"
    )

    if (file.exists(result_file)) {

      result <- readRDS(result_file)

      need(
        identical(
          result$signature,
          run_signature
        ),
        "Completed deterministic run has a different signature."
      )

      message(
        "Matching within-between run is already complete: ",
        result$status
      )

      job <- list(
        process = NULL,
        out = out,
        log = NULL,
        completed = TRUE
      )

      state$job <- job

      options(
        phenoimpact.mundlak.job = job
      )

      if (monitor) {
        status(job)
      }

      return(invisible(job))
    }

    lock <- file.path(
      out_root,
      paste0(
        "queue_",
        substr(run_signature, 1, 12),
        ".lock"
      )
    )

    if (dir.exists(lock)) {

      alive <- lock_owner_alive(lock)

      need(
        !alive,
        paste(
          "A matching within-between workflow appears to be running:",
          lock
        )
      )

      unlink(
        lock,
        recursive = TRUE
      )
    }

    need(
      dir.create(
        lock,
        recursive = TRUE
      ),
      "Cannot acquire within-between workflow lock."
    )

    lock_token <- paste(
      Sys.getpid(),
      format(
        Sys.time(),
        "%Y%m%d%H%M%OS6"
      ),
      sep = "_"
    )

    parent_handle <- ps::ps_handle()

    atomic_save(
      list(
        pid = Sys.getpid(),
        create_time =
          ps::ps_create_time(parent_handle),
        token = lock_token
      ),
      file.path(
        lock,
        "owner.rds"
      )
    )

    cfg <- list(
      version = VERSION,
      run_signature = run_signature,
      project_root = project_root,
      bridge_run = sources$bridge_run,
      source_run = sources$source_run,
      prepared = sources$prepared,
      current_fits = sources$fits,
      source_metadata = source_metadata,
      out = out,
      progress = file.path(
        out,
        "progress.rds"
      ),
      lock = lock,
      lock_token = lock_token,
      min_disk = minimum_free_GiB,
      iter = as.integer(
        optimizer_iterations
      ),
      eval = as.integer(
        optimizer_evaluations
      )
    )

    cfg$engine <- dump_engine(out)

    atomic_save(
      cfg,
      file.path(
        out,
        "run_config.rds"
      )
    )

    csv(
      source_metadata,
      file.path(
        out,
        "source_manifest.csv"
      )
    )

    phase(
      cfg,
      "QUEUED",
      paste(
        "Two within-between fits:",
        paste(
          GROUPS,
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

    child <- NULL
    success <- FALSE

    on.exit(
      if (!success) {

        if (
          !is.null(child) &&
          child$is_alive()
        ) {
          child$kill()
        }

        release_lock(cfg)
      },
      add = TRUE
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

        e$worker_entry(cfg)
      },
      args = list(
        engine = cfg$engine,
        cfg = cfg
      ),
      libpath = .libPaths(),
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

    child_handle <- ps::ps_handle(
      child$get_pid()
    )

    atomic_save(
      list(
        pid = child$get_pid(),
        create_time =
          ps::ps_create_time(child_handle),
        token = lock_token
      ),
      file.path(
        lock,
        "owner.rds"
      )
    )

    job <- list(
      process = child,
      out = out,
      log = log_file,
      started = Sys.time(),
      completed = FALSE
    )

    state$job <- job

    options(
      phenoimpact.mundlak.job = job
    )

    success <- TRUE

    message(
      "Within-between workflow launched. Keep R open."
    )

    message(
      "Deterministic output/checkpoint folder: ",
      out
    )

    message(
      "Monitor: pheno_mundlak$status()"
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
        "phenoimpact.mundlak.job"
      )
    }

    need(
      !is.null(job),
      "No registered within-between job."
    )

    progress_file <- file.path(
      job$out,
      "progress.rds"
    )

    if (file.exists(progress_file)) {

      z <- tryCatch(
        readRDS(progress_file),
        error = function(e) NULL
      )

      if (!is.null(z)) {
        print(z)
      }
    }

    alive <- if (
      !is.null(job$process)
    ) {
      job$process$is_alive()
    } else {
      FALSE
    }

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

      cat("\nCompleted fits:\n")

      print(
        utils::read.csv(
          partial,
          stringsAsFactors = FALSE
        ),
        row.names = FALSE
      )
    }

    if (
      !is.null(job$process) &&
      !alive
    ) {

      cat(
        "Exit status:",
        job$process$get_exit_status(),
        "\n"
      )
    }

    key <- file.path(
      job$out,
      "MUNDLAK_key_effects.csv"
    )

    if (
      file.exists(key) &&
      !alive
    ) {

      cat(
        "\nWithin-between temporal effects:\n"
      )

      print(
        utils::read.csv(
          key,
          stringsAsFactors = FALSE,
          check.names = FALSE
        ),
        row.names = FALSE
      )
    }

    differences <- file.path(
      job$out,
      "MUNDLAK_within_between_differences.csv"
    )

    if (
      file.exists(differences) &&
      !alive
    ) {

      cat(
        "\nWithin-minus-between contrasts:\n"
      )

      print(
        utils::read.csv(
          differences,
          stringsAsFactors = FALSE,
          check.names = FALSE
        ),
        row.names = FALSE
      )
    }

    joint <- file.path(
      job$out,
      "MUNDLAK_joint_Wald_tests.csv"
    )

    if (
      file.exists(joint) &&
      !alive
    ) {

      cat(
        "\nJoint Wald tests:\n"
      )

      print(
        utils::read.csv(
          joint,
          stringsAsFactors = FALSE,
          check.names = FALSE
        ),
        row.names = FALSE
      )
    }

    balance <- file.path(
      job$out,
      "MUNDLAK_onset_offset_balance.csv"
    )

    if (
      file.exists(balance) &&
      !alive
    ) {

      cat(
        "\nOnset-offset balance within and between species:\n"
      )

      print(
        utils::read.csv(
          balance,
          stringsAsFactors = FALSE,
          check.names = FALSE
        ),
        row.names = FALSE
      )
    }

    zip_file <- file.path(
      job$out,
      "mundlak_results_to_review.zip"
    )

    if (
      file.exists(zip_file) &&
      !alive
    ) {

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

    invisible(job)
  }


  list(
    run = run,
    status = status
  )
})


run_mundlak_geography_plasticity <- function(...) {
  pheno_mundlak$run(...)
}


status_mundlak_geography_plasticity <- function(...) {
  pheno_mundlak$status(...)
}
