# =============================================================================
# phenoIMPACT | BRIDGE: geography + phenological plasticity -> abundance trends
# Version 1.0 | 2026-10-08
#
# SCIENTIFIC QUESTION
#   Does the association between phenological plasticity and population trend
#   persist after explicitly controlling the within-species latitude/longitude
#   gradient? And does latitude retain an association after plasticity enters?
#
# TWO MODELS, fitted sequentially:
#   1) multivoltine
#   2) univoltine
#
# MODEL
#   ABUND_INDEX ~ year_decade * (
#       latitude_within_species_10deg +
#       longitude_within_species_10deg +
#       onset_plasticity_z +
#       offset_plasticity_bio_z
#   )
#   + species intercept + year slope
#   + site intercept + year slope
#   + population AR(1)
#   + annual IID Matern spatial field
#
# IMPORTANT
#   - Same annual rows as the FINAL onset-adjusted plasticity-abundance models.
#   - Same Gamma/log likelihood, random effects, AR(1), mesh and annual field.
#   - Same onset_plasticity_z and offset_plasticity_bio_z used in those models.
#   - Latitude/longitude are centred within species using equal population
#     weights WITHIN EACH fitted voltinism subset.
#   - Current fitted spatial models are WARM STARTS only. All parameters are
#     re-estimated.
#   - This bridge asks whether the EXISTING plasticity-demography association
#     survives geographic adjustment. It is not a causal mediation analysis.
#   - Plasticity is not yet decomposed into within- vs between-species effects.
#     Species-centred plasticity is exported for a possible stricter sensitivity.
#
# RUN
#   source(file.choose(), encoding = "UTF-8")
#   bridge_job <- pheno_bridge$run(monitor = FALSE)
#
# MONITOR
#   pheno_bridge$status()
#
# Keep R open. Only ONE large fit runs at a time.
# =============================================================================


pheno_bridge <- local({

  VERSION <- "bridge_geography_plasticity_spatial_15km_1.0"

  ROOT <- "E:/phenoIMPACT project/code/phenoIMPACT"

  SOURCE_RUN <- paste0(
    ROOT,
    "/output/population_trends/spatial_abundance_onset_adjusted_15km/",
    "run_20261004_214012_44136"
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
  # Current accepted source fits
  # ---------------------------------------------------------------------------

  current_fits <- function(source_run) {

    multi <- file.path(
      source_run,
      "multivoltine__full",
      "fit.rds"
    )

    need(
      file.exists(multi),
      "Current multivoltine fit missing."
    )

    refinement_root <- file.path(
      source_run,
      "refinement_univoltine"
    )

    need(
      dir.exists(refinement_root),
      "Univoltine refinement directory missing."
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
        "univoltine__full",
        "fit.rds"
      )

      if (!all(file.exists(c(result_file, fit_file)))) {
        return(FALSE)
      }

      z <- tryCatch(
        readRDS(result_file),
        error = function(e) NULL
      )

      !is.null(z) && identical(z$status, "FIT_OK")

    }, logical(1))]

    need(
      length(good) == 1L,
      paste(
        "Expected exactly one FIT_OK univoltine refinement; found",
        length(good)
      )
    )

    list(
      multivoltine = normalizePath(
        multi,
        winslash = "/",
        mustWork = TRUE
      ),
      univoltine = normalizePath(
        file.path(
          good,
          "univoltine__full",
          "fit.rds"
        ),
        winslash = "/",
        mustWork = TRUE
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
  # Prepare ONE group
  # ---------------------------------------------------------------------------

  prepare_group <- function(cfg, group) {

    phase(
      cfg,
      "PREPARING",
      "Reading final spatial abundance input and accepted fit",
      group
    )

    input_file <- file.path(
      cfg$source_run,
      paste0("input_", group, ".rds")
    )

    need(
      file.exists(input_file),
      paste("Missing source input:", input_file)
    )

    input <- readRDS(input_file)
    source_fit <- readRDS(cfg$current_fits[[group]])

    need(
      isTRUE(source_fit$checks$numerical_checks_ok),
      paste(group, "source fit requires numerical review.")
    )

    need(
      identical(
        source_fit$data_signature,
        input$data_signature
      ),
      paste(
        group,
        "accepted source fit and source input have different data signatures."
      )
    )

    need(
      all(c(
        "SPECIES",
        "SITE_ID",
        "pop_id",
        "YEAR",
        "year_decade",
        "ABUND_INDEX",
        "onset_plasticity_z",
        "offset_plasticity_bio_z"
      ) %in% names(input$frame)),
      paste(group, "source frame has unexpected columns.")
    )

    d <- input$frame

    d$SPECIES <- as.character(d$SPECIES)
    d$SITE_ID <- as.character(d$SITE_ID)
    d$pop_id <- as.character(d$pop_id)

    sites <- input$sites

    need(
      all(c(
        "SITE_ID", "x_km", "y_km"
      ) %in% names(sites)),
      paste(group, "source site table is incomplete.")
    )

    sites$SITE_ID <- as.character(sites$SITE_ID)

    need(
      !anyDuplicated(sites$SITE_ID),
      paste(group, "source site table has duplicate SITE_ID.")
    )

    # -----------------------------------------------------------------------
    # Convert exact source coordinates EPSG:3035 -> lon/lat
    # -----------------------------------------------------------------------

    site_sf <- sites
    site_sf$x_m <- site_sf$x_km * 1000
    site_sf$y_m <- site_sf$y_km * 1000

    site_sf <- sf::st_as_sf(
      site_sf,
      coords = c("x_m", "y_m"),
      crs = 3035
    )

    ll <- sf::st_coordinates(
      sf::st_transform(site_sf, 4326)
    )

    sites$longitude_deg <- ll[, 1]
    sites$latitude_deg <- ll[, 2]

    # -----------------------------------------------------------------------
    # Population-level table
    # -----------------------------------------------------------------------

    pop <- unique(
      d[c(
        "pop_id",
        "SPECIES",
        "SITE_ID",
        "onset_plasticity_z",
        "offset_plasticity_bio_z"
      )]
    )

    need(
      !anyDuplicated(pop$pop_id),
      paste(group, "population table has duplicate pop_id.")
    )

    m_site <- match(
      pop$SITE_ID,
      sites$SITE_ID
    )

    need(
      !anyNA(m_site),
      paste(group, "population-site coordinate join failed.")
    )

    pop$longitude_deg <- sites$longitude_deg[m_site]
    pop$latitude_deg <- sites$latitude_deg[m_site]

    # -----------------------------------------------------------------------
    # Geographic centring: equal POPULATION weight within species
    # -----------------------------------------------------------------------

    lat_mean <- tapply(
      pop$latitude_deg,
      pop$SPECIES,
      mean
    )

    lon_mean <- tapply(
      pop$longitude_deg,
      pop$SPECIES,
      mean
    )

    pop$latitude_species_mean <-
      unname(lat_mean[pop$SPECIES])

    pop$longitude_species_mean <-
      unname(lon_mean[pop$SPECIES])

    pop$latitude_within_species_deg <-
      pop$latitude_deg -
      pop$latitude_species_mean

    pop$longitude_within_species_deg <-
      pop$longitude_deg -
      pop$longitude_species_mean

    pop$latitude_within_species_10deg <-
      pop$latitude_within_species_deg / 10

    pop$longitude_within_species_10deg <-
      pop$longitude_within_species_deg / 10

    # -----------------------------------------------------------------------
    # Export species-centred plasticity for later sensitivity.
    # NOT used in this primary bridge fit.
    # -----------------------------------------------------------------------

    onset_sp_mean <- tapply(
      pop$onset_plasticity_z,
      pop$SPECIES,
      mean
    )

    offset_sp_mean <- tapply(
      pop$offset_plasticity_bio_z,
      pop$SPECIES,
      mean
    )

    pop$onset_plasticity_species_mean_z <-
      unname(onset_sp_mean[pop$SPECIES])

    pop$offset_plasticity_species_mean_z <-
      unname(offset_sp_mean[pop$SPECIES])

    pop$onset_plasticity_within_species_z <-
      pop$onset_plasticity_z -
      pop$onset_plasticity_species_mean_z

    pop$offset_plasticity_within_species_z <-
      pop$offset_plasticity_bio_z -
      pop$offset_plasticity_species_mean_z

    # Centring checks.
    lat_sum <- rowsum(
      pop$latitude_within_species_deg,
      pop$SPECIES,
      reorder = FALSE
    )

    lon_sum <- rowsum(
      pop$longitude_within_species_deg,
      pop$SPECIES,
      reorder = FALSE
    )

    need(
      max(abs(lat_sum)) < 1e-8,
      paste(group, "latitude centring failed.")
    )

    need(
      max(abs(lon_sum)) < 1e-8,
      paste(group, "longitude centring failed.")
    )

    # -----------------------------------------------------------------------
    # Join population geography back to exact annual rows
    # -----------------------------------------------------------------------

    m_pop <- match(
      d$pop_id,
      pop$pop_id
    )

    need(
      !anyNA(m_pop),
      paste(group, "population geography join failed.")
    )

    add_vars <- c(
      "latitude_deg",
      "longitude_deg",
      "latitude_species_mean",
      "longitude_species_mean",
      "latitude_within_species_deg",
      "longitude_within_species_deg",
      "latitude_within_species_10deg",
      "longitude_within_species_10deg",
      "onset_plasticity_species_mean_z",
      "offset_plasticity_species_mean_z",
      "onset_plasticity_within_species_z",
      "offset_plasticity_within_species_z"
    )

    for (v in add_vars) {
      d[[v]] <- pop[[v]][m_pop]
    }

    # -----------------------------------------------------------------------
    # BRIDGE fixed-effects design
    # -----------------------------------------------------------------------

    X <- stats::model.matrix(
      ~ year_decade * (
        latitude_within_species_10deg +
        longitude_within_species_10deg +
        onset_plasticity_z +
        offset_plasticity_bio_z
      ),
      data = d
    )

    expected <- c(
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
      identical(colnames(X), expected),
      paste(
        group,
        "unexpected bridge fixed-effects design:\n",
        paste(colnames(X), collapse = "\n")
      )
    )

    need(
      qr(X)$rank == ncol(X),
      paste(group, "bridge fixed-effects design is rank deficient.")
    )

    # -----------------------------------------------------------------------
    # Reuse exact source mesh
    # -----------------------------------------------------------------------

    mesh_file <- file.path(
      cfg$source_run,
      "mesh.rds"
    )

    need(
      file.exists(mesh_file),
      "Source mesh.rds is missing."
    )

    mesh <- readRDS(mesh_file)

    A <- methods::as(
      fmesher::fm_basis(
        mesh$mesh,
        loc = as.matrix(
          sites[c("x_km", "y_km")]
        )
      ),
      "CsparseMatrix"
    )

    need(
      all(abs(Matrix::rowSums(A) - 1) < 1e-8),
      paste(group, "some sites lie outside the source mesh.")
    )

    # -----------------------------------------------------------------------
    # TMB data: identical response/random structure, new X only
    # -----------------------------------------------------------------------

    data <- input$data

    data$X <- unname(X)
    data$A <- A
    data$M0 <- mesh$M0
    data$M1 <- mesh$M1
    data$M2 <- mesh$M2
    data$use_spatial <- 1L

    need(
      length(data$y) == nrow(d),
      paste(group, "response/frame length mismatch.")
    )

    # -----------------------------------------------------------------------
    # Warm start from accepted FINAL spatial plasticity model
    # -----------------------------------------------------------------------

    parameters <- source_fit$parameters

    source_beta <- source_fit$beta

    beta0 <- setNames(
      rep(0, ncol(X)),
      colnames(X)
    )

    shared <- intersect(
      names(source_beta),
      names(beta0)
    )

    beta0[shared] <- source_beta[shared]

    parameters$beta <- unname(beta0)

    # Random/hyperparameter dimensions must still match exactly.
    need(
      nrow(parameters$b_species) ==
        length(unique(d$SPECIES)),
      paste(group, "species RE dimension mismatch.")
    )

    need(
      nrow(parameters$b_site) ==
        nrow(sites),
      paste(group, "site RE dimension mismatch.")
    )

    need(
      length(parameters$ar_obs) == nrow(d),
      paste(group, "AR1 latent dimension mismatch.")
    )

    need(
      nrow(parameters$field) ==
        nrow(mesh$M0) &&
        ncol(parameters$field) ==
        length(input$years),
      paste(group, "spatial-field dimension mismatch.")
    )

    # -----------------------------------------------------------------------
    # Predictor correlation / VIF, equal population weight
    # -----------------------------------------------------------------------

    predictors <- c(
      "latitude_within_species_10deg",
      "longitude_within_species_10deg",
      "onset_plasticity_z",
      "offset_plasticity_bio_z"
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

    # Within-species plasticity correlations with geography.
    within_vars <- c(
      "latitude_within_species_10deg",
      "longitude_within_species_10deg",
      "onset_plasticity_within_species_z",
      "offset_plasticity_within_species_z"
    )

    within_cor <- stats::cor(
      pop[within_vars],
      use = "complete.obs"
    )

    within_cor_out <- data.frame(
      predictor = rownames(within_cor),
      within_cor,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )

    # -----------------------------------------------------------------------
    # Signature and export
    # -----------------------------------------------------------------------

    signature <- digest::digest(
      list(
        source_data_signature = input$data_signature,
        frame = d[c(
          "SPECIES",
          "SITE_ID",
          "pop_id",
          "YEAR",
          "ABUND_INDEX",
          "year_decade",
          "latitude_within_species_10deg",
          "longitude_within_species_10deg",
          "onset_plasticity_z",
          "offset_plasticity_bio_z"
        )],
        X = X
      ),
      algo = "sha256"
    )

    prepared <- list(
      version = VERSION,
      group = group,
      source_data_signature = input$data_signature,
      signature = signature,
      frame = d,
      populations = pop,
      sites = sites,
      years = input$years,
      X = X,
      data = data,
      parameters = parameters,
      source_fit_path = cfg$current_fits[[group]],
      source_fit_beta = source_fit$beta,
      source_fit_V = source_fit$V,
      source_fit_checks = source_fit$checks
    )

    prepared_file <- file.path(
      cfg$out,
      paste0("prepared_bridge_", group, ".rds")
    )

    atomic_save(
      prepared,
      prepared_file
    )

    csv(
      pop,
      file.path(
        cfg$out,
        paste0("population_predictors_", group, ".csv")
      )
    )

    csv(
      cor_out,
      file.path(
        cfg$out,
        paste0("predictor_correlations_", group, ".csv")
      )
    )

    csv(
      within_cor_out,
      file.path(
        cfg$out,
        paste0(
          "within_species_predictor_correlations_",
          group,
          ".csv"
        )
      )
    )

    csv(
      vif_out,
      file.path(
        cfg$out,
        paste0("predictor_VIF_", group, ".csv")
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
      source_fit = cfg$current_fits[[group]],
      source_data_signature = input$data_signature,
      bridge_signature = signature,
      stringsAsFactors = FALSE
    )
  }


  # ---------------------------------------------------------------------------
  # Fit ONE bridge model
  # ---------------------------------------------------------------------------

  fit_group <- function(cfg, group) {

    load_dll(cfg$source_run)

    prepared_file <- file.path(
      cfg$out,
      paste0("prepared_bridge_", group, ".rds")
    )

    need(
      file.exists(prepared_file),
      paste("Prepared bridge input missing for", group)
    )

    p <- readRDS(prepared_file)

    out_group <- file.path(
      cfg$out,
      paste0(group, "__bridge")
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
        paste(group, "cached bridge fit mismatch.")
      )

      phase(
        cfg,
        "REUSING",
        "Completed bridge fit found",
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
      ) z else NULL

    } else NULL

    optimized <- if (file.exists(opt_file)) {

      z <- readRDS(opt_file)

      if (
        identical(z$version, VERSION) &&
        identical(z$signature, p$signature)
      ) z else NULL

    } else NULL

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
      try(TMB::FreeADFun(obj), silent = TRUE),
      add = TRUE
    )

    if (!is.null(warm) && is.null(optimized)) {
      obj$par <- warm$par
    }

    if (!is.null(optimized)) {
      obj$par <- optimized$optimizer$par
    }

    warnings <- character()

    capture_warning <- function(w) {

      warnings <<- unique(
        c(warnings, conditionMessage(w))
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
          "species"
        ),
        group
      )

      set.seed(
        if (group == "multivoltine") {
          202610081L
        } else {
          202610082L
        }
      )

      opt <- withCallingHandlers(
        nlminb(
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
            "-> BFGS"
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
          second$value <= opt$objective + 1e-6
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

      val <- obj$fn(opt$par)

      need(
        is.finite(val) &&
          abs(val - opt$objective) < 0.01,
        paste(group, "final objective mismatch.")
      )

      atomic_save(
        list(
          version = VERSION,
          signature = p$signature,
          parameters = current_parameters(opt$par),
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

      val <- obj$fn(opt$par)

      need(
        is.finite(val) &&
          abs(val - opt$objective) < 0.01,
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

    g <- as.numeric(
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
      as.numeric(opt$par[beta_ids]),
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
      all(is.finite(g)) &&
      all(is.finite(V_full))
    ) {
      sqrt(
        max(
          0,
          as.numeric(
            crossprod(
              g,
              V_full %*% g
            )
          )
        )
      )
    } else Inf

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
      all(is.finite(V)) &&
      all(diag(V) > 0) &&
      is.finite(scaled_gradient) &&
      scaled_gradient < GRAD_TOL &&
      is.finite(opt$objective)

    checks <- data.frame(
      group = group,
      convergence_code = opt$convergence,
      pdHess = isTRUE(sd$pdHess),
      max_abs_gradient = max(abs(g)),
      gradient_covariance_norm = scaled_gradient,
      numerical_checks_ok = good,
      logLik = -opt$objective,
      parameter_df = length(opt$par),
      AIC = 2 * opt$objective + 2 * length(opt$par),
      nobs = nrow(p$frame),
      spatial_range_km = as.numeric(report$range),
      spatial_sd = as.numeric(report$spatial_sd),
      range_greater_than_1_5_extent =
        as.numeric(report$range) > 1.5 * extent,
      spatial_sd_below_0_001 =
        as.numeric(report$spatial_sd) < 0.001,
      AR1_rho = as.numeric(report$rho),
      Gamma_shape = as.numeric(report$shape),
      stringsAsFactors = FALSE
    )

    fit <- list(
      version = VERSION,
      signature = p$signature,
      source_data_signature = p$source_data_signature,
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
    # Key bridge effects
    # -----------------------------------------------------------------------

    key_terms <- c(
      "year_decade:latitude_within_species_10deg",
      "year_decade:longitude_within_species_10deg",
      "year_decade:onset_plasticity_z",
      "year_decade:offset_plasticity_bio_z"
    )

    key <- fx[
      match(key_terms, fx$term),
      ,
      drop = FALSE
    ]

    need(
      !anyNA(key$term),
      paste(group, "one or more key bridge interactions are missing.")
    )

    key$group <- group

    key$unit <- c(
      "+10 degrees latitude",
      "+10 degrees longitude",
      "+1 SD onset advancement plasticity",
      if (group == "univoltine") {
        "+1 SD stronger offset advancement"
      } else {
        "+1 SD stronger offset delay/extension"
      }
    )

    key$ratio_of_decadal_growth_factors <-
      exp(key$estimate)

    key$ratio_conf_low <-
      exp(key$conf.low)

    key$ratio_conf_high <-
      exp(key$conf.high)

    key$percent_difference_in_decadal_growth_factor <-
      100 * expm1(key$estimate)

    key$percent_low <-
      100 * expm1(key$conf.low)

    key$percent_high <-
      100 * expm1(key$conf.high)

    key <- key[
      c(
        "group",
        "term",
        "unit",
        "estimate",
        "std.error",
        "conf.low",
        "conf.high",
        "p_Wald",
        "ratio_of_decadal_growth_factors",
        "ratio_conf_low",
        "ratio_conf_high",
        "percent_difference_in_decadal_growth_factor",
        "percent_low",
        "percent_high"
      )
    ]

    csv(
      key,
      file.path(
        out_group,
        "key_bridge_effects.csv"
      )
    )

    # -----------------------------------------------------------------------
    # Compare plasticity effects to FINAL plasticity-only spatial model
    # -----------------------------------------------------------------------

    source_fx <- fixed_table(
      p$source_fit_beta,
      p$source_fit_V
    )

    compare_terms <- c(
      "year_decade:onset_plasticity_z",
      "year_decade:offset_plasticity_bio_z"
    )

    source_sub <- source_fx[
      match(compare_terms, source_fx$term),
      ,
      drop = FALSE
    ]

    bridge_sub <- fx[
      match(compare_terms, fx$term),
      ,
      drop = FALSE
    ]

    need(
      !anyNA(source_sub$term) &&
        !anyNA(bridge_sub$term),
      paste(group, "plasticity comparison terms are missing.")
    )

    comparison <- data.frame(
      group = group,
      term = compare_terms,
      plasticity_only_estimate = source_sub$estimate,
      plasticity_only_low = source_sub$conf.low,
      plasticity_only_high = source_sub$conf.high,
      plasticity_only_p_Wald = source_sub$p_Wald,
      bridge_estimate = bridge_sub$estimate,
      bridge_low = bridge_sub$conf.low,
      bridge_high = bridge_sub$conf.high,
      bridge_p_Wald = bridge_sub$p_Wald,
      absolute_change_bridge_minus_source =
        bridge_sub$estimate -
        source_sub$estimate,
      proportional_change_in_magnitude =
        ifelse(
          abs(source_sub$estimate) > 1e-12,
          abs(bridge_sub$estimate) /
            abs(source_sub$estimate),
          NA_real_
        ),
      same_sign =
        sign(bridge_sub$estimate) ==
        sign(source_sub$estimate),
      stringsAsFactors = FALSE
    )

    csv(
      comparison,
      file.path(
        out_group,
        "plasticity_effect_comparison.csv"
      )
    )

    # -----------------------------------------------------------------------
    # Nested comparison: does geography add information to plasticity-only?
    # 4 added fixed terms = lat, lon, year:lat, year:lon.
    # -----------------------------------------------------------------------

    source_checks <- p$source_fit_checks

    need(
      nrow(source_checks) == 1L,
      paste(group, "unexpected source fit-check table.")
    )

    LR <- 2 * (
      checks$logLik -
      source_checks$logLik
    )

    df_add <- checks$parameter_df -
      source_checks$parameter_df

    nested <- data.frame(
      group = group,
      source_model = "plasticity_only_spatial",
      bridge_model = "plasticity_plus_geography_spatial",
      logLik_source = source_checks$logLik,
      logLik_bridge = checks$logLik,
      AIC_source = source_checks$AIC,
      AIC_bridge = checks$AIC,
      delta_AIC_bridge_minus_source =
        checks$AIC - source_checks$AIC,
      LR_statistic = max(0, LR),
      added_parameters = df_add,
      p_LRT = if (
        is.finite(LR) &&
        LR >= -1e-5 &&
        df_add > 0
      ) {
        stats::pchisq(
          max(0, LR),
          df = df_add,
          lower.tail = FALSE
        )
      } else NA_real_,
      note = paste0(
        "Nested test of all four added geographic fixed terms together: ",
        "latitude, longitude, year:latitude, year:longitude."
      ),
      stringsAsFactors = FALSE
    )

    csv(
      nested,
      file.path(
        out_group,
        "nested_geography_addition_test.csv"
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
        "lat beta",
        signif(
          beta[
            "year_decade:latitude_within_species_10deg"
          ],
          5
        ),
        "| onset beta",
        signif(
          beta[
            "year_decade:onset_plasticity_z"
          ],
          5
        ),
        "| offset beta",
        signif(
          beta[
            "year_decade:offset_plasticity_bio_z"
          ],
          5
        ),
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

  combine_outputs <- function(cfg, statuses) {

    key_all <- list()
    comp_all <- list()
    nested_all <- list()
    checks_all <- list()

    for (group in GROUPS) {

      folder <- file.path(
        cfg$out,
        paste0(group, "__bridge")
      )

      files <- list(
        key = file.path(
          folder,
          "key_bridge_effects.csv"
        ),
        comp = file.path(
          folder,
          "plasticity_effect_comparison.csv"
        ),
        nested = file.path(
          folder,
          "nested_geography_addition_test.csv"
        ),
        checks = file.path(
          folder,
          "fit_checks.csv"
        )
      )

      if (file.exists(files$key)) {
        key_all[[group]] <-
          utils::read.csv(
            files$key,
            stringsAsFactors = FALSE
          )
      }

      if (file.exists(files$comp)) {
        comp_all[[group]] <-
          utils::read.csv(
            files$comp,
            stringsAsFactors = FALSE
          )
      }

      if (file.exists(files$nested)) {
        nested_all[[group]] <-
          utils::read.csv(
            files$nested,
            stringsAsFactors = FALSE
          )
      }

      if (file.exists(files$checks)) {
        checks_all[[group]] <-
          utils::read.csv(
            files$checks,
            stringsAsFactors = FALSE
          )
      }
    }

    if (length(key_all)) {
      csv(
        do.call(rbind, key_all),
        file.path(
          cfg$out,
          "BRIDGE_key_effects.csv"
        )
      )
    }

    if (length(comp_all)) {
      csv(
        do.call(rbind, comp_all),
        file.path(
          cfg$out,
          "BRIDGE_plasticity_effect_comparison.csv"
        )
      )
    }

    if (length(nested_all)) {
      csv(
        do.call(rbind, nested_all),
        file.path(
          cfg$out,
          "BRIDGE_nested_geography_tests.csv"
        )
      )
    }

    if (length(checks_all)) {
      csv(
        do.call(rbind, checks_all),
        file.path(
          cfg$out,
          "BRIDGE_fit_checks.csv"
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

    check_packages()
    disk_check(cfg$out, cfg$min_disk)
    load_dll(cfg$source_run)

    tryCatch({

      phase(
        cfg,
        "PREPARING",
        "Preparing both bridge inputs"
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
          "bridge_data_summary.csv"
        )
      )

      statuses <- data.frame(
        group = character(),
        status = character(),
        stringsAsFactors = FALSE
      )

      # One fit at a time.
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

        gc()
      }

      combine_outputs(
        cfg,
        statuses
      )

      writeLines(
        c(
          "phenoIMPACT | GEOGRAPHY + PLASTICITY BRIDGE MODELS",
          "",
          paste("Source run:", cfg$source_run),
          "",
          "Primary question:",
          "Does the existing plasticity-demography association persist after",
          "controlling within-species latitude and longitude?",
          "",
          "Model:",
          "ABUND_INDEX ~ year_decade * (latitude_within_species_10deg +",
          "longitude_within_species_10deg + onset_plasticity_z +",
          "offset_plasticity_bio_z), plus the SAME species/site temporal",
          "random effects, population AR(1), Gamma/log likelihood and annual",
          "IID Matern spatial field as the final abundance models.",
          "",
          "Latitude/longitude are population-weighted and centred within species",
          "inside each fitted voltinism subset.",
          "",
          "Plasticity variables retain the exact scaling from the final",
          "plasticity-abundance models. Thus this bridge is directly comparable",
          "to those models.",
          "",
          "The bridge is NOT a causal mediation analysis.",
          "The primary interpretation is statistical overlap / independent",
          "association after mutual adjustment.",
          "",
          "Species-centred plasticity variables are exported but are NOT used",
          "in this primary fit. They can support a stricter within-species",
          "plasticity sensitivity after the bridge results are reviewed.",
          "",
          "BRIDGE_key_effects.csv contains the four focal year interactions.",
          "BRIDGE_plasticity_effect_comparison.csv compares onset/offset effects",
          "before vs after adding geography.",
          "BRIDGE_nested_geography_tests.csv tests all four added geography fixed",
          "terms jointly against the final plasticity-only spatial model."
        ),
        file.path(
          cfg$out,
          "README_bridge.txt"
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
        "bridge_results_to_review.zip"
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

      phase(
        cfg,
        terminal_status,
        zip_path,
        terminal = TRUE
      )

      list(
        status = terminal_status,
        output = cfg$out,
        zip = zip_path
      )

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
      "bridge_workflow_engine.R"
    )

    dump(
      c(
        "VERSION",
        "ROOT",
        "SOURCE_RUN",
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
    source_run = SOURCE_RUN,
    monitor = FALSE,
    minimum_free_GiB = 25,
    optimizer_iterations = 1200L,
    optimizer_evaluations = 1800L
  ) {

    check_packages()

    current <- state$job

    if (
      !is.null(current) &&
      current$process$is_alive()
    ) {
      message("Bridge workflow already running.")
      return(invisible(current))
    }

    project_root <- normalizePath(
      project_root,
      winslash = "/",
      mustWork = TRUE
    )

    source_run <- normalizePath(
      source_run,
      winslash = "/",
      mustWork = TRUE
    )

    fits <- current_fits(
      source_run
    )

    out_root <- file.path(
      project_root,
      "output",
      "population_trends",
      "bridge_geography_plasticity_spatial_15km"
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
        format(
          Sys.time(),
          "%Y%m%d_%H%M%S"
        ),
        "_",
        Sys.getpid()
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
      project_root = project_root,
      source_run = source_run,
      current_fits = fits,
      out = out,
      progress = file.path(
        out,
        "progress.rds"
      ),
      min_disk = minimum_free_GiB,
      iter = as.integer(
        optimizer_iterations
      ),
      eval = as.integer(
        optimizer_evaluations
      )
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

    phase(
      cfg,
      "QUEUED",
      paste(
        "Two bridge fits:",
        paste(GROUPS, collapse = " -> ")
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

        e$worker(cfg)
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

    job <- list(
      process = child,
      out = out,
      log = log_file,
      started = Sys.time()
    )

    state$job <- job

    options(
      phenoimpact.bridge.job = job
    )

    message(
      "Bridge workflow launched. Keep R open."
    )

    message(
      "Output: ",
      out
    )

    message(
      "Monitor: pheno_bridge$status()"
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
        "phenoimpact.bridge.job"
      )
    }

    need(
      !is.null(job),
      "No registered bridge job."
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

    cat(
      "\nAlive:",
      job$process$is_alive(),
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

    if (!job$process$is_alive()) {

      cat(
        "Exit status:",
        job$process$get_exit_status(),
        "\n"
      )

      key <- file.path(
        job$out,
        "BRIDGE_key_effects.csv"
      )

      if (file.exists(key)) {

        cat(
          "\nKey bridge effects:\n"
        )

        print(
          utils::read.csv(
            key,
            stringsAsFactors = FALSE
          ),
          row.names = FALSE
        )
      }

      comp <- file.path(
        job$out,
        "BRIDGE_plasticity_effect_comparison.csv"
      )

      if (file.exists(comp)) {

        cat(
          "\nPlasticity before vs after geography:\n"
        )

        print(
          utils::read.csv(
            comp,
            stringsAsFactors = FALSE
          ),
          row.names = FALSE
        )
      }

      nested <- file.path(
        job$out,
        "BRIDGE_nested_geography_tests.csv"
      )

      if (file.exists(nested)) {

        cat(
          "\nDoes geography add information to plasticity-only model?\n"
        )

        print(
          utils::read.csv(
            nested,
            stringsAsFactors = FALSE
          ),
          row.names = FALSE
        )
      }

      zip_file <- file.path(
        job$out,
        "bridge_results_to_review.zip"
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
    run = run,
    status = status
  )
})


run_bridge_geography_plasticity <- function(...) {
  pheno_bridge$run(...)
}

status_bridge_geography_plasticity <- function(...) {
  pheno_bridge$status(...)
}
