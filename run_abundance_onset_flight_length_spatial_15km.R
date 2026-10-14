# =============================================================================
# phenoIMPACT | Population trends: ONSET + direct FLIGHT-LENGTH plasticity
# Version 1.0 | 2026-10-12
#
# SCIENTIFIC QUESTION
# At equal onset plasticity, do populations that lengthen their flight season
# more strongly in warm years show more favourable abundance trends?
#
# Model fitted separately for univoltines and multivoltines:
#
#   Gamma/log abundance ~
#     year_decade * onset_plasticity_z +
#     year_decade * flight_length_plasticity_z
#
# with the same:
#   - abundance response and rows as the final QC abundance models
#   - species intercepts + temporal slopes
#   - site intercepts + temporal slopes
#   - stationary population AR(1) with true calendar gaps
#   - annual IID Matern spatial field
#   - reviewed 15-km mesh and unchanged C++ likelihood
#
# Only the second phenological predictor changes:
#   FINAL offset plasticity -> FINAL spatial direct flight-length plasticity.
#
# Positive flight_length_plasticity means a longer flight period in warmer
# years. Because the group mean slope is negative, a higher value may also mean
# "shortens less" rather than literally lengthens for every population.
#
# The key coefficient is:
#   year_decade:flight_length_plasticity_z
#
# It is conditional on onset_plasticity_z because both interactions are fitted
# simultaneously.
#
# Sources are explicit:
#   abundance rows / onset predictor / AR1 / mesh:
#     spatial_abundance_qc_offset_final_15km/prep_6c08f872b268
#   reviewed warm-start abundance fits:
#     spatial_abundance_qc_offset_final_15km/run_eee755f6aac4
#   spatial direct flight-length plasticity:
#     flight_length_direct_spatial_15km/run_bc5b3c03313a/
#       all_population_flight_length_plasticity_spatial.csv
#
# All abundance parameters are re-estimated. First-stage plasticity uncertainty
# is not propagated, consistently with the primary two-stage analysis.
#
# Run:
#   source(file.choose(), encoding = "UTF-8")
#   onset_flight_job <- pheno_onset_flight_abundance$run(monitor = FALSE)
#   pheno_onset_flight_abundance$status()
#
# Keep R open. Do not run another large TMB job concurrently.
# =============================================================================

pheno_onset_flight_abundance <- local({

  VERSION <- "onset_flight_length_spatial_abundance_1.0"
  ROOT <- "E:/phenoIMPACT project/code/phenoIMPACT"
  SOURCE_PREP <- file.path(
    ROOT, "output", "population_trends",
    "spatial_abundance_qc_offset_final_15km", "prep_6c08f872b268"
  )
  SOURCE_ABUNDANCE_RUN <- file.path(
    ROOT, "output", "population_trends",
    "spatial_abundance_qc_offset_final_15km", "run_eee755f6aac4"
  )
  FLIGHT_LENGTH_FILE <- file.path(
    ROOT, "output", "phenology_plasticity",
    "flight_length_direct_spatial_15km", "run_bc5b3c03313a",
    "all_population_flight_length_plasticity_spatial.csv"
  )
  OUTPUT_ROOT <- file.path(
    ROOT, "output", "population_trends",
    "spatial_abundance_onset_flight_length_15km"
  )
  GROUPS <- c("univoltine", "multivoltine")
  CPP_MD5 <- "88e15b8506c2e85a2e1e0d3839101530"
  GRADIENT_TOL <- 0.05
  CORE <- c("TMB", "Matrix", "fmesher")
  state <- new.env(parent = emptyenv())

  need <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)
  same_num <- function(a, b, tol = 1e-8) {
    is.numeric(a) && is.numeric(b) &&
      length(a) == length(b) &&
      all(is.finite(a)) && all(is.finite(b)) &&
      (length(a) == 0L ||
         max(abs(as.numeric(a) - as.numeric(b))) <= tol)
  }
  csv <- function(x, path) {
    utils::write.csv(x, path, row.names = FALSE, na = "")
  }
  stamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

  atomic_rds <- function(x, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    tmp <- tempfile(".writing_", tmpdir = dirname(path), fileext = ".rds")
    on.exit(unlink(tmp), add = TRUE)
    saveRDS(x, tmp, compress = "gzip")

    bak <- paste0(path, ".previous")
    if (file.exists(bak)) unlink(bak)
    had <- file.exists(path)
    if (had) {
      need(file.rename(path, bak), paste("Cannot retain previous checkpoint:", path))
    }
    if (!file.rename(tmp, path)) {
      if (had) file.rename(bak, path)
      stop("Cannot finalize checkpoint: ", path, call. = FALSE)
    }
    if (file.exists(bak)) unlink(bak)
    invisible(path)
  }

  phase <- function(cfg, stage, group = "", detail = "", terminal = FALSE) {
    x <- list(
      time = Sys.time(),
      stage = stage,
      group = group,
      detail = detail,
      terminal = terminal
    )
    atomic_rds(x, file.path(cfg$out, "progress.rds"))
    message("[", stamp(), "] ", group, " | ", stage, " | ", detail)
  }

  check_packages <- function() {
    pkgs <- c(CORE, "callr", "ps", "digest")
    miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
    need(!length(miss),
         paste("Missing packages:", paste(miss, collapse = ", ")))
    invisible(TRUE)
  }

  package_versions <- function() {
    data.frame(
      package = CORE,
      version = vapply(
        CORE,
        function(p) as.character(utils::packageVersion(p)),
        character(1)
      ),
      stringsAsFactors = FALSE
    )
  }

  disk_check <- function(path, minimum_GiB = 50, verbose = TRUE) {
    d <- ps::ps_disk_usage(path)
    need(nrow(d) == 1L && is.finite(d$available),
         "Could not determine available disk space.")
    free <- as.numeric(d$available) / 1024^3
    if (verbose) message("Disk preflight: ", round(free, 1), " GiB available.")
    need(free >= minimum_GiB,
         paste("Only", round(free, 1), "GiB available; require at least",
               minimum_GiB, "GiB."))
    invisible(free)
  }

  fixed_table <- function(beta, V) {
    se <- sqrt(diag(V))
    need(all(is.finite(beta)) &&
           all(is.finite(se)) &&
           all(se > 0),
         "Fixed estimates/SE are not finite and positive.")
    z <- beta / se
    data.frame(
      term = names(beta),
      estimate = as.numeric(beta),
      std.error = as.numeric(se),
      conf.low = as.numeric(beta - 1.96 * se),
      conf.high = as.numeric(beta + 1.96 * se),
      p_Wald = as.numeric(2 * stats::pnorm(-abs(z))),
      stringsAsFactors = FALSE
    )
  }

  predictions <- function(beta, V, years) {
    rows <- list()
    for (v in c("onset_plasticity_z", "flight_length_plasticity_z")) {
      for (z in c(-1, 0, 1)) {
        it <- paste0("year_decade:", v)
        need(all(c("year_decade", it) %in% names(beta)),
             paste("Missing prediction term:", it))
        tt <- (years - min(years)) / 10
        slope <- beta["year_decade"] + beta[it] * z
        vv <- V["year_decade", "year_decade"] +
          z^2 * V[it, it] +
          2 * z * V["year_decade", it]
        need(is.finite(vv) && vv >= -1e-10,
             "Invalid prediction variance.")
        se <- sqrt(max(vv, 0)) * tt
        eta <- as.numeric(slope) * tt
        rows[[length(rows) + 1L]] <- data.frame(
          variable = v,
          plasticity_z = z,
          YEAR = years,
          percent_change = 100 * expm1(eta),
          low = 100 * expm1(eta - 1.96 * se),
          high = 100 * expm1(eta + 1.96 * se)
        )
      }
    }
    do.call(rbind, rows)
  }

  checked_covariance <- function(V, n) {
    need(is.matrix(V) && identical(dim(V), as.integer(c(n, n))) &&
           all(is.finite(V)),
         "Marginal covariance is missing/non-finite.")
    sym <- max(abs(V - t(V))) / max(1, max(abs(V)))
    need(sym < 1e-8, "Marginal covariance is not numerically symmetric.")
    V <- (V + t(V)) / 2
    need(all(diag(V) > 0), "Marginal covariance has non-positive diagonal.")
    ch <- tryCatch(chol(V), error = function(e) NULL)
    need(!is.null(ch), "Marginal covariance is not positive definite.")
    V
  }

  scaled_gradient <- function(g, V) {
    V <- checked_covariance(V, length(g))
    x <- as.numeric(crossprod(g, V %*% g))
    need(is.finite(x) && x >= 0, "Invalid covariance-scaled gradient.")
    sqrt(x)
  }

  newton_direction <- function(g, V) {
    V <- checked_covariance(V, length(g))
    norm <- scaled_gradient(g, V)
    shrink <- if (norm > 1) 1 / norm else 1
    direction <- -as.numeric(V %*% g) * shrink
    dot <- sum(g * direction)
    need(all(is.finite(direction)) && is.finite(dot) && dot < 0,
         "Could not obtain a finite descent direction.")
    list(direction = direction, dot = dot, shrink = shrink)
  }

  line_search <- function(
    par, objective, gradient, covariance, fn,
    max_backtracks = 6L
  ) {
    step <- newton_direction(gradient, covariance)
    rec <- vector("list", max_backtracks + 1L)

    for (j in 0:max_backtracks) {
      alpha <- 0.5^j
      p <- par + alpha * step$direction
      value <- fn(p)
      armijo <- objective + 1e-4 * alpha * step$dot
      accept <- is.finite(value) && value < objective && value <= armijo

      rec[[j + 1L]] <- data.frame(
        backtrack = j,
        alpha = alpha,
        direction_shrink = step$shrink,
        objective_before = objective,
        objective_trial = value,
        Armijo_bound = armijo,
        accepted = accept
      )

      if (accept) {
        return(list(
          accepted = TRUE,
          par = p,
          objective = value,
          alpha = alpha,
          history = do.call(rbind, rec[seq_len(j + 1L)])
        ))
      }
    }

    list(
      accepted = FALSE,
      par = par,
      objective = objective,
      alpha = 0,
      history = do.call(rbind, rec)
    )
  }

  self_test <- function() {
    H <- matrix(c(4, 1, 1, 3), 2, 2)
    V <- solve(H)
    p <- c(.02, -.03)
    fn <- function(x) as.numeric(crossprod(x, H %*% x)) / 2
    gr <- function(x) as.numeric(H %*% x)
    s <- line_search(p, fn(p), gr(p), V, fn, max_backtracks = 6L)
    need(s$accepted && max(abs(s$par)) < 1e-10,
         "Newton/backtracking self-test failed.")
    invisible(TRUE)
  }

  # Exact reviewed likelihood construction.
  with_mesh <- function(input, mesh) {
    data <- input$data
    data$use_spatial <- 1L

    data$A <- methods::as(
      fmesher::fm_basis(
        mesh$mesh,
        loc = as.matrix(input$sites[c("x_km", "y_km")])
      ),
      "CsparseMatrix"
    )
    need(all(abs(Matrix::rowSums(data$A) - 1) < 1e-8),
         "Some retained site coordinates fall outside the saved mesh.")

    data$M0 <- mesh$M0
    data$M1 <- mesh$M1
    data$M2 <- mesh$M2

    list(data = data, par = input$par)
  }

  make_objective <- function(x) {
    data <- x$data
    for (nm in c("A", "M0", "M1", "M2")) {
      data[[nm]] <- methods::as(
        methods::as(
          methods::as(data[[nm]], "dMatrix"),
          "generalMatrix"
        ),
        "TsparseMatrix"
      )
    }

    obj <- TMB::MakeADFun(
      data = data,
      parameters = x$par,
      random = c("b_species", "b_site", "ar_obs", "field"),
      map = NULL,
      DLL = "trend_spde",
      silent = TRUE,
      inner.control = list(maxit = 1000, trace = FALSE)
    )

    list(obj = obj, data = data)
  }

  current_parameters <- function(obj, fixed_par) {
    # Ensure latent optimum corresponds to this fixed point.
    value <- obj$fn(fixed_par)
    need(is.finite(value), "Non-finite objective while extracting parameters.")
    ans <- obj$env$parList(par = obj$env$last.par)

    fixed_names <- c(
      "beta", "log_sd", "rho_raw", "log_shape",
      "log_range", "log_spatial_sd"
    )
    for (nm in fixed_names) {
      need(
        same_num(
          as.numeric(ans[[nm]]),
          as.numeric(fixed_par[names(fixed_par) == nm]),
          1e-7
        ),
        paste("Fixed parameter mapping mismatch:", nm)
      )
    }

    need(all(vapply(
      ans,
      function(z) is.numeric(z) && all(is.finite(z)),
      logical(1)
    )), "Non-finite current parameter array.")

    ans
  }

  assess <- function(obj, par) {
    warnings <- character()
    catcher <- function(w) {
      warnings <<- unique(c(warnings, conditionMessage(w)))
      message("WARNING: ", conditionMessage(w))
      invokeRestart("muffleWarning")
    }

    value <- withCallingHandlers(obj$fn(par), warning = catcher)
    gradient <- withCallingHandlers(as.numeric(obj$gr(par)), warning = catcher)

    sdr <- withCallingHandlers(
      TMB::sdreport(
        obj,
        par.fixed = par,
        getJointPrecision = FALSE,
        getReportCovariance = FALSE
      ),
      warning = catcher
    )

    V <- sdr$cov.fixed
    covariance_ok <- tryCatch({
      checked_covariance(V, length(par))
      TRUE
    }, error = function(e) FALSE)

    norm <- if (covariance_ok && all(is.finite(gradient))) {
      scaled_gradient(gradient, V)
    } else {
      Inf
    }

    pass <- isTRUE(sdr$pdHess) &&
      covariance_ok &&
      all(is.finite(gradient)) &&
      is.finite(norm) &&
      norm < GRADIENT_TOL &&
      is.finite(value)

    list(
      objective = as.numeric(value),
      gradient = gradient,
      covariance = V,
      pdHess = isTRUE(sdr$pdHess),
      scaled_gradient = norm,
      covariance_ok = covariance_ok,
      pass = pass,
      warnings = warnings
    )
  }

  prepare_onset_flight_inputs <- function(cfg) {
    check_packages()
    dir.create(cfg$prep, recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(cfg$prep, "compiled"), showWarnings = FALSE)

    flight <- utils::read.csv(
      cfg$flight_length_file,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    required_flight <- c(
      "pop_id", "SPECIES", "SITE_ID", "voltinism",
      "flight_length_plasticity_days_per_C",
      "flight_length_plasticity_z"
    )
    need(
      all(required_flight %in% names(flight)),
      "Spatial flight-length plasticity table lacks required columns."
    )
    for (nm in c("pop_id", "SPECIES", "SITE_ID", "voltinism")) {
      flight[[nm]] <- as.character(flight[[nm]])
    }
    need(
      !anyDuplicated(flight$pop_id) &&
        all(flight$voltinism %in% GROUPS) &&
        all(is.finite(flight$flight_length_plasticity_days_per_C)) &&
        all(is.finite(flight$flight_length_plasticity_z)),
      "Invalid/duplicated spatial flight-length plasticities."
    )

    mesh_source <- file.path(cfg$source_prep, "mesh.rds")
    cpp_source <- file.path(
      cfg$source_abundance_run, "compiled", "trend_spde.cpp"
    )
    need(file.exists(mesh_source), "Source mesh.rds is missing.")
    need(file.exists(cpp_source), "Source trend_spde.cpp is missing.")
    need(
      identical(unname(tools::md5sum(cpp_source)), CPP_MD5),
      "Source C++ likelihood differs from the reviewed model."
    )
    need(
      file.copy(mesh_source, file.path(cfg$prep, "mesh.rds"), overwrite = TRUE),
      "Could not copy the reviewed mesh."
    )
    need(
      file.copy(
        cpp_source,
        file.path(cfg$prep, "compiled", "trend_spde.cpp"),
        overwrite = TRUE
      ),
      "Could not copy the reviewed C++ likelihood."
    )

    summary_rows <- list()
    scaling_rows <- list()

    for (group in GROUPS) {
      phase(
        cfg, "PREPARING_INPUT", group,
        "Replacing only offset plasticity with spatial direct flight-length plasticity"
      )

      source_input_path <- file.path(
        cfg$source_prep, paste0("input_", group, ".rds")
      )
      source_fit_path <- file.path(
        cfg$source_abundance_run, paste0(group, "__full"), "fit.rds"
      )
      need(file.exists(source_input_path), paste("Missing:", source_input_path))
      need(file.exists(source_fit_path), paste("Missing:", source_fit_path))

      old <- readRDS(source_input_path)
      source_fit <- readRDS(source_fit_path)

      need(
        identical(old$group, group) &&
          identical(source_fit$group, group) &&
          identical(source_fit$variant, "full") &&
          isTRUE(source_fit$checks$numerical_checks_ok) &&
          identical(source_fit$data_signature, old$data_signature),
        paste("Source abundance input/fit is inconsistent:", group)
      )
      need(
        identical(
          old$data_signature,
          digest::digest(list(old$frame, old$data$X), algo = "sha256")
        ),
        paste("Source abundance data signature cannot be reproduced:", group)
      )

      p <- flight[flight$voltinism == group, , drop = FALSE]
      need(nrow(p) > 1L, paste("No flight-length populations:", group))

      frame <- old$frame
      frame$pop_id <- as.character(frame$pop_id)
      frame$SPECIES <- as.character(frame$SPECIES)
      frame$SITE_ID <- as.character(frame$SITE_ID)

      ix <- match(frame$pop_id, p$pop_id)
      need(
        !anyNA(ix),
        paste(
          "Some final abundance populations lack spatial flight-length plasticity:",
          group
        )
      )
      need(
        identical(frame$SPECIES, p$SPECIES[ix]) &&
          identical(frame$SITE_ID, p$SITE_ID[ix]),
        paste("Population identity mismatch:", group)
      )

      frame$flight_length_plasticity_days_per_C <-
        p$flight_length_plasticity_days_per_C[ix]
      frame$flight_length_plasticity_z <-
        p$flight_length_plasticity_z[ix]

      need(
        all(is.finite(frame$onset_plasticity_z)) &&
          all(is.finite(frame$flight_length_plasticity_z)),
        paste("Non-finite phenological predictors:", group)
      )

      X <- stats::model.matrix(
        ~ year_decade * onset_plasticity_z +
          year_decade * flight_length_plasticity_z,
        frame
      )
      need(
        nrow(X) == nrow(frame) &&
          all(is.finite(X)) &&
          qr(X)$rank == ncol(X),
        paste("New abundance fixed design is invalid:", group)
      )

      expected <- c(
        "(Intercept)",
        "year_decade",
        "onset_plasticity_z",
        "flight_length_plasticity_z",
        "year_decade:onset_plasticity_z",
        "year_decade:flight_length_plasticity_z"
      )
      need(
        identical(colnames(X), expected),
        paste(
          "Unexpected new design-column order:",
          paste(colnames(X), collapse = ", ")
        )
      )

      data <- old$data
      old_X <- data$X
      data$X <- X

      # Everything except X remains bit-for-bit identical.
      other <- setdiff(names(data), "X")
      need(
        identical(data[other], old$data[other]),
        paste("Response/group/time/AR1 structure changed:", group)
      )
      need(
        same_num(data$y, old$data$y, 0) &&
          same_num(data$time, old$data$time, 0),
        paste("Response or year scale changed:", group)
      )

      # Warm start from the reviewed FINAL onset+offset abundance fit.
      pars <- source_fit$parameters
      old_beta <- source_fit$beta
      need(
        all(c(
          "(Intercept)", "year_decade", "onset_plasticity_z",
          "offset_plasticity_bio_z",
          "year_decade:onset_plasticity_z",
          "year_decade:offset_plasticity_bio_z"
        ) %in% names(old_beta)),
        paste("Unexpected source abundance beta names:", group)
      )

      beta_start <- stats::setNames(rep(0, ncol(X)), colnames(X))
      common <- intersect(names(beta_start), names(old_beta))
      beta_start[common] <- old_beta[common]

      # The old offset coefficients are reasonable standardized warm starts
      # for the new direct flight-length terms; all coefficients remain free.
      beta_start["flight_length_plasticity_z"] <-
        old_beta["offset_plasticity_bio_z"]
      beta_start["year_decade:flight_length_plasticity_z"] <-
        old_beta["year_decade:offset_plasticity_bio_z"]

      need(all(is.finite(beta_start)), "Non-finite beta warm start.")
      pars$beta <- unname(beta_start)

      input <- list(
        signature = cfg$signature,
        schema = VERSION,
        group = group,
        data = data,
        frame = frame,
        sites = old$sites,
        years = old$years,
        par = pars,
        data_signature = digest::digest(list(frame, X), algo = "sha256"),
        source_data_signature = old$data_signature,
        source_fit = normalizePath(
          source_fit_path, winslash = "/", mustWork = TRUE
        ),
        source_input = normalizePath(
          source_input_path, winslash = "/", mustWork = TRUE
        ),
        onset_predictor_preserved_exactly = TRUE,
        flight_length_predictor_spatial = TRUE,
        flight_length_units = "days per degree C; standardized within voltinism before abundance join",
        abundance_rows_unchanged = TRUE,
        AR1_structure_unchanged = TRUE,
        mesh_unchanged = TRUE,
        diagnostics_status = "PENDING"
      )

      atomic_rds(
        input,
        file.path(cfg$prep, paste0("input_", group, ".rds"))
      )

      pop <- unique(frame[, c(
        "pop_id", "onset_plasticity_z",
        "offset_plasticity_bio_z",
        "flight_length_plasticity_z",
        "flight_length_plasticity_days_per_C"
      )])
      summary_rows[[group]] <- data.frame(
        group = group,
        n_observations = nrow(frame),
        n_populations = nrow(pop),
        all_abundance_populations_matched = TRUE,
        onset_predictor_unchanged = same_num(
          frame$onset_plasticity_z,
          old$frame$onset_plasticity_z,
          0
        ),
        response_rows_AR1_mesh_unchanged = TRUE,
        cor_onset_vs_flight_length = stats::cor(
          pop$onset_plasticity_z,
          pop$flight_length_plasticity_z
        ),
        cor_old_offset_vs_flight_length = stats::cor(
          pop$offset_plasticity_bio_z,
          pop$flight_length_plasticity_z
        ),
        min_flight_length_z = min(pop$flight_length_plasticity_z),
        max_flight_length_z = max(pop$flight_length_plasticity_z),
        stringsAsFactors = FALSE
      )

      scaling_rows[[group]] <- data.frame(
        group = group,
        n_spatial_flight_length_populations = nrow(p),
        mean_days_per_C = mean(p$flight_length_plasticity_days_per_C),
        SD_days_per_C = stats::sd(p$flight_length_plasticity_days_per_C),
        scaling_universe =
          "all populations in accepted spatial flight-length fit; before abundance join",
        stringsAsFactors = FALSE
      )

      rm(old, source_fit, frame, data, X, pars, input, p, ix)
      invisible(gc())
    }

    csv(
      do.call(rbind, summary_rows),
      file.path(cfg$prep, "preparation_summary.csv")
    )
    csv(
      do.call(rbind, scaling_rows),
      file.path(cfg$prep, "flight_length_scaling.csv")
    )
    atomic_rds(
      list(
        signature = cfg$signature,
        status = "PREPARED",
        time = Sys.time()
      ),
      file.path(cfg$prep, "preparation_complete.rds")
    )

    invisible(TRUE)
  }

  compile_worker <- function(cfg) {
    check_packages()
    source_cpp <- file.path(cfg$prep, "compiled", "trend_spde.cpp")
    target_dir <- file.path(cfg$out, "compiled")
    target_cpp <- file.path(target_dir, "trend_spde.cpp")
    marker <- file.path(target_dir, "compile_complete.rds")
    dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)

    need(
      identical(unname(tools::md5sum(source_cpp)), CPP_MD5),
      "Prepared C++ likelihood differs from the reviewed template."
    )

    if (file.exists(marker) &&
        file.exists(TMB::dynlib(file.path(target_dir, "trend_spde")))) {
      x <- readRDS(marker)
      need(identical(x$signature, cfg$signature),
           "Compiled cache signature mismatch.")
      return("COMPILED_CACHED")
    }

    need(file.copy(source_cpp, target_cpp, overwrite = TRUE),
         "Could not copy C++ likelihood.")
    need(identical(unname(tools::md5sum(target_cpp)), CPP_MD5),
         "Copied C++ likelihood MD5 mismatch.")

    phase(cfg, "COMPILING", "", "Compiling unchanged reviewed TMB likelihood")

    oldwd <- getwd()
    on.exit(setwd(oldwd), add = TRUE)
    setwd(target_dir)

    code <- TMB::compile("trend_spde.cpp", flags = "-O1 -g0")
    need(isTRUE(code == 0L) &&
           file.exists(TMB::dynlib("trend_spde")),
         "TMB compilation failed.")

    atomic_rds(
      list(
        signature = cfg$signature,
        cpp_md5 = CPP_MD5,
        package_versions = package_versions()
      ),
      marker
    )

    "COMPILED"
  }

  load_dll <- function(cfg) {
    dll <- TMB::dynlib(file.path(cfg$out, "compiled", "trend_spde"))
    need(file.exists(dll), "Compiled trend_spde DLL missing.")
    dyn.load(dll)
    TMB::openmp(n = 1L, DLL = "trend_spde")
    invisible(dll)
  }

  source_fit_for_comparison <- function(input) {
    need(is.character(input$source_fit) &&
           length(input$source_fit) == 1L &&
           file.exists(input$source_fit),
         "Prepared input source_fit is missing.")
    fit <- readRDS(input$source_fit)
    need(identical(fit$group, input$group) &&
           identical(fit$variant, "full") &&
           isTRUE(fit$checks$numerical_checks_ok),
         "Prepared warm-start source fit is not the reviewed FIT_OK full model.")
    fit
  }

  export_fit <- function(fit, input, cfg, group) {
    out <- file.path(cfg$out, paste0(group, "__full"))
    dir.create(out, recursive = TRUE, showWarnings = FALSE)

    csv(fit$checks, file.path(out, "fit_checks.csv"))
    csv(fixed_table(fit$beta, fit$V), file.path(out, "fixed_effects.csv"))
    csv(
      data.frame(term = rownames(fit$V), fit$V, check.names = FALSE),
      file.path(out, "fixed_effect_covariance.csv")
    )
    csv(
      predictions(fit$beta, fit$V, input$years),
      file.path(out, "predictions_spatial.csv")
    )

    src <- source_fit_for_comparison(input)
    before <- fixed_table(src$beta, src$V)
    after <- fixed_table(fit$beta, fit$V)

    comp <- merge(
      before, after,
      by = "term",
      all = TRUE,
      sort = FALSE,
      suffixes = c("_before_final_QC", "_after_final_QC")
    )
    comp$estimate_change <-
      comp$estimate_after_final_QC - comp$estimate_before_final_QC
    comp$change_in_previous_SE <-
      comp$estimate_change / comp$std.error_before_final_QC

    csv(comp, file.path(out, "coefficient_comparison.csv"))

    focal <- after[
      after$term %in% c(
        "year_decade:onset_plasticity_z",
        "year_decade:flight_length_plasticity_z"
      ),
      ,
      drop = FALSE
    ]
    focal$ratio_of_decadal_factors_per_SD <- exp(focal$estimate)
    focal$ratio_conf.low <- exp(focal$conf.low)
    focal$ratio_conf.high <- exp(focal$conf.high)
    focal$question <- c(
      "Effect of stronger onset advancement plasticity, conditional on direct flight-length plasticity.",
      "Effect of stronger flight-length extension plasticity, conditional on onset plasticity."
    )[match(
      focal$term,
      c(
        "year_decade:onset_plasticity_z",
        "year_decade:flight_length_plasticity_z"
      )
    )]
    csv(focal, file.path(out, "focal_interactions.csv"))

    csv(
      data.frame(
        SITE_ID = rep(input$sites$SITE_ID, times = length(input$years)),
        YEAR = rep(input$years, each = nrow(input$sites)),
        field = as.vector(fit$field_at_sites)
      ),
      file.path(out, "annual_spatial_field.csv")
    )

    writeLines(unique(fit$warnings), file.path(out, "warnings.txt"))

    writeLines(c(
      paste("ONSET + direct FLIGHT-LENGTH abundance model:", group),
      paste("Observations:", nrow(input$frame)),
      paste("Populations:", length(unique(input$frame$pop_id))),
      "ONSET predictor: exactly preserved from the final QC abundance input.",
      "FLIGHT-LENGTH predictor: spatially corrected direct seasonal-duration plasticity.",
      "Positive flight-length values mean a longer season in warmer years;",
      "relative values can also represent less shortening when the group mean is negative.",
      "Response rows, year scale, population AR1 links/gaps and 15-km mesh are unchanged.",
      "Gamma/log; species/site intercepts and temporal slopes; population AR1;",
      "annual IID Matern spatial field.",
      "Both phenological predictors and both year interactions are fitted simultaneously.",
      "All parameters are re-estimated.",
      "Key test: year_decade:flight_length_plasticity_z, conditional on onset plasticity.",
      "No reduced LRTs here. p-values are normal-Wald.",
      "No propagation of first-stage plasticity uncertainty.",
      "Residual diagnostics remain PENDING until focal coefficients are reviewed."
    ), file.path(out, "README.txt"))

    invisible(comp)
  }

  fit_worker <- function(cfg, group) {
    check_packages()
    disk_check(cfg$out, cfg$minimum_free_GiB)

    out <- file.path(cfg$out, paste0(group, "__full"))
    dir.create(out, recursive = TRUE, showWarnings = FALSE)

    input <- readRDS(file.path(cfg$prep, paste0("input_", group, ".rds")))
    need(
      identical(input$schema, VERSION) &&
        identical(input$group, group) &&
        isTRUE(input$onset_predictor_preserved_exactly) &&
        isTRUE(input$flight_length_predictor_spatial) &&
        isTRUE(input$abundance_rows_unchanged) &&
        isTRUE(input$AR1_structure_unchanged) &&
        isTRUE(input$mesh_unchanged),
      paste("Unexpected onset+flight-length input specification:", group)
    )
    need(
      identical(
        input$data_signature,
        digest::digest(list(input$frame, input$data$X), algo = "sha256")
      ),
      paste("Prepared input signature cannot be reproduced:", group)
    )

    final_file <- file.path(out, "fit.rds")
    done_file <- file.path(out, "completed.rds")

    if (file.exists(done_file) && file.exists(final_file)) {
      done <- readRDS(done_file)
      fit <- readRDS(final_file)
      need(identical(done$signature, cfg$signature) &&
             identical(fit$signature, cfg$signature) &&
             identical(fit$data_signature, input$data_signature),
           "Completed fit cache mismatch.")
      export_fit(fit, input, cfg, group)
      return(done$status)
    }

    load_dll(cfg)
    mesh <- readRDS(file.path(cfg$prep, "mesh.rds"))

    phase(
      cfg, "BUILDING_MODEL", group,
      paste(nrow(input$frame), "observations;",
            length(unique(input$frame$pop_id)), "populations")
    )

    built <- make_objective(with_mesh(input, mesh))
    obj <- built$obj
    on.exit(try(TMB::FreeADFun(obj), silent = TRUE), add = TRUE)

    # Prepared warm start: verify the FIXED marginal parameters exactly and
    # verify dimensions/indices of the RANDOM arrays. Do NOT require random
    # effects to remain numerically identical after obj$fn(): with Laplace
    # random effects, TMB internally optimizes b_species, b_site, ar_obs and
    # field conditional on the supplied fixed parameters. Their post-evaluation
    # modes are therefore expected to differ from their starting values.
    start_parameters <- input$par
    start_fixed <- obj$par
    expected_names <- c(
      rep("beta", ncol(input$data$X)),
      rep("log_sd", 5L),
      "rho_raw", "log_shape", "log_range", "log_spatial_sd"
    )
    need(identical(names(start_fixed), expected_names),
         "Unexpected marginal parameter order.")

    fixed_start_names <- c(
      "beta", "log_sd", "rho_raw", "log_shape",
      "log_range", "log_spatial_sd"
    )
    for (nm in fixed_start_names) {
      need(
        same_num(
          as.numeric(start_fixed[names(start_fixed) == nm]),
          as.numeric(start_parameters[[nm]]),
          1e-8
        ),
        paste("Prepared FIXED warm start was not reproduced:", nm)
      )
    }

    need(
      is.matrix(start_parameters$b_species) &&
        nrow(start_parameters$b_species) ==
          (max(input$data$species) + 1L) &&
        ncol(start_parameters$b_species) == 2L,
      "Prepared b_species dimensions do not match the retained species index."
    )
    need(
      is.matrix(start_parameters$b_site) &&
        nrow(start_parameters$b_site) ==
          (max(input$data$site) + 1L) &&
        ncol(start_parameters$b_site) == 2L,
      "Prepared b_site dimensions do not match the retained site index."
    )
    need(
      length(start_parameters$ar_obs) == length(input$data$y),
      "Prepared AR1 random-effect vector does not match retained observations."
    )
    need(
      is.matrix(start_parameters$field) &&
        ncol(start_parameters$field) == length(input$years) &&
        nrow(start_parameters$field) == nrow(mesh$M0),
      "Prepared spatial-field dimensions do not match the saved mesh/year universe."
    )

    # One finite likelihood replay verifies that the supplied random warm starts
    # can be used by the rebuilt objective. The resulting random modes are
    # allowed to move, as required by the Laplace inner optimization.
    warm_start_objective <- withCallingHandlers(
      obj$fn(start_fixed),
      warning = function(w) {
        message("Warm-start replay warning: ", conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    need(
      is.finite(warm_start_objective),
      "Prepared warm start produced a non-finite initial likelihood."
    )

    warnings <- character()
    catcher <- function(w) {
      warnings <<- unique(c(warnings, conditionMessage(w)))
      message("WARNING: ", conditionMessage(w))
      invokeRestart("muffleWarning")
    }

    checkpoint <- file.path(out, "last_evaluation.rds")
    optimized_file <- file.path(out, "optimized_parameters.rds")

    optimized <- if (file.exists(optimized_file)) {
      x <- readRDS(optimized_file)
      need(identical(x$signature, cfg$signature) &&
             identical(x$data_signature, input$data_signature),
           "Optimized checkpoint mismatch.")
      x
    } else NULL

    warm <- if (is.null(optimized) && file.exists(checkpoint)) {
      x <- readRDS(checkpoint)
      need(identical(x$signature, cfg$signature) &&
             identical(x$data_signature, input$data_signature),
           "Evaluation checkpoint mismatch.")
      x
    } else NULL

    count <- 0L
    ng <- 0L
    last_save <- Sys.time() - 3600

    extract_arrays <- function(par) {
      current_parameters(obj, par)
    }

    objective <- function(par) {
      value <- withCallingHandlers(obj$fn(par), warning = catcher)
      count <<- count + 1L

      if (is.finite(value) &&
          (count == 1L ||
             as.numeric(difftime(Sys.time(), last_save, units = "secs")) >= 300)) {
        disk_check(cfg$out, 10, verbose = FALSE)
        pars <- extract_arrays(par)
        atomic_rds(
          list(
            signature = cfg$signature,
            data_signature = input$data_signature,
            par = par,
            parameters = pars,
            objective = value,
            time = Sys.time()
          ),
          checkpoint
        )
        phase(
          cfg, "FITTING", group,
          paste("function", count, "gradient", ng,
                "objective", signif(value, 10),
                "| portable checkpoint saved")
        )
        last_save <<- Sys.time()
      }
      value
    }

    gradient <- function(par) {
      ng <<- ng + 1L
      withCallingHandlers(obj$gr(par), warning = catcher)
    }

    started <- Sys.time()

    if (!is.null(optimized)) {
      phase(
        cfg, "RECOVERING_OPTIMIZED_POINT", group,
        "Reusing saved optimized fixed parameters; no optimizer repeat"
      )
      opt <- optimized$optimizer
      need(identical(names(opt$par), names(obj$par)),
           "Saved optimizer parameter names changed.")
      replay <- obj$fn(opt$par)
      need(is.finite(replay) &&
             abs(replay - opt$objective) < .01,
           "Saved optimized objective cannot be reproduced.")

    } else {
      if (!is.null(warm)) {
        need(identical(names(warm$par), names(obj$par)),
             "Warm checkpoint parameter names changed.")
        obj$par <- warm$par
        phase(
          cfg, "FITTING", group,
          "Restarting from last portable evaluation checkpoint"
        )
      } else {
        phase(
          cfg, "FITTING", group,
          "Warm start from reviewed current spatial abundance fit; all parameters free"
        )
      }

      set.seed(20261007L + match(group, GROUPS))

      opt <- withCallingHandlers(
        stats::nlminb(
          obj$par,
          objective,
          gradient,
          control = list(
            iter.max = cfg$optimizer_iterations,
            eval.max = cfg$optimizer_evaluations,
            rel.tol = 1e-10
          )
        ),
        warning = catcher
      )

      if (opt$convergence != 0L) {
        phase(
          cfg, "BFGS_REFINEMENT", group,
          paste("nlminb convergence code", opt$convergence)
        )
        second <- withCallingHandlers(
          stats::optim(
            opt$par,
            objective,
            gradient,
            method = "BFGS",
            control = list(maxit = 300L, reltol = 1e-10)
          ),
          warning = catcher
        )
        if (is.finite(second$value) &&
            second$value <= opt$objective + 1e-6) {
          opt <- list(
            par = second$par,
            objective = second$value,
            convergence = second$convergence,
            message = "BFGS refinement after nlminb",
            first_optimizer = opt,
            counts = second$counts
          )
        }
      }

      replay <- obj$fn(opt$par)
      need(is.finite(replay) &&
             abs(replay - opt$objective) < .01,
           "Final optimizer objective mismatch.")

      arrays <- extract_arrays(opt$par)
      atomic_rds(
        list(
          signature = cfg$signature,
          data_signature = input$data_signature,
          optimizer = opt,
          parameters = arrays
        ),
        optimized_file
      )
    }

    # Initial uncertainty/numerical assessment.
    phase(
      cfg, "UNCERTAINTY", group,
      "Optimized parameters saved; computing marginal Hessian/covariance"
    )

    a <- assess(obj, opt$par)
    warnings <- unique(c(warnings, a$warnings))

    history <- data.frame(
      step = 0L,
      objective = a$objective,
      max_abs_gradient = max(abs(a$gradient)),
      gradient_covariance_norm = a$scaled_gradient,
      pdHess = a$pdHess,
      covariance_ok = a$covariance_ok,
      numerical_checks_ok =
        isTRUE(opt$convergence == 0L) && isTRUE(a$pass),
      note = "post-optimizer assessment",
      stringsAsFactors = FALSE
    )
    csv(history, file.path(out, "optimization_history.csv"))

    # Bounded local Newton refinement if needed and covariance is usable.
    line_hist <- data.frame()
    par_current <- opt$par
    objective_current <- a$objective
    assessment <- a
    custom_steps <- 0L
    custom_ok <- isTRUE(opt$convergence == 0L) && isTRUE(a$pass)

    if (!custom_ok &&
        isTRUE(a$pdHess) &&
        isTRUE(a$covariance_ok) &&
        is.finite(a$scaled_gradient)) {

      for (step in seq_len(cfg$max_newton_steps)) {
        phase(
          cfg, "NEWTON_REFINEMENT", group,
          paste("Step", step, "of at most", cfg$max_newton_steps,
                "| scaled gradient", signif(assessment$scaled_gradient, 6))
        )

        trial <- line_search(
          par = par_current,
          objective = objective_current,
          gradient = assessment$gradient,
          covariance = assessment$covariance,
          fn = function(p) {
            withCallingHandlers(obj$fn(p), warning = catcher)
          },
          max_backtracks = cfg$max_backtracks
        )

        h <- cbind(step = step, trial$history)
        line_hist <- rbind(line_hist, h)
        csv(line_hist, file.path(out, "line_search_history.csv"))

        if (!trial$accepted) {
          obj$fn(par_current)
          break
        }

        par_current <- trial$par
        objective_current <- trial$objective
        custom_steps <- step

        # Save accepted arrays BEFORE expensive new Hessian/covariance.
        accepted_arrays <- extract_arrays(par_current)
        atomic_rds(
          list(
            signature = cfg$signature,
            data_signature = input$data_signature,
            step = step,
            par = par_current,
            parameters = accepted_arrays,
            objective = objective_current,
            alpha = trial$alpha
          ),
          file.path(out, sprintf("accepted_newton_step_%02d.rds", step))
        )

        phase(
          cfg, "UNCERTAINTY", group,
          paste("Accepted Newton step", step,
                "| recomputing fresh marginal covariance")
        )

        assessment <- assess(obj, par_current)
        warnings <- unique(c(warnings, assessment$warnings))

        row <- data.frame(
          step = step,
          objective = assessment$objective,
          max_abs_gradient = max(abs(assessment$gradient)),
          gradient_covariance_norm = assessment$scaled_gradient,
          pdHess = assessment$pdHess,
          covariance_ok = assessment$covariance_ok,
          numerical_checks_ok = isTRUE(assessment$pass),
          note = paste("accepted damped Newton; alpha", trial$alpha),
          stringsAsFactors = FALSE
        )
        history <- rbind(history, row)
        csv(history, file.path(out, "optimization_history.csv"))

        if (isTRUE(assessment$pass)) {
          custom_ok <- TRUE
          break
        }

        if (!isTRUE(assessment$pdHess) ||
            !isTRUE(assessment$covariance_ok)) {
          break
        }
      }
    }

    # Final point.
    final_par <- if (custom_steps > 0L) par_current else opt$par
    final_assessment <- if (custom_steps > 0L) assessment else a

    # Criterion after local refinement: if optimizer had non-zero code but
    # Newton refinement reached the documented gradient/Hessian criterion,
    # record convergence_code 0 for THIS bounded final solver and preserve the
    # original optimizer code separately.
    numerical_ok <- isTRUE(final_assessment$pass) &&
      is.finite(final_assessment$objective)

    obj$fn(final_par)
    pars <- extract_arrays(final_par)
    report <- obj$report(obj$env$last.par)

    ids <- which(names(final_par) == "beta")
    need(length(ids) == ncol(input$data$X),
         "Final beta mapping failed.")

    beta <- stats::setNames(
      as.numeric(final_par[ids]),
      colnames(input$data$X)
    )
    Vbeta <- final_assessment$covariance[ids, ids, drop = FALSE]
    dimnames(Vbeta) <- list(names(beta), names(beta))

    extent <- max(
      diff(range(input$sites$x_km)),
      diff(range(input$sites$y_km))
    )

    checks <- data.frame(
      convergence_code = if (numerical_ok) 0L else opt$convergence,
      source_optimizer_convergence = opt$convergence,
      accepted_newton_steps = custom_steps,
      pdHess = final_assessment$pdHess,
      max_abs_gradient = max(abs(final_assessment$gradient)),
      gradient_covariance_norm = final_assessment$scaled_gradient,
      numerical_checks_ok = numerical_ok,
      logLik = -final_assessment$objective,
      parameter_df = length(final_par),
      AIC = 2 * final_assessment$objective + 2 * length(final_par),
      nobs = nrow(input$frame),
      spatial_range_km = as.numeric(report$range),
      spatial_sd = as.numeric(report$spatial_sd),
      range_greater_than_1_5_extent =
        as.numeric(report$range) > 1.5 * extent,
      spatial_sd_below_0_001 =
        as.numeric(report$spatial_sd) < .001,
      AR1_rho = as.numeric(report$rho),
      Gamma_shape = as.numeric(report$shape),
      stringsAsFactors = FALSE
    )

    fit <- list(
      signature = cfg$signature,
      data_signature = input$data_signature,
      parameters = pars,
      optimizer = list(
        par = final_par,
        objective = final_assessment$objective,
        convergence = checks$convergence_code,
        original_optimizer = opt,
        accepted_newton_steps = custom_steps,
        solver = if (custom_steps > 0L)
          "nlminb/BFGS + bounded damped Newton"
        else
          "nlminb/BFGS"
      ),
      beta = beta,
      V = Vbeta,
      cov_fixed = final_assessment$covariance,
      checks = checks,
      eta = as.numeric(report$eta),
      field_at_sites = report$spatial,
      group = group,
      variant = "full",
      fit_minutes = as.numeric(difftime(Sys.time(), started, units = "mins")),
      warnings = unique(warnings),
      onset_predictor_preserved = TRUE,
      flight_length_predictor_spatial = TRUE,
      residual_diagnostics = "PENDING"
    )

    # Save final fit BEFORE exports.
    atomic_rds(fit, final_file)
    export_fit(fit, input, cfg, group)

    status <- if (numerical_ok) "FIT_OK" else "NUMERICAL_REVIEW_REQUIRED"

    atomic_rds(
      list(
        signature = cfg$signature,
        status = status,
        data_signature = input$data_signature,
        time = Sys.time()
      ),
      done_file
    )

    phase(
      cfg,
      status,
      group,
      paste(
        "scaled gradient",
        signif(checks$gradient_covariance_norm, 6),
        "| Newton steps", custom_steps
      )
    )

    status
  }

  collect <- function(cfg) {
    for (nm in c(
      "fit_checks",
      "fixed_effects",
      "coefficient_comparison",
      "focal_interactions"
    )) {
      rows <- lapply(GROUPS, function(g) {
        p <- file.path(
          cfg$out,
          paste0(g, "__full"),
          paste0(nm, ".csv")
        )
        if (!file.exists(p)) return(NULL)
        cbind(
          group = g,
          utils::read.csv(
            p,
            stringsAsFactors = FALSE,
            check.names = FALSE
          )
        )
      })
      rows <- Filter(Negate(is.null), rows)
      if (length(rows)) {
        csv(
          do.call(rbind, rows),
          file.path(cfg$out, paste0("all_", nm, ".csv"))
        )
      }
    }

    if (requireNamespace("zip", quietly = TRUE)) {
      tryCatch({
        rel <- list.files(
          cfg$out,
          recursive = TRUE,
          full.names = FALSE
        )
        rel <- rel[
          grepl("\\.(csv|txt|log|R|cpp)$", rel, ignore.case = TRUE)
        ]
        zip::zipr(
          file.path(cfg$out, "results_to_review.zip"),
          files = rel,
          root = cfg$out,
          mode = "mirror",
          include_directories = FALSE
        )
      }, error = function(e)
        message("Optional review ZIP failed: ", conditionMessage(e)))
    }

    invisible(NULL)
  }

  invoke_worker <- function(cfg, fun, group = "") {
    log <- file.path(
      cfg$out, "logs",
      paste0(
        fun, "_", group, "_",
        format(Sys.time(), "%Y%m%d_%H%M%S"),
        ".log"
      )
    )

    p <- callr::r_bg(
      function(engine, cfg, fun, group) {
        e <- new.env(parent = globalenv())
        sys.source(engine, envir = e)
        args <- if (fun == "fit_worker") {
          list(cfg = cfg, group = group)
        } else {
          list(cfg = cfg)
        }
        do.call(e[[fun]], args)
      },
      args = list(
        engine = cfg$engine,
        cfg = cfg,
        fun = fun,
        group = group
      ),
      libpath = cfg$libpath,
      stdout = log,
      stderr = "2>&1",
      supervise = TRUE,
      user_profile = FALSE,
      system_profile = FALSE,
      wd = cfg$out,
      env = c(
        callr::rcmd_safe_env(),
        OMP_NUM_THREADS = "1",
        OPENBLAS_NUM_THREADS = "1",
        MKL_NUM_THREADS = "1"
      )
    )

    on.exit(if (p$is_alive()) p$kill(), add = TRUE)

    atomic_rds(
      list(
        pid = p$get_pid(),
        fun = fun,
        group = group,
        log = log
      ),
      file.path(cfg$out, "active_worker.rds")
    )

    p$wait()
    p$get_result()
  }

  controller <- function(cfg) {
    on.exit(unlink(cfg$lock, recursive = TRUE), add = TRUE)

    status <- data.frame(
      task = c("compile", paste0("abundance_", GROUPS)),
      status = "PENDING",
      message = "",
      stringsAsFactors = FALSE
    )
    save_status <- function()
      csv(status, file.path(cfg$out, "workflow_status.csv"))
    save_status()

    current <- NA_integer_

    tryCatch({
      check_packages()
      self_test()
      disk_check(cfg$out, cfg$minimum_free_GiB)

      phase(
        cfg, "PREFLIGHT_COMPLETE", "",
        "Onset+flight-length inputs, preserved abundance structure and numerical self-tests validated"
      )

      current <- 1L
      status$status[1] <- "RUNNING"
      save_status()
      status$status[1] <- invoke_worker(cfg, "compile_worker")
      save_status()

      for (i in seq_along(GROUPS)) {
        row <- i + 1L
        current <- row
        status$status[row] <- "RUNNING"
        save_status()

        disk_check(cfg$out, cfg$minimum_free_GiB)

        ans <- tryCatch(
          invoke_worker(cfg, "fit_worker", GROUPS[i]),
          error = function(e) {
            status$message[row] <<- conditionMessage(e)
            "ERROR"
          }
        )

        status$status[row] <- ans
        save_status()
        collect(cfg)
      }

      ok <- all(
        status$status[match(paste0("abundance_", GROUPS), status$task)] ==
          "FIT_OK"
      )

      phase(
        cfg,
        if (ok) "FITS_COMPLETE" else "FINISHED_WITH_ISSUES",
        "",
        if (ok)
          "Both onset+flight-length abundance fits are numerically complete; review conditional flight-length effects"
        else
          "Inspect workflow_status.csv and per-group numerical history",
        terminal = TRUE
      )

      invisible(status)

    }, error = function(e) {
      if (!is.na(current)) {
        status$status[current] <- "ERROR"
        status$message[current] <- conditionMessage(e)
      }
      save_status()
      writeLines(
        conditionMessage(e),
        file.path(cfg$out, "controller_error.txt")
      )
      try(
        phase(
          cfg, "STOPPED_WITH_ERROR", "",
          conditionMessage(e),
          terminal = TRUE
        ),
        silent = TRUE
      )
      invisible(status)

    }, finally = {
      save_status()
      try(collect(cfg), silent = TRUE)
    })
  }

  dump_engine <- function(out) {
    env <- environment(dump_engine)
    n <- ls(env, all.names = TRUE)
    fun <- n[vapply(
      n,
      function(x) is.function(get(x, envir = env)),
      logical(1)
    )]

    path <- file.path(out, "workflow_engine.R")
    dump(
      c(
        "VERSION", "ROOT", "SOURCE_PREP", "SOURCE_ABUNDANCE_RUN",
        "FLIGHT_LENGTH_FILE", "OUTPUT_ROOT", "GROUPS",
        "CPP_MD5", "GRADIENT_TOL", "CORE", fun
      ),
      file = path,
      envir = env,
      control = "all"
    )
    normalizePath(path, winslash = "/", mustWork = TRUE)
  }

  get_job <- function() {
    j <- state$job
    if (is.null(j)) {
      j <- getOption("phenoimpact.onset_flight_abundance.job")
    }
    j
  }

  status <- function(job = get_job()) {
    need(!is.null(job), "No onset+flight-length abundance job in this R session.")

    alive <- job$process$is_alive()
    message(
      "Elapsed: ",
      round(
        as.numeric(difftime(Sys.time(), job$started, units = "hours")),
        2
      ),
      " h | controller running: ", alive
    )

    p <- file.path(job$out, "progress.rds")
    if (file.exists(p)) {
      x <- tryCatch(readRDS(p), error = function(e) NULL)
      if (!is.null(x)) {
        message(x$stage, " | ", x$group, " | ", x$detail)
      }
    }

    f <- file.path(job$out, "workflow_status.csv")
    if (file.exists(f)) {
      try(
        print(
          utils::read.csv(f, stringsAsFactors = FALSE),
          row.names = FALSE
        ),
        silent = TRUE
      )
    }

    message("Output: ", job$out)

    invisible(list(
      alive = alive,
      output = job$out
    ))
  }

  watch <- function(job = get_job(), every = 60) {
    need(!is.null(job), "No onset+flight-length abundance job.")
    tryCatch({
      repeat {
        status(job)
        if (!job$process$is_alive()) break
        Sys.sleep(every)
      }
    }, interrupt = function(e)
      message(
        "Monitoring paused; the background queue continues. ",
        "Use pheno_onset_flight_abundance$watch()."
      ))
    invisible(job)
  }

  show_log <- function(n = 40L, job = get_job()) {
    need(!is.null(job), "No onset+flight-length abundance job.")
    a <- file.path(job$out, "active_worker.rds")
    w <- if (file.exists(a)) tryCatch(readRDS(a), error = function(e) NULL) else NULL
    log <- if (!is.null(w) && file.exists(w$log)) w$log else job$log
    message("Log: ", log)
    cat(utils::tail(readLines(log, warn = FALSE), n), sep = "\n")
    invisible(log)
  }

  run <- function(
    source_prep = SOURCE_PREP,
    source_abundance_run = SOURCE_ABUNDANCE_RUN,
    flight_length_file = FLIGHT_LENGTH_FILE,
    minimum_free_GiB = 50,
    optimizer_iterations = NULL,
    optimizer_evaluations = NULL,
    max_newton_steps = 2L,
    max_backtracks = 6L,
    monitor = FALSE
  ) {
    check_packages()
    self_test()

    source_prep <- normalizePath(
      source_prep, winslash = "/", mustWork = TRUE
    )
    source_abundance_run <- normalizePath(
      source_abundance_run, winslash = "/", mustWork = TRUE
    )
    flight_length_file <- normalizePath(
      flight_length_file, winslash = "/", mustWork = TRUE
    )

    required <- c(
      file.path(source_prep, paste0("input_", GROUPS, ".rds")),
      file.path(source_prep, "mesh.rds"),
      file.path(
        source_abundance_run,
        paste0(GROUPS, "__full"),
        "fit.rds"
      ),
      file.path(source_abundance_run, "run_config.rds"),
      file.path(source_abundance_run, "compiled", "trend_spde.cpp"),
      flight_length_file
    )
    need(
      all(file.exists(required)),
      paste(
        "Required source file missing:\n",
        paste(required[!file.exists(required)], collapse = "\n")
      )
    )
    need(
      identical(
        unname(tools::md5sum(
          file.path(source_abundance_run, "compiled", "trend_spde.cpp")
        )),
        CPP_MD5
      ),
      "Reviewed abundance C++ template has changed."
    )

    src_cfg <- readRDS(
      file.path(source_abundance_run, "run_config.rds")
    )
    if (is.null(optimizer_iterations)) {
      optimizer_iterations <- src_cfg$optimizer_iterations
    }
    if (is.null(optimizer_evaluations)) {
      optimizer_evaluations <- src_cfg$optimizer_evaluations
    }

    need(
      is.numeric(optimizer_iterations) &&
        length(optimizer_iterations) == 1L &&
        is.finite(optimizer_iterations) &&
        optimizer_iterations > 0,
      "Invalid optimizer_iterations."
    )
    need(
      is.numeric(optimizer_evaluations) &&
        length(optimizer_evaluations) == 1L &&
        is.finite(optimizer_evaluations) &&
        optimizer_evaluations > 0,
      "Invalid optimizer_evaluations."
    )
    need(
      is.numeric(max_newton_steps) &&
        length(max_newton_steps) == 1L &&
        max_newton_steps >= 0 &&
        max_newton_steps <= 3 &&
        max_newton_steps == floor(max_newton_steps),
      "max_newton_steps must be 0..3."
    )
    need(
      is.numeric(max_backtracks) &&
        length(max_backtracks) == 1L &&
        max_backtracks >= 0 &&
        max_backtracks <= 8 &&
        max_backtracks == floor(max_backtracks),
      "max_backtracks must be 0..8."
    )

    current <- get_job()
    if (!is.null(current) && current$process$is_alive()) {
      message("Onset+flight-length abundance queue already running.")
      return(invisible(current))
    }

    disk_check(source_prep, minimum_free_GiB)

    versions <- package_versions()
    hashes <- unname(tools::md5sum(required))

    signature <- digest::digest(
      list(
        version = VERSION,
        files = normalizePath(required, winslash = "/", mustWork = TRUE),
        hashes = hashes,
        packages = versions,
        formula =
          "~ year_decade * onset_plasticity_z + year_decade * flight_length_plasticity_z",
        optimizer_iterations = optimizer_iterations,
        optimizer_evaluations = optimizer_evaluations,
        max_newton_steps = max_newton_steps,
        max_backtracks = max_backtracks,
        gradient_tolerance = GRADIENT_TOL
      ),
      algo = "sha256"
    )

    dir.create(OUTPUT_ROOT, recursive = TRUE, showWarnings = FALSE)
    out <- file.path(
      OUTPUT_ROOT,
      paste0("run_", substr(signature, 1, 12))
    )
    dir.create(file.path(out, "logs"), recursive = TRUE, showWarnings = FALSE)
    out <- normalizePath(out, winslash = "/", mustWork = TRUE)

    lock <- file.path(OUTPUT_ROOT, "queue.lock")
    if (dir.exists(lock)) {
      owner_file <- file.path(lock, "owner.rds")
      owner <- if (file.exists(owner_file)) {
        tryCatch(readRDS(owner_file), error = function(e) NULL)
      } else NULL
      alive <- if (is.null(owner)) FALSE else tryCatch(
        ps::ps_is_running(
          ps::ps_handle(owner$pid, time = owner$create_time)
        ),
        error = function(e) FALSE
      )
      need(
        !alive,
        paste("Another onset+flight-length abundance queue may be running:", lock)
      )
      unlink(lock, recursive = TRUE)
    }
    need(dir.create(lock), "Cannot acquire abundance queue lock.")

    launched <- FALSE
    on.exit(if (!launched) unlink(lock, recursive = TRUE), add = TRUE)

    cfg <- list(
      version = VERSION,
      signature = signature,
      root = ROOT,
      source_prep = source_prep,
      source_abundance_run = source_abundance_run,
      flight_length_file = flight_length_file,
      prep = file.path(out, "prepared_inputs"),
      out = out,
      lock = lock,
      libpath = .libPaths(),
      minimum_free_GiB = minimum_free_GiB,
      optimizer_iterations = as.integer(optimizer_iterations),
      optimizer_evaluations = as.integer(optimizer_evaluations),
      max_newton_steps = as.integer(max_newton_steps),
      max_backtracks = as.integer(max_backtracks)
    )

    # Lightweight, deterministic input replacement happens before the background
    # optimizer queue. It neither evaluates nor fits the abundance likelihood.
    prep_done <- file.path(cfg$prep, "preparation_complete.rds")
    if (file.exists(prep_done)) {
      done <- readRDS(prep_done)
      need(
        identical(done$signature, signature),
        "Prepared onset+flight-length input cache mismatch."
      )
      message("Reusing prepared onset+flight-length abundance inputs.")
    } else {
      prepare_onset_flight_inputs(cfg)
    }

    cfg$engine <- dump_engine(out)
    atomic_rds(cfg, file.path(out, "run_config.rds"))
    csv(
      data.frame(
        file = normalizePath(required, winslash = "/", mustWork = TRUE),
        md5 = hashes
      ),
      file.path(out, "source_manifest.csv")
    )
    csv(versions, file.path(out, "package_versions.csv"))

    writeLines(c(
      "SPATIAL ABUNDANCE MODELS: ONSET + DIRECT FLIGHT-LENGTH PLASTICITY",
      "",
      paste("Workflow:", VERSION),
      paste("Signature:", signature),
      paste("Final abundance source:", source_abundance_run),
      paste("Flight-length plasticity source:", flight_length_file),
      "",
      "Question: at equal onset plasticity, do populations with greater direct",
      "flight-length extension in warm years show more favourable trends?",
      "",
      "The response rows, year scale, species/site indices, population AR1 links",
      "and annual 15-km spatial field structure are unchanged.",
      "ONSET plasticity is preserved exactly.",
      "The old OFFSET predictor is replaced by spatial direct flight-length",
      "plasticity standardized over all accepted phenology populations within",
      "each voltinism group before the abundance join.",
      "",
      "Gamma/log likelihood; species/site intercept+slope; stationary population",
      "AR1 with true calendar gaps; annual IID Matern field.",
      "All parameters are re-estimated from reviewed warm starts.",
      "",
      "Primary test: year_decade:flight_length_plasticity_z.",
      "It is conditional on onset plasticity because both are in the same model.",
      "No first-stage plasticity uncertainty is propagated.",
      "No reduced LRTs or residual diagnostics are run before focal review."
    ), file.path(out, "README.txt"))

    log <- file.path(
      out, "logs",
      paste0("controller_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log")
    )

    proc <- callr::r_bg(
      function(engine, cfg) {
        e <- new.env(parent = globalenv())
        sys.source(engine, envir = e)
        e$controller(cfg)
      },
      args = list(engine = cfg$engine, cfg = cfg),
      libpath = cfg$libpath,
      stdout = log,
      stderr = "2>&1",
      supervise = TRUE,
      user_profile = FALSE,
      system_profile = FALSE,
      wd = out,
      env = c(
        callr::rcmd_safe_env(),
        OMP_NUM_THREADS = "1",
        OPENBLAS_NUM_THREADS = "1",
        MKL_NUM_THREADS = "1"
      )
    )

    h <- ps::ps_handle(proc$get_pid())
    atomic_rds(
      list(
        pid = proc$get_pid(),
        create_time = ps::ps_create_time(h)
      ),
      file.path(lock, "owner.rds")
    )

    job <- list(
      process = proc,
      out = out,
      started = Sys.time(),
      log = log
    )
    state$job <- job
    options(phenoimpact.onset_flight_abundance.job = job)
    launched <- TRUE

    message("ONSET + FLIGHT-LENGTH abundance queue launched.")
    message("Output: ", out)
    message(
      "Keep R open. Check with pheno_onset_flight_abundance$status()."
    )

    if (monitor) watch(job)
    invisible(job)
  }

  stop_job <- function(job = get_job()) {
    need(!is.null(job), "No onset+flight-length abundance job.")
    if (job$process$is_alive()) job$process$kill()
    message(
      "Onset+flight-length abundance queue stopped. ",
      "Saved portable checkpoints/fits remain."
    )
    invisible(job)
  }

  list(
    run = run,
    status = status,
    watch = watch,
    log = show_log,
    stop = stop_job,
    self_test = self_test
  )
})

message("Functions loaded. Start with:")
message("  onset_flight_job <- pheno_onset_flight_abundance$run(monitor = FALSE)")
message("  pheno_onset_flight_abundance$status()")
