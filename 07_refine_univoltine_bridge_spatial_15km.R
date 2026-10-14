# =============================================================================
# phenoIMPACT | Refine/adjudicate ONLY the UNIVOLTINE bridge model
# Version 1.0 | 2026-10-08
#
# SOURCE MODEL
#   bridge_geography_plasticity_spatial_15km/
#   run_20261008_100948_5724/univoltine__bridge/fit.rds
#
# WHY
#   The source bridge fit has:
#     nlminb/BFGS convergence code = 1
#     positive-definite Hessian     = TRUE
#     covariance-scaled gradient   = 0.034 < 0.05
#     finite fixed-effect covariance
#
#   Therefore this script does NOT restart the full fit. It reconstructs the
#   exact saved TMB likelihood and:
#
#     1) verifies the saved objective, coefficients and fitted eta;
#     2) recalculates a fresh marginal Hessian/covariance at the exact source
#        parameter point;
#     3) declares numerical completion there if the prespecified criterion
#        remains satisfied;
#     4) otherwise allows at most TWO damped Newton updates, with backtracking,
#        using the full marginal covariance:
#
#           direction = -V %*% gradient
#
#   Every accepted parameter point is checkpointed BEFORE recalculating
#   uncertainty. The multivoltine bridge fit and all source files are READ ONLY.
#
# PASS CRITERION
#   - positive Hessian
#   - positive-definite finite marginal covariance
#   - covariance-scaled gradient < 0.05
#   - finite objective no worse than the source objective
#
# A convergence code of 0 in the refined output means the documented custom
# criterion above was met. It is NOT a new nlminb/BFGS return code. The source
# optimizer and its code are preserved in the output.
#
# RUN
#   source(file.choose(), encoding = "UTF-8")
#   bridge_refinement_job <- pheno_bridge_refine$run(monitor = FALSE)
#
# MONITOR
#   pheno_bridge_refine$status()
#
# If R or Windows stops, run pheno_bridge_refine$run() again with the same
# settings. The deterministic output folder and saved checkpoint are reused.
# Keep R open and do not run another large model simultaneously.
# =============================================================================


pheno_bridge_refine <- local({

  VERSION <- "univoltine_bridge_refinement_1.0"
  SOURCE_BRIDGE_VERSION <- "bridge_geography_plasticity_spatial_15km_1.0"

  DEFAULT_BRIDGE_RUN <- paste0(
    "E:/phenoIMPACT project/code/phenoIMPACT/output/population_trends/",
    "bridge_geography_plasticity_spatial_15km/",
    "run_20261008_100948_5724"
  )

  CORE_VERSIONS <- c(
    TMB = "1.9.17",
    Matrix = "1.7.3",
    fmesher = "0.5.0"
  )

  GRADIENT_TOL <- 0.05

  state <- new.env(parent = emptyenv())


  # ---------------------------------------------------------------------------
  # Utilities
  # ---------------------------------------------------------------------------

  need <- function(ok, msg) {
    if (!isTRUE(ok)) stop(msg, call. = FALSE)
  }


  same_num <- function(x, y, tol = 1e-8) {

    is.numeric(x) &&
      is.numeric(y) &&
      length(x) == length(y) &&
      all(is.finite(x)) &&
      all(is.finite(y)) &&
      (
        length(x) == 0L ||
          max(abs(as.numeric(x) - as.numeric(y))) <= tol
      )
  }


  csv <- function(x, path) {
    utils::write.csv(x, path, row.names = FALSE, na = "")
  }


  atomic_save <- function(x, path) {

    dir.create(
      dirname(path),
      recursive = TRUE,
      showWarnings = FALSE
    )

    tmp <- tempfile(
      ".writing_",
      tmpdir = dirname(path),
      fileext = ".rds"
    )

    on.exit(unlink(tmp), add = TRUE)

    saveRDS(
      x,
      tmp,
      compress = "gzip"
    )

    backup <- paste0(path, ".previous")

    if (file.exists(backup)) {
      need(
        unlink(backup) == 0L,
        paste("Cannot remove old backup:", backup)
      )
    }

    had <- file.exists(path)

    if (had) {
      need(
        file.rename(path, backup),
        paste("Cannot retain previous checkpoint:", path)
      )
    }

    if (!file.rename(tmp, path)) {

      if (had && file.exists(backup)) {
        file.rename(backup, path)
      }

      stop(
        "Cannot finalize checkpoint: ",
        path,
        call. = FALSE
      )
    }

    invisible(path)
  }


  phase <- function(
    cfg,
    stage,
    detail = "",
    terminal = FALSE
  ) {

    now <- Sys.time()

    atomic_save(
      list(
        stage = stage,
        detail = detail,
        time = now,
        terminal = terminal,
        started = cfg$launched_at,
        finished = if (terminal) now else NULL
      ),
      file.path(cfg$out, "progress.rds")
    )

    message(
      "[",
      format(now, "%Y-%m-%d %H:%M:%S"),
      "] ",
      stage,
      if (nzchar(detail)) paste0(" | ", detail) else ""
    )
  }


  check_packages <- function() {

    packages <- c(
      names(CORE_VERSIONS),
      "callr",
      "ps",
      "digest"
    )

    missing <- packages[
      !vapply(
        packages,
        requireNamespace,
        logical(1),
        quietly = TRUE
      )
    ]

    need(
      !length(missing),
      paste(
        "Missing packages:",
        paste(missing, collapse = ", "),
        "\nNo packages are installed automatically."
      )
    )

    current <- vapply(
      names(CORE_VERSIONS),
      function(p) {
        as.character(utils::packageVersion(p))
      },
      character(1)
    )

    need(
      identical(current, CORE_VERSIONS),
      paste(
        "Core package versions must match the reviewed fit:",
        paste(
          names(CORE_VERSIONS),
          CORE_VERSIONS,
          collapse = "; "
        )
      )
    )
  }


  disk_check <- function(
    path,
    minimum_GiB,
    verbose = TRUE
  ) {

    d <- ps::ps_disk_usage(path)

    need(
      nrow(d) == 1L &&
        is.finite(d$available),
      "Cannot check free disk space."
    )

    free <- as.numeric(d$available) / 1024^3

    if (verbose) {
      message(
        "Disk preflight: ",
        round(free, 1),
        " GiB available."
      )
    }

    need(
      free >= minimum_GiB,
      paste(
        "Only",
        round(free, 1),
        "GiB free; require at least",
        minimum_GiB,
        "GiB."
      )
    )

    invisible(free)
  }


  checked_covariance <- function(V, n) {

    need(
      is.matrix(V) &&
        identical(dim(V), as.integer(c(n, n))) &&
        all(is.finite(V)),
      "Marginal parameter covariance is missing or non-finite."
    )

    symmetry <-
      max(abs(V - t(V))) /
      max(1, max(abs(V)))

    need(
      symmetry < 1e-8,
      "Marginal covariance is not numerically symmetric."
    )

    V <- (V + t(V)) / 2

    need(
      all(diag(V) > 0),
      "Marginal covariance has a non-positive diagonal."
    )

    ch <- tryCatch(
      chol(V),
      error = function(e) NULL
    )

    need(
      !is.null(ch),
      paste0(
        "Marginal covariance is not positive definite; ",
        "no ridge is added."
      )
    )

    V
  }


  gradient_norm <- function(g, V) {

    need(
      is.numeric(g) &&
        all(is.finite(g)),
      "Gradient contains non-finite values."
    )

    V <- checked_covariance(
      V,
      length(g)
    )

    value <- as.numeric(
      crossprod(
        g,
        V %*% g
      )
    )

    need(
      is.finite(value) &&
        value >= 0,
      "Invalid covariance-scaled gradient."
    )

    sqrt(value)
  }


  newton_direction <- function(g, V) {

    V <- checked_covariance(
      V,
      length(g)
    )

    norm <- gradient_norm(g, V)

    # Trust safeguard: at most one local covariance-metric unit.
    shrink <- if (norm > 1) {
      1 / norm
    } else {
      1
    }

    direction <-
      -as.numeric(V %*% g) *
      shrink

    dot <- sum(g * direction)

    need(
      all(is.finite(direction)) &&
        is.finite(dot) &&
        dot < 0,
      paste0(
        "No finite descent direction; ",
        "optimization is stopped for review."
      )
    )

    list(
      direction = direction,
      g_dot_direction = dot,
      shrink = shrink
    )
  }


  backtrack <- function(
    par,
    value,
    gradient,
    covariance,
    evaluate,
    max_backtracks = 6L
  ) {

    step <- newton_direction(
      gradient,
      covariance
    )

    records <- vector(
      "list",
      max_backtracks + 1L
    )

    for (j in 0:max_backtracks) {

      alpha <- 0.5^j

      trial <- par +
        alpha * step$direction

      trial_value <- evaluate(
        trial,
        j,
        alpha
      )

      armijo_bound <-
        value +
        1e-4 *
        alpha *
        step$g_dot_direction

      accepted <-
        is.finite(trial_value) &&
        trial_value < value &&
        trial_value <= armijo_bound

      records[[j + 1L]] <- data.frame(
        backtrack = j,
        alpha = alpha,
        direction_shrink = step$shrink,
        objective_before = value,
        objective_trial = trial_value,
        Armijo_bound = armijo_bound,
        accepted = accepted
      )

      if (accepted) {

        return(
          list(
            accepted = TRUE,
            par = trial,
            objective = trial_value,
            alpha = alpha,
            history = do.call(
              rbind,
              records[seq_len(j + 1L)]
            )
          )
        )
      }
    }

    list(
      accepted = FALSE,
      par = par,
      objective = value,
      alpha = 0,
      history = do.call(rbind, records)
    )
  }


  self_test <- function() {

    H <- matrix(
      c(4, 1, 1, 3),
      2,
      2
    )

    V <- solve(H)

    par <- c(a = 0.02, b = -0.03)

    fn <- function(x) {
      as.numeric(
        crossprod(
          x,
          H %*% x
        )
      ) / 2
    }

    gr <- function(x) {
      as.numeric(H %*% x)
    }

    one <- backtrack(
      par,
      fn(par),
      gr(par),
      V,
      function(x, j, a) fn(x)
    )

    need(
      one$accepted &&
        one$alpha == 1 &&
        max(abs(one$par)) < 1e-12,
      "Newton self-test failed."
    )

    larger <- backtrack(
      par,
      fn(par),
      gr(par),
      8 * V,
      function(x, j, a) fn(x)
    )

    need(
      larger$accepted &&
        larger$alpha < 1 &&
        larger$objective < fn(par),
      "Backtracking self-test failed."
    )

    reject <- backtrack(
      par,
      fn(par),
      gr(par),
      V,
      function(x, j, a) Inf,
      2L
    )

    need(
      !reject$accepted &&
        identical(reject$par, par) &&
        nrow(reject$history) == 3,
      "Rejected-step self-test failed."
    )

    invisible(TRUE)
  }


  # ---------------------------------------------------------------------------
  # Exact bridge objective reconstruction
  # ---------------------------------------------------------------------------

  sparse_TMB <- function(x) {

    methods::as(
      methods::as(
        methods::as(x, "dMatrix"),
        "generalMatrix"
      ),
      "TsparseMatrix"
    )
  }


  make_objective <- function(
    data,
    parameters
  ) {

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
      map = NULL,
      DLL = "trend_spde",
      silent = TRUE,
      inner.control = list(
        maxit = 1000,
        trace = FALSE
      )
    )
  }


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
      paste("Source DLL is missing:", dll)
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


  fixed_table <- function(beta, V) {

    diagV <- diag(V)

    se <- ifelse(
      is.finite(diagV) &
        diagV > 0,
      sqrt(pmax(0, diagV)),
      NA_real_
    )

    z <- beta / se

    data.frame(
      term = names(beta),
      estimate = as.numeric(beta),
      std.error = as.numeric(se),
      conf.low = as.numeric(beta - 1.96 * se),
      conf.high = as.numeric(beta + 1.96 * se),
      z = as.numeric(z),
      p_Wald = as.numeric(
        2 * stats::pnorm(-abs(z))
      ),
      stringsAsFactors = FALSE
    )
  }


  linear_contrast <- function(
    beta,
    V,
    weights,
    label,
    unit,
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
      unit = unit,
      estimate = estimate,
      std.error = se,
      conf.low = estimate - 1.96 * se,
      conf.high = estimate + 1.96 * se,
      p_Wald = 2 * stats::pnorm(
        -abs(estimate / se)
      ),
      ratio_of_decadal_factors = exp(estimate),
      ratio_conf_low = exp(estimate - 1.96 * se),
      ratio_conf_high = exp(estimate + 1.96 * se),
      percent_difference = 100 * expm1(estimate),
      percent_low = 100 * expm1(
        estimate - 1.96 * se
      ),
      percent_high = 100 * expm1(
        estimate + 1.96 * se
      ),
      interpretation = interpretation,
      stringsAsFactors = FALSE
    )
  }


  # ---------------------------------------------------------------------------
  # Deterministic configuration and source manifest
  # ---------------------------------------------------------------------------

  configure <- function(
    bridge_run,
    max_steps,
    max_backtracks,
    minimum_free_GiB
  ) {

    bridge_run <- normalizePath(
      bridge_run,
      winslash = "/",
      mustWork = TRUE
    )

    bridge_cfg_file <- file.path(
      bridge_run,
      "run_config.rds"
    )

    prepared_file <- file.path(
      bridge_run,
      "prepared_bridge_univoltine.rds"
    )

    source_fit_file <- file.path(
      bridge_run,
      "univoltine__bridge",
      "fit.rds"
    )

    source_optimized_file <- file.path(
      bridge_run,
      "univoltine__bridge",
      "optimized_parameters.rds"
    )

    required <- c(
      bridge_cfg_file,
      prepared_file,
      source_fit_file,
      source_optimized_file
    )

    need(
      all(file.exists(required)),
      paste(
        "Missing bridge source file(s):",
        paste(
          required[!file.exists(required)],
          collapse = "\n"
        )
      )
    )

    bridge_cfg <- readRDS(bridge_cfg_file)

    need(
      identical(
        bridge_cfg$version,
        SOURCE_BRIDGE_VERSION
      ),
      "Unexpected bridge workflow version."
    )

    source_run <- normalizePath(
      bridge_cfg$source_run,
      winslash = "/",
      mustWork = TRUE
    )

    cpp <- file.path(
      source_run,
      "compiled",
      "trend_spde.cpp"
    )

    dll <- TMB::dynlib(
      file.path(
        source_run,
        "compiled",
        "trend_spde"
      )
    )

    mesh <- file.path(
      source_run,
      "mesh.rds"
    )

    extra <- c(cpp, dll, mesh)

    need(
      all(file.exists(extra)),
      paste(
        "Missing likelihood/mesh source file(s):",
        paste(
          extra[!file.exists(extra)],
          collapse = "\n"
        )
      )
    )

    manifest_paths <- c(
      bridge_cfg_file,
      prepared_file,
      source_fit_file,
      source_optimized_file,
      cpp,
      dll,
      mesh
    )

    message(
      "Fingerprinting ",
      length(manifest_paths),
      " source files..."
    )

    manifest <- data.frame(
      role = c(
        "bridge_run_config",
        "prepared_bridge_univoltine",
        "source_bridge_fit",
        "source_bridge_optimized_parameters",
        "likelihood_CPP",
        "likelihood_DLL",
        "source_mesh"
      ),
      path = normalizePath(
        manifest_paths,
        winslash = "/",
        mustWork = TRUE
      ),
      bytes = file.info(manifest_paths)$size,
      md5 = unname(tools::md5sum(manifest_paths)),
      stringsAsFactors = FALSE
    )

    need(
      !anyNA(manifest$md5),
      "Could not fingerprint all source files."
    )

    signature <- digest::digest(
      list(
        version = VERSION,
        manifest = manifest,
        core_versions = CORE_VERSIONS,
        max_steps = as.integer(max_steps),
        max_backtracks = as.integer(max_backtracks),
        gradient_tolerance = GRADIENT_TOL
      ),
      algo = "sha256"
    )

    root <- file.path(
      bridge_run,
      "refinement_univoltine_bridge"
    )

    dir.create(
      root,
      recursive = TRUE,
      showWarnings = FALSE
    )

    out <- file.path(
      root,
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
      source_bridge_version = SOURCE_BRIDGE_VERSION,
      signature = signature,
      bridge_run = bridge_run,
      source_run = source_run,
      prepared_file = prepared_file,
      source_fit_file = source_fit_file,
      source_optimized_file = source_optimized_file,
      out = out,
      root = normalizePath(
        root,
        winslash = "/",
        mustWork = TRUE
      ),
      lock = file.path(root, "queue.lock"),
      max_steps = as.integer(max_steps),
      max_backtracks = as.integer(max_backtracks),
      minimum_free_GiB = minimum_free_GiB,
      manifest = manifest
    )

    cfg
  }


  # ---------------------------------------------------------------------------
  # Review archive
  # ---------------------------------------------------------------------------

  archive_review <- function(cfg) {

    if (file.exists(cfg$log_file)) {

      try(
        file.copy(
          cfg$log_file,
          file.path(
            cfg$out,
            "log_snapshot.txt"
          ),
          overwrite = TRUE
        ),
        silent = TRUE
      )
    }

    rel <- list.files(
      cfg$out,
      recursive = TRUE,
      full.names = FALSE
    )

    rel <- rel[
      grepl(
        "\\.(csv|txt|R)$",
        rel,
        ignore.case = TRUE
      )
    ]

    zip_path <- file.path(
      cfg$out,
      "refinement_results_to_review.zip"
    )

    if (requireNamespace("zip", quietly = TRUE)) {

      tryCatch(
        zip::zipr(
          zip_path,
          files = rel,
          root = cfg$out,
          mode = "mirror",
          include_directories = FALSE
        ),
        error = function(e) {
          message(
            "Optional review ZIP failed; ",
            "CSVs remain safe: ",
            conditionMessage(e)
          )
        }
      )

    } else {

      old <- getwd()
      on.exit(setwd(old), add = TRUE)

      setwd(cfg$out)

      try(
        utils::zip(
          zipfile = zip_path,
          files = rel
        ),
        silent = TRUE
      )
    }

    invisible(zip_path)
  }


  # ---------------------------------------------------------------------------
  # Core worker
  # ---------------------------------------------------------------------------

  refine_worker <- function(cfg) {

    check_packages()
    self_test()
    disk_check(
      cfg$out,
      cfg$minimum_free_GiB
    )

    # Check that sources still match the signature created by the parent.
    need(
      identical(
        unname(tools::md5sum(cfg$manifest$path)),
        as.character(cfg$manifest$md5)
      ),
      paste0(
        "One or more source files changed after configuration. ",
        "No refinement was attempted."
      )
    )

    csv(
      cfg$manifest,
      file.path(
        cfg$out,
        "source_manifest.csv"
      )
    )

    source_fit <- readRDS(
      cfg$source_fit_file
    )

    source_optimized <- readRDS(
      cfg$source_optimized_file
    )

    prepared <- readRDS(
      cfg$prepared_file
    )

    need(
      identical(source_fit$version, SOURCE_BRIDGE_VERSION),
      "Unexpected source bridge fit version."
    )

    need(
      identical(prepared$version, SOURCE_BRIDGE_VERSION),
      "Unexpected prepared bridge version."
    )

    need(
      identical(source_fit$group, "univoltine") &&
        identical(prepared$group, "univoltine"),
      "Wrong source group."
    )

    need(
      identical(
        source_fit$signature,
        prepared$signature
      ),
      "Source fit and prepared bridge signatures disagree."
    )

    need(
      identical(
        source_optimized$signature,
        prepared$signature
      ),
      "Optimized checkpoint and prepared bridge signatures disagree."
    )

    need(
      identical(
        source_fit$source_data_signature,
        prepared$source_data_signature
      ),
      "Source abundance data signature mismatch."
    )

    need(
      isTRUE(source_fit$checks$pdHess) &&
        source_fit$checks$convergence_code == 1L,
      paste0(
        "This refinement expects the reviewed bridge fit with ",
        "source convergence code 1 and positive Hessian."
      )
    )

    need(
      isTRUE(
        source_fit$checks$gradient_covariance_norm < GRADIENT_TOL
      ),
      paste0(
        "Source scaled gradient is not already below the documented ",
        "threshold; this bounded adjudication is not the correct workflow."
      )
    )

    need(
      same_num(
        source_fit$optimizer$par,
        source_optimized$optimizer$par
      ),
      "Source fit and optimized checkpoint have different outer parameters."
    )

    need(
      same_num(
        source_fit$beta,
        source_fit$optimizer$par[
          names(source_fit$optimizer$par) == "beta"
        ]
      ),
      "Source beta mapping is inconsistent."
    )

    npar <- length(
      source_fit$optimizer$par
    )

    expected_names <- c(
      rep("beta", ncol(prepared$X)),
      rep("log_sd", 5L),
      "rho_raw",
      "log_shape",
      "log_range",
      "log_spatial_sd"
    )

    need(
      identical(
        names(source_fit$optimizer$par),
        expected_names
      ),
      "Marginal parameter order is not the expected bridge order."
    )

    need(
      identical(
        colnames(prepared$X),
        names(source_fit$beta)
      ),
      "Prepared fixed-effect columns and source beta names differ."
    )

    source_cov <- checked_covariance(
      source_fit$cov_fixed,
      npar
    )

    checkpoint_file <- file.path(
      cfg$out,
      "checkpoint.rds"
    )

    if (file.exists(checkpoint_file)) {

      cp <- readRDS(checkpoint_file)

      need(
        identical(cp$signature, cfg$signature) &&
          identical(
            cp$source_bridge_signature,
            prepared$signature
          ),
        "Refinement checkpoint does not match this source/settings."
      )

    } else {

      cp <- list(
        signature = cfg$signature,
        source_bridge_signature = prepared$signature,
        step = 0L,
        par = source_fit$optimizer$par,
        parameters = source_fit$parameters,
        objective = source_fit$optimizer$objective,
        covariance = source_cov,
        pdHess = TRUE,
        fresh_covariance = FALSE,
        needs_uncertainty = FALSE,
        gradient = NULL,
        history = data.frame(),
        line_search = data.frame()
      )
    }

    need(
      cp$step <= cfg$max_steps,
      "Saved accepted-step count exceeds configured maximum."
    )

    phase(
      cfg,
      "BUILDING_TMB",
      paste0(
        "Reconstructing the exact univoltine bridge likelihood at step ",
        cp$step
      )
    )

    load_dll(cfg$source_run)

    obj <- make_objective(
      prepared$data,
      cp$parameters
    )

    on.exit(
      try(TMB::FreeADFun(obj), silent = TRUE),
      add = TRUE
    )

    need(
      identical(
        names(obj$par),
        names(cp$par)
      ) &&
        length(obj$par) == length(cp$par),
      "Rebuilt marginal parameter vector has changed."
    )

    obj$par <- cp$par

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

      z <- obj$env$parList(
        par = obj$env$last.par
      )

      for (nm in c(
        "beta",
        "log_sd",
        "rho_raw",
        "log_shape",
        "log_range",
        "log_spatial_sd"
      )) {

        need(
          same_num(
            as.numeric(z[[nm]]),
            as.numeric(
              par[names(par) == nm]
            )
          ),
          paste(
            "Current parameter arrays are misaligned:",
            nm
          )
        )
      }

      need(
        all(vapply(
          z,
          function(x) {
            is.numeric(x) &&
              all(is.finite(x))
          },
          logical(1)
        )),
        "Non-finite current parameter array."
      )

      z
    }


    phase(
      cfg,
      "VERIFYING_RECONSTRUCTION",
      "Replaying objective, fixed effects and fitted eta"
    )

    replay <- withCallingHandlers(
      obj$fn(cp$par),
      warning = capture_warning
    )

    need(
      is.finite(replay) &&
        abs(replay - cp$objective) < 0.001,
      paste0(
        "Saved negative log-likelihood was not reproduced within 0.001; ",
        "no refinement attempted."
      )
    )

    report <- obj$report(
      obj$env$last.par
    )

    eta_error <- if (cp$step == 0L) {
      max(
        abs(
          as.numeric(report$eta) -
            as.numeric(source_fit$eta)
        )
      )
    } else {
      NA_real_
    }

    if (cp$step == 0L) {

      need(
        is.finite(eta_error) &&
          eta_error < 1e-5,
        paste0(
          "Saved fitted eta was not reproduced; ",
          "no refinement attempted."
        )
      )
    }

    csv(
      data.frame(
        step = cp$step,
        n_observations = nrow(prepared$frame),
        source_objective = cp$objective,
        replayed_objective = replay,
        absolute_objective_error =
          abs(replay - cp$objective),
        eta_max_abs_error = eta_error,
        unchanged_likelihood = TRUE,
        unchanged_rows_design = TRUE
      ),
      file.path(
        cfg$out,
        paste0(
          "reconstruction_checks_step_",
          cp$step,
          ".csv"
        )
      )
    )

    cp$objective <- as.numeric(replay)
    cp$parameters <- current_parameters(cp$par)

    cp$gradient <- withCallingHandlers(
      as.numeric(obj$gr(cp$par)),
      warning = capture_warning
    )

    need(
      length(cp$gradient) == npar &&
        all(is.finite(cp$gradient)),
      "Invalid reconstructed gradient."
    )


    save_checkpoint <- function() {

      disk_check(
        cfg$out,
        5,
        verbose = FALSE
      )

      atomic_save(
        cp,
        checkpoint_file
      )
    }


    assess_uncertainty <- function() {

      phase(
        cfg,
        "UNCERTAINTY",
        paste(
          "Step",
          cp$step,
          "| computing a fresh marginal Hessian/covariance"
        )
      )

      disk_check(
        cfg$out,
        cfg$minimum_free_GiB,
        verbose = FALSE
      )

      sd <- withCallingHandlers(
        TMB::sdreport(
          obj,
          par.fixed = cp$par,
          getJointPrecision = FALSE,
          getReportCovariance = FALSE
        ),
        warning = capture_warning
      )

      need(
        is.matrix(sd$cov.fixed) &&
          identical(
            dim(sd$cov.fixed),
            as.integer(c(npar, npar))
          ),
        paste0(
          "sdreport returned an unexpected covariance shape. ",
          "Accepted parameters remain checkpointed."
        )
      )

      cp$covariance <<- sd$cov.fixed
      cp$pdHess <<- isTRUE(sd$pdHess)

      cp$gradient <<- withCallingHandlers(
        as.numeric(obj$gr(cp$par)),
        warning = capture_warning
      )

      value <- withCallingHandlers(
        obj$fn(cp$par),
        warning = capture_warning
      )

      need(
        is.finite(value) &&
          abs(value - cp$objective) < 1e-5,
        paste0(
          "Objective changed during uncertainty calculation. ",
          "Parameter checkpoint remains safe."
        )
      )

      cp$objective <<- as.numeric(value)
      cp$parameters <<- current_parameters(cp$par)
      cp$needs_uncertainty <<- FALSE
      cp$fresh_covariance <<- TRUE

      rm(sd)
      invisible(gc())

      save_checkpoint()
    }


    assessment <- function() {

      covariance_ok <- tryCatch(
        {
          checked_covariance(
            cp$covariance,
            npar
          )
          TRUE
        },
        error = function(e) FALSE
      )

      norm <- if (
        covariance_ok &&
        all(is.finite(cp$gradient))
      ) {
        gradient_norm(
          cp$gradient,
          cp$covariance
        )
      } else {
        Inf
      }

      pass <-
        isTRUE(cp$pdHess) &&
        covariance_ok &&
        isTRUE(cp$fresh_covariance) &&
        !isTRUE(cp$needs_uncertainty) &&
        is.finite(norm) &&
        norm < GRADIENT_TOL &&
        is.finite(cp$objective) &&
        cp$objective <=
          source_fit$optimizer$objective + 0.001

      list(
        covariance_ok = covariance_ok,
        norm = norm,
        pass = pass
      )
    }


    add_history <- function(note) {

      a <- assessment()

      row <- data.frame(
        step = cp$step,
        objective = cp$objective,
        change_from_source =
          cp$objective -
          source_fit$optimizer$objective,
        max_abs_gradient =
          max(abs(cp$gradient)),
        gradient_covariance_norm = a$norm,
        pdHess = isTRUE(cp$pdHess),
        fresh_covariance =
          isTRUE(cp$fresh_covariance),
        numerical_checks_ok = a$pass,
        note = note,
        time = format(
          Sys.time(),
          "%Y-%m-%d %H:%M:%S"
        ),
        stringsAsFactors = FALSE
      )

      cp$history <<- rbind(
        cp$history,
        row
      )

      csv(
        cp$history,
        file.path(
          cfg$out,
          "optimization_history.csv"
        )
      )

      save_checkpoint()

      phase(
        cfg,
        "NUMERICAL_CHECK",
        paste(
          "Step",
          cp$step,
          "| scaled gradient",
          signif(a$norm, 7),
          "| criterion <",
          GRADIENT_TOL,
          "|",
          note
        )
      )

      a
    }


    if (isTRUE(cp$needs_uncertainty)) {
      assess_uncertainty()
    }

    a <- add_history(
      if (cp$step == 0L) {
        "Replayed source bridge fit"
      } else {
        "Resumed accepted refinement checkpoint"
      }
    )

    # Even though the source covariance already suggested norm < 0.05,
    # require one fresh covariance calculation at the exact parameter point.
    if (
      a$norm < GRADIENT_TOL &&
      !isTRUE(cp$fresh_covariance)
    ) {

      cp$needs_uncertainty <- TRUE
      save_checkpoint()

      assess_uncertainty()

      a <- add_history(
        "Fresh covariance at current parameter point"
      )
    }

    stop_reason <- "STEP_LIMIT_REACHED"

    while (
      !a$pass &&
      cp$step < cfg$max_steps
    ) {

      if (
        !a$covariance_ok ||
        !isTRUE(cp$pdHess)
      ) {

        stop_reason <- "COVARIANCE_REVIEW"
        break
      }

      step_id <- cp$step + 1L

      disk_check(
        cfg$out,
        cfg$minimum_free_GiB,
        verbose = FALSE
      )

      phase(
        cfg,
        "NEWTON_DIRECTION",
        paste(
          "Step",
          step_id,
          "of at most",
          cfg$max_steps,
          "| all 19 marginal parameters remain free"
        )
      )

      trial <- backtrack(
        cp$par,
        cp$objective,
        cp$gradient,
        cp$covariance,
        evaluate = function(par, j, alpha) {

          phase(
            cfg,
            "LINE_SEARCH",
            paste(
              "Step",
              step_id,
              "| trial",
              j + 1L,
              "of",
              cfg$max_backtracks + 1L,
              "| alpha",
              alpha
            )
          )

          as.numeric(
            withCallingHandlers(
              obj$fn(par),
              warning = capture_warning
            )
          )
        },
        max_backtracks = cfg$max_backtracks
      )

      line <- cbind(
        step = step_id,
        trial$history
      )

      cp$line_search <- rbind(
        cp$line_search,
        line
      )

      csv(
        cp$line_search,
        file.path(
          cfg$out,
          "line_search_history.csv"
        )
      )

      if (!trial$accepted) {

        # Restore last accepted point after rejected trials.
        obj$fn(cp$par)

        cp$parameters <- current_parameters(cp$par)

        save_checkpoint()

        stop_reason <- "NO_ACCEPTABLE_DESCENT"
        break
      }

      cp$par <- trial$par
      cp$objective <- trial$objective
      cp$step <- step_id
      cp$parameters <- current_parameters(cp$par)

      cp$needs_uncertainty <- TRUE
      cp$fresh_covariance <- FALSE
      cp$covariance <- NULL
      cp$gradient <- NULL
      cp$pdHess <- NA

      save_checkpoint()

      atomic_save(
        cp,
        file.path(
          cfg$out,
          sprintf(
            "accepted_step_%02d_parameters.rds",
            cp$step
          )
        )
      )

      phase(
        cfg,
        "ACCEPTED_CHECKPOINT_SAVED",
        paste(
          "Step",
          cp$step,
          "| objective",
          format(
            cp$objective,
            digits = 15
          ),
          "| alpha",
          trial$alpha
        )
      )

      assess_uncertainty()

      a <- add_history(
        paste(
          "Accepted damped-Newton step; alpha",
          trial$alpha
        )
      )
    }

    a <- assessment()

    if (a$pass) {

      stop_reason <- if (cp$step == 0L) {
        "SCALED_GRADIENT_CRITERION_MET_AT_SOURCE_POINT"
      } else {
        "SCALED_GRADIENT_CRITERION_MET_AFTER_REFINEMENT"
      }
    }

    status <- if (a$pass) {
      "FIT_OK"
    } else {
      "NUMERICAL_REVIEW_REQUIRED"
    }

    phase(
      cfg,
      "EXPORTING",
      paste(
        status,
        "|",
        stop_reason,
        "| source bridge fit remains untouched"
      )
    )

    obj$fn(cp$par)
    cp$parameters <- current_parameters(cp$par)

    report <- obj$report(
      obj$env$last.par
    )

    beta_ids <- which(
      names(cp$par) == "beta"
    )

    beta <- setNames(
      as.numeric(cp$par[beta_ids]),
      colnames(prepared$X)
    )

    Vbeta <- cp$covariance[
      beta_ids,
      beta_ids,
      drop = FALSE
    ]

    dimnames(Vbeta) <- list(
      names(beta),
      names(beta)
    )

    extent <- max(
      diff(range(prepared$sites$x_km)),
      diff(range(prepared$sites$y_km))
    )

    custom_code <- if (a$pass) {
      0L
    } else {
      1L
    }

    checks <- data.frame(
      group = "univoltine",
      convergence_code = custom_code,
      source_optimizer_convergence_code =
        source_fit$optimizer$convergence,
      pdHess = isTRUE(cp$pdHess),
      max_abs_gradient =
        max(abs(cp$gradient)),
      gradient_covariance_norm = a$norm,
      numerical_checks_ok = a$pass,
      logLik = -cp$objective,
      parameter_df = npar,
      AIC = 2 * cp$objective + 2 * npar,
      nobs = nrow(prepared$frame),
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
      accepted_Newton_steps = cp$step,
      refinement_reason = stop_reason,
      stringsAsFactors = FALSE
    )

    optimizer <- list(
      par = cp$par,
      objective = cp$objective,
      convergence = custom_code,
      message = stop_reason,
      method = paste0(
        "custom numerical adjudication / damped Newton; ",
        "covariance-scaled gradient < ",
        GRADIENT_TOL
      ),
      iterations = cp$step,
      source_optimizer = source_fit$optimizer
    )

    refined_fit <- list(
      version = VERSION,
      source_bridge_version = SOURCE_BRIDGE_VERSION,
      source_bridge_signature = prepared$signature,
      refinement_signature = cfg$signature,
      parameters = cp$parameters,
      optimizer = optimizer,
      beta = beta,
      V = Vbeta,
      cov_fixed = cp$covariance,
      checks = checks,
      eta = as.numeric(report$eta),
      field_at_sites = report$spatial,
      group = "univoltine",
      variant = "bridge",
      warnings = unique(warnings),
      source_fit = cfg$source_fit_file,
      source_optimizer_code =
        source_fit$optimizer$convergence
    )

    dest <- file.path(
      cfg$out,
      "univoltine__bridge"
    )

    dir.create(
      dest,
      recursive = TRUE,
      showWarnings = FALSE
    )

    atomic_save(
      refined_fit,
      file.path(dest, "fit.rds")
    )

    atomic_save(
      list(
        version = VERSION,
        signature = cfg$signature,
        source_bridge_signature = prepared$signature,
        parameters = cp$parameters,
        optimizer = optimizer
      ),
      file.path(
        dest,
        "optimized_parameters.rds"
      )
    )

    csv(
      checks,
      file.path(
        dest,
        "fit_checks.csv"
      )
    )

    writeLines(
      unique(warnings),
      file.path(
        dest,
        "warnings.txt"
      )
    )

    refined_fx <- fixed_table(
      beta,
      Vbeta
    )

    csv(
      refined_fx,
      file.path(
        dest,
        "fixed_effects.csv"
      )
    )

    csv(
      data.frame(
        term = rownames(Vbeta),
        Vbeta,
        check.names = FALSE
      ),
      file.path(
        dest,
        "fixed_effect_covariance.csv"
      )
    )

    # -----------------------------------------------------------------------
    # Before/after numerical comparison
    # -----------------------------------------------------------------------

    before_fx <- fixed_table(
      source_fit$beta,
      source_fit$V
    )

    comparison <- merge(
      before_fx,
      refined_fx,
      by = "term",
      sort = FALSE,
      all = TRUE,
      suffixes = c("_before", "_after")
    )

    comparison$estimate_change <-
      comparison$estimate_after -
      comparison$estimate_before

    comparison$change_in_before_SE <-
      comparison$estimate_change /
      comparison$std.error_before

    comparison$sign_changed <-
      sign(comparison$estimate_before) !=
      sign(comparison$estimate_after)

    csv(
      comparison,
      file.path(
        cfg$out,
        "coefficient_comparison_before_after.csv"
      )
    )

    common_checks <- intersect(
      names(source_fit$checks),
      names(checks)
    )

    csv(
      rbind(
        cbind(
          stage = "source_bridge",
          source_fit$checks[common_checks]
        ),
        cbind(
          stage = "refined_bridge",
          checks[common_checks]
        )
      ),
      file.path(
        cfg$out,
        "fit_checks_before_after.csv"
      )
    )

    csv(
      data.frame(
        index = seq_along(cp$par),
        term = names(cp$par),
        before =
          as.numeric(source_fit$optimizer$par),
        after = as.numeric(cp$par),
        change =
          as.numeric(
            cp$par -
            source_fit$optimizer$par
          ),
        gradient_after = cp$gradient,
        stringsAsFactors = FALSE
      ),
      file.path(
        cfg$out,
        "all_marginal_parameter_changes.csv"
      )
    )

    # -----------------------------------------------------------------------
    # Refined key bridge effects
    # -----------------------------------------------------------------------

    key_terms <- c(
      "year_decade:latitude_within_species_10deg",
      "year_decade:longitude_within_species_10deg",
      "year_decade:onset_plasticity_z",
      "year_decade:offset_plasticity_bio_z"
    )

    key <- refined_fx[
      match(key_terms, refined_fx$term),
      ,
      drop = FALSE
    ]

    need(
      !anyNA(key$term),
      "One or more refined key bridge terms are missing."
    )

    key$group <- "univoltine"

    key$unit <- c(
      "+10 degrees latitude",
      "+10 degrees longitude",
      "+1 SD onset advancement plasticity",
      "+1 SD stronger offset advancement"
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

    key <- key[c(
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
    )]

    csv(
      key,
      file.path(
        dest,
        "key_bridge_effects.csv"
      )
    )

    # -----------------------------------------------------------------------
    # Plasticity-only source vs refined bridge
    # -----------------------------------------------------------------------

    plasticity_source_fx <- fixed_table(
      prepared$source_fit_beta,
      prepared$source_fit_V
    )

    plasticity_terms <- c(
      "year_decade:onset_plasticity_z",
      "year_decade:offset_plasticity_bio_z"
    )

    source_sub <- plasticity_source_fx[
      match(
        plasticity_terms,
        plasticity_source_fx$term
      ),
      ,
      drop = FALSE
    ]

    refined_sub <- refined_fx[
      match(
        plasticity_terms,
        refined_fx$term
      ),
      ,
      drop = FALSE
    ]

    plasticity_comparison <- data.frame(
      group = "univoltine",
      term = plasticity_terms,
      plasticity_only_estimate =
        source_sub$estimate,
      plasticity_only_low =
        source_sub$conf.low,
      plasticity_only_high =
        source_sub$conf.high,
      plasticity_only_p_Wald =
        source_sub$p_Wald,
      refined_bridge_estimate =
        refined_sub$estimate,
      refined_bridge_low =
        refined_sub$conf.low,
      refined_bridge_high =
        refined_sub$conf.high,
      refined_bridge_p_Wald =
        refined_sub$p_Wald,
      absolute_change_bridge_minus_source =
        refined_sub$estimate -
        source_sub$estimate,
      proportional_change_in_magnitude =
        abs(refined_sub$estimate) /
        abs(source_sub$estimate),
      same_sign =
        sign(refined_sub$estimate) ==
        sign(source_sub$estimate),
      stringsAsFactors = FALSE
    )

    csv(
      plasticity_comparison,
      file.path(
        dest,
        "plasticity_effect_comparison.csv"
      )
    )

    # -----------------------------------------------------------------------
    # Nested geography addition LRT
    # -----------------------------------------------------------------------

    source_checks <- prepared$source_fit_checks

    LR <- 2 * (
      checks$logLik -
      source_checks$logLik
    )

    df_added <-
      checks$parameter_df -
      source_checks$parameter_df

    nested <- data.frame(
      group = "univoltine",
      source_model =
        "plasticity_only_spatial",
      bridge_model =
        "plasticity_plus_geography_spatial_refined",
      logLik_source =
        source_checks$logLik,
      logLik_bridge =
        checks$logLik,
      AIC_source =
        source_checks$AIC,
      AIC_bridge =
        checks$AIC,
      delta_AIC_bridge_minus_source =
        checks$AIC -
        source_checks$AIC,
      LR_statistic = max(0, LR),
      added_parameters = df_added,
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
      note = paste0(
        "Joint test of four added geographic fixed terms: ",
        "latitude, longitude, year:latitude, year:longitude."
      ),
      stringsAsFactors = FALSE
    )

    csv(
      nested,
      file.path(
        dest,
        "nested_geography_addition_test.csv"
      )
    )

    # -----------------------------------------------------------------------
    # Joint Wald test of TEMPORAL geography only (2 df)
    # -----------------------------------------------------------------------

    temporal_geo_terms <- c(
      "year_decade:latitude_within_species_10deg",
      "year_decade:longitude_within_species_10deg"
    )

    q <- beta[temporal_geo_terms]

    Vq <- Vbeta[
      temporal_geo_terms,
      temporal_geo_terms,
      drop = FALSE
    ]

    wald <- as.numeric(
      crossprod(
        q,
        solve(Vq, q)
      )
    )

    joint_geo <- data.frame(
      group = "univoltine",
      test =
        "year:latitude + year:longitude",
      Wald_chisq = wald,
      df = 2L,
      p_Wald_joint =
        stats::pchisq(
          wald,
          df = 2,
          lower.tail = FALSE
        ),
      interpretation = paste0(
        "Joint 2-df test of geographic effects on abundance trend only; ",
        "main geographic effects are not included in this contrast."
      ),
      stringsAsFactors = FALSE
    )

    csv(
      joint_geo,
      file.path(
        dest,
        "joint_temporal_geography_Wald.csv"
      )
    )

    # -----------------------------------------------------------------------
    # Onset-offset reparameterization, controlling geography
    # -----------------------------------------------------------------------

    onset_term <-
      "year_decade:onset_plasticity_z"

    offset_term <-
      "year_decade:offset_plasticity_bio_z"

    balance <- rbind(
      linear_contrast(
        beta,
        Vbeta,
        weights = c(
          setNames(1, onset_term),
          setNames(1, offset_term)
        ),
        label = "coordinated_advancement",
        unit = paste0(
          "+1 SD onset advancement AND ",
          "+1 SD offset advancement"
        ),
        interpretation = paste0(
          "Whole-season coordinated advancement, ",
          "controlling season-extension axis and geography."
        )
      ),
      linear_contrast(
        beta,
        Vbeta,
        weights = c(
          setNames(1, onset_term),
          setNames(-1, offset_term)
        ),
        label = "season_extension",
        unit = paste0(
          "+1 SD onset advancement AND ",
          "1 SD relatively later termination"
        ),
        interpretation = paste0(
          "Earlier onset combined with relatively later offset, ",
          "controlling coordinated-advancement axis and geography."
        )
      )
    )

    balance$group <- "univoltine"

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
        dest,
        "onset_offset_balance_controlling_geography.csv"
      )
    )

    # -----------------------------------------------------------------------
    # Final result object
    # -----------------------------------------------------------------------

    result <- list(
      signature = cfg$signature,
      source_bridge_signature = prepared$signature,
      status = status,
      reason = stop_reason,
      started_at = cfg$launched_at,
      finished_at = Sys.time(),
      output = cfg$out,
      fit_file = file.path(
        dest,
        "fit.rds"
      ),
      accepted_steps = cp$step,
      source_optimizer_code =
        source_fit$optimizer$convergence,
      scaled_gradient_before =
        source_fit$checks$gradient_covariance_norm,
      scaled_gradient_after = a$norm
    )

    atomic_save(
      result,
      file.path(
        cfg$out,
        "result.rds"
      )
    )

    csv(
      data.frame(
        status = status,
        reason = stop_reason,
        accepted_steps = cp$step,
        gradient_threshold = GRADIENT_TOL,
        source_optimizer_code =
          source_fit$optimizer$convergence,
        gradient_before =
          result$scaled_gradient_before,
        gradient_after =
          result$scaled_gradient_after,
        objective_change =
          cp$objective -
          source_fit$optimizer$objective,
        latitude_before =
          source_fit$beta[
            "year_decade:latitude_within_species_10deg"
          ],
        latitude_after =
          beta[
            "year_decade:latitude_within_species_10deg"
          ],
        onset_before =
          source_fit$beta[
            "year_decade:onset_plasticity_z"
          ],
        onset_after =
          beta[
            "year_decade:onset_plasticity_z"
          ],
        offset_before =
          source_fit$beta[
            "year_decade:offset_plasticity_bio_z"
          ],
        offset_after =
          beta[
            "year_decade:offset_plasticity_bio_z"
          ],
        fit_file = result$fit_file,
        stringsAsFactors = FALSE
      ),
      file.path(
        cfg$out,
        "refinement_summary.csv"
      )
    )

    writeLines(
      c(
        "REFINEMENT / NUMERICAL ADJUDICATION OF THE UNIVOLTINE BRIDGE MODEL",
        "",
        paste("Source bridge run:", cfg$bridge_run),
        paste("Source fit:", cfg$source_fit_file),
        paste("Workflow version:", VERSION),
        "",
        "The model, observations, predictors, scaling, likelihood, random",
        "effects, AR1 process, source mesh and annual Matern field are unchanged.",
        "",
        paste0(
          "Pass criterion: positive Hessian/covariance and covariance-scaled ",
          "gradient < ", GRADIENT_TOL, "."
        ),
        "",
        paste(
          "Accepted damped-Newton steps:",
          cp$step
        ),
        paste(
          "Final reason:",
          stop_reason
        ),
        "",
        "A refined convergence code of 0 denotes the documented custom",
        "criterion, not a new nlminb/BFGS return code.",
        "",
        "Source files and the multivoltine bridge model remain untouched.",
        "",
        "Additional exports include:",
        "- joint 2-df Wald test of year:latitude + year:longitude",
        "- onset/offset reparameterization into coordinated advancement and",
        "  season extension while controlling geography",
        "- source-vs-refined coefficient and numerical comparisons"
      ),
      file.path(
        cfg$out,
        "README_refinement.txt"
      )
    )

    writeLines(
      capture.output(
        utils::sessionInfo()
      ),
      file.path(
        cfg$out,
        "sessionInfo.txt"
      )
    )

    phase(
      cfg,
      status,
      paste(
        "Scaled gradient",
        signif(a$norm, 7),
        "|",
        stop_reason
      ),
      terminal = TRUE
    )

    archive_review(cfg)

    result
  }


  # ---------------------------------------------------------------------------
  # Lock / background process helpers
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

    need(
      !is.null(owner) &&
        !is.null(owner$pid) &&
        !is.null(owner$create_time),
      paste(
        "Unidentified refinement lock; inspect:",
        lock
      )
    )

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

    owner_file <- file.path(
      cfg$lock,
      "owner.rds"
    )

    owner <- if (file.exists(owner_file)) {
      tryCatch(
        readRDS(owner_file),
        error = function(e) NULL
      )
    } else {
      NULL
    }

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

    invisible(NULL)
  }


  worker_entry <- function(cfg) {

    on.exit(
      release_lock(cfg),
      add = TRUE
    )

    tryCatch(
      refine_worker(cfg),
      error = function(e) {

        writeLines(
          conditionMessage(e),
          file.path(
            cfg$out,
            "error.txt"
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

        try(
          archive_review(cfg),
          silent = TRUE
        )

        list(
          status = "STOPPED_WITH_ERROR",
          error = conditionMessage(e),
          output = cfg$out
        )
      }
    )
  }


  dump_engine <- function(out) {

    e <- environment(dump_engine)

    all_names <- ls(
      e,
      all.names = TRUE
    )

    functions <- all_names[
      vapply(
        all_names,
        function(nm) {
          is.function(
            get(nm, envir = e)
          )
        },
        logical(1)
      )
    ]

    path <- file.path(
      out,
      "refinement_engine.R"
    )

    dump(
      c(
        "VERSION",
        "SOURCE_BRIDGE_VERSION",
        "DEFAULT_BRIDGE_RUN",
        "CORE_VERSIONS",
        "GRADIENT_TOL",
        functions
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
  # Public run/status
  # ---------------------------------------------------------------------------

  run <- function(
    bridge_run = DEFAULT_BRIDGE_RUN,
    max_steps = 2L,
    max_backtracks = 6L,
    minimum_free_GiB = 25,
    monitor = FALSE
  ) {

    check_packages()
    self_test()

    need(
      is.numeric(max_steps) &&
        length(max_steps) == 1L &&
        is.finite(max_steps) &&
        max_steps == floor(max_steps) &&
        max_steps >= 0L &&
        max_steps <= 3L,
      "max_steps must be 0, 1, 2 or 3."
    )

    need(
      is.numeric(max_backtracks) &&
        length(max_backtracks) == 1L &&
        is.finite(max_backtracks) &&
        max_backtracks == floor(max_backtracks) &&
        max_backtracks >= 0L &&
        max_backtracks <= 8L,
      "max_backtracks must be 0 to 8."
    )

    need(
      is.numeric(minimum_free_GiB) &&
        length(minimum_free_GiB) == 1L &&
        is.finite(minimum_free_GiB) &&
        minimum_free_GiB >= 10,
      "minimum_free_GiB must be at least 10."
    )

    current <- state$job

    if (
      !is.null(current) &&
      !is.null(current$process) &&
      current$process$is_alive()
    ) {

      message(
        "The univoltine bridge refinement is already running."
      )

      return(invisible(current))
    }

    cfg <- configure(
      bridge_run = bridge_run,
      max_steps = max_steps,
      max_backtracks = max_backtracks,
      minimum_free_GiB = minimum_free_GiB
    )

    disk_check(
      cfg$out,
      cfg$minimum_free_GiB
    )

    completed_file <- file.path(
      cfg$out,
      "result.rds"
    )

    if (file.exists(completed_file)) {

      result <- readRDS(completed_file)

      need(
        identical(
          result$signature,
          cfg$signature
        ),
        "Completed result signature mismatch."
      )

      message(
        "Matching refinement is already complete: ",
        result$status
      )

      job <- list(
        process = NULL,
        out = cfg$out,
        log = NULL,
        completed = TRUE
      )

      state$job <- job

      options(
        phenoimpact.bridge.refinement.job = job
      )

      if (monitor) {
        status(job)
      }

      return(invisible(job))
    }

    if (dir.exists(cfg$lock)) {

      alive <- lock_owner_alive(cfg$lock)

      need(
        !alive,
        paste(
          "A matching bridge refinement appears to be running:",
          cfg$lock
        )
      )

      unlink(
        cfg$lock,
        recursive = TRUE
      )
    }

    need(
      dir.create(
        cfg$lock,
        recursive = TRUE
      ),
      "Cannot acquire the refinement lock."
    )

    cfg$lock_token <- paste(
      Sys.getpid(),
      format(
        Sys.time(),
        "%Y%m%d%H%M%OS6"
      ),
      sep = "_"
    )

    cfg$launched_at <- Sys.time()
    cfg$libpath <- .libPaths()

    parent_handle <- ps::ps_handle()

    atomic_save(
      list(
        pid = Sys.getpid(),
        create_time =
          ps::ps_create_time(parent_handle),
        token = cfg$lock_token
      ),
      file.path(
        cfg$lock,
        "owner.rds"
      )
    )

    cfg$engine <- dump_engine(
      cfg$out
    )

    cfg$log_file <- file.path(
      cfg$out,
      paste0(
        "refinement_",
        format(
          Sys.time(),
          "%Y%m%d_%H%M%S"
        ),
        ".log"
      )
    )

    csv(
      cfg$manifest,
      file.path(
        cfg$out,
        "source_manifest.csv"
      )
    )

    atomic_save(
      cfg,
      file.path(
        cfg$out,
        "refinement_config.rds"
      )
    )

    phase(
      cfg,
      "QUEUED",
      paste0(
        "Only univoltine bridge numerical adjudication; ",
        "source and multivoltine fits are read only"
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
      libpath = cfg$libpath,
      stdout = cfg$log_file,
      stderr = "2>&1",
      supervise = TRUE,
      user_profile = FALSE,
      system_profile = FALSE,
      wd = cfg$bridge_run,
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
        token = cfg$lock_token
      ),
      file.path(
        cfg$lock,
        "owner.rds"
      )
    )

    job <- list(
      process = child,
      out = cfg$out,
      log = cfg$log_file,
      started = cfg$launched_at,
      completed = FALSE
    )

    state$job <- job

    options(
      phenoimpact.bridge.refinement.job = job
    )

    success <- TRUE

    message(
      "Univoltine bridge refinement launched."
    )

    message(
      "Deterministic output/checkpoint folder: ",
      cfg$out
    )

    message(
      "At most ",
      cfg$max_steps,
      " accepted Newton steps."
    )

    message(
      "Monitor with: pheno_bridge_refine$status()"
    )

    if (monitor) {
      status(job)
    }

    invisible(job)
  }


  status <- function(job = NULL) {

    if (is.null(job)) {
      job <- state$job
    }

    if (is.null(job)) {
      job <- getOption(
        "phenoimpact.bridge.refinement.job"
      )
    }

    need(
      !is.null(job),
      "No univoltine bridge refinement is registered."
    )

    progress_file <- file.path(
      job$out,
      "progress.rds"
    )

    if (file.exists(progress_file)) {

      p <- tryCatch(
        readRDS(progress_file),
        error = function(e) NULL
      )

      if (!is.null(p)) {
        print(p)
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

    history_file <- file.path(
      job$out,
      "optimization_history.csv"
    )

    if (file.exists(history_file)) {

      cat("\nNumerical history:\n")

      print(
        utils::read.csv(
          history_file,
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

    result_file <- file.path(
      job$out,
      "result.rds"
    )

    if (file.exists(result_file)) {

      result <- readRDS(result_file)

      cat(
        "\nResult:\n"
      )

      print(result)

      summary_file <- file.path(
        job$out,
        "refinement_summary.csv"
      )

      if (file.exists(summary_file)) {

        cat(
          "\nRefinement summary:\n"
        )

        print(
          utils::read.csv(
            summary_file,
            stringsAsFactors = FALSE
          ),
          row.names = FALSE
        )
      }

      balance_file <- file.path(
        job$out,
        "univoltine__bridge",
        "onset_offset_balance_controlling_geography.csv"
      )

      if (file.exists(balance_file)) {

        cat(
          "\nOnset-offset balance controlling geography:\n"
        )

        print(
          utils::read.csv(
            balance_file,
            stringsAsFactors = FALSE
          ),
          row.names = FALSE
        )
      }

      zip_file <- file.path(
        job$out,
        "refinement_results_to_review.zip"
      )

      if (file.exists(zip_file)) {

        cat(
          "\nFILE TO SEND:\n",
          zip_file,
          "\n",
          sep = ""
        )
      }
    }

    error_file <- file.path(
      job$out,
      "error.txt"
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


run_univoltine_bridge_refinement <- function(...) {
  pheno_bridge_refine$run(...)
}


status_univoltine_bridge_refinement <- function(...) {
  pheno_bridge_refine$status(...)
}
