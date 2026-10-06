# ============================================================================
# phenoIMPACT | Refine ONLY the saved onset-adjusted UNIVOLTINE abundance model
# Version 1.0 | 2026-10-05
#
# source(file.choose(), encoding = "UTF-8")
# refinement_job <- pheno_refine_univoltine$run(monitor = FALSE)
# pheno_refine_univoltine$status()
#
# Source: spatial_abundance_onset_adjusted_15km/
#         run_20261004_214012_44136/univoltine__full/fit.rds
#
# The model, observations, plasticities, scaling, AR1 links, random effects and
# 15-km mesh are UNCHANGED. No phenology/extraction/multivoltine job is run.
# Rebuild native TMB state from the saved numeric arrays. Verify the original
# objective, coefficients and predictions BEFORE further optimization.
#
# At most TWO damped Newton updates by default. The direction is -V %*% g,
# using the full marginal parameter covariance (not only the beta block).
# Backtracking accepts only a sufficiently LOWER negative log-likelihood.
# After an accepted step, recalculate TMB::sdreport at that exact point,
# then reassess sqrt(t(g) %*% V %*% g) < 0.05, the original threshold.
# A fresh covariance is required before declaring success. No ridge, changed
# priors, fixed variance parameters, relaxed threshold or automatic BFGS loop.
#
# The optimizer used here is explicitly a custom damped-Newton solver; its
# convergence=0 means its documented numerical criterion has been met. It is
# NOT a newly returned nlminb/BFGS convergence code. Original optimizer status
# is retained separately. Numerical completion does not validate residuals.
#
# Source files are READ ONLY. Results/checkpoints are written in a NEW folder.
# Accepted parameter arrays are saved BEFORE covariance calculation. The
# entire native objective and the full joint precision are NEVER serialized.
# Keep R open; do not run another large model at the same time.
#
# R/TMB references (for the API; not a claim of full-data execution testing):
# https://search.r-project.org/CRAN/refmans/TMB/html/MakeADFun.html
# https://search.r-project.org/CRAN/refmans/TMB/html/sdreport.html
# ============================================================================

pheno_refine_univoltine <- local({
  VERSION <- "univoltine_spatial_abundance_refinement_1.0"
  SOURCE_VERSION <- "updated_plasticity_spatial_abundance_1.0"
  DEFAULT_SOURCE <- paste0(
    "E:/phenoIMPACT project/code/phenoIMPACT/output/population_trends/",
    "spatial_abundance_onset_adjusted_15km/run_20261004_214012_44136")
  CPP_MD5 <- "88e15b8506c2e85a2e1e0d3839101530"
  CORE_VERSIONS <- c(TMB = "1.9.17", Matrix = "1.7.3", fmesher = "0.5.0")
  GRADIENT_TOL <- 0.05
  state <- new.env(parent = emptyenv())

  need <- function(ok, text) if (!isTRUE(ok)) stop(text, call. = FALSE)
  same_num <- function(x, y, tol = 1e-8) {
    is.numeric(x) && is.numeric(y) && length(x) == length(y) &&
      all(is.finite(x)) && all(is.finite(y)) &&
      (length(x) == 0L || max(abs(as.numeric(x) - as.numeric(y))) <= tol)
  }
  write_csv <- function(x, path) utils::write.csv(x, path, row.names = FALSE, na = "")
  atomic_save <- function(x, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    tmp <- tempfile(".writing_", tmpdir = dirname(path), fileext = ".rds")
    on.exit(unlink(tmp), add = TRUE)
    saveRDS(x, tmp, compress = "gzip")
    # Only derived refinement files use this helper, never source inputs.
    bak <- paste0(path, ".previous")
    if (file.exists(bak)) need(unlink(bak) == 0L, paste("Cannot remove old backup:", bak))
    had <- file.exists(path)
    if (had) need(file.rename(path, bak), paste("Cannot retain previous checkpoint:", path))
    if (!file.rename(tmp, path)) {
      if (had) file.rename(bak, path)
      stop("Cannot finalize checkpoint: ", path, call. = FALSE)
    }
    invisible(path)
  }
  phase <- function(cfg, stage, detail = "", terminal = FALSE) {
    now <- Sys.time()
    atomic_save(list(stage = stage, detail = detail, time = now,
      terminal = terminal, started = cfg$launched_at,
      finished = if (terminal) now else NULL), file.path(cfg$out, "progress.rds"))
    message("[", format(now, "%Y-%m-%d %H:%M:%S"), "] ", stage, " | ", detail)
  }
  check_packages <- function() {
    packages <- c(names(CORE_VERSIONS), "callr", "ps", "digest")
    missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
    need(!length(missing), paste("Missing packages:", paste(missing, collapse = ", "),
      "Use the original R library; no packages are installed by this script."))
    current <- vapply(names(CORE_VERSIONS), function(p)
      as.character(utils::packageVersion(p)), character(1))
    need(identical(current, CORE_VERSIONS), paste(
      "Core package versions must match the reviewed fit:",
      paste(names(CORE_VERSIONS), CORE_VERSIONS, collapse = "; ")))
  }
  disk_check <- function(path, min_GiB, verbose = TRUE) {
    d <- ps::ps_disk_usage(path)
    need(nrow(d) == 1L && is.finite(d$available), "Cannot check free disk space.")
    free <- as.numeric(d$available) / 1024^3
    if (verbose) message("Disk preflight: ", round(free, 1), " GiB available.")
    need(free >= min_GiB, paste("Only", round(free, 1), "GiB free; reserve is", min_GiB,
      "GiB. No further expensive step started."))
    invisible(free)
  }
  checked_covariance <- function(V, n) {
    need(is.matrix(V) && identical(dim(V), as.integer(c(n, n))) &&
      all(is.finite(V)), "Marginal parameter covariance is missing or non-finite.")
    symmetry <- max(abs(V - t(V))) / max(1, max(abs(V)))
    need(symmetry < 1e-8, "Marginal covariance is not numerically symmetric.")
    V <- (V + t(V)) / 2
    need(all(diag(V) > 0), "Marginal covariance has a non-positive diagonal.")
    ch <- tryCatch(chol(V), error = function(e) NULL)
    need(!is.null(ch), "Marginal covariance is not positive definite; no ridge is added.")
    V
  }
  gradient_norm <- function(g, V) {
    need(is.numeric(g) && all(is.finite(g)), "Gradient contains non-finite values.")
    V <- checked_covariance(V, length(g))
    value <- as.numeric(crossprod(g, V %*% g))
    need(is.finite(value) && value >= 0, "Invalid covariance-scaled gradient.")
    sqrt(value)
  }
  newton_direction <- function(g, V) {
    V <- checked_covariance(V, length(g))
    norm <- gradient_norm(g, V)
    # Trust safeguard: at most one local covariance-metric unit per trial.
    # In this source fit norm is ~0.0825, so this cap should be inactive.
    shrink <- if (norm > 1) 1 / norm else 1
    direction <- -as.numeric(V %*% g) * shrink
    dot <- sum(g * direction)
    need(all(is.finite(direction)) && is.finite(dot) && dot < 0,
      "No finite descent direction; optimization is stopped for review.")
    list(direction = direction, g_dot_direction = dot, shrink = shrink)
  }
  backtrack <- function(p, value, g, V, evaluate, max_backtracks = 6L) {
    step <- newton_direction(g, V)
    records <- vector("list", max_backtracks + 1L)
    for (j in 0:max_backtracks) {
      alpha <- 0.5^j
      trial <- p + alpha * step$direction
      f_trial <- evaluate(trial, j, alpha)
      bound <- value + 1e-4 * alpha * step$g_dot_direction
      accept <- is.finite(f_trial) && f_trial < value && f_trial <= bound
      records[[j + 1L]] <- data.frame(backtrack = j, alpha = alpha,
        direction_shrink = step$shrink, objective_before = value,
        objective_trial = f_trial, Armijo_bound = bound, accepted = accept)
      if (accept) return(list(accepted = TRUE, par = trial, objective = f_trial,
        alpha = alpha, history = do.call(rbind, records[seq_len(j + 1L)])))
    }
    list(accepted = FALSE, par = p, objective = value, alpha = 0,
      history = do.call(rbind, records))
  }
  self_test <- function() {
    # Deterministic checks of Newton orientation, covariance scaling, rejection,
    # and real backtracking. Pure R; no large data or TMB objective is required.
    H <- matrix(c(4, 1, 1, 3), 2, 2)
    V <- solve(H); p <- c(a = 0.02, b = -0.03)
    fn <- function(x) as.numeric(crossprod(x, H %*% x)) / 2
    gr <- function(x) as.numeric(H %*% x)
    s <- backtrack(p, fn(p), gr(p), V, function(x, j, a) fn(x))
    need(s$accepted && s$alpha == 1 && max(abs(s$par)) < 1e-12,
      "Newton self-test failed.")
    large <- backtrack(p, fn(p), gr(p), 8 * V, function(x, j, a) fn(x))
    need(large$accepted && large$alpha < 1 && large$objective < fn(p),
      "Backtracking self-test failed.")
    reject <- backtrack(p, fn(p), gr(p), V, function(x, j, a) Inf, 2L)
    need(!reject$accepted && identical(reject$par, p) && nrow(reject$history) == 3,
      "Rejected-step self-test failed.")
    need(abs(gradient_norm(gr(p), V)^2 - 2 * fn(p)) < 1e-12,
      "Gradient scaling self-test failed.")
    bad <- tryCatch({checked_covariance(diag(c(1, -1)), 2L); FALSE},
      error = function(e) TRUE)
    need(bad, "Indefinite covariance was not rejected.")
    invisible(TRUE)
  }

  # These two builders preserve the exact reference model construction.
  with_mesh <- function(input, mesh) {
    d <- input$data
    d$use_spatial <- 1L
    d$A <- methods::as(fmesher::fm_basis(mesh$mesh,
      loc = as.matrix(input$sites[c("x_km", "y_km")])), "CsparseMatrix")
    need(all(abs(Matrix::rowSums(d$A) - 1) < 1e-8), "Coordinates lie outside the saved mesh.")
    d$M0 <- mesh$M0; d$M1 <- mesh$M1; d$M2 <- mesh$M2
    d
  }
  make_objective <- function(data, parameters) {
    for (nm in c("A", "M0", "M1", "M2")) data[[nm]] <- methods::as(
      methods::as(methods::as(data[[nm]], "dMatrix"), "generalMatrix"), "TsparseMatrix")
    TMB::MakeADFun(data = data, parameters = parameters,
      random = c("b_species", "b_site", "ar_obs", "field"), map = NULL,
      DLL = "trend_spde", silent = TRUE,
      inner.control = list(maxit = 1000, trace = FALSE))
  }
  fixed_table <- function(beta, V) {
    diagV <- diag(V)
    se <- ifelse(is.finite(diagV) & diagV > 0, sqrt(pmax(0, diagV)), NA_real_)
    data.frame(term = names(beta), estimate = as.numeric(beta), std.error = as.numeric(se),
      conf.low = as.numeric(beta - 1.96 * se), conf.high = as.numeric(beta + 1.96 * se),
      p_Wald = as.numeric(2 * stats::pnorm(-abs(beta / se))))
  }
  predictions <- function(beta, V, years) {
    ans <- list()
    for (v in c("onset_plasticity_z", "offset_plasticity_bio_z")) for (z in c(-1, 0, 1)) {
      it <- paste0("year_decade:", v); tt <- (years - min(years)) / 10
      slope <- beta["year_decade"] + z * beta[it]
      variance <- V["year_decade", "year_decade"] + z^2 * V[it, it] +
        2 * z * V["year_decade", it]
      need(is.finite(variance) && variance >= -1e-10, "Invalid prediction variance.")
      eta <- as.numeric(slope) * tt; se <- sqrt(max(0, variance)) * tt
      ans[[length(ans) + 1L]] <- data.frame(variable = v, plasticity_z = z, YEAR = years,
        percent_change = 100 * expm1(eta), low = 100 * expm1(eta - 1.96 * se),
        high = 100 * expm1(eta + 1.96 * se))
    }
    do.call(rbind, ans)
  }

  fingerprint <- function(cfg) {
    dll <- TMB::dynlib(file.path(cfg$source_run, "compiled", "trend_spde"))
    paths <- c(file.path(cfg$source_run, c("run_config.rds", "input_univoltine.rds",
      "mesh.rds", "compiled/compile_complete.rds", "compiled/trend_spde.cpp",
      "univoltine__full/fit.rds")), dll)
    need(all(file.exists(paths)), paste("Missing required source file(s):",
      paste(paths[!file.exists(paths)], collapse = "\n")))
    answer <- lapply(seq_along(paths), function(i) {
      phase(cfg, "FINGERPRINTING", paste(i, "/", length(paths), basename(paths[i])))
      h <- unname(tools::md5sum(paths[i]))
      need(!is.na(h), paste("Cannot read:", paths[i]))
      data.frame(path = normalizePath(paths[i], winslash = "/", mustWork = TRUE),
        bytes = file.info(paths[i])$size, md5 = h, stringsAsFactors = FALSE)
    })
    manifest <- do.call(rbind, answer)
    signature <- digest::digest(list(version = VERSION, manifest = manifest,
      core = CORE_VERSIONS, engine_md5 = unname(tools::md5sum(cfg$engine)),
      max_steps = cfg$max_steps, max_backtracks = cfg$max_backtracks, tolerance = GRADIENT_TOL), algo = "sha256")
    if (!is.null(cfg$signature)) need(identical(signature, cfg$signature),
      "Source files/settings changed. Resume refused; previous refinement files are untouched.")
    cfg$signature <- signature
    write_csv(manifest, file.path(cfg$out, "source_manifest.csv"))
    atomic_save(cfg, file.path(cfg$out, "refinement_config.rds"))
    cfg
  }
  archive_review <- function(cfg) {
    if (!requireNamespace("zip", quietly = TRUE)) {
      message("Optional zip package unavailable; review CSVs/log remain in the output folder.")
      return(invisible(NULL))
    }
    # Logs are copied to avoid archiving an open Windows log handle.
    if (file.exists(cfg$log_file)) try(file.copy(cfg$log_file,
      file.path(cfg$out, "log_snapshot.txt"), overwrite = TRUE), silent = TRUE)
    rel <- list.files(cfg$out, recursive = TRUE, full.names = FALSE)
    rel <- rel[grepl("\\.(csv|txt|R)$", rel, ignore.case = TRUE)]
    tryCatch(zip::zipr(file.path(cfg$out, "results_to_review.zip"), files = rel,
      root = cfg$out, mode = "mirror", include_directories = FALSE),
      error = function(e) message("Optional review ZIP failed; CSVs are safe: ", conditionMessage(e)))
  }

  refine_worker <- function(cfg) {
    check_packages(); self_test(); disk_check(cfg$out, cfg$minimum_free_GiB)
    cfg <- fingerprint(cfg)
    if (file.exists(file.path(cfg$out, "result.rds"))) {
      done <- readRDS(file.path(cfg$out, "result.rds"))
      need(identical(done$signature, cfg$signature), "Completed refinement signature mismatch.")
      phase(cfg, done$status, "This bounded attempt is already complete; no optimization repeated.", TRUE)
      return(invisible(done))
    }
    phase(cfg, "READING_SAVED_FIT", "Reading UNIVOLTINE arrays only; original files are read-only")
    src_cfg <- readRDS(file.path(cfg$source_run, "run_config.rds"))
    input <- readRDS(file.path(cfg$source_run, "input_univoltine.rds"))
    original <- readRDS(file.path(cfg$source_run, "univoltine__full", "fit.rds"))
    marker <- readRDS(file.path(cfg$source_run, "compiled", "compile_complete.rds"))
    need(identical(src_cfg$version, SOURCE_VERSION) &&
      identical(input$schema, SOURCE_VERSION), "Unexpected source workflow version.")
    need(identical(original$signature, src_cfg$signature) &&
      identical(input$signature, src_cfg$signature) &&
      identical(marker$signature, src_cfg$signature), "Source fit/input/compiled-engine signatures disagree.")
    need(identical(input$group, "univoltine") && identical(original$group, "univoltine") &&
      identical(original$variant, "full") && isTRUE(original$offset_adjusted_for_onset) &&
      isTRUE(input$offset_adjusted_for_onset), "Wrong source model/group or old plasticities.")
    need(identical(original$data_signature, input$data_signature) &&
      identical(input$data_signature, digest::digest(list(input$frame, input$data$X), algo = "sha256")),
      "Saved abundance rows/design do not match the fit.")
    need(isTRUE(original$checks$pdHess) && original$optimizer$convergence == 0L,
      "This refinement expects the reviewed fit with code 0 and positive Hessian.")
    need(identical(marker$cpp_md5, CPP_MD5) &&
      identical(unname(tools::md5sum(file.path(cfg$source_run, "compiled", "trend_spde.cpp"))), CPP_MD5),
      "The likelihood template is not the reviewed Gamma/log spatial model.")
    npar <- length(original$optimizer$par)
    V0 <- checked_covariance(original$cov_fixed, npar)
    beta_ids <- which(names(original$optimizer$par) == "beta")
    need(length(beta_ids) == ncol(input$data$X) &&
      identical(names(original$beta), colnames(input$data$X)) &&
      same_num(original$beta, original$optimizer$par[beta_ids]) &&
      same_num(as.numeric(original$V), as.numeric(V0[beta_ids, beta_ids, drop = FALSE])),
      "Source beta/covariance mapping is inconsistent.")
    need(nrow(input$frame) == length(input$data$y) &&
      all(is.finite(input$data$y) & input$data$y > 0), "Invalid Gamma observation input.")
    checkpoint_path <- file.path(cfg$out, "checkpoint.rds")
    if (file.exists(checkpoint_path)) {
      cp <- readRDS(checkpoint_path)
      need(identical(cp$signature, cfg$signature) &&
        identical(cp$data_signature, input$data_signature), "Refinement checkpoint mismatch.")
    } else {
      cp <- list(signature = cfg$signature, data_signature = input$data_signature,
        step = 0L, par = original$optimizer$par, parameters = original$parameters,
        objective = original$optimizer$objective, covariance = V0,
        pdHess = TRUE, fresh_covariance = FALSE, needs_uncertainty = FALSE,
        gradient = NULL, history = data.frame(), line_search = data.frame())
    }
    need(cp$step <= cfg$max_steps, "Saved step count exceeds the configured limit.")
    mesh <- readRDS(file.path(cfg$source_run, "mesh.rds"))
    phase(cfg, "BUILDING_TMB", "Reconstructing the SAME spatial model at the saved parameter values")
    data <- with_mesh(input, mesh)
    need(identical(dim(cp$parameters$field), as.integer(c(nrow(mesh$M0), length(input$years)))),
      "Saved annual field has unexpected dimensions.")
    dll <- TMB::dynlib(file.path(cfg$source_run, "compiled", "trend_spde"))
    dyn.load(dll); TMB::openmp(n = 1L, DLL = "trend_spde")
    obj <- make_objective(data, cp$parameters)
    on.exit(try(TMB::FreeADFun(obj), silent = TRUE), add = TRUE)
    need(identical(names(obj$par), names(cp$par)) && same_num(obj$par, cp$par),
      "Rebuilt marginal parameters do not match the saved parameter vector.")
    expected_names <- c(rep("beta", ncol(data$X)), rep("log_sd", 5L),
      "rho_raw", "log_shape", "log_range", "log_spatial_sd")
    need(identical(names(cp$par), expected_names), "Marginal parameter order changed.")
    warnings <- character()
    capture_warning <- function(w) {
      warnings <<- unique(c(warnings, conditionMessage(w)))
      message("WARNING: ", conditionMessage(w)); invokeRestart("muffleWarning")
    }
    current_parameters <- function(p) {
      a <- obj$env$parList(par = obj$env$last.par)
      for (nm in c("beta", "log_sd", "rho_raw", "log_shape", "log_range", "log_spatial_sd"))
        need(same_num(as.numeric(a[[nm]]), as.numeric(p[names(p) == nm])),
          paste("Current parameter arrays are misaligned:", nm))
      need(all(vapply(a, function(x) is.numeric(x) && all(is.finite(x)), logical(1))),
        "Non-finite parameter array; checkpoint refused.")
      a
    }
    phase(cfg, "VERIFYING_RECONSTRUCTION", "Replaying likelihood and checking predictions before refinement")
    replay <- withCallingHandlers(obj$fn(cp$par), warning = capture_warning)
    need(is.finite(replay) && abs(replay - cp$objective) < 0.001,
      "Saved negative log-likelihood was not reproduced within 0.001; no refinement attempted.")
    report <- obj$report(obj$env$last.par)
    eta_error <- if (cp$step == 0L) max(abs(as.numeric(report$eta) - original$eta)) else NA_real_
    if (cp$step == 0L) need(is.finite(eta_error) && eta_error < 1e-5,
      "Saved linear predictions were not reproduced; no refinement attempted.")
    write_csv(data.frame(step = cp$step, n_observations = nrow(input$frame),
      source_objective = cp$objective, replayed_objective = replay,
      absolute_objective_error = abs(replay - cp$objective), eta_max_abs_error = eta_error,
      identical_input_signature = TRUE, unchanged_CPP = TRUE),
      file.path(cfg$out, paste0("reconstruction_checks_step_", cp$step, ".csv")))
    cp$objective <- as.numeric(replay); cp$parameters <- current_parameters(cp$par)
    cp$gradient <- withCallingHandlers(as.numeric(obj$gr(cp$par)), warning = capture_warning)
    need(length(cp$gradient) == npar && all(is.finite(cp$gradient)), "Invalid reconstructed gradient.")
    save_cp <- function() {
      disk_check(cfg$out, 5, verbose = FALSE)
      atomic_save(cp, checkpoint_path)
    }
    assess_uncertainty <- function() {
      # Checkpoint deliberately has needs_uncertainty=TRUE until this completes.
      phase(cfg, "UNCERTAINTY", paste("Step", cp$step,
        "| accepted parameters saved; computing a NEW marginal Hessian/covariance"))
      disk_check(cfg$out, cfg$minimum_free_GiB, verbose = FALSE)
      sd <- withCallingHandlers(TMB::sdreport(obj, par.fixed = cp$par,
        getJointPrecision = FALSE, getReportCovariance = FALSE), warning = capture_warning)
      need(is.matrix(sd$cov.fixed) && identical(dim(sd$cov.fixed), as.integer(c(npar, npar))),
        "sdreport returned an unexpected covariance shape. Accepted parameters are retained.")
      cp$covariance <<- sd$cov.fixed
      cp$pdHess <<- isTRUE(sd$pdHess)
      cp$gradient <<- withCallingHandlers(as.numeric(obj$gr(cp$par)), warning = capture_warning)
      val <- withCallingHandlers(obj$fn(cp$par), warning = capture_warning)
      need(is.finite(val) && abs(val - cp$objective) < 1e-5,
        "Objective changed during uncertainty calculation. Parameter checkpoint retained.")
      cp$objective <<- as.numeric(val)
      cp$parameters <<- current_parameters(cp$par)
      cp$needs_uncertainty <<- FALSE
      cp$fresh_covariance <<- TRUE
      # Store arrays only: not sdreport's random covariance or native objective.
      rm(sd); invisible(gc())
      save_cp()
    }
    assessment <- function() {
      cov_ok <- tryCatch({checked_covariance(cp$covariance, npar); TRUE}, error = function(e) FALSE)
      norm <- if (cov_ok && all(is.finite(cp$gradient)))
        gradient_norm(cp$gradient, cp$covariance) else Inf
      list(covariance_ok = cov_ok, norm = norm,
        pass = isTRUE(cp$pdHess) && cov_ok && isTRUE(cp$fresh_covariance) &&
          !isTRUE(cp$needs_uncertainty) && is.finite(norm) && norm < GRADIENT_TOL &&
          is.finite(cp$objective) && cp$objective <= original$optimizer$objective + 0.001)
    }
    add_history <- function(note) {
      a <- assessment()
      row <- data.frame(step = cp$step, objective = cp$objective,
        change_from_source = cp$objective - original$optimizer$objective,
        max_abs_gradient = max(abs(cp$gradient)), gradient_covariance_norm = a$norm,
        pdHess = isTRUE(cp$pdHess), fresh_covariance = isTRUE(cp$fresh_covariance),
        numerical_checks_ok = a$pass, note = note, time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
      cp$history <<- rbind(cp$history, row)
      write_csv(cp$history, file.path(cfg$out, "optimization_history.csv"))
      save_cp()
      phase(cfg, "NUMERICAL_CHECK", paste("Step", cp$step, "| scaled gradient",
        signif(a$norm, 6), "| criterion <", GRADIENT_TOL, "|", note))
      a
    }
    if (cp$needs_uncertainty) assess_uncertainty()
    a <- add_history(if (cp$step == 0L) "Replayed source fit" else "Resumed accepted checkpoint")
    if (a$norm < GRADIENT_TOL && !cp$fresh_covariance) {
      cp$needs_uncertainty <- TRUE; save_cp(); assess_uncertainty()
      a <- add_history("Fresh covariance at the original parameter point")
    }
    stop_reason <- "STEP_LIMIT_REACHED"
    while (!a$pass && cp$step < cfg$max_steps) {
      if (!a$covariance_ok || !isTRUE(cp$pdHess)) {stop_reason <- "COVARIANCE_REVIEW"; break}
      step_id <- cp$step + 1L
      disk_check(cfg$out, cfg$minimum_free_GiB, verbose = FALSE)
      phase(cfg, "NEWTON_DIRECTION", paste("Step", step_id, "of at most", cfg$max_steps,
        "| all marginal parameters remain free"))
      trial <- backtrack(cp$par, cp$objective, cp$gradient, cp$covariance,
        evaluate = function(p, j, alpha) {
          phase(cfg, "LINE_SEARCH", paste("Step", step_id, "| trial", j + 1L,
            "of", cfg$max_backtracks + 1L, "| alpha", alpha))
          # Non-finite trial values are rejected. An actual evaluation error is
          # not swallowed (e.g. out-of-memory); it stops with the checkpoint safe.
          val <- withCallingHandlers(obj$fn(p), warning = capture_warning)
          as.numeric(val)
        }, max_backtracks = cfg$max_backtracks)
      h <- cbind(step = step_id, trial$history)
      cp$line_search <- rbind(cp$line_search, h)
      write_csv(cp$line_search, file.path(cfg$out, "line_search_history.csv"))
      if (!trial$accepted) {
        obj$fn(cp$par) # Restore the accepted base point after rejected trials.
        cp$parameters <- current_parameters(cp$par); save_cp()
        stop_reason <- "NO_ACCEPTABLE_DESCENT"; break
      }
      cp$par <- trial$par; cp$objective <- trial$objective; cp$step <- step_id
      cp$parameters <- current_parameters(cp$par)
      cp$needs_uncertainty <- TRUE; cp$fresh_covariance <- FALSE
      cp$covariance <- NULL; cp$gradient <- NULL; cp$pdHess <- NA
      save_cp()
      atomic_save(cp, file.path(cfg$out, sprintf("accepted_step_%02d_parameters.rds", cp$step)))
      phase(cfg, "ACCEPTED_CHECKPOINT_SAVED", paste("Step", cp$step,
        "| objective", format(cp$objective, digits = 15), "| alpha", trial$alpha))
      assess_uncertainty()
      a <- add_history(paste("Accepted damped-Newton step; alpha", trial$alpha))
    }
    a <- assessment()
    if (a$pass) stop_reason <- "SCALED_GRADIENT_CRITERION_MET"
    status <- if (a$pass) "FIT_OK" else "NUMERICAL_REVIEW_REQUIRED"
    phase(cfg, "EXPORTING", paste(status, "|", stop_reason, "| source fit remains untouched"))
    obj$fn(cp$par); cp$parameters <- current_parameters(cp$par)
    report <- obj$report(obj$env$last.par)
    beta <- stats::setNames(as.numeric(cp$par[beta_ids]), colnames(data$X))
    Vbeta <- cp$covariance[beta_ids, beta_ids, drop = FALSE]
    dimnames(Vbeta) <- list(names(beta), names(beta))
    extent <- max(diff(range(input$sites$x_km)), diff(range(input$sites$y_km)))
    code <- if (a$pass) 0L else 1L
    checks <- data.frame(convergence_code = code, pdHess = isTRUE(cp$pdHess),
      max_abs_gradient = max(abs(cp$gradient)), gradient_covariance_norm = a$norm,
      numerical_checks_ok = a$pass, logLik = -cp$objective, parameter_df = npar,
      AIC = 2 * cp$objective + 2 * npar, nobs = nrow(input$frame),
      spatial_range_km = as.numeric(report$range), spatial_sd = as.numeric(report$spatial_sd),
      range_greater_than_1_5_extent = as.numeric(report$range) > 1.5 * extent,
      spatial_sd_below_0_001 = as.numeric(report$spatial_sd) < 0.001,
      AR1_rho = as.numeric(report$rho), Gamma_shape = as.numeric(report$shape))
    optimizer <- list(par = cp$par, objective = cp$objective, convergence = code,
      message = stop_reason, method = "custom damped Newton, covariance-scaled gradient < 0.05",
      iterations = cp$step, source_optimizer = original$optimizer)
    fit <- list(signature = original$signature, data_signature = input$data_signature,
      parameters = cp$parameters, optimizer = optimizer, beta = beta, V = Vbeta,
      cov_fixed = cp$covariance, checks = checks, eta = as.numeric(report$eta),
      field_at_sites = report$spatial, group = "univoltine", variant = "full",
      fit_minutes = as.numeric(difftime(Sys.time(), cfg$launched_at, units = "mins")),
      warnings = unique(warnings), offset_adjusted_for_onset = TRUE,
      residual_diagnostics = "PENDING", refinement_signature = cfg$signature,
      source_fit = file.path(cfg$source_run, "univoltine__full", "fit.rds"),
      source_fit_minutes = original$fit_minutes, refinement_reason = stop_reason)
    dest <- file.path(cfg$out, "univoltine__full")
    dir.create(dest, showWarnings = FALSE)
    atomic_save(fit, file.path(dest, "fit.rds"))
    atomic_save(list(signature = cfg$signature, data_signature = input$data_signature,
      parameters = cp$parameters, optimizer = optimizer), file.path(dest, "optimized_parameters.rds"))
    write_csv(checks, file.path(dest, "fit_checks.csv"))
    write_csv(fixed_table(beta, Vbeta), file.path(dest, "fixed_effects.csv"))
    write_csv(data.frame(term = names(beta), Vbeta, check.names = FALSE),
      file.path(dest, "fixed_effect_covariance.csv"))
    before <- fixed_table(original$beta, original$V); after <- fixed_table(beta, Vbeta)
    comparison <- merge(before, after, by = "term", sort = FALSE, all = TRUE,
      suffixes = c("_before", "_after"))
    comparison$estimate_change <- comparison$estimate_after - comparison$estimate_before
    comparison$change_in_before_SE <- comparison$estimate_change / comparison$std.error_before
    write_csv(comparison, file.path(cfg$out, "coefficient_comparison.csv"))
    cols <- intersect(names(original$checks), names(checks))
    write_csv(rbind(cbind(stage = "source", original$checks[cols]),
                    cbind(stage = "refined", checks[cols])),
      file.path(cfg$out, "fit_checks_before_after.csv"))
    write_csv(data.frame(index = seq_along(cp$par), term = names(cp$par),
      before = as.numeric(original$optimizer$par), after = as.numeric(cp$par),
      change = as.numeric(cp$par - original$optimizer$par),
      gradient_after = cp$gradient), file.path(cfg$out, "all_marginal_parameter_changes.csv"))
    if (a$pass) write_csv(predictions(beta, Vbeta, input$years),
      file.path(dest, "predictions_spatial.csv"))
    write_csv(data.frame(SITE_ID = rep(input$sites$SITE_ID, times = length(input$years)),
      YEAR = rep(input$years, each = nrow(input$sites)), field = as.vector(report$spatial)),
      file.path(dest, "annual_spatial_field.csv"))
    writeLines(unique(warnings), file.path(dest, "warnings.txt"))
    used <- utils::read.csv(file.path(cfg$out, "source_manifest.csv"), stringsAsFactors = FALSE)
    need(identical(unname(tools::md5sum(used$path)), as.character(used$md5)),
      "Source files changed during refinement. Do not use these exports without review.")
    result <- list(signature = cfg$signature, status = status, reason = stop_reason,
      started_at = cfg$launched_at, finished_at = Sys.time(), output = cfg$out,
      fit_file = file.path(dest, "fit.rds"), accepted_steps = cp$step,
      scaled_gradient_before = original$checks$gradient_covariance_norm,
      scaled_gradient_after = a$norm, residual_diagnostics = "PENDING")
    atomic_save(result, file.path(cfg$out, "result.rds"))
    write_csv(data.frame(status = status, reason = stop_reason, accepted_steps = cp$step,
      threshold = GRADIENT_TOL, gradient_before = result$scaled_gradient_before,
      gradient_after = a$norm, objective_change = cp$objective - original$optimizer$objective,
      elapsed_minutes_this_execution = fit$fit_minutes, residual_diagnostics = "PENDING",
      fit_file = result$fit_file), file.path(cfg$out, "refinement_summary.csv"))
    writeLines(capture.output(utils::sessionInfo()), file.path(cfg$out, "sessionInfo.txt"))
    phase(cfg, status, paste("Scaled gradient", signif(a$norm, 6), "|", stop_reason,
      "| residual diagnostics still pending"), terminal = TRUE)
    archive_review(cfg)
    result
  }

  get_job <- function() {
    j <- state$job
    if (is.null(j)) j <- getOption("phenoimpact.univoltine_refinement.job")
    j
  }
  lock_owner_alive <- function(lock) {
    path <- file.path(lock, "owner.rds")
    owner <- tryCatch(readRDS(path), error = function(e) NULL)
    need(!is.null(owner) && !is.null(owner$pid) && !is.null(owner$create_time),
      paste("Unidentified queue lock; inspect before rerunning:", lock))
    tryCatch(ps::ps_is_running(ps::ps_handle(owner$pid, time = owner$create_time)),
      error = function(e) FALSE)
  }
  release_lock <- function(cfg) {
    path <- file.path(cfg$lock, "owner.rds")
    owner <- if (file.exists(path)) tryCatch(readRDS(path), error = function(e) NULL) else NULL
    if (!is.null(owner) && identical(owner$token, cfg$lock_token)) unlink(cfg$lock, recursive = TRUE)
    invisible(NULL)
  }
  worker_entry <- function(cfg) {
    on.exit(release_lock(cfg), add = TRUE)
    # Parent finishes registering the child's PID before any expensive work.
    ready <- file.path(cfg$out, "launch_ready.rds")
    start <- Sys.time()
    repeat {
      r <- if (file.exists(ready)) tryCatch(readRDS(ready), error = function(e) NULL) else NULL
      if (!is.null(r) && identical(r$token, cfg$lock_token)) break
      need(as.numeric(difftime(Sys.time(), start, units = "secs")) < 60,
        "Parent did not finalize background launch.")
      Sys.sleep(0.1)
    }
    tryCatch(refine_worker(cfg), error = function(e) {
      writeLines(conditionMessage(e), file.path(cfg$out, "error.txt"))
      try(phase(cfg, "STOPPED_WITH_ERROR", conditionMessage(e), terminal = TRUE), silent = TRUE)
      try(archive_review(cfg), silent = TRUE)
      message("STOPPED: ", conditionMessage(e), "\nSaved source and accepted checkpoints are retained.")
      list(status = "STOPPED_WITH_ERROR", error = conditionMessage(e), output = cfg$out)
    })
  }
  dump_engine <- function(out) {
    e <- environment(dump_engine)
    all <- ls(e, all.names = TRUE)
    fun <- all[vapply(all, function(n) is.function(get(n, envir = e)), logical(1))]
    path <- file.path(out, "refinement_engine.R")
    dump(c("VERSION", "SOURCE_VERSION", "DEFAULT_SOURCE", "CPP_MD5", "CORE_VERSIONS",
      "GRADIENT_TOL", fun), file = path, envir = e, control = "all")
    normalizePath(path, winslash = "/", mustWork = TRUE)
  }
  launch <- function(cfg, monitor = FALSE) {
    check_packages(); disk_check(cfg$out, cfg$minimum_free_GiB)
    if (dir.exists(cfg$lock)) {
      need(!lock_owner_alive(cfg$lock), paste("An abundance/refinement queue is already running:", cfg$lock))
      unlink(cfg$lock, recursive = TRUE)
    }
    need(dir.create(cfg$lock, recursive = TRUE), "Cannot acquire the abundance queue lock.")
    cfg$lock_token <- paste(Sys.getpid(), format(Sys.time(), "%Y%m%d%H%M%OS6"), sep = "_")
    cfg$launched_at <- Sys.time(); cfg$libpath <- .libPaths()
    h <- ps::ps_handle()
    atomic_save(list(pid = Sys.getpid(), create_time = ps::ps_create_time(h), token = cfg$lock_token),
      file.path(cfg$lock, "owner.rds"))
    success <- FALSE; child <- NULL
    on.exit(if (!success) {
      if (!is.null(child) && child$is_alive()) child$kill()
      release_lock(cfg)
    }, add = TRUE)
    cfg$log_file <- file.path(cfg$out, paste0("refinement_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
    atomic_save(cfg, file.path(cfg$out, "refinement_config.rds"))
    phase(cfg, "QUEUED", "Only univoltine refinement; source model and multivoltine untouched")
    child <- callr::r_bg(function(engine, cfg) {
      e <- new.env(parent = globalenv()); sys.source(engine, envir = e)
      e$worker_entry(cfg)
    }, args = list(engine = cfg$engine, cfg = cfg), libpath = cfg$libpath,
      stdout = cfg$log_file, stderr = "2>&1", supervise = TRUE,
      user_profile = FALSE, system_profile = FALSE, wd = cfg$source_run,
      env = c(callr::rcmd_safe_env(), OMP_NUM_THREADS = "1",
        OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1"))
    h <- ps::ps_handle(child$get_pid())
    atomic_save(list(pid = child$get_pid(), create_time = ps::ps_create_time(h), token = cfg$lock_token),
      file.path(cfg$lock, "owner.rds"))
    atomic_save(list(token = cfg$lock_token), file.path(cfg$out, "launch_ready.rds"))
    job <- list(process = child, out = cfg$out, started = cfg$launched_at, log = cfg$log_file)
    state$job <- job; options(phenoimpact.univoltine_refinement.job = job)
    success <- TRUE
    message("Univoltine refinement launched. At most ", cfg$max_steps, " accepted Newton steps.")
    message("Original fit/extraction/multivoltine are read-only. Output: ", cfg$out)
    message("Keep R open. Progress: pheno_refine_univoltine$status()")
    if (monitor) watch(job)
    invisible(job)
  }
  run <- function(source_run = DEFAULT_SOURCE, max_steps = 2L, max_backtracks = 6L,
                  minimum_free_GiB = 25, monitor = FALSE) {
    current <- get_job()
    if (!is.null(current) && current$process$is_alive()) {
      message("This refinement is already running."); return(invisible(current))
    }
    check_packages(); self_test()
    source_run <- normalizePath(source_run, winslash = "/", mustWork = TRUE)
    need(dir.exists(source_run), "Source run directory not found.")
    for (pair in list(c(max_steps, 1, 3), c(max_backtracks, 0, 8)))
      need(is.numeric(pair) && length(pair) == 3L && all(is.finite(pair)) &&
        pair[1] == floor(pair[1]) && pair[1] >= pair[2] && pair[1] <= pair[3],
        "Use max_steps=1..3 and max_backtracks=0..8 (default 2 and 6).")
    need(length(minimum_free_GiB) == 1L && is.numeric(minimum_free_GiB) &&
      is.finite(minimum_free_GiB) && minimum_free_GiB >= 10, "minimum_free_GiB must be >=10.")
    need(is.logical(monitor) && length(monitor) == 1L && !is.na(monitor), "Invalid monitor.")
    src_cfg <- readRDS(file.path(source_run, "run_config.rds"))
    need(identical(src_cfg$version, SOURCE_VERSION) &&
      identical(normalizePath(src_cfg$out, winslash = "/", mustWork = TRUE), source_run),
      "Unexpected or moved source run. No automatic search for other fits is used.")
    need(is.character(src_cfg$lock) && length(src_cfg$lock) == 1L, "Source queue lock is missing.")
    if (dir.exists(src_cfg$lock)) need(!lock_owner_alive(src_cfg$lock),
      "The original abundance queue is active; do not run concurrently.")
    out <- file.path(source_run, "refinement_univoltine",
      paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid()))
    need(!dir.exists(out), "Refinement output already exists; use resume() instead.")
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    cfg <- list(version = VERSION, signature = NULL, source_run = source_run,
      out = normalizePath(out, winslash = "/", mustWork = TRUE),
      lock = src_cfg$lock, max_steps = as.integer(max_steps),
      max_backtracks = as.integer(max_backtracks), minimum_free_GiB = minimum_free_GiB)
    cfg$engine <- dump_engine(cfg$out)
    writeLines(c("BOUNDED REFINEMENT OF THE UNIVOLTINE SPATIAL ABUNDANCE MODEL",
      paste("Source run:", cfg$source_run), paste("Code version:", VERSION),
      "Only the univoltine model is reconstructed and refined; no other fit is loaded.",
      "Same data, covariates, likelihood, priors, AR1 structure, mesh and parameterization.",
      "Direction = -V %*% g using all 15 marginal parameters, not the beta-only covariance.",
      "Default maximum: 2 accepted steps, each with up to 6 halvings of its step length.",
      "Line search requires strict objective decrease and the Armijo condition (c1=1e-4).",
      "A local covariance-metric step-length cap of 1 is a numerical safeguard only.",
      "The full marginal covariance is recomputed after each accepted step.",
      "Pass criterion: scaled gradient <0.05, positive Hessian/covariance, finite estimates and non-worse objective.",
      "Convergence code 0 is from this custom Newton solver, NOT a new nlminb/BFGS return code.",
      "Source optimizer status is retained in optimizer$source_optimizer.",
      "No ridge, changed variance constraints, loosened criterion or automatic extra optimizers.",
      "A bounded number of steps is NOT a wall-time guarantee; Hessian calculations can be slow.",
      "Accepted portable parameter arrays are saved BEFORE uncertainty calculation.",
      "Resume continues a saved accepted point, recalculating pending covariance before new steps.",
      "Original fit.rds and workflow_status.csv are NOT overwritten. Downstream code must explicitly select the reviewed new fit.",
      "No full joint precision or native TMB environment is serialized.",
      "A numerical pass is not residual validation. NEW residual diagnostics remain PENDING.",
      "The point-estimate plasticity approach is unchanged; its uncertainty is not propagated.",
      "results_to_review.zip contains CSV/text/code exports, not the fitted RDS arrays.",
      "API references:", "https://search.r-project.org/CRAN/refmans/TMB/html/MakeADFun.html",
      "https://search.r-project.org/CRAN/refmans/TMB/html/sdreport.html"), file.path(out, "README.txt"))
    launch(cfg, monitor)
  }
  status <- function(job = get_job()) {
    need(!is.null(job), "No refinement job in this session.")
    p <- tryCatch(readRDS(file.path(job$out, "progress.rds")), error = function(e) NULL)
    alive <- job$process$is_alive()
    end <- if (!alive && !is.null(p$finished)) p$finished else Sys.time()
    label <- if (!alive && !is.null(p$finished)) "Recorded duration" else "Elapsed since launch"
    message(label, ": ", round(as.numeric(difftime(end, job$started, units = "hours")), 2),
      " h | worker running: ", alive)
    if (!is.null(p)) message(p$stage, " | ", p$detail)
    f <- file.path(job$out, "optimization_history.csv")
    if (file.exists(f)) try(print(utils::tail(utils::read.csv(f), 4), row.names = FALSE), silent = TRUE)
    message("Output: ", job$out)
    if (!alive && !identical(job$process$get_exit_status(), 0L))
      message("Worker exit code: ", job$process$get_exit_status(), ". Inspect log: ", job$log)
    invisible(list(alive = alive, progress = p, output = job$out))
  }
  watch <- function(job = get_job(), every = 60) {
    need(!is.null(job), "No refinement job in this session.")
    tryCatch(repeat {
      status(job); if (!job$process$is_alive()) break; Sys.sleep(every)
    }, interrupt = function(e) message("Display paused; the background refinement continues."))
    invisible(job)
  }
  show_log <- function(n = 30L, job = get_job()) {
    need(!is.null(job), "No refinement job in this session.")
    message("Log: ", job$log)
    cat(utils::tail(readLines(job$log, warn = FALSE), n), sep = "\n")
    invisible(job$log)
  }
  resume <- function(out_dir = NULL, monitor = FALSE) {
    job <- get_job()
    if (!is.null(job) && job$process$is_alive()) {
      message("Refinement already running."); return(invisible(job))
    }
    if (is.null(out_dir)) {
      need(!is.null(job), "Supply out_dir when resuming from a new R session.")
      out_dir <- job$out
    }
    cfg <- readRDS(file.path(out_dir, "refinement_config.rds"))
    need(identical(cfg$version, VERSION), "Use the same refinement script version.")
    need(identical(normalizePath(out_dir, winslash = "/", mustWork = TRUE), cfg$out),
      "Refinement directory was moved; do not silently change source paths.")
    need(file.exists(cfg$engine), "Saved refinement_engine.R is missing.")
    launch(cfg, monitor)
  }
  stop_job <- function(job = get_job()) {
    need(!is.null(job), "No refinement job in this session.")
    if (job$process$is_alive()) job$process$kill()
    message("Worker stopped. Source fit and accepted checkpoints retained; unfinished calculations may be lost.")
    invisible(job)
  }
  list(run = run, status = status, watch = watch, log = show_log,
       resume = resume, stop = stop_job, self_test = self_test)
})

message("Functions loaded. Start: refinement_job <- pheno_refine_univoltine$run(monitor = FALSE)")
