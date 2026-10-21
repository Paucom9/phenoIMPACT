# phenoIMPACT: explicit, isolated rescue of offset_multivoltine from final v1.
# source() only defines functions. Rscript without --execute only audits metadata.
# Heavy modes require a fresh Rscript session and refuse active independent jobs.
pheno_env_offset_rescue <- local({
  VERSION <- "offset_multivoltine_rescue_1.0.0"
  SOURCE_VERSION <- "environment_mundlak_by_voltinism_qc_final_15km_2.1_resilient"
  ANALYSIS <- "offset_multivoltine"
  ROOT <- "E:/phenoIMPACT project/code/phenoIMPACT"
  need <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
  stamp <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  elapsed <- function() unname(proc.time()[["elapsed"]])
  path <- function(x) normalizePath(x, winslash = "/", mustWork = TRUE)
  small <- function(file) {
    s <- file.info(file)$size
    need(length(s) == 1L && is.finite(s) && s > 0 && s < 1024^2,
      paste("Missing/oversized metadata:", file))
    readRDS(file)
  }
  settings <- function(overrides = list()) {
    defaults <- list(max_passes = 2L, iter.max = 100L, eval.max = 200L,
      rel.tol = 1e-12, x.tol = 1e-10, repeat_evaluations = 3L,
      gradient_repeat_abs = 1e-6, gradient_repeat_rel = 1e-3,
      objective_repeat_abs = 1e-7, objective_progress_abs = 1e-7,
      objective_progress_rel = 1e-12, gradient_progress_abs = 1e-5,
      gradient_progress_rel = .01, stagnation_patience = 1L,
      maximum_pass_seconds = 3600, maximum_total_seconds = 21600,
      checkpoint_seconds = 120, minimum_free_GiB = 40, minimum_disk_GiB = 40,
      rebuild_objective_tolerance = 1e-5)
    need(is.list(overrides) && (!length(overrides) ||
      (!is.null(names(overrides)) && !anyDuplicated(names(overrides)) &&
       all(names(overrides) %in% names(defaults)))), "Unknown/unnamed settings.")
    z <- utils::modifyList(defaults, overrides)
    for (name in names(z)) need(is.numeric(z[[name]]) && length(z[[name]]) == 1L &&
      is.finite(z[[name]]) && z[[name]] > 0, paste("Invalid setting:", name))
    for (name in c("max_passes", "iter.max", "eval.max", "repeat_evaluations", "stagnation_patience"))
      need(z[[name]] == floor(z[[name]]), paste("Integer required:", name))
    need(z$repeat_evaluations >= 2 && z$max_passes <= 10 &&
      z$gradient_progress_rel < 1 && z$gradient_repeat_rel < 1,
      "Invalid repetition/pass/relative-gradient setting.")
    z
  }
  geo <- function(root) {
    e <- new.env(parent = globalenv())
    sys.source(file.path(root, "geography_mundlak_by_voltinism_qc_final.R"), envir = e)
    g <- environment(e$pheno_geo_mundlak$run)
    need(identical(g$VERSION, "geography_mundlak_by_voltinism_qc_final_15km_1.1.0"),
      "Geographic worker checker 1.1.0 is required.")
    g
  }
  manifest <- function(files) {
    files <- sort(unique(vapply(files, path, character(1))))
    info <- file.info(files)
    need(all(!info$isdir & is.finite(info$size) & info$size > 0), "Invalid source files.")
    hash <- rep(NA_character_, length(files))
    tiny <- info$size < 1024^2
    hash[tiny] <- unname(tools::md5sum(files[tiny]))
    data.frame(path = files, bytes = info$size, modified = as.numeric(info$mtime),
      small_file_md5 = hash, stringsAsFactors = FALSE)
  }
  resources <- function(root) {
    packages <- c("sdmTMB", "TMB", "Matrix", "fmesher", "ps")
    versions <- vapply(packages, function(p) as.character(utils::packageVersion(p)), character(1))
    need(identical(unname(versions["sdmTMB"]), "1.1.0"), "sdmTMB 1.1.0 is required.")
    need(requireNamespace("ps", quietly = TRUE), "ps is required.")
    list(packages = versions, available_ram_GiB = ps::ps_system_memory()$avail / 1024^3,
      available_disk_GiB = ps::ps_disk_usage(root)$available / 1024^3)
  }
  audit <- function(project_root = ROOT, source_run = NULL, control = list(), quiet = FALSE) {
    ctl <- settings(control)
    root <- path(project_root)
    if (is.null(source_run)) source_run <- file.path(root, "output/phenology_plasticity",
      "environment_mundlak_by_voltinism_qc_final_15km/run_de0890c85d4a")
    errors <- character()
    result <- list(version = VERSION, observed_at = stamp(), project_root = root,
      source_run = source_run, analysis = ANALYSIS, mode = "metadata_only",
      ready = FALSE, control = ctl)
    record <- function(expr) tryCatch(force(expr), error = function(e) {
      errors <<- c(errors, conditionMessage(e)); NULL
    })
    record({
      run <- path(source_run)
      g <- geo(root)
      roots <- g$refinement_roots(run)
      d <- file.path(roots$v1, ANALYSIS)
      cfg <- small(file.path(run, "run_config.rds"))
      done <- small(file.path(run, ANALYSIS, "completed.rds"))
      refined <- small(file.path(d, "refinement_completed.rds"))
      need(identical(cfg$version, SOURCE_VERSION) &&
        is.character(cfg$signature) && length(cfg$signature) == 1L &&
        grepl("^[a-f0-9]{64}$", cfg$signature) && identical(done$signature, cfg$signature),
        "Source run/completion identity mismatch.")
      summary <- refined$summary
      need(is.data.frame(summary) && nrow(summary) == 1L &&
        identical(summary$analysis, ANALYSIS) &&
        identical(refined$source_status, done$status) &&
        summary$refinement_status %in% c("NUMERICAL_REVIEW_REQUIRED_AFTER_MAX_PASSES",
          "NUMERICAL_REVIEW_REQUIRED", "FIT_OK_REFINED", "FIT_OK_AFTER_REBUILD"),
        "A completed v1 offset_multivoltine source is required.")
      files <- c(fit = file.path(d, "fit_refined_final.rds"),
        prepared = file.path(run, ANALYSIS, "prepared_input.rds"),
        checks = file.path(d, "fit_checks_refined.csv"),
        fixed_effects = file.path(run, ANALYSIS, "fixed_effects.csv"),
        parameters = file.path(d, "parameters_refined.rds"),
        completed = file.path(d, "refinement_completed.rds"),
        original_completed = file.path(run, ANALYSIS, "completed.rds"),
        run_config = file.path(run, "run_config.rds"))
      result$source_run <- run
      result$files <- files
      result$signature <- cfg$signature
      result$source_status <- summary$refinement_status
      result$manifest <- manifest(c(files, roots$config_files))
      need(file.info(files[["checks"]])$size < 1024^2, "Oversized source checks.")
      checks <- utils::read.csv(files[["checks"]])
      need(nrow(checks) == 1L && is.finite(checks$objective) &&
        is.finite(checks$max_abs_gradient) && is.finite(checks$nobs), "Invalid source checks.")
      result$source_checks <- checks
      result$active_workers <- g$active_jobs(root, run)
      need(!length(result$active_workers), paste("Active independent/refinement workers:",
        paste(result$active_workers, collapse = "; ")))
    })
    record({
      result <- utils::modifyList(result, resources(root))
      need(result$available_ram_GiB >= ctl$minimum_free_GiB, "Insufficient available RAM.")
      need(result$available_disk_GiB >= ctl$minimum_disk_GiB, "Insufficient available disk.")
    })
    record(need(!dir.exists(file.path(root,
      "output/phenology_plasticity/offset_multivoltine_rescue/rescue.lock")),
      "Rescue lock exists; inspect its owner and output before another attempt."))
    result$errors <- unique(errors)
    result$ready <- !length(errors)
    if (!quiet) {
      cat("Offset multivoltine rescue: metadata only; no fit loaded or started.\n")
      cat("Ready:", result$ready, "\n")
      if (!is.null(result$source_checks)) print(result$source_checks)
      if (length(errors)) cat(paste0("- ", unique(errors)), sep = "\n")
    }
    invisible(result)
  }
  write_rds <- function(x, file) {
    tmp <- tempfile(".checkpoint-", dirname(file))
    saveRDS(x, tmp, compress = "gzip")
    on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
    backup <- paste0(file, ".previous")
    if (file.exists(backup)) need(unlink(backup) == 0L, "Cannot rotate old checkpoint.")
    if (file.exists(file)) need(file.rename(file, backup), "Cannot preserve previous checkpoint.")
    if (!file.rename(tmp, file)) {
      if (file.exists(backup)) file.rename(backup, file)
      stop("Cannot finalize checkpoint: ", file, call. = FALSE)
    }
    invisible(file)
  }
  csv <- function(x, file) utils::write.csv(x, file, row.names = FALSE, na = "")
  stop_reason <- function(reason) {
    stop(structure(list(message = reason, call = NULL),
      class = c("offset_rescue_stop", "error", "condition")))
  }
  guard <- function(a, owner = NULL) {
    need(identical(a$manifest, manifest(a$manifest$path)), "Source changed after audit.")
    active <- geo(a$project_root)$active_jobs(a$project_root, a$source_run)
    if (!is.null(owner)) active <- setdiff(active, owner)
    need(!length(active), paste("Active independent/refinement workers:", paste(active, collapse = "; ")))
    invisible(TRUE)
  }
  phase <- function(out, stage, detail = "") {
    write_rds(list(stage = stage, detail = detail, observed_at = stamp()),
      file.path(out, "progress.rds"))
    message(stamp(), " | ", stage, " | ", detail)
  }
  budget <- function(deadline, out) {
    if (file.exists(file.path(out, "STOP_REQUESTED"))) stop_reason("STOP_REQUESTED")
    if (elapsed() >= deadline) stop_reason("TIME_BUDGET_REACHED")
  }
  parameter_table <- function(theta, gradient, lower, upper) {
    labels <- names(theta)
    if (is.null(labels)) labels <- rep("parameter", length(theta))
    data.frame(index = seq_along(theta), parameter = labels,
      id = make.unique(labels), estimate = as.numeric(theta),
      gradient = as.numeric(gradient), absolute_gradient = abs(gradient),
      lower = lower, upper = upper,
      distance_to_lower = as.numeric(theta) - lower,
      distance_to_upper = upper - as.numeric(theta), stringsAsFactors = FALSE)
  }
  repeat_diagnostics <- function(obj, theta, lower, upper, ctl, tick = function() NULL) {
    n <- length(theta)
    need(n > 0L && all(is.finite(theta)) && length(lower) == n && length(upper) == n &&
      !anyNA(lower) && !anyNA(upper) && all(theta >= lower & theta <= upper),
      "Invalid outer parameter vector/bounds.")
    values <- numeric(ctl$repeat_evaluations)
    gradients <- matrix(NA_real_, ctl$repeat_evaluations, n)
    timing <- numeric(ctl$repeat_evaluations)
    for (i in seq_len(ctl$repeat_evaluations)) {
      tick(); began <- elapsed()
      values[[i]] <- obj$fn(theta)
      tick()
      gradients[i, ] <- obj$gr(theta)
      timing[[i]] <- elapsed() - began
      tick()
    }
    need(all(is.finite(values)) && all(is.finite(gradients)), "Non-finite objective/gradient.")
    spread <- apply(gradients, 2L, function(x) diff(range(x)))
    noise <- diff(range(values))
    g <- gradients[nrow(gradients), ]
    tol_f <- ctl$objective_repeat_abs + 128 * .Machine$double.eps * max(1, abs(values))
    stable <- noise <= tol_f &&
      all(spread <= ctl$gradient_repeat_abs + ctl$gradient_repeat_rel * abs(g))
    tab <- parameter_table(theta, g, lower, upper)
    tab$term <- if (length(obj$parameter_labels) == n) obj$parameter_labels else tab$id
    tab$gradient_spread <- spread
    list(theta = theta, objective = tail(values, 1), gradient = g,
      max_gradient = max(abs(g)), stable = stable, objective_noise = noise,
      parameters = tab, evaluations = data.frame(repetition = seq_along(values),
        objective = values, max_abs_gradient = apply(abs(gradients), 1, max),
        seconds = timing), gradients = gradients)
  }
  decision <- function(before, after, convergence, ctl) {
    tolerance <- max(ctl$objective_repeat_abs, before$objective_noise, after$objective_noise,
      128 * .Machine$double.eps * max(1, abs(before$objective)))
    f_gain <- before$objective - after$objective
    g_gain <- before$max_gradient - after$max_gradient
    f_progress <- f_gain > max(tolerance, ctl$objective_progress_abs,
      ctl$objective_progress_rel * max(1, abs(before$objective)))
    g_progress <- g_gain > max(ctl$gradient_progress_abs,
      ctl$gradient_progress_rel * before$max_gradient)
    acceptable_objective <- f_gain >= -tolerance
    candidate <- after$stable && acceptable_objective &&
      (f_progress || g_progress || after$max_gradient <= before$max_gradient)
    qualifies <- candidate && isTRUE(convergence == 0L) && after$max_gradient < .001
    reason <- if (!after$stable) "UNSTABLE_GRADIENT" else if (!acceptable_objective)
      "OBJECTIVE_WORSE" else if (!candidate) "GRADIENT_WORSE_WITHOUT_OBJECTIVE_PROGRESS" else
      if (qualifies) "GRADIENT_TARGET_REACHED" else if (!f_progress && !g_progress)
      "STAGNATION" else "PROGRESS"
    list(reason = reason, retain = candidate, qualifies = qualifies,
      stagnated = !f_progress && !g_progress, objective_gain = f_gain,
      gradient_gain = g_gain, objective_tolerance = tolerance)
  }
  save_diagnostics <- function(d, out, name) {
    csv(d$parameters, file.path(out, paste0(name, "_gradient.csv")))
    csv(d$evaluations, file.path(out, paste0(name, "_repeatability.csv")))
    write_rds(d, file.path(out, paste0(name, "_diagnostics.rds")))
  }
  checkpoint <- function(obj, theta, a, label, diagnostics = NULL, optimizer = NULL) {
    joint <- obj$env$last.par
    random <- obj$env$random
    fixed <- if (length(random)) joint[-random] else joint
    need(length(fixed) == length(theta) && all(is.finite(joint)) &&
      max(abs(as.numeric(fixed) - theta)) < 1e-8, "TMB current state differs from requested parameters.")
    list(version = VERSION, analysis = ANALYSIS, source_signature = a$signature,
      source_manifest = a$manifest, label = label, observed_at = stamp(),
      theta = theta, joint_parameters = joint,
      parlist = obj$env$parList(par = joint), diagnostics = diagnostics, optimizer = optimizer,
      numerically_accepted = FALSE)
  }
  validate_source <- function(saved, p, a) {
    need(identical(saved$version, "environment_mundlak_wait_refine_1.0") &&
      identical(saved$refinement_status, a$source_status) &&
      identical(saved$analysis, ANALYSIS) && identical(saved$source_signature, a$signature) &&
      identical(path(saved$source_run), a$source_run) &&
      identical(p$version, SOURCE_VERSION) && identical(p$analysis, ANALYSIS) &&
      identical(p$spec$voltinism, "multivoltine") &&
      identical(p$spec$anomaly, "clim_anomaly_tw90") &&
      identical(p$spec$control, "ONSET_mean_z"), "Saved/prepared source identity mismatch.")
    fit <- saved$fit
    need(inherits(fit, "sdmTMB") && as.character(fit$version) == "1.1.0" &&
      isFALSE(fit$reml) && identical(fit$time, "YEAR") && all(fit$spatial == "off") &&
      all(fit$spatiotemporal == "iid") && identical(fit$family$family, "gaussian") &&
      identical(fit$family$link, "identity"), "Unexpected source model.")
    need(identical(fit$priors, sdmTMB::sdmTMBpriors()) &&
      (is.null(fit$offset) || all(fit$offset == 0)) &&
      (is.null(fit$tmb_data$weights) || all(fit$tmb_data$weights == 1)) &&
      is.null(fit$time_varying) && is.null(fit$spatial_varying) &&
      is.null(fit$nonlocal_formula) && !isTRUE(as.logical(fit$tmb_data$has_smooths)),
      "Unexpected source priors/weights/varying effects.")
    f <- if (inherits(fit$formula, "formula")) fit$formula else fit$formula[[1L]]
    columns <- unique(c(all.vars(f), "YEAR", "x_km", "y_km",
      if ("source_row" %in% names(p$data)) "source_row"))
    same_rows <- nrow(fit$data) == nrow(p$data) &&
      all(columns %in% names(fit$data)) && all(columns %in% names(p$data)) &&
      all(vapply(columns, function(name) {
        x <- fit$data[[name]]; y <- p$data[[name]]
        if (is.factor(x) || is.factor(y)) identical(as.character(x), as.character(y)) else identical(x, y)
      }, logical(1)))
    need(identical(paste(deparse(f), collapse = " "), paste(deparse(p$formula), collapse = " ")) &&
      same_rows &&
      isTRUE(all.equal(fit$spde, p$mesh, tolerance = 0, check.attributes = TRUE)),
      "Source rows/formula/mesh mismatch.")
    need(nrow(fit$data) == a$source_checks$nobs &&
      identical(p$classification_audit$group, "multivoltine") &&
      isTRUE(p$classification_audit$retained_rows == nrow(p$data)),
      "Source sample/voltinism audit mismatch.")
    invisible(TRUE)
  }
  load_source <- function(a, out, tick) {
    # This function is deliberately unreachable from audit mode.
    phase(out, "LOADING_V1", "Reading the final v1 fit; no source file will be written.")
    tick()
    saved <- readRDS(a$files[["fit"]])
    # Drop the very large report covariance before taking another fit reference.
    saved$fit$sd_report[["cov"]] <- NULL
    saved$fit$sd_report$env <- NULL
    invisible(gc())
    tick()
    p <- readRDS(a$files[["prepared"]])
    validate_source(saved, p, a)
    fit <- saved$fit
    old <- list(model = fit$model, gradients = fit$gradients,
      sdreport_parameters = fit$sd_report$par.fixed,
      sdreport_gradient = fit$sd_report$gradient.fixed,
      sdreport_covariance = fit$sd_report$cov.fixed,
      pdHess = fit$sd_report$pdHess, nlminb_control = fit$nlminb_control)
    write_rds(old, file.path(out, "source_optimizer_and_covariance.rds"))
    saved$fit <- NULL
    phase(out, "REBUILDING", "Same scientific model; do_fit = FALSE.")
    tick()
    rebuilt <- sdmTMB::sdmTMB(formula = p$formula, data = p$data, mesh = p$mesh,
      time = "YEAR", family = stats::gaussian(link = "identity"), spatial = "off",
      spatiotemporal = "iid", reml = FALSE, previous_fit = fit, do_fit = FALSE,
      silent = TRUE, control = sdmTMB::sdmTMBcontrol(multiphase = FALSE,
        newton_loops = 0L, nlminb_loops = 1L, parallel = 1L,
        get_joint_precision = FALSE, collapse_spatial_variance = FALSE))
    need(identical(names(rebuilt$tmb_obj$par), names(fit$model$par)) &&
      length(rebuilt$tmb_obj$par) == length(fit$model$par) &&
      identical(rebuilt$tmb_map, fit$tmb_map) &&
      isTRUE(all.equal(rebuilt$lower, fit$lower, tolerance = 0)) &&
      isTRUE(all.equal(rebuilt$upper, fit$upper, tolerance = 0)),
      "Rebuilt outer parameters/map/bounds differ from v1.")
    fit$tmb_obj <- rebuilt$tmb_obj
    fixed <- utils::read.csv(a$files[["fixed_effects"]], stringsAsFactors = FALSE)
    labels <- make.unique(names(fit$model$par))
    beta <- which(names(fit$model$par) == "b_j")
    need(length(beta) == nrow(fixed), "Cannot label fixed-effect gradients.")
    labels[beta] <- fixed$term
    fit$tmb_obj$parameter_labels <- labels
    fit$lower <- rep_len(rebuilt$lower, length(fit$model$par))
    fit$upper <- rep_len(rebuilt$upper, length(fit$model$par))
    rm(rebuilt, p); invisible(gc())
    tick()
    list(fit = fit, saved = saved, original = old)
  }
  optimize_once <- function(obj, current, bounds, a, out, pass, ctl, tick,
                            minimize = stats::nlminb) {
    started <- elapsed()
    deadline <- started + ctl$maximum_pass_seconds
    counts <- c(fn = 0L, gr = 0L)
    last_save <- -Inf
    trace_file <- file.path(out, sprintf("pass_%02d_evaluations.csv", pass))
    writeLines("kind,count,seconds,value", trace_file)
    check_time <- function() { tick(); if (elapsed() >= deadline) stop_reason("PASS_TIME_BUDGET_REACHED") }
    fn <- function(theta) {
      check_time(); began <- elapsed()
      value <- obj$fn(theta); counts[["fn"]] <<- counts[["fn"]] + 1L
      cat(sprintf("fn,%d,%.9g,%.17g\n", counts[["fn"]], elapsed() - began, value),
        file = trace_file, append = TRUE)
      need(is.finite(value), "Non-finite optimizer objective.")
      if (elapsed() - last_save >= ctl$checkpoint_seconds) {
        write_rds(checkpoint(obj, theta, a, "LAST_EVALUATED_NOT_ACCEPTED"),
          file.path(out, "evaluation_checkpoint.rds"))
        last_save <<- elapsed()
      }
      check_time(); value
    }
    gr <- function(theta) {
      check_time(); began <- elapsed()
      g <- obj$gr(theta); counts[["gr"]] <<- counts[["gr"]] + 1L
      need(length(g) == length(theta) && all(is.finite(g)), "Invalid optimizer gradient.")
      cat(sprintf("gr,%d,%.9g,%.17g\n", counts[["gr"]], elapsed() - began, max(abs(g))),
        file = trace_file, append = TRUE)
      check_time(); g
    }
    warnings <- character()
    opt <- tryCatch(withCallingHandlers(minimize(current$theta, fn, gr,
        lower = bounds$lower, upper = bounds$upper,
        control = ctl[c("iter.max", "eval.max", "rel.tol", "x.tol")]),
        warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }),
      error = function(err) {
        write_rds(list(error = conditionMessage(err), condition_class = class(err),
          callback_counts = counts, warnings = unique(warnings), seconds = elapsed() - started,
          controls = ctl[c("iter.max", "eval.max", "rel.tol", "x.tol")]),
          file.path(out, sprintf("pass_%02d_optimizer_error.rds", pass)))
        stop(err)
      })
    record <- list(optimizer = opt, warnings = unique(warnings),
      callback_counts = counts, seconds = elapsed() - started,
      controls = ctl[c("iter.max", "eval.max", "rel.tol", "x.tol")])
    write_rds(record, file.path(out, sprintf("pass_%02d_optimizer.rds", pass)))
    record
  }
  rescue_loop <- function(obj, current, bounds, a, out, ctl, tick, optimizer = optimize_once) {
    history <- list()
    opt <- list(par = current$theta, convergence = NA_integer_, objective = current$objective)
    stagnant <- 0L
    status <- "MAX_PASSES_REACHED"
    for (pass in seq_len(ctl$max_passes)) {
      tick(); phase(out, "OPTIMIZING", paste("Pass", pass, "of", ctl$max_passes))
      record <- optimizer(obj, current, bounds, a, out, pass, ctl, tick)
      candidate_opt <- record$optimizer
      after <- repeat_diagnostics(obj, candidate_opt$par, bounds$lower, bounds$upper, ctl, tick)
      save_diagnostics(after, out, sprintf("pass_%02d", pass))
      d <- decision(current, after, candidate_opt$convergence, ctl)
      history[[pass]] <- data.frame(pass = pass, objective = after$objective,
        max_abs_gradient = after$max_gradient, convergence = candidate_opt$convergence,
        optimizer_message = if (is.null(candidate_opt$message)) "" else candidate_opt$message,
        objective_gain = d$objective_gain, gradient_gain = d$gradient_gain,
        max_parameter_change = max(abs(after$theta - current$theta)),
        stable = after$stable, retained = d$retain, reason = d$reason,
        optimization_seconds = record$seconds, stringsAsFactors = FALSE)
      csv(do.call(rbind, history), file.path(out, "rescue_history.csv"))
      write_rds(checkpoint(obj, after$theta, a, d$reason, after, candidate_opt),
        file.path(out, sprintf("pass_%02d_checkpoint.rds", pass)))
      if (!d$retain) { status <- d$reason; break }
      current <- after; opt <- candidate_opt
      write_rds(checkpoint(obj, current$theta, a, "RETAINED_NOT_ACCEPTED", current, opt),
        file.path(out, "retained_checkpoint.rds"))
      if (d$qualifies) { status <- "GRADIENT_TARGET_REACHED"; break }
      stagnant <- if (d$stagnated) stagnant + 1L else 0L
      if (stagnant >= ctl$stagnation_patience) { status <- "STAGNATION"; break }
    }
    list(status = status, current = current, optimizer = opt, history = history)
  }
  final_gate <- function(d, opt, sd, initial_objective, ctl) {
    tol <- max(ctl$objective_repeat_abs, d$objective_noise,
      128 * .Machine$double.eps * max(1, abs(initial_objective)))
    n <- length(d$theta)
    covariance <- sd$cov.fixed
    covariance_ok <- is.matrix(covariance) && identical(dim(covariance), c(n, n)) &&
      all(is.finite(covariance)) &&
      isTRUE(tryCatch({ chol(covariance); TRUE }, error = function(e) FALSE))
    sd_gradient <- as.numeric(sd$gradient.fixed)
    same_point <- is.numeric(sd$par.fixed) && length(sd$par.fixed) == n &&
      identical(names(sd$par.fixed), names(d$theta)) &&
      all(is.finite(sd$par.fixed)) && max(abs(sd$par.fixed - d$theta)) < 1e-10
    consistent <- length(sd_gradient) == n && all(is.finite(sd_gradient)) &&
      all(abs(sd_gradient - d$gradient) <= ctl$gradient_repeat_abs +
        ctl$gradient_repeat_rel * abs(d$gradient))
    okay <- isTRUE(opt$convergence == 0L) && d$stable && d$max_gradient < .001 &&
      same_point && consistent && max(abs(sd_gradient)) < .001 &&
      isTRUE(sd$pdHess) && covariance_ok && d$objective <= initial_objective + tol
    scaled <- if (covariance_ok) sqrt(drop(crossprod(d$gradient, covariance %*% d$gradient))) else NA_real_
    list(accepted = okay, covariance_ok = covariance_ok, same_point = same_point,
      gradients_consistent = consistent, gradient_covariance_norm = scaled,
      raw_gradient_threshold = .001)
  }
  run <- function(project_root = ROOT, source_run = NULL, mode = "audit", control = list()) {
    need(length(mode) == 1L && mode %in% c("audit", "diagnose", "rescue"), "Unknown mode.")
    a <- audit(project_root, source_run, control)
    if (mode == "audit") return(invisible(a))
    need(a$ready, paste("Rescue blocked before loading any fit:", paste(a$errors, collapse = "; ")))
    need(!interactive(), "Use a fresh Rscript --vanilla process for heavy modes.")
    guard(a)
    ctl <- a$control
    parent <- file.path(a$project_root, "output/phenology_plasticity/offset_multivoltine_rescue")
    dir.create(parent, recursive = TRUE, showWarnings = FALSE)
    parent <- path(parent)
    lock <- file.path(parent, "rescue.lock")
    need(dir.create(lock, showWarnings = FALSE), "Rescue lock exists; inspect it before another launch.")
    lock_pid <- Sys.getpid()
    lock_created <- as.numeric(ps::ps_create_time(ps::ps_handle()))
    on.exit({
      # Only delete this attempt's lock, never a source/output tree.
      if (dir.exists(lock) && identical(dirname(path(lock)), parent)) {
        lock_owner <- tryCatch(small(file.path(lock, "owner.rds")), error = function(e) NULL)
        if ((is.null(lock_owner) && !length(list.files(lock, all.files = TRUE, no.. = TRUE))) ||
            (!is.null(lock_owner) && identical(lock_owner$pid, lock_pid) &&
             identical(as.numeric(lock_owner$create_time), lock_created))) unlink(lock, recursive = TRUE)
      }
    }, add = TRUE)
    out <- tempfile(paste0("attempt_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_"), parent)
    need(dir.create(out), "Cannot create a fresh rescue directory.")
    out <- path(out)
    owner <- file.path(out, "worker_owner.rds")
    ownership <- list(pid = Sys.getpid(), create_time = ps::ps_create_time(ps::ps_handle()),
      started = Sys.time(), out = out)
    write_rds(ownership, owner)
    # The directory lock is atomic; the worker owner is recorded only once to avoid self-blocking.
    write_rds(list(pid = ownership$pid, create_time = ownership$create_time, output = out),
      file.path(lock, "owner.rds"))
    write_rds(list(audit = a, mode = mode, version = VERSION), file.path(out, "rescue_config.rds"))
    deadline <- elapsed() + ctl$maximum_total_seconds
    tick <- function() budget(deadline, out)
    result <- tryCatch({
      guard(a, owner)
      need(requireNamespace("sdmTMB", quietly = TRUE) && requireNamespace("TMB", quietly = TRUE),
        "Required native packages unavailable.")
      options(sdmTMB.cores = 1L)
      TMB::openmp(n = 1L, DLL = "sdmTMB")
      data <- load_source(a, out, tick)
      fit <- data$fit
      obj <- fit$tmb_obj
      bounds <- list(lower = fit$lower, upper = fit$upper)
      phase(out, "DIAGNOSING", "Repeated objective/gradient at the same outer parameters.")
      current <- repeat_diagnostics(obj, fit$model$par, bounds$lower, bounds$upper, ctl, tick)
      save_diagnostics(current, out, "initial")
      need(abs(current$objective - fit$model$objective) <= ctl$rebuild_objective_tolerance &&
        abs(current$objective - a$source_checks$objective) <= ctl$rebuild_objective_tolerance,
        "Rebuilt likelihood differs from saved v1.")
      write_rds(checkpoint(obj, current$theta, a, "INITIAL_NOT_ACCEPTED", current, fit$model),
        file.path(out, "initial_checkpoint.rds"))
      if (mode == "diagnose") {
        list(status = if (current$stable) "DIAGNOSED" else "UNSTABLE_GRADIENT",
          accepted = FALSE, max_abs_gradient = current$max_gradient)
      } else if (!current$stable) {
        list(status = "UNSTABLE_GRADIENT", accepted = FALSE)
      } else {
        guard(a, owner)
        initial_objective <- current$objective
        if (current$max_gradient < .001 && isTRUE(fit$model$convergence == 0L)) {
          loop <- list(status = "GRADIENT_TARGET_REACHED", current = current, optimizer = fit$model)
        } else {
          # Recheck all independent workers and source files before each pass.
          checked_optimizer <- function(...) { guard(a, owner); optimize_once(...) }
          loop <- rescue_loop(obj, current, bounds, a, out, ctl, tick, checked_optimizer)
        }
        if (loop$status != "GRADIENT_TARGET_REACHED") {
          list(status = loop$status, accepted = FALSE, max_abs_gradient = loop$current$max_gradient)
        } else {
          tick(); guard(a, owner)
          phase(out, "FINAL_COVARIANCE", "Explicit parameter point; no covariance recomputation during rescue passes.")
          began <- elapsed()
          sd <- TMB::sdreport(obj, par.fixed = loop$current$theta,
            getJointPrecision = FALSE, getReportCovariance = FALSE)
          covariance_seconds <- elapsed() - began
          tick()
          final <- repeat_diagnostics(obj, loop$current$theta, bounds$lower, bounds$upper, ctl, tick)
          save_diagnostics(final, out, "final")
          gate <- final_gate(final, loop$optimizer, sd, initial_objective, ctl)
          write_rds(list(gate = gate, covariance_seconds = covariance_seconds,
            parameters = sd$par.fixed, gradient = sd$gradient.fixed, covariance = sd$cov.fixed,
            pdHess = sd$pdHess), file.path(out, "final_covariance_diagnostics.rds"))
          if (!gate$accepted) list(status = "FINAL_CHECKS_FAILED", accepted = FALSE, checks = gate) else {
            guard(a, owner); tick()
            cp <- checkpoint(obj, final$theta, a, "VALIDATED", final, loop$optimizer)
            fit$model <- loop$optimizer
            fit$model$objective <- final$objective
            fit$sd_report <- sd
            fit$sd_report$env <- NULL
            fit$gradients <- final$gradient
            fit$parlist <- cp$parlist
            fit$last.par.best <- cp$joint_parameters
            fit$pos_def_hessian <- TRUE
            fit$bad_eig <- FALSE
            fit$nlminb_control <- ctl[c("iter.max", "eval.max", "rel.tol", "x.tol")]
            # Do not serialize a live native objective or its potentially stale best-point cache.
            fit$tmb_obj <- NULL
            payload <- list(version = VERSION, analysis = ANALYSIS, source_run = a$source_run,
              source_signature = a$signature, source_manifest = a$manifest,
              source_checkpoint = a$files[["fit"]], refinement_status = "FIT_OK_RESCUED",
              checks = gate, fit = fit)
            write_rds(payload, file.path(out, "fit_rescued_final.rds"))
            csv(data.frame(convergence_code = fit$model$convergence, pdHess = TRUE,
              max_abs_gradient = final$max_gradient, numerical_checks_ok = TRUE,
              objective = final$objective, logLik = -final$objective,
              AIC = 2 * final$objective + 2 * length(final$theta), nobs = nrow(fit$data)),
              file.path(out, "fit_checks_rescued.csv"))
            list(status = "FIT_OK_RESCUED_REVIEW_REQUIRED", accepted = TRUE,
              automatic_downstream_registration = FALSE, checks = gate)
          }
        }
      }
    }, offset_rescue_stop = function(e) list(status = conditionMessage(e), accepted = FALSE),
       error = function(e) list(status = "ERROR", accepted = FALSE, error = conditionMessage(e)))
    result$observed_at <- stamp()
    result$output <- out
    write_rds(result, file.path(out, "result.rds"))
    phase(out, result$status, "Rescue ended; v1/v2 and downstream selection were not modified.")
    print(result)
    invisible(result)
  }
  list(audit = audit, run = run)
})

audit_offset_multivoltine_rescue <- function(...) pheno_env_offset_rescue$audit(...)
run_offset_multivoltine_rescue <- function(...) pheno_env_offset_rescue$run(...)

if (sys.nframe() == 0L && !interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  execute <- "--execute" %in% args
  args <- args[args != "--execute"]
  allowed <- c("mode", "project-root", "source-run", "control-file")
  values <- list()
  for (arg in args) {
    pieces <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    if (!startsWith(arg, "--") || length(pieces) < 2L || !pieces[[1]] %in% allowed ||
        pieces[[1]] %in% names(values)) stop("Invalid/repeated CLI argument: ", arg)
    values[[pieces[[1]]]] <- paste(pieces[-1L], collapse = "=")
  }
  mode <- if (is.null(values$mode)) "audit" else values$mode
  if (mode != "audit" && !execute) stop("Heavy modes require explicit --execute.")
  call <- list(mode = mode)
  if (!is.null(values[["project-root"]])) call$project_root <- values[["project-root"]]
  if (!is.null(values[["source-run"]])) call$source_run <- values[["source-run"]]
  if (!is.null(values[["control-file"]])) {
    file <- values[["control-file"]]
    if (!file.exists(file) || file.info(file)$size >= 1024^2) stop("Invalid control RDS.")
    call$control <- readRDS(file)
  }
  answer <- do.call(run_offset_multivoltine_rescue, call)
  if (mode != "audit" && !answer$status %in% c("DIAGNOSED", "FIT_OK_RESCUED_REVIEW_REQUIRED"))
    quit(status = 2L)
}
