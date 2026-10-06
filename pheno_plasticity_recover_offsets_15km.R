# pheno_plasticity_recover_offsets_15km.R -- 2026-09-27
# In the SAME R session, press Escape once to leave the sensitivity queue.
# The current worker remains alive. Source this file to wait for that worker,
# diagnose the univoltine offset, then recover the multivoltine offset.
# To stop the current sensitivity calculation immediately (losing unsaved work),
# use stop_spatial() AFTER leaving its queue, before sourcing this file.
# Results: output/phenology_plasticity/spatial_models/offset_recovery_15km
# Original model files are never overwritten. No packages are installed/updated.
# This workflow attempts numerical recovery; it cannot guarantee convergence.
# options(pheno.offset.recovery.functions_only = TRUE) disables automatic start.
pheno_offset_recovery <- local({
write_table <- function(x, path) {
    utils::write.csv(x, path, row.names = FALSE, na = "")
}

normal_ci <- function(x) {
    required <- c("term", "estimate", "std.error")
    if (!all(required %in% names(x))) 
        stop("Unexpected fixed-effect table.")
    x <- x[, required]
    names(x)[names(x) == "std.error"] <- "SE"
    z <- stats::qnorm(0.975)
    x$lower_95 <- x$estimate - z * x$SE
    x$upper_95 <- x$estimate + z * x$SE
    x$p_Wald <- 2 * stats::pnorm(-abs(x$estimate/x$SE))
    x$CI_excludes_zero <- x$lower_95 > 0 | x$upper_95 < 0
    x$interval_method <- "Normal Wald 95%; not an LRT"
    tibble::as_tibble(x)
}

fit_checks <- function(fit, cfg) {
    gradients <- fit$gradients
    grad <- if (length(gradients) && all(is.finite(gradients))) {
        max(abs(gradients))
    }
    else NA_real_
    code <- fit$model$convergence
    if (length(code) != 1L) 
        code <- NA_integer_
    pd <- isTRUE(fit$sd_report$pdHess)
    fixed <- tryCatch(sdmTMB::tidy(fit, effects = "fixed", conf.int = FALSE), error = function(e) NULL)
    se_ok <- !is.null(fixed) && all(c("estimate", "std.error") %in% names(fixed)) && nrow(fixed) > 0L && all(is.finite(fixed$estimate)) && 
        all(is.finite(fixed$std.error) & fixed$std.error > 0)
    sanity <- tryCatch(sdmTMB::sanity(fit, gradient_thresh = cfg$gradient_threshold, silent = TRUE), error = function(e) NULL)
    all_sanity <- !is.null(sanity) && isTRUE(all(unlist(sanity)))
    data.frame(convergence_code = code, positive_definite_Hessian = pd, max_abs_gradient = grad, finite_fixed_SE = se_ok, 
        all_sanity_checks_pass = all_sanity, fit_ok = isTRUE(code == 0L) && pd && is.finite(grad) && grad < cfg$gradient_threshold && 
            se_ok)
}

formula_text <- function(f) paste(deparse(f, width.cutoff = 500L), collapse = " ")

atomic_save <- function(object, path) {
    tmp <- tempfile(pattern = "checkpoint_", tmpdir = dirname(path), fileext = ".rds")
    on.exit(unlink(tmp), add = TRUE)
    saveRDS(object, tmp, compress = "gzip")
    backup <- paste0(path, ".previous")
    existed <- file.exists(path)
    if (existed && !file.copy(path, backup, overwrite = TRUE)) 
        stop("Cannot back up ", path)
    if (existed && !file.remove(path)) 
        stop("Cannot replace ", path)
    if (!file.rename(tmp, path)) {
        if (existed) 
            file.copy(backup, path, overwrite = TRUE)
        stop("Cannot finalize checkpoint: ", path)
    }
    if (existed) 
        unlink(backup)
    invisible(path)
}

moderators <- function(w) {
    c(photoperiod = paste0("photo_tw", w), background = paste0("clim_background_tw", w), predictability = paste0("clim_predictability_tw", 
        w), trend = paste0("clim_trend_tw", w))
}

interaction_term <- function(terms, a, b) {
    hit <- intersect(c(paste(a, b, sep = ":"), paste(b, a, sep = ":")), terms)
    if (length(hit) != 1L) 
        stop("Cannot identify interaction: ", a, " x ", b)
    hit
}

fixed_covariance <- function(fit, terms) {
    raw <- stats::vcov(fit)
    candidates <- if (is.matrix(raw) || inherits(raw, "Matrix")) 
        list(raw)
    else raw
    for (v in candidates) {
        if ((is.matrix(v) || inherits(v, "Matrix")) && all(terms %in% rownames(v)) && all(terms %in% colnames(v))) {
            V <- as.matrix(v)[terms, terms, drop = FALSE]
            if (all(is.finite(V))) 
                return(V)
        }
    }
    stop("Cannot align vcov(fit) with fixed-effect coefficient names.")
}

fit_statistics <- function(fit) {
    ll <- stats::logLik(fit)
    data.frame(n_observations = nrow(fit$data), logLik = as.numeric(ll), n_parameters = attr(ll, "df"), AIC = stats::AIC(fit), 
        BIC = stats::BIC(fit))
}

write_primary_results <- function(saved, p, out_dir) {
    fixed <- normal_ci(sdmTMB::tidy(saved$fit, effects = "fixed", conf.int = FALSE))
    if (!setequal(fixed$term, p$original_coefficients$term)) {
        stop("Spatial fixed-effect terms differ from the original model.")
    }
    mods <- moderators(p$best_window)
    interactions <- fixed[match(vapply(mods, function(m) interaction_term(fixed$term, p$anomaly, m), character(1)), 
        fixed$term), ]
    interactions$moderator <- names(mods)
    lrt_path <- file.path(out_dir, "interaction_LRT.csv")
    if (file.exists(lrt_path)) 
        interactions <- dplyr::left_join(interactions, utils::read.csv(lrt_path), by = "term")
    else {
        interactions$p_LRT <- NA_real_
        interactions$LRT_status <- "not_run"
    }
    interactions$LRT_status[is.na(interactions$LRT_status)] <- "not_run"
    random <- sdmTMB::tidy(saved$fit, effects = "ran_pars", conf.int = TRUE)
    comparison <- dplyr::bind_rows(dplyr::mutate(p$original_coefficients, model = "original_lmer"), dplyr::mutate(fixed, 
        model = "annual_spatial"))
    focal <- comparison[comparison$term %in% c(p$anomaly, interactions$term), ]
    changes <- dplyr::mutate(dplyr::inner_join(p$original_coefficients, fixed, by = "term", suffix = c("_original", 
        "_spatial")), estimate_change = estimate_spatial - estimate_original, SE_ratio = SE_spatial/SE_original)
    write_table(fixed, file.path(out_dir, "fixed_effects.csv"))
    write_table(random, file.path(out_dir, "spatial_parameters.csv"))
    write_table(interactions, file.path(out_dir, "main_interactions.csv"))
    write_table(changes, file.path(out_dir, "coefficient_changes.csv"))
    graph <- ggplot2::ggplot(focal, ggplot2::aes(term, estimate, colour = model)) + ggplot2::geom_hline(yintercept = 0, 
        linetype = 2, colour = "grey50") + ggplot2::geom_pointrange(ggplot2::aes(ymin = lower_95, ymax = upper_95), 
        position = ggplot2::position_dodge(width = 0.55)) + ggplot2::coord_flip() + ggplot2::theme_bw() + ggplot2::labs(x = NULL, 
        y = "Estimate and normal Wald 95% CI", colour = "Model")
    ggplot2::ggsave(file.path(out_dir, "coefficient_comparison.png"), graph, width = 12, height = 6, dpi = 250)
    b <- stats::setNames(fixed$estimate, fixed$term)
    V <- fixed_covariance(saved$fit, names(b))
    if (!isTRUE(all.equal(unname(sqrt(diag(V))), fixed$SE, tolerance = 1e-05))) {
        stop("Fixed-effect covariance diagonal and reported standard errors disagree.")
    }
    curves <- lapply(seq_along(mods), function(i) {
        mod <- mods[[i]]
        term <- interaction_term(names(b), p$anomaly, mod)
        limits <- stats::quantile(p$data[[mod]], c(0.025, 0.975), names = FALSE)
        x <- seq(limits[1], limits[2], length.out = 150L)
        slope <- b[[p$anomaly]] + b[[term]] * x
        variance <- V[p$anomaly, p$anomaly] + x^2 * V[term, term] + 2 * x * V[p$anomaly, term]
        if (any(variance < -1e-08)) 
            stop("Negative slope variance.")
        se <- sqrt(pmax(0, variance))
        data.frame(moderator = names(mods)[i], moderator_value = x, slope = slope, SE = se, lower_95 = slope - 
            stats::qnorm(0.975) * se, upper_95 = slope + stats::qnorm(0.975) * se)
    })
    curves <- dplyr::bind_rows(curves)
    write_table(curves, file.path(out_dir, "plasticity_slopes.csv"))
    graph <- ggplot2::ggplot(curves, ggplot2::aes(moderator_value, slope)) + ggplot2::geom_hline(yintercept = 0, 
        linetype = 2, colour = "grey50") + ggplot2::geom_ribbon(ggplot2::aes(ymin = lower_95, ymax = upper_95), 
        fill = "#2166ac", alpha = 0.2) + ggplot2::geom_line(colour = "#2166ac") + ggplot2::facet_wrap(~moderator, 
        scales = "free_x") + ggplot2::theme_bw() + ggplot2::labs(x = "Moderator (original model scale)", y = paste(p$response, 
        "slope (days per anomaly model unit)"), caption = "Population fixed effects; other moderators = 0. Pointwise normal Wald 95% CI.")
    ggplot2::ggsave(file.path(out_dir, "plasticity_slopes.png"), graph, width = 11, height = 7, dpi = 250)
    readme <- data.frame(topic = c("Window", "Model", "Mesh", "Intervals", "Tests", "Slopes"), details = c(paste("Preserved original window:", 
        p$best_window, "days"), "Gaussian ML; site/species intercepts, independent species anomaly slope, site/year intercept, annual IID Matern field.", 
        "15 km construction cutoff; not a correlation range or an optimized resolution; original reference mesh.", 
        "Normal Wald 95% confidence intervals; inference conditional on the selected window and model.", "p_Wald is a normal Wald test; p_LRT is a matched spatial reduced-model LRT, when run.", 
        "Partial anomaly slopes with other moderators at zero; central observed moderator range; covariance included."))
    writexl::write_xlsx(list(README = readme, main_interactions = interactions, fixed_effects = fixed, coefficient_changes = changes, 
        spatial_parameters = random, fit_checks = cbind(saved$checks, fit_statistics(saved$fit)), plasticity_slopes = curves), 
        file.path(out_dir, "results.xlsx"))
    invisible(interactions)
}

# Recovery logic; bundled with the unchanged primary scientific export helpers.

recovery_packages <- function() {
  pkgs <- c("sdmTMB", "TMB", "Matrix", "fmesher", "lme4", "dplyr", "tibble", "ggplot2", "writexl")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))
  invisible(pkgs)
}

set_phase <- function(path, phase, gradient = NA_real_) {
  atomic_save(list(phase = phase, gradient = gradient, updated = Sys.time()), path)
  message("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", phase)
}

same_parameters <- function(a, b) {
  identical(names(a), names(b)) && length(a) == length(b) &&
    isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = 1e-10))
}

recovery_status <- function(checks, details) {
  if (!isTRUE(checks$fit_ok)) return("NUMERICAL_CHECKS_FAILED")
  if (!is.data.frame(details) || !nrow(details) || !isTRUE(all(details$passed))) return("ADDITIONAL_CHECKS_REQUIRE_REVIEW")
  "NUMERICAL_AND_SANITY_CHECKS_OK"
}

diagnose_fit <- function(saved, directory, label, cfg) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  fit <- saved$fit
  checks <- fit_checks(fit, cfg)
  sanity_error <- NULL
  sanity <- tryCatch(sdmTMB::sanity(fit, gradient_thresh = cfg$gradient_threshold, silent = TRUE),
    error = function(e) { sanity_error <<- conditionMessage(e); NULL })
  flat <- unlist(sanity)
  if (!length(flat)) {
    details <- data.frame(check = "sanity_evaluation", passed = FALSE)
  } else {
    details <- data.frame(check = names(flat), passed = vapply(as.list(flat), isTRUE, logical(1)))
  }
  stats <- fit_statistics(fit)
  write_table(cbind(checks, stats), file.path(directory, paste0("fit_checks_", label, ".csv")))
  write_table(details, file.path(directory, paste0("sanity_checks_", label, ".csv")))
  writeLines(c(capture.output(print(details, row.names = FALSE)), sanity_error),
    file.path(directory, paste0("sanity_", label, ".txt")))
  g <- as.numeric(fit$gradients)
  nm <- names(fit$model$par)
  if (length(nm) != length(g)) nm <- paste0("parameter_", seq_along(g))
  terms <- tryCatch(sdmTMB::tidy(fit, effects = "fixed")$term, error = function(e) character())
  if (sum(nm == "b_j") == length(terms)) nm[nm == "b_j"] <- paste0("b_j: ", terms)
  gradients <- data.frame(index = seq_along(g), parameter = nm, gradient = g, abs_gradient = abs(g))
  write_table(gradients[order(-gradients$abs_gradient), ],
    file.path(directory, paste0("gradients_", label, ".csv")))
  pars <- sdmTMB::tidy(fit, effects = "ran_pars", conf.int = TRUE)
  write_table(pars, file.path(directory, paste0("spatial_parameters_", label, ".csv")))
  xy <- fit$spde$xy_cols
  if (length(xy) == 2L && all(xy %in% names(fit$data))) {
    spans <- vapply(fit$data[, xy, drop = FALSE], function(x) diff(range(x)), numeric(1))
    ranges <- pars$estimate[pars$term == "range"]
    if (length(ranges)) write_table(data.frame(
      range_km = ranges, x_span_km = spans[1], y_span_km = spans[2],
      range_to_greatest_axis_span = ranges / max(spans),
      exceeds_sanity_range_rule = ranges > 1.5 * max(spans)),
      file.path(directory, paste0("range_context_", label, ".csv")))
  }
  list(checks = checks, details = details, statistics = stats,
    status = recovery_status(checks, details))
}

verify_inputs <- function(saved, p, meta) {
  if (!inherits(saved$fit, "sdmTMB") || !identical(saved$label, "spatial_field") ||
      !identical(saved$run_id, meta$run_id) || !nzchar(saved$fingerprint)) {
    stop("Checkpoint identity does not match the original metadata.")
  }
  if (!isTRUE(meta$settings$mesh_cutoff_km == 15) ||
      !isTRUE(meta$settings$gradient_threshold == 0.001)) {
    stop("Expected the original 15-km model and 0.001 gradient threshold.")
  }
  if (!identical(formula_text(p$formula), formula_text(meta$formula)) ||
      !all(names(p$data) %in% names(saved$fit$data)) ||
      !isTRUE(all.equal(saved$fit$data[, names(p$data), drop = FALSE], p$data, check.attributes = TRUE))) {
    stop("Saved data or formula differ from the prepared original model.")
  }
  packages <- meta$package_versions
  critical <- c("sdmTMB", "TMB", "Matrix", "fmesher", "lme4")
  versions <- packages[packages$package %in% critical, , drop = FALSE]
  if (!all(critical %in% versions$package)) stop("Original package-version record is incomplete.")
  current <- vapply(versions$package, function(x) as.character(utils::packageVersion(x)), character(1))
  if (!identical(unname(current), as.character(versions$version))) {
    stop("Package versions changed since the saved fit. Restore the original versions before recovery.")
  }
  invisible(TRUE)
}

# Rebuild native TMB state in this fresh process. Never call a serialized native
# function pointer to continue optimization. previous_fit supplies saved estimates.
rebuild_fit <- function(fit, p) {
  rebuilt <- sdmTMB::sdmTMB(formula = p$formula, data = p$data, mesh = fit$spde,
    time = "YEAR", family = stats::gaussian(link = "identity"), spatial = "off",
    spatiotemporal = "iid", reml = FALSE, previous_fit = fit, do_fit = FALSE,
    control = sdmTMB::sdmTMBcontrol(multiphase = FALSE), silent = FALSE)
  if (!identical(names(rebuilt$tmb_obj$par), names(fit$model$par)) ||
      length(rebuilt$tmb_obj$par) != length(fit$model$par) ||
      !isTRUE(all.equal(rebuilt$tmb_data, fit$tmb_data, check.attributes = TRUE)) ||
      !isTRUE(all.equal(rebuilt$tmb_map, fit$tmb_map, check.attributes = TRUE))) {
    stop("Reconstructed TMB model does not match the saved model.")
  }
  value <- rebuilt$tmb_obj$fn(fit$model$par)
  old_value <- fit$model$objective
  if (length(value) != 1L || !is.finite(value) || length(old_value) != 1L ||
      !is.finite(old_value) || abs(value - old_value) > 0.01) {
    stop("Saved and reconstructed likelihoods disagree; no optimization attempted.")
  }
  fit$tmb_obj <- rebuilt$tmb_obj
  fit$lower <- rebuilt$lower
  fit$upper <- rebuilt$upper
  if (any(is.finite(fit$lower)) || any(is.finite(fit$upper))) {
    stop("Finite parameter bounds require a separate recovery strategy; Newton steps were not run.")
  }
  fit$model$objective <- value
  fit
}

newton_step <- function(fit) {
  fit <- sdmTMB::run_extra_optimization(fit, nlminb_loops = 0L, newton_loops = 1L)
  # run_extra_optimization updates sd_report but some releases retain these caches.
  if (!same_parameters(fit$sd_report$par.fixed, fit$model$par)) {
    stop("Optimization and uncertainty calculations refer to different parameter values.")
  }
  fit$pos_def_hessian <- isTRUE(fit$sd_report$pdHess)
  fit$tmb_obj$fn(fit$model$par)
  fit$parlist <- fit$tmb_obj$env$parList(par = fit$tmb_obj$env$last.par.best)
  fit$last.par.best <- fit$tmb_obj$env$last.par.best
  fit
}

recover_worker <- function(source_dir, out_dir, phase_file, analysis, max_steps) {
  recovery_packages()
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  set_phase(phase_file, paste(analysis, "reading saved fit"))
  started <- Sys.time()
  source_path <- file.path(source_dir, "fit_spatial_field.rds")
  saved <- readRDS(source_path)
  meta <- readRDS(file.path(source_dir, "metadata.rds"))
  p <- readRDS(file.path(source_dir, "prepared_input.rds"))
  verify_inputs(saved, p, meta)
  cfg <- meta$settings
  writeLines(c(paste("Analysis:", analysis), paste("Original checkpoint:", source_path),
    paste("Original run:", meta$run_id), "Mesh construction cutoff: 15 km.",
    "Original observations, scales, formula, annual IID fields and ML retained.",
    "Gradient threshold: 0.001. Original files are not overwritten.",
    "Univoltine: detailed diagnosis only. Multivoltine: up to three Newton updates from saved estimates.",
    "Passing numerical checks does not establish absence of residual spatial autocorrelation."),
    file.path(out_dir, "README_recovery.txt"))
  writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
  write_table(meta$validation, file.path(out_dir, "validation.csv"))
  set_phase(phase_file, paste(analysis, "checking original fit"))
  before <- diagnose_fit(saved, out_dir, "original", cfg)
  report <- before
  messages <- character()
  capture_warning <- function(w) {
    messages <<- c(messages, conditionMessage(w)); invokeRestart("muffleWarning")
  }
  recovery_path <- file.path(out_dir, "fit_spatial_field.rds")
  attempts <- 0L
  history_path <- file.path(out_dir, "optimization_history.csv")
  history <- if (file.exists(history_path)) utils::read.csv(history_path) else NULL
  if (analysis == "offset_multivoltine") {
    if (file.exists(recovery_path)) {
      resumed <- readRDS(recovery_path)
      verify_inputs(resumed, p, meta)
      if (!identical(resumed$fingerprint, saved$fingerprint)) stop("Recovery fingerprint mismatch.")
      saved <- resumed
      report <- diagnose_fit(saved, out_dir, "resumed", cfg)
    }
    if (!isTRUE(report$checks$fit_ok)) {
      set_phase(phase_file, paste(analysis, "rebuilding native model from saved estimates"), report$checks$max_abs_gradient)
      saved$fit <- withCallingHandlers(rebuild_fit(saved$fit, p), warning = capture_warning)
      for (i in seq_len(max_steps)) {
        set_phase(phase_file, paste(analysis, "Newton step", i, "of at most", max_steps), report$checks$max_abs_gradient)
        old_par <- saved$fit$model$par
        old_objective <- saved$fit$model$objective
        candidate <- withCallingHandlers(newton_step(saved$fit), warning = capture_warning)
        candidate_objective <- candidate$model$objective
        if (!is.finite(candidate_objective) || candidate_objective > old_objective + 1e-6) {
          stop("Newton update worsened the likelihood; the previous saved checkpoint is retained.")
        }
        changed <- !same_parameters(old_par, candidate$model$par)
        saved$fit <- candidate
        saved$checks <- fit_checks(candidate, cfg)
        saved$warnings <- unique(c(saved$warnings, messages))
        saved$recovery_source <- source_path
        saved$recovery_newton_steps <- if (is.null(saved$recovery_newton_steps)) 1L else saved$recovery_newton_steps + 1L
        saved$recovery_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
        set_phase(phase_file, paste(analysis, "saving recovery checkpoint"), saved$checks$max_abs_gradient)
        atomic_save(saved, recovery_path)
        attempts <- i
        report <- diagnose_fit(saved, out_dir, paste0("newton_", saved$recovery_newton_steps), cfg)
        row <- data.frame(timestamp = as.character(Sys.time()), step = saved$recovery_newton_steps,
          objective_before = old_objective, objective_after = candidate_objective,
          parameters_changed = changed, max_abs_gradient = report$checks$max_abs_gradient,
          fit_ok = report$checks$fit_ok, all_sanity_checks_pass = report$checks$all_sanity_checks_pass)
        history <- rbind(history, row)
        write_table(history, history_path)
        rm(candidate); invisible(gc())
        if (isTRUE(report$checks$fit_ok) || !changed) break
      }
    }
  }
  report <- diagnose_fit(saved, out_dir, "final", cfg)
  status <- report$status
  # Export only after both numerical and additional checks pass. A large-range
  # warning is retained for scientific review; it is not silenced by a refit.
  if (analysis == "offset_multivoltine" && status == "NUMERICAL_AND_SANITY_CHECKS_OK") {
    set_phase(phase_file, paste(analysis, "exporting recovered results"), report$checks$max_abs_gradient)
    saved$checks <- report$checks
    write_primary_results(saved, p, out_dir)
  }
  failed <- report$details$check[!report$details$passed & report$details$check != "all_ok"]
  result <- data.frame(analysis = analysis, status = status, max_abs_gradient = report$checks$max_abs_gradient,
    failed_checks = paste(failed, collapse = "; "), newton_steps_this_run = attempts,
    elapsed_minutes = as.numeric(difftime(Sys.time(), started, units = "mins")),
    source_directory = source_dir, output_directory = out_dir)
  write_table(result, file.path(out_dir, "recovery_status.csv"))
  writeLines(unique(c(saved$warnings, messages)), file.path(out_dir, "warnings_recovery.txt"))
  set_phase(phase_file, paste(analysis, status), report$checks$max_abs_gradient)
  result
}

find_offset <- function(root, analysis) {
  dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[startsWith(basename(dirs), paste0("primary__", analysis, "_"))]
  needed <- c("fit_spatial_field.rds", "metadata.rds", "prepared_input.rds")
  dirs <- dirs[vapply(dirs, function(d) all(file.exists(file.path(d, needed))), logical(1))]
  if (length(dirs) != 1L) stop("Expected exactly one original saved fit for ", analysis,
    "; found ", length(dirs), ". No model started.")
  normalizePath(dirs, winslash = "/", mustWork = TRUE)
}

monitor_job <- function(job, interval = 60) {
  tryCatch({
    repeat {
      state <- tryCatch(readRDS(job$phase_file), error = function(e) NULL)
      phase <- if (is.null(state)) "Starting" else state$phase
      gradient <- if (is.null(state) || !is.finite(state$gradient)) "not yet checked" else format(state$gradient, digits = 4)
      mins <- as.numeric(difftime(Sys.time(), job$started, units = "mins"))
      cat(sprintf("[%s] %s | elapsed %.0f min | last gradient %s | worker running: %s\n",
        format(Sys.time(), "%Y-%m-%d %H:%M:%S"), phase, mins, gradient, job$process$is_alive()))
      if (!job$process$is_alive()) break
      deadline <- Sys.time() + interval
      while (job$process$is_alive() && Sys.time() < deadline) Sys.sleep(1)
    }
    job$process$get_result()
  }, interrupt = function(e) {
    message("Display paused. The worker continues. No following recovery task will start.")
    invisible(NULL)
  })
}

wait_current <- function(job) {
  if (is.null(job) || is.null(job$process)) return(invisible(TRUE))
  tryCatch({
    while (job$process$is_alive()) {
      message("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        "] Waiting for the current worker. No other 30-km model will be launched by this script.")
      deadline <- Sys.time() + 60
      while (job$process$is_alive() && Sys.time() < deadline) Sys.sleep(1)
    }
    tryCatch(job$process$get_result(), error = function(e) message("Previous worker: ", conditionMessage(e)))
  }, interrupt = function(e) stop("Recovery queue paused. The current worker continues.", call. = FALSE))
  invisible(TRUE)
}

run_recovery <- function(root = here::here("output", "phenology_plasticity", "spatial_models"), max_steps = 3L) {
  if (!requireNamespace("callr", quietly = TRUE)) stop("Package callr is required.")
  if (length(max_steps) != 1L || !is.finite(max_steps) || max_steps < 1 || max_steps != as.integer(max_steps)) stop("Invalid max_steps.")
  analyses <- c("offset_univoltine", "offset_multivoltine")
  sources <- setNames(vapply(analyses, function(a) find_offset(root, a), character(1)), analyses)
  # Called only after Escape has unwound the old sensitivity queue.
  wait_current(getOption("pheno.spatial.active_job"))
  wait_current(getOption("pheno.offset.recovery.active_job"))
  out_root <- file.path(root, "offset_recovery_15km")
  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
  out_root <- normalizePath(out_root, winslash = "/", mustWork = TRUE)
  engine <- file.path(out_root, "recovery_engine.R")
  engine_names <- c("write_table", "normal_ci", "fit_checks", "formula_text", "atomic_save", "moderators",
    "interaction_term", "fixed_covariance", "fit_statistics", "write_primary_results", "recovery_packages",
    "set_phase", "same_parameters", "recovery_status", "diagnose_fit", "verify_inputs", "rebuild_fit", "newton_step", "recover_worker")
  dump(engine_names, file = engine, envir = environment(run_recovery), control = "all")
  statuses <- list()
  jobs <- list()
  for (analysis in analyses) {
    out <- file.path(out_root, basename(sources[[analysis]]))
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    monitor <- tempfile(pattern = paste0("monitor_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_"), tmpdir = out)
    dir.create(monitor)
    phase_file <- file.path(monitor, "progress.rds")
    log <- file.path(monitor, "worker.log")
    process <- callr::r_bg(function(engine, args) {
      e <- new.env(parent = globalenv())
      sys.source(engine, envir = e)
      do.call(e$recover_worker, args)
    }, args = list(engine = engine, args = list(source_dir = sources[[analysis]], out_dir = out,
      phase_file = phase_file, analysis = analysis, max_steps = as.integer(max_steps))),
      stdout = log, stderr = "2>&1", supervise = TRUE, wd = getwd())
    job <- list(process = process, phase_file = phase_file, started = Sys.time(), log_file = log, output_root = out_root)
    options(pheno.offset.recovery.active_job = job)
    jobs[[analysis]] <- job
    message("Recovery ", analysis, ". Detailed output: ", log)
    result <- tryCatch(monitor_job(job), error = function(e) {
      data.frame(analysis = analysis, status = "WORKER_ERROR", max_abs_gradient = NA_real_,
        failed_checks = conditionMessage(e), newton_steps_this_run = NA_integer_, elapsed_minutes = NA_real_,
        source_directory = sources[[analysis]], output_directory = out)
    })
    if (process$is_alive()) stop("Recovery queue paused; the active worker continues. No next task started.")
    if (is.null(result)) stop("Worker stopped without a confirmed result; inspect ", log)
    statuses[[analysis]] <- result
    status <- do.call(rbind, statuses)
    write_table(status, file.path(out_root, "recovery_summary.csv"))
    print(result[, c("analysis", "status", "max_abs_gradient", "failed_checks")], row.names = FALSE)
  }
  message("Recovery tasks finished. Review recovery_summary.csv in: ", out_root)
  message("Original fits are preserved. Review the reports before adopting any recovered fit.")
  invisible(list(status = do.call(rbind, statuses), jobs = jobs, output_root = out_root))
}

list(run = run_recovery, watch = function() monitor_job(getOption("pheno.offset.recovery.active_job")))

})
if (!isTRUE(getOption("pheno.offset.recovery.functions_only", FALSE))) {
  offset_recovery_job <- pheno_offset_recovery$run()
}
