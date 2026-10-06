# ============================================================================
# phenoIMPACT | Refit the TWO PRIMARY spatial OFFSET models with onset control
# Version 1.0 | 2026-10-01
#
# source(file.choose(), encoding = "UTF-8")
#
# ONE scientific change: add the annual ONSET_mean as an additive covariate.
# Internally ONSET_mean_z = (ONSET_mean - mean)/SD, separately in each fitted
# subset. With an intercept this is the same model as adding raw ONSET_mean.
# All existing covariates/scales, the selected 90-day window, rows, random terms,
# Gaussian observation family and the ORIGINAL 15-km annual IID field are kept.
#
# This script DOES refit both offset models. It does NOT refit onset/first peak,
# reselect temperature windows, run reduced-model LRTs, run error-distribution
# sensitivity models, re-extract population plasticity, or refit abundance.
# Those downstream steps must use the reviewed, onset-adjusted fits later.
#
# Old fits are INPUTS only. In particular, the recovered multivoltine checkpoint
# is mandatory; the failed original multivoltine fit is never substituted.
# Matched source-row/offset checks precede any new full-data fit. Missing onset
# stops preparation: no imputation, extrapolation or silent removal of rows.
#
# Existing estimated parameters initialize optimization; ALL parameters of the
# new model are re-estimated. A start list is used, NOT previous_fit across
# changed formulas. Public API: sdmTMB v1.1.0 / sdmTMBcontrol(start = ...).
# https://sdmtmb.github.io/sdmTMB/reference/sdmTMBcontrol.html
# https://raw.githubusercontent.com/sdmTMB/sdmTMB/v1.1.0/R/fit.R
#
# The controller and each job are separate R processes; jobs run sequentially.
# Escape pauses monitoring only. Keep R open and the computer awake.
#   pheno_offset_onset$status()
#   pheno_offset_onset$watch()
#   pheno_offset_onset$stop()
# Re-sourcing resumes the SAME run when inputs/configuration are unchanged.
# Completed fits are reused. An incomplete numerical fit is retained for review;
# it is NOT silently refitted on re-sourcing. Interrupted fits with no completed
# fit file must restart from their stored initial parameters.
# ============================================================================

pheno_offset_onset <- local({
  VERSION <- "offset_onset_spatial_15km_1.0"
  ROOT <- "E:/phenoIMPACT project/code/phenoIMPACT"
  ANALYSES <- c("offset_univoltine", "offset_multivoltine")
  CORE <- c("sdmTMB", "TMB", "Matrix", "fmesher", "lme4")
  START_KEYS <- c("b_j", "ln_kappa", "ln_tau_E", "ln_phi",
                  "re_cov_pars", "re_b_pars", "epsilon_st")
  state <- new.env(parent = emptyenv())

  need <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
  text_formula <- function(x) paste(deparse(x, width.cutoff = 500L), collapse = " ")
  csv <- function(x, path) utils::write.csv(x, path, row.names = FALSE, na = "")
  md5 <- function(paths) unname(tools::md5sum(paths))
  same_numbers <- function(x, y, tolerance = 1e-8) {
    length(x) == length(y) && all(is.finite(x)) && all(is.finite(y)) &&
      (length(x) == 0L || max(abs(as.numeric(x) - as.numeric(y))) <= tolerance)
  }
  stamp_time <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  atomic_rds <- function(x, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    tmp <- tempfile(pattern = ".writing_", tmpdir = dirname(path), fileext = ".rds")
    on.exit(unlink(tmp), add = TRUE)
    saveRDS(x, tmp, compress = "gzip")
    # This function is only used inside the NEW output run, never on old inputs.
    bak <- paste0(path, ".previous")
    had <- file.exists(path)
    if (had) {
      need(file.copy(path, bak, overwrite = TRUE), paste("Cannot back up", path))
      need(file.remove(path), paste("Cannot replace", path))
    }
    if (!file.rename(tmp, path)) {
      if (had) file.copy(bak, path, overwrite = TRUE)
      stop("Cannot finalize checkpoint: ", path, call. = FALSE)
    }
    if (had) unlink(bak)
    invisible(path)
  }
  phase <- function(cfg, stage, analysis = "", detail = "") {
    rec <- list(time = Sys.time(), stage = stage, analysis = analysis, detail = detail)
    atomic_rds(rec, file.path(cfg$run, "progress.rds"))
    message("[", stamp_time(), "] ", analysis, " | ", stage, " | ", detail)
  }
  versions <- function() data.frame(package = CORE,
    version = vapply(CORE, function(p) as.character(utils::packageVersion(p)), character(1)))
  check_packages <- function() {
    packages <- c(CORE, "data.table", "digest", "callr", "ps")
    missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
    need(!length(missing), paste("Missing packages:", paste(missing, collapse = ", "),
      "\nUse the original R library. Nothing is installed or updated by this script."))
    need(as.character(utils::packageVersion("sdmTMB")) == "1.1.0",
      "Use the original sdmTMB 1.1.0 library, not the separate error-sensitivity library.")
    invisible(TRUE)
  }
  verify_versions <- function(meta) {
    old <- meta$package_versions
    now <- versions()
    need(all(c("package", "version") %in% names(old)) &&
           all(CORE %in% old$package), "Incomplete original package record.")
    want <- as.character(old$version[match(now$package, old$package)])
    need(identical(now$version, want), paste(
      "Core package versions differ from the original fit. Restore the original library:",
      paste(paste(now$package, want, sep = " = "), collapse = "; ")))
  }
  locate <- function(source_root, analysis) {
    dirs <- list.dirs(source_root, full.names = TRUE, recursive = FALSE)
    dirs <- dirs[startsWith(basename(dirs), paste0("primary__", analysis, "_"))]
    required <- c("prepared_input.rds", "metadata.rds", "fit_spatial_field.rds")
    dirs <- dirs[vapply(dirs, function(d) all(file.exists(file.path(d, required))), logical(1))]
    need(length(dirs) == 1L, paste("Expected exactly one PRIMARY source for", analysis,
      "; found", length(dirs), ". No automatic choice between ambiguous versions."))
    cp <- if (analysis == "offset_multivoltine") {
      file.path(source_root, "offset_recovery_15km", basename(dirs), "fit_spatial_field.rds")
    } else file.path(dirs, "fit_spatial_field.rds")
    need(file.exists(cp), paste("Required checkpoint missing:", cp))
    list(source = normalizePath(dirs, winslash = "/", mustWork = TRUE),
         checkpoint = normalizePath(cp, winslash = "/", mustWork = TRUE))
  }
  fixed_table <- function(fit, formula, data) {
    Xnames <- colnames(stats::model.matrix(lme4::nobars(formula),
      data[seq_len(min(3L, nrow(data))), , drop = FALSE]))
    ix <- which(names(fit$sd_report$par.fixed) == "b_j")
    need(length(ix) == length(Xnames) && !anyDuplicated(Xnames),
      "Cannot align fixed coefficients with the saved formula.")
    need(identical(names(fit$model$par), names(fit$sd_report$par.fixed)) &&
           same_numbers(fit$model$par, fit$sd_report$par.fixed),
      "Optimization and covariance refer to different fixed parameter values.")
    b <- as.numeric(fit$sd_report$par.fixed[ix])
    V <- as.matrix(fit$sd_report$cov.fixed[ix, ix, drop = FALSE])
    need(all(is.finite(b)) && all(is.finite(V)) && all(diag(V) > 0),
      "Missing/non-finite coefficient estimates or covariance.")
    dimnames(V) <- list(Xnames, Xnames)
    se <- sqrt(diag(V)); z <- stats::qnorm(.975)
    list(table = data.frame(term = Xnames, estimate = b, std.error = unname(se),
      conf.low = b - z * se, conf.high = b + z * se,
      p_Wald = 2 * stats::pnorm(-abs(b/se)), row.names = NULL), V = V)
  }
  numerical_checks <- function(fit, fixed = NULL) {
    g <- fit$gradients
    maxg <- if (length(g) && all(is.finite(g))) max(abs(g)) else NA_real_
    code <- if (length(fit$model$convergence) == 1L) fit$model$convergence else NA_integer_
    se_ok <- !is.null(fixed) && nrow(fixed) > 0L &&
      all(is.finite(fixed$estimate)) && all(is.finite(fixed$std.error) & fixed$std.error > 0)
    pd <- isTRUE(fit$sd_report$pdHess)
    data.frame(convergence_code = code, positive_definite_Hessian = pd,
      max_abs_gradient = maxg, finite_fixed_SE = se_ok,
      numerical_ok = isTRUE(code == 0L) && pd && is.finite(maxg) && maxg < .001 && se_ok)
  }
  read_parameters <- function(fit) {
    # parList is an R-level unpacker, NOT a call to a saved native likelihood.
    pars <- fit$tmb_obj$env$parList(par = fit$tmb_obj$env$last.par.best)
    need(all(START_KEYS %in% names(pars)), "Saved starting parameters are incomplete.")
    need(all(vapply(pars[START_KEYS], function(x) is.numeric(x) &&
      all(is.finite(x)), logical(1))), "Non-finite saved starting parameters.")
    pars[START_KEYS]
  }
  add_covariate <- function(old_formula, d, old_b) {
    need("ONSET_mean" %in% names(d) && is.numeric(d$ONSET_mean) &&
           all(is.finite(d$ONSET_mean)), "Every fitted row needs a finite ONSET_mean.")
    centre <- mean(d$ONSET_mean); spread <- stats::sd(d$ONSET_mean)
    need(is.finite(spread) && spread > 0, "ONSET_mean is constant or invalid.")
    d$ONSET_mean_z <- (d$ONSET_mean - centre) / spread
    f <- stats::update.formula(old_formula, . ~ . + ONSET_mean_z)
    environment(f) <- baseenv()
    old_terms <- attr(stats::terms(lme4::nobars(old_formula)), "term.labels")
    new_terms <- attr(stats::terms(lme4::nobars(f)), "term.labels")
    need(setequal(new_terms, c(old_terms, "ONSET_mean_z")), "Unexpected new fixed terms.")
    need(identical(vapply(lme4::findbars(old_formula), text_formula, character(1)),
                   vapply(lme4::findbars(f), text_formula, character(1))),
      "Random-effects specification changed unexpectedly.")
    X0 <- stats::model.matrix(lme4::nobars(old_formula), d)
    X1 <- stats::model.matrix(lme4::nobars(f), d)
    need(nrow(X1) == nrow(d) && all(is.finite(X1)) && qr(X1)$rank == ncol(X1),
      "Incomplete or rank-deficient onset-adjusted model matrix.")
    need(setequal(names(old_b), colnames(X0)) && all(colnames(X0) %in% colnames(X1)),
      "Original fixed coefficients/design cannot be aligned.")
    b1 <- stats::setNames(numeric(ncol(X1)), colnames(X1))
    b1[names(old_b)] <- old_b
    err <- max(abs(as.numeric(X0 %*% old_b[colnames(X0)]) - as.numeric(X1 %*% b1)))
    need(is.finite(err) && err < 1e-8, "Warm-start fixed predictions do not reproduce the old fit.")
    list(data = d, formula = f, beta_start = unname(b1),
      scaling = data.frame(variable = "ONSET_mean", fitted_variable = "ONSET_mean_z",
        centre = centre, SD = spread, raw_units = "day of year", n = nrow(d)),
      design_check = data.frame(n_rows = nrow(d), n_old_terms = ncol(X0),
        n_new_terms = ncol(X1), matrix_rank = qr(X1)$rank,
        max_initial_prediction_difference = err))
  }
  attach_onset <- function(d, table, out) {
    keys <- c("SPECIES", "SITE_ID", "YEAR")
    required <- c(keys, "ONSET_mean", "OFFSET_mean")
    need(all(required %in% names(table)), "Phenology table lacks onset/offset or key columns.")
    tab <- as.data.frame(table[, required, drop = FALSE])
    for (k in c("SPECIES", "SITE_ID")) tab[[k]] <- as.character(tab[[k]])
    yr <- suppressWarnings(as.numeric(as.character(tab$YEAR)))
    need(all(is.finite(yr) & yr == floor(yr)) && !anyNA(tab[keys]),
      "Invalid keys/years in phenology input.")
    tab$YEAR <- as.integer(yr)
    # Identical repeated records are harmless; conflicting records are not.
    tab <- unique(tab)
    dup <- duplicated(tab[keys]) | duplicated(tab[keys], fromLast = TRUE)
    if (any(dup)) {
      csv(tab[dup, , drop = FALSE], file.path(out, "conflicting_phenology_records.csv"))
      stop("Conflicting onset/offset values for a species-site-year; see exported audit.", call. = FALSE)
    }
    req <- data.frame(SPECIES = as.character(d$SPECIES), SITE_ID = as.character(d$SITE_ID),
                      YEAR = as.integer(d$YEAR), original_row = seq_len(nrow(d)))
    joined <- merge(req, tab, by = keys, all.x = TRUE, sort = FALSE)
    joined <- joined[order(joined$original_row), , drop = FALSE]
    need(nrow(joined) == nrow(d) && identical(joined$original_row, seq_len(nrow(d))),
      "Phenology join altered the number/order of fitted rows.")
    need(is.numeric(joined$ONSET_mean) && is.numeric(joined$OFFSET_mean),
      "Onset and offset columns must be numeric.")
    bad <- !is.finite(joined$ONSET_mean) | !is.finite(joined$OFFSET_mean)
    if (any(bad)) {
      csv(joined[bad, , drop = FALSE], file.path(out, "missing_onset_or_offset.csv"))
      stop("Missing onset/offset in matched fitted rows. No rows were dropped and no fit was started.", call. = FALSE)
    }
    delta <- abs(joined$OFFSET_mean - d$OFFSET_mean)
    if (any(delta > 1e-7)) {
      audit <- cbind(joined, original_OFFSET_mean = d$OFFSET_mean, abs_difference = delta)
      csv(audit[delta > 1e-7, , drop = FALSE], file.path(out, "offset_source_mismatch.csv"))
      stop("Phenology CSV differs from the original fitted offset data; review its version.", call. = FALSE)
    }
    if ("ONSET_mean" %in% names(d)) need(same_numbers(d$ONSET_mean, joined$ONSET_mean),
      "Saved onset and phenology CSV disagree.")
    d$ONSET_mean <- joined$ONSET_mean
    csv(data.frame(n_matched = nrow(d), n_missing_onset = 0L,
      max_offset_difference = max(delta), n_onset_after_offset = sum(d$ONSET_mean > d$OFFSET_mean)),
      file.path(out, "onset_join_checks.csv"))
    if (any(d$ONSET_mean > d$OFFSET_mean)) {
      csv(d[d$ONSET_mean > d$OFFSET_mean, c(keys, "ONSET_mean", "OFFSET_mean")],
          file.path(out, "onset_after_offset_records.csv"))
      stop("Some onset dates are after offset. Review the exported records before fitting.", call. = FALSE)
    }
    d
  }
  prepare_worker <- function(cfg, analysis) {
    check_packages()
    out <- file.path(cfg$run, analysis); dir.create(out, recursive = TRUE, showWarnings = FALSE)
    complete <- file.path(out, "preparation_complete.rds")
    if (file.exists(complete)) {
      done <- readRDS(complete); need(identical(done$signature, cfg$signature), "Preparation cache mismatch.")
      need(file.exists(file.path(out, "input.rds")), "Prepared input is missing.")
      return("PREPARED_CACHED")
    }
    phase(cfg, "PREPARING", analysis, "Reading original fitted rows and checkpoint")
    job <- cfg$jobs[[analysis]]
    p <- readRDS(file.path(job$source, "prepared_input.rds"))
    meta <- readRDS(file.path(job$source, "metadata.rds")); verify_versions(meta)
    saved <- readRDS(job$checkpoint); fit <- saved$fit
    need(inherits(fit, "sdmTMB") && identical(saved$run_id, meta$run_id) &&
      identical(saved$label, "spatial_field"), "Original checkpoint identity mismatch.")
    need(isTRUE(meta$settings$mesh_cutoff_km == 15) &&
      isTRUE(p$best_window == 90) && identical(p$response, "OFFSET_mean") &&
      identical(p$anomaly, "clim_anomaly_tw90"), "Unexpected mesh/window/response; nothing was changed.")
    need(identical(text_formula(p$formula), text_formula(meta$formula)) &&
      isTRUE(all.equal(fit$data[, names(p$data), drop = FALSE], p$data)),
      "Original checkpoint and prepared input disagree.")
    expected <- stats::as.formula(paste0("OFFSET_mean ~ clim_anomaly_tw90 * (",
      "photo_tw90 + clim_background_tw90 + clim_predictability_tw90 + clim_trend_tw90)",
      " + (1 | SITE_ID) + (1 | SPECIES) + (0 + clim_anomaly_tw90 | SPECIES_slope)",
      " + (1 | site_year_id)"))
    need(setequal(attr(stats::terms(lme4::nobars(p$formula)), "term.labels"),
                  attr(stats::terms(lme4::nobars(expected)), "term.labels")),
      "Unexpected original fixed terms; this script must only add onset.")
    need(setequal(vapply(lme4::findbars(p$formula), text_formula, character(1)),
                  vapply(lme4::findbars(expected), text_formula, character(1))),
      "Unexpected original random effects.")
    need(identical(fit$family$family, "gaussian") && identical(fit$family$link, "identity") &&
      !isTRUE(fit$reml) && !any(fit$offset != 0), "Unexpected family/link/REML/model offset.")
    need(length(fit$spatial) == 1L && length(fit$spatiotemporal) == 1L &&
      all(fit$spatial %in% c("off", FALSE)) && all(fit$spatiotemporal == "iid"),
      "Expected annual IID spatial fields, not a static field or AR1 field.")
    need(!anyNA(p$data) && !anyDuplicated(p$data$source_row) &&
      !anyDuplicated(p$data[c("SPECIES", "SITE_ID", "YEAR")]), "Invalid original fitted rows.")
    old_fixed <- fixed_table(fit, p$formula, p$data)
    need(numerical_checks(fit, old_fixed$table)$numerical_ok,
      "Original checkpoint fails numerical checks. The recovered multivoltine fit is required.")
    pars <- read_parameters(fit)
    need(same_numbers(pars$b_j, old_fixed$table$estimate, 1e-6), "Stale original parameter list.")
    mesh <- fit$spde
    need(inherits(mesh, "sdmTMBmesh") && nrow(mesh$loc_xy) == nrow(p$data) &&
      same_numbers(as.numeric(as.matrix(mesh$loc_xy)),
                   as.numeric(as.matrix(p$data[c("x_km", "y_km")]))),
      "Saved mesh projection does not match fitted coordinates.")
    old_objective <- fit$model$objective
    rm(saved, fit); invisible(gc())
    phase(cfg, "PREPARING", analysis, "Joining annual onset; no rows may be removed")
    tab <- data.table::fread(cfg$phenology_file,
      select = c("SPECIES", "SITE_ID", "YEAR", "ONSET_mean", "OFFSET_mean"),
      colClasses = list(character = c("SPECIES", "SITE_ID")), data.table = FALSE,
      nThread = 1L, showProgress = FALSE)
    d <- attach_onset(p$data, tab, out); rm(tab)
    added <- add_covariate(p$formula, d,
      stats::setNames(old_fixed$table$estimate, old_fixed$table$term))
    pars$b_j <- added$beta_start
    input <- list(signature = cfg$signature, analysis = analysis, data = added$data,
      formula = added$formula, old_formula = p$formula, mesh = mesh, start = pars,
      onset_scaling = added$scaling, old_fixed = old_fixed$table, old_V = old_fixed$V,
      old_objective = old_objective, source = job, original_run_id = meta$run_id,
      source_rows = p$data$source_row, settings = list(mesh_cutoff_km = 15,
      best_window = 90, spatial = "off", spatiotemporal = "iid", time = "YEAR", reml = FALSE))
    csv(old_fixed$table, file.path(out, "old_fixed_effects.csv"))
    csv(added$scaling, file.path(out, "onset_scaling.csv"))
    csv(added$design_check, file.path(out, "design_checks.csv"))
    csv(versions(), file.path(out, "package_versions.csv"))
    writeLines(c("ORIGINAL:", text_formula(p$formula), "", "NEW:", text_formula(added$formula),
      "", "ONSET_mean_z is annual onset, centred/scaled over these EXACT fitted rows.",
      "All other predictors retain their saved scale. The 90-day window is held fixed.",
      "No observations were removed. The original triangulation and projection are reused."),
      file.path(out, "model_specification.txt"))
    atomic_rds(input, file.path(out, "input.rds"))
    atomic_rds(list(signature = cfg$signature, time = Sys.time()), complete)
    "PREPARED"
  }
  export_tables <- function(fit, input, out) {
    fx <- fixed_table(fit, input$formula, input$data)
    ch <- numerical_checks(fit, fx$table)
    csv(ch, file.path(out, "fit_checks.csv"))
    csv(fx$table, file.path(out, "fixed_effects.csv"))
    csv(data.frame(term = rownames(fx$V), fx$V, check.names = FALSE),
        file.path(out, "fixed_effect_covariance.csv"))
    a <- fx$table[grepl("clim_anomaly_tw90", fx$table$term, fixed = TRUE), , drop = FALSE]
    csv(a, file.path(out, "plasticity_coefficients.csv"))
    comp <- merge(input$old_fixed, fx$table, by = "term", all = TRUE,
      suffixes = c("_without_onset", "_with_onset"), sort = FALSE)
    comp$estimate_change <- comp$estimate_with_onset - comp$estimate_without_onset
    comp$SE_ratio <- comp$std.error_with_onset / comp$std.error_without_onset
    csv(comp, file.path(out, "coefficient_comparison.csv"))
    old_k <- nrow(input$old_fixed); new_k <- nrow(fx$table)
    need(new_k == old_k + 1L, "Expected one additional fixed coefficient.")
    csv(data.frame(n_rows = nrow(input$data), extra_fixed_parameters = new_k - old_k,
      old_objective = input$old_objective, new_objective = fit$model$objective,
      delta_AIC_with_minus_without = 2 * (fit$model$objective - input$old_objective) + 2,
      numerical_ok = ch$numerical_ok), file.path(out, "model_comparison.csv"))
    # Explicit back-transformation of the onset coefficient for interpretation.
    row <- fx$table[fx$table$term == "ONSET_mean_z", , drop = FALSE]
    need(nrow(row) == 1L, "Onset covariate is missing from the fitted model.")
    for (k in c("estimate", "std.error", "conf.low", "conf.high"))
      row[[k]] <- row[[k]] / input$onset_scaling$SD
    row$term <- "ONSET_mean (per day; equivalent coefficient)"
    csv(row, file.path(out, "onset_covariate_per_day.csv"))
    ch
  }
  fit_worker <- function(cfg, analysis) {
    check_packages()
    options(sdmTMB.cores = 1L)
    # v1.1.0's parallel control is ignored; set native threads explicitly.
    TMB::openmp(n = 1L, DLL = "sdmTMB")
    out <- file.path(cfg$run, analysis)
    input <- readRDS(file.path(out, "input.rds"))
    need(identical(input$signature, cfg$signature), "Input signature mismatch.")
    done_path <- file.path(out, "completed.rds")
    final_path <- file.path(out, "fit_with_onset.rds")
    if (file.exists(done_path)) {
      done <- readRDS(done_path)
      need(identical(done$signature, cfg$signature) && file.exists(final_path), "Completed cache mismatch.")
      return(done$status)
    }
    # A fit saved just before an export failure can be re-exported without any refit.
    if (file.exists(final_path)) {
      saved <- readRDS(final_path)
      need(identical(saved$signature, cfg$signature), "Saved fit signature mismatch.")
      phase(cfg, "EXPORTING_CACHED", analysis, "No new optimization")
      ch <- export_tables(saved$fit, input, out)
      need(ch$numerical_ok, "Saved new fit needs numerical review; it was not silently refitted.")
      atomic_rds(list(signature = cfg$signature, status = "FIT_OK", time = Sys.time()), done_path)
      return("FIT_OK")
    }
    need(!file.exists(file.path(out, "fit_initial.rds")),
      "An initial fit is retained without a final checkpoint; review it before another optimization.")
    phase(cfg, "FITTING", analysis, paste(nrow(input$data), "rows; annual onset now included"))
    warnings <- character()
    catcher <- function(w) {
      warnings <<- unique(c(warnings, conditionMessage(w)))
      message("WARNING: ", conditionMessage(w)); invokeRestart("muffleWarning")
    }
    set.seed(cfg$fit_seed)
    t0 <- Sys.time()
    fit <- withCallingHandlers(sdmTMB::sdmTMB(
      formula = input$formula, data = input$data, mesh = input$mesh,
      time = "YEAR", family = stats::gaussian(link = "identity"), spatial = "off",
      spatiotemporal = "iid", reml = FALSE, silent = FALSE,
      control = sdmTMB::sdmTMBcontrol(start = input$start, multiphase = FALSE,
        nlminb_loops = 1L, newton_loops = 0L, parallel = 1L,
        get_joint_precision = TRUE, collapse_spatial_variance = FALSE)), warning = catcher)
    # Save BEFORE exports/checks/extra optimization. Original inputs are untouched.
    record <- function(f, attempt) list(signature = cfg$signature, analysis = analysis,
      formula = input$formula, onset_scaling = input$onset_scaling,
      original_source = input$source, fit = f, attempt = attempt, warnings = warnings,
      elapsed_minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")))
    atomic_rds(record(fit, 0L), file.path(out, "fit_initial.rds"))
    csv(data.frame(warning = warnings), file.path(out, "warnings.csv"))
    need(identical(fit$data$source_row, input$data$source_row) &&
      same_numbers(fit$data$ONSET_mean_z, input$data$ONSET_mean_z), "Fitted rows/covariate changed.")
    ch <- tryCatch(numerical_checks(fit, fixed_table(fit, input$formula, input$data)$table),
      error = function(e) numerical_checks(fit))
    history <- cbind(attempt = 0L, objective = fit$model$objective, ch)
    csv(history, file.path(out, "optimization_history.csv"))
    for (attempt in seq_len(cfg$max_extra_newton)) {
      if (isTRUE(ch$numerical_ok)) break
      need(!any(is.finite(fit$lower)) && !any(is.finite(fit$upper)),
        "Finite optimizer bounds require separate recovery; no unbounded Newton step was attempted.")
      phase(cfg, "EXTRA_OPTIMIZATION", analysis, paste("Newton step", attempt,
        "of at most", cfg$max_extra_newton, "; prior checkpoints retained"))
      next_fit <- tryCatch(withCallingHandlers(sdmTMB::run_extra_optimization(
        fit, nlminb_loops = 0L, newton_loops = 1L), warning = catcher),
        error = function(e) { warnings <<- c(warnings, conditionMessage(e)); NULL })
      if (is.null(next_fit)) break
      # Save the attempt even if worse; keep the better likelihood as the working fit.
      atomic_rds(record(next_fit, attempt), file.path(out, sprintf("fit_newton_%02d.rds", attempt)))
      if (!is.finite(next_fit$model$objective) ||
          next_fit$model$objective > fit$model$objective + .001) break
      fit <- next_fit
      ch <- tryCatch(numerical_checks(fit, fixed_table(fit, input$formula, input$data)$table),
        error = function(e) numerical_checks(fit))
      history <- rbind(history, cbind(attempt = attempt, objective = fit$model$objective, ch))
      csv(history, file.path(out, "optimization_history.csv"))
    }
    # Synchronize point-estimate caches after any extra optimization.
    fit$tmb_obj$fn(fit$model$par)
    fit$parlist <- fit$tmb_obj$env$parList(par = fit$tmb_obj$env$last.par.best)
    fit$last.par.best <- fit$tmb_obj$env$last.par.best
    atomic_rds(record(fit, tail(history$attempt, 1L)), final_path)
    writeLines(unique(warnings), file.path(out, "warnings.txt"))
    phase(cfg, "EXPORTING", analysis, "Coefficients, covariance and comparison with old fit")
    ch <- export_tables(fit, input, out)
    need(ch$numerical_ok, "New fit requires numerical review. All checkpoints are retained; no downstream models started.")
    # Numerical sanity warnings (e.g. very long range) are retained, not hidden.
    sanity <- tryCatch(capture.output(sdmTMB::sanity(fit, gradient_thresh = .001)),
      error = function(e) paste("Sanity export failed:", conditionMessage(e)))
    writeLines(sanity, file.path(out, "sanity.txt"))
    tryCatch(csv(as.data.frame(sdmTMB::tidy(fit, effects = "ran_pars", conf.int = TRUE)),
      file.path(out, "random_parameters.csv")), error = function(e)
        writeLines(conditionMessage(e), file.path(out, "random_parameter_export_error.txt")))
    # Save portable parameter estimates too; no native calls are needed to read these.
    atomic_rds(fit$parlist, file.path(out, "parameters.rds"))
    writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
    atomic_rds(list(signature = cfg$signature, status = "FIT_OK", time = Sys.time()), done_path)
    "FIT_OK"
  }
  self_test <- function() {
    # Pure R checks; not a substitute for executing a full sdmTMB fit.
    d <- data.frame(OFFSET_mean = seq_len(40) + 180,
      ONSET_mean = 90 + (seq_len(40) %% 11), x = sin(seq_len(40)),
      w = cos(seq_len(40)/3), SITE_ID = factor(rep(1:5, 8)))
    f <- OFFSET_mean ~ x * w + (1 | SITE_ID)
    b <- c("(Intercept)" = 200, x = -3, w = 2, "x:w" = .2)
    a <- add_covariate(f, d, b)
    need(ncol(stats::model.matrix(lme4::nobars(a$formula), a$data)) == 5L,
      "Onset addition self-test failed.")
    raw <- stats::model.matrix(~ x*w + ONSET_mean, d)
    std <- stats::model.matrix(lme4::nobars(a$formula), a$data)
    bs <- stats::setNames(a$beta_start, colnames(std)); bs["ONSET_mean_z"] <- 7
    br <- stats::setNames(numeric(ncol(raw)), colnames(raw)); br[names(b)] <- b
    br["ONSET_mean"] <- 7/a$scaling$SD
    br["(Intercept)"] <- b["(Intercept)"] - 7*a$scaling$centre/a$scaling$SD
    need(same_numbers(as.numeric(raw %*% br), as.numeric(std %*% bs)),
      "Raw-vs-standardized onset self-test failed.")
    invisible(TRUE)
  }
  collect <- function(cfg) {
    for (name in c("fit_checks", "coefficient_comparison", "plasticity_coefficients", "onset_covariate_per_day")) {
      rows <- lapply(ANALYSES, function(a) {
        p <- file.path(cfg$run, a, paste0(name, ".csv"))
        if (!file.exists(p)) return(NULL)
        cbind(analysis = a, utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE))
      })
      rows <- Filter(Negate(is.null), rows)
      if (length(rows)) csv(do.call(rbind, rows), file.path(cfg$run, paste0("all_", name, ".csv")))
    }
    # Text/CSV review bundle only. No multi-GB fits are put into the ZIP.
    if (requireNamespace("zip", quietly = TRUE)) tryCatch({
      all <- list.files(cfg$run, recursive = TRUE, full.names = FALSE)
      keep <- grepl("\\.(csv|txt|R|log)$", all, ignore.case = TRUE)
      zip::zipr(file.path(cfg$run, "results_to_review.zip"), all[keep], root = cfg$run,
                include_directories = FALSE)
    }, error = function(e) message("Optional review ZIP failed: ", conditionMessage(e)))
    invisible(cfg$run)
  }
  invoke_worker <- function(cfg, fun, analysis) {
    log <- file.path(cfg$run, analysis, paste0(fun, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
    dir.create(dirname(log), recursive = TRUE, showWarnings = FALSE)
    child <- callr::r_bg(function(engine, cfg, fun, analysis) {
      e <- new.env(parent = globalenv()); sys.source(engine, envir = e)
      do.call(e[[fun]], list(cfg = cfg, analysis = analysis))
    }, args = list(cfg$engine, cfg, fun, analysis), libpath = cfg$libpath,
      stdout = log, stderr = "2>&1", supervise = TRUE,
      user_profile = FALSE, system_profile = FALSE,
      env = c(callr::rcmd_safe_env(), OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1"))
    on.exit(if (child$is_alive()) child$kill(), add = TRUE)
    atomic_rds(list(pid = child$get_pid(), stage = fun, analysis = analysis, log = log),
      file.path(cfg$run, "active_worker.rds"))
    child$wait(); child$get_result()
  }
  controller <- function(cfg) {
    on.exit(unlink(cfg$lock, recursive = TRUE), add = TRUE)
    status <- data.frame(analysis = ANALYSES, preparation = "PENDING", fit = "PENDING", message = "")
    save_status <- function() csv(status, file.path(cfg$run, "workflow_status.csv"))
    save_status()
    tryCatch({
      check_packages(); self_test()
      phase(cfg, "PREFLIGHT", detail = "Specification self-tests passed")
      # Prepare BOTH inputs before starting either long fit.
      for (i in seq_along(ANALYSES)) {
        status$preparation[i] <- "RUNNING"; save_status()
        ans <- tryCatch(invoke_worker(cfg, "prepare_worker", ANALYSES[i]),
          error = function(e) {
            status$preparation[i] <<- "PREPARATION_FAILED"
            status$message[i] <<- conditionMessage(e)
            status$fit[status$fit == "PENDING"] <<- "NOT_STARTED"
            save_status()
            stop(conditionMessage(e), call. = FALSE)
          })
        status$preparation[i] <- ans; save_status()
      }
      if (cfg$prepare_only) {
        phase(cfg, "PREPARED_ONLY", detail = "No models fitted")
        return(invisible(status))
      }
      for (i in seq_along(ANALYSES)) {
        status$fit[i] <- "RUNNING"; save_status()
        ans <- tryCatch(invoke_worker(cfg, "fit_worker", ANALYSES[i]), error = function(e) {
          status$message[i] <<- conditionMessage(e); "REVIEW_OR_ERROR"
        })
        status$fit[i] <- ans; save_status(); collect(cfg)
      }
      ok <- all(status$fit == "FIT_OK")
      phase(cfg, if (ok) "FITS_COMPLETE" else "FINISHED_WITH_ISSUES", detail =
        if (ok) "Review coefficients and redo residual diagnostics before downstream extraction/abundance fits"
        else "Inspect workflow_status.csv and per-model logs; completed fits retained")
      invisible(status)
    }, error = function(e) {
      writeLines(conditionMessage(e), file.path(cfg$run, "controller_error.txt"))
      phase(cfg, "STOPPED", detail = conditionMessage(e))
      message("Stopped before remaining fits: ", conditionMessage(e))
      invisible(status)
    }, finally = {save_status(); collect(cfg)})
  }
  engine_file <- function(run) {
    e <- environment(engine_file)
    n <- ls(e, all.names = TRUE)
    fun <- n[vapply(n, function(x) is.function(get(x, envir = e)), logical(1))]
    # Functions needing state are exported but only run/watch/status use it.
    path <- file.path(run, "workflow_engine.R")
    dump(c("VERSION", "ROOT", "ANALYSES", "CORE", "START_KEYS", fun),
         file = path, envir = e, control = "all")
    normalizePath(path, winslash = "/", mustWork = TRUE)
  }
  status <- function(job = state$job) {
    need(!is.null(job), "No job in this R session.")
    path <- file.path(job$run, "progress.rds")
    p <- if (file.exists(path)) tryCatch(readRDS(path), error = function(e) NULL) else NULL
    message("Elapsed: ", round(as.numeric(difftime(Sys.time(), job$started, units = "hours")), 2),
      " h | controller running: ", job$process$is_alive())
    if (!is.null(p)) message(p$stage, " | ", p$analysis, " | ", p$detail)
    q <- file.path(job$run, "workflow_status.csv")
    if (file.exists(q)) print(utils::read.csv(q, stringsAsFactors = FALSE), row.names = FALSE)
    message("Output: ", job$run)
    invisible(p)
  }
  watch <- function(job = state$job, every = 60) {
    need(!is.null(job), "No job in this R session.")
    tryCatch({
      repeat {
        status(job)
        if (!job$process$is_alive()) break
        Sys.sleep(every)
      }
      exit <- job$process$get_exit_status()
      if (!identical(exit, 0L)) message("Controller exit status: ", exit,
        ". See controller log: ", job$log)
    }, interrupt = function(e) message("Monitoring paused; the queue continues. Use pheno_offset_onset$watch()."))
    invisible(job)
  }
  stop_job <- function(job = state$job) {
    need(!is.null(job), "No job in this R session.")
    if (job$process$is_alive()) job$process$kill()
    message("Queue stopped. Completed checkpoints retained; unsaved optimizer work is lost.")
    invisible(job)
  }
  run <- function(root = ROOT, phenology_file = NULL, prepare_only = FALSE,
                  max_extra_newton = 3L, monitor = TRUE) {
    check_packages()
    need(dir.exists(root), paste("Project directory not found:", root))
    need(is.logical(prepare_only) && length(prepare_only) == 1L && !is.na(prepare_only), "Invalid prepare_only.")
    need(is.logical(monitor) && length(monitor) == 1L && !is.na(monitor), "Invalid monitor.")
    need(is.numeric(max_extra_newton) && length(max_extra_newton) == 1L &&
      is.finite(max_extra_newton) && max_extra_newton >= 0 && max_extra_newton <= 3 &&
      max_extra_newton == floor(max_extra_newton), "max_extra_newton must be 0, 1, 2 or 3.")
    if (!is.null(state$job) && state$job$process$is_alive()) {
      message("This queue is already running."); if (monitor) watch(state$job); return(invisible(state$job))
    }
    root <- normalizePath(root, winslash = "/", mustWork = TRUE)
    source_root <- file.path(root, "output", "phenology_plasticity", "spatial_models")
    need(dir.exists(source_root), "Original spatial_models directory is missing.")
    if (is.null(phenology_file)) phenology_file <- file.path(root, "output", "pheno_estimates_allspp.csv")
    need(file.exists(phenology_file), paste("Phenology file missing:", phenology_file))
    phenology_file <- normalizePath(phenology_file, winslash = "/", mustWork = TRUE)
    jobs <- stats::setNames(lapply(ANALYSES, function(a) locate(source_root, a)), ANALYSES)
    inputs <- unique(c(phenology_file, unlist(lapply(jobs, function(j)
      c(j$checkpoint, file.path(j$source, c("metadata.rds", "prepared_input.rds")))))))
    hashes <- md5(inputs); need(!anyNA(hashes), "Cannot fingerprint input files.")
    signature <- digest::digest(list(version = VERSION, sources = inputs, md5 = hashes,
      versions = versions(), max_extra_newton = max_extra_newton), algo = "sha256")
    output <- file.path(root, "output", "phenology_plasticity", "offset_with_onset_15km")
    dir.create(output, recursive = TRUE, showWarnings = FALSE)
    lock <- file.path(output, "queue.lock")
    if (dir.exists(lock)) {
      who <- file.path(lock, "owner.rds")
      owner <- if (file.exists(who)) tryCatch(readRDS(who), error = function(e) NULL) else NULL
      alive <- if (is.null(owner)) TRUE else tryCatch(
        ps::ps_is_running(ps::ps_handle(owner$pid, time = owner$create_time)), error = function(e) FALSE)
      need(!alive, paste("An onset-adjusted queue may already be running. Inspect:", lock))
      unlink(lock, recursive = TRUE)
    }
    need(dir.create(lock), "Cannot acquire the new queue lock.")
    launched <- FALSE
    on.exit(if (!launched) unlink(lock, recursive = TRUE), add = TRUE)
    run_dir <- file.path(output, paste0("run_", substr(signature, 1, 12)))
    dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
    cfg <- list(signature = signature, root = root, run = run_dir, lock = lock,
      jobs = jobs, phenology_file = phenology_file, prepare_only = prepare_only,
      max_extra_newton = as.integer(max_extra_newton), fit_seed = 20260921L,
      libpath = .libPaths(), input_paths = inputs, input_md5 = hashes)
    cfg$engine <- engine_file(run_dir)
    atomic_rds(cfg, file.path(run_dir, "run_config.rds"))
    csv(data.frame(file = inputs, md5 = hashes), file.path(run_dir, "source_manifest.csv"))
    writeLines(c("TWO OFFSET MODELS WITH ANNUAL ONSET CONTROL", "",
      "Only ONSET_mean_z is added to the original primary Gaussian spatial offset formulas.",
      "The raw onset and its centring/scaling are saved. No onset interaction is added.",
      "Source: original fitted rows; multivoltine uses the recovered original checkpoint.",
      "90-day temperature window, all moderator interactions, random structures and mesh retained.",
      "The selected window is held fixed to isolate the added covariate, NOT reselected by AIC.",
      "Every row must match a finite annual onset, with identical offset to the old fit.",
      "The initial coefficient of ONSET_mean_z is zero; all model parameters are re-estimated.",
      "New full fits are saved before extra Newton steps and before exports.",
      "A numerical pass is not complete validation; new residual diagnostics are still required.",
      "No LRTs for moderator interactions are run; exported p-values are normal Wald tests.",
      "No 30-km, error-distribution, sampling-support or geographic/phylogenetic sensitivity is run here.",
      "Existing sensitivity results do not automatically validate this changed specification.",
      "No onset/first-peak, plasticity extraction or population abundance model is refitted here.",
      "DO NOT use the old extractor on these fits: it deliberately rejects the extra onset term.",
      "Update extraction and refit BOTH abundance models after reviewing these new offset fits.",
      "Both onset and offset abundance-interaction estimates may change in that downstream refit.",
      "Onset-adjusted offset plasticity is a conditional response, not a total seasonal shift.",
      "Measurement uncertainty in the onset/offset estimates is not propagated by this change.",
      "Core package versions must match the original 1.1.0 environment. No packages are installed.",
      "Results are new files; original models and inputs remain untouched."), file.path(run_dir, "README.txt"))
    log <- file.path(run_dir, paste0("controller_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
    proc <- callr::r_bg(function(engine, cfg) {
      e <- new.env(parent = globalenv()); sys.source(engine, envir = e)
      # Refuse changes after the job was configured.
      e$need(identical(e$md5(cfg$input_paths), cfg$input_md5), "Inputs changed before execution.")
      e$controller(cfg)
    }, args = list(cfg$engine, cfg), libpath = cfg$libpath, stdout = log,
      stderr = "2>&1", supervise = TRUE, user_profile = FALSE, system_profile = FALSE)
    handle <- ps::ps_handle(proc$get_pid())
    atomic_rds(list(pid = proc$get_pid(), create_time = ps::ps_create_time(handle)), file.path(lock, "owner.rds"))
    job <- list(process = proc, run = run_dir, started = Sys.time(), log = log)
    state$job <- job; launched <- TRUE
    message("Two primary OFFSET refits with onset control. Output: ", run_dir)
    message("No changes to onset/first-peak fits or abundance fits. Keep this R session open.")
    if (monitor) watch(job)
    invisible(job)
  }
  list(run = run, status = status, watch = watch, stop = stop_job, self_test = self_test)
})

# Define functions only for custom settings:
# options(phenoimpact.offset_onset_functions_only = TRUE)
# source(file.choose(), encoding = "UTF-8")
# pheno_offset_onset$run(prepare_only = TRUE)  # validates both inputs; no fitting
# pheno_offset_onset$run()                    # runs/reuses the two full offset fits
if (!isTRUE(getOption("phenoimpact.offset_onset_functions_only", FALSE))) {
  offset_onset_job <- pheno_offset_onset$run()
}
