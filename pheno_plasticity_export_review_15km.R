# pheno_plasticity_export_review_15km.R
# Read-only export of the four saved 15-km fits (2026-09-27).
# Multivoltine MUST use the recovered checkpoint. Original models are never edited.
# No optimizer is called. Native TMB state is rebuilt at the saved parameters.
# Run from the phenoIMPACT R project: source(file.choose()).
# Keep this R session open. A new timestamped export folder and ZIP are created.
# To load only: options(pheno.review.functions_only = TRUE).
# API: pheno_spatial_review$run(root = <spatial_models directory>)
#      pheno_spatial_review$watch()  # resume the current worker display
# Escape pauses display/queue; the current worker continues, as in the main workflow.

pheno_spatial_review <- local({
write_table <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE, na = "")
}

formula_text <- function(f) paste(deparse(f, width.cutoff = 500L), collapse = " ")

atomic_save <- function(object, path) {
  tmp <- tempfile(pattern = "checkpoint_", tmpdir = dirname(path), fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(object, tmp, compress = "gzip")
  # Windows cannot rename over an existing file. Keep a backup until success.
  backup <- paste0(path, ".previous")
  existed <- file.exists(path)
  if (existed && !file.copy(path, backup, overwrite = TRUE)) stop("Cannot back up ", path)
  if (existed && !file.remove(path)) stop("Cannot replace ", path)
  if (!file.rename(tmp, path)) {
    if (existed) file.copy(backup, path, overwrite = TRUE)
    stop("Cannot finalize checkpoint: ", path)
  }
  if (existed) unlink(backup)
  invisible(path)
}

fit_checks <- function(fit, cfg) {
  gradients <- fit$gradients
  grad <- if (length(gradients) && all(is.finite(gradients))) {
    max(abs(gradients))
  } else NA_real_
  code <- fit$model$convergence
  if (length(code) != 1L) code <- NA_integer_
  pd <- isTRUE(fit$sd_report$pdHess)
  fixed <- tryCatch(sdmTMB::tidy(fit, effects = "fixed", conf.int = FALSE),
                    error = function(e) NULL)
  se_ok <- !is.null(fixed) && all(c("estimate", "std.error") %in% names(fixed)) &&
    nrow(fixed) > 0L && all(is.finite(fixed$estimate)) &&
    all(is.finite(fixed$std.error) & fixed$std.error > 0)
  sanity <- tryCatch(sdmTMB::sanity(fit, gradient_thresh = cfg$gradient_threshold,
                                  silent = TRUE), error = function(e) NULL)
  all_sanity <- !is.null(sanity) && isTRUE(all(unlist(sanity)))
  data.frame(convergence_code = code, positive_definite_Hessian = pd,
             max_abs_gradient = grad, finite_fixed_SE = se_ok,
             all_sanity_checks_pass = all_sanity,
             fit_ok = isTRUE(code == 0L) && pd && is.finite(grad) &&
               grad < cfg$gradient_threshold && se_ok)
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

set_phase <- function(path, phase) {
  atomic_save(list(phase = phase, updated = Sys.time()), path)
  message("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", phase)
}

locate_source <- function(root, analysis) {
  dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[startsWith(basename(dirs), paste0("primary__", analysis, "_"))]
  needed <- c("metadata.rds", "prepared_input.rds", "fit_spatial_field.rds")
  dirs <- dirs[vapply(dirs, function(d) all(file.exists(file.path(d, needed))), logical(1))]
  if (length(dirs) != 1L) stop("Expected one saved primary model for ", analysis,
    "; found ", length(dirs), ". No model selected automatically.")
  normalizePath(dirs, winslash = "/", mustWork = TRUE)
}

select_checkpoint <- function(root, source, analysis) {
  if (analysis == "offset_multivoltine") {
    path <- file.path(root, "offset_recovery_15km", basename(source), "fit_spatial_field.rds")
    if (!file.exists(path)) stop("Recovered multivoltine checkpoint missing: ", path,
      ". The original failed fit will not be substituted.")
  } else path <- file.path(source, "fit_spatial_field.rds")
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

write_gzip_csv <- function(x, path) {
  con <- gzfile(path, open = "wt", compression = 6)
  on.exit(close(con))
  utils::write.csv(x, con, row.names = FALSE, na = "")
}

align_predictions <- function(pred, d) {
  if (!all(c("source_row", "est", "epsilon_st") %in% names(pred)) ||
      nrow(pred) != nrow(d) || anyDuplicated(pred$source_row) || anyDuplicated(d$source_row)) {
    stop("Prediction columns or row identifiers are invalid.")
  }
  index <- match(d$source_row, pred$source_row)
  if (anyNA(index)) stop("Predictions do not align with prepared data.")
  out <- pred[index, , drop = FALSE]
  if (any(!is.finite(out$est)) || any(!is.finite(out$epsilon_st))) {
    stop("Non-finite predictions or annual field values.")
  }
  out
}

export_worker <- function(source, checkpoint, out, phase_file, analysis, seed) {
  pkgs <- c("sdmTMB", "TMB", "Matrix", "fmesher", "lme4", "dplyr")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))
  started <- Sys.time()
  set_phase(phase_file, "Reading saved fit and input")
  p <- readRDS(file.path(source, "prepared_input.rds"))
  meta <- readRDS(file.path(source, "metadata.rds"))
  saved <- readRDS(checkpoint)
  verify_inputs(saved, p, meta)
  if (!identical(saved$fit$family$family, "gaussian") ||
      !identical(saved$fit$family$link, "identity")) stop("Expected Gaussian identity model.")
  if (!isTRUE(saved$fit$model$convergence == 0L) ||
      !isTRUE(saved$fit$sd_report$pdHess) || !length(saved$fit$gradients) ||
      any(!is.finite(saved$fit$gradients)) || max(abs(saved$fit$gradients)) >= 0.001) {
    stop("The selected checkpoint fails numerical convergence checks.")
  }
  original_parameters <- saved$fit$model$par
  set_phase(phase_file, "Rebuilding native state at saved estimates (no refit)")
  fit <- rebuild_fit(saved$fit, p)
  rm(saved)
  invisible(gc())
  if (!identical(original_parameters, fit$model$par)) stop("Saved parameters changed.")
  checks <- fit_checks(fit, list(gradient_threshold = 0.001))
  if (!isTRUE(checks$fit_ok)) stop("Selected checkpoint fails fixed-effect checks.")
  write_table(checks, file.path(out, "fit_checks.csv"))
  details <- sdmTMB::sanity(fit, gradient_thresh = 0.001, silent = TRUE)
  write_table(data.frame(check = names(details), passed = as.logical(unlist(details))),
    file.path(out, "sanity_checks.csv"))
  parameters <- sdmTMB::tidy(fit, effects = "ran_pars", conf.int = TRUE)
  write_table(parameters, file.path(out, "spatial_parameters.csv"))
  write_table(sdmTMB::tidy(fit, effects = "fixed", conf.int = TRUE),
    file.path(out, "fixed_effects.csv"))
  writeLines(capture.output(utils::sessionInfo()), file.path(out, "sessionInfo.txt"))
  write_table(data.frame(analysis = analysis, source_directory = source,
    checkpoint = checkpoint, run_id = meta$run_id, cutoff_km = meta$settings$mesh_cutoff_km,
    residual_seed = seed, no_refit = TRUE), file.path(out, "provenance.csv"))

  set_phase(phase_file, "Extracting fitted values and annual spatial fields")
  d <- p$data
  pred <- align_predictions(stats::predict(fit), d)
  keep <- unique(c("source_row", "SITE_ID", "SPECIES", "bms_id", "YEAR", "x_km", "y_km",
    grep("^(clim_|photo_)", names(d), value = TRUE)))
  if (!all(keep %in% names(d))) stop("Required row or coordinate columns missing.")
  rows <- d[, keep, drop = FALSE]
  rows$observed <- d[[p$response]]
  rows$fitted <- pred$est
  rows$original_response_residual <- p$original_residuals
  rows$spatial_response_residual <- rows$observed - rows$fitted
  rows$annual_field <- pred$epsilon_st
  if (length(p$original_residuals) != nrow(d) ||
      any(!is.finite(as.matrix(rows[, c("observed", "fitted", "original_response_residual",
        "spatial_response_residual", "annual_field")])))) stop("Invalid residual alignment or values.")
  rm(pred)

  fields <- rows |>
    dplyr::group_by(bms_id, YEAR, SITE_ID, x_km, y_km) |>
    dplyr::summarise(field_spread = diff(range(annual_field)),
      annual_field = mean(annual_field),
      n_observations = dplyr::n(), .groups = "drop")
  if (any(fields$field_spread > 1e-6)) stop("Annual field differs within site/year.")
  fields$field_spread <- NULL
  # Collapse co-located site IDs before calculating equally weighted annual means.
  locations <- fields |>
    dplyr::group_by(YEAR, x_km, y_km) |>
    dplyr::summarise(annual_field = mean(annual_field), .groups = "drop")
  years <- locations |>
    dplyr::group_by(YEAR) |>
    dplyr::summarise(n_locations = dplyr::n(),
      mean_field = mean(annual_field), sd_field_across_locations = stats::sd(annual_field),
      minimum_field = min(annual_field), maximum_field = max(annual_field),
      q05 = as.numeric(stats::quantile(annual_field, 0.05)),
      q95 = as.numeric(stats::quantile(annual_field, 0.95)), .groups = "drop")
  fields <- dplyr::left_join(fields, years[, c("YEAR", "mean_field")], by = "YEAR")
  fields$field_centered_within_year <- fields$annual_field - fields$mean_field
  write_table(fields, file.path(out, "annual_fields.csv"))
  write_table(years, file.path(out, "annual_field_summary.csv"))
  write_table(data.frame(response = p$response, n_observations = nrow(d),
    x_span_km = diff(range(d$x_km)), y_span_km = diff(range(d$y_km))),
    file.path(out, "data_summary.csv"))
  # Save useful response residuals even if the MVN calculation later fails.
  write_gzip_csv(rows, file.path(out, "residuals_response.csv.gz"))

  set_phase(phase_file, "Calculating MVN residuals (one reproducible draw)")
  mvn_error <- NULL
  mvn <- tryCatch({
    set.seed(seed)
    r <- as.numeric(stats::residuals(fit, type = "mle-mvn"))
    index <- match(d$source_row, fit$data$source_row)
    if (length(r) != nrow(fit$data) || anyNA(index) || any(!is.finite(r))) {
      stop("Invalid or misaligned MVN residuals.")
    }
    r[index]
  }, error = function(e) { mvn_error <<- conditionMessage(e); NULL })
  if (!identical(original_parameters, fit$model$par)) stop("Parameter estimates changed.")
  if (!is.null(mvn)) {
    set_phase(phase_file, "Writing residual tables")
    write_gzip_csv(data.frame(source_row = d$source_row, mle_mvn_residual = mvn),
      file.path(out, "residuals_mle_mvn.csv.gz"))
  } else writeLines(mvn_error, file.path(out, "mvn_error.txt"))
  result <- data.frame(analysis = analysis,
    status = if (is.null(mvn_error)) "EXPORT_OK" else "PARTIAL_MVN_FAILED",
    n_observations = nrow(d), max_abs_gradient = checks$max_abs_gradient,
    all_sanity_checks_pass = checks$all_sanity_checks_pass,
    elapsed_minutes = as.numeric(difftime(Sys.time(), started, units = "mins")),
    notes = if (is.null(mvn_error)) "" else mvn_error)
  write_table(result, file.path(out, "export_status.csv"))
  set_phase(phase_file, result$status)
  result
}

monitor_job <- function(job, interval = 60) {
  if (is.null(job)) stop("No export worker is registered.")
  tryCatch({
    repeat {
      state <- if (file.exists(job$phase_file)) {
        suppressWarnings(tryCatch(readRDS(job$phase_file), error = function(e) NULL))
      } else NULL
      phase <- if (is.null(state)) "Starting" else state$phase
      cat(sprintf("[%s] %s (%d/4) | %s | elapsed %.0f min | worker running: %s\n",
        format(Sys.time(), "%Y-%m-%d %H:%M:%S"), job$analysis, job$index, phase,
        as.numeric(difftime(Sys.time(), job$started, units = "mins")), job$process$is_alive()))
      if (!job$process$is_alive()) break
      deadline <- Sys.time() + interval
      while (job$process$is_alive() && Sys.time() < deadline) Sys.sleep(1)
    }
    job$process$get_result()
  }, interrupt = function(e) {
    message("Display paused. The current export continues; later models have not started.")
    invisible(NULL)
  })
}

wait_active <- function(job) {
  if (is.null(job) || is.null(job$process)) return(invisible(TRUE))
  while (job$process$is_alive()) {
    message("Waiting for the active worker before starting this export.")
    deadline <- Sys.time() + 60
    while (job$process$is_alive() && Sys.time() < deadline) Sys.sleep(1)
  }
  invisible(TRUE)
}

make_zip <- function(out) {
  destination <- paste0(out, ".zip")
  files <- list.files(out, recursive = TRUE, full.names = FALSE)
  # Native TMB objects are deliberately never included in this export.
  files <- files[!grepl("\\.rds$", files, ignore.case = TRUE)]
  if (requireNamespace("zip", quietly = TRUE)) {
    zip::zipr(destination, files = files, root = out, mode = "mirror")
  } else if (.Platform$OS.type == "windows") {
    ps_quote <- function(x) paste0("'", gsub("'", "''", x, fixed = TRUE), "'")
    # ZipFile avoids wildcard interpretation of project paths.
    command <- paste(c("$ErrorActionPreference = 'Stop'",
      "Add-Type -AssemblyName System.IO.Compression.FileSystem",
      paste0("[System.IO.Compression.ZipFile]::CreateFromDirectory(",
        ps_quote(out), ", ", ps_quote(destination), ")")), collapse = "; ")
    status <- system2("powershell.exe", c("-NoProfile", "-NonInteractive", "-Command", shQuote(command)))
    if (status != 0L) stop("ZIP creation failed; compress the export folder manually.")
  } else {
    previous <- setwd(out)
    on.exit(setwd(previous), add = TRUE)
    utils::zip(destination, files = files)
  }
  if (!file.exists(destination) || file.info(destination)$size == 0) stop("ZIP was not created.")
  destination
}

run_review <- function(root = here::here("output", "phenology_plasticity", "spatial_models")) {
  if (!requireNamespace("callr", quietly = TRUE)) stop("Package callr is required.")
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  analyses <- c("onset", "first_peak", "offset_univoltine", "offset_multivoltine")
  sources <- setNames(vapply(analyses, function(a) locate_source(root, a), character(1)), analyses)
  checkpoints <- setNames(vapply(analyses, function(a) select_checkpoint(root, sources[[a]], a), character(1)), analyses)
  for (option in c("pheno.spatial.active_job", "pheno.offset.recovery.active_job", "pheno.review.active_job")) {
    wait_active(getOption(option))
  }
  out_root <- tempfile(paste0("review_15km_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_"), tmpdir = root)
  dir.create(out_root)
  out_root <- normalizePath(out_root, winslash = "/", mustWork = TRUE)
  engine <- file.path(out_root, "export_engine.R")
  engine_names <- c("write_table", "formula_text", "atomic_save", "fit_checks", "verify_inputs",
    "rebuild_fit", "set_phase", "write_gzip_csv", "align_predictions", "export_worker")
  dump(engine_names, file = engine, envir = environment(run_review), control = "all")
  writeLines(c(
    "Read-only diagnostic export of the four primary 15-km fits.",
    "Multivoltine uses the recovered checkpoint; see each provenance.csv.",
    "No optimizer is called; native state is rebuilt at saved parameter values.",
    "EXPORT_OK indicates successful extraction, not model adequacy.",
    "Annual fields are conditional estimates at sampled site/year locations, in days.",
    "field_centered_within_year subtracts the equally weighted mean over distinct coordinates that year.",
    "Annual means and dispersions also reflect changing spatial sampling coverage.",
    "Response residuals are observed minus fitted; original and spatial fits use matched rows.",
    "MVN residuals use one seeded approximate posterior draw, as recommended by sdmTMB.",
    "MVN residuals join response residuals by source_row within each analysis.",
    "All random effects are included in conditional fitted values.",
    "Coordinates are EPSG:3035 in kilometres. Sampled locations only; no interpolation/extrapolation.",
    "These files support maps, QQ plots, and residual spatial checks; no adequacy decision is automated.",
    "Raw model checkpoints are not included. A failed/partial export must be reviewed in export_summary.csv.",
    "Sources: https://sdmtmb.github.io/sdmTMB/reference/predict.sdmTMB.html",
    "https://sdmtmb.github.io/sdmTMB/reference/residuals.sdmTMB.html"
  ), file.path(out_root, "README.txt"))
  results <- list()
  jobs <- list()
  for (i in seq_along(analyses)) {
    analysis <- analyses[i]
    out <- file.path(out_root, analysis)
    dir.create(out)
    phase_file <- file.path(out, "progress.rds")
    log <- file.path(out, "worker.log")
    process <- callr::r_bg(function(engine, args) {
      e <- new.env(parent = globalenv())
      sys.source(engine, envir = e)
      do.call(e$export_worker, args)
    }, args = list(engine = engine, args = list(source = sources[[analysis]],
      checkpoint = checkpoints[[analysis]], out = out, phase_file = phase_file,
      analysis = analysis, seed = 20260922L)), stdout = log, stderr = "2>&1",
      supervise = TRUE, wd = getwd())
    job <- list(process = process, phase_file = phase_file, started = Sys.time(),
      analysis = analysis, index = i, log_file = log, output_root = out_root)
    options(pheno.review.active_job = job)
    jobs[[analysis]] <- job
    message("Export ", i, "/4: ", analysis, ". Log: ", log)
    result <- tryCatch(monitor_job(job), error = function(e) data.frame(analysis = analysis,
      status = "WORKER_ERROR", n_observations = NA_integer_, max_abs_gradient = NA_real_,
      all_sanity_checks_pass = NA, elapsed_minutes = NA_real_, notes = conditionMessage(e)))
    if (process$is_alive()) stop("Export queue paused. Keep R open; the current worker continues.")
    if (is.null(result)) stop("Export stopped without a confirmed result. See ", log)
    results[[analysis]] <- result
    write_table(do.call(rbind, results), file.path(out_root, "export_summary.csv"))
    print(result[, c("analysis", "status", "elapsed_minutes")], row.names = FALSE)
    # Only transient monitor state; no model checkpoint is ever removed.
    if (file.exists(phase_file)) unlink(phase_file)
  }
  archive <- tryCatch(make_zip(out_root), error = function(e) {
    message(conditionMessage(e)); NA_character_
  })
  message("Export finished. Review export_summary.csv for failed or partial tasks.")
  if (!is.na(archive)) message("SEND THIS ZIP: ", archive) else message("Compress and send this folder: ", out_root)
  invisible(list(status = do.call(rbind, results), jobs = jobs, output_root = out_root, zip_file = archive))
}

list(run = run_review, watch = function() monitor_job(getOption("pheno.review.active_job")))
})

if (!isTRUE(getOption("pheno.review.functions_only", FALSE))) {
  spatial_review_job <- pheno_spatial_review$run()
}
