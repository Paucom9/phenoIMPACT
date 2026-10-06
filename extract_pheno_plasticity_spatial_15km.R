# Contextual population plasticity from EXISTING spatial fits (15 km).
# Run from the phenoIMPACT R project: source(file.choose()).
# No refitting, prediction, native-state reconstruction or model-file changes.
# Sources: onset, offset_univoltine, RECOVERED offset_multivoltine.
# The original extraction script refitted different models: correlated species
# intercept/slopes and ONSET_mean in the offsets. Those terms are NOT added here.
# Set options(pheno.extract.functions_only = TRUE) to load without running.
# API: pheno_spatial_extract$run(); pheno_spatial_extract$watch()

pheno_spatial_extract <- local({
write_table <- function(x, path) utils::write.csv(x, path, row.names = FALSE, na = "")
formula_text <- function(x) paste(deparse(x, width.cutoff = 500L), collapse = " ")
assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)

atomic_save <- function(x, path) {
  tmp <- tempfile("writing_", tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(x, tmp)
  # These are disposable progress/output files, never model checkpoints.
  if (file.exists(path) && !file.remove(path)) stop("Cannot replace: ", path)
  if (!file.rename(tmp, path)) stop("Cannot finalize: ", path)
  invisible(path)
}
set_phase <- function(path, phase) {
  atomic_save(list(phase = phase, updated = Sys.time()), path)
  message("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", phase)
}
specification <- function(analysis) {
  assert(analysis %in% c("onset", "offset_univoltine", "offset_multivoltine"), "Unknown analysis.")
  window <- if (analysis == "onset") 60L else 90L
  list(response = if (analysis == "onset") "ONSET_mean" else "OFFSET_mean",
    anomaly = paste0("clim_anomaly_tw", window),
    moderators = paste0(c("photo_tw", "clim_background_tw", "clim_predictability_tw", "clim_trend_tw"), window),
    name = switch(analysis,
      onset = "onset_advancement_plasticity_contextual",
      offset_univoltine = "offset_univoltine_termination_plasticity_contextual",
      offset_multivoltine = "offset_multivoltine_delay_plasticity_contextual"),
    sign = if (analysis == "onset") -1 else 1)
}
locate_source <- function(root, analysis) {
  dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[startsWith(basename(dirs), paste0("primary__", analysis, "_"))]
  required <- c("metadata.rds", "prepared_input.rds", "fit_spatial_field.rds")
  dirs <- dirs[vapply(dirs, function(d) all(file.exists(file.path(d, required))), logical(1))]
  assert(length(dirs) == 1L, paste("Expected one primary model directory for", analysis,
    "but found", length(dirs), "-- no automatic selection."))
  normalizePath(dirs, winslash = "/", mustWork = TRUE)
}
select_checkpoint <- function(root, source, analysis) {
  path <- if (analysis == "offset_multivoltine") {
    file.path(root, "offset_recovery_15km", basename(source), "fit_spatial_field.rds")
  } else file.path(source, "fit_spatial_field.rds")
  assert(file.exists(path), paste("Required checkpoint missing:", path,
    "The original failed multivoltine fit is never substituted."))
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

verify_inputs <- function(saved, p, meta, spec) {
  fit <- saved$fit
  assert(inherits(fit, "sdmTMB") && identical(saved$label, "spatial_field") &&
    identical(saved$run_id, meta$run_id) && length(saved$fingerprint) == 1L &&
    nzchar(saved$fingerprint), "Checkpoint identity does not match metadata.")
  assert(isTRUE(meta$settings$mesh_cutoff_km == 15), "Expected the 15-km mesh.")
  assert(identical(p$response, spec$response) && identical(p$anomaly, spec$anomaly),
    "Response or temperature window differs from the expected model.")
  assert(identical(formula_text(p$formula), formula_text(meta$formula)) &&
    all(names(p$data) %in% names(fit$data)) &&
    isTRUE(all.equal(fit$data[, names(p$data), drop = FALSE], p$data, check.attributes = TRUE)),
    "Checkpoint and prepared input disagree on formula, data or row order.")
  versions <- meta$package_versions
  assert(any(versions$package == "sdmTMB" & versions$version == "1.1.0"),
    "This static extractor is validated for saved sdmTMB 1.1.0 objects only.")
  assert(identical(fit$family$family, "gaussian") && identical(fit$family$link, "identity"),
    "Expected a Gaussian identity-link model.")
  assert(length(fit$split_formula) == 1L, "Expected a single model component.")
  expected <- stats::as.formula(paste(spec$response, "~", spec$anomaly, "* (",
    paste(spec$moderators, collapse = " + "), ")"))
  actual <- fit$split_formula[[1]]$form_no_bars
  assert(setequal(attr(stats::terms(actual), "term.labels"),
    attr(stats::terms(expected), "term.labels")), "Unexpected fixed effects; extraction stopped.")
  cnms <- fit$split_formula[[1]]$re_cov_terms$cnms
  assert(length(cnms) == 4L && !anyDuplicated(names(cnms)) &&
    setequal(names(cnms), c("SITE_ID", "SPECIES", "SPECIES_slope", "site_year_id")) &&
    identical(cnms$SPECIES_slope, spec$anomaly) &&
    all(vapply(cnms[setdiff(names(cnms), "SPECIES_slope")],
      function(x) identical(x, "(Intercept)"), logical(1))), "Unexpected random-effect structure.")
  assert(identical(as.character(p$data$SPECIES_slope), as.character(p$data$SPECIES)),
    "Species-slope grouping does not match species identifiers.")
  assert(isTRUE(fit$model$convergence == 0L) && isTRUE(fit$sd_report$pdHess) &&
    length(fit$gradients) > 0L && all(is.finite(fit$gradients)) &&
    max(abs(fit$gradients)) < 0.001, "Selected checkpoint fails numerical convergence checks.")
  assert(identical(names(fit$model$par), names(fit$sd_report$par.fixed)) &&
    isTRUE(all.equal(as.numeric(fit$model$par), as.numeric(fit$sd_report$par.fixed), tolerance = 1e-8)),
    "Saved estimates and uncertainty report refer to different parameter values.")
  invisible(TRUE)
}

# Read plain saved arrays; DO NOT call tidy(), predict(), sanity(), or tmb_obj.
# tidy.sdmTMB() calls reinitialize(), which can rebuild the native model after
# readRDS(). Mapping below follows sdmTMB v1.1.0 R/tidy.R and R/parsing.R.
# https://github.com/sdmTMB/sdmTMB/blob/v1.1.0/R/tidy.R
extract_coefficients <- function(fit, spec) {
  fixed_names <- colnames(stats::model.matrix(fit$split_formula[[1]]$form_no_bars,
    fit$data[seq_len(min(2L, nrow(fit$data))), , drop = FALSE]))
  ix <- which(names(fit$sd_report$par.fixed) == "b_j")
  assert(length(ix) == length(fixed_names) && !anyDuplicated(fixed_names), "Fixed-coefficient mapping failed.")
  beta <- stats::setNames(as.numeric(fit$sd_report$par.fixed[ix]), fixed_names)
  covariance <- fit$sd_report$cov.fixed[ix, ix, drop = FALSE]
  assert(all(is.finite(beta)) && all(is.finite(covariance)) && all(diag(covariance) > 0),
    "Invalid fixed estimates or covariance.")
  dimnames(covariance) <- list(fixed_names, fixed_names)

  # re_b_df records zero-based parameter blocks and the actual group level IDs.
  # Use those IDs, never alphabetical assumptions or factor-position joins.
  sf <- fit$split_formula[[1]]
  blocks <- sf$re_cov_terms$re_b_df
  assert(all(c("start", "end", "group_indices", "level") %in% names(blocks)),
    "Saved random-effect mapping is missing.")
  all_ix <- unlist(Map(function(a, b) seq.int(a, b), blocks$start, blocks$end), use.names = FALSE)
  values <- fit$sd_report$value[grepl("^re_b_pars($|\\[)", names(fit$sd_report$value))]
  errors <- fit$sd_report$sd[grepl("^re_b_pars($|\\[)", names(fit$sd_report$value))]
  active <- !is.na(fit$tmb_map$re_b_pars)
  assert(length(values) == length(active) && length(errors) == length(active) &&
    identical(as.integer(all_ix + 1L), which(as.vector(active))),
    "Random-effect parameter blocks do not match the saved report.")
  # Unlike start/end, group_indices in re_b_df are ONE-based.
  group <- match("SPECIES_slope", names(sf$re_cov_terms$cnms))
  b <- blocks[blocks$group_indices == group, , drop = FALSE]
  assert(nrow(b) > 0L && all(b$start == b$end) && !anyDuplicated(as.character(b$level)),
    "Expected exactly one random slope per species.")
  slopes <- data.frame(SPECIES = as.character(b$level), sp_slope_dev = as.numeric(values[b$start + 1L]),
    sp_slope_dev_marginal_SE = as.numeric(errors[b$start + 1L]))
  assert(setequal(slopes$SPECIES, unique(as.character(fit$data$SPECIES))) &&
    all(is.finite(slopes$sp_slope_dev)) && all(is.finite(slopes$sp_slope_dev_marginal_SE)) &&
    all(slopes$sp_slope_dev_marginal_SE >= 0), "Species slope estimates or IDs are invalid.")
  list(beta = beta, covariance = covariance, slopes = slopes)
}

read_voltinism <- function(path) {
  x <- utils::read.csv(path, sep = ";", stringsAsFactors = FALSE)
  assert("SPECIES" %in% names(x) && !anyDuplicated(x$SPECIES), "Voltinism file has missing or duplicated species keys.")
  pieces <- lapply(setdiff(names(x), "SPECIES"), function(bms) {
    name <- if (bms == "ES.CTBMS") "ES-CTBMS" else if (bms == "ES.ZEBMS") "ES-ZEBMS" else bms
    data.frame(SPECIES = as.character(x$SPECIES), bms_id = name,
      voltinism_lookup = trimws(as.character(x[[bms]])), stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, pieces)
  out$voltinism_lookup[out$voltinism_lookup == ""] <- NA_character_
  assert(!anyDuplicated(out[c("SPECIES", "bms_id")]), "Duplicated species/network voltinism keys.")
  out
}

population_slopes <- function(d, coef, spec, analysis, lookup) {
  keep <- c("SPECIES", "SITE_ID", "bms_id", "YEAR", spec$anomaly, spec$moderators)
  assert(all(keep %in% names(d)), "Required fitted-data columns are missing.")
  d <- d[, keep, drop = FALSE]
  for (key in c("SPECIES", "SITE_ID", "bms_id")) d[[key]] <- as.character(d[[key]])
  assert(!anyNA(d) && all(is.finite(as.matrix(d[, c(spec$anomaly, spec$moderators), drop = FALSE]))),
    "Missing or non-finite fitted rows; no additional filtering is permitted.")
  assert(!anyDuplicated(d[c("SPECIES", "SITE_ID", "YEAR")]), "Repeated population/year rows require review before averaging.")
  networks <- unique(d[c("SPECIES", "SITE_ID", "bms_id")])
  assert(!anyDuplicated(networks[c("SPECIES", "SITE_ID")]), "A population occurs in more than one network.")
  context <- d |>
    dplyr::group_by(SPECIES, SITE_ID, bms_id) |>
    dplyr::summarise(n_years = dplyr::n_distinct(YEAR), n_observations = dplyr::n(),
      first_year = min(YEAR), last_year = max(YEAR),
      anomaly_SD_within_population = stats::sd(.data[[spec$anomaly]]),
      dplyr::across(dplyr::all_of(spec$moderators), mean), .groups = "drop") |>
    as.data.frame()
  context$pop_id <- paste(context$SPECIES, context$SITE_ID, sep = "_")
  assert(!anyDuplicated(context$pop_id), "The legacy population ID contains a collision.")
  context <- dplyr::left_join(context, lookup, by = c("SPECIES", "bms_id"))
  context$voltinism <- context$voltinism_lookup
  if (analysis != "onset") {
    expected <- sub("offset_", "", analysis, fixed = TRUE)
    assert(all(is.na(context$voltinism_lookup) | context$voltinism_lookup == expected),
      "Current voltinism lookup conflicts with the fitted offset subset.")
    context$voltinism <- expected # Known from the saved offset analysis.
  }
  context$voltinism_lookup_missing <- is.na(context$voltinism_lookup)
  index <- match(context$SPECIES, coef$slopes$SPECIES)
  assert(!anyNA(index), "Species random slopes did not match population IDs.")
  context$sp_slope_dev <- coef$slopes$sp_slope_dev[index]
  context$sp_slope_dev_marginal_SE <- coef$slopes$sp_slope_dev_marginal_SE[index]
  context$pop_slope_dev <- 0 # These models contain NO population random slopes.
  assert(spec$anomaly %in% names(coef$beta), "Fixed anomaly slope not found.")
  context$fixed_anomaly_slope <- unname(coef$beta[spec$anomaly])
  context$interaction_contribution <- 0
  for (moderator in spec$moderators) {
    term <- intersect(c(paste(spec$anomaly, moderator, sep = ":"),
      paste(moderator, spec$anomaly, sep = ":")), names(coef$beta))
    assert(length(term) == 1L, paste("Missing/ambiguous interaction:", moderator))
    contribution <- as.numeric(coef$beta[term]) * context[[moderator]]
    context[[paste0("contribution__", moderator)]] <- contribution
    context$interaction_contribution <- context$interaction_contribution + contribution
  }
  context$plasticity_raw <- context$fixed_anomaly_slope + context$sp_slope_dev + context$interaction_contribution
  context[[paste0(spec$name, "_raw")]] <- context$plasticity_raw
  context[[spec$name]] <- spec$sign * context$plasticity_raw
  assert(all(is.finite(context$plasticity_raw)), "Non-finite contextual slopes.")
  context$analysis <- analysis
  context[order(context$SPECIES, context$SITE_ID), , drop = FALSE]
}

extract_worker <- function(source, checkpoint, out, phase_file, analysis, lookup) {
  started <- Sys.time()
  assert(requireNamespace("dplyr", quietly = TRUE), "Package dplyr is required.")
  set_phase(phase_file, "Reading saved fit and exact fitted rows (no fitting)")
  p <- readRDS(file.path(source, "prepared_input.rds"))
  meta <- readRDS(file.path(source, "metadata.rds"))
  saved <- readRDS(checkpoint)
  spec <- specification(analysis)
  set_phase(phase_file, "Checking checkpoint, formula and saved numerical diagnostics")
  verify_inputs(saved, p, meta, spec)
  set_phase(phase_file, "Reading saved fixed coefficients and species slopes")
  coef <- extract_coefficients(saved$fit, spec)
  set_phase(phase_file, "Averaging fitted population contexts and calculating plasticity")
  populations <- population_slopes(p$data, coef, spec, analysis, lookup)
  write_table(populations, file.path(out, "population_plasticity_components.csv"))
  write_table(coef$slopes, file.path(out, "species_random_slopes.csv"))
  write_table(data.frame(term = names(coef$beta), estimate = unname(coef$beta),
    std.error = sqrt(diag(coef$covariance))), file.path(out, "fixed_effects.csv"))
  write_table(data.frame(term = rownames(coef$covariance), coef$covariance, check.names = FALSE),
    file.path(out, "fixed_effect_covariance.csv"))
  write_table(unique(p$data[c("source_row", "SPECIES", "SITE_ID", "YEAR")]),
    file.path(out, "fitted_population_years.csv"))
  write_table(saved$checks, file.path(out, "saved_fit_checks.csv"))
  checks_ok <- isTRUE(saved$checks$all_sanity_checks_pass)
  provenance <- data.frame(analysis = analysis, source_directory = source, checkpoint = checkpoint,
    run_id = meta$run_id, checkpoint_fingerprint = saved$fingerprint, cutoff_km = 15,
    formula = formula_text(p$formula), no_refit = TRUE, native_state_rebuilt = FALSE,
    contextual_slope_uncertainty_propagated = FALSE,
    offset_adjusted_for_onset = FALSE, population_random_slope = FALSE)
  write_table(provenance, file.path(out, "provenance.csv"))
  writeLines(capture.output(utils::sessionInfo()), file.path(out, "extraction_sessionInfo.txt"))
  write_table(meta$package_versions, file.path(out, "fitting_package_versions.csv"))
  result <- list(populations = populations, provenance = provenance,
    summary = data.frame(analysis = analysis, status = "EXTRACTED", n_observations = nrow(p$data),
      n_populations = nrow(populations), n_species = nrow(coef$slopes),
      max_abs_gradient = max(abs(saved$fit$gradients)),
      saved_all_sanity_checks_pass = checks_ok, n_missing_voltinism_lookup = sum(populations$voltinism_lookup_missing),
      elapsed_minutes = as.numeric(difftime(Sys.time(), started, units = "mins")),
      note = if (checks_ok) "" else "Retained diagnostic warning: review saved_fit_checks.csv; univoltine range warning is known."))
  atomic_save(result, file.path(out, "extraction.rds"))
  set_phase(phase_file, "Extraction saved")
  result$summary
}

combine_results <- function(results) {
  keys <- c("pop_id", "SPECIES", "SITE_ID", "bms_id")
  identities <- unique(do.call(rbind, lapply(results, function(x) x$populations[keys])))
  assert(!anyDuplicated(identities$pop_id), "Population identifiers conflict between analyses.")
  all <- identities
  metadata <- dplyr::bind_rows(lapply(results, function(x) x$populations[c(keys, "voltinism")]))
  known <- unique(metadata[!is.na(metadata$voltinism), c("pop_id", "voltinism")])
  assert(!anyDuplicated(known$pop_id), "Voltinism differs across saved model subsets.")
  all <- dplyr::left_join(all, known, by = "pop_id")
  for (analysis in names(results)) {
    spec <- specification(analysis)
    p <- results[[analysis]]$populations
    x <- p[, c("pop_id", "n_years", paste0(spec$name, "_raw"), spec$name), drop = FALSE]
    names(x)[names(x) == "n_years"] <- switch(analysis,
      onset = "n_years_onset", offset_univoltine = "n_years_offset_uni", offset_multivoltine = "n_years_offset_multi")
    all <- dplyr::left_join(all, x, by = "pop_id")
  }
  all$offset_termination_plasticity_contextual <- ifelse(all$voltinism == "univoltine",
    all$offset_univoltine_termination_plasticity_contextual,
    ifelse(all$voltinism == "multivoltine", all$offset_multivoltine_delay_plasticity_contextual, NA_real_))
  all$offset_termination_plasticity_contextual_raw <- all$offset_termination_plasticity_contextual
  all <- all[order(all$SPECIES, all$SITE_ID), , drop = FALSE]
  # Legacy table was anchored on onset. Also save the full union, so no
  # offset-only populations disappear without an explicit audit.
  legacy <- all[all$pop_id %in% results$onset$populations$pop_id, , drop = FALSE]
  legacy$bms_id <- NULL
  list(legacy = legacy, all = all,
    offset_only = all[!all$pop_id %in% legacy$pop_id, , drop = FALSE])
}

monitor_job <- function(job, interval = 30) {
  assert(!is.null(job$process), "No extraction job is available.")
  tryCatch({
    repeat {
      phase <- "Starting"
      if (file.exists(job$phase_file)) {
        state <- tryCatch(suppressWarnings(readRDS(job$phase_file)), error = function(e) NULL)
        if (!is.null(state)) phase <- state$phase
      }
      alive <- job$process$is_alive()
      message("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", job$analysis,
        " (", job$index, "/3) | ", phase, " | elapsed ",
        round(as.numeric(difftime(Sys.time(), job$started, units = "mins")), 1),
        " min | worker running: ", alive)
      if (!alive) break
      Sys.sleep(interval)
    }
    job$process$get_result()
  }, interrupt = function(e) {
    message("Display/queue paused; extraction worker continues. Keep R open. Resume display with pheno_spatial_extract$watch().")
    NULL
  })
}

wait_active <- function(job) {
  if (is.null(job$process)) return(invisible(NULL))
  while (isTRUE(tryCatch(job$process$is_alive(), error = function(e) FALSE))) {
    message("[", format(Sys.time(), "%H:%M:%S"), "] Waiting for the existing worker to finish before loading another large fit.")
    Sys.sleep(30)
  }
}

run_extract <- function(root = here::here("output", "phenology_plasticity", "spatial_models"),
                        voltinism_file = here::here("data", "voltinism", "species_country_voltinism.csv")) {
  assert(requireNamespace("callr", quietly = TRUE), "Package callr is required.")
  assert(requireNamespace("dplyr", quietly = TRUE), "Package dplyr is required.")
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  lookup <- read_voltinism(voltinism_file)
  analyses <- c("onset", "offset_univoltine", "offset_multivoltine")
  sources <- stats::setNames(vapply(analyses, function(a) locate_source(root, a), character(1)), analyses)
  checkpoints <- stats::setNames(vapply(analyses, function(a) select_checkpoint(root, sources[[a]], a), character(1)), analyses)
  for (option in c("pheno.spatial.active_job", "pheno.offset.recovery.active_job",
                   "pheno.review.active_job", "pheno.extract.active_job")) wait_active(getOption(option))
  if (exists("spatial_job", envir = .GlobalEnv, inherits = FALSE)) wait_active(get("spatial_job", envir = .GlobalEnv))
  out_root <- tempfile(paste0("population_plasticity_15km_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_"), tmpdir = root)
  dir.create(out_root)
  out_root <- normalizePath(out_root, winslash = "/", mustWork = TRUE)
  engine <- file.path(out_root, "extraction_engine.R")
  dump(c("write_table", "formula_text", "assert", "atomic_save", "set_phase", "specification",
    "verify_inputs", "extract_coefficients", "population_slopes", "extract_worker"),
    file = engine, envir = environment(run_extract), control = "all")
  write_table(data.frame(path = normalizePath(voltinism_file, winslash = "/"),
    md5 = unname(tools::md5sum(voltinism_file))), file.path(out_root, "voltinism_lookup_source.csv"))
  writeLines(c(
    "Contextual population plasticity from the existing 15-km spatial models; no refit.",
    "Raw slope = fixed anomaly coefficient + species slope deviation + sum(interaction coefficient * population mean moderator).",
    "Population contexts average the exact fitted years, without restandardization or additional filtering.",
    "All four interactions are retained, regardless of their significance.",
    "Positive onset_advancement = earlier onset under increasing thermal anomaly (minus raw slope).",
    "Both offset columns retain the raw sign: positive = later termination, negative = earlier termination.",
    "Units are days per one unit of the FITTED anomaly variable. Do not assume days/degree C without checking its original scaling.",
    "Species deviations are partially pooled. No independent population random slope was fitted.",
    "Additive spatial fields and random intercepts are not added to a temperature derivative.",
    "Unlike the old extraction script, offsets are not conditioned on ONSET_mean and species intercept/slope effects are independent.",
    "Comparisons with that old extraction are therefore not an isolated test of spatial correction.",
    "Estimates only: full contextual-slope uncertainty and covariance are NOT propagated.",
    "sp_slope_dev_marginal_SE describes only the species deviation, NOT the complete population slope.",
    "The main CSV preserves the original onset-anchored population set and legacy column names.",
    "The all_populations CSV retains offset-only populations; excluded_from_onset_anchored_table.csv records them.",
    "Missing voltinism lookup values are flagged; fitted offset subset membership supplies their known classification.",
    "Recovered multivoltine checkpoint is mandatory. The univoltine range warning remains; extraction does not validate model adequacy.",
    "First peak is not extracted because the original population-plasticity script used onset and the two offsets only.",
    "Progress reports stage and elapsed minutes; no unreliable completion-time estimate is supplied."
  ), file.path(out_root, "README.txt"))
  status <- list()
  for (i in seq_along(analyses)) {
    analysis <- analyses[i]
    out <- file.path(out_root, analysis)
    dir.create(out)
    phase_file <- file.path(out, "progress.rds")
    log <- file.path(out, "worker.log")
    process <- callr::r_bg(function(engine, args) {
      e <- new.env(parent = globalenv())
      sys.source(engine, envir = e)
      do.call(e$extract_worker, args)
    }, args = list(engine = engine, args = list(source = sources[[analysis]],
      checkpoint = checkpoints[[analysis]], out = out, phase_file = phase_file,
      analysis = analysis, lookup = lookup)), stdout = log, stderr = "2>&1", supervise = TRUE, wd = getwd())
    job <- list(process = process, phase_file = phase_file, analysis = analysis,
      index = i, started = Sys.time(), log_file = log, output_root = out_root)
    options(pheno.extract.active_job = job)
    message("Extracting ", analysis, ". Log: ", log)
    result <- tryCatch(monitor_job(job), error = function(e) data.frame(analysis = analysis,
      status = "ERROR", note = conditionMessage(e)))
    if (process$is_alive()) stop("Queue paused; the current extraction continues. Output: ", out_root)
    if (is.null(result)) stop("Extraction did not return a confirmed result. Inspect: ", log)
    status[[analysis]] <- result
    write_table(dplyr::bind_rows(status), file.path(out_root, "extraction_summary.csv"))
  }
  assert(all(vapply(status, function(x) identical(x$status, "EXTRACTED"), logical(1))),
    paste("Some extractions failed. No combined table was published. Inspect:", out_root))
  results <- stats::setNames(lapply(analyses, function(a) readRDS(file.path(out_root, a, "extraction.rds"))), analyses)
  tables <- combine_results(results)
  stem <- "phenological_population_plasticity_contextual_spatial_15km"
  csv <- file.path(out_root, paste0(stem, ".csv"))
  write_table(tables$legacy, csv)
  write_table(tables$all, file.path(out_root, paste0(stem, "_all_populations.csv")))
  write_table(tables$offset_only, file.path(out_root, "excluded_from_onset_anchored_table.csv"))
  atomic_save(tables, file.path(out_root, paste0(stem, ".rds")))
  if (requireNamespace("writexl", quietly = TRUE)) {
    tryCatch(writexl::write_xlsx(list(plasticity = tables$legacy,
      all_populations = tables$all, extraction_summary = dplyr::bind_rows(status)),
      file.path(out_root, paste0(stem, ".xlsx"))), error = function(e) warning("CSV/RDS saved; XLSX failed: ", conditionMessage(e)))
  }
  message("Extraction complete. Main table: ", csv,
    "\nReview extraction_summary.csv and README.txt before fitting population-trend models.")
  invisible(list(output_root = out_root, csv = csv, status = dplyr::bind_rows(status)))
}
list(run = run_extract, watch = function() monitor_job(getOption("pheno.extract.active_job")))
})

if (!isTRUE(getOption("pheno.extract.functions_only", FALSE))) {
  spatial_plasticity_extraction <- pheno_spatial_extract$run()
}
