# ============================================================================
# phenoIMPACT | Updated offset plasticity -> spatial abundance-trend models
# Version 1.0 | 2026-10-04
#
# DEFINES FUNCTIONS ONLY. In RStudio:
# source(file.choose(), encoding = "UTF-8")
# updated_abundance_job <- pheno_updated_abundance$run(monitor = FALSE)
# pheno_updated_abundance$status()
#
# Sources are explicit; no "latest" model is selected:
#   offsets:   offset_with_onset_15km/run_3911d7bc87ac
#   abundance: spatial_abundance_15km/run_650c41f15949
#
# ONSET plasticity is UNCHANGED: preserve the previous preparation_plasticity.csv
# (all extracted populations used for scaling, BEFORE the abundance join).
# Extract NEW contextual OFFSET slopes from the two completed onset-adjusted
# fits. These are partial anomaly derivatives HOLDING ANNUAL ONSET FIXED.
# Additive ONSET_mean_z does not enter that derivative; its inclusion changes
# the estimated anomaly, interaction and species-slope coefficients.
#
# Keep the EXACT abundance observations, their ordering, time scale, group
# membership, species/site random terms, population AR1 links/gaps, and saved
# 15-km mesh of run_650c41f15949. Replace only the offset-plasticity columns of X.
# Standardize the NEW offsets over the SAME population universe within group.
# Re-estimate ALL abundance coefficients, random effects and variance parameters.
# Warm starts are from the previous SPATIAL abundance fits, NOT fixed values.
#
# Gamma/log; independent species/site intercepts and year slopes; stationary
# population AR1 (irregular integer gaps); IID annual Matern fields shared by
# species, with range/variance shared among years. Original C++ template reused.
# Only TWO abundance fits. No new non-spatial baseline or reduced-model LRTs.
# No phenology models are refitted. No onset/first-peak model is loaded.
#
# Sequential processes. Fits are saved as portable arrays, NOT native TMB
# environments. Successful evaluations/optimized parameters are checkpointed
# before uncertainty estimation. TMB::sdreport(getJointPrecision = FALSE).
#
# Residual diagnostics of the NEW offsets and abundance models remain PENDING.
# Passing numerical checks is not full model validation. Plasticity-estimation
# uncertainty is NOT propagated (same plug-in approach as the previous run).
# No prediction of overnight completion time is made.
#
# All new outputs go into a NEW directory. Original inputs are read-only.
# No packages or fonts are installed/updated. Keep the parent R session open.
# ============================================================================

pheno_updated_abundance <- local({
  VERSION <- "updated_plasticity_spatial_abundance_1.0"
  ROOT <- "E:/phenoIMPACT project/code/phenoIMPACT"
  GROUPS <- c("univoltine", "multivoltine")
  ANALYSES <- paste0("offset_", GROUPS)
  CORE <- c("TMB", "Matrix", "fmesher")
  state <- new.env(parent = emptyenv())
  assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
  formula_text <- function(f) paste(deparse(f, width.cutoff = 500L), collapse = " ")
  write_table <- function(x, path) utils::write.csv(x, path, row.names = FALSE, na = "")
  stamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  same <- function(a, b, tolerance = 1e-8) isTRUE(all.equal(a, b, tolerance = tolerance))
  same_num <- function(a, b, tolerance = 1e-8) {
    length(a) == length(b) && all(is.finite(a)) && all(is.finite(b)) &&
      (length(a) == 0L || max(abs(as.numeric(a) - as.numeric(b))) <= tolerance)
  }
  atomic_save <- function(x, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    tmp <- tempfile(".writing_", tmpdir = dirname(path), fileext = ".rds")
    on.exit(unlink(tmp), add = TRUE)
    saveRDS(x, tmp, compress = "gzip")
    bak <- paste0(path, ".previous")
    # Rename, do not duplicate multi-GB files. This only touches this NEW run.
    if (file.exists(bak)) assert(unlink(bak) == 0L, paste("Cannot replace backup:", bak))
    had <- file.exists(path)
    if (had) assert(file.rename(path, bak), paste("Cannot retain previous checkpoint:", path))
    if (!file.rename(tmp, path)) {
      if (had) file.rename(bak, path)
      stop("Cannot finalize checkpoint: ", path, call. = FALSE)
    }
    invisible(path)
  }
  phase <- function(cfg, stage, task = "", detail = "") {
    atomic_save(list(stage = stage, task = task, detail = detail, time = Sys.time()),
                file.path(cfg$out, "progress.rds"))
    message("[", stamp(), "] ", task, " | ", stage, " | ", detail)
  }
  disk_check <- function(path, minimum_GiB = 50, verbose = TRUE) {
    x <- ps::ps_disk_usage(path)
    assert(nrow(x) == 1L && is.finite(x$available), "Could not determine free disk space.")
    free <- as.numeric(x$available) / 1024^3
    if (verbose) message("Disk preflight: ", round(free, 1), " GiB available.")
    assert(free >= minimum_GiB, paste("Only", round(free, 1), "GiB free; need at least",
      minimum_GiB, "GiB. No further model is started."))
    invisible(free)
  }
  check_packages <- function() {
    p <- c("callr", "ps", "digest", "dplyr", "lme4", "sdmTMB", CORE)
    missing <- p[!vapply(p, requireNamespace, logical(1), quietly = TRUE)]
    assert(!length(missing), paste("Missing packages:", paste(missing, collapse = ", "),
      "-- use the original R library. Nothing was installed or refitted."))
    assert(as.character(utils::packageVersion("sdmTMB")) == "1.1.0",
           "Use the original sdmTMB 1.1.0 library, not the error-sensitivity library.")
    invisible(TRUE)
  }
  verify_versions <- function(file, packages) {
    old <- utils::read.csv(file, stringsAsFactors = FALSE)
    assert(all(c("package", "version") %in% names(old)) && all(packages %in% old$package),
           paste("Incomplete package-version file:", file))
    now <- vapply(packages, function(p) as.character(utils::packageVersion(p)), character(1))
    wanted <- as.character(old$version[match(packages, old$package)])
    assert(identical(unname(now), wanted), paste("Package versions differ. Required:",
      paste(paste(packages, wanted, sep = "="), collapse = ", ")))
  }
  check_cache <- function(path, cfg) {
    if (!file.exists(path)) return(NULL)
    x <- readRDS(path)
    assert(identical(x$signature, cfg$signature), paste("Cache/input mismatch:", path))
    x
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

  extract_worker <- function(cfg, analysis) {
    check_packages()
    out <- file.path(cfg$out, "extraction", analysis)
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    dest <- file.path(out, "extraction.rds")
    if (!is.null(check_cache(dest, cfg))) return("EXTRACTED_CACHED")
    src <- file.path(cfg$offset_run, analysis)
    verify_versions(file.path(src, "package_versions.csv"), c("sdmTMB", CORE, "lme4"))
    phase(cfg, "EXTRACTING", analysis, "Reading the completed onset-adjusted fit; no refitting")
    inp <- readRDS(file.path(src, "input.rds"))
    done <- readRDS(file.path(src, "completed.rds"))
    saved <- readRDS(file.path(src, "fit_with_onset.rds"))
    f <- saved$fit
    assert(identical(saved$analysis, analysis) && identical(inp$analysis, analysis) &&
      identical(done$signature, saved$signature) && identical(inp$signature, saved$signature) &&
      identical(done$status, "FIT_OK"), "Offset checkpoint/completion identities do not match.")
    assert(inherits(f, "sdmTMB") && identical(f$family$family, "gaussian") &&
      identical(f$family$link, "identity") && !isTRUE(f$reml), "Unexpected offset family/link/REML.")
    assert(isTRUE(f$model$convergence == 0L) && isTRUE(f$sd_report$pdHess) &&
      length(f$gradients) > 0L && all(is.finite(f$gradients)) && max(abs(f$gradients)) < .001,
      "Offset fit does not pass the saved numerical checks.")
    assert(identical(names(f$model$par), names(f$sd_report$par.fixed)) &&
      same_num(f$model$par, f$sd_report$par.fixed), "Stale covariance or fitted parameters.")
    assert(identical(formula_text(inp$formula), formula_text(saved$formula)) &&
      same(f$data[, names(inp$data), drop = FALSE], inp$data), "Offset data/formula mismatch.")
    assert(isTRUE(inp$settings$mesh_cutoff_km == 15) && isTRUE(inp$settings$best_window == 90) &&
      identical(inp$settings$spatial, "off") && identical(inp$settings$spatiotemporal, "iid"),
      "Expected 15-km annual IID spatial offset model at the 90-day window.")
    spec <- specification(analysis)
    sf <- f$split_formula[[1]]
    expected <- stats::as.formula(paste(spec$response, "~ ONSET_mean_z +", spec$anomaly,
      "* (", paste(spec$moderators, collapse = "+"), ")"))
    assert(length(f$split_formula) == 1L &&
      setequal(attr(stats::terms(sf$form_no_bars), "term.labels"),
               attr(stats::terms(expected), "term.labels")),
      "Offset fixed formula must include additive ONSET_mean_z and the four original interactions.")
    cn <- sf$re_cov_terms$cnms
    assert(length(cn) == 4L && setequal(names(cn), c("SITE_ID", "SPECIES", "SPECIES_slope", "site_year_id")) &&
      identical(cn$SPECIES_slope, spec$anomaly) &&
      all(vapply(cn[setdiff(names(cn), "SPECIES_slope")], function(x) identical(x, "(Intercept)"), logical(1))),
      "Unexpected offset random-effects mapping.")
    assert(identical(as.character(inp$data$SPECIES_slope), as.character(inp$data$SPECIES)),
      "Species slope grouping is inconsistent.")
    assert(all(is.finite(inp$data$ONSET_mean)) && all(is.finite(inp$data$ONSET_mean_z)) &&
      same_num(inp$data$ONSET_mean_z,
        (inp$data$ONSET_mean - inp$onset_scaling$centre) / inp$onset_scaling$SD),
      "Annual onset covariate/scaling is inconsistent.")
    phase(cfg, "EXTRACTING", analysis, "Reading saved coefficients and species slopes, without native TMB calls")
    coef <- extract_coefficients(f, spec)
    # Additional check after run_extra_optimization: compare the independent
    # portable parameter file with the ADREPORT values used for extraction.
    pars <- readRDS(file.path(src, "parameters.rds"))
    assert(same_num(pars$b_j, coef$beta, 1e-6), "Portable fixed parameters differ from the saved report.")
    values <- f$sd_report$value[grepl("^re_b_pars($|\\[)", names(f$sd_report$value))]
    assert(same_num(as.vector(pars$re_b_pars), values, 1e-6),
      "Random-effect report is stale relative to parameters.rds; extraction stopped.")
    group <- sub("offset_", "", analysis, fixed = TRUE)
    lookup <- unique(data.frame(SPECIES = as.character(inp$data$SPECIES),
      bms_id = as.character(inp$data$bms_id), stringsAsFactors = FALSE))
    # Group membership comes from the verified fitted subset, not a new lookup.
    lookup$voltinism_lookup <- group
    p <- population_slopes(inp$data, coef, spec, analysis, lookup)
    p$offset_adjusted_for_onset <- TRUE
    p$offset_conditional_raw <- p$plasticity_raw
    write_table(p, file.path(out, "population_plasticity_components.csv"))
    write_table(coef$slopes, file.path(out, "species_random_slopes.csv"))
    write_table(data.frame(term = names(coef$beta), estimate = unname(coef$beta),
      std.error = sqrt(diag(coef$covariance))), file.path(out, "fixed_effects.csv"))
    write_table(data.frame(term = rownames(coef$covariance), coef$covariance, check.names = FALSE),
      file.path(out, "fixed_effect_covariance.csv"))
    prov <- data.frame(analysis = analysis, source_run = cfg$offset_run,
      source_checkpoint = file.path(src, "fit_with_onset.rds"), source_signature = saved$signature,
      formula = formula_text(inp$formula), n_observations = nrow(inp$data),
      n_populations = nrow(p), n_species = nrow(coef$slopes), max_abs_gradient = max(abs(f$gradients)),
      offset_adjusted_for_onset = TRUE, new_phenology_fit = FALSE,
      plasticity_uncertainty_propagated = FALSE, residual_review = "PENDING")
    write_table(prov, file.path(out, "provenance.csv"))
    atomic_save(list(signature = cfg$signature, populations = p, provenance = prov), dest)
    "EXTRACTED"
  }

  update_plasticity_table <- function(old, offset_tables) {
    required <- c("pop_id", "SPECIES", "SITE_ID", "voltinism",
      "onset_advancement_plasticity_contextual", "offset_termination_plasticity_contextual",
      "onset_plasticity_z", "offset_plasticity_bio_z")
    assert(all(required %in% names(old)) && nrow(old) > 0L,
      "The previous preparation_plasticity.csv lacks required columns.")
    for (k in c("pop_id", "SPECIES", "SITE_ID", "voltinism")) old[[k]] <- as.character(old[[k]])
    assert(!anyNA(old[required]) && !anyDuplicated(old$pop_id) &&
      !anyDuplicated(old[c("SPECIES", "SITE_ID")]) && all(old$voltinism %in% GROUPS) &&
      identical(old$pop_id, paste(old$SPECIES, old$SITE_ID, sep = "_")),
      "Invalid previous population identities/values. No populations were filtered.")
    p <- old
    p$previous_offset_termination_plasticity_contextual <- old$offset_termination_plasticity_contextual
    p$previous_offset_plasticity_bio_z <- old$offset_plasticity_bio_z
    scaling <- audits <- list()
    for (g in GROUPS) {
      i <- which(p$voltinism == g); a <- offset_tables[[g]]
      assert(length(i) > 1L && !anyDuplicated(a$pop_id), "Invalid offset population table.")
      ix <- match(p$pop_id[i], a$pop_id)
      assert(!anyNA(ix), paste("New offset estimates missing for existing", g,
        "scaling populations. No changed sample will be silently fitted."))
      assert(identical(as.character(a$SPECIES[ix]), p$SPECIES[i]) &&
        identical(as.character(a$SITE_ID[ix]), p$SITE_ID[i]) &&
        all(as.character(a$voltinism[ix]) == g), "Offset population-key/group mismatch.")
      raw <- a$offset_conditional_raw[ix]
      assert(all(is.finite(raw)), "Non-finite new offset plasticity.")
      on <- p$onset_advancement_plasticity_contextual[i]
      on_sd <- stats::sd(on); on_mean <- mean(on)
      assert(is.finite(on_sd) && on_sd > 0 &&
        same_num((on - on_mean)/on_sd, p$onset_plasticity_z[i], 1e-6),
        "Onset standardization is not reproduced over the previous population universe.")
      old_bio <- if (g == "univoltine") -old$offset_termination_plasticity_contextual[i] else old$offset_termination_plasticity_contextual[i]
      assert(same_num(as.numeric(scale(old_bio)), old$offset_plasticity_bio_z[i], 1e-6),
        "Previous offset standardization not reproduced; inspect the preparation CSV.")
      bio <- if (g == "univoltine") -raw else raw
      sd_bio <- stats::sd(bio); mu_bio <- mean(bio)
      assert(is.finite(sd_bio) && sd_bio > 0, "New offset plasticity is constant/invalid.")
      p$offset_termination_plasticity_contextual[i] <- raw
      p$offset_plasticity_bio[i] <- bio
      p$offset_plasticity_bio_z[i] <- (bio - mu_bio)/sd_bio
      scaling[[g]] <- data.frame(group = g, n_populations = length(i), onset_mean = on_mean,
        onset_SD = on_sd, previous_offset_bio_mean = mean(old_bio), previous_offset_bio_SD = stats::sd(old_bio),
        new_offset_bio_mean = mu_bio, new_offset_bio_SD = sd_bio,
        scaling_population = "identical previous extracted populations; before abundance matching")
      audits[[g]] <- data.frame(group = g, n_populations = length(i),
        old_new_offset_raw_cor = stats::cor(old$offset_termination_plasticity_contextual[i], raw),
        old_new_offset_z_cor = stats::cor(old$offset_plasticity_bio_z[i], p$offset_plasticity_bio_z[i]),
        median_abs_offset_z_change = stats::median(abs(p$offset_plasticity_bio_z[i] - old$offset_plasticity_bio_z[i])),
        max_abs_offset_z_change = max(abs(p$offset_plasticity_bio_z[i] - old$offset_plasticity_bio_z[i])))
    }
    p$offset_adjusted_for_onset <- TRUE
    # Main/legacy names intentionally retain their established meaning.
    list(populations = p, scaling = do.call(rbind, scaling), audit = do.call(rbind, audits))
  }

  prepare_abundance_worker <- function(cfg) {
    check_packages()
    done_path <- file.path(cfg$out, "preparation_complete.rds")
    if (!is.null(check_cache(done_path, cfg))) {
      assert(all(file.exists(file.path(cfg$out, paste0("input_", GROUPS, ".rds")))),
        "Missing prepared abundance inputs.")
      return("PREPARED_CACHED")
    }
    phase(cfg, "PREPARING_ABUNDANCE", "", "Preserving the previous population universe and ONSET estimates")
    oldp <- utils::read.csv(cfg$previous_plasticity_file, stringsAsFactors = FALSE, check.names = FALSE)
    offsets <- lapply(ANALYSES, function(a) {
      x <- check_cache(file.path(cfg$out, "extraction", a, "extraction.rds"), cfg)
      assert(!is.null(x), "Both new offset extractions must be complete before abundance preparation.")
      x$populations
    }); names(offsets) <- GROUPS
    updated <- update_plasticity_table(oldp, offsets); p <- updated$populations
    write_table(p, file.path(cfg$out, "phenological_population_plasticity_onset_adjusted_15km.csv"))
    write_table(updated$scaling, file.path(cfg$out, "plasticity_scaling.csv"))
    write_table(updated$audit, file.path(cfg$out, "plasticity_comparison.csv"))
    # The old abundance mesh is reused bit-for-bit, never rebuilt at a new cutoff.
    assert(file.copy(file.path(cfg$abundance_reference, "mesh.rds"),
      file.path(cfg$out, "mesh.rds"), overwrite = TRUE), "Cannot copy the reference mesh.")
    assert(identical(unname(tools::md5sum(file.path(cfg$abundance_reference, "mesh.rds"))),
                     unname(tools::md5sum(file.path(cfg$out, "mesh.rds")))), "Mesh-copy mismatch.")
    rows <- list()
    for (g in GROUPS) {
      phase(cfg, "PREPARING_ABUNDANCE", g, "Replacing offset predictors; identical response, rows, times and AR1 links")
      old <- readRDS(file.path(cfg$abundance_reference, paste0("input_", g, ".rds")))
      fit <- readRDS(file.path(cfg$abundance_reference, paste0(g, "__full"), "fit.rds"))
      assert(identical(old$group, g) && identical(fit$group, g) && identical(fit$variant, "full") &&
        identical(fit$data_signature, old$data_signature) && isTRUE(fit$checks$numerical_checks_ok),
        "The reference spatial abundance fit/input is inconsistent or needs numerical review.")
      d <- old$frame; ix <- match(d$pop_id, p$pop_id)
      assert(!anyNA(ix) && all(p$voltinism[ix] == g) &&
        identical(as.character(d$SPECIES), p$SPECIES[ix]) && identical(as.character(d$SITE_ID), p$SITE_ID[ix]),
        "Updated plasticities do not match every original abundance row.")
      assert(same_num(d$onset_plasticity_z, p$onset_plasticity_z[ix], 1e-6) &&
        same_num(d$offset_plasticity_bio_z, p$previous_offset_plasticity_bio_z[ix], 1e-6),
        "Reference abundance predictors do not match the previous plasticity preparation file.")
      assert(identical(old$data_signature, digest::digest(list(d, old$data$X), algo = "sha256")),
        "Original abundance frame/design signature mismatch.")
      # Verify portable spatial parameters and all original random-effect indices
      # without replaying the very expensive likelihood or calling an optimizer.
      op <- fit$parameters
      eta0 <- as.numeric(old$data$X %*% op$beta) +
        op$b_species[old$data$species + 1L, 1L] + old$data$time * op$b_species[old$data$species + 1L, 2L] +
        op$b_site[old$data$site + 1L, 1L] + old$data$time * op$b_site[old$data$site + 1L, 2L] + op$ar_obs +
        fit$field_at_sites[cbind(old$data$site + 1L, old$data$year + 1L)]
      err <- max(abs(eta0 - fit$eta))
      assert(is.finite(err) && err < 1e-5, "Previous abundance predictor cannot be reconstructed.")
      new_frame <- d
      # ONSET z is retained EXACTLY as used by the old abundance model.
      new_frame$offset_plasticity_bio_z <- p$offset_plasticity_bio_z[ix]
      X <- stats::model.matrix(~ year_decade * onset_plasticity_z +
                                year_decade * offset_plasticity_bio_z, new_frame)
      assert(setequal(colnames(X), colnames(old$data$X)), "Changed abundance design terms.")
      X <- X[, colnames(old$data$X), drop = FALSE]
      changing <- c("offset_plasticity_bio_z", "year_decade:offset_plasticity_bio_z")
      keep <- setdiff(colnames(X), changing)
      assert(same_num(X[, keep, drop = FALSE], old$data$X[, keep, drop = FALSE], 1e-10) &&
        nrow(X) == nrow(d) && all(is.finite(X)) && qr(X)$rank == ncol(X),
        "Unexpected change or rank-deficient new abundance design.")
      newdata <- old$data; newdata$X <- X
      assert(identical(newdata[setdiff(names(newdata), "X")], old$data[setdiff(names(old$data), "X")]),
        "Response, temporal scale, grouping or AR1 links were changed.")
      inp <- list(signature = cfg$signature, schema = VERSION, group = g, data = newdata,
        frame = new_frame, sites = old$sites, years = old$years, par = op,
        data_signature = digest::digest(list(new_frame, X), algo = "sha256"),
        previous_data_signature = old$data_signature, previous_spatial_beta = fit$beta,
        previous_spatial_V = fit$V, previous_spatial_checks = fit$checks,
        offset_adjusted_for_onset = TRUE, diagnostics_status = "PENDING")
      atomic_save(inp, file.path(cfg$out, paste0("input_", g, ".rds")))
      rows[[g]] <- data.frame(group = g, n_observations = nrow(d), n_populations = length(unique(d$pop_id)),
        n_species = nrow(op$b_species), n_sites = nrow(op$b_site), first_year = min(d$YEAR), last_year = max(d$YEAR),
        same_abundance_rows = TRUE, same_response = TRUE, same_time_AR1 = TRUE, same_onset_predictor = TRUE,
        old_eta_max_abs_error = err, max_offset_z_change = max(abs(new_frame$offset_plasticity_bio_z - d$offset_plasticity_bio_z)),
        data_signature = inp$data_signature)
      rm(old, fit, op, eta0, inp, newdata, new_frame, X, d); invisible(gc())
    }
    write_table(do.call(rbind, rows), file.path(cfg$out, "abundance_preparation_checks.csv"))
    atomic_save(list(signature = cfg$signature, status = "PREPARED", time = Sys.time()), done_path)
    "PREPARED"
  }

export_predictions <-
function(beta,V,years) {
  rows<-list()
  for(v in c('onset_plasticity_z','offset_plasticity_bio_z'))for(z in c(-1,0,1)) {
    it<-paste0('year_decade:',v);t<-(years-min(years))/10
    slope<-beta['year_decade']+beta[it]*z
    variance<-V['year_decade','year_decade']+z^2*V[it,it]+2*z*V['year_decade',it]
    assert(variance>=-1e-10,'Invalid fixed-effect prediction variance.')
    se<-sqrt(max(variance,0))*t;eta<-as.numeric(slope)*t
    rows[[length(rows)+1L]]<-data.frame(variable=v,plasticity_z=z,YEAR=years,
      percent_change=100*expm1(eta),low=100*expm1(eta-1.96*se),high=100*expm1(eta+1.96*se))
  }
  do.call(rbind,rows)
}

fixed_table <-
function(beta,V) {
  se<-sqrt(diag(V));z<-beta/se
  data.frame(term=names(beta),estimate=as.numeric(beta),std.error=as.numeric(se),
    conf.low=as.numeric(beta-1.96*se),conf.high=as.numeric(beta+1.96*se),
    p_Wald=as.numeric(2*pnorm(-abs(z))))
}

make_objective <-
function(input, spatial=TRUE, variant='full') {
  data<-input$data;par<-input$par
  for(nm in c('A','M0','M1','M2')) data[[nm]]<-methods::as(methods::as(methods::as(data[[nm]],'dMatrix'),'generalMatrix'),'TsparseMatrix')
  if(variant!='full') {
    variable<-switch(variant,without_onset_interaction='onset_plasticity_z',
      without_offset_interaction='offset_plasticity_bio_z',stop('Unknown variant.'))
    drop<-match(paste0('year_decade:',variable),colnames(data$X))
    assert(!is.na(drop),'Missing interaction.')
    data$X<-data$X[,-drop,drop=FALSE];par$beta<-par$beta[-drop]
  }
  random<-c('b_species','b_site','ar_obs')
  map<-NULL
  if(spatial) random<-c(random,'field') else map<-list(log_range=factor(NA),log_spatial_sd=factor(NA))
  obj<-TMB::MakeADFun(data=data,parameters=par,random=random,map=map,DLL='trend_spde',silent=TRUE,
    inner.control=list(maxit=1000,trace=FALSE))
  list(obj=obj,data=data,par=par)
}

with_mesh <-
function(input, mesh=NULL, spatial=TRUE) {
  data<-input$data;par<-input$par
  zero<-Matrix::sparseMatrix(i=integer(),j=integer(),dims=c(0L,0L))
  data$use_spatial<-as.integer(spatial)
  if(spatial) {
    assert(!is.null(mesh),'A mesh is required.')
    data$A<-methods::as(fmesher::fm_basis(mesh$mesh,loc=as.matrix(input$sites[c('x_km','y_km')])), 'CsparseMatrix')
    assert(all(abs(Matrix::rowSums(data$A)-1)<1e-8),'Some coordinates lie outside the mesh.')
    data$M0<-mesh$M0;data$M1<-mesh$M1;data$M2<-mesh$M2
    par$field<-matrix(0,nrow(mesh$M0),length(input$years))
  } else {
    data$A<-Matrix::sparseMatrix(i=integer(),j=integer(),dims=c(nrow(input$sites),0L))
    data$M0<-zero;data$M1<-zero;data$M2<-zero
    par$field<-matrix(numeric(),0,0)
  }
  list(data=data,par=par)
}

  compile_worker <- function(cfg) {
    check_packages()
    folder <- file.path(cfg$out, "compiled")
    dir.create(folder, recursive = TRUE, showWarnings = FALSE)
    source <- file.path(cfg$abundance_reference, "compiled", "trend_spde.cpp")
    target <- file.path(folder, "trend_spde.cpp")
    mark <- file.path(folder, "compile_complete.rds")
    ready <- check_cache(mark, cfg)
    if (!is.null(ready) && file.exists(TMB::dynlib(file.path(folder, "trend_spde"))) &&
        identical(unname(tools::md5sum(target)), unname(tools::md5sum(source)))) return("COMPILED_CACHED")
    assert(file.copy(source, target, overwrite = TRUE), "Cannot copy the original C++ template.")
    phase(cfg, "COMPILING", "abundance", "Compiling the unchanged reference model in the NEW run directory")
    oldwd <- getwd(); on.exit(setwd(oldwd), add = TRUE); setwd(folder)
    code <- TMB::compile("trend_spde.cpp", flags = "-O1 -g0")
    assert(isTRUE(code == 0L) && file.exists(TMB::dynlib("trend_spde")),
      "TMB compilation failed. Use the same Rtools/toolchain as the original abundance run.")
    dyn.load(TMB::dynlib("trend_spde"))
    TMB::openmp(n = 1L, DLL = "trend_spde")
    err <- template_self_test()
    write_table(data.frame(test = "small_joint_likelihood_R_vs_CPP", absolute_error = err,
      passed = TRUE), file.path(folder, "template_self_test.csv"))
    atomic_save(list(signature = cfg$signature, cpp_md5 = unname(tools::md5sum(source))), mark)
    "COMPILED"
  }
  template_self_test <- function() {
    # Deterministic synthetic small example. NO model fitting; compare the exact
    # joint Gamma + Gaussian + stationary AR1 + Matern prior likelihood with R.
    n <- 12L; ns <- 2L; nt <- 2L; ny <- 3L; nv <- 2L
    tm <- rep(c(-.1, 0, .1), 4L)
    year <- rep(0:2, 4L); sp <- rep(c(0L, 1L), each = 6L)
    site <- rep(rep(0:1, each = 3L), 2L)
    prev <- rep(c(-1L, 0L, 1L), 4L)
    for (k in 0:3) prev[k*3L + 2:3] <- k*3L + 0:1
    d <- list(y = seq(1, 2.1, length.out = n), X = cbind(1, tm), time = tm,
      species = sp, site = site, year = year, previous = prev, gap = rep(1L, n),
      A = Matrix::Diagonal(2L), M0 = Matrix::Diagonal(2L),
      M1 = Matrix::Diagonal(2L, x = rep(.4,2L)), M2 = Matrix::Diagonal(2L, x = rep(.7,2L)), use_spatial = 1L)
    p <- list(beta = c(.5, -.2), log_sd = log(c(.3,.2,.25,.15,.35)),
      rho_raw = .4, log_shape = log(2.5), b_species = matrix(c(.1,-.1,.02,-.02), ns, 2),
      b_site = matrix(c(.03,-.03,.01,-.01), nt, 2), ar_obs = sin(seq_len(n))*.04,
      log_range = log(3), log_spatial_sd = log(.25), field = matrix(seq_len(6)*.01, nv, ny))
    for (nm in c("A", "M0", "M1", "M2")) d[[nm]] <- methods::as(
      methods::as(methods::as(d[[nm]], "dMatrix"), "generalMatrix"), "TsparseMatrix")
    obj <- TMB::MakeADFun(data = d, parameters = p, DLL = "trend_spde", silent = TRUE)
    on.exit(TMB::FreeADFun(obj), add = TRUE)
    got <- obj$fn(obj$par)
    ss <- exp(p$log_sd); rho <- p$rho_raw/sqrt(1+p$rho_raw^2)
    expected <- 0
    for (j in 1:2) expected <- expected - sum(stats::dnorm(p$b_species[,j],0,ss[j],log=TRUE)) -
      sum(stats::dnorm(p$b_site[,j],0,ss[j+2],log=TRUE))
    for (i in seq_len(n)) {
      mu <- if (prev[i]<0) 0 else rho*p$ar_obs[prev[i]+1L]
      sd <- if (prev[i]<0) ss[5] else ss[5]*sqrt(1-rho^2)
      expected <- expected - stats::dnorm(p$ar_obs[i],mu,sd,log=TRUE)
    }
    kappa <- sqrt(8)/exp(p$log_range)
    tau <- 1/(sqrt(4*pi)*kappa*exp(p$log_spatial_sd))
    Q <- as.matrix(tau^2*(kappa^4*d$M0+2*kappa^2*d$M1+d$M2))
    for (j in seq_len(ny)) expected <- expected + .5*(nv*log(2*pi)-as.numeric(determinant(Q, logarithm=TRUE)$modulus)+
      as.numeric(crossprod(p$field[,j], Q%*%p$field[,j])))
    eta <- as.numeric(d$X%*%p$beta) + p$b_species[sp+1L,1] + tm*p$b_species[sp+1L,2] +
      p$b_site[site+1L,1] + tm*p$b_site[site+1L,2] + p$ar_obs + p$field[cbind(site+1L,year+1L)]
    expected <- expected - sum(stats::dgamma(d$y,shape=exp(p$log_shape),scale=exp(eta)/exp(p$log_shape),log=TRUE))
    err <- abs(got-expected)
    assert(is.finite(err) && err < 1e-8, "Small joint-likelihood self-test failed; no abundance model was fitted.")
    err
  }
  load_dll <- function(cfg) {
    path <- TMB::dynlib(file.path(cfg$out, "compiled", "trend_spde"))
    assert(file.exists(path), "Newly compiled abundance engine is missing.")
    dyn.load(path); TMB::openmp(n = 1L, DLL = "trend_spde")
  }
  export_abundance <- function(fit, input, cfg, group) {
    out <- file.path(cfg$out, paste0(group, "__full"))
    write_table(fit$checks, file.path(out, "fit_checks.csv"))
    writeLines(fit$warnings, file.path(out, "warnings.txt"))
    new <- fixed_table(fit$beta, fit$V)
    old <- fixed_table(input$previous_spatial_beta, input$previous_spatial_V)
    write_table(new, file.path(out, "fixed_effects.csv"))
    write_table(data.frame(term = rownames(fit$V), fit$V, check.names = FALSE),
      file.path(out, "fixed_effect_covariance.csv"))
    write_table(merge(old, new, by = "term", all = TRUE,
      suffixes = c("_previous_spatial", "_updated_spatial")), file.path(out, "coefficient_comparison.csv"))
    if (isTRUE(fit$checks$numerical_checks_ok)) {
      write_table(export_predictions(fit$beta, fit$V, input$years), file.path(out, "predictions_spatial.csv"))
      write_table(export_predictions(input$previous_spatial_beta, input$previous_spatial_V, input$years),
        file.path(out, "predictions_previous_spatial.csv"))
    }
    write_table(data.frame(SITE_ID = rep(input$sites$SITE_ID, times = length(input$years)),
      YEAR = rep(input$years, each = nrow(input$sites)), field = as.vector(fit$field_at_sites)),
      file.path(out, "annual_spatial_field.csv"))
    writeLines(c("MODEL: Gamma/log abundance with the same observations and mesh as the reference run.",
      "Only the offset-plasticity predictors were updated; all fitted parameters were re-estimated.",
      "ONSET predictor and population scaling universe are unchanged.",
      "New OFFSET plasticity is conditional on annual ONSET_mean.",
      "predictions_spatial.csv: fixed-effect relative trajectories; other plasticity at mean.",
      "Species/site, population AR1 and annual spatial effects set to zero for these curves.",
      "Ribbons: pointwise normal-Wald 95% confidence intervals; not prediction intervals.",
      "Previous comparison also uses spatial abundance models, not non-spatial baselines.",
      "The offset SD has been recomputed, using the identical previous population universe.",
      "p_Wald, NOT reduced-model LRTs. No plasticity-estimation uncertainty propagation.",
      "NUMERICAL CHECKS DO NOT REPLACE RESIDUAL DIAGNOSTICS: new residual review PENDING."),
      file.path(out, "README.txt"))
  }
  fit_abundance_worker <- function(cfg, group) {
    check_packages(); disk_check(cfg$out, cfg$minimum_free_GiB)
    out <- file.path(cfg$out, paste0(group, "__full"))
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    input <- readRDS(file.path(cfg$out, paste0("input_", group, ".rds")))
    assert(identical(input$signature, cfg$signature), "New abundance input signature mismatch.")
    cache <- file.path(out, "fit.rds")
    existing <- check_cache(cache, cfg)
    if (!is.null(existing)) {
      assert(identical(existing$data_signature, input$data_signature), "Saved fit/data mismatch.")
      export_abundance(existing, input, cfg, group)
      return(if (isTRUE(existing$checks$numerical_checks_ok)) "FIT_OK_CACHED" else "NUMERICAL_REVIEW_REQUIRED")
    }
    load_dll(cfg)
    phase(cfg, "BUILDING_SPATIAL_MODEL", group, paste(nrow(input$frame), "rows; same 15-km mesh; new offset plasticity"))
    mesh <- readRDS(file.path(cfg$out, "mesh.rds"))
    a <- with_mesh(input, mesh, spatial = TRUE)
    a$par <- input$par # with_mesh initially zeros fields; restore old SPATIAL starts.
    assert(identical(dim(a$par$field), as.integer(c(nrow(mesh$M0), length(input$years)))),
      "Spatial starting-field dimensions changed.")
    optimized_path <- file.path(out, "optimized_parameters.rds")
    checkpoint_path <- file.path(out, "last_evaluation.rds")
    optimized <- check_cache(optimized_path, cfg)
    warm <- if (is.null(optimized)) check_cache(checkpoint_path, cfg) else NULL
    if (!is.null(optimized)) {
      assert(identical(optimized$data_signature, input$data_signature), "Optimized checkpoint/data mismatch.")
      a$par <- optimized$parameters
    } else if (!is.null(warm)) {
      assert(identical(warm$data_signature, input$data_signature), "Evaluation checkpoint/data mismatch.")
      a$par <- warm$parameters
    }
    built <- make_objective(a, spatial = TRUE, variant = "full"); obj <- built$obj
    on.exit(try(TMB::FreeADFun(obj), silent = TRUE), add = TRUE)
    warnings <- character()
    catch_warning <- function(w) {
      warnings <<- unique(c(warnings, conditionMessage(w)))
      message("WARNING: ", conditionMessage(w)); invokeRestart("muffleWarning")
    }
    started <- Sys.time(); count <- 0L; ng <- 0L; last <- Sys.time()-3600
    extract_current <- function(par) {
      # last.par is the successful CURRENT evaluation, not an earlier best point.
      params <- obj$env$parList(par = obj$env$last.par)
      assert(same_num(params$beta, par[names(par) == "beta"], 1e-8),
        "Current random/fixed parameter checkpoint is misaligned.")
      params
    }
    objective <- function(par) {
      value <- obj$fn(par); count <<- count + 1L
      if (is.finite(value) && (count == 1L || as.numeric(difftime(Sys.time(),last,units="secs")) >= 60)) {
        disk_check(cfg$out, 10, verbose = FALSE)
        params <- extract_current(par)
        atomic_save(list(signature = cfg$signature, data_signature = input$data_signature,
          par = par, parameters = params, objective = value, time = Sys.time()), checkpoint_path)
        phase(cfg, "FITTING_ABUNDANCE", group, paste("function",count,"gradient",ng,
          "objective",signif(value,9),"| checkpoint saved"))
        last <<- Sys.time()
      }
      value
    }
    gradient <- function(par) {ng <<- ng+1L; obj$gr(par)}
    if (!is.null(optimized)) {
      phase(cfg, "RECOVERING_UNCERTAINTY", group, "Reusing optimized parameters; no repeated optimizer run")
      opt <- optimized$optimizer
      assert(identical(names(obj$par), names(opt$par)), "Optimizer parameter names changed.")
      v <- obj$fn(opt$par)
      assert(is.finite(v) && abs(v-opt$objective) < .01,
        "Saved optimized likelihood could not be reproduced.")
    } else {
      phase(cfg, "FITTING_ABUNDANCE", group, if (is.null(warm))
        "Warm start from the previous spatial model; all parameters free" else
        "Warm restart from saved evaluation; not an exact optimizer continuation")
      set.seed(20261004L + match(group, GROUPS))
      opt <- withCallingHandlers(stats::nlminb(obj$par, objective, gradient,
        control = list(iter.max = cfg$optimizer_iterations, eval.max = cfg$optimizer_evaluations,
          rel.tol = 1e-10)), warning = catch_warning)
      if (opt$convergence != 0L) {
        phase(cfg, "REFINING_ABUNDANCE", group, "BFGS refinement after non-zero nlminb convergence code")
        second <- withCallingHandlers(stats::optim(opt$par, objective, gradient, method = "BFGS",
          control = list(maxit = 300L, reltol = 1e-10)), warning = catch_warning)
        if (is.finite(second$value) && second$value <= opt$objective + 1e-6) opt <- list(
          par = second$par, objective = second$value, convergence = second$convergence,
          message = "BFGS refinement", first_optimizer = opt, counts = second$counts)
      }
      value <- obj$fn(opt$par)
      assert(is.finite(value) && abs(value-opt$objective) < .01, "Final objective mismatch.")
      parameters <- extract_current(opt$par)
      optimized <- list(signature = cfg$signature, data_signature = input$data_signature,
        parameters = parameters, optimizer = opt)
      atomic_save(optimized, optimized_path)
    }
    phase(cfg, "UNCERTAINTY", group, "Optimized parameters saved; computing fixed-effect Hessian and covariance")
    sdr <- withCallingHandlers(TMB::sdreport(obj, par.fixed = opt$par,
      getJointPrecision = FALSE, getReportCovariance = FALSE), warning = catch_warning)
    g <- as.numeric(obj$gr(opt$par)); Vfull <- sdr$cov.fixed
    ids <- which(names(opt$par) == "beta")
    beta <- stats::setNames(as.numeric(opt$par[ids]), colnames(built$data$X))
    V <- Vfull[ids,ids,drop=FALSE]; dimnames(V) <- list(names(beta),names(beta))
    scaled <- if (all(is.finite(g)) && all(is.finite(Vfull)))
      sqrt(max(0,as.numeric(crossprod(g,Vfull%*%g)))) else Inf
    good <- opt$convergence == 0L && isTRUE(sdr$pdHess) && all(is.finite(V)) &&
      all(diag(V)>0) && is.finite(scaled) && scaled<.05 && is.finite(opt$objective)
    obj$fn(opt$par); parameters <- extract_current(opt$par)
    report <- obj$report(obj$env$last.par)
    extent <- max(diff(range(input$sites$x_km)),diff(range(input$sites$y_km)))
    checks <- data.frame(convergence_code = opt$convergence, pdHess = isTRUE(sdr$pdHess),
      max_abs_gradient = max(abs(g)), gradient_covariance_norm = scaled, numerical_checks_ok = good,
      logLik = -opt$objective, parameter_df = length(opt$par), AIC = 2*opt$objective+2*length(opt$par),
      nobs = nrow(input$frame), spatial_range_km = as.numeric(report$range), spatial_sd = as.numeric(report$spatial_sd),
      range_greater_than_1_5_extent = as.numeric(report$range)>1.5*extent,
      spatial_sd_below_0_001 = as.numeric(report$spatial_sd)<.001,
      AR1_rho = as.numeric(report$rho), Gamma_shape = as.numeric(report$shape))
    fit <- list(signature = cfg$signature, data_signature = input$data_signature,
      parameters = parameters, optimizer = opt, beta = beta, V = V, cov_fixed = Vfull,
      checks = checks, eta = as.numeric(report$eta), field_at_sites = report$spatial,
      group = group, variant = "full", fit_minutes = as.numeric(difftime(Sys.time(),started,units="mins")),
      warnings = unique(warnings), offset_adjusted_for_onset = TRUE, residual_diagnostics = "PENDING")
    # Preserve the fit BEFORE any CSV/prediction export or review-archive step.
    atomic_save(fit, cache)
    export_abundance(fit, input, cfg, group)
    if (isTRUE(good)) {
      atomic_save(list(signature=cfg$signature,status="FIT_OK",time=Sys.time()),file.path(out,"completed.rds"))
      "FIT_OK"
    } else "NUMERICAL_REVIEW_REQUIRED"
  }

  collect <- function(cfg) {
    for (name in c("fit_checks", "fixed_effects", "coefficient_comparison")) {
      rows <- lapply(GROUPS, function(g) {
        path <- file.path(cfg$out,paste0(g,"__full"),paste0(name,".csv"))
        if (!file.exists(path)) return(NULL)
        cbind(group=g,utils::read.csv(path,stringsAsFactors=FALSE,check.names=FALSE))
      })
      rows <- Filter(Negate(is.null),rows)
      if (length(rows)) write_table(do.call(rbind,rows),file.path(cfg$out,paste0("all_",name,".csv")))
    }
    if (requireNamespace("zip",quietly=TRUE)) tryCatch({
      relative <- list.files(cfg$out,recursive=TRUE,full.names=FALSE)
      relative <- relative[grepl("\\.(csv|txt|log|R|cpp)$",relative,ignore.case=TRUE)]
      # Preserve subdirectories: avoid duplicate flattened filenames in the ZIP.
      zip::zipr(file.path(cfg$out,"results_to_review.zip"),files=relative,
        root=cfg$out,mode="mirror",include_directories=FALSE)
    },error=function(e)message("Optional review ZIP failed (saved fits retained): ",conditionMessage(e)))
    invisible(NULL)
  }
  invoke_worker <- function(cfg, fun, task = "") {
    log <- file.path(cfg$out, "logs", paste0(fun,"_",task,"_",format(Sys.time(),"%Y%m%d_%H%M%S"),".log"))
    args <- if (fun == "extract_worker") list(cfg=cfg,analysis=task) else if (fun == "fit_abundance_worker")
      list(cfg=cfg,group=task) else list(cfg=cfg)
    p <- callr::r_bg(function(engine,fun,args) {
      e <- new.env(parent=globalenv()); sys.source(engine,envir=e)
      do.call(e[[fun]],args)
    },args=list(cfg$engine,fun,args),libpath=cfg$libpath,stdout=log,stderr="2>&1",
      supervise=TRUE,user_profile=FALSE,system_profile=FALSE,wd=cfg$root,
      env=c(callr::rcmd_safe_env(),OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1"))
    on.exit(if(p$is_alive())p$kill(),add=TRUE)
    atomic_save(list(pid=p$get_pid(),fun=fun,task=task,log=log),file.path(cfg$out,"active_worker.rds"))
    p$wait(); p$get_result()
  }
  scientific_self_test <- function() {
    # Unit test population-universe preservation, onset invariance, offset sign
    # orientation, restandardization, and missing-population rejection.
    old <- data.frame(pop_id=paste0("s",1:6,"_a"),SPECIES=paste0("s",1:6),SITE_ID="a",
      voltinism=rep(GROUPS,each=3),onset_advancement_plasticity_contextual=rep(c(1,2,4),2),
      offset_termination_plasticity_contextual=c(-1,-2,-4,1,2,4),
      onset_plasticity_z=rep(as.numeric(scale(c(1,2,4))),2),
      offset_plasticity_bio_z=rep(as.numeric(scale(c(1,2,4))),2),stringsAsFactors=FALSE)
    tabs <- lapply(GROUPS,function(g) {
      a <- old[old$voltinism==g,c("pop_id","SPECIES","SITE_ID","voltinism")]
      a$offset_conditional_raw <- if(g=="univoltine")c(-1,-3,-5) else c(1,3,5)
      a
    });names(tabs)<-GROUPS
    out<-update_plasticity_table(old,tabs)$populations
    assert(identical(out$pop_id,old$pop_id) && identical(out$onset_plasticity_z,old$onset_plasticity_z) &&
      same_num(out$offset_plasticity_bio_z,rep(c(-1,0,1),2)), "Plasticity-update unit test failed.")
    broken<-tabs;broken$univoltine<-broken$univoltine[-1,,drop=FALSE]
    refused<-tryCatch({update_plasticity_table(old,broken);FALSE},error=function(e)TRUE)
    assert(refused,"Missing-population unit test failed.")
    # Annual onset is held fixed when taking the partial anomaly derivative.
    f<-function(T,onset,m) 5 + .03*onset + T*(-2 + .4*m)
    assert(abs((f(1,100,.5)-f(-1,100,.5))/2 - (-2+.4*.5))<1e-12,
      "Conditional-offset derivative test failed.")
    invisible(TRUE)
  }
  resolve_sources <- function(cfg) {
    assert(dir.exists(cfg$offset_run) && dir.exists(cfg$abundance_reference),
      "The explicit offset/reference abundance run directory is missing.")
    status <- utils::read.csv(file.path(cfg$offset_run,"workflow_status.csv"),stringsAsFactors=FALSE)
    assert(all(c("analysis","fit")%in%names(status)) &&
      all(status$fit[match(ANALYSES,status$analysis)]=="FIT_OK"),
      "Both requested onset-adjusted offset fits must be FIT_OK.")
    base <- utils::read.csv(file.path(cfg$abundance_reference,"baseline_sources.csv"),stringsAsFactors=FALSE)
    assert(all(c("group","baseline")%in%names(base)) && all(GROUPS%in%base$group),
      "Missing original abundance source manifest.")
    if(is.null(cfg$previous_plasticity_file)) {
      origin<-unique(dirname(dirname(base$baseline[match(GROUPS,base$group)])))
      assert(length(origin)==1L,"Multiple previous plasticity preparations; specify previous_plasticity_file explicitly.")
      cfg$previous_plasticity_file<-file.path(origin,"preparation_plasticity.csv")
    }
    required <- unique(c(cfg$previous_plasticity_file,
      file.path(cfg$offset_run,"workflow_status.csv"),
      unlist(lapply(ANALYSES,function(a)file.path(cfg$offset_run,a,
        c("fit_with_onset.rds","input.rds","completed.rds","parameters.rds","package_versions.csv")))),
      file.path(cfg$abundance_reference,c("run_config.rds","baseline_sources.csv","package_versions.csv","mesh.rds","compiled/trend_spde.cpp")),
      file.path(cfg$abundance_reference,paste0("input_",GROUPS,".rds")),
      file.path(cfg$abundance_reference,paste0(GROUPS,"__full"),"fit.rds")))
    absent<-required[!file.exists(required)]
    assert(!length(absent),paste("Required source file(s) missing:\n",paste(absent,collapse="\n"),
      "\nNo fallback to other runs, missing populations or non-spatial models is permitted."))
    verify_versions(file.path(cfg$abundance_reference,"package_versions.csv"),CORE)
    original_cfg<-readRDS(file.path(cfg$abundance_reference,"run_config.rds"))
    assert(isTRUE(original_cfg$mesh_cutoff_km==15),"Reference abundance mesh is not recorded as 15 km.")
    for(nm in c("optimizer_iterations","optimizer_evaluations")) {
      v<-original_cfg[[nm]]
      assert(is.numeric(v)&&length(v)==1L&&is.finite(v)&&v>0,paste("Missing reference optimizer setting:",nm))
      cfg[[nm]]<-v
    }
    # Pin the actual C++ likelihood used by run_650c41f15949.
    assert(identical(unname(tools::md5sum(file.path(cfg$abundance_reference,"compiled/trend_spde.cpp"))),
      "88e15b8506c2e85a2e1e0d3839101530"),"The reference C++ template differs from the reviewed model.")
    rows<-lapply(seq_along(required),function(i) {
      phase(cfg,"FINGERPRINTING","",paste(i,"/",length(required),basename(required[i])))
      message("Reading MD5: ",required[i])
      h<-unname(tools::md5sum(required[i]))
      assert(!is.na(h),paste("Could not fingerprint:",required[i]))
      data.frame(path=normalizePath(required[i],winslash="/",mustWork=TRUE),
        bytes=file.info(required[i])$size,md5=h,stringsAsFactors=FALSE)
    })
    manifest<-do.call(rbind,rows)
    signature<-digest::digest(list(version=VERSION,manifest=manifest[c("path","md5")],
      engine_md5=unname(tools::md5sum(cfg$engine)),
      optimizer_iterations=cfg$optimizer_iterations,optimizer_evaluations=cfg$optimizer_evaluations),algo="sha256")
    if(!is.null(cfg$signature))assert(identical(cfg$signature,signature),
      "Inputs or workflow code changed since this run. Resume refused; saved results remain untouched.")
    cfg$signature<-signature
    write_table(manifest,file.path(cfg$out,"source_manifest.csv"))
    atomic_save(cfg,file.path(cfg$out,"run_config.rds"))
    cfg
  }
  controller <- function(cfg) {
    on.exit(unlink(cfg$lock,recursive=TRUE),add=TRUE)
    tasks<-data.frame(task=c(paste0("extract_",ANALYSES),"prepare_abundance","compile",
      paste0("abundance_",GROUPS)),status="PENDING",message="",stringsAsFactors=FALSE)
    save_status<-function()write_table(tasks,file.path(cfg$out,"workflow_status.csv"))
    save_status()
    current<-NA_integer_
    result<-tryCatch({
      check_packages();scientific_self_test()
      disk_check(cfg$out,cfg$minimum_free_GiB)
      cfg<-resolve_sources(cfg)
      phase(cfg,"PREFLIGHT_COMPLETE","","Source identities and unit tests passed. No phenology refits will run.")
      steps<-list(list(fun="extract_worker",task=ANALYSES[1]),
        list(fun="extract_worker",task=ANALYSES[2]),list(fun="prepare_abundance_worker",task=""),
        list(fun="compile_worker",task=""),list(fun="fit_abundance_worker",task=GROUPS[1]),
        list(fun="fit_abundance_worker",task=GROUPS[2]))
      for(i in seq_along(steps)) {
        if(cfg$prepare_only && i>=4L) {tasks$status[i]<-"NOT_REQUESTED";next}
        current<-i; tasks$status[i]<-"RUNNING";save_status()
        disk_check(cfg$out,cfg$minimum_free_GiB)
        s<-steps[[i]]
        ans<-invoke_worker(cfg,s$fun,s$task)
        tasks$status[i]<-ans;save_status()
        tryCatch(collect(cfg),error=function(e)message("Review export issue: ",conditionMessage(e)))
      }
      ok<-all(tasks$status[5:6]%in%c("FIT_OK","FIT_OK_CACHED"))
      final<-if(cfg$prepare_only)"PREPARED_ONLY" else if(ok)"FITS_COMPLETE" else "FINISHED_WITH_NUMERICAL_REVIEW"
      phase(cfg,final,"",if(cfg$prepare_only)"New plasticity and abundance inputs saved; no abundance fit started" else
        "Review NEW offset and abundance residuals before final inference. No old diagnostics were reused.")
      list(complete=ok,output=cfg$out,status=tasks)
    },error=function(e) {
      if(!is.na(current)) {tasks$status[current]<-"ERROR";tasks$message[current]<-conditionMessage(e)} else {tasks$status[tasks$status=="PENDING"]<-"NOT_STARTED"}
      save_status()
      writeLines(conditionMessage(e),file.path(cfg$out,"controller_error.txt"))
      try(phase(cfg,"STOPPED","",conditionMessage(e)),silent=TRUE)
      message("STOPPED. Checkpoints retained; no further jobs started. ",conditionMessage(e))
      list(complete=FALSE,output=cfg$out,status=tasks,error=conditionMessage(e))
    },finally={
      try(save_status(),silent=TRUE)
      try(collect(cfg),silent=TRUE)
    })
    result
  }
  dump_engine <- function(out) {
    e<-environment(dump_engine)
    n<-ls(e,all.names=TRUE)
    f<-n[vapply(n,function(x)is.function(get(x,envir=e)),logical(1))]
    path<-file.path(out,"workflow_engine.R")
    dump(c("VERSION","ROOT","GROUPS","ANALYSES","CORE",f),file=path,envir=e,control="all")
    normalizePath(path,winslash="/",mustWork=TRUE)
  }
  get_job <- function() {
    j<-state$job
    if(is.null(j))j<-getOption("phenoimpact.updated_abundance.job")
    j
  }
  status <- function(job=get_job()) {
    assert(!is.null(job),"No updated-abundance job in this session.")
    alive<-job$process$is_alive()
    message("Elapsed: ",round(as.numeric(difftime(Sys.time(),job$started,units="hours")),2),
      " h | controller running: ",alive)
    path<-file.path(job$out,"progress.rds")
    p<-if(file.exists(path))tryCatch(readRDS(path),error=function(e)NULL)else NULL
    if(!is.null(p))message(p$stage," | ",p$task," | ",p$detail)
    path<-file.path(job$out,"workflow_status.csv")
    if(file.exists(path))try(print(utils::read.csv(path,stringsAsFactors=FALSE),row.names=FALSE),silent=TRUE)
    message("Output: ",job$out)
    if(!alive && !identical(job$process$get_exit_status(),0L))message(
      "Controller exit code: ",job$process$get_exit_status(),". Log: ",job$log)
    invisible(list(alive=alive,progress=p,output=job$out))
  }
  watch <- function(job=get_job(),every=60) {
    assert(!is.null(job),"No updated-abundance job.")
    tryCatch({repeat{status(job);if(!job$process$is_alive())break;Sys.sleep(every)}},
      interrupt=function(e)message("Monitor paused; background job CONTINUES. Keep R open."))
    invisible(job)
  }
  show_log <- function(n=30L,job=get_job()) {
    assert(!is.null(job),"No updated-abundance job.")
    path<-file.path(job$out,"active_worker.rds")
    a<-if(file.exists(path))tryCatch(readRDS(path),error=function(e)NULL)else NULL
    log<-if(!is.null(a)&&file.exists(a$log))a$log else job$log
    message("Log: ",log)
    cat(tail(readLines(log,warn=FALSE),n),sep="\n")
    invisible(log)
  }
  acquire_lock <- function(lock) {
    if(dir.exists(lock)) {
      f<-file.path(lock,"owner.rds")
      owner<-if(file.exists(f))tryCatch(readRDS(f),error=function(e)NULL)else NULL
      assert(!is.null(owner),paste("Unidentified queue lock; inspect before proceeding:",lock))
      alive<-tryCatch(ps::ps_is_running(ps::ps_handle(owner$pid,time=owner$create_time)),error=function(e)FALSE)
      assert(!alive,paste("A queue is already running; do not start a second copy:",lock))
      unlink(lock,recursive=TRUE)
    }
    assert(dir.create(lock),paste("Cannot acquire queue lock:",lock))
    h<-ps::ps_handle()
    atomic_save(list(pid=Sys.getpid(),create_time=ps::ps_create_time(h)),file.path(lock,"owner.rds"))
  }
  launch <- function(cfg,monitor=FALSE) {
    check_packages()
    disk_check(cfg$out,cfg$minimum_free_GiB)
    acquire_lock(cfg$lock)
    launched<-FALSE
    on.exit(if(!launched)unlink(cfg$lock,recursive=TRUE),add=TRUE)
    log<-file.path(cfg$out,"logs",paste0("controller_",format(Sys.time(),"%Y%m%d_%H%M%S"),".log"))
    p<-callr::r_bg(function(engine,cfg){
      e<-new.env(parent=globalenv());sys.source(engine,envir=e);e$controller(cfg)
    },args=list(cfg$engine,cfg),libpath=cfg$libpath,stdout=log,stderr="2>&1",
      supervise=TRUE,user_profile=FALSE,system_profile=FALSE,wd=cfg$root)
    h<-ps::ps_handle(p$get_pid())
    atomic_save(list(pid=p$get_pid(),create_time=ps::ps_create_time(h)),file.path(cfg$lock,"owner.rds"))
    job<-list(process=p,out=cfg$out,started=Sys.time(),log=log)
    state$job<-job;options(phenoimpact.updated_abundance.job=job);launched<-TRUE
    message("Queue launched: updated OFFSET plasticity -> two SPATIAL abundance models.")
    message("Output: ",cfg$out)
    message("Hashing/extraction run in the background; inspect pheno_updated_abundance$status().")
    message("Keep R open. Do not run another large model concurrently.")
    if(monitor)watch(job)
    invisible(job)
  }
  run <- function(root=ROOT,offset_run=NULL,abundance_reference=NULL,
                  previous_plasticity_file=NULL,out_dir=NULL,prepare_only=FALSE,
                  minimum_free_GiB=50,monitor=FALSE) {
    assert(is.logical(prepare_only)&&length(prepare_only)==1L&&!is.na(prepare_only),"Invalid prepare_only.")
    assert(is.logical(monitor)&&length(monitor)==1L&&!is.na(monitor),"Invalid monitor.")
    assert(is.numeric(minimum_free_GiB)&&length(minimum_free_GiB)==1L&&is.finite(minimum_free_GiB)&&minimum_free_GiB>=10,
      "minimum_free_GiB must be >= 10 (default 50). This is a safety reserve, not a measured disk requirement.")
    current<-get_job()
    if(!is.null(current)&&current$process$is_alive()) {
      message("This queue is already running.");if(monitor)watch(current);return(invisible(current))
    }
    check_packages()
    root<-normalizePath(root,winslash="/",mustWork=TRUE)
    if(is.null(offset_run))offset_run<-file.path(root,"output","phenology_plasticity","offset_with_onset_15km","run_3911d7bc87ac")
    if(is.null(abundance_reference))abundance_reference<-file.path(root,"output","population_trends","spatial_abundance_15km","run_650c41f15949")
    for(d in c(offset_run,abundance_reference))assert(dir.exists(d),paste("Missing source directory:",d))
    output_root<-file.path(root,"output","population_trends","spatial_abundance_onset_adjusted_15km")
    dir.create(output_root,recursive=TRUE,showWarnings=FALSE)
    if(is.null(out_dir))out_dir<-file.path(output_root,paste0("run_",format(Sys.time(),"%Y%m%d_%H%M%S"),"_",Sys.getpid()))
    assert(!file.exists(file.path(out_dir,"run_config.rds")),
      "That run already exists. Use pheno_updated_abundance$resume(run_dir=...) rather than overwriting it.")
    # Refuse output locations within either read-only source run.
    candidate<-tolower(normalizePath(out_dir,winslash="/",mustWork=FALSE))
    for(d in c(offset_run,abundance_reference)) {
      d<-tolower(normalizePath(d,winslash="/",mustWork=TRUE))
      assert(candidate!=d&&!startsWith(paste0(candidate,"/"),paste0(d,"/")),"Output cannot be inside a source run.")
    }
    dir.create(file.path(out_dir,"logs"),recursive=TRUE,showWarnings=FALSE)
    cfg<-list(version=VERSION,signature=NULL,root=root,out=normalizePath(out_dir,winslash="/",mustWork=TRUE),
      offset_run=normalizePath(offset_run,winslash="/",mustWork=TRUE),
      abundance_reference=normalizePath(abundance_reference,winslash="/",mustWork=TRUE),
      previous_plasticity_file=previous_plasticity_file,prepare_only=prepare_only,
      minimum_free_GiB=minimum_free_GiB,libpath=.libPaths(),lock=file.path(output_root,"queue.lock"))
    cfg$engine<-dump_engine(cfg$out)
    atomic_save(cfg,file.path(cfg$out,"run_config.rds"))
    writeLines(c("UPDATED PLASTICITY + SPATIAL ABUNDANCE (15 KM)",
      paste("Workflow version:",VERSION),paste("Offset source:",cfg$offset_run),
      paste("Previous SPATIAL abundance reference:",cfg$abundance_reference),
      "ONSET estimates preserved from the original preparation_plasticity.csv.",
      "NEW OFFSET slopes extracted from completed fits that include annual ONSET_mean_z.",
      "All original abundance observations, time/group/AR1 structures and the mesh are retained.",
      "New offsets standardized over the SAME previous population universe, before abundance matching.",
      "No propagation of plasticity-estimation uncertainty; no phenology refits; no LRTs.",
      "Two full spatial Gamma/log fits, all parameters re-estimated, run sequentially.",
      "TMB native environments and full joint precision matrices are not serialized.",
      "Source manifests/hashes, checks and portable optimizer checkpoints are retained.",
      "Automatic diagnostic status: PENDING for NEW offsets and NEW abundance models.",
      "FITS_COMPLETE means numerical completion, NOT final scientific validation.",
      "The observed 15/30-km sensitivity and old residual diagnostics do not validate this changed fit.",
      "Resume with pheno_updated_abundance$resume(run_dir=<this folder>, monitor=FALSE).",
      "The review ZIP uses mirror mode, preserving subdirectories and distinct filenames."),
      file.path(cfg$out,"README.txt"))
    phase(cfg,"QUEUED","","Preparing the background controller")
    launch(cfg,monitor)
  }
  resume <- function(run_dir=NULL,monitor=FALSE,prepare_only=NULL) {
    current<-get_job()
    if(!is.null(current)&&current$process$is_alive()) {
      message("This queue is already running.");return(invisible(current))
    }
    if(is.null(run_dir)) {
      assert(!is.null(current),"Supply run_dir when resuming from a new R session.")
      run_dir<-current$out
    }
    cfg<-readRDS(file.path(run_dir,"run_config.rds"))
    assert(identical(cfg$version,VERSION),"Resume with the same workflow version.")
    assert(identical(normalizePath(run_dir,winslash="/",mustWork=TRUE),cfg$out),
      "Do not move a partly completed run; its internal paths are fixed.")
    if(!is.null(prepare_only)) {
      assert(is.logical(prepare_only)&&length(prepare_only)==1L&&!is.na(prepare_only),"Invalid prepare_only.")
      cfg$prepare_only<-prepare_only
    }
    cfg$libpath<-.libPaths()
    launch(cfg,monitor)
  }
  stop_job <- function(job=get_job()) {
    assert(!is.null(job),"No updated-abundance job.")
    if(job$process$is_alive())job$process$kill()
    message("Controller stopped. Completed checkpoints retained; unsaved work can be lost.")
    invisible(job)
  }
  list(run=run,resume=resume,status=status,watch=watch,log=show_log,stop=stop_job,self_test=scientific_self_test)
})
# No autorun: call pheno_updated_abundance$run(monitor = FALSE) explicitly.

message("Functions loaded. Start with pheno_updated_abundance$run(monitor = FALSE).")
