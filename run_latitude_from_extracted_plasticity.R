# ============================================================================
# phenoIMPACT | Latitude patterns in EXTRACTED contextual plasticity
# 2026-10-05 | v1.0.0 | Base R >= 4.1; no package installation
#
# source(file.choose(), encoding = "UTF-8")
# latitud <- run_latitude_from_plasticities()
#
# READS CSVs ONLY. No readRDS(), sdmTMB(), TMB(), likelihoods or large checkpoints.
# Sources: original spatial ONSET; updated, onset-adjusted spatial OFFSETs.
# The original anomaly x latitude models are NOT fitted or used as responses.
# Their site_latitudes.csv files supply coordinates only.
#
# SECOND-STAGE ESTIMAND: linear geographical summaries of model-predicted
# contextual sensitivities over the SAMPLED population contexts.
# Global and between-species total slopes are DESCRIPTIVE (no naive SE/p-value).
# Within-species slope uncertainty is propagated analytically from the saved
# fixed-effect covariance: species slope deviations cancel on demeaning.
# These are not independent observations, causal latitude effects, or new
# estimates of population random slopes. No second spatial field is fitted.
#
# All component populations are retained for each response, not just the subset
# with usable abundance. The combined updated CSV is a consistency check.
# No anomaly back-transformation is guessed: units remain days/fitted anomaly unit.
# No multiplicity-adjusted tests or new response/residual uncertainty are claimed.
# Output folders are new and input files are never overwritten.
# ============================================================================

pheno_latitude_from_plasticities <- local({
  VERSION <- "latitude_extracted_plasticity_1.0.0"
  ROOT <- "E:/phenoIMPACT project/code/phenoIMPACT"
  ANALYSES <- c("onset", "offset_univoltine", "offset_multivoltine")
  MAX_CSV_MIB <- 256
  need <- function(ok, text) if (!isTRUE(ok)) stop(text, call. = FALSE)
  csv <- function(x, path) utils::write.csv(x, path, row.names = FALSE, na = "")
  msg <- function(...) message(format(Sys.time(), "%H:%M:%S"), " | ", ...)
  pathnorm <- function(x, must = TRUE) normalizePath(x, winslash = "/", mustWork = must)
  numeric_equal <- function(a, b, tol = 1e-6) {
    length(a) == length(b) && all(is.finite(a)) && all(is.finite(b)) &&
      max(abs(a - b), 0) <= tol * max(1, abs(a), abs(b))
  }
  true_values <- function(x) {
    length(x) > 0L && !anyNA(x) && all(toupper(as.character(x)) %in% c("TRUE", "1"))
  }
  read_csv <- function(path) {
    need(file.exists(path), paste("Required CSV not found:", path))
    z <- file.info(path)
    need(!is.na(z$size) && z$size > 0 && z$size <= MAX_CSV_MIB * 1024^2,
      paste("Refusing empty/oversized CSV:", path))
    header <- utils::read.csv(path, nrows = 0, check.names = FALSE)
    cc <- rep(NA_character_, ncol(header))
    cc[names(header) %in% c("SPECIES", "SITE_ID", "bms_id", "pop_id", "term")] <- "character"
    utils::read.csv(path, colClasses = cc, stringsAsFactors = FALSE,
      check.names = FALSE, na.strings = c("NA", ""))
  }
  require_columns <- function(x, cols, where) {
    miss <- setdiff(cols, names(x))
    need(!length(miss), paste(where, "is missing:", paste(miss, collapse = ", ")))
  }
  remap_project_path <- function(path, root) {
    path <- gsub("\\\\", "/", as.character(path))
    if (file.exists(path)) return(pathnorm(path))
    # Remapping preserves the exact run/file beneath output, not 'latest'.
    k <- regexpr("/output/", path, fixed = TRUE)
    if (length(path) == 1L && k > 0L) {
      candidate <- file.path(root, substring(path, k + 1L))
      if (file.exists(candidate)) return(pathnorm(candidate))
    }
    stop("Pinned source path does not exist: ", path, call. = FALSE)
  }

  locate <- function(root = ROOT, updated_run = NULL, onset_dir = NULL,
                     coordinates_files = NULL) {
    root <- pathnorm(root)
    if (is.null(updated_run)) updated_run <- file.path(root, "output", "population_trends",
      "spatial_abundance_onset_adjusted_15km", "run_20261004_214012_44136")
    need(dir.exists(updated_run), paste("Pinned updated run not found:", updated_run,
      "\nUse updated_run= only to name an explicitly verified replacement; no latest-run fallback."))
    updated_run <- pathnorm(updated_run)
    main <- file.path(updated_run, "phenological_population_plasticity_onset_adjusted_15km.csv")
    need(file.exists(main), paste("Updated consistency table not found:", main))
    manifest_path <- NULL
    if (is.null(onset_dir)) {
      manifest_path <- file.path(root, "output", "population_trends", "spatial_plasticity_15km",
        "run_e30730052267", "input_sources.csv")
      m <- read_csv(manifest_path)
      require_columns(m, c("input", "path"), manifest_path)
      i <- which(m$input == "spatial_plasticity")
      need(length(i) == 1L, "Original source manifest does not identify exactly one spatial-plasticity CSV.")
      original_csv <- remap_project_path(m$path[i], root)
      onset_dir <- file.path(dirname(original_csv), "onset")
    }
    need(dir.exists(onset_dir), paste("Original onset extraction directory not found:", onset_dir,
      "\nSpecify onset_dir= with the exact original extraction/onset folder, not a phenology model folder."))
    dirs <- c(onset = pathnorm(onset_dir),
      offset_univoltine = file.path(updated_run, "extraction", "offset_univoltine"),
      offset_multivoltine = file.path(updated_run, "extraction", "offset_multivoltine"))
    files <- list()
    for (a in ANALYSES) {
      files[[a]] <- setNames(file.path(dirs[[a]], c("population_plasticity_components.csv",
        "fixed_effects.csv", "fixed_effect_covariance.csv", "provenance.csv")),
        c("components", "fixed", "covariance", "provenance"))
    }
    if (is.null(coordinates_files)) {
      cr <- file.path(root, "output", "phenology_plasticity", "latitude_models_15km")
      coordinates_files <- list.files(cr, pattern = "^site_latitudes\\.csv$",
        full.names = TRUE, recursive = TRUE)
      coordinates_files <- coordinates_files[grepl("^primary__", basename(dirname(coordinates_files)))]
    }
    need(length(coordinates_files) > 0L,
      "No site_latitudes.csv found. Supply coordinates_files= with CSVs containing SITE_ID, bms_id and latitude_deg.")
    all_paths <- c(main, manifest_path, unlist(files, use.names = FALSE), coordinates_files)
    missing <- all_paths[!file.exists(all_paths)]
    need(!length(missing), paste("Required CSVs are missing:\n", paste(missing, collapse = "\n")))
    sizes <- file.info(all_paths)$size
    need(all(is.finite(sizes) & sizes > 0 & sizes <= MAX_CSV_MIB * 1024^2) &&
      sum(sizes) <= 768 * 1024^2, "Input CSV size guard triggered; no files were loaded.")
    list(root = root, updated_run = updated_run, main = main, original_manifest = manifest_path,
      files = files, coordinates = sort(unique(pathnorm(coordinates_files))),
      all_paths = unique(pathnorm(all_paths)))
  }

  key <- function(site, network) {
    need(!anyNA(site) && !anyNA(network) && all(nzchar(site)) && all(nzchar(network)),
      "Missing/empty site or network identities.")
    need(!any(grepl("\r", site, fixed = TRUE)) && !any(grepl("\r", network, fixed = TRUE)),
      "Unexpected separator in an identity.")
    paste(site, network, sep = "\r")
  }
  read_coordinates <- function(paths) {
    rows <- lapply(paths, function(f) {
      d <- read_csv(f)
      require_columns(d, c("SITE_ID", "bms_id", "latitude_deg"), f)
      d <- d[c("SITE_ID", "bms_id", "latitude_deg")]
      need(is.numeric(d$latitude_deg) && all(is.finite(d$latitude_deg)) &&
        all(abs(d$latitude_deg) <= 90), paste("Invalid latitude:", f))
      d
    })
    all <- do.call(rbind, rows)
    k <- key(all$SITE_ID, all$bms_id)
    first <- !duplicated(k)
    # Match into the unique reference, not the first rows of the full table.
    ref <- all$latitude_deg[first][match(k, k[first])]
    need(all(abs(all$latitude_deg - ref) <= 1e-7),
      "Coordinate exports disagree for the same site/network. No coordinates were averaged.")
    all <- all[first, , drop = FALSE]
    rownames(all) <- NULL
    all
  }

  read_analysis <- function(a, paths, main, coords, out) {
    p <- read_csv(paths[["components"]]); fe <- read_csv(paths[["fixed"]])
    vc <- read_csv(paths[["covariance"]]); prov <- read_csv(paths[["provenance"]])
    window <- if (a == "onset") 60L else 90L
    anomaly <- paste0("clim_anomaly_tw", window)
    mods <- paste0(c("photo_tw", "clim_background_tw", "clim_predictability_tw", "clim_trend_tw"), window)
    require_columns(p, c("pop_id", "SPECIES", "SITE_ID", "bms_id", "n_years",
      "plasticity_raw", "sp_slope_dev", "pop_slope_dev", mods), paths[["components"]])
    require_columns(fe, c("term", "estimate", "std.error"), paths[["fixed"]])
    require_columns(vc, "term", paths[["covariance"]])
    require_columns(prov, c("analysis", "formula", "offset_adjusted_for_onset"), paths[["provenance"]])
    need(nrow(prov) == 1L && identical(as.character(prov$analysis), a), "Analysis provenance mismatch.")
    need(!is.na(prov$formula) && !grepl("latitude", prov$formula, ignore.case = TRUE),
      "These inputs appear to come from the old latitude model, not the environmental model.")
    if (a == "onset") {
      require_columns(prov, "cutoff_km", "Onset provenance")
      need(isTRUE(prov$cutoff_km == 15), "Onset extraction is not marked as 15 km.")
    } else {
      require_columns(p, c("offset_adjusted_for_onset", "offset_conditional_raw"), "Updated offset table")
      require_columns(prov, "source_run", "Updated offset provenance")
      need(true_values(prov$offset_adjusted_for_onset) && true_values(p$offset_adjusted_for_onset),
        "Offset inputs are NOT explicitly marked as adjusted for annual onset.")
      need(basename(gsub("\\\\", "/", prov$source_run)) == "run_3911d7bc87ac",
        "Offset source is not the pinned run_3911d7bc87ac. Review an intended source change explicitly.")
      need(numeric_equal(p$plasticity_raw, p$offset_conditional_raw), "Offset raw/conditional columns disagree.")
      need("ONSET_mean_z" %in% fe$term &&
        !any(grepl("ONSET_mean_z:|:ONSET_mean_z", fe$term)),
        "Expected an additive annual ONSET_mean_z covariate, without onset interactions.")
    }
    need(nrow(p) > 1L && !anyNA(p[c("pop_id", "SPECIES", "SITE_ID", "bms_id")]) &&
      !anyDuplicated(p$pop_id) && !anyDuplicated(p[c("SPECIES", "SITE_ID")]),
      "Population rows must be unique species/site estimates, not annual observations.")
    need(all(nzchar(p$SPECIES)) && all(p$pop_id == paste(p$SPECIES, p$SITE_ID, sep = "_")),
      "Unexpected population identities.")
    vals <- c("plasticity_raw", "sp_slope_dev", "pop_slope_dev", "n_years", mods)
    need(all(vapply(p[vals], is.numeric, logical(1))) &&
      all(is.finite(as.matrix(p[vals]))) && all(p$n_years >= 1), "Invalid numeric population data.")
    need(all(p$pop_slope_dev == 0), "Population random slopes are present; this analytical covariance propagation no longer applies.")
    dev <- tapply(p$sp_slope_dev, p$SPECIES, function(z) diff(range(z)))
    need(all(dev < 1e-8), "Species slope deviations are not constant within species.")
    need(!anyNA(fe$term) && !anyDuplicated(fe$term) && all(is.finite(fe$estimate)) &&
      all(is.finite(fe$std.error) & fe$std.error > 0), "Invalid fixed-effect table.")
    beta <- stats::setNames(fe$estimate, fe$term)
    need(!anyDuplicated(vc$term) && setequal(vc$term, names(beta)) &&
      all(names(beta) %in% names(vc)), "Covariance coefficient names do not match the fixed effects.")
    V <- as.matrix(vc[match(names(beta), vc$term), names(beta), drop = FALSE])
    storage.mode(V) <- "double"; dimnames(V) <- list(names(beta), names(beta))
    need(all(is.finite(V)) && max(abs(V - t(V))) <= 1e-7 * max(1, abs(V)),
      "Non-finite/asymmetric covariance matrix.")
    V <- (V + t(V)) / 2
    eig <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
    need(min(eig) >= -1e-9 * max(1, abs(eig)), "Fixed-effect covariance is not positive semidefinite.")
    need(numeric_equal(sqrt(pmax(diag(V), 0)), fe$std.error), "Covariance diagonal does not reproduce fixed-effect SEs.")
    need(anomaly %in% names(beta), "Expected anomaly coefficient missing.")
    D <- matrix(0, nrow(p), length(beta), dimnames = list(NULL, names(beta)))
    D[, anomaly] <- 1
    for (m in mods) {
      term <- intersect(c(paste(anomaly, m, sep = ":"), paste(m, anomaly, sep = ":")), names(beta))
      need(length(term) == 1L, paste("Missing/ambiguous moderator interaction:", m))
      D[, term] <- p[[m]]
    }
    reconstructed <- as.numeric(D %*% beta) + p$sp_slope_dev
    need(numeric_equal(reconstructed, p$plasticity_raw),
      "Extracted plasticities do not reproduce fixed + species + environmental components.")
    main_rows <- if (a == "onset") seq_len(nrow(main)) else which(main$voltinism == sub("offset_", "", a, fixed = TRUE))
    i <- match(main$pop_id[main_rows], p$pop_id)
    need(length(i) > 0 && !anyNA(i), "Some populations in the updated consistency CSV are missing from this extraction.")
    target <- if (a == "onset") -main$onset_advancement_plasticity_contextual[main_rows] else
      main$offset_termination_plasticity_contextual[main_rows]
    need(identical(p$SPECIES[i], main$SPECIES[main_rows]) && identical(p$SITE_ID[i], main$SITE_ID[main_rows]) &&
      numeric_equal(p$plasticity_raw[i], target), "Extraction and current combined plasticity CSV disagree.")
    j <- match(key(p$SITE_ID, p$bms_id), key(coords$SITE_ID, coords$bms_id))
    if (anyNA(j)) {
      csv(p[is.na(j), c("pop_id", "SPECIES", "SITE_ID", "bms_id")], file.path(out, "missing_coordinates.csv"))
      stop("Coordinates are missing for some populations; see missing_coordinates.csv. No rows were dropped.", call. = FALSE)
    }
    p$latitude_deg <- coords$latitude_deg[j]
    p$analysis <- a
    p$plasticity_oriented <- (if (a == "offset_multivoltine") 1 else -1) * p$plasticity_raw
    list(analysis = a, p = p, D = D, beta = beta, V = V, provenance = prov,
      n_combined_check = length(i), n_extra_vs_combined = nrow(p) - length(i),
      reconstruction_error = max(abs(reconstructed - p$plasticity_raw)))
  }

  grouped_means <- function(M, group) {
    if (is.null(dim(M))) M <- matrix(M, ncol = 1L)
    sums <- rowsum(M, group = group, reorder = FALSE)
    index <- match(group, rownames(sums))
    n <- tabulate(index, nbins = nrow(sums))
    means <- sweep(sums, 1L, n, "/")
    list(mean = means, index = index, n = n, levels = rownames(sums),
      centred = M - means[index, , drop = FALSE])
  }
  propagation <- function(C, V) {
    variance <- as.numeric(crossprod(C, V %*% C))
    need(is.finite(variance) && variance >= -1e-9, "Negative/non-finite propagated slope variance.")
    sqrt(max(0, variance))
  }
  projection <- function(x, y, D, beta, V, w = rep(1, length(x)), within = FALSE) {
    need(length(x) == length(y) && nrow(D) == length(x) && length(w) == length(x) &&
      all(is.finite(x)) && all(is.finite(y)) && all(is.finite(w) & w > 0), "Invalid projection inputs.")
    xc <- if (within) x else x - stats::weighted.mean(x, w)
    denom <- sum(w * xc^2)
    need(is.finite(denom) && denom > 1e-12, "Insufficient latitudinal variation for this projection.")
    u <- w * xc / denom
    C <- colSums(D * u)
    estimate <- sum(u * y)
    fixed <- sum(C * beta)
    if (within) need(numeric_equal(estimate, fixed), "Within-species cancellation check failed.")
    se <- if (within) propagation(C, V) else NA_real_
    list(estimate = estimate, fixed_component = fixed, species_component = estimate - fixed,
      SE = se, lower = estimate - stats::qnorm(.975) * se,
      upper = estimate + stats::qnorm(.975) * se, contrast = C,
      intercept = if (within) 0 else stats::weighted.mean(y, w) - estimate * stats::weighted.mean(x, w),
      denominator = denom, weights = w, u = u)
  }
  row_result <- function(a, name, weighting, fit, n_pop, n_species, n_units, inference) {
    sg <- if (a == "offset_multivoltine") 1 else -1
    lo <- if (sg == 1) fit$lower else -fit$upper
    hi <- if (sg == 1) fit$upper else -fit$lower
    data.frame(analysis = a, contrast = name, weighting = weighting,
      n_populations = n_pop, n_species = n_species, n_units_in_projection = n_units,
      estimate_raw_per_10deg = fit$estimate, SE_source_model = fit$SE,
      lower_95_raw = fit$lower, upper_95_raw = fit$upper,
      estimate_oriented_per_10deg = sg * fit$estimate,
      lower_95_oriented = lo, upper_95_oriented = hi,
      fixed_context_component_raw = fit$fixed_component,
      species_deviation_component_raw = fit$species_component,
      CI_excludes_zero = if (is.finite(fit$SE)) fit$lower > 0 || fit$upper < 0 else NA,
      interval_scope = inference,
      units = "change in days/fitted-anomaly-unit per 10 degrees latitude",
      stringsAsFactors = FALSE)
  }
  save_plot <- function(path, f) {
    grDevices::png(path, width = 1750, height = 1350, res = 210)
    on.exit(grDevices::dev.off(), add = TRUE)
    graphics::par(mar = c(4.8, 5.0, 1.2, 1.0), las = 1, family = "serif", bty = "l")
    f()
  }
  make_plots <- function(p, species, within_x, within_y, fits, out) {
    rawlab <- "Predicted sensitivity (days / fitted anomaly unit)"
    save_plot(file.path(out, "global_latitude.png"), function() {
      graphics::plot(p$latitude_deg, p$plasticity_raw, pch = 16, cex = .45,
        col = grDevices::adjustcolor("black", alpha.f = .14), xlab = "Latitude (degrees N)", ylab = rawlab)
      xx <- range(p$latitude_deg)
      graphics::lines(xx, fits$global$intercept + fits$global$estimate * xx / 10, lwd = 2)
    })
    save_plot(file.path(out, "within_species_latitude.png"), function() {
      graphics::plot(within_x * 10, within_y, pch = 16, cex = .45,
        col = grDevices::adjustcolor("black", alpha.f = .14),
        xlab = "Latitude minus sampled species mean (degrees)",
        ylab = "Predicted sensitivity minus species mean")
      xx <- seq(min(within_x), max(within_x), length.out = 201)
      lo <- pmin(xx * fits$within$lower, xx * fits$within$upper)
      hi <- pmax(xx * fits$within$lower, xx * fits$within$upper)
      graphics::polygon(c(xx * 10, rev(xx * 10)), c(lo, rev(hi)),
        col = grDevices::adjustcolor("grey50", alpha.f = .25), border = NA)
      graphics::abline(h = 0, v = 0, lty = 3)
      graphics::lines(xx * 10, xx * fits$within$estimate, lwd = 2)
    })
    save_plot(file.path(out, "between_species_latitude.png"), function() {
      graphics::plot(species$mean_latitude_deg, species$mean_plasticity_raw,
        pch = 16, cex = .7, xlab = "Mean latitude of sampled populations (degrees N)", ylab = rawlab)
      xx <- range(species$mean_latitude_deg)
      graphics::lines(xx, fits$between$intercept + fits$between$estimate * xx / 10, lwd = 2)
    })
  }

  analyse <- function(obj, out, make_figures = TRUE) {
    p <- obj$p; D <- obj$D; beta <- obj$beta; V <- obj$V; a <- obj$analysis
    means <- grouped_means(cbind(latitude = p$latitude_deg / 10, plasticity = p$plasticity_raw, D), p$SPECIES)
    xw <- means$centred[, 1]; yw <- means$centred[, 2]
    Dw <- means$centred[, -(1:2), drop = FALSE]
    weights_balanced <- 1 / means$n[means$index]
    ns <- length(means$n); np <- nrow(p)
    fits <- list(
      global = projection(p$latitude_deg / 10, p$plasticity_raw, D, beta, V),
      global_balanced = projection(p$latitude_deg / 10, p$plasticity_raw, D, beta, V, weights_balanced),
      between = projection(means$mean[, 1], means$mean[, 2], means$mean[, -(1:2), drop = FALSE], beta, V),
      within = projection(xw, yw, Dw, beta, V, within = TRUE),
      within_balanced = projection(xw, yw, Dw, beta, V, weights_balanced, within = TRUE))
    # Uncertainty from species deviations, including covariance with beta, cancels
    # exactly because the within contrast sums to zero for EACH species.
    need(max(abs(rowsum(fits$within$u, p$SPECIES, reorder = FALSE))) < 1e-8 &&
      max(abs(rowsum(fits$within_balanced$u, p$SPECIES, reorder = FALSE))) < 1e-8,
      "A within contrast does not sum to zero within each species.")
    desc <- "DESCRIPTIVE ONLY: total-slope joint beta/species covariance unavailable; SE/CI intentionally omitted"
    inf <- "Source fixed-effect covariance; conditional on observed contexts and fitted model; NOT independent new latitude evidence"
    result <- rbind(
      row_result(a, "global", "equal_population", fits$global, np, ns, np, desc),
      row_result(a, "global", "equal_species_total_weight", fits$global_balanced, np, ns, np, desc),
      row_result(a, "between_species", "equal_species", fits$between, np, ns, ns, desc),
      row_result(a, "within_species", "equal_population", fits$within, np, ns, np, inf),
      row_result(a, "within_species", "equal_species_total_weight", fits$within_balanced, np, ns, np, inf))
    contrasts <- do.call(rbind, lapply(names(fits), function(n) data.frame(
      analysis = a, projection = n, term = names(beta), contrast_loading = fits[[n]]$contrast,
      fixed_estimate = beta, contribution_raw_per_10deg = fits[[n]]$contrast * beta)))
    p$species_mean_latitude_deg <- means$mean[means$index, 1] * 10
    p$latitude_within_species_deg <- xw * 10
    p$species_mean_plasticity_raw <- means$mean[means$index, 2]
    p$plasticity_within_species_raw <- yw
    p$species_balanced_weight <- weights_balanced
    species <- do.call(rbind, lapply(seq_len(ns), function(k) {
      ix <- which(means$index == k); den <- sum(xw[ix]^2)
      slope <- se <- NA_real_
      if (den > 1e-12) {
        C <- colSums(Dw[ix, , drop = FALSE] * xw[ix]) / den
        slope <- sum(C * beta); se <- propagation(C, V)
      }
      data.frame(analysis = a, SPECIES = means$levels[k], n_populations = length(ix),
        n_sites = length(unique(p$SITE_ID[ix])), min_latitude_deg = min(p$latitude_deg[ix]),
        max_latitude_deg = max(p$latitude_deg[ix]), latitude_span_deg = diff(range(p$latitude_deg[ix])),
        mean_latitude_deg = means$mean[k, 1] * 10, mean_plasticity_raw = means$mean[k, 2],
        median_years_per_population = stats::median(p$n_years[ix]), min_years_per_population = min(p$n_years[ix]),
        supplies_within_latitude_information = den > 1e-12,
        within_design_information_fraction = den / fits$within$denominator,
        predicted_within_slope_raw_per_10deg = slope, SE_source_model = se,
        lower_95_raw = slope - stats::qnorm(.975) * se, upper_95_raw = slope + stats::qnorm(.975) * se)
    }))
    # Leave-one-species-out changes the geographical SUMMARY only. No source fit.
    total_num <- colSums(Dw * xw); total_den <- sum(xw^2)
    loo <- do.call(rbind, lapply(seq_len(ns), function(k) {
      ix <- which(means$index == k); den <- total_den - sum(xw[ix]^2)
      estimate <- se <- NA_real_
      if (den > 1e-12) {
        C <- (total_num - colSums(Dw[ix, , drop = FALSE] * xw[ix])) / den
        estimate <- sum(C * beta); se <- propagation(C, V)
      }
      data.frame(analysis = a, omitted_species = means$levels[k], removed_populations = length(ix),
        estimate_raw_per_10deg = estimate, SE_source_model = se,
        lower_95_raw = estimate - stats::qnorm(.975) * se, upper_95_raw = estimate + stats::qnorm(.975) * se,
        source_model_refitted = FALSE)
    }))
    bins <- split(seq_len(np), 5 * floor(p$latitude_deg / 5))
    support <- do.call(rbind, lapply(names(bins), function(n) {
      ix <- bins[[n]]
      data.frame(analysis = a, band_start_degrees = as.numeric(n), band_end_degrees = as.numeric(n) + 5,
        n_populations = length(ix), n_species = length(unique(p$SPECIES[ix])), n_sites = length(unique(p$SITE_ID[ix])),
        mean_predicted_sensitivity_raw = mean(p$plasticity_raw[ix]), median_years = stats::median(p$n_years[ix]))
    }))
    audit <- data.frame(analysis = a, n_populations = np, n_species = ns,
      n_sites = length(unique(p$SITE_ID)), min_latitude_deg = min(p$latitude_deg), max_latitude_deg = max(p$latitude_deg),
      n_species_within_support = sum(species$supplies_within_latitude_information),
      n_populations_in_combined_consistency_check = obj$n_combined_check,
      n_additional_component_populations_retained = obj$n_extra_vs_combined,
      maximum_reconstruction_error = obj$reconstruction_error,
      offset_adjusted_for_onset = a != "onset", population_rows_dropped = 0L,
      within_total_equals_context_only = TRUE, no_source_model_refit = TRUE)
    csv(result, file.path(out, "latitude_summary.csv"))
    csv(contrasts, file.path(out, "linear_contrast_coefficients.csv"))
    csv(p, file.path(out, "population_plasticity_with_latitude.csv"))
    csv(species, file.path(out, "species_latitude_support_and_slopes.csv"))
    csv(support, file.path(out, "latitude_bands_support.csv"))
    csv(loo, file.path(out, "leave_one_species_out_summary.csv"))
    csv(audit, file.path(out, "validation.csv"))
    if (make_figures) tryCatch(make_plots(p, species, xw, yw, fits, out), error = function(e) {
      writeLines(conditionMessage(e), file.path(out, "figure_error.txt"))
      warning("Numerical summaries saved; figure failed: ", conditionMessage(e), call. = FALSE)
    })
    list(summary = result, validation = audit)
  }

  self_test <- function() {
    # Exact constructed case: species shifts cannot affect centred slopes.
    s <- rep(c("A", "B", "C"), c(3, 4, 5))
    x <- c(3, 4, 5, 4, 5, 6, 7, 5, 5.5, 6, 6.5, 7)
    si <- match(s, c("A", "B", "C"))
    m <- 2 * x + c(-4, 0, 9)[si]
    D <- cbind(anomaly = 1, interaction = m)
    beta <- c(anomaly = -2, interaction = .3)
    V <- matrix(c(.04, .006, .006, .01), 2, dimnames = list(names(beta), names(beta)))
    y <- as.numeric(D %*% beta) + c(5, -3, 7)[si]
    g <- grouped_means(cbind(x, y, D), s)
    f <- projection(g$centred[, 1], g$centred[, 2], g$centred[, -(1:2), drop = FALSE], beta, V, within = TRUE)
    fb <- projection(g$centred[, 1], g$centred[, 2], g$centred[, -(1:2), drop = FALSE], beta, V,
      1 / g$n[g$index], within = TRUE)
    need(abs(f$estimate - .6) < 1e-10 && abs(f$SE - .2) < 1e-10 &&
      abs(fb$estimate - .6) < 1e-10 && abs(fb$SE - .2) < 1e-10,
      "Internal analytical self-test failed.")
    # Compare only the point coefficient with an explicit species fixed-effect fit.
    ordinary <- stats::lm(y ~ x + factor(s))
    need(abs(unname(stats::coef(ordinary)["x"]) - f$estimate) < 1e-10,
      "Demeaned and species fixed-effect projections differ in the internal self-test.")
    # A large additional species-only component must leave both within results unchanged.
    y2 <- y + c(900, -450, 1300)[si]
    g2 <- grouped_means(cbind(x, y2, D), s)
    f2 <- projection(g2$centred[, 1], g2$centred[, 2], g2$centred[, -(1:2), drop = FALSE], beta, V, within = TRUE)
    need(abs(f2$estimate - f$estimate) < 1e-10, "Species-deviation cancellation self-test failed.")
    invisible(TRUE)
  }

  zip_output <- function(out) {
    target <- file.path(out, "latitude_patterns_to_review.zip")
    paths <- list.files(out, recursive = TRUE, full.names = TRUE, all.files = FALSE)
    paths <- paths[!dir.exists(paths) & !grepl("\\.zip$", paths, ignore.case = TRUE)]
    relative <- substring(paths, nchar(out) + 2L)
    # Mirror mode preserves directories and repeated basenames.
    ok <- tryCatch({
      if (requireNamespace("zip", quietly = TRUE)) {
        zip::zip(target, files = relative, root = out, mode = "mirror")
      } else if (.Platform$OS.type == "windows") {
        # The temporary destination is outside the source directory.
        tmp <- paste0(out, "_review.zip")
        quote_ps <- function(x) paste0("'", gsub("'", "''", x, fixed = TRUE), "'")
        command <- paste0("Add-Type -AssemblyName System.IO.Compression.FileSystem; ",
          "[System.IO.Compression.ZipFile]::CreateFromDirectory(", quote_ps(pathnorm(out)), ", ",
          quote_ps(pathnorm(tmp, FALSE)), ", [System.IO.Compression.CompressionLevel]::Optimal, $false)")
        status <- system2("powershell", c("-NoProfile", "-NonInteractive", "-Command", shQuote(command)), stdout = TRUE, stderr = TRUE)
        code <- attr(status, "status")
        need(is.null(code) || code == 0L, paste("Archive command failed:", paste(status, collapse = "\n")))
        need(file.exists(tmp) && file.rename(tmp, target), "Cannot finalize review ZIP.")
      } else {
        need(nzchar(Sys.which("zip")), "No ZIP utility available.")
        old <- getwd(); on.exit(setwd(old), add = TRUE); setwd(out)
        relative <- substring(paths, nchar(out) + 2L)
        code <- utils::zip(target, files = relative, flags = "-q")
        need(is.null(code) || code == 0L, "System zip failed.")
      }
      need(file.exists(target) && file.info(target)$size > 0, "Review ZIP was not created.")
      contents <- utils::unzip(target, list = TRUE)$Name
      contents <- contents[!grepl("/$", contents)]
      need(!anyDuplicated(contents) && setequal(contents, relative),
        "ZIP entries differ from expected relative paths; do not use this archive.")
      TRUE
    }, error = function(e) {
      warning("Results are saved, but ZIP creation failed: ", conditionMessage(e), call. = FALSE)
      FALSE
    })
    if (ok) pathnorm(target) else NA_character_
  }

  run <- function(root = ROOT, updated_run = NULL, onset_dir = NULL,
                  coordinates_files = NULL, out_dir = NULL, collect_only = FALSE,
                  make_figures = TRUE) {
    need(is.logical(collect_only) && length(collect_only) == 1L && !is.na(collect_only), "Invalid collect_only.")
    need(is.logical(make_figures) && length(make_figures) == 1L && !is.na(make_figures), "Invalid make_figures.")
    self_test()
    msg("Locating the pinned CSV inputs. No phenology model is loaded or refitted.")
    src <- locate(root, updated_run, onset_dir, coordinates_files)
    if (is.null(out_dir)) {
      base <- file.path(src$root, "output", "phenology_plasticity", "latitude_from_extracted_plasticity_15km")
      dir.create(base, recursive = TRUE, showWarnings = FALSE)
      out_dir <- tempfile(paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_"), tmpdir = base)
    }
    need(!dir.exists(out_dir) && !file.exists(out_dir), "Output already exists. Use a new directory; nothing will be overwritten.")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    need(dir.exists(out_dir), "Could not create output directory.")
    out <- pathnorm(out_dir)
    complete <- FALSE
    on.exit(if (!complete) writeLines("STOPPED before successful completion. Read the error; do not interpret partial outputs as complete.",
      file.path(out, "WORKFLOW_STOPPED.txt")), add = TRUE)
    initial_hashes <- unname(tools::md5sum(src$all_paths))
    need(!anyNA(initial_hashes), "Could not fingerprint an input CSV.")
    manifest <- data.frame(file = src$all_paths, bytes = file.info(src$all_paths)$size,
      modified = as.character(file.info(src$all_paths)$mtime), md5 = initial_hashes)
    csv(manifest, file.path(out, "input_manifest.csv"))
    dirs <- c(ANALYSES, "inputs/coordinates", paste0("inputs/", ANALYSES), "inputs/consistency")
    for (d in dirs) dir.create(file.path(out, d), recursive = TRUE, showWarnings = FALSE)
    copied <- list()
    cp <- function(source, destination) {
      need(file.copy(source, destination, overwrite = FALSE), paste("Could not collect CSV:", source))
      need(identical(unname(tools::md5sum(source)), unname(tools::md5sum(destination))), "Collected CSV hash mismatch.")
      copied[[length(copied) + 1L]] <<- data.frame(source = source,
        review_relative_path = substring(destination, nchar(out) + 2L))
    }
    for (a in ANALYSES) for (f in src$files[[a]]) cp(f, file.path(out, "inputs", a, basename(f)))
    cp(src$main, file.path(out, "inputs", "consistency", basename(src$main)))
    if (!is.null(src$original_manifest)) cp(src$original_manifest,
      file.path(out, "inputs", "consistency", "original_input_sources.csv"))
    for (i in seq_along(src$coordinates)) cp(src$coordinates[i],
      file.path(out, "inputs", "coordinates", paste0(sprintf("%02d_", i), "site_latitudes.csv")))
    csv(do.call(rbind, copied), file.path(out, "collected_input_files.csv"))
    main <- read_csv(src$main)
    require_columns(main, c("pop_id", "SPECIES", "SITE_ID", "voltinism",
      "onset_advancement_plasticity_contextual", "offset_termination_plasticity_contextual",
      "offset_adjusted_for_onset"), "Current combined plasticity table")
    need(!anyNA(main[c("pop_id", "SPECIES", "SITE_ID", "voltinism")]) &&
      !anyDuplicated(main$pop_id) && all(main$voltinism %in% c("univoltine", "multivoltine")) &&
      true_values(main$offset_adjusted_for_onset), "Combined updated table fails identity/offset-adjustment checks.")
    coords <- read_coordinates(src$coordinates)
    csv(coords, file.path(out, "verified_site_latitudes.csv"))
    summaries <- validations <- list()
    for (a in ANALYSES) {
      msg(a, ": validating components, current estimates, coordinates and covariance.")
      target <- file.path(out, a)
      obj <- read_analysis(a, src$files[[a]], main, coords, target)
      if (!collect_only) {
        result <- analyse(obj, target, make_figures)
        summaries[[a]] <- result$summary; validations[[a]] <- result$validation
        msg(a, ": geographical projections saved (", nrow(obj$p), " populations).")
      }
      rm(obj); invisible(gc())
    }
    summary <- if (length(summaries)) do.call(rbind, summaries) else data.frame()
    if (length(summaries)) {
      csv(summary, file.path(out, "all_latitude_summaries.csv"))
      csv(do.call(rbind, validations), file.path(out, "all_validation_checks.csv"))
    }
    notes <- c(
      "phenoIMPACT: latitude patterns of model-predicted contextual plasticity", paste("Version:", VERSION), "",
      "Inputs are CSVs only. No checkpoint is read, copied, modified or refitted.",
      paste("Pinned updated run:", src$updated_run),
      "ONSET: original 15-km spatial environmental extraction, 60-day anomaly window.",
      "OFFSET: both updated environmental fits with additive annual ONSET_mean_z, 90-day window.",
      "All populations in each component extraction are used, even when absent from abundance models.",
      "The updated combined CSV checks consistency; it does not define the geographical-analysis sample.",
      "Means give each sampled population one vote within a species, irrespective of years followed.",
      "No years/abundance-based weighting, new filters, or causal interpretation is applied.", "",
      "MODELS / PROJECTIONS", "x = latitude in degrees / 10; y = raw contextual sensitivity.",
      "Global: y ~ 1 + x (equal population; sensitivity with 1/n_population_species weights).",
      "Between species: species mean y ~ 1 + species mean x (one row/vote per species).",
      "Within species: (y - species mean y) ~ 0 + (x - species mean x).",
      "Within sensitivity: each species has total weight 1 via population weights 1/n_species_populations.",
      "Equal total species weight does not remove the larger slope leverage of wider latitude ranges.",
      "Species with no within-species latitude variation do not inform the within slope.",
      "Within coefficients match a regression with species-specific intercepts for these weights.", "",
      "UNCERTAINTY", "Each extracted slope is y = D beta + species_slope_deviation.",
      "D contains 1 for the anomaly term and population means for its moderator interactions.",
      "Within-species centring removes the constant anomaly slope and species deviations exactly.",
      "For the within projection, g = C beta, C = sum_i(w_i x_within_i D_within_i)/sum_i(w_i x_within_i^2).",
      "Var(g) = C V_beta C'; 95% intervals use a normal approximation and the FULL saved fixed-effect covariance.",
      "No independent noise is drawn for each population, and no naive regression SE/p-value is reported.",
      "For the specified within contrast, joint covariance with species deviations is not needed: their coefficients are zero.",
      "These intervals are conditional on sampled locations, population contexts, fixed model specification and fitted covariance parameters.",
      "They do not include uncertainty in phenology extraction, moderator measurements, window/model selection or unsampled populations.",
      "Global/between TOTAL estimates include species deviations; complete joint covariance is not exported, so their SE/CI are intentionally blank.",
      "Do not interpret the within projection as new independent evidence that latitude causes plasticity.",
      "It geographically summarizes relationships already contained in the environmental phenology model.",
      "Leave-one-species-out removes species from the SUMMARY, not from the original phenology fit.",
      "No tests compare onset and offsets; those fits have no exported cross-model covariance.", "",
      "SIGNS / FIGURES", "Raw sensitivity < 0 means earlier phenology in warmer-anomaly years; > 0 means later.",
      "Oriented exports multiply onset/univoltine offset by -1 and multivoltine offset by +1.",
      "Thus positive oriented values mean stronger advancement (onset/uni offset) or delay (multi offset), not necessarily adaptive responses.",
      "No absolute values are used. Offset sensitivities hold annual onset fixed.",
      "Temperature units are inherited. Do not relabel as days/degree C without verifying original scaling.",
      "Figures use raw sensitivity. The within figure is species-centred in BOTH axes; its ribbon is a model-derived projection interval.",
      "Global and between figures intentionally have no uncertainty ribbons.", "",
      "QUALITY CONTROL", "Internal analytical self-test passed before reading project data.",
      "Input identities, sign conventions, coefficient reconstruction, covariance and coordinates are checked; no silent row drops.",
      "These checks do not replace residual diagnostics/validation of the source phenology models.",
      "References: van de Pol & Wright (2009), Animal Behaviour 77:753-758, doi:10.1016/j.anbehav.2008.11.006.",
      "Houslay & Wilson (2017), Behavioral Ecology 28:948-952, doi:10.1093/beheco/arx023.")
    writeLines(notes, file.path(out, "README_METHODS.txt"))
    writeLines(capture.output(utils::sessionInfo()), file.path(out, "sessionInfo.txt"))
    need(identical(initial_hashes, unname(tools::md5sum(src$all_paths))), "An input changed during execution. Outputs require review.")
    writeLines(if (collect_only) "INPUTS_COLLECTED_AND_VALIDATED; NO GEOGRAPHICAL PROJECTIONS RUN" else
      "GEOGRAPHICAL_PROJECTIONS_COMPLETE; NOT A NEW VALIDATION OF SOURCE PHENOLOGY FITS", file.path(out, "STATUS.txt"))
    complete <- TRUE
    archive <- zip_output(out)
    cat("\nOutput directory:\n", out, "\n\n", sep = "")
    if (nrow(summary)) print(summary[c("analysis", "contrast", "weighting", "estimate_raw_per_10deg",
      "lower_95_raw", "upper_95_raw")], row.names = FALSE)
    if (!is.na(archive)) cat("\nFILE TO SEND:\n", archive, "\n", sep = "") else
      cat("\nCompress the output directory above and send that ZIP. Numerical outputs remain saved.\n")
    invisible(list(output_directory = out, zip_to_send = archive, summary = summary))
  }
  list(run = run, locate = locate, self_test = self_test)
})

run_latitude_from_plasticities <- pheno_latitude_from_plasticities$run
