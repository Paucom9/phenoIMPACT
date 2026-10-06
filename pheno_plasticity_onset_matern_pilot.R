# phenoIMPACT: onset / annual Matern field pilot
# Prepared 2026-09-21. Run from your existing R project with source(this_file).
#
# INPUT: output/phenology_plasticity/plasticity_main_models.rds
# Set input_rds below if the file is elsewhere. All fitted rows are retained.
# The saved onset formula, climatic window and predictor scales are preserved.
#
# COMPARISON (both fitted by ML in sdmTMB):
#   site_year: original fixed/random terms + site x year random intercept
#   spatial_field: the same + an independent Matern field for each year
# The annual field is shared by species, with a common range and variance over
# years; it is not a spatially varying anomaly slope or a latitude analysis.
#
# OUTPUT: diagnostics/spatial_field_pilot/onset_mesh15km_<fingerprint>/
#   comparison.xlsx              main tables (see README sheet)
#   fit_site_year.rds             checkpoint saved before fitting the next model
#   fit_spatial_field.rds         checkpoint, also saved if convergence fails
#   moran_all_tests.csv           all tests, including skipped groups
#   focal_coefficients.csv        estimates, SEs and normal Wald 95% CIs
#   residuals_site_year.rds       response + mle-mvn residuals, in fitted row order
#   residuals_spatial_field.rds   same for spatial model
#   coefficient_comparison.png, moran_comparison.png, residual_checks.pdf
#   mesh.png, metadata.rds, sessionInfo.txt, validation.csv, fit_checks.csv
#
# RERUN: identical data, model settings and relevant package versions reuse
# checkpoints. Diagnostic settings can change without refitting models.
# Use check_only = TRUE for input/mesh checks before starting lengthy fits.
# A finer mesh (e.g. mesh_cutoff_km = 10) creates a separate output directory.
# A 15-km cutoff is a pilot choice, not a demonstrated adequate resolution.
#
# INTERPRETATION:
# - Check fit_checks and spatial_parameters before interpreting coefficients.
# - Compare effect estimates/CIs, not only changes in significance or AIC.
# - Response residuals permit comparison with July diagnostics. mle-mvn uses
#   one prespecified posterior draw of random effects, as recommended by
#   sdmTMB; its diagnostic result is stochastic. Do not select a favourable seed.
# - Moran screens use species-mean residuals per site/year, then average sites
#   with identical coordinates. Unequal species counts can affect variances;
#   permutation p-values are exploratory, not a calibrated model-level test.
# - BH is calculated over valid BMS/year tests per variant/residual type/band,
#   and also over all bands within each variant/residual type.
# - Aggregation can conceal species-specific spatial dependence. This pilot
#   evaluates the shared signal; it does not establish adequacy for every species.
# - Lower residual Moran I alone does not validate the spatial model. Inspect
#   residual plots, coefficient changes, variance/range uncertainty and mesh
#   sensitivity before extending to other responses. No naive chi-square LRT
#   of zero field variance is used (variance lies on a boundary).
#
# API references checked when preparing this script:
# https://sdmtmb.github.io/sdmTMB/reference/make_mesh.html
# https://sdmtmb.github.io/sdmTMB/reference/residuals.sdmTMB.html
# https://sdmtmb.github.io/sdmTMB/reference/tidy.sdmTMB.html
# https://sdmtmb.github.io/sdmTMB/reference/sanity.html
# https://r-spatial.github.io/spdep/reference/moran.mc.html

pilot_config <- list(
  input_rds = NULL,                 # NULL = the project path shown above
  output_root = NULL,               # NULL = project diagnostics directory
  mesh_cutoff_km = 15,
  check_only = FALSE,
  force_refit = FALSE,
  gradient_threshold = 0.001,
  extra_optimization_attempts = 2L,
  minimum_connected_sites = 10L,
  n_permutations = 1999L,
  fit_seed = 20260921L,
  residual_seed = 20260922L,
  permutation_seed = 20260923L,
  distance_bands = data.frame(
    distance_band = c("0-25 km", "25-50 km", "50-100 km", "100-200 km"),
    d_min_km = c(0, 25, 50, 100), d_max_km = c(25, 50, 100, 200)
  )
)

# Dependencies are checked, never installed automatically.
check_packages <- function() {
  pkgs <- c("lme4", "lmerTest", "sdmTMB", "spdep", "dplyr", "tidyr",
            "tibble", "ggplot2", "writexl", "here", "digest")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Missing R packages. Install them first:\ninstall.packages(c(",
         paste(sprintf('"%s"', missing), collapse = ", "), "))")
  }
  invisible(pkgs)
}

write_table <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE, na = "")
}

normal_ci <- function(x) {
  required <- c("term", "estimate", "std.error")
  if (!all(required %in% names(x))) stop("Unexpected fixed-effect table.")
  x <- x[, required]
  names(x)[names(x) == "std.error"] <- "SE"
  z <- stats::qnorm(0.975)
  x$lower_95 <- x$estimate - z * x$SE
  x$upper_95 <- x$estimate + z * x$SE
  x$CI_excludes_zero <- x$lower_95 > 0 | x$upper_95 < 0
  x$interval_method <- "Normal Wald 95%; not an LRT"
  tibble::as_tibble(x)
}

same_values <- function(a, b) {
  if (is.numeric(a) && is.numeric(b)) {
    isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = 1e-10))
  } else identical(as.character(a), as.character(b))
}

prepare_onset <- function(bundle) {
  entry <- bundle$primary$onset
  if (is.null(entry) || !inherits(entry$model, "merMod")) {
    stop("Expected a fitted lme4 model at bundle$primary$onset$model.")
  }
  original <- entry$model
  fixed_formula <- lme4::nobars(stats::formula(original))
  mf <- stats::model.frame(original)
  raw <- as.data.frame(entry$data)
  required <- unique(c(all.vars(stats::formula(original)), "bms_id", "YEAR"))
  if (!all(required %in% names(raw))) stop("Saved data lack required columns.")
  if (!identical(all.vars(fixed_formula[[2]]), "ONSET_mean")) {
    stop("The primary onset response is not ONSET_mean.")
  }

  index <- match(rownames(mf), rownames(raw))
  if (anyNA(index)) {
    if (nrow(raw) != nrow(mf)) stop("Cannot align saved and fitted rows.")
    index <- seq_len(nrow(raw))
  }
  d <- raw[index, required, drop = FALSE]
  # Equal row counts are not enough: verify response, predictors AND group IDs.
  shared <- intersect(names(mf), names(d))
  bad <- shared[!vapply(shared, function(nm) same_values(mf[[nm]], d[[nm]]),
                       logical(1))]
  if (length(bad)) stop("Saved/fitted rows disagree for: ", paste(bad, collapse = ", "))
  old_X <- lme4::getME(original, "X")
  new_X <- stats::model.matrix(fixed_formula, data = d)
  if (!identical(colnames(old_X), colnames(new_X)) ||
      !isTRUE(all.equal(unname(old_X), unname(new_X), check.attributes = FALSE,
                       tolerance = 1e-10))) {
    stop("Fixed-effect design matrix differs from the fitted model.")
  }
  d$source_row <- index
  d$SITE_ID <- as.character(d$SITE_ID)
  d$SPECIES <- as.character(d$SPECIES)
  d$bms_id <- as.character(d$bms_id)
  yr <- suppressWarnings(as.numeric(as.character(d$YEAR)))
  if (any(!is.finite(yr)) || any(yr != floor(yr))) stop("YEAR must contain integer years.")
  d$YEAR <- as.integer(yr)
  if (anyNA(d) || any(!nzchar(d$SITE_ID)) || any(!nzchar(d$bms_id))) {
    stop("Missing data or empty IDs in fitted rows; do not silently drop rows.")
  }
  if (anyDuplicated(d[, c("SITE_ID", "SPECIES", "YEAR")])) {
    stop("Duplicate species/site/year rows: inspect the input before fitting.")
  }

  coords <- as.data.frame(bundle$site_coordinates)
  if (!all(c("SITE_ID", "bms_id", "x_3035", "y_3035") %in% names(coords))) {
    stop("Expected site_coordinates with EPSG:3035 coordinates in metres.")
  }
  coords <- coords |>
    dplyr::transmute(SITE_ID = as.character(SITE_ID), bms_id = as.character(bms_id),
                     x_3035 = as.numeric(x_3035), y_3035 = as.numeric(y_3035)) |>
    dplyr::distinct()
  if (anyDuplicated(coords[, c("SITE_ID", "bms_id")])) {
    stop("Conflicting coordinate records for the same site/network.")
  }
  sites <- dplyr::distinct(d, SITE_ID, bms_id)
  if (anyDuplicated(sites$SITE_ID)) {
    stop("SITE_ID occurs in multiple networks; review the original site grouping.")
  }
  d <- dplyr::left_join(d, coords, by = c("SITE_ID", "bms_id"))
  if (nrow(d) != nrow(mf) || any(!is.finite(d$x_3035)) ||
      any(!is.finite(d$y_3035))) stop("Coordinate join failed or changed fitted rows.")
  if (max(abs(d$x_3035)) <= 180 && max(abs(d$y_3035)) <= 90) {
    stop("Coordinates look like longitude/latitude, not EPSG:3035 metres.")
  }
  d$SITE_ID <- factor(d$SITE_ID)
  d$SPECIES <- factor(d$SPECIES)
  # A second grouping-column name explicitly preserves independent species
  # intercepts/slopes even across sdmTMB versions with different term handling.
  d$SPECIES_slope <- d$SPECIES
  d$site_year_id <- interaction(d$SITE_ID, d$YEAR, drop = TRUE, lex.order = TRUE)
  d$x_km <- d$x_3035 / 1000
  d$y_km <- d$y_3035 / 1000
  anomaly <- paste0("clim_anomaly_tw", entry$best_window)
  if (!anomaly %in% names(d)) stop("Best-window anomaly column not found.")

  # Verify the random structure rather than silently replacing an unknown one.
  cnms <- lme4::getME(original, "cnms")
  actual <- sort(paste(names(cnms), vapply(cnms, paste, character(1), collapse = ","),
                       sep = ":"))
  expected <- sort(c("SITE_ID:(Intercept)", "SPECIES:(Intercept)",
                     paste0("SPECIES:", anomaly)))
  if (!identical(actual, expected)) {
    stop("Unexpected original random-effects structure: ", paste(actual, collapse = "; "))
  }
  rhs <- paste(deparse(fixed_formula[[3]], width.cutoff = 500L), collapse = " ")
  model_formula <- stats::as.formula(paste0(
    "ONSET_mean ~ ", rhs, " + (1 | SITE_ID) + (1 | SPECIES) + (0 + ",
    anomaly, " | SPECIES_slope) + (1 | site_year_id)"
  ), env = baseenv())
  beta <- lme4::fixef(original)
  se <- sqrt(diag(as.matrix(stats::vcov(original))))
  original_coefficients <- normal_ci(data.frame(
    term = names(beta), estimate = unname(beta), std.error = unname(se)))
  residual <- as.numeric(stats::residuals(original))
  if (length(residual) != nrow(d) || any(!is.finite(residual))) {
    stop("Original residuals do not align with fitted rows.")
  }
  list(data = as.data.frame(d), formula = model_formula,
       original_formula = paste(deparse(stats::formula(original)), collapse = " "),
       best_window = entry$best_window, anomaly = anomaly,
       original_coefficients = original_coefficients, original_residuals = residual)
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

fit_or_resume <- function(variant, d, formula, mesh, run_id, out_dir, cfg) {
  path <- file.path(out_dir, paste0("fit_", variant, ".rds"))
  if (file.exists(path) && !cfg$force_refit) {
    message("Reusing saved fit: ", variant)
    saved <- readRDS(path)
    if (!identical(saved$run_id, run_id) || !identical(saved$variant, variant)) {
      stop("Checkpoint does not match current data/settings.")
    }
    saved$checks <- fit_checks(saved$fit, cfg)
    return(saved)
  }
  messages <- character()
  capture_warning <- function(w) {
    messages <<- c(messages, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
  start <- proc.time()[["elapsed"]]
  set.seed(cfg$fit_seed)
  message("Fitting ", variant, " with ", nrow(d), " observations...")
  fit <- withCallingHandlers(sdmTMB::sdmTMB(
    formula = formula, data = d, mesh = mesh, time = "YEAR",
    family = stats::gaussian(link = "identity"), spatial = "off",
    spatiotemporal = if (variant == "spatial_field") "iid" else "off",
    reml = FALSE, silent = FALSE
  ), warning = capture_warning)
  checks <- fit_checks(fit, cfg)
  attempts <- 0L
  while (!checks$fit_ok && attempts < cfg$extra_optimization_attempts) {
    attempts <- attempts + 1L
    message("Extra optimization ", attempts, " for ", variant)
    # Avoid dense Newton Hessians for this large analysis.
    improved <- tryCatch(withCallingHandlers(sdmTMB::run_extra_optimization(
      fit, nlminb_loops = 1, newton_loops = 0), warning = capture_warning),
      error = function(e) {
        messages <<- c(messages, paste("Extra optimization failed:", conditionMessage(e)))
        NULL
      })
    # Retain and save the available fit if an extra optimization step errors.
    if (is.null(improved)) break
    fit <- improved
    checks <- fit_checks(fit, cfg)
  }
  saved <- list(run_id = run_id, variant = variant, fit = fit, checks = checks,
                warnings = unique(messages), extra_optimization = attempts,
                elapsed_minutes = (proc.time()[["elapsed"]] - start) / 60)
  saveRDS(saved, path, compress = "gzip")
  writeLines(unique(messages), file.path(out_dir, paste0("warnings_", variant, ".txt")))
  saved
}

aligned_residuals <- function(fit, d, seed) {
  rows <- fit$data$source_row
  if (is.null(rows) || anyDuplicated(rows) || length(rows) != nrow(d)) {
    stop("The fitted object lacks unique source-row identifiers.")
  }
  index <- match(d$source_row, rows)
  if (anyNA(index)) stop("Fitted residual rows cannot be matched to the input.")
  response <- as.numeric(stats::residuals(fit, type = "response"))
  set.seed(seed)
  mvn <- as.numeric(stats::residuals(fit, type = "mle-mvn"))
  if (length(response) != length(rows) || length(mvn) != length(rows)) {
    stop("Residual lengths differ from the fitted data.")
  }
  list(response = response[index], `mle-mvn` = mvn[index])
}

aggregate_residuals <- function(d, residual, variant, residual_type) {
  if (length(residual) != nrow(d) || any(!is.finite(residual))) {
    stop("Non-finite or misaligned residuals for ", variant, " / ", residual_type)
  }
  r <- d[, c("SITE_ID", "bms_id", "YEAR", "x_km", "y_km")]
  r$residual <- as.numeric(residual)
  r |>
    dplyr::group_by(bms_id, YEAR, SITE_ID, x_km, y_km) |>
    dplyr::summarise(residual = mean(residual), n_observations = dplyr::n(),
                     .groups = "drop") |>
    # Keep site means equally weighted when multiple IDs share coordinates.
    dplyr::group_by(bms_id, YEAR, x_km, y_km) |>
    dplyr::summarise(residual = mean(residual), n_sites = dplyr::n(),
                     n_observations = sum(n_observations), .groups = "drop") |>
    dplyr::mutate(model_variant = variant, residual_type = residual_type)
}

bh_valid <- function(p) {
  out <- rep(NA_real_, length(p))
  valid <- is.finite(p)
  out[valid] <- stats::p.adjust(p[valid], method = "BH")
  out
}

moran_one <- function(d, lo, hi, seed, cfg) {
  d <- d[order(d$x_km, d$y_km), ]
  result <- data.frame(n_locations_total = nrow(d), n_locations_used = NA_integer_,
                       n_links = NA_integer_, n_components = NA_integer_,
                       expected_I = NA_real_, moran_I = NA_real_,
                       excess_I = NA_real_, p_permutation = NA_real_,
                       status = "not_run", notes = "")
  if (nrow(d) < cfg$minimum_connected_sites) {
    result$status <- "too_few_locations"
    return(result)
  }
  notes <- character()
  out <- tryCatch(withCallingHandlers({
    xy <- as.matrix(d[, c("x_km", "y_km")])
    # First band [0,25]; later bands (25,50], etc.: no overlapping boundaries.
    bounds <- c(if (lo == 0) "GE" else "GT", "LE")
    nb <- spdep::dnearneigh(xy, d1 = lo, d2 = hi, longlat = FALSE, bounds = bounds)
    keep <- spdep::card(nb) > 0L
    result$n_locations_used <- sum(keep)
    if (sum(keep) < cfg$minimum_connected_sites) {
      result$status <- "too_few_connected_locations"
      return(result)
    }
    xy <- xy[keep, , drop = FALSE]
    values <- d$residual[keep]
    nb <- spdep::dnearneigh(xy, d1 = lo, d2 = hi, longlat = FALSE, bounds = bounds)
    result$n_links <- as.integer(sum(spdep::card(nb)) / 2)
    result$n_components <- spdep::n.comp.nb(nb)$nc
    if (result$n_links < cfg$minimum_connected_sites) {
      result$status <- "too_few_links"
      return(result)
    }
    if (!is.finite(stats::sd(values)) || stats::sd(values) == 0) {
      result$status <- "zero_residual_variance"
      return(result)
    }
    weights <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
    set.seed(seed)
    test <- spdep::moran.mc(values, listw = weights, nsim = cfg$n_permutations,
                          alternative = "greater", zero.policy = TRUE)
    result$expected_I <- -1 / (length(values) - 1)
    result$moran_I <- as.numeric(test$statistic)
    result$excess_I <- result$moran_I - result$expected_I
    result$p_permutation <- test$p.value
    result$status <- "OK"
    result
  }, warning = function(w) {
    notes <<- c(notes, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) {
    result$status <- paste0("ERROR: ", conditionMessage(e))
    result
  })
  out$notes <- paste(unique(notes), collapse = " | ")
  out
}

run_moran <- function(site_residuals, cfg) {
  group_cols <- c("model_variant", "residual_type", "bms_id", "YEAR")
  groups <- site_residuals |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    tidyr::nest() |>
    dplyr::ungroup()
  output <- vector("list", nrow(groups) * nrow(cfg$distance_bands))
  k <- 0L
  for (i in seq_len(nrow(groups))) {
    for (j in seq_len(nrow(cfg$distance_bands))) {
      k <- k + 1L
      band <- cfg$distance_bands[j, ]
      # Same permutations for matched network/year/band, independent of loop order.
      key <- paste(cfg$permutation_seed, groups$bms_id[i], groups$YEAR[i],
                    band$distance_band, sep = "|")
      seed <- strtoi(substr(digest::digest(key, algo = "xxhash32"), 1, 7), base = 16L)
      output[[k]] <- dplyr::bind_cols(groups[i, group_cols], band,
        moran_one(groups$data[[i]], band$d_min_km, band$d_max_km, seed, cfg))
    }
    if (i %% 25L == 0L) message("Moran groups: ", i, " / ", nrow(groups))
  }
  dplyr::bind_rows(output) |>
    dplyr::group_by(model_variant, residual_type, distance_band) |>
    dplyr::mutate(p_FDR_band = bh_valid(p_permutation)) |>
    dplyr::ungroup() |>
    dplyr::group_by(model_variant, residual_type) |>
    dplyr::mutate(p_FDR_all_bands = bh_valid(p_permutation)) |>
    dplyr::ungroup()
}

summarise_moran <- function(tests) {
  summary <- tests |>
    dplyr::filter(status == "OK") |>
    dplyr::group_by(model_variant, residual_type, distance_band) |>
    dplyr::summarise(n_tests = dplyr::n(), median_I = stats::median(moran_I),
      median_excess_I = stats::median(excess_I),
      n_positive_FDR_band = sum(excess_I > 0 & p_FDR_band < 0.05),
      n_positive_FDR_all_bands = sum(excess_I > 0 & p_FDR_all_bands < 0.05),
      .groups = "drop")
  keep <- c("residual_type", "bms_id", "YEAR", "distance_band",
            "n_locations_used", "excess_I", "p_FDR_band", "p_FDR_all_bands")
  base <- tests |>
    dplyr::filter(model_variant == "site_year", status == "OK") |>
    dplyr::select(dplyr::all_of(keep))
  spatial <- tests |>
    dplyr::filter(model_variant == "spatial_field", status == "OK") |>
    dplyr::select(dplyr::all_of(keep))
  paired <- dplyr::inner_join(base, spatial,
    by = c("residual_type", "bms_id", "YEAR", "distance_band"),
    suffix = c("_site_year", "_spatial")) |>
    dplyr::mutate(reduction_in_excess_I = excess_I_site_year - excess_I_spatial)
  if (any(paired$n_locations_used_site_year != paired$n_locations_used_spatial)) {
    stop("Paired Moran tests used different locations.")
  }
  paired_summary <- paired |>
    dplyr::group_by(residual_type, distance_band) |>
    dplyr::summarise(n_paired_tests = dplyr::n(),
      median_excess_I_site_year = stats::median(excess_I_site_year),
      median_excess_I_spatial = stats::median(excess_I_spatial),
      median_reduction = stats::median(reduction_in_excess_I),
      proportion_reduced = mean(reduction_in_excess_I > 0),
      n_positive_FDR_site_year = sum(excess_I_site_year > 0 & p_FDR_band_site_year < 0.05),
      n_positive_FDR_spatial = sum(excess_I_spatial > 0 & p_FDR_band_spatial < 0.05),
      .groups = "drop")
  list(moran_summary = summary, paired_BMS_year = paired, paired_summary = paired_summary)
}

plot_residual_checks <- function(d, vectors, variant) {
  for (type in names(vectors)) {
    r <- vectors[[type]]
    # Quantiles avoid plotting hundreds of thousands of overlapping QQ points.
    p <- seq(0.001, 0.999, length.out = 999)
    graphics::plot(stats::qnorm(p), stats::quantile(r, p), pch = 16, cex = 0.5,
                   main = paste(variant, type), xlab = "Standard normal quantile",
                   ylab = "Residual quantile")
    if (type == "mle-mvn") graphics::abline(0, 1, col = "red") else {
      qs <- stats::quantile(r, c(0.25, 0.75))
      slope <- diff(qs) / diff(stats::qnorm(c(0.25, 0.75)))
      graphics::abline(qs[1] - slope * stats::qnorm(0.25), slope, col = "red")
    }
  }
  fitted <- d$ONSET_mean - vectors$response
  i <- unique(round(seq(1, nrow(d), length.out = min(nrow(d), 30000L))))
  graphics::plot(fitted[i], vectors$response[i], pch = 16, cex = 0.3,
                 col = grDevices::adjustcolor("black", alpha.f = 0.15),
                 main = paste(variant, "response residuals"),
                 xlab = "Fitted onset (day of year)", ylab = "Observed - fitted (days)")
  graphics::abline(h = 0, col = "red")
  graphics::plot(d$YEAR[i], vectors[["mle-mvn"]][i], pch = 16, cex = 0.3,
                 col = grDevices::adjustcolor("black", alpha.f = 0.15),
                 main = paste(variant, "mle-mvn by year"),
                 xlab = "Year", ylab = "mle-mvn residual")
  graphics::abline(h = 0, col = "red")
}

run_onset_pilot <- function(cfg = pilot_config) {
  pkgs <- check_packages()
  if (!is.numeric(cfg$mesh_cutoff_km) || length(cfg$mesh_cutoff_km) != 1L ||
      !is.finite(cfg$mesh_cutoff_km) || cfg$mesh_cutoff_km <= 0) stop("Invalid mesh cutoff.")
  if (cfg$n_permutations < 99 || cfg$n_permutations != as.integer(cfg$n_permutations)) {
    stop("n_permutations must be an integer >= 99.")
  }
  path <- cfg$input_rds
  if (is.null(path)) path <- here::here("output", "phenology_plasticity",
                                       "plasticity_main_models.rds")
  if (!file.exists(path)) stop("Input not found: ", path, "\nSet pilot_config$input_rds.")
  message("Reading the saved model bundle; this may require several GB of RAM.")
  bundle <- readRDS(path)
  prepared <- prepare_onset(bundle)
  bundle_created <- bundle$created_at
  rm(bundle)
  invisible(gc())
  d <- prepared$data
  package_versions <- data.frame(package = pkgs,
    version = vapply(pkgs, function(p) as.character(utils::packageVersion(p)), character(1)))
  # Deliberately exclude diagnostic settings so models can be reused.
  signature <- list(code_version = "onset_matern_pilot_1.0", data = d,
    formula = paste(deparse(prepared$formula), collapse = " "),
    cutoff = cfg$mesh_cutoff_km, seed = cfg$fit_seed,
    gradient_threshold = cfg$gradient_threshold,
    extra_optimization_attempts = cfg$extra_optimization_attempts,
    packages = package_versions[package_versions$package %in% c("sdmTMB", "lme4"), ],
    TMB_version = as.character(utils::packageVersion("TMB")),
    fmesher_version = as.character(utils::packageVersion("fmesher")))
  run_id <- digest::digest(signature, algo = "sha256")
  rm(signature)
  root <- cfg$output_root
  if (is.null(root)) root <- here::here("output", "phenology_plasticity", "diagnostics",
                                       "spatial_field_pilot")
  out_dir <- file.path(root, paste0("onset_mesh", cfg$mesh_cutoff_km, "km_", substr(run_id, 1, 10)))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
  unique_xy <- dplyr::distinct(d, x_km, y_km)
  set.seed(cfg$fit_seed)
  triangulation <- sdmTMB::make_mesh(unique_xy, xy_cols = c("x_km", "y_km"),
                                    cutoff = cfg$mesh_cutoff_km, type = "cutoff")
  # CRITICAL FIX: the observation-to-mesh projection must contain EVERY fitted
  # row in its actual order, not just one row per unique location.
  mesh <- sdmTMB::make_mesh(d, xy_cols = c("x_km", "y_km"), mesh = triangulation$mesh)
  if (nrow(mesh$loc_xy) != nrow(d) ||
      !isTRUE(all.equal(unname(as.matrix(mesh$loc_xy)),
                       unname(as.matrix(d[, c("x_km", "y_km")])),
                       check.attributes = FALSE))) stop("Mesh/data alignment failed.")
  grDevices::png(file.path(out_dir, "mesh.png"), width = 1800, height = 1500, res = 180)
  tryCatch(plot(triangulation), finally = grDevices::dev.off())
  validation <- data.frame(model = "primary__onset", best_window = prepared$best_window,
    n_observations = nrow(d), n_species = nlevels(d$SPECIES), n_sites = nlevels(d$SITE_ID),
    n_site_years = nlevels(d$site_year_id), n_unique_coordinates = nrow(unique_xy),
    n_mesh_vertices = mesh$mesh$n, mesh_cutoff_km = cfg$mesh_cutoff_km,
    first_year = min(d$YEAR), last_year = max(d$YEAR), coordinates = "EPSG:3035 / km",
    fitted_rows_and_design_verified = TRUE)
  write_table(validation, file.path(out_dir, "validation.csv"))
  saveRDS(list(run_id = run_id, input = normalizePath(path), bundle_created = bundle_created,
    settings = cfg, validation = validation, original_formula = prepared$original_formula,
    formula = prepared$formula, package_versions = package_versions,
    mesh = triangulation, source_rows = d$source_row),
    file.path(out_dir, "metadata.rds"), compress = "gzip")
  print(validation)
  if (cfg$check_only) {
    message("Input and mesh checks completed. Set check_only = FALSE to fit.\n", out_dir)
    return(invisible(out_dir))
  }

  coefficients <- list(original_lmer = dplyr::mutate(prepared$original_coefficients,
                                                     model_variant = "original_lmer"))
  aggregated <- list(aggregate_residuals(d, prepared$original_residuals,
                                        "original_lmer", "response"))
  checks_all <- list()
  random_parameters <- list()
  grDevices::pdf(file.path(out_dir, "residual_checks.pdf"), width = 10, height = 8)
  pdf_device <- grDevices::dev.cur()
  on.exit({
    if (pdf_device %in% grDevices::dev.list()) grDevices::dev.off(pdf_device)
  }, add = TRUE)
  graphics::par(mfrow = c(2, 2))
  for (variant in c("site_year", "spatial_field")) {
    saved <- fit_or_resume(variant, d, prepared$formula, mesh, run_id, out_dir, cfg)
    fit <- saved$fit
    check <- saved$checks
    check$model_variant <- variant
    check$elapsed_minutes <- saved$elapsed_minutes
    check$extra_optimization <- saved$extra_optimization
    check$logLik <- as.numeric(stats::logLik(fit))
    check$AIC <- stats::AIC(fit)
    checks_all[[variant]] <- check
    write_table(dplyr::bind_rows(checks_all), file.path(out_dir, "fit_checks.csv"))
    if (!check$fit_ok) {
      stop("Fit did not pass numerical checks: ", variant,
           ". Saved its checkpoint and fit_checks.csv in ", out_dir,
           ". Inspect these before interpreting results or refitting.")
    }
    if (!check$all_sanity_checks_pass) {
      warning(variant, ": some additional sanity checks failed; inspect sanity output and parameters.")
    }
    writeLines(capture.output(sdmTMB::sanity(fit, gradient_thresh = cfg$gradient_threshold)),
               file.path(out_dir, paste0("sanity_", variant, ".txt")))
    coefficients[[variant]] <- normal_ci(sdmTMB::tidy(fit, effects = "fixed", conf.int = FALSE)) |>
      dplyr::mutate(model_variant = variant)
    if (!setequal(coefficients[[variant]]$term, coefficients$original_lmer$term)) {
      stop("Fitted fixed-effect terms differ from the saved original model.")
    }
    random_parameters[[variant]] <- sdmTMB::tidy(fit, effects = "ran_pars", conf.int = TRUE) |>
      dplyr::mutate(model_variant = variant)
    vectors <- aligned_residuals(fit, d, cfg$residual_seed)
    for (type in names(vectors)) {
      aggregated[[length(aggregated) + 1L]] <- aggregate_residuals(d, vectors[[type]], variant, type)
    }
    saveRDS(list(run_id = run_id, seed = cfg$residual_seed, source_rows = d$source_row,
                  residuals = vectors), file.path(out_dir, paste0("residuals_", variant, ".rds")),
            compress = "gzip")
    plot_residual_checks(d, vectors, variant)
    rm(fit, saved, vectors)
    invisible(gc())
  }
  grDevices::dev.off(pdf_device)
  model_fit <- dplyr::bind_rows(checks_all)
  model_fit$delta_AIC_from_best <- model_fit$AIC - min(model_fit$AIC)
  model_fit$delta_AIC_spatial_minus_site_year <-
    model_fit$AIC[model_fit$model_variant == "spatial_field"] -
    model_fit$AIC[model_fit$model_variant == "site_year"]
  write_table(model_fit, file.path(out_dir, "fit_checks.csv"))
  coefficients <- dplyr::bind_rows(coefficients)
  coefficients$focal <- coefficients$term == prepared$anomaly |
    (grepl(prepared$anomaly, coefficients$term, fixed = TRUE) &
       grepl(":", coefficients$term, fixed = TRUE))
  focal <- dplyr::filter(coefficients, focal)
  base <- focal |>
    dplyr::filter(model_variant == "site_year") |>
    dplyr::select(-model_variant, -focal, -interval_method)
  spatial <- focal |>
    dplyr::filter(model_variant == "spatial_field") |>
    dplyr::select(-model_variant, -focal, -interval_method)
  changes <- dplyr::inner_join(base, spatial, by = "term", suffix = c("_site_year", "_spatial")) |>
    dplyr::mutate(estimate_change = estimate_spatial - estimate_site_year,
                  SE_ratio = SE_spatial / SE_site_year,
                  sign_retained = sign(estimate_spatial) == sign(estimate_site_year))
  write_table(focal, file.path(out_dir, "focal_coefficients.csv"))
  site_residuals <- dplyr::bind_rows(aggregated)
  saveRDS(site_residuals, file.path(out_dir, "aggregated_residuals.rds"), compress = "gzip")
  message("Computing spatial diagnostic screens...")
  tests <- run_moran(site_residuals, cfg)
  write_table(tests, file.path(out_dir, "moran_all_tests.csv"))
  summaries <- summarise_moran(tests)
  if (any(grepl("^ERROR", tests$status))) warning("Some Moran tests failed; inspect moran_all_tests.")
  if (!nrow(summaries$paired_summary)) warning("No valid paired Moran tests; spatial adequacy is unresolved.")

  # Include the original estimates as context; only the two sdmTMB fits are
  # matched comparisons of covariance structures within one fitting engine.
  variant_levels <- c("original_lmer", "site_year", "spatial_field")
  focal$model_variant <- factor(focal$model_variant, levels = variant_levels)
  p <- ggplot2::ggplot(focal, ggplot2::aes(x = term, y = estimate, colour = model_variant)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
    ggplot2::geom_pointrange(ggplot2::aes(ymin = lower_95, ymax = upper_95),
                            position = ggplot2::position_dodge(width = 0.55)) +
    ggplot2::coord_flip() + ggplot2::theme_bw() +
    ggplot2::labs(x = NULL, y = "Estimate and normal Wald 95% CI", colour = "Model")
  ggplot2::ggsave(file.path(out_dir, "coefficient_comparison.png"), p,
                  width = 12, height = 6, dpi = 250)
  plot_tests <- tests |>
    dplyr::filter(status == "OK") |>
    dplyr::mutate(model_variant = factor(model_variant, levels = variant_levels),
      distance_band = factor(distance_band, levels = cfg$distance_bands$distance_band))
  if (nrow(plot_tests)) {
    p <- ggplot2::ggplot(plot_tests, ggplot2::aes(model_variant, excess_I, fill = model_variant)) +
      ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
      ggplot2::geom_boxplot(outlier.size = 0.5) +
      ggplot2::facet_grid(residual_type ~ distance_band) + ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 35, hjust = 1)) +
      ggplot2::labs(x = NULL, y = "Moran I minus its null expectation",
        caption = "Distributions across valid BMS/year groups; boxplots are not confidence intervals.")
    ggplot2::ggsave(file.path(out_dir, "moran_comparison.png"), p, width = 13, height = 7, dpi = 250)
  }
  readme <- data.frame(topic = c("Comparison", "Convergence", "Fixed effects", "Field",
    "Residuals", "Moran", "FDR", "Limits", "Mesh", "Next decision"),
    details = c(
      "Same data/formula in sdmTMB; ML; site/year versus site/year plus annual IID field.",
      "fit_ok checks optimizer, Hessian, gradient and fixed SEs; inspect additional sanity checks too.",
      "Normal Wald 95% CIs; no missing tidy p-values interpreted as non-significance; these are not LRTs.",
      "Isotropic Matern; one range/variance shared across years and species. Range units: km.",
      paste0("Response residuals plus mle-mvn with prespecified seed ", cfg$residual_seed, "."),
      "Exploratory one-sided positive-autocorrelation screens; species/site means and duplicate-coordinate means.",
      "BH over valid network/year tests within each variant/type/band; second column pools all distance bands.",
      "Unequal averaging variances and species-specific structure limit aggregate permutation-test interpretation.",
      "Pilot cutoff is not validated; examine mesh and repeat with finer resolution if needed.",
      "Evaluate residual structure, effect sizes/CIs and parameter uncertainty together before extending models."
    ))
  writexl::write_xlsx(c(list(README = readme, model_fit = model_fit,
    focal_coefficients = as.data.frame(focal), coefficient_changes = changes,
    all_coefficients = coefficients, spatial_parameters = dplyr::bind_rows(random_parameters)),
    summaries, list(moran_all_tests = tests, validation = validation,
                    distance_bands = cfg$distance_bands, package_versions = package_versions)),
    file.path(out_dir, "comparison.xlsx"))
  message("Completed. Inspect comparison.xlsx: model_fit, coefficient_changes, paired_summary.\n", out_dir)
  invisible(out_dir)
}

# Set options(pheno.pilot.functions_only = TRUE) to load helpers without running.
if (!isTRUE(getOption("pheno.pilot.functions_only", FALSE))) run_onset_pilot(pilot_config)
