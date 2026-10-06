# phenoIMPACT -- weekend mesh sensitivity queue, 2026-09-25
# Default: four additional primary full models at 30 km (coarser than 15 km).
# The 15-km primary workflow and its saved fits are not modified.
# Run inside the existing R project. No package updates are performed.
#
# 1. Leave the active fitting worker running. If watch_spatial/espera_model is
#    occupying the console, press Esc once. The queue guard may report that the
#    current model is still active; this does NOT kill that worker.
# 2. In the SAME RStudio session run: source(file.choose()) and select this file.
# 3. This script waits for the active worker and its exports, then runs onset,
#    first_peak, offset_univoltine, and offset_multivoltine, each in a fresh R
#    subprocess. Keep RStudio open and the computer awake throughout.
#
# Output: output/phenology_plasticity/spatial_models_mesh_sensitivity/
# Each cutoff/analysis has its own manifests, checkpoints and exports.
# queue_status.csv tracks all four jobs, including crashes and stage failures.
# The initial ETA reuses the 15-km pilot only as a rough reference; completion
# before Monday is not guaranteed. Mesh cutoff is NOT a correlation range.
#
# Scientific specification inherited from pheno_plasticity_spatial.R version 4:
# same verified observations, climatic windows, scales and fixed effects;
# Gaussian ML, site/species intercepts, independent species anomaly slope,
# site-year intercept, IID annual Matern field; spatial=off, reml=FALSE.
# Only mesh construction cutoff changes. No additional LRTs or window selection.
# Coarser-mesh sensitivity alone does not establish numerical convergence as
# mesh resolution increases; a finer mesh would answer that separate question.
#
# Advanced, optional: to define functions without starting the queue:
# options(pheno.mesh.sensitivity.functions_only = TRUE)
# source(file.choose())
# pheno_mesh_sensitivity$run(cutoffs_km = 30)
#
# All workflow functions below live in a private environment, so loading this
# file does not replace functions used by an already-running primary fit.

pheno_mesh_sensitivity <- local({
spatial_config <- list(
  input_rds = NULL,       # NULL = default project path above
  output_root = NULL,     # NULL = output/phenology_plasticity/spatial_models
  pilot_root = NULL,      # NULL = output/phenology_plasticity/diagnostics/spatial_field_pilot
  analysis_groups = "primary",
  analyses = NULL,        # or c("onset", "first_peak", "offset_univoltine",
                          #      "offset_multivoltine")
  stages = "fit",
  reference_fit_hours = 707.3215 / 60, # observed onset spatial runtime; fallback ETA only
  mesh_cutoff_km = 15,
  check_only = FALSE,     # validates every selected input and builds meshes
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
  # readRDS() restores a cached mesh but does not register fmesher's plot methods.
  # make_mesh() loads fmesher on a fresh run; resume must load it explicitly too.
  # Keep the returned pkgs vector unchanged: the onset pilot fingerprint uses it.
  dependencies <- unique(c(pkgs, "fmesher"))
  missing <- dependencies[!vapply(dependencies, requireNamespace, logical(1), quietly = TRUE)]
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
  x$p_Wald <- 2 * stats::pnorm(-abs(x$estimate / x$SE))
  x$CI_excludes_zero <- x$lower_95 > 0 | x$upper_95 < 0
  x$interval_method <- "Normal Wald 95%; not an LRT"
  tibble::as_tibble(x)
}

same_values <- function(a, b) {
  if (is.numeric(a) && is.numeric(b)) {
    isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = 1e-10))
  } else identical(as.character(a), as.character(b))
}

prepare_analysis <- function(entry, coordinates, expected_response) {
  if (is.null(entry) || !inherits(entry$model, "merMod")) {
    stop("Expected a fitted lme4 model in the selected bundle entry.")
  }
  original <- entry$model
  if (length(entry$best_window) != 1L || !entry$best_window %in% c(30, 60, 90)) {
    stop("Expected one saved climatic window: 30, 60 or 90 days.")
  }
  fixed_formula <- lme4::nobars(stats::formula(original))
  mf <- stats::model.frame(original)
  raw <- as.data.frame(entry$data)
  required <- unique(c(all.vars(stats::formula(original)), "bms_id", "YEAR"))
  if (!all(required %in% names(raw))) stop("Saved data lack required columns.")
  if (!identical(all.vars(fixed_formula[[2]]), expected_response)) {
    stop("The saved model response does not match this analysis.")
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
  numeric_cols <- vapply(d, is.numeric, logical(1))
  if (any(vapply(d[numeric_cols], function(x) any(!is.finite(x)), logical(1)))) {
    stop("Non-finite response/predictor values in original fitted rows.")
  }
  if (anyNA(d) || any(!nzchar(d$SITE_ID)) || any(!nzchar(d$bms_id))) {
    stop("Missing data or empty IDs in fitted rows; do not silently drop rows.")
  }
  if (anyDuplicated(d[, c("SITE_ID", "SPECIES", "YEAR")])) {
    stop("Duplicate species/site/year rows: inspect the input before fitting.")
  }

  coords <- as.data.frame(coordinates)
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
    expected_response, " ~ ", rhs, " + (1 | SITE_ID) + (1 | SPECIES) + (0 + ",
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
  # Retain unused climatic windows separately: their missing values must not
  # remove observations from the primary fit.
  extra_names <- grep("^(clim_|photo).*tw[0-9]+$", names(raw), value = TRUE)
  extras <- raw[index, setdiff(extra_names, names(d)), drop = FALSE]
  validate_fixed_design(model_formula, d)
  list(data = as.data.frame(d), formula = model_formula,
       extras = extras, response = expected_response,
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

plot_residual_checks <- function(d, vectors, variant, response) {
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
  fitted <- d[[response]] - vectors$response
  i <- unique(round(seq(1, nrow(d), length.out = min(nrow(d), 30000L))))
  graphics::plot(fitted[i], vectors$response[i], pch = 16, cex = 0.3,
                 col = grDevices::adjustcolor("black", alpha.f = 0.15),
                 main = paste(variant, "response residuals"),
                 xlab = paste("Fitted", response, "(day of year)"), ylab = "Observed - fitted (days)")
  graphics::abline(h = 0, col = "red")
  graphics::plot(d$YEAR[i], vectors[["mle-mvn"]][i], pch = 16, cex = 0.3,
                 col = grDevices::adjustcolor("black", alpha.f = 0.15),
                 main = paste(variant, "mle-mvn by year"),
                 xlab = "Year", ylab = "mle-mvn residual")
  graphics::abline(h = 0, col = "red")
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

validate_fixed_design <- function(formula, d) {
  X <- stats::model.matrix(lme4::nobars(formula), d)
  if (nrow(X) != nrow(d) || any(!is.finite(X)) || qr(X)$rank < ncol(X)) {
    stop("Non-finite, incomplete or rank-deficient fixed-effect design.")
  }
  invisible(colnames(X))
}

moderators <- function(w) {
  c(photoperiod = paste0("photo_tw", w),
    background = paste0("clim_background_tw", w),
    predictability = paste0("clim_predictability_tw", w),
    trend = paste0("clim_trend_tw", w))
}

interaction_term <- function(terms, a, b) {
  hit <- intersect(c(paste(a, b, sep = ":"), paste(b, a, sep = ":")), terms)
  if (length(hit) != 1L) stop("Cannot identify interaction: ", a, " x ", b)
  hit
}

make_formula <- function(response, w, mods = moderators(w)) {
  anomaly <- paste0("clim_anomaly_tw", w)
  fixed <- if (length(mods)) paste0(anomaly, " * (", paste(mods, collapse = " + "), ")") else anomaly
  stats::as.formula(paste0(response, " ~ ", fixed,
    " + (1 | SITE_ID) + (1 | SPECIES) + (0 + ", anomaly,
    " | SPECIES_slope) + (1 | site_year_id)"), env = baseenv())
}

fixed_covariance <- function(fit, terms) {
  raw <- stats::vcov(fit)
  # Gaussian sdmTMB versions may expose a matrix or a list of model matrices.
  candidates <- if (is.matrix(raw) || inherits(raw, "Matrix")) list(raw) else raw
  for (v in candidates) {
    if ((is.matrix(v) || inherits(v, "Matrix")) &&
        all(terms %in% rownames(v)) && all(terms %in% colnames(v))) {
      V <- as.matrix(v)[terms, terms, drop = FALSE]
      if (all(is.finite(V))) return(V)
    }
  }
  stop("Cannot align vcov(fit) with fixed-effect coefficient names.")
}

fit_statistics <- function(fit) {
  ll <- stats::logLik(fit)
  data.frame(n_observations = nrow(fit$data), logLik = as.numeric(ll),
    n_parameters = attr(ll, "df"), AIC = stats::AIC(fit), BIC = stats::BIC(fit))
}

fit_fingerprint <- function(d, formula, mesh, run_id, cfg) {
  digest::digest(list(version = "spatial_fit_1.0", run_id = run_id,
    formula = formula_text(formula), data = d, mesh = mesh$mesh,
    spatial = "off", spatiotemporal = "iid", time = "YEAR", reml = FALSE,
    seed = cfg$fit_seed), algo = "sha256")
}

reuse_onset_pilot <- function(p, mesh, run_id, out_dir, cfg) {
  # The primary 15-km pilot cannot be reused for another triangulation.
  if (cfg$mesh_cutoff_km != 15) return(invisible(FALSE))
  destination <- file.path(out_dir, "fit_spatial_field.rds")
  if (cfg$force_refit || file.exists(destination)) return(invisible(FALSE))
  pkgs <- check_packages()
  versions <- data.frame(package = pkgs,
    version = vapply(pkgs, function(x) as.character(utils::packageVersion(x)), character(1)))
  # Reconstruct the exact signature used by onset_matern_pilot_1.0.
  signature <- list(code_version = "onset_matern_pilot_1.0", data = p$data,
    formula = paste(deparse(p$formula), collapse = " "),
    cutoff = cfg$mesh_cutoff_km, seed = cfg$fit_seed,
    gradient_threshold = cfg$gradient_threshold,
    extra_optimization_attempts = cfg$extra_optimization_attempts,
    packages = versions[versions$package %in% c("sdmTMB", "lme4"), ],
    TMB_version = as.character(utils::packageVersion("TMB")),
    fmesher_version = as.character(utils::packageVersion("fmesher")))
  pilot_id <- digest::digest(signature, algo = "sha256")
  root <- cfg$pilot_root
  if (is.null(root)) root <- here::here("output", "phenology_plasticity", "diagnostics", "spatial_field_pilot")
  directory <- file.path(root, paste0("onset_mesh15km_", substr(pilot_id, 1, 10)))
  path <- file.path(directory, "fit_spatial_field.rds")
  metadata_path <- file.path(directory, "metadata.rds")
  if (!file.exists(path) || !file.exists(metadata_path)) return(invisible(FALSE))
  meta <- readRDS(metadata_path)
  if (!identical(meta$run_id, pilot_id) ||
      !isTRUE(all.equal(meta$mesh$mesh, mesh$mesh, check.attributes = TRUE))) {
    message("Pilot triangulation differs; fitting a new onset model.")
    return(invisible(FALSE))
  }
  old <- readRDS(path)
  if (!identical(old$run_id, pilot_id) || !identical(old$variant, "spatial_field") ||
      !inherits(old$fit, "sdmTMB")) stop("Unexpected pilot checkpoint contents: ", path)
  old$checks <- fit_checks(old$fit, cfg)
  if (!old$checks$fit_ok || !all(names(p$data) %in% names(old$fit$data)) ||
      !isTRUE(all.equal(old$fit$data[, names(p$data), drop = FALSE], p$data,
                       check.attributes = TRUE))) {
    message("Pilot fit/data verification failed; fitting a new onset model.")
    return(invisible(FALSE))
  }
  old$fingerprint <- fit_fingerprint(p$data, p$formula, mesh, run_id, cfg)
  old$run_id <- run_id
  old$label <- "spatial_field"
  old$reused_from <- normalizePath(path, winslash = "/")
  atomic_save(old, destination)
  message("Reused the verified onset pilot: ", path)
  invisible(TRUE)
}

fit_model <- function(label, d, formula, mesh, run_id, out_dir, cfg) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(out_dir, paste0("fit_", label, ".rds"))
  progress_fit(path, "checking")
  progress_finished <- FALSE
  on.exit(if (!progress_finished) progress_fit(path, "failed"), add = TRUE)
  validate_fixed_design(formula, d)
  fingerprint <- fit_fingerprint(d, formula, mesh, run_id, cfg)
  reused <- file.exists(path) && !cfg$force_refit
  if (file.exists(path) && !cfg$force_refit) {
    saved <- readRDS(path)
    if (!identical(saved$fingerprint, fingerprint)) {
      stop("Checkpoint mismatch: ", path, ". Do not mix changed models/data.")
    }
    message("Reusing ", label)
    saved$checks <- fit_checks(saved$fit, cfg)
  } else {
    messages <- character()
    capture_warning <- function(w) {
      messages <<- c(messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
    start <- proc.time()[["elapsed"]]
    progress_fit(path, "running")
    set.seed(cfg$fit_seed)
    message("Fitting ", label, ": ", nrow(d), " rows")
    fit <- withCallingHandlers(sdmTMB::sdmTMB(
      formula = formula, data = d, mesh = mesh, time = "YEAR",
      family = stats::gaussian(link = "identity"), spatial = "off",
      spatiotemporal = "iid", reml = FALSE, silent = FALSE
    ), warning = capture_warning)
    checks <- fit_checks(fit, cfg)
    saved <- list(fingerprint = fingerprint, run_id = run_id, label = label,
      fit = fit, checks = checks, warnings = unique(messages),
      extra_optimization = 0L, elapsed_minutes = (proc.time()[["elapsed"]] - start) / 60)
    # Save the initial fit BEFORE any additional optimization.
    atomic_save(saved, path)
    for (attempt in seq_len(cfg$extra_optimization_attempts)) {
      if (checks$fit_ok) break
      progress_phase(paste("Extra optimization", attempt, "for", label))
      improved <- tryCatch(withCallingHandlers(sdmTMB::run_extra_optimization(
        fit, nlminb_loops = 1, newton_loops = 0), warning = capture_warning),
        error = function(e) {
          messages <<- c(messages, paste("Extra optimization failed:", conditionMessage(e)))
          NULL
        })
      if (is.null(improved)) break
      fit <- improved
      checks <- fit_checks(fit, cfg)
      saved$fit <- fit
      saved$checks <- checks
      saved$extra_optimization <- attempt
      saved$elapsed_minutes <- (proc.time()[["elapsed"]] - start) / 60
      saved$warnings <- unique(messages)
      atomic_save(saved, path)
    }
    saved$warnings <- unique(messages)
    atomic_save(saved, path)
  }
  # Check data alignment even for a resumed fit before inference or diagnostics.
  if (!identical(as.integer(saved$fit$data$source_row), as.integer(d$source_row))) {
    stop("Fit/input row order differs: ", path)
  }
  check <- cbind(saved$checks, fit_statistics(saved$fit),
    elapsed_minutes = saved$elapsed_minutes, extra_optimization = saved$extra_optimization)
  write_table(check, file.path(out_dir, paste0("fit_checks_", label, ".csv")))
  writeLines(saved$warnings, file.path(out_dir, paste0("warnings_", label, ".txt")))
  if (!saved$checks$fit_ok) stop("Numerical checks failed; checkpoint retained: ", path)
  writeLines(capture.output(sdmTMB::sanity(saved$fit,
    gradient_thresh = cfg$gradient_threshold)), file.path(out_dir, paste0("sanity_", label, ".txt")))
  progress_fit(path, if (reused) "cached" else "done", saved$elapsed_minutes * 60)
  progress_finished <- TRUE
  saved
}

build_mesh <- function(d, cfg) {
  xy <- dplyr::distinct(d, x_km, y_km)
  set.seed(cfg$fit_seed)
  tri <- sdmTMB::make_mesh(xy, xy_cols = c("x_km", "y_km"),
    cutoff = cfg$mesh_cutoff_km, type = "cutoff")
  # Project ALL observations onto the triangulation, in their fitted order.
  mesh <- project_mesh(d, tri)
  list(triangulation = tri, mesh = mesh)
}

project_mesh <- function(d, tri) {
  mesh <- sdmTMB::make_mesh(d, xy_cols = c("x_km", "y_km"), mesh = tri$mesh)
  if (nrow(mesh$loc_xy) != nrow(d) ||
      !isTRUE(all.equal(unname(as.matrix(mesh$loc_xy)),
        unname(as.matrix(d[, c("x_km", "y_km")])), check.attributes = FALSE))) {
    stop("Mesh/data projection is misaligned.")
  }
  mesh
}

write_primary_results <- function(saved, p, out_dir, mesh_cutoff_km) {
  fixed <- normal_ci(sdmTMB::tidy(saved$fit, effects = "fixed", conf.int = FALSE))
  if (!setequal(fixed$term, p$original_coefficients$term)) {
    stop("Spatial fixed-effect terms differ from the original model.")
  }
  mods <- moderators(p$best_window)
  interactions <- fixed[match(vapply(mods, function(m)
    interaction_term(fixed$term, p$anomaly, m), character(1)), fixed$term), ]
  interactions$moderator <- names(mods)
  # Preserve completed LRTs on an identical rerun of the fit stage.
  lrt_path <- file.path(out_dir, "interaction_LRT.csv")
  if (file.exists(lrt_path)) interactions <- dplyr::left_join(interactions,
    utils::read.csv(lrt_path), by = "term") else {
    interactions$p_LRT <- NA_real_
    interactions$LRT_status <- "not_run"
  }
  interactions$LRT_status[is.na(interactions$LRT_status)] <- "not_run"
  random <- sdmTMB::tidy(saved$fit, effects = "ran_pars", conf.int = TRUE)
  comparison <- dplyr::bind_rows(
    dplyr::mutate(p$original_coefficients, model = "original_lmer"),
    dplyr::mutate(fixed, model = "annual_spatial"))
  focal <- comparison[comparison$term %in% c(p$anomaly, interactions$term), ]
  changes <- dplyr::inner_join(p$original_coefficients, fixed, by = "term",
    suffix = c("_original", "_spatial")) |>
    dplyr::mutate(estimate_change = estimate_spatial - estimate_original,
                  SE_ratio = SE_spatial / SE_original)
  write_table(fixed, file.path(out_dir, "fixed_effects.csv"))
  write_table(random, file.path(out_dir, "spatial_parameters.csv"))
  write_table(interactions, file.path(out_dir, "main_interactions.csv"))
  write_table(changes, file.path(out_dir, "coefficient_changes.csv"))
  graph <- ggplot2::ggplot(focal, ggplot2::aes(term, estimate, colour = model)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
    ggplot2::geom_pointrange(ggplot2::aes(ymin = lower_95, ymax = upper_95),
      position = ggplot2::position_dodge(width = 0.55)) +
    ggplot2::coord_flip() + ggplot2::theme_bw() +
    ggplot2::labs(x = NULL, y = "Estimate and normal Wald 95% CI", colour = "Model")
  ggplot2::ggsave(file.path(out_dir, "coefficient_comparison.png"), graph,
    width = 12, height = 6, dpi = 250)
  # Slope = b_anomaly + b_interaction * moderator, other moderators at zero.
  b <- stats::setNames(fixed$estimate, fixed$term)
  V <- fixed_covariance(saved$fit, names(b))
  if (!isTRUE(all.equal(unname(sqrt(diag(V))), fixed$SE, tolerance = 1e-5))) {
    stop("Fixed-effect covariance diagonal and reported standard errors disagree.")
  }
  curves <- lapply(seq_along(mods), function(i) {
    mod <- mods[[i]]
    term <- interaction_term(names(b), p$anomaly, mod)
    # Central observed range avoids imposing +/-2 on unstandardized predictors.
    limits <- stats::quantile(p$data[[mod]], c(0.025, 0.975), names = FALSE)
    x <- seq(limits[1], limits[2], length.out = 150L)
    slope <- b[[p$anomaly]] + b[[term]] * x
    variance <- V[p$anomaly, p$anomaly] + x^2 * V[term, term] +
      2 * x * V[p$anomaly, term]
    if (any(variance < -1e-8)) stop("Negative slope variance.")
    se <- sqrt(pmax(0, variance))
    data.frame(moderator = names(mods)[i], moderator_value = x, slope = slope,
      SE = se, lower_95 = slope - stats::qnorm(.975) * se,
      upper_95 = slope + stats::qnorm(.975) * se)
  })
  curves <- dplyr::bind_rows(curves)
  write_table(curves, file.path(out_dir, "plasticity_slopes.csv"))
  graph <- ggplot2::ggplot(curves, ggplot2::aes(moderator_value, slope)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower_95, ymax = upper_95),
      fill = "#2166ac", alpha = .2) + ggplot2::geom_line(colour = "#2166ac") +
    ggplot2::facet_wrap(~moderator, scales = "free_x") + ggplot2::theme_bw() +
    ggplot2::labs(x = "Moderator (original model scale)",
      y = paste(p$response, "slope (days per anomaly model unit)"),
      caption = "Population fixed effects; other moderators = 0. Pointwise normal Wald 95% CI.")
  ggplot2::ggsave(file.path(out_dir, "plasticity_slopes.png"), graph,
    width = 11, height = 7, dpi = 250)
  readme <- data.frame(topic = c("Window", "Model", "Mesh", "Intervals", "Tests", "Slopes"),
    details = c(paste("Preserved original window:", p$best_window, "days"),
      "Gaussian ML; site/species intercepts, independent species anomaly slope, site/year intercept, annual IID Matern field.",
      paste(mesh_cutoff_km, "km construction cutoff; mesh sensitivity fit; primary reference is 15 km. Cutoff is not a correlation range or an optimized resolution."),
      "Normal Wald 95% confidence intervals; inference conditional on the selected window and model.",
      "p_Wald is a normal Wald test; p_LRT is a matched spatial reduced-model LRT, when run.",
      "Partial anomaly slopes with other moderators at zero; central observed moderator range; covariance included."))
  writexl::write_xlsx(list(README = readme, main_interactions = interactions,
    fixed_effects = fixed, coefficient_changes = changes,
    spatial_parameters = random, fit_checks = cbind(saved$checks, fit_statistics(saved$fit)),
    plasticity_slopes = curves), file.path(out_dir, "results.xlsx"))
  invisible(interactions)
}

run_interactions <- function(saved, p, mesh, run_id, out_dir, cfg) {
  full <- saved$fit
  labels <- attr(stats::terms(lme4::nobars(p$formula)), "term.labels")
  mods <- moderators(p$best_window)
  rows <- list()
  for (i in seq_along(mods)) {
    term <- interaction_term(labels, p$anomaly, mods[[i]])
    reduced_formula <- stats::update.formula(p$formula, paste(". ~ . -", term))
    reduced_labels <- attr(stats::terms(lme4::nobars(reduced_formula)), "term.labels")
    if (!setequal(reduced_labels, setdiff(labels, term))) stop("Incorrect reduced-model terms.")
    row <- data.frame(term = term, Chisq = NA_real_, LRT_df = 1L, p_LRT = NA_real_,
      delta_AIC_reduced_minus_full = NA_real_, LRT_status = "not_run")
    row <- tryCatch({
      reduced <- fit_model(paste0("drop_", names(mods)[i]), p$data,
        reduced_formula, mesh, run_id, file.path(out_dir, "interactions"), cfg)$fit
      if (!identical(full$data$source_row, reduced$data$source_row)) stop("LRT rows differ.")
      ll_full <- stats::logLik(full)
      ll_reduced <- stats::logLik(reduced)
      df <- attr(ll_full, "df") - attr(ll_reduced, "df")
      chi <- 2 * (as.numeric(ll_full) - as.numeric(ll_reduced))
      if (!is.finite(chi) || df != 1L) stop("Invalid fixed-interaction LRT degrees of freedom.")
      if (chi < -1e-3) stop("Reduced likelihood exceeds full likelihood; recheck optimization.")
      row$Chisq <- max(0, chi)
      row$LRT_df <- df
      row$p_LRT <- stats::pchisq(row$Chisq, df = df, lower.tail = FALSE)
      row$delta_AIC_reduced_minus_full <- stats::AIC(reduced) - stats::AIC(full)
      row$LRT_status <- "OK"
      rm(reduced)
      row
    }, error = function(e) {
      row$LRT_status <- paste("ERROR:", conditionMessage(e))
      warning(row$LRT_status, call. = FALSE)
      row
    })
    rows[[i]] <- row
    write_table(dplyr::bind_rows(rows), file.path(out_dir, "interaction_LRT.csv"))
    invisible(gc())
  }
  write_primary_results(saved, p, out_dir, cfg$mesh_cutoff_km)
  if (any(vapply(rows, function(r) r$LRT_status != "OK", logical(1)))) {
    stop("Some interaction LRTs failed; see interaction_LRT.csv. Successful fits are retained.")
  }
  invisible(NULL)
}

run_diagnostics <- function(saved, p, out_dir, cfg) {
  vectors <- aligned_residuals(saved$fit, p$data, cfg$residual_seed)
  atomic_save(list(seed = cfg$residual_seed, source_rows = p$data$source_row,
    residuals = vectors), file.path(out_dir, "residuals_spatial_field.rds"))
  grDevices::pdf(file.path(out_dir, "residual_checks.pdf"), width = 10, height = 8)
  device <- grDevices::dev.cur()
  on.exit(if (device %in% grDevices::dev.list()) grDevices::dev.off(device), add = TRUE)
  graphics::par(mfrow = c(2, 2))
  plot_residual_checks(p$data, vectors, "annual_spatial", p$response)
  grDevices::dev.off(device)
  aggregate <- list(aggregate_residuals(p$data, p$original_residuals, "original_lmer", "response"))
  for (type in names(vectors)) aggregate[[length(aggregate) + 1L]] <-
    aggregate_residuals(p$data, vectors[[type]], "annual_spatial", type)
  aggregate <- dplyr::bind_rows(aggregate)
  atomic_save(aggregate, file.path(out_dir, "aggregated_residuals.rds"))
  tests <- run_moran(aggregate, cfg)
  write_table(tests, file.path(out_dir, "moran_all_tests.csv"))
  summary <- tests |>
    dplyr::group_by(model_variant, residual_type, distance_band) |>
    dplyr::summarise(n_tests = dplyr::n(), n_valid = sum(status == "OK"),
      n_positive_FDR = sum(status == "OK" & excess_I > 0 & p_FDR_band < .05, na.rm = TRUE),
      median_excess_I = if (any(status == "OK")) stats::median(excess_I[status == "OK"]) else NA_real_,
      .groups = "drop")
  write_table(summary, file.path(out_dir, "moran_summary.csv"))
  writeLines(c(paste("Residual seed:", cfg$residual_seed),
    paste("Permutation seed:", cfg$permutation_seed), paste("Permutations:", cfg$n_permutations),
    "One-sided exploratory positive autocorrelation screens of aggregated residuals.",
    "BH separately by model/residual type/band; second FDR column pools bands.",
    "Unequal aggregation variances and species-specific patterns limit interpretation."),
    file.path(out_dir, "diagnostic_notes.txt"))
  ok <- tests[tests$status == "OK", ]
  if (nrow(ok)) {
    ok$distance_band <- factor(ok$distance_band, levels = cfg$distance_bands$distance_band)
    graph <- ggplot2::ggplot(ok, ggplot2::aes(model_variant, excess_I, fill = model_variant)) +
      ggplot2::geom_hline(yintercept = 0, linetype = 2) + ggplot2::geom_boxplot(outlier.size = .4) +
      ggplot2::facet_grid(residual_type ~ distance_band) + ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 35, hjust = 1)) +
      ggplot2::labs(x = NULL, y = "Moran I minus null expectation")
    ggplot2::ggsave(file.path(out_dir, "moran_comparison.png"), graph, width = 12, height = 7, dpi = 250)
  }
  if (any(grepl("^ERROR", tests$status))) warning("Some Moran screens failed; see moran_all_tests.csv.")
  invisible(NULL)
}

add_extras <- function(p) {
  d <- p$data
  for (nm in names(p$extras)) d[[nm]] <- p$extras[[nm]]
  d
}

run_windows <- function(p, tri, run_id, out_dir, cfg) {
  d <- add_extras(p)
  needed <- unlist(lapply(c(30, 60, 90), function(w) c(paste0("clim_anomaly_tw", w), moderators(w))))
  if (!all(needed %in% names(d))) stop("Saved data lack predictors for the three-window comparison.")
  if (anyNA(d[, needed]) || any(!is.finite(as.matrix(d[, needed])))) {
    stop("Three-window comparison requires the same complete fitted rows; no silent filtering.")
  }
  mesh <- project_mesh(d, tri)
  tables <- list()
  coefs <- list()
  path <- file.path(out_dir, "windows")
  for (w in c(30, 60, 90)) for (kind in c("anomaly_only", "full")) {
    formula <- make_formula(p$response, w, if (kind == "full") moderators(w) else character())
    label <- paste(kind, w, sep = "_")
    fit <- fit_model(label, d, formula, mesh, run_id, path, cfg)$fit
    tables[[label]] <- cbind(window = w, model = kind, fit_statistics(fit))
    if (kind == "full") coefs[[label]] <- dplyr::mutate(
      normal_ci(sdmTMB::tidy(fit, effects = "fixed", conf.int = FALSE)), window = w)
    write_table(dplyr::bind_rows(tables), file.path(path, "window_AIC.csv"))
    rm(fit)
    invisible(gc())
  }
  table <- dplyr::bind_rows(tables) |>
    dplyr::group_by(model) |>
    dplyr::mutate(delta_AIC = AIC - min(AIC),
      weight = exp(-.5 * delta_AIC) / sum(exp(-.5 * delta_AIC))) |>
    dplyr::ungroup()
  write_table(table, file.path(path, "window_AIC.csv"))
  write_table(dplyr::bind_rows(coefs), file.path(path, "window_fixed_effects.csv"))
  writeLines("Exploratory spatial window comparison. Primary models retain the original selected window.",
    file.path(path, "README.txt"))
}

run_alternatives <- function(p, tri, run_id, out_dir, cfg) {
  d <- add_extras(p)
  mods <- moderators(p$best_window)
  autocorr <- paste0("clim_autocorr_tw", p$best_window)
  needed <- c(p$response, p$anomaly, mods, autocorr)
  if (!all(needed %in% names(d))) stop("Saved data lack climatic autocorrelation for alternatives.")
  keep <- stats::complete.cases(d[, needed]) & apply(d[, needed], 1, function(x) all(is.finite(x)))
  d <- droplevels(d[keep, , drop = FALSE])
  if (!nrow(d)) stop("No common complete rows for alternative predictor sets.")
  mesh <- project_mesh(d, tri)
  sets <- list(main = mods, main_plus_autocorr = c(mods, autocorr),
    autocorr_instead_predictability = c(mods[names(mods) != "predictability"], autocorr),
    no_trend = mods[names(mods) != "trend"])
  path <- file.path(out_dir, "alternatives")
  rows <- list()
  coefs <- list()
  for (nm in names(sets)) {
    fit <- fit_model(nm, d, make_formula(p$response, p$best_window, sets[[nm]]),
      mesh, run_id, path, cfg)$fit
    rows[[nm]] <- cbind(model = nm, fit_statistics(fit), n_original = nrow(p$data))
    coefs[[nm]] <- dplyr::mutate(normal_ci(sdmTMB::tidy(fit,
      effects = "fixed", conf.int = FALSE)), model = nm)
    write_table(dplyr::bind_rows(rows), file.path(path, "alternative_AIC.csv"))
    rm(fit)
    invisible(gc())
  }
  table <- dplyr::bind_rows(rows)
  table$delta_AIC <- table$AIC - min(table$AIC)
  write_table(table, file.path(path, "alternative_AIC.csv"))
  write_table(dplyr::bind_rows(coefs), file.path(path, "alternative_fixed_effects.csv"))
  writeLines("All alternatives use the same complete subset of the original primary fitted rows and same triangulation. Do not compare AIC to models fitted on different rows.",
    file.path(path, "README.txt"))
}

run_spatial_plasticity_worker <- function(cfg = spatial_config) {
  pkgs <- check_packages()
  allowed_stages <- c("fit", "interactions", "diagnostics", "windows", "alternatives")
  if (!length(cfg$stages) || any(!cfg$stages %in% allowed_stages)) stop("Unknown or empty stages.")
  if (!length(cfg$analysis_groups) || any(!cfg$analysis_groups %in% c("primary", "sensitivity_1zero"))) {
    stop("Unknown or empty analysis_groups.")
  }
  if (length(cfg$mesh_cutoff_km) != 1L || !is.finite(cfg$mesh_cutoff_km) ||
      cfg$mesh_cutoff_km <= 0 || cfg$mesh_cutoff_km == 15) {
    stop("This separate sensitivity workflow requires a positive cutoff other than 15 km.")
  }
  if (length(cfg$n_permutations) != 1L || !is.finite(cfg$n_permutations) ||
      cfg$n_permutations < 99 || cfg$n_permutations != as.integer(cfg$n_permutations)) {
    stop("n_permutations must be an integer >=99.")
  }
  input <- cfg$input_rds
  if (is.null(input)) input <- here::here("output", "phenology_plasticity", "plasticity_main_models.rds")
  if (!file.exists(input)) stop("Input not found: ", input, ". Set spatial_config$input_rds.")
  root <- cfg$output_root
  if (is.null(root)) root <- here::here("output", "phenology_plasticity", "spatial_models")
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  versions <- data.frame(package = unique(c(pkgs, "TMB", "fmesher", "Matrix")), stringsAsFactors = FALSE)
  versions$version <- vapply(versions$package, function(x) as.character(utils::packageVersion(x)), character(1))
  writeLines(capture.output(sessionInfo()), file.path(root, "sessionInfo.txt"))
  message("Reading original bundle. Large models may require several GB of memory.")
  progress_phase("Reading the original model bundle")
  bundle <- readRDS(input)
  responses <- c(onset = "ONSET_mean", first_peak = "FIRST_PEAK",
    offset_univoltine = "OFFSET_mean", offset_multivoltine = "OFFSET_mean")
  jobs <- list()
  manifest <- list()
  for (group in unique(cfg$analysis_groups)) {
    entries <- bundle[[group]]
    if (is.null(entries)) stop("Bundle lacks group: ", group)
    selected <- if (is.null(cfg$analyses)) names(entries) else intersect(cfg$analyses, names(entries))
    for (analysis in selected) {
      if (!analysis %in% names(responses)) stop("Unrecognized analysis: ", analysis)
      key <- paste(group, analysis, sep = "__")
      message("Validating ", key)
      progress_phase(paste("Preparing data and mesh:", key))
      p <- prepare_analysis(entries[[analysis]], bundle$site_coordinates, responses[[analysis]])
      expected_terms <- attr(stats::terms(lme4::nobars(make_formula(p$response, p$best_window))), "term.labels")
      actual_terms <- attr(stats::terms(lme4::nobars(p$formula)), "term.labels")
      if (!setequal(expected_terms, actual_terms)) stop("Unexpected fixed formula for ", key)
      signature <- list(version = "pheno_spatial_1.0", data = p$data,
        formula = formula_text(p$formula), cutoff = cfg$mesh_cutoff_km,
        seed = cfg$fit_seed, gradient_threshold = cfg$gradient_threshold,
        extra_optimization_attempts = cfg$extra_optimization_attempts,
        packages = versions[versions$package %in% c("sdmTMB", "lme4", "TMB", "fmesher", "Matrix"), ])
      run_id <- digest::digest(signature, algo = "sha256")
      rm(signature)
      out_dir <- file.path(root, paste0(key, "_", substr(run_id, 1, 12)))
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      metadata_path <- file.path(out_dir, "metadata.rds")
      if (file.exists(metadata_path)) {
        old <- readRDS(metadata_path)
        if (!identical(old$run_id, run_id)) stop("Metadata hash mismatch.")
        tri <- old$triangulation
      } else {
        built <- build_mesh(p$data, cfg)
        tri <- built$triangulation
        rm(built)
      }
      validation <- data.frame(analysis_group = group, analysis = analysis,
        response = p$response, selected_window_days = p$best_window,
        n_observations = nrow(p$data), n_species = nlevels(p$data$SPECIES),
        n_sites = nlevels(p$data$SITE_ID), n_site_years = nlevels(p$data$site_year_id),
        n_mesh_vertices = tri$mesh$n, cutoff_km = cfg$mesh_cutoff_km,
        fitted_rows_and_design_verified = TRUE, window_selection = "preserved_original")
      write_table(validation, file.path(out_dir, "validation.csv"))
      message("Writing mesh figure: ", key)
      progress_phase(paste("Writing mesh figure:", key))
      grDevices::png(file.path(out_dir, "mesh.png"), width = 1800, height = 1500, res = 180)
      tryCatch(graphics::plot(tri), finally = grDevices::dev.off())
      atomic_save(list(run_id = run_id, input = normalizePath(input, winslash = "/"),
        bundle_created = bundle$created_at, settings = cfg, validation = validation,
        formula = p$formula, original_formula = p$original_formula,
        package_versions = versions, triangulation = tri, source_rows = p$data$source_row), metadata_path)
      # Keep compact prepared inputs in files; do not hold all fits in memory.
      atomic_save(p, file.path(out_dir, "prepared_input.rds"))
      jobs[[key]] <- list(run_id = run_id, out_dir = out_dir)
      manifest[[key]] <- cbind(validation, run_id = run_id, output_directory = out_dir,
        checkpoint = file.path(out_dir, "fit_spatial_field.rds"))
      rm(p, tri)
      invisible(gc())
    }
  }
  rm(bundle, entries)
  invisible(gc())
  if (!length(jobs)) stop("No analyses selected.")
  if (!is.null(cfg$analyses)) {
    found <- unique(dplyr::bind_rows(manifest)$analysis)
    if (!all(cfg$analyses %in% found)) stop("Requested analyses missing: ", paste(setdiff(cfg$analyses, found), collapse = ", "))
  }
  manifest <- dplyr::bind_rows(manifest)
  write_table(manifest, file.path(root, "model_manifest.csv"))
  if (cfg$check_only) {
    message("Input/design/mesh checks complete for ", length(jobs), " analyses. No models fitted.\n", root)
    return(invisible(manifest))
  }
  stages <- allowed_stages[allowed_stages %in% cfg$stages]
  if (!is.null(.pheno_progress$state)) {
    .pheno_progress$state$models <- progress_plan(manifest, stages, cfg)
    .pheno_progress$state$stages_total <- length(stages) * length(jobs)
    progress_save()
  }
  status <- list()
  for (stage in stages) for (key in names(jobs)) {
    job <- jobs[[key]]
    out_dir <- job$out_dir
    message("[", stage, "] ", key)
    progress_phase(paste(stage, key))
    result <- tryCatch({
      p <- readRDS(file.path(out_dir, "prepared_input.rds"))
      meta <- readRDS(file.path(out_dir, "metadata.rds"))
      mesh <- project_mesh(p$data, meta$triangulation)
      if (stage %in% c("fit", "interactions", "diagnostics")) {
        checkpoint <- file.path(out_dir, "fit_spatial_field.rds")
        if (stage == "fit" && key == "primary__onset") {
          reuse_onset_pilot(p, mesh, job$run_id, out_dir, cfg)
        }
        if (stage != "fit" && !file.exists(checkpoint)) stop("Run stage 'fit' first for this analysis.")
        main_cfg <- cfg
        # A forced refit is performed once in the fit stage, not again per stage.
        if (stage != "fit") main_cfg$force_refit <- FALSE
        saved <- fit_model("spatial_field", p$data, p$formula, mesh, job$run_id, out_dir, main_cfg)
        if (stage == "fit") {
          progress_phase(paste("Writing tables and figures:", key))
          # Previous LRT numbers cannot accompany a newly forced full fit.
          if (cfg$force_refit) unlink(file.path(out_dir, "interaction_LRT.csv"))
          write_primary_results(saved, p, out_dir, cfg$mesh_cutoff_km)
        }
        if (stage == "interactions") run_interactions(saved, p, mesh, job$run_id, out_dir, cfg)
        if (stage == "diagnostics") run_diagnostics(saved, p, out_dir, cfg)
      }
      if (stage == "windows") run_windows(p, meta$triangulation, job$run_id, out_dir, cfg)
      if (stage == "alternatives") run_alternatives(p, meta$triangulation, job$run_id, out_dir, cfg)
      "OK"
    }, error = function(e) {
      message("Failed ", stage, " / ", key, ": ", conditionMessage(e))
      paste("ERROR:", conditionMessage(e))
    })
    status[[paste(key, stage)]] <- data.frame(analysis = key, stage = stage,
      status = result, timestamp = as.character(Sys.time()), output_directory = out_dir)
    write_table(dplyr::bind_rows(status), file.path(root, "workflow_status.csv"))
    progress_stage_end(key, stage, result)
    if (exists("saved", inherits = FALSE)) rm(saved)
    if (exists("p", inherits = FALSE)) rm(p)
    if (exists("meta", inherits = FALSE)) rm(meta)
    if (exists("mesh", inherits = FALSE)) rm(mesh)
    invisible(gc())
  }
  all_interactions <- lapply(seq_len(nrow(manifest)), function(i) {
    path <- file.path(manifest$output_directory[i], "main_interactions.csv")
    if (!file.exists(path)) return(NULL)
    cbind(analysis_group = manifest$analysis_group[i], analysis = manifest$analysis[i],
      window = manifest$selected_window_days[i], utils::read.csv(path))
  })
  all_interactions <- dplyr::bind_rows(all_interactions)
  if (nrow(all_interactions)) write_table(all_interactions, file.path(root, "all_main_interactions.csv"))
  status <- dplyr::bind_rows(status)
  if (any(status$status != "OK")) warning("Some stages failed. Inspect workflow_status.csv; rerun resumes completed fits.", call. = FALSE)
  message("Workflow finished. Results and status: ", root)
  invisible(list(manifest = manifest, status = status, interactions = all_interactions))
}

# Live monitoring runs in the parent R session; model fitting runs in one child
# R session. A timer in the fitting session itself would freeze inside TMB.
.pheno_worker_functions <- c(
  "check_packages", "write_table", "normal_ci", "same_values", "prepare_analysis",
  "fit_checks", "aligned_residuals", "aggregate_residuals", "bh_valid", "moran_one",
  "run_moran", "plot_residual_checks", "formula_text", "atomic_save",
  "validate_fixed_design", "moderators", "interaction_term", "make_formula",
  "fixed_covariance", "fit_statistics", "fit_fingerprint", "reuse_onset_pilot",
  "fit_model", "build_mesh", "project_mesh", "write_primary_results",
  "run_interactions", "run_diagnostics", "add_extras", "run_windows",
  "run_alternatives", "run_spatial_plasticity_worker", "progress_save",
  "progress_phase", "progress_plan", "progress_fit", "progress_stage_end"
)

write_worker_snapshot <- function(path, cfg, scope = environment(run_spatial_plasticity)) {
  missing <- .pheno_worker_functions[!vapply(.pheno_worker_functions, function(nm)
    exists(nm, envir = scope, mode = "function", inherits = TRUE), logical(1))]
  if (length(missing)) stop("Run the whole script first. Missing functions: ", paste(missing, collapse = ", "))
  definitions <- new.env(parent = baseenv())
  definitions$spatial_config <- cfg
  for (nm in .pheno_worker_functions) definitions[[nm]] <- get(nm, envir = scope, inherits = TRUE)
  # Export only this workflow's functions and settings. No filename discovery,
  # editor API, or serialization of unrelated objects in the user's workspace.
  dump(c("spatial_config", .pheno_worker_functions), file = path, envir = definitions)
  cat("\n.pheno_progress <- new.env(parent = emptyenv())\n.pheno_progress$state <- NULL\n",
      file = path, append = TRUE)
  invisible(path)
}

.pheno_progress <- new.env(parent = emptyenv())
.pheno_progress$state <- NULL

progress_save <- function() {
  s <- .pheno_progress$state
  if (!is.null(s)) {
    s$updated_at <- as.numeric(Sys.time())
    .pheno_progress$state <- s
    atomic_save(s, .pheno_progress$path)
  }
  invisible(NULL)
}

progress_phase <- function(text) {
  if (!is.null(.pheno_progress$state)) {
    .pheno_progress$state$phase <- text
    progress_save()
  }
}

progress_plan <- function(manifest, stages, cfg) {
  rows <- list()
  for (stage in stages) for (i in seq_len(nrow(manifest))) {
    directory <- manifest$output_directory[i]
    labels <- switch(stage,
      fit = "spatial_field",
      interactions = paste0("drop_", names(moderators(manifest$selected_window_days[i]))),
      windows = as.vector(rbind(paste0("anomaly_only_", c(30, 60, 90)),
                               paste0("full_", c(30, 60, 90)))),
      alternatives = c("main", "main_plus_autocorr", "autocorr_instead_predictability", "no_trend"),
      character())
    if (stage != "fit") directory <- file.path(directory, stage)
    for (label in labels) {
      reference <- NA_real_
      old <- file.path(directory, paste0("fit_checks_", label, ".csv"))
      if (file.exists(old)) reference <- tryCatch({
        tab <- utils::read.csv(old)
        as.numeric(tab$elapsed_minutes[1]) * 60
      }, error = function(e) NA_real_)
      rows[[length(rows) + 1L]] <- data.frame(
        id = file.path(directory, paste0("fit_", label, ".rds")),
        analysis = paste(manifest$analysis_group[i], manifest$analysis[i], sep = "__"),
        stage = stage, model = label, status = "pending", started_at = NA_real_,
        ended_at = NA_real_, reference_seconds = reference, duration_seconds = NA_real_,
        stringsAsFactors = FALSE)
    }
  }
  if (length(rows)) do.call(rbind, rows) else data.frame(
    id = character(), analysis = character(), stage = character(), model = character(),
    status = character(), started_at = numeric(), ended_at = numeric(),
    reference_seconds = numeric(), duration_seconds = numeric())
}

progress_fit <- function(path, status, duration_seconds = NA_real_) {
  s <- .pheno_progress$state
  if (is.null(s) || is.null(s$models)) return(invisible(NULL))
  i <- match(path, s$models$id)
  if (is.na(i)) return(invisible(NULL)) # Reading a full fit for later stages.
  # A full fit reused for LRT/diagnostics must not be counted a second time.
  if (s$models$status[i] %in% c("done", "cached", "failed", "skipped")) return(invisible(NULL))
  s$models$status[i] <- status
  if (status == "checking") s$models$started_at[i] <- as.numeric(Sys.time())
  if (status == "running") s$models$started_at[i] <- as.numeric(Sys.time())
  if (status %in% c("done", "cached", "failed")) {
    s$models$ended_at[i] <- as.numeric(Sys.time())
    s$models$duration_seconds[i] <- duration_seconds
  }
  .pheno_progress$state <- s
  progress_save()
}

progress_stage_end <- function(key, stage, result) {
  s <- .pheno_progress$state
  if (is.null(s)) return(invisible(NULL))
  pending <- s$models$analysis == key & s$models$stage == stage &
    s$models$status %in% c("pending", "checking", "running")
  s$models$status[pending] <- "skipped"
  s$stages_finished <- s$stages_finished + 1L
  if (result != "OK") s$stage_errors <- s$stage_errors + 1L
  .pheno_progress$state <- s
  progress_save()
}

format_elapsed <- function(seconds) {
  if (length(seconds) != 1L || !is.finite(seconds)) return("unknown")
  minutes <- floor(max(0, seconds) / 60)
  sprintf("%dh %02dm", minutes %/% 60, minutes %% 60)
}

progress_summary <- function(s, now = as.numeric(Sys.time())) {
  m <- s$models
  n <- if (is.null(m)) 0L else nrow(m)
  done <- if (n) sum(m$status %in% c("done", "cached")) else 0L
  failed <- if (n) sum(m$status %in% c("failed", "skipped")) else 0L
  active <- if (n) which(m$status %in% c("checking", "running")) else integer()
  remaining <- if (n) which(m$status %in% c("pending", "checking", "running")) else integer()
  history <- if (n) which(is.finite(m$duration_seconds) & m$duration_seconds > 0 &
    m$status %in% c("done", "cached")) else integer()
  predictions <- vapply(seq_len(n), function(i) {
    if (is.finite(m$reference_seconds[i]) && m$reference_seconds[i] > 0) return(m$reference_seconds[i])
    same_analysis <- history[m$analysis[history] == m$analysis[i]]
    if (length(same_analysis)) return(stats::median(m$duration_seconds[same_analysis]))
    if (length(history)) return(stats::median(m$duration_seconds[history]))
    s$fallback_seconds
  }, numeric(1))
  eta <- if (n) sum(predictions[remaining]) else NA_real_
  overdue <- FALSE
  current <- s$phase
  active_elapsed <- NA_real_
  if (length(active)) {
    i <- active[1]
    active_elapsed <- max(0, (if (s$finished) s$ended_at else now) - m$started_at[i])
    current <- paste(m$analysis[i], m$model[i], paste0("[", m$status[i], "]"))
    if (startsWith(s$phase, "Extra optimization")) current <- paste(current, s$phase, sep = " | ")
    if (m$status[i] == "running") {
      overdue <- is.finite(predictions[i]) && active_elapsed >= predictions[i]
      eta <- if (overdue) NA_real_ else eta - active_elapsed
    }
  }
  if (s$finished && !s$aborted && !failed && !s$stage_errors) eta <- 0
  if (s$aborted || failed || s$stage_errors) eta <- NA_real_
  resolved <- done + failed
  data.frame(timestamp = format(as.POSIXct(now, origin = "1970-01-01"), "%Y-%m-%d %H:%M:%S"),
    elapsed_seconds = max(0, (if (s$finished) s$ended_at else now) - s$started_at),
    current = current, current_elapsed_seconds = active_elapsed,
    fits_total = n, fits_successful = done, fits_cached = if (n) sum(m$status == "cached") else 0L,
    fits_failed_or_skipped = failed, fits_resolved = resolved,
    stages_finished = s$stages_finished, stages_total = s$stages_total,
    stage_errors = s$stage_errors, approximate_fit_seconds_remaining = eta,
    current_fit_exceeds_reference = overdue, finished = s$finished, aborted = s$aborted,
    eta_scope = "Fitting only; preparation, checkpoint verification, exports and diagnostics excluded.",
    stringsAsFactors = FALSE)
}

spatial_status <- function(job = getOption("pheno.spatial.active_job"), quiet = FALSE) {
  if (is.null(job)) stop("No running job. Start with run_spatial_plasticity().")
  s <- tryCatch(readRDS(job$progress_file), error = function(e) NULL)
  if (is.null(s)) {
    if (!quiet) message("Progress state is being updated; try again shortly.")
    return(invisible(NULL))
  }
  alive <- job$process$is_alive()
  if (!alive && !s$finished) {
    # A killed worker cannot run its R on.exit() handler.
    ended <- tryCatch(as.numeric(job$process$get_end_time()), error = function(e) NA_real_)
    s$ended_at <- if (length(ended) == 1L && is.finite(ended)) ended else as.numeric(Sys.time())
    s$finished <- TRUE
    s$aborted <- TRUE
    s$phase <- "Worker stopped before completion"
    if (!is.null(s$models)) s$models$status[s$models$status %in% c("checking", "running")] <- "failed"
    atomic_save(s, job$progress_file)
  }
  x <- progress_summary(s)
  x$worker_alive <- alive
  # Explicit progress is completed tasks, never an invented within-fit percent.
  percent <- if (x$fits_total) round(100 * x$fits_resolved / x$fits_total) else NA_real_
  eta <- if (is.finite(x$approximate_fit_seconds_remaining)) {
    if (x$approximate_fit_seconds_remaining > 0) paste0("~", ceiling(x$approximate_fit_seconds_remaining / 3600), " h") else "no fits pending"
  } else if (x$aborted || x$stage_errors || x$fits_failed_or_skipped) {
    "unavailable: inspect failed/stopped tasks"
  } else if (x$current_fit_exceeds_reference) {
    "unknown: current fit exceeds its reference duration"
  } else "unavailable during preparation or diagnostics"
  lines <- c(paste0("[", x$timestamp, "] Fits resolved: ", x$fits_resolved, "/", x$fits_total,
      if (is.finite(percent)) paste0(" (", percent, "%)") else "",
      " | cached: ", x$fits_cached, " | failed/skipped: ", x$fits_failed_or_skipped),
    paste0("Current: ", x$current, " | current fit: ",
      if (is.finite(x$current_elapsed_seconds)) format_elapsed(x$current_elapsed_seconds) else "not fitting"),
    paste0("Elapsed: ", format_elapsed(x$elapsed_seconds), " | approximate fitting time remaining: ", eta),
    paste0("Stages: ", x$stages_finished, "/", x$stages_total,
      " | stage errors: ", x$stage_errors, " | worker running: ", x$worker_alive),
    if (x$aborted) "STOPPED WITH ERROR: inspect worker.log and restart from saved fits." else x$eta_scope)
  writeLines(lines, file.path(job$monitor_directory, "progress.txt"))
  write_table(x, file.path(job$monitor_directory, "progress.csv"))
  if (!quiet) cat(paste(lines, collapse = "\n"), "\n\n")
  invisible(x)
}

watch_spatial <- function(job = getOption("pheno.spatial.active_job"),
                          interval_seconds = 60) {
  if (is.null(job)) stop("No job to monitor.")
  if (length(interval_seconds) != 1L || !is.finite(interval_seconds) || interval_seconds < 1) {
    stop("interval_seconds must be >=1.")
  }
  tryCatch({
    repeat {
      spatial_status(job)
      if (!job$process$is_alive()) break
      # This parent session stays responsive while TMB runs in the child.
      deadline <- Sys.time() + interval_seconds
      while (Sys.time() < deadline && job$process$is_alive()) Sys.sleep(1)
    }
    result <- job$process$get_result() # Surface worker errors; never report false success.
    message("Finished. See results and workflow_status.csv in: ", job$output_root)
    invisible(result)
  }, interrupt = function(e) {
    message("Display paused; FITTING CONTINUES. Use spatial_status() for one update,\n",
      "watch_spatial() to resume, or stop_spatial() to stop the calculation. Keep this R session open.")
    invisible(job)
  })
}

stop_spatial <- function(job = getOption("pheno.spatial.active_job")) {
  if (is.null(job)) return(invisible(FALSE))
  if (job$process$is_alive()) job$process$kill_tree()
  message("Worker stopped. Completed checkpoints are retained; the current unsaved fit must restart.")
  invisible(TRUE)
}

pilot_reference_seconds <- function(cfg) {
  root <- cfg$pilot_root
  files <- if (!is.null(root) && dir.exists(root)) list.files(root,
    pattern = "^fit_checks\\.csv$", recursive = TRUE, full.names = TRUE) else character()
  files <- files[grepl("onset_mesh15km_", files, fixed = TRUE)]
  if (length(files)) files <- files[order(file.info(files)$mtime, decreasing = TRUE)]
  for (path in files) {
    value <- tryCatch({
      x <- utils::read.csv(path)
      as.numeric(x$elapsed_minutes[x$model_variant == "spatial_field" & x$fit_ok][1]) * 60
    }, error = function(e) NA_real_)
    if (length(value) == 1L && is.finite(value) && value > 0) return(value)
  }
  cfg$reference_fit_hours * 3600
}

run_spatial_plasticity <- function(cfg = spatial_config, script_path = NULL) {
  # script_path is retained for compatibility; loaded functions are snapshotted.
  check_packages()
  if (!requireNamespace("callr", quietly = TRUE)) {
    stop('Live monitoring needs callr. Install it once with install.packages("callr").')
  }
  previous <- getOption("pheno.spatial.active_job")
  if (!is.null(previous) && previous$process$is_alive()) {
    stop("A spatial job is already running. Use spatial_status() or watch_spatial().")
  }
  if (!is.numeric(cfg$reference_fit_hours) || length(cfg$reference_fit_hours) != 1L ||
      !is.finite(cfg$reference_fit_hours) || cfg$reference_fit_hours <= 0) stop("Invalid reference_fit_hours.")
  if (is.null(cfg$input_rds)) cfg$input_rds <- here::here("output", "phenology_plasticity", "plasticity_main_models.rds")
  if (!file.exists(cfg$input_rds)) stop("Input not found: ", cfg$input_rds)
  cfg$input_rds <- normalizePath(cfg$input_rds, winslash = "/", mustWork = TRUE)
  if (is.null(cfg$output_root)) cfg$output_root <- here::here("output", "phenology_plasticity", "spatial_models")
  dir.create(cfg$output_root, recursive = TRUE, showWarnings = FALSE)
  cfg$output_root <- normalizePath(cfg$output_root, winslash = "/", mustWork = TRUE)
  if (is.null(cfg$pilot_root)) cfg$pilot_root <- here::here("output", "phenology_plasticity", "diagnostics", "spatial_field_pilot")
  cfg$pilot_root <- normalizePath(cfg$pilot_root, winslash = "/", mustWork = FALSE)
  directory <- tempfile(pattern = paste0("monitor_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_"),
                        tmpdir = cfg$output_root)
  dir.create(directory)
  snapshot <- file.path(directory, "analysis_script.R")
  cfg$progress_file <- file.path(directory, "progress_state.rds")
  write_worker_snapshot(snapshot, cfg)
  s <- list(started_at = as.numeric(Sys.time()), updated_at = as.numeric(Sys.time()),
    ended_at = NA_real_, phase = "Starting; reading and validating data", models = NULL,
    stages_finished = 0L, stages_total = 0L, stage_errors = 0L,
    fallback_seconds = pilot_reference_seconds(cfg), finished = FALSE, aborted = FALSE)
  atomic_save(s, cfg$progress_file)
  log <- file.path(directory, "worker.log")
  process <- callr::r_bg(function(script, config, project_directory, r_options) {
    setwd(project_directory)
    options(r_options)
    options(pheno.spatial.functions_only = TRUE)
    source(script, local = .GlobalEnv)
    .pheno_progress$path <- config$progress_file
    .pheno_progress$state <- readRDS(config$progress_file)
    completed <- FALSE
    on.exit({
      .pheno_progress$state$finished <- TRUE
      .pheno_progress$state$aborted <- !completed
      .pheno_progress$state$ended_at <- as.numeric(Sys.time())
      .pheno_progress$state$phase <- if (completed) "Finished; inspect workflow_status.csv" else "Stopped with error"
      if (!completed && !is.null(.pheno_progress$state$models)) {
        m <- .pheno_progress$state$models
        m$status[m$status %in% c("checking", "running")] <- "failed"
        .pheno_progress$state$models <- m
      }
      progress_save()
    }, add = TRUE)
    result <- run_spatial_plasticity_worker(config)
    completed <- TRUE
    result
  }, args = list(script = snapshot, config = cfg, project_directory = getwd(),
    r_options = options()[intersect(names(options()), c("contrasts", "na.action", "warn"))]),
    libpath = .libPaths(), stdout = log, stderr = log, supervise = TRUE,
    user_profile = FALSE, system_profile = FALSE, wd = getwd())
  job <- list(process = process, progress_file = cfg$progress_file,
    monitor_directory = directory, output_root = cfg$output_root, log_file = log)
  options(pheno.spatial.active_job = job)
  writeLines(c(paste("Worker PID:", process$get_pid()), paste("Log:", log),
    "Use spatial_status() for elapsed time and a rough fitting ETA; watch_spatial() refreshes every minute.",
    "Ctrl+C during watch_spatial() pauses only the display. stop_spatial() terminates the fitting worker.",
    "Keep the parent R session open and the computer awake; closing R terminates this worker."),
    file.path(directory, "README.txt"))
  message("Started one fitting worker. Monitor: watch_spatial() or spatial_status().\n",
    "Keep R open. Detailed output: ", log)
  invisible(job)
}


  run_queue <- function(cutoffs_km = 30,
                        base_config = get0("peak_config", envir = .GlobalEnv,
                                           inherits = FALSE, ifnotfound = spatial_config),
                        output_root = NULL) {
    cutoffs_km <- unique(as.numeric(cutoffs_km))
    if (!length(cutoffs_km) || any(!is.finite(cutoffs_km)) ||
        any(cutoffs_km <= 0 | cutoffs_km == 15)) {
      stop("Supply positive sensitivity cutoffs other than the 15-km reference.")
    }
    cfg <- base_config
    cfg$analysis_groups <- "primary"
    cfg$stages <- "fit"
    cfg$check_only <- FALSE
    cfg$force_refit <- FALSE
    active <- getOption("pheno.spatial.active_job")
    primary_root <- cfg$output_root
    if (is.null(primary_root) && !is.null(active)) primary_root <- active$output_root
    if (is.null(primary_root)) {
      primary_root <- here::here("output", "phenology_plasticity", "spatial_models")
    }
    if (is.null(output_root)) output_root <- paste0(primary_root, "_mesh_sensitivity")
    primary_root <- normalizePath(primary_root, winslash = "/", mustWork = FALSE)
    proposed_root <- normalizePath(output_root, winslash = "/", mustWork = FALSE)
    if (tolower(proposed_root) == tolower(primary_root) ||
        startsWith(tolower(paste0(primary_root, "/")), tolower(paste0(proposed_root, "/")))) {
      stop("Sensitivity output must be separate from the primary output directory.")
    }
    dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
    output_root <- normalizePath(output_root, winslash = "/", mustWork = TRUE)
    tag <- paste0("queue_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid())
    archive <- file.path(output_root, paste0(tag, ".csv"))
    log_path <- file.path(output_root, paste0(tag, ".log"))
    note <- function(...) {
      line <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste0(..., collapse = ""))
      message(line)
      cat(line, "\n", file = log_path, append = TRUE)
    }
    analyses <- c("onset", "first_peak", "offset_univoltine", "offset_multivoltine")
    queue <- expand.grid(analysis = analyses, cutoff_km = cutoffs_km,
                         KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    queue$status <- "pending"
    queue$started_at <- queue$ended_at <- queue$message <- ""
    queue$output_root <- queue$monitor_directory <- queue$fit_directory <- ""
    queue$exit_status <- NA_real_
    flush_queue <- function() {
      write_table(queue, archive)
      write_table(queue, file.path(output_root, "queue_status.csv"))
    }
    flush_queue()
    writeLines(c(
      "Mesh sensitivity: primary reference is 15 km; sensitivity fits are stored separately.",
      paste("Requested cutoffs (km):", paste(cutoffs_km, collapse = ", ")),
      "Same original observations, climatic window, scales, fixed effects and random effects.",
      "Gaussian ML; annual IID Matern fields; no extra LRT, window or predictor fits.",
      "One new R subprocess per analysis; wait for the active process before starting.",
      "Matching saved sensitivity fits are reused. Failed jobs are logged; later jobs continue.",
      "The first ETA uses the 15-km pilot and is only a rough reference for another mesh.",
      "Compare coefficient estimates and uncertainty with the 15-km results; do not select by p-values.",
      "Keep the parent R session open and the computer awake.",
      paste("Reference output:", primary_root)
    ), file.path(output_root, "README.txt"))
    note("SENSITIVITY QUEUED: ", nrow(queue), " fits; cutoffs ", paste(cutoffs_km, collapse = ", "), " km.")
    note("Results: ", output_root)

    # Never start a sensitivity fit while the current worker is still alive.
    if (!is.null(active)) {
      if (active$process$is_alive()) note("Waiting for the current primary fit and its exports. It continues unchanged.")
      tryCatch(watch_spatial(active), error = function(e) {
        note("Previous worker reported: ", conditionMessage(e))
      })
      if (active$process$is_alive()) {
        stop("Queue paused while the current worker remains alive. No new fit started.")
      }
    }

    collected <- list()
    jobs <- list()
    for (i in seq_len(nrow(queue))) {
      analysis <- queue$analysis[i]
      cutoff <- queue$cutoff_km[i]
      cutoff_label <- gsub(".", "p", format(cutoff, trim = TRUE, scientific = FALSE), fixed = TRUE)
      cfg$analyses <- analysis
      cfg$mesh_cutoff_km <- cutoff
      # Per-analysis roots also preserve each job's manifest and workflow status.
      cfg$output_root <- file.path(output_root, paste0("cutoff_", cutoff_label, "km"), analysis)
      queue$output_root[i] <- cfg$output_root
      queue$status[i] <- "starting"
      queue$started_at[i] <- as.character(Sys.time())
      flush_queue()
      note("Starting ", i, "/", nrow(queue), ": ", analysis, " at ", cutoff, " km.")
      job <- NULL
      failure <- NULL
      result <- tryCatch({
        job <- run_spatial_plasticity(cfg)
        jobs[[paste(cutoff, analysis, sep = "__")]] <- job
        queue$monitor_directory[i] <- job$monitor_directory
        queue$status[i] <- "running"
        flush_queue()
        watch_spatial(job)
      }, error = function(e) {
        failure <<- conditionMessage(e)
        NULL
      })
      if (!is.null(job) && job$process$is_alive()) {
        queue$status[i] <- "running_queue_paused"
        queue$message[i] <- "Monitor interrupted or failed while worker remains alive. No following job started."
        flush_queue()
        stop(queue$message[i])
      }
      queue$ended_at[i] <- as.character(Sys.time())
      if (!is.null(job)) {
        code <- job$process$get_exit_status()
        if (length(code) == 1L) queue$exit_status[i] <- as.numeric(code)
      }
      successful <- is.null(failure) && is.list(result) && is.data.frame(result$status) &&
        nrow(result$status) == 1L && isTRUE(all(result$status$status == "OK"))
      if (successful) {
        queue$status[i] <- "OK"
        queue$fit_directory[i] <- result$manifest$output_directory[1]
        path <- file.path(queue$fit_directory[i], "main_interactions.csv")
        tryCatch({
          tab <- utils::read.csv(path, stringsAsFactors = FALSE)
          collected[[length(collected) + 1L]] <- cbind(
            sensitivity_cutoff_km = cutoff, sensitivity_analysis = analysis, tab)
          write_table(do.call(rbind, collected), file.path(output_root, "all_main_interactions_sensitivity.csv"))
        }, error = function(e) {
          queue$message[i] <<- paste("Fit/export stage OK; aggregate table warning:", conditionMessage(e))
        })
      } else {
        queue$status[i] <- "FAILED"
        if (!is.null(failure)) queue$message[i] <- failure
        else if (is.list(result) && is.data.frame(result$status)) {
          queue$message[i] <- paste(result$status$status, collapse = "; ")
        } else queue$message[i] <- "Worker returned no confirmed successful fit/export stage."
      }
      flush_queue()
      note(analysis, " at ", cutoff, " km: ", queue$status[i],
           if (nzchar(queue$message[i])) paste0(" -- ", queue$message[i]) else "")
    }
    note("QUEUE FINISHED: ", sum(queue$status == "OK"), "/", nrow(queue),
         " successful. Review queue_status.csv and per-model checks.")
    invisible(list(status = queue, jobs = jobs, output_root = output_root))
  }
  list(run = run_queue)
})

# This option is separate from pheno.spatial.functions_only.
# Sourcing this file queues the default four 30-km models automatically.
if (!isTRUE(getOption("pheno.mesh.sensitivity.functions_only", FALSE))) {
  pheno_mesh_sensitivity$run()
}
