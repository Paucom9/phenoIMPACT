# phenoIMPACT -- geographic patterns in phenological plasticity
# 2026-09-28 | startup fix v1.0.1 | R >= 4.1
#
# Run from the existing phenoIMPACT R project:
#   source("pheno_plasticity_latitude_15km.R")
#   latitude_job <- run_latitude()
#   watch_latitude()
# Sourcing only defines functions. run_latitude() starts ONE background worker.
# Keep R and the computer running. Escape/Ctrl+C pauses the monitor only.
#   latitude_status()   # current status
#   stop_latitude()     # stop worker; retain completed checkpoints
#
# Optional preflight (no fits):
#   cfg <- latitude_config
#   cfg$check_only <- TRUE
#   latitude_job <- run_latitude(cfg)
#   watch_latitude()
# To include the complementary first peak:
#   cfg <- latitude_config
#   cfg$analyses <- c(cfg$analyses, "first_peak")
#   latitude_job <- run_latitude(cfg)
# To begin with onset alone, use cfg$analyses <- "onset".
# Later rerun with all three; verified completed fits are reused.
#
# Input: output/phenology_plasticity/plasticity_main_models.rds
# Output: output/phenology_plasticity/latitude_models_15km/
# No changes to the input bundle or existing environmental-model outputs.
# Self-contained: no need to source the old spatial script.
#
# Full:    date ~ anomaly * latitude_z + random effects + annual spatial field
# Reduced: date ~ anomaly + latitude_z + identical random/spatial structure
# Random: (1|SITE_ID) + (1|SPECIES) + (0+anomaly|SPECIES_slope)
#         + (1|site_year_id)
# Gaussian ML; spatial="off"; spatiotemporal="iid"; time="YEAR".
# IID annual Matern fields share covariance parameters across years/species.
# Mesh cutoff 15 km is a mesh setting, NOT the estimated correlation range.
# Preserves original fitted rows, temperature window and anomaly scale.
# Latitude is reconstructed from verified EPSG:3035 metres, then standardized
# across distinct sites in each analysis. Degrees north remain in exports.
# Full and reduced models use identical rows and the same mesh.
# No environmental moderators are included: this tests the overall linear
# north-south association, combining within- and between-species information.
# It does not establish an intraspecific or causal geographic effect.
#
# Per analysis: fit_full/reduced.rds, convergence checks, fixed_effects.csv,
# interaction_test.csv (LRT and Wald kept distinct), latitude_slopes.csv,
# latitude_plasticity.png/pdf (no titles), latitude_support.csv and metadata.
# Root: geographic_tests.csv, model_manifest.csv, workflow_status.csv.
# Curves: fixed-effect temperature sensitivity; pointwise Wald 95% CIs from
# the full coefficient covariance matrix; not prediction intervals.
# Negative raw slopes mean advancement in warm years; positive mean delay.
# A sign-oriented column is also exported (-raw for onset/first peak/uni offset,
# +raw for multi offset); no absolute values are used.
# Anomaly units are inherited. Do not label days/degree C without recovering
# the original temperature scaling. The script does not rescale anomalies.
# No phylogenetic or within/between-species decomposition is performed here.
#
# Reused validation/fitting helpers originate from pheno_plasticity_spatial.R
# version 4 (2026-09-24). Source metadata and current package versions are saved.

latitude_config <- list(
  input_rds = NULL,
  output_root = NULL,
  analyses = c("onset", "offset_univoltine", "offset_multivoltine"),
  analysis_group = "primary",
  mesh_cutoff_km = 15,
  check_only = FALSE,
  force_refit = FALSE,
  gradient_threshold = 0.001,
  extra_optimization_attempts = 2L,
  fit_seed = 20260921L
)


pheno_latitude <- local({

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

check_packages <- function() {
  pkgs <- c("lme4", "lmerTest", "sdmTMB", "TMB", "fmesher", "Matrix", "sf",
            "dplyr", "tibble", "ggplot2", "here", "digest", "callr")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "),
    ". Install them before starting; this script does not install packages.")
  invisible(pkgs)
}

geographic_formula <- function(response, anomaly, interaction = TRUE) {
  operator <- if (interaction) " * " else " + "
  stats::as.formula(paste0(response, " ~ ", anomaly, operator, "latitude_z",
    " + (1 | SITE_ID) + (1 | SPECIES) + (0 + ", anomaly,
    " | SPECIES_slope) + (1 | site_year_id)"), env = baseenv())
}

prepare_geography <- function(entry, coordinates, response) {
  # Reuse the original alignment/design/random-effect validation BEFORE
  # replacing the environmental fixed effects with latitude.
  original <- prepare_analysis(entry, coordinates, response)
  keep <- unique(c(response, original$anomaly, "SITE_ID", "SPECIES", "SPECIES_slope",
    "YEAR", "bms_id", "source_row", "site_year_id", "x_3035", "y_3035", "x_km", "y_km"))
  d <- original$data[, keep, drop = FALSE]
  sites <- unique(d[, c("SITE_ID", "bms_id", "x_3035", "y_3035")])
  if (anyDuplicated(sites$SITE_ID)) stop("Conflicting site coordinates.")
  geo <- sf::st_as_sf(sites, coords = c("x_3035", "y_3035"), crs = 3035)
  ll <- sf::st_coordinates(sf::st_transform(geo, 4326))
  sites$longitude_deg <- ll[, 1]
  sites$latitude_deg <- ll[, 2]
  if (any(!is.finite(ll)) || any(abs(sites$latitude_deg) > 90) ||
      any(abs(sites$longitude_deg) > 180)) stop("Invalid coordinate transformation.")
  center <- mean(sites$latitude_deg)
  spread <- stats::sd(sites$latitude_deg)
  if (!is.finite(spread) || spread <= 0) stop("No usable latitudinal gradient.")
  index <- match(as.character(d$SITE_ID), as.character(sites$SITE_ID))
  if (anyNA(index)) stop("Latitude join failed.")
  d$latitude_deg <- sites$latitude_deg[index]
  d$latitude_z <- (d$latitude_deg - center) / spread
  full <- geographic_formula(response, original$anomaly, TRUE)
  reduced <- geographic_formula(response, original$anomaly, FALSE)
  validate_fixed_design(full, d)
  validate_fixed_design(reduced, d)
  tf <- attr(stats::terms(lme4::nobars(full)), "term.labels")
  tr <- attr(stats::terms(lme4::nobars(reduced)), "term.labels")
  interaction <- paste(original$anomaly, "latitude_z", sep = ":")
  if (!setequal(setdiff(tf, tr), interaction) || !all(tr %in% tf)) {
    stop("Full and reduced models differ by more than the latitude interaction.")
  }
  list(data = d, sites = sites, formula = full, reduced_formula = reduced,
    response = response, anomaly = original$anomaly, interaction = interaction,
    best_window = original$best_window, latitude_center = center,
    latitude_sd = spread, original_formula = original$original_formula)
}

.monitor <- new.env(parent = emptyenv())
.monitor$path <- NULL
.monitor$state <- NULL

progress_save <- function() {
  if (is.null(.monitor$path) || is.null(.monitor$state)) return(invisible(NULL))
  .monitor$state$updated_at <- as.numeric(Sys.time())
  atomic_save(.monitor$state, .monitor$path)
}

progress_phase <- function(text) {
  if (is.null(.monitor$state)) return(invisible(NULL))
  .monitor$state$phase <- text
  progress_save()
}

progress_fit <- function(path, status, duration_seconds = NA_real_) {
  if (is.null(.monitor$state)) return(invisible(NULL))
  m <- .monitor$state$models
  i <- match(path, m$path)
  if (!is.na(i)) {
    m$status[i] <- status
    if (status == "running") m$started_at[i] <- as.numeric(Sys.time())
    if (is.finite(duration_seconds)) m$duration_seconds[i] <- duration_seconds
    .monitor$state$models <- m
  }
  .monitor$state$phase <- paste(status, basename(dirname(path)), basename(path))
  progress_save()
}

empty_test <- function(analysis, p) {
  data.frame(analysis = analysis, window_days = p$best_window,
    n_observations = nrow(p$data), n_sites = nlevels(p$data$SITE_ID),
    n_species = nlevels(p$data$SPECIES), term = p$interaction,
    estimate_per_SD_latitude = NA_real_, SE = NA_real_,
    lower_95 = NA_real_, upper_95 = NA_real_, p_Wald = NA_real_,
    estimate_per_10_degrees = NA_real_, lower_95_per_10_degrees = NA_real_,
    upper_95_per_10_degrees = NA_real_,
    Chisq = NA_real_, LRT_df = NA_integer_, p_LRT = NA_real_,
    delta_AIC_reduced_minus_full = NA_real_,
    full_status = "not_run", LRT_status = "not_run")
}

write_full_results <- function(saved, p, analysis, out_dir, result) {
  fixed <- normal_ci(sdmTMB::tidy(saved$fit, effects = "fixed", conf.int = FALSE))
  expected <- colnames(stats::model.matrix(lme4::nobars(p$formula), p$data))
  if (!setequal(fixed$term, expected)) stop("Unexpected fitted fixed-effect terms.")
  row <- match(p$interaction, fixed$term)
  if (is.na(row)) stop("Latitude interaction missing from fitted model.")
  z <- fixed[row, ]
  result$estimate_per_SD_latitude <- z$estimate
  result$SE <- z$SE
  result$lower_95 <- z$lower_95
  result$upper_95 <- z$upper_95
  result$p_Wald <- z$p_Wald
  result$estimate_per_10_degrees <- z$estimate * 10 / p$latitude_sd
  result$lower_95_per_10_degrees <- z$lower_95 * 10 / p$latitude_sd
  result$upper_95_per_10_degrees <- z$upper_95 * 10 / p$latitude_sd
  result$full_status <- "OK"
  write_table(fixed, file.path(out_dir, "fixed_effects.csv"))
  write_table(result, file.path(out_dir, "interaction_test.csv"))

  b <- stats::setNames(fixed$estimate, fixed$term)
  V <- fixed_covariance(saved$fit, names(b))
  if (!isTRUE(all.equal(unname(sqrt(diag(V))), fixed$SE, tolerance = 1e-5))) {
    stop("Coefficient covariance does not match fitted standard errors.")
  }
  lat <- seq(min(p$sites$latitude_deg), max(p$sites$latitude_deg), length.out = 150)
  lat_z <- (lat - p$latitude_center) / p$latitude_sd
  L <- matrix(0, length(lat), length(b), dimnames = list(NULL, names(b)))
  L[, p$anomaly] <- 1
  L[, p$interaction] <- lat_z
  estimate <- as.numeric(L %*% b)
  variance <- rowSums((L %*% V) * L)
  if (any(variance < -1e-8)) stop("Negative slope variance; inspect covariance.")
  se <- sqrt(pmax(variance, 0))
  q <- stats::qnorm(0.975)
  direction <- if (analysis == "offset_multivoltine") 1 else -1
  slopes <- data.frame(analysis = analysis, latitude_deg = lat, latitude_z = lat_z,
    sensitivity = estimate, SE = se, lower_95 = estimate - q * se,
    upper_95 = estimate + q * se, sign_oriented_sensitivity = direction * estimate,
    oriented_lower_95 = direction * estimate - q * se,
    oriented_upper_95 = direction * estimate + q * se,
    anomaly_units = "saved input anomaly unit; scale unchanged",
    interval = "pointwise fixed-effect normal Wald 95%")
  write_table(slopes, file.path(out_dir, "latitude_slopes.csv"))
  graph <- ggplot2::ggplot(slopes, ggplot2::aes(latitude_deg, sensitivity)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey55") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower_95, ymax = upper_95),
      fill = "#3B718F", alpha = 0.22) +
    ggplot2::geom_line(colour = "#24546F", linewidth = 0.8) +
    ggplot2::geom_rug(data = p$sites, ggplot2::aes(x = latitude_deg),
      inherit.aes = FALSE, sides = "b", alpha = 0.2) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::labs(x = "Latitude (degrees N)",
      y = "Temperature sensitivity (days per anomaly unit)")
  ggplot2::ggsave(file.path(out_dir, "latitude_plasticity.png"), graph,
    width = 6.5, height = 4.5, dpi = 300)
  ggplot2::ggsave(file.path(out_dir, "latitude_plasticity.pdf"), graph,
    width = 6.5, height = 4.5)
  result
}

compare_models <- function(full, reduced, result) {
  if (!identical(full$data$source_row, reduced$data$source_row)) stop("LRT rows differ.")
  lf <- stats::logLik(full)
  lr <- stats::logLik(reduced)
  df <- attr(lf, "df") - attr(lr, "df")
  chi <- 2 * (as.numeric(lf) - as.numeric(lr))
  if (!is.finite(chi) || length(df) != 1L || df != 1L) stop("Invalid interaction LRT.")
  if (chi < -1e-3) stop("Reduced likelihood exceeds full; recheck optimization.")
  result$Chisq <- max(0, chi)
  result$LRT_df <- df
  result$p_LRT <- stats::pchisq(result$Chisq, df, lower.tail = FALSE)
  result$delta_AIC_reduced_minus_full <- stats::AIC(reduced) - stats::AIC(full)
  result$LRT_status <- "OK"
  result
}

write_support <- function(p, out_dir) {
  support <- p$data |>
    dplyr::mutate(latitude_band_start = 5 * floor(latitude_deg / 5)) |>
    dplyr::group_by(latitude_band_start) |>
    dplyr::summarise(n_observations = dplyr::n(),
      n_sites = dplyr::n_distinct(SITE_ID), n_species = dplyr::n_distinct(SPECIES),
      .groups = "drop")
  write_table(support, file.path(out_dir, "latitude_support.csv"))
  write_table(p$sites, file.path(out_dir, "site_latitudes.csv"))
}

worker <- function(cfg) {
  .monitor$path <- cfg$progress_file
  .monitor$state <- readRDS(cfg$progress_file)
  completed <- FALSE
  on.exit({
    .monitor$state$finished <- TRUE
    .monitor$state$aborted <- !completed
    .monitor$state$ended_at <- as.numeric(Sys.time())
    .monitor$state$phase <- if (!completed) "Stopped with error; inspect worker.log" else
      if (.monitor$state$stage_errors > 0) "Finished with errors; inspect workflow_status.csv" else
      if (cfg$check_only) "Preflight complete; no models fitted" else "Finished"
    progress_save()
  }, add = TRUE)
  pkgs <- check_packages()
  versions <- data.frame(package = pkgs, version = vapply(pkgs, function(x)
    as.character(utils::packageVersion(x)), character(1)))
  writeLines(capture.output(sessionInfo()), file.path(cfg$monitor_directory, "sessionInfo.txt"))
  progress_phase("Reading original model bundle")
  bundle <- readRDS(cfg$input_rds)
  entries <- bundle[[cfg$analysis_group]]
  if (is.null(entries) || !all(cfg$analyses %in% names(entries))) {
    stop("Requested analyses are missing from the input bundle.")
  }
  responses <- c(onset = "ONSET_mean", first_peak = "FIRST_PEAK",
    offset_univoltine = "OFFSET_mean", offset_multivoltine = "OFFSET_mean")
  jobs <- list()
  manifest <- list()
  model_progress <- list()
  initial_tests <- list()
  for (analysis in cfg$analyses) {
    progress_phase(paste("Validating data, latitude and mesh:", analysis))
    p <- prepare_geography(entries[[analysis]], bundle$site_coordinates, responses[[analysis]])
    run_id <- digest::digest(list(version = "latitude_1.0", data = p$data,
      formula = formula_text(p$formula), reduced = formula_text(p$reduced_formula),
      cutoff = cfg$mesh_cutoff_km, seed = cfg$fit_seed, latitude_center = p$latitude_center,
      latitude_sd = p$latitude_sd, packages = versions,
      threshold = cfg$gradient_threshold, extra_optimization = cfg$extra_optimization_attempts),
      algo = "sha256")
    out_dir <- file.path(cfg$output_root, paste0(cfg$analysis_group, "__", analysis,
      "_", substr(run_id, 1, 12)))
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    metadata_path <- file.path(out_dir, "metadata.rds")
    if (file.exists(metadata_path)) {
      meta <- readRDS(metadata_path)
      if (!identical(meta$run_id, run_id)) stop("Metadata mismatch.")
      tri <- meta$triangulation
      checked_mesh <- project_mesh(p$data, tri)
      rm(checked_mesh)
    } else {
      built <- build_mesh(p$data, cfg)
      tri <- built$triangulation
      rm(built)
    }
    validation <- data.frame(analysis_group = cfg$analysis_group, analysis = analysis,
      response = p$response, selected_window_days = p$best_window,
      n_observations = nrow(p$data), n_sites = nrow(p$sites),
      n_species = nlevels(p$data$SPECIES), n_site_years = nlevels(p$data$site_year_id),
      latitude_min = min(p$sites$latitude_deg), latitude_max = max(p$sites$latitude_deg),
      latitude_center = p$latitude_center, latitude_sd = p$latitude_sd,
      mesh_cutoff_km = cfg$mesh_cutoff_km, n_mesh_vertices = tri$mesh$n,
      window_selection = "preserved_original", fitted_rows_verified = TRUE,
      latitude_scaling = "equal weight per distinct site", anomaly_scaling = "inherited")
    write_table(validation, file.path(out_dir, "validation.csv"))
    write_support(p, out_dir)
    atomic_save(list(run_id = run_id, input_rds = cfg$input_rds,
      input_bundle_created = bundle$created_at, settings = cfg, validation = validation,
      formula = p$formula, reduced_formula = p$reduced_formula,
      original_formula = p$original_formula, package_versions = versions,
      triangulation = tri, helper_source = "pheno_plasticity_spatial.R v4 2026-09-24"),
      metadata_path)
    atomic_save(p, file.path(out_dir, "prepared_input.rds"))
    jobs[[analysis]] <- list(out_dir = out_dir, run_id = run_id)
    initial_tests[[analysis]] <- empty_test(analysis, p)
    manifest[[analysis]] <- cbind(validation, run_id = run_id, output_directory = out_dir)
    model_progress[[analysis]] <- data.frame(analysis = analysis,
      model = c("full", "reduced"), path = file.path(out_dir, c("fit_full.rds", "fit_reduced.rds")),
      status = "pending", started_at = NA_real_, duration_seconds = NA_real_)
    rm(p, tri)
    invisible(gc())
  }
  rm(bundle, entries)
  invisible(gc())
  manifest <- dplyr::bind_rows(manifest)
  write_table(manifest, file.path(cfg$output_root, "model_manifest.csv"))
  .monitor$state$models <- dplyr::bind_rows(model_progress)
  progress_save()
  if (cfg$check_only) {
    .monitor$state$models$status <- "not_run_preflight"
    completed <- TRUE
    return(manifest)
  }

  # Finish one phenophase at a time: a full result is available after each pair.
  workflow <- list()
  all_tests <- initial_tests
  write_table(dplyr::bind_rows(all_tests), file.path(cfg$output_root, "geographic_tests.csv"))
  write_table(data.frame(analysis = cfg$analyses, stage = "pending", status = "not_run",
    message = ""), file.path(cfg$output_root, "workflow_status.csv"))
  for (analysis in names(jobs)) {
    job <- jobs[[analysis]]
    out_dir <- job$out_dir
    p <- readRDS(file.path(out_dir, "prepared_input.rds"))
    meta <- readRDS(file.path(out_dir, "metadata.rds"))
    mesh <- project_mesh(p$data, meta$triangulation)
    result <- empty_test(analysis, p)
    # Clear any old reported LRT before retrying; stored model fits remain.
    write_table(result, file.path(out_dir, "interaction_test.csv"))
    if (cfg$force_refit) {
      # These are derived files belonging only to this geographic analysis.
      # Remove them before a forced refit so a failed rerun cannot show old curves.
      unlink(file.path(out_dir, c("fixed_effects.csv", "latitude_slopes.csv",
        "latitude_plasticity.png", "latitude_plasticity.pdf")))
    }
    full <- NULL
    reduced <- NULL
    full_error <- NULL
    full <- tryCatch({
      fit_model("full", p$data, p$formula, mesh, job$run_id, out_dir, cfg)
    }, error = function(e) {
      full_error <<- conditionMessage(e)
      NULL
    })
    if (is.null(full)) {
      result$full_status <- paste("ERROR:", full_error)
      result$LRT_status <- "skipped_full_fit_failed"
      progress_fit(file.path(out_dir, "fit_reduced.rds"), "skipped")
      workflow[[paste0(analysis, "_full")]] <- data.frame(analysis = analysis,
        stage = "full", status = "ERROR", message = full_error)
      .monitor$state$stage_errors <- .monitor$state$stage_errors + 1L
    } else {
      export_error <- NULL
      result <- tryCatch(write_full_results(full, p, analysis, out_dir, result),
        error = function(e) {
          export_error <<- conditionMessage(e)
          result$full_status <- paste("fit_OK_export_ERROR:", export_error)
          result
        })
      workflow[[paste0(analysis, "_full")]] <- data.frame(analysis = analysis,
        stage = "full", status = if (is.null(export_error)) "OK" else "ERROR",
        message = if (is.null(export_error)) "" else export_error)
      if (!is.null(export_error)) .monitor$state$stage_errors <- .monitor$state$stage_errors + 1L
      lrt_error <- NULL
      result <- tryCatch({
        reduced <- fit_model("reduced", p$data, p$reduced_formula, mesh,
          job$run_id, out_dir, cfg)
        compare_models(full$fit, reduced$fit, result)
      }, error = function(e) {
        lrt_error <<- conditionMessage(e)
        result$LRT_status <- paste("ERROR:", lrt_error)
        result
      })
      workflow[[paste0(analysis, "_lrt")]] <- data.frame(analysis = analysis,
        stage = "LRT", status = if (is.null(lrt_error)) "OK" else "ERROR",
        message = if (is.null(lrt_error)) "" else lrt_error)
      if (!is.null(lrt_error)) .monitor$state$stage_errors <- .monitor$state$stage_errors + 1L
    }
    all_tests[[analysis]] <- result
    write_table(result, file.path(out_dir, "interaction_test.csv"))
    write_table(dplyr::bind_rows(all_tests), file.path(cfg$output_root, "geographic_tests.csv"))
    write_table(dplyr::bind_rows(workflow), file.path(cfg$output_root, "workflow_status.csv"))
    progress_save()
    rm(p, meta, mesh, full, reduced)
    invisible(gc())
  }
  completed <- TRUE
  dplyr::bind_rows(all_tests)
}

start <- function(cfg) {
  check_packages()
  allowed <- c("onset", "first_peak", "offset_univoltine", "offset_multivoltine")
  if (!length(cfg$analyses) || anyNA(cfg$analyses) || any(!cfg$analyses %in% allowed)) {
    stop("Unknown or empty analyses.")
  }
  cfg$analyses <- unique(cfg$analyses)
  if (length(cfg$analysis_group) != 1L || is.na(cfg$analysis_group) ||
      !cfg$analysis_group %in% c("primary", "sensitivity_1zero")) stop("Invalid analysis_group.")
  if (!identical(as.numeric(cfg$mesh_cutoff_km), 15)) stop("This workflow uses a 15-km mesh cutoff.")
  for (nm in c("check_only", "force_refit")) {
    if (!is.logical(cfg[[nm]]) || length(cfg[[nm]]) != 1L || is.na(cfg[[nm]])) stop("Invalid ", nm)
  }
  if (length(cfg$gradient_threshold) != 1L || !is.finite(cfg$gradient_threshold) ||
      cfg$gradient_threshold <= 0) stop("Invalid gradient_threshold.")
  if (length(cfg$extra_optimization_attempts) != 1L || !is.finite(cfg$extra_optimization_attempts) ||
      cfg$extra_optimization_attempts < 0 || cfg$extra_optimization_attempts %% 1 != 0) {
    stop("Invalid extra_optimization_attempts.")
  }
  if (length(cfg$fit_seed) != 1L || !is.finite(cfg$fit_seed) || cfg$fit_seed %% 1 != 0 ||
      cfg$fit_seed < 0 || cfg$fit_seed > .Machine$integer.max) stop("Invalid fit_seed.")
  for (key in c("pheno.latitude.active_job", "pheno.spatial.active_job")) {
    previous <- getOption(key)
    if (!is.null(previous) && isTRUE(tryCatch(previous$process$is_alive(), error = function(e) FALSE))) {
      stop("An analysis worker is already running (", key, "). Let it finish before starting another.")
    }
  }
  if (is.null(cfg$input_rds)) cfg$input_rds <- here::here("output", "phenology_plasticity", "plasticity_main_models.rds")
  if (!file.exists(cfg$input_rds)) stop("Input not found: ", cfg$input_rds)
  cfg$input_rds <- normalizePath(cfg$input_rds, winslash = "/", mustWork = TRUE)
  if (is.null(cfg$output_root)) cfg$output_root <- here::here("output", "phenology_plasticity", "latitude_models_15km")
  dir.create(cfg$output_root, recursive = TRUE, showWarnings = FALSE)
  cfg$output_root <- normalizePath(cfg$output_root, winslash = "/", mustWork = TRUE)
  monitor <- tempfile(pattern = paste0("monitor_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_"),
    tmpdir = cfg$output_root)
  dir.create(monitor)
  cfg$monitor_directory <- monitor
  cfg$progress_file <- file.path(monitor, "progress_state.rds")
  log <- file.path(monitor, "worker.log")
  state <- list(started_at = as.numeric(Sys.time()), updated_at = as.numeric(Sys.time()),
    ended_at = NA_real_, phase = "Starting worker", models = NULL,
    stage_errors = 0L, finished = FALSE, aborted = FALSE)
  atomic_save(state, cfg$progress_file)
  # Preserve the toolkit's lexical environment: callr otherwise resets it to
  # .GlobalEnv and loses .monitor and every locally defined helper.
  # This local environment contains only the toolkit, not the user's data.
  process <- callr::r_bg(worker, args = list(cfg = cfg), libpath = .libPaths(),
    stdout = log, stderr = log, supervise = TRUE, user_profile = FALSE,
    system_profile = FALSE, package = TRUE, wd = getwd())
  job <- list(process = process, progress_file = cfg$progress_file,
    monitor_directory = monitor, log_file = log, output_root = cfg$output_root)
  options(pheno.latitude.active_job = job)
  message("Started one latitude worker (", length(cfg$analyses), " analyses; ",
    if (cfg$check_only) "preflight only" else paste(2 * length(cfg$analyses), "fits"), ").\n",
    "Use watch_latitude() or latitude_status(). Keep R open and the computer awake.\n",
    "Log: ", log)
  invisible(job)
}

status <- function(job = getOption("pheno.latitude.active_job"), quiet = FALSE) {
  if (is.null(job)) stop("No latitude job. Run run_latitude() first.")
  s <- tryCatch(readRDS(job$progress_file), error = function(e) NULL)
  if (is.null(s)) {
    if (!quiet) message("Progress is being updated; try again shortly.")
    return(invisible(NULL))
  }
  alive <- job$process$is_alive()
  if (!alive && !s$finished) {
    s$finished <- TRUE
    s$aborted <- TRUE
    ended <- tryCatch(as.numeric(job$process$get_end_time()), error = function(e) NA_real_)
    s$ended_at <- if (length(ended) == 1L && is.finite(ended)) ended else as.numeric(Sys.time())
    s$phase <- "Worker stopped; inspect worker.log"
    if (!is.null(s$models)) s$models$status[s$models$status %in% c("running", "checking")] <- "failed"
    atomic_save(s, job$progress_file)
  }
  now <- if (isTRUE(s$finished) && is.finite(s$ended_at)) s$ended_at else as.numeric(Sys.time())
  m <- s$models
  total <- if (is.null(m)) 0L else nrow(m)
  done <- if (is.null(m)) 0L else sum(m$status %in% c("done", "cached"))
  failed <- if (is.null(m)) 0L else sum(m$status %in% c("failed", "skipped"))
  if (!quiet) {
    cat(sprintf("Elapsed: %.1f h | completed/cached fits: %d/%d | failed/skipped: %d\n",
      (now - s$started_at) / 3600, done, total, failed))
    cat("Current:", s$phase, "| worker running:", alive, "\n")
    if (!is.null(m)) {
      running <- which(m$status == "running")
      if (length(running)) cat(sprintf("Current fit elapsed: %.1f h\n",
        (now - m$started_at[running[1]]) / 3600))
    }
    if (s$stage_errors > 0) cat("Stage errors:", s$stage_errors, "(see workflow_status.csv)\n")
  }
  invisible(s)
}

watch <- function(job = getOption("pheno.latitude.active_job"), interval_seconds = 60) {
  if (is.null(job)) stop("No latitude job. Run run_latitude() first.")
  if (length(interval_seconds) != 1L || !is.finite(interval_seconds) || interval_seconds < 1) {
    stop("interval_seconds must be >= 1.")
  }
  tryCatch({
    repeat {
      status(job)
      if (!job$process$is_alive()) break
      deadline <- Sys.time() + interval_seconds
      while (Sys.time() < deadline && job$process$is_alive()) Sys.sleep(1)
    }
    result <- job$process$get_result()
    message("Worker ended. Results and workflow_status.csv: ", job$output_root)
    invisible(result)
  }, interrupt = function(e) {
    message("Monitor paused; FITTING CONTINUES. Use latitude_status(), watch_latitude(), or stop_latitude().")
    invisible(job)
  })
}

stop_job <- function(job = getOption("pheno.latitude.active_job")) {
  if (is.null(job)) return(invisible(FALSE))
  if (job$process$is_alive()) job$process$kill_tree()
  message("Latitude worker stopped. Saved fits are retained; an unfinished fit must restart.")
  invisible(TRUE)
}

list(start = start, status = status, watch = watch, stop = stop_job)
})

run_latitude <- function(cfg = latitude_config) pheno_latitude$start(cfg)
latitude_status <- function(job = getOption("pheno.latitude.active_job"), quiet = FALSE) {
  pheno_latitude$status(job, quiet)
}
watch_latitude <- function(job = getOption("pheno.latitude.active_job"), interval_seconds = 60) {
  pheno_latitude$watch(job, interval_seconds)
}
stop_latitude <- function(job = getOption("pheno.latitude.active_job")) pheno_latitude$stop(job)

message("Latitude functions loaded. Run latitude_job <- run_latitude(); then watch_latitude().")
