# ============================================================================================ #
# pheno_plasticity_spatial_field_sensitivity.R
#
# Test a continuous spatiotemporal correlation structure for the selected
# phenology-plasticity models.
#
# The spatial model adds an independent Matérn spatial field for every year,
# with one spatial range shared across years. This is the continuous-covariance
# analogue of the earlier 50-km x year block sensitivity, but without arbitrary
# grid boundaries.
#
# Main comparison:
#   site-year sdmTMB model
#     vs.
#   the same model + IID spatiotemporal Matérn field
#
# The model retains:
#   (1 | SITE_ID)
#   (1 | SPECIES)
#   (0 + anomaly | SPECIES)
#   (1 | site_year_id)
#
# Requires:
#   output/phenology_plasticity/plasticity_main_models.rds
#
# Outputs:
#   output/phenology_plasticity/diagnostics/spatial_field/
#     plasticity_spatial_field_models.rds
#     plasticity_spatial_field_comparison.xlsx
#     spatial_field_moran_comparison.png
#
# Coordinates in the saved model bundle are EPSG:3035 metres.
# ============================================================================================ #


# ---------------------------------------------------------------------------- #
# Settings
# ---------------------------------------------------------------------------- #

analysis_group <- "primary"

# Start with onset as a pilot. After it runs successfully, use:
# models_to_run <- "all"
models_to_run <- c("onset")

# The previous correlogram showed dependence mainly within 25-50 km.
# A 15-km mesh cutoff gives enough resolution for that scale without placing a
# mesh vertex at every transect.
mesh_cutoff_km <- 15

distance_bands <- tibble::tribble(
  ~distance_band, ~d_min_km, ~d_max_km,
  "0-25 km",               0,         25,
  "25-50 km",             25,         50
)

n_permutations <- 499L
minimum_connected_sites <- 10L

# Re-run one extra optimization cycle only if the first fit has a bad Hessian
# or a maximum absolute gradient above this threshold.
gradient_threshold <- 0.001

set.seed(20260724)


# ---------------------------------------------------------------------------- #
# Packages and paths
# ---------------------------------------------------------------------------- #

required_packages <- c(
  "dplyr", "tidyr", "purrr", "tibble", "ggplot2",
  "lme4", "sdmTMB", "spdep", "writexl", "here"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Install the missing packages first: ",
    paste(missing_packages, collapse = ", "),
    "\nFor sdmTMB use: install.packages(\"sdmTMB\")"
  )
}

rds_path <- here::here(
  "output",
  "phenology_plasticity",
  "plasticity_main_models.rds"
)

if (!file.exists(rds_path)) {
  stop("Model bundle not found: ", rds_path)
}

out_dir <- here::here(
  "output",
  "phenology_plasticity",
  "diagnostics",
  "spatial_field"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

model_bundle <- readRDS(rds_path)


# ---------------------------------------------------------------------------- #
# Recover entries, exact fitted rows, and coordinates
# ---------------------------------------------------------------------------- #

match_entry_data <- function(entry) {
  fitted_frame <- stats::model.frame(entry$model)
  original_data <- entry$data

  row_index <- match(rownames(fitted_frame), rownames(original_data))

  if (anyNA(row_index)) {
    if (nrow(fitted_frame) == nrow(original_data)) {
      warning(
        "Row names did not match; using row order because row counts agree."
      )
      row_index <- seq_len(nrow(original_data))
    } else {
      stop("Could not match the saved rows to the fitted model.")
    }
  }

  original_data[row_index, , drop = FALSE]
}

model_entries <- model_bundle[[analysis_group]]

if (is.null(model_entries) || length(model_entries) == 0) {
  stop("No models found for analysis group: ", analysis_group)
}

if (!identical(models_to_run, "all")) {
  unknown_models <- setdiff(models_to_run, names(model_entries))

  if (length(unknown_models) > 0) {
    stop(
      "Unknown models in models_to_run: ",
      paste(unknown_models, collapse = ", "),
      "\nAvailable models: ",
      paste(names(model_entries), collapse = ", ")
    )
  }

  model_entries <- model_entries[models_to_run]
}

names(model_entries) <- paste(
  analysis_group,
  names(model_entries),
  sep = "__"
)

site_coordinates <- model_bundle$site_coordinates |>
  dplyr::transmute(
    SITE_ID = as.character(SITE_ID),
    bms_id = as.character(bms_id),
    x_3035 = as.numeric(x_3035),
    y_3035 = as.numeric(y_3035)
  ) |>
  dplyr::filter(
    !is.na(SITE_ID),
    !is.na(bms_id),
    is.finite(x_3035),
    is.finite(y_3035)
  ) |>
  dplyr::distinct(SITE_ID, bms_id, .keep_all = TRUE)

prepare_entry_data <- function(entry, model_name) {
  d <- match_entry_data(entry) |>
    dplyr::mutate(
      SITE_ID = as.character(SITE_ID),
      SPECIES = as.character(SPECIES),
      bms_id = as.character(bms_id),
      YEAR = as.integer(YEAR)
    ) |>
    dplyr::select(-dplyr::any_of(c("x_3035", "y_3035"))) |>
    dplyr::left_join(
      site_coordinates,
      by = c("SITE_ID", "bms_id")
    ) |>
    dplyr::mutate(
      SITE_ID = factor(SITE_ID),
      SPECIES = factor(SPECIES),
      bms_id = factor(bms_id),
      site_year_id = interaction(
        SITE_ID,
        YEAR,
        drop = TRUE,
        lex.order = TRUE
      ),
      x_km = x_3035 / 1000,
      y_km = y_3035 / 1000
    )

  if (any(!is.finite(d$x_km)) || any(!is.finite(d$y_km))) {
    stop(model_name, ": some fitted observations lack EPSG:3035 coordinates.")
  }

  if (anyNA(d$YEAR)) {
    stop(model_name, ": YEAR contains missing values.")
  }

  d
}

prepared_data <- purrr::imap(
  model_entries,
  prepare_entry_data
)


# ---------------------------------------------------------------------------- #
# Build one shared mesh from all sites used by the selected models
# ---------------------------------------------------------------------------- #

mesh_data <- purrr::map_dfr(
  prepared_data,
  ~ dplyr::distinct(.x, x_km, y_km)
) |>
  dplyr::distinct(x_km, y_km)

message(
  "Building a shared spatial mesh from ",
  nrow(mesh_data),
  " unique coordinates..."
)

spatial_mesh <- sdmTMB::make_mesh(
  data = mesh_data,
  xy_cols = c("x_km", "y_km"),
  cutoff = mesh_cutoff_km,
  type = "cutoff"
)

n_mesh_vertices <- spatial_mesh$mesh$n

message(
  "Mesh built with ",
  n_mesh_vertices,
  " vertices and a ",
  mesh_cutoff_km,
  "-km cutoff."
)


# ---------------------------------------------------------------------------- #
# Recreate the exact fixed effects and known random structure
# ---------------------------------------------------------------------------- #

make_sdm_formula <- function(entry) {
  original_formula <- stats::formula(entry$model)
  fixed_formula <- lme4::nobars(original_formula)

  response_text <- paste(
    deparse(fixed_formula[[2]], width.cutoff = 500L),
    collapse = ""
  )

  fixed_rhs_text <- paste(
    deparse(fixed_formula[[3]], width.cutoff = 500L),
    collapse = ""
  )

  anomaly <- paste0("clim_anomaly_tw", entry$best_window)

  # The two SPECIES terms reproduce:
  #   (1 + anomaly || SPECIES)
  # without requiring || support in a particular sdmTMB version.
  formula_text <- paste0(
    response_text,
    " ~ ",
    fixed_rhs_text,
    " + (1 | SITE_ID)",
    " + (1 | SPECIES)",
    " + (0 + ",
    anomaly,
    " | SPECIES)",
    " + (1 | site_year_id)"
  )

  stats::as.formula(
    formula_text,
    env = environment(original_formula)
  )
}

sdm_formulas <- purrr::map(model_entries, make_sdm_formula)


# ---------------------------------------------------------------------------- #
# Fit matched non-spatial and spatiotemporal-field models
# ---------------------------------------------------------------------------- #

max_abs_gradient <- function(model) {
  if (is.null(model$gradients) || length(model$gradients) == 0) {
    return(NA_real_)
  }

  max(abs(model$gradients), na.rm = TRUE)
}

pd_hessian <- function(model) {
  value <- tryCatch(
    model$sd_report$pdHess,
    error = function(e) NA
  )

  if (length(value) == 0) NA else isTRUE(value)
}

fit_one_sdm_model <- function(
    formula,
    data,
    mesh,
    include_spatial_field,
    model_name) {

  variant <- if (include_spatial_field) {
    "spatial_field"
  } else {
    "site_year"
  }

  message("Fitting ", model_name, " [", variant, "]...")

  start_time <- proc.time()[["elapsed"]]

  fitted_model <- sdmTMB::sdmTMB(
    formula = formula,
    data = data,
    mesh = mesh,
    time = "YEAR",
    family = stats::gaussian(link = "identity"),
    spatial = "off",
    spatiotemporal = if (include_spatial_field) "iid" else "off",
    reml = FALSE,
    silent = FALSE
  )

  first_gradient <- max_abs_gradient(fitted_model)
  first_pd_hessian <- pd_hessian(fitted_model)

  needs_extra_optimization <- (
    !isTRUE(first_pd_hessian) ||
      !is.finite(first_gradient) ||
      first_gradient > gradient_threshold
  )

  if (needs_extra_optimization) {
    message(
      model_name,
      " [",
      variant,
      "] needs one extra optimization cycle."
    )

    fitted_model <- sdmTMB::run_extra_optimization(
      fitted_model,
      nlminb_loops = 1,
      newton_loops = 1
    )
  }

  elapsed_seconds <- proc.time()[["elapsed"]] - start_time

  list(
    model = fitted_model,
    elapsed_seconds = elapsed_seconds,
    extra_optimization = needs_extra_optimization
  )
}

fit_one_model_pair <- function(entry, model_name) {
  d <- prepared_data[[model_name]]
  model_formula <- sdm_formulas[[model_name]]

  nonspatial <- fit_one_sdm_model(
    formula = model_formula,
    data = d,
    mesh = spatial_mesh,
    include_spatial_field = FALSE,
    model_name = model_name
  )

  spatial <- fit_one_sdm_model(
    formula = model_formula,
    data = d,
    mesh = spatial_mesh,
    include_spatial_field = TRUE,
    model_name = model_name
  )

  list(
    original_lmer = entry$model,
    data = d,
    formula = model_formula,
    site_year = nonspatial$model,
    spatial_field = spatial$model,
    site_year_elapsed_seconds = nonspatial$elapsed_seconds,
    spatial_field_elapsed_seconds = spatial$elapsed_seconds,
    site_year_extra_optimization = nonspatial$extra_optimization,
    spatial_field_extra_optimization = spatial$extra_optimization
  )
}

model_pairs <- purrr::imap(
  model_entries,
  fit_one_model_pair
)

models_rds <- file.path(
  out_dir,
  "plasticity_spatial_field_models.rds"
)

saveRDS(
  list(
    created_at = as.character(Sys.time()),
    settings = list(
      analysis_group = analysis_group,
      models_to_run = names(model_entries),
      mesh_cutoff_km = mesh_cutoff_km,
      n_mesh_vertices = n_mesh_vertices
    ),
    mesh = spatial_mesh,
    models = model_pairs
  ),
  file = models_rds,
  compress = "gzip"
)


# ---------------------------------------------------------------------------- #
# Model fit, convergence, and spatial parameter summaries
# ---------------------------------------------------------------------------- #

safe_sanity <- function(model) {
  checks <- tryCatch(
    sdmTMB::sanity(
      model,
      gradient_thresh = gradient_threshold,
      silent = TRUE
    ),
    error = function(e) NULL
  )

  if (is.null(checks)) {
    return(NA)
  }

  isTRUE(all(unlist(checks)))
}

safe_tidy_random <- function(model) {
  out <- tryCatch(
    sdmTMB::tidy(
      model,
      effects = "ran_pars",
      conf.int = TRUE
    ),
    error = function(e) NULL
  )

  if (is.null(out)) {
    out <- tryCatch(
      sdmTMB::tidy(
        model,
        effects = "ran_par",
        conf.int = TRUE
      ),
      error = function(e) NULL
    )
  }

  if (is.null(out)) {
    return(tibble::tibble())
  }

  tibble::as_tibble(out)
}

model_fit <- purrr::imap_dfr(
  model_pairs,
  function(x, model_name) {
    dplyr::bind_rows(
      tibble::tibble(
        model = model_name,
        model_variant = "site_year",
        n_observations = nrow(x$data),
        n_sites = dplyr::n_distinct(x$data$SITE_ID),
        n_site_years = dplyr::n_distinct(x$data$site_year_id),
        n_mesh_vertices = n_mesh_vertices,
        logLik = as.numeric(stats::logLik(x$site_year)),
        AIC = stats::AIC(x$site_year),
        max_abs_gradient = max_abs_gradient(x$site_year),
        positive_definite_Hessian = pd_hessian(x$site_year),
        all_sanity_checks_pass = safe_sanity(x$site_year),
        elapsed_minutes = x$site_year_elapsed_seconds / 60,
        extra_optimization = x$site_year_extra_optimization
      ),
      tibble::tibble(
        model = model_name,
        model_variant = "spatial_field",
        n_observations = nrow(x$data),
        n_sites = dplyr::n_distinct(x$data$SITE_ID),
        n_site_years = dplyr::n_distinct(x$data$site_year_id),
        n_mesh_vertices = n_mesh_vertices,
        logLik = as.numeric(stats::logLik(x$spatial_field)),
        AIC = stats::AIC(x$spatial_field),
        max_abs_gradient = max_abs_gradient(x$spatial_field),
        positive_definite_Hessian = pd_hessian(x$spatial_field),
        all_sanity_checks_pass = safe_sanity(x$spatial_field),
        elapsed_minutes = x$spatial_field_elapsed_seconds / 60,
        extra_optimization = x$spatial_field_extra_optimization
      )
    )
  }
) |>
  dplyr::group_by(model) |>
  dplyr::mutate(
    delta_AIC_from_best = AIC - min(AIC),
    delta_AIC_spatial_minus_site_year =
      AIC[model_variant == "spatial_field"] -
      AIC[model_variant == "site_year"]
  ) |>
  dplyr::ungroup()

spatial_parameters <- purrr::imap_dfr(
  model_pairs,
  function(x, model_name) {
    safe_tidy_random(x$spatial_field) |>
      dplyr::mutate(
        model = model_name,
        .before = 1
      )
  }
)


# ---------------------------------------------------------------------------- #
# Compare focal coefficients within the same sdmTMB engine
# ---------------------------------------------------------------------------- #

safe_tidy_fixed <- function(model) {
  out <- tryCatch(
    sdmTMB::tidy(
      model,
      effects = "fixed",
      conf.int = TRUE
    ),
    error = function(e) NULL
  )

  if (is.null(out)) {
    out <- tryCatch(
      sdmTMB::tidy(
        model,
        conf.int = TRUE
      ),
      error = function(e) NULL
    )
  }

  if (is.null(out)) {
    stop("Could not extract fixed-effect coefficients from an sdmTMB model.")
  }

  out <- tibble::as_tibble(out)

  standard_names <- c(
    "std.error" = "SE",
    "std_error" = "SE",
    "p.value" = "p",
    "p_value" = "p",
    "conf.low" = "lower_95",
    "conf_low" = "lower_95",
    "conf.high" = "upper_95",
    "conf_high" = "upper_95"
  )

  for (old_name in names(standard_names)) {
    if (old_name %in% names(out)) {
      names(out)[names(out) == old_name] <- standard_names[[old_name]]
    }
  }

  required_columns <- c("term", "estimate", "SE")
  missing_columns <- setdiff(required_columns, names(out))

  if (length(missing_columns) > 0) {
    stop(
      "Unexpected sdmTMB tidy output; missing: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  if (!"p" %in% names(out)) {
    out$p <- NA_real_
  }

  if (!"lower_95" %in% names(out)) {
    out$lower_95 <- out$estimate - 1.96 * out$SE
  }

  if (!"upper_95" %in% names(out)) {
    out$upper_95 <- out$estimate + 1.96 * out$SE
  }

  out |>
    dplyr::select(
      term,
      estimate,
      SE,
      lower_95,
      upper_95,
      p
    )
}

focal_coefficients <- purrr::imap_dfr(
  model_pairs,
  function(x, model_name) {
    entry <- model_entries[[model_name]]
    anomaly <- paste0("clim_anomaly_tw", entry$best_window)

    focal_terms <- names(lme4::fixef(entry$model))
    focal_terms <- focal_terms[
      focal_terms == anomaly |
        (
          grepl(anomaly, focal_terms, fixed = TRUE) &
            grepl(":", focal_terms, fixed = TRUE)
        )
    ]

    site_year_coefficients <- safe_tidy_fixed(x$site_year) |>
      dplyr::filter(term %in% focal_terms) |>
      dplyr::rename_with(
        ~ paste0(.x, "_site_year"),
        -term
      )

    spatial_coefficients <- safe_tidy_fixed(x$spatial_field) |>
      dplyr::filter(term %in% focal_terms) |>
      dplyr::rename_with(
        ~ paste0(.x, "_spatial"),
        -term
      )

    site_year_coefficients |>
      dplyr::full_join(
        spatial_coefficients,
        by = "term"
      ) |>
      dplyr::mutate(
        model = model_name,
        estimate_change =
          estimate_spatial - estimate_site_year,
        relative_estimate_change = dplyr::if_else(
          estimate_site_year == 0,
          NA_real_,
          estimate_change / abs(estimate_site_year)
        ),
        SE_ratio = SE_spatial / SE_site_year,
        sign_retained =
          sign(estimate_spatial) == sign(estimate_site_year),
        significant_site_year = p_site_year < 0.05,
        significant_spatial = p_spatial < 0.05,
        .before = 1
      )
  }
)


# ---------------------------------------------------------------------------- #
# Conditional residuals
# ---------------------------------------------------------------------------- #

extract_sdm_residuals <- function(
    fitted_model,
    d,
    model_name,
    model_variant) {

  raw_residual <- stats::residuals(
    fitted_model,
    type = "response"
  )

  if (length(raw_residual) != nrow(d)) {
    stop(
      model_name,
      " [",
      model_variant,
      "]: residual and data lengths differ."
    )
  }

  tibble::tibble(
    model = model_name,
    model_variant = model_variant,
    bms_id = as.character(d$bms_id),
    YEAR = d$YEAR,
    SITE_ID = as.character(d$SITE_ID),
    SPECIES = as.character(d$SPECIES),
    x_3035 = d$x_3035,
    y_3035 = d$y_3035,
    residual = as.numeric(raw_residual)
  )
}

extract_lmer_residuals <- function(
    fitted_model,
    d,
    model_name) {

  raw_residual <- stats::residuals(fitted_model)

  if (length(raw_residual) != nrow(d)) {
    stop(
      model_name,
      " [original_lmer]: residual and data lengths differ."
    )
  }

  tibble::tibble(
    model = model_name,
    model_variant = "original_lmer",
    bms_id = as.character(d$bms_id),
    YEAR = d$YEAR,
    SITE_ID = as.character(d$SITE_ID),
    SPECIES = as.character(d$SPECIES),
    x_3035 = d$x_3035,
    y_3035 = d$y_3035,
    residual = as.numeric(raw_residual)
  )
}

residual_data <- purrr::imap_dfr(
  model_pairs,
  function(x, model_name) {
    dplyr::bind_rows(
      extract_lmer_residuals(
        fitted_model = x$original_lmer,
        d = x$data,
        model_name = model_name
      ),
      extract_sdm_residuals(
        fitted_model = x$site_year,
        d = x$data,
        model_name = model_name,
        model_variant = "site_year"
      ),
      extract_sdm_residuals(
        fitted_model = x$spatial_field,
        d = x$data,
        model_name = model_name,
        model_variant = "spatial_field"
      )
    )
  }
)

site_year_residuals <- residual_data |>
  dplyr::group_by(
    model,
    model_variant,
    bms_id,
    YEAR,
    SITE_ID,
    x_3035,
    y_3035
  ) |>
  dplyr::summarise(
    residual = mean(residual, na.rm = TRUE),
    n_species = dplyr::n_distinct(SPECIES),
    .groups = "drop"
  ) |>
  dplyr::filter(
    is.finite(x_3035),
    is.finite(y_3035),
    is.finite(residual)
  )


# ---------------------------------------------------------------------------- #
# Moran's I within every BMS-year and distance band
# ---------------------------------------------------------------------------- #

empty_moran_result <- function(
    status,
    n_total,
    n_used = NA_integer_) {

  tibble::tibble(
    n_sites_total = as.integer(n_total),
    n_sites_used = as.integer(n_used),
    n_links = NA_integer_,
    expected_I = NA_real_,
    moran_I = NA_real_,
    excess_I = NA_real_,
    p_permutation = NA_real_,
    status = status
  )
}

moran_distance_one <- function(
    d,
    d_min_km,
    d_max_km,
    minimum_sites,
    nsim) {

  d <- d |>
    dplyr::arrange(SITE_ID)

  n_total <- nrow(d)

  if (n_total < minimum_sites) {
    return(empty_moran_result("too_few_sites", n_total))
  }

  coordinates <- as.matrix(d[, c("x_3035", "y_3035")])

  tryCatch({
    initial_neighbours <- spdep::dnearneigh(
      coordinates,
      d1 = d_min_km * 1000,
      d2 = d_max_km * 1000,
      row.names = seq_len(n_total),
      longlat = FALSE
    )

    connected <- spdep::card(initial_neighbours) > 0
    n_used <- sum(connected)

    if (n_used < minimum_sites) {
      return(empty_moran_result(
        "too_few_connected_sites",
        n_total,
        n_used
      ))
    }

    d_used <- d[connected, , drop = FALSE]
    coordinates_used <- coordinates[connected, , drop = FALSE]

    neighbours_used <- spdep::dnearneigh(
      coordinates_used,
      d1 = d_min_km * 1000,
      d2 = d_max_km * 1000,
      row.names = seq_len(n_used),
      longlat = FALSE
    )

    n_links <- as.integer(sum(spdep::card(neighbours_used)) / 2)

    if (n_links < minimum_sites) {
      return(empty_moran_result(
        "too_few_links",
        n_total,
        n_used
      ))
    }

    if (
      !is.finite(stats::sd(d_used$residual)) ||
        stats::sd(d_used$residual) == 0
    ) {
      return(empty_moran_result(
        "zero_residual_variance",
        n_total,
        n_used
      ))
    }

    weights <- spdep::nb2listw(
      neighbours_used,
      style = "W",
      zero.policy = TRUE
    )

    asymptotic_test <- spdep::moran.test(
      d_used$residual,
      listw = weights,
      randomisation = TRUE,
      alternative = "greater",
      zero.policy = TRUE
    )

    permutation_test <- spdep::moran.mc(
      d_used$residual,
      listw = weights,
      nsim = nsim,
      alternative = "greater",
      zero.policy = TRUE
    )

    observed_I <- unname(
      as.numeric(
        asymptotic_test$estimate[["Moran I statistic"]]
      )
    )

    expected_I <- -1 / (n_used - 1)

    tibble::tibble(
      n_sites_total = as.integer(n_total),
      n_sites_used = as.integer(n_used),
      n_links = as.integer(n_links),
      expected_I = expected_I,
      moran_I = observed_I,
      excess_I = observed_I - expected_I,
      p_permutation = permutation_test$p.value,
      status = "OK"
    )
  }, error = function(e) {
    empty_moran_result(
      paste0("ERROR: ", conditionMessage(e)),
      n_total
    )
  })
}

nested_groups <- site_year_residuals |>
  dplyr::group_by(
    model,
    model_variant,
    bms_id,
    YEAR
  ) |>
  tidyr::nest() |>
  dplyr::ungroup() |>
  dplyr::mutate(group_index = dplyr::row_number())

test_grid <- dplyr::cross_join(
  tibble::tibble(group_index = nested_groups$group_index),
  distance_bands
)

message(
  "Running ",
  nrow(test_grid),
  " Moran permutation tests..."
)

moran_by_bms_year <- purrr::pmap_dfr(
  test_grid,
  function(group_index, distance_band, d_min_km, d_max_km) {
    group <- nested_groups[group_index, ]

    dplyr::bind_cols(
      group |>
        dplyr::select(
          model,
          model_variant,
          bms_id,
          YEAR
        ),
      tibble::tibble(
        distance_band = distance_band,
        d_min_km = d_min_km,
        d_max_km = d_max_km
      ),
      moran_distance_one(
        d = group$data[[1]],
        d_min_km = d_min_km,
        d_max_km = d_max_km,
        minimum_sites = minimum_connected_sites,
        nsim = n_permutations
      )
    )
  }
) |>
  dplyr::group_by(
    model,
    model_variant,
    distance_band
  ) |>
  dplyr::mutate(
    p_FDR = stats::p.adjust(
      p_permutation,
      method = "BH"
    )
  ) |>
  dplyr::ungroup()


# ---------------------------------------------------------------------------- #
# Summarise the residual spatial change
# ---------------------------------------------------------------------------- #

moran_summary <- moran_by_bms_year |>
  dplyr::filter(status == "OK") |>
  dplyr::group_by(
    model,
    model_variant,
    distance_band
  ) |>
  dplyr::summarise(
    n_tests = dplyr::n(),
    mean_moran_I = mean(moran_I),
    median_moran_I = stats::median(moran_I),
    mean_excess_I = mean(excess_I),
    median_excess_I = stats::median(excess_I),
    proportion_positive = mean(excess_I > 0),
    n_positive_raw_p_lt_0_05 =
      sum(excess_I > 0 & p_permutation < 0.05),
    n_positive_FDR_lt_0_05 =
      sum(excess_I > 0 & p_FDR < 0.05),
    .groups = "drop"
  )

paired_moran <- moran_by_bms_year |>
  dplyr::filter(
    status == "OK",
    model_variant %in% c(
      "site_year",
      "spatial_field"
    )
  ) |>
  dplyr::select(
    model,
    bms_id,
    YEAR,
    distance_band,
    model_variant,
    n_sites_used,
    moran_I,
    excess_I,
    p_permutation,
    p_FDR
  ) |>
  tidyr::pivot_wider(
    names_from = model_variant,
    values_from = c(
      n_sites_used,
      moran_I,
      excess_I,
      p_permutation,
      p_FDR
    )
  ) |>
  dplyr::filter(
    is.finite(excess_I_site_year),
    is.finite(excess_I_spatial_field)
  ) |>
  dplyr::mutate(
    change_in_excess_I =
      excess_I_spatial_field - excess_I_site_year,
    reduction_in_excess_I =
      excess_I_site_year - excess_I_spatial_field,
    spatial_signal_reduced =
      excess_I_spatial_field < excess_I_site_year
  )

paired_summary <- paired_moran |>
  dplyr::group_by(model, distance_band) |>
  dplyr::summarise(
    n_paired_tests = dplyr::n(),
    median_excess_I_site_year =
      stats::median(excess_I_site_year),
    median_excess_I_spatial_field =
      stats::median(excess_I_spatial_field),
    median_reduction_in_excess_I =
      stats::median(reduction_in_excess_I),
    mean_reduction_in_excess_I =
      mean(reduction_in_excess_I),
    proportion_tests_reduced =
      mean(spatial_signal_reduced),
    n_FDR_significant_site_year =
      sum(
        excess_I_site_year > 0 &
          p_FDR_site_year < 0.05
      ),
    n_FDR_significant_spatial_field =
      sum(
        excess_I_spatial_field > 0 &
          p_FDR_spatial_field < 0.05
      ),
    .groups = "drop"
  )


# ---------------------------------------------------------------------------- #
# Figure and Excel export
# ---------------------------------------------------------------------------- #

plot_summary <- moran_by_bms_year |>
  dplyr::filter(status == "OK") |>
  dplyr::group_by(
    model,
    model_variant,
    distance_band
  ) |>
  dplyr::summarise(
    median_excess_I = stats::median(excess_I),
    q25 = stats::quantile(excess_I, 0.25),
    q75 = stats::quantile(excess_I, 0.75),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    model_variant = factor(
      model_variant,
      levels = c(
        "original_lmer",
        "site_year",
        "spatial_field"
      ),
      labels = c(
        "Original lmer",
        "Site x year",
        "Spatial field"
      )
    )
  )

p_comparison <- ggplot2::ggplot(
  plot_summary,
  ggplot2::aes(
    x = model_variant,
    y = median_excess_I,
    group = distance_band,
    colour = distance_band
  )
) +
  ggplot2::geom_hline(
    yintercept = 0,
    linetype = 2,
    colour = "grey55"
  ) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(size = 2.2) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = q25, ymax = q75),
    width = 0.08
  ) +
  ggplot2::facet_wrap(~ model, scales = "free_y") +
  ggplot2::labs(
    x = NULL,
    y = "Median excess Moran's I",
    colour = "Distance band"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    legend.position = "bottom",
    axis.text.x = ggplot2::element_text(
      angle = 15,
      hjust = 1
    )
  )

figure_path <- file.path(
  out_dir,
  "spatial_field_moran_comparison.png"
)

ggplot2::ggsave(
  filename = figure_path,
  plot = p_comparison,
  width = 10,
  height = 7,
  dpi = 300
)

mesh_summary <- tibble::tibble(
  mesh_cutoff_km = mesh_cutoff_km,
  n_unique_coordinates = nrow(mesh_data),
  n_mesh_vertices = n_mesh_vertices,
  spatial_structure =
    "Independent Matérn spatial field per year; shared range",
  range_interpretation =
    "Distance in km where correlation is approximately 0.13"
)

output_xlsx <- file.path(
  out_dir,
  "plasticity_spatial_field_comparison.xlsx"
)

writexl::write_xlsx(
  list(
    paired_summary = paired_summary,
    model_fit = model_fit,
    focal_coefficients = focal_coefficients,
    spatial_parameters = spatial_parameters,
    moran_summary = moran_summary,
    paired_BMS_year = paired_moran,
    moran_all_tests = moran_by_bms_year,
    mesh_summary = mesh_summary
  ),
  path = output_xlsx
)

message("Saved comparison tables: ", output_xlsx)
message("Saved fitted models: ", models_rds)
message("Saved comparison figure: ", figure_path)
message(
  "Key sheets: paired_summary, model_fit, focal_coefficients, ",
  "and spatial_parameters."
)
