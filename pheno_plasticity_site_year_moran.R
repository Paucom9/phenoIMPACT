# ============================================================================================ #
# pheno_plasticity_site_year_moran.R
#
# Directly compare residual spatial autocorrelation between:
#   1. the original phenology-plasticity models; and
#   2. the same models plus a SITE_ID x YEAR random intercept.
#
# Requires:
#   output/phenology_plasticity/plasticity_main_models.rds
#
# Outputs:
#   output/phenology_plasticity/diagnostics/site_year_comparison/
#     plasticity_site_year_models.rds
#     plasticity_site_year_moran_comparison.xlsx
#     site_year_moran_comparison.png
#
# Coordinates saved in the model bundle are EPSG:3035 metres.
# ============================================================================================ #


# ---------------------------------------------------------------------------- #
# Settings
# ---------------------------------------------------------------------------- #

analysis_group <- "primary"

# These are the distances over which the previous correlogram showed positive
# residual spatial autocorrelation.
distance_bands <- tibble::tribble(
  ~distance_band, ~d_min_km, ~d_max_km,
  "0-25 km",               0,         25,
  "25-50 km",             25,         50
)

n_permutations <- 499L
minimum_connected_sites <- 10L

set.seed(20260724)


# ---------------------------------------------------------------------------- #
# Packages and paths
# ---------------------------------------------------------------------------- #

required_packages <- c(
  "dplyr", "tidyr", "purrr", "tibble", "ggplot2",
  "lme4", "lmerTest", "spdep", "writexl", "here"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Install the missing packages first: ",
    paste(missing_packages, collapse = ", ")
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
  "site_year_comparison"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

model_bundle <- readRDS(rds_path)


# ---------------------------------------------------------------------------- #
# Recover the exact rows used by each original model
# ---------------------------------------------------------------------------- #

match_entry_data <- function(entry) {
  fitted_frame <- stats::model.frame(entry$model)
  original_data <- entry$data

  row_index <- match(rownames(fitted_frame), rownames(original_data))

  if (anyNA(row_index)) {
    if (nrow(fitted_frame) == nrow(original_data)) {
      warning("Row names did not match; using row order because row counts agree.")
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

names(model_entries) <- paste(analysis_group, names(model_entries), sep = "__")

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


# ---------------------------------------------------------------------------- #
# Fit the site-year models
# ---------------------------------------------------------------------------- #

model_convergence_message <- function(fitted_model) {
  messages <- fitted_model@optinfo$conv$lme4$messages
  if (is.null(messages)) "OK" else paste(messages, collapse = "; ")
}

coefficient_table <- function(fitted_model) {
  coefficient_matrix <- as.data.frame(summary(fitted_model)$coefficients)
  p_column <- grep("^Pr\\(", names(coefficient_matrix), value = TRUE)
  df_column <- intersect("df", names(coefficient_matrix))

  tibble::tibble(
    term = rownames(coefficient_matrix),
    estimate = coefficient_matrix[["Estimate"]],
    SE = coefficient_matrix[["Std. Error"]],
    df = if (length(df_column) == 1) {
      coefficient_matrix[[df_column]]
    } else {
      NA_real_
    },
    p = if (length(p_column) >= 1) {
      coefficient_matrix[[p_column[1]]]
    } else {
      NA_real_
    }
  )
}

fit_site_year_model <- function(entry, model_name) {
  original_model <- entry$model

  d <- match_entry_data(entry) |>
    dplyr::mutate(
      SITE_ID = as.character(SITE_ID),
      SPECIES = as.character(SPECIES),
      YEAR = as.integer(YEAR),
      bms_id = as.character(bms_id),
      site_year_id = interaction(
        SITE_ID,
        YEAR,
        drop = TRUE,
        lex.order = TRUE
      )
    ) |>
    dplyr::left_join(
      site_coordinates,
      by = c("SITE_ID", "bms_id")
    )

  if (any(!is.finite(d$x_3035)) || any(!is.finite(d$y_3035))) {
    stop(model_name, ": some fitted observations lack EPSG:3035 coordinates.")
  }

  site_year_formula <- stats::update(
    stats::formula(original_model),
    . ~ . + (1 | site_year_id)
  )

  message("Fitting site-year model: ", model_name)

  site_year_model <- lmerTest::lmer(
    site_year_formula,
    data = d,
    REML = lme4::isREML(original_model),
    control = lme4::lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 2e6)
    )
  )

  list(
    original_model = original_model,
    site_year_model = site_year_model,
    data = d
  )
}

refitted_models <- purrr::imap(
  model_entries,
  ~ fit_site_year_model(.x, .y)
)

saveRDS(
  refitted_models,
  file = file.path(out_dir, "plasticity_site_year_models.rds"),
  compress = "gzip"
)


# ---------------------------------------------------------------------------- #
# Compare model fit and focal coefficients
# ---------------------------------------------------------------------------- #

extract_site_year_sd <- function(fitted_model) {
  variance_components <- as.data.frame(lme4::VarCorr(fitted_model))
  value <- variance_components$sdcor[
    match("site_year_id", variance_components$grp)
  ]

  if (length(value) == 0) NA_real_ else value
}

fit_comparison <- purrr::imap_dfr(
  refitted_models,
  function(x, model_name) {
    tibble::tibble(
      model = model_name,
      n_observations = stats::nobs(x$original_model),
      original_AIC = stats::AIC(x$original_model),
      site_year_AIC = stats::AIC(x$site_year_model),
      delta_AIC_site_year_minus_original =
        stats::AIC(x$site_year_model) - stats::AIC(x$original_model),
      original_sigma = stats::sigma(x$original_model),
      site_year_sigma = stats::sigma(x$site_year_model),
      site_year_SD = extract_site_year_sd(x$site_year_model),
      site_year_singular =
        lme4::isSingular(x$site_year_model, tol = 1e-4),
      site_year_convergence =
        model_convergence_message(x$site_year_model)
    )
  }
)

coefficient_comparison <- purrr::imap_dfr(
  refitted_models,
  function(x, model_name) {
    anomaly <- paste0(
      "clim_anomaly_tw",
      model_entries[[model_name]]$best_window
    )

    focal_terms <- names(lme4::fixef(x$original_model))
    focal_terms <- focal_terms[
      focal_terms == anomaly |
        (
          grepl(anomaly, focal_terms, fixed = TRUE) &
            grepl(":", focal_terms, fixed = TRUE)
        )
    ]

    original <- coefficient_table(x$original_model) |>
      dplyr::filter(term %in% focal_terms) |>
      dplyr::rename_with(~ paste0(.x, "_original"), -term)

    site_year <- coefficient_table(x$site_year_model) |>
      dplyr::filter(term %in% focal_terms) |>
      dplyr::rename_with(~ paste0(.x, "_site_year"), -term)

    original |>
      dplyr::left_join(site_year, by = "term") |>
      dplyr::mutate(
        model = model_name,
        estimate_change = estimate_site_year - estimate_original,
        relative_estimate_change = dplyr::if_else(
          estimate_original == 0,
          NA_real_,
          estimate_change / abs(estimate_original)
        ),
        SE_ratio = SE_site_year / SE_original,
        .before = 1
      )
  }
)


# ---------------------------------------------------------------------------- #
# Extract conditional residuals from both model structures
# ---------------------------------------------------------------------------- #

extract_residuals <- function(x, model_name, model_variant) {
  fitted_model <- if (model_variant == "original") {
    x$original_model
  } else {
    x$site_year_model
  }

  raw_residual <- stats::residuals(fitted_model)

  if (length(raw_residual) != nrow(x$data)) {
    stop(model_name, ": residual and data lengths differ.")
  }

  tibble::tibble(
    model = model_name,
    model_variant = model_variant,
    bms_id = x$data$bms_id,
    YEAR = x$data$YEAR,
    SITE_ID = x$data$SITE_ID,
    SPECIES = x$data$SPECIES,
    x_3035 = x$data$x_3035,
    y_3035 = x$data$y_3035,
    scaled_residual =
      as.numeric(raw_residual / stats::sigma(fitted_model))
  )
}

residual_data <- purrr::imap_dfr(
  refitted_models,
  function(x, model_name) {
    dplyr::bind_rows(
      extract_residuals(x, model_name, "original"),
      extract_residuals(x, model_name, "site_year")
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
    scaled_residual = mean(scaled_residual, na.rm = TRUE),
    n_species = dplyr::n_distinct(SPECIES),
    .groups = "drop"
  ) |>
  dplyr::filter(
    is.finite(x_3035),
    is.finite(y_3035),
    is.finite(scaled_residual)
  )


# ---------------------------------------------------------------------------- #
# Moran's I within each BMS-year and distance band
# ---------------------------------------------------------------------------- #

empty_moran_result <- function(status, n_total, n_used = NA_integer_) {
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

  result <- tryCatch({
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
      return(empty_moran_result("too_few_links", n_total, n_used))
    }

    if (!is.finite(stats::sd(d_used$scaled_residual)) ||
        stats::sd(d_used$scaled_residual) == 0) {
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
      d_used$scaled_residual,
      listw = weights,
      randomisation = TRUE,
      alternative = "greater",
      zero.policy = TRUE
    )

    permutation_test <- spdep::moran.mc(
      d_used$scaled_residual,
      listw = weights,
      nsim = nsim,
      alternative = "greater",
      zero.policy = TRUE
    )

    observed_I <- unname(
      as.numeric(asymptotic_test$estimate[["Moran I statistic"]])
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

  result
}

nested_groups <- site_year_residuals |>
  dplyr::group_by(model, model_variant, bms_id, YEAR) |>
  tidyr::nest() |>
  dplyr::ungroup() |>
  dplyr::mutate(group_index = dplyr::row_number())

test_grid <- dplyr::cross_join(
  tibble::tibble(group_index = nested_groups$group_index),
  distance_bands
)

message(
  "Running ", nrow(test_grid),
  " Moran permutation tests..."
)

moran_by_bms_year <- purrr::pmap_dfr(
  test_grid,
  function(group_index, distance_band, d_min_km, d_max_km) {
    group <- nested_groups[group_index, ]

    dplyr::bind_cols(
      group |>
        dplyr::select(model, model_variant, bms_id, YEAR),
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
  dplyr::group_by(model, model_variant, distance_band) |>
  dplyr::mutate(
    p_FDR = stats::p.adjust(p_permutation, method = "BH")
  ) |>
  dplyr::ungroup()


# ---------------------------------------------------------------------------- #
# Summarise whether site-year reduces the spatial signal
# ---------------------------------------------------------------------------- #

moran_summary <- moran_by_bms_year |>
  dplyr::filter(status == "OK") |>
  dplyr::group_by(model, model_variant, distance_band) |>
  dplyr::summarise(
    n_tests = dplyr::n(),
    mean_moran_I = mean(moran_I),
    median_moran_I = stats::median(moran_I),
    mean_excess_I = mean(excess_I),
    median_excess_I = stats::median(excess_I),
    proportion_positive = mean(excess_I > 0),
    n_positive_raw_p_lt_0.05 =
      sum(excess_I > 0 & p_permutation < 0.05),
    n_positive_FDR_lt_0.05 =
      sum(excess_I > 0 & p_FDR < 0.05),
    .groups = "drop"
  )

paired_moran <- moran_by_bms_year |>
  dplyr::filter(status == "OK") |>
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
    is.finite(excess_I_original),
    is.finite(excess_I_site_year)
  ) |>
  dplyr::mutate(
    change_in_excess_I =
      excess_I_site_year - excess_I_original,
    reduction_in_excess_I =
      excess_I_original - excess_I_site_year,
    spatial_signal_reduced =
      excess_I_site_year < excess_I_original
  )

paired_summary <- paired_moran |>
  dplyr::group_by(model, distance_band) |>
  dplyr::summarise(
    n_paired_tests = dplyr::n(),
    median_excess_I_original =
      stats::median(excess_I_original),
    median_excess_I_site_year =
      stats::median(excess_I_site_year),
    median_reduction_in_excess_I =
      stats::median(reduction_in_excess_I),
    mean_reduction_in_excess_I =
      mean(reduction_in_excess_I),
    proportion_tests_reduced =
      mean(spatial_signal_reduced),
    n_FDR_significant_original =
      sum(excess_I_original > 0 & p_FDR_original < 0.05),
    n_FDR_significant_site_year =
      sum(excess_I_site_year > 0 & p_FDR_site_year < 0.05),
    .groups = "drop"
  )


# ---------------------------------------------------------------------------- #
# Plot and export
# ---------------------------------------------------------------------------- #

plot_summary <- moran_by_bms_year |>
  dplyr::filter(status == "OK") |>
  dplyr::group_by(model, model_variant, distance_band) |>
  dplyr::summarise(
    median_excess_I = stats::median(excess_I),
    q25 = stats::quantile(excess_I, 0.25),
    q75 = stats::quantile(excess_I, 0.75),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    model_variant = factor(
      model_variant,
      levels = c("original", "site_year"),
      labels = c("Original", "Site x year")
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
    axis.text.x = ggplot2::element_text(angle = 15, hjust = 1)
  )

ggplot2::ggsave(
  filename = file.path(out_dir, "site_year_moran_comparison.png"),
  plot = p_comparison,
  width = 10,
  height = 7,
  dpi = 300
)

output_xlsx <- file.path(
  out_dir,
  "plasticity_site_year_moran_comparison.xlsx"
)

writexl::write_xlsx(
  list(
    paired_summary = paired_summary,
    moran_summary = moran_summary,
    paired_BMS_year = paired_moran,
    moran_all_tests = moran_by_bms_year,
    model_fit = fit_comparison,
    focal_coefficients = coefficient_comparison
  ),
  path = output_xlsx
)

message("Saved comparison tables: ", output_xlsx)
message(
  "Saved fitted site-year models: ",
  file.path(out_dir, "plasticity_site_year_models.rds")
)
message(
  "Saved comparison figure: ",
  file.path(out_dir, "site_year_moran_comparison.png")
)
