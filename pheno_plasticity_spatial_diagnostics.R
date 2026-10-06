# ============================================================================================ #
# pheno_plasticity_spatial_diagnostics.R
#
# Deeper spatial diagnostics and sensitivity analyses for the selected
# phenological temperature-sensitivity models.
#
# Requires:
#   output/phenology_plasticity/plasticity_main_models.rds
#
# Main analyses:
#   1. Moran's I within BMS-year using several k-nearest-neighbour graphs.
#   2. Distance-band correlograms within BMS-year.
#   3. Optional Moran's I within BMS-species-year.
#   4. Optional coefficient sensitivity to site-year and spatial-block-year
#      random intercepts.
#
# Coordinates in plasticity_main_models.rds are already EPSG:3035 metres.
# ============================================================================================ #


# ---------------------------------------------------------------------------- #
# Settings
# ---------------------------------------------------------------------------- #

analysis_groups <- "primary"

# Neighbour definitions used to check that the result is not specific to k = 5.
k_values <- c(3L, 5L, 10L)

# Absolute distance bands for the correlogram.
distance_bands <- tibble::tribble(
  ~distance_band, ~d_min_km, ~d_max_km,
  "0-25 km",               0,         25,
  "25-50 km",             25,         50,
  "50-100 km",            50,        100,
  "100-200 km",          100,        200,
  "200-400 km",          200,        400
)

n_permutations <- 499L
minimum_sites <- 10L

# This is informative but potentially slow because it creates many
# BMS-species-year tests. Start with FALSE.
run_species_year_tests <- FALSE
n_permutations_species <- 199L

# This refits the four primary models several times. Run the diagnostics first,
# then turn it on. Ideally choose block sizes slightly larger than the distance
# over which the correlogram remains positive.
run_spatial_block_sensitivity <- TRUE
block_sizes_km <- c(50, 100)
block_origin_shifts <- c(0, 0.5)

set.seed(20260723)


# ---------------------------------------------------------------------------- #
# Packages and paths
# ---------------------------------------------------------------------------- #

required_packages <- c(
  "dplyr", "tidyr", "purrr", "tibble", "ggplot2",
  "lme4", "lmerTest", "spdep", "writexl", "here", "scales"
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

model_bundle <- readRDS(rds_path)

out_dir <- here::here(
  "output",
  "phenology_plasticity",
  "diagnostics",
  "spatial_deep"
)

figure_dir <- file.path(out_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)


# ---------------------------------------------------------------------------- #
# Recover the exact fitted rows and their conditional residuals
# ---------------------------------------------------------------------------- #

flatten_model_entries <- function(bundle, groups) {
  entries <- list()

  for (group_name in intersect(groups, names(bundle))) {
    for (model_name in names(bundle[[group_name]])) {
      entries[[paste(group_name, model_name, sep = "__")]] <-
        bundle[[group_name]][[model_name]]
    }
  }

  entries
}

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

extract_residual_data <- function(entry, model_name) {
  fitted_model <- entry$model
  d <- match_entry_data(entry)

  required <- c("YEAR", "SITE_ID", "SPECIES", "bms_id")
  missing <- setdiff(required, names(d))

  if (length(missing) > 0) {
    stop(
      model_name, " lacks required columns: ",
      paste(missing, collapse = ", ")
    )
  }

  raw_residual <- stats::residuals(fitted_model)

  if (length(raw_residual) != nrow(d)) {
    stop(model_name, ": residuals and matched model rows have different lengths.")
  }

  tibble::tibble(
    model = model_name,
    YEAR = as.integer(d$YEAR),
    SITE_ID = as.character(d$SITE_ID),
    SPECIES = as.character(d$SPECIES),
    bms_id = as.character(d$bms_id),
    scaled_residual = as.numeric(raw_residual / stats::sigma(fitted_model))
  )
}

model_entries <- flatten_model_entries(model_bundle, analysis_groups)

if (length(model_entries) == 0) {
  stop("No model entries found for: ", paste(analysis_groups, collapse = ", "))
}

residual_data <- purrr::imap_dfr(
  model_entries,
  ~ extract_residual_data(.x, .y)
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

# Geographic longitude/latitude would normally be below 180/90. Values in the
# EPSG:3035 coordinate table should instead be in metres and much larger.
if (
  max(abs(site_coordinates$x_3035), na.rm = TRUE) < 1000 ||
  max(abs(site_coordinates$y_3035), na.rm = TRUE) < 1000
) {
  stop(
    "Coordinates do not look like EPSG:3035 metres. ",
    "Check how site_coordinates was created."
  )
}


# ---------------------------------------------------------------------------- #
# Datasets for the two spatial questions
# ---------------------------------------------------------------------------- #

# Common spatial signal in each BMS-year, averaging across species at each site.
site_year_residuals <- residual_data |>
  dplyr::group_by(model, bms_id, YEAR, SITE_ID) |>
  dplyr::summarise(
    scaled_residual = mean(scaled_residual, na.rm = TRUE),
    n_species = dplyr::n_distinct(SPECIES),
    .groups = "drop"
  ) |>
  dplyr::left_join(
    site_coordinates,
    by = c("SITE_ID", "bms_id")
  ) |>
  dplyr::filter(
    is.finite(scaled_residual),
    is.finite(x_3035),
    is.finite(y_3035)
  )

# Species-specific spatial signal in each BMS-year.
species_year_residuals <- residual_data |>
  dplyr::group_by(model, bms_id, YEAR, SPECIES, SITE_ID) |>
  dplyr::summarise(
    scaled_residual = mean(scaled_residual, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::left_join(
    site_coordinates,
    by = c("SITE_ID", "bms_id")
  ) |>
  dplyr::filter(
    is.finite(scaled_residual),
    is.finite(x_3035),
    is.finite(y_3035)
  )


# ---------------------------------------------------------------------------- #
# Moran's I helpers
# ---------------------------------------------------------------------------- #

empty_moran_result <- function(
    status,
    n_sites_total = NA_integer_,
    n_sites_used = NA_integer_) {

  tibble::tibble(
    n_sites_total = as.integer(n_sites_total),
    n_sites_used = as.integer(n_sites_used),
    n_links = NA_integer_,
    n_components = NA_integer_,
    median_link_distance_km = NA_real_,
    moran_I = NA_real_,
    expected_I = NA_real_,
    excess_I = NA_real_,
    p_asymptotic = NA_real_,
    p_permutation = NA_real_,
    status = status
  )
}

collapse_coincident_sites <- function(d) {
  d |>
    dplyr::group_by(x_3035, y_3035) |>
    dplyr::summarise(
      scaled_residual = mean(scaled_residual, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::filter(is.finite(scaled_residual))
}

moran_knn_one <- function(
    d,
    k,
    nsim = 499L,
    minimum_sites = 10L) {

  d <- collapse_coincident_sites(d)
  n_sites <- nrow(d)

  if (n_sites < minimum_sites) {
    return(empty_moran_result("too_few_sites", n_sites, n_sites))
  }

  if (n_sites <= k) {
    return(empty_moran_result("k_not_smaller_than_n", n_sites, n_sites))
  }

  if (!is.finite(stats::sd(d$scaled_residual)) ||
      stats::sd(d$scaled_residual) == 0) {
    return(empty_moran_result("zero_residual_variance", n_sites, n_sites))
  }

  coordinates <- as.matrix(d[, c("x_3035", "y_3035")])

  tryCatch({
    neighbours <- spdep::knn2nb(
      spdep::knearneigh(coordinates, k = as.integer(k), longlat = FALSE),
      sym = TRUE
    )

    weights <- spdep::nb2listw(
      neighbours,
      style = "W",
      zero.policy = TRUE
    )

    asymptotic_test <- spdep::moran.test(
      d$scaled_residual,
      listw = weights,
      randomisation = TRUE,
      alternative = "greater",
      zero.policy = TRUE
    )

    permutation_test <- spdep::moran.mc(
      d$scaled_residual,
      listw = weights,
      nsim = as.integer(nsim),
      alternative = "greater",
      zero.policy = TRUE
    )

    link_distances <- unlist(
      spdep::nbdists(neighbours, coordinates, longlat = FALSE),
      use.names = FALSE
    ) / 1000

    observed_I <- unname(
      as.numeric(asymptotic_test$estimate[["Moran I statistic"]])
    )
    expected_I <- -1 / (n_sites - 1)

    tibble::tibble(
      n_sites_total = n_sites,
      n_sites_used = n_sites,
      n_links = as.integer(sum(spdep::card(neighbours)) / 2),
      n_components = as.integer(spdep::n.comp.nb(neighbours)$nc),
      median_link_distance_km = stats::median(link_distances),
      moran_I = observed_I,
      expected_I = expected_I,
      excess_I = observed_I - expected_I,
      p_asymptotic = asymptotic_test$p.value,
      p_permutation = permutation_test$p.value,
      status = "OK"
    )
  }, error = function(e) {
    empty_moran_result(
      paste0("error: ", conditionMessage(e)),
      n_sites,
      n_sites
    )
  })
}

moran_distance_one <- function(
    d,
    d_min_km,
    d_max_km,
    nsim = 499L,
    minimum_sites = 10L) {

  d <- collapse_coincident_sites(d)
  n_total <- nrow(d)

  if (n_total < minimum_sites) {
    return(empty_moran_result("too_few_sites", n_total, n_total))
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

    # Remove sites with no neighbour in this particular distance band.
    keep <- spdep::card(initial_neighbours) > 0
    d_used <- d[keep, , drop = FALSE]
    coordinates_used <- coordinates[keep, , drop = FALSE]
    n_used <- nrow(d_used)

    if (n_used < minimum_sites) {
      return(empty_moran_result("too_few_connected_sites", n_total, n_used))
    }

    neighbours <- spdep::dnearneigh(
      coordinates_used,
      d1 = d_min_km * 1000,
      d2 = d_max_km * 1000,
      row.names = seq_len(n_used),
      longlat = FALSE
    )

    n_links <- as.integer(sum(spdep::card(neighbours)) / 2)

    if (n_links < minimum_sites) {
      return(empty_moran_result("too_few_links", n_total, n_used))
    }

    if (!is.finite(stats::sd(d_used$scaled_residual)) ||
        stats::sd(d_used$scaled_residual) == 0) {
      return(empty_moran_result("zero_residual_variance", n_total, n_used))
    }

    weights <- spdep::nb2listw(
      neighbours,
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
      nsim = as.integer(nsim),
      alternative = "greater",
      zero.policy = TRUE
    )

    link_distances <- unlist(
      spdep::nbdists(neighbours, coordinates_used, longlat = FALSE),
      use.names = FALSE
    ) / 1000

    observed_I <- unname(
      as.numeric(asymptotic_test$estimate[["Moran I statistic"]])
    )
    expected_I <- -1 / (n_used - 1)

    tibble::tibble(
      n_sites_total = n_total,
      n_sites_used = n_used,
      n_links = n_links,
      n_components = as.integer(spdep::n.comp.nb(neighbours)$nc),
      median_link_distance_km = stats::median(link_distances),
      moran_I = observed_I,
      expected_I = expected_I,
      excess_I = observed_I - expected_I,
      p_asymptotic = asymptotic_test$p.value,
      p_permutation = permutation_test$p.value,
      status = "OK"
    )
  }, error = function(e) {
    empty_moran_result(
      paste0("error: ", conditionMessage(e)),
      n_total,
      NA_integer_
    )
  })
}

add_multiple_test_corrections <- function(results, family_var) {
  # p_BH_within_family:
  #   BMS-year tests corrected separately for each k or distance band.
  # p_BH_all_families:
  #   additionally treats all k values or all distance bands as tests.
  # BH uses the analytical randomisation p-value. The permutation p-value is
  # retained as a robustness check; with 499 permutations it is too coarse for
  # a large family of multiplicity-adjusted tests.
  results |>
    dplyr::group_by(model, .data[[family_var]]) |>
    dplyr::mutate(
      p_BH_within_family = stats::p.adjust(
        p_asymptotic,
        method = "BH"
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::group_by(model) |>
    dplyr::mutate(
      p_BH_all_families = stats::p.adjust(
        p_asymptotic,
        method = "BH"
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::group_by(model, bms_id, .data[[family_var]]) |>
    dplyr::mutate(
      p_BH_within_BMS = stats::p.adjust(
        p_asymptotic,
        method = "BH"
      )
    ) |>
    dplyr::ungroup()
}


# ---------------------------------------------------------------------------- #
# 1. k-nearest-neighbour sensitivity
# ---------------------------------------------------------------------------- #

run_knn_tests <- function(
    d,
    grouping_variables,
    k_values,
    nsim,
    minimum_sites) {

  nested <- d |>
    dplyr::group_by(dplyr::across(dplyr::all_of(grouping_variables))) |>
    tidyr::nest() |>
    dplyr::ungroup() |>
    dplyr::mutate(k = list(as.integer(k_values))) |>
    tidyr::unnest(k)

  results <- nested |>
    dplyr::mutate(
      result = purrr::map2(
        data,
        k,
        ~ moran_knn_one(
          .x,
          k = .y,
          nsim = nsim,
          minimum_sites = minimum_sites
        )
      )
    ) |>
    dplyr::select(-data) |>
    tidyr::unnest(result)

  add_multiple_test_corrections(results, family_var = "k")
}

message("Running BMS-year Moran tests across k values...")

knn_site_year <- run_knn_tests(
  d = site_year_residuals,
  grouping_variables = c("model", "bms_id", "YEAR"),
  k_values = k_values,
  nsim = n_permutations,
  minimum_sites = minimum_sites
)

knn_summary <- knn_site_year |>
  dplyr::filter(status == "OK") |>
  dplyr::group_by(model, k) |>
  dplyr::summarise(
    n_tests = dplyr::n(),
    mean_moran_I = mean(moran_I),
    median_moran_I = stats::median(moran_I),
    mean_excess_I = mean(excess_I),
    median_excess_I = stats::median(excess_I),
    proportion_positive = mean(excess_I > 0),
    n_asymptotic_p_lt_0_05 = sum(
      excess_I > 0 & p_asymptotic < 0.05,
      na.rm = TRUE
    ),
    n_permutation_p_lt_0_05 = sum(
      excess_I > 0 & p_permutation < 0.05,
      na.rm = TRUE
    ),
    n_BH_within_k_lt_0_05 = sum(
      excess_I > 0 & p_BH_within_family < 0.05,
      na.rm = TRUE
    ),
    n_BH_all_k_lt_0_05 = sum(
      excess_I > 0 & p_BH_all_families < 0.05,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

knn_summary_by_BMS <- knn_site_year |>
  dplyr::filter(status == "OK") |>
  dplyr::group_by(model, k, bms_id) |>
  dplyr::summarise(
    n_tests = dplyr::n(),
    median_n_sites = stats::median(n_sites_used),
    median_link_distance_km = stats::median(median_link_distance_km),
    mean_moran_I = mean(moran_I),
    median_moran_I = stats::median(moran_I),
    mean_excess_I = mean(excess_I),
    proportion_positive = mean(excess_I > 0),
    n_asymptotic_p_lt_0_05 = sum(
      excess_I > 0 & p_asymptotic < 0.05,
      na.rm = TRUE
    ),
    n_permutation_p_lt_0_05 = sum(
      excess_I > 0 & p_permutation < 0.05,
      na.rm = TRUE
    ),
    n_BH_within_BMS_lt_0_05 = sum(
      excess_I > 0 & p_BH_within_BMS < 0.05,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


# ---------------------------------------------------------------------------- #
# 2. Distance-band correlograms
# ---------------------------------------------------------------------------- #

run_distance_tests <- function(
    d,
    grouping_variables,
    bands,
    nsim,
    minimum_sites) {

  nested <- d |>
    dplyr::group_by(dplyr::across(dplyr::all_of(grouping_variables))) |>
    tidyr::nest() |>
    dplyr::ungroup() |>
    dplyr::mutate(band_data = list(bands)) |>
    tidyr::unnest(band_data)

  results <- nested |>
    dplyr::mutate(
      result = purrr::pmap(
        list(data, d_min_km, d_max_km),
        ~ moran_distance_one(
          d = ..1,
          d_min_km = ..2,
          d_max_km = ..3,
          nsim = nsim,
          minimum_sites = minimum_sites
        )
      )
    ) |>
    dplyr::select(-data) |>
    tidyr::unnest(result)

  add_multiple_test_corrections(results, family_var = "distance_band")
}

message("Running distance-band correlograms...")

distance_site_year <- run_distance_tests(
  d = site_year_residuals,
  grouping_variables = c("model", "bms_id", "YEAR"),
  bands = distance_bands,
  nsim = n_permutations,
  minimum_sites = minimum_sites
)

distance_summary <- distance_site_year |>
  dplyr::filter(status == "OK") |>
  dplyr::group_by(model, distance_band, d_min_km, d_max_km) |>
  dplyr::summarise(
    n_tests = dplyr::n(),
    mean_moran_I = mean(moran_I),
    median_moran_I = stats::median(moran_I),
    mean_excess_I = mean(excess_I),
    median_excess_I = stats::median(excess_I),
    q25_excess_I = stats::quantile(excess_I, 0.25),
    q75_excess_I = stats::quantile(excess_I, 0.75),
    proportion_positive = mean(excess_I > 0),
    n_asymptotic_p_lt_0_05 = sum(
      excess_I > 0 & p_asymptotic < 0.05,
      na.rm = TRUE
    ),
    n_permutation_p_lt_0_05 = sum(
      excess_I > 0 & p_permutation < 0.05,
      na.rm = TRUE
    ),
    n_BH_within_band_lt_0_05 = sum(
      excess_I > 0 & p_BH_within_family < 0.05,
      na.rm = TRUE
    ),
    n_BH_all_bands_lt_0_05 = sum(
      excess_I > 0 & p_BH_all_families < 0.05,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) |>
  dplyr::mutate(distance_mid_km = (d_min_km + d_max_km) / 2)

distance_summary_by_BMS <- distance_site_year |>
  dplyr::filter(status == "OK") |>
  dplyr::group_by(
    model,
    bms_id,
    distance_band,
    d_min_km,
    d_max_km
  ) |>
  dplyr::summarise(
    n_tests = dplyr::n(),
    median_n_sites = stats::median(n_sites_used),
    median_n_links = stats::median(n_links),
    mean_moran_I = mean(moran_I),
    median_moran_I = stats::median(moran_I),
    mean_excess_I = mean(excess_I),
    proportion_positive = mean(excess_I > 0),
    n_asymptotic_p_lt_0_05 = sum(
      excess_I > 0 & p_asymptotic < 0.05,
      na.rm = TRUE
    ),
    n_permutation_p_lt_0_05 = sum(
      excess_I > 0 & p_permutation < 0.05,
      na.rm = TRUE
    ),
    n_BH_within_BMS_lt_0_05 = sum(
      excess_I > 0 & p_BH_within_BMS < 0.05,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


# ---------------------------------------------------------------------------- #
# 3. Optional species-specific spatial tests
# ---------------------------------------------------------------------------- #

if (run_species_year_tests) {
  message("Running BMS-species-year Moran tests...")

  knn_species_year <- run_knn_tests(
    d = species_year_residuals,
    grouping_variables = c("model", "bms_id", "SPECIES", "YEAR"),
    k_values = 5L,
    nsim = n_permutations_species,
    minimum_sites = minimum_sites
  )

  knn_species_summary <- knn_species_year |>
    dplyr::filter(status == "OK") |>
    dplyr::group_by(model, bms_id, k) |>
    dplyr::summarise(
      n_species_year_tests = dplyr::n(),
      n_species = dplyr::n_distinct(SPECIES),
      mean_moran_I = mean(moran_I),
      median_moran_I = stats::median(moran_I),
      mean_excess_I = mean(excess_I),
      proportion_positive = mean(excess_I > 0),
      n_asymptotic_p_lt_0_05 = sum(
        excess_I > 0 & p_asymptotic < 0.05,
        na.rm = TRUE
      ),
      n_permutation_p_lt_0_05 = sum(
        excess_I > 0 & p_permutation < 0.05,
        na.rm = TRUE
      ),
      n_BH_model_lt_0_05 = sum(
        excess_I > 0 & p_BH_all_families < 0.05,
        na.rm = TRUE
      ),
      .groups = "drop"
    )
} else {
  knn_species_year <- tibble::tibble()
  knn_species_summary <- tibble::tibble()
}


# ---------------------------------------------------------------------------- #
# Figures
# ---------------------------------------------------------------------------- #

valid_knn <- knn_site_year |>
  dplyr::filter(status == "OK")

p_knn <- ggplot2::ggplot(
  valid_knn,
  ggplot2::aes(
    x = factor(k),
    y = excess_I,
    colour = bms_id
  )
) +
  ggplot2::geom_hline(
    yintercept = 0,
    linetype = 2,
    colour = "grey40"
  ) +
  ggplot2::geom_jitter(
    width = 0.15,
    height = 0,
    alpha = 0.55,
    size = 1.3
  ) +
  ggplot2::stat_summary(
    ggplot2::aes(group = k),
    fun = stats::median,
    geom = "crossbar",
    width = 0.55,
    colour = "black",
    linewidth = 0.45
  ) +
  ggplot2::facet_wrap(~ model, scales = "free_y") +
  ggplot2::labs(
    x = "Number of nearest neighbours (k)",
    y = "Moran's I minus its null expectation",
    colour = "BMS"
  ) +
  ggplot2::theme_bw()

ggplot2::ggsave(
  file.path(figure_dir, "moran_knn_sensitivity.png"),
  p_knn,
  width = 10,
  height = 7,
  dpi = 300
)

p_distance <- ggplot2::ggplot(
  distance_summary,
  ggplot2::aes(x = distance_mid_km, y = median_excess_I)
) +
  ggplot2::geom_hline(
    yintercept = 0,
    linetype = 2,
    colour = "grey40"
  ) +
  ggplot2::geom_ribbon(
    ggplot2::aes(
      ymin = q25_excess_I,
      ymax = q75_excess_I
    ),
    fill = "#2166AC",
    alpha = 0.20
  ) +
  ggplot2::geom_line(colour = "#2166AC", linewidth = 0.8) +
  ggplot2::geom_point(colour = "#2166AC", size = 2) +
  ggplot2::facet_wrap(~ model, scales = "free_y") +
  ggplot2::scale_x_continuous(
    breaks = (distance_bands$d_min_km + distance_bands$d_max_km) / 2,
    labels = distance_bands$distance_band
  ) +
  ggplot2::labs(
    x = "Distance band",
    y = "Median Moran's I minus its null expectation"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 35, hjust = 1)
  )

ggplot2::ggsave(
  file.path(figure_dir, "moran_distance_correlogram.png"),
  p_distance,
  width = 10,
  height = 7,
  dpi = 300
)

# Map the four strongest BMS-year patterns for each model at k = 5.
map_k <- if (5L %in% valid_knn$k) 5L else min(valid_knn$k)

strongest_groups <- valid_knn |>
  dplyr::filter(k == map_k) |>
  dplyr::group_by(model) |>
  dplyr::slice_max(excess_I, n = 4, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::select(model, bms_id, YEAR, moran_I, excess_I)

map_data <- site_year_residuals |>
  dplyr::inner_join(
    strongest_groups,
    by = c("model", "bms_id", "YEAR")
  ) |>
  dplyr::mutate(
    x_km = x_3035 / 1000,
    y_km = y_3035 / 1000,
    panel = paste0(
      bms_id,
      " – ",
      YEAR,
      " | I = ",
      round(moran_I, 3)
    )
  )

if (nrow(map_data) > 0) {
  residual_limit <- stats::quantile(
    abs(map_data$scaled_residual),
    0.98,
    na.rm = TRUE
  )

  p_maps <- ggplot2::ggplot(
    map_data,
    ggplot2::aes(
      x = x_km,
      y = y_km,
      colour = scaled_residual
    )
  ) +
    ggplot2::geom_point(size = 1.8, alpha = 0.85) +
    ggplot2::scale_colour_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-residual_limit, residual_limit),
      oob = scales::squish
    ) +
    ggplot2::facet_wrap(
      ~ model + panel,
      scales = "free"
    ) +
    ggplot2::theme(aspect.ratio = 1) +
    ggplot2::labs(
      x = "EPSG:3035 x (km)",
      y = "EPSG:3035 y (km)",
      colour = "Mean scaled\nresidual"
    ) +
    ggplot2::theme_bw()

  ggplot2::ggsave(
    file.path(figure_dir, "strongest_spatial_residual_patterns.png"),
    p_maps,
    width = 13,
    height = 10,
    dpi = 300
  )
}


# ---------------------------------------------------------------------------- #
# 4. Optional spatial-block-year sensitivity of focal coefficients
# ---------------------------------------------------------------------------- #

coefficient_table <- function(model) {
  coefficient_matrix <- as.data.frame(summary(model)$coefficients)
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

model_convergence_message <- function(model) {
  messages <- model@optinfo$conv$lme4$messages
  if (is.null(messages)) "OK" else paste(messages, collapse = "; ")
}

fit_one_block_sensitivity <- function(
    entry,
    model_name,
    coordinates,
    block_km = NA_real_,
    shift_fraction = NA_real_) {

  original_model <- entry$model
  d <- match_entry_data(entry) |>
    dplyr::mutate(
      SITE_ID = as.character(SITE_ID),
      bms_id = as.character(bms_id),
      YEAR = as.integer(YEAR)
    ) |>
    dplyr::left_join(
      coordinates,
      by = c("SITE_ID", "bms_id")
    )

  if (any(!is.finite(d$x_3035)) || any(!is.finite(d$y_3035))) {
    stop(model_name, ": some fitted observations lack EPSG:3035 coordinates.")
  }

  d <- d |>
    dplyr::mutate(
      site_year_id = interaction(
        SITE_ID,
        YEAR,
        drop = TRUE,
        lex.order = TRUE
      )
    )

  include_spatial_block <- is.finite(block_km)

  if (include_spatial_block) {
    block_size_m <- block_km * 1000
    shift_m <- shift_fraction * block_size_m

    d <- d |>
      dplyr::mutate(
        block_x = floor((x_3035 + shift_m) / block_size_m),
        block_y = floor((y_3035 + shift_m) / block_size_m),
        spatial_block = interaction(
          bms_id,
          block_x,
          block_y,
          drop = TRUE,
          lex.order = TRUE
        ),
        spacetime_block = interaction(
          spatial_block,
          YEAR,
          drop = TRUE,
          lex.order = TRUE
        )
      )

    block_sizes <- d |>
      dplyr::distinct(spacetime_block, SITE_ID) |>
      dplyr::count(spacetime_block, name = "n_sites")

    sensitivity_formula <- stats::update(
      stats::formula(original_model),
      . ~ . + (1 | site_year_id) + (1 | spacetime_block)
    )
  } else {
    block_sizes <- tibble::tibble(n_sites = integer())

    sensitivity_formula <- stats::update(
      stats::formula(original_model),
      . ~ . + (1 | site_year_id)
    )
  }

  fitted_model <- lmerTest::lmer(
    sensitivity_formula,
    data = d,
    REML = FALSE,
    control = lme4::lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 2e6)
    )
  )

  anomaly <- paste0("clim_anomaly_tw", entry$best_window)
  focal_terms <- names(lme4::fixef(original_model))
  focal_terms <- focal_terms[
    focal_terms == anomaly |
      (
        grepl(anomaly, focal_terms, fixed = TRUE) &
          grepl(":", focal_terms, fixed = TRUE)
      )
  ]

  original_coefficients <- coefficient_table(original_model) |>
    dplyr::filter(term %in% focal_terms) |>
    dplyr::rename_with(
      ~ paste0(.x, "_original"),
      -term
    )

  sensitivity_coefficients <- coefficient_table(fitted_model) |>
    dplyr::filter(term %in% focal_terms) |>
    dplyr::rename_with(
      ~ paste0(.x, "_spatial"),
      -term
    )

  coefficient_comparison <- original_coefficients |>
    dplyr::left_join(sensitivity_coefficients, by = "term") |>
    dplyr::mutate(
      model = model_name,
      model_variant = if (include_spatial_block) {
        "site_year_plus_spatial_block"
      } else {
        "site_year_only"
      },
      block_km = block_km,
      shift_fraction = shift_fraction,
      estimate_change = estimate_spatial - estimate_original,
      relative_estimate_change = dplyr::if_else(
        estimate_original == 0,
        NA_real_,
        estimate_change / abs(estimate_original)
      ),
      SE_ratio = SE_spatial / SE_original,
      .before = 1
    )

  variance_components <- as.data.frame(lme4::VarCorr(fitted_model))

  fit_summary <- tibble::tibble(
    model = model_name,
    model_variant = if (include_spatial_block) {
      "site_year_plus_spatial_block"
    } else {
      "site_year_only"
    },
    block_km = block_km,
    shift_fraction = shift_fraction,
    n_spacetime_blocks = if (include_spatial_block) {
      nrow(block_sizes)
    } else {
      NA_integer_
    },
    median_sites_per_spacetime_block = if (include_spatial_block) {
      stats::median(block_sizes$n_sites)
    } else {
      NA_real_
    },
    proportion_multisite_blocks = if (include_spatial_block) {
      mean(block_sizes$n_sites >= 2)
    } else {
      NA_real_
    },
    site_year_SD = variance_components$sdcor[
      match("site_year_id", variance_components$grp)
    ],
    spacetime_block_SD = if (include_spatial_block) {
      variance_components$sdcor[
        match("spacetime_block", variance_components$grp)
      ]
    } else {
      NA_real_
    },
    singular_fit = lme4::isSingular(fitted_model, tol = 1e-4),
    convergence = model_convergence_message(fitted_model),
    AIC = stats::AIC(fitted_model)
  )

  list(
    coefficients = coefficient_comparison,
    fit = fit_summary
  )
}

if (run_spatial_block_sensitivity) {
  message("Refitting spatial-block-year sensitivity models...")

  block_specifications <- tidyr::expand_grid(
    block_km = block_sizes_km,
    shift_fraction = block_origin_shifts
  )

  block_coefficients <- list()
  block_fit_summary <- list()
  result_index <- 1L

  for (model_name in names(model_entries)) {
    message("Spatial sensitivity: ", model_name, " | site-year only")

    site_year_result <- tryCatch(
      fit_one_block_sensitivity(
        entry = model_entries[[model_name]],
        model_name = model_name,
        coordinates = site_coordinates
      ),
      error = function(e) {
        list(
          coefficients = tibble::tibble(),
          fit = tibble::tibble(
            model = model_name,
            model_variant = "site_year_only",
            block_km = NA_real_,
            shift_fraction = NA_real_,
            n_spacetime_blocks = NA_integer_,
            median_sites_per_spacetime_block = NA_real_,
            proportion_multisite_blocks = NA_real_,
            site_year_SD = NA_real_,
            spacetime_block_SD = NA_real_,
            singular_fit = NA,
            convergence = paste0("ERROR: ", conditionMessage(e)),
            AIC = NA_real_
          )
        )
      }
    )

    block_coefficients[[result_index]] <- site_year_result$coefficients
    block_fit_summary[[result_index]] <- site_year_result$fit
    result_index <- result_index + 1L

    for (specification in seq_len(nrow(block_specifications))) {
      block_km <- block_specifications$block_km[specification]
      shift_fraction <- block_specifications$shift_fraction[specification]

      message(
        "Spatial sensitivity: ", model_name,
        " | block = ", block_km, " km",
        " | shift = ", shift_fraction
      )

      result <- tryCatch(
        fit_one_block_sensitivity(
          entry = model_entries[[model_name]],
          model_name = model_name,
          coordinates = site_coordinates,
          block_km = block_km,
          shift_fraction = shift_fraction
        ),
        error = function(e) {
          list(
            coefficients = tibble::tibble(),
            fit = tibble::tibble(
              model = model_name,
              model_variant = "site_year_plus_spatial_block",
              block_km = block_km,
              shift_fraction = shift_fraction,
              n_spacetime_blocks = NA_integer_,
              median_sites_per_spacetime_block = NA_real_,
              proportion_multisite_blocks = NA_real_,
              site_year_SD = NA_real_,
              spacetime_block_SD = NA_real_,
              singular_fit = NA,
              convergence = paste0("ERROR: ", conditionMessage(e)),
              AIC = NA_real_
            )
          )
        }
      )

      block_coefficients[[result_index]] <- result$coefficients
      block_fit_summary[[result_index]] <- result$fit
      result_index <- result_index + 1L
    }
  }

  block_coefficients <- dplyr::bind_rows(block_coefficients)
  block_fit_summary <- dplyr::bind_rows(block_fit_summary)
} else {
  block_coefficients <- tibble::tibble()
  block_fit_summary <- tibble::tibble()
}


# ---------------------------------------------------------------------------- #
# Export
# ---------------------------------------------------------------------------- #

output_tables <- list(
  knn_site_year = knn_site_year,
  knn_summary = knn_summary,
  knn_summary_BMS = knn_summary_by_BMS,
  distance_site_year = distance_site_year,
  distance_summary = distance_summary,
  distance_summary_BMS = distance_summary_by_BMS,
  species_year = knn_species_year,
  species_summary = knn_species_summary,
  block_coefficients = block_coefficients,
  block_fit_summary = block_fit_summary
)

output_tables <- output_tables[
  vapply(output_tables, ncol, integer(1)) > 0
]

output_xlsx <- file.path(
  out_dir,
  "plasticity_spatial_diagnostics_deep.xlsx"
)

writexl::write_xlsx(
  output_tables,
  path = output_xlsx
)

message("Saved spatial diagnostic tables: ", output_xlsx)
message("Saved spatial diagnostic figures: ", figure_dir)

