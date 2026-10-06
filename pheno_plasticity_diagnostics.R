# ============================================================================================ #
# pheno_plasticity_diagnostics.R
#
# Residual diagnostics for the selected phenological plasticity models:
#   3. residual temporal and spatial dependence
#   7. residual distribution, heteroscedasticity and extreme populations
#
# Run pheno_plasticity.R once before this script. That script creates:
# output/phenology_plasticity/plasticity_main_models.rds
# ============================================================================================ #

required_packages <- c(
  "dplyr", "tidyr", "purrr", "tibble", "ggplot2",
  "lme4", "lmerTest", "spdep", "writexl", "here"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Install the missing packages before running diagnostics: ",
    paste(missing_packages, collapse = ", ")
  )
}

library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(ggplot2)
library(lme4)
library(lmerTest)

# Use c("primary") here if only the four primary models should be checked.
analysis_groups <- c("primary", "sensitivity_1zero")

rds_path <- here::here(
  "output",
  "phenology_plasticity",
  "plasticity_main_models.rds"
)

if (!file.exists(rds_path)) {
  stop(
    "Model bundle not found: ", rds_path, "\n",
    "Run the updated pheno_plasticity.R once to fit and save the models."
  )
}

model_bundle <- readRDS(rds_path)

out_dir <- here::here(
  "output",
  "phenology_plasticity",
  "diagnostics"
)

figure_dir <- file.path(out_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)


# ---------------------------------------------------------------------------- #
# Helpers
# ---------------------------------------------------------------------------- #

safe_file_name <- function(x) {
  gsub("[^A-Za-z0-9_]+", "_", x)
}

flatten_model_entries <- function(bundle, groups) {
  entries <- list()

  for (group_name in intersect(groups, names(bundle))) {
    group_entries <- bundle[[group_name]]

    for (model_name in names(group_entries)) {
      entries[[paste(group_name, model_name, sep = "__")]] <-
        group_entries[[model_name]]
    }
  }

  entries
}

extract_residual_data <- function(entry) {
  model <- entry$model
  original_data <- entry$data
  model_data <- model.frame(model)

  # The saved data retain their original row names. Match them to the model
  # frame so residuals are joined to the correct year and population.
  row_index <- match(rownames(model_data), rownames(original_data))

  if (anyNA(row_index)) {
    if (nrow(model_data) == nrow(original_data)) {
      warning(
        "Could not match row names; using row order because row counts agree."
      )
      row_index <- seq_len(nrow(original_data))
    } else {
      stop("Could not match the saved data rows to the fitted model.")
    }
  }

  d <- original_data[row_index, , drop = FALSE]

  required_columns <- c("YEAR", "SITE_ID", "SPECIES")
  missing_columns <- setdiff(required_columns, names(d))

  if (length(missing_columns) > 0) {
    stop(
      "Saved model data lack required columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  sp_site <- if ("sp_site" %in% names(d)) {
    as.character(d$sp_site)
  } else {
    paste(d$SPECIES, d$SITE_ID, sep = "__")
  }

  bms_id <- if ("bms_id" %in% names(d)) {
    as.character(d$bms_id)
  } else {
    rep(NA_character_, nrow(d))
  }

  raw_residual <- residuals(model)
  model_sigma <- sigma(model)

  tibble(
    YEAR = as.integer(d$YEAR),
    SITE_ID = as.character(d$SITE_ID),
    SPECIES = as.character(d$SPECIES),
    sp_site = sp_site,
    bms_id = bms_id,
    fitted = as.numeric(fitted(model)),
    residual = as.numeric(raw_residual),
    scaled_residual = as.numeric(raw_residual / model_sigma)
  )
}

safe_skewness <- function(x) {
  x <- x[is.finite(x)]
  s <- stats::sd(x)
  if (length(x) < 3 || !is.finite(s) || s == 0) return(NA_real_)
  mean(((x - mean(x)) / s)^3)
}

safe_excess_kurtosis <- function(x) {
  x <- x[is.finite(x)]
  s <- stats::sd(x)
  if (length(x) < 4 || !is.finite(s) || s == 0) return(NA_real_)
  mean(((x - mean(x)) / s)^4) - 3
}

model_check_summary <- function(
    fitted_model,
    residual_data,
    model_name,
    entry) {
  
  convergence_message <- fitted_model@optinfo$conv$lme4$messages
  convergence_message <- if (is.null(convergence_message)) {
    "OK"
  } else {
    paste(convergence_message, collapse = "; ")
  }
  
  gradient <- tryCatch(
    fitted_model@optinfo$derivs$gradient,
    error = function(e) NA_real_
  )
  
  max_abs_gradient <- if (all(is.na(gradient))) {
    NA_real_
  } else {
    max(abs(gradient), na.rm = TRUE)
  }
  
  r <- residual_data$scaled_residual
  
  tibble(
    model = model_name,
    analysis = entry$analysis,
    response = entry$response,
    best_window = entry$best_window,
    n_observations = nrow(stats::model.frame(fitted_model)),
    residual_sigma = stats::sigma(fitted_model),
    residual_skewness = safe_skewness(r),
    residual_excess_kurtosis = safe_excess_kurtosis(r),
    proportion_abs_residual_gt_2 = mean(abs(r) > 2, na.rm = TRUE),
    proportion_abs_residual_gt_3 = mean(abs(r) > 3, na.rm = TRUE),
    singular_fit = lme4::isSingular(fitted_model, tol = 1e-4),
    max_abs_gradient = max_abs_gradient,
    convergence = convergence_message
  )
}

save_residual_plots <- function(residual_data, model_name, file_name) {
  # Plot at most 30,000 points to keep large diagnostics responsive.
  set.seed(42)
  plot_data <- if (nrow(residual_data) > 30000) {
    dplyr::slice_sample(residual_data, n = 30000)
  } else {
    residual_data
  }

  grDevices::png(
    filename = file_name,
    width = 2400,
    height = 1800,
    res = 250
  )

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old_par)
    grDevices::dev.off()
  }, add = TRUE)

  graphics::par(mfrow = c(2, 2), mar = c(4.3, 4.3, 3.0, 1.0))

  graphics::plot(
    plot_data$fitted,
    plot_data$scaled_residual,
    pch = 16,
    cex = 0.35,
    col = grDevices::adjustcolor("black", alpha.f = 0.15),
    xlab = "Conditional fitted value",
    ylab = "Scaled conditional residual",
    main = "Residuals vs fitted"
  )
  graphics::abline(h = 0, lty = 2, col = "grey40")
  if (length(unique(plot_data$fitted)) > 2) {
    graphics::lines(
      stats::lowess(plot_data$fitted, plot_data$scaled_residual),
      col = "red3",
      lwd = 2
    )
  }

  stats::qqnorm(
    plot_data$scaled_residual,
    pch = 16,
    cex = 0.35,
    col = grDevices::adjustcolor("black", alpha.f = 0.2),
    main = "Normal Q-Q",
    ylab = "Scaled conditional residual"
  )
  stats::qqline(plot_data$scaled_residual, col = "red3", lwd = 2)

  graphics::hist(
    plot_data$scaled_residual,
    breaks = 40,
    col = "grey80",
    border = "white",
    main = "Residual distribution",
    xlab = "Scaled conditional residual"
  )

  graphics::plot(
    plot_data$fitted,
    sqrt(abs(plot_data$scaled_residual)),
    pch = 16,
    cex = 0.35,
    col = grDevices::adjustcolor("black", alpha.f = 0.15),
    xlab = "Conditional fitted value",
    ylab = expression(sqrt("|scaled residual|")),
    main = "Scale-location"
  )
  if (length(unique(plot_data$fitted)) > 2) {
    graphics::lines(
      stats::lowess(
        plot_data$fitted,
        sqrt(abs(plot_data$scaled_residual))
      ),
      col = "red3",
      lwd = 2
    )
  }

  graphics::mtext(
    model_name,
    side = 3,
    outer = TRUE,
    line = -1.2,
    font = 2
  )
}

save_random_effect_plots <- function(model, model_name, file_name) {
  random_effects <- lme4::ranef(model)

  plot_index <- do.call(
    rbind,
    lapply(names(random_effects), function(group_name) {
      data.frame(
        group = group_name,
        term = names(random_effects[[group_name]]),
        stringsAsFactors = FALSE
      )
    })
  )

  if (is.null(plot_index) || nrow(plot_index) == 0) return(invisible(NULL))

  n_rows <- ceiling(nrow(plot_index) / 2)

  grDevices::png(
    filename = file_name,
    width = 2200,
    height = max(1200, 850 * n_rows),
    res = 250
  )

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old_par)
    grDevices::dev.off()
  }, add = TRUE)

  graphics::par(
    mfrow = c(n_rows, 2),
    mar = c(4.3, 4.3, 3.0, 1.0),
    oma = c(0, 0, 2, 0)
  )

  for (i in seq_len(nrow(plot_index))) {
    group_name <- plot_index$group[i]
    term_name <- plot_index$term[i]
    values <- random_effects[[group_name]][[term_name]]

    stats::qqnorm(
      values,
      pch = 16,
      cex = 0.5,
      col = grDevices::adjustcolor("black", alpha.f = 0.35),
      main = paste(group_name, term_name, sep = ": "),
      ylab = "Estimated random effect"
    )
    stats::qqline(values, col = "red3", lwd = 2)
  }

  if (nrow(plot_index) %% 2 == 1) graphics::plot.new()

  graphics::mtext(
    paste(model_name, "- random-effect Q-Q plots"),
    side = 3,
    outer = TRUE,
    line = 0.3,
    font = 2
  )
}


# ---------------------------------------------------------------------------- #
# Temporal residual dependence
# ---------------------------------------------------------------------------- #

temporal_diagnostics <- function(residual_data, model_name, figure_dir) {
  # Collapse any accidental duplicates to one residual per population-year.
  pop_year <- residual_data |>
    group_by(sp_site, YEAR) |>
    summarise(
      scaled_residual = mean(scaled_residual, na.rm = TRUE),
      .groups = "drop"
    ) |>
    arrange(sp_site, YEAR) |>
    group_by(sp_site) |>
    mutate(
      previous_year = lag(YEAR),
      lag_residual = lag(scaled_residual),
      consecutive = YEAR - previous_year == 1L
    ) |>
    ungroup()

  consecutive_pairs <- pop_year |>
    filter(
      consecutive,
      is.finite(lag_residual),
      is.finite(scaled_residual)
    )

  n_pairs <- nrow(consecutive_pairs)
  pooled_rho <- if (
    n_pairs >= 4 &&
      stats::sd(consecutive_pairs$lag_residual) > 0 &&
      stats::sd(consecutive_pairs$scaled_residual) > 0
  ) {
    stats::cor(
      consecutive_pairs$scaled_residual,
      consecutive_pairs$lag_residual
    )
  } else {
    NA_real_
  }

  # Fisher intervals are descriptive here: consecutive pairs within a
  # population are not completely independent.
  fisher_ci <- if (n_pairs > 3 && is.finite(pooled_rho)) {
    z <- atanh(max(min(pooled_rho, 0.999999), -0.999999))
    se <- 1 / sqrt(n_pairs - 3)
    tanh(z + c(-1, 1) * 1.96 * se)
  } else {
    c(NA_real_, NA_real_)
  }

  population_rho <- consecutive_pairs |>
    group_by(sp_site) |>
    summarise(
      n_consecutive_pairs = n(),
      lag1_rho = if (
        n_consecutive_pairs >= 4 &&
          stats::sd(lag_residual) > 0 &&
          stats::sd(scaled_residual) > 0
      ) {
        stats::cor(scaled_residual, lag_residual)
      } else {
        NA_real_
      },
      .groups = "drop"
    ) |>
    filter(is.finite(lag1_rho)) |>
    mutate(model = model_name, .before = 1)

  year_summary <- residual_data |>
    group_by(YEAR) |>
    summarise(
      n = n(),
      mean_residual = mean(scaled_residual, na.rm = TRUE),
      se_residual = stats::sd(scaled_residual, na.rm = TRUE) / sqrt(n),
      .groups = "drop"
    )

  temporal_summary <- tibble(
    model = model_name,
    n_consecutive_pairs = n_pairs,
    pooled_lag1_rho = pooled_rho,
    pooled_lag1_lower_95 = fisher_ci[1],
    pooled_lag1_upper_95 = fisher_ci[2],
    n_populations_with_rho = nrow(population_rho),
    median_population_rho = if (nrow(population_rho) > 0) {
      stats::median(population_rho$lag1_rho)
    } else {
      NA_real_
    },
    population_rho_q25 = if (nrow(population_rho) > 0) {
      stats::quantile(population_rho$lag1_rho, 0.25, names = FALSE)
    } else {
      NA_real_
    },
    population_rho_q75 = if (nrow(population_rho) > 0) {
      stats::quantile(population_rho$lag1_rho, 0.75, names = FALSE)
    } else {
      NA_real_
    },
    max_abs_annual_mean_residual = if (nrow(year_summary) > 0) {
      max(abs(year_summary$mean_residual), na.rm = TRUE)
    } else {
      NA_real_
    }
  )

  if (nrow(consecutive_pairs) > 0) {
    set.seed(42)
    pair_plot_data <- if (nrow(consecutive_pairs) > 30000) {
      dplyr::slice_sample(consecutive_pairs, n = 30000)
    } else {
      consecutive_pairs
    }

    p_lag <- ggplot(
      pair_plot_data,
      aes(x = lag_residual, y = scaled_residual)
    ) +
      geom_hline(yintercept = 0, colour = "grey75") +
      geom_vline(xintercept = 0, colour = "grey75") +
      geom_point(alpha = 0.12, size = 0.7) +
      geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "red3") +
      labs(
        title = model_name,
        subtitle = paste0(
          "Consecutive years only; pooled lag-1 r = ",
          round(pooled_rho, 3)
        ),
        x = "Residual in year t - 1",
        y = "Residual in year t"
      ) +
      theme_bw()

    ggsave(
      file.path(
        figure_dir,
        paste0(safe_file_name(model_name), "_temporal_lag1.png")
      ),
      p_lag,
      width = 6.5,
      height = 5.2,
      dpi = 300
    )
  }

  p_year <- ggplot(year_summary, aes(x = YEAR, y = mean_residual)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
    geom_ribbon(
      aes(
        ymin = mean_residual - 1.96 * se_residual,
        ymax = mean_residual + 1.96 * se_residual
      ),
      fill = "grey75",
      alpha = 0.5
    ) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 1.4) +
    labs(
      title = model_name,
      subtitle = "Mean conditional residual shared across observations",
      x = "Year",
      y = "Mean scaled residual"
    ) +
    theme_bw()

  ggsave(
    file.path(
      figure_dir,
      paste0(safe_file_name(model_name), "_annual_residuals.png")
    ),
    p_year,
    width = 7.2,
    height = 4.8,
    dpi = 300
  )

  list(
    summary = temporal_summary,
    population_rho = population_rho,
    year_summary = year_summary |>
      mutate(model = model_name, .before = 1)
  )
}


# ---------------------------------------------------------------------------- #
# Spatial residual dependence
# ---------------------------------------------------------------------------- #

moran_one_group <- function(d, k_max = 5L, minimum_sites = 10L) {
  # If two sites share exactly the same coordinates, average them so the
  # nearest-neighbour graph remains well defined.
  d <- d |>
    group_by(x_3035, y_3035) |>
    summarise(
      scaled_residual = mean(scaled_residual, na.rm = TRUE),
      .groups = "drop"
    ) |>
    filter(is.finite(scaled_residual))

  n_sites <- nrow(d)

  empty_result <- function(reason) {
    tibble(
      n_sites = n_sites,
      k = NA_integer_,
      moran_I = NA_real_,
      expected_I = NA_real_,
      variance_I = NA_real_,
      p_value = NA_real_,
      status = reason
    )
  }

  if (n_sites < minimum_sites) return(empty_result("too_few_sites"))
  if (stats::sd(d$scaled_residual) == 0) {
    return(empty_result("zero_residual_variance"))
  }

  k_use <- min(as.integer(k_max), n_sites - 1L)
  coordinates <- as.matrix(d[, c("x_3035", "y_3035")])

  tryCatch({
    neighbour_graph <- spdep::knearneigh(coordinates, k = k_use)
    neighbours <- spdep::knn2nb(neighbour_graph, sym = TRUE)
    weights <- spdep::nb2listw(
      neighbours,
      style = "W",
      zero.policy = TRUE
    )

    test <- spdep::moran.test(
      d$scaled_residual,
      listw = weights,
      randomisation = TRUE,
      alternative = "two.sided",
      zero.policy = TRUE
    )

    tibble(
      n_sites = n_sites,
      k = k_use,
      moran_I = unname(test$estimate[["Moran I statistic"]]),
      expected_I = unname(test$estimate[["Expectation"]]),
      variance_I = unname(test$estimate[["Variance"]]),
      p_value = test$p.value,
      status = "OK"
    )
  }, error = function(e) {
    empty_result(paste0("error: ", conditionMessage(e)))
  })
}

spatial_diagnostics <- function(
    residual_data,
    site_coordinates,
    model_name,
    figure_dir) {

  if (all(is.na(residual_data$bms_id))) {
    warning(model_name, ": bms_id is unavailable; skipping spatial test.")
    return(list(summary = tibble(), by_year = tibble()))
  }

  site_year <- residual_data |>
    group_by(bms_id, YEAR, SITE_ID) |>
    summarise(
      scaled_residual = mean(scaled_residual, na.rm = TRUE),
      .groups = "drop"
    ) |>
    left_join(
      site_coordinates,
      by = c("SITE_ID", "bms_id")
    ) |>
    filter(
      is.finite(x_3035),
      is.finite(y_3035),
      is.finite(scaled_residual)
    )

  if (nrow(site_year) == 0) {
    warning(model_name, ": no residuals could be joined to coordinates.")
    return(list(summary = tibble(), by_year = tibble()))
  }

  moran_by_year <- site_year |>
    group_by(bms_id, YEAR) |>
    tidyr::nest() |>
    mutate(result = map(data, moran_one_group)) |>
    select(-data) |>
    tidyr::unnest(result) |>
    mutate(
      model = model_name,
      p_FDR = ifelse(
        status == "OK",
        stats::p.adjust(p_value, method = "BH"),
        NA_real_
      ),
      .before = 1
    )

  valid_moran <- moran_by_year |>
    filter(status == "OK", is.finite(moran_I))

  spatial_summary <- tibble(
    model = model_name,
    n_bms_year_tests = nrow(valid_moran),
    mean_moran_I = if (nrow(valid_moran) > 0) {
      mean(valid_moran$moran_I)
    } else {
      NA_real_
    },
    median_moran_I = if (nrow(valid_moran) > 0) {
      stats::median(valid_moran$moran_I)
    } else {
      NA_real_
    },
    proportion_positive = if (nrow(valid_moran) > 0) {
      mean(valid_moran$moran_I > 0)
    } else {
      NA_real_
    },
    n_positive_p_lt_0_05 = sum(
      valid_moran$moran_I > 0 & valid_moran$p_value < 0.05,
      na.rm = TRUE
    ),
    n_positive_FDR_lt_0_05 = sum(
      valid_moran$moran_I > 0 & valid_moran$p_FDR < 0.05,
      na.rm = TRUE
    )
  )

  if (nrow(valid_moran) > 0) {
    p_spatial <- ggplot(
      valid_moran,
      aes(x = YEAR, y = moran_I, colour = bms_id)
    ) +
      geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
      geom_point(aes(shape = p_FDR < 0.05), size = 2, alpha = 0.8) +
      labs(
        title = model_name,
        subtitle = "Moran's I on site-year mean residuals, within BMS and year",
        x = "Year",
        y = "Moran's I",
        colour = "BMS",
        shape = "FDR < 0.05"
      ) +
      theme_bw() +
      theme(legend.position = "right")

    ggsave(
      file.path(
        figure_dir,
        paste0(safe_file_name(model_name), "_spatial_moran.png")
      ),
      p_spatial,
      width = 8.5,
      height = 5.2,
      dpi = 300
    )
  }

  list(summary = spatial_summary, by_year = moran_by_year)
}


# ---------------------------------------------------------------------------- #
# Run all diagnostics
# ---------------------------------------------------------------------------- #

model_entries <- flatten_model_entries(model_bundle, analysis_groups)

if (length(model_entries) == 0) {
  stop("No model entries found for: ", paste(analysis_groups, collapse = ", "))
}

diagnostic_results <- purrr::imap(model_entries, function(entry, model_name) {
  message("Running diagnostics: ", model_name)

  residual_data <- extract_residual_data(entry)

  checks <- model_check_summary(
    entry$model,
    residual_data,
    model_name,
    entry
  )

  save_residual_plots(
    residual_data,
    model_name,
    file.path(
      figure_dir,
      paste0(safe_file_name(model_name), "_residual_checks.png")
    )
  )

  save_random_effect_plots(
    entry$model,
    model_name,
    file.path(
      figure_dir,
      paste0(safe_file_name(model_name), "_random_effect_qq.png")
    )
  )

  influential_populations <- residual_data |>
    group_by(sp_site, SPECIES, SITE_ID) |>
    summarise(
      n = n(),
      rms_scaled_residual = sqrt(mean(scaled_residual^2, na.rm = TRUE)),
      max_abs_scaled_residual = max(abs(scaled_residual), na.rm = TRUE),
      .groups = "drop"
    ) |>
    arrange(desc(max_abs_scaled_residual)) |>
    slice_head(n = 25) |>
    mutate(model = model_name, .before = 1)

  temporal <- temporal_diagnostics(
    residual_data,
    model_name,
    figure_dir
  )

  spatial <- spatial_diagnostics(
    residual_data,
    model_bundle$site_coordinates,
    model_name,
    figure_dir
  )

  list(
    checks = checks,
    influential_populations = influential_populations,
    temporal = temporal,
    spatial = spatial
  )
})

model_checks <- map_dfr(diagnostic_results, "checks")
influential_populations <- map_dfr(
  diagnostic_results,
  "influential_populations"
)
temporal_summary <- map_dfr(
  diagnostic_results,
  c("temporal", "summary")
)
temporal_by_population <- map_dfr(
  diagnostic_results,
  c("temporal", "population_rho")
)
annual_residuals <- map_dfr(
  diagnostic_results,
  c("temporal", "year_summary")
)
spatial_summary <- map_dfr(
  diagnostic_results,
  c("spatial", "summary")
)
spatial_moran_by_year <- map_dfr(
  diagnostic_results,
  c("spatial", "by_year")
)

diagnostic_workbook <- file.path(
  out_dir,
  "plasticity_residual_diagnostics.xlsx"
)

writexl::write_xlsx(
  list(
    model_checks = model_checks,
    temporal_summary = temporal_summary,
    temporal_by_population = temporal_by_population,
    annual_residuals = annual_residuals,
    spatial_summary = spatial_summary,
    spatial_moran_by_year = spatial_moran_by_year,
    influential_populations = influential_populations
  ),
  path = diagnostic_workbook
)

message("Diagnostics written to: ", diagnostic_workbook)
message("Diagnostic figures written to: ", figure_dir)
