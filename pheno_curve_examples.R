# ============================================================================================ #
# pheno_curve_examples.R
#
# Author: Pau Colom
# Date: 2026-06-01
#
# Description: This script prepares panel B of Figure 1. It shows phenological curves
# for multiple years from one univoltine and one multivoltine population. Fitted GAM
# curves are coloured by thermal anomaly, illustrating how seasonal phenology shifts
# among colder and warmer years.
#
# ============================================================================================ #

# Clean session
rm(list = ls())

#### Load required libraries ####
# ---
library(dplyr)
library(data.table)
library(mgcv)
library(lubridate)
library(ggplot2)
library(here)
library(extrafont)
library(purrr)
# ---

#### Data import ####
# ---

here::here() # Check current working directory

ebms_count_df <- read.csv(
  here::here("data", "ebms_count.csv"),
  sep = ",",
  dec = "."
)

ebms_visit_df <- read.csv(
  here::here("data", "ebms_visit.csv"),
  sep = ",",
  dec = "."
)

phenology_estimates <- read.csv(
  here::here("output", "pheno_estimates_allspp.csv"),
  sep = ",",
  dec = "."
)

clim_vars <- read.csv(
  here::here("output", "climate", "climate_variables_all_phenophases.csv"),
  sep = ",",
  dec = "."
)
# ---

# ---------------------------------------------------------------------------- #
#### User options ####
# ---------------------------------------------------------------------------- #

# Thermal anomaly used to colour curves
anomaly_var <- "clim_anomaly_tw90"

# Climate row used to extract the anomaly.
# For the figure, ONSET_mean is a clear option because it represents the
# pre-onset thermal window.
climate_pheno_type <- "ONSET_mean"

# Number of years to show per population
n_years_to_plot <- 20

# Minimum data quality for candidate years
min_visits <- 15
min_positive_dates <- 3
min_total_count <- 5

# Examples: one univoltine and one multivoltine population.
# Replace these if better examples are found.
examples <- tibble::tribble(
  ~example_label,  ~SPECIES,                  ~SITE_ID,
  "Univoltine",    "Anthocharis cardamines", "ES-CTBMS.5",
  "Multivoltine",  "Lysandra bellargus",    "UKBMS.2010"
)

# ---------------------------------------------------------------------------- #
#### Prepare data ####
# ---------------------------------------------------------------------------- #

m_count <- data.table(ebms_count_df)
m_visit <- data.table(ebms_visit_df)

setnames(
  m_count,
  c("transect_id", "visit_date", "species_name", "count"),
  c("SITE_ID", "DATE", "SPECIES", "COUNT")
)

setnames(
  m_visit,
  c("transect_id", "visit_date"),
  c("SITE_ID", "DATE")
)

m_count[, DATE := as.Date(DATE)]
m_visit[, DATE := as.Date(DATE)]

m_count[, year := as.integer(year)]
m_visit[, year := as.integer(year)]

phenology_estimates <- phenology_estimates |>
  dplyr::mutate(
    SITE_ID = as.character(SITE_ID),
    SPECIES = as.character(SPECIES),
    YEAR = as.integer(YEAR)
  )

clim_vars <- clim_vars |>
  dplyr::mutate(
    SITE_ID = as.character(SITE_ID),
    SPECIES = as.character(SPECIES),
    YEAR = as.integer(YEAR)
  )

# Keep one climatic anomaly per species x site x year
if ("pheno_type" %in% names(clim_vars)) {
  clim_anomaly_df <- clim_vars |>
    dplyr::filter(pheno_type == climate_pheno_type)
} else {
  clim_anomaly_df <- clim_vars
}

clim_anomaly_df <- clim_anomaly_df |>
  dplyr::select(
    SPECIES,
    SITE_ID,
    YEAR,
    anomaly = dplyr::all_of(anomaly_var)
  ) |>
  dplyr::distinct()

# ---------------------------------------------------------------------------- #
#### Data summaries for selecting years ####
# ---------------------------------------------------------------------------- #

count_stats <- m_count |>
  as.data.frame() |>
  dplyr::mutate(
    SITE_ID = as.character(SITE_ID),
    SPECIES = as.character(SPECIES),
    YEAR = as.integer(year),
    DATE = as.Date(DATE)
  ) |>
  dplyr::group_by(SPECIES, SITE_ID, YEAR) |>
  dplyr::summarise(
    n_positive_dates = dplyr::n_distinct(DATE[COUNT > 0]),
    total_count = sum(COUNT, na.rm = TRUE),
    max_count = max(COUNT, na.rm = TRUE),
    .groups = "drop"
  )

visit_stats <- m_visit |>
  as.data.frame() |>
  dplyr::mutate(
    SITE_ID = as.character(SITE_ID),
    YEAR = as.integer(year),
    DATE = as.Date(DATE)
  ) |>
  dplyr::group_by(SITE_ID, YEAR) |>
  dplyr::summarise(
    n_visits = dplyr::n_distinct(visit_id),
    .groups = "drop"
  )

# ---------------------------------------------------------------------------- #
#### Select representative years across the thermal-anomaly gradient ####
# ---------------------------------------------------------------------------- #

select_years_for_population <- function(species_pick,
                                        site_pick,
                                        example_label,
                                        phenology_estimates,
                                        clim_anomaly_df,
                                        count_stats,
                                        visit_stats,
                                        n_years_to_plot = 7,
                                        min_visits = 10,
                                        min_positive_dates = 4,
                                        min_total_count = 10) {
  
  candidates <- phenology_estimates |>
    dplyr::filter(
      SPECIES == species_pick,
      SITE_ID == site_pick,
      !is.na(ONSET_mean),
      !is.na(OFFSET_mean),
      ONSET_mean >= 60,
      OFFSET_mean <= 304
    ) |>
    dplyr::left_join(
      clim_anomaly_df,
      by = c("SPECIES", "SITE_ID", "YEAR")
    ) |>
    dplyr::left_join(
      count_stats,
      by = c("SPECIES", "SITE_ID", "YEAR")
    ) |>
    dplyr::left_join(
      visit_stats,
      by = c("SITE_ID", "YEAR")
    ) |>
    dplyr::filter(
      !is.na(anomaly),
      n_visits >= min_visits,
      n_positive_dates >= min_positive_dates,
      total_count >= min_total_count
    ) |>
    dplyr::arrange(anomaly)
  
  if (nrow(candidates) == 0) {
    stop(paste("No suitable years found for:", species_pick, site_pick))
  }
  
  idx <- unique(round(seq(
    from = 1,
    to = nrow(candidates),
    length.out = min(n_years_to_plot, nrow(candidates))
  )))
  
  candidates[idx, ] |>
    dplyr::mutate(
      example_label = example_label
    )
}

selected_years <- purrr::pmap_dfr(
  examples,
  ~ select_years_for_population(
    species_pick = ..2,
    site_pick = ..3,
    example_label = ..1,
    phenology_estimates = phenology_estimates,
    clim_anomaly_df = clim_anomaly_df,
    count_stats = count_stats,
    visit_stats = visit_stats,
    n_years_to_plot = n_years_to_plot,
    min_visits = min_visits,
    min_positive_dates = min_positive_dates,
    min_total_count = min_total_count
  )
)

selected_years |>
  dplyr::select(
    example_label,
    SPECIES,
    SITE_ID,
    YEAR,
    anomaly,
    n_visits,
    n_positive_dates,
    total_count,
    ONSET_mean,
    OFFSET_mean
  ) |>
  dplyr::arrange(example_label, anomaly)

# ---------------------------------------------------------------------------- #
#### Function to prepare one curve for one year ####
# ---------------------------------------------------------------------------- #

prepare_curve_year <- function(species_pick,
                               site_pick,
                               year_pick,
                               example_label,
                               anomaly,
                               onset_day,
                               offset_day,
                               m_count,
                               m_visit) {
  
  sub_count <- m_count[
    SPECIES == species_pick &
      SITE_ID == site_pick &
      year == year_pick
  ]
  
  sub_visit <- m_visit[
    SITE_ID == site_pick &
      year == year_pick
  ]
  
  if (nrow(sub_count) == 0) {
    stop(paste("No count data for:", species_pick, site_pick, year_pick))
  }
  
  if (nrow(sub_visit) == 0) {
    stop(paste("No visit data for:", site_pick, year_pick))
  }
  
  sub_count[, julian_day := lubridate::yday(DATE)]
  sub_visit[, julian_day := lubridate::yday(DATE)]
  
  # Add zeros for visited dates where the species was not recorded
  missing_dates <- sub_visit[!DATE %in% sub_count$DATE]
  
  missing_rows <- missing_dates |>
    as.data.frame() |>
    dplyr::mutate(
      SPECIES = species_pick,
      COUNT = 0,
      julian_day = lubridate::yday(DATE)
    ) |>
    dplyr::select(
      SITE_ID,
      DATE,
      year,
      SPECIES,
      COUNT,
      julian_day
    )
  
  observed_counts <- sub_count |>
    as.data.frame() |>
    dplyr::select(
      SITE_ID,
      DATE,
      year,
      SPECIES,
      COUNT,
      julian_day
    ) |>
    dplyr::bind_rows(missing_rows) |>
    dplyr::arrange(julian_day) |>
    dplyr::mutate(
      example_label = example_label,
      anomaly = anomaly
    )
  
  # Structural zeros outside the monitoring season
  anchor_zeros <- data.frame(
    SITE_ID = site_pick,
    DATE = as.Date(NA),
    year = year_pick,
    SPECIES = species_pick,
    COUNT = 0,
    julian_day = c(1:30, 335:365)
  )
  
  all_counts <- observed_counts |>
    dplyr::select(SITE_ID, DATE, year, SPECIES, COUNT, julian_day) |>
    dplyr::bind_rows(anchor_zeros)
  
  # Fit GAM
  gam_model <- mgcv::gam(
    COUNT ~ s(julian_day),
    data = all_counts,
    family = nb()
  )
  
  julian_days <- 1:365
  
  predicted_curve <- data.frame(
    julian_day = julian_days,
    predicted_count = predict(
      gam_model,
      newdata = data.frame(julian_day = julian_days),
      type = "response"
    ),
    example_label = example_label,
    YEAR = year_pick,
    anomaly = anomaly
  )
  
  boundary_df <- data.frame(
    example_label = example_label,
    YEAR = year_pick,
    anomaly = anomaly,
    boundary = c("Onset", "Offset"),
    x = c(onset_day, offset_day)
  )
  
  list(
    observed = observed_counts,
    predicted = predicted_curve,
    boundaries = boundary_df
  )
}

# ---------------------------------------------------------------------------- #
#### Prepare all curves ####
# ---------------------------------------------------------------------------- #

curve_list <- lapply(seq_len(nrow(selected_years)), function(i) {
  
  row_i <- selected_years[i, ]
  
  prepare_curve_year(
    species_pick = row_i$SPECIES,
    site_pick = row_i$SITE_ID,
    year_pick = row_i$YEAR,
    example_label = row_i$example_label,
    anomaly = row_i$anomaly,
    onset_day = row_i$ONSET_mean,
    offset_day = row_i$OFFSET_mean,
    m_count = m_count,
    m_visit = m_visit
  )
})

observed_df <- dplyr::bind_rows(lapply(curve_list, `[[`, "observed"))
predicted_df <- dplyr::bind_rows(lapply(curve_list, `[[`, "predicted"))
boundaries_df <- dplyr::bind_rows(lapply(curve_list, `[[`, "boundaries"))

# ---------------------------------------------------------------------------- #
#### Plot panel B ####
# ---------------------------------------------------------------------------- #

facet_order <- c("Univoltine", "Multivoltine")

predicted_df <- predicted_df |>
  dplyr::mutate(
    example_label = factor(example_label, levels = facet_order)
  )

boundaries_df <- boundaries_df |>
  dplyr::mutate(
    example_label = factor(example_label, levels = facet_order)
  )

fig1b_pheno_anomaly_curves <- ggplot() +
  # Onset and offset for each selected year
  geom_vline(
    data = boundaries_df,
    aes(
      xintercept = x,
      colour = anomaly
    ),
    linetype = "dashed",
    linewidth = 0.5,
    alpha = 0.4,
    show.legend = FALSE
  ) +
  # Fitted GAM curves coloured by thermal anomaly
  geom_line(
    data = predicted_df,
    aes(
      x = julian_day,
      y = predicted_count,
      group = interaction(example_label, YEAR),
      colour = anomaly
    ),
    linewidth = 1,
    alpha = 0.8
  ) +
  facet_wrap(
    ~ example_label,
    ncol = 1,
    scales = "free_y"
  ) +
  scale_colour_gradient2(
    low = "#2166AC",
    mid = "grey85",
    high = "#B2182B",
    midpoint = 0,
    name = "Temperature anomaly (°C)"
  ) +
  scale_x_continuous(
    limits = c(60, 304),
    breaks = c(60, 91, 121, 152, 182, 213, 244, 274),
    labels = c("Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct"),
    expand = c(0.01, 0.01)
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0.12))
  ) +
  labs(
    x = NULL,
    y = "Butterfly count"
  ) +
  theme_classic(
    base_family = "Garamond",
    base_size = 14
  ) +
  theme(
    strip.background = element_blank(),
    strip.text = element_blank(),
    axis.text = element_text(colour = "black"),
    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.3
    ),
    axis.line = element_line(
      colour = "black",
      linewidth = 0.4
    ),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.5
    ),
    panel.spacing = grid::unit(0.2, "cm"),
    
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.justification = "center",
    
    legend.title = element_text(
      face = "bold",
      size = 10,
      hjust = 0.5
    ),
    legend.text = element_text(
      size = 10
    ),
    legend.margin = margin(t = 2, r = 2, b = 2, l = 2),
    
    plot.margin = margin(5, 5, 5, 5)
  ) +
  guides(
    colour = guide_colourbar(
      direction = "horizontal",
      barheight = grid::unit(0.35, "cm"),
      barwidth  = grid::unit(4.0, "cm"),
      title.position = "top",
      title.hjust = 0.5
    )
  )

fig1b_pheno_anomaly_curves
# ---------------------------------------------------------------------------- #
#### Save figure ####
# ---------------------------------------------------------------------------- #

out_fig_dir <- here::here("output", "figures")
dir.create(out_fig_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(
  filename = file.path(out_fig_dir, "fig1b_pheno_curves_temperature_anomaly.png"),
  plot = fig1b_pheno_anomaly_curves,
  width = 5.5,
  height = 7.2,
  dpi = 600
)


# ============================================================================================ #
# End of script
# ============================================================================================ #