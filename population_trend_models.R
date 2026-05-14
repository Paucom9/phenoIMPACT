# ============================================================================================ #
# population_trend_models.R
#
# Author: Pau Colom
# Date: 2026-05-13
#
# Description:
# Test whether phenological plasticity predicts long-term population trends.
#
# Population trends are modelled directly from GAM-derived annual abundance indices:
#
# ABUND_INDEX ~ year * phenological plasticity
#
# Models use Gamma GLMMs with log link, appropriate for positive continuous abundance indices.
#
# Main models:
#   1. Onset advancement plasticity
#   2. Offset termination plasticity in univoltine populations
#   3. Offset delay plasticity in multivoltine populations
# ============================================================================================ #

rm(list = ls())

#### Load required libraries ####

library(dplyr)
library(tidyr)
library(glmmTMB)
library(DHARMa)
library(ggplot2)
library(here)
library(writexl)

#### Data import ####

phenology_estimates <- read.csv(
  here("output", "pheno_abund_estimates_allspp.csv"),
  sep = ",",
  dec = "."
)

phenological_plasticity <- read.csv(
  here("output", "phenological_population_plasticity_contextual.csv"),
  sep = ",",
  dec = "."
)

#### Helper functions ####

zscale <- function(x) {
  as.numeric(scale(x))
}

#### 0. Filter aberrant GAM-derived abundance indices ####

# Exclude extreme GAM-derived abundance values
# These likely correspond to unstable GAM fits producing unrealistic abundance indices

abund_threshold <- quantile(
  phenology_estimates$ABUND_INDEX,
  probs = 0.999,
  na.rm = TRUE
)

abund_threshold

phenology_estimates_clean <- phenology_estimates |>
  dplyr::filter(
    !is.na(ABUND_INDEX),
    ABUND_INDEX > 0,
    ABUND_INDEX <= abund_threshold
  )

# Check how many rows are removed

abundance_filter_summary <- phenology_estimates |>
  dplyr::summarise(
    n_total = dplyr::n(),
    n_removed = sum(
      is.na(ABUND_INDEX) |
        ABUND_INDEX <= 0 |
        ABUND_INDEX > abund_threshold
    ),
    prop_removed = n_removed / n_total
  )

abundance_filter_summary
summary(phenology_estimates_clean$ABUND_INDEX)

#### 1. Prepare annual abundance dataset ####

# For population trends we use all available annual abundance estimates.
# We do not require consecutive years, because we are not calculating lambda.

df_abund <- phenology_estimates_clean |>
  dplyr::distinct(
    SPECIES,
    SITE_ID,
    YEAR,
    ABUND_INDEX
  ) |>
  dplyr::mutate(
    YEAR = as.integer(YEAR),
    SPECIES = as.character(SPECIES),
    SITE_ID = as.character(SITE_ID),
    pop_id = paste(SPECIES, SITE_ID, sep = "_")
  ) |>
  dplyr::filter(
    !is.na(SPECIES),
    !is.na(SITE_ID),
    !is.na(YEAR),
    !is.na(pop_id),
    !is.na(ABUND_INDEX),
    ABUND_INDEX > 0
  ) |>
  dplyr::mutate(
    log_abund_index = log(ABUND_INDEX),
    year_z = zscale(YEAR),
    SPECIES = as.factor(SPECIES),
    SITE_ID = as.factor(SITE_ID),
    YEAR = as.factor(YEAR),
    pop_id = as.factor(pop_id)
  )

summary(df_abund$ABUND_INDEX)
summary(df_abund$log_abund_index)

#### 2. Onset advancement plasticity and population trends ####

df_trend_onset <- df_abund |>
  dplyr::left_join(
    phenological_plasticity |>
      dplyr::select(
        pop_id,
        voltinism,
        onset_advancement_plasticity_contextual
      ),
    by = "pop_id"
  ) |>
  dplyr::filter(
    !is.na(ABUND_INDEX),
    !is.na(year_z),
    !is.na(onset_advancement_plasticity_contextual),
    !is.na(voltinism)
  ) |>
  dplyr::mutate(
    voltinism = as.factor(voltinism),
    onset_plasticity_z = zscale(onset_advancement_plasticity_contextual)
  )

summary(df_trend_onset$ABUND_INDEX)
summary(df_trend_onset$onset_advancement_plasticity_contextual)
summary(df_trend_onset$onset_plasticity_z)

mod_trend_onset_gamma <- glmmTMB(
  ABUND_INDEX ~ year_z * onset_plasticity_z +
    (1 + year_z || SPECIES) +
    (1 + year_z || SITE_ID),
  data = df_trend_onset,
  family = Gamma(link = "log")
)

summary(mod_trend_onset_gamma)

#### 3. Offset plasticity and population trends ####

# Offset plasticity has different biological meanings depending on voltinism:
#
# Univoltine populations:
#   warmer years may advance offset / terminate the flight period earlier.
#   Therefore, higher biologically oriented plasticity = stronger advancement of offset.
#
# Multivoltine populations:
#   warmer years may delay offset / extend the flight period.
#   Therefore, higher biologically oriented plasticity = stronger delay of offset.
#
# To make "high plasticity" biologically meaningful in both groups:
#   univoltine   = -offset_termination_plasticity_contextual
#   multivoltine =  offset_termination_plasticity_contextual

df_trend_offset <- df_abund |>
  dplyr::left_join(
    phenological_plasticity |>
      dplyr::select(
        pop_id,
        voltinism,
        offset_termination_plasticity_contextual
      ),
    by = "pop_id"
  ) |>
  dplyr::filter(
    !is.na(ABUND_INDEX),
    !is.na(year_z),
    !is.na(offset_termination_plasticity_contextual),
    !is.na(voltinism)
  ) |>
  dplyr::mutate(
    voltinism = as.factor(voltinism),
    offset_plasticity_bio = dplyr::case_when(
      voltinism == "univoltine" ~ -offset_termination_plasticity_contextual,
      voltinism == "multivoltine" ~  offset_termination_plasticity_contextual,
      TRUE ~ NA_real_
    )
  ) |>
  dplyr::filter(
    !is.na(offset_plasticity_bio)
  ) |>
  dplyr::group_by(voltinism) |>
  dplyr::mutate(
    offset_plasticity_bio_z = zscale(offset_plasticity_bio)
  ) |>
  dplyr::ungroup()

offset_plasticity_summary <- df_trend_offset |>
  dplyr::group_by(voltinism) |>
  dplyr::summarise(
    n = dplyr::n(),
    n_species = dplyr::n_distinct(SPECIES),
    n_sites = dplyr::n_distinct(SITE_ID),
    min_raw = min(offset_termination_plasticity_contextual, na.rm = TRUE),
    mean_raw = mean(offset_termination_plasticity_contextual, na.rm = TRUE),
    max_raw = max(offset_termination_plasticity_contextual, na.rm = TRUE),
    min_bio = min(offset_plasticity_bio, na.rm = TRUE),
    mean_bio = mean(offset_plasticity_bio, na.rm = TRUE),
    max_bio = max(offset_plasticity_bio, na.rm = TRUE),
    .groups = "drop"
  )

offset_plasticity_summary

#### 4. Univoltine offset model ####

df_trend_offset_uni <- df_trend_offset |>
  dplyr::filter(voltinism == "univoltine") |>
  droplevels()

summary(df_trend_offset_uni$ABUND_INDEX)
summary(df_trend_offset_uni$offset_plasticity_bio)
summary(df_trend_offset_uni$offset_plasticity_bio_z)

mod_trend_offset_uni_gamma <- glmmTMB(
  ABUND_INDEX ~ year_z * offset_plasticity_bio_z +
    (1 + year_z || SPECIES) +
    (1 + year_z || SITE_ID),
  data = df_trend_offset_uni,
  family = Gamma(link = "log")
)

summary(mod_trend_offset_uni_gamma)

#### 5. Multivoltine offset model ####

df_trend_offset_multi <- df_trend_offset |>
  dplyr::filter(voltinism == "multivoltine") |>
  droplevels()

summary(df_trend_offset_multi$ABUND_INDEX)
summary(df_trend_offset_multi$offset_plasticity_bio)
summary(df_trend_offset_multi$offset_plasticity_bio_z)

mod_trend_offset_multi_gamma <- glmmTMB(
  ABUND_INDEX ~ year_z * offset_plasticity_bio_z +
    (1 + year_z || SPECIES) +
    (1 + year_z || SITE_ID),
  data = df_trend_offset_multi,
  family = Gamma(link = "log")
)

summary(mod_trend_offset_multi_gamma)


#### 6. Function to plot relative abundance change from Gamma models ####

plot_gamma_percent_change <- function(model,
                                      plasticity_var,
                                      x_min = -2,
                                      x_max = 2,
                                      x_by = 0.1,
                                      title = NULL) {
  
  b <- glmmTMB::fixef(model)$cond
  V <- as.matrix(vcov(model)$cond)
  
  interaction_name <- paste("year_z", plasticity_var, sep = ":")
  
  if (!interaction_name %in% names(b)) {
    interaction_name <- paste(plasticity_var, "year_z", sep = ":")
  }
  
  if (!interaction_name %in% names(b)) {
    stop("Could not find interaction term between year_z and plasticity variable.")
  }
  
  x_vals <- seq(x_min, x_max, by = x_by)
  plast_vals <- c(1, 0, -1)
  x0 <- min(x_vals)
  
  pred <- expand.grid(
    year_z = x_vals,
    plasticity_z = plast_vals
  ) |>
    dplyr::mutate(
      group = factor(
        plasticity_z,
        levels = c(1, 0, -1),
        labels = c("High", "Mean", "Low")
      ),
      delta_year = year_z - x0,
      
      # Slope on the log scale
      slope = b["year_z"] +
        b[interaction_name] * plasticity_z,
      
      # Relative change on the log scale
      log_rel = delta_year * slope,
      
      # Uncertainty of the slope
      var_slope =
        V["year_z", "year_z"] +
        plasticity_z^2 * V[interaction_name, interaction_name] +
        2 * plasticity_z * V["year_z", interaction_name],
      
      se_log_rel = abs(delta_year) * sqrt(var_slope),
      
      low_log_rel = log_rel - 1.96 * se_log_rel,
      high_log_rel = log_rel + 1.96 * se_log_rel,
      
      # Convert to percent change
      percent_change = (exp(log_rel) - 1) * 100,
      low_percent = (exp(low_log_rel) - 1) * 100,
      high_percent = (exp(high_log_rel) - 1) * 100
    )
  
  p <- ggplot(
    pred,
    aes(
      x = year_z,
      y = percent_change,
      colour = group,
      fill = group
    )
  ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      colour = "grey40"
    ) +
    geom_ribbon(
      aes(
        ymin = low_percent,
        ymax = high_percent
      ),
      alpha = 0.10,
      colour = NA
    ) +
    geom_line(linewidth = 1.4) +
    scale_colour_manual(
      breaks = c("High", "Mean", "Low"),
      values = c(
        "High" = "#0072B2",
        "Mean" = "#009E73",
        "Low" = "#D55E00"
      )
    ) +
    scale_fill_manual(
      breaks = c("High", "Mean", "Low"),
      values = c(
        "High" = "#0072B2",
        "Mean" = "#009E73",
        "Low" = "#D55E00"
      )
    ) +
    labs(
      x = "Year (standardized)",
      y = "Relative change in abundance (%)",
      colour = "Plasticity",
      fill = "Plasticity",
      title = title
    ) +
    theme_classic(base_size = 15, base_family = "Garamond")
  
  return(p)
}

#### 7. Plot model effects ####

plot_trend_onset <- plot_gamma_percent_change(
  model = mod_trend_onset_gamma,
  plasticity_var = "onset_plasticity_z",
  title = "Onset advancement plasticity"
)

plot_trend_offset_uni <- plot_gamma_percent_change(
  model = mod_trend_offset_uni_gamma,
  plasticity_var = "offset_plasticity_bio_z",
  title = "Offset termination plasticity: univoltine populations"
)

plot_trend_offset_multi <- plot_gamma_percent_change(
  model = mod_trend_offset_multi_gamma,
  plasticity_var = "offset_plasticity_bio_z",
  title = "Offset delay plasticity: multivoltine populations"
)

plot_trend_onset
plot_trend_offset_uni
plot_trend_offset_multi

#### 8. Save outputs ####

saveRDS(
  mod_trend_onset_gamma,
  here("output", "model_trend_onset_gamma.rds")
)

saveRDS(
  mod_trend_offset_uni_gamma,
  here("output", "model_trend_offset_univoltine_gamma.rds")
)

saveRDS(
  mod_trend_offset_multi_gamma,
  here("output", "model_trend_offset_multivoltine_gamma.rds")
)

saveRDS(
  res_trend_onset_gamma,
  here("output", "diagnostics_trend_onset_gamma_DHARMa.rds")
)

saveRDS(
  res_trend_offset_uni_gamma,
  here("output", "diagnostics_trend_offset_univoltine_gamma_DHARMa.rds")
)

saveRDS(
  res_trend_offset_multi_gamma,
  here("output", "diagnostics_trend_offset_multivoltine_gamma_DHARMa.rds")
)

ggsave(
  filename = here("output", "figures", "plot_trend_onset_gamma.png"),
  plot = plot_trend_onset,
  width = 7,
  height = 5,
  dpi = 300
)

ggsave(
  filename = here("output", "figures", "plot_trend_offset_univoltine_gamma.png"),
  plot = plot_trend_offset_uni,
  width = 7,
  height = 5,
  dpi = 300
)

ggsave(
  filename = here("output", "figures", "plot_trend_offset_multivoltine_gamma.png"),
  plot = plot_trend_offset_multi,
  width = 7,
  height = 5,
  dpi = 300
)

write.csv(
  df_trend_onset,
  here("output", "df_population_trend_onset.csv"),
  row.names = FALSE
)

write.csv(
  df_trend_offset,
  here("output", "df_population_trend_offset.csv"),
  row.names = FALSE
)

write.csv(
  offset_plasticity_summary,
  here("output", "summary_offset_plasticity_by_voltinism.csv"),
  row.names = FALSE
)

writexl::write_xlsx(
  list(
    abundance_filter_summary = abundance_filter_summary,
    df_abund = df_abund,
    df_trend_onset = df_trend_onset,
    df_trend_offset = df_trend_offset,
    df_trend_offset_uni = df_trend_offset_uni,
    df_trend_offset_multi = df_trend_offset_multi,
    offset_plasticity_summary = offset_plasticity_summary
  ),
  here("output", "population_trend_model_datasets.xlsx")
)
