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
# ABUND_INDEX ~ year_decade * onset plasticity + year_decade * offset plasticity
#
# Models use Gamma GLMMs with log link, appropriate for positive continuous abundance indices.
#
# Main models:
#   1. Combined onset advancement + offset termination plasticity in univoltine populations
#   2. Combined onset advancement + offset delay plasticity in multivoltine populations
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
library(extrafont)
library(performance)
library(tibble)

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

abund_threshold <- quantile(
  phenology_estimates$ABUND_INDEX,
  probs = 0.999,
  na.rm = TRUE
)

phenology_estimates_clean <- phenology_estimates |>
  dplyr::filter(
    !is.na(ABUND_INDEX),
    ABUND_INDEX > 0,
    ABUND_INDEX <= abund_threshold
  )

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

df_abund_raw <- phenology_estimates_clean |>
  dplyr::distinct(
    SPECIES,
    SITE_ID,
    YEAR,
    ABUND_INDEX
  ) |>
  dplyr::mutate(
    year_num = as.integer(YEAR),
    SPECIES = as.character(SPECIES),
    SITE_ID = as.character(SITE_ID),
    pop_id = paste(SPECIES, SITE_ID, sep = "_")
  ) |>
  dplyr::filter(
    !is.na(SPECIES),
    !is.na(SITE_ID),
    !is.na(year_num),
    !is.na(pop_id),
    !is.na(ABUND_INDEX),
    ABUND_INDEX > 0
  )

year_center <- mean(df_abund_raw$year_num, na.rm = TRUE)

df_abund <- df_abund_raw |>
  dplyr::mutate(
    log_abund_index = log(ABUND_INDEX),
    year_decade = (year_num - year_center) / 10,
    SPECIES = as.factor(SPECIES),
    SITE_ID = as.factor(SITE_ID),
    YEAR = as.factor(year_num),
    pop_id = as.character(pop_id)
  )

summary(df_abund$ABUND_INDEX)
summary(df_abund$log_abund_index)

#### 2. Prepare population-level plasticity dataset ####

# Offset plasticity has different biological meanings depending on voltinism:
#
# Univoltine populations:
#   higher biological offset plasticity = stronger advancement / earlier termination.
#
# Multivoltine populations:
#   higher biological offset plasticity = stronger delay / late-season extension.
#
# To make "high offset plasticity" biologically meaningful in both groups:
#   univoltine   = -offset_termination_plasticity_contextual
#   multivoltine =  offset_termination_plasticity_contextual

plasticity_pop <- phenological_plasticity |>
  dplyr::select(
    pop_id,
    voltinism,
    onset_advancement_plasticity_contextual,
    offset_termination_plasticity_contextual
  ) |>
  dplyr::mutate(
    pop_id = as.character(pop_id),
    voltinism = as.factor(voltinism),
    offset_plasticity_bio = dplyr::case_when(
      voltinism == "univoltine" ~ -offset_termination_plasticity_contextual,
      voltinism == "multivoltine" ~  offset_termination_plasticity_contextual,
      TRUE ~ NA_real_
    )
  ) |>
  dplyr::filter(
    !is.na(pop_id),
    !is.na(voltinism),
    !is.na(onset_advancement_plasticity_contextual),
    !is.na(offset_termination_plasticity_contextual),
    !is.na(offset_plasticity_bio)
  ) |>
  dplyr::distinct(
    pop_id,
    voltinism,
    onset_advancement_plasticity_contextual,
    offset_termination_plasticity_contextual,
    offset_plasticity_bio
  ) |>
  dplyr::group_by(voltinism) |>
  dplyr::mutate(
    onset_plasticity_z = zscale(onset_advancement_plasticity_contextual),
    offset_plasticity_bio_z = zscale(offset_plasticity_bio)
  ) |>
  dplyr::ungroup()

#### 3. Check correlation between onset and offset plasticity ####

plasticity_correlations <- plasticity_pop |>
  dplyr::group_by(voltinism) |>
  dplyr::summarise(
    n_populations = dplyr::n(),
    
    pearson_raw = cor(
      onset_advancement_plasticity_contextual,
      offset_termination_plasticity_contextual,
      use = "complete.obs",
      method = "pearson"
    ),
    
    spearman_raw = cor(
      onset_advancement_plasticity_contextual,
      offset_termination_plasticity_contextual,
      use = "complete.obs",
      method = "spearman"
    ),
    
    pearson_bio = cor(
      onset_advancement_plasticity_contextual,
      offset_plasticity_bio,
      use = "complete.obs",
      method = "pearson"
    ),
    
    spearman_bio = cor(
      onset_advancement_plasticity_contextual,
      offset_plasticity_bio,
      use = "complete.obs",
      method = "spearman"
    ),
    
    .groups = "drop"
  )

plasticity_correlations

plot_plasticity_correlation <- ggplot(
  plasticity_pop,
  aes(
    x = onset_plasticity_z,
    y = offset_plasticity_bio_z
  )
) +
  geom_point(alpha = 0.30) +
  geom_smooth(method = "lm", se = TRUE) +
  facet_wrap(~ voltinism) +
  theme_classic(base_size = 15, base_family = "Garamond") +
  labs(
    x = "Onset advancement plasticity",
    y = "Offset plasticity",
    title = "Correlation between onset and offset plasticity"
  )

plot_plasticity_correlation

#### 4. Create combined trend dataset ####

df_trend_combined <- df_abund |>
  dplyr::left_join(
    plasticity_pop,
    by = "pop_id"
  ) |>
  dplyr::filter(
    !is.na(ABUND_INDEX),
    !is.na(year_decade),
    !is.na(onset_plasticity_z),
    !is.na(offset_plasticity_bio_z),
    !is.na(voltinism)
  ) |>
  dplyr::mutate(
    SPECIES = as.factor(SPECIES),
    SITE_ID = as.factor(SITE_ID),
    pop_id = as.factor(pop_id),
    voltinism = as.factor(voltinism),
        year_fac = factor(
      year_num,
      levels = min(year_num, na.rm = TRUE):max(year_num, na.rm = TRUE)
    )
  )

combined_data_summary <- df_trend_combined |>
  dplyr::group_by(voltinism) |>
  dplyr::summarise(
    n_observations = dplyr::n(),
    n_populations = dplyr::n_distinct(pop_id),
    n_species = dplyr::n_distinct(SPECIES),
    n_sites = dplyr::n_distinct(SITE_ID),
    min_year = min(year_num, na.rm = TRUE),
    max_year = max(year_num, na.rm = TRUE),
    .groups = "drop"
  )

combined_data_summary

#### 5. Split by voltinism ####

df_trend_combined_uni <- df_trend_combined |>
  dplyr::filter(voltinism == "univoltine") |>
  droplevels()

df_trend_combined_multi <- df_trend_combined |>
  dplyr::filter(voltinism == "multivoltine") |>
  droplevels()

#### 6. Main model: univoltine populations ####

mod_trend_combined_uni_gamma <- glmmTMB(
  ABUND_INDEX ~
    year_decade * onset_plasticity_z +
    year_decade * offset_plasticity_bio_z +
    (1 + year_decade || SPECIES) +
    (1 + year_decade || SITE_ID) +
    ar1(year_fac + 0 | pop_id),
  data = df_trend_combined_uni,
  family = Gamma(link = "log")
)

summary(mod_trend_combined_uni_gamma)

#### 7. Main model: multivoltine populations ####

mod_trend_combined_multi_gamma <- glmmTMB(
  ABUND_INDEX ~
    year_decade * onset_plasticity_z +
    year_decade * offset_plasticity_bio_z +
    (1 + year_decade || SPECIES) +
    (1 + year_decade || SITE_ID) +
    ar1(year_fac + 0 | pop_id),
  data = df_trend_combined_multi,
  family = Gamma(link = "log")
)

summary(mod_trend_combined_multi_gamma)

#### 8. Check multicollinearity ####

collinearity_uni <- performance::check_collinearity(mod_trend_combined_uni_gamma)
collinearity_multi <- performance::check_collinearity(mod_trend_combined_multi_gamma)

collinearity_uni
collinearity_multi

#### 9. Extract key model terms ####

extract_key_terms <- function(model) {
  
  coef_tab <- as.data.frame(summary(model)$coefficients$cond)
  
  coef_tab |>
    tibble::rownames_to_column("term") |>
    dplyr::filter(
      term %in% c(
        "year_decade",
        "onset_plasticity_z",
        "offset_plasticity_bio_z",
        "year_decade:onset_plasticity_z",
        "year_decade:offset_plasticity_bio_z",
        "onset_plasticity_z:year_decade",
        "offset_plasticity_bio_z:year_decade"
      )
    ) |>
    dplyr::mutate(
      percent_change_per_decade = dplyr::if_else(
        grepl("year_decade", term),
        (exp(Estimate) - 1) * 100,
        NA_real_
      )
    )
}

key_terms_uni <- extract_key_terms(mod_trend_combined_uni_gamma) |>
  dplyr::mutate(model = "univoltine")

key_terms_multi <- extract_key_terms(mod_trend_combined_multi_gamma) |>
  dplyr::mutate(model = "multivoltine")

key_terms_combined <- dplyr::bind_rows(
  key_terms_uni,
  key_terms_multi
)

key_terms_combined

#### 10. Function to plot relative abundance change from Gamma models ####

# This plots the effect of one plasticity variable while holding the other plasticity variable at its mean,
# because both plasticity variables are z-scaled and therefore mean = 0.

plot_gamma_percent_change <- function(model,
                                      data,
                                      plasticity_var,
                                      year_var = "year_decade",
                                      year_num_var = "year_num",
                                      year_center,
                                      title = NULL) {
  
  b <- glmmTMB::fixef(model)$cond
  V <- as.matrix(vcov(model)$cond)
  
  interaction_name <- paste(year_var, plasticity_var, sep = ":")
  
  if (!interaction_name %in% names(b)) {
    interaction_name <- paste(plasticity_var, year_var, sep = ":")
  }
  
  if (!interaction_name %in% names(b)) {
    stop("Could not find interaction term.")
  }
  
  year_vals <- seq(
    min(data[[year_num_var]], na.rm = TRUE),
    max(data[[year_num_var]], na.rm = TRUE),
    by = 1
  )
  
  year_df <- tibble::tibble(
    year_num = year_vals
  )
  
  year_df[[year_var]] <- (year_df$year_num - year_center) / 10
  
  x0 <- min(year_df[[year_var]], na.rm = TRUE)
  
  pred <- tidyr::expand_grid(
    year_df,
    plasticity_z = c(1, 0, -1)
  ) |>
    dplyr::mutate(
      group = factor(
        plasticity_z,
        levels = c(1, 0, -1),
        labels = c("High (+1 SD)", "Mean", "Low (-1 SD)")
      ),
      
      delta_year = .data[[year_var]] - x0,
      
      slope = b[year_var] +
        b[interaction_name] * plasticity_z,
      
      log_rel = delta_year * slope,
      
      var_slope =
        V[year_var, year_var] +
        plasticity_z^2 * V[interaction_name, interaction_name] +
        2 * plasticity_z * V[year_var, interaction_name],
      
      se_log_rel = abs(delta_year) * sqrt(var_slope),
      
      low_log_rel  = log_rel - 1.96 * se_log_rel,
      high_log_rel = log_rel + 1.96 * se_log_rel,
      
      percent_change = (exp(log_rel) - 1) * 100,
      low_percent    = (exp(low_log_rel) - 1) * 100,
      high_percent   = (exp(high_log_rel) - 1) * 100
    ) |>
    dplyr::arrange(group, year_num)
  
  ggplot(
    pred,
    aes(
      x = year_num,
      y = percent_change,
      colour = group,
      fill = group,
      group = group
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
      values = c(
        "High (+1 SD)" = "#009E73",
        "Mean" = "#0072B2",
        "Low (-1 SD)" = "#D55E00"
      )
    ) +
    scale_fill_manual(
      values = c(
        "High (+1 SD)" = "#009E73",
        "Mean" = "#0072B2",
        "Low (-1 SD)" = "#D55E00"
      )
    ) +
    labs(
      x = "Year",
      y = "Relative change in abundance (%)",
      colour = "Plasticity",
      fill = "Plasticity",
      title = title
    ) +
    theme_classic(base_size = 15, base_family = "Garamond")
}

#### 11. Plot model effects ####

plot_trend_onset_uni <- plot_gamma_percent_change(
  model = mod_trend_combined_uni_gamma,
  data = df_trend_combined_uni,
  plasticity_var = "onset_plasticity_z",
  year_center = year_center,
  title = "Onset advancement plasticity: univoltine populations"
)

plot_trend_offset_uni <- plot_gamma_percent_change(
  model = mod_trend_combined_uni_gamma,
  data = df_trend_combined_uni,
  plasticity_var = "offset_plasticity_bio_z",
  year_center = year_center,
  title = "Offset termination plasticity: univoltine populations"
)

plot_trend_onset_multi <- plot_gamma_percent_change(
  model = mod_trend_combined_multi_gamma,
  data = df_trend_combined_multi,
  plasticity_var = "onset_plasticity_z",
  year_center = year_center,
  title = "Onset advancement plasticity: multivoltine populations"
)

plot_trend_offset_multi <- plot_gamma_percent_change(
  model = mod_trend_combined_multi_gamma,
  data = df_trend_combined_multi,
  plasticity_var = "offset_plasticity_bio_z",
  year_center = year_center,
  title = "Offset delay plasticity: multivoltine populations"
)

plot_trend_onset_uni
plot_trend_offset_uni
plot_trend_onset_multi
plot_trend_offset_multi


#### 12. Plot distribution population trends ####

get_model_based_population_trends <- function(model, data, group_name) {
  
  b <- glmmTMB::fixef(model)$cond
  
  int_onset <- "year_decade:onset_plasticity_z"
  if (!int_onset %in% names(b)) {
    int_onset <- "onset_plasticity_z:year_decade"
  }
  
  int_offset <- "year_decade:offset_plasticity_bio_z"
  if (!int_offset %in% names(b)) {
    int_offset <- "offset_plasticity_bio_z:year_decade"
  }
  
  sp_re <- as.data.frame(glmmTMB::ranef(model)$cond$SPECIES) |>
    tibble::rownames_to_column("SPECIES")
  
  site_re <- as.data.frame(glmmTMB::ranef(model)$cond$SITE_ID) |>
    tibble::rownames_to_column("SITE_ID")
  
  sp_year_col <- grep("year_decade", names(sp_re), value = TRUE)[1]
  site_year_col <- grep("year_decade", names(site_re), value = TRUE)[1]
  
  sp_re <- sp_re |>
    dplyr::transmute(
      SPECIES = as.character(SPECIES),
      species_year_re = .data[[sp_year_col]]
    )
  
  site_re <- site_re |>
    dplyr::transmute(
      SITE_ID = as.character(SITE_ID),
      site_year_re = .data[[site_year_col]]
    )
  
  data |>
    dplyr::distinct(
      pop_id,
      SPECIES,
      SITE_ID,
      onset_plasticity_z,
      offset_plasticity_bio_z
    ) |>
    dplyr::mutate(
      SPECIES = as.character(SPECIES),
      SITE_ID = as.character(SITE_ID)
    ) |>
    dplyr::left_join(sp_re, by = "SPECIES") |>
    dplyr::left_join(site_re, by = "SITE_ID") |>
    dplyr::mutate(
      voltinism = group_name,
      
      beta_log_per_decade =
        b["year_decade"] +
        b[int_onset] * onset_plasticity_z +
        b[int_offset] * offset_plasticity_bio_z +
        dplyr::coalesce(species_year_re, 0) +
        dplyr::coalesce(site_year_re, 0),
      
      trend_percent_decade =
        (exp(beta_log_per_decade) - 1) * 100
    )
}

model_trends_uni <- get_model_based_population_trends(
  model = mod_trend_combined_uni_gamma,
  data = df_trend_combined_uni,
  group_name = "univoltine"
)

model_trends_multi <- get_model_based_population_trends(
  model = mod_trend_combined_multi_gamma,
  data = df_trend_combined_multi,
  group_name = "multivoltine"
)

model_trends <- dplyr::bind_rows(
  model_trends_uni,
  model_trends_multi
)

model_trends |>
  dplyr::group_by(voltinism) |>
  dplyr::summarise(
    n_populations = dplyr::n(),
    mean_trend = mean(trend_percent_decade, na.rm = TRUE),
    median_trend = median(trend_percent_decade, na.rm = TRUE),
    q25 = quantile(trend_percent_decade, 0.25, na.rm = TRUE),
    q75 = quantile(trend_percent_decade, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

model_trends <- model_trends |>
  dplyr::mutate(
    voltinism = dplyr::case_when(
      voltinism == "univoltine" ~ "Univoltine",
      voltinism == "multivoltine" ~ "Multivoltine",
      TRUE ~ as.character(voltinism)
    )
  )

trend_summary <- model_trends |>
  dplyr::group_by(voltinism) |>
  dplyr::summarise(
    median_trend = median(trend_percent_decade, na.rm = TRUE),
    prop_positive = mean(trend_percent_decade > 0, na.rm = TRUE) * 100,
    .groups = "drop"
  )

xlims <- quantile(
  model_trends$trend_percent_decade,
  probs = c(0.01, 0.99),
  na.rm = TRUE
)

plot_model_trend_hist_overlap <- ggplot(
  model_trends,
  aes(
    x = trend_percent_decade,
    fill = voltinism,
    colour = voltinism
  )
) +
  geom_histogram(
    aes(y = after_stat(density)),
    binwidth = 5,
    position = "identity",
    alpha = 0.35,
    colour = "white"
  ) +
  geom_density(
    linewidth = 1.1,
    alpha = 0
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = "grey30"
  ) +
  geom_vline(
    data = trend_summary,
    aes(xintercept = median_trend, colour = voltinism),
    linewidth = 1.1,
    linetype = "dashed"
  ) +
  coord_cartesian(xlim = xlims) +
  scale_colour_manual(
    breaks = c("Univoltine", "Multivoltine"),
    values = c(
      "Univoltine" = "#2F5D62",
      "Multivoltine" = "#B85C38"
    )
  ) +
  scale_fill_manual(
    breaks = c("Univoltine", "Multivoltine"),
    values = c(
      "Univoltine" = "#2F5D62",
      "Multivoltine" = "#B85C38"
    )
  ) +
  theme_classic(base_size = 15, base_family = "Garamond") +
  labs(
    x = "Abundance change (% per decade)",
    y = "Density",
    colour = "",
    fill = "",
    title = ""
  ) +
  theme(
    legend.position = "top"
  )
  

plot_model_trend_hist_overlap

#### 12. Save outputs ####

dir.create(
  here("output", "figures"),
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(
  mod_trend_combined_uni_gamma,
  here("output", "model_trend_combined_univoltine_gamma.rds")
)

saveRDS(
  mod_trend_combined_multi_gamma,
  here("output", "model_trend_combined_multivoltine_gamma.rds")
)

saveRDS(
  collinearity_uni,
  here("output", "collinearity_combined_univoltine.rds")
)

saveRDS(
  collinearity_multi,
  here("output", "collinearity_combined_multivoltine.rds")
)


ggsave(
  filename = here("output", "figures", "plot_density_pop_trends_univol_multivol.png"),
  plot = plot_model_trend_hist_overlap,
  width = 7,
  height = 5,
  dpi = 300
)

ggsave(
  filename = here("output", "figures", "plot_trend_onset_univoltine_gamma.png"),
  plot = plot_trend_onset_uni,
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
  filename = here("output", "figures", "plot_trend_onset_multivoltine_gamma.png"),
  plot = plot_trend_onset_multi,
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

ggsave(
  filename = here("output", "figures", "plot_plasticity_correlation.png"),
  plot = plot_plasticity_correlation,
  width = 7,
  height = 5,
  dpi = 300
)

write.csv(
  df_trend_combined,
  here("output", "df_population_trend_combined.csv"),
  row.names = FALSE
)

write.csv(
  plasticity_correlations,
  here("output", "summary_plasticity_correlations.csv"),
  row.names = FALSE
)

write.csv(
  combined_data_summary,
  here("output", "summary_combined_population_trend_dataset.csv"),
  row.names = FALSE
)

write.csv(
  key_terms_combined,
  here("output", "summary_combined_population_trend_model_terms.csv"),
  row.names = FALSE
)

writexl::write_xlsx(
  list(
    abundance_filter_summary = abundance_filter_summary,
    combined_data_summary = combined_data_summary,
    plasticity_correlations = plasticity_correlations,
    key_terms_combined = key_terms_combined,
    df_abund = df_abund,
    plasticity_pop = plasticity_pop,
    df_trend_combined = df_trend_combined,
    df_trend_combined_uni = df_trend_combined_uni,
    df_trend_combined_multi = df_trend_combined_multi
  ),
  here("output", "population_trend_combined_model_datasets.xlsx")
)
