# ============================================================================================ #
# pheno_trend_models.R
#
# Author: Pau Colom
# Date: 2026-02-20
#
# Description: This script analyzes phenological trends across different sites and species 
# by fitting linear mixed models. It explores the distribution of these trends, their 
# correlations, and how they relate to climatic regions, latitude, and temperature. 
# The script also includes visualizations to illustrate the findings.
#
# ============================================================================================ #

#### Load required libraries ####
# ---
library(dplyr) # For data manipulation
library(tidyr) # For data tidying
library(broom) # For tidying model outputs
library(purrr) # For functional programming
library(nlme) # For generalized least squares models
library(here) # For file path management
library(reshape2)
library(ggplot2)
library(lme4)
library(glmmTMB)
library(broom.mixed)
library(emmeans)
library(ggeffects)
library(patchwork)
library(mgcv)
library(sf)
library(rnaturalearth)
library(MuMIn)
library(tibble)
library(lmerTest)
library(extrafont)
# ---

#### Data Import and Preparation ####
# ---
here::here() # Check the current working directory
# ---
pheno_trends_site  <- read.csv(here("output", "pheno_temporal_trends_allspp.csv"), sep = ",", dec = ".")
clim_vars <- read.csv(here::here("output", "climate", "climate_variables.csv"), sep = ",", dec = ".")
ebms_clim_df   <- read.csv(here("data", "ebms_transect_climate.csv"), sep = ",", dec = ".")
ebms_coord_df  <- read.csv(here("data", "ebms_transect_coord.csv"), sep = ",", dec = ".")
species_traits <- read.csv(here::here("data", "species_trait_table.csv"), sep = ";", dec = ".")

str(pheno_trends_site)
str(clim_vars)
str(ebms_clim_df)
str(ebms_coord_df)


# Merge datasets

voltinism_spp <- species_traits %>%
  dplyr::transmute(
    SPECIES = gsub("_", " ", Taxon),
    voltinism = dplyr::case_when(
      Vol_max <= 1.5 ~ "univoltine",
      Vol_max >= 2   ~ "multivoltine",
      TRUE ~ NA_character_
    )
  )

# --- Climate (1 row per SITE × SPECIES) ---
clim_vars_sp <- clim_vars |>
  distinct(SITE_ID, SPECIES, .keep_all = TRUE)

# --- Merge everything ---
df <- pheno_trends_site |>
  
  # climate
  left_join(
    clim_vars_sp |>
      dplyr::select(
        SITE_ID,
        SPECIES,
        clim_background_temp_90_sc,
        clim_trend_temp_90_sc,
        clim_stability_temp_90_sc,
        clim_autocorr_temp_90_sc,
        photo_90_sc
      ),
    by = c("SITE_ID", "SPECIES")
  ) |>
  
  # climate region
  left_join(
    ebms_clim_df,
    by = c("SITE_ID" = "transect_id")
  ) |>
  
  # coordinates
  left_join(
    ebms_coord_df,
    by = c("SITE_ID" = "transect_id")
  ) |>
  
  # voltinism
  left_join(
    voltinism_spp,
    by = "SPECIES"
  ) |>
  
  mutate(
    voltinism = factor(voltinism,
                       levels = c("univoltine", "multivoltine"))
  )

str(df)
# ---


#### Models ####

# Fit a predefined set of hypothesis-driven models
# H1: warming only
# H2: warming × background climate
# H3: warming × voltinism
# H4: warming × background + warming × voltinism

model_set <- function(d) {
  
  # ---- remove missing values (required for lmer) ----
  d <- d %>%
    dplyr::filter(
      !is.na(estimate),
      !is.na(std.error),
      !is.na(clim_background_temp_90_sc),
      !is.na(clim_trend_temp_90_sc),
      !is.na(voltinism)
    )
  
  list(
    
    # H1: warming only
    m1 = lmer(
      estimate ~ clim_trend_temp_90_sc +
        (1 |  SPECIES) +  # random slope for warming
        (1 | SITE_ID),                    # spatial grouping
      data = d,
      weights = 1 / (std.error^2),        # meta-analytic weighting
      REML = FALSE,
      control = lme4::lmerControl(optimizer = "bobyqa")
    ),
    
    
    # H2: warming × background climate
    m2 = lmer(
      estimate ~ clim_trend_temp_90_sc * clim_background_temp_90_sc +
        (1 |  SPECIES) +  # random slope for warming
        (1 | SITE_ID),
      data = d,
      weights = 1 / (std.error^2),
      REML = FALSE,
      control = lme4::lmerControl(optimizer = "bobyqa")
    ),
    
    # H3: warming × voltinism
    m3 = lmer(
      estimate ~ clim_trend_temp_90_sc * voltinism +
        (1 |  SPECIES) +  # random slope for warming
        (1 | SITE_ID),
      data = d,
      weights = 1 / (std.error^2),
      REML = FALSE,
      control = lme4::lmerControl(optimizer = "bobyqa")
    ),
    
    # H4: full model (both moderators)
    m4 = lmer(
      estimate ~ clim_trend_temp_90_sc * clim_background_temp_90_sc +
        clim_trend_temp_90_sc * voltinism +
        (1 |  SPECIES) +  # random slope for warming
        (1 | SITE_ID),
      data = d,
      weights = 1 / (std.error^2),
      REML = FALSE,
      control = lme4::lmerControl(optimizer = "bobyqa")
    )
  )
}


#### Model selection (ΔAIC ≤ 2 retained) ####

d_onset <- df %>%
  dplyr::filter(phenovar == "ONSET_mean")

mods <- model_set(d_onset)

sel <- MuMIn::model.sel(mods)

best_models <- as.data.frame(sel) %>%
  tibble::rownames_to_column("model") %>%
  dplyr::filter(delta <= 2)


#### Extract coefficients from selected models ####

table_final <- purrr::map_dfr(seq_len(nrow(best_models)), function(i) {
  
  mod_name <- best_models$model[i]
  mod      <- mods[[mod_name]]
  
  coefs <- broom.mixed::tidy(mod, effects = "fixed")
  
  coefs %>%
    dplyr::select(term, estimate, std.error, p.value) %>%
    dplyr::mutate(
      model = mod_name,
      delta = best_models$delta[i],
      weight = best_models$weight[i]
    ) %>%
    tidyr::pivot_wider(
      names_from = term,
      values_from = c(estimate, std.error, p.value)
    )
})


#### Rename predictors (clean table labels) ####

table_final <- table_final %>%
  dplyr::rename(
    Warming = estimate_clim_trend_temp_90_sc,
    Background = estimate_clim_background_temp_90_sc,
    `Warming×Background` = `estimate_clim_trend_temp_90_sc:clim_background_temp_90_sc`,
    Voltinism = estimate_voltinismmultivoltine,
    `Warming×Voltinism` = `estimate_clim_trend_temp_90_sc:voltinismmultivoltine`
  )


#### Significance stars ####

stars <- function(p) {
  dplyr::case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE ~ ""
  )
}

#### Final clean table (publication-ready) ####

table_clean <- table_final %>%
  dplyr::transmute(
    model,
    delta = round(delta, 2),
    weight = round(weight, 2),
    
    Warming = paste0(
      round(Warming, 3),
      stars(p.value_clim_trend_temp_90_sc)
    ),
    
    Background = paste0(
      round(Background, 3),
      stars(p.value_clim_background_temp_90_sc)
    ),
    
    `Warming×Background` = paste0(
      round(`Warming×Background`, 3),
      stars(`p.value_clim_trend_temp_90_sc:clim_background_temp_90_sc`)
    ),
    
    Voltinism = paste0(
      round(Voltinism, 3),
      stars(p.value_voltinismmultivoltine)
    ),
    
    `Warming×Voltinism` = paste0(
      round(`Warming×Voltinism`, 3),
      stars(`p.value_clim_trend_temp_90_sc:voltinismmultivoltine`)
    )
  ) %>%
  dplyr::arrange(delta)

table_clean
  
#### Plots ####

# Onset 
mod_onset <- best_models[["ONSET_mean"]]$models[["m2"]]

pred <- ggpredict(mod_onset, terms = "clim_background_sc [-2:2]")

ggplot(pred, aes(x, predicted)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(
    x = "Background climate",
    y = "Onset trend"
  ) +
  theme_classic(base_family = "Garamond") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

pred <- ggpredict(mod_onset, terms = "clim_trend_sc [-2:2]")

ggplot(pred, aes(x, predicted)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(
    x = "Warming trend",
    y = "Onset trend"
  ) +
  theme_classic(base_family = "Garamond") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

# Offset (mean)
mod_offset <- best_models[["OFFSET_mean"]]$models[["m4"]]

pred <- ggpredict(mod_offset, terms = "clim_background_sc [-2:2]")

ggplot(pred, aes(x, predicted)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(
    x = "Background climate",
    y = "Offset trend"
  ) +
  theme_classic(base_family = "Garamond") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

# Offset (var)
mod_offset_var <- best_models[["OFFSET_var"]]$models[["m3"]]

pred <- ggpredict(
  mod_offset_var,
  terms = c("clim_trend_sc [-2:2]", "voltinism")
)

ggplot(pred, aes(x, predicted, color = group)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group),
              alpha = 0.2, color = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(
    x = "Warming",
    y = "Offset trend",
    color = "Voltinism",
    fill = "Voltinism"
  ) +
  theme_classic(base_family = "Garamond") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

# Flight length 

mod_fl <- best_models[["FLIGHT_LENGTH_mean"]]$models[["m2"]]

pred <- ggpredict(mod, terms = "clim_background_sc [-2:2]")

ggplot(pred, aes(x, predicted)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(
    x = "Background climate",
    y = "Flight length trend"
  ) +
  theme_classic(base_family = "Garamond") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )


################### Alternative model selection #######################################

#### Model selection ####

df$YEAR_c <- scale(df$YEAR, scale = FALSE)

model_set_year <- function(d) {
  
  d <- d %>%
    dplyr::filter(
      !is.na(ONSET_mean),
      !is.na(YEAR_c),
      !is.na(clim_background_temp_90_sc),
      !is.na(clim_trend_temp_90_sc),
      !is.na(photo_90_sc)
    )
  
  ctrl <- lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 1e6)
  )
  
  f_base <- ONSET_mean ~ YEAR_c +
    (1 | SPECIES) +
    (1 | SITE_ID)
  
  list(
    
    # H1
    m1 = lmer(f_base, data = d, REML = FALSE, control = ctrl),
    
    # H2
    m2 = lmer(update(f_base,
                     . ~ . +
                       YEAR_c:clim_trend_temp_90_sc),
              data = d, REML = FALSE, control = ctrl),
    
    # H3
    m3 = lmer(update(f_base,
                     . ~ . +
                       YEAR_c:clim_background_temp_90_sc),
              data = d, REML = FALSE, control = ctrl),
    
    # H4
    m4 = lmer(update(f_base,
                     . ~ . +
                       YEAR_c:photo_90_sc),
              data = d, REML = FALSE, control = ctrl),
    
    # H5
    m5 = lmer(update(f_base,
                     . ~ . + 
                       YEAR_c:clim_trend_temp_90_sc +
                       YEAR_c:clim_background_temp_90_sc),
              data = d, REML = FALSE, control = ctrl),
    
    # H6
    m6 = lmer(update(f_base,
                     . ~ . +
                       YEAR_c:clim_trend_temp_90_sc +
                       YEAR_c:photo_90_sc),
              data = d, REML = FALSE, control = ctrl),
    # H7
    
    m7 = lmer(update(f_base,
                     . ~ . +
                       YEAR_c:clim_background_temp_90_sc +
                       YEAR_c:photo_90_sc),
              data = d, REML = FALSE, control = ctrl),
    
    # H8
    m8 = lmer(update(f_base,
                     . ~ . +
                       YEAR_c:clim_background_temp_90_sc +
                       YEAR_c:clim_trend_temp_90_sc + 
                       YEAR_c:photo_90_sc),
              data = d, REML = FALSE, control = ctrl),
    
    m9 = lmer(update(f_base,
                     . ~ . +
                       YEAR_c:clim_trend_temp_90_sc:clim_background_temp_90_sc +
                       YEAR_c:photo_90_sc),
              data = d, REML = FALSE, control = ctrl),
    
    m10 = lmer(update(f_base,
                      . ~ . +
                        YEAR_c:clim_trend_temp_90_sc:photo_90_sc +
                        YEAR_c:clim_background_temp_90_sc),
               data = d, REML = FALSE, control = ctrl),
    
    m11 = lmer(update(f_base,
                      . ~ . + YEAR_c +
                        YEAR_c:clim_trend_temp_90_sc +
                        YEAR_c:clim_background_temp_90_sc +
                        YEAR_c:photo_90_sc +
                        YEAR_c:clim_trend_temp_90_sc:clim_background_temp_90_sc +
                        YEAR_c:clim_trend_temp_90_sc:photo_90_sc),
               data = d, REML = FALSE, control = ctrl)
    
    
    
  )
}

mods <- model_set_year(df)

sel <- MuMIn::model.sel(mods)

sel_df <- as.data.frame(sel) %>%
  tibble::rownames_to_column("model") %>%
  dplyr::filter(delta <= 2)

sel_df

sel_full <- as.data.frame(sel) %>%
  tibble::rownames_to_column("model") %>%
  dplyr::arrange(delta)

# Check collinearity
performance::check_collinearity(mods$m11)


table_final <- purrr::map_dfr(seq_len(nrow(sel_full)), function(i) {
  
  mod_name <- sel_full$model[i]
  mod      <- mods[[mod_name]]
  
  coefs <- broom.mixed::tidy(mod, effects = "fixed")
  
  if (!"p.value" %in% names(coefs)) {
    coefs$p.value <- NA
  }
  
  coefs %>%
    dplyr::select(term, estimate, std.error, p.value) %>%
    dplyr::mutate(
      model = mod_name,
      delta = sel_full$delta[i],
      weight = sel_full$weight[i]
    ) %>%
    tidyr::pivot_wider(
      names_from = term,
      values_from = c(estimate, std.error, p.value)
    )
})

table_final <- table_final %>%
  dplyr::rename(
    year = estimate_YEAR_c,
    Background = estimate_clim_background_temp_90_sc,
    Predictability = estimate_clim_trend_temp_90_sc,
    Photoperiod = estimate_photo_90_sc,
    
    `year×Background` =
      `estimate_YEAR_c:clim_background_temp_90_sc`,
    
    `year×Predictability` =
      `estimate_YEAR_c:clim_trend_temp_90_sc`,
    
    `year×Photoperiod` =
      `estimate_YEAR_c:photo_90_sc`
  )


table_final %>%
  dplyr::select(model, delta, weight,
                year, Background, Predictability, Photoperiod,
                `year×Background`,
                `year×Predictability`,
                `year×Photoperiod`)



stars <- function(p) {
  dplyr::case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE ~ ""
  )
}

table_clean <- table_final %>%
  dplyr::transmute(
    model,
    delta = round(delta, 1),
    
    weight = dplyr::case_when(
      weight == 0 ~ "<1e-50",
      weight < 1e-3 ~ formatC(weight, format = "e", digits = 2),
      TRUE ~ formatC(weight, format = "f", digits = 3)
    ),
    
    year = paste0(
      round(year, 2),
      stars(p.value_YEAR_c)
    ),
    
    Background = paste0(
      round(Background, 2),
      stars(p.value_clim_background_temp_90_sc)
    ),
    
    Predictability = paste0(
      round(Predictability, 2),
      stars(p.value_clim_trend_temp_90_sc)
    ),
    
    Photoperiod = paste0(
      round(Photoperiod, 2),
      stars(p.value_photo_90_sc)
    ),
    
    `year×Background` = paste0(
      round(`year×Background`, 2),
      stars(`p.value_YEAR_c:clim_background_temp_90_sc`)
    ),
    
    `year×Predictability` = paste0(
      round(`year×Predictability`, 2),
      stars(`p.value_YEAR_c:clim_trend_temp_90_sc`)
    ),
    
    `year×Photoperiod` = paste0(
      round(`year×Photoperiod`, 2),
      stars(`p.value_YEAR_c:photo_90_sc`)
    )
  ) %>%
  dplyr::arrange(delta) %>%
  dplyr::mutate(
    across(-c(model, delta, weight),
           ~ ifelse(. %in% c("NA", NA), "—", .))
  )

table_clean



#### Plot effects ####

#1. Forest plot function

make_forest <- function(model){
  
  plot_df <- broom.mixed::tidy(model, effects = "fixed") %>%
    dplyr::filter(term != "(Intercept)") %>%
    dplyr::mutate(
      lower = estimate - 1.96 * std.error,
      upper = estimate + 1.96 * std.error,
      
      term = dplyr::recode(
        term,
        clim_anomaly_temp_90_sc = "Plasticity",
        clim_background_temp_90_sc = "Background",
        clim_trend_temp_90_sc = "Predictability",
        photo_90_sc = "Photoperiod",
        `clim_anomaly_temp_90_sc:clim_background_temp_90_sc` = "Plasticity × Background",
        `clim_anomaly_temp_90_sc:clim_trend_temp_90_sc` = "Plasticity × Predictability",
        `clim_anomaly_temp_90_sc:photo_90_sc` = "Plasticity × Photoperiod"
      ),
      
      term = factor(term, levels = rev(c(
        "Plasticity",
        "Background",
        "Predictability",
        "Photoperiod",
        "Plasticity × Background",
        "Plasticity × Predictability",
        "Plasticity × Photoperiod"
      )))
    )
  
  p <- ggplot(plot_df,
              aes(x = term,
                  y = estimate,
                  ymin = lower,
                  ymax = upper)) +
    geom_pointrange(size = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    coord_flip() +
    theme_classic() +
    labs(
      x = "",
      y = "Effect size (estimate ± 95% CI)"
    )
  
  return(p)
}

p <- make_forest(mods$m8)
p

ggsave(
  filename = here::here("output", "figures", "forest_plasticity_m8B.png"),
  plot = p,
  width = 6,
  height = 4,
  dpi = 300
)


#2. Plot prediction effects.

final_mod <- mods$m11

# Year x background

pred <- ggpredict(
  final_mod,
  terms = c("YEAR_c", "clim_background_temp_90_sc[-1:1]")
)

pred$group <- factor(pred$group, levels = c("1", "0", "-1"))

pred$group <- factor(
  pred$group,
  levels = c("1", "0", "-1"),
  labels = c("Warm", "Mean", "Cold")
)

ggplot(pred, aes(x, predicted, color = group)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group),
              alpha = 0.2, color = NA) +
  
  scale_color_manual(values = c("Warm" = "#F28E2B", "Mean" = "#2CA58D", "Cold" = "#6C6BD1")) +
  scale_fill_manual(values = c("Warm" = "#F28E2B", "Mean" = "#2CA58D", "Cold" = "#6C6BD1")) +
  
  labs(
    x = "Year (centered)",
    y = "Onset (day of the year)",
    color = "Background climate",
    fill = "Background climate"
  ) +
  theme_classic(base_family = "Garamond") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8))


# Year x background

pred <- ggpredict(
  final_mod,
  terms = c("YEAR_c", "clim_trend_temp_90_sc[-1:1]")
)

pred$group <- factor(pred$group, levels = c("1", "0", "-1"))

pred$group <- factor(
  pred$group,
  levels = c("1", "0", "-1"),
  labels = c("Strong", "Mean", "Weak")
)

ggplot(pred, aes(x, predicted, color = group)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group),
              alpha = 0.2, color = NA) +
  
  scale_color_manual(values = c("Strong" = "#F28E2B", "Mean" = "#2CA58D", "Weak" = "#6C6BD1")) +
  scale_fill_manual(values = c("Strong" = "#F28E2B", "Mean" = "#2CA58D", "Weak" = "#6C6BD1")) +
  
  labs(
    x = "Year (centered)",
    y = "Onset (day of the year)",
    color = "Warming trend",
    fill = "Warming trend"
  ) +
  theme_classic(base_family = "Garamond") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8))




# Year x photoperiod


pred <- ggpredict(
  final_mod,
  terms = c("YEAR_c", "photo_90_sc[-1:1]")
)

pred$group <- factor(pred$group, levels = c("1", "0", "-1"))

pred$group <- factor(
  pred$group,
  levels = c("1", "0", "-1"),
  labels = c("High", "Mean", "Low")
)

ggplot(pred, aes(x, predicted, color = group)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group),
              alpha = 0.2, color = NA) +
  
  scale_color_manual(values = c("High" = "#F28E2B", "Mean" = "#2CA58D", "Low" = "#6C6BD1")) +
  scale_fill_manual(values = c("High" = "#F28E2B", "Mean" = "#2CA58D", "Low" = "#6C6BD1")) +
  
  labs(
    x = "Year (centered)",
    y = "Onset (day of the year)",
    color = "Photoperiod",
    fill = "Photoperiod"
  ) +
  theme_classic(base_family = "Garamond") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8))


# Triple interacció

pred <- ggpredict(
  mods$m11,
  terms = c(
    "YEAR_c",
    "clim_trend_temp_90_sc[1,0,-1]",
    "clim_background_temp_90_sc[-2,0,2]"
  )
)

# arreglar labels
pred$group <- factor(pred$group, levels = c("1", "0", "-1"),
                     labels = c("Stong warming", "Intermediate", "Weak warming"))

pred$facet <- factor(pred$facet, levels = c("-2", "0", "2"),
                     labels = c("Cold climate", "Intermediate", "Warm climate"))

ggplot(pred, aes(x = x, y = predicted, color = group)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group),
              alpha = 0.2, color = NA) +
  
  facet_wrap(~ facet) +
  
  scale_color_manual(values = c(
    "Stong warming" = "#F28E2B",
    "Intermediate" = "#2CA58D",
    "Weak warming" = "#6C6BD1"
  )) +
  scale_fill_manual(values = c(
    "Stong warming" = "#F28E2B",
    "Intermediate" = "#2CA58D",
    "Weak warming" = "#6C6BD1"
  )) +
  
  labs(
    x = "Year (centered)",
    y = "Onset (day of the year)",
    color = "Warming trend",
    fill = "Warming trend"
  ) +
  
  theme_classic(base_family = "Garamond") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

