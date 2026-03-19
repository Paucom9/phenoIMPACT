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

df <- pheno_trends_site |>
  
  # climate
  left_join(
    clim_vars |>
      distinct(
        transect_id,
        clim_background,
        clim_var,
        clim_pred_lag,
        clim_trend,
        clim_background_sc,
        clim_var_sc,
        clim_pred_lag_sc,
        clim_trend_sc
      ),
    by = c("SITE_ID" = "transect_id")
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
    dplyr::mutate(
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
      !is.na(clim_background_sc),
      !is.na(clim_trend_sc),
      !is.na(voltinism)
    )
  
  list(
    
    # H1: warming only
    m1 = lmer(
      estimate ~ clim_trend_sc +
        (1 + clim_trend_sc || SPECIES) +  # random slope for warming
        (1 | SITE_ID),                    # spatial grouping
      data = d,
      weights = 1 / (std.error^2),        # meta-analytic weighting
      REML = FALSE,
      control = lme4::lmerControl(optimizer = "bobyqa")
    ),
    
    # H2: warming × background climate
    m2 = lmer(
      estimate ~ clim_trend_sc * clim_background_sc +
        (1 + clim_trend_sc || SPECIES) +
        (1 | SITE_ID),
      data = d,
      weights = 1 / (std.error^2),
      REML = FALSE,
      control = lme4::lmerControl(optimizer = "bobyqa")
    ),
    
    # H3: warming × voltinism
    m3 = lmer(
      estimate ~ clim_trend_sc * voltinism +
        (1 + clim_trend_sc || SPECIES) +
        (1 | SITE_ID),
      data = d,
      weights = 1 / (std.error^2),
      REML = FALSE,
      control = lme4::lmerControl(optimizer = "bobyqa")
    ),
    
    # H4: full model (both moderators)
    m4 = lmer(
      estimate ~ clim_trend_sc * clim_background_sc +
        clim_trend_sc * voltinism +
        (1 + clim_trend_sc || SPECIES) +
        (1 | SITE_ID),
      data = d,
      weights = 1 / (std.error^2),
      REML = FALSE,
      control = lme4::lmerControl(optimizer = "bobyqa")
    )
  )
}


#### Model selection (ΔAIC ≤ 2 retained) ####

best_models <- df %>%
  split(.$phenovar) %>%
  purrr::map(function(d) {
    
    mods <- model_set(d)
    sel  <- MuMIn::model.sel(mods)
    
    sel_df <- as.data.frame(sel) %>%
      tibble::rownames_to_column("model") %>%
      dplyr::filter(delta <= 2)   # retain competing models
    
    list(
      models = mods,
      selection = sel_df
    )
  })


#### Extract coefficients from selected models ####

table_final <- purrr::imap_dfr(best_models, function(x, phen) {
  
  purrr::map_dfr(seq_len(nrow(x$selection)), function(i) {
    
    mod_name <- x$selection$model[i]
    mod      <- x$models[[mod_name]]
    
    coefs <- broom.mixed::tidy(mod, effects = "fixed")
    
    coefs %>%
      dplyr::select(term, estimate, std.error, p.value) %>%
      dplyr::mutate(
        phenovar = phen,
        model = mod_name,
        delta = x$selection$delta[i],
        weight = x$selection$weight[i]
      ) %>%
      tidyr::pivot_wider(
        names_from = term,
        values_from = c(estimate, std.error, p.value)
      )
  })
})


#### Rename predictors (clean table labels) ####

table_final <- table_final %>%
  dplyr::rename(
    Warming = estimate_clim_trend_sc,
    Background = estimate_clim_background_sc,
    `Warming×Background` = `estimate_clim_trend_sc:clim_background_sc`,
    Voltinism = estimate_voltinismmultivoltine,
    `Warming×Voltinism` = `estimate_clim_trend_sc:voltinismmultivoltine`
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
    phenovar,
    model,
    delta = round(delta, 2),
    weight = round(weight, 2),
    
    Warming = paste0(round(Warming, 3),
                     stars(p.value_clim_trend_sc)),
    
    Background = paste0(round(Background, 3),
                        stars(p.value_clim_background_sc)),
    
    `Warming×Background` = paste0(
      round(`Warming×Background`, 3),
      stars(`p.value_clim_trend_sc:clim_background_sc`)
    ),
    
    Voltinism = paste0(
      round(Voltinism, 3),
      stars(p.value_voltinismmultivoltine)
    ),
    
    `Warming×Voltinism` = paste0(
      round(`Warming×Voltinism`, 3),
      stars(`p.value_clim_trend_sc:voltinismmultivoltine`)
    )
  ) %>%
  dplyr::arrange(phenovar, delta)  # order by support

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

