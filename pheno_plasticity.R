# ============================================================================================ #
# pheno_plasticity.R
#
# Author: Pau Colom
# Date: 2026-02-20
#
# Desciription: This script analyzes phenological plasticity in response to temperature anomalies, 
# and how this plasticity is influenced by photoperiod and other climatic variables. 
#
# ============================================================================================ #

# Clean session
rm(list = ls())

#### Load required libraries ####
# ---
library(data.table)  # For efficient data handling
library(lme4)
library(dplyr)       # For data manipulation
library(tidyr)
library(purrr)
library(broom.mixed)
library(performance)
library(ggplot2)
library(sf)
library(stringr)
library(ggeffects)
library(Matrix)
library(lmerTest)
library(extrafont)
library(writexl)
library(MuMIn)
library(here)


# ---

#### Data Import and Preparation ####
# ---
here::here() # Check the current working directory
# ---
phenology_estimates  <- read.csv(here::here("output", "pheno_estimates_allspp.csv"), sep = ",", dec = ".")
clim_vars <- read.csv(here::here("output", "climate", "climate_variables_all_phenophases.csv"), sep = ",", dec = ".")
ebms_transect_coord <- read.csv(here::here("data", "ebms_transect_coord.csv"), sep = ",", dec = ".")
voltinism <- read.csv(here::here("data", "voltinism", "species_country_voltinism.csv"), sep = ";", dec = ".")


# ---

# Merge datasets #
# Site coordinates + BMS identity
coords_sf <- ebms_transect_coord |>
  dplyr::filter(!is.na(transect_lon),
                !is.na(transect_lat)) |>
  dplyr::distinct(bms_id, transect_id, transect_lon, transect_lat) |>
  sf::st_as_sf(
    coords = c("transect_lon", "transect_lat"),
    crs = 3035
  ) |>
  sf::st_transform(4326)

coord_site <- coords_sf |>
  dplyr::mutate(latitude = sf::st_coordinates(coords_sf)[, 2]) |>
  sf::st_drop_geometry() |>
  dplyr::select(
    SITE_ID = transect_id,
    bms_id,
    latitude
  )

# Add latitude and BMS identity to climate variables
clim_vars <- clim_vars |>
  dplyr::select(-latitude, -dplyr::any_of("bms_id")) |>
  dplyr::left_join(coord_site, by = "SITE_ID") |>
  dplyr::mutate(latitude = as.numeric(scale(latitude)))

# Voltinism in long format
voltinism_long <- voltinism |>
  tidyr::pivot_longer(
    cols = -SPECIES,
    names_to = "bms_id",
    values_to = "voltinism"
  ) |>
  dplyr::mutate(
    bms_id = dplyr::recode(
      bms_id,
      "ES.CTBMS" = "ES-CTBMS",
      "ES.ZEBMS" = "ES-ZEBMS"
    )
  )

# Merge phenology + climate + voltinism
df <- phenology_estimates |>
  dplyr::left_join(
    clim_vars,
    by = c("SPECIES", "SITE_ID", "YEAR")
  ) |>
  dplyr::left_join(
    voltinism_long,
    by = c("SPECIES", "bms_id")
  ) |>
  dplyr::mutate(
    SPECIES = factor(SPECIES),
    SITE_ID = factor(SITE_ID),
    bms_id = factor(bms_id),
    voltinism = factor(voltinism),
    sp_site = interaction(SPECIES, SITE_ID, drop = TRUE)
  )


# ------------------------------------------------------------------------------------- #
#### Data exploration #### 

# Inspect climatic anomalies # 

df_unique <- df %>%
  group_by(SITE_ID, SPECIES, YEAR) %>%
  slice(1) %>%
  ungroup()

med <- median(df_unique$clim_anomaly_tw30, na.rm = TRUE)

p <- ggplot(df_unique, aes(x = clim_anomaly_tw30)) +
  geom_histogram(fill = "grey", color = "white", bins = 30) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1) +
  geom_vline(xintercept = med, color = "red", linetype = "dashed", size = 1) +
  annotate("text",
           x = med,
           y = Inf,
           label = paste("Median =", round(med, 2)),
           vjust = 2,
           hjust = 0,
           color = "red") +
  labs(
    x = "Temperature anomaly (°C)",
    y = "Count"
  ) +
  theme_minimal()

ggsave(
  filename = here::here("output", "figures", "distribution_clim_anomalies_27_4_26.png"),
  plot = p,
  width = 8,
  height = 6
)

# Inspect distribution of predictor variables #

# variables to plot (explicit → evita errors)
vars_plot <- c(
  "clim_background_tw90",
  "clim_trend_tw90",
  "clim_autocorr_tw90",
  "clim_stability_tw90",
  "clim_predictability_tw90",
  "photo_tw90",
  "latitude"
)

# remove pseudo-replication (one row per site × year × pheno_type)
df_plot <- df |>
  dplyr::distinct(
    SITE_ID,
    YEAR,
    pheno_type,
    dplyr::across(dplyr::all_of(vars_plot))
  )

# long format
df_long <- df_plot |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(vars_plot),
    names_to = "variable",
    values_to = "value"
  )

# labels
labels <- c(
  clim_background_tw90      = "Background temperature\n(°C)",
  clim_trend_tw90           = "Temperature trend\n(°C per decade)",
  clim_autocorr_tw90        = "Autocorrelation\n(lag-1)",
  clim_stability_tw90       = "Stability\n(-SD)",
  clim_predictability_tw90  = "Predictability\n(-Var residuals)",
  photo_tw90                = "Photoperiod\n(daily hours)",
  latitude                  = "Latitude\n(degrees)"
)

# apply labels safely
df_long$variable <- factor(
  df_long$variable,
  levels = names(labels),
  labels = labels
)

# plot
p <- ggplot(df_long, aes(x = value)) +
  geom_histogram(bins = 30, fill = "grey70", color = "black") +
  facet_wrap(~ variable, scales = "free", ncol = 3) +
  theme_classic(base_size = 12) +
  labs(x = NULL, y = "Frequency") +
  theme(
    strip.text = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
  )

# save
ggsave(
  filename = here::here("output", "figures", "predictor_histograms_26_3_26B.png"),
  plot = p,
  width = 9,
  height = 7,
  dpi = 300
)

# ---

# Correlation among predictors # 

cor_mat <- df_cor |>
  select(clim_background_temp_90,
         clim_trend_temp_90,
         clim_autocorr_temp_90,
         clim_stability_temp_90,
         clim_predictability_temp_90,
         photo_winter_fixed_fixed,
         latitude) |>
  cor(use = "complete.obs")


new_names <- c(
  "Background",
  "Trend",
  "Autocorr",
  "Stability",
  "Predictability",
  "Photoperiod",
  "Latitude"
)

rownames(cor_mat) <- new_names
colnames(cor_mat) <- new_names

cor_df <- as.data.frame(as.table(cor_mat))

# keep only lower triangle
cor_df <- cor_df |>
  filter(as.numeric(Var1) > as.numeric(Var2))

corr_clim_pred_sds <- ggplot(cor_df, aes(Var1, Var2, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = round(Freq, 2)), size = 4) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, limits = c(-1,1)) +
  theme_minimal() +
  labs(x = NULL, y = NULL, fill = "r") +
  coord_equal()

corr_clim_pred_sds

# Save plot
ggsave(
  filename = here::here("output", "figures", "corr_site_vars_26_3_26B.png"),
  plot = corr_clim_pred_sds,
  width = 6,
  height = 5,
  dpi = 300
)

# ---------------------------------------------------------------------------- #
#### Onset (mean) #### 

df_onset_mean <- subset(df, df$pheno_type == "ONSET_mean")

# Anomaly time-window selection # 

mod_30 <- lmer(
  ONSET_mean ~ clim_anomaly_tw30 +
    (1 | SITE_ID) + (1 + clim_anomaly_tw30 | SPECIES),
  data = df_onset_mean, REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e6)
  )
)

mod_60 <- lmer(
  ONSET_mean ~ clim_anomaly_tw60 +
    (1 | SITE_ID) + (1 + clim_anomaly_tw60 | SPECIES),
  data = df_onset_mean, REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e6)
  )
)

mod_90 <- lmer(
  ONSET_mean ~ clim_anomaly_tw90 +
    (1 | SITE_ID) + (1 + clim_anomaly_tw90 | SPECIES),
  data = df_onset_mean, REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e6)
  )
)

AIC(mod_30, mod_60, mod_90) ### Best time-window is 60 days before onset
summary(mod_30)
summary(mod_60)
summary(mod_90)
fixef(mod_30)["clim_anomaly_tw30"]
fixef(mod_60)["clim_anomaly_tw60"]
fixef(mod_90)["clim_anomaly_tw90"]


# Plot effects

extract_effect <- function(mod, var, label) {
  tidy(mod) %>%
    filter(term == var) %>%
    mutate(window = label)
}

df_eff <- bind_rows(
  extract_effect(mod_30, "clim_anomaly_tw30", "30 days"),
  extract_effect(mod_60, "clim_anomaly_tw60", "60 days"),
  extract_effect(mod_90, "clim_anomaly_tw90", "90 days")
) %>%
  mutate(
    lower = estimate - 1.96 * std.error,
    upper = estimate + 1.96 * std.error
  )

ggplot(df_eff, aes(x = window, y = estimate)) +
  
  geom_point(size = 3) +
  
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0.1) +
  
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  theme_classic(base_family = "Garamond") +
  
  labs(
    x = "Time window",
    y = "Phenological sensitivity\n(slope of onset vs temperature anomaly)"
  )

make_pred <- function(mod, var, label) {
  
  newdat <- data.frame(
    anomaly = seq(-2, 2, length.out = 100)
  )
  
  names(newdat) <- var
  
  pred <- predict(
    mod,
    newdata = newdat,
    re.form = NA,   # important!
    se.fit = TRUE
  )
  
  data.frame(
    anomaly = newdat[[var]],
    fit = pred$fit,
    se = pred$se.fit,
    window = label
  )
}

df_pred <- dplyr::bind_rows(
  make_pred(mod_30, "clim_anomaly_tw30", "30 days"),
  make_pred(mod_60, "clim_anomaly_tw60", "60 days"),
  make_pred(mod_90, "clim_anomaly_tw90", "90 days")
)

df_pred <- df_pred %>%
  dplyr::mutate(
    lower = fit - 1.96 * se,
    upper = fit + 1.96 * se
  )

ggplot(df_pred, aes(anomaly, fit,
                    color = window,
                    fill = window)) +
  
  geom_line(size = 1.2) +
  
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.2, color = NA) +
  
  theme_classic(base_family = "Garamond") +
  
  labs(
    x = "Temperature anomaly",
    y = "Predicted onset",
    color = "Time window",
    fill = "Time window"
  )


# ---------------------------------------------------------------------------- #
# Model selection #

d_o_m <- df_onset_mean %>%
  dplyr::filter(
    !is.na(ONSET_mean),
    !is.na(clim_anomaly_tw60),
    !is.na(photo_tw60),
    !is.na(clim_background_tw60),
    !is.na(clim_predictability_tw60),
    !is.na(clim_autocorr_tw60),
    !is.na(clim_trend_tw60)
  )

global_mod_o_m <- lmer(
  ONSET_mean ~ clim_anomaly_tw60 *
    (photo_tw60 +
       clim_background_tw60 +
       clim_predictability_tw60 +
       clim_autocorr_tw60 +
       clim_trend_tw60) +
    (1 | SITE_ID) +
    (1 + clim_anomaly_tw60 | SPECIES),
  data = d_o_m, REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e6)
  )
)

options(na.action = "na.fail")

dd_o_m <- dredge(
  global_mod_o_m,
  fixed = "clim_anomaly_tw60",
  subset =
    dc("photo_tw60", "clim_anomaly_tw60:photo_tw60") &
    dc("clim_background_tw60", "clim_anomaly_tw60:clim_background_tw60") &
    dc("clim_predictability_tw60", "clim_anomaly_tw60:clim_predictability_tw60") &
    dc("clim_autocorr_tw60", "clim_anomaly_tw60:clim_autocorr_tw60") &
    dc("clim_trend_tw60", "clim_anomaly_tw60:clim_trend_tw60")
)

model.sel(dd_o_m)

best_mod_o_m <- get.models(dd_o_m, subset = 1)[[1]]
summary(best_mod_o_m)
check_collinearity(best_mod_o_m)

# --- all models ---
ms_o_m <- model.sel(dd_o_m)

all_models_o_m <- as.data.frame(ms_o_m)

# --- stars function ---
stars <- function(p) {
  ifelse(p < 0.001, "***",
         ifelse(p < 0.01, "**",
                ifelse(p < 0.05, "*", "")))
}

# --- p-values from best model ---
pvals <- summary(best_mod_o_m)$coefficients[, "Pr(>|t|)"]

# --- clean table ---
fmt <- function(x, pvals, name) {
  
  p <- pvals[name]
  
  if (is.na(p) && grepl(":", name)) {
    parts <- strsplit(name, ":")[[1]]
    alt_name <- paste0(parts[2], ":", parts[1])
    if (alt_name %in% names(pvals)) {
      p <- pvals[alt_name]
    }
  }
  
  ifelse(
    is.na(x),
    "—",
    paste0(
      round(x, 2),
      ifelse(!is.na(p), stars(p), "")
    )
  )
}

table_clean_o_m <- all_models_o_m %>%
  transmute(
    delta = round(delta, 2),
    weight = signif(weight, 2),
    
    Plasticity = fmt(clim_anomaly_tw60, pvals, "clim_anomaly_tw60"),
    Photoperiod = fmt(photo_tw60, pvals, "photo_tw60"),
    Background = fmt(clim_background_tw60, pvals, "clim_background_tw60"),
    Predictability = fmt(clim_predictability_tw60, pvals, "clim_predictability_tw60"),
    Autocorrelation = fmt(clim_autocorr_tw60, pvals, "clim_autocorr_tw60"),
    Trend = fmt(clim_trend_tw60, pvals, "clim_trend_tw60"),
    
    `Plast×Photo` = fmt(`clim_anomaly_tw60:photo_tw60`, pvals, "clim_anomaly_tw60:photo_tw60"),
    `Plast×Background` = fmt(`clim_anomaly_tw60:clim_background_tw60`, pvals, "clim_anomaly_tw60:clim_background_tw60"),
    `Plast×Predictability` = fmt(`clim_anomaly_tw60:clim_predictability_tw60`, pvals, "clim_anomaly_tw60:clim_predictability_tw60"),
    `Plast×Autocorr` = fmt(`clim_anomaly_tw60:clim_autocorr_tw60`, pvals, "clim_anomaly_tw60:clim_autocorr_tw60"),
    `Plast×Trend` = fmt(`clim_anomaly_tw60:clim_trend_tw60`, pvals, "clim_anomaly_tw60:clim_trend_tw60")
  )

table_clean_o_m


write_xlsx(
  table_clean_o_m,
  here::here("output", "onset_mean_plasticity_table_result.xlsx")
)

# ---------------------------------------------------------------------------- #
# Plot interactions of the best model #

mod_final_o_m <- lmer(
  ONSET_mean ~ clim_anomaly_tw60 *
    (photo_tw60 +
       clim_background_tw60 +
       clim_predictability_tw60 + clim_autocorr_tw60) +
    (1 | SITE_ID) +
    (1 + clim_anomaly_tw60 | SPECIES),
  
  data = d_o_m, REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e6)
  )
)

summary(mod_final_o_m)
performance::check_collinearity(mod_final_o_m)

# sequences for predictor variables

photo_seq <- seq(-2, 2, length.out = 100)
trend_seq    <- seq(-2, 2, length.out = 100)
pred_seq  <- seq(-2, 2, length.out = 100)
auto_seq  <- seq(-2, 2, length.out = 100)
back_seq  <- seq(-2, 2, length.out = 100)


# get slopes of plasticity across gradients

get_slope_df <- function(mod, var_name, var_seq, label,
                         anomaly_name = "clim_anomaly_tw60") {
  
  b <- fixef(mod)
  V <- vcov(mod)
  
  coef_main <- anomaly_name
  
  coef_int <- paste0(coef_main, ":", var_name)
  if (!coef_int %in% names(b)) {
    coef_int <- paste0(var_name, ":", coef_main)
  }
  
  if (!all(c(coef_main, coef_int) %in% names(b))) {
    stop(paste("No s'han trobat coeficients per:", var_name))
  }
  
  coef_names <- c(coef_main, coef_int)
  
  get_slope <- function(x) {
    
    g <- c(1, x)
    
    slope <- sum(b[coef_names] * g)
    
    se <- as.numeric(sqrt(t(g) %*% V[coef_names, coef_names] %*% g))
    
    c(slope = slope, se = se)
  }
  
  out <- do.call(rbind, lapply(var_seq, get_slope))
  
  data.frame(
    x = var_seq,
    slope = out[, "slope"],
    lower = out[, "slope"] - 1.96 * out[, "se"],
    upper = out[, "slope"] + 1.96 * out[, "se"],
    variable = label
  )
}

df_photo <- get_slope_df(
  mod_final_o_m,
  "photo_tw60",
  photo_seq,
  "Photoperiod"
)


df_pred <- get_slope_df(
  mod_final_o_m,
  "clim_predictability_tw60",
  pred_seq,
  "Temperature stability"
)

df_auto <- get_slope_df(
  mod_final_o_m,
  "clim_autocorr_tw60",
  auto_seq,
  "Temperature autocorrelation"
)

df_back <- get_slope_df(
  mod_final_o_m,
  "clim_background_tw60",
  back_seq,
  "Mean temperature"
)

# combine data frames

df_all <- dplyr::bind_rows(
  df_photo,
  df_pred,
  df_auto,
  df_back
)

# order facets

df_all$variable <- factor(
  df_all$variable,
  levels = c(
    "Photoperiod",
    "Temperature trend",
    "Temperature stability",
    "Temperature autocorrelation",
    "Mean temperature"
  )
)

# plot

plasticity_interactions_plot <- ggplot(df_all, aes(x, slope, color = variable, fill = variable)) +
  
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.15, color = NA) +
  
  geom_line(size = 0.8, alpha = 0.6) +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  
  scale_color_manual(values = c(
    "Photoperiod" = "#1b9e77",
    "Temperature stability" = "#7570b3",
    "Temperature autocorrelation" = "#E6AB02",
    "Mean temperature" = "#d95f02"
    
  )) +
  
  scale_fill_manual(values = c(
    "Photoperiod" = "#1b9e77",
    "Temperature stability" = "#7570b3",
    "Temperature autocorrelation" = "#E6AB02",
    "Mean temperature" = "#d95f02"
  )) +
  
  coord_cartesian(ylim = c(-9, -1)) +
  
  theme_classic(base_family = "Garamond", base_size = 14) +
  
  labs(
    x = "Environmental gradient",
    y = "Phenological plasticity\n(slope of onset vs temperature anomaly)",
    color = "Environmental gradient",
    fill = "Environmental gradient"
  )
plasticity_interactions_plot


ggsave(
  filename = here("output", "figures", "onset_mean_plasticity_interactions_plot.png"),
  plot = plasticity_interactions_plot,
  width = 7,
  height = 5,
  dpi = 300
)


# ---------------------------------------------------------------------------- #
#### Onset (variance) #### 

df_onset_var <- subset(df, df$pheno_type == "ONSET_var")

# Anomaly time-window selection # 

options(na.action = "na.omit")

mod_30_o_v <- lmer(
  ONSET_var ~ clim_anomaly_tw30 +
    (1 | SITE_ID) + (1 + clim_anomaly_tw30 | SPECIES),
  data = df_onset_var, REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e6)
  )
)

mod_60_o_v <- lmer(
  ONSET_var ~ clim_anomaly_tw60 +
    (1 | SITE_ID) + (1 + clim_anomaly_tw60 | SPECIES),
  data = df_onset_var, REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e6)
  )
)

mod_90_o_v <- lmer(
  ONSET_var ~ clim_anomaly_tw90 +
    (1 | SITE_ID) + (1 + clim_anomaly_tw90 | SPECIES),
  data = df_onset_var, REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e6)
  )
)

AIC(mod_30_o_v, mod_60_o_v, mod_90_o_v) ### Best time-window is 60 days before onset
summary(mod_30_o_v)
summary(mod_60_o_v)
summary(mod_90_o_v)
fixef(mod_30_o_v)["clim_anomaly_tw30"]
fixef(mod_60_o_v)["clim_anomaly_tw60"]
fixef(mod_90_o_v)["clim_anomaly_tw90"]


# Plot effects

extract_effect <- function(mod, var, label) {
  tidy(mod) %>%
    filter(term == var) %>%
    mutate(window = label)
}

df_eff <- bind_rows(
  extract_effect(mod_30_o_v, "clim_anomaly_tw30", "30 days"),
  extract_effect(mod_60_o_v, "clim_anomaly_tw60", "60 days"),
  extract_effect(mod_90_o_v, "clim_anomaly_tw90", "90 days")
) %>%
  mutate(
    lower = estimate - 1.96 * std.error,
    upper = estimate + 1.96 * std.error
  )

ggplot(df_eff, aes(x = window, y = estimate)) +
  
  geom_point(size = 3) +
  
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0.1) +
  
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  theme_classic(base_family = "Garamond") +
  
  labs(
    x = "Time window",
    y = "Phenological sensitivity\n(slope of onset vs temperature anomaly)"
  )

make_pred <- function(mod, var, label) {
  
  newdat <- data.frame(
    anomaly = seq(-2, 2, length.out = 100)
  )
  
  names(newdat) <- var
  
  pred <- predict(
    mod,
    newdata = newdat,
    re.form = NA,   # important!
    se.fit = TRUE
  )
  
  data.frame(
    anomaly = newdat[[var]],
    fit = pred$fit,
    se = pred$se.fit,
    window = label
  )
}

df_pred <- dplyr::bind_rows(
  make_pred(mod_30_o_v, "clim_anomaly_tw30", "30 days"),
  make_pred(mod_60_o_v, "clim_anomaly_tw60", "60 days"),
  make_pred(mod_90_o_v, "clim_anomaly_tw90", "90 days")
)

df_pred <- df_pred %>%
  dplyr::mutate(
    lower = fit - 1.96 * se,
    upper = fit + 1.96 * se
  )

ggplot(df_pred, aes(anomaly, fit,
                    color = window,
                    fill = window)) +
  
  geom_line(size = 1.2) +
  
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.2, color = NA) +
  
  theme_classic(base_family = "Garamond") +
  
  labs(
    x = "Temperature anomaly",
    y = "Predicted onset",
    color = "Time window",
    fill = "Time window"
  )


# ---------------------------------------------------------------------------- #
# Model selection #

d_o_v <- df_onset_var %>%
  dplyr::filter(
    !is.na(ONSET_var),
    !is.na(clim_anomaly_tw60),
    !is.na(photo_tw60),
    !is.na(clim_background_tw60),
    !is.na(clim_predictability_tw60),
    !is.na(clim_autocorr_tw60),
    !is.na(clim_trend_tw60)
  )

global_mod_o_v <- lmer(
  ONSET_var ~ clim_anomaly_tw60 *
    (photo_tw60 +
       clim_background_tw60 +
       clim_predictability_tw60 +
       clim_autocorr_tw60 +
       clim_trend_tw60) +
    (1 | SITE_ID) +
    (1 + clim_anomaly_tw60 | SPECIES),
  data = d_o_v, REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e6)
  )
)

options(na.action = "na.fail")

dd_o_v <- dredge(
  global_mod_o_v,
  fixed = "clim_anomaly_tw60",
  subset =
    dc("photo_tw60", "clim_anomaly_tw60:photo_tw60") &
    dc("clim_background_tw60", "clim_anomaly_tw60:clim_background_tw60") &
    dc("clim_predictability_tw60", "clim_anomaly_tw60:clim_predictability_tw60") &
    dc("clim_autocorr_tw60", "clim_anomaly_tw60:clim_autocorr_tw60") &
    dc("clim_trend_tw60", "clim_anomaly_tw60:clim_trend_tw60")
)

model.sel(dd_o_v)

best_mod_o_v <- get.models(dd_o_v, subset = 1)[[1]]
summary(best_mod_o_v)
check_collinearity(best_mod_o_v)

# --- all models ---
ms_o_v <- model.sel(dd_o_v)

all_models_o_v <- as.data.frame(ms_o_v)

# --- stars function ---
stars <- function(p) {
  ifelse(p < 0.001, "***",
         ifelse(p < 0.01, "**",
                ifelse(p < 0.05, "*", "")))
}

# --- p-values from best model ---
pvals <- summary(best_mod_o_v)$coefficients[, "Pr(>|t|)"]

# --- clean table ---
fmt <- function(x, pvals, name) {
  
  p <- pvals[name]
  
  if (is.na(p) && grepl(":", name)) {
    parts <- strsplit(name, ":")[[1]]
    alt_name <- paste0(parts[2], ":", parts[1])
    if (alt_name %in% names(pvals)) {
      p <- pvals[alt_name]
    }
  }
  
  ifelse(
    is.na(x),
    "—",
    paste0(
      round(x, 2),
      ifelse(!is.na(p), stars(p), "")
    )
  )
}

table_clean_o_v <- all_models_o_v %>%
  transmute(
    delta = round(delta, 2),
    weight = signif(weight, 2),
    
    Plasticity = fmt(clim_anomaly_tw60, pvals, "clim_anomaly_tw60"),
    Photoperiod = fmt(photo_tw60, pvals, "photo_tw60"),
    Background = fmt(clim_background_tw60, pvals, "clim_background_tw60"),
    Predictability = fmt(clim_predictability_tw60, pvals, "clim_predictability_tw60"),
    Autocorrelation = fmt(clim_autocorr_tw60, pvals, "clim_autocorr_tw60"),
    Trend = fmt(clim_trend_tw60, pvals, "clim_trend_tw60"),
    
    `Plast×Photo` = fmt(`clim_anomaly_tw60:photo_tw60`, pvals, "clim_anomaly_tw60:photo_tw60"),
    `Plast×Background` = fmt(`clim_anomaly_tw60:clim_background_tw60`, pvals, "clim_anomaly_tw60:clim_background_tw60"),
    `Plast×Predictability` = fmt(`clim_anomaly_tw60:clim_predictability_tw60`, pvals, "clim_anomaly_tw60:clim_predictability_tw60"),
    `Plast×Autocorr` = fmt(`clim_anomaly_tw60:clim_autocorr_tw60`, pvals, "clim_anomaly_tw60:clim_autocorr_tw60"),
    `Plast×Trend` = fmt(`clim_anomaly_tw60:clim_trend_tw60`, pvals, "clim_anomaly_tw60:clim_trend_tw60")
  )

table_clean_o_v


write_xlsx(
  table_clean_o_v,
  here::here("output", "onset_var_plasticity_table_result.xlsx")
)

# ---------------------------------------------------------------------------- #
# Plot interactions of the best model #

mod_final_o_v <- lmer(
  ONSET_var ~ clim_anomaly_tw60 *
    (photo_tw60 +
       clim_background_tw60 +
       clim_predictability_tw60 + clim_autocorr_tw60 + clim_trend_tw60) +
    (1 | SITE_ID) +
    (1 + clim_anomaly_tw60 | SPECIES),
  
  data = d_o_v, REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e6)
  )
)

summary(mod_final_o_v)
performance::check_collinearity(mod_final_o_v)

# sequences for predictor variables

photo_seq <- seq(-2, 2, length.out = 100)
trend_seq    <- seq(-2, 2, length.out = 100)
pred_seq  <- seq(-2, 2, length.out = 100)
auto_seq  <- seq(-2, 2, length.out = 100)
back_seq  <- seq(-2, 2, length.out = 100)
df_trend <- seq(-2, 2, length.out = 100)

# get slopes of plasticity across gradients

get_slope_df <- function(mod, var_name, var_seq, label,
                         anomaly_name = "clim_anomaly_tw60") {
  
  b <- fixef(mod)
  V <- vcov(mod)
  
  coef_main <- anomaly_name
  
  coef_int <- paste0(coef_main, ":", var_name)
  if (!coef_int %in% names(b)) {
    coef_int <- paste0(var_name, ":", coef_main)
  }
  
  if (!all(c(coef_main, coef_int) %in% names(b))) {
    stop(paste("No s'han trobat coeficients per:", var_name))
  }
  
  coef_names <- c(coef_main, coef_int)
  
  get_slope <- function(x) {
    
    g <- c(1, x)
    
    slope <- sum(b[coef_names] * g)
    
    se <- as.numeric(sqrt(t(g) %*% V[coef_names, coef_names] %*% g))
    
    c(slope = slope, se = se)
  }
  
  out <- do.call(rbind, lapply(var_seq, get_slope))
  
  data.frame(
    x = var_seq,
    slope = out[, "slope"],
    lower = out[, "slope"] - 1.96 * out[, "se"],
    upper = out[, "slope"] + 1.96 * out[, "se"],
    variable = label
  )
}

df_photo <- get_slope_df(
  mod_final_o_v,
  "photo_tw60",
  photo_seq,
  "Photoperiod"
)


df_pred <- get_slope_df(
  mod_final_o_v,
  "clim_predictability_tw60",
  pred_seq,
  "Temperature stability"
)

df_auto <- get_slope_df(
  mod_final_o_v,
  "clim_autocorr_tw60",
  auto_seq,
  "Temperature autocorrelation"
)

df_back <- get_slope_df(
  mod_final_o_v,
  "clim_background_tw60",
  back_seq,
  "Mean temperature"
)

df_trend <- get_slope_df(
  mod_final_o_v,
  "clim_trend_tw60",
  trend_seq,
  "Temperature trend"
)

# combine data frames

df_all <- dplyr::bind_rows(
  df_photo,
  df_pred,
  df_auto,
  df_back,
  df_trend
)

# order facets

df_all$variable <- factor(
  df_all$variable,
  levels = c(
    "Photoperiod",
    "Temperature trend",
    "Temperature stability",
    "Temperature autocorrelation",
    "Mean temperature"
  )
)

# plot

plasticity_interactions_plot <- ggplot(df_all, aes(x, slope, color = variable, fill = variable)) +
  
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.15, color = NA) +
  
  geom_line(size = 0.8, alpha = 0.6) +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  
  scale_color_manual(values = c(
    "Photoperiod" = "#1b9e77",
    "Temperature stability" = "#7570b3",
    "Temperature autocorrelation" = "#E6AB02",
    "Mean temperature" = "#d95f02",
    "Temperature trend" = "#66a34e"
    
  )) +
  
  scale_fill_manual(values = c(
    "Photoperiod" = "#1b9e77",
    "Temperature stability" = "#7570b3",
    "Temperature autocorrelation" = "#E6AB02",
    "Mean temperature" = "#d95f02",
    "Temperature trend" = "#66a34e"
  )) +
  
  coord_cartesian(ylim = c(-9.5, -2)) +
  
  theme_classic(base_family = "Garamond", base_size = 14) +
  
  labs(
    x = "Environmental gradient",
    y = "Phenological plasticity\n(slope of onset vs temperature anomaly)",
    color = "Environmental gradient",
    fill = "Environmental gradient"
  )
plasticity_interactions_plot


ggsave(
  filename = here("output", "figures", "onset_var_plasticity_interactions_plot.png"),
  plot = plasticity_interactions_plot,
  width = 7,
  height = 5,
  dpi = 300
)


# ---------------------------------------------------------------------------- #
#### Peak (first peak) #### 

df_first_peak <- subset(df, df$pheno_type == "FIRST_PEAK")

# Anomaly time-window selection # 

options(na.action = "na.omit")

mod_30_f_p <- lmer(
  FIRST_PEAK ~ clim_anomaly_tw30 +
    (1 | SITE_ID) + (1 + clim_anomaly_tw30 | SPECIES),
  data = df_first_peak, REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e6)
  )
)

mod_60_f_p <- lmer(
  FIRST_PEAK ~ clim_anomaly_tw60 +
    (1 | SITE_ID) + (1 + clim_anomaly_tw60 | SPECIES),
  data = df_first_peak, REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e6)
  )
)

mod_90_f_p <- lmer(
  FIRST_PEAK ~ clim_anomaly_tw90 +
    (1 | SITE_ID) + (1 + clim_anomaly_tw90 | SPECIES),
  data = df_first_peak, REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e6)
  )
)

AIC(mod_30_f_p, mod_60_f_p, mod_90_f_p) ### Best time-window is 60 days before onset
summary(mod_30_f_p)
summary(mod_60_f_p)
summary(mod_90_f_p)
fixef(mod_30_f_p)["clim_anomaly_tw30"]
fixef(mod_60_f_p)["clim_anomaly_tw60"]
fixef(mod_90_f_p)["clim_anomaly_tw90"]


# Plot effects

extract_effect <- function(mod, var, label) {
  tidy(mod) %>%
    filter(term == var) %>%
    mutate(window = label)
}

df_eff <- bind_rows(
  extract_effect(mod_30_f_p, "clim_anomaly_tw30", "30 days"),
  extract_effect(mod_60_f_p, "clim_anomaly_tw60", "60 days"),
  extract_effect(mod_90_f_p, "clim_anomaly_tw90", "90 days")
) %>%
  mutate(
    lower = estimate - 1.96 * std.error,
    upper = estimate + 1.96 * std.error
  )

ggplot(df_eff, aes(x = window, y = estimate)) +
  
  geom_point(size = 3) +
  
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0.1) +
  
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  theme_classic(base_family = "Garamond") +
  
  labs(
    x = "Time window",
    y = "Phenological sensitivity\n(slope of onset vs temperature anomaly)"
  )

make_pred <- function(mod, var, label) {
  
  newdat <- data.frame(
    anomaly = seq(-2, 2, length.out = 100)
  )
  
  names(newdat) <- var
  
  pred <- predict(
    mod,
    newdata = newdat,
    re.form = NA,   # important!
    se.fit = TRUE
  )
  
  data.frame(
    anomaly = newdat[[var]],
    fit = pred$fit,
    se = pred$se.fit,
    window = label
  )
}

df_pred <- dplyr::bind_rows(
  make_pred(mod_30_f_p, "clim_anomaly_tw30", "30 days"),
  make_pred(mod_60_f_p, "clim_anomaly_tw60", "60 days"),
  make_pred(mod_90_f_p, "clim_anomaly_tw90", "90 days")
)

df_pred <- df_pred %>%
  dplyr::mutate(
    lower = fit - 1.96 * se,
    upper = fit + 1.96 * se
  )

ggplot(df_pred, aes(anomaly, fit,
                    color = window,
                    fill = window)) +
  
  geom_line(size = 1.2) +
  
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.2, color = NA) +
  
  theme_classic(base_family = "Garamond") +
  
  labs(
    x = "Temperature anomaly",
    y = "Predicted first peak day",
    color = "Time window",
    fill = "Time window"
  )


# ---------------------------------------------------------------------------- #
# Model selection #

d_f_p <- df_first_peak %>%
  dplyr::filter(
    !is.na(FIRST_PEAK),
    !is.na(clim_anomaly_tw90),
    !is.na(photo_tw90),
    !is.na(clim_background_tw90),
    !is.na(clim_predictability_tw90),
    !is.na(clim_autocorr_tw90),
    !is.na(clim_trend_tw90)
  )

global_mod_f_p <- lmer(
  FIRST_PEAK ~ clim_anomaly_tw90 *
    (photo_tw90 +
       clim_background_tw90 +
       clim_predictability_tw90 +
       clim_autocorr_tw90 +
       clim_trend_tw90) +
    (1 | SITE_ID) +
    (1 + clim_anomaly_tw90 | SPECIES),
  data = d_f_p, REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e6)
  )
)

options(na.action = "na.fail")

dd_f_p <- dredge(
  global_mod_f_p,
  fixed = "clim_anomaly_tw90",
  subset =
    dc("photo_tw90", "clim_anomaly_tw90:photo_tw90") &
    dc("clim_background_tw90", "clim_anomaly_tw90:clim_background_tw90") &
    dc("clim_predictability_tw90", "clim_anomaly_tw90:clim_predictability_tw90") &
    dc("clim_autocorr_tw90", "clim_anomaly_tw90:clim_autocorr_tw90") &
    dc("clim_trend_tw90", "clim_anomaly_tw90:clim_trend_tw90")
)

model.sel(dd_f_p)

best_mod_f_p <- get.models(dd_f_p, subset = 1)[[1]]
summary(best_mod_f_p)
check_collinearity(best_mod_f_p)

# --- all models ---
ms_f_p <- model.sel(dd_f_p)

all_models_f_p <- as.data.frame(ms_f_p)

# --- stars function ---
stars <- function(p) {
  ifelse(p < 0.001, "***",
         ifelse(p < 0.01, "**",
                ifelse(p < 0.05, "*", "")))
}

# --- p-values from best model ---
pvals <- summary(best_mod_f_p)$coefficients[, "Pr(>|t|)"]

# --- clean table ---
fmt <- function(x, pvals, name) {
  
  p <- pvals[name]
  
  if (is.na(p) && grepl(":", name)) {
    parts <- strsplit(name, ":")[[1]]
    alt_name <- paste0(parts[2], ":", parts[1])
    if (alt_name %in% names(pvals)) {
      p <- pvals[alt_name]
    }
  }
  
  ifelse(
    is.na(x),
    "—",
    paste0(
      round(x, 2),
      ifelse(!is.na(p), stars(p), "")
    )
  )
}

table_clean_f_p <- all_models_f_p %>%
  transmute(
    delta = round(delta, 2),
    weight = signif(weight, 2),
    
    Plasticity = fmt(clim_anomaly_tw90, pvals, "clim_anomaly_tw90"),
    Photoperiod = fmt(photo_tw90, pvals, "photo_tw90"),
    Background = fmt(clim_background_tw90, pvals, "clim_background_tw90"),
    Predictability = fmt(clim_predictability_tw90, pvals, "clim_predictability_tw90"),
    Autocorrelation = fmt(clim_autocorr_tw90, pvals, "clim_autocorr_tw90"),
    Trend = fmt(clim_trend_tw90, pvals, "clim_trend_tw90"),
    
    `Plast×Photo` = fmt(`clim_anomaly_tw90:photo_tw90`, pvals, "clim_anomaly_tw90:photo_tw90"),
    `Plast×Background` = fmt(`clim_anomaly_tw90:clim_background_tw90`, pvals, "clim_anomaly_tw90:clim_background_tw90"),
    `Plast×Predictability` = fmt(`clim_anomaly_tw90:clim_predictability_tw90`, pvals, "clim_anomaly_tw90:clim_predictability_tw90"),
    `Plast×Autocorr` = fmt(`clim_anomaly_tw90:clim_autocorr_tw90`, pvals, "clim_anomaly_tw90:clim_autocorr_tw90"),
    `Plast×Trend` = fmt(`clim_anomaly_tw90:clim_trend_tw90`, pvals, "clim_anomaly_tw90:clim_trend_tw90")
  )

table_clean_f_p


write_xlsx(
  table_clean_f_p,
  here::here("output", "onset_first_peak_plasticity_table_result.xlsx")
)

# ---------------------------------------------------------------------------- #
# Plot interactions of the best model #

mod_final_f_p <- lmer(
  FIRST_PEAK ~ clim_anomaly_tw90 *
    (photo_tw90 +
       clim_background_tw90 +
       clim_predictability_tw90 + clim_autocorr_tw90) +
    (1 | SITE_ID) +
    (1 + clim_anomaly_tw90 | SPECIES),
  
  data = d_f_p, REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e6)
  )
)

summary(mod_final_f_p)
performance::check_collinearity(mod_final_f_p)

# sequences for predictor variables

photo_seq <- seq(-2, 2, length.out = 100)
trend_seq    <- seq(-2, 2, length.out = 100)
pred_seq  <- seq(-2, 2, length.out = 100)
auto_seq  <- seq(-2, 2, length.out = 100)
back_seq  <- seq(-2, 2, length.out = 100)


# get slopes of plasticity across gradients

get_slope_df <- function(mod, var_name, var_seq, label,
                         anomaly_name = "clim_anomaly_tw90") {
  
  b <- fixef(mod)
  V <- vcov(mod)
  
  coef_main <- anomaly_name
  
  coef_int <- paste0(coef_main, ":", var_name)
  if (!coef_int %in% names(b)) {
    coef_int <- paste0(var_name, ":", coef_main)
  }
  
  if (!all(c(coef_main, coef_int) %in% names(b))) {
    stop(paste("No s'han trobat coeficients per:", var_name))
  }
  
  coef_names <- c(coef_main, coef_int)
  
  get_slope <- function(x) {
    
    g <- c(1, x)
    
    slope <- sum(b[coef_names] * g)
    
    se <- as.numeric(sqrt(t(g) %*% V[coef_names, coef_names] %*% g))
    
    c(slope = slope, se = se)
  }
  
  out <- do.call(rbind, lapply(var_seq, get_slope))
  
  data.frame(
    x = var_seq,
    slope = out[, "slope"],
    lower = out[, "slope"] - 1.96 * out[, "se"],
    upper = out[, "slope"] + 1.96 * out[, "se"],
    variable = label
  )
}

df_photo <- get_slope_df(
  mod_final_f_p,
  "photo_tw90",
  photo_seq,
  "Photoperiod"
)


df_pred <- get_slope_df(
  mod_final_f_p,
  "clim_predictability_tw90",
  pred_seq,
  "Temperature stability"
)

df_auto <- get_slope_df(
  mod_final_f_p,
  "clim_autocorr_tw90",
  auto_seq,
  "Temperature autocorrelation"
)

df_back <- get_slope_df(
  mod_final_f_p,
  "clim_background_tw90",
  back_seq,
  "Mean temperature"
)

# combine data frames

df_all <- dplyr::bind_rows(
  df_photo,
  df_pred,
  df_auto,
  df_back
)

# order facets

df_all$variable <- factor(
  df_all$variable,
  levels = c(
    "Photoperiod",
    "Temperature trend",
    "Temperature stability",
    "Temperature autocorrelation",
    "Mean temperature"
  )
)

# plot

first_peak_plasticity_interactions_plot <- ggplot(df_all, aes(x, slope, color = variable, fill = variable)) +
  
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.15, color = NA) +
  
  geom_line(size = 0.8, alpha = 0.6) +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  
  scale_color_manual(values = c(
    "Photoperiod" = "#1b9e77",
    "Temperature stability" = "#7570b3",
    "Temperature autocorrelation" = "#E6AB02",
    "Mean temperature" = "#d95f02"
    
  )) +
  
  scale_fill_manual(values = c(
    "Photoperiod" = "#1b9e77",
    "Temperature stability" = "#7570b3",
    "Temperature autocorrelation" = "#E6AB02",
    "Mean temperature" = "#d95f02"
  )) +
  
  coord_cartesian(ylim = c(-7, -2.5)) +
  
  theme_classic(base_family = "Garamond", base_size = 14) +
  
  labs(
    x = "Environmental gradient",
    y = "Phenological plasticity\n(slope of first peak day vs temperature anomaly)",
    color = "Environmental gradient",
    fill = "Environmental gradient"
  )
first_peak_plasticity_interactions_plot


ggsave(
  filename = here("output", "figures", "first_peak_plasticity_interactions_plot.png"),
  plot = first_peak_plasticity_interactions_plot,
  width = 7,
  height = 5,
  dpi = 300
)


# ---------------------------------------------------------------------------- #
#### Peak (max) #### 

df_peakday <- subset(df, df$pheno_type == "PEAKDAY")

# Anomaly time-window selection # 

options(na.action = "na.omit")

mod_30_p_d <- lmer(
  PEAKDAY ~ clim_anomaly_tw30 +
    (1 | SITE_ID) + (1 + clim_anomaly_tw30 | SPECIES),
  data = df_peakday, REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e6)
  )
)

mod_60_p_d <- lmer(
  PEAKDAY ~ clim_anomaly_tw60 +
    (1 | SITE_ID) + (1 + clim_anomaly_tw60 | SPECIES),
  data = df_peakday, REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e6)
  )
)

mod_90_p_d <- lmer(
  PEAKDAY ~ clim_anomaly_tw90 +
    (1 | SITE_ID) + (1 + clim_anomaly_tw90 | SPECIES),
  data = df_peakday, REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e6)
  )
)

AIC(mod_30_p_d, mod_60_p_d, mod_90_p_d) ### Best time-window is 60 days before onset
summary(mod_30_p_d)
summary(mod_60_p_d)
summary(mod_90_p_d)
fixef(mod_30_p_d)["clim_anomaly_tw30"]
fixef(mod_60_p_d)["clim_anomaly_tw60"]
fixef(mod_90_p_d)["clim_anomaly_tw90"]


# Plot effects

extract_effect <- function(mod, var, label) {
  tidy(mod) %>%
    filter(term == var) %>%
    mutate(window = label)
}

df_eff <- bind_rows(
  extract_effect(mod_30_p_d, "clim_anomaly_tw30", "30 days"),
  extract_effect(mod_60_p_d, "clim_anomaly_tw60", "60 days"),
  extract_effect(mod_90_p_d, "clim_anomaly_tw90", "90 days")
) %>%
  mutate(
    lower = estimate - 1.96 * std.error,
    upper = estimate + 1.96 * std.error
  )

ggplot(df_eff, aes(x = window, y = estimate)) +
  
  geom_point(size = 3) +
  
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0.1) +
  
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  theme_classic(base_family = "Garamond") +
  
  labs(
    x = "Time window",
    y = "Phenological sensitivity\n(slope of onset vs temperature anomaly)"
  )

make_pred <- function(mod, var, label) {
  
  newdat <- data.frame(
    anomaly = seq(-2, 2, length.out = 100)
  )
  
  names(newdat) <- var
  
  pred <- predict(
    mod,
    newdata = newdat,
    re.form = NA,   # important!
    se.fit = TRUE
  )
  
  data.frame(
    anomaly = newdat[[var]],
    fit = pred$fit,
    se = pred$se.fit,
    window = label
  )
}

df_pred <- dplyr::bind_rows(
  make_pred(mod_30_p_d, "clim_anomaly_tw30", "30 days"),
  make_pred(mod_60_p_d, "clim_anomaly_tw60", "60 days"),
  make_pred(mod_90_p_d, "clim_anomaly_tw90", "90 days")
)

df_pred <- df_pred %>%
  dplyr::mutate(
    lower = fit - 1.96 * se,
    upper = fit + 1.96 * se
  )

ggplot(df_pred, aes(anomaly, fit,
                    color = window,
                    fill = window)) +
  
  geom_line(size = 1.2) +
  
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.2, color = NA) +
  
  theme_classic(base_family = "Garamond") +
  
  labs(
    x = "Temperature anomaly",
    y = "Predicted first peak day",
    color = "Time window",
    fill = "Time window"
  )


# ---------------------------------------------------------------------------- #
# Model selection #

d_p_d <- df_peakday %>%
  dplyr::filter(
    !is.na(PEAKDAY),
    !is.na(clim_anomaly_tw90),
    !is.na(photo_tw90),
    !is.na(clim_background_tw90),
    !is.na(clim_predictability_tw90),
    !is.na(clim_autocorr_tw90),
    !is.na(clim_trend_tw90)
  )

global_mod_p_d <- lmer(
  PEAKDAY ~ clim_anomaly_tw90 *
    (photo_tw90 +
       clim_background_tw90 +
       clim_predictability_tw90 +
       clim_autocorr_tw90 +
       clim_trend_tw90) +
    (1 | SITE_ID) +
    (1 + clim_anomaly_tw90 | SPECIES),
  data = d_p_d, REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e6)
  )
)

options(na.action = "na.fail")

dd_p_d <- dredge(
  global_mod_p_d,
  fixed = "clim_anomaly_tw90",
  subset =
    dc("photo_tw90", "clim_anomaly_tw90:photo_tw90") &
    dc("clim_background_tw90", "clim_anomaly_tw90:clim_background_tw90") &
    dc("clim_predictability_tw90", "clim_anomaly_tw90:clim_predictability_tw90") &
    dc("clim_autocorr_tw90", "clim_anomaly_tw90:clim_autocorr_tw90") &
    dc("clim_trend_tw90", "clim_anomaly_tw90:clim_trend_tw90")
)

model.sel(dd_p_d)

best_mod_p_d <- get.models(dd_p_d, subset = 1)[[1]]
summary(best_mod_p_d)
check_collinearity(best_mod_p_d)

# --- all models ---
ms_p_d <- model.sel(dd_p_d)

all_models_p_d <- as.data.frame(ms_p_d)

# --- stars function ---
stars <- function(p) {
  ifelse(p < 0.001, "***",
         ifelse(p < 0.01, "**",
                ifelse(p < 0.05, "*", "")))
}

# --- p-values from best model ---
pvals <- summary(best_mod_p_d)$coefficients[, "Pr(>|t|)"]

# --- clean table ---
fmt <- function(x, pvals, name) {
  
  p <- pvals[name]
  
  if (is.na(p) && grepl(":", name)) {
    parts <- strsplit(name, ":")[[1]]
    alt_name <- paste0(parts[2], ":", parts[1])
    if (alt_name %in% names(pvals)) {
      p <- pvals[alt_name]
    }
  }
  
  ifelse(
    is.na(x),
    "—",
    paste0(
      round(x, 2),
      ifelse(!is.na(p), stars(p), "")
    )
  )
}

table_clean_p_d <- all_models_p_d %>%
  transmute(
    delta = round(delta, 2),
    weight = signif(weight, 2),
    
    Plasticity = fmt(clim_anomaly_tw90, pvals, "clim_anomaly_tw90"),
    Photoperiod = fmt(photo_tw90, pvals, "photo_tw90"),
    Background = fmt(clim_background_tw90, pvals, "clim_background_tw90"),
    Predictability = fmt(clim_predictability_tw90, pvals, "clim_predictability_tw90"),
    Autocorrelation = fmt(clim_autocorr_tw90, pvals, "clim_autocorr_tw90"),
    Trend = fmt(clim_trend_tw90, pvals, "clim_trend_tw90"),
    
    `Plast×Photo` = fmt(`clim_anomaly_tw90:photo_tw90`, pvals, "clim_anomaly_tw90:photo_tw90"),
    `Plast×Background` = fmt(`clim_anomaly_tw90:clim_background_tw90`, pvals, "clim_anomaly_tw90:clim_background_tw90"),
    `Plast×Predictability` = fmt(`clim_anomaly_tw90:clim_predictability_tw90`, pvals, "clim_anomaly_tw90:clim_predictability_tw90"),
    `Plast×Autocorr` = fmt(`clim_anomaly_tw90:clim_autocorr_tw90`, pvals, "clim_anomaly_tw90:clim_autocorr_tw90"),
    `Plast×Trend` = fmt(`clim_anomaly_tw90:clim_trend_tw90`, pvals, "clim_anomaly_tw90:clim_trend_tw90")
  )

table_clean_p_d


write_xlsx(
  table_clean_p_d,
  here::here("output", "peakday_plasticity_table_result.xlsx")
)

# ---------------------------------------------------------------------------- #
# Plot interactions of the best model #

mod_final_p_d <- lmer(
  PEAKDAY ~ clim_anomaly_tw90 *
    (photo_tw90 + clim_autocorr_tw90) +
    (1 | SITE_ID) +
    (1 + clim_anomaly_tw90 | SPECIES),
  
  data = d_p_d, REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e6)
  )
)

summary(mod_final_p_d)
performance::check_collinearity(mod_final_p_d)

# sequences for predictor variables

photo_seq <- seq(-2, 2, length.out = 100)
auto_seq  <- seq(-2, 2, length.out = 100)

# get slopes of plasticity across gradients

get_slope_df <- function(mod, var_name, var_seq, label,
                         anomaly_name = "clim_anomaly_tw90") {
  
  b <- fixef(mod)
  V <- vcov(mod)
  
  coef_main <- anomaly_name
  
  coef_int <- paste0(coef_main, ":", var_name)
  if (!coef_int %in% names(b)) {
    coef_int <- paste0(var_name, ":", coef_main)
  }
  
  if (!all(c(coef_main, coef_int) %in% names(b))) {
    stop(paste("No s'han trobat coeficients per:", var_name))
  }
  
  coef_names <- c(coef_main, coef_int)
  
  get_slope <- function(x) {
    
    g <- c(1, x)
    
    slope <- sum(b[coef_names] * g)
    
    se <- as.numeric(sqrt(t(g) %*% V[coef_names, coef_names] %*% g))
    
    c(slope = slope, se = se)
  }
  
  out <- do.call(rbind, lapply(var_seq, get_slope))
  
  data.frame(
    x = var_seq,
    slope = out[, "slope"],
    lower = out[, "slope"] - 1.96 * out[, "se"],
    upper = out[, "slope"] + 1.96 * out[, "se"],
    variable = label
  )
}

df_photo <- get_slope_df(
  mod_final_p_d,
  "photo_tw90",
  photo_seq,
  "Photoperiod"
)



df_auto <- get_slope_df(
  mod_final_p_d,
  "clim_autocorr_tw90",
  auto_seq,
  "Temperature autocorrelation"
)


# combine data frames

df_all <- dplyr::bind_rows(
  df_photo,
  df_auto

)

# order facets

df_all$variable <- factor(
  df_all$variable,
  levels = c(
    "Photoperiod",
    "Temperature trend",
    "Temperature stability",
    "Temperature autocorrelation",
    "Mean temperature"
  )
)

# plot

peakday_plasticity_interactions_plot <- ggplot(df_all, aes(x, slope, color = variable, fill = variable)) +
  
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.15, color = NA) +
  
  geom_line(size = 0.8, alpha = 0.6) +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  
  scale_color_manual(values = c(
    "Photoperiod" = "#1b9e77",
    "Temperature stability" = "#7570b3",
    "Temperature autocorrelation" = "#E6AB02",
    "Mean temperature" = "#d95f02"
    
  )) +
  
  scale_fill_manual(values = c(
    "Photoperiod" = "#1b9e77",
    "Temperature stability" = "#7570b3",
    "Temperature autocorrelation" = "#E6AB02",
    "Mean temperature" = "#d95f02"
  )) +
  
  coord_cartesian(ylim = c(-6.5, -2)) +
  
  theme_classic(base_family = "Garamond", base_size = 14) +
  
  labs(
    x = "Environmental gradient",
    y = "Phenological plasticity\n(slope of onset vs temperature anomaly)",
    color = "Environmental gradient",
    fill = "Environmental gradient"
  )
peakday_plasticity_interactions_plot


ggsave(
  filename = here("output", "figures", "peakday_plasticity_interactions_plot.png"),
  plot = peakday_plasticity_interactions_plot,
  width = 7,
  height = 5,
  dpi = 300
)


# ---------------------------------------------------------------------------- #
#### Offset (mean) #### 

df_offset_mean <- subset(df, df$pheno_type == "OFFSET_mean")

# Anomaly time-window selection # 

#1. Univoltine species 

options(na.action = "na.omit")

mod_30_off_uni <- lmer(
  OFFSET_mean ~ clim_anomaly_tw30 +
    (1 | SITE_ID) + (1 + clim_anomaly_tw30 | SPECIES),
  data = df_offset_mean |> dplyr::filter(voltinism == "univoltine"), 
  REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e6)
  )
)

mod_60_off_uni <- lmer(
  OFFSET_mean ~ clim_anomaly_tw60 +
    (1 | SITE_ID) + (1 + clim_anomaly_tw60 | SPECIES),
  data = df_offset_mean |> dplyr::filter(voltinism == "univoltine"), 
  REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e6)
  )
)

mod_90_off_uni <- lmer(
  OFFSET_mean ~ clim_anomaly_tw90 +
    (1 | SITE_ID) + (1 + clim_anomaly_tw90 | SPECIES),
  data = df_offset_mean |> dplyr::filter(voltinism == "univoltine"), 
  REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e6)
  )
)

AIC(mod_30_off_uni, mod_60_off_uni, mod_90_off_uni) ### Best time-window is 60 days before offset
summary(mod_30_off_uni)
summary(mod_60_off_uni)
summary(mod_90_off_uni)
fixef(mod_30_off_uni)["clim_anomaly_tw30"]
fixef(mod_60_off_uni)["clim_anomaly_tw60"]
fixef(mod_90_off_uni)["clim_anomaly_tw90"]


# Plot effects

extract_effect <- function(mod, var, label) {
  tidy(mod) %>%
    filter(term == var) %>%
    mutate(window = label)
}

df_eff <- bind_rows(
  extract_effect(mod_30_off_uni, "clim_anomaly_tw30", "30 days"),
  extract_effect(mod_60_off_uni, "clim_anomaly_tw60", "60 days"),
  extract_effect(mod_90_off_uni, "clim_anomaly_tw90", "90 days")
) %>%
  mutate(
    lower = estimate - 1.96 * std.error,
    upper = estimate + 1.96 * std.error
  )

ggplot(df_eff, aes(x = window, y = estimate)) +
  
  geom_point(size = 3) +
  
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0.1) +
  
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  theme_classic(base_family = "Garamond") +
  
  labs(
    x = "Time window",
    y = "Phenological sensitivity\n(slope of offset vs temperature anomaly)"
  )

make_pred <- function(mod, var, label) {
  
  newdat <- data.frame(
    anomaly = seq(-2, 2, length.out = 100)
  )
  
  names(newdat) <- var
  
  pred <- predict(
    mod,
    newdata = newdat,
    re.form = NA,   # important!
    se.fit = TRUE
  )
  
  data.frame(
    anomaly = newdat[[var]],
    fit = pred$fit,
    se = pred$se.fit,
    window = label
  )
}

df_pred <- dplyr::bind_rows(
  make_pred(mod_30_off_uni, "clim_anomaly_tw30", "30 days"),
  make_pred(mod_60_off_uni, "clim_anomaly_tw60", "60 days"),
  make_pred(mod_90_off_uni, "clim_anomaly_tw90", "90 days")
)

df_pred <- df_pred %>%
  dplyr::mutate(
    lower = fit - 1.96 * se,
    upper = fit + 1.96 * se
  )

ggplot(df_pred, aes(anomaly, fit,
                    color = window,
                    fill = window)) +
  
  geom_line(size = 1.2) +
  
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.2, color = NA) +
  
  theme_classic(base_family = "Garamond") +
  
  labs(
    x = "Temperature anomaly",
    y = "Predicted offset",
    color = "Time window",
    fill = "Time window"
  )


#1. Multivoltine species 

options(na.action = "na.omit")

mod_30_off_multi <- lmer(
  OFFSET_mean ~ ONSET_mean + clim_anomaly_tw30 +
    (1 | SITE_ID) + (1 + clim_anomaly_tw30 | SPECIES),
  data = df_offset_mean |> dplyr::filter(voltinism == "multivoltine"), 
  REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e6)
  )
)

mod_60_off_multi <- lmer(
  OFFSET_mean ~ ONSET_mean + clim_anomaly_tw60 +
    (1 | SITE_ID) + (1 + clim_anomaly_tw60 | SPECIES),
  data = df_offset_mean |> dplyr::filter(voltinism == "multivoltine"), 
  REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e6)
  )
)

mod_90_off_multi <- lmer(
  OFFSET_mean ~ ONSET_mean + clim_anomaly_tw90 +
    (1 | SITE_ID) + (1 + clim_anomaly_tw90 | SPECIES),
  data = df_offset_mean |> dplyr::filter(voltinism == "multivoltine"), 
  REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e6)
  )
)

AIC(mod_30_off_multi, mod_60_off_multi, mod_90_off_multi) ### Best time-window is 60 days before offset
summary(mod_30_off_multi)
summary(mod_60_off_multi)
summary(mod_90_off_multi)
fixef(mod_30_off_multi)["clim_anomaly_tw30"]
fixef(mod_60_off_multi)["clim_anomaly_tw60"]
fixef(mod_90_off_multi)["clim_anomaly_tw90"]


# Plot effects

df_eff <- bind_rows(
  extract_effect(mod_30_off_multi, "clim_anomaly_tw30", "30 days"),
  extract_effect(mod_60_off_multi, "clim_anomaly_tw60", "60 days"),
  extract_effect(mod_90_off_multi, "clim_anomaly_tw90", "90 days")
) %>%
  mutate(
    lower = estimate - 1.96 * std.error,
    upper = estimate + 1.96 * std.error
  )

ggplot(df_eff, aes(x = window, y = estimate)) +
  
  geom_point(size = 3) +
  
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0.1) +
  
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  theme_classic(base_family = "Garamond") +
  
  labs(
    x = "Time window",
    y = "Phenological sensitivity\n(slope of offset vs temperature anomaly)"
  )

make_pred <- function(model, anomaly_var, label) {
  
  xseq <- seq(
    min(model.frame(model)[[anomaly_var]], na.rm = TRUE),
    max(model.frame(model)[[anomaly_var]], na.rm = TRUE),
    length.out = 100
  )
  
  nd <- data.frame(
    ONSET_mean = mean(model.frame(model)$ONSET_mean, na.rm = TRUE),
    photo_tw90 = 0,
    clim_background_tw90 = 0,
    clim_predictability_tw90 = 0,
    clim_autocorr_tw90 = 0,
    clim_trend_tw90 = 0
  )
  
  nd <- nd[rep(1, length(xseq)), ]
  nd[[anomaly_var]] <- xseq
  
  pred <- predict(model, newdata = nd, re.form = NA, se.fit = TRUE)
  
  tibble::tibble(
    x = xseq,
    fit = pred$fit,
    se = pred$se.fit,
    window = label
  )
}

df_pred <- dplyr::bind_rows(
  make_pred(mod_30_off_multi, "clim_anomaly_tw30", "30 days"),
  make_pred(mod_60_off_multi, "clim_anomaly_tw60", "60 days"),
  make_pred(mod_90_off_multi, "clim_anomaly_tw90", "90 days")
) |>
  dplyr::mutate(
    lower = fit - 1.96 * se,
    upper = fit + 1.96 * se
  )


ggplot(df_pred, aes(x, fit,
                    color = window,
                    fill = window)) +
  
  geom_line(size = 1.2) +
  
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.2, color = NA) +
  
  theme_classic(base_family = "Garamond") +
  
  labs(
    x = "Temperature anomaly",
    y = "Predicted offset",
    color = "Time window",
    fill = "Time window"
  )



# ---------------------------------------------------------------------------- #
# Model selection #

d_off_m <- df_offset_mean %>%
  dplyr::filter(
    !is.na(OFFSET_mean),
    !is.na(clim_anomaly_tw90),
    !is.na(photo_tw90),
    !is.na(clim_background_tw90),
    !is.na(clim_predictability_tw90),
    !is.na(clim_autocorr_tw90),
    !is.na(clim_trend_tw90)
  )


# Univoltine model
global_mod_off_m_uni <- lmer(
  OFFSET_mean ~ ONSET_mean + clim_anomaly_tw90 *
    (photo_tw90 +
       clim_background_tw90 +
       clim_predictability_tw90 +
       clim_autocorr_tw90 +
       clim_trend_tw90) +
    (1 | SITE_ID) +
    (1 + clim_anomaly_tw90 | SPECIES),
  data = d_off_m |> dplyr::filter(voltinism == "univoltine"),
  REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e6)
  )
)

# Multivoltine model
global_mod_off_m_multi <- lmer(
  OFFSET_mean ~ ONSET_mean + clim_anomaly_tw90 *
    (photo_tw90 +
       clim_background_tw90 +
       clim_predictability_tw90 +
       clim_autocorr_tw90 +
       clim_trend_tw90) +
    (1 | SITE_ID) +
    (1 + clim_anomaly_tw90 | SPECIES),
  data = d_off_m |> dplyr::filter(voltinism == "multivoltine"),
  REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e6)
  )
)
options(na.action = "na.fail")


# Model selection univoltine
dd_off_m_uni <- dredge(
  global_mod_off_m_uni,
  fixed = "clim_anomaly_tw90",
  subset =
    dc("photo_tw90", "clim_anomaly_tw90:photo_tw90") &
    dc("clim_background_tw90", "clim_anomaly_tw90:clim_background_tw90") &
    dc("clim_predictability_tw90", "clim_anomaly_tw90:clim_predictability_tw90") &
    dc("clim_autocorr_tw90", "clim_anomaly_tw90:clim_autocorr_tw90") &
    dc("clim_trend_tw90", "clim_anomaly_tw90:clim_trend_tw90")
)

model.sel(dd_off_m_uni)

best_mod_off_m_uni <- get.models(dd_off_m_uni, subset = 1)[[1]]
summary(best_mod_off_m_uni)
check_collinearity(best_mod_off_m_uni)

# --- all models ---
ms_off_m_uni <- model.sel(dd_off_m_uni)

all_models_off_m_uni <- as.data.frame(ms_off_m_uni)

# --- stars function ---
stars <- function(p) {
  ifelse(p < 0.001, "***",
         ifelse(p < 0.01, "**",
                ifelse(p < 0.05, "*", "")))
}

# --- p-values from best model ---
pvals <- summary(best_mod_off_m_uni)$coefficients[, "Pr(>|t|)"]

# --- clean table ---
fmt <- function(x, pvals, name) {
  
  p <- pvals[name]
  
  if (is.na(p) && grepl(":", name)) {
    parts <- strsplit(name, ":")[[1]]
    alt_name <- paste0(parts[2], ":", parts[1])
    if (alt_name %in% names(pvals)) {
      p <- pvals[alt_name]
    }
  }
  
  ifelse(
    is.na(x),
    "—",
    paste0(
      round(x, 2),
      ifelse(!is.na(p), stars(p), "")
    )
  )
}

table_clean_off_m_uni <- all_models_off_m_uni %>%
  transmute(
    delta = round(delta, 2),
    weight = signif(weight, 2),
    
    Plasticity = fmt(clim_anomaly_tw90, pvals, "clim_anomaly_tw90"),
    Photoperiod = fmt(photo_tw90, pvals, "photo_tw90"),
    Background = fmt(clim_background_tw90, pvals, "clim_background_tw90"),
    Predictability = fmt(clim_predictability_tw90, pvals, "clim_predictability_tw90"),
    Autocorrelation = fmt(clim_autocorr_tw90, pvals, "clim_autocorr_tw90"),
    Trend = fmt(clim_trend_tw90, pvals, "clim_trend_tw90"),
    
    `Plast×Photo` = fmt(`clim_anomaly_tw90:photo_tw90`, pvals, "clim_anomaly_tw90:photo_tw90"),
    `Plast×Background` = fmt(`clim_anomaly_tw90:clim_background_tw90`, pvals, "clim_anomaly_tw90:clim_background_tw90"),
    `Plast×Predictability` = fmt(`clim_anomaly_tw90:clim_predictability_tw90`, pvals, "clim_anomaly_tw90:clim_predictability_tw90"),
    `Plast×Autocorr` = fmt(`clim_anomaly_tw90:clim_autocorr_tw90`, pvals, "clim_anomaly_tw90:clim_autocorr_tw90"),
    `Plast×Trend` = fmt(`clim_anomaly_tw90:clim_trend_tw90`, pvals, "clim_anomaly_tw90:clim_trend_tw90")
  )

table_clean_off_m_uni


write_xlsx(
  table_clean_off_m_uni,
  here::here("output", "uni_offset_mean_plasticity_table_result.xlsx")
)


# Model selection multivoltine
dd_off_m_multi <- dredge(
  global_mod_off_m_multi,
  fixed = "clim_anomaly_tw90",
  subset =
    dc("photo_tw90", "clim_anomaly_tw90:photo_tw90") &
    dc("clim_background_tw90", "clim_anomaly_tw90:clim_background_tw90") &
    dc("clim_predictability_tw90", "clim_anomaly_tw90:clim_predictability_tw90") &
    dc("clim_autocorr_tw90", "clim_anomaly_tw90:clim_autocorr_tw90") &
    dc("clim_trend_tw90", "clim_anomaly_tw90:clim_trend_tw90")
)

model.sel(dd_off_m_multi)

best_mod_off_m_multi <- get.models(dd_off_m_multi, subset = 1)[[1]]
summary(best_mod_off_m_multi)
check_collinearity(best_mod_off_m_multi)

# --- all models ---
ms_off_m_multi <- model.sel(dd_off_m_multi)

all_models_off_m_multi <- as.data.frame(ms_off_m_multi)

# --- stars function ---
stars <- function(p) {
  ifelse(p < 0.001, "***",
         ifelse(p < 0.01, "**",
                ifelse(p < 0.05, "*", "")))
}

# --- p-values from best model ---
pvals <- summary(best_mod_off_m_multi)$coefficients[, "Pr(>|t|)"]

# --- clean table ---
fmt <- function(x, pvals, name) {
  
  p <- pvals[name]
  
  if (is.na(p) && grepl(":", name)) {
    parts <- strsplit(name, ":")[[1]]
    alt_name <- paste0(parts[2], ":", parts[1])
    if (alt_name %in% names(pvals)) {
      p <- pvals[alt_name]
    }
  }
  
  ifelse(
    is.na(x),
    "—",
    paste0(
      round(x, 2),
      ifelse(!is.na(p), stars(p), "")
    )
  )
}

table_clean_off_m_multi <- all_models_off_m_multi %>%
  transmute(
    delta = round(delta, 2),
    weight = signif(weight, 2),
    
    Plasticity = fmt(clim_anomaly_tw90, pvals, "clim_anomaly_tw90"),
    Photoperiod = fmt(photo_tw90, pvals, "photo_tw90"),
    Background = fmt(clim_background_tw90, pvals, "clim_background_tw90"),
    Predictability = fmt(clim_predictability_tw90, pvals, "clim_predictability_tw90"),
    Autocorrelation = fmt(clim_autocorr_tw90, pvals, "clim_autocorr_tw90"),
    Trend = fmt(clim_trend_tw90, pvals, "clim_trend_tw90"),
    
    `Plast×Photo` = fmt(`clim_anomaly_tw90:photo_tw90`, pvals, "clim_anomaly_tw90:photo_tw90"),
    `Plast×Background` = fmt(`clim_anomaly_tw90:clim_background_tw90`, pvals, "clim_anomaly_tw90:clim_background_tw90"),
    `Plast×Predictability` = fmt(`clim_anomaly_tw90:clim_predictability_tw90`, pvals, "clim_anomaly_tw90:clim_predictability_tw90"),
    `Plast×Autocorr` = fmt(`clim_anomaly_tw90:clim_autocorr_tw90`, pvals, "clim_anomaly_tw90:clim_autocorr_tw90"),
    `Plast×Trend` = fmt(`clim_anomaly_tw90:clim_trend_tw90`, pvals, "clim_anomaly_tw90:clim_trend_tw90")
  )

table_clean_off_m_multi


write_xlsx(
  table_clean_off_m_multi,
  here::here("output", "multi_offset_mean_plasticity_table_result.xlsx")
)


# ---------------------------------------------------------------------------- #
# Plot interactions of the best model #

# ------------ Univoltine species -------------- #

mod_final_off_m_uni <- lmer(
  OFFSET_mean ~ ONSET_mean + clim_anomaly_tw90 *
    (photo_tw90 +
       clim_background_tw90 +
       clim_predictability_tw90) +
    (1 | SITE_ID) +
    (1 + clim_anomaly_tw90 | SPECIES),
  data = d_off_m |> dplyr::filter(voltinism == "univoltine"),
  REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e6)
  )
)

summary(mod_final_off_m_uni)
performance::check_collinearity(mod_final_off_m_uni)

# sequences for predictor variables

photo_seq <- seq(-2, 2, length.out = 100)
trend_seq    <- seq(-2, 2, length.out = 100)
pred_seq  <- seq(-2, 2, length.out = 100)
auto_seq  <- seq(-2, 2, length.out = 100)
back_seq  <- seq(-2, 2, length.out = 100)


# get slopes of plasticity across gradients

get_slope_df <- function(mod, var_name, var_seq, label,
                         anomaly_name = "clim_anomaly_tw90") {
  
  b <- fixef(mod)
  V <- vcov(mod)
  
  coef_main <- anomaly_name
  
  coef_int <- paste0(coef_main, ":", var_name)
  if (!coef_int %in% names(b)) {
    coef_int <- paste0(var_name, ":", coef_main)
  }
  
  if (!all(c(coef_main, coef_int) %in% names(b))) {
    stop(paste("No s'han trobat coeficients per:", var_name))
  }
  
  coef_names <- c(coef_main, coef_int)
  
  get_slope <- function(x) {
    
    g <- c(1, x)
    
    slope <- sum(b[coef_names] * g)
    
    se <- as.numeric(sqrt(t(g) %*% V[coef_names, coef_names] %*% g))
    
    c(slope = slope, se = se)
  }
  
  out <- do.call(rbind, lapply(var_seq, get_slope))
  
  data.frame(
    x = var_seq,
    slope = out[, "slope"],
    lower = out[, "slope"] - 1.96 * out[, "se"],
    upper = out[, "slope"] + 1.96 * out[, "se"],
    variable = label
  )
}

df_photo <- get_slope_df(
  mod_final_off_m_uni,
  "photo_tw90",
  photo_seq,
  "Photoperiod"
)


df_pred <- get_slope_df(
  mod_final_off_m_uni,
  "clim_predictability_tw90",
  pred_seq,
  "Temperature stability"
)

df_auto <- get_slope_df(
  mod_final_off_m_uni,
  "clim_autocorr_tw90",
  auto_seq,
  "Temperature autocorrelation"
)

df_back <- get_slope_df(
  mod_final_off_m_uni,
  "clim_background_tw90",
  back_seq,
  "Mean temperature"
)

# combine data frames

df_all <- dplyr::bind_rows(
  df_photo,
  df_pred,
  df_back
)

# order facets

df_all$variable <- factor(
  df_all$variable,
  levels = c(
    "Photoperiod",
    "Temperature trend",
    "Temperature stability",
    "Temperature autocorrelation",
    "Mean temperature"
  )
)

# plot

univoltine_offset_plasticity_interactions_plot <- ggplot(df_all, aes(x, slope, color = variable, fill = variable)) +
  
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.15, color = NA) +
  
  geom_line(size = 0.8, alpha = 0.6) +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  
  scale_color_manual(values = c(
    "Photoperiod" = "#1b9e77",
    "Temperature stability" = "#7570b3",
    "Temperature autocorrelation" = "#E6AB02",
    "Mean temperature" = "#d95f02"
    
  )) +
  
  scale_fill_manual(values = c(
    "Photoperiod" = "#1b9e77",
    "Temperature stability" = "#7570b3",
    "Temperature autocorrelation" = "#E6AB02",
    "Mean temperature" = "#d95f02"
  )) +
  
  coord_cartesian(ylim = c(-5.5, -1.5)) +
  
  theme_classic(base_family = "Garamond", base_size = 14) +
  
  labs(
    x = "Environmental gradient",
    y = "Phenological plasticity\n(slope of onset vs temperature anomaly)",
    color = "Environmental gradient",
    fill = "Environmental gradient"
  )
univoltine_offset_plasticity_interactions_plot


ggsave(
  filename = here("output", "figures", "univoltine_offset_plasticity_interactions_plot.png"),
  plot = univoltine_offset_plasticity_interactions_plot,
  width = 7,
  height = 5,
  dpi = 300
)


# ------------ Multivoltine species -------------- #

mod_final_off_m_multi <- lmer(
  OFFSET_mean ~ ONSET_mean + clim_anomaly_tw90 *
    (photo_tw90 +
       clim_background_tw90 +
       clim_predictability_tw90
     + clim_autocorr_tw90) +
    (1 | SITE_ID) +
    (1 + clim_anomaly_tw90 | SPECIES),
  data = d_off_m |> dplyr::filter(voltinism == "multivoltine"),
  REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e6)
  )
)

summary(mod_final_off_m_multi)
performance::check_collinearity(mod_final_off_m_multi)

# sequences for predictor variables

photo_seq <- seq(-2, 2, length.out = 100)
trend_seq    <- seq(-2, 2, length.out = 100)
pred_seq  <- seq(-2, 2, length.out = 100)
auto_seq  <- seq(-2, 2, length.out = 100)
back_seq  <- seq(-2, 2, length.out = 100)


# get slopes of plasticity across gradients

get_slope_df <- function(mod, var_name, var_seq, label,
                         anomaly_name = "clim_anomaly_tw90") {
  
  b <- fixef(mod)
  V <- vcov(mod)
  
  coef_main <- anomaly_name
  
  coef_int <- paste0(coef_main, ":", var_name)
  if (!coef_int %in% names(b)) {
    coef_int <- paste0(var_name, ":", coef_main)
  }
  
  if (!all(c(coef_main, coef_int) %in% names(b))) {
    stop(paste("No s'han trobat coeficients per:", var_name))
  }
  
  coef_names <- c(coef_main, coef_int)
  
  get_slope <- function(x) {
    
    g <- c(1, x)
    
    slope <- sum(b[coef_names] * g)
    
    se <- as.numeric(sqrt(t(g) %*% V[coef_names, coef_names] %*% g))
    
    c(slope = slope, se = se)
  }
  
  out <- do.call(rbind, lapply(var_seq, get_slope))
  
  data.frame(
    x = var_seq,
    slope = out[, "slope"],
    lower = out[, "slope"] - 1.96 * out[, "se"],
    upper = out[, "slope"] + 1.96 * out[, "se"],
    variable = label
  )
}

df_photo <- get_slope_df(
  mod_final_off_m_multi,
  "photo_tw90",
  photo_seq,
  "Photoperiod"
)


df_pred <- get_slope_df(
  mod_final_off_m_multi,
  "clim_predictability_tw90",
  pred_seq,
  "Temperature stability"
)

df_auto <- get_slope_df(
  mod_final_off_m_multi,
  "clim_autocorr_tw90",
  auto_seq,
  "Temperature autocorrelation"
)

df_back <- get_slope_df(
  mod_final_off_m_multi,
  "clim_background_tw90",
  back_seq,
  "Mean temperature"
)

# combine data frames

df_all <- dplyr::bind_rows(
  df_photo,
  df_pred,
  df_back,
  df_auto
)

# order facets

df_all$variable <- factor(
  df_all$variable,
  levels = c(
    "Photoperiod",
    "Temperature trend",
    "Temperature stability",
    "Temperature autocorrelation",
    "Mean temperature"
  )
)

# plot

multivoltine_offset_plasticity_interactions_plot <- ggplot(df_all, aes(x, slope, color = variable, fill = variable)) +
  
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.15, color = NA) +
  
  geom_line(size = 0.8, alpha = 0.6) +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  
  scale_color_manual(values = c(
    "Photoperiod" = "#1b9e77",
    "Temperature stability" = "#7570b3",
    "Temperature autocorrelation" = "#E6AB02",
    "Mean temperature" = "#d95f02"
    
  )) +
  
  scale_fill_manual(values = c(
    "Photoperiod" = "#1b9e77",
    "Temperature stability" = "#7570b3",
    "Temperature autocorrelation" = "#E6AB02",
    "Mean temperature" = "#d95f02"
  )) +
  
  coord_cartesian(ylim = c(-4.5, 9)) +
  
  theme_classic(base_family = "Garamond", base_size = 14) +
  
  labs(
    x = "Environmental gradient",
    y = "Phenological plasticity\n(slope of onset vs temperature anomaly)",
    color = "Environmental gradient",
    fill = "Environmental gradient"
  )
multivoltine_offset_plasticity_interactions_plot


ggsave(
  filename = here("output", "figures", "multivoltine_offset_plasticity_interactions_plot.png"),
  plot = multivoltine_offset_plasticity_interactions_plot,
  width = 7,
  height = 5,
  dpi = 300
)
