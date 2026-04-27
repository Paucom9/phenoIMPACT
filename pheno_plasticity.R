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


# ---

#### Data Import and Preparation ####
# ---
here::here() # Check the current working directory
# ---
phenology_estimates  <- read.csv(here::here("output", "pheno_estimates_allspp.csv"), sep = ",", dec = ".")
clim_vars <- read.csv(here::here("output", "climate", "climate_variables.csv"), sep = ",", dec = ".")
ebms_transect_coord <- read.csv(here::here("data", "ebms_transect_coord.csv"), sep = ",", dec = ".")
species_traits <- read.csv(here::here("data", "species_trait_table.csv"), sep = ";", dec = ".")

# ---

# Merge datasets #
# Filter and convert to sf
coords_sf <- ebms_transect_coord |>
  filter(!is.na(transect_lon),
         !is.na(transect_lat)) |>
  distinct(transect_id, transect_lon, transect_lat) |>
  st_as_sf(
    coords = c("transect_lon", "transect_lat"),
    crs = 3035
  )

# Transform to WGS84
coords_wgs84 <- st_transform(coords_sf, 4326)

# Extract latitude
coord_site <- coords_wgs84 |>
  mutate(latitude = st_coordinates(coords_wgs84)[,2]) |>
  st_drop_geometry() |>
  select(transect_id, latitude)

clim_vars <- clim_vars |>
  left_join(coord_site, by = c("SITE_ID" = "transect_id")) |>
  mutate(latitude = scale(latitude)[,1])

df <- phenology_estimates |>
  left_join(
    clim_vars,
    by = c("SPECIES","SITE_ID", "YEAR")
  )

df <- df |>
  mutate(
    SPECIES = factor(SPECIES),
    SITE_ID = factor(SITE_ID)
  )

df$sp_site <- interaction(df$SPECIES, df$SITE_ID, drop = TRUE)


# ------------------------------------------------------------------------------------- #
#### Data exploration #### 

# Inspect climatic anomalies # 

med <- median(df$clim_anomaly_temp_30, na.rm = TRUE)

p <- ggplot(df, aes(x = clim_anomaly_temp_30)) +
  geom_histogram(fill = "grey", color = "white", bins = 30) +
  
  # line at zero
  geom_vline(xintercept = 0, linetype = "dashed", size = 1) +
  
  # line at median
  geom_vline(xintercept = med, color = "red", linetype = "dashed", size = 1) +
  
  # annotation
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
  filename = here::here("output", "figures", "distribution_clim_anomalies_26_3_26.png"),
  plot = p,
  width = 8,
  height = 6
)

# Inspect predictor variables #

df_cor <- df |>
  distinct(SITE_ID,
           clim_background_temp_90,
           clim_trend_temp_90,
           clim_autocorr_temp_90,
           clim_stability_temp_90,
           clim_predictability_temp_90,
           photo_winter_fixed_fixed,
           latitude)

# long format 
df_long <- df_cor |>
  pivot_longer(
    cols = -SITE_ID,
    names_to = "variable",
    values_to = "value"
  )

# labels
labels <- c(
  clim_background_temp_90   = "Background temperature\n(°C)",
  clim_trend_temp_90        = "Temperature trend\n(°C per decade)",
  clim_autocorr_temp_90     = "Autocorrelation\n(lag-1)",
  clim_stability_temp_90    = "Stability\n(-SD)",
  clim_predictability_temp_90 = "Predictability\n(-Var residuals)",
  photo_winter_fixed_fixed = "Photoperiod\n(daily hours)",
  latitude                  = "Latitude\n(degrees)"
)

df_long$variable <- factor(df_long$variable,
                           levels = names(labels),
                           labels = labels)

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
  filename = here::here("output", "figures", "climate_histograms_26_3_26B.png"),
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
#### Anomaly time-window selection #### 

mod_30 <- lmer(
  ONSET_mean ~ clim_anomaly_temp_30 +
    (1 | SITE_ID) + (1 + clim_anomaly_temp_30 | SPECIES),
  data = df, REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e6)
  )
)

mod_60 <- lmer(
  ONSET_mean ~ clim_anomaly_temp_60 +
    (1 | SITE_ID) + (1 + clim_anomaly_temp_60 | SPECIES),
  data = df, REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e6)
  )
)

mod_90 <- lmer(
  ONSET_mean ~ clim_anomaly_temp_90 +
    (1 | SITE_ID) + (1 + clim_anomaly_temp_90 | SPECIES),
  data = df, REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e6)
  )
)

AIC(mod_30, mod_60, mod_90) ### Best time-window is 60 days before onset
summary(mod_30)
summary(mod_60)
summary(mod_90)
fixef(mod_30)["clim_anomaly_temp_30"]
fixef(mod_60)["clim_anomaly_temp_60"]
fixef(mod_90)["clim_anomaly_temp_90"]


# Plot effects

library(broom.mixed)
library(dplyr)

extract_effect <- function(mod, var, label) {
  tidy(mod) %>%
    filter(term == var) %>%
    mutate(window = label)
}

df_eff <- bind_rows(
  extract_effect(mod_30, "clim_anomaly_temp_30", "30 days"),
  extract_effect(mod_60, "clim_anomaly_temp_60", "60 days"),
  extract_effect(mod_90, "clim_anomaly_temp_90", "90 days")
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
  make_pred(mod_30, "clim_anomaly_temp_30", "30 days"),
  make_pred(mod_60, "clim_anomaly_temp_60", "60 days"),
  make_pred(mod_90, "clim_anomaly_temp_90", "90 days")
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
#### photoperiod time-window selection ####

mod_photo30 <- lmer(
  ONSET_mean ~ clim_anomaly_temp_60 * photo_30 +
    (1 | SITE_ID) + (1 + clim_anomaly_temp_60 | SPECIES),
  data = df, REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e6)
  )
)

mod_photo60 <- lmer(
  ONSET_mean ~ clim_anomaly_temp_60 * photo_60 +
    (1 | SITE_ID) + (1 + clim_anomaly_temp_60 | SPECIES),
  data = df, REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e6)
  )
)

mod_photo90 <- lmer(
  ONSET_mean ~ clim_anomaly_temp_60 * photo_90 +
    (1 | SITE_ID) + (1 + clim_anomaly_temp_60 | SPECIES),
  data = df, REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e6)
  )
)

AIC(mod_photo30, mod_photo60, mod_photo90) # best time-window is 90 days pre-onset
summary(mod_photo30)
summary(mod_photo60)
summary(mod_photo90)

extract_interaction <- function(mod, var, photo, label) {
  
  term_name <- paste0(var, ":", photo)
  
  tidy(mod) %>%
    filter(term == term_name) %>%
    mutate(window = label)
}

df_eff_photo <- bind_rows(
  extract_interaction(mod_photo30, "clim_anomaly_temp_60", "photo_30", "30 days"),
  extract_interaction(mod_photo60, "clim_anomaly_temp_60", "photo_60", "60 days"),
  extract_interaction(mod_photo90, "clim_anomaly_temp_60", "photo_90", "90 days")
) %>%
  mutate(
    lower = estimate - 1.96 * std.error,
    upper = estimate + 1.96 * std.error
  )

ggplot(df_eff_photo, aes(x = window, y = estimate)) +
  
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  theme_classic(base_family = "Garamond") +
  
  labs(
    x = "Photoperiod time window",
    y = "Effect of photoperiod on plasticity\n(interaction with temperature anomaly)"
  )

photo_vars <- c("photo_30", "photo_60", "photo_90")

get_slope_df <- function(photo_var) {
  
  form <- as.formula(
    paste0("ONSET_mean ~ clim_anomaly_temp_60 * (", photo_var, ") + ",
           "(1 | SITE_ID) + (1 + clim_anomaly_temp_60 | SPECIES)")
  )
  
  mod <- lmer(form, data = df, REML = FALSE,
              control = lmerControl(optimizer = "bobyqa",
                                    optCtrl = list(maxfun = 2e6)))
  
  b <- fixef(mod)
  V <- vcov(mod)
  
  coef_names <- c(
    "clim_anomaly_temp_60",
    paste0("clim_anomaly_temp_60:", photo_var)
  )
  
  x_seq <- seq(-2.5, 2.5, length.out = 200)
  
  get_slope <- function(x) {
    
    g <- c(1, x)
    
    slope <- as.numeric(b[coef_names] %*% g)
    
    se <- sqrt(as.numeric(t(g) %*% V[coef_names, coef_names] %*% g))
    
    c(slope = slope, se = se)
  }
  
  out <- t(sapply(x_seq, get_slope))
  
  data.frame(
    photo = x_seq,
    slope = out[, "slope"],
    lower = out[, "slope"] - 1.96 * out[, "se"],
    upper = out[, "slope"] + 1.96 * out[, "se"],
    window = photo_var
  )
}


# ---------------------------------------------------------------------------- #
#### Model selection ####

d <- df %>%
  dplyr::filter(
    !is.na(ONSET_mean),
    !is.na(clim_anomaly_temp_60),
    !is.na(photo_60),
    !is.na(clim_background_temp_60),
    !is.na(clim_predictability_temp_60),
    !is.na(clim_autocorr_temp_60),
    !is.na(clim_trend_temp_60)
  )

global_mod <- lmer(
  ONSET_mean ~ clim_anomaly_temp_60 *
    (photo_60 +
       clim_background_temp_60 +
       clim_predictability_temp_60 +
       clim_autocorr_temp_60 +
       clim_trend_temp_60) +
    (1 | SITE_ID) +
    (1 + clim_anomaly_temp_60 | SPECIES),
  data = d, REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e6)
  )
)


dd <- dredge(
  global_mod,
  fixed = "clim_anomaly_temp_60",
  subset =
    dc("photo_60", "clim_anomaly_temp_60:photo_60") &
    dc("clim_background_temp_60", "clim_anomaly_temp_60:clim_background_temp_60") &
    dc("clim_predictability_temp_60", "clim_anomaly_temp_60:clim_predictability_temp_60") &
    dc("clim_autocorr_temp_60", "clim_anomaly_temp_60:clim_autocorr_temp_60") &
    dc("clim_trend_temp_60", "clim_anomaly_temp_60:clim_trend_temp_60")
)

model.sel(dd)

best_mod <- get.models(dd, subset = 1)[[1]]
summary(best_mod)
check_collinearity(best_mod)

# --- all models ---
ms <- model.sel(dd)

all_models <- as.data.frame(ms)

# --- stars function ---
stars <- function(p) {
  ifelse(p < 0.001, "***",
         ifelse(p < 0.01, "**",
                ifelse(p < 0.05, "*", "")))
}

# --- p-values from best model ---
pvals <- summary(best_mod)$coefficients[, "Pr(>|t|)"]

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

table_clean <- all_models %>%
  transmute(
    delta = round(delta, 2),
    weight = signif(weight, 2),
    
    Plasticity = fmt(clim_anomaly_temp_60, pvals, "clim_anomaly_temp_60"),
    Photoperiod = fmt(photo_60, pvals, "photo_60"),
    Background = fmt(clim_background_temp_60, pvals, "clim_background_temp_60"),
    Predictability = fmt(clim_predictability_temp_60, pvals, "clim_predictability_temp_60"),
    Autocorrelation = fmt(clim_autocorr_temp_60, pvals, "clim_autocorr_temp_60"),
    Trend = fmt(clim_trend_temp_60, pvals, "clim_trend_temp_60"),
    
    `Plast×Photo` = fmt(`clim_anomaly_temp_60:photo_60`, pvals, "clim_anomaly_temp_60:photo_60"),
    `Plast×Background` = fmt(`clim_anomaly_temp_60:clim_background_temp_60`, pvals, "clim_anomaly_temp_60:clim_background_temp_60"),
    `Plast×Predictability` = fmt(`clim_anomaly_temp_60:clim_predictability_temp_60`, pvals, "clim_anomaly_temp_60:clim_predictability_temp_60"),
    `Plast×Autocorr` = fmt(`clim_anomaly_temp_60:clim_autocorr_temp_60`, pvals, "clim_anomaly_temp_60:clim_autocorr_temp_60"),
    `Plast×Trend` = fmt(`clim_anomaly_temp_60:clim_trend_temp_60`, pvals, "clim_anomaly_temp_60:clim_trend_temp_60")
  )

table_clean


write_xlsx(
  table_clean,
  here::here("output", "onset_plasticity_table_result_v2.xlsx")
)

# ---------------------------------------------------------------------------- #
#### Plot interactions of the best model ####

mod_final <- lmer(
  ONSET_mean ~ clim_anomaly_temp_60 *
    (photo_60 +
       clim_background_temp_60 +
       clim_predictability_temp_60 + clim_autocorr_temp_60) +
    (1 | SITE_ID) +
    (1 + clim_anomaly_temp_60 | SPECIES),
  
  data = df, REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e6)
  )
)

summary(mod_final)
performance::check_collinearity(mod_final)



# sequences for predictor variables

photo_seq <- seq(-2, 2, length.out = 100)
trend_seq    <- seq(-2, 2, length.out = 100)
pred_seq  <- seq(-2, 2, length.out = 100)
auto_seq  <- seq(-2, 2, length.out = 100)
back_seq  <- seq(-2, 2, length.out = 100)


# get slopes of plasticity across gradients

get_slope_df <- function(mod, var_name, var_seq, label,
                         anomaly_name = "clim_anomaly_temp_60") {
  
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
  mod_final,
  "photo_60",
  photo_seq,
  "Photoperiod"
)


df_pred <- get_slope_df(
  mod_final,
  "clim_predictability_temp_60",
  pred_seq,
  "Temperature stability"
)

df_auto <- get_slope_df(
  mod_final,
  "clim_autocorr_temp_60",
  auto_seq,
  "Temperature autocorrelation"
)

df_back <- get_slope_df(
  mod_final,
  "clim_background_temp_60",
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
  filename = here("output", "figures", "plasticity_interactions_plot.png"),
  plot = plasticity_interactions_plot,
  width = 7,
  height = 5,
  dpi = 300
)

