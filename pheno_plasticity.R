# ============================================================================================ #
# pheno_plasticity.R
#
# Author: Pau Colom
# Date: 2026-02-20
#
# Desciription: 
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

#### Merge datasets ####
# 1. Filter and convert to sf
coords_sf <- ebms_transect_coord |>
  filter(!is.na(transect_lon),
         !is.na(transect_lat)) |>
  distinct(transect_id, transect_lon, transect_lat) |>
  st_as_sf(
    coords = c("transect_lon", "transect_lat"),
    crs = 3035
  )

# 2. Transform to WGS84
coords_wgs84 <- st_transform(coords_sf, 4326)

# 3. Extract latitude
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

str(df)


df <- df |>
  mutate(
    SPECIES = factor(SPECIES),
    SITE_ID = factor(SITE_ID)
  )

df$sp_site <- interaction(df$SPECIES, df$SITE_ID, drop = TRUE)

# Inspect climatic anomalies

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

#### Inspect predictor variables ####

df_cor <- df |>
  distinct(SITE_ID,
           clim_background_temp_90,
           clim_trend_temp_90,
           clim_autocorr_temp_90,
           clim_stability_temp_90,
           clim_predictability_temp_90,
           photo_winter_fixed_fixed,
           latitude)

# --- long format ---
df_long <- df_cor |>
  pivot_longer(
    cols = -SITE_ID,
    names_to = "variable",
    values_to = "value"
  )

# --- labels bonics ---
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

# --- plot ---
p <- ggplot(df_long, aes(x = value)) +
  geom_histogram(bins = 30, fill = "grey70", color = "black") +
  facet_wrap(~ variable, scales = "free", ncol = 3) +
  theme_classic(base_size = 12) +
  labs(x = NULL, y = "Frequency") +
  theme(
    strip.text = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
  )

# --- save ---
ggsave(
  filename = here::here("output", "figures", "climate_histograms_26_3_26B.png"),
  plot = p,
  width = 9,
  height = 7,
  dpi = 300
)

# ---

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

ggsave(
  filename = here::here("output", "figures", "corr_site_vars_26_3_26B.png"),
  plot = corr_clim_pred_sds,
  width = 6,
  height = 5,
  dpi = 300
)


#### Models: testing non-linear response of phenological sensitivity to temperature experienced pre-onset ####

# models to test the non-linear relationship between phenological sensitivity (slope of onset ~ temperature anomaly) and background temperature (mean pre-onset temperature), 
# for each combination of anomaly and background variables (using different time-windows)

anomaly_vars   <- c("clim_anomaly_temp_30", 
                    "clim_anomaly_temp_60", 
                    "clim_anomaly_temp_90")

background_vars <- c("clim_background_temp_30", 
                     "clim_background_temp_60", 
                     "clim_background_temp_90")

grid <- expand.grid(
  anomaly = anomaly_vars,
  background = background_vars,
  stringsAsFactors = FALSE
)

fit_models <- function(anom, back) {
  
  # fórmules
  f_lin <- as.formula(
    paste0("ONSET_mean ~ ", anom, " * ", back,
           " + (1 | SITE_ID) + (1 + ", anom, " | SPECIES)")
  )
  
  f_quad <- as.formula(
    paste0("ONSET_mean ~ ", anom, " * (", back, " + I(", back, "^2))",
           " + (1 | SITE_ID) + (1 + ", anom, " | SPECIES)")
  )
  
  # models
  mod_lin <- lmer(f_lin, data = df, REML = FALSE,
                  control = lmerControl(optimizer = "bobyqa",
                                        optCtrl = list(maxfun = 2e6)))
  
  mod_quad <- lmer(f_quad, data = df, REML = FALSE,
                   control = lmerControl(optimizer = "bobyqa",
                                         optCtrl = list(maxfun = 2e6)))
  
  # coeficients
  coefs <- fixef(mod_quad)
  
  beta_lin  <- coefs[paste0(anom, ":", back)]
  beta_quad <- coefs[paste0(anom, ":I(", back, "^2)")]
  
  # òptim (evitar divisió per 0)
  optimum <- ifelse(abs(beta_quad) > 1e-6,
                    -beta_lin / (2 * beta_quad),
                    NA)
  
  tibble(
    anomaly = anom,
    background = back,
    
    AIC_lin = AIC(mod_lin),
    AIC_quad = AIC(mod_quad),
    deltaAIC = AIC(mod_lin) - AIC(mod_quad),
    
    beta_linear = beta_lin,
    beta_quad = beta_quad,
    
    quad_sign = case_when(
      beta_quad > 0 ~ "U",
      beta_quad < 0 ~ "inverted_U",
      TRUE ~ "none"
    ),
    
    optimum = optimum
  )
}

results_full <- pmap_dfr(grid, ~fit_models(..1, ..2)) %>%
  arrange(AIC_quad)

results_full

#  Plotting the slope of the relationship between onset and temperature anomaly across the range of background temperatures, with confidence intervals.

get_slope_df_ci <- function(anom, back) {
  
  formula <- as.formula(
    paste0("ONSET_mean ~ ", anom, " * (", back, " + I(", back, "^2)) + ",
           "(1 | SITE_ID) + (1 + ", anom, " | SPECIES)")
  )
  
  mod <- lmer(
    formula,
    data = df,
    REML = FALSE,
    control = lmerControl(optimizer = "bobyqa",
                          optCtrl = list(maxfun = 2e6))
  )
  
  b <- fixef(mod)
  V <- vcov(mod)
  
  bg_seq <- seq(-2, 2, length.out = 200)
  
  # noms coeficients
  cn <- c(
    anom,
    paste0(anom, ":", back),
    paste0(anom, ":I(", back, "^2)")
  )
  
  # matriu X
  X <- cbind(
    1,
    bg_seq,
    bg_seq^2
  )
  
  colnames(X) <- cn
  
  # slope
  slope <- as.numeric(X %*% b[cn])
  
  # SE
  V_sub <- V[cn, cn]
  se <- sqrt(diag(X %*% V_sub %*% t(X)))
  
  tibble(
    bg = bg_seq,
    slope = slope,
    lower = slope - 1.96 * se,
    upper = slope + 1.96 * se,
    anomaly = anom,
    background = back
  )
}

slopes_ci_df <- pmap_dfr(grid, ~get_slope_df_ci(..1, ..2))

scale_color_manual(
  values = c("#E64B35", "#00A087", "#4DBBD5"),
  labels = c("30 days", "60 days", "90 days"),
  name = "Temperature anomaly\n(window)"
) +
  scale_fill_manual(
    values = c("#E64B35", "#00A087", "#4DBBD5"),
    labels = c("30 days", "60 days", "90 days"),
    name = "Temperature anomaly\n(window)"
  )

labeller = labeller(
  background = c(
    clim_background_temp_30 = "Pre-onset temperature (30 days)",
    clim_background_temp_60 = "Pre-onset temperature (60 days)",
    clim_background_temp_90 = "Pre-onset temperature (90 days)"
  )
)


ggplot(slopes_ci_df, aes(bg, slope, color = anomaly, fill = anomaly)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4) +
  geom_line(size = 1.2) +
  
  facet_wrap(
    ~ background,
    nrow = 1,
    labeller = as_labeller(c(
      clim_background_temp_30 = "30-days site mean TW",
      clim_background_temp_60 = "60-days site mean TW",
      clim_background_temp_90 = "90-days site mean TW"
    ))
  ) +
  
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  scale_color_manual(
    values = c("#2CA58D", "#F28E2B", "#6C6BD1"),
    labels = c("30 days", "60 days", "90 days"),
    name = "Temperature anomaly\n(time window)"
  ) +
  
  scale_fill_manual(
    values = c("#2CA58D", "#F28E2B", "#6C6BD1"),
    labels = c("30 days", "60 days", "90 days"),
    name = "Temperature anomaly\n(time window)"
  ) +
  
  labs(
    x = "Site mean pre-onset temperature (scaled)",
    y = "Phenological sensitivity\n(slope onset ~ temperature anomaly)"
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )



#### General Models on phenotypic plasticity ####

#### General model selection ####

model_set_plasticity <- function(d) {
  
  d <- d %>%
    dplyr::filter(
      !is.na(ONSET_mean),
      !is.na(clim_anomaly_temp_90),
      !is.na(clim_background_temp_90),
      !is.na(clim_predictability_temp_90),
      !is.na(photo_winter_fixed_fixed)
    )
  
  ctrl <- lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e6)
  )
  
  f_base <- ONSET_mean ~ clim_anomaly_temp_90 +
    (1 + clim_anomaly_temp_90 || SPECIES) +
    (1 | SITE_ID)
  
  list(
    
    # H1
    m1 = lmer(f_base, data = d, REML = FALSE, control = ctrl),
    
    # H2
    m2 = lmer(update(f_base,
                     . ~ . + clim_background_temp_90 +
                       clim_anomaly_temp_90:clim_background_temp_90),
              data = d, REML = FALSE, control = ctrl),
    
    # H3
    m3 = lmer(update(f_base,
                     . ~ . + clim_predictability_temp_90 +
                       clim_anomaly_temp_90:clim_predictability_temp_90),
              data = d, REML = FALSE, control = ctrl),
    
    # H4
    m4 = lmer(update(f_base,
                     . ~ . + photo_winter_fixed_fixed +
                       clim_anomaly_temp_90:photo_winter_fixed_fixed),
              data = d, REML = FALSE, control = ctrl),
    
    # H5
    m5 = lmer(update(f_base,
                     . ~ . +
                       clim_background_temp_90 +
                       clim_predictability_temp_90 +
                       clim_anomaly_temp_90:clim_background_temp_90 +
                       clim_anomaly_temp_90:clim_predictability_temp_90),
              data = d, REML = FALSE, control = ctrl),
    
    # H6
    m6 = lmer(update(f_base,
                     . ~ . +
                       clim_background_temp_90 +
                       photo_winter_fixed_fixed +
                       clim_anomaly_temp_90:clim_background_temp_90 +
                       clim_anomaly_temp_90:photo_winter_fixed_fixed),
              data = d, REML = FALSE, control = ctrl),
    
    # H7
    m7 = lmer(update(f_base,
                     . ~ . +
                       clim_predictability_temp_90 +
                       photo_winter_fixed_fixed +
                       clim_anomaly_temp_90:clim_predictability_temp_90 +
                       clim_anomaly_temp_90:photo_winter_fixed_fixed),
              data = d, REML = FALSE, control = ctrl),
    
    # H8
    m8 = lmer(update(f_base,
                     . ~ . +
                       clim_background_temp_90 +
                       clim_predictability_temp_90 +
                       photo_winter_fixed_fixed +
                       clim_anomaly_temp_90:clim_background_temp_90 +
                       clim_anomaly_temp_90:clim_predictability_temp_90 +
                       clim_anomaly_temp_90:photo_winter_fixed_fixed),
              data = d, REML = FALSE, control = ctrl)
  )
}

mods <- model_set_plasticity(df)

sel <- MuMIn::model.sel(mods)

sel_df <- as.data.frame(sel) %>%
  tibble::rownames_to_column("model") %>%
  dplyr::filter(delta <= 2)

sel_df

sel_full <- as.data.frame(sel) %>%
  tibble::rownames_to_column("model") %>%
  dplyr::arrange(delta)

# Check collinearity
performance::check_collinearity(mods$m8)


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
    Plasticity = estimate_clim_anomaly_temp_90,
    Background = estimate_clim_background_temp_90,
    Predictability = estimate_clim_predictability_temp_90,
    Photoperiod = estimate_photo_winter_fixed_fixed,
    
    `Plasticity×Background` =
      `estimate_clim_anomaly_temp_90:clim_background_temp_90`,
    
    `Plasticity×Predictability` =
      `estimate_clim_anomaly_temp_90:clim_predictability_temp_90`,
    
    `Plasticity×Photoperiod` =
      `estimate_clim_anomaly_temp_90:photo_winter_fixed_fixed`
  )


table_final %>%
  dplyr::select(model, delta, weight,
                Plasticity, Background, Predictability, Photoperiod,
                `Plasticity×Background`,
                `Plasticity×Predictability`,
                `Plasticity×Photoperiod`)



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
    
    Plasticity = paste0(
      round(Plasticity, 2),
      stars(p.value_clim_anomaly_temp_90)
    ),
    
    Background = paste0(
      round(Background, 2),
      stars(p.value_clim_background_temp_90)
    ),
    
    Predictability = paste0(
      round(Predictability, 2),
      stars(p.value_clim_predictability_temp_90)
    ),
    
    Photoperiod = paste0(
      round(Photoperiod, 2),
      stars(p.value_photo_winter_fixed_fixed)
    ),
    
    `Plasticity×Background` = paste0(
      round(`Plasticity×Background`, 2),
      stars(`p.value_clim_anomaly_temp_90:clim_background_temp_90`)
    ),
    
    `Plasticity×Predictability` = paste0(
      round(`Plasticity×Predictability`, 2),
      stars(`p.value_clim_anomaly_temp_90:clim_predictability_temp_90`)
    ),
    
    `Plasticity×Photoperiod` = paste0(
      round(`Plasticity×Photoperiod`, 2),
      stars(`p.value_clim_anomaly_temp_90:photo_winter_fixed_fixed`)
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
        clim_anomaly_temp_90 = "Plasticity",
        clim_background_temp_90 = "Background",
        clim_predictability_temp_90 = "Predictability",
        photo_winter_fixed_fixed = "Photoperiod",
        `clim_anomaly_temp_90:clim_background_temp_90` = "Plasticity × Background",
        `clim_anomaly_temp_90:clim_predictability_temp_90` = "Plasticity × Predictability",
        `clim_anomaly_temp_90:photo_winter_fixed_fixed` = "Plasticity × Photoperiod"
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

final_mod <- mods$m8


# Plasticity x background

pred <- ggpredict(
  final_mod,
  terms = c("clim_anomaly_temp_90", "clim_background_temp_90[-1:1]")
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
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  
  scale_color_manual(
    values = c(
      "High" = "#6C6BD1",  # blau-lila
      "Mean" = "#F28E2B",  # taronja
      "Low"  = "#2CA58D"   # verd-turquesa
    )
  ) +
  scale_fill_manual(
    values = c(
      "High" = "#6C6BD1",
      "Mean" = "#F28E2B",
      "Low"  = "#2CA58D"
    )
  ) +
  
  labs(
    x = "Temperature anomaly",
    y = "Onset (day of the year)",
    color = "Mean temperature",
    fill = "Mean temperature"
  ) +
  theme_classic(base_family = "Garamond") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

# Plasticity x photoperiod

pred <- ggpredict(
  final_mod,
  terms = c("clim_anomaly_temp_90", "photo_winter_fixed_fixed[-1:1]")
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
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  
  scale_color_manual(
    values = c(
      "High" = "#6C6BD1",  # blau-lila
      "Mean" = "#F28E2B",  # taronja
      "Low"  = "#2CA58D"   # verd-turquesa
    )
  ) +
  scale_fill_manual(
    values = c(
      "High" = "#6C6BD1",
      "Mean" = "#F28E2B",
      "Low"  = "#2CA58D"
    )
  ) +
  
  labs(
    x = "Temperature anomaly",
    y = "Onset (day of the year)",
    color = "Photoperiod",
    fill = "Photoperiod"
  ) +
  theme_classic(base_family = "Garamond") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )



# Plasticity x predictability

pred <- ggpredict(
  final_mod,
  terms = c("clim_anomaly_temp_90", "clim_predictability_temp_90[-1:1]")
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
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  
  scale_color_manual(
    values = c(
      "High" = "#6C6BD1",  # blau-lila
      "Mean" = "#F28E2B",  # taronja
      "Low"  = "#2CA58D"   # verd-turquesa
    )
  ) +
  scale_fill_manual(
    values = c(
      "High" = "#6C6BD1",
      "Mean" = "#F28E2B",
      "Low"  = "#2CA58D"
    )
  ) +
  
  labs(
    x = "Temperature anomaly",
    y = "Onset (day of the year)",
    color = "Predictability",
    fill = "Predictability"
  ) +
  theme_classic(base_family = "Garamond") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )
