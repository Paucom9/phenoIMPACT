# ============================================================================================ #
# pheno_plasticity.R
#
# Author: Pau Colom
# Date: 2026-02-20
#
# Desciription: 
#
# ============================================================================================ #

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
coords_sf <- ebms_coord_df |>
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
  mutate(latitude_sc = scale(latitude)[,1])

df <- phenology_estimates |>
  left_join(
    clim_vars,
    by = c("SITE_ID", "YEAR")
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
  filename = here::here("output", "figures", "distribution_clim_anomalies.png"),
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
           photo_90,
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
  photo_90 = "Photoperiod\n(daily hours)",
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
  filename = here::here("output", "figures", "climate_histograms_26_3_26.png"),
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
         photo_90,
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
  filename = here::here("output", "figures", "corr_site_vars_26_3_26.png"),
  plot = corr_clim_pred_sds,
  width = 6,
  height = 5,
  dpi = 300
)


#### Models phenotypic plasticity ####

#### Model selection ####

model_set_plasticity <- function(d) {
  
  d <- d %>%
    dplyr::filter(
      !is.na(ONSET_mean),
      !is.na(clim_anomaly_temp_90_sc),
      !is.na(clim_background_temp_90_sc),
      !is.na(clim_predictability_temp_90_sc),
      !is.na(photo_90_sc)
    )
  
  ctrl <- lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e6)
  )
  
  f_base <- ONSET_mean ~ clim_anomaly_temp_90_sc +
    (1 + clim_anomaly_temp_90_sc || SPECIES) +
    (1 | SITE_ID)
  
  list(
    
    # H1
    m1 = lmer(f_base, data = d, REML = FALSE, control = ctrl),
    
    # H2
    m2 = lmer(update(f_base,
                     . ~ . + clim_background_temp_90_sc +
                       clim_anomaly_temp_90_sc:clim_background_temp_90_sc),
              data = d, REML = FALSE, control = ctrl),
    
    # H3
    m3 = lmer(update(f_base,
                     . ~ . + clim_predictability_temp_90_sc +
                       clim_anomaly_temp_90_sc:clim_predictability_temp_90_sc),
              data = d, REML = FALSE, control = ctrl),
    
    # H4
    m4 = lmer(update(f_base,
                     . ~ . + photo_90_sc +
                       clim_anomaly_temp_90_sc:photo_90_sc),
              data = d, REML = FALSE, control = ctrl),
    
    # H5
    m5 = lmer(update(f_base,
                     . ~ . +
                       clim_background_temp_90_sc +
                       clim_predictability_temp_90_sc +
                       clim_anomaly_temp_90_sc:clim_background_temp_90_sc +
                       clim_anomaly_temp_90_sc:clim_predictability_temp_90_sc),
              data = d, REML = FALSE, control = ctrl),
    
    # H6
    m6 = lmer(update(f_base,
                     . ~ . +
                       clim_background_temp_90_sc +
                       photo_90_sc +
                       clim_anomaly_temp_90_sc:clim_background_temp_90_sc +
                       clim_anomaly_temp_90_sc:photo_90_sc),
              data = d, REML = FALSE, control = ctrl),
    
    # H7
    m7 = lmer(update(f_base,
                     . ~ . +
                       clim_predictability_temp_90_sc +
                       photo_90_sc +
                       clim_anomaly_temp_90_sc:clim_predictability_temp_90_sc +
                       clim_anomaly_temp_90_sc:photo_90_sc),
              data = d, REML = FALSE, control = ctrl),
    
    # H8
    m8 = lmer(update(f_base,
                     . ~ . +
                       clim_background_temp_90_sc +
                       clim_predictability_temp_90_sc +
                       photo_90_sc +
                       clim_anomaly_temp_90_sc:clim_background_temp_90_sc +
                       clim_anomaly_temp_90_sc:clim_predictability_temp_90_sc +
                       clim_anomaly_temp_90_sc:photo_90_sc),
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
    Plasticity = estimate_clim_anomaly_temp_90_sc,
    Background = estimate_clim_background_temp_90_sc,
    Predictability = estimate_clim_predictability_temp_90_sc,
    Photoperiod = estimate_photo_90_sc,
    
    `Plasticity×Background` =
      `estimate_clim_anomaly_temp_90_sc:clim_background_temp_90_sc`,
    
    `Plasticity×Predictability` =
      `estimate_clim_anomaly_temp_90_sc:clim_predictability_temp_90_sc`,
    
    `Plasticity×Photoperiod` =
      `estimate_clim_anomaly_temp_90_sc:photo_90_sc`
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
      stars(p.value_clim_anomaly_temp_90_sc)
    ),
    
    Background = paste0(
      round(Background, 2),
      stars(p.value_clim_background_temp_90_sc)
    ),
    
    Predictability = paste0(
      round(Predictability, 2),
      stars(p.value_clim_predictability_temp_90_sc)
    ),
    
    Photoperiod = paste0(
      round(Photoperiod, 2),
      stars(p.value_photo_90_sc)
    ),
    
    `Plasticity×Background` = paste0(
      round(`Plasticity×Background`, 2),
      stars(`p.value_clim_anomaly_temp_90_sc:clim_background_temp_90_sc`)
    ),
    
    `Plasticity×Predictability` = paste0(
      round(`Plasticity×Predictability`, 2),
      stars(`p.value_clim_anomaly_temp_90_sc:clim_predictability_temp_90_sc`)
    ),
    
    `Plasticity×Photoperiod` = paste0(
      round(`Plasticity×Photoperiod`, 2),
      stars(`p.value_clim_anomaly_temp_90_sc:photo_90_sc`)
    )
  ) %>%
  dplyr::arrange(delta) %>%
  dplyr::mutate(
    across(-c(model, delta, weight),
           ~ ifelse(. %in% c("NA", NA), "—", .))
  )

table_clean





# 2. Forest plot function

make_forest_plot <- function(results, label){
  
  plot_df <- results %>%
    dplyr::filter(effect == "fixed",
                  term != "(Intercept)") %>%
    dplyr::mutate(
      lower = estimate - 1.96 * std.error,
      upper = estimate + 1.96 * std.error,
      
      # ---- Rename phenological variables ----
      phenovar = dplyr::recode(
        phenovar,
        ONSET_mean = "On.(mean)",
        ONSET_var = "On.(var.)",
        PEAKDAY = "PD",
        OFFSET_mean = "Off.(mean)",
        OFFSET_var = "Off.(var.)",
        FLIGHT_LENGTH_mean = "FL(mean)",
        FLIGHT_LENGTH_var = "FL(var.)"
      ),
      
      # ---- Order phenological variables ----
      phenovar = factor(
        phenovar,
        levels = rev(c(
          "On.(mean)",
          "On.(var.)",
          "PD",
          "Off.(mean)",
          "Off.(var.)",
          "FL(mean)",
          "FL(var.)"
        ))
      ),
      
      # ---- Rename predictors ----
      term = dplyr::recode(
        term,
        clim_anomaly_sc = "Clim. Anomaly (sensitivity)",
        clim_background_sc = "Clim. Background",
        clim_pred_sd_sc = "Clim. Predictability (1/SD)",
        clim_pred_lag_sc = "Clim. Predictability (lag-1)",
        latitude_sc = "Latitude",
        `clim_anomaly_sc:clim_background_sc` = "Sensitivity × Clim. Back.",
        `clim_anomaly_sc:clim_pred_sd_sc` = "Sensitivity × Clim. 1/SD",
        `clim_anomaly_sc:clim_pred_lag_sc` = "Sensitivity × Clim. lag-1",
        `clim_anomaly_sc:latitude_sc` = "Sensitivity × Latitude"
      ),
      
      # ---- Order predictors ----
      term = factor(
        term,
        levels = c(
          "Clim. Anomaly (sensitivity)",
          "Clim. Background",
          "Clim. Predictability (1/SD)",
          "Clim. Predictability (lag-1)",
          "Latitude",
          "Sensitivity × Clim. Back.",
          "Sensitivity × Clim. 1/SD",
          "Sensitivity × Clim. lag-1",
          "Sensitivity × Latitude"
        )
      )
    )
  
  p <- ggplot(
    plot_df,
    aes(x = phenovar,
        y = estimate,
        ymin = lower,
        ymax = upper)
  ) +
    geom_pointrange(size = 0.4) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_wrap(~ term, scales = "free_x") +
    coord_flip() +
    theme_bw() +
    theme(
      strip.text = element_text(size = 9),
      axis.text.y = element_text(size = 8)
    ) +
    labs(
      x = "Phenological metric",
      y = "Effect size (estimate ± 95% CI)"
    )
  
  ggsave(
    here::here("output", "figures",
               paste0("forest_", label, ".png")),
    p,
    width = 8,
    height = 5,
    dpi = 300
  )
  
  return(p)
}

# 3. Interaction settings

mods <- c("clim_background_sc",
          "clim_pred_sd_sc",
          "clim_pred_lag_sc",
          "latitude_sc")

pretty_mods <- c(
  clim_background_sc = "Climate\nbackground",
  clim_pred_sd_sc    = "Climate\npredictability\n(SD)",
  clim_pred_lag_sc   = "Climate\npredictability\n(lag-1)",
  latitude_sc        = "Latitude"
)

moderator_ranges <- setNames(
  lapply(mods, function(v)
    seq(min(df[[v]], na.rm = TRUE),
        max(df[[v]], na.rm = TRUE),
        length.out = 150)
  ),
  mods
)


# 4. Prediction grid for interaction effects

# Moderators to test
mods <- c(
  "clim_background_sc",
  "clim_pred_sd_sc",
  "clim_pred_lag_sc",
  "latitude_sc"
)

# Labels for colour legend
pretty_mods <- c(
  clim_background_sc = "Climate\nbackground",
  clim_pred_sd_sc    = "Climate\npredictability\n(SD)",
  clim_pred_lag_sc   = "Climate\npredictability\n(lag-1)",
  latitude_sc        = "Latitude"
)


# Create prediction grid for a given model and moderator
make_pred_interaction <- function(model, name, moderator, moderator_ranges){
  
  beta <- lme4::fixef(model)
  mod_seq <- moderator_ranges[[moderator]]
  
  grid <- expand.grid(
    clim_anomaly_sc = seq(-2, 3.5, length.out = 150),
    moderator_value = mod_seq
  )
  
  # set all moderators to 0
  grid[mods] <- 0
  
  # assign focal moderator
  grid[[moderator]] <- grid$moderator_value
  
  # predicted response
  grid$pred <-
    beta["(Intercept)"] +
    beta["clim_anomaly_sc"] * grid$clim_anomaly_sc +
    beta[moderator] * grid[[moderator]] +
    beta[paste0("clim_anomaly_sc:", moderator)] *
    grid$clim_anomaly_sc * grid[[moderator]]
  
  grid$phenovar <- name
  grid$moderator_name <- moderator
  
  return(grid)
}

# 5. Interaction plots

plot_interaction <- function(data, moderator_label){
  
  ggplot(
    dplyr::filter(data, moderator_name == moderator_label),
    aes(
      x = clim_anomaly_sc,
      y = pred,
      group = moderator_value,
      colour = moderator_value
    )
  ) +
    geom_line(alpha = 0.6) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_colour_viridis_c(
      option = "D",
      name = pretty_mods[[moderator_label]]
    ) +
    facet_wrap(~ phenovar, scales = "free_y") +
    theme_bw() +
    labs(
      x = "Temperature anomaly",
      y = "Day of the year or number of days"
    )
}

make_interaction_plots <- function(models, label, data){
  
  moderator_ranges <- setNames(
    lapply(mods, function(v)
      seq(
        min(data[[v]], na.rm = TRUE),
        max(data[[v]], na.rm = TRUE),
        length.out = 150
      )
    ),
    mods
  )
  
  grid_all <- purrr::map_dfr(
    mods,
    function(mod){
      purrr::imap_dfr(
        models,
        ~ make_pred_interaction(.x, .y, mod, moderator_ranges)
      )
    }
  )
  
  grid_all <- grid_all %>%
    dplyr::mutate(
      phenovar = dplyr::recode(
        phenovar,
        ONSET_mean = "Onset (mean)",
        ONSET_var = "Onset (variance)",
        PEAKDAY = "Peak day",
        OFFSET_mean = "Offset (mean)",
        OFFSET_var = "Offset (variance)",
        FLIGHT_LENGTH_mean = "Flight length (mean)",
        FLIGHT_LENGTH_var = "Flight length (variance)"
      ),
      phenovar = factor(
        phenovar,
        levels = c(
          "Onset (mean)",
          "Onset (variance)",
          "Peak day",
          "Offset (mean)",
          "Offset (variance)",
          "Flight length (mean)",
          "Flight length (variance)"
        )
      )
    )
  
  plots <- list()
  
  for(mod in mods){
    
    p <- plot_interaction(grid_all, mod)
    
    file_path <- here::here(
      "output","figures",
      paste0(mod,"_",label,".png")
    )
    
    print(file_path)
    
    ggsave(
      filename = file_path,
      plot = p,
      width = 6,
      height = 5,
      dpi = 300
    )
    
    plots[[mod]] <- p
  }
  
  return(plots)
}


# 6. Add voltinism trait

voltinism_spp <- species_traits %>%
  dplyr::transmute(
    Taxon = gsub("_", " ", Taxon),
    univoltine = Vol_max <= 1.5,
    multivoltine = Vol_max >= 2,
    strict_multivoltine = Vol_min >= 2
  )

df_long <- df_long %>%
  mutate(SPECIES = as.character(SPECIES)) %>%
  left_join(voltinism_spp, by = c("SPECIES" = "Taxon"))

# 7. Run functions

# GLOBAL
fit_global <- fit_pheno_models(df_long)
str(df_long)

make_forest_plot(fit_global$results, "global")

make_interaction_plots(
  fit_global$models,
  "global",
  df_long
)

# UNIVOLTINE
df_uni <- df_long %>%
  filter(univoltine)
str(df_uni)
length(unique(df_long$SPECIES))
length(unique(df_long$SITE_ID))

fit_uni <- fit_pheno_models(df_uni)

make_forest_plot(fit_uni$results, "univoltine")

make_interaction_plots(
  fit_uni$models,
  "univoltine",
  df_uni
)

# MULTIVOLTINE
df_multi <- df_long %>%
  filter(multivoltine)
str(df_multi)
length(unique(df_multi$SPECIES))
length(unique(df_multi$SITE_ID))

fit_multi <- fit_pheno_models(df_multi)

make_forest_plot(fit_multi$results, "multivoltine")

make_interaction_plots(
  fit_multi$models,
  "multivoltine",
  df_multi
)

# Strict MULTIVOLTINE
df_s_multi <- df_long %>%
  filter(strict_multivoltine)
str(df_s_multi)
unique(df_s_multi$SPECIES)
length(unique(df_s_multi$SPECIES))
length(unique(df_s_multi$SITE_ID))

fit_s_multi <- fit_pheno_models(df_s_multi)

make_forest_plot(fit_s_multi$results, "strict_multivoltine")

make_interaction_plots(
  fit_s_multi$models,
  "strict_multivoltine",
  df_s_multi
)

