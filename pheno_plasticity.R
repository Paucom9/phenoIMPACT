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
# ---

#### Data Import and Preparation ####
# ---
here::here() # Check the current working directory
# ---
phenology_estimates  <- read.csv(here::here("output", "pheno_estimates_allspp.csv"), sep = ",", dec = ".")
clim_vars <- read.csv(here::here("output", "climate", "climate_variables.csv"), sep = ",", dec = ".")
ebms_transect_coord <- read.csv(here::here("data", "ebms_transect_coord.csv"), sep = ",", dec = ".")
species_traits <- read.csv(here::here("data", "species_trait_table.csv"), sep = ";", dec = ".")

str(phenology_estimates)
str(clim_vars)
str(ebms_transect_coord)
str(species_traits)

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
  left_join(coord_site, by = "transect_id") |>
  mutate(latitude_sc = scale(latitude)[,1])

df <- phenology_estimates |>
  left_join(
    clim_vars,
    by = c("SITE_ID" = "transect_id",
           "YEAR"    = "year")
  )

str(df)


df <- df |>
  mutate(
    SPECIES = factor(SPECIES),
    SITE_ID = factor(SITE_ID)
  )

df$sp_site <- interaction(df$SPECIES, df$SITE_ID, drop = TRUE)

# Inspect climatic anomalies

png(filename = here::here("output", "figures", "distribution_clim_anomalies.png"),
 width = 800, height = 600)

med <- median(df$clim_anomaly, na.rm = TRUE)

hist(df$clim_anomaly,
     main = "",
     xlab = "Temperature anomaly (°C)",
     col = "grey",
     border = "white")

# Line at zero
abline(v = 0, col = "black", lwd = 2, lty = 2)

# Line at median
abline(v = med, col = "red", lwd = 2, lty = 2)

text(med, par("usr")[4]*0.9,
     labels = paste("Median =", round(med, 2)),
     pos = 4,
     col = "red")

dev.off()

#### Inspect predictor variables ####

df_cor <- df |>
  distinct(SITE_ID,
           clim_background,
           clim_pred_sd,
           clim_pred_lag,
           latitude)

png(
  filename = here::here("output", "figures", "climate_histograms.png"),
  width = 9,
  height = 7,
  units = "in",
  res = 300
)

par(mfrow = c(2,2))


hist(df_cor$clim_background,
     main = "Mean annual temperature",
     xlab = "Mean temp.")

hist(df_cor$clim_pred_sd,
     main = "Temperature predictability (1/SD)",
     xlab = "Inverse Temp. SD")

hist(df_cor$clim_pred_lag,
     main = "Temperature predictability (lag-1)",
     xlab = "Lag-1 autocorrelation")

hist(df_cor$latitude,
     main = "Latitude",
     xlab = "Degrees")

dev.off()



# ---

cor_mat <- df_cor |>
  select(clim_background,
         clim_pred_sd,
         clim_pred_lag,
         latitude) |>
  cor(use = "complete.obs")


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
  filename = here::here("output", "figures", "corr_site_vars.png"),
  plot = corr_clim_pred_sds,
  width = 5,
  height = 5,
  dpi = 300
)


#### Model phenotypic plasticity: Do climate and latitude explain variation in plasticity? ####

# 0. Prepara the dataset

df_long <- df |>
  pivot_longer(
    cols = c(
      ONSET_mean,
      ONSET_var,
      PEAKDAY,
      OFFSET_mean,
      OFFSET_var,
      FLIGHT_LENGTH_mean,
      FLIGHT_LENGTH_var
    ),
    names_to = "phenovar",
    values_to = "pheno_value"
  ) |>
  dplyr::mutate(
    phenovar = factor(
      phenovar,
      levels = c(
        "ONSET_mean",
        "ONSET_var",
        "PEAKDAY",
        "OFFSET_mean",
        "OFFSET_var",
        "FLIGHT_LENGTH_mean",
        "FLIGHT_LENGTH_var"
      )
    )
  )


# 1. Model function

fit_pheno_models <- function(data) {
  
  models <- data %>%
    split(.$phenovar) %>%
    purrr::map(~ lme4::lmer(
      pheno_value ~
        clim_anomaly_sc * clim_background_sc +
        clim_anomaly_sc * clim_pred_sd_sc +
        clim_anomaly_sc * clim_pred_lag_sc +
        clim_anomaly_sc * latitude_sc +
        (1 | SPECIES) +
        (1 | SPECIES:SITE_ID) +
        (0 + clim_anomaly_sc | SPECIES) +
        (0 + clim_anomaly_sc | SPECIES:SITE_ID),
      data = .x,
      REML = FALSE,
      control = lmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = 1e6)
      )
    ))
  
  results <- purrr::map_df(
    models,
    ~ broom.mixed::tidy(.x, effects = "fixed"),
    .id = "phenovar"
  )
  
  list(models = models, results = results)
}


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

