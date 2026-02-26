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
mean_temperature <- read.csv(here::here("output", "climate", "mean_annual_temperature.csv"), sep = ",", dec = ".")
ebms_transect_coord <- read.csv(here::here("data", "ebms_transect_coord.csv"), sep = ",", dec = ".")
str(phenology_estimates)
str(mean_temperature)
str(ebms_transect_coord)
# ---

#### Calculate within-temperature(1) and climate predictability(2) ####
#(1) i.e., temporal anomaly around each transect’s long-term mean
#(2) i.e., temporal autocorrelation

clim_vars <- mean_temperature |>
  arrange(transect_id, year) |>
  group_by(transect_id) |>
  mutate(
    clim_background = mean(avg_temp, na.rm = TRUE),
    clim_var      = sd(avg_temp, na.rm = TRUE),
    clim_anomaly  = avg_temp - clim_background,
    n_years      = sum(!is.na(avg_temp)),
    
    clim_pred = if (first(n_years) >= 10) {
      cor(avg_temp[-n()], avg_temp[-1], use = "complete.obs")
    } else {
      NA_real_
    }
  ) |>
  ungroup() |>
  mutate(
    clim_background_sc = scale(clim_background)[,1],
    clim_var_sc      = scale(clim_var)[,1],
    clim_anomaly_sc  = scale(clim_anomaly)[,1],
    clim_pred_sc    = scale(clim_pred)[,1]
  )

str(clim_vars)

#### Inspect climate variables ####

clim_site <- clim_vars |>
  distinct(transect_id,
           clim_background,
           clim_anomaly,
           clim_var,
           clim_pred,
           clim_anomaly_sc,
           clim_background_sc,
           clim_var_sc,
           clim_pred_sc)

png(
  filename = here::here("output", "figures", "climate_histograms.png"),
  width = 9,
  height = 7,
  units = "in",
  res = 300
)

par(mfrow = c(2,2))

hist(clim_site$clim_anomaly,
     main = "Year temperature anomaly (within-site)",
     xlab = "Temp. anomaly")

hist(clim_site$clim_backgroud,
     main = "Site mean annual temperature",
     xlab = "Mean temp.")

hist(clim_site$clim_var,
     main = "Site interannual variability",
     xlab = "Temp. SD")

hist(clim_site$clim_pred,
     main = "Site temperature predictability",
     xlab = "Lag-1 autocorrelation")

dev.off()



# ---

cor_mat <- clim_site |>
  select(clim_anomaly,
         clim_background,
         clim_var,
         clim_pred) |>
  cor(use = "complete.obs")

cor_df <- as.data.frame(as.table(cor_mat))

# keep only lower triangle
cor_df <- cor_df |>
  filter(as.numeric(Var1) > as.numeric(Var2))

corr_clim_vars <- ggplot(cor_df, aes(Var1, Var2, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = round(Freq, 2)), size = 4) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, limits = c(-1,1)) +
  theme_minimal() +
  labs(x = NULL, y = NULL, fill = "r") +
  coord_equal()


ggsave(
  filename = here::here("output", "figures", "corr_climate_vars.png"),
  plot = corr_clim_vars,
  width = 5,
  height = 5,
  dpi = 300
)

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

#### Model phenotypic plasticity: Do climate and latitude explain variation in plasticity? ####

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
  )

models_pheno_sensitivity <- df_long %>%
  split(.$phenovar) %>%                           # fit a separate model for each phenology variable
  map(~ lmer(
    pheno_value ~ 
      clim_anomaly_sc * clim_background_sc +      # sensitivity interacting with background climate
      clim_anomaly_sc * clim_var_sc +             # sensitivity interacting with climate variability
      clim_anomaly_sc * clim_pred_sc +            # sensitivity interacting with climate predictability
      clim_anomaly_sc * latitude_sc +             # sensitivity interacting with latitude
      (1 | SPECIES) +                             # random intercept for species
      (1 | SPECIES:SITE_ID) +                     # random intercept for population (species × site)
      (0 + clim_anomaly_sc | SPECIES) +           # random slope of sensitivity across species
      (0 + clim_anomaly_sc | SPECIES:SITE_ID),    # random slope of sensitivity across populations
    data = .x,  
    REML = FALSE,
    control = lmerControl(                        # robust optimizer for complex models
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 1e6)
    )
  ))


results <- map_df(
  models_pheno_sensitivity,
  ~ broom.mixed::tidy(.x, effects = "fixed"),
  .id = "phenovar"
)

print(results, n = Inf)


#### Check multicollinearity ####

m_vif <- lmer(
  ONSET_mean ~ clim_anomaly_sc * clim_background_sc + 
    clim_anomaly_sc * clim_var_sc +
    clim_anomaly_sc * clim_pred_sc + 
    clim_anomaly_sc * latitude_sc + 
    (0 + clim_anomaly_sc | SPECIES) +
    (0 + clim_anomaly_sc | SPECIES:SITE_ID) +
    (1 | SPECIES) + 
    (1 | SPECIES:SITE_ID), 
  data = df,
  REML = FALSE
)

performance::check_collinearity(m_vif)


#### Plot effects ####

plot_df <- results %>%
  filter(effect == "fixed",
         term != "(Intercept)") %>%
  mutate(
    lower = estimate - 1.96 * std.error,
    upper = estimate + 1.96 * std.error,
    
    # ---- Rename phenological variables ----
    phenovar = recode(phenovar,
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
    term = recode(term,
                  clim_anomaly_sc = "Clim. Anomaly (sensitivity)",
                  clim_background_sc = "Clim. Background",
                  clim_var_sc = "Clim. Variability",
                  clim_pred_sc = "Clim. Predictability",
                  latitude_sc = "Latitude",
                  `clim_anomaly_sc:clim_background_sc` = "Sensitivity × Clim. Back.",
                  `clim_anomaly_sc:clim_var_sc` = "Sensitivity × Clim. Var.",
                  `clim_anomaly_sc:clim_pred_sc` = "Sensitivity × Clim. Pred.",
                  `clim_anomaly_sc:latitude_sc` = "Sensitivity × Latitude"
    ),
    
    # ---- Order predictors (facet order) ----
    term = factor(
      term,
      levels = c(
        "Clim. Anomaly (sensitivity)",
        "Clim. Background",
        "Clim. Variability",
        "Clim. Predictability",
        "Latitude",
        "Sensitivity × Clim. Back.",
        "Sensitivity × Clim. Var.",
        "Sensitivity × Clim. Pred.",
        "Sensitivity × Latitude"
      )
    )
  )


interaction_effects <- ggplot(plot_df,
       aes(x = phenovar,
           y = estimate,
           ymin = lower,
           ymax = upper)) +
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

interaction_effects

ggsave(
  filename = here::here("output", "figures", "sensitivity_interaction_effects.png"),
  plot = interaction_effects,
  width = 8,
  height = 5,
  dpi = 300
)

#### Plot interactions effects ####

mods <- c("clim_background_sc",
          "clim_var_sc",
          "clim_pred_sc",
          "latitude_sc")

pretty_mods <- c(
  clim_background_sc = "Climate 
back",
  clim_var_sc      = "Climate 
var.",
  clim_pred_sc    = "Climate 
pred.",
  latitude_sc     = "Latitude"
)

moderator_ranges <- setNames(
  lapply(mods, function(v)
    seq(min(df[[v]], na.rm = TRUE),
        max(df[[v]], na.rm = TRUE),
        length.out = 150)
  ),
  mods
)

make_pred_interaction <- function(model, name, moderator) {
  
  beta <- fixef(model)
  mod_seq <- moderator_ranges[[moderator]]
  
  grid <- expand.grid(
    clim_anomaly_sc = seq(-2, 3.5, length.out = 150),
    moderator_value = mod_seq
  )
  
  # initialize all moderators at 0
  grid[mods] <- 0
  
  # assign focal moderator
  grid[[moderator]] <- grid$moderator_value
  
  # prediction (only relevant terms)
  grid$pred <-
    beta["(Intercept)"] +
    beta["clim_anomaly_sc"] * grid$clim_anomaly_sc +
    beta[moderator] * grid[[moderator]] +
    beta[paste0("clim_anomaly_sc:", moderator)] *
    grid$clim_anomaly_sc * grid[[moderator]]
  
  grid$phenovar <- name
  grid$moderator_name <- moderator
  
  grid
}

grid_all <- map_dfr(
  mods,
  function(mod) {
    imap_dfr(models_pheno_sensitivity,
             ~ make_pred_interaction(.x, .y, mod))
  }
)


grid_all <- grid_all %>%
  mutate(
    phenovar = recode(
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

plot_interaction <- function(data, moderator_label) {
  
  ggplot(
    dplyr::filter(data, moderator_name == moderator_label),
    aes(x = clim_anomaly_sc,
        y = pred,
        group = moderator_value,
        colour = moderator_value)
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

clim_var_fig <- plot_interaction(grid_all, "clim_var_sc")
clim_background_fig <- plot_interaction(grid_all, "clim_background_sc")
clim_pred_fig <- plot_interaction(grid_all, "clim_pred_sc")
clim_lat_fig <- plot_interaction(grid_all, "latitude_sc")

ggsave(
  filename = here::here("output", "figures", "clim_var_interactions.png"),
  plot = clim_var_fig,
  width = 6,
  height = 5,
  dpi = 300)

ggsave(
  filename = here::here("output", "figures", "clim_background_interactions.png"),
  plot = clim_background_fig,
  width = 6,
  height = 5,
  dpi = 300)

ggsave(
  filename = here::here("output", "figures", "clim_pred_interactions.png"),
  plot = clim_pred_fig,
  width = 6,
  height = 5,
  dpi = 300)

ggsave(
  filename = here::here("output", "figures", "clim_lat_interactions.png"),
  plot = clim_lat_fig,
  width = 6,
  height = 5,
  dpi = 300)
