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
    temp_between = mean(avg_temp, na.rm = TRUE),
    temp_sd      = sd(avg_temp, na.rm = TRUE),
    temp_within  = avg_temp - temp_between,
    n_years      = sum(!is.na(avg_temp)),
    
    temp_acf1 = if (first(n_years) >= 10) {
      cor(avg_temp[-n()], avg_temp[-1], use = "complete.obs")
    } else {
      NA_real_
    }
  ) |>
  ungroup() |>
  mutate(
    temp_between_sc = scale(temp_between)[,1],
    temp_sd_sc      = scale(temp_sd)[,1],
    temp_within_sc  = scale(temp_within)[,1],
    temp_acf1_sc    = scale(temp_acf1)[,1]
  )

str(clim_vars)

#### Inspect climate variables ####

clim_site <- clim_vars |>
  distinct(transect_id,
           temp_between,
           temp_sd,
           temp_acf1,
           temp_between_sc,
           temp_sd_sc,
           temp_acf1_sc)

par(mfrow = c(1,3))

hist(clim_site$temp_between,
     main = "Mean annual temp",
     xlab = "temp_between")

hist(clim_site$temp_sd,
     main = "Interannual SD",
     xlab = "temp_sd")

hist(clim_site$temp_acf1,
     main = "Lag-1 autocorrelation",
     xlab = "temp_acf1")
# ---

clim_site <- clim_vars |>
  distinct(transect_id,
           temp_between,
           temp_sd,
           temp_acf1)

round(
  cor(clim_site[, c("temp_between",
                    "temp_sd",
                    "temp_acf1")],
      use = "complete.obs"),
  2
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

m_explain <- lmer(
  ONSET_mean ~ 
    
  # ---
  # Fixed effects
  # ---
    
  # Plastic response to interannual temperature anomaly
  temp_within_sc *
    
    # Long-term mean temperature (spatial climate)
    temp_between_sc +
    
    # Plasticity moderated by interannual temperature variability
    temp_within_sc * temp_sd_sc +
    
    # Plasticity moderated by climatic predictability (lag-1 autocorrelation)
    temp_within_sc * temp_acf1_sc +
    
    # Plasticity moderated by latitude (photoperiod proxy)
    temp_within_sc * latitude_sc +
    
    
    # ---
  # Random effects
  # ---
  
  # Species-specific plasticity (reaction norm slope differs among species)
  (0 + temp_within_sc | SPECIES) +
    
    # Species-specific baseline onset timing
    (1 | SPECIES) +
    
    # Species-site specific plasticity 
    # (local populations may differ in slope)
    (0 + temp_within_sc | SPECIES:SITE_ID) +
    
    # Species-site specific baseline onset
    # (local populations differ in mean phenology)
    (1 | SPECIES:SITE_ID),
  
  
  data = df,
  REML = FALSE,
  
  # Increase optimizer iterations for complex random structure
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 1e6)
  )
)
summary(m_explain)


# Plasticity is reduced in warmer climates, with populations in colder regions showing stronger advances in onset during warm years.
# Plasticity is reduced in more climatically variable (noisy) environments, suggesting that high interannual variability weakens temperature sensitivity.
# Plasticity is stronger in more climatically predictable environments, indicating that reliable year-to-year temperature structure enhances phenological responsiveness.
# Plasticity tends to be weaker at higher latitudes, consistent with a greater role of photoperiodic constraints limiting temperature-driven shifts (though this effect is relatively weak).

#### Check consistency among phenovars ####

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

models_pheno <- df_long %>%
  split(.$phenovar) %>%
  map(~ lmer(
    pheno_value ~ 
      temp_within_sc * temp_between_sc +
      temp_within_sc * temp_sd_sc +
      temp_within_sc * temp_acf1_sc +
      temp_within_sc * latitude_sc +
      (0 + temp_within_sc | SPECIES) +
      (0 + temp_within_sc | SPECIES:SITE_ID) +
      (1 | SPECIES) +
      (1 | SPECIES:SITE_ID),
    data = .x,  
    REML = FALSE,
    control = lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 1e6)
    )
  ))

results <- map_df(
  models_pheno,
  ~ broom.mixed::tidy(.x, effects = "fixed"),
  .id = "phenovar"
)

print(results, n = Inf)

#### Check multicollinearity ####

m_vif <- lmer(
  ONSET_mean ~ temp_within_sc * temp_between_sc + 
    temp_within_sc * temp_sd_sc +
    temp_within_sc * temp_acf1_sc + 
    temp_within_sc * latitude_sc + 
    (0 + temp_within | SPECIES) +
    (0 + temp_within | SPECIES:SITE_ID) +
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
                  temp_within_sc = "Temp. anomaly (plasticity)",
                  temp_between_sc = "Mean temperature",
                  temp_sd_sc = "Temp. variability",
                  temp_acf1_sc = "Temp. predictability",
                  latitude_sc = "Latitude",
                  `temp_within_sc:temp_between_sc` = "Plasticity × Mean climate",
                  `temp_within_sc:temp_sd_sc` = "Plasticity × Variability",
                  `temp_within_sc:temp_acf1_sc` = "Plasticity × Predictability",
                  `temp_within_sc:latitude_sc` = "Plasticity × Latitude"
    ),
    
    # ---- Order predictors (facet order) ----
    term = factor(
      term,
      levels = c(
        "Temp. anomaly (plasticity)",
        "Mean temperature",
        "Temp. variability",
        "Temp. predictability",
        "Latitude",
        "Plasticity × Mean climate",
        "Plasticity × Variability",
        "Plasticity × Predictability",
        "Plasticity × Latitude"
      )
    )
  )


ggplot(plot_df,
       aes(x = phenovar,
           y = estimate,
           ymin = lower,
           ymax = upper)) +
  geom_pointrange(size = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ term, scales = "free_y") +
  coord_flip() +
  theme_bw() +
  theme(
    strip.text = element_text(size = 9),
    axis.text.y = element_text(size = 8)
  )


#### Plot interactions effects ####

mods <- c("temp_between_sc",
          "temp_sd_sc",
          "temp_acf1_sc",
          "latitude_sc")

pretty_mods <- c(
  temp_between_sc = "Mean clim.",
  temp_sd_sc      = "Clim. var.",
  temp_acf1_sc    = "Clim. pred.",
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
    temp_within_sc = seq(-2, 3.5, length.out = 150),
    moderator_value = mod_seq
  )
  
  # initialize all moderators at 0
  grid[mods] <- 0
  
  # assign focal moderator
  grid[[moderator]] <- grid$moderator_value
  
  # prediction (only relevant terms)
  grid$pred <-
    beta["(Intercept)"] +
    beta["temp_within_sc"] * grid$temp_within_sc +
    beta[moderator] * grid[[moderator]] +
    beta[paste0("temp_within_sc:", moderator)] *
    grid$temp_within_sc * grid[[moderator]]
  
  grid$phenovar <- name
  grid$moderator_name <- moderator
  
  grid
}

grid_all <- map_dfr(
  mods,
  function(mod) {
    imap_dfr(models_pheno,
             ~ make_pred_interaction(.x, .y, mod))
  }
)

library(dplyr)

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
    aes(x = temp_within_sc,
        y = pred,
        group = moderator_value,
        colour = moderator_value)
  ) +
    geom_line(alpha = 0.6) +
    scale_colour_viridis_c(
      option = "D",
      name = pretty_mods[[moderator_label]]
    ) +
    facet_wrap(~ phenovar, scales = "free_y") +
    theme_bw() +
    labs(
      x = "Temperature anomaly",
      y = "Phenology"
    )
}

plot_interaction(grid_all, "temp_sd_sc")
plot_interaction(grid_all, "temp_between_sc")
plot_interaction(grid_all, "temp_acf1_sc")
plot_interaction(grid_all, "latitude_sc")

