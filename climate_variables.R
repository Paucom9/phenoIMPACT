# ============================================================================================ #
# climate_variables.R
#
# Author: Pau Colom
# Date: 2026-03-19
#
# Desciription: 
#
# ============================================================================================ #

library(dplyr)

#### Data Import and Preparation ####
# ---
here::here() # Check the current working directory
# ---
mean_temperature <- read.csv(here::here("output", "climate", "mean_annual_temperature.csv"), sep = ",", dec = ".")

str(mean_temperature)

# ---

#### XXX ####


# --- 1. Warming trend per transect ---
clim_trend_df <- mean_temperature |>
  group_by(transect_id) |>
  summarise(
    clim_trend = if (sum(!is.na(avg_temp)) >= 10) {
      coef(lm(avg_temp ~ year))[["year"]] * 10  # °C per decade
    } else {
      NA_real_
    },
    .groups = "drop"
  )

# --- 2. Main dataset ---
clim_vars <- mean_temperature |>
  arrange(transect_id, year) |>
  group_by(transect_id) |>
  mutate(
    clim_background = mean(avg_temp, na.rm = TRUE),
    clim_anomaly    = avg_temp - clim_background,
    
    # variability (better name)
    clim_var = sd(avg_temp, na.rm = TRUE),
    
    # lag-1 autocorrelation (predictability)
    n_years = sum(!is.na(avg_temp)),
    clim_pred_lag = if (first(n_years) >= 10) {
      acf(avg_temp, plot = FALSE, lag.max = 1,
          na.action = na.pass)$acf[2]
    } else {
      NA_real_
    }
  ) |>
  ungroup() |>
  
  # --- 3. add warming ---
  left_join(clim_trend_df, by = "transect_id") |>
  
  # --- 4. scale ---
  mutate(
    clim_background_sc = scale(clim_background)[,1],
    clim_anomaly_sc    = scale(clim_anomaly)[,1],
    clim_var_sc        = scale(clim_var)[,1],
    clim_pred_lag_sc   = scale(clim_pred_lag)[,1],
    clim_trend_sc      = scale(clim_trend)[,1]
  )

str(clim_vars)


write.csv(
  clim_vars,
  here::here("output", "climate", "climate_variables.csv"),
  row.names = FALSE
)
