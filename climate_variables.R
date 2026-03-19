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

clim_vars <- mean_temperature |>
  arrange(transect_id, year) |>
  group_by(transect_id) |>
  mutate(
    clim_background = mean(avg_temp, na.rm = TRUE),
    clim_anomaly    = avg_temp - clim_background,
    n_years         = sum(!is.na(avg_temp)),
    
    # Predictability based on variability (higher = more predictable)
    clim_pred_sd = -sd(avg_temp, na.rm = TRUE),
    
    # Predictability based on temporal structure (lag-1 autocorrelation)
    clim_pred_lag = if (first(n_years) >= 10) {
      acf(avg_temp, plot = FALSE, lag.max = 1,
          na.action = na.pass)$acf[2]
    } else {
      NA_real_
    }
  ) |>
  ungroup() |>
  mutate(
    clim_background_sc = scale(clim_background)[,1],
    clim_anomaly_sc    = scale(clim_anomaly)[,1],
    clim_pred_sd_sc    = scale(clim_pred_sd)[,1],
    clim_pred_lag_sc   = scale(clim_pred_lag)[,1]
  )

str(clim_vars)


write.csv(
  clim_vars,
  file.path("output", "climate", "climate_variables.csv"),
  row.names = FALSE
)
