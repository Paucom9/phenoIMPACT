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
library(here)
library(data.table)

#### Data Import and Preparation ####
# ---
here::here() # Check the current working directory
# ---
phenology_estimates  <- read.csv(here::here("output", "pheno_estimates_allspp.csv"), sep = ",", dec = ".")
daily_temperature <- read.csv(here::here("output", "climate", "daily_temperature.csv"), sep = ",", dec = ".")


# ---


#### Calculate mean temperatures for relevant time windows based on mean onset ####

# Select site–year combinations with phenology data

sites_years <- phenology_estimates %>%
  distinct(SITE_ID, YEAR) %>%            # unique site–year combinations
  mutate(transect_id = SITE_ID)          # match naming with temperature dataset


# Add previous year (important for pre-onset windows crossing years)

sites_years_ext <- sites_years %>%
  bind_rows(
    sites_years %>%
      mutate(YEAR = YEAR - 1)                 # include previous year to capture pre-onset days
  ) %>%
  distinct()


# Use annual onset (SPECIES × SITE_ID × YEAR)

onset_dt <- as.data.table(phenology_estimates)[
  , .(SPECIES, SITE_ID, YEAR, ONSET_mean)
]


# Create continuous time variable (to cross years properly)

dt_filtered[, day_global := year * 365 + julian_day]
onset_dt[, onset_global := YEAR * 365 + ONSET_mean]

# Compute temperature in windows preceding onset

temp_windows <- dt_filtered[
  onset_dt,
  on = .(SITE_ID),
  by = .EACHI,
  .(
    SPECIES = i.SPECIES,
    SITE_ID = i.SITE_ID,
    YEAR    = i.YEAR,
    
    # mean temperature in the last 30 days before onset
    temp_30 = mean(
      temp[day_global <= i.onset_global &
             day_global > (i.onset_global - 30)],
      na.rm = TRUE
    ),
    
    # mean temperature in the last 60 days before onset
    temp_60 = mean(
      temp[day_global <= i.onset_global &
             day_global > (i.onset_global - 60)],
      na.rm = TRUE
    ),
    
    # mean temperature in the last 90 days before onset
    temp_90 = mean(
      temp[day_global <= i.onset_global &
             day_global > (i.onset_global - 90)],
      na.rm = TRUE
    )
  )
]

temp_windows <- temp_windows[, !duplicated(names(temp_windows)), with = FALSE]


#### Calculate climate variables based on relevant time windows ####

temp_site_year <- temp_windows %>%
  group_by(SITE_ID, YEAR) %>%
  summarise(
    temp_30 = mean(temp_30, na.rm = TRUE),
    temp_60 = mean(temp_60, na.rm = TRUE),
    temp_90 = mean(temp_90, na.rm = TRUE),
    .groups = "drop"
  )

clim_trend_df <- temp_site_year |>
  group_by(SITE_ID) |>
  summarise(
    clim_trend = if (sum(!is.na(temp_90)) >= 10) {
      coef(lm(temp_60 ~ YEAR))[["YEAR"]] * 10   # °C per decade
    } else {
      NA_real_
    },
    .groups = "drop"
  )

clim_vars <- temp_site_year |>
  arrange(SITE_ID, YEAR) |>
  group_by(SITE_ID) |>
  mutate(
    clim_background = mean(temp_60, na.rm = TRUE),
    clim_anomaly    = temp_60 - clim_background,
    
    # interannual variability
    clim_var = sd(temp_60, na.rm = TRUE),
    
    # number of years
    n_years = sum(!is.na(temp_60)),
    
    # lag-1 autocorrelation (predictability)
    clim_pred_lag = if (first(n_years) >= 10) {
      acf(temp_60, plot = FALSE, lag.max = 1,
          na.action = na.pass)$acf[2]
    } else {
      NA_real_
    }
  ) |>
  ungroup() |>
  
  # --- add warming ---
  left_join(clim_trend_df, by = "SITE_ID") |>
  
  # --- scale ---
  mutate(
    clim_background_sc = scale(clim_background)[,1],
    clim_anomaly_sc    = scale(clim_anomaly)[,1],
    clim_var_sc        = scale(clim_var)[,1],
    clim_pred_lag_sc   = scale(clim_pred_lag)[,1],
    clim_trend_sc      = scale(clim_trend)[,1]
  )


write.csv(
  clim_vars,
  here::here("output", "climate", "climate_variables.csv"),
  row.names = FALSE
)
