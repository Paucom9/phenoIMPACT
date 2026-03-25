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
phenology_estimates  <- fread(here::here("output", "pheno_estimates_allspp.csv"))
daily_temperature <- fread(here::here("output", "climate", "daily_temperature.csv"))

# ---


#### Calculate mean temperatures for relevant time windows based on mean onset ####

# --- Ensure naming ---
setnames(daily_temperature, tolower(names(daily_temperature)))
daily_temperature[, SITE_ID := transect_id]

# --- Keep only sites with phenology ---
valid_sites <- unique(phenology_estimates$SITE_ID)

daily_temperature <- daily_temperature[
  SITE_ID %in% valid_sites
]

# --- Mean onset per site ---
onset_site <- phenology_estimates %>%
  group_by(SITE_ID) %>%
  summarise(
    onset_mean_site = mean(ONSET_mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  as.data.table()

# --- Merge onset into daily data ---
dt <- merge(daily_temperature, onset_site, by = "SITE_ID", all.x = FALSE)

# --- Relative day to onset (within year) ---
dt[, rel_day := julian_day - onset_mean_site]

# --- Compute windows per SITE × YEAR ---
temp_windows <- dt[
  rel_day <= 0 & rel_day > -90,
  .(
    temp_30 = mean(temp[rel_day > -30], na.rm = TRUE),
    temp_60 = mean(temp[rel_day > -60], na.rm = TRUE),
    temp_90 = mean(temp, na.rm = TRUE)
  ),
  by = .(SITE_ID, year)
]

# --- Rename ---
setnames(temp_windows, "year", "YEAR")


#### calculate temperature-based variables ####

#### --- 1. Predictability (Cuchot-style, daily data) ----

dt_pred <- copy(dt)

# --- assign windows ---
dt_pred[, window :=
          fifelse(rel_day > -30 & rel_day <= 0, "temp_30",
                  fifelse(rel_day > -60 & rel_day <= 0, "temp_60",
                          fifelse(rel_day > -90 & rel_day <= 0, "temp_90", NA_character_)))
]

dt_pred <- dt_pred[!is.na(window)]

# --- climatology (mean per DOY within window) ---
mean_profile <- dt_pred[
  , .(temp_mean = mean(temp, na.rm = TRUE)),
  by = .(SITE_ID, window, julian_day)
]

# --- residuals ---
dt_pred <- mean_profile[dt_pred, on = .(SITE_ID, window, julian_day)]
dt_pred[, residual := temp - temp_mean]

# --- predictability per site × window ---
pred_window <- dt_pred[
  , {
    n_years <- uniqueN(year)
    
    list(
      clim_predictability = if (n_years >= 10) {
        -var(residual, na.rm = TRUE)
      } else NA_real_
    )
  },
  by = .(SITE_ID, window)
]

#---

#### --- 2. Annual temperature (windows already computed) ----

temp_long <- melt(
  temp_windows,
  id.vars = c("SITE_ID", "YEAR"),
  measure.vars = c("temp_30", "temp_60", "temp_90"),
  variable.name = "window",
  value.name = "temp"
)

#---

#### --- 3. Climate metrics per site × window ----

clim_site_window <- temp_long[
  , {
    n_years <- sum(!is.na(temp))
    
    list(
      n_years = n_years,
      clim_background = mean(temp, na.rm = TRUE),
      
      clim_trend = if (n_years >= 10)
        coef(lm(temp ~ YEAR))[["YEAR"]] * 10 else NA_real_,
      
      clim_stability = if (n_years >= 10)
        -sd(temp, na.rm = TRUE) else NA_real_,
      
      clim_autocorr = if (n_years >= 10)
        acf(temp, plot = FALSE, lag.max = 1,
            na.action = na.pass)$acf[2] else NA_real_
    )
  },
  by = .(SITE_ID, window)
]

#---

#### --- 4. ADD predictability (the correct one) ----

clim_site_window <- pred_window[
  clim_site_window,
  on = .(SITE_ID, window)
]

#---

#### --- 5. Merge back + anomalies ----

temp_long <- clim_site_window[
  temp_long,
  on = .(SITE_ID, window)
]

temp_long[, clim_anomaly := temp - clim_background]

#---

#### --- 6. Wide format ----

clim_vars <- dcast(
  temp_long,
  SITE_ID + YEAR ~ window,
  value.var = c("clim_anomaly",
                "clim_background",
                "clim_trend",
                "clim_stability",
                "clim_autocorr",
                "clim_predictability")
)

#---

#### --- 7. Scaling ----

cols_to_scale <- setdiff(names(clim_vars), c("SITE_ID", "YEAR"))

clim_vars[
  , (paste0(cols_to_scale, "_sc")) :=
    lapply(.SD, function(x) as.numeric(scale(x))),
  .SDcols = cols_to_scale
]
#### Save ####

write.csv(
  clim_vars,
  here::here("output", "climate", "climate_variables.csv"),
  row.names = FALSE
)