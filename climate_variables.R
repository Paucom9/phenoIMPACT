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
library(sf)

#### Data Import and Preparation ####
# ---
here::here() # Check the current working directory
# ---
phenology_estimates  <- fread(here::here("output", "pheno_estimates_allspp.csv"))
daily_temperature <- fread(here::here("output", "climate", "daily_temperature.csv"))
ebms_coord_df <- fread(here::here("data", "ebms_transect_coord.csv"))

coords_sf <- ebms_coord_df |> filter(!is.na(transect_lon), !is.na(transect_lat)) |> distinct(transect_id, transect_lon, transect_lat) |> st_as_sf( coords = c("transect_lon", "transect_lat"), crs = 3035 )
coords_wgs84 <- st_transform(coords_sf, 4326)
coord_site <- coords_wgs84 |> mutate(latitude = st_coordinates(coords_wgs84)[,2]) |> st_drop_geometry() |> select(transect_id, latitude)
coords_dt <- data.table(coord_site)

# ---

#### Identify relevant time windows based on mean onset ####

# --- Ensure naming ---
daily_temperature[, SITE_ID := transect_id]

# --- Keep only sites with phenology ---
valid_sites <- unique(phenology_estimates$SITE_ID)

daily_temperature <- daily_temperature[
  SITE_ID %in% valid_sites
]

# --- Mean onset per site ---
onset_mean_sp <- phenology_estimates %>%
  group_by(SPECIES, SITE_ID) %>%
  summarise(onset_mean = mean(ONSET_mean, na.rm = TRUE),
            .groups = "drop") %>%
  as.data.table()

onset_mean_sp <- onset_mean_sp[
  SITE_ID %in% daily_temperature$SITE_ID
]

# --- Create dt (BASE DATASET) ---
dt <- daily_temperature[
  onset_mean_sp,
  on = .(SITE_ID),
  allow.cartesian = TRUE
]

dt[, rel_day := julian_day - onset_mean]

dt <- dt[rel_day <= 0 & rel_day > -90]

# --- ADD COORDS  ---
dt <- merge(dt, coords_dt, by = ("SITE_ID" = "transect_id"), all.x = TRUE)


# --- Photoperiod ---
daylength <- function(lat, doy) {
  lat_rad <- lat * pi / 180
  decl <- 23.44 * pi / 180 * sin(2 * pi * (doy - 80) / 365)
  ha <- acos(-tan(lat_rad) * tan(decl))
  24 * ha / pi
}

dt[, photoperiod := daylength(LATITUDE, julian_day)]

# --- Temperature windows ---
temp_windows <- dt[
  rel_day <= 0 & rel_day > -90,
  .(
    temp_30 = mean(temp[rel_day > -30], na.rm = TRUE),
    temp_60 = mean(temp[rel_day > -60], na.rm = TRUE),
    temp_90 = mean(temp, na.rm = TRUE)
  ),
  by = .(SITE_ID, SPECIES, year)
]

# --- Photoperiod windows ---
photo_windows <- dt[
  rel_day <= 0 & rel_day > -90,
  .(
    photo_30 = mean(photoperiod[rel_day > -30], na.rm = TRUE),
    photo_60 = mean(photoperiod[rel_day > -60], na.rm = TRUE),
    photo_90 = mean(photoperiod, na.rm = TRUE)
  ),
  by = .(SITE_ID, SPECIES, year)
]

setnames(temp_windows, "year", "YEAR")
setnames(photo_windows, "year", "YEAR")




#### calculate temperature-based variables ####

# Predictability

dt_pred <- copy(dt)

# --- assign windows ---
dt_pred[, window :=
          fifelse(rel_day > -30 & rel_day <= 0, "temp_30",
                  fifelse(rel_day > -60 & rel_day <= 0, "temp_60",
                          fifelse(rel_day > -90 & rel_day <= 0, "temp_90", NA_character_)))
]

dt_pred <- dt_pred[!is.na(window)]

# --- climatology (species-specific) ---
mean_profile <- dt_pred[
  , .(temp_mean = mean(temp, na.rm = TRUE)),
  by = .(SITE_ID, SPECIES, window, julian_day)
]

# --- join + residuals ---
dt_pred <- mean_profile[
  dt_pred,
  on = .(SITE_ID, SPECIES, window, julian_day)
]

dt_pred[, residual := temp - temp_mean]

# --- predictability ---
pred_window <- dt_pred[
  , {
    n_years <- uniqueN(year)
    
    list(
      clim_predictability = if (n_years >= 10) {
        -var(residual, na.rm = TRUE)
      } else NA_real_
    )
  },
  by = .(SITE_ID, SPECIES, window)
]

#---


# Climate metrics per site × window

temp_long <- melt(
  temp_windows,
  id.vars = c("SITE_ID", "SPECIES", "YEAR"),
  measure.vars = c("temp_30", "temp_60", "temp_90"),
  variable.name = "window",
  value.name = "temp"
)

clim_site_window <- temp_long[
  , {
    
    valid <- !is.na(temp)
    n_years <- uniqueN(YEAR[valid])
    
    if (n_years >= 10 && sum(valid) > 0) {
      
      list(
        clim_background = mean(temp, na.rm = TRUE),
        clim_trend = coef(lm(temp ~ YEAR))[["YEAR"]] * 10,
        clim_stability = -sd(temp, na.rm = TRUE),
        clim_autocorr = acf(temp, plot = FALSE, lag.max = 1,
                            na.action = na.pass)$acf[2]
      )
      
    } else {
      
      list(
        clim_background = NA_real_,
        clim_trend = NA_real_,
        clim_stability = NA_real_,
        clim_autocorr = NA_real_
      )
    }
    
  },
  by = .(SITE_ID, SPECIES, window)
]

#---

# ADD predictability ----

clim_site_window <- pred_window[
  clim_site_window,
  on = .(SITE_ID, SPECIES, window)
]

#---

# Merge back + anomalies ----

temp_long <- clim_site_window[
  temp_long,
  on = .(SITE_ID, SPECIES, window)
]

temp_long[, clim_anomaly := temp - clim_background]

#---

# Wide format ----

clim_vars <- dcast(
  temp_long,
  SITE_ID + SPECIES + YEAR ~ window,
  value.var = c("clim_anomaly",
                "clim_background",
                "clim_trend",
                "clim_stability",
                "clim_autocorr",
                "clim_predictability")
)

clim_vars <- merge(
  clim_vars,
  photo_windows,
  by = c("SITE_ID", "SPECIES", "YEAR"),
  all.x = TRUE
)

#---

# Scaling ----

cols_to_scale <- setdiff(names(clim_vars), c("SITE_ID", "SPECIES", "YEAR"))

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



####################################################################################

#--- Mean offset per site ---
  offset_mean_sp <- phenology_estimates %>%
  group_by(SPECIES, SITE_ID) %>%
  summarise(offset_mean = mean(OFFSET_mean, na.rm = TRUE),
            .groups = "drop") %>%
  as.data.table()

offset_mean_sp <- offset_mean_sp[
  SITE_ID %in% daily_temperature$SITE_ID
]

# --- Create dt (BASE DATASET) ---
dt <- daily_temperature[
  offset_mean_sp,
  on = .(SITE_ID),
  allow.cartesian = TRUE
]

dt[, rel_day := julian_day - offset_mean]

dt <- dt[rel_day <= 0 & rel_day > -90]

# --- ADD COORDS  ---
dt <- merge(dt, coords_dt, by = ("SITE_ID" = "transect_id"), all.x = TRUE)


# --- Photoperiod ---
daylength <- function(lat, doy) {
  lat_rad <- lat * pi / 180
  decl <- 23.44 * pi / 180 * sin(2 * pi * (doy - 80) / 365)
  ha <- acos(-tan(lat_rad) * tan(decl))
  24 * ha / pi
}

dt[, photoperiod := daylength(latitude, julian_day)]

# --- Temperature windows ---
temp_windows <- dt[
  rel_day <= 0 & rel_day > -90,
  .(
    temp_30 = mean(temp[rel_day > -30], na.rm = TRUE),
    temp_60 = mean(temp[rel_day > -60], na.rm = TRUE),
    temp_90 = mean(temp, na.rm = TRUE)
  ),
  by = .(SITE_ID, SPECIES, year)
]

# --- Photoperiod windows ---
photo_windows <- dt[
  rel_day <= 0 & rel_day > -90,
  .(
    photo_30 = mean(photoperiod[rel_day > -30], na.rm = TRUE),
    photo_60 = mean(photoperiod[rel_day > -60], na.rm = TRUE),
    photo_90 = mean(photoperiod, na.rm = TRUE)
  ),
  by = .(SITE_ID, SPECIES, year)
]

setnames(temp_windows, "year", "YEAR")
setnames(photo_windows, "year", "YEAR")




#### calculate temperature-based variables ####

# Predictability

dt_pred <- copy(dt)

# --- assign windows ---
dt_pred[, window :=
          fifelse(rel_day > -30 & rel_day <= 0, "temp_30",
                  fifelse(rel_day > -60 & rel_day <= 0, "temp_60",
                          fifelse(rel_day > -90 & rel_day <= 0, "temp_90", NA_character_)))
]

dt_pred <- dt_pred[!is.na(window)]

# --- climatology (species-specific) ---
mean_profile <- dt_pred[
  , .(temp_mean = mean(temp, na.rm = TRUE)),
  by = .(SITE_ID, SPECIES, window, julian_day)
]

# --- join + residuals ---
dt_pred <- mean_profile[
  dt_pred,
  on = .(SITE_ID, SPECIES, window, julian_day)
]

dt_pred[, residual := temp - temp_mean]

# --- predictability ---
pred_window <- dt_pred[
  , {
    n_years <- uniqueN(year)
    
    list(
      clim_predictability = if (n_years >= 10) {
        -var(residual, na.rm = TRUE)
      } else NA_real_
    )
  },
  by = .(SITE_ID, SPECIES, window)
]

#---


# Climate metrics per site × window

temp_long <- melt(
  temp_windows,
  id.vars = c("SITE_ID", "SPECIES", "YEAR"),
  measure.vars = c("temp_30", "temp_60", "temp_90"),
  variable.name = "window",
  value.name = "temp"
)

clim_site_window <- temp_long[
  , {
    
    valid <- !is.na(temp)
    n_years <- uniqueN(YEAR[valid])
    
    if (n_years >= 10 && sum(valid) > 0) {
      
      list(
        clim_background = mean(temp, na.rm = TRUE),
        clim_trend = coef(lm(temp ~ YEAR))[["YEAR"]] * 10,
        clim_stability = -sd(temp, na.rm = TRUE),
        clim_autocorr = acf(temp, plot = FALSE, lag.max = 1,
                            na.action = na.pass)$acf[2]
      )
      
    } else {
      
      list(
        clim_background = NA_real_,
        clim_trend = NA_real_,
        clim_stability = NA_real_,
        clim_autocorr = NA_real_
      )
    }
    
  },
  by = .(SITE_ID, SPECIES, window)
]

#---

# ADD predictability ----

clim_site_window <- pred_window[
  clim_site_window,
  on = .(SITE_ID, SPECIES, window)
]

#---

# Merge back + anomalies ----

temp_long <- clim_site_window[
  temp_long,
  on = .(SITE_ID, SPECIES, window)
]

temp_long[, clim_anomaly := temp - clim_background]

#---

# Wide format ----

clim_vars <- dcast(
  temp_long,
  SITE_ID + SPECIES + YEAR ~ window,
  value.var = c("clim_anomaly",
                "clim_background",
                "clim_trend",
                "clim_stability",
                "clim_autocorr",
                "clim_predictability")
)

clim_vars <- merge(
  clim_vars,
  photo_windows,
  by = c("SITE_ID", "SPECIES", "YEAR"),
  all.x = TRUE
)

#---

# Scaling ----

cols_to_scale <- setdiff(names(clim_vars), c("SITE_ID", "SPECIES", "YEAR"))

clim_vars[
  , (paste0(cols_to_scale, "_sc")) :=
    lapply(.SD, function(x) as.numeric(scale(x))),
  .SDcols = cols_to_scale
]


#### Save ####

write.csv(
  clim_vars,
  here::here("output", "climate", "climate_variables_offset.csv"),
  row.names = FALSE
)
