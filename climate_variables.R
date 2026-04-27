# ============================================================================================ #
# climate_variables.R
#
# Author: Pau Colom
# Date: 2026-03-19
#
# Desciription: 
#
# ============================================================================================ #

# Clean session
rm(list = ls())

# Libraries

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

# Get photoperiod data #

# --- Ensure naming ---
daily_temperature[, SITE_ID := transect_id]
daily_temperature[, YEAR := year]

# --- Keep only sites with phenology ---
valid_sites <- unique(phenology_estimates$SITE_ID)

daily_temperature <- daily_temperature[
  SITE_ID %in% valid_sites
]

# --- Add coordinates (needed for photoperiod) ---
daily_temperature <- merge(
  daily_temperature,
  coords_dt,
  by.x = "SITE_ID",
  by.y = "transect_id",
  all.x = TRUE
)

# --- Photoperiod ---
daylength <- function(lat, doy) {
  lat_rad <- lat * pi / 180
  decl <- 23.44 * pi / 180 * sin(2 * pi * (doy - 80) / 365)
  
  x <- -tan(lat_rad) * tan(decl)
  x <- pmin(pmax(x, -1), 1)
  
  ha <- acos(x)
  24 * ha / pi
}

daily_temperature[, photoperiod := daylength(latitude, julian_day)]


#### Identify relevant time windows for each phenophase and calculate relative clim vars ####

compute_climate_windows <- function(pheno_var, dt_climate, pheno_dt) {
  
  message(">>> Processing: ", pheno_var)
  
  # phenology
  pheno_mean <- pheno_dt[
    , .(pheno = mean(get(pheno_var), na.rm = TRUE)),
    by = .(SPECIES, SITE_ID)
  ]
  
  # split climate per site (CLAU)
  dt_split <- split(dt_climate, by = "SITE_ID", keep.by = TRUE)
  
  res <- vector("list", length(dt_split))
  names(res) <- names(dt_split)
  
  for (site in names(dt_split)) {
    
    dt_site <- dt_split[[site]]
    pheno_site <- pheno_mean[SITE_ID == site]
    
    if (nrow(pheno_site) == 0) next
    
    site_res <- vector("list", nrow(pheno_site))
    
    for (i in seq_len(nrow(pheno_site))) {
      
      ph <- pheno_site$pheno[i]
      sp <- pheno_site$SPECIES[i]
      
      dt_site[, rel_day := julian_day - ph]
      dt_sub <- dt_site[rel_day <= 0 & rel_day > -90]
      
      if (nrow(dt_sub) == 0) next
      
      out <- dt_sub[
        , .(
          tw30 = mean(temp[rel_day > -30], na.rm = TRUE),
          tw60 = mean(temp[rel_day > -60], na.rm = TRUE),
          tw90 = mean(temp, na.rm = TRUE),
          
          photo_tw30 = mean(photoperiod[rel_day > -30], na.rm = TRUE),
          photo_tw60 = mean(photoperiod[rel_day > -60], na.rm = TRUE),
          photo_tw90 = mean(photoperiod, na.rm = TRUE)
        ),
        by = .(YEAR = year)
      ]
      
      out[, `:=`(
        SITE_ID = site,
        SPECIES = sp,
        pheno_type = pheno_var
      )]
      
      site_res[[i]] <- out
    }
    
    res[[site]] <- rbindlist(site_res, fill = TRUE)
    
    message("   site done: ", site)
  }
  
  result <- rbindlist(res, fill = TRUE)
  
  message("<<< Done: ", pheno_var)
  
  return(result)
}

# Run the function for all phenophases #

pheno_vars <- c(
  "ONSET_mean",
  "ONSET_var",
  "PEAKDAY",
  "FIRST_PEAK",
  "OFFSET_mean",
  "OFFSET_var"
)

clim_all <- rbindlist(
  lapply(pheno_vars, compute_climate_windows,
         dt_climate = daily_temperature,
         pheno_dt = phenology_estimates),
  fill = TRUE
)

str(clim_all)


# Scale variables #

cols_to_scale <- setdiff(
  names(clim_all),
  c("SITE_ID", "SPECIES", "YEAR", "pheno_type")
)

clim_all[
  , (cols_to_scale) := lapply(.SD, function(x) as.numeric(scale(x))),
  .SDcols = cols_to_scale
]

# Calculate climate metrics for each site, species, and phenophase #

clim_vars <- clim_all[
  , `:=`(
    clim_background_tw30 = mean(tw30, na.rm = TRUE),
    clim_background_tw60 = mean(tw60, na.rm = TRUE),
    clim_background_tw90 = mean(tw90, na.rm = TRUE),
    
    clim_stability_tw30 = -sd(tw30, na.rm = TRUE),
    clim_stability_tw60 = -sd(tw60, na.rm = TRUE),
    clim_stability_tw90 = -sd(tw90, na.rm = TRUE),
    
    clim_trend_tw30 = if (sum(!is.na(tw30)) > 5)
      coef(lm(tw30 ~ YEAR))[["YEAR"]] * 10 else NA_real_,
    clim_trend_tw60 = if (sum(!is.na(tw60)) > 5)
      coef(lm(tw60 ~ YEAR))[["YEAR"]] * 10 else NA_real_,
    clim_trend_tw90 = if (sum(!is.na(tw90)) > 5)
      coef(lm(tw90 ~ YEAR))[["YEAR"]] * 10 else NA_real_,
    
    clim_autocorr_tw30 = if (sum(!is.na(tw30)) > 5)
      acf(tw30[!is.na(tw30)], plot = FALSE, lag.max = 1)$acf[2] else NA_real_,
    clim_autocorr_tw60 = if (sum(!is.na(tw60)) > 5)
      acf(tw60[!is.na(tw60)], plot = FALSE, lag.max = 1)$acf[2] else NA_real_,
    clim_autocorr_tw90 = if (sum(!is.na(tw90)) > 5)
      acf(tw90[!is.na(tw90)], plot = FALSE, lag.max = 1)$acf[2] else NA_real_
  ),
  by = .(SITE_ID, SPECIES, pheno_type)
]

# Calculate climate anomalies for each observation #

clim_vars[
  , `:=`(
    clim_anomaly_tw30 = tw30 - clim_background_tw30,
    clim_anomaly_tw60 = tw60 - clim_background_tw60,
    clim_anomaly_tw90 = tw90 - clim_background_tw90
  )
]

# Calculate predictability as the negative variance of the residuals from a climatology model #

compute_predictability <- function(pheno_var, dt_climate, pheno_dt) {
  
  pheno_mean <- pheno_dt[
    , .(pheno = mean(get(pheno_var), na.rm = TRUE)),
    by = .(SPECIES, SITE_ID)
  ]
  
  dt_split <- split(dt_climate, by = "SITE_ID", keep.by = TRUE)
  
  res <- list()
  
  for (site in names(dt_split)) {
    
    dt_site <- dt_split[[site]]
    pheno_site <- pheno_mean[SITE_ID == site]
    
    for (i in seq_len(nrow(pheno_site))) {
      
      sp <- pheno_site$SPECIES[i]
      ph <- pheno_site$pheno[i]
      
      dt_site[, rel_day := julian_day - ph]
      dt_sub <- dt_site[rel_day <= 0 & rel_day > -90]
      
      if (nrow(dt_sub) == 0) next
      
      dt_sub[, twindow :=
               fifelse(rel_day > -30, "tw30",
                       fifelse(rel_day > -60, "tw60", "tw90"))]
      
      # climatology
      clim <- dt_sub[
        , .(temp_mean = mean(temp)),
        by = .(twindow, julian_day)
      ]
      
      dt_sub <- clim[dt_sub, on = .(twindow, julian_day)]
      dt_sub[, residual := temp - temp_mean]
      
      pred <- dt_sub[
        , .(clim_predictability = -var(residual, na.rm = TRUE)),
        by = twindow
      ]
      
      pred[, `:=`(SITE_ID = site, SPECIES = sp, pheno_type = pheno_var)]
      
      res[[length(res) + 1]] <- pred
    }
  }
  
  rbindlist(res)
}

pred_all <- rbindlist(
  lapply(pheno_vars, compute_predictability,
         dt_climate = daily_temperature,
         pheno_dt = phenology_estimates)
)

pred_wide <- dcast(
  pred_all,
  SITE_ID + SPECIES + pheno_type ~ twindow,
  value.var = "clim_predictability"
)

setnames(pred_wide,
         c("tw30", "tw60", "tw90"),
         c("clim_predictability_tw30",
           "clim_predictability_tw60",
           "clim_predictability_tw90"))


# Merge predictability metrics back into the main climate dataset #

clim_vars_final <- merge(
  clim_vars,
  pred_wide,
  by = c("SITE_ID", "SPECIES", "pheno_type"),
  all.x = TRUE
)

# Scale predictability #

pred_cols <- grep("^clim_predictability", names(clim_vars_final), value = TRUE)

clim_vars_final[
  , (pred_cols) := lapply(.SD, function(x) as.numeric(scale(x))),
  .SDcols = pred_cols
]

clim_vars_final[, c("tw30", "tw60", "tw90") := NULL]

str(clim_vars_final)

# Save #

write.csv(
  clim_vars_final,
  here::here("output", "climate", "climate_variables_all_phenophases.csv"),
  row.names = FALSE
)
