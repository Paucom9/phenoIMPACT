# ============================================================================================
# annual_mean_temperature.R
#
# Author: Pau Colom
# Date: 2026-02-19
#
# Desciription: This script extracts mean annual temperature data for each transect site 
# across four defined time periods. It uses the `climateExtract` package to retrieve 
# climate data, processes it to calculate annual averages, and compiles the results into a 
# comprehensive dataset. The final outputs include both the full annual dataset and a 
# summary of mean temperatures per site across the entire time series.
#
# ============================================================================================


#### Load required libraries ####
# ---
library(sf)
library(raster)
library(climateExtract)
library(terra)
library(magrittr)
library(dplyr)
library(data.table)
library(geodata)
library(tidyr)
library(progress)
library(tidyr)
library(lubridate)
library(here)
# ---


#### Data Import and Preparation ####
# ---
here::here() # Check the current working directory
# ---
# --- eBMS data
# Import transect coordinates
ebms_coord_df  <- read.csv(here("data", "ebms_transect_coord.csv"), sep = ",", dec = ".")
# Import country codes
country_codes  <- read.csv(here("data", "country_codes.csv"), sep = ";", dec = ".")

# --- Merge data
head(ebms_coord_df)
head(country_codes)

m_coord_c <- merge(ebms_coord_df, country_codes, by = "bms_id", all.x = TRUE)
m_coord_c<- na.omit(m_coord_c)

# --- Define the main directory for saving outputs
# Permanent output folder
output_dir <- here("output", "climate")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Temporary working folder (auto-clean when session ends)
temp_dir <- file.path(tempdir(), "climate_extract")

if (!dir.exists(temp_dir)) {
  dir.create(temp_dir, recursive = TRUE)
}

# Define ECAD periods for estracting mean temperature 
Period      <- c("1965-1979", "1980-1994", "1995-2010", "2011-2025")
First_year  <- c(1965, 1980, 1995, 2011)
Last_year   <- c(1979, 1994, 2010, 2025)

# Store everything here
combined_results <- data.frame()

#### Loop to extract mean temperatures across sites ####

# Loop across countries
for (country in unique(m_coord_c$country_code)) {
  
  country_coord <- subset(m_coord_c, country_code == country)
  
  cat("Processing country:", unique(country_coord$Country.Name), "\n")
  
  sf_point <- st_as_sf(
    country_coord[!is.na(country_coord$transect_lon), ],
    coords = c("transect_lon", "transect_lat"),
    crs = 3035
  )
  
  m_coord_sf_transformed <- st_transform(sf_point, crs = 4326)
  country_coord$lon_wgs84 <- st_coordinates(m_coord_sf_transformed)[, 1]
  country_coord$lat_wgs84 <- st_coordinates(m_coord_sf_transformed)[, 2]
  
  grid_points <- st_make_grid(sf_point, cellsize = 100000, square = TRUE, what = "polygons")
  grid_points <- grid_points[sf_point]
  bms_bbsf <- st_as_sf(st_as_sfc(st_bbox(st_transform(grid_points, crs = 4326))))
  
  country_results <- data.frame()
  
  for (i in seq_along(Period)) {
    
    Period_chunk      <- Period[i]
    First_year_chunk  <- First_year[i]
    Last_year_chunk   <- Last_year[i]
    
    cat("   Processing period:", Period_chunk, "\n")
    
    mean_output_filename <- file.path(
      temp_dir,
      paste0("raster_mean_temp_", country, "_", Period_chunk, ".tif")
    )
    
    extract_nc_value(
      first_year     = First_year_chunk,
      last_year      = Last_year_chunk,
      local_file     = FALSE,
      file_path      = temp_dir,
      sml_chunk      = Period_chunk,
      spatial_extent = bms_bbsf,
      clim_variable  = "mean temp",
      statistic      = "mean",
      grid_size      = 0.1,
      ecad_v         = NULL,
      write_raster   = TRUE,
      out            = mean_output_filename,
      return_data    = TRUE
    )
    
    rbk_mean <- terra::rast(mean_output_filename)
    
    annual_avg_temp <- temporal_aggregate(
      x             = rbk_mean,
      agg_function  = "mean",
      variable_name = "average temp",
      time_step     = "annual"
    )
    
    for (j in seq_len(nlyr(annual_avg_temp))) {
      
      year  <- names(annual_avg_temp)[j]
      coords <- country_coord[, c("lon_wgs84", "lat_wgs84")]
      
      temp_values <- terra::extract(
        annual_avg_temp[[j]],
        coords,
        method = "simple",
        ID = FALSE
      )
      
      temp_df <- data.frame(
        transect_id = country_coord$transect_id,
        year        = as.numeric(year),
        avg_temp    = temp_values[, 1],
      )
      
      country_results <- bind_rows(country_results, temp_df)
    }
  }
  
  combined_results <- bind_rows(combined_results, country_results)
  
  cat("Completed country:", unique(country_coord$Country.Name), "\n")
}

cat("Processing complete for all countries.\n")

#### Save the data ####

# View the combined data
head(combined_results)
str(combined_results)

# Save mean annual temperature
write.csv(
  combined_results,
  file.path(output_dir, "mean_annual_temperature.csv"),
  row.names = FALSE
)









