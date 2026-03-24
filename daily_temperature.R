# ============================================================================================
# daily_temperature.R
#
# Author: Pau Colom
# Date: 2026-03-23
#
# Desciription: 
#
# ============================================================================================

#### Libraries ####
library(sf)
library(climateExtract)
library(terra)
library(dplyr)
library(data.table)
library(lubridate)
library(here)
library(progress)

#### Data ####
here::here()
ebms_coord_df  <- read.csv(here("data", "ebms_transect_coord.csv"))
country_codes  <- read.csv(here("data", "country_codes.csv"), sep = ";")

m_coord_c <- merge(ebms_coord_df, country_codes, by = "bms_id", all.x = TRUE) %>%
  na.omit()

output_dir <- here("output", "climate")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

temp_dir <- file.path(tempdir(), "climate_extract")
dir.create(temp_dir, showWarnings = FALSE)

Period      <- c("1965-1979", "1980-1994", "1995-2010", "2011-2025")
First_year  <- c(1965, 1980, 1995, 2011)
Last_year   <- c(1979, 1994, 2010, 2025)

combined_results <- list()

data_dir <- here("data", "climate_raw")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

list.files(here(), pattern = ".nc", full.names = TRUE) |> file.remove()

options(timeout = 10000)

#### LOOP ####

library(progress)

countries <- unique(m_coord_c$country_code)

total_steps <- length(countries) * length(Period)

pb <- progress_bar$new(
  format = "  downloading [:bar] :percent eta: :eta",
  total = total_steps,
  clear = FALSE,
  width = 60
)

combined_results <- list()

for (country in countries) {
  
  country_coord <- subset(m_coord_c, country_code == country)
  
  cat("\nProcessing:", unique(country_coord$Country.Name), "\n")
  
  # ---- coords ----
  sf_point <- st_as_sf(
    country_coord,
    coords = c("transect_lon", "transect_lat"),
    crs = 3035
  )
  
  sf_wgs84 <- st_transform(sf_point, 4326)
  
  coords_mat <- st_coordinates(sf_wgs84)
  
  country_coord$lon <- coords_mat[,1]
  country_coord$lat <- coords_mat[,2]
  
  bbox <- st_bbox(sf_wgs84) %>%
    st_as_sfc() %>%
    st_as_sf()
  
  coords <- country_coord[, c("lon", "lat")]
  
  country_results <- list()
  
  for (i in seq_along(Period)) {
    
    cat("  Period:", Period[i], "\n")
    
    out_file <- file.path(
      temp_dir,
      paste0("temp_", country, "_", Period[i], ".tif")
    )
    
    # ---- download & raster ----
    extract_nc_value(
      first_year     = First_year[i],
      last_year      = Last_year[i],
      local_file     = FALSE,
      file_path      = data_dir,
      sml_chunk      = Period[i],
      spatial_extent = bbox,
      clim_variable  = "mean temp",
      statistic      = "mean",
      grid_size      = 0.1,
      ecad_v         = NULL,
      write_raster   = TRUE,
      out            = out_file,
      return_data    = TRUE
    )
    
    r <- terra::rast(out_file)
    
    layer_names <- names(r)
    
    # ---- extract daily temps ----
    for (j in seq_len(nlyr(r))) {
      
      temp_vals <- terra::extract(
        r[[j]],
        coords,
        method = "simple",
        ID = FALSE
      )
      
      # ⚠️ IMPORTANT: adapta segons format real
      date_j <- as.Date(layer_names[j])
      
      df_tmp <- data.frame(
        transect_id = country_coord$transect_id,
        date        = date_j,
        temp        = temp_vals[,1]
      )
      
      country_results[[length(country_results) + 1]] <- df_tmp
    }
    
    pb$tick()   # 🔥 tick quan TOT el període està fet
  }
  
  combined_results[[country]] <- bind_rows(country_results)
  
  cat("Done:", country, "\n")
}

#### Combine final ####

climate_daily <- bind_rows(combined_results) %>%
  mutate(
    year = lubridate::year(date),
    julian_day = lubridate::yday(date)
  )

#### Save ####
write.csv(
  climate_daily,
  file.path(output_dir, "daily_temperature.csv"),
  row.names = FALSE
)
