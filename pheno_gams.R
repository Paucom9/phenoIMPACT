# Clean environment and set the working directory to the location of the data files
rm(list = ls())  # Remove all objects from the current R session to ensure a clean working environment
setwd("D:/phenoIMPACT project/ebms data/eBMS dataset 2025")  

# Load required libraries
# ----
library(data.table)  # For efficient data handling
library(mgcv)        # For generalized additive models
library(dplyr)       # For data manipulation
library(tidyr)       # For data tidying
library(broom)       # For converting statistical analysis objects into tidy data frames
library(stringr)     # For string manipulation
library(lubridate)   # For easy and intuitive work with dates and times
library(doParallel)  # For increasing loop performance
library(changepoint) # For change point analyses
library(sf)
# ----


# ---- Data Import and Preparation ---- #
#-----
# eBMS data

# Import butterfly count data
ebms_count_df <- read.csv("ebms_count.csv", sep = ",", dec = ".")
# Import visit data
ebms_visit_df <- read.csv("ebms_visit.csv", sep = ",", dec = ".")
# Import climate region data
ebms_clim_df <- read.csv("ebms_transect_climate.csv", sep = ",", dec = ".")
# Import transect coordinates
ebms_coord_df <- read.csv("ebms_transect_coord.csv", sep = ",", dec = ".")
# Import country codes
country_codes <- read.csv("country_codes.csv", sep = ";", dec = ".")

# Extract bms_id from transect_id and select relevant columns
ebms_clim_df <- ebms_clim_df %>%
  mutate(bms_id = str_extract(transect_id, "^[^.]*")) %>%
  dplyr::select(bms_id, transect_id, genzname)

## Transform data frames to data tables

m_count <- data.table(ebms_count_df)
m_visit <- data.table(ebms_visit_df)
m_clim <- data.table(ebms_clim_df)
dt_country_cod <- data.table(country_codes)

## Change column names

setnames(m_visit, c('transect_id', 'visit_date'), c('SITE_ID', 'DATE'))
setnames(m_count, c('transect_id', 'visit_date','species_name', 'count'),
         c('SITE_ID', 'DATE', 'SPECIES', 'COUNT'))
setnames(m_clim, c('transect_id', 'genzname'),
         c('SITE_ID', 'RCLIM'))

# Perform a left join to add RCLIM from m_clim to m_visit based on bms_id and SITE_ID
m_visit <- m_visit[m_clim, on = .(bms_id, SITE_ID), nomatch = 0]

## Perform a left join to merge m_clim into m_count

m_count <- left_join(m_count, m_clim, by = c("SITE_ID", "bms_id"))

# Merge m_count with dt_country_cod to include country_code
m_count <- merge(m_count, dt_country_cod, by = "bms_id", all.x = TRUE)

# Create and ID factor
m_count$ID <- paste(m_count$SPECIES, m_count$SITE_ID, sep = "_")

# Year as factor
m_count$year <- as.factor(m_count$year)
m_visit$year <- as.factor(m_visit$year)

# Check count and visit data
head(m_count)
head(m_visit)



# --- Calculate pheno estimates  --- #

# Create an anchor argument to add zeros (0) before and after the monitoring season 
# This ensures that the flight curve starts and end at zero

anchor <- data.table(
  bms_id = rep(NA_character_, 151),
  visit_id = rep(NA_integer_, 151),
  SITE_ID = rep(NA_character_, 151),
  DATE = rep(as.Date(NA), 151),
  year = rep(NA_integer_, 151),
  month = rep(NA_integer_, 151),
  day = rep(NA_integer_, 151),
  SPECIES = rep(NA_character_, 151),
  COUNT = rep(0, 151),
  RCLIM = rep(NA_character_, 151),
  country_code = rep(NA_character_, 151),
  Country.Name = rep(NA_character_, 151),
  ID = rep(NA_character_, 151),
  ebms_partner = rep(NA, 151),
  week = rep(NA_integer_, 151),
  # Assign a single sequence for julian_day without repetition
  julian_day = c(1:59, 274:365)
)

# Function to find local maxima (i.e. peak abundance)

find_peaks <- function(x = pred$pred,
                       ignore_threshold = 0,
                       span = 3,
                       strict = TRUE) {
  range_x <- range(x, finite = TRUE)
  min_x <- range_x[1]
  max_x <- range_x[2]
  x <- ifelse(!is.finite(x), min_x, x)
  # the next two lines cater for the case when max_x < 0, which is quite
  # common with logs
  delta <- max_x - min_x
  top_flag <- ignore_threshold > 0.0
  scaled_threshold <- delta * abs(ignore_threshold)
  pks <- splus2R::peaks(x = x, span = span, strict = strict)
  if (abs(ignore_threshold) < 1e-5)
    return(pks)
  if (top_flag) {
    return(ifelse(x - min_x > scaled_threshold, pks, FALSE))
  } else {
    return(ifelse(max_x - x > scaled_threshold, pks, FALSE))
  }
}

# Initialize an empty data frame to store phenology estimates
phenology_estimates <- data.frame(ID = numeric(), 
                                  YEAR = numeric(), 
                                  SPECIES = character(), 
                                  SITE_ID = character(), 
                                  ONSET_mean = numeric(),
                                  ONSET_var = numeric(), 
                                  OFFSET_mean = numeric(), 
                                  OFFSET_var = numeric(), 
                                  PEAKDAY = numeric(), 
                                  FIRST_PEAK = numeric(),
                                  LAST_PEAK = numeric(),
                                  FLIGHT_LENGTH_mean = numeric(),
                                  FLIGHT_LENGTH_var = numeric(),
                                  n_peaks = numeric(),
                                  stringsAsFactors = FALSE)


# Create a progress bar with an estimated total count
pb <- progress_estimated(length(unique(m_count$ID)))

for(id in unique(m_count$ID)){
  
  sub_count <- m_count[ID == id]
  site <- unique(sub_count$SITE_ID)
  
  # Update the progress bar
  pb$tick()$print()
  
  for(YEAR in unique(sub_count$year)){
    
    sub_count_year <- sub_count[year == YEAR]
    sub_visit <- m_visit[year == YEAR & SITE_ID == site]
    
    if(nrow(sub_count_year) >= 3 & nrow(sub_visit) >= 10){
      
      # Identify dates in sub_visit not present in sub_count_year
      missing_dates <- sub_visit[!sub_visit$DATE %in% sub_count_year$DATE, ]
      
      # Create new rows in sub_count_year for missing dates with COUNT = 0
      missing_rows <- missing_dates %>%
        mutate(COUNT = 0)
      
      # Include zero counts in the species*site df
      all_counts <- bind_rows(sub_count_year, missing_rows)
      
      # Convert DATE column to Date class
      all_counts <- all_counts %>%
        mutate(DATE = as.Date(DATE))
      
      # Calculate Julian day
      all_counts <- all_counts %>%
        mutate(julian_day = yday(DATE))
      
      # Bind the anchor rows with the counts
      all_counts <- rbind(all_counts, anchor)
      
      tryCatch({
        
        # Fit a GAM model
        gam_model <- gam(COUNT ~ s(julian_day), data = all_counts, family = nb)
        
        # Initialize phenology estimates as NA by default
        onset_mean <- NA
        offset_mean <- NA
        flight_length_mean <- NA
        onset_var <- NA
        offset_var <- NA
        flight_length_var <- NA
        peak_day <- NA
        first_peak <- NA
        last_peak <- NA
        n_peaks <- NA
        
        # Check if the GAM model fit was successful
        if (!inherits(gam_model, "try-error")) {
          
          # Create a sequence of Julian days from 1 to 365
          julian_days <- 1:365
          new_data <- data.frame(julian_day = julian_days)
          # Predict COUNT using the GAM model
          predict_count <- predict(gam_model, newdata = new_data, type = "response")
          
          plot(predict_count)
          
          pheno <- data.frame(Julian_Day = julian_days, Predicted_Count = predict_count)
          
          # --- Phenology estimates --- #
          peaks <- which(unlist(find_peaks(pheno$Predicted_Count, ignore_threshold = 0.1, span = 11)))
          n_peaks <- length(peaks)
          
          # First and last peaks
          if (n_peaks > 0) {
            first_peak <- peaks[1]         # First peak
            last_peak <- peaks[n_peaks]    # Last peak
          } else {
            first_peak <- NA               # No peaks detected
            last_peak <- NA                # No peaks detected
          }
          
          peak_day <- which.max(pheno$Predicted_Count) # Day of the flight curve in which the abundance is maximum
          
          # Identify onset and offset of butterfly activity using change point analysis
          cp_mean <- cpt.mean(pheno$Predicted_Count, method = "PELT", penalty = "Manual", pen.value = 2)
          cp_var <- cpt.var(pheno$Predicted_Count, method = "PELT", penalty = "Manual", pen.value = 0.05)
          
          # Extract the change points
          change_points_mean <- cpts(cp_mean)
          change_points_var <- cpts(cp_var)
          
          # Check if n_peaks == 0, or if change_points_mean is empty or contains Inf
          # The detected change points, extracted with cpts(cp_mean), indicate shifts in butterfly activity levels
          # The first and last detected change points are used to define the activity period
            if (n_peaks > 0 && length(change_points_mean) > 0 && !any(is.infinite(change_points_mean))) {
            onset_mean <- min(change_points_mean) # First day of appearance with cp_mean method
            offset_mean <- max(change_points_mean) # Last day of appearance with cp_mean method
            flight_length_mean <- length(onset_mean:offset_mean) # Length of the flight period
            
            onset_var <- min(change_points_var) # First day of appearance with cp_var method
            offset_var <- max(change_points_var) # Last day of appearance with cp_var method
            flight_length_var <- length(onset_var:offset_var) # Length of the flight period
          }
        }
        
      }, error = function(e) {
        # Print the error message for debugging
        print(paste("Error for ", id, " in ", YEAR, ":", e))
      })
      
      # Save the results into the data frame
      phenology_estimates <- rbind(phenology_estimates, 
                                   data.frame(ID = id, 
                                              YEAR = YEAR, 
                                              SPECIES = unique(sub_count_year$SPECIES), 
                                              SITE_ID = unique(sub_count_year$SITE_ID), 
                                              ONSET_mean = onset_mean,
                                              ONSET_var = onset_var, 
                                              OFFSET_mean = offset_mean, 
                                              OFFSET_var = offset_var, 
                                              PEAKDAY = peak_day, 
                                              FIRST_PEAK = first_peak,
                                              LAST_PEAK = last_peak,
                                              FLIGHT_LENGTH_mean = flight_length_mean,
                                              FLIGHT_LENGTH_var = flight_length_var,
                                              N_PEAKS = n_peaks))
      
    }
  }
}

#----

head(phenology_estimates)

# Specify the file path and name
file_path <- "D:/phenoIMPACT project/outputs/pheno_estimates.csv"

# Save as a CSV file
write.csv(phenology_estimates, file = file_path, row.names = FALSE)

# ----


