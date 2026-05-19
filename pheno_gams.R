# ============================================================================================ #
# pheno_gams.R
#
# Author: Pau Colom
# Date: 2026-02-19
#
# Desciription: This script fits generalized additive models (GAMs) to butterfly count data 
# from the eBMS to estimate phenological parameters such as onset, offset, peak day, 
# and flight length for each species × site × year combination. 
# The script includes data import and preparation, site and species filtering, GAM fitting, 
# phenology estimation, and exploratory analyses of the resulting phenology estimates.
#
# ============================================================================================ #


#### Load required libraries ####
# ---
library(data.table)  # For efficient data handling
library(mgcv)        # For generalized additive models
library(dplyr)       # For data manipulation
library(tidyr)       # For data tidying
library(broom)       # For converting statistical analysis objects into tidy data frames
library(stringr)     # For string manipulation
library(lubridate)   # For easy and intuitive work with dates and times
library(doParallel)  # For increasing loop performance
library(changepoint) # For change point analyses
library(sf)          # For handling spatial data
library(ggplot2)   # For data visualization
library(reshape2)     # For reshaping data frames
library(here)       # For constructing file paths in a way that is independent of the operating system
# ---

#### Data Import and Preparation ####
# ---
here::here() # Check the current working directory
# ---

# --- eBMS data
# Import butterfly count data
ebms_count_df  <- read.csv(here::here("data", "ebms_count.csv"), sep = ",", dec = ".")
# Import visit data
ebms_visit_df  <- read.csv(here::here("data", "ebms_visit.csv"), sep = ",", dec = ".")
# Import climate region data
ebms_clim_df   <- read.csv(here::here("data", "ebms_transect_climate.csv"), sep = ",", dec = ".")
# Import transect coordinates
ebms_coord_df  <- read.csv(here::here("data", "ebms_transect_coord.csv"), sep = ",", dec = ".")
# Import country codes
country_codes  <- read.csv(here::here("data", "country_codes.csv"), sep = ";", dec = ".")

# --- Extract bms_id from transect_id and select relevant columns
ebms_clim_df <- ebms_clim_df %>%
  mutate(bms_id = str_extract(transect_id, "^[^.]*")) %>%
  dplyr::select(bms_id, transect_id, genzname)

# --- Transform data frames to data tables
m_count <- data.table(ebms_count_df)
m_visit <- data.table(ebms_visit_df)
m_clim <- data.table(ebms_clim_df)
dt_country_cod <- data.table(country_codes)

# --- Change column names
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
# ---

##### Site filtering #####

# Explore mean number of visits x year in each BMS_ID #
vis_sy <- m_visit[
  , .(n_visits = uniqueN(visit_id)),
  by = .(bms_id, SITE_ID, year)
]

vis_sy[
  , .(
    mean_visits_year = mean(n_visits),
    sd_visits_year   = sd(n_visits)
  ),
  by = bms_id
][order(-mean_visits_year)]


# Count visits per site x year #
vis_year <- m_visit[
  , .(n_visits = uniqueN(visit_id)), 
  by = .(SITE_ID, year)
]

# Keep only site x year combinations with at least 10 visits
vis_year_ok <- vis_year[n_visits >= 10]

# Kepp sites with at least 10 years of adequate visits 
sites_ok <- vis_year_ok[
  , .N, by = SITE_ID
][N >= 10, SITE_ID]

# Subsets datasets to keep only the selected sites
m_visit_filt <- m_visit[SITE_ID %in% sites_ok]
m_count_filt <- m_count[SITE_ID %in% sites_ok]
head(m_visit_filt)
head(m_count_filt)

# Explore number of sites per BMS_ID after filtering #
m_count_filt[, uniqueN(SITE_ID), by = bms_id][order(-V1)]
# ---

##### Sampling support / data quality after site filtering #####

# Prepare dates
m_visit_filt[, DATE := as.Date(DATE)]
m_count_filt[, DATE := as.Date(DATE)]

m_visit_filt[, visit_doy := lubridate::yday(DATE)]
m_count_filt[, obs_doy := lubridate::yday(DATE)]

# Unique site-year visit dates
visit_dates_sy <- unique(
  m_visit_filt[, .(
    SITE_ID,
    year,
    bms_id,
    visit_id,
    DATE,
    visit_doy
  )]
)

# Site-year visit summary
visit_summary_sy <- visit_dates_sy[
  , .(
    n_visits_site_year = uniqueN(visit_id),
    first_visit_doy = min(visit_doy, na.rm = TRUE),
    last_visit_doy  = max(visit_doy, na.rm = TRUE)
  ),
  by = .(SITE_ID, year, bms_id)
]

# Species-site-year positive observation summary
obs_summary_ssy <- m_count_filt[
  !is.na(COUNT) & COUNT > 0,
  .(
    n_positive_dates = uniqueN(DATE),
    first_obs_doy = min(obs_doy, na.rm = TRUE),
    last_obs_doy  = max(obs_doy, na.rm = TRUE)
  ),
  by = .(ID, SPECIES, SITE_ID, year, bms_id)
]

# Keep the same species-site-year candidates that enter the GAM loop:
# ≥3 positive dates and ≥10 site-year visits
sampling_support_ssy <- obs_summary_ssy |>
  as.data.frame() |>
  dplyr::left_join(
    visit_summary_sy |>
      as.data.frame() |>
      dplyr::select(
        SITE_ID,
        year,
        n_visits_site_year,
        first_visit_doy,
        last_visit_doy
      ),
    by = c("SITE_ID", "year")
  ) |>
  dplyr::filter(
    n_positive_dates >= 3,
    n_visits_site_year >= 10
  )

# Count real zero visits before first observation and after last observation
zero_support_ssy <- sampling_support_ssy |>
  dplyr::select(
    ID,
    SPECIES,
    SITE_ID,
    year,
    first_obs_doy,
    last_obs_doy
  ) |>
  dplyr::left_join(
    visit_dates_sy |>
      as.data.frame() |>
      dplyr::select(
        SITE_ID,
        year,
        visit_id,
        visit_doy
      ),
    by = c("SITE_ID", "year")
  ) |>
  dplyr::group_by(ID, SPECIES, SITE_ID, year) |>
  dplyr::summarise(
    n_zero_visits_before_first_obs =
      dplyr::n_distinct(visit_id[visit_doy < first_obs_doy]),
    
    n_zero_visits_after_last_obs =
      dplyr::n_distinct(visit_id[visit_doy > last_obs_doy]),
    
    .groups = "drop"
  )

sampling_support_ssy <- sampling_support_ssy |>
  dplyr::left_join(
    zero_support_ssy,
    by = c("ID", "SPECIES", "SITE_ID", "year")
  ) |>
  dplyr::mutate(
    gap_first_visit_first_obs = first_obs_doy - first_visit_doy,
    gap_last_obs_last_visit = last_visit_doy - last_obs_doy,
    
    has_zero_before_first_obs =
      n_zero_visits_before_first_obs >= 1,
    
    has_zero_after_last_obs =
      n_zero_visits_after_last_obs >= 1,
    
    has_zero_before_and_after =
      has_zero_before_first_obs & has_zero_after_last_obs
  )

# Quick inspection
sampling_quality_summary <- sampling_support_ssy |>
  dplyr::summarise(
    n_species_site_year = dplyr::n(),
    
    mean_n_visits_site_year =
      mean(n_visits_site_year, na.rm = TRUE),
    
    mean_n_positive_dates =
      mean(n_positive_dates, na.rm = TRUE),
    
    mean_zero_visits_before_first_obs =
      mean(n_zero_visits_before_first_obs, na.rm = TRUE),
    
    mean_zero_visits_after_last_obs =
      mean(n_zero_visits_after_last_obs, na.rm = TRUE),
    
    median_gap_first_visit_first_obs =
      median(gap_first_visit_first_obs, na.rm = TRUE),
    
    median_gap_last_obs_last_visit =
      median(gap_last_obs_last_visit, na.rm = TRUE),
    
    prop_has_zero_before_first_obs =
      mean(has_zero_before_first_obs, na.rm = TRUE),
    
    prop_has_zero_after_last_obs =
      mean(has_zero_after_last_obs, na.rm = TRUE),
    
    prop_has_zero_before_and_after =
      mean(has_zero_before_and_after, na.rm = TRUE)
  )

sampling_quality_by_bms <- sampling_support_ssy |>
  dplyr::group_by(bms_id) |>
  dplyr::summarise(
    n_species_site_year = dplyr::n(),
    
    mean_n_visits_site_year =
      mean(n_visits_site_year, na.rm = TRUE),
    
    mean_n_positive_dates =
      mean(n_positive_dates, na.rm = TRUE),
    
    prop_has_zero_before_first_obs =
      mean(has_zero_before_first_obs, na.rm = TRUE),
    
    prop_has_zero_after_last_obs =
      mean(has_zero_after_last_obs, na.rm = TRUE),
    
    prop_has_zero_before_and_after =
      mean(has_zero_before_and_after, na.rm = TRUE),
    
    median_gap_first_visit_first_obs =
      median(gap_first_visit_first_obs, na.rm = TRUE),
    
    median_gap_last_obs_last_visit =
      median(gap_last_obs_last_visit, na.rm = TRUE),
    
    .groups = "drop"
  ) |>
  dplyr::arrange(desc(n_species_site_year))

sampling_quality_by_year <- sampling_support_ssy |>
  dplyr::group_by(year) |>
  dplyr::summarise(
    n_species_site_year = dplyr::n(),
    
    prop_has_zero_before_first_obs =
      mean(has_zero_before_first_obs, na.rm = TRUE),
    
    prop_has_zero_after_last_obs =
      mean(has_zero_after_last_obs, na.rm = TRUE),
    
    prop_has_zero_before_and_after =
      mean(has_zero_before_and_after, na.rm = TRUE),
    
    .groups = "drop"
  )

sampling_quality_summary |>
  tibble::as_tibble() |>
  print(n = Inf, width = Inf)
sampling_quality_by_bms |>
  tibble::as_tibble() |>
  print(n = Inf, width = Inf)
sampling_quality_by_year |>
  tibble::as_tibble() |>
  print(n = Inf, width = Inf)

##### Save sampling support in phenology_estimates format #####

phenology_sampling_support <- sampling_support_ssy |>
  dplyr::mutate(
    YEAR = as.integer(as.character(year)),
    ID = as.character(ID),
    SPECIES = as.character(SPECIES),
    SITE_ID = as.character(SITE_ID)
  ) |>
  dplyr::select(
    ID,
    YEAR,
    SPECIES,
    SITE_ID,
    bms_id,
    n_visits_site_year,
    n_positive_dates,
    first_visit_doy,
    last_visit_doy,
    first_obs_doy,
    last_obs_doy,
    n_zero_visits_before_first_obs,
    n_zero_visits_after_last_obs,
    gap_first_visit_first_obs,
    gap_last_obs_last_visit,
    has_zero_before_first_obs,
    has_zero_after_last_obs,
    has_zero_before_and_after
  )

head(phenology_sampling_support)
str(phenology_sampling_support)

write.csv(
  phenology_sampling_support,
  file = here("output", "phenology_sampling_support_allspp.csv"),
  row.names = FALSE
)


##### Functions #####


# --- Calculate pheno estimates

# Function to find local maxima (i.e. peak abundance)
find_peaks <- function(x,
                       ignore_threshold = 0,
                       span = 3,
                       strict = TRUE) {
  
  # Replace non-finite values
  x[!is.finite(x)] <- min(x, na.rm = TRUE)
  
  # Detect local maxima
  pks <- splus2R::peaks(x = x, span = span, strict = strict)
  
  # If no threshold → return raw peaks
  if (ignore_threshold <= 0) return(pks)
  
  # Threshold relative to maximum (🔥 key change)
  max_x <- max(x, na.rm = TRUE)
  
  return(ifelse(x > ignore_threshold * max_x, pks, FALSE))
}


##### Loop gams and extract pheno and abundance metrics #####

phenology_estimates <- data.frame(
  ID = character(),
  YEAR = integer(),
  SPECIES = character(),
  SITE_ID = character(),
  ABUND_INDEX = numeric(),
  ONSET_mean = numeric(),
  ONSET_var = numeric(),
  OFFSET_mean = numeric(),
  OFFSET_var = numeric(),
  PEAKDAY = numeric(),
  FIRST_PEAK = numeric(),
  LAST_PEAK = numeric(),
  FLIGHT_LENGTH_mean = numeric(),
  FLIGHT_LENGTH_var = numeric(),
  N_PEAKS = numeric(),
  stringsAsFactors = FALSE
)

pb <- progress_estimated(length(unique(m_count_filt$ID)))

for (id in unique(m_count_filt$ID)) {
  
  sub_count <- m_count_filt[ID == id]
  site <- unique(sub_count$SITE_ID)
  
  pb$tick()$print()
  
  for (sp in unique(sub_count$SPECIES)) {
    
    sub_count_species <- sub_count[SPECIES == sp]
    
    for (yr in unique(sub_count_species$year)) {
      
      sub_count_year <- sub_count_species[year == yr]
      sub_visit <- m_visit_filt[year == yr & SITE_ID == site]
      
      onset_mean <- offset_mean <- flight_length_mean <- NA
      onset_var  <- offset_var  <- flight_length_var  <- NA
      peak_day <- first_peak <- last_peak <- n_peaks <- NA
      abund_index <- NA
      
      # Require ≥3 weeks with records and ≥10 visits
      if (
        uniqueN(as.Date(sub_count_year$DATE)) >= 3 &
        data.table::uniqueN(sub_visit$visit_id) >= 10
      ) {
        
        # Zero-fill missing visits within the monitoring season
        missing_dates <- sub_visit[!DATE %in% sub_count_year$DATE]
        
        missing_rows <- missing_dates %>%
          mutate(
            SPECIES = sp,
            COUNT = 0
          )
        
        all_counts <- bind_rows(sub_count_year, missing_rows) %>%
          mutate(
            DATE = as.Date(DATE),
            julian_day = yday(DATE)
          )
        
        # Structural zeros outside main monitoring season
        anchor_year <- data.table(
          SITE_ID = site,
          year = yr,
          SPECIES = sp,
          COUNT = 0,
          julian_day = c(1:30, 335:365)
        )
        
        all_counts <- rbind(all_counts, anchor_year, fill = TRUE)
        
        tryCatch({
          
          gam_model <- gam(
            COUNT ~ s(julian_day),
            data = all_counts,
            family = nb
          )
          
          julian_days <- 1:365
          
          predict_count <- predict(
            gam_model,
            newdata = data.frame(julian_day = julian_days),
            type = "response"
          )
          
          abund_index <- sum(predict_count, na.rm = TRUE)
          
          pheno <- data.frame(
            Julian_Day = julian_days,
            Predicted_Count = predict_count
          )
          
          peaks <- which(unlist(find_peaks(
            pheno$Predicted_Count,
            ignore_threshold = 0.2,
            span = 11
          )))
          
          n_peaks <- length(peaks)
          
          if (n_peaks > 0) {
            first_peak <- peaks[1]
            last_peak  <- peaks[n_peaks]
          }
          
          peak_day <- which.max(pheno$Predicted_Count)
          
          cp_mean <- cpt.mean(
            pheno$Predicted_Count,
            method = "PELT",
            penalty = "Manual",
            pen.value = 2
          )
          
          cp_var <- cpt.var(
            pheno$Predicted_Count,
            method = "PELT",
            penalty = "Manual",
            pen.value = 0.05
          )
          
          cps_mean <- cpts(cp_mean)
          cps_var  <- cpts(cp_var)
          
          if (n_peaks > 0 &&
              length(cps_mean) > 0 &&
              length(cps_var) > 0 &&
              !any(is.infinite(cps_mean)) &&
              !any(is.infinite(cps_var))) {
            
            onset_mean <- min(cps_mean)
            offset_mean <- max(cps_mean)
            flight_length_mean <- offset_mean - onset_mean + 1
            
            onset_var <- min(cps_var)
            offset_var <- max(cps_var)
            flight_length_var <- offset_var - onset_var + 1
          }
          
        }, error = function(e) {
          message("Error for ", id, " / ", sp, " / ", yr, ": ", e$message)
        })
        
        phenology_estimates <- rbind(
          phenology_estimates,
          data.frame(
            ID = id,
            YEAR = yr,
            SPECIES = sp,
            SITE_ID = site,
            ABUND_INDEX = abund_index,
            ONSET_mean = onset_mean,
            ONSET_var = onset_var,
            OFFSET_mean = offset_mean,
            OFFSET_var = offset_var,
            PEAKDAY = peak_day,
            FIRST_PEAK = first_peak,
            LAST_PEAK = last_peak,
            FLIGHT_LENGTH_mean = flight_length_mean,
            FLIGHT_LENGTH_var = flight_length_var,
            N_PEAKS = n_peaks
          )
        )
      }
    }
  }
}
#---


##### Save the data #####
head(phenology_estimates)
str(phenology_estimates)

# Save as CSV inside project
write.csv(
  phenology_estimates,
  file = here("output", "pheno_abund_estimates_allspp.csv"),
  row.names = FALSE
)
#---


#### Plot pheno-curve examples ####

plot_pheno_curve <- function(species_pick,
                             site_pick,
                             year_pick,
                             m_count_filt,
                             m_visit_filt,
                             pen_mean = 2,
                             pen_var  = 0.05,
                             save_plot = FALSE) {
  
  
  # Filter data
  sub_count <- m_count_filt %>%
    filter(SPECIES == species_pick,
           SITE_ID == site_pick,
           year == year_pick)
  
  sub_visit <- m_visit_filt %>%
    filter(SITE_ID == site_pick,
           year == year_pick)
  
  if(nrow(sub_count) == 0)
    stop("No data for this species-site-year combination.")
  
  # Dates
  sub_count$DATE <- as.Date(sub_count$DATE)
  sub_visit$DATE <- as.Date(sub_visit$DATE)
  
  sub_count$julian_day <- yday(sub_count$DATE)
  sub_visit$julian_day <- yday(sub_visit$DATE)
  
  # Zero-fill
  missing_dates <- sub_visit %>%
    filter(!DATE %in% sub_count$DATE)
  
  missing_rows <- missing_dates %>%
    mutate(COUNT = 0) %>%
    select(julian_day, COUNT)
  
  all_counts <- bind_rows(
    sub_count %>% select(julian_day, COUNT),
    missing_rows
  )
  
  # Structural zeros
  anchor <- data.frame(
    julian_day = c(1:30, 335:365),
    COUNT = 0
  )
  
  all_counts <- bind_rows(all_counts, anchor)
  
  # GAM
  gam_model <- gam(COUNT ~ s(julian_day),
                   data = all_counts,
                   family = nb())
  
  julian_days <- 1:365
  pred <- predict(gam_model,
                  newdata = data.frame(julian_day = julian_days),
                  type = "response")
  
  pheno <- data.frame(
    Julian_Day = julian_days,
    Predicted_Count = pred
  )
  
  # Change points
  cp_mean <- cpt.mean(pheno$Predicted_Count,
                      method = "PELT",
                      penalty = "Manual",
                      pen.value = pen_mean)
  
  cp_var <- cpt.var(pheno$Predicted_Count,
                    method = "PELT",
                    penalty = "Manual",
                    pen.value = pen_var)
  
  cps_mean <- cpts(cp_mean)
  cps_var  <- cpts(cp_var)
  
  onset_mean  <- min(cps_mean)
  offset_mean <- max(cps_mean)
  onset_var   <- min(cps_var)
  offset_var  <- max(cps_var)
  
  # Peak
  peak_day   <- which.max(pheno$Predicted_Count)
  peak_value <- max(pheno$Predicted_Count)
  
  # Plot object
  plot_obj <- ggplot() +
    geom_line(data = pheno,
              aes(x = Julian_Day, y = Predicted_Count),
              colour = "#1B9E77",
              linewidth = 1.4) +
    geom_point(data = all_counts,
               aes(x = julian_day, y = COUNT),
               colour = "grey40",
               alpha = 0.6,
               size = 2) +
    geom_vline(xintercept = onset_mean,
               colour = "#D95F02",
               linetype = "dashed",
               linewidth = 1.4) +
    geom_vline(xintercept = offset_mean,
               colour = "#D95F02",
               linetype = "dashed",
               linewidth = 1.4) +
    
    geom_vline(xintercept = onset_var,
               colour = "#7570B3",
               linetype = "dotted",
               linewidth = 1.2
    ) +
    geom_vline(xintercept = offset_var,
               colour = "#7570B3",
               linetype = "dotted",
               linewidth = 1.2) +
    
    geom_point(aes(x = peak_day, y = peak_value),
               colour = "red",
               size = 4) +
    
    geom_segment(aes(x = peak_day, xend = peak_day,
                     y = 0, yend = peak_value),
                 colour = "red",
                 linewidth = 1.2
    ) +
    
    labs(
      title = paste0("Species: ", species_pick),
      subtitle = paste0("Site: ", site_pick,
                        "   |   Year: ", year_pick),
      x = "Day of the year",
      y = "Relative abundance"
    ) +
    
    theme_classic(base_size = 14) +
    theme(plot.title = element_text(face = "bold"))
  
  # Save if requested
  if(save_plot){
    
    dir.create(here::here("output", "figures"), showWarnings = FALSE)
    
    filename <- paste0("pheno_",
                       gsub(" ", "_", species_pick), "_",
                       site_pick, "_",
                       year_pick, ".png")
    
    ggsave(
      filename = here::here("output", "figures", filename),
      plot = plot_obj,
      width = 8,
      height = 5,
      dpi = 300
    )
  }
  
  return(plot_obj)
}


plot_pheno_curve(
  species_pick = "Boloria selene",
  site_pick    = "FIBMS.246",
  year_pick    = "2023",
  m_count_filt = m_count_filt,
  m_visit_filt = m_visit_filt,
  save_plot    = TRUE
)
#---
