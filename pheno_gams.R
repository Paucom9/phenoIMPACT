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
library(sf)          # For handling spatial data
library(ggplot2)   # For data visualization
library(reshape2)

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

<<<<<<< HEAD


# --- Calculate pheno estimates  --- #
=======
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
>>>>>>> d378979 (site and species filtering + loop fixed issues + ouput exploration)

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


#####

##### Species filtering #####

# --- Explore species with most site x year combinations (minimum 3 positive counts) --- #

# Keep only positive counts
dt_pos <- m_count_filt[COUNT > 0]

# Count distinc weeeks per species x site x year
spec_sy_weeks <- dt_pos[
  , .(n_weeks = uniqueN(paste(year, month, day))), 
  by = .(SPECIES, SITE_ID, year)
]

# Keep only species x site x year combinations with at least 3 weeks with positive counts
spec_sy_ok <- spec_sy_weeks[n_weeks >= 3]

# Count valid site × year combos per species
spec_cov <- spec_sy_ok[
  , .(n_site_years_3w = .N),
  by = SPECIES
][order(-n_site_years_3w)]

# --- Subset Maniola jurtina (species with more data) --- #
mj_keys <- m_count_filt[
  SPECIES == "Maniola jurtina",
  unique(paste(SITE_ID, DATE))
]

mj_visit <- m_visit_filt[
  paste(SITE_ID, DATE) %in% mj_keys
]

mj_count <- m_count_filt[SPECIES == "Maniola jurtina"]



#####


# --- Calculate pheno estimates  --- #

# Function to find local maxima (i.e. peak abundance)
<<<<<<< HEAD

find_peaks <- function(x = pred$pred,
=======
find_peaks <- function(x,
>>>>>>> d378979 (site and species filtering + loop fixed issues + ouput exploration)
                       ignore_threshold = 0,
                       span = 3,
                       strict = TRUE) {
  range_x <- range(x, finite = TRUE)
  min_x <- range_x[1]
  max_x <- range_x[2]
  x <- ifelse(!is.finite(x), min_x, x)
  delta <- max_x - min_x
  top_flag <- ignore_threshold > 0.0
  scaled_threshold <- delta * abs(ignore_threshold)
  pks <- splus2R::peaks(x = x, span = span, strict = strict)
  if (abs(ignore_threshold) < 1e-5) return(pks)
  if (top_flag) {
    return(ifelse(x - min_x > scaled_threshold, pks, FALSE))
  } else {
    return(ifelse(max_x - x > scaled_threshold, pks, FALSE))
  }
}

# Initialize output
phenology_estimates <- data.frame(
  ID = character(),
  YEAR = integer(),
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
  N_PEAKS = numeric(),
  stringsAsFactors = FALSE
)

pb <- progress_estimated(length(unique(mj_count$ID)))

<<<<<<< HEAD
# Create a progress bar with an estimated total count
pb <- progress_estimated(length(unique(m_count$ID)))

for(id in unique(m_count$ID)){
  
  sub_count <- m_count[ID == id]
=======
for (id in unique(mj_count$ID)) {
  
  sub_count <- mj_count[ID == id]
>>>>>>> d378979 (site and species filtering + loop fixed issues + ouput exploration)
  site <- unique(sub_count$SITE_ID)
  
  pb$tick()$print()
  
  for (YEAR in unique(sub_count$year)) {
    
    sub_count_year <- sub_count[year == YEAR]
    sub_visit <- mj_visit[year == YEAR & SITE_ID == site]
    
    # Require ≥3 weeks with records and ≥10 visits
    if (uniqueN(as.Date(sub_count_year$DATE)) >= 3 & nrow(sub_visit) >= 10) {
      
      # Zero-fill missing visits within the monitoring season
      missing_dates <- sub_visit[!DATE %in% sub_count_year$DATE]
      missing_rows <- missing_dates %>% mutate(COUNT = 0)
      
      all_counts <- bind_rows(sub_count_year, missing_rows) %>%
        mutate(
          DATE = as.Date(DATE),
          julian_day = yday(DATE)
        )
      
      # ---- ANCHOR (structural zeros outside March–September) ----
      anchor_year <- data.table(
        SITE_ID = site,
        year = YEAR,
        SPECIES = unique(sub_count_year$SPECIES),
        COUNT = 0,
        julian_day = c(1:30, 305:365)  # Jan–early Feb & Oct–Dec
      )
      
      all_counts <- rbind(all_counts, anchor_year, fill = TRUE)
      
      tryCatch({
        
        gam_model <- gam(COUNT ~ s(julian_day), data = all_counts, family = nb)
        
        onset_mean <- offset_mean <- flight_length_mean <- NA
        onset_var  <- offset_var  <- flight_length_var  <- NA
        peak_day <- first_peak <- last_peak <- n_peaks <- NA
        
        if (!inherits(gam_model, "try-error")) {
          
          julian_days <- 1:365
          predict_count <- predict(
            gam_model,
            newdata = data.frame(julian_day = julian_days),
            type = "response"
          )
          
          pheno <- data.frame(
            Julian_Day = julian_days,
            Predicted_Count = predict_count
          )
          
          peaks <- which(unlist(find_peaks(
            pheno$Predicted_Count,
            ignore_threshold = 0.1,
            span = 11
          )))
          
          n_peaks <- length(peaks)
          if (n_peaks > 0) {
            first_peak <- peaks[1]
            last_peak  <- peaks[n_peaks]
          }
          
          peak_day <- which.max(pheno$Predicted_Count)
          
          cp_mean <- cpt.mean(pheno$Predicted_Count, method = "PELT",
                              penalty = "Manual", pen.value = 2)
          cp_var  <- cpt.var(pheno$Predicted_Count, method = "PELT",
                             penalty = "Manual", pen.value = 0.05)
          
          cps_mean <- cpts(cp_mean)
          cps_var  <- cpts(cp_var)
          
          if (n_peaks > 0 && length(cps_mean) > 0 && !any(is.infinite(cps_mean))) {
            onset_mean <- min(cps_mean)
            offset_mean <- max(cps_mean)
            flight_length_mean <- offset_mean - onset_mean + 1
            
            onset_var <- min(cps_var)
            offset_var <- max(cps_var)
            flight_length_var <- offset_var - onset_var + 1
          }
        }
        
      }, error = function(e) {
        message("Error for ", id, " in ", YEAR, ": ", e$message)
      })
      
      phenology_estimates <- rbind(
        phenology_estimates,
        data.frame(
          ID = id,
          YEAR = YEAR,
          SPECIES = unique(sub_count_year$SPECIES),
          SITE_ID = site,
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

#----

head(phenology_estimates)

# Specify the file path and name
file_path <- "D:/phenoIMPACT project/outputs/pheno_estimates.csv"

# Save as a CSV file
write.csv(phenology_estimates, file = file_path, row.names = FALSE)

# ----
<<<<<<< HEAD
=======


##### Output exploration #####

head(phenology_estimates)
str(phenology_estimates)
summary(phenology_estimates)


# --- Plot distributions of phenology variables --- #####
vars <- c(
  "N_PEAKS",
  "ONSET_mean", "ONSET_var",
  "PEAKDAY",
  "OFFSET_mean", "OFFSET_var",
  "FLIGHT_LENGTH_mean", "FLIGHT_LENGTH_var"
)

ph_long <- phenology_estimates |>
  pivot_longer(cols = all_of(vars),
               names_to = "var",
               values_to = "value") |>
  mutate(var = factor(var, levels = vars))

ph_median <- ph_long |>
  group_by(var) |>
  summarise(
    median = median(value, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(ph_long, aes(x = value)) +
  geom_histogram(bins = 30, fill = "grey70", color = "black") +
  geom_vline(data = ph_median,
             aes(xintercept = median),
             linewidth = 0.7) +
  facet_wrap(~var, scales = "free") +
  theme_minimal()

# --- Correlation matrix of phenology variables --- #####
vars_pheno <- c(
  "ONSET_mean",
  "ONSET_var",
  "PEAKDAY",
  "OFFSET_mean",
  "OFFSET_var",
  "FLIGHT_LENGTH_mean",
  "FLIGHT_LENGTH_var"
)

cor_mat <- cor(
  phenology_estimates[, vars_pheno],
  use = "pairwise.complete.obs"
)

cor_df <- melt(cor_mat)

ggplot(cor_df, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", value)),
            size = 3) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0, limits = c(-1, 1)
  ) +
  coord_fixed() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank()
  )


# --- Test relationships between phenology variables and latitude/longitude --- #####

# Merge phenology estimates with transect coordinates
pheno_geo <- phenology_estimates |>
  left_join(
    ebms_coord_df,
    by = c("SITE_ID" = "transect_id")
  )
head(pheno_geo)


vars_pheno <- c(
  "ONSET_mean",
  "ONSET_var",
  "PEAKDAY",
  "OFFSET_mean",
  "OFFSET_var",
  "FLIGHT_LENGTH_mean",
  "FLIGHT_LENGTH_var"
)

ph_long_geo <- pheno_geo |>
  select(transect_lat, all_of(vars_pheno)) |>
  pivot_longer(
    cols = all_of(vars_pheno),
    names_to = "var",
    values_to = "value"
  )

ggplot(ph_long_geo, aes(transect_lat, value)) +
  geom_point(alpha = 0.2, size = 0.6) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  facet_wrap(~var, scales = "free_y") +
  theme_minimal() +
  labs(
    x = "North–south gradient (EPSG:3035)",
    y = "Phenological estimate"
  )


>>>>>>> d378979 (site and species filtering + loop fixed issues + ouput exploration)


