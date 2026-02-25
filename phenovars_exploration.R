# ============================================================================================ #
# phenovars_exploration.R
#
# Author: Pau Colom
# Date: 2026-02-20
#
# Desciription: This script explores the phenological variables estimated for each species × 
# site combination. It includes visualizations of the distributions of these variables, 
# a correlation matrix to examine relationships between them, and scatter plots to investigate 
# potential associations with latitude. The goal is to gain insights into the patterns and 
# relationships in the phenological data before conducting further analyses.
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
library(rnaturalearth) # For obtaining natural earth map data
library(ggplot2)   # For data visualization
library(reshape2)     # For reshaping data frames
library(here)       # For constructing file paths in a way that is independent of the operating system
# ---

#### Data Import and Preparation ####
# ---
here::here() # Check the current working directory
# ---
phenology_estimates  <- read.csv(here("output", "pheno_estimates_allspp.csv"), sep = ",", dec = ".")
str(phenology_estimates)

head(phenology_estimates)
str(phenology_estimates)
summary(phenology_estimates)


#### Plot distributions of phenology variables ####

# Exclude rows with NAs
pheno_estimates_clean <- na.exclude(phenology_estimates)
str(pheno_estimates_clean)

vars <- c(
  "N_PEAKS",
  "ONSET_mean", "ONSET_var",
  "PEAKDAY",
  "OFFSET_mean", "OFFSET_var",
  "FLIGHT_LENGTH_mean", "FLIGHT_LENGTH_var"
)

pretty_labels <- c(
  N_PEAKS = "Number of detected peaks",
  ONSET_mean = "Onset (mean)",
  ONSET_var = "Onset (variance)",
  PEAKDAY = "Peak day",
  OFFSET_mean = "Offset (mean)",
  OFFSET_var = "Offset (variance)",
  FLIGHT_LENGTH_mean = "Flight length (mean)",
  FLIGHT_LENGTH_var = "Flight length (variance)"
)

ph_long <- pheno_estimates_clean |>
  pivot_longer(cols = all_of(vars),
               names_to = "var",
               values_to = "value") |>
  mutate(
    var = factor(var, levels = vars),
    var_label = factor(pretty_labels[var],
                       levels = pretty_labels[vars])
  )

ph_median <- ph_long |>
  group_by(var_label) |>
  summarise(
    median = median(value, na.rm = TRUE),
    .groups = "drop"
  )

distr_plot <-  ggplot(ph_long, aes(x = value)) +
  geom_histogram(bins = 30, fill = "grey70", color = "black") +
  geom_vline(data = ph_median,
             aes(xintercept = median),
             linewidth = 0.7) +
  facet_wrap(~var_label, scales = "free") +
  theme_minimal()

distr_plot

ggsave(
  filename = here::here("output", "figures", "distribution_phenovars.png"),
  plot = distr_plot,
  width = 8,
  height = 5,
  dpi = 300
)

#### Correlation matrix of phenology variables ####

vars_pheno <- c(
  "ONSET_mean",
  "ONSET_var",
  "PEAKDAY",
  "OFFSET_mean",
  "OFFSET_var",
  "FLIGHT_LENGTH_mean",
  "FLIGHT_LENGTH_var"
)

pretty_labels <- c(
  ONSET_mean = "Onset (mean)",
  ONSET_var = "Onset (variance)",
  PEAKDAY = "Peak day",
  OFFSET_mean = "Offset (mean)",
  OFFSET_var = "Offset (variance)",
  FLIGHT_LENGTH_mean = "Flight length (mean)",
  FLIGHT_LENGTH_var = "Flight length (variance)"
)

# Correlation matrix
cor_mat <- cor(
  pheno_estimates_clean[, vars_pheno],
  use = "pairwise.complete.obs"
)

# Rename and reorder
cor_mat <- cor_mat[vars_pheno, vars_pheno]
colnames(cor_mat) <- pretty_labels[vars_pheno]
rownames(cor_mat) <- pretty_labels[vars_pheno]

cor_df <- melt(cor_mat)

# Keep only lower triangle
cor_df <- cor_df %>%
  filter(as.numeric(Var1) > as.numeric(Var2))

cor_phenovars<-ggplot(cor_df, aes(Var2, Var1, fill = value)) +
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

cor_phenovars

ggsave(
  filename = here::here("output", "figures", "correlogram_phenovars.png"),
  plot = cor_phenovars,
  width = 9,
  height = 6,
  dpi = 300
)


#### Test relationships between phenology variables and latitude/longitude ####

# Merge phenology estimates with transect coordinates

str(pheno_estimates_clean)
str(ebms_coord_df)

pheno_geo <- pheno_estimates_clean |>
  left_join(
    ebms_coord_df,
    by = c("SITE_ID" = "transect_id")
  )
head(pheno_geo)


pheno_sf <- pheno_geo %>%
  st_as_sf(coords = c("transect_lon", "transect_lat"),
           crs = 3035) %>%          # current CRS
  st_transform(crs = 4326)         # convert to lat–lon

# Extract lon/lat in degrees
coords <- st_coordinates(pheno_sf)

pheno_geo <- pheno_sf %>%
  mutate(
    lon = coords[,1],
    lat = coords[,2]
  ) %>%
  st_drop_geometry()


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
  select(lat, all_of(vars_pheno)) |>
  pivot_longer(
    cols = all_of(vars_pheno),
    names_to = "var",
    values_to = "value"
  )


europe <- ne_countries(
  scale = "medium",
  continent = "Europe",
  returnclass = "sf"
)

site_pheno <- pheno_geo %>%
  filter(!is.na(lon), !is.na(lat)) %>%
  group_by(SITE_ID, lon, lat) %>%
  summarise(
    ONSET_mean = mean(ONSET_mean, na.rm = TRUE),
    PEAKDAY = mean(PEAKDAY, na.rm = TRUE),
    OFFSET_mean = mean(OFFSET_mean, na.rm = TRUE),
    FLIGHT_LENGTH_mean = mean(FLIGHT_LENGTH_mean, na.rm = TRUE),
    .groups = "drop"
  )

# Convert sites to sf
sites_sf <- st_as_sf(site_pheno,
                     coords = c("lon", "lat"),
                     crs = 4326)

bbox_sites <- st_bbox(sites_sf)

buffer <- 1

bbox_expanded <- bbox_sites
bbox_expanded["xmin"] <- bbox_sites["xmin"] - buffer
bbox_expanded["xmax"] <- bbox_sites["xmax"] + buffer
bbox_expanded["ymin"] <- bbox_sites["ymin"] - buffer
bbox_expanded["ymax"] <- bbox_sites["ymax"] + buffer

europe_crop <- st_crop(europe, bbox_expanded)


ggplot() +
  geom_sf(data = europe_crop, fill = "grey95", color = "grey70") +
  geom_sf(data = sites_sf,
          aes(color = ONSET_mean),
          size = 2) +
  scale_color_viridis_c() +
  theme_minimal()
