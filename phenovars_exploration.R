# ============================================================================================ #
# phenovars_exploration.R
#
# Author: Pau Colom
# Date: 2026-02-20
#
# Desciription:
#
# ============================================================================================ #


#### Load required libraries ####
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
library(reshape2)     # For reshaping data frames
library(here)       # For constructing file paths in a way that is independent of the operating system


# ---- Data Import and Preparation ---- #

here::here() # Check the current working directory

#-----

phenology_estimates  <- read.csv(here("output", "pheno_estimates_allspp.csv"), sep = ",", dec = ".")

str(phenology_estimates)


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
