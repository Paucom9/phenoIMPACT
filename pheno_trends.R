# ============================================================================================ #
# 02_pheno_trends.R
#
# Author: Pau Colom
# Date: 2026-02-19
#
# Description: This script calculates phenological trends for each species × site combination 
# using linear models with an AR(1) correlation structure. It then visualizes the distribution 
# of these trends and explores correlations between different phenological variables. 
#
# ============================================================================================ #

#### Load required libraries ####
# ---
library(dplyr) # For data manipulation
library(tidyr) # For data tidying
library(broom) # For tidying model outputs
library(purrr) # For functional programming
library(nlme) # For generalized least squares models
library(here) # For file path management
library(reshape2)
library(ggplot2)
library(lme4)
library(glmmTMB)
library(broom.mixed)
library(emmeans)
library(ggeffects)
library(patchwork)
library(mgcv)
library(sf)
library(rnaturalearth)


# ----

# ---- Data Import and Preparation ---- #

here::here() # Check the current working directory

#-----

df  <- read.csv(here("output", "pheno_estimates_allspp.csv"), sep = ",", dec = ".")

str(df)

#-----


# ---- Calculate phenological trends ---- #

# Select variables
vars_pheno <- c(
  "ONSET_mean",
  "ONSET_var",
  "PEAKDAY",
  "OFFSET_mean",
  "OFFSET_var",
  "FLIGHT_LENGTH_mean",
  "FLIGHT_LENGTH_var"
)

# 1) number of years per species × site
n_years_df <- df |>
  mutate(YEAR_num = as.numeric(as.character(YEAR))) |>
  group_by(SPECIES, SITE_ID) |>
  summarise(N_years = n(), .groups = "drop")

# 2) trends for all phenology variables
pheno_trends_site <- map_dfr(vars_pheno, function(v) {
  
  df %>%
    mutate(YEAR_num = as.numeric(as.character(YEAR))) %>%
    group_by(SPECIES, SITE_ID) %>%
    filter(n() >= 10) %>%
    group_modify(~{
      
      model <- try(
        gls(
          reformulate("YEAR_num", v),
          data = .x,
          correlation = corAR1(form = ~ YEAR_num) # 
        ),
        silent = TRUE
      )
      
      if(inherits(model, "try-error")) {
        return(tibble())
      }
      
      coef_tab <- summary(model)$tTable
      
      if(!"YEAR_num" %in% rownames(coef_tab)) {
        return(tibble())
      }
      
      tibble(
        term      = "YEAR_num",
        estimate  = coef_tab["YEAR_num", "Value"],
        std.error = coef_tab["YEAR_num", "Std.Error"],
        statistic = coef_tab["YEAR_num", "t-value"],
        p.value   = coef_tab["YEAR_num", "p-value"]
      )
      
    }) %>%
    mutate(phenovar = v) %>%
    ungroup()
  
}) %>%
  left_join(n_years_df, by = c("SPECIES", "SITE_ID"))

head(pheno_trends_site)
str(pheno_trends_site)

# Save as CSV inside project
write.csv(
  pheno_trends_site,
  file = here("output", "pheno_temporal_trends_allspp.csv"),
  row.names = FALSE
)