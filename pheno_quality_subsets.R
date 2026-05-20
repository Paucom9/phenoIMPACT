# ============================================================================================ #
# pheno_quality_subsets.R
#
# Author: Pau Colom
# Date: 2026-05-20
#
# Description: Join phenology estimates with sampling-support variables and create 
# sensitivity subsets
# ============================================================================================ #

rm(list = ls())

#### Load required libraries ####

library(dplyr)
library(here)
library(writexl)

# ---------------------------------------------------------------------------- #
#### Load data ####
# ---------------------------------------------------------------------------- #

phenology_estimates <- read.csv(
  here("output", "pheno_abund_estimates_allspp.csv"),
  sep = ",",
  dec = "."
)

phenology_sampling_support <- read.csv(
  here("output", "phenology_sampling_support_allspp.csv"),
  sep = ",",
  dec = "."
)

# ---------------------------------------------------------------------------- #
#### Prepare keys ####
# ---------------------------------------------------------------------------- #

phenology_estimates <- phenology_estimates |>
  dplyr::mutate(
    ID = as.character(ID),
    YEAR = as.integer(YEAR),
    SPECIES = as.character(SPECIES),
    SITE_ID = as.character(SITE_ID)
  )

phenology_sampling_support <- phenology_sampling_support |>
  dplyr::mutate(
    ID = as.character(ID),
    YEAR = as.integer(YEAR),
    SPECIES = as.character(SPECIES),
    SITE_ID = as.character(SITE_ID)
  )

# ---------------------------------------------------------------------------- #
#### Join phenology estimates + sampling support ####
# ---------------------------------------------------------------------------- #

phenology_estimates_quality <- phenology_estimates |>
  dplyr::left_join(
    phenology_sampling_support,
    by = c("ID", "YEAR", "SPECIES", "SITE_ID")
  ) |>
  dplyr::mutate(
    
    # Boundary estimates outside the real sampling window
    onset_before_first_visit = dplyr::case_when(
      !is.na(ONSET_mean) & !is.na(first_visit_doy) ~
        ONSET_mean < first_visit_doy,
      TRUE ~ NA
    ),
    
    offset_after_last_visit = dplyr::case_when(
      !is.na(OFFSET_mean) & !is.na(last_visit_doy) ~
        OFFSET_mean > last_visit_doy,
      TRUE ~ NA
    ),
    
    # Gaps between estimated boundaries and observed occurrences
    gap_onset_first_obs = dplyr::case_when(
      !is.na(ONSET_mean) & !is.na(first_obs_doy) ~
        first_obs_doy - ONSET_mean,
      TRUE ~ NA_real_
    ),
    
    gap_last_obs_offset = dplyr::case_when(
      !is.na(OFFSET_mean) & !is.na(last_obs_doy) ~
        OFFSET_mean - last_obs_doy,
      TRUE ~ NA_real_
    ),
    
    # Relaxed quality criteria: at least one observed zero
    onset_supported_1zero =
      !is.na(ONSET_mean) &
      n_zero_visits_before_first_obs >= 1 &
      onset_before_first_visit == FALSE,
    
    offset_supported_1zero =
      !is.na(OFFSET_mean) &
      n_zero_visits_after_last_obs >= 1 &
      offset_after_last_visit == FALSE,
    
    boundaries_supported_1zero =
      onset_supported_1zero & offset_supported_1zero,
    
    # Strict quality criteria: at least two observed zeros
    onset_supported_2zero =
      !is.na(ONSET_mean) &
      n_zero_visits_before_first_obs >= 2 &
      onset_before_first_visit == FALSE,
    
    offset_supported_2zero =
      !is.na(OFFSET_mean) &
      n_zero_visits_after_last_obs >= 2 &
      offset_after_last_visit == FALSE,
    
    boundaries_supported_2zero =
      onset_supported_2zero & offset_supported_2zero
  )

# ---------------------------------------------------------------------------- #
#### Check join success ####
# ---------------------------------------------------------------------------- #

join_check <- phenology_estimates_quality |>
  dplyr::summarise(
    n_rows = dplyr::n(),
    n_with_sampling_support = sum(!is.na(n_visits_site_year)),
    prop_with_sampling_support = n_with_sampling_support / n_rows
  )

join_check

# ---------------------------------------------------------------------------- #
#### Quality summaries ####
# ---------------------------------------------------------------------------- #

quality_summary_overall <- phenology_estimates_quality |>
  dplyr::summarise(
    n_cases = dplyr::n(),
    
    prop_onset_before_first_visit =
      mean(onset_before_first_visit, na.rm = TRUE),
    
    prop_offset_after_last_visit =
      mean(offset_after_last_visit, na.rm = TRUE),
    
    mean_zero_visits_before_first_obs =
      mean(n_zero_visits_before_first_obs, na.rm = TRUE),
    
    mean_zero_visits_after_last_obs =
      mean(n_zero_visits_after_last_obs, na.rm = TRUE),
    
    median_gap_onset_first_obs =
      median(gap_onset_first_obs, na.rm = TRUE),
    
    median_gap_last_obs_offset =
      median(gap_last_obs_offset, na.rm = TRUE),
    
    prop_onset_supported_1zero =
      mean(onset_supported_1zero, na.rm = TRUE),
    
    prop_offset_supported_1zero =
      mean(offset_supported_1zero, na.rm = TRUE),
    
    prop_boundaries_supported_1zero =
      mean(boundaries_supported_1zero, na.rm = TRUE),
    
    prop_onset_supported_2zero =
      mean(onset_supported_2zero, na.rm = TRUE),
    
    prop_offset_supported_2zero =
      mean(offset_supported_2zero, na.rm = TRUE),
    
    prop_boundaries_supported_2zero =
      mean(boundaries_supported_2zero, na.rm = TRUE)
  )

quality_summary_by_bms <- phenology_estimates_quality |>
  dplyr::group_by(bms_id) |>
  dplyr::summarise(
    n_cases = dplyr::n(),
    prop_onset_supported_1zero = mean(onset_supported_1zero, na.rm = TRUE),
    prop_offset_supported_1zero = mean(offset_supported_1zero, na.rm = TRUE),
    prop_boundaries_supported_1zero = mean(boundaries_supported_1zero, na.rm = TRUE),
    prop_onset_supported_2zero = mean(onset_supported_2zero, na.rm = TRUE),
    prop_offset_supported_2zero = mean(offset_supported_2zero, na.rm = TRUE),
    prop_boundaries_supported_2zero = mean(boundaries_supported_2zero, na.rm = TRUE),
    median_gap_onset_first_obs = median(gap_onset_first_obs, na.rm = TRUE),
    median_gap_last_obs_offset = median(gap_last_obs_offset, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::arrange(desc(n_cases))

quality_summary_by_year <- phenology_estimates_quality |>
  dplyr::group_by(YEAR) |>
  dplyr::summarise(
    n_cases = dplyr::n(),
    prop_onset_supported_1zero = mean(onset_supported_1zero, na.rm = TRUE),
    prop_offset_supported_1zero = mean(offset_supported_1zero, na.rm = TRUE),
    prop_boundaries_supported_1zero = mean(boundaries_supported_1zero, na.rm = TRUE),
    .groups = "drop"
  )

quality_summary_overall
quality_summary_by_bms
quality_summary_by_year

# ---------------------------------------------------------------------------- #
#### Create sensitivity subsets ####
# ---------------------------------------------------------------------------- #

pheno_all <- phenology_estimates_quality

pheno_onset_1zero <- phenology_estimates_quality |>
  dplyr::filter(onset_supported_1zero)

pheno_onset_2zero <- phenology_estimates_quality |>
  dplyr::filter(onset_supported_2zero)

pheno_offset_1zero <- phenology_estimates_quality |>
  dplyr::filter(offset_supported_1zero)

pheno_offset_2zero <- phenology_estimates_quality |>
  dplyr::filter(offset_supported_2zero)

pheno_boundaries_1zero <- phenology_estimates_quality |>
  dplyr::filter(boundaries_supported_1zero)

pheno_boundaries_2zero <- phenology_estimates_quality |>
  dplyr::filter(boundaries_supported_2zero)

subset_summary <- tibble::tibble(
  subset = c(
    "all",
    "onset_1zero",
    "onset_2zero",
    "offset_1zero",
    "offset_2zero",
    "boundaries_1zero",
    "boundaries_2zero"
  ),
  n_rows = c(
    nrow(pheno_all),
    nrow(pheno_onset_1zero),
    nrow(pheno_onset_2zero),
    nrow(pheno_offset_1zero),
    nrow(pheno_offset_2zero),
    nrow(pheno_boundaries_1zero),
    nrow(pheno_boundaries_2zero)
  )
) |>
  dplyr::mutate(
    prop_retained = n_rows / nrow(pheno_all)
  )

subset_summary

# ---------------------------------------------------------------------------- #
#### Save outputs ####
# ---------------------------------------------------------------------------- #

out_dir <- here("output", "phenology_quality_subsets")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(
  phenology_estimates_quality,
  file = file.path(out_dir, "pheno_abund_estimates_allspp_with_quality.csv"),
  row.names = FALSE
)

write.csv(
  pheno_onset_1zero,
  file = file.path(out_dir, "pheno_subset_onset_1zero.csv"),
  row.names = FALSE
)

write.csv(
  pheno_onset_2zero,
  file = file.path(out_dir, "pheno_subset_onset_2zero.csv"),
  row.names = FALSE
)

write.csv(
  pheno_offset_1zero,
  file = file.path(out_dir, "pheno_subset_offset_1zero.csv"),
  row.names = FALSE
)

write.csv(
  pheno_offset_2zero,
  file = file.path(out_dir, "pheno_subset_offset_2zero.csv"),
  row.names = FALSE
)

write.csv(
  pheno_boundaries_1zero,
  file = file.path(out_dir, "pheno_subset_boundaries_1zero.csv"),
  row.names = FALSE
)

write.csv(
  pheno_boundaries_2zero,
  file = file.path(out_dir, "pheno_subset_boundaries_2zero.csv"),
  row.names = FALSE
)

writexl::write_xlsx(
  list(
    join_check = join_check,
    subset_summary = subset_summary,
    quality_summary_overall = quality_summary_overall,
    quality_summary_by_bms = quality_summary_by_bms,
    quality_summary_by_year = quality_summary_by_year
  ),
  file.path(out_dir, "phenology_quality_subset_summaries.xlsx")
)
