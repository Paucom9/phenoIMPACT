# ============================================================================================ #
# extract_pheno_plasticity.R
#
# Author: Pau Colom
# Date: 2026-05-12
#
# Description:
# Estimate population-level phenological plasticity using partial pooling.
# Plasticity is estimated as the population-level slope of phenological timing
# against thermal anomalies.
# ============================================================================================ #

rm(list = ls())

#### Load required libraries ####
library(data.table)
library(lme4)
library(lmerTest)
library(dplyr)
library(tidyr)
library(tibble)
library(performance)
library(here)
library(sf)
library(writexl)

#### Data Import ####

phenology_estimates <- read.csv(
  here("output", "pheno_estimates_allspp.csv"),
  sep = ",",
  dec = "."
)

clim_vars <- read.csv(
  here("output", "climate", "climate_variables_all_phenophases.csv"),
  sep = ",",
  dec = "."
)

ebms_transect_coord <- read.csv(
  here("data", "ebms_transect_coord.csv"),
  sep = ",",
  dec = "."
)

voltinism <- read.csv(
  here("data", "voltinism", "species_country_voltinism.csv"),
  sep = ";",
  dec = "."
)

#### Data preparation ####

# Site coordinates + BMS identity
coords_sf <- ebms_transect_coord |>
  dplyr::filter(
    !is.na(transect_lon),
    !is.na(transect_lat)
  ) |>
  dplyr::distinct(
    bms_id,
    transect_id,
    transect_lon,
    transect_lat
  ) |>
  sf::st_as_sf(
    coords = c("transect_lon", "transect_lat"),
    crs = 3035
  ) |>
  sf::st_transform(4326)

coord_site <- coords_sf |>
  dplyr::mutate(
    latitude = sf::st_coordinates(coords_sf)[, 2]
  ) |>
  sf::st_drop_geometry() |>
  dplyr::select(
    SITE_ID = transect_id,
    bms_id,
    latitude
  )

# Add latitude and BMS identity to climate variables
clim_vars <- clim_vars |>
  dplyr::select(-dplyr::any_of("bms_id")) |>
  dplyr::left_join(coord_site, by = "SITE_ID") |>
  dplyr::mutate(
    latitude = as.numeric(scale(latitude))
  )

# Voltinism in long format
voltinism_long <- voltinism |>
  tidyr::pivot_longer(
    cols = -SPECIES,
    names_to = "bms_id",
    values_to = "voltinism"
  ) |>
  dplyr::mutate(
    bms_id = dplyr::recode(
      bms_id,
      "ES.CTBMS" = "ES-CTBMS",
      "ES.ZEBMS" = "ES-ZEBMS"
    )
  )

# Merge phenology + climate + voltinism
df <- phenology_estimates |>
  dplyr::left_join(
    clim_vars,
    by = c("SPECIES", "SITE_ID", "YEAR")
  ) |>
  dplyr::left_join(
    voltinism_long,
    by = c("SPECIES", "bms_id")
  ) |>
  dplyr::mutate(
    SPECIES = factor(SPECIES),
    SITE_ID = factor(SITE_ID),
    YEAR = as.integer(YEAR),
    bms_id = factor(bms_id),
    voltinism = factor(voltinism),
    pop_id = paste(SPECIES, SITE_ID, sep = "_")
  )

#### Function to estimate and extract contextual population-level plasticity ####

fit_extract_contextual_plasticity <- function(data,
                                              response_var,
                                              anomaly_var,
                                              moderators,
                                              plasticity_name,
                                              control_vars = NULL,
                                              invert_slope = FALSE,
                                              include_pop_slope = FALSE) {
  
  message("\n============================================================")
  message("Fitting contextual plasticity model for: ", plasticity_name)
  message("============================================================")
  
  required_vars <- c(
    response_var,
    anomaly_var,
    moderators,
    control_vars,
    "SPECIES",
    "SITE_ID",
    "YEAR",
    "pop_id",
    "voltinism"
  )
  
  d <- data |>
    dplyr::filter(
      dplyr::if_all(
        dplyr::all_of(required_vars),
        ~ !is.na(.)
      )
    ) |>
    dplyr::mutate(
      SPECIES = factor(SPECIES),
      SITE_ID = factor(SITE_ID),
      YEAR = factor(YEAR),
      pop_id = factor(pop_id),
      voltinism = factor(voltinism)
    )
  
  message("Number of observations: ", nrow(d))
  message("Number of species: ", dplyr::n_distinct(d$SPECIES))
  message("Number of sites: ", dplyr::n_distinct(d$SITE_ID))
  message("Number of populations: ", dplyr::n_distinct(d$pop_id))
  message("Number of years: ", dplyr::n_distinct(d$YEAR))
  
  if (nrow(d) == 0) {
    stop("No data left after filtering complete cases.")
  }
  
  control_formula <- if (!is.null(control_vars)) {
    paste0(paste(control_vars, collapse = " + "), " + ")
  } else {
    ""
  }
  
  moderator_formula <- paste(moderators, collapse = " + ")
  
  random_formula <- if (include_pop_slope) {
    paste0(
      " + (1 | SITE_ID)",
      " + (1 + ", anomaly_var, " | SPECIES)",
      " + (1 + ", anomaly_var, " || pop_id)"
    )
  } else {
    paste0(
      " + (1 | SITE_ID)",
      " + (1 + ", anomaly_var, " | SPECIES)"
    )
  }
  
  form <- as.formula(
    paste0(
      response_var, " ~ ",
      control_formula,
      anomaly_var, " * (", moderator_formula, ")",
      random_formula
    )
  )
  
  print(form)
  
  mod <- lmer(
    form,
    data = d,
    REML = FALSE,
    control = lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 2e6)
    )
  )
  
  print(summary(mod))
  
  singular <- isSingular(mod)
  convergence <- performance::check_convergence(mod)
  
  message("Singular fit: ", singular)
  print(convergence)
  
  gc()
  
  b <- fixef(mod)
  re <- ranef(mod, condVar = FALSE)
  
  message("Random-effect components:")
  print(names(re))
  
  # Helper to get interaction name safely
  get_interaction <- function(beta, x1, x2) {
    
    term1 <- paste0(x1, ":", x2)
    term2 <- paste0(x2, ":", x1)
    
    if (term1 %in% names(beta)) {
      return(as.numeric(beta[term1]))
    }
    
    if (term2 %in% names(beta)) {
      return(as.numeric(beta[term2]))
    }
    
    stop("Interaction not found: ", x1, " x ", x2)
  }
  
  fixed_slope <- as.numeric(b[anomaly_var])
  
  interaction_betas <- purrr::map_dbl(
    moderators,
    ~ get_interaction(b, anomaly_var, .x)
  )
  
  names(interaction_betas) <- moderators
  
  message("Fixed anomaly slope: ", fixed_slope)
  message("Interaction slopes:")
  print(interaction_betas)
  
  # Species-level random slope
  sp_re <- as.data.frame(re$SPECIES)
  
  if (!anomaly_var %in% colnames(sp_re)) {
    stop("Species-level random slope not found for ", anomaly_var)
  }
  
  sp_slopes <- data.frame(
    SPECIES = rownames(sp_re),
    sp_slope_dev = as.numeric(sp_re[[anomaly_var]]),
    stringsAsFactors = FALSE
  )
  
  message("SD species slope dev: ", sd(sp_slopes$sp_slope_dev, na.rm = TRUE))
  
  # Optional population-level random slope
  if (include_pop_slope) {
    
    pop_re_name <- names(re)[
      grepl("^pop_id", names(re)) &
        vapply(
          re,
          function(x) anomaly_var %in% colnames(as.data.frame(x)),
          logical(1)
        )
    ][1]
    
    if (is.na(pop_re_name)) {
      stop("Population-level random slope not found for ", anomaly_var)
    }
    
    pop_re <- as.data.frame(re[[pop_re_name]])
    
    pop_slopes <- data.frame(
      pop_id = rownames(pop_re),
      pop_slope_dev = as.numeric(pop_re[[anomaly_var]]),
      stringsAsFactors = FALSE
    )
    
    message("SD population slope dev: ", sd(pop_slopes$pop_slope_dev, na.rm = TRUE))
    
  } else {
    
    pop_slopes <- d |>
      dplyr::mutate(pop_id = as.character(pop_id)) |>
      dplyr::distinct(pop_id) |>
      dplyr::mutate(pop_slope_dev = 0)
  }
  
  # Mean moderator values per population
  pop_context <- d |>
    dplyr::mutate(
      SPECIES = as.character(SPECIES),
      SITE_ID = as.character(SITE_ID),
      pop_id = as.character(pop_id),
      voltinism = as.character(voltinism)
    ) |>
    dplyr::group_by(pop_id, SPECIES, SITE_ID, voltinism) |>
    dplyr::summarise(
      n_years = dplyr::n_distinct(YEAR),
      dplyr::across(
        dplyr::all_of(moderators),
        ~ mean(.x, na.rm = TRUE)
      ),
      .groups = "drop"
    )
  
  plast <- pop_context |>
    dplyr::left_join(sp_slopes, by = "SPECIES") |>
    dplyr::left_join(pop_slopes, by = "pop_id")
  
  if (any(is.na(plast$sp_slope_dev))) {
    stop("Some species did not match species-level random slopes.")
  }
  
  if (any(is.na(plast$pop_slope_dev))) {
    stop("Some populations did not match population-level random slopes.")
  }
  
  # Contextual slope
  plast$plasticity_raw <- fixed_slope +
    plast$sp_slope_dev +
    plast$pop_slope_dev
  
  for (m in moderators) {
    plast$plasticity_raw <- plast$plasticity_raw +
      interaction_betas[m] * plast[[m]]
  }
  
  if (invert_slope) {
    plast$plasticity_final <- -plast$plasticity_raw
  } else {
    plast$plasticity_final <- plast$plasticity_raw
  }
  
  message("Summary raw contextual plasticity:")
  print(summary(plast$plasticity_raw))
  
  message("Summary final contextual plasticity:")
  print(summary(plast$plasticity_final))
  
  message("SD final plasticity: ", sd(plast$plasticity_final, na.rm = TRUE))
  
  if (sd(plast$plasticity_final, na.rm = TRUE) < 1e-8) {
    stop("Extracted contextual plasticity is constant. Something went wrong.")
  }
  
  plast <- plast |>
    dplyr::rename(
      !!paste0(plasticity_name, "_raw") := plasticity_raw,
      !!plasticity_name := plasticity_final
    )
  
  return(
    list(
      model = mod,
      plasticity = plast,
      diagnostics = list(
        singular = singular,
        convergence = convergence,
        n_obs = nrow(d),
        n_species = dplyr::n_distinct(d$SPECIES),
        n_pop = dplyr::n_distinct(d$pop_id),
        n_years = dplyr::n_distinct(d$YEAR),
        formula = form,
        fixed_slope = fixed_slope,
        interaction_betas = interaction_betas,
        sd_species_slope_dev = sd(sp_slopes$sp_slope_dev, na.rm = TRUE),
        include_pop_slope = include_pop_slope
      )
    )
  )
}

#### 1. Onset contextual plasticity ####

df_onset <- df |>
  dplyr::filter(pheno_type == "ONSET_mean")

onset_res <- fit_extract_contextual_plasticity(
  data = df_onset,
  response_var = "ONSET_mean",
  anomaly_var = "clim_anomaly_tw60",
  moderators = c(
    "photo_tw60",
    "clim_background_tw60",
    "clim_predictability_tw60",
    "clim_autocorr_tw60"
  ),
  plasticity_name = "onset_advancement_plasticity_contextual",
  invert_slope = TRUE,
  include_pop_slope = FALSE
)

onset_plasticity <- onset_res$plasticity

summary(onset_plasticity$onset_advancement_plasticity_contextual)

#### 2. Offset contextual plasticity: univoltine ####

df_offset_uni <- df |>
  dplyr::filter(
    pheno_type == "OFFSET_mean",
    voltinism == "univoltine"
  )

offset_uni_res <- fit_extract_contextual_plasticity(
  data = df_offset_uni,
  response_var = "OFFSET_mean",
  anomaly_var = "clim_anomaly_tw90",
  moderators = c(
    "photo_tw90",
    "clim_background_tw90",
    "clim_predictability_tw90",
    "clim_autocorr_tw90"
  ),
  control_vars = "ONSET_mean",
  plasticity_name = "offset_univoltine_termination_plasticity_contextual",
  invert_slope = FALSE,
  include_pop_slope = FALSE
)

offset_uni_plasticity <- offset_uni_res$plasticity

summary(offset_uni_plasticity$offset_univoltine_termination_plasticity_contextual)

#### 3. Offset contextual plasticity: multivoltine ####

df_offset_multi <- df |>
  dplyr::filter(
    pheno_type == "OFFSET_mean",
    voltinism == "multivoltine"
  )

offset_multi_res <- fit_extract_contextual_plasticity(
  data = df_offset_multi,
  response_var = "OFFSET_mean",
  anomaly_var = "clim_anomaly_tw90",
  moderators = c(
    "photo_tw90",
    "clim_background_tw90",
    "clim_predictability_tw90",
    "clim_autocorr_tw90"
  ),
  control_vars = "ONSET_mean",
  plasticity_name = "offset_multivoltine_delay_plasticity_contextual",
  invert_slope = FALSE,
  include_pop_slope = FALSE
)

offset_multi_plasticity <- offset_multi_res$plasticity

summary(offset_multi_plasticity$offset_multivoltine_delay_plasticity_contextual)

#### 4. Combine contextual plasticity outputs ####

phenological_plasticity_contextual <- onset_plasticity |>
  dplyr::select(
    pop_id,
    SPECIES,
    SITE_ID,
    voltinism,
    n_years_onset = n_years,
    onset_advancement_plasticity_contextual_raw,
    onset_advancement_plasticity_contextual
  ) |>
  dplyr::left_join(
    offset_uni_plasticity |>
      dplyr::select(
        pop_id,
        n_years_offset_uni = n_years,
        offset_univoltine_termination_plasticity_contextual_raw,
        offset_univoltine_termination_plasticity_contextual
      ),
    by = "pop_id"
  ) |>
  dplyr::left_join(
    offset_multi_plasticity |>
      dplyr::select(
        pop_id,
        n_years_offset_multi = n_years,
        offset_multivoltine_delay_plasticity_contextual_raw,
        offset_multivoltine_delay_plasticity_contextual
      ),
    by = "pop_id"
  ) |>
  dplyr::mutate(
    offset_termination_plasticity_contextual = dplyr::case_when(
      voltinism == "univoltine" ~ offset_univoltine_termination_plasticity_contextual,
      voltinism == "multivoltine" ~ offset_multivoltine_delay_plasticity_contextual,
      TRUE ~ NA_real_
    ),
    offset_termination_plasticity_contextual_raw = dplyr::case_when(
      voltinism == "univoltine" ~ offset_univoltine_termination_plasticity_contextual_raw,
      voltinism == "multivoltine" ~ offset_multivoltine_delay_plasticity_contextual_raw,
      TRUE ~ NA_real_
    )
  )

summary(phenological_plasticity_contextual$onset_advancement_plasticity_contextual)
summary(phenological_plasticity_contextual$offset_termination_plasticity_contextual)

#### 5. Save outputs ####

write.csv(
  phenological_plasticity_contextual,
  here("output", "phenological_population_plasticity_contextual.csv"),
  row.names = FALSE
)

writexl::write_xlsx(
  phenological_plasticity_contextual,
  here("output", "phenological_population_plasticity_contextual.xlsx")
)

saveRDS(
  onset_res$model,
  here("output", "model_onset_contextual_plasticity.rds")
)

saveRDS(
  offset_uni_res$model,
  here("output", "model_offset_univoltine_contextual_plasticity.rds")
)

saveRDS(
  offset_multi_res$model,
  here("output", "model_offset_multivoltine_contextual_plasticity.rds")
)

diagnostics_plasticity_contextual <- list(
  onset = onset_res$diagnostics,
  offset_univoltine = offset_uni_res$diagnostics,
  offset_multivoltine = offset_multi_res$diagnostics
)

saveRDS(
  diagnostics_plasticity_contextual,
  here("output", "diagnostics_phenological_population_plasticity_contextual.rds")
)