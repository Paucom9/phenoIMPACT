# ============================================================================================ #
# population_growth_models.R
#
# Author: Pau Colom
# Date: 2026-05-12
#
# Description:
# Test whether phenological plasticity predicts annual population growth
# under thermal anomalies.
#
# Annual population growth is calculated as:
# lambda = log(N_t / N_t-1)
# ============================================================================================ #

rm(list = ls())

#### Load required libraries ####

library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)
library(performance)
library(ggeffects)
library(ggplot2)
library(here)
library(writexl)

#### Data import ####

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

phenological_plasticity <- read.csv(
  here("output", "phenological_population_plasticity_contextual.csv"),
  sep = ",",
  dec = "."
)

#### Helper functions ####

zscale <- function(x) {
  as.numeric(scale(x))
}

#### 1. Calculate annual population growth ####

df_lambda <- phenology_estimates |>
  dplyr::distinct(
    SPECIES,
    SITE_ID,
    YEAR,
    ABUND_INDEX
  ) |>
  dplyr::filter(
    !is.na(ABUND_INDEX),
    ABUND_INDEX > 0
  ) |>
  dplyr::mutate(
    YEAR = as.integer(YEAR),
    SPECIES = as.factor(SPECIES),
    SITE_ID = as.factor(SITE_ID),
    pop_id = paste(SPECIES, SITE_ID, sep = "_")
  ) |>
  dplyr::arrange(pop_id, YEAR) |>
  dplyr::group_by(pop_id) |>
  dplyr::mutate(
    log_N = log(ABUND_INDEX),
    log_N_lag = dplyr::lag(log_N),
    year_lag = dplyr::lag(YEAR),
    lambda = log_N - log_N_lag,
    consecutive_year = YEAR == year_lag + 1
  ) |>
  dplyr::ungroup() |>
  dplyr::filter(
    consecutive_year,
    !is.na(lambda),
    !is.na(log_N_lag)
  )

summary(df_lambda$lambda)

#### 2. Prepare climate data ####

# Climate variables for onset models
clim_onset <- clim_vars |>
  dplyr::filter(pheno_type == "ONSET_mean") |>
  dplyr::select(
    SPECIES,
    SITE_ID,
    YEAR,
    clim_anomaly_tw60
  ) |>
  dplyr::distinct()

# Climate variables for offset models
clim_offset <- clim_vars |>
  dplyr::filter(pheno_type == "OFFSET_mean") |>
  dplyr::select(
    SPECIES,
    SITE_ID,
    YEAR,
    clim_anomaly_tw90
  ) |>
  dplyr::distinct()

#### 3. Onset plasticity and annual population growth ####

df_demo_onset <- df_lambda |>
  dplyr::left_join(
    clim_onset,
    by = c("SPECIES", "SITE_ID", "YEAR")
  ) |>
  dplyr::left_join(
    phenological_plasticity |>
      dplyr::select(
        pop_id,
        voltinism,
        onset_advancement_plasticity_contextual
      ),
    by = "pop_id"
  ) |>
  dplyr::filter(
    !is.na(lambda),
    !is.na(log_N_lag),
    !is.na(clim_anomaly_tw60),
    !is.na(onset_advancement_plasticity_contextual),
    !is.na(voltinism)
  ) |>
  dplyr::mutate(
    SPECIES = as.factor(SPECIES),
    SITE_ID = as.factor(SITE_ID),
    YEAR = as.factor(YEAR),
    pop_id = as.factor(pop_id),
    voltinism = as.factor(voltinism),
    lambda_z = zscale(lambda),
    log_N_lag_z = zscale(log_N_lag),
    clim_anomaly_z = zscale(clim_anomaly_tw60),
    onset_plasticity_z = zscale(onset_advancement_plasticity_contextual)
  )

summary(df_demo_onset$lambda)
summary(df_demo_onset$onset_advancement_plasticity_contextual)

mod_lambda_onset <- lmer(
  lambda ~ log_N_lag_z +
    clim_anomaly_z * onset_plasticity_z +
    I(clim_anomaly_z^2) +
    (1 | SPECIES) +
    (1 | SITE_ID) +
    (1 | YEAR),
  data = df_demo_onset,
  REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e6)
  )
)

summary(mod_lambda_onset)
isSingular(mod_lambda_onset)
performance::check_convergence(mod_lambda_onset)

#### 4. Offset termination plasticity and annual population growth ####

df_demo_offset <- df_lambda |>
  dplyr::left_join(
    clim_offset,
    by = c("SPECIES", "SITE_ID", "YEAR")
  ) |>
  dplyr::left_join(
    phenological_plasticity |>
      dplyr::select(
        pop_id,
        voltinism,
        offset_termination_plasticity_contextual
      ),
    by = "pop_id"
  ) |>
  dplyr::filter(
    !is.na(lambda),
    !is.na(log_N_lag),
    !is.na(clim_anomaly_tw90),
    !is.na(offset_termination_plasticity_contextual),
    !is.na(voltinism)
  ) |>
  dplyr::mutate(
    SPECIES = as.factor(SPECIES),
    SITE_ID = as.factor(SITE_ID),
    YEAR = as.factor(YEAR),
    pop_id = as.factor(pop_id),
    voltinism = as.factor(voltinism),
    log_N_lag_z = zscale(log_N_lag),
    clim_anomaly_z = zscale(clim_anomaly_tw90),
    offset_plasticity_z = zscale(offset_termination_plasticity_contextual)
  )

summary(df_demo_offset$lambda)
summary(df_demo_offset$offset_termination_plasticity_contextual)

mod_lambda_offset <- lmer(
  lambda ~ log_N_lag_z +
    clim_anomaly_z * offset_plasticity_z +
    I(clim_anomaly_z^2) +
    (1 | SPECIES) +
    (1 | SITE_ID) +
    (1 | YEAR),
  data = df_demo_offset,
  REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e6)
  )
)

summary(mod_lambda_offset)
isSingular(mod_lambda_offset)
performance::check_convergence(mod_lambda_offset)

#### 5. Offset model with voltinism interaction ####

mod_lambda_offset_voltinism <- lmer(
  lambda ~ log_N_lag_z +
    clim_anomaly_z * offset_plasticity_z * voltinism +
    I(clim_anomaly_z^2) +
    (1 | SPECIES) +
    (1 | SITE_ID) +
    (1 | YEAR),
  data = df_demo_offset,
  REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e6)
  )
)

summary(mod_lambda_offset_voltinism)
isSingular(mod_lambda_offset_voltinism)
performance::check_convergence(mod_lambda_offset_voltinism)

#### 6. Multivoltine-only offset model ####

df_demo_offset_multi <- df_demo_offset |>
  dplyr::filter(voltinism == "multivoltine") |>
  dplyr::mutate(
    offset_delay_plasticity_z = zscale(offset_termination_plasticity_contextual)
  )

mod_lambda_offset_multi <- lmer(
  lambda ~ log_N_lag_z +
    clim_anomaly_z * offset_delay_plasticity_z +
    I(clim_anomaly_z^2) +
    (1 | SPECIES) +
    (1 | SITE_ID) +
    (1 | YEAR),
  data = df_demo_offset_multi,
  REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e6)
  )
)

summary(mod_lambda_offset_multi)
isSingular(mod_lambda_offset_multi)
performance::check_convergence(mod_lambda_offset_multi)

#### 7. Full model: onset + offset ####

df_demo_full <- df_lambda |>
  dplyr::left_join(
    clim_onset,
    by = c("SPECIES", "SITE_ID", "YEAR")
  ) |>
  dplyr::left_join(
    clim_offset,
    by = c("SPECIES", "SITE_ID", "YEAR")
  ) |>
  dplyr::left_join(
    phenological_plasticity |>
      dplyr::select(
        pop_id,
        voltinism,
        onset_advancement_plasticity_contextual,
        offset_termination_plasticity_contextual
      ),
    by = "pop_id"
  ) |>
  dplyr::filter(
    !is.na(lambda),
    !is.na(log_N_lag),
    !is.na(clim_anomaly_tw60),
    !is.na(clim_anomaly_tw90),
    !is.na(onset_advancement_plasticity_contextual),
    !is.na(offset_termination_plasticity_contextual),
    !is.na(voltinism)
  ) |>
  dplyr::mutate(
    SPECIES = as.factor(SPECIES),
    SITE_ID = as.factor(SITE_ID),
    YEAR = as.factor(YEAR),
    pop_id = as.factor(pop_id),
    voltinism = as.factor(voltinism),
    log_N_lag_z = zscale(log_N_lag),
    clim_anomaly_onset_z = zscale(clim_anomaly_tw60),
    clim_anomaly_offset_z = zscale(clim_anomaly_tw90),
    onset_plasticity_z = zscale(onset_advancement_plasticity_contextual),
    offset_plasticity_z = zscale(offset_termination_plasticity_contextual)
  )

mod_lambda_full <- lmer(
  lambda ~ log_N_lag_z +
    clim_anomaly_onset_z * onset_plasticity_z +
    clim_anomaly_offset_z * offset_plasticity_z +
    I(clim_anomaly_onset_z^2) +
    I(clim_anomaly_offset_z^2) +
    (1 | SPECIES) +
    (1 | SITE_ID) +
    (1 | YEAR),
  data = df_demo_full,
  REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e6)
  )
)

summary(mod_lambda_full)
isSingular(mod_lambda_full)
performance::check_convergence(mod_lambda_full)

#### 8. Save outputs ####

saveRDS(
  mod_lambda_onset,
  here("output", "model_lambda_onset_plasticity.rds")
)

saveRDS(
  mod_lambda_offset,
  here("output", "model_lambda_offset_plasticity.rds")
)

saveRDS(
  mod_lambda_offset_voltinism,
  here("output", "model_lambda_offset_plasticity_voltinism.rds")
)

saveRDS(
  mod_lambda_offset_multi,
  here("output", "model_lambda_offset_multivoltine.rds")
)

saveRDS(
  mod_lambda_full,
  here("output", "model_lambda_full_onset_offset_plasticity.rds")
)

write.csv(
  df_demo_full,
  here("output", "df_demographic_models_full.csv"),
  row.names = FALSE
)

writexl::write_xlsx(
  list(
    df_lambda = df_lambda,
    df_demo_onset = df_demo_onset,
    df_demo_offset = df_demo_offset,
    df_demo_full = df_demo_full
  ),
  here("output", "demographic_model_datasets.xlsx")
)