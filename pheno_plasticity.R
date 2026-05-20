# ============================================================================================ #
# pheno_plasticity.R
#
# Author: Pau Colom
# Date: 2026-02-20
#
# Desciription: This script analyzes phenological plasticity in response to temperature anomalies, 
# and how this plasticity is influenced by photoperiod and other climatic variables. 
#
# ============================================================================================ #

# Clean session
rm(list = ls())

#### Load required libraries ####
# ---
library(data.table)  # For efficient data handling
library(lme4)
library(dplyr)       # For data manipulation
library(tidyr)
library(purrr)
library(broom.mixed)
library(performance)
library(ggplot2)
library(sf)
library(stringr)
library(ggeffects)
library(Matrix)
library(lmerTest)
library(extrafont)
library(writexl)
library(MuMIn)
library(here)


# ---

#### Data Import and Preparation ####
# ---

here::here() # Check the current working directory

phenology_estimates <- read.csv(
  here::here("output", "pheno_estimates_allspp.csv"),
  sep = ",",
  dec = "."
)

phenology_sampling_support <- read.csv(
  here::here("output", "phenology_sampling_support_allspp.csv"),
  sep = ",",
  dec = "."
)

clim_vars <- read.csv(
  here::here("output", "climate", "climate_variables_all_phenophases.csv"),
  sep = ",",
  dec = "."
)

ebms_transect_coord <- read.csv(
  here::here("data", "ebms_transect_coord.csv"),
  sep = ",",
  dec = "."
)

voltinism <- read.csv(
  here::here("data", "voltinism", "species_country_voltinism.csv"),
  sep = ";",
  dec = "."
)

# ---------------------------------------------------------------------------- #
#### Prepare phenology estimates + sampling support ####
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
  ) |>
  # Drop bms_id here to avoid duplicated bms_id columns after joining climate data
  dplyr::select(-dplyr::any_of("bms_id"))

phenology_estimates <- phenology_estimates |>
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

# Check that sampling support joined correctly
sampling_join_check <- phenology_estimates |>
  dplyr::summarise(
    n_rows = dplyr::n(),
    n_with_sampling_support = sum(!is.na(n_visits_site_year)),
    prop_with_sampling_support = n_with_sampling_support / n_rows,
    prop_onset_supported_1zero = mean(onset_supported_1zero, na.rm = TRUE),
    prop_offset_supported_1zero = mean(offset_supported_1zero, na.rm = TRUE),
    prop_boundaries_supported_1zero = mean(boundaries_supported_1zero, na.rm = TRUE)
  )

sampling_join_check

# ---------------------------------------------------------------------------- #
#### Merge datasets ####
# ---------------------------------------------------------------------------- #

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

# Merge phenology + sampling support + climate + voltinism
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
    bms_id = factor(bms_id),
    voltinism = factor(voltinism),
    sp_site = interaction(SPECIES, SITE_ID, drop = TRUE)
  )

# Optional check
df |>
  dplyr::summarise(
    n_rows = dplyr::n(),
    prop_onset_supported_1zero = mean(onset_supported_1zero, na.rm = TRUE),
    prop_offset_supported_1zero = mean(offset_supported_1zero, na.rm = TRUE),
    prop_boundaries_supported_1zero = mean(boundaries_supported_1zero, na.rm = TRUE)
  )

# ------------------------------------------------------------------------------------- #
#### Data exploration #### 

# Inspect climatic anomalies # 

df_unique <- df %>%
  group_by(SITE_ID, SPECIES, YEAR) %>%
  slice(1) %>%
  ungroup()

med <- median(df_unique$clim_anomaly_tw30, na.rm = TRUE)

p <- ggplot(df_unique, aes(x = clim_anomaly_tw30)) +
  geom_histogram(fill = "grey", color = "white", bins = 30) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1) +
  geom_vline(xintercept = med, color = "red", linetype = "dashed", size = 1) +
  annotate("text",
           x = med,
           y = Inf,
           label = paste("Median =", round(med, 2)),
           vjust = 2,
           hjust = 0,
           color = "red") +
  labs(
    x = "Temperature anomaly (°C)",
    y = "Count"
  ) +
  theme_minimal()

ggsave(
  filename = here::here("output", "figures", "distribution_clim_anomalies_27_4_26.png"),
  plot = p,
  width = 8,
  height = 6
)

# Inspect distribution of predictor variables #

# variables to plot
vars_plot <- c(
  "clim_background_tw60",
  "clim_trend_tw60",
  "clim_autocorr_tw60",
  "clim_stability_tw60",
  "clim_predictability_tw60",
  "photo_tw60"
)

# remove pseudo-replication (one row per site × year × pheno_type)
df_plot <- df |>
  dplyr::filter(pheno_type == "ONSET_mean") |>
  dplyr::distinct(
    SITE_ID,
    dplyr::across(dplyr::all_of(vars_plot))
  )

# long format
df_long <- df_plot |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(vars_plot),
    names_to = "variable",
    values_to = "value"
  )

# labels
labels <- c(
  clim_background_tw90      = "Background temperature\n(°C)",
  clim_trend_tw90           = "Temperature trend\n(°C per decade)",
  clim_autocorr_tw90        = "Autocorrelation\n(lag-1)",
  clim_stability_tw90       = "Stability\n(-SD)",
  clim_predictability_tw90  = "Predictability\n(-Var residuals)",
  photo_tw90                = "Photoperiod\n(daily hours)",
  latitude                  = "Latitude\n(degrees)"
)

# apply labels safely
df_long$variable <- factor(
  df_long$variable,
  levels = names(labels),
  labels = labels
)

# plot
p <- ggplot(df_long, aes(x = value)) +
  geom_histogram(bins = 30, fill = "grey70", color = "black") +
  facet_wrap(~ variable, scales = "free", ncol = 3) +
  theme_classic(base_size = 12) +
  labs(x = NULL, y = "Frequency") +
  theme(
    strip.text = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
  )

# save
ggsave(
  filename = here::here("output", "figures", "predictor_histograms_26_3_26B.png"),
  plot = p,
  width = 9,
  height = 7,
  dpi = 300
)

# ---

# Correlation among predictors # 

cor_mat <- df_plot |>
  select(clim_background_tw60,
         clim_trend_tw60,
         clim_autocorr_tw60,
         clim_predictability_tw60,
         photo_tw60) |>
  cor(use = "complete.obs")


new_names <- c(
  "Mean temp.",
  "Temp. trend",
  "Temp. autocorr.",
  "Temp. stab.",
  "Photoperiod"
)

rownames(cor_mat) <- new_names
colnames(cor_mat) <- new_names

cor_df <- as.data.frame(as.table(cor_mat))

# keep only lower triangle
cor_df <- cor_df |>
  filter(as.numeric(Var1) > as.numeric(Var2))

corr_clim_pred_sds <- ggplot(cor_df, aes(Var1, Var2, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = round(Freq, 2)), size = 4) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, limits = c(-1,1)) +
  theme_minimal() +
  labs(x = NULL, y = NULL, fill = "r") +
  coord_equal()

corr_clim_pred_sds

# Save plot
ggsave(
  filename = here::here("output", "figures", "corr_site_vars_4_5_26.png"),
  plot = corr_clim_pred_sds,
  width = 6,
  height = 5,
  dpi = 300
)


#### Functions for phenological plasticity analysis ####

# ---------------------------------------------------------------------------- #
# Helper functions
# ---------------------------------------------------------------------------- #

get_vars_window <- function(window, include_autocorr = FALSE) {
  
  vars <- c(
    paste0("clim_anomaly_tw", window),
    paste0("photo_tw", window),
    paste0("clim_background_tw", window),
    paste0("clim_predictability_tw", window),
    paste0("clim_trend_tw", window)
  )
  
  if (include_autocorr) {
    vars <- c(vars, paste0("clim_autocorr_tw", window))
  }
  
  vars
}

get_moderators_window <- function(window, include_autocorr = FALSE) {
  
  mods <- c(
    paste0("photo_tw", window),
    paste0("clim_background_tw", window),
    paste0("clim_predictability_tw", window),
    paste0("clim_trend_tw", window)
  )
  
  if (include_autocorr) {
    mods <- c(mods, paste0("clim_autocorr_tw", window))
  }
  
  mods
}

fit_lmer_model <- function(formula, data) {
  
  lmer(
    formula,
    data = data,
    REML = FALSE,
    control = lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 2e6)
    )
  )
}

make_complete_data <- function(data, response, windows,
                               include_autocorr = FALSE) {
  
  vars_needed <- unique(c(
    response,
    "SITE_ID",
    "SPECIES",
    unlist(lapply(windows, get_vars_window, include_autocorr = include_autocorr))
  ))
  
  vars_needed <- vars_needed[vars_needed %in% names(data)]
  
  data[stats::complete.cases(data[, vars_needed]), ]
}

make_manual_aic_table <- function(models) {
  
  tibble(
    window = as.numeric(names(models)),
    logLik = map_dbl(models, ~ as.numeric(logLik(.x))),
    df = map_dbl(models, ~ attr(logLik(.x), "df")),
    nobs = map_int(models, nobs)
  ) |>
    mutate(
      AIC = -2 * logLik + 2 * df,
      BIC = -2 * logLik + log(nobs) * df,
      delta_AIC = AIC - min(AIC),
      weight = exp(-0.5 * delta_AIC) / sum(exp(-0.5 * delta_AIC))
    ) |>
    arrange(AIC)
}

make_plasticity_formula <- function(response, window, random_cor = FALSE) {
  
  anomaly <- paste0("clim_anomaly_tw", window)
  
  random_species <- if (random_cor) {
    paste0("(1 + ", anomaly, " | SPECIES)")
  } else {
    paste0("(1 + ", anomaly, " || SPECIES)")
  }
  
  as.formula(paste0(
    response, " ~ ",
    anomaly, " + ",
    "(1 | SITE_ID) + ",
    random_species
  ))
}

make_full_formula <- function(response, window,
                              include_autocorr = FALSE,
                              random_cor = FALSE) {
  
  anomaly <- paste0("clim_anomaly_tw", window)
  moderators <- get_moderators_window(window, include_autocorr)
  
  random_species <- if (random_cor) {
    paste0("(1 + ", anomaly, " | SPECIES)")
  } else {
    paste0("(1 + ", anomaly, " || SPECIES)")
  }
  
  as.formula(paste0(
    response, " ~ ",
    anomaly, " * (", paste(moderators, collapse = " + "), ") + ",
    "(1 | SITE_ID) + ",
    random_species
  ))
}

find_interaction_name <- function(coef_names, anomaly, moderator) {
  
  candidates <- c(
    paste0(anomaly, ":", moderator),
    paste0(moderator, ":", anomaly)
  )
  
  out <- candidates[candidates %in% coef_names][1]
  
  if (is.na(out)) {
    stop(paste("Interaction not found:", anomaly, "x", moderator))
  }
  
  out
}

get_fixed_effects_table <- function(model) {
  
  coefs <- summary(model)$coefficients |>
    as.data.frame() |>
    tibble::rownames_to_column("term")
  
  names(coefs) <- gsub("Std. Error", "SE", names(coefs), fixed = TRUE)
  names(coefs) <- gsub("Pr\\(>\\|t\\|\\)", "p", names(coefs))
  
  coefs |>
    mutate(
      lower_95 = Estimate - 1.96 * SE,
      upper_95 = Estimate + 1.96 * SE
    )
}

# ---------------------------------------------------------------------------- #
# Window selection and full models
# ---------------------------------------------------------------------------- #

fit_window_selection_models <- function(data, response,
                                        windows = c(30, 60, 90),
                                        random_cor = FALSE) {
  
  d_tw <- make_complete_data(
    data = data,
    response = response,
    windows = windows,
    include_autocorr = FALSE
  )
  
  models <- setNames(
    lapply(windows, function(w) {
      fit_lmer_model(
        make_plasticity_formula(
          response = response,
          window = w,
          random_cor = random_cor
        ),
        d_tw
      )
    }),
    as.character(windows)
  )
  
  list(
    data = d_tw,
    models = models,
    table = make_manual_aic_table(models)
  )
}

fit_full_models_all_windows <- function(data, response,
                                        windows = c(30, 60, 90),
                                        random_cor = FALSE) {
  
  d_tw <- make_complete_data(
    data = data,
    response = response,
    windows = windows,
    include_autocorr = FALSE
  )
  
  models <- setNames(
    lapply(windows, function(w) {
      fit_lmer_model(
        make_full_formula(
          response = response,
          window = w,
          include_autocorr = FALSE,
          random_cor = random_cor
        ),
        d_tw
      )
    }),
    as.character(windows)
  )
  
  list(
    data = d_tw,
    models = models,
    table = make_manual_aic_table(models)
  )
}

# ---------------------------------------------------------------------------- #
# LRTs and interaction tables
# ---------------------------------------------------------------------------- #

get_lrt_for_interaction <- function(full_model, moderator, window, data, response,
                                    random_cor = FALSE) {
  
  anomaly <- paste0("clim_anomaly_tw", window)
  
  all_moderators <- c(
    paste0("photo_tw", window),
    paste0("clim_background_tw", window),
    paste0("clim_predictability_tw", window),
    paste0("clim_trend_tw", window)
  )
  
  reduced_moderators <- setdiff(all_moderators, moderator)
  
  random_species <- if (random_cor) {
    paste0("(1 + ", anomaly, " | SPECIES)")
  } else {
    paste0("(1 + ", anomaly, " || SPECIES)")
  }
  
  reduced_formula <- as.formula(paste0(
    response, " ~ ",
    anomaly, " * (", paste(reduced_moderators, collapse = " + "), ") + ",
    anomaly, " + ", moderator, " + ",
    "(1 | SITE_ID) + ",
    random_species
  ))
  
  reduced_model <- fit_lmer_model(reduced_formula, data)
  
  lrt <- anova(reduced_model, full_model)
  
  logLik_reduced <- as.numeric(logLik(reduced_model))
  logLik_full <- as.numeric(logLik(full_model))
  df_reduced <- attr(logLik(reduced_model), "df")
  df_full <- attr(logLik(full_model), "df")
  
  AIC_reduced <- -2 * logLik_reduced + 2 * df_reduced
  AIC_full <- -2 * logLik_full + 2 * df_full
  
  BIC_reduced <- -2 * logLik_reduced + log(nobs(reduced_model)) * df_reduced
  BIC_full <- -2 * logLik_full + log(nobs(full_model)) * df_full
  
  data.frame(
    moderator_var = moderator,
    interaction = paste0(anomaly, ":", moderator),
    Chisq = lrt$Chisq[2],
    LRT_df = lrt$Df[2],
    p_LRT = lrt$`Pr(>Chisq)`[2],
    AIC_reduced = AIC_reduced,
    AIC_full = AIC_full,
    delta_AIC = AIC_reduced - AIC_full,
    BIC_reduced = BIC_reduced,
    BIC_full = BIC_full,
    delta_BIC = BIC_reduced - BIC_full
  )
}

get_main_interaction_table <- function(model, window, data, response,
                                       random_cor = FALSE) {
  
  anomaly <- paste0("clim_anomaly_tw", window)
  
  moderator_lookup <- tibble(
    moderator_var = c(
      paste0("photo_tw", window),
      paste0("clim_background_tw", window),
      paste0("clim_predictability_tw", window),
      paste0("clim_trend_tw", window)
    ),
    moderator = c(
      "Photoperiod",
      "Mean temperature",
      "Temperature predictability",
      "Temperature trend"
    ),
    interpretation = c(
      "Stronger advancement under longer photoperiods",
      "Weaker advancement in warmer baseline climates",
      "Stronger advancement in more predictable thermal environments",
      "Stronger advancement in faster-warming environments"
    )
  )
  
  fixed_table <- get_fixed_effects_table(model)
  
  lrt_table <- map_dfr(
    moderator_lookup$moderator_var,
    ~ get_lrt_for_interaction(
      full_model = model,
      moderator = .x,
      window = window,
      data = data,
      response = response,
      random_cor = random_cor
    )
  )
  
  moderator_lookup |>
    rowwise() |>
    mutate(
      term = find_interaction_name(
        coef_names = fixed_table$term,
        anomaly = anomaly,
        moderator = moderator_var
      )
    ) |>
    ungroup() |>
    left_join(fixed_table, by = "term") |>
    left_join(lrt_table, by = "moderator_var") |>
    select(
      all_of(c(
        "moderator",
        "term",
        "Estimate",
        "SE",
        "lower_95",
        "upper_95",
        "Chisq",
        "LRT_df",
        "p_LRT",
        "delta_AIC",
        "delta_BIC",
        "interpretation"
      ))
    )
}

# ---------------------------------------------------------------------------- #
# Consistency table across full 30/60/90 models - simplified wide version
# ---------------------------------------------------------------------------- #

get_consistency_table_windows <- function(full_models, response) {
  
  consistency_long <- purrr::map_dfr(names(full_models), function(w) {
    
    model <- full_models[[w]]
    window <- as.numeric(w)
    fixed_table <- get_fixed_effects_table(model)
    anomaly <- paste0("clim_anomaly_tw", window)
    
    moderators <- tibble::tibble(
      moderator_var = c(
        paste0("photo_tw", window),
        paste0("clim_background_tw", window),
        paste0("clim_predictability_tw", window),
        paste0("clim_trend_tw", window)
      ),
      moderator = c(
        "Photoperiod",
        "Mean temperature",
        "Temperature predictability",
        "Temperature trend"
      )
    )
    
    moderators |>
      dplyr::rowwise() |>
      dplyr::mutate(
        term = find_interaction_name(
          coef_names = fixed_table$term,
          anomaly = anomaly,
          moderator = moderator_var
        )
      ) |>
      dplyr::ungroup() |>
      dplyr::left_join(fixed_table, by = "term") |>
      dplyr::mutate(
        response = response,
        window = window,
        sig = dplyr::case_when(
          p < 0.001 ~ "***",
          p < 0.01  ~ "**",
          p < 0.05  ~ "*",
          TRUE ~ ""
        ),
        estimate_ci = paste0(
          round(Estimate, 2),
          " [",
          round(lower_95, 2),
          ", ",
          round(upper_95, 2),
          "]",
          sig
        )
      ) |>
      dplyr::select(
        response,
        window,
        moderator,
        Estimate,
        SE,
        lower_95,
        upper_95,
        p,
        sig,
        estimate_ci
      )
  })
  
  consistency_long |>
    dplyr::select(window, moderator, estimate_ci) |>
    tidyr::pivot_wider(
      names_from = moderator,
      values_from = estimate_ci
    ) |>
    dplyr::arrange(window)
}
# ---------------------------------------------------------------------------- #
# Plot functions
# ---------------------------------------------------------------------------- #

get_slope_df <- function(model, var_name, var_seq, label, anomaly_name) {
  
  b <- lme4::fixef(model)
  V <- as.matrix(vcov(model))
  
  coef_main <- anomaly_name
  coef_int <- paste0(coef_main, ":", var_name)
  
  if (!coef_int %in% names(b)) {
    coef_int <- paste0(var_name, ":", coef_main)
  }
  
  coef_names <- c(coef_main, coef_int)
  
  get_slope <- function(x) {
    
    g <- c(1, x)
    
    slope <- sum(b[coef_names] * g)
    
    se <- as.numeric(
      sqrt(t(g) %*% V[coef_names, coef_names] %*% g)
    )
    
    c(slope = slope, se = se)
  }
  
  out <- do.call(rbind, lapply(var_seq, get_slope))
  
  data.frame(
    x = var_seq,
    slope = out[, "slope"],
    lower = out[, "slope"] - 1.96 * out[, "se"],
    upper = out[, "slope"] + 1.96 * out[, "se"],
    variable = label
  )
}

make_plasticity_plot <- function(model, window,
                                 phenology_label = "phenology",
                                 ylim = NULL,
                                 facet = FALSE) {
  
  anomaly <- paste0("clim_anomaly_tw", window)
  seq_x <- seq(-2, 2, length.out = 100)
  
  df_plot <- bind_rows(
    get_slope_df(model, paste0("photo_tw", window), seq_x, "Photoperiod", anomaly),
    get_slope_df(model, paste0("clim_background_tw", window), seq_x, "Mean temperature", anomaly),
    get_slope_df(model, paste0("clim_predictability_tw", window), seq_x, "Temperature predictability", anomaly),
    get_slope_df(model, paste0("clim_trend_tw", window), seq_x, "Temperature trend", anomaly)
  )
  
  df_plot$variable <- factor(
    df_plot$variable,
    levels = c(
      "Photoperiod",
      "Mean temperature",
      "Temperature predictability",
      "Temperature trend"
    )
  )
  
  p <- ggplot(
    df_plot,
    aes(x = x, y = slope, color = variable, fill = variable)
  ) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.8, alpha = 0.75) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    scale_color_manual(values = c(
      "Photoperiod" = "#1b9e77",
      "Mean temperature" = "#E6AB02",
      "Temperature predictability" = "#7570b3",
      "Temperature trend" = "#d95f02"
    )) +
    scale_fill_manual(values = c(
      "Photoperiod" = "#1b9e77",
      "Mean temperature" = "#E6AB02",
      "Temperature predictability" = "#7570b3",
      "Temperature trend" = "#d95f02"
    )) +
    theme_classic(base_family = "Garamond", base_size = 16) +
    labs(
      x = "Environmental gradient",
      y = paste0("Slope of ", phenology_label, " vs. temperature anomaly"),
      color = "",
      fill = ""
    )
  
  if (!is.null(ylim)) {
    p <- p + coord_cartesian(ylim = ylim)
  }
  
  if (facet) {
    p <- p +
      facet_wrap(~ variable, nrow = 1) +
      guides(color = "none", fill = "none")
  }
  
  p
}

# ---------------------------------------------------------------------------- #
# Save and return figures for all full models
# ---------------------------------------------------------------------------- #

save_full_model_figures_all_windows <- function(full_models,
                                                output_prefix,
                                                out_dir,
                                                phenology_label = "phenology",
                                                plot_ylim = NULL,
                                                facet_plot = FALSE) {
  
  plots <- purrr::imap(full_models, function(mod, w) {
    
    p <- make_plasticity_plot(
      model = mod,
      window = as.numeric(w),
      phenology_label = phenology_label,
      ylim = plot_ylim,
      facet = facet_plot
    )
    
    ggsave(
      filename = file.path(
        out_dir,
        "figures",
        paste0(output_prefix, "_full_model_tw", w, "_plasticity_interactions.png")
      ),
      plot = p,
      width = 7,
      height = 5,
      dpi = 300
    )
    
    p
  })
  
  plots
}

# ---------------------------------------------------------------------------- #
# Alternative predictor models
# ---------------------------------------------------------------------------- #

safe_max_vif <- function(model) {
  
  out <- tryCatch(
    performance::check_collinearity(model),
    error = function(e) NULL
  )
  
  if (is.null(out)) return(NA_real_)
  
  suppressWarnings(max(out$VIF, na.rm = TRUE))
}

make_alt_aic_table <- function(models) {
  
  tibble(
    model = names(models),
    logLik = map_dbl(models, ~ as.numeric(logLik(.x))),
    df = map_dbl(models, ~ attr(logLik(.x), "df")),
    nobs = map_int(models, nobs),
    max_VIF = map_dbl(models, safe_max_vif)
  ) |>
    mutate(
      AIC = -2 * logLik + 2 * df,
      BIC = -2 * logLik + log(nobs) * df,
      delta_AIC = AIC - min(AIC),
      weight = exp(-0.5 * delta_AIC) / sum(exp(-0.5 * delta_AIC))
    ) |>
    arrange(AIC)
}

fit_alternative_models <- function(data, response, window,
                                   random_cor = FALSE) {
  
  anomaly <- paste0("clim_anomaly_tw", window)
  
  build_custom_formula <- function(moderators) {
    
    random_species <- if (random_cor) {
      paste0("(1 + ", anomaly, " | SPECIES)")
    } else {
      paste0("(1 + ", anomaly, " || SPECIES)")
    }
    
    as.formula(paste0(
      response, " ~ ",
      anomaly, " * (", paste(moderators, collapse = " + "), ") + ",
      "(1 | SITE_ID) + ",
      random_species
    ))
  }
  
  mods_main <- c(
    paste0("photo_tw", window),
    paste0("clim_background_tw", window),
    paste0("clim_predictability_tw", window),
    paste0("clim_trend_tw", window)
  )
  
  mods_plus_autocorr <- c(
    mods_main,
    paste0("clim_autocorr_tw", window)
  )
  
  mods_autocorr_instead <- c(
    paste0("photo_tw", window),
    paste0("clim_background_tw", window),
    paste0("clim_autocorr_tw", window),
    paste0("clim_trend_tw", window)
  )
  
  mods_no_trend <- c(
    paste0("photo_tw", window),
    paste0("clim_background_tw", window),
    paste0("clim_predictability_tw", window)
  )
  
  models <- list(
    main = fit_lmer_model(build_custom_formula(mods_main), data),
    main_plus_autocorr = fit_lmer_model(build_custom_formula(mods_plus_autocorr), data),
    autocorr_instead_predictability = fit_lmer_model(build_custom_formula(mods_autocorr_instead), data),
    no_trend = fit_lmer_model(build_custom_formula(mods_no_trend), data)
  )
  
  list(
    models = models,
    table = make_alt_aic_table(models)
  )
}

# ---------------------------------------------------------------------------- #
# Main wrapper
# ---------------------------------------------------------------------------- #

run_plasticity_protocol <- function(data,
                                    response,
                                    output_prefix,
                                    windows = c(30, 60, 90),
                                    out_dir = here("output", "phenology_plasticity"),
                                    random_cor = FALSE,
                                    phenology_label = "phenology",
                                    plot_ylim = NULL,
                                    facet_plot = FALSE) {
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
  
  window_selection <- fit_window_selection_models(
    data = data,
    response = response,
    windows = windows,
    random_cor = random_cor
  )
  
  tw_selection_table <- window_selection$table
  best_window <- tw_selection_table$window[1]
  
  full_window_models <- fit_full_models_all_windows(
    data = data,
    response = response,
    windows = windows,
    random_cor = random_cor
  )
  
  full_models_by_window <- full_window_models$models
  full_window_table <- full_window_models$table
  
  mod_main <- full_models_by_window[[as.character(best_window)]]
  d_main <- full_window_models$data
  
  fixed_table <- get_fixed_effects_table(mod_main)
  
  interaction_table <- get_main_interaction_table(
    model = mod_main,
    window = best_window,
    data = d_main,
    response = response,
    random_cor = random_cor
  )
  
  consistency_windows <- get_consistency_table_windows(
    full_models = full_models_by_window,
    response = response
  )
  
  plasticity_plot <- make_plasticity_plot(
    model = mod_main,
    window = best_window,
    phenology_label = phenology_label,
    ylim = plot_ylim,
    facet = facet_plot
  )
  
  ggsave(
    filename = file.path(
      out_dir,
      "figures",
      paste0(output_prefix, "_main_tw", best_window, "_plasticity_interactions.png")
    ),
    plot = plasticity_plot,
    width = 7,
    height = 5,
    dpi = 300
  )
  
  plots_all_windows <- save_full_model_figures_all_windows(
    full_models = full_models_by_window,
    output_prefix = output_prefix,
    out_dir = out_dir,
    phenology_label = phenology_label,
    plot_ylim = plot_ylim,
    facet_plot = facet_plot
  )
  
  vars_alt <- c(
    response,
    "SITE_ID",
    "SPECIES",
    get_vars_window(best_window, include_autocorr = TRUE)
  )
  
  vars_alt <- vars_alt[vars_alt %in% names(data)]
  d_alt <- data[stats::complete.cases(data[, vars_alt]), ]
  
  alt_models <- fit_alternative_models(
    data = d_alt,
    response = response,
    window = best_window,
    random_cor = random_cor
  )
  
  write_xlsx(
    list(
      window_selection_anomaly_only = tw_selection_table,
      full_model_AIC_all_windows = full_window_table,
      main_interactions = interaction_table,
      main_fixed_effects = fixed_table,
      consistency_across_windows = consistency_windows,
      alternative_predictor_sets = alt_models$table
    ),
    path = file.path(
      out_dir,
      "tables",
      paste0(output_prefix, "_plasticity_protocol_tables.xlsx")
    )
  )
  
  list(
    best_window = best_window,
    window_selection_models = window_selection$models,
    window_selection_table = tw_selection_table,
    full_models_by_window = full_models_by_window,
    full_window_table = full_window_table,
    main_model = mod_main,
    main_interaction_table = interaction_table,
    main_fixed_effects = fixed_table,
    consistency_windows = consistency_windows,
    alternative_models = alt_models,
    plot = plasticity_plot,
    plots_all_windows = plots_all_windows
  )
}

# ---------------------------------------------------------------------------- #
#### Onset mean ####
# ---------------------------------------------------------------------------- #

df_onset_mean <- df |>
  filter(pheno_type == "ONSET_mean")

res_onset_mean <- run_plasticity_protocol(
  data = df_onset_mean,
  response = "ONSET_mean",
  output_prefix = "onset_mean",
  windows = c(30, 60, 90),
  out_dir = here("output", "phenology_plasticity"),
  random_cor = FALSE,
  phenology_label = "onset",
  plot_ylim = NULL,
  facet_plot = FALSE
)

res_onset_mean$best_window
res_onset_mean$window_selection_table
res_onset_mean$full_window_table
res_onset_mean$main_interaction_table
res_onset_mean$consistency_windows
summary(res_onset_mean$main_model)
performance::check_collinearity(res_onset_mean$main_model)
res_onset_mean$plots_all_windows[["30"]]
res_onset_mean$plots_all_windows[["60"]]
res_onset_mean$plots_all_windows[["90"]]

# ---------------------------------------------------------------------------- #
#### First peak ####
# ---------------------------------------------------------------------------- #

df_first_peak <- df |>
  dplyr::filter(pheno_type == "FIRST_PEAK")

res_first_peak <- run_plasticity_protocol(
  data = df_first_peak,
  response = "FIRST_PEAK",
  output_prefix = "first_peak",
  windows = c(30, 60, 90),
  out_dir = here("output", "phenology_plasticity"),
  random_cor = FALSE,
  phenology_label = "first peak",
  plot_ylim = NULL,
  facet_plot = FALSE
)

res_first_peak$best_window
res_first_peak$window_selection_table
res_first_peak$main_interaction_table
res_first_peak$consistency_windows
res_first_peak$plots_all_windows[["30"]]
res_first_peak$plots_all_windows[["60"]]
res_first_peak$plots_all_windows[["90"]]


# ---------------------------------------------------------------------------- #
#### Offset mean - univoltine ####
# ---------------------------------------------------------------------------- #

df_offset_uni <- df |>
  dplyr::filter(
    pheno_type == "OFFSET_mean",
    voltinism == "univoltine"
  )

res_offset_uni <- run_plasticity_protocol(
  data = df_offset_uni,
  response = "OFFSET_mean",
  output_prefix = "offset_mean_univoltine",
  windows = c(30, 60, 90),
  out_dir = here("output", "phenology_plasticity"),
  random_cor = FALSE,
  phenology_label = "offset",
  plot_ylim = NULL,
  facet_plot = FALSE
)

res_offset_uni$best_window
res_offset_uni$window_selection_table
res_offset_uni$main_interaction_table
res_offset_uni$consistency_windows
res_offset_uni$plots_all_windows[["30"]]
res_offset_uni$plots_all_windows[["60"]]
res_offset_uni$plots_all_windows[["90"]]


# ---------------------------------------------------------------------------- #
#### Offset mean - multivoltine ####
# ---------------------------------------------------------------------------- #

df_offset_multi <- df |>
  dplyr::filter(
    pheno_type == "OFFSET_mean",
    voltinism == "multivoltine"
  )

res_offset_multi <- run_plasticity_protocol(
  data = df_offset_multi,
  response = "OFFSET_mean",
  output_prefix = "offset_mean_multivoltine",
  windows = c(30, 60, 90),
  out_dir = here("output", "phenology_plasticity"),
  random_cor = FALSE,
  phenology_label = "offset",
  plot_ylim = NULL,
  facet_plot = FALSE
)

res_offset_multi$best_window
res_offset_multi$window_selection_table
res_offset_multi$main_interaction_table
res_offset_multi$consistency_windows
res_offset_multi$plots_all_windows[["30"]]
res_offset_multi$plots_all_windows[["60"]]
res_offset_multi$plots_all_windows[["90"]]


# ---------------------------------------------------------------------------- #
#### Sensitivity analysis: sampling support 1 zero ####
# ---------------------------------------------------------------------------- #

out_dir_sens_1zero <- here("output", "phenology_plasticity_sensitivity_1zero")


# ---------------------------------------------------------------------------- #
#### Onset mean - 1 zero before first observation ####
# ---------------------------------------------------------------------------- #

df_onset_mean_1zero <- df |>
  dplyr::filter(
    pheno_type == "ONSET_mean",
    onset_supported_1zero
  )

res_onset_mean_1zero <- run_plasticity_protocol(
  data = df_onset_mean_1zero,
  response = "ONSET_mean",
  output_prefix = "onset_mean_1zero",
  windows = c(30, 60, 90),
  out_dir = out_dir_sens_1zero,
  random_cor = FALSE,
  phenology_label = "onset",
  plot_ylim = NULL,
  facet_plot = FALSE
)

res_onset_mean_1zero$best_window
res_onset_mean_1zero$window_selection_table
res_onset_mean_1zero$main_interaction_table
res_onset_mean_1zero$consistency_windows
res_onset_mean_1zero$plots_all_windows[["30"]]
res_onset_mean_1zero$plots_all_windows[["60"]]
res_onset_mean_1zero$plots_all_windows[["90"]]

# ---------------------------------------------------------------------------- #
#### Offset mean - univoltine - 1 zero after last observation ####
# ---------------------------------------------------------------------------- #

df_offset_uni_1zero <- df |>
  dplyr::filter(
    pheno_type == "OFFSET_mean",
    voltinism == "univoltine",
    offset_supported_1zero
  )

res_offset_uni_1zero <- run_plasticity_protocol(
  data = df_offset_uni_1zero,
  response = "OFFSET_mean",
  output_prefix = "offset_mean_univoltine_1zero",
  windows = c(30, 60, 90),
  out_dir = out_dir_sens_1zero,
  random_cor = FALSE,
  phenology_label = "offset",
  plot_ylim = NULL,
  facet_plot = FALSE
)

res_offset_uni_1zero$best_window
res_offset_uni_1zero$window_selection_table
res_offset_uni_1zero$main_interaction_table
res_offset_uni_1zero$consistency_windows
res_offset_uni_1zero$plots_all_windows[["30"]]
res_offset_uni_1zero$plots_all_windows[["60"]]
res_offset_uni_1zero$plots_all_windows[["90"]]


# ---------------------------------------------------------------------------- #
#### Offset mean - multivoltine - 1 zero after last observation ####
# ---------------------------------------------------------------------------- #

df_offset_multi_1zero <- df |>
  dplyr::filter(
    pheno_type == "OFFSET_mean",
    voltinism == "multivoltine",
    offset_supported_1zero
  )

res_offset_multi_1zero <- run_plasticity_protocol(
  data = df_offset_multi_1zero,
  response = "OFFSET_mean",
  output_prefix = "offset_mean_multivoltine_1zero",
  windows = c(30, 60, 90),
  out_dir = out_dir_sens_1zero,
  random_cor = FALSE,
  phenology_label = "offset",
  plot_ylim = NULL,
  facet_plot = FALSE
)

res_offset_multi_1zero$best_window
res_offset_multi_1zero$window_selection_table
res_offset_multi_1zero$main_interaction_table
res_offset_multi_1zero$consistency_windows
res_offset_multi_1zero$plots_all_windows[["30"]]
res_offset_multi_1zero$plots_all_windows[["60"]]
res_offset_multi_1zero$plots_all_windows[["90"]]
