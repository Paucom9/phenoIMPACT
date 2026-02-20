# ============================================================================================ #
# pheno_plasticity.R
#
# Author: Pau Colom
# Date: 2026-02-20
#
# Desciription:
#
# ============================================================================================ #

#### Load required libraries ####
library(data.table)  # For efficient data handling
library(lme4)
library(dplyr)       # For data manipulation
library(tidyr)
library(purrr)
library(broom.mixed)
library(performance)
library(ggplot2)





# ---- Data Import and Preparation ---- #

here::here() # Check the current working directory

#---

phenology_estimates  <- read.csv(here("output", "pheno_estimates_allspp.csv"), sep = ",", dec = ".")
mean_temperature <- read.csv(here("output", "climate", "mean_annual_temperature.csv"), sep = ",", dec = ".")

str(phenology_estimates)
str(mean_temperature)

#### Calculate within-temperature ####
#(i.e., temporal anomaly around each transect’s long-term mean) ####

temp_partitioned <- mean_temperature |>
  group_by(transect_id) |>
  mutate(
    temp_between = mean(avg_temp, na.rm = TRUE),
    temp_sd      = sd(avg_temp, na.rm = TRUE),
    temp_within  = avg_temp - temp_between
  ) |>
  ungroup() |>
  mutate(
    temp_between_sc = scale(temp_between)[,1],
    temp_sd_sc      = scale(temp_sd)[,1],
    temp_within_sc  = scale(temp_within)[,1]   # optional but recommended
  )

#### Merge datasets ####
df <- phenology_estimates |>
  left_join(
    temp_partitioned,
    by = c("SITE_ID" = "transect_id",
           "YEAR" = "year")
  )
str(df)

df <- df |>
  mutate(
    SPECIES = factor(SPECIES),
    SITE_ID = factor(SITE_ID)
  )

df$sp_site <- interaction(df$SPECIES, df$SITE_ID, drop = TRUE)

#### Model plasticity ####

m_pop <- lmer(
  ONSET_mean ~ temp_within +
    (1 | SPECIES) +
    (1 | SPECIES:SITE_ID) +
    (0 + within_temp | SPECIES:SITE_ID),
  data = df,
  REML = FALSE
)

m_no_pop <- lmer(
  ONSET_mean ~ temp_within +
    (1 | SPECIES) +
    (1 | SPECIES:SITE_ID),
  data = df,
  REML = FALSE
)

m_plast_partition <- lmer(
  ONSET_mean ~ temp_within +
    (1 | SPECIES) +
    (1 | SPECIES:SITE_ID) +
    (0 + within_temp | SPECIES) +
    (0 + within_temp | SPECIES:SITE_ID),
  data = df,
  REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)

summary(m_no_pop)
summary(m_pop)

anova(m_no_pop, m_pop) # Adding population-level random slopes significantly improves the model.

summary(m_plast_partition) 

anova(m_pop, m_plast_partition)


#### Do climate explain variation in plasticity? ####

m_explain <- lmer(
  ONSET_mean ~ temp_within * temp_between +
    temp_within * temp_sd +
    (0 + temp_within | SPECIES) +
    (0 + temp_within | SPECIES:SITE_ID) +
    (1 | SPECIES) +
    (1 | SPECIES:SITE_ID),
  data = df,
  REML = FALSE
)
summary(m_explain)


m_explain_sc <- lmer(
  ONSET_mean ~ temp_within_sc * temp_between_sc +
    temp_within_sc * temp_sd_sc +
    (0 + temp_within_sc | SPECIES) +
    (0 + temp_within_sc | SPECIES:SITE_ID) +
    (1 | SPECIES) +
    (1 | SPECIES:SITE_ID),
  data = df,
  REML = FALSE
)

summary(m_explain_sc)
# Cold climates → stronger response to anomalies; Warm climates → buffered response
# Stable climates → strong tracking of anomalies; Variable climates → weaker sensitivity

# Check consistency among phenovars #



df_long <- df |>
  pivot_longer(
    cols = c(
      ONSET_mean,
      ONSET_var,
      PEAKDAY,
      OFFSET_mean,
      OFFSET_var,
      FLIGHT_LENGTH_mean,
      FLIGHT_LENGTH_var
    ),
    names_to = "phenovar",
    values_to = "pheno_value"
  )

str(df_long)

models_pheno <- df_long %>%
  split(df_long$phenovar) %>%
  map(~ lmer(
    pheno_value ~ temp_within_sc * temp_between_sc +
      temp_within_sc * temp_sd_sc +
      (0 + temp_within_sc | SPECIES) +
      (0 + temp_within_sc | SPECIES:SITE_ID) +
      (1 | SPECIES) +
      (1 | SPECIES:SITE_ID),
    data = .x,
    REML = FALSE,
    control = lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 2e5)
    )
  ))


results <- map_df(
  models_pheno,
  ~ tidy(.x, effects = "fixed"),
  .id = "phenovar"
)


print(results, n = Inf)


# Check multicollinearity #

m_vif <- lmer(
  pheno_value ~ temp_within_sc * temp_between_sc +
    temp_within_sc * temp_sd_sc +
    (1 | SPECIES) +
    (1 | SPECIES:SITE_ID),
  data = df_long[df_long$phenovar == "ONSET_mean",],
  REML = FALSE
)

performance::check_collinearity(m_vif)


# Plot effects #


plot_df <- results %>%
  filter(term %in% c(
    "temp_between_sc",
    "temp_within_sc",
    "temp_within_sc:temp_between_sc",
    "temp_within_sc:temp_sd_sc"
  )) %>%
  mutate(
    lower = estimate - 1.96 * std.error,
    upper = estimate + 1.96 * std.error
  )


ggplot(plot_df,
       aes(x = phenovar, y = estimate, ymin = lower, ymax = upper)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ term, scales = "free_y") +
  coord_flip() +
  theme_bw()
