# ============================================================================================ #
# pheno_trends.R
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
# ---



#### Data Import and Preparation ####
# ---
here::here() # Check the current working directory
# ---
df  <- read.csv(here("output", "pheno_estimates_allspp.csv"), sep = ",", dec = ".")
str(df)
# ---


#### Calculate phenological trends ####
# --- Select variables
vars_pheno <- c(
  "ONSET_mean",
  "ONSET_var",
  "PEAKDAY",
  "OFFSET_mean",
  "OFFSET_var",
  "FLIGHT_LENGTH_mean",
  "FLIGHT_LENGTH_var"
)

# ---  number of years per species × site
n_years_df <- df |>
  mutate(YEAR_num = as.numeric(as.character(YEAR))) |>
  group_by(SPECIES, SITE_ID) |>
  summarise(N_years = n(), .groups = "drop")

# --- trends for all phenology variables
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

# --- Save as CSV inside project
write.csv(
  pheno_trends_site,
  file = here("output", "pheno_temporal_trends_allspp.csv"),
  row.names = FALSE
)


#### Data exploration ####

#### Plot histograms of trends for each phenology variable ####
# Define order
vars_order <- c(
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

# Prepare data
pheno_trends_plot <- pheno_trends_site |>
  mutate(
    trend_decade = estimate,
    phenovar = factor(phenovar, levels = vars_order),
    phenovar_label = factor(pretty_labels[phenovar],
                            levels = pretty_labels[vars_order])
  )

# Compute medians
median_df <- pheno_trends_plot |>
  group_by(phenovar_label) |>
  summarise(
    med = median(trend_decade, na.rm = TRUE),
    .groups = "drop"
  )

# Plot
distr_pheno_trends <-ggplot(pheno_trends_plot, aes(x = trend_decade)) +
  geom_histogram(bins = 30, fill = "grey70", color = "black") +
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "grey40") +
  geom_vline(
    data = median_df,
    aes(xintercept = med),
    linewidth = 0.8
  ) +
  facet_wrap(~ phenovar_label, scales = "free") +
  theme_minimal() +
  labs(
    x = "Trend (days per year)",
    y = "Number of species-sites"
  )

distr_pheno_trends

ggsave(
  filename = here::here("output", "figures", "distribution_pheno_trends.png"),
  plot = distr_pheno_trends,
  width = 7,
  height = 5,
  dpi = 300
)

#### Correlation between phenovar trends ####

# Wide format + convert to days per decade
pheno_trends_wide <- pheno_trends_site |>
  select(SPECIES, SITE_ID, phenovar, estimate) |>
  pivot_wider(
    names_from = phenovar,
    values_from = estimate
  ) |>
  mutate(across(all_of(vars_order), ~ .x * 10))

# Correlation matrix
cor_mat <- pheno_trends_wide |>
  select(all_of(vars_order)) |>
  cor(use = "complete.obs")

# Reorder + rename
cor_mat <- cor_mat[vars_order, vars_order]
colnames(cor_mat) <- pretty_labels[vars_order]
rownames(cor_mat) <- pretty_labels[vars_order]

cor_df <- melt(cor_mat)

# Keep lower triangle only
cor_df <- cor_df |>
  filter(as.numeric(Var1) > as.numeric(Var2))

# Plot
correologram_pheno_trends<- ggplot(cor_df, aes(Var2, Var1, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", value)), size = 3) +
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

correologram_pheno_trends

ggsave(
  filename = here::here("output", "figures", "correlogram_phenotrends.png"),
  plot = correologram_pheno_trends,
  width = 9,
  height = 6,
  dpi = 300
)