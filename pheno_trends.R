# ==================================================================================================================
# 02_pheno_trends.R
#
# Author: Pau Colom
# Date: 2026-02-19
#
# Description: This script calculates phenological trends for each species × site combination using linear models
# with an AR(1) correlation structure. It then visualizes the distribution of these trends and explores correlations
# between different phenological variables. 
#
# ==================================================================================================================

# Load required libraries
# ----
library(dplyr) # For data manipulation
library(broom) # For tidying model outputs
library(purrr) # For functional programming
library(nlme) # For generalized least squares models
library(here) # For file path management

# ----

# ---- Data Import and Preparation ---- #

here::here() # Check the current working directory

#-----

df  <- read.csv(here("output", "pheno_estimates.csv"), sep = ",", dec = ".")

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
  file = here("output", "pheno_temporal_trends.csv"),
  row.names = FALSE
)


# --- Plot histograms of trends for each phenology variable --- #
# compute median per phenovar
median_df <- pheno_trends_site |>
  mutate(trend_decade = estimate * 10) |>
  group_by(phenovar) |>
  summarise(
    med = median(trend_decade, na.rm = TRUE),
    .groups = "drop"
  )

pheno_trends_site |>
  mutate(trend_decade = estimate * 10) |>
  ggplot(aes(x = trend_decade)) +
  geom_histogram(bins = 30, fill = "grey70", color = "black") +
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "grey40") +
  geom_vline(
    data = median_df,
    aes(xintercept = med),
    linewidth = 0.8
  ) +
  facet_wrap(~ phenovar, scales = "free") +
  theme_minimal() +
  labs(
    x = "Trend (days per decade)",
    y = "Number of sites"
  )

# --- Correlation between phenovars --- #


library(reshape2)
library(ggplot2)

cor_df <- melt(cor_mat)

ggplot(cor_df, aes(Var1, Var2, fill = value)) +
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


pheno_trends_wide <- pheno_trends_site |>
  select(SPECIES, SITE_ID, phenovar, estimate) |>
  pivot_wider(
    names_from = phenovar,
    values_from = estimate
  ) |>
  mutate(across(starts_with("FLIGHT_LENGTH"), ~ .x * 10),
         across(c("ONSET_mean", "OFFSET_mean", "PEAKDAY"), ~ .x * 10))



plot_df <- pheno_trends_wide |>
  select(
    FLIGHT_LENGTH_mean,
    ONSET_mean,
    OFFSET_mean,
    PEAKDAY
  ) |>
  pivot_longer(
    cols = c(ONSET_mean, OFFSET_mean, PEAKDAY),
    names_to = "phenovar",
    values_to = "trend"
  )

# correlation per panel
cor_df <- plot_df |>
  group_by(phenovar) |>
  summarise(
    r = cor(trend, FLIGHT_LENGTH_mean, use = "complete.obs"),
    .groups = "drop"
  )

ggplot(plot_df, aes(trend, FLIGHT_LENGTH_mean)) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.9) +
  geom_text(
    data = cor_df,
    aes(
      x = -Inf, y = Inf,
      label = paste0("r = ", round(r, 2))
    ),
    hjust = -0.1, vjust = 1.2,
    inherit.aes = FALSE,
    size = 4
  ) +
  facet_wrap(~ phenovar, scales = "free_x") +
  theme_minimal() +
  labs(
    x = "Phenology trend (days per decade)",
    y = "Flight length trend (days per decade)"
  )

