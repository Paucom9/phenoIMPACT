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
library(tidyr) # For data tidying
library(broom) # For tidying model outputs
library(purrr) # For functional programming
library(nlme) # For generalized least squares models
library(here) # For file path management
library(reshape2)
library(ggplot2)
library(lme4)
library(broom.mixed)
library(emmeans)
library(ggeffects)

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

cor_mat <- pheno_trends_site |>
  select(SPECIES, SITE_ID, phenovar, estimate) |>
  pivot_wider(
    names_from = phenovar,
    values_from = estimate
  ) |>
  select(-SPECIES, -SITE_ID) |>
  cor(use = "complete.obs")

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


# --- Pheno trends vs. climatic regions --- #

ebms_clim_df   <- read.csv(here("data", "ebms_transect_climate.csv"), sep = ",", dec = ".")

pheno_clim_df <- pheno_trends_site |>
  left_join(ebms_clim_df, by = c("SITE_ID" = "transect_id"))

head(pheno_clim_df)


models_clim <- pheno_clim_df %>%
  split(.$phenovar) %>%
  map(~ lm(
    estimate ~ genzname,
    data = .x,
    weights = 1 / std.error^2
  ))


results_fixed <- map_df(
  models_clim,
  ~ tidy(.x, effects = "fixed"),
  .id = "phenovar"
)

print(results_fixed, n = Inf)

anova_results <- map_df(
  models_clim,
  ~ as.data.frame(anova(.x)),
  .id = "phenovar"
)

anova_results

# Plot estimated marginal means for Flight length mean (p = 0.09)

m_fl <- models_clim[["FLIGHT_LENGTH_mean"]]
emm_fl <- emmeans(m_fl, ~ genzname)
emm_df <- as.data.frame(emm_fl)

ggplot(emm_df, aes(x = genzname, y = emmean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = "Climatic region",
    y = "Trend in flight length (days/year)",
    title = "Flight length temporal trend across climatic regions"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


# --- Pheno trends vs. latitude --- #

ebms_coord_df  <- read.csv(here("data", "ebms_transect_coord.csv"), sep = ",", dec = ".")

pheno_coord_df <- pheno_trends_site |>
  left_join(ebms_coord_df, by = c("SITE_ID" = "transect_id"))

head(pheno_coord_df)
str(pheno_coord_df)

models_latitude <- pheno_coord_df %>%
  split(.$phenovar) %>%
  map(~ lm(
    estimate ~ scale(transect_lat) + scale(transect_lon) + scale(N_years),
    data = .x,
    weights = 1 / std.error^2
  ))

results_latitude <- map_df(
  models_latitude,
  broom::tidy,
  .id = "phenovar"
)

results_latitude

anova_results_lat <- map_df(
  models_latitude,
  ~ as.data.frame(anova(.x)),
  .id = "phenovar"
)

anova_results_lat


# Plot latitudinal effects

pred_latitude <- map_df(
  models_latitude,
  ~ as.data.frame(ggpredict(.x, terms = "transect_lat [all]")),
  .id = "phenovar"
)

pred_lon <- map_df(
  models_latitude,
  ~ as.data.frame(ggpredict(.x, terms = "transect_lon [all]")),
  .id = "phenovar"
)

ggplot(pred_latitude, aes(x = x, y = predicted)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ phenovar, scales = "free_y") +
  theme_bw() +
  labs(
    x = "Latitude (projected meters)",
    y = "Trend (days per year)",
    title = "Latitudinal variation in phenological trends"
  )

ggplot(pred_lon, aes(x = x, y = predicted)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ phenovar, scales = "free_y") +
  theme_bw() +
  labs(
    x = "Longitude (projected meters)",
    y = "Trend (days per year)",
    title = "Longitudinal variation in phenological trends"
  )


# Gam spatial approach

models_gam <- pheno_coord_df %>%
  split(.$phenovar) %>%
  map(~ gam(
    estimate ~ s(transect_lon, transect_lat, k = 20),
    data = .x,
    weights = 1 / std.error^2,
    method = "REML"
  ))

gam_results <- map_df(
  models_gam,
  ~ as.data.frame(summary(.x)$s.table),
  .id = "phenovar"
)

gam_results


# Europe map
europe_3035 <- ne_countries(
  continent = "Europe",
  scale = "medium",
  returnclass = "sf"
) %>%
  st_transform(3035)

plot_gam_surface <- function(pv) {
  
  mod <- models_gam[[pv]]
  
  dat <- pheno_coord_df %>%
    filter(phenovar == pv)
  
  # Prediction grid
  grid <- expand.grid(
    transect_lon = seq(min(dat$transect_lon),
                       max(dat$transect_lon),
                       length.out = 80),
    transect_lat = seq(min(dat$transect_lat),
                       max(dat$transect_lat),
                       length.out = 80)
  )
  
  grid$pred <- predict(mod, newdata = grid)
  
  # Partial residuals
  dat$partial_resid <- residuals(mod, type = "response")
  
  sites_sf <- st_as_sf(
    dat,
    coords = c("transect_lon", "transect_lat"),
    crs = 3035
  )
  
  # Bounding box
  expand_factor <- 300000
  
  bbox <- st_bbox(sites_sf)
  
  bbox_exp <- bbox
  bbox_exp["xmin"] <- bbox["xmin"] - expand_factor
  bbox_exp["xmax"] <- bbox["xmax"] + expand_factor
  bbox_exp["ymin"] <- bbox["ymin"] - expand_factor
  bbox_exp["ymax"] <- bbox["ymax"] + expand_factor
  
  europe_crop <- st_crop(europe_3035, bbox_exp)
  
  grid_crop <- grid %>%
    filter(
      between(transect_lon, bbox_exp["xmin"], bbox_exp["xmax"]),
      between(transect_lat, bbox_exp["ymin"], bbox_exp["ymax"])
    )
  
  # Plot
  ggplot() +
    geom_sf(data = europe_crop,
            fill = "grey95",
            color = "grey70") +
    geom_tile(data = grid_crop,
              aes(transect_lon, transect_lat, fill = pred),
              alpha = 0.85) +
    geom_sf(data = sites_sf,
            size = 0.8,
            alpha = 0.6,
            color = "black") +
    scale_fill_viridis_c(option = "C",
                         name = "Trend") +
    coord_sf(
      xlim = c(bbox_exp["xmin"], bbox_exp["xmax"]),
      ylim = c(bbox_exp["ymin"], bbox_exp["ymax"]),
      expand = FALSE
    ) +
    theme_void() +
    labs(title = pv)
}

plots <- lapply(names(models_gam), plot_gam_surface)

wrap_plots(plots, ncol = 3)
