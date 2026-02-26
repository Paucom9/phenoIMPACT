# ============================================================================================ #
# pheno_trend_models.R
#
# Author: Pau Colom
# Date: 2026-02-20
#
# Description: This script analyzes phenological trends across different sites and species 
# by fitting linear mixed models. It explores the distribution of these trends, their 
# correlations, and how they relate to climatic regions, latitude, and temperature. 
# The script also includes visualizations to illustrate the findings.
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
pheno_trends_site  <- read.csv(here("output", "pheno_temporal_trends_allspp.csv"), sep = ",", dec = ".")
str(pheno_trends_site)
# ---

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



#### Pheno trends vs. climatic regions ####
# Climate region data
ebms_clim_df   <- read.csv(here("data", "ebms_transect_climate.csv"), sep = ",", dec = ".")

pheno_clim_df <- pheno_trends_site |>
  left_join(ebms_clim_df, by = c("SITE_ID" = "transect_id"))

head(pheno_clim_df)

# Count observations per phenovar and climatic region
clim_counts <- pheno_clim_df %>%
  group_by(genzname) %>%
  summarise(
    N_obs = n(),
    N_species = n_distinct(SPECIES),
    N_sites = n_distinct(SITE_ID),
    .groups = "drop"
  ) %>%
  arrange(genzname)

clim_counts

# Model

models_clim <- pheno_clim_df %>%
  split(.$phenovar) %>%           # fit one model per phenological metric
  map(~ lmer(
    scale(estimate) ~             # response: standardized phenological trend (slope)
      genzname +                  # fixed effect: climatic region
      (1 | SPECIES) +             # random intercept: species
      (1 | SITE_ID),              # random intercept: site
    data = .x,
    weights = 1 / (std.error^2),  # Precision weights (inverse of trend variance)
    REML = FALSE
  ))

results_fixed <- map_df(
  models_clim,
  ~ tidy(.x, effects = "fixed"),
  .id = "phenovar"
)

print(results_fixed, n = Inf)

anova_results <- imap_dfr(
  models_clim,
  function(m1, name) {
    
    m0 <- update(m1, . ~ . - genzname)
    
    a <- anova(m0, m1)  # LRT
    
    data.frame(
      phenovar = name,
      Df = a$Df[2],
      Chisq = a$Chisq[2],
      p.value = a$`Pr(>Chisq)`[2]
    )
  }
)

anova_results

# Plot estimated marginal means for each phenovar

# Estimated marginal means
emm_all <- emm_all %>%
  mutate(
    phenovar = factor(phenovar, levels = vars_order),
    phenovar_label = factor(pretty_names[phenovar],
                            levels = pretty_names[vars_order])
  )


emm_all$genzname <- factor(
  emm_all$genzname,
  levels = c(
    "E. Cold and wet",
    "F. Extremely cold and mesic",
    "G. Cold and mesic",
    "H. Cool temperate and dry",
    "J. Cool temperate and moist",
    "K. Warm temperate and mesic",
    "L. Warm temperate and xeric"
  )
)

emm_all$genz_letter <- factor(
  sub("\\..*", "", emm_all$genzname),
  levels = c("E","F","G","H","J","K","L")
)

trends_by_rclim <- ggplot(emm_all, aes(x = genz_letter, y = emmean)) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.15) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ phenovar_label, scales = "free_y")+
  labs(
    x = "Climatic region",
    y = "Estimated temporal trend (days/year)"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold")
  )

trends_by_rclim

ggsave(
  filename = here::here("output", "figures", "trends_by_rclims.png"),
  plot = trends_by_rclim,
  width = 9,
  height = 5,
  dpi = 300
)


#### Pheno trends vs. latitude ####

ebms_coord_df  <- read.csv(here("data", "ebms_transect_coord.csv"), sep = ",", dec = ".")

pheno_coord_df <- pheno_trends_site |>
  left_join(ebms_coord_df, by = c("SITE_ID" = "transect_id"))

head(pheno_coord_df)
str(pheno_coord_df)

models_latlon <- pheno_coord_df %>%
  split(.$phenovar) %>%
  map(~ lmer(
    scale(estimate) ~ 
      scale(transect_lat) +
      scale(transect_lon) +
      scale(N_years) +
      (1 | SPECIES) +
      (1 | SITE_ID),
    data = .x,
    weights = 1 / (std.error^2),
    REML = FALSE
  ))

results_latlon <- map_df(
  models_latlon,
  broom::tidy,
  .id = "phenovar"
)

results_latlon

lrt_results_latlon <- imap_dfr(
  models_latlon,
  function(m_full, name) {
    
    # Remove latitude
    m_no_lat <- update(m_full, . ~ . - scale(transect_lat))
    a_lat <- anova(m_no_lat, m_full)
    
    # Remove longitude
    m_no_lon <- update(m_full, . ~ . - scale(transect_lon))
    a_lon <- anova(m_no_lon, m_full)
    
    # Remove N_years
    m_no_years <- update(m_full, . ~ . - scale(N_years))
    a_years <- anova(m_no_years, m_full)
    
    data.frame(
      phenovar = name,
      
      lat_Chisq  = a_lat$Chisq[2],
      lat_df     = a_lat$Df[2],
      lat_p      = a_lat$`Pr(>Chisq)`[2],
      
      lon_Chisq  = a_lon$Chisq[2],
      lon_df     = a_lon$Df[2],
      lon_p      = a_lon$`Pr(>Chisq)`[2],
      
      years_Chisq = a_years$Chisq[2],
      years_df    = a_years$Df[2],
      years_p     = a_years$`Pr(>Chisq)`[2]
    )
  }
)

lrt_results_latlon

# Plot latitudinal effects

pred_latitude <- map_df(
  models_latlon,
  ~ as.data.frame(ggpredict(.x, terms = "transect_lat [all]")),
  .id = "phenovar"
)

pred_lon <- map_df(
  models_latlon,
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


#### Gam spatial approach ####

# ---- Prepare data ----

pheno_coord_df <- pheno_coord_df %>%
  mutate(
    SPECIES = factor(SPECIES),
    SITE_ID = factor(SITE_ID)
  )

# Global scaling (store attributes)
lon_sc_obj <- scale(pheno_coord_df$transect_lon)
lat_sc_obj <- scale(pheno_coord_df$transect_lat)

pheno_coord_df$lon_sc <- as.numeric(lon_sc_obj)
pheno_coord_df$lat_sc <- as.numeric(lat_sc_obj)

lon_center <- attr(lon_sc_obj, "scaled:center")
lon_scale  <- attr(lon_sc_obj, "scaled:scale")
lat_center <- attr(lat_sc_obj, "scaled:center")
lat_scale  <- attr(lat_sc_obj, "scaled:scale")

# ---- Fit models ----

models_gam <- pheno_coord_df %>%
  split(.$phenovar) %>%             # fit one model per phenological metric
  map(~ bam(
    estimate ~                      # response: phenological trend (slope)
      s(lon_sc, lat_sc, k = 50) +   # 2D smooth for spatial trend (longitude x latitude)
      s(SPECIES, bs = "re") +       # random effect: species
      s(SITE_ID, bs = "re"),        # random effect: site
    data = .x,
    weights = 1 / (std.error^2),    # Precision weights (inverse of trend variance)
    method = "fREML",
    discrete = TRUE
  ))

# ---- Extract smooth table ----

gam_results <- map_df(
  models_gam,
  ~ as.data.frame(summary(.x)$s.table),
  .id = "phenovar"
)

# ---- Europe map ----

europe_3035 <- ne_countries(
  continent = "Europe",
  scale = "medium",
  returnclass = "sf"
) %>%
  st_transform(3035)

# Desired order
desired_order <- c(
  "ONSET_mean", "ONSET_var", "PEAKDAY",
  "OFFSET_mean", "OFFSET_var",
  "FLIGHT_LENGTH_mean", "FLIGHT_LENGTH_var"
)

models_gam <- models_gam[desired_order]

pretty_names <- c(
  ONSET_mean = "Onset (mean)",
  ONSET_var = "Onset (var)",
  PEAKDAY = "Peak day",
  OFFSET_mean = "Offset (mean)",
  OFFSET_var = "Offset (var)",
  FLIGHT_LENGTH_mean = "Flight length (mean)",
  FLIGHT_LENGTH_var = "Flight length (var)"
)

# ---- Plot function ----

plot_gam_surface <- function(pv) {
  
  mod <- models_gam[[pv]]
  dat <- pheno_coord_df %>% filter(phenovar == pv)
  
  # Prediction grid
  grid <- expand.grid(
    transect_lon = seq(min(dat$transect_lon),
                       max(dat$transect_lon),
                       length.out = 80),
    transect_lat = seq(min(dat$transect_lat),
                       max(dat$transect_lat),
                       length.out = 80)
  )
  
  # Use global scaling
  grid$lon_sc <- (grid$transect_lon - lon_center) / lon_scale
  grid$lat_sc <- (grid$transect_lat - lat_center) / lat_scale
  
  # Dummy RE levels
  grid$SPECIES <- dat$SPECIES[1]
  grid$SITE_ID <- dat$SITE_ID[1]
  
  # Extract spatial smooth explicitly
  terms_pred <- predict(mod, newdata = grid, type = "terms")
  grid$pred  <- terms_pred[, "s(lon_sc,lat_sc)"]
  
  # Convert sites to sf
  sites_sf <- st_as_sf(
    dat,
    coords = c("transect_lon", "transect_lat"),
    crs = 3035
  )
  
  # Bounding box
  bbox <- st_bbox(sites_sf)
  expand_factor <- 3e5
  
  bbox_exp <- bbox
  bbox_exp[c("xmin","xmax")] <- bbox[c("xmin","xmax")] + c(-expand_factor, expand_factor)
  bbox_exp[c("ymin","ymax")] <- bbox[c("ymin","ymax")] + c(-expand_factor, expand_factor)
  
  europe_crop <- st_crop(europe_3035, bbox_exp)
  
  # Plot
  ggplot() +
    geom_sf(data = europe_crop,
            fill = "grey95",
            color = "grey70",
            linewidth = 0.2) +
    geom_tile(data = grid,
              aes(transect_lon, transect_lat, fill = pred),
              alpha = 0.9) +
    geom_sf(data = sites_sf,
            size = 0.3,
            alpha = 0.5,
            color = "black") +
    scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      name = "Spatial\ntrend"
    ) +
    coord_sf(
      xlim = c(bbox_exp["xmin"], bbox_exp["xmax"]),
      ylim = c(bbox_exp["ymin"], bbox_exp["ymax"]),
      expand = FALSE
    ) +
    theme_void() +
    labs(title = pretty_names[pv])
}

# ---- Generate all plots ----

plot_maps <- map(names(models_gam), plot_gam_surface)

combined_map_fig <- wrap_plots(plot_maps, ncol = 3)


ggsave(
  filename = here::here("output", "figures", "pheno_trends_lat_lon.png"),
  plot = combined_map_fig,
  width = 7,
  height = 9,
  dpi = 300
)

#### Pheno trends vs. temperature ####
temp_df  <- read.csv(here("output", "climate", "mean_temperature_site.csv"), sep = ",", dec = ".")

pheno_temp_df <- pheno_trends_site |>
  left_join(temp_df, by = c("SITE_ID" = "transect_id"))

head(pheno_temp_df)
str(pheno_temp_df)

pheno_temp_df <- pheno_temp_df %>%
  mutate(mean_temp_sc = scale(mean_temp)[,1])

models_temp <- pheno_temp_df %>%
  split(.$phenovar) %>%            # fit one model per phenological metric
  map(~ lmer(
    estimate ~                     # response: phenological trend (slope)
      mean_temp_sc +               # fixed effect: standardized mean temperature      
      (1 | SPECIES) +              # random intercept: species
      (1 | SITE_ID),               # random intercept: site
    data = .x,
    weights = 1 / (std.error^2),   # Precision weights (inverse of trend variance)
    REML = FALSE
  ))


results_temp <- models_temp %>%
  purrr::map_df(
    ~ broom::tidy(.x, conf.int = TRUE),
    .id = "phenovar"
  ) %>%
  dplyr::filter(term == "mean_temp_sc")

results_temp

anova_temp <- pheno_temp_df %>%
  split(.$phenovar) %>%
  imap_dfr(function(df, name) {
    
    # Keep only complete cases for ALL variables used in full model
    df2 <- df %>%
      dplyr::filter(
        !is.na(estimate),
        !is.na(mean_temp_sc),
        !is.na(std.error),
        !is.na(SPECIES),
        !is.na(SITE_ID)
      )
    
    # Now both models use EXACT same df2
    m_full <- lmer(
      estimate ~ mean_temp_sc +
        (1 | SPECIES) +
        (1 | SITE_ID),
      data = df2,
      weights = 1 / (std.error^2),
      REML = FALSE
    )
    
    m_null <- lmer(
      estimate ~ 
        (1 | SPECIES) +
        (1 | SITE_ID),
      data = df2,
      weights = 1 / (std.error^2),
      REML = FALSE
    )
    
    a <- anova(m_null, m_full)
    
    data.frame(
      phenovar = name,
      Chisq = a$Chisq[2],
      df = a$Df[2],
      p.value = a$`Pr(>Chisq)`[2]
    )
  })

anova_temp

# Plot effects

desired_order <- c(
  "ONSET_mean",
  "ONSET_var",
  "PEAKDAY",
  "OFFSET_mean",
  "OFFSET_var",
  "FLIGHT_LENGTH_mean",
  "FLIGHT_LENGTH_var"
)

pretty_names <- c(
  ONSET_mean = "Onset (mean)",
  ONSET_var = "Onset (variance)",
  PEAKDAY = "Peak day",
  OFFSET_mean = "Offset (mean)",
  OFFSET_var = "Offset (variance)",
  FLIGHT_LENGTH_mean = "Flight length (mean)",
  FLIGHT_LENGTH_var = "Flight length (variance)"
)

pred_temp <- purrr::map_df(
  models_temp,
  ~ as.data.frame(ggeffects::ggpredict(.x, terms = "mean_temp_sc [all]")),
  .id = "phenovar"
)

pred_temp$phenovar <- factor(
  pred_temp$phenovar,
  levels = desired_order
)

pred_temp$phenovar <- factor(
  pred_temp$phenovar,
  labels = pretty_names[desired_order]
)

pheno_trends_temp <- ggplot(pred_temp, aes(x = x, y = predicted)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ phenovar, scales = "free_y") +
  theme_bw() +
  labs(
    x = "Standardized mean annual temperature",
    y = "Trend (days per year)",
    title = ""
  )

ggsave(
  filename = here::here("output", "figures", "pheno_trends_meantemp.png"),
  plot = pheno_trends_temp,
  width = 9,
  height = 5,
  dpi = 300
)


# Interpretation: FL decrease at warmer region because of stronger advance in the offset.
# Are we capturing the actual offset or there is a methodological issue here???

