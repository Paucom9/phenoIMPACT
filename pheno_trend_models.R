# ============================================================================================ #
# 03_pheno_trend_models.R
#
# Author: Pau Colom
# Date: 2026-02-20
#
# Description:
#
# ============================================================================================ #

#### Load required libraries ####
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

pheno_trends_site  <- read.csv(here("output", "pheno_temporal_trends_allspp.csv"), sep = ",", dec = ".")

str(pheno_trends_site)

#-----


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
  map(~ lmer(
    scale(estimate) ~ genzname +
      (1 | SPECIES) +
      (1 | SITE_ID),
    data = .x,
    weights = 1 / (std.error^2),
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
emm_all <- imap_dfr(
  models_clim,
  function(model, name) {
    
    emm <- emmeans(model, ~ genzname)
    df  <- as.data.frame(emm)
    
    df$phenovar <- name
    df
  }
)

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

ggplot(emm_all, aes(x = genz_letter, y = emmean)) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.15) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ phenovar, scales = "free_y") +
  labs(
    x = "Climatic region",
    y = "Estimated temporal trend (days/year)"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold")
  )








# --- Pheno trends vs. latitude --- #

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


# Gam spatial approach

pheno_coord_df$SPECIES <- factor(pheno_coord_df$SPECIES)
pheno_coord_df$SITE_ID <- factor(pheno_coord_df$SITE_ID)

pheno_coord_df <- pheno_coord_df %>%
  mutate(
    lon_sc = scale(transect_lon),
    lat_sc = scale(transect_lat)
  )


models_gam <- pheno_coord_df %>%
  split(.$phenovar) %>%
  map(~ bam(
    estimate ~ 
      s(lon_sc, lat_sc, k = 50) +
      s(SPECIES, bs = "re") +
      s(SITE_ID, bs = "re"),
    data = .x,
    weights = 1 / (std.error^2),
    method = "fREML",
    discrete = TRUE
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

# Order
desired_order <- c(
  "ONSET_mean",
  "ONSET_var",
  "PEAKDAY",
  "OFFSET_mean",
  "OFFSET_var",
  "FLIGHT_LENGTH_mean",
  "FLIGHT_LENGTH_var"
)

models_gam <- models_gam[desired_order]

names(models_gam)

# Cleaner short labels
pretty_names <- c(
  ONSET_mean = "Onset (mean)",
  ONSET_var = "Onset (var)",
  PEAKDAY = "Peak day",
  OFFSET_mean = "Offset (mean)",
  OFFSET_var = "Offset (var)",
  FLIGHT_LENGTH_mean = "Flight length (mean)",
  FLIGHT_LENGTH_var = "Flight length (var)"
)

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
  
  # Scale as in model
  grid$lon_sc <- (grid$transect_lon - mean(dat$transect_lon)) / sd(dat$transect_lon)
  grid$lat_sc <- (grid$transect_lat - mean(dat$transect_lat)) / sd(dat$transect_lat)
  
  # Dummy RE
  grid$SPECIES <- dat$SPECIES[1]
  grid$SITE_ID <- dat$SITE_ID[1]
  
  # Spatial smooth only
  grid$pred <- predict(mod, newdata = grid, type = "terms")[,1]
  
  # Convert to sf
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
  
  ggplot() +
    geom_sf(data = europe_crop,
            fill = "grey95",
            color = "grey70",
            linewidth = 0.2) +
    geom_tile(data = grid_crop,
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

plots <- lapply(names(models_gam), plot_gam_surface)

wrap_plots(plots, ncol = 3, guides = "collect") &
  theme(legend.position = "right")

# --- Pheno trends vs. temperature --- #

temp_df  <- read.csv(here("output", "climate", "mean_temperature_site.csv"), sep = ",", dec = ".")

pheno_temp_df <- pheno_trends_site |>
  left_join(temp_df, by = c("SITE_ID" = "transect_id"))

head(pheno_temp_df)
str(pheno_temp_df)

pheno_temp_df <- pheno_temp_df %>%
  mutate(mean_temp_sc = scale(mean_temp)[,1])

models_temp <- pheno_temp_df %>%
  split(.$phenovar) %>%
  map(~ lmer(
    estimate ~ mean_temp_sc +
      (1 | SPECIES) +
      (1 | SITE_ID),
    data = .x,
    weights = 1 / (std.error^2),
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

ggplot(pred_temp, aes(x = x, y = predicted)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ phenovar, scales = "free_y") +
  theme_bw() +
  labs(
    x = "Mean annual temperature (°C)",
    y = "Trend (days per year)",
    title = "Temperature variation in phenological trends"
  )

# Interpretation: FL decrease at warmer region because of stronger advance in the offset.
# Are we capturing the actual offset or there is a methodological issue here???

