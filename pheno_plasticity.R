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
# ---

#### Data Import and Preparation ####
# ---
here::here() # Check the current working directory
# ---
phenology_estimates  <- read.csv(here::here("output", "pheno_estimates_allspp.csv"), sep = ",", dec = ".")
clim_vars <- read.csv(here::here("output", "climate", "climate_variables.csv"), sep = ",", dec = ".")
ebms_transect_coord <- read.csv(here::here("data", "ebms_transect_coord.csv"), sep = ",", dec = ".")
species_traits <- read.csv(here::here("data", "species_trait_table.csv"), sep = ";", dec = ".")

# ---

#### Merge datasets ####
# 1. Filter and convert to sf
coords_sf <- ebms_coord_df |>
  filter(!is.na(transect_lon),
         !is.na(transect_lat)) |>
  distinct(transect_id, transect_lon, transect_lat) |>
  st_as_sf(
    coords = c("transect_lon", "transect_lat"),
    crs = 3035
  )

# 2. Transform to WGS84
coords_wgs84 <- st_transform(coords_sf, 4326)

# 3. Extract latitude
coord_site <- coords_wgs84 |>
  mutate(latitude = st_coordinates(coords_wgs84)[,2]) |>
  st_drop_geometry() |>
  select(transect_id, latitude)

clim_vars <- clim_vars |>
  left_join(coord_site, by = c("SITE_ID" = "transect_id")) |>
  mutate(latitude_sc = scale(latitude)[,1])

df <- phenology_estimates |>
  left_join(
    clim_vars,
    by = c("SPECIES","SITE_ID", "YEAR")
  )

str(df)


df <- df |>
  mutate(
    SPECIES = factor(SPECIES),
    SITE_ID = factor(SITE_ID)
  )

df$sp_site <- interaction(df$SPECIES, df$SITE_ID, drop = TRUE)

# Inspect climatic anomalies

med <- median(df$clim_anomaly_temp_30, na.rm = TRUE)

p <- ggplot(df, aes(x = clim_anomaly_temp_30)) +
  geom_histogram(fill = "grey", color = "white", bins = 30) +
  
  # line at zero
  geom_vline(xintercept = 0, linetype = "dashed", size = 1) +
  
  # line at median
  geom_vline(xintercept = med, color = "red", linetype = "dashed", size = 1) +
  
  # annotation
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
  filename = here::here("output", "figures", "distribution_clim_anomalies_26_3_26.png"),
  plot = p,
  width = 8,
  height = 6
)

#### Inspect predictor variables ####

df_cor <- df |>
  distinct(SITE_ID,
           clim_background_temp_90,
           clim_trend_temp_90,
           clim_autocorr_temp_90,
           clim_stability_temp_90,
           clim_predictability_temp_90,
           photo_90,
           latitude)

# --- long format ---
df_long <- df_cor |>
  pivot_longer(
    cols = -SITE_ID,
    names_to = "variable",
    values_to = "value"
  )

# --- labels bonics ---
labels <- c(
  clim_background_temp_90   = "Background temperature\n(°C)",
  clim_trend_temp_90        = "Temperature trend\n(°C per decade)",
  clim_autocorr_temp_90     = "Autocorrelation\n(lag-1)",
  clim_stability_temp_90    = "Stability\n(-SD)",
  clim_predictability_temp_90 = "Predictability\n(-Var residuals)",
  photo_90 = "Photoperiod\n(daily hours)",
  latitude                  = "Latitude\n(degrees)"
)

df_long$variable <- factor(df_long$variable,
                           levels = names(labels),
                           labels = labels)

# --- plot ---
p <- ggplot(df_long, aes(x = value)) +
  geom_histogram(bins = 30, fill = "grey70", color = "black") +
  facet_wrap(~ variable, scales = "free", ncol = 3) +
  theme_classic(base_size = 12) +
  labs(x = NULL, y = "Frequency") +
  theme(
    strip.text = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
  )

# --- save ---
ggsave(
  filename = here::here("output", "figures", "climate_histograms_26_3_26B.png"),
  plot = p,
  width = 9,
  height = 7,
  dpi = 300
)

# ---

cor_mat <- df_cor |>
  select(clim_background_temp_90,
         clim_trend_temp_90,
         clim_autocorr_temp_90,
         clim_stability_temp_90,
         clim_predictability_temp_90,
         photo_90,
         latitude) |>
  cor(use = "complete.obs")


new_names <- c(
  "Background",
  "Trend",
  "Autocorr",
  "Stability",
  "Predictability",
  "Photoperiod",
  "Latitude"
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

ggsave(
  filename = here::here("output", "figures", "corr_site_vars_26_3_26B.png"),
  plot = corr_clim_pred_sds,
  width = 6,
  height = 5,
  dpi = 300
)


#### Models phenotypic plasticity ####

#### Model selection ####

model_set_plasticity <- function(d) {
  
  d <- d %>%
    dplyr::filter(
      !is.na(ONSET_mean),
      !is.na(clim_anomaly_temp_90_sc),
      !is.na(clim_background_temp_90_sc),
      !is.na(clim_predictability_temp_90_sc),
      !is.na(photo_90_sc)
    )
  
  ctrl <- lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e6)
  )
  
  f_base <- ONSET_mean ~ clim_anomaly_temp_90_sc +
    (1 + clim_anomaly_temp_90_sc || SPECIES) +
    (1 | SITE_ID)
  
  list(
    
    # H1
    m1 = lmer(f_base, data = d, REML = FALSE, control = ctrl),
    
    # H2
    m2 = lmer(update(f_base,
                     . ~ . + clim_background_temp_90_sc +
                       clim_anomaly_temp_90_sc:clim_background_temp_90_sc),
              data = d, REML = FALSE, control = ctrl),
    
    # H3
    m3 = lmer(update(f_base,
                     . ~ . + clim_predictability_temp_90_sc +
                       clim_anomaly_temp_90_sc:clim_predictability_temp_90_sc),
              data = d, REML = FALSE, control = ctrl),
    
    # H4
    m4 = lmer(update(f_base,
                     . ~ . + photo_90_sc +
                       clim_anomaly_temp_90_sc:photo_90_sc),
              data = d, REML = FALSE, control = ctrl),
    
    # H5
    m5 = lmer(update(f_base,
                     . ~ . +
                       clim_background_temp_90_sc +
                       clim_predictability_temp_90_sc +
                       clim_anomaly_temp_90_sc:clim_background_temp_90_sc +
                       clim_anomaly_temp_90_sc:clim_predictability_temp_90_sc),
              data = d, REML = FALSE, control = ctrl),
    
    # H6
    m6 = lmer(update(f_base,
                     . ~ . +
                       clim_background_temp_90_sc +
                       photo_90_sc +
                       clim_anomaly_temp_90_sc:clim_background_temp_90_sc +
                       clim_anomaly_temp_90_sc:photo_90_sc),
              data = d, REML = FALSE, control = ctrl),
    
    # H7
    m7 = lmer(update(f_base,
                     . ~ . +
                       clim_predictability_temp_90_sc +
                       photo_90_sc +
                       clim_anomaly_temp_90_sc:clim_predictability_temp_90_sc +
                       clim_anomaly_temp_90_sc:photo_90_sc),
              data = d, REML = FALSE, control = ctrl),
    
    # H8
    m8 = lmer(update(f_base,
                     . ~ . +
                       clim_background_temp_90_sc +
                       clim_predictability_temp_90_sc +
                       photo_90_sc +
                       clim_anomaly_temp_90_sc:clim_background_temp_90_sc +
                       clim_anomaly_temp_90_sc:clim_predictability_temp_90_sc +
                       clim_anomaly_temp_90_sc:photo_90_sc),
              data = d, REML = FALSE, control = ctrl)
  )
}

mods <- model_set_plasticity(df)

sel <- MuMIn::model.sel(mods)

sel_df <- as.data.frame(sel) %>%
  tibble::rownames_to_column("model") %>%
  dplyr::filter(delta <= 2)

sel_df

sel_full <- as.data.frame(sel) %>%
  tibble::rownames_to_column("model") %>%
  dplyr::arrange(delta)

# Check collinearity
performance::check_collinearity(mods$m8)


table_final <- purrr::map_dfr(seq_len(nrow(sel_full)), function(i) {
  
  mod_name <- sel_full$model[i]
  mod      <- mods[[mod_name]]
  
  coefs <- broom.mixed::tidy(mod, effects = "fixed")
  
  if (!"p.value" %in% names(coefs)) {
    coefs$p.value <- NA
  }
  
  coefs %>%
    dplyr::select(term, estimate, std.error, p.value) %>%
    dplyr::mutate(
      model = mod_name,
      delta = sel_full$delta[i],
      weight = sel_full$weight[i]
    ) %>%
    tidyr::pivot_wider(
      names_from = term,
      values_from = c(estimate, std.error, p.value)
    )
})

table_final <- table_final %>%
  dplyr::rename(
    Plasticity = estimate_clim_anomaly_temp_90_sc,
    Background = estimate_clim_background_temp_90_sc,
    Predictability = estimate_clim_predictability_temp_90_sc,
    Photoperiod = estimate_photo_90_sc,
    
    `Plasticity×Background` =
      `estimate_clim_anomaly_temp_90_sc:clim_background_temp_90_sc`,
    
    `Plasticity×Predictability` =
      `estimate_clim_anomaly_temp_90_sc:clim_predictability_temp_90_sc`,
    
    `Plasticity×Photoperiod` =
      `estimate_clim_anomaly_temp_90_sc:photo_90_sc`
  )


table_final %>%
  dplyr::select(model, delta, weight,
                Plasticity, Background, Predictability, Photoperiod,
                `Plasticity×Background`,
                `Plasticity×Predictability`,
                `Plasticity×Photoperiod`)



stars <- function(p) {
  dplyr::case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE ~ ""
  )
}

table_clean <- table_final %>%
  dplyr::transmute(
    model,
    delta = round(delta, 1),
    
    weight = dplyr::case_when(
      weight == 0 ~ "<1e-50",
      weight < 1e-3 ~ formatC(weight, format = "e", digits = 2),
      TRUE ~ formatC(weight, format = "f", digits = 3)
    ),
    
    Plasticity = paste0(
      round(Plasticity, 2),
      stars(p.value_clim_anomaly_temp_90_sc)
    ),
    
    Background = paste0(
      round(Background, 2),
      stars(p.value_clim_background_temp_90_sc)
    ),
    
    Predictability = paste0(
      round(Predictability, 2),
      stars(p.value_clim_predictability_temp_90_sc)
    ),
    
    Photoperiod = paste0(
      round(Photoperiod, 2),
      stars(p.value_photo_90_sc)
    ),
    
    `Plasticity×Background` = paste0(
      round(`Plasticity×Background`, 2),
      stars(`p.value_clim_anomaly_temp_90_sc:clim_background_temp_90_sc`)
    ),
    
    `Plasticity×Predictability` = paste0(
      round(`Plasticity×Predictability`, 2),
      stars(`p.value_clim_anomaly_temp_90_sc:clim_predictability_temp_90_sc`)
    ),
    
    `Plasticity×Photoperiod` = paste0(
      round(`Plasticity×Photoperiod`, 2),
      stars(`p.value_clim_anomaly_temp_90_sc:photo_90_sc`)
    )
  ) %>%
  dplyr::arrange(delta) %>%
  dplyr::mutate(
    across(-c(model, delta, weight),
           ~ ifelse(. %in% c("NA", NA), "—", .))
  )

table_clean



#### Plot effects ####

#1. Forest plot function

make_forest <- function(model){
  
  plot_df <- broom.mixed::tidy(model, effects = "fixed") %>%
    dplyr::filter(term != "(Intercept)") %>%
    dplyr::mutate(
      lower = estimate - 1.96 * std.error,
      upper = estimate + 1.96 * std.error,
      
      term = dplyr::recode(
        term,
        clim_anomaly_temp_90_sc = "Plasticity",
        clim_background_temp_90_sc = "Background",
        clim_predictability_temp_90_sc = "Predictability",
        photo_90_sc = "Photoperiod",
        `clim_anomaly_temp_90_sc:clim_background_temp_90_sc` = "Plasticity × Background",
        `clim_anomaly_temp_90_sc:clim_predictability_temp_90_sc` = "Plasticity × Predictability",
        `clim_anomaly_temp_90_sc:photo_90_sc` = "Plasticity × Photoperiod"
      ),
      
      term = factor(term, levels = rev(c(
        "Plasticity",
        "Background",
        "Predictability",
        "Photoperiod",
        "Plasticity × Background",
        "Plasticity × Predictability",
        "Plasticity × Photoperiod"
      )))
    )
  
  p <- ggplot(plot_df,
              aes(x = term,
                  y = estimate,
                  ymin = lower,
                  ymax = upper)) +
    geom_pointrange(size = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    coord_flip() +
    theme_classic() +
    labs(
      x = "",
      y = "Effect size (estimate ± 95% CI)"
    )
  
  return(p)
}

p <- make_forest(mods$m8)
p

ggsave(
  filename = here::here("output", "figures", "forest_plasticity_m8B.png"),
  plot = p,
  width = 6,
  height = 4,
  dpi = 300
)


#2. Plot prediction effects.

final_mod <- mods$m8


# Plasticity x background

pred <- ggpredict(
  final_mod,
  terms = c("clim_anomaly_temp_90_sc", "clim_background_temp_90_sc[-1:1]")
)

pred$group <- factor(pred$group, levels = c("1", "0", "-1"))

pred$group <- factor(
  pred$group,
  levels = c("1", "0", "-1"),
  labels = c("High", "Mean", "Low")
)

ggplot(pred, aes(x, predicted, color = group)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group),
              alpha = 0.2, color = NA) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  
  scale_color_manual(
    values = c(
      "High" = "#6C6BD1",  # blau-lila
      "Mean" = "#F28E2B",  # taronja
      "Low"  = "#2CA58D"   # verd-turquesa
    )
  ) +
  scale_fill_manual(
    values = c(
      "High" = "#6C6BD1",
      "Mean" = "#F28E2B",
      "Low"  = "#2CA58D"
    )
  ) +
  
  labs(
    x = "Temperature anomaly",
    y = "Onset (day of the year)",
    color = "Mean temperature",
    fill = "Mean temperature"
  ) +
  theme_classic(base_family = "Garamond") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

# Plasticity x photoperiod

pred <- ggpredict(
  final_mod,
  terms = c("clim_anomaly_temp_90_sc", "photo_90_sc[-1:1]")
)

pred$group <- factor(pred$group, levels = c("1", "0", "-1"))

pred$group <- factor(
  pred$group,
  levels = c("1", "0", "-1"),
  labels = c("High", "Mean", "Low")
)

ggplot(pred, aes(x, predicted, color = group)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group),
              alpha = 0.2, color = NA) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  
  scale_color_manual(
    values = c(
      "High" = "#6C6BD1",  # blau-lila
      "Mean" = "#F28E2B",  # taronja
      "Low"  = "#2CA58D"   # verd-turquesa
    )
  ) +
  scale_fill_manual(
    values = c(
      "High" = "#6C6BD1",
      "Mean" = "#F28E2B",
      "Low"  = "#2CA58D"
    )
  ) +
  
  labs(
    x = "Temperature anomaly",
    y = "Onset (day of the year)",
    color = "Photoperiod",
    fill = "Photoperiod"
  ) +
  theme_classic(base_family = "Garamond") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )



# Plasticity x predictability

pred <- ggpredict(
  final_mod,
  terms = c("clim_anomaly_temp_90_sc", "clim_predictability_temp_90_sc[-1:1]")
)

pred$group <- factor(pred$group, levels = c("1", "0", "-1"))

pred$group <- factor(
  pred$group,
  levels = c("1", "0", "-1"),
  labels = c("High", "Mean", "Low")
)


ggplot(pred, aes(x, predicted, color = group)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group),
              alpha = 0.2, color = NA) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  
  scale_color_manual(
    values = c(
      "High" = "#6C6BD1",  # blau-lila
      "Mean" = "#F28E2B",  # taronja
      "Low"  = "#2CA58D"   # verd-turquesa
    )
  ) +
  scale_fill_manual(
    values = c(
      "High" = "#6C6BD1",
      "Mean" = "#F28E2B",
      "Low"  = "#2CA58D"
    )
  ) +
  
  labs(
    x = "Temperature anomaly",
    y = "Onset (day of the year)",
    color = "Predictability",
    fill = "Predictability"
  ) +
  theme_classic(base_family = "Garamond") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )
