# ============================================================================================ #
# plot_conceptual_plasticity.R
#
# Author: Pau Colom
# Date: 2026-06-01
#
# Description: Conceptual panel for Figure 1 showing how environmental context
# modulates phenological plasticity to temperature anomalies.
#
# Onset is shown in the first row:
#   - warmer years -> earlier onset
#
# Offset is shown in the second row:
#   - univoltine populations: warmer years -> earlier offset
#   - multivoltine populations: warmer years -> later offset
#
# Environmental contexts are shown as Low / Average / High.
#
# ============================================================================================ #

# Clean session
rm(list = ls())

#### Load required libraries ####
library(dplyr)
library(tidyr)
library(ggplot2)
library(here)
library(extrafont)

# ---------------------------------------------------------------------------- #
#### Create conceptual data ####
# ---------------------------------------------------------------------------- #

x_seq <- seq(-2, 2, length.out = 200)

# IMPORTANT:
# level_class = Low / Average / High value of the environmental variable itself
#
# Hypotheses represented:
# - Daylength: higher daylength -> stronger plasticity
# - Mean temperature: lower mean temperature -> stronger plasticity
# - Temperature predictability: higher predictability -> stronger plasticity
# - Temperature trend: lower warming trend -> stronger plasticity

moderator_info <- tibble::tribble(
  ~moderator,                   ~level_class, ~level_rank, ~slope_mag,
  
  "Daylength",                  "Low",        1,           0.7,
  "Daylength",                  "Average",    2,           1.7,
  "Daylength",                  "High",       3,           2.8,
  
  "Mean temperature",           "Low",        1,           2.8,
  "Mean temperature",           "Average",    2,           1.7,
  "Mean temperature",           "High",       3,           0.7,
  
  "Temperature predictability", "Low",        1,           0.7,
  "Temperature predictability", "Average",    2,           1.6,
  "Temperature predictability", "High",       3,           2.6,
  
  "Temperature trend",          "Low",        1,           2.6,
  "Temperature trend",          "Average",    2,           1.5,
  "Temperature trend",          "High",       3,           0.6
)

moderator_info <- moderator_info |>
  mutate(
    moderator = factor(
      moderator,
      levels = c(
        "Daylength",
        "Mean temperature",
        "Temperature predictability",
        "Temperature trend"
      )
    ),
    level_class = factor(
      level_class,
      levels = c("Low", "Average", "High")
    )
  )

# ---------------------------------------------------------------------------- #
#### Build conceptual line data ####
# ---------------------------------------------------------------------------- #

# Onset:
# warmer years -> earlier onset
onset_df <- moderator_info |>
  crossing(anomaly = x_seq) |>
  mutate(
    panel_row = "Onset",
    group = "Onset",
    timing = -slope_mag * anomaly
  )

# Offset in univoltines:
# warmer years -> earlier offset
offset_uni_df <- moderator_info |>
  crossing(anomaly = x_seq) |>
  mutate(
    panel_row = "Offset",
    group = "Univoltine",
    timing = -slope_mag * anomaly
  )

# Offset in multivoltines:
# warmer years -> later offset
offset_multi_df <- moderator_info |>
  crossing(anomaly = x_seq) |>
  mutate(
    panel_row = "Offset",
    group = "Multivoltine",
    timing = slope_mag * anomaly
  )

concept_df <- bind_rows(
  onset_df,
  offset_uni_df,
  offset_multi_df
) |>
  mutate(
    panel_row = factor(
      panel_row,
      levels = c("Onset", "Offset")
    ),
    group = factor(
      group,
      levels = c("Onset", "Univoltine", "Multivoltine")
    )
  )

# ---------------------------------------------------------------------------- #
#### Plot conceptual panel C ####
# ---------------------------------------------------------------------------- #

fig1c_conceptual <- ggplot() +
  
  # Onset lines: solid, coloured by environmental level
  geom_line(
    data = concept_df |>
      filter(panel_row == "Onset"),
    aes(
      x = anomaly,
      y = timing,
      group = interaction(moderator, level_class),
      colour = level_class
    ),
    linewidth = 1.15,
    lineend = "round"
  ) +
  
  # Offset lines: univoltine vs multivoltine distinguished by linetype
  geom_line(
    data = concept_df |>
      filter(panel_row == "Offset"),
    aes(
      x = anomaly,
      y = timing,
      group = interaction(moderator, group, level_class),
      colour = level_class,
      linetype = group
    ),
    linewidth = 1.15,
    lineend = "round"
  ) +
  
  facet_grid(
    panel_row ~ moderator,
    switch = "y"
  ) +
  
  # Blue - green - orange palette, similar to map
  scale_colour_manual(
    values = c(
      "Low" = "#3B528B",
      "Average" = "#22A884",
      "High" = "#F28E2B"
    ),
    breaks = c("High", "Average", "Low"),
    name = NULL
  ) +
  
  scale_linetype_manual(
    values = c(
      "Univoltine" = "solid",
      "Multivoltine" = "22"
    ),
    breaks = c("Univoltine", "Multivoltine"),
    name = NULL
  ) +
  
  scale_x_continuous(
    limits = c(-2.2, 2.2),
    breaks = c(-1.75, 0, 1.75),
    labels = c("Cold", "Average", "Warm"),
    expand = c(0.01, 0.01)
  ) +
  
  scale_y_continuous(
    limits = c(-6, 6),
    breaks = NULL,
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  
  labs(
    x = "Temperature anomaly",
    y = "Phenological timing"
  ) +
  
  theme_classic(
    base_family = "Garamond",
    base_size = 13
  ) +
  
  theme(
    strip.background = element_rect(
      fill = "grey95",
      colour = "black",
      linewidth = 0.4
    ),
    strip.placement = "outside",
    strip.text.x = element_text(
      face = "bold",
      size = 11
    ),
    # rotate row labels the other way for easier reading
    strip.text.y.left = element_text(
      angle = 90,
      face = "bold",
      size = 11
    ),
    
    axis.text = element_text(
      colour = "black",
      size = 10
    ),
    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.3
    ),
    axis.line = element_line(
      colour = "black",
      linewidth = 0.4
    ),
    
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.4
    ),
    
    panel.spacing.x = grid::unit(0.45, "cm"),
    panel.spacing.y = grid::unit(0.35, "cm"),
    
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.text = element_text(size = 10),
    legend.key.width = grid::unit(1.6, "cm"),
    
    plot.margin = margin(5, 5, 5, 5)
  ) +
  
  guides(
    colour = guide_legend(
      order = 1,
      override.aes = list(
        linetype = "solid",
        linewidth = 1.3
      )
    ),
    linetype = guide_legend(
      order = 2,
      override.aes = list(
        colour = "black",
        linewidth = 1.3
      )
    )
  )

fig1c_conceptual

# ---------------------------------------------------------------------------- #
#### Save figure ####
# ---------------------------------------------------------------------------- #

out_fig_dir <- here::here("output", "figures")
dir.create(out_fig_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(
  filename = file.path(out_fig_dir, "fig1c_conceptual_plasticity.png"),
  plot = fig1c_conceptual,
  width = 9.5,
  height = 5.2,
  dpi = 600
)

ggsave(
  filename = file.path(out_fig_dir, "fig1c_conceptual_plasticity.pdf"),
  plot = fig1c_conceptual,
  width = 9.5,
  height = 5.2
)

# ============================================================================================ #
# End of script
# ============================================================================================ #
