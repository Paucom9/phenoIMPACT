# ============================================================================================ #
# map_study_sites.R
#
# Author: Pau Colom
# Date: 2026-06-01
#
# Description: This script maps the eBMS transects/sites included in the phenological analyses 
# by combining the final phenology estimates with the transect coordinate file. 
# Transect coordinates are kept in EPSG:3035, the original projection of the eBMS coordinate file.
#
# ============================================================================================ #

# Clean session
rm(list = ls())

#### Load required libraries ####
# ---
library(dplyr)         # For data manipulation
library(stringr)       # For string manipulation
library(sf)            # For handling spatial data
library(ggplot2)       # For data visualization
library(rnaturalearth) # For country boundaries
library(here)          # For constructing project-independent file paths
library(extrafont)     # For publication fonts
# ---

#### Data Import and Preparation ####
# ---

here::here() # Check the current working directory

# --- Phenology estimates used in the analyses
phenology_estimates <- read.csv(
  here::here("output", "pheno_estimates_allspp.csv"),
  sep = ",",
  dec = "."
)

# --- Transect coordinates
# Note: despite the column names, transect_lon and transect_lat are projected
# coordinates in EPSG:3035, not longitude/latitude in WGS84.
ebms_transect_coord <- read.csv(
  here::here("data", "ebms_transect_coord.csv"),
  sep = ",",
  dec = "."
)
# ---

# ---------------------------------------------------------------------------- #
#### Prepare sites included in the phenological analyses ####
# ---------------------------------------------------------------------------- #

analysis_sites <- phenology_estimates |>
  dplyr::mutate(
    SITE_ID = as.character(SITE_ID),
    SPECIES = as.character(SPECIES),
    YEAR = as.integer(YEAR)
  ) |>
  dplyr::filter(!is.na(SITE_ID)) |>
  dplyr::distinct(SITE_ID)

analysis_sites |>
  dplyr::summarise(n_sites = dplyr::n())

# ---------------------------------------------------------------------------- #
#### Prepare coordinates in EPSG:3035 ####
# ---------------------------------------------------------------------------- #

coord_site <- ebms_transect_coord |>
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
  dplyr::rename(
    SITE_ID = transect_id,
    x_3035 = transect_lon,
    y_3035 = transect_lat
  )

# ---------------------------------------------------------------------------- #
#### Number of species retained per site ####
# ---------------------------------------------------------------------------- #

site_species <- phenology_estimates |>
  dplyr::mutate(
    SITE_ID = as.character(SITE_ID),
    SPECIES = as.character(SPECIES)
  ) |>
  dplyr::filter(
    !is.na(SITE_ID),
    !is.na(SPECIES)
  ) |>
  dplyr::group_by(SITE_ID) |>
  dplyr::summarise(
    n_species = dplyr::n_distinct(SPECIES),
    .groups = "drop"
  )

# ---------------------------------------------------------------------------- #
#### Prepare map data ####
# ---------------------------------------------------------------------------- #

map_sites <- analysis_sites |>
  dplyr::left_join(coord_site, by = "SITE_ID") |>
  dplyr::left_join(site_species, by = "SITE_ID") |>
  dplyr::mutate(
    bms_id = dplyr::case_when(
      !is.na(bms_id) ~ as.character(bms_id),
      TRUE ~ stringr::str_extract(SITE_ID, "^[^.]*")
    )
  )

# Check missing coordinates
missing_coords <- map_sites |>
  dplyr::filter(
    is.na(x_3035) |
      is.na(y_3035)
  )

missing_coords

# Keep only sites with valid coordinates
map_sites <- map_sites |>
  dplyr::filter(
    !is.na(x_3035),
    !is.na(y_3035)
  )

map_sites_sf <- map_sites |>
  sf::st_as_sf(
    coords = c("x_3035", "y_3035"),
    crs = 3035,
    remove = FALSE
  )

# ---------------------------------------------------------------------------- #
#### Summary tables ####
# ---------------------------------------------------------------------------- #

map_summary_by_bms <- map_sites |>
  dplyr::group_by(bms_id) |>
  dplyr::summarise(
    n_sites = dplyr::n_distinct(SITE_ID),
    mean_n_species = mean(n_species, na.rm = TRUE),
    median_n_species = median(n_species, na.rm = TRUE),
    min_n_species = min(n_species, na.rm = TRUE),
    max_n_species = max(n_species, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::arrange(dplyr::desc(n_sites))

map_summary_by_bms

phenology_summary_by_bms <- phenology_estimates |>
  dplyr::mutate(
    SITE_ID = as.character(SITE_ID),
    SPECIES = as.character(SPECIES),
    YEAR = as.integer(YEAR)
  ) |>
  dplyr::left_join(
    coord_site |>
      dplyr::select(SITE_ID, bms_id),
    by = "SITE_ID"
  ) |>
  dplyr::mutate(
    bms_id = dplyr::case_when(
      !is.na(bms_id) ~ as.character(bms_id),
      TRUE ~ stringr::str_extract(SITE_ID, "^[^.]*")
    )
  ) |>
  dplyr::group_by(bms_id) |>
  dplyr::summarise(
    n_sites = dplyr::n_distinct(SITE_ID),
    n_species = dplyr::n_distinct(SPECIES),
    first_year = min(YEAR, na.rm = TRUE),
    last_year = max(YEAR, na.rm = TRUE),
    n_years = dplyr::n_distinct(YEAR),
    n_pheno_estimates = dplyr::n(),
    .groups = "drop"
  ) |>
  dplyr::arrange(dplyr::desc(n_sites))

phenology_summary_by_bms

# Save summary tables
out_table_dir <- here::here("output", "figures", "fig1_map_tables")
dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(
  map_summary_by_bms,
  file = file.path(out_table_dir, "fig1_map_sites_by_bms.csv"),
  row.names = FALSE
)

write.csv(
  phenology_summary_by_bms,
  file = file.path(out_table_dir, "fig1_map_phenology_summary_by_bms.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------------------------- #
#### Prepare map background in EPSG:3035 ####
# ---------------------------------------------------------------------------- #

map_background <- rnaturalearth::ne_countries(
  scale = "medium",
  returnclass = "sf"
) |>
  dplyr::filter(
    continent == "Europe" |
      admin %in% c("Turkey", "Morocco", "Algeria", "Tunisia")
  ) |>
  sf::st_transform(3035)

# Dynamic map limits based on site coordinates, in metres
bbox_sites <- sf::st_bbox(map_sites_sf)

map_padding <- 250000 # 250 km

xlim_map <- c(
  as.numeric(bbox_sites["xmin"]) - map_padding,
  as.numeric(bbox_sites["xmax"]) + map_padding
)

ylim_map <- c(
  as.numeric(bbox_sites["ymin"]) - map_padding,
  as.numeric(bbox_sites["ymax"]) + map_padding
)

# ---------------------------------------------------------------------------- #
#### Plot map: clean manuscript version ####
# ---------------------------------------------------------------------------- #

fig1a_map_clean <- ggplot() +
  geom_sf(
    data = map_background,
    fill = "grey92",
    colour = "grey70",
    linewidth = 0.25
  ) +
  geom_sf(
    data = map_sites_sf,
    aes(fill = n_species),
    shape = 21,
    colour = "black",
    stroke = 0.2,
    size = 1.7,
    alpha = 0.9
  ) +
  scale_fill_viridis_c(
    option = "H",
    name = "No. species",
    guide = guide_colourbar(
      barheight = grid::unit(2.2, "cm"),
      barwidth  = grid::unit(0.5, "cm"),
      title.position = "top",
      title.hjust = 0.5
    )
  ) +
  coord_sf(
    crs = sf::st_crs(3035),
    datum = sf::st_crs(4326),
    xlim = xlim_map,
    ylim = ylim_map,
    expand = FALSE
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_classic(
    base_family = "Garamond",
    base_size = 14
  ) +
  theme(
    axis.text = element_text(
      colour = "black",
      size = 11
    ),
    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.3
    ),
    axis.line = element_blank(),
    
    panel.grid.major = element_line(
      colour = "grey85",
      linewidth = 0.35
    ),
    panel.background = element_rect(
      fill = "white",
      colour = NA
    ),
    plot.background = element_rect(
      fill = "white",
      colour = NA
    ),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.5
    ),
    
    legend.position = c(0.93, 0.05),
    legend.justification = c(1, 0),
    legend.background = element_rect(
      fill = "white",
      colour = "black",
      linewidth = 0.35
    ),
    legend.title = element_text(
      face = "bold",
      size = 10
    ),
    legend.text = element_text(
      size = 10
    ),
    legend.margin = margin(4, 5, 4, 5),
    legend.spacing.y = grid::unit(0.1, "cm"),
    
    plot.margin = margin(5, 5, 5, 5)
  )

fig1a_map_clean

# ---------------------------------------------------------------------------- #
#### Save figures ####
# ---------------------------------------------------------------------------- #

out_fig_dir <- here::here("output", "figures")
dir.create(out_fig_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(
  filename = file.path(out_fig_dir, "fig1a_map_sites_europe_clean_3035.png"),
  plot = fig1a_map_clean,
  width = 5.8,
  height = 6.2,
  dpi = 600
)


# ============================================================================================ #
# End of script
# ============================================================================================ #