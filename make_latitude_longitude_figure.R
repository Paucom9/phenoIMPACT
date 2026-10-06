
# =============================================================================
# phenoIMPACT | Figure: geographic patterns in predicted phenological plasticity
# Reproducible plotting script for the latitude + longitude analysis
#
# PURPOSE
#   Produce the 3 x 2 manuscript figure:
#     top:    maps of within-species relative responses
#     bottom: partial within-species latitude relationships
#
# INTERPRETATION OF THE MAP COLOURS
#   Values are oriented biologically and centred within species.
#   Positive values mean:
#     onset                -> MORE onset advancement than the species mean
#     offset_univoltine    -> MORE offset advancement than the species mean
#     offset_multivoltine  -> MORE offset delay than the species mean
#
#   The map response also removes the fitted within-species longitude component.
#   It is therefore NOT absolute plasticity.
#
# INPUTS
#   1) The reviewed latitude-from-extracted-plasticity run:
#      output/phenology_plasticity/latitude_from_extracted_plasticity_15km/
#      run_20261005_124834_165c6dd85b83
#
#   2) longitude_sensitivity_results.csv supplied with this script.
#
# OUTPUTS
#   manuscript_latitude_longitude_patterns.png
#   manuscript_latitude_longitude_patterns.pdf
#   figure_site_values.csv
#   figure_partial_values.csv
#
# NOTES
#   - No phenology model is loaded or refitted.
#   - No spatial field is refitted.
#   - The latitude slope/CI labels are the reviewed source-model covariance
#     projections from the latitude + longitude sensitivity analysis.
#   - Font is explicitly requested as Garamond. If Garamond is unavailable,
#     R/graphics may substitute another serif font; check the exported PDF.
# =============================================================================


make_latitude_longitude_figure <- function(
  project_root = "E:/phenoIMPACT project/code/phenoIMPACT",
  latitude_run = file.path(
    project_root,
    "output", "phenology_plasticity",
    "latitude_from_extracted_plasticity_15km",
    "run_20261005_124834_165c6dd85b83"
  ),
  longitude_results = NULL,
  output_dir = file.path(latitude_run, "manuscript_figure_latitude_longitude"),
  font_family = "Garamond",
  width_in = 15,
  height_in = 9.5,
  dpi = 400
) {

  required_pkgs <- c("ggplot2", "dplyr", "sf", "patchwork", "rnaturalearth")
  missing <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      "Missing package(s): ", paste(missing, collapse = ", "),
      ". Install them manually before running this script."
    )
  }

  `%>%` <- dplyr::`%>%`

  if (!dir.exists(latitude_run)) stop("Latitude run not found: ", latitude_run)

  if (is.null(longitude_results)) {
    longitude_results <- file.path(
      project_root,
      "data",
      "longitude_sensitivity_results.csv"
    )
  }

  if (!file.exists(longitude_results)) stop("Longitude results not found: ", longitude_results)

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Register the requested Windows font name when possible.
  if (.Platform$OS.type == "windows") {
    try(grDevices::windowsFonts(Garamond = grDevices::windowsFont("Garamond")), silent = TRUE)
  }

  analyses <- c("onset", "offset_univoltine", "offset_multivoltine")

  titles <- c(
    onset = "Onset",
    offset_univoltine = "Offset — univoltine",
    offset_multivoltine = "Offset — multivoltine"
  )

  ylabels <- c(
    onset = "Relative onset advancement",
    offset_univoltine = "Relative offset advancement",
    offset_multivoltine = "Relative offset delay"
  )

  # Biological orientation:
  # onset/univoltine: raw negative = advancement -> multiply by -1
  # multivoltine: raw positive = delay -> retain sign
  orient <- c(
    onset = -1,
    offset_univoltine = -1,
    offset_multivoltine = 1
  )

  # ---------------------------------------------------------------------------
  # Coordinates
  # ---------------------------------------------------------------------------
  coord_files <- list.files(
    file.path(latitude_run, "inputs", "coordinates"),
    pattern = "\\.csv$",
    full.names = TRUE
  )
  if (!length(coord_files)) stop("No coordinate CSVs found in the latitude run.")

  coords_all <- dplyr::bind_rows(lapply(coord_files, function(f) {
    x <- utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
    need <- c("SITE_ID", "bms_id", "longitude_deg", "latitude_deg")
    if (!all(need %in% names(x))) stop("Coordinate file lacks required columns: ", f)
    x[, need]
  }))

  coords_all$SITE_ID <- as.character(coords_all$SITE_ID)
  coords_all$bms_id <- as.character(coords_all$bms_id)

  coord_check <- coords_all %>%
    dplyr::group_by(SITE_ID, bms_id) %>%
    dplyr::summarise(
      lon_range = max(longitude_deg) - min(longitude_deg),
      lat_range = max(latitude_deg) - min(latitude_deg),
      .groups = "drop"
    )

  if (any(coord_check$lon_range > 1e-7 | coord_check$lat_range > 1e-7)) {
    stop("Coordinate exports disagree for at least one site.")
  }

  coords <- coords_all %>%
    dplyr::distinct(SITE_ID, bms_id, longitude_deg, latitude_deg)

  # ---------------------------------------------------------------------------
  # Reviewed latitude + longitude coefficients
  # ---------------------------------------------------------------------------
  lr <- utils::read.csv(longitude_results, stringsAsFactors = FALSE, check.names = FALSE)

  lr <- lr %>%
    dplyr::filter(
      weighting == "equal_population",
      specification == "latitude_plus_longitude_within_species"
    )

  needed_lr <- c(
    "analysis", "term", "estimate_raw_per_10deg",
    "lower_95_raw", "upper_95_raw"
  )
  if (!all(needed_lr %in% names(lr))) stop("Unexpected longitude-results schema.")

  # ---------------------------------------------------------------------------
  # Basemap
  # ---------------------------------------------------------------------------
  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  world <- sf::st_crop(world, xmin = -12, xmax = 33, ymin = 38, ymax = 66)

  map_panels <- list()
  scatter_panels <- list()
  all_site_values <- list()
  all_partial_values <- list()

  for (a in analyses) {

    pop_file <- file.path(latitude_run, a, "population_plasticity_with_latitude.csv")
    if (!file.exists(pop_file)) stop("Missing: ", pop_file)

    p <- utils::read.csv(pop_file, stringsAsFactors = FALSE, check.names = FALSE)

    need <- c(
      "SPECIES", "SITE_ID", "bms_id",
      "plasticity_within_species_raw",
      "latitude_within_species_deg"
    )
    if (!all(need %in% names(p))) stop("Population table lacks required columns for ", a)

    p$SPECIES <- as.character(p$SPECIES)
    p$SITE_ID <- as.character(p$SITE_ID)
    p$bms_id <- as.character(p$bms_id)

    p <- dplyr::left_join(
      p, coords,
      by = c("SITE_ID", "bms_id"),
      relationship = "many-to-one"
    )

    if (anyNA(p$longitude_deg) || anyNA(p$latitude_deg)) {
      stop("Missing coordinates after join for ", a)
    }

    # Species-centred longitude.
    p <- p %>%
      dplyr::group_by(SPECIES) %>%
      dplyr::mutate(
        longitude_within_species_deg = longitude_deg - mean(longitude_deg)
      ) %>%
      dplyr::ungroup()

    x_lat <- p$latitude_within_species_deg / 10
    z_lon <- p$longitude_within_species_deg / 10
    y <- orient[[a]] * p$plasticity_within_species_raw

    this_lr <- lr[lr$analysis == a, , drop = FALSE]
    if (!all(c("latitude", "longitude") %in% this_lr$term)) {
      stop("Latitude/longitude coefficients missing for ", a)
    }

    get_row <- function(term) this_lr[this_lr$term == term, , drop = FALSE]

    lat_row <- get_row("latitude")
    lon_row <- get_row("longitude")

    # Orient coefficients in the same biological direction as y.
    b_lat <- orient[[a]] * lat_row$estimate_raw_per_10deg
    lat_ci <- sort(orient[[a]] * c(lat_row$lower_95_raw, lat_row$upper_95_raw))
    b_lon <- orient[[a]] * lon_row$estimate_raw_per_10deg

    # -------------------------------------------------------------------------
    # MAP
    #
    # Site-level colour is the average, across species present at the site, of:
    # biologically oriented within-species response minus the fitted longitude
    # component. Thus:
    #   positive = more advancement (onset/uni) or more delay (multi)
    #              than the species' sampled mean, conditional on longitude.
    # -------------------------------------------------------------------------
    p$map_response <- y - b_lon * z_lon

    site <- p %>%
      dplyr::group_by(SITE_ID, bms_id, longitude_deg, latitude_deg) %>%
      dplyr::summarise(
        relative_response = mean(map_response),
        n_species = dplyr::n_distinct(SPECIES),
        .groups = "drop"
      )

    site$analysis <- a
    all_site_values[[a]] <- site

    q <- as.numeric(stats::quantile(abs(site$relative_response), 0.98, na.rm = TRUE))
    if (!is.finite(q) || q <= 0) q <- max(abs(site$relative_response), na.rm = TRUE)

    p_map <- ggplot2::ggplot() +
      ggplot2::geom_sf(data = world, fill = "grey97", linewidth = 0.22) +
      ggplot2::geom_point(
        data = site,
        ggplot2::aes(
          x = longitude_deg,
          y = latitude_deg,
          colour = relative_response,
          size = n_species
        ),
        alpha = 0.80
      ) +
      ggplot2::scale_colour_gradient2(
        midpoint = 0,
        limits = c(-q, q),
        oob = scales::squish,
        name = paste0(ylabels[[a]], "\n(0 = species mean)")
      ) +
      ggplot2::scale_size_continuous(range = c(0.7, 2.4), guide = "none") +
      ggplot2::coord_sf(
        xlim = c(-12, 33),
        ylim = c(38, 66),
        expand = FALSE
      ) +
      ggplot2::labs(
        title = titles[[a]],
        x = NULL,
        y = NULL
      ) +
      ggplot2::theme_classic(base_family = font_family, base_size = 11) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          hjust = 0.5, face = "bold", size = 14
        ),
        legend.title = ggplot2::element_text(size = 9.5),
        legend.text = ggplot2::element_text(size = 8.5),
        legend.key.height = grid::unit(26, "mm"),
        axis.text = ggplot2::element_text(size = 8.5)
      )

    map_panels[[a]] <- p_map

    # -------------------------------------------------------------------------
    # BOTTOM PANEL: partial within-species latitude relationship
    #
    # Frisch-Waugh-Lovell visualization:
    # residualize both latitude and oriented response on within-species longitude.
    # The resulting line has the same latitude coefficient as the simultaneous
    # latitude + longitude projection.
    # -------------------------------------------------------------------------
    denom_z <- sum(z_lon^2)
    if (!is.finite(denom_z) || denom_z <= 0) stop("No longitude variation for ", a)

    lat_res <- x_lat - z_lon * sum(x_lat * z_lon) / denom_z
    y_res <- y - z_lon * sum(y * z_lon) / denom_z

    partial <- data.frame(
      analysis = a,
      SPECIES = p$SPECIES,
      SITE_ID = p$SITE_ID,
      latitude_partial_deg = lat_res * 10,
      relative_response_partial = y_res
    )
    all_partial_values[[a]] <- partial

    xseq <- seq(
      stats::quantile(partial$latitude_partial_deg, 0.005, na.rm = TRUE),
      stats::quantile(partial$latitude_partial_deg, 0.995, na.rm = TRUE),
      length.out = 250
    )

    line <- data.frame(
      x = xseq,
      fit = b_lat * xseq / 10,
      low = lat_ci[1] * xseq / 10,
      high = lat_ci[2] * xseq / 10
    )
    # CI lines cross at x=0; keep ymin/ymax ordered.
    line$ymin <- pmin(line$low, line$high)
    line$ymax <- pmax(line$low, line$high)

    label <- sprintf(
      "+%.2f [%.2f, %.2f] per 10°",
      b_lat, lat_ci[1], lat_ci[2]
    )

    p_scatter <- ggplot2::ggplot(
      partial,
      ggplot2::aes(x = latitude_partial_deg, y = relative_response_partial)
    ) +
      ggplot2::geom_hline(yintercept = 0, linetype = 3, linewidth = 0.35) +
      ggplot2::geom_vline(xintercept = 0, linetype = 3, linewidth = 0.35) +
      ggplot2::geom_point(alpha = 0.12, size = 0.55) +
      ggplot2::geom_ribbon(
        data = line,
        ggplot2::aes(x = x, ymin = ymin, ymax = ymax),
        inherit.aes = FALSE,
        alpha = 0.16
      ) +
      ggplot2::geom_line(
        data = line,
        ggplot2::aes(x = x, y = fit),
        inherit.aes = FALSE,
        linewidth = 0.9
      ) +
      ggplot2::annotate(
        "label",
        x = -Inf, y = Inf,
        label = label,
        hjust = -0.05, vjust = 1.15,
        family = font_family,
        size = 3.5,
        label.size = 0.25
      ) +
      ggplot2::labs(
        title = titles[[a]],
        x = "Latitude within species (°)",
        y = ylabels[[a]]
      ) +
      ggplot2::theme_classic(base_family = font_family, base_size = 11) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          hjust = 0.5, face = "bold", size = 14
        ),
        axis.title = ggplot2::element_text(size = 11),
        axis.text = ggplot2::element_text(size = 9)
      )

    scatter_panels[[a]] <- p_scatter
  }

  # ---------------------------------------------------------------------------
  # Assemble figure
  # ---------------------------------------------------------------------------
  top <- map_panels[[1]] + map_panels[[2]] + map_panels[[3]]
  bottom <- scatter_panels[[1]] + scatter_panels[[2]] + scatter_panels[[3]]

  fig <- (top / bottom) +
    patchwork::plot_annotation(
      tag_levels = "a",
      theme = ggplot2::theme(
        plot.tag = ggplot2::element_text(
          family = font_family, face = "bold", size = 12
        )
      )
    )

  png_file <- file.path(output_dir, "manuscript_latitude_longitude_patterns.png")
  pdf_file <- file.path(output_dir, "manuscript_latitude_longitude_patterns.pdf")

  ggplot2::ggsave(
    png_file, fig,
    width = width_in, height = height_in,
    units = "in", dpi = dpi, bg = "white"
  )

  ggplot2::ggsave(
    pdf_file, fig,
    width = width_in, height = height_in,
    units = "in", device = grDevices::cairo_pdf, bg = "white"
  )

  utils::write.csv(
    dplyr::bind_rows(all_site_values),
    file.path(output_dir, "figure_site_values.csv"),
    row.names = FALSE
  )

  utils::write.csv(
    dplyr::bind_rows(all_partial_values),
    file.path(output_dir, "figure_partial_values.csv"),
    row.names = FALSE
  )

  caption_note <- paste(
    "Map colours represent population-level deviations from the sampled species mean",
    "in biologically oriented predicted thermal responses after accounting for",
    "within-species longitude. Positive values indicate greater onset advancement",
    "for onset, greater offset advancement for univoltines, and greater offset delay",
    "for multivoltines. Bottom panels show the corresponding partial within-species",
    "latitudinal relationships. Latitude effects and 95% intervals come from the",
    "reviewed source-model covariance projections; no phenology model or spatial",
    "field is refitted here."
  )

  writeLines(
    c(
      paste("PNG:", normalizePath(png_file, winslash = "/", mustWork = TRUE)),
      paste("PDF:", normalizePath(pdf_file, winslash = "/", mustWork = TRUE)),
      "",
      "Suggested caption note:",
      caption_note
    ),
    file.path(output_dir, "figure_notes.txt")
  )

  message("Figure saved:")
  message("  ", png_file)
  message("  ", pdf_file)

  invisible(list(
    figure = fig,
    png = png_file,
    pdf = pdf_file,
    site_values = dplyr::bind_rows(all_site_values),
    partial_values = dplyr::bind_rows(all_partial_values),
    caption_note = caption_note
  ))
}

# Example:
# fig_lat <- make_latitude_longitude_figure()
