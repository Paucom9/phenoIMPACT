# ============================================================================ #
# phenoIMPACT | Population trends with the 15-km spatial correction
# Reproduce Pau's four reference figures in R / ggplot2.
#
# Input:  univoltine__full/predictions_spatial.csv
#         multivoltine__full/predictions_spatial.csv
# Output: four PNG figures (7 x 5 inches, 300 dpi), editable ggplot objects,
#         and copies of the plotted predictions / input provenance.
#
# NO MODEL IS LOADED OR REFITTED. The original scripts are NOT sourced.
# The saved predictions and their confidence limits are used without alteration.
# Required: ggplot2 >= 3.4.0. Optional: ragg for PNG rendering.
# Garamond must already be installed (as for the reference figures).
#
# These are fixed-effect relative trajectories from the spatially fitted models:
#   100 * [exp((year - 1990)/10 * (beta_year + z * beta_year:plasticity)) - 1].
# The other plasticity predictor is held at zero (its standardized mean).
# Species/site, population AR(1), and annual spatial effects are set to zero;
# they are NOT plotted as realized site-specific trajectories or marginalized.
# Ribbons are the exported pointwise 95% confidence intervals based on the
# fixed-effect covariance, NOT prediction intervals for new populations.
# ============================================================================ #

plot_population_trends_spatial <- function(
    run_dir = paste0(
      "E:/phenoIMPACT project/code/phenoIMPACT/output/population_trends/",
      "spatial_abundance_15km/run_650c41f15949"),
    out_dir = file.path(run_dir, "figures_R_original_style"),
    font_family = "Garamond",
    base_size = 16,
    width = 7,
    height = 5,
    dpi = 300,
    show_titles = TRUE,
    show_plots = interactive()) {

  # ---- Checks and fonts ---------------------------------------------------- #
  fail <- function(...) stop(..., call. = FALSE)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    fail("Package 'ggplot2' is required. Install it, then source this script again.")
  }
  if (utils::packageVersion("ggplot2") < "3.4.0") {
    fail("This script requires ggplot2 >= 3.4.0 (linewidth support).")
  }
  if (!dir.exists(run_dir)) fail("Run directory not found: ", run_dir)
  if (!is.character(font_family) || length(font_family) != 1L ||
      is.na(font_family) || !nzchar(font_family)) fail("Invalid font_family.")
  sizes <- c(base_size, width, height, dpi)
  if (length(sizes) != 4L || any(!is.finite(sizes)) || any(sizes <= 0)) {
    fail("base_size, width, height and dpi must each be a positive number.")
  }
  if (!is.logical(show_titles) || length(show_titles) != 1L || is.na(show_titles) ||
      !is.logical(show_plots) || length(show_plots) != 1L || is.na(show_plots)) {
    fail("show_titles and show_plots must be TRUE or FALSE.")
  }

  # Check the actual installed family when systemfonts is available; do not
  # silently substitute another font. No fonts or packages are installed here.
  if (requireNamespace("systemfonts", quietly = TRUE)) {
    installed <- unique(systemfonts::system_fonts()$family)
    if (!tolower(font_family) %in% tolower(installed)) {
      fail("Font '", font_family, "' is not installed. Set font_family to an ",
           "installed family, or install your original Garamond before rerunning.")
    }
  }
  if (.Platform$OS.type == "windows") {
    # This also supports the Windows graphics device when ragg is unavailable.
    do.call(grDevices::windowsFonts,
            stats::setNames(list(grDevices::windowsFont(font_family)), font_family))
  }

  # ---- Read ONLY the spatial prediction exports ---------------------------- #
  groups <- c("univoltine", "multivoltine")
  input_files <- stats::setNames(
    file.path(run_dir, paste0(groups, "__full"), "predictions_spatial.csv"), groups)
  if (any(!file.exists(input_files))) {
    fail("Missing spatial prediction CSV(s):\n",
         paste(input_files[!file.exists(input_files)], collapse = "\n"))
  }

  read_predictions <- function(path) {
    d <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    numeric_cols <- c("plasticity_z", "YEAR", "percent_change", "low", "high")
    required <- c("variable", numeric_cols)
    if (!all(required %in% names(d)) || !nrow(d)) {
      fail("Unexpected columns or empty prediction export: ", path)
    }
    if (anyNA(d[required]) ||
        !all(vapply(d[numeric_cols], function(x) {
          is.numeric(x) && all(is.finite(x))
        }, logical(1)))) {
      fail("Non-numeric, missing or non-finite prediction values: ", path)
    }
    variables <- c("onset_plasticity_z", "offset_plasticity_bio_z")
    if (!setequal(unique(d$variable), variables) ||
        !setequal(unique(d$plasticity_z), c(-1, 0, 1)) ||
        anyDuplicated(d[c("variable", "plasticity_z", "YEAR")])) {
      fail("Unexpected predictor levels or duplicate prediction rows: ", path)
    }
    if (any(d$YEAR != floor(d$YEAR)) ||
        any(d$low > d$percent_change + 1e-8) ||
        any(d$high < d$percent_change - 1e-8)) {
      fail("Invalid years or confidence limits: ", path)
    }
    years <- sort(unique(d$YEAR))
    if (length(years) < 2L) fail("At least two years are required: ", path)
    for (v in variables) for (z in c(-1, 0, 1)) {
      a <- d[d$variable == v & d$plasticity_z == z, , drop = FALSE]
      if (!identical(sort(a$YEAR), years)) {
        fail("Prediction trajectories do not have the same years: ", path)
      }
      first <- a[a$YEAR == min(years), c("percent_change", "low", "high")]
      if (any(abs(as.matrix(first)) > 1e-8)) {
        fail("The prediction export is not relative to its first year: ", path)
      }
    }
    d[order(d$variable, -d$plasticity_z, d$YEAR), , drop = FALSE]
  }

  predictions <- lapply(input_files, read_predictions)
  if (!identical(sort(unique(predictions[[1]]$YEAR)),
                 sort(unique(predictions[[2]]$YEAR)))) {
    fail("The two groups have different year ranges; review before plotting.")
  }

  # ---- The same titles, colours and legend as the supplied images ---------- #
  jobs <- data.frame(
    key = c("onset_univoltine", "offset_univoltine",
            "onset_multivoltine", "offset_multivoltine"),
    group = c("univoltine", "univoltine", "multivoltine", "multivoltine"),
    variable = c("onset_plasticity_z", "offset_plasticity_bio_z",
                 "onset_plasticity_z", "offset_plasticity_bio_z"),
    title = c(
      "Emergence plasticity: univoltine populations",
      "Offset advancement plasticity: univoltine populations",
      "Onset advancement plasticity: multivoltine populations",
      "Offset delay plasticity: multivoltine populations"),
    stringsAsFactors = FALSE
  )
  legend_levels <- c("High (+1 SD)", "Mean", "Low (-1 SD)")
  palette <- stats::setNames(c("#009E73", "#0072B2", "#D55E00"), legend_levels)

  make_plot <- function(d, variable, title) {
    d <- d[d$variable == variable, , drop = FALSE]
    # offset_plasticity_bio_z is ALREADY biologically oriented. Do not flip it:
    # +1 SD = stronger advancement for univoltines, more delay for multivoltines.
    d$Plasticity <- factor(d$plasticity_z,
                          levels = c(1, 0, -1), labels = legend_levels)
    ggplot2::ggplot(
      d, ggplot2::aes(x = YEAR, y = percent_change,
                     colour = Plasticity, fill = Plasticity, group = Plasticity)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = low, ymax = high),
                           alpha = 0.10, colour = NA) +
      ggplot2::geom_line(linewidth = 1.4) +
      ggplot2::scale_colour_manual(name = "Plasticity", values = palette,
                                   breaks = legend_levels, drop = FALSE) +
      ggplot2::scale_fill_manual(name = "Plasticity", values = palette,
                                 breaks = legend_levels, drop = FALSE) +
      ggplot2::scale_x_continuous(
        breaks = function(limits) {
          seq(ceiling(limits[1] / 10) * 10, floor(limits[2] / 10) * 10, by = 10)
        }, expand = ggplot2::expansion(mult = 0.05)) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.05)) +
      ggplot2::labs(
        title = if (show_titles) title else NULL,
        x = "Year", y = "Relative change in abundance (%)") +
      ggplot2::theme_classic(base_family = font_family, base_size = base_size) +
      ggplot2::theme(
        legend.position = "right",
        legend.direction = "vertical",
        plot.title = ggplot2::element_text(hjust = 0, face = "plain"))
    # Intentionally no ylim(): preserve the complete spatial-model intervals.
    # In particular, do not reuse the negative-only scales of the old figures.
  }

  # ---- Export figures; leave input files and fitted models untouched -------- #
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(out_dir)) fail("Cannot create output folder: ", out_dir)
  use_ragg <- requireNamespace("ragg", quietly = TRUE)
  png_device <- if (use_ragg) ragg::agg_png else grDevices::png
  device_args <- if (use_ragg) list() else if (.Platform$OS.type == "windows") {
    list(type = "windows")
  } else if (capabilities("cairo")) {
    list(type = "cairo")
  } else list()

  plots <- stats::setNames(vector("list", nrow(jobs)), jobs$key)
  for (i in seq_len(nrow(jobs))) {
    job <- jobs[i, ]
    p <- make_plot(predictions[[job$group]], job$variable, job$title)
    filename <- file.path(out_dir,
                          paste0("plot_trend_", job$key, "_gamma_spatial_15km.png"))
    do.call(ggplot2::ggsave, c(list(
      filename = filename, plot = p, device = png_device,
      width = width, height = height, units = "in", dpi = dpi, bg = "white"),
      device_args))
    plots[[job$key]] <- p
    if (show_plots) print(p)
    message("Saved: ", basename(filename))
  }

  saveRDS(plots, file.path(out_dir, "plots_spatial_15km.rds"))
  plotted <- do.call(rbind, lapply(groups, function(g) {
    data.frame(group = g, predictions[[g]], row.names = NULL)
  }))
  utils::write.csv(plotted, file.path(out_dir, "predictions_used.csv"), row.names = FALSE)
  utils::write.csv(data.frame(
    group = groups, source_file = unname(input_files),
    md5 = unname(tools::md5sum(input_files)), font_family = font_family,
    base_size = base_size, width_in = width, height_in = height, dpi = dpi,
    png_device = if (use_ragg) "ragg::agg_png" else "grDevices::png",
    ggplot2_version = as.character(utils::packageVersion("ggplot2"))),
    file.path(out_dir, "figure_manifest.csv"), row.names = FALSE)

  message("DONE. Four figures saved. No fitted models were loaded or refitted.")
  message("Output: ", normalizePath(out_dir, winslash = "/", mustWork = TRUE))
  invisible(list(plots = plots, predictions = plotted, output_dir = out_dir))
}

# ---- Run with the defaults -------------------------------------------------- #
# In RStudio: source(file.choose(), encoding = "UTF-8")
#
# To define the function WITHOUT executing it, first run:
# options(phenoimpact.spatial_plots_functions_only = TRUE)
# Then source this file and call, for example:
# population_trends_spatial_figures <- plot_population_trends_spatial(
#   font_family = "Garamond", show_titles = FALSE)
#
# Plot objects remain available for editing, for example:
# population_trends_spatial_figures$plots$onset_univoltine
#
if (!isTRUE(getOption("phenoimpact.spatial_plots_functions_only", FALSE))) {
  population_trends_spatial_figures <- plot_population_trends_spatial()
}
