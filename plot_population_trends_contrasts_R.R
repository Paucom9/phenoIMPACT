# ============================================================================ #
# phenoIMPACT | Spatial population-trend models (15-km mesh)
# High (+1 SD) minus low (-1 SD) plasticity: contrasts of temporal slopes.
#
# Reads ONLY: <group>__full/fixed_effects.csv, for each voltinism group.
# Does NOT load/refit models, source other scripts, or install packages/fonts.
#
# For the linear predictor used by the fitted models:
#   slope(z) = beta_year + z * beta_year:plasticity
#   slope(+1) - slope(-1) = 2 * beta_year:plasticity
# Therefore estimate, standard error and 95% CI limits are multiplied by 2.
# The shared beta_year cancels; its variance must NOT be added to this contrast.
# This is NOT a subtraction of the two trajectory confidence bands.
#
# Units: difference in log-abundance slope per decade; NOT percentage points.
# Other plasticity is held constant. High/low use within-group standardized units.
# High onset = greater advancement in both groups.
# High offset = greater advancement (univoltines), more delay (multivoltines).
# The exported intervals are individual 95% Wald CIs, not simultaneous CIs.
# This does not test differences between voltinism groups or quantify explained
# variance. It plots the fitted associations, not causal effects.
#
# Run in RStudio: source(file.choose(), encoding = "UTF-8")
# Graphics reference: https://ggplot2.tidyverse.org/reference/geom_linerange.html
# Export reference:   https://ggplot2.tidyverse.org/reference/ggsave.html
# ============================================================================ #

plot_population_trends_contrasts <- function(
    run_dir = paste0(
      "E:/phenoIMPACT project/code/phenoIMPACT/output/population_trends/",
      "spatial_abundance_15km/run_650c41f15949"),
    out_dir = file.path(run_dir, "figures_plasticity_contrasts"),
    font_family = "Garamond",
    base_size = 16,
    width = 7,
    height = 5,
    dpi = 300,
    group_colours = c(univoltine = "#0072B2", multivoltine = "#D55E00"),
    show_plot = interactive()) {

  fail <- function(...) stop(..., call. = FALSE)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    fail("Package 'ggplot2' is required. No model has been touched.")
  }
  if (utils::packageVersion("ggplot2") < "3.4.0") {
    fail("ggplot2 >= 3.4.0 is required for linewidth support.")
  }
  if (!is.character(run_dir) || length(run_dir) != 1L || is.na(run_dir) ||
      !dir.exists(run_dir)) fail("Run directory not found: ", run_dir)
  if (!is.character(out_dir) || length(out_dir) != 1L || is.na(out_dir) ||
      !nzchar(out_dir)) fail("out_dir must be a non-empty path.")
  if (!is.character(font_family) || length(font_family) != 1L ||
      is.na(font_family) || !nzchar(font_family)) fail("Invalid font_family.")
  for (nm in c("base_size", "width", "height", "dpi")) {
    v <- get(nm)
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v <= 0) {
      fail(nm, " must be one finite positive number.")
    }
  }
  if (!is.logical(show_plot) || length(show_plot) != 1L || is.na(show_plot)) {
    fail("show_plot must be TRUE or FALSE.")
  }
  groups <- c("univoltine", "multivoltine")
  if (!is.character(group_colours) || length(group_colours) != 2L ||
      anyNA(group_colours) || anyDuplicated(names(group_colours)) ||
      !setequal(names(group_colours), groups)) {
    fail("group_colours must contain named colours for univoltine and multivoltine.")
  }
  tryCatch(grDevices::col2rgb(group_colours), error = function(e) {
    fail("Invalid group colour: ", conditionMessage(e))
  })

  # Do not silently replace Garamond if a font inventory is available.
  if (requireNamespace("systemfonts", quietly = TRUE)) {
    installed <- unique(systemfonts::system_fonts()$family)
    if (!tolower(font_family) %in% tolower(installed)) {
      fail("Font '", font_family, "' is not installed. Use your original Garamond ",
           "installation or set font_family explicitly to an installed family.")
    }
  }
  if (.Platform$OS.type == "windows") {
    do.call(grDevices::windowsFonts,
      stats::setNames(list(grDevices::windowsFont(font_family)), font_family))
  }

  # Source files are explicitly spatial-model exports, not the baseline tables.
  input_files <- stats::setNames(
    file.path(run_dir, paste0(groups, "__full"), "fixed_effects.csv"), groups)
  if (any(!file.exists(input_files))) {
    fail("Missing spatial fixed-effect export(s):\n",
         paste(input_files[!file.exists(input_files)], collapse = "\n"))
  }
  source_md5 <- unname(tools::md5sum(input_files))

  read_effects <- function(path, group) {
    d <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    needed <- c("term", "estimate", "std.error", "conf.low", "conf.high", "p_Wald")
    if (!nrow(d) || !all(needed %in% names(d))) {
      fail("Unexpected fixed-effect columns or empty export: ", path)
    }
    if (anyNA(d$term) || anyDuplicated(d$term)) {
      fail("Missing or duplicated coefficient names: ", path)
    }
    variables <- c(onset = "onset_plasticity_z", offset = "offset_plasticity_bio_z")
    ans <- lapply(names(variables), function(component) {
      v <- unname(variables[component])
      candidates <- c(paste0("year_decade:", v), paste0(v, ":year_decade"))
      a <- d[d$term %in% candidates, needed, drop = FALSE]
      if (nrow(a) != 1L) fail("Expected exactly one year interaction for ", v, " in ", path)
      nums <- needed[needed != "term"]
      if (!all(vapply(a[nums], function(x) {
        is.numeric(x) && length(x) == 1L && is.finite(x)
      }, logical(1)))) fail("Missing or non-finite interaction values: ", path)
      if (a$std.error <= 0 || a$conf.low > a$estimate ||
          a$conf.high < a$estimate || a$p_Wald < 0 || a$p_Wald > 1) {
        fail("Invalid interaction estimate, interval, SE or p-value: ", path)
      }
      # Confirm a 95% normal Wald interval, accepting the source's rounded 1.96.
      # Preserve the exported limits exactly rather than recomputing them.
      z_limits <- c(a$estimate - a$conf.low, a$conf.high - a$estimate) / a$std.error
      if (any(abs(z_limits - stats::qnorm(0.975)) > 1e-4)) {
        fail("Input interval is not the expected 95% Wald interval: ", path,
             ". Review the source before plotting; no interval was replaced.")
      }
      data.frame(
        group = group, component = component, source_term = a$term,
        beta_interaction = a$estimate, SE_interaction = a$std.error,
        interaction_CI_low = a$conf.low, interaction_CI_high = a$conf.high,
        plasticity_high_z = 1, plasticity_low_z = -1,
        contrast = 2 * a$estimate,
        contrast_SE = 2 * a$std.error,
        contrast_CI_low = 2 * a$conf.low,
        contrast_CI_high = 2 * a$conf.high,
        p_Wald = a$p_Wald,
        source_file = normalizePath(path, winslash = "/", mustWork = TRUE),
        stringsAsFactors = FALSE)
    })
    do.call(rbind, ans)
  }
  contrasts <- do.call(rbind, lapply(groups, function(g) {
    read_effects(input_files[[g]], g)
  }))
  # Four separate rows: keep each offset's biological direction explicit.
  order_keys <- c("onset_univoltine", "onset_multivoltine",
                  "offset_univoltine", "offset_multivoltine")
  idx <- match(order_keys, paste(contrasts$component, contrasts$group, sep = "_"))
  if (anyNA(idx) || anyDuplicated(idx)) fail("Unexpected contrast layout.")
  contrasts <- contrasts[idx, , drop = FALSE]
  rownames(contrasts) <- NULL
  contrasts$row_y <- c(4, 3, 1.75, 0.75)
  contrasts$row_label <- c(
    "Onset advancement\nUnivoltine populations",
    "Onset advancement\nMultivoltine populations",
    "Offset advancement\nUnivoltine populations",
    "Offset delay\nMultivoltine populations")

  # A symmetric axis keeps positive and negative contrasts visually comparable.
  # The limit is derived from the full intervals; none are cropped.
  limit <- 1.10 * max(abs(c(contrasts$contrast_CI_low, contrasts$contrast_CI_high)))
  p <- ggplot2::ggplot(contrasts, ggplot2::aes(
    x = contrast, y = row_y, colour = group, shape = group)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                       colour = "grey40", linewidth = 0.5) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = contrast_CI_low, xmax = contrast_CI_high),
      orientation = "y", width = 0.14, linewidth = 0.85) +
    ggplot2::geom_point(size = 3.6) +
    ggplot2::scale_colour_manual(values = group_colours, breaks = groups) +
    ggplot2::scale_shape_manual(values = c(univoltine = 16, multivoltine = 17),
                                breaks = groups) +
    ggplot2::scale_x_continuous(
      limits = c(-limit, limit),
      breaks = function(limits) pretty(limits, n = 6),
      labels = function(x) formatC(x, format = "f", digits = 2),
      expand = ggplot2::expansion(mult = 0)) +
    ggplot2::scale_y_continuous(
      breaks = contrasts$row_y, labels = contrasts$row_label,
      expand = ggplot2::expansion(add = 0.55)) +
    ggplot2::labs(
      x = "Difference in trend: high - low plasticity\n(log-abundance per decade)",
      y = NULL) +
    ggplot2::theme_classic(base_family = font_family, base_size = base_size) +
    ggplot2::theme(
      legend.position = "none",  # Each row already states its voltinism group.
      axis.text.y = ggplot2::element_text(size = base_size * 0.81,
                                        colour = "black", lineheight = 0.95),
      axis.text.x = ggplot2::element_text(size = base_size * 0.81),
      axis.title.x = ggplot2::element_text(size = base_size * 0.90,
                                         margin = ggplot2::margin(t = 10)),
      axis.ticks.y = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(10, 12, 10, 10))

  # Export only to this plotting subfolder. Input files and fits remain untouched.
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(out_dir)) fail("Cannot create output folder: ", out_dir)
  png_file <- file.path(out_dir, "plot_plasticity_trend_contrasts_spatial_15km.png")
  csv_file <- file.path(out_dir, "plasticity_trend_contrasts_spatial_15km.csv")
  rds_file <- file.path(out_dir, "plot_plasticity_trend_contrasts_spatial_15km.rds")
  use_ragg <- requireNamespace("ragg", quietly = TRUE)
  device <- if (use_ragg) ragg::agg_png else grDevices::png
  extra <- if (use_ragg) list() else if (.Platform$OS.type == "windows") {
    list(type = "windows")
  } else if (capabilities("cairo")) list(type = "cairo") else list()
  do.call(ggplot2::ggsave, c(list(filename = png_file, plot = p, device = device,
    width = width, height = height, units = "in", dpi = dpi, bg = "white"), extra))
  utils::write.csv(contrasts, csv_file, row.names = FALSE)
  saveRDS(p, rds_file)
  utils::write.csv(data.frame(
    group = groups, source_file = unname(input_files), md5 = source_md5,
    font_family = font_family, base_size = base_size,
    width_in = width, height_in = height, dpi = dpi,
    ggplot2_version = as.character(utils::packageVersion("ggplot2"))),
    file.path(out_dir, "figure_manifest.csv"), row.names = FALSE)
  writeLines(c(
    "SPATIAL POPULATION-TREND CONTRASTS: HIGH MINUS LOW PLASTICITY",
    "High = +1 SD; low = -1 SD, using the model's within-group standardization.",
    "Contrast = 2 * beta(year_decade:plasticity); SE and CI limits also multiplied by 2.",
    "Units: difference in temporal log-abundance slope per decade, not percentage points.",
    "The other plasticity is held constant; the shared temporal coefficient cancels.",
    "Positive values: a more positive slope at high plasticity; negative values: the reverse.",
    "High onset means stronger advancement in both voltinism groups.",
    "High offset means stronger advancement in univoltines, more delay in multivoltines.",
    "Intervals are the transformed individual 95% Wald CIs from the spatial fit exports.",
    "p_Wald is unchanged by the positive multiplier 2; no multiplicity adjustment is added.",
    "No between-voltinism test or explained-variance comparison is performed.",
    "No coefficients are hard-coded. No trajectories, random effects or models are loaded.",
    "No model is refitted. No confidence interval is cropped. No causal effect is estimated.",
    "Rows and groups are identified in the y labels; colour and shape are redundant cues.",
    "The PNG has no title/caption. The RDS contains the editable ggplot object."
  ), file.path(out_dir, "README.txt"))
  if (!identical(unname(tools::md5sum(input_files)), source_md5)) {
    fail("Input files changed during plotting. Review the input files and rerun.")
  }
  if (show_plot) print(p)
  message("DONE. Four high-minus-low contrasts; no model was loaded or refitted.")
  message("Output: ", normalizePath(out_dir, winslash = "/", mustWork = TRUE))
  print(contrasts[c("group", "component", "contrast", "contrast_CI_low",
                    "contrast_CI_high", "p_Wald")], row.names = FALSE)
  invisible(list(plot = p, contrasts = contrasts,
                 files = c(png = png_file, csv = csv_file, rds = rds_file),
                 output_dir = out_dir))
}

# Source to run with defaults. For custom settings, define without running first:
# options(phenoimpact.contrast_plots_functions_only = TRUE)
# source(file.choose(), encoding = "UTF-8")
# population_trends_contrasts <- plot_population_trends_contrasts(width = 7.5)
# The returned plot is editable: population_trends_contrasts$plot
if (!isTRUE(getOption("phenoimpact.contrast_plots_functions_only", FALSE))) {
  population_trends_contrasts <- plot_population_trends_contrasts()
}
