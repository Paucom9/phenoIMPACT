# Phenological plasticity: spatial models, 15-km mesh, 28 September 2026
# Reproduces the four principal curves in slides 29-32 of the original deck.
# Uses the verified saved estimates; this script does NOT fit any model.
# The multivoltine result comes from the numerically recovered checkpoint.
#
# Usage in RStudio: source(file.choose())
# Default output: output/phenology_plasticity/spatial_models/figures_15km
# Optional output override, set BEFORE sourcing:
# options(pheno.spatial.figures.output_dir = "E:/my_folder")
# Original Windows font: Garamond (as in pheno_plasticity.R).
# Optional font override: options(pheno.spatial.figures.font = "EB Garamond")
# Load functions only: options(pheno.spatial.figures.functions_only = TRUE)
# Then: pheno_spatial_figures_15km$run(output_dir = "E:/my_folder")
#
# Conditional fixed-effect slope: b_anomaly + x * b_interaction.
# Other moderators are fixed at zero in their inherited model scales.
# 95% pointwise normal Wald intervals include coefficient covariance.
# Solid/dashed lines encode the two-sided, unadjusted normal Wald test of
# the anomaly x moderator interaction (p < 0.05 / p >= 0.05).
# Curves and ribbons stop at the observed marginal predictor range.
# Curves are not species-specific random slopes or population averages.
# The 15-km mesh cutoff is not a spatial correlation range.
# Only the selected windows are used: onset 60 d; all others 90 d.
# See README.md for source, covariance verification and interpretation.

pheno_spatial_figures_15km <- local({
  parameters <- utils::read.csv(text = 'analysis,window,phenology_label,moderator,b_anomaly,b_interaction,var_anomaly,var_interaction,cov_anomaly_interaction,observed_min,observed_max,central_min,central_max
onset,60,onset,photoperiod,-2.50677783497258,-1.41776076904681,0.024379709803687925,0.008504659668175098,-0.0012364564272870475,-3.24099396995394,2.88259718512563,-2.33712121816527,1.04888098752902
onset,60,onset,background,-2.50677783497258,0.785294191956072,0.024379709803687925,0.013126696094430991,0.003978141293057231,-4.95632010142018,2.87475452837388,-2.16068761430988,0.656741759629387
onset,60,onset,predictability,-2.50677783497258,-0.331845979864984,0.024379709803687925,0.00241977498949158,-0.0004200689577234689,-13.1344471113885,1.65105553314022,-2.28521805135259,0.973786484627663
onset,60,onset,trend,-2.50677783497258,-0.0289865226028226,0.024379709803687925,0.0028813543540538375,-0.0005670446222822166,-5.62914985896075,6.82798394012212,-1.66924381256895,1.84769020180517
first_peak,90,first peak,photoperiod,-2.34860583413151,-1.64829582327738,0.03993283013342476,0.026529331885822682,-0.005336603128731513,-2.58950311386539,2.74048378695558,-1.41123030357036,1.21018342165281
first_peak,90,first peak,background,-2.34860583413151,0.272100271870992,0.03993283013342476,0.02657706852207545,0.0026089573583493527,-3.70360765462454,2.9136819826796,-1.41139637967442,1.19869427659318
first_peak,90,first peak,predictability,-2.34860583413151,-0.538349633840222,0.03993283013342476,0.01253663957846276,0.0007298071055016481,-6.8458710821459,1.56894167795788,-1.65305683379563,1.08561610050657
first_peak,90,first peak,trend,-2.34860583413151,-0.158637102155111,0.03993283013342476,0.007685902136066169,-0.00267125194601246,-4.97545087805525,5.63520004994828,-1.76186975647399,2.04518424848898
offset_univoltine,90,offset,photoperiod,-1.97827888469861,-0.615770630439595,0.06628199442095603,0.03026261308505221,-0.012984982137103335,-2.16920630988667,2.76037165283533,-0.509314514079218,1.57306833977492
offset_univoltine,90,offset,background,-1.97827888469861,-0.271058116649808,0.06628199442095603,0.03002901954402982,-0.012037379288294103,-3.03739979552026,2.95344772355418,-0.761330429499775,1.62346144547838
offset_univoltine,90,offset,predictability,-1.97827888469861,-0.751726397291479,0.06628199442095603,0.019900178246035977,-0.002221954435901757,-4.09909585764332,1.58562488147312,-1.2147523727541,1.26455497723192
offset_univoltine,90,offset,trend,-1.97827888469861,0.180781876571304,0.06628199442095603,0.006919514709789343,-0.0017431307724377633,-4.54372427975082,5.8027018157789,-2.02685110993499,2.06330402234141
offset_multivoltine,90,offset,photoperiod,0.3818912170617485,1.290535986969837,0.34505269999166965,0.06645673345936694,0.002687711415519701,-2.04043900864663,2.75468670533953,-0.737814265194715,1.00724149254686
offset_multivoltine,90,offset,background,0.3818912170617485,-1.352521928733278,0.34505269999166965,0.08386950984872761,-0.0856756152803348,-2.51202363276738,2.98055982746189,0.290881797789844,2.44470614481174
offset_multivoltine,90,offset,predictability,0.3818912170617485,-0.1465918260733355,0.34505269999166965,0.08604205819322927,0.038766569055628176,-3.26180841471071,1.5843811748522,-1.09895225698758,1.24303642592404
offset_multivoltine,90,offset,trend,0.3818912170617485,-0.373923710434463,0.34505269999166965,0.015179886666848314,0.003692760908067836,-3.28808120708399,5.8027018157789,-1.8260804307235,2.864745114857839
',
                                stringsAsFactors = FALSE, check.names = FALSE)
  model_order <- c("onset", "first_peak", "offset_univoltine", "offset_multivoltine")
  moderator_order <- c("photoperiod", "background", "predictability", "trend")
  labels <- c(photoperiod = "Photoperiod", background = "Mean temperature",
              predictability = "Temperature predictability", trend = "Temperature trend")
  colours <- c("Photoperiod" = "#1b9e77", "Mean temperature" = "#E6AB02",
               "Temperature predictability" = "#7570b3", "Temperature trend" = "#d95f02")
  display_names <- c(onset = "Onset", first_peak = "First peak",
                    offset_univoltine = "Offset: univoltine",
                    offset_multivoltine = "Offset: multivoltine")

  validate_parameters <- function(p) {
    stopifnot(nrow(p) == 16L, !anyDuplicated(p[c("analysis", "moderator")]),
              setequal(p$analysis, model_order), setequal(p$moderator, moderator_order))
    numerics <- c("b_anomaly", "b_interaction", "var_anomaly", "var_interaction",
                  "cov_anomaly_interaction", "observed_min", "observed_max",
                  "central_min", "central_max")
    stopifnot(all(vapply(p[numerics], function(x) all(is.finite(x)), logical(1))),
              all(p$var_anomaly > 0), all(p$var_interaction > 0),
              all(p$var_anomaly * p$var_interaction > p$cov_anomaly_interaction^2),
              all(p$observed_min < p$observed_max))
    for (m in model_order) {
      z <- p[p$analysis == m, ]
      stopifnot(nrow(z) == 4L, length(unique(z$b_anomaly)) == 1L,
                length(unique(z$var_anomaly)) == 1L,
                all(z$window == if (m == "onset") 60L else 90L))
    }
    invisible(TRUE)
  }

  make_curves <- function(p = parameters, range = c("original", "central")) {
    validate_parameters(p)
    range <- match.arg(range)
    rows <- lapply(seq_len(nrow(p)), function(i) {
      r <- p[i, ]
      limits <- if (range == "original") c(-2, 2) else c(r$central_min, r$central_max)
      limits <- c(max(limits[1], r$observed_min), min(limits[2], r$observed_max))
      if (limits[1] >= limits[2]) stop("No observed range to plot: ", r$analysis, "/", r$moderator)
      extra <- c(0, r$observed_min, r$observed_max)
      x <- sort(unique(c(seq(limits[1], limits[2], length.out = 100L),
                         extra[extra >= limits[1] & extra <= limits[2]])))
      slope <- r$b_anomaly + r$b_interaction * x
      variance <- r$var_anomaly + x^2 * r$var_interaction +
        2 * x * r$cov_anomaly_interaction
      if (any(variance < -1e-10)) stop("Negative slope variance: ", r$analysis)
      se <- sqrt(pmax(0, variance))
      p_interaction <- 2 * stats::pnorm(-abs(r$b_interaction / sqrt(r$var_interaction)))
      data.frame(analysis = r$analysis, window = r$window,
        moderator = r$moderator, x = x, slope = slope, SE = se,
        lower = slope - stats::qnorm(.975) * se,
        upper = slope + stats::qnorm(.975) * se,
        observed_min = r$observed_min, observed_max = r$observed_max,
        supported = x >= r$observed_min & x <= r$observed_max,
        p_interaction_Wald = p_interaction,
        significance = if (p_interaction < .05) "p < 0.05" else "p >= 0.05",
        stringsAsFactors = FALSE)
    })
    d <- do.call(rbind, rows)
    d$variable <- factor(unname(labels[d$moderator]), levels = unname(labels))
    d$significance <- factor(d$significance, levels = c("p < 0.05", "p >= 0.05"))
    rownames(d) <- NULL
    d
  }

  resolve_font <- function() {
    requested <- getOption("pheno.spatial.figures.font", NULL)
    if (!is.null(requested)) return(requested)
    if (.Platform$OS.type == "windows") {
      if (requireNamespace("extrafont", quietly = TRUE)) {
        try(extrafont::loadfonts(device = "win", quiet = TRUE), silent = TRUE)
      }
      grDevices::windowsFonts(Garamond = grDevices::windowsFont("Garamond"))
      return("Garamond")
    }
    available <- character()
    if (requireNamespace("systemfonts", quietly = TRUE)) {
      available <- unique(systemfonts::system_fonts()$family)
    } else if (nzchar(Sys.which("fc-list"))) {
      available <- tryCatch(system2("fc-list", c(":", "family"), stdout = TRUE),
                            error = function(e) character())
    }
    available <- trimws(unlist(strsplit(available, ",", fixed = TRUE)))
    for (f in c("Garamond", "EB Garamond")) {
      if (f %in% available) return(f)
    }
    warning("Garamond was not found; using the system serif font. ",
            "Set options(pheno.spatial.figures.font = ...) to select a font.")
    "serif"
  }

  make_plot <- function(d, model, font, ylim = NULL, legend = TRUE) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package ggplot2 is required.")
    z <- d[d$analysis == model, ]
    if (!nrow(z)) stop("No curves for ", model)
    stopifnot(all(z$supported))
    key_rows <- z[match(unname(labels), as.character(z$variable)), ]
    key_types <- ifelse(key_rows$p_interaction_Wald < .05, "solid", "dashed")
    p <- ggplot2::ggplot(z, ggplot2::aes(x = x, y = slope, colour = variable,
                                       fill = variable)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper),
                           alpha = .15, colour = NA, show.legend = FALSE) +
      ggplot2::geom_line(ggplot2::aes(linetype = significance),
                         linewidth = .8, alpha = .75) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
      ggplot2::scale_colour_manual(values = colours, drop = FALSE) +
      ggplot2::scale_fill_manual(values = colours, drop = FALSE, guide = "none") +
      ggplot2::scale_linetype_manual(
        values = c("p < 0.05" = "solid", "p >= 0.05" = "dashed"),
        labels = c("p < 0.05", "p \u2265 0.05"), drop = FALSE) +
      ggplot2::theme_classic(base_family = font, base_size = 16) +
      ggplot2::theme(legend.key.width = grid::unit(1.6, "lines"),
                     legend.title = ggplot2::element_text(size = 13),
                     legend.spacing.y = grid::unit(.2, "cm")) +
      ggplot2::guides(
        colour = ggplot2::guide_legend(order = 1,
          override.aes = list(linetype = key_types)),
        linetype = ggplot2::guide_legend(order = 2,
          override.aes = list(colour = "grey20", alpha = 1))) +
      ggplot2::labs(x = "Environmental gradient",
        y = paste0("Slope of ", parameters$phenology_label[match(model, parameters$analysis)],
                   " vs. temperature anomaly"), colour = "", fill = "",
        linetype = "Anomaly \u00d7 moderator\n(Wald test)")
    if (!is.null(ylim)) p <- p + ggplot2::coord_cartesian(ylim = ylim)
    if (!legend) p <- p + ggplot2::theme(legend.position = "none")
    p
  }

  run <- function(output_dir = getOption("pheno.spatial.figures.output_dir", NULL),
                  range = c("original", "central"), font = resolve_font()) {
    range <- match.arg(range)
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package ggplot2 is required.")
    if (is.null(output_dir)) {
      root <- if (requireNamespace("here", quietly = TRUE)) here::here() else getwd()
      output_dir <- file.path(root, "output", "phenology_plasticity", "spatial_models",
                              "figures_15km")
    }
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(output_dir)) stop("Cannot create output directory: ", output_dir)
    start <- Sys.time()
    d <- make_curves(range = range)
    utils::write.csv(d, file.path(output_dir, paste0("plotted_slopes_15km_", range, ".csv")),
                     row.names = FALSE)
    utils::write.csv(parameters, file.path(output_dir, "slope_parameters_15km.csv"),
                     row.names = FALSE)
    tests <- parameters[c("analysis", "window", "moderator", "b_interaction", "var_interaction")]
    tests$SE_interaction <- sqrt(tests$var_interaction)
    tests$p_Wald <- 2 * stats::pnorm(-abs(tests$b_interaction / tests$SE_interaction))
    tests$line_type <- ifelse(tests$p_Wald < .05, "solid", "dashed")
    tests$test <- "Two-sided normal Wald; no multiplicity correction"
    utils::write.csv(tests, file.path(output_dir, "interaction_tests_15km.csv"), row.names = FALSE)
    plots <- stats::setNames(vector("list", length(model_order)), model_order)
    files <- character()
    for (i in seq_along(model_order)) {
      m <- model_order[i]
      plots[[m]] <- make_plot(d, m, font)
      base <- file.path(output_dir, paste0(m, "_spatial_15km_plasticity_", range))
      ggplot2::ggsave(paste0(base, ".png"), plots[[m]], width = 7, height = 5,
                      units = "in", dpi = 300, bg = "white")
      files <- c(files, paste0(base, ".png"))
      if (capabilities("cairo")) {
        ggplot2::ggsave(paste0(base, ".pdf"), plots[[m]], width = 7, height = 5,
                        units = "in", device = grDevices::cairo_pdf, bg = "white")
        files <- c(files, paste0(base, ".pdf"))
      } else {
        warning("Cairo PDF is unavailable; PNG was saved for ", m)
      }
      message(sprintf("Figure %d/4 saved: %s (%.0f seconds elapsed)",
                      i, display_names[[m]], as.numeric(difftime(Sys.time(), start, units = "secs"))))
    }
    notes <- c(
      "Spatial plasticity curves: 15-km mesh, selected climatic windows only.",
      "Onset: 60 d. First peak and both offsets: 90 d.",
      "Multivoltine estimates are from the recovered numerical fit.",
      "Conditional fixed-effect slopes; all other moderators are held at zero.",
      "Bands: pointwise normal Wald 95% confidence intervals, including coefficient covariance.",
      "Solid coloured lines: anomaly x moderator interaction p < 0.05; dashed: p >= 0.05.",
      "Tests are two-sided normal Wald tests without multiplicity correction, not LRTs.",
      "Line style refers to the interaction coefficient, not the conditional slope at each x.",
      "Curves and ribbons are restricted to the observed marginal predictor range; no extrapolation.",
      "The horizontal grey dashed line marks a zero slope.",
      "Predictor scales are inherited unchanged; no conversion to days per degree C is made.",
      "No refitting, window selection, likelihood-ratio tests or population-trend modelling.",
      paste("Font:", font), paste("Gradient range mode:", range),
      "Original design: pheno_plasticity(1).R, make_plasticity_plot().")
    writeLines(notes, file.path(output_dir, "figure_notes.txt"))
    writeLines(capture.output(utils::sessionInfo()), file.path(output_dir, "figure_sessionInfo.txt"))
    message("Done. Figures: ", normalizePath(output_dir, winslash = "/"))
    invisible(list(plots = plots, curves = d, parameters = parameters,
                   files = files, output_dir = output_dir))
  }

  list(run = run, curves = make_curves, plot = make_plot, parameters = parameters,
       model_order = model_order, display_names = display_names)
})

if (!isTRUE(getOption("pheno.spatial.figures.functions_only", FALSE))) {
  pheno_spatial_figures_15km$run()
}
