# ============================================================================ #
# phenoIMPACT | Spatial population-trend models, 15-km mesh
# Two panels: relative effects of onset / offset plasticity on decadal change.
#
# Source in RStudio: source(file.choose(), encoding = "UTF-8")
#
# Main figure reads <group>__full/fixed_effects.csv ONLY.
# Optional joint contrasts also read the two SAVED fit.rds lists for beta and V.
# NO fitting, optimization, DLL loading, simulation, or package/font installation.
# No input file or previous figure is changed. All new outputs have their own folder.
#
# Main estimand, at each z in [-1, +1]:
#   delta_log_slope = z * beta(year_decade:plasticity)
#   percent = 100 * (exp(delta_log_slope) - 1)
# This is the percentage difference in the DECADAL MULTIPLICATIVE CHANGE FACTOR,
# relative to the same group's mean plasticity, other plasticity held constant.
# It is NOT an abundance trend, a percentage-point difference between trends,
# a variance-explained measure, or a causal effect.
#
# High onset = stronger advancement in both groups.
# High offset = stronger advancement (univoltine), more delay (multivoltine).
# -1 SD means below the group mean, NOT necessarily a delay or no plasticity.
#
# Bands: transformed pointwise 95% Wald CIs for the contrast, not prediction
# intervals or simultaneous bands. At z=0, estimate and CI are zero BY DEFINITION.
# For z<0 the lower / upper limits are reversed before exponentiation.
#
# Optional JOINT contrasts: changes from (0,0) of exactly 1 SD in each component.
# These are EXPLORATORY combinations in the existing additive trend model.
# 1 SD onset and 1 SD offset do NOT represent an equal shift in days or days/degree.
# Standardization remains that of the original fitted models, separately by group.
#
# Graphics: https://ggplot2.tidyverse.org/reference/geom_ribbon.html
#           https://ggplot2.tidyverse.org/reference/facet_wrap.html
# Contrasts: https://rvlenth.github.io/emmeans/articles/comparisons.html
# ============================================================================ #

plot_population_trends_relative_effects <- function(
    run_dir = paste0(
      "E:/phenoIMPACT project/code/phenoIMPACT/output/population_trends/",
      "spatial_abundance_15km/run_650c41f15949"),
    out_dir = file.path(run_dir, "figures_relative_plasticity_effects"),
    font_family = "Garamond",
    base_size = 15,
    width = 9.5,
    height = 4.3,
    dpi = 300,
    group_colours = c(univoltine = "#0072B2", multivoltine = "#D55E00"),
    compute_joint = TRUE,
    export_svg = TRUE,
    show_plot = interactive()) {

  fail <- function(...) stop(..., call. = FALSE)
  require_true <- function(ok, message) if (!isTRUE(ok)) fail(message)
  scalar_positive <- function(x) is.numeric(x) && length(x) == 1L &&
    is.finite(x) && x > 0
  groups <- c("univoltine", "multivoltine")
  variables <- c(onset = "onset_plasticity_z", offset = "offset_plasticity_bio_z")

  require_true(requireNamespace("ggplot2", quietly = TRUE),
               "Package ggplot2 is required; no model has been touched.")
  require_true(utils::packageVersion("ggplot2") >= "3.4.0",
               "This script requires ggplot2 >= 3.4.0.")
  require_true(is.character(run_dir) && length(run_dir) == 1L &&
                 !is.na(run_dir) && dir.exists(run_dir), "Run directory not found.")
  require_true(is.character(out_dir) && length(out_dir) == 1L &&
                 !is.na(out_dir) && nzchar(out_dir), "Invalid out_dir.")
  require_true(is.character(font_family) && length(font_family) == 1L &&
                 !is.na(font_family) && nzchar(font_family), "Invalid font_family.")
  for (value in list(base_size, width, height, dpi)) {
    require_true(scalar_positive(value), "Figure sizes and dpi must be positive scalars.")
  }
  for (value in list(compute_joint, export_svg, show_plot)) {
    require_true(is.logical(value) && length(value) == 1L && !is.na(value),
                 "compute_joint, export_svg and show_plot must be TRUE/FALSE.")
  }
  require_true(is.character(group_colours) && length(group_colours) == 2L &&
                 !anyNA(group_colours) && !anyDuplicated(names(group_colours)) &&
                 setequal(names(group_colours), groups),
               "Supply two named group_colours: univoltine and multivoltine.")
  tryCatch(grDevices::col2rgb(group_colours), error = function(e) {
    fail("Invalid group colours: ", conditionMessage(e))
  })
  if (requireNamespace("systemfonts", quietly = TRUE)) {
    families <- unique(systemfonts::system_fonts()$family)
    require_true(tolower(font_family) %in% tolower(families),
      paste0("Font '", font_family, "' is not installed. Choose an installed font explicitly."))
  }
  if (.Platform$OS.type == "windows") {
    do.call(grDevices::windowsFonts,
      stats::setNames(list(grDevices::windowsFont(font_family)), font_family))
  }

  csv_paths <- stats::setNames(
    file.path(run_dir, paste0(groups, "__full"), "fixed_effects.csv"), groups)
  require_true(all(file.exists(csv_paths)), "A spatial fixed_effects.csv is missing.")
  input_manifest <- data.frame(group = groups, input_type = "fixed_effects",
    source_file = unname(csv_paths), md5 = unname(tools::md5sum(csv_paths)))
  all_fixed <- stats::setNames(vector("list", 2L), groups)
  interaction_rows <- list()
  find_term <- function(term_names, variable) {
    candidate <- c(paste0("year_decade:", variable), paste0(variable, ":year_decade"))
    found <- candidate[candidate %in% term_names]
    require_true(length(found) == 1L, paste("Missing/ambiguous interaction:", variable))
    found[[1L]]
  }
  for (g in groups) {
    a <- utils::read.csv(csv_paths[[g]], stringsAsFactors = FALSE, check.names = FALSE)
    needed <- c("term", "estimate", "std.error", "conf.low", "conf.high", "p_Wald")
    require_true(nrow(a) > 0L && all(needed %in% names(a)) &&
                   !anyNA(a$term) && !anyDuplicated(a$term),
                 paste("Unexpected fixed-effect table for", g))
    require_true(all(vapply(a[setdiff(needed, "term")], function(v) {
      is.numeric(v) && all(is.finite(v))
    }, logical(1))), paste("Invalid numeric coefficients for", g))
    require_true(all(a$std.error > 0 & a$conf.low <= a$estimate &
                       a$estimate <= a$conf.high & a$p_Wald >= 0 & a$p_Wald <= 1),
                 paste("Invalid uncertainty estimates for", g))
    all_fixed[[g]] <- a
    for (component in names(variables)) {
      term <- find_term(a$term, variables[[component]])
      r <- a[a$term == term, needed, drop = FALSE]
      zci <- c(r$estimate - r$conf.low, r$conf.high - r$estimate) / r$std.error
      require_true(all(abs(zci - stats::qnorm(0.975)) < 1e-4),
                   "Expected pointwise 95% normal Wald intervals in the CSV.")
      interaction_rows[[paste(g, component)]] <- data.frame(
        group = g, component = component, term = term,
        estimate = r$estimate, SE = r$std.error,
        CI_low = r$conf.low, CI_high = r$conf.high, p_Wald = r$p_Wald)
    }
  }
  interactions <- do.call(rbind, interaction_rows)
  rownames(interactions) <- NULL

  # Keep EXACTLY the original group-specific biological orientation.
  z_grid <- sort(unique(c(seq(-1, 1, length.out = 201L), 0)))
  curves <- do.call(rbind, lapply(seq_len(nrow(interactions)), function(i) {
    a <- interactions[i, ]
    log_est <- z_grid * a$estimate
    log_lo <- pmin(z_grid * a$CI_low, z_grid * a$CI_high)
    log_hi <- pmax(z_grid * a$CI_low, z_grid * a$CI_high)
    data.frame(group = a$group, component = a$component, plasticity_z = z_grid,
      delta_log_per_decade = log_est, delta_log_CI_low = log_lo, delta_log_CI_high = log_hi,
      effect_percent = 100 * expm1(log_est),
      low_percent = 100 * expm1(log_lo), high_percent = 100 * expm1(log_hi))
  }))
  rownames(curves) <- NULL
  require_true(all(is.finite(as.matrix(curves[c("effect_percent", "low_percent", "high_percent")]))),
               "Non-finite transformed contrasts.")
  require_true(all(curves$low_percent <= curves$effect_percent &
                     curves$effect_percent <= curves$high_percent), "CI order is incorrect.")
  require_true(all(as.matrix(curves[curves$plasticity_z == 0,
    c("effect_percent", "low_percent", "high_percent")]) == 0), "The reference contrast must be zero.")
  curves$panel <- factor(curves$component, levels = c("onset", "offset"),
                        labels = c("a   Onset", "b   Offset"))
  curves$group <- factor(curves$group, levels = groups)
  endpoints <- curves[curves$plasticity_z == 1, , drop = FALSE]
  endpoints$label <- ifelse(endpoints$group == "univoltine", "Univoltine", "Multivoltine")
  endpoints$label[endpoints$component == "offset" & endpoints$group == "univoltine"] <-
    "Univoltine\n(advancement)"
  endpoints$label[endpoints$component == "offset" & endpoints$group == "multivoltine"] <-
    "Multivoltine\n(delay)"
  lim <- 1.10 * max(abs(c(curves$low_percent, curves$high_percent)))
  percent_ticks <- function(x) {
    ifelse(abs(x) < 1e-10, "0", paste0(ifelse(x > 0, "+", ""),
      format(x, trim = TRUE, scientific = FALSE), "%"))
  }

  p <- ggplot2::ggplot(curves, ggplot2::aes(x = plasticity_z, y = effect_percent,
      colour = group, fill = group, group = group)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey50", linetype = "dashed", linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.3) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = low_percent, ymax = high_percent),
                        alpha = 0.13, colour = NA) +
    ggplot2::geom_line(ggplot2::aes(linetype = group), linewidth = 1.05, lineend = "round") +
    ggplot2::geom_segment(data = endpoints, inherit.aes = FALSE,
      ggplot2::aes(x = 1, xend = 1.07, y = effect_percent, yend = effect_percent, colour = group),
      linewidth = 0.5, show.legend = FALSE) +
    ggplot2::geom_text(data = endpoints, inherit.aes = FALSE,
      ggplot2::aes(x = 1.12, y = effect_percent, label = label, colour = group),
      family = font_family, size = base_size * 0.245, hjust = 0, vjust = 0.5,
      lineheight = 0.95, show.legend = FALSE) +
    ggplot2::facet_wrap(~panel, nrow = 1, scales = "fixed") +
    ggplot2::scale_colour_manual(values = group_colours, drop = FALSE) +
    ggplot2::scale_fill_manual(values = group_colours, drop = FALSE) +
    ggplot2::scale_linetype_manual(values = c(univoltine = "solid", multivoltine = "22")) +
    ggplot2::scale_x_continuous(breaks = c(-1, 0, 1), labels = c("-1", "Mean", "+1"),
                               expand = ggplot2::expansion(mult = 0)) +
    ggplot2::scale_y_continuous(labels = percent_ticks,
                               expand = ggplot2::expansion(mult = 0)) +
    ggplot2::coord_cartesian(xlim = c(-1.08, 1.94), ylim = c(-lim, lim), clip = "off") +
    ggplot2::labs(x = "Plasticity (SD from the group mean)",
      y = "Decadal change factor\n(% difference from mean plasticity)") +
    ggplot2::theme_classic(base_family = font_family, base_size = base_size) +
    ggplot2::theme(legend.position = "none",
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(hjust = 0, face = "bold", size = base_size),
      axis.text = ggplot2::element_text(colour = "grey25", size = base_size * 0.85),
      axis.title = ggplot2::element_text(size = base_size * 0.90),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 9)),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 9)),
      axis.line = ggplot2::element_line(linewidth = 0.45),
      axis.ticks = ggplot2::element_line(linewidth = 0.4),
      panel.spacing.x = grid::unit(1.0, "lines"),
      plot.margin = ggplot2::margin(8, 14, 8, 8))

  # Read a saved parameter list only for a 6 x 6 covariance, never call a model.
  joint <- covariances <- checks <- list()
  if (compute_joint) {
    for (g in groups) {
      path <- file.path(run_dir, paste0(g, "__full"), "fit.rds")
      require_true(file.exists(path), paste("Missing fit.rds for joint contrast:", g))
      message("Reading saved beta and V for joint contrasts: ", g, " (no refit).")
      saved <- readRDS(path)
      require_true(is.list(saved) && all(c("beta", "V", "group", "variant") %in% names(saved)),
                   "Expected the saved spatial-fit parameter list.")
      require_true(identical(as.character(saved$group), g) &&
                     identical(as.character(saved$variant), "full"), "Incorrect fit group/variant.")
      if (!is.null(saved$checks$numerical_checks_ok)) {
        require_true(isTRUE(saved$checks$numerical_checks_ok), "The saved fit requires numerical review.")
      }
      b <- saved$beta; V <- as.matrix(saved$V)
      require_true(is.numeric(b) && !is.null(names(b)) && !anyDuplicated(names(b)) &&
        all(is.finite(b)) && !is.null(rownames(V)) && !is.null(colnames(V)) &&
        !anyDuplicated(rownames(V)) && !anyDuplicated(colnames(V)) &&
        setequal(rownames(V), names(b)) && setequal(colnames(V), names(b)),
        "Coefficients or covariance names are incomplete/ambiguous.")
      V <- V[names(b), names(b), drop = FALSE]
      require_true(all(is.finite(V)) && max(abs(V - t(V))) < 1e-8 &&
        min(eigen(V, symmetric = TRUE, only.values = TRUE)$values) > 0,
        "Fixed-effect covariance is not finite symmetric positive definite.")
      a <- all_fixed[[g]]
      require_true(setequal(a$term, names(b)), "Fit and CSV coefficient names differ.")
      a <- a[match(names(b), a$term), , drop = FALSE]
      beta_error <- max(abs(b - a$estimate))
      SE_error <- max(abs(sqrt(diag(V)) - a$std.error))
      require_true(beta_error < 1e-8 && SE_error < 1e-8,
                   "Fit and CSV estimates/SEs disagree; no joint contrast is produced.")
      onset_term <- find_term(names(b), variables[["onset"]])
      offset_term <- find_term(names(b), variables[["offset"]])
      offset_direction <- if (g == "univoltine") 1 else -1
      scenarios <- data.frame(
        scenario = c("advance_both_equal_SD", "advance_onset_shift_offset_towards_delay"),
        delta_onset_z = c(1, 1), delta_offset_bio_z = c(offset_direction, -offset_direction))
      for (j in seq_len(nrow(scenarios))) {
        s <- scenarios[j, ]
        L <- stats::setNames(numeric(length(b)), names(b))
        L[onset_term] <- s$delta_onset_z
        L[offset_term] <- s$delta_offset_bio_z
        est <- as.numeric(sum(L * b))
        variance <- as.numeric(t(L) %*% V %*% L)
        require_true(is.finite(variance) && variance > 0, "Invalid joint contrast variance.")
        SE <- sqrt(variance); lo <- est - stats::qnorm(0.975) * SE
        hi <- est + stats::qnorm(0.975) * SE
        joint[[paste(g, j)]] <- data.frame(group = g, s,
          contrast_log_per_decade = est, SE = SE, CI_low = lo, CI_high = hi,
          p_Wald = 2 * stats::pnorm(-abs(est / SE)),
          relative_decadal_factor_percent = 100 * expm1(est),
          percent_low = 100 * expm1(lo), percent_high = 100 * expm1(hi),
          source_fit = normalizePath(path, winslash = "/", mustWork = TRUE))
      }
      covariances[[g]] <- V
      checks[[g]] <- data.frame(group = g, max_beta_error = beta_error, max_SE_error = SE_error,
        cov_onset_offset = V[onset_term, offset_term],
        cor_onset_offset = V[onset_term, offset_term] /
          sqrt(V[onset_term, onset_term] * V[offset_term, offset_term]))
      input_manifest <- rbind(input_manifest, data.frame(group = g, input_type = "saved_fit",
        source_file = path, md5 = unname(tools::md5sum(path))))
      rm(saved); invisible(gc())
    }
    joint <- do.call(rbind, joint); rownames(joint) <- NULL
    # Exploratory tests: also report a Holm adjustment across these FOUR contrasts.
    # Confidence limits remain pointwise 95% intervals, not adjusted intervals.
    joint$p_Holm_four_exploratory <- stats::p.adjust(joint$p_Wald, method = "holm")
    checks <- do.call(rbind, checks); rownames(checks) <- NULL
  } else {
    joint <- checks <- NULL
  }

  # All outputs are derived; original model files remain read-only.
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  require_true(dir.exists(out_dir), "Cannot create output folder.")
  png_path <- file.path(out_dir, "plasticity_relative_effects_spatial_15km.png")
  use_ragg <- requireNamespace("ragg", quietly = TRUE)
  device <- if (use_ragg) ragg::agg_png else grDevices::png
  extra <- if (use_ragg) list() else if (.Platform$OS.type == "windows") {
    list(type = "windows")
  } else if (capabilities("cairo")) list(type = "cairo") else list()
  do.call(ggplot2::ggsave, c(list(filename = png_path, plot = p, device = device,
    width = width, height = height, units = "in", dpi = dpi, bg = "white"), extra))
  svg_path <- NA_character_
  if (export_svg) {
    svg_device <- if (requireNamespace("svglite", quietly = TRUE)) {
      svglite::svglite
    } else if (capabilities("cairo")) grDevices::svg else NULL
    if (!is.null(svg_device)) {
      svg_path <- file.path(out_dir, "plasticity_relative_effects_spatial_15km.svg")
      ggplot2::ggsave(filename = svg_path, plot = p, device = svg_device,
                     width = width, height = height, units = "in", bg = "white")
    } else warning("No SVG device available; PNG and editable ggplot are still saved.")
  }
  utils::write.csv(curves, file.path(out_dir, "relative_effect_curves.csv"), row.names = FALSE)
  utils::write.csv(interactions, file.path(out_dir, "source_interactions.csv"), row.names = FALSE)
  utils::write.csv(input_manifest, file.path(out_dir, "input_manifest.csv"), row.names = FALSE)
  if (compute_joint) {
    utils::write.csv(joint, file.path(out_dir, "joint_plasticity_contrasts.csv"), row.names = FALSE)
    utils::write.csv(checks, file.path(out_dir, "joint_reconstruction_checks.csv"), row.names = FALSE)
    for (g in groups) utils::write.csv(data.frame(term = rownames(covariances[[g]]),
      covariances[[g]], check.names = FALSE),
      file.path(out_dir, paste0(g, "__fixed_effect_covariance.csv")), row.names = FALSE)
  }
  saveRDS(p, file.path(out_dir, "plasticity_relative_effects_spatial_15km.rds"))
  writeLines(capture.output(utils::sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
  writeLines(c(
    "RELATIVE PLASTICITY EFFECTS | TWO-PANEL SPATIAL TREND FIGURE",
    "Main figure: 100 * [exp(z * beta_year:plasticity) - 1].",
    "Reference: z=0 of the SAME group. Other plasticity unchanged; fixed-effect contrast.",
    "Y is a percentage difference in decadal multiplicative change factor, NOT percentage points.",
    "Positive values do not necessarily imply growing populations; zero is a comparison reference.",
    "Band limits use the exported 95% Wald coefficient limits, reordered when z < 0.",
    "At z=0, the estimate and its interval vanish by construction (a self-comparison).",
    "No simulation or new information is supplied by plotting a denser z grid.",
    "Both groups are separately standardized: 1 SD is NOT a common day/degree unit.",
    "Onset +z means stronger advancement in both groups; -z means below average, not necessarily delay.",
    "Offset +z: stronger advancement for univoltines, more delay for multivoltines.",
    "No between-group hypothesis test, variance-explained measure, or causal effect is estimated.",
    "JOINT EXPLORATORY CONTRASTS, when requested:",
    "advance_both_equal_SD: (delta_onset_z,delta_offset_bio_z) = (1,1) uni; (1,-1) multi.",
    "advance_onset_shift_offset_towards_delay: (1,-1) uni; (1,1) multi.",
    "Reference for all joint contrasts: group-specific mean of both standardized predictors.",
    "Delta = L' beta; Var(Delta) = L' V L, including onset-offset coefficient covariance.",
    "Pointwise 95% Wald intervals, unadjusted p-values and Holm p-values for the 4 exploratory tests.",
    "These are additive-model contrasts, not tests of an onset-by-offset interaction.",
    "Equal SD shifts do NOT demonstrate equal date shifts or unchanged flight-season length.",
    "They do not establish that extending flight duration causes better population trends.",
    "A rigid seasonal shift requires weights based on original, comparable plasticity units.",
    "The joint combinations have not been screened for bivariate observational support here.",
    "Both coefficients and their uncertainty condition on the existing estimated plasticity inputs.",
    "NO MODEL WAS REFITTED. No model parameters, input files or old figures were changed.",
    "The PNG/SVG have no overall title or caption; panel labels identify onset and offset."
  ), file.path(out_dir, "README.txt"))
  require_true(identical(unname(tools::md5sum(input_manifest$source_file)), input_manifest$md5),
               "An input file changed during this export; review and rerun.")
  if (show_plot) print(p)
  message("DONE. Relative-effect figure exported; no models were refitted.")
  message("Output: ", normalizePath(out_dir, winslash = "/", mustWork = TRUE))
  if (compute_joint) {
    message("EXPLORATORY JOINT CONTRASTS, 1 SD each (not equal changes in days):")
    print(joint[c("group", "scenario", "contrast_log_per_decade", "CI_low", "CI_high",
                  "p_Wald", "p_Holm_four_exploratory")], row.names = FALSE)
  }
  invisible(list(plot = p, curves = curves, interactions = interactions, joint = joint,
    files = c(png = png_path, svg = svg_path), output_dir = out_dir))
}

# Define without executing:
# options(phenoimpact.relative_effects_functions_only = TRUE)
# source(file.choose(), encoding = "UTF-8")
# population_trends_relative_effects <- plot_population_trends_relative_effects()
#
# Plot only, without reading the saved fit lists:
# population_trends_relative_effects <- plot_population_trends_relative_effects(compute_joint = FALSE)
#
# Editable plot: population_trends_relative_effects$plot
if (!isTRUE(getOption("phenoimpact.relative_effects_functions_only", FALSE))) {
  population_trends_relative_effects <- plot_population_trends_relative_effects()
}
