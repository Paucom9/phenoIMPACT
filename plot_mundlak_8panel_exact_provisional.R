
# ============================================================================
# phenoIMPACT | Mundlak figure (8 panels: onset + offset only)
# Provisional outputs, exact fixed-effect predictions + exact 95% CIs
#
# Layout:
#   Columns 1-2: Univoltine (Within / Between species)
#   Columns 3-4: Multivoltine (Within / Between species)
#   Row 1: Onset
#   Row 2: Offset
#
# Significance codes are placed bottom-left inside each panel.
#
# No model is refitted.
# Predictions use:
#   slope = d' beta
#   Var(slope) = d' V d
# and then:
#   Relative change (%) = 100 * [exp(slope * decades_since_1990) - 1]
#
# High / Low = +1 SD / -1 SD of the relevant Mundlak component.
# ============================================================================

plot_mundlak_8panel <- function(
  zip_file = file.choose(),
  out_dir = dirname(zip_file),
  prefix = "mundlak_8panel_exact_provisional",
  font_family = "Arial",
  width_mm = 180,
  height_mm = 108,
  dpi = 600
) {

  pkgs <- c("ggplot2", "patchwork", "dplyr", "grid")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      "Missing packages: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  read_zip_csv <- function(path) {
    con <- unz(zip_file, path)
    on.exit(close(con), add = TRUE)
    utils::read.csv(
      con,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }

  groups <- c("univoltine", "multivoltine")

  fixed <- stats::setNames(
    lapply(groups, function(g) {
      x <- read_zip_csv(
        file.path(paste0(g, "__mundlak"), "fixed_effects.csv")
      )
      rownames(x) <- x$term
      x
    }),
    groups
  )

  vcovs <- stats::setNames(
    lapply(groups, function(g) {
      x <- read_zip_csv(
        file.path(
          paste0(g, "__mundlak"),
          "fixed_effect_covariance.csv"
        )
      )
      rownames(x) <- x$term
      out <- as.matrix(x[, -1, drop = FALSE])
      storage.mode(out) <- "double"
      out
    }),
    groups
  )

  scales <- stats::setNames(
    lapply(groups, function(g) {
      read_zip_csv(paste0("component_scales_", g, ".csv"))
    }),
    groups
  )

  component_sd <- function(group, component) {
    x <- scales[[group]]
    value <- x$SD_for_descriptive_rescaling[
      x$component == component
    ]
    if (length(value) != 1L || !is.finite(value)) {
      stop(
        "Could not resolve SD for ",
        group, " / ", component,
        call. = FALSE
      )
    }
    value
  }

  sig_symbol <- function(p) {
    if (p < 0.001) return("***")
    if (p < 0.01)  return("**")
    if (p < 0.05)  return("*")
    "n.s."
  }

  panel_spec <- data.frame(
    panel = letters[1:8],
    row = rep(c("Onset", "Offset"), each = 4),
    group = rep(
      c(
        "univoltine", "univoltine",
        "multivoltine", "multivoltine"
      ),
      2
    ),
    scale = rep(c("within", "between", "within", "between"), 2),
    stringsAsFactors = FALSE
  )

  get_component <- function(group, row, scale) {

    if (row == "Onset") {
      component <- paste0("onset_", scale, "_species")
      if (scale == "within") {
        variable <- "onset_plasticity_within_species_z"
      } else {
        variable <- "onset_plasticity_species_mean_z"
      }
    } else {
      component <- paste0("offset_", scale, "_species")
      if (scale == "within") {
        variable <- "offset_plasticity_within_species_z"
      } else {
        variable <- "offset_plasticity_species_mean_z"
      }
    }

    list(
      variable = variable,
      component = component,
      SD = component_sd(group, component)
    )
  }

  years <- 1990:2023

  build_curve <- function(
    group,
    row,
    scale,
    plasticity_z
  ) {

    comp <- get_component(group, row, scale)

    fe <- fixed[[group]]
    V <- vcovs[[group]]

    terms <- rownames(V)
    beta <- stats::setNames(fe$estimate, fe$term)

    int_term <- paste0("year_decade:", comp$variable)

    needed <- c("year_decade", int_term)
    if (!all(needed %in% terms)) {
      stop(
        "Missing term(s) for ",
        group, " / ", row, " / ", scale,
        call. = FALSE
      )
    }

    d <- stats::setNames(rep(0, length(terms)), terms)
    d["year_decade"] <- 1
    d[int_term] <- plasticity_z * comp$SD

    slope <- sum(d * beta[names(d)])
    slope_se <- sqrt(
      drop(
        t(d) %*%
          V[names(d), names(d), drop = FALSE] %*%
          d
      )
    )

    decades <- (years - min(years)) / 10

    data.frame(
      YEAR = years,
      estimate = 100 * expm1(slope * decades),
      low = 100 * expm1(
        (slope - 1.96 * slope_se) * decades
      ),
      high = 100 * expm1(
        (slope + 1.96 * slope_se) * decades
      ),
      slope_log_per_decade = slope,
      slope_SE = slope_se,
      component_SD = comp$SD,
      variable = comp$variable,
      stringsAsFactors = FALSE
    )
  }

  predictions <- do.call(
    rbind,
    lapply(seq_len(nrow(panel_spec)), function(i) {

      spec <- panel_spec[i, ]

      do.call(
        rbind,
        lapply(
          c(1, 0, -1),
          function(z) {

            d <- build_curve(
              spec$group,
              spec$row,
              spec$scale,
              z
            )

            d$panel <- spec$panel
            d$row <- spec$row
            d$group <- spec$group
            d$scale <- spec$scale
            d$plasticity_z <- z
            d$Plasticity <- factor(
              z,
              levels = c(1, 0, -1),
              labels = c(
                "High (+1 SD)",
                "Mean",
                "Low (-1 SD)"
              )
            )
            d
          }
        )
      )
    })
  )

  # Add exact p-values for the focal interaction.
  verification <- do.call(
    rbind,
    lapply(seq_len(nrow(panel_spec)), function(i) {

      spec <- panel_spec[i, ]
      comp <- get_component(
        spec$group,
        spec$row,
        spec$scale
      )

      term <- paste0(
        "year_decade:",
        comp$variable
      )

      fe <- fixed[[spec$group]]

      data.frame(
        panel = spec$panel,
        row = spec$row,
        group = spec$group,
        scale = spec$scale,
        term = term,
        estimate = fe[term, "estimate"],
        std.error = fe[term, "std.error"],
        conf.low = fe[term, "conf.low"],
        conf.high = fe[term, "conf.high"],
        p_Wald = fe[term, "p_Wald"],
        significance = sig_symbol(fe[term, "p_Wald"]),
        component_SD = comp$SD,
        stringsAsFactors = FALSE
      )
    })
  )

  predictions <- dplyr::left_join(
    predictions,
    verification[
      ,
      c(
        "panel",
        "p_Wald",
        "significance"
      )
    ],
    by = "panel"
  )

  palette <- c(
    "High (+1 SD)" = "#009E73",
    "Mean" = "#0072B2",
    "Low (-1 SD)" = "#D55E00"
  )

  nature_theme <- function() {
    ggplot2::theme_classic(
      base_family = font_family,
      base_size = 9
    ) +
      ggplot2::theme(
        panel.background = ggplot2::element_rect(
          fill = "white",
          colour = NA
        ),
        plot.background = ggplot2::element_rect(
          fill = "white",
          colour = NA
        ),
        axis.line = ggplot2::element_line(
          linewidth = 0.4,
          colour = "black"
        ),
        axis.ticks = ggplot2::element_line(
          linewidth = 0.4,
          colour = "black"
        ),
        plot.margin = ggplot2::margin(
          2, 2, 2, 2,
          unit = "pt"
        ),
        legend.position = "bottom",
        legend.title = ggplot2::element_text(
          face = "bold"
        ),
        legend.key.width = grid::unit(
          12,
          "pt"
        ),
        legend.spacing.x = grid::unit(
          4,
          "pt"
        )
      )
  }

  make_panel <- function(
    panel_letter,
    show_y = TRUE,
    show_x = TRUE
  ) {

    d <- predictions[
      predictions$panel == panel_letter,
      ,
      drop = FALSE
    ]

    meta <- verification[
      verification$panel == panel_letter,
      ,
      drop = FALSE
    ]

    group <- meta$group

    if (group == "univoltine") {
      ylim <- c(-50, 5)
      ybreaks <- c(-50, -40, -30, -20, -10, 0)
    } else {
      ylim <- c(-30, 30)
      ybreaks <- c(-30, -20, -10, 0, 10, 20, 30)
    }

    p <- ggplot2::ggplot() +
      ggplot2::geom_hline(
        yintercept = 0,
        linetype = "dashed",
        colour = "grey45",
        linewidth = 0.4
      ) +

      # Exact CIs from d' V d.
      ggplot2::geom_ribbon(
        data = d,
        ggplot2::aes(
          x = YEAR,
          ymin = low,
          ymax = high,
          fill = Plasticity,
          group = Plasticity
        ),
        alpha = 0.10,
        colour = NA
      ) +

      ggplot2::geom_line(
        data = d,
        ggplot2::aes(
          x = YEAR,
          y = estimate,
          colour = Plasticity,
          group = Plasticity
        ),
        linewidth = 0.7
      ) +

      ggplot2::scale_colour_manual(
        name = "Plasticity",
        values = palette,
        breaks = names(palette),
        drop = FALSE
      ) +

      ggplot2::scale_fill_manual(
        name = "Plasticity",
        values = palette,
        breaks = names(palette),
        drop = FALSE
      ) +

      ggplot2::scale_x_continuous(
        breaks = c(1990, 2000, 2010, 2020),
        expand = ggplot2::expansion(
          mult = c(0.02, 0.03)
        )
      ) +

      ggplot2::scale_y_continuous(
        breaks = ybreaks
      ) +

      ggplot2::coord_cartesian(
        ylim = ylim,
        clip = "on"
      ) +

      ggplot2::labs(
        tag = panel_letter,
        x = NULL,
        y = NULL
      ) +

      # Significance bottom-left.
      ggplot2::annotate(
        "text",
        x = 1990,
        y = ylim[1] + 0.04 * diff(ylim),
        label = meta$significance,
        hjust = 0,
        vjust = 0,
        family = font_family,
        fontface = ifelse(
          meta$significance == "n.s.",
          "plain",
          "bold"
        ),
        size = 3.2
      ) +

      nature_theme() +

      ggplot2::theme(
        plot.tag = ggplot2::element_text(
          face = "bold",
          size = 9
        ),
        plot.tag.position = c(
          -0.08,
          1.02
        ),
        axis.text.y = if (show_y) {
          ggplot2::element_text()
        } else {
          ggplot2::element_blank()
        },
        axis.ticks.y = if (show_y) {
          ggplot2::element_line(
            linewidth = 0.4
          )
        } else {
          ggplot2::element_blank()
        },
        axis.text.x = if (show_x) {
          ggplot2::element_text()
        } else {
          ggplot2::element_blank()
        },
        axis.ticks.x = if (show_x) {
          ggplot2::element_line(
            linewidth = 0.4
          )
        } else {
          ggplot2::element_blank()
        }
      )

    p
  }

  p_a <- make_panel("a", TRUE, TRUE)
  p_b <- make_panel("b", TRUE, TRUE)
  p_c <- make_panel("c", TRUE, TRUE)
  p_d <- make_panel("d", TRUE, TRUE)

  p_e <- make_panel("e", TRUE, TRUE)
  p_f <- make_panel("f", TRUE, TRUE)
  p_g <- make_panel("g", TRUE, TRUE)
  p_h <- make_panel("h", TRUE, TRUE)

  body <- (
    p_a | p_b | p_c | p_d
  ) / (
    p_e | p_f | p_g | p_h
  ) +
    patchwork::plot_layout(
      guides = "collect",
      heights = c(1, 1)
    )

  # Header matching the agreed mock-up structure.
  header_df <- data.frame(
    x = c(1, 2, 3, 4),
    label = c(
      "Within species",
      "Between species",
      "Within species",
      "Between species"
    )
  )

  header <- ggplot2::ggplot() +
    ggplot2::annotate(
      "text",
      x = 1.5,
      y = 0.90,
      label = "Univoltine",
      family = font_family,
      fontface = "bold",
      size = 4.0
    ) +
    ggplot2::annotate(
      "text",
      x = 3.5,
      y = 0.90,
      label = "Multivoltine",
      family = font_family,
      fontface = "bold",
      size = 4.0
    ) +
    ggplot2::annotate(
      "segment",
      x = 0.65, xend = 2.35,
      y = 0.73, yend = 0.73,
      linewidth = 0.4
    ) +
    ggplot2::annotate(
      "segment",
      x = 2.65, xend = 4.35,
      y = 0.73, yend = 0.73,
      linewidth = 0.4
    ) +
    ggplot2::annotate(
      "segment",
      x = 0.65, xend = 0.65,
      y = 0.73, yend = 0.60,
      linewidth = 0.4
    ) +
    ggplot2::annotate(
      "segment",
      x = 2.35, xend = 2.35,
      y = 0.73, yend = 0.60,
      linewidth = 0.4
    ) +
    ggplot2::annotate(
      "segment",
      x = 2.65, xend = 2.65,
      y = 0.73, yend = 0.60,
      linewidth = 0.4
    ) +
    ggplot2::annotate(
      "segment",
      x = 4.35, xend = 4.35,
      y = 0.73, yend = 0.60,
      linewidth = 0.4
    ) +
    ggplot2::geom_text(
      data = header_df,
      ggplot2::aes(
        x = x,
        y = 0.40,
        label = label
      ),
      family = font_family,
      fontface = "bold",
      size = 3.2
    ) +
    ggplot2::coord_cartesian(
      xlim = c(0.5, 4.5),
      ylim = c(0, 1)
    ) +
    ggplot2::theme_void()

  # Left strip: one global y-axis title + row labels.
  left_strip <- ggplot2::ggplot() +
    ggplot2::annotate(
      "text",
      x = 0.15,
      y = 1,
      label = "Relative change in abundance (%)",
      angle = 90,
      family = font_family,
      size = 3.0
    ) +
    ggplot2::annotate(
      "text",
      x = 0.65,
      y = 1.50,
      label = "Onset",
      angle = 90,
      family = font_family,
      fontface = "bold",
      size = 3.7
    ) +
    ggplot2::annotate(
      "text",
      x = 0.65,
      y = 0.50,
      label = "Offset",
      angle = 90,
      family = font_family,
      fontface = "bold",
      size = 3.7
    ) +
    ggplot2::coord_cartesian(
      xlim = c(0, 1),
      ylim = c(0, 2)
    ) +
    ggplot2::theme_void()

  figure <- (
    patchwork::plot_spacer() | header
  ) /
    (
      left_strip | body
    ) +
    patchwork::plot_layout(
      heights = c(0.14, 1),
      widths = c(0.075, 1),
      guides = "collect"
    ) &
    ggplot2::theme(
      legend.position = "bottom",
      legend.direction = "horizontal"
    )

  dir.create(
    out_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  source_path <- file.path(
    out_dir,
    paste0(prefix, "_source_data.csv")
  )

  verification_path <- file.path(
    out_dir,
    paste0(prefix, "_verification.csv")
  )

  utils::write.csv(
    predictions,
    source_path,
    row.names = FALSE
  )

  utils::write.csv(
    verification,
    verification_path,
    row.names = FALSE
  )

  png_path <- file.path(
    out_dir,
    paste0(prefix, ".png")
  )

  pdf_path <- file.path(
    out_dir,
    paste0(prefix, ".pdf")
  )

  ggplot2::ggsave(
    png_path,
    figure,
    width = width_mm,
    height = height_mm,
    units = "mm",
    dpi = dpi,
    bg = "white"
  )

  ggplot2::ggsave(
    pdf_path,
    figure,
    width = width_mm,
    height = height_mm,
    units = "mm",
    device = grDevices::cairo_pdf,
    bg = "white"
  )

  message("Saved:")
  message("  ", png_path)
  message("  ", pdf_path)
  message("  ", source_path)
  message("  ", verification_path)

  invisible(
    list(
      figure = figure,
      predictions = predictions,
      verification = verification,
      png = png_path,
      pdf = pdf_path,
      source_data = source_path,
      verification_file = verification_path
    )
  )
}

# Run interactively:
mundlak_8panel_figure <- plot_mundlak_8panel()
