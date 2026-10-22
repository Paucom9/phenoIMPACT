# phenoIMPACT: final environmental models, within-species plasticity curves.
# No fitting, optimization or likelihood evaluation. Load ONE saved model at a time.
# Uses the LOCAL post_environment_mundlak.R selector, including documented reviews.
# Source this file, then run: env_figure <- plot_environment_mundlak_nature()

plot_environment_mundlak_nature <- function(
    project_root = "E:/phenoIMPACT project/code/phenoIMPACT",
    environment_run = file.path(project_root, "output", "phenology_plasticity",
      "environment_mundlak_by_voltinism_qc_final_15km", "run_de0890c85d4a"),
    out_dir = file.path(environment_run, "figures_environment_nature"),
    font_family = "Arial", width_mm = 180, height_mm = 155,
    base_size = 9, dpi = 600, gradient_limit = 2,
    range_probs = c(.025, .975), show_plot = interactive()) {

  need <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)
  pkgs <- c("ggplot2", "patchwork", "lme4")
  need(all(vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)),
       paste("Required packages:", paste(pkgs, collapse = ", ")))
  need(length(range_probs) == 2 && all(is.finite(range_probs)) &&
       0 <= range_probs[1] && range_probs[1] < range_probs[2] && range_probs[2] <= 1,
       "range_probs must be two increasing probabilities between 0 and 1.")
  need(is.finite(gradient_limit) && gradient_limit > 0, "Invalid gradient_limit.")
  need(capabilities("cairo"), "Cairo is required for vector PDF export.")
  selector <- file.path(project_root, "post_environment_mundlak.R")
  need(file.exists(selector), "Local post_environment_mundlak.R not found.")
  e <- new.env(parent = globalenv())
  sys.source(selector, envir = e)
  pf <- e$pheno_post_env$preflight(project_root = project_root,
                                  environment_run = environment_run)
  need(!isTRUE(pf$active), "An environmental worker is active; export after it finishes.")
  sources <- pf$environment
  analyses <- c("onset_univoltine", "onset_multivoltine",
                "offset_univoltine", "offset_multivoltine")
  need(all(analyses %in% sources$analysis), "Incomplete source selection.")
  sources <- sources[match(analyses, sources$analysis), , drop = FALSE]
  need(all(sources$ready), "Some selected environmental sources are not accepted.")
  print(sources[c("analysis", "kind", "refinement_version", "status", "fit")], row.names = FALSE)
  mods <- c("photoperiod", "background", "predictability", "trend")
  labels <- c("Photoperiod", "Background temperature", "Temperature predictability", "Warming trend")
  colours <- setNames(c("#009E73", "#E69F00", "#7570B3", "#D55E00"), labels)
  curves <- tests <- references <- list()

  for (i in seq_along(analyses)) {
    a <- analyses[i]
    message("Reading selected source: ", a, " (one model in memory)")
    p <- readRDS(sources$prepared[i])
    saved <- readRDS(sources$fit[i])
    fit <- if (inherits(saved, "sdmTMB")) saved else saved$fit
    need(inherits(fit, "sdmTMB"), paste("Invalid fit:", a))
    if (!is.null(saved$analysis)) need(identical(saved$analysis, a), "Analysis identity mismatch.")
    original <- readRDS(sources$original_marker[i])
    signature <- if (sources$kind[i] == "refined") saved$source_signature else saved$signature
    need(!is.null(signature) && identical(signature, original$signature),
         paste("Source provenance mismatch:", a))
    form <- p$formula
    need(nrow(fit$data) == nrow(p$data), "Fit/prepared row counts differ.")
    vars <- unique(all.vars(lme4::nobars(form)))
    for (v in vars) need(isTRUE(all.equal(fit$data[[v]], p$data[[v]],
                                         check.attributes = FALSE)),
                         paste("Fit/prepared data mismatch:", v))
    terms <- colnames(model.matrix(lme4::nobars(form), p$data[1, , drop = FALSE]))
    ix <- which(names(fit$sd_report$par.fixed) == "b_j")
    need(length(ix) == length(terms), "Cannot align fixed coefficients.")
    need(identical(names(fit$model$par), names(fit$sd_report$par.fixed)) &&
         isTRUE(all.equal(as.numeric(fit$model$par),
                          as.numeric(fit$sd_report$par.fixed), tolerance = 1e-6)),
         "Covariance and point estimates are not synchronized.")
    beta <- setNames(as.numeric(fit$sd_report$par.fixed[ix]), terms)
    V <- as.matrix(fit$sd_report$cov.fixed[ix, ix, drop = FALSE])
    dimnames(V) <- list(terms, terms)
    need(all(is.finite(beta)) && all(is.finite(V)) && all(diag(V) > 0),
         "Invalid fixed-effect covariance.")
    need(min(eigen((V + t(V))/2, symmetric = TRUE, only.values = TRUE)$values) >
           -1e-8 * max(diag(V)), "Fixed covariance is not positive semidefinite.")
    anomaly <- if (startsWith(a, "onset_")) "clim_anomaly_tw60" else "clim_anomaly_tw90"
    interaction <- function(v) {
      hit <- intersect(c(paste(anomaly, v, sep = ":"), paste(v, anomaly, sep = ":")), terms)
      need(length(hit) == 1, paste("Missing/ambiguous interaction:", v))
      hit
    }
    envvars <- as.vector(outer(mods, c("within_species", "species_mean"), paste, sep = "_"))
    need(all(envvars %in% names(p$data)), "Missing environmental components.")
    pop <- unique(p$data[c("SPECIES", "SITE_ID", envvars)])
    need(!anyDuplicated(pop[c("SPECIES", "SITE_ID")]),
         "Environmental components vary within a population.")
    # Reference: within-species deviations = 0; between-species context =
    # equal-population mean in this fitted group. Avoids arbitrary inherited zeros.
    d0 <- setNames(rep(0, length(beta)), terms)
    d0[anomaly] <- 1
    for (m in mods) {
      v <- paste0(m, "_species_mean")
      ref <- mean(pop[[v]])
      d0[interaction(v)] <- ref
      references[[length(references)+1L]] <- data.frame(analysis = a, variable = v, value = ref)
    }
    for (j in seq_along(mods)) {
      m <- mods[j]; v <- paste0(m, "_within_species"); term <- interaction(v)
      sc <- sd(pop[[v]])
      need(is.finite(sc) && sc > 0, "Zero within-species component SD.")
      limits <- as.numeric(quantile(pop[[v]] / sc, range_probs, names = FALSE))
      limits <- c(max(-gradient_limit, limits[1]), min(gradient_limit, limits[2]))
      need(limits[1] < limits[2], "No supported gradient to plot.")
      x <- sort(unique(c(seq(limits[1], limits[2], length.out = 201),
                          if (limits[1] <= 0 && limits[2] >= 0) 0)))
      D <- matrix(rep(d0, each = length(x)), nrow = length(x), dimnames = list(NULL, terms))
      D[, term] <- x * sc
      est <- as.vector(D %*% beta)
      variance <- rowSums((D %*% V) * D)
      need(all(variance >= -1e-8), "Negative prediction variance.")
      se <- sqrt(pmax(variance, 0))
      pv <- 2 * pnorm(-abs(beta[term] / sqrt(V[term, term])))
      symbol <- if (pv < .001) "***" else if (pv < .01) "**" else if (pv < .05) "*" else "n.s."
      curves[[length(curves)+1L]] <- data.frame(analysis = a, moderator = labels[j],
        x = x, slope = est, SE = se, lower = est - 1.96*se, upper = est + 1.96*se,
        p_Wald = pv, line_type = if (pv < .05) "Supported" else "Not supported",
        component_SD = sc)
      tests[[length(tests)+1L]] <- data.frame(analysis = a, moderator = labels[j],
        term = term, estimate = beta[term], SE = sqrt(V[term, term]),
        p_Wald = pv, symbol = symbol, component_SD = sc)
    }
    rm(p, saved, fit, V, pop, D); invisible(gc())
  }
  dat <- do.call(rbind, curves); tst <- do.call(rbind, tests)
  dat$moderator <- factor(dat$moderator, levels = labels)
  dat$line_type <- factor(dat$line_type, levels = c("Supported", "Not supported"))
  plots <- list()
  rowlimits <- lapply(c("onset_", "offset_"), function(prefix) {
    d <- dat[startsWith(dat$analysis, prefix), ]
    lim <- range(c(d$lower, d$upper, 0)); lim + c(-1,1)*diff(lim)*.07
  })
  for (i in seq_along(analyses)) {
    a <- analyses[i]; z <- dat[dat$analysis == a, ]; tt <- tst[tst$analysis == a, ]
    lim <- rowlimits[[if (i <= 2) 1 else 2]]
    tt$x <- min(z$x); tt$y <- lim[1] + diff(lim) * (.06 + .065 * (3:0))
    tt$annotation <- paste(tt$moderator, tt$symbol)
    plots[[i]] <- ggplot2::ggplot(z, ggplot2::aes(x, slope, colour = moderator, fill = moderator)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = .36) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha = .10,
                           colour = NA, show.legend = FALSE) +
      ggplot2::geom_line(ggplot2::aes(linetype = line_type), linewidth = .6) +
      ggplot2::geom_text(data = tt, ggplot2::aes(x = x, y = y, label = annotation,
         colour = moderator), inherit.aes = FALSE, hjust = 0, family = font_family,
         size = 2.8, show.legend = FALSE) +
      ggplot2::scale_colour_manual(values = colours, drop = FALSE) +
      ggplot2::scale_fill_manual(values = colours, drop = FALSE, guide = "none") +
      ggplot2::scale_linetype_manual(values = c(Supported = "solid", `Not supported` = "dashed"),
                                    guide = "none") +
      ggplot2::scale_x_continuous(breaks = -2:2) +
      ggplot2::coord_cartesian(xlim = c(-gradient_limit, gradient_limit), ylim = lim) +
      ggplot2::labs(tag = letters[i], colour = NULL,
        title = if (i <= 2) c("Univoltine", "Multivoltine")[i] else NULL,
        x = if (i > 2) "Within-species environmental gradient (SD)" else NULL,
        y = if (i %% 2 == 1) paste0(if (i == 1) "Onset" else "Offset",
          " response\n(days per SD of temperature anomaly)") else NULL) +
      ggplot2::theme_classic(base_family = font_family, base_size = base_size) +
      ggplot2::theme(axis.line = ggplot2::element_line(linewidth = .36),
        axis.ticks = ggplot2::element_line(linewidth = .36),
        plot.title = ggplot2::element_text(face = "bold", hjust = .5, size = base_size+1),
        plot.tag = ggplot2::element_text(face = "bold", size = base_size),
        legend.position = "bottom", legend.key.width = grid::unit(12,"mm"),
        plot.margin = ggplot2::margin(5,7,5,5)) +
      ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2,
                                                   override.aes = list(linetype = "solid")))
  }
  figure <- patchwork::wrap_plots(plots, ncol = 2, guides = "collect") &
    ggplot2::theme(legend.position = "bottom")
  # A new subfolder on each execution preserves previous exports.
  out <- file.path(out_dir, paste0("export_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  need(!dir.exists(out), "Export folder already exists; retry in a second.")
  dir.create(out, recursive = TRUE)
  csv <- function(x, name) write.csv(x, file.path(out, name), row.names = FALSE)
  csv(dat, "figure_source_data.csv"); csv(tst, "interaction_tests.csv")
  csv(do.call(rbind, references), "reference_environment.csv")
  csv(sources, "selected_environmental_sources.csv")
  ggplot2::ggsave(file.path(out, "environmental_plasticity.pdf"), figure,
    width = width_mm, height = height_mm, units = "mm", device = grDevices::cairo_pdf, bg = "white")
  ggplot2::ggsave(file.path(out, "environmental_plasticity.png"), figure,
    width = width_mm, height = height_mm, units = "mm", dpi = dpi, bg = "white")
  writeLines(c(
    "Conditional fixed-effect temperature-anomaly slopes; within-species moderators only.",
    "X: population deviation from species mean divided by its equal-population SD within this analysis.",
    "Each curve varies one moderator; other within-species components = 0.",
    "Between-species components held at equal-population means in each analysis.",
    "Y: days per ONE UNIT of the originally standardized anomaly, NOT days per degree C.",
    "Inherited anomaly standardization is not necessarily identical across windows/groups.",
    "Bands: pointwise 95% normal Wald CIs, using the FULL fixed coefficient covariance.",
    "These are not prediction intervals or population-specific random-slope intervals.",
    "Curves restricted to central observed marginal ranges and gradient_limit; joint support is not guaranteed.",
    "Solid/dashed: unadjusted two-sided Wald interaction p < .05 / >= .05.",
    "Stars refer to interactions, not conditional slopes at individual X values.",
    "All curves share the same reference slope at X=0 within each model.",
    "Offset response is conditional on annual onset; its additive term does not enter the anomaly derivative.",
    "Local source selector acceptance is respected, including documented numerical reviews.",
    "No model fitted, optimized, or evaluated; originals unchanged.",
    "Sensitivity analyses and propagation of plasticity uncertainty remain separate work."
  ), file.path(out, "figure_notes.txt"))
  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
  if (show_plot) print(figure)
  message("Figure saved to: ", normalizePath(out, winslash = "/"))
  invisible(list(figure = figure, curves = dat, tests = tst, sources = sources, output_dir = out))
}
