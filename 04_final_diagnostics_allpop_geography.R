# =============================================================================
# phenoIMPACT | Final residual diagnostics for ALL-population geography model
# Version 1.0 | 2026-10-08
#
# NO REFITTING.
#
# Source model:
#   output/population_trends/independent_geography_all_populations_15km/
#
# The script automatically selects the most recent completed FIT_OK run with:
#   fit.rds
#   prepared_all_populations.rds
#
# Diagnostics:
#   1) numerical-fit audit (saved checks; no arbitrary extra optimization)
#   2) randomized/PIT-normal residual distribution for Gamma/log model
#   3) residual behaviour across fitted-value deciles
#   4) temporal residual correlation by calendar lag after fitted AR(1)
#   5) annual spatial Moran's I after fitted annual Matern field
#
# Residuals are CONDITIONAL on the fitted random effects, AR(1), and annual
# spatial field. They are diagnostics, not new inferential tests.
#
# RUN
#   source(file.choose(), encoding = "UTF-8")
#   final_diag <- run_allpop_geography_final_diagnostics()
#
# Then send:
#   final_diagnostics_to_review.zip
# =============================================================================

run_allpop_geography_final_diagnostics <- function(
  project_root = "E:/phenoIMPACT project/code/phenoIMPACT",
  analysis_run = NULL,
  spatial_k = 8L,
  spatial_permutations = 499L,
  seed = 20261008L,
  make_figures = TRUE
) {

  need <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)
  csv <- function(x, p) utils::write.csv(x, p, row.names = FALSE, na = "")

  required <- c("spdep", "ggplot2")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  need(!length(missing),
       paste("Missing package(s):", paste(missing, collapse = ", "),
             "\nInstall them manually before running this diagnostic script."))

  project_root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)

  root <- file.path(
    project_root, "output", "population_trends",
    "independent_geography_all_populations_15km"
  )
  need(dir.exists(root), paste("Analysis root not found:", root))

  # ---------------------------------------------------------------------------
  # Locate the reviewed completed run
  # ---------------------------------------------------------------------------

  if (is.null(analysis_run)) {
    dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)

    ok <- dirs[vapply(dirs, function(d) {
      req <- file.path(d, c(
        "fit.rds",
        "prepared_all_populations.rds",
        "workflow_status.csv"
      ))
      if (!all(file.exists(req))) return(FALSE)

      st <- tryCatch(
        utils::read.csv(file.path(d, "workflow_status.csv"),
                        stringsAsFactors = FALSE),
        error = function(e) NULL
      )
      fit <- tryCatch(readRDS(file.path(d, "fit.rds")),
                      error = function(e) NULL)

      !is.null(st) && nrow(st) == 1L &&
        identical(as.character(st$status[1]), "FIT_OK") &&
        !is.null(fit) &&
        isTRUE(fit$checks$numerical_checks_ok)
    }, logical(1))]

    need(length(ok) > 0L, "No completed FIT_OK all-population geography run found.")

    mt <- file.info(file.path(ok, "fit.rds"))$mtime
    analysis_run <- ok[which.max(mt)]

  } else {
    analysis_run <- normalizePath(analysis_run, winslash = "/", mustWork = TRUE)
  }

  analysis_run <- normalizePath(analysis_run, winslash = "/", mustWork = TRUE)
  message("Diagnostic source run: ", analysis_run)

  fit_path <- file.path(analysis_run, "fit.rds")
  prep_path <- file.path(analysis_run, "prepared_all_populations.rds")

  need(file.exists(fit_path), "fit.rds missing.")
  need(file.exists(prep_path), "prepared_all_populations.rds missing.")

  fit <- readRDS(fit_path)
  p <- readRDS(prep_path)

  need(identical(fit$signature, p$signature),
       "Fit and prepared input signatures disagree.")
  need(isTRUE(fit$checks$numerical_checks_ok),
       "Saved fit does not pass its numerical checks.")
  need(length(fit$eta) == nrow(p$frame),
       "Saved fitted linear predictor and input rows differ.")
  need(length(p$data$y) == nrow(p$frame),
       "Saved response and frame rows differ.")

  # ---------------------------------------------------------------------------
  # Output
  # ---------------------------------------------------------------------------

  out <- file.path(analysis_run, "final_residual_diagnostics")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  out <- normalizePath(out, winslash = "/", mustWork = TRUE)

  # ---------------------------------------------------------------------------
  # 1. Numerical audit
  # ---------------------------------------------------------------------------

  checks <- fit$checks
  numerical <- data.frame(
    check = c(
      "optimizer_convergence_code_zero",
      "positive_definite_Hessian",
      "finite_fixed_effect_covariance",
      "positive_fixed_effect_variances",
      "covariance_scaled_gradient_below_0.05",
      "saved_numerical_checks_ok"
    ),
    value = c(
      as.character(checks$convergence_code),
      as.character(checks$pdHess),
      as.character(all(is.finite(fit$V))),
      as.character(all(diag(fit$V) > 0)),
      format(checks$gradient_covariance_norm, digits = 8),
      as.character(checks$numerical_checks_ok)
    ),
    passed = c(
      checks$convergence_code == 0,
      isTRUE(checks$pdHess),
      all(is.finite(fit$V)),
      all(diag(fit$V) > 0),
      is.finite(checks$gradient_covariance_norm) &&
        checks$gradient_covariance_norm < 0.05,
      isTRUE(checks$numerical_checks_ok)
    ),
    stringsAsFactors = FALSE
  )

  csv(numerical, file.path(out, "numerical_audit.csv"))
  csv(checks, file.path(out, "source_fit_checks.csv"))

  # Important: raw max gradient is exported but is not used as a standalone
  # convergence threshold because gradients are parameter-scale dependent.
  csv(
    data.frame(
      max_abs_gradient = checks$max_abs_gradient,
      covariance_scaled_gradient = checks$gradient_covariance_norm,
      predefined_scaled_gradient_threshold = 0.05,
      forced_refinement_performed = FALSE,
      reason = paste0(
        "No forced optimization: convergence=0, pdHess=TRUE and ",
        "covariance-scaled gradient already satisfies the predefined criterion."
      )
    ),
    file.path(out, "refinement_decision.csv")
  )

  # ---------------------------------------------------------------------------
  # 2. Conditional Gamma PIT-normal residuals
  # ---------------------------------------------------------------------------

  y <- as.numeric(p$data$y)
  eta <- as.numeric(fit$eta)
  shape <- as.numeric(checks$Gamma_shape)

  need(all(is.finite(y)) && all(y > 0), "Invalid Gamma responses.")
  need(all(is.finite(eta)), "Non-finite fitted eta.")
  need(length(shape) == 1L && is.finite(shape) && shape > 0,
       "Invalid fitted Gamma shape.")

  mu <- exp(eta)
  need(all(is.finite(mu)) && all(mu > 0), "Invalid fitted means.")

  # Continuous PIT: no randomization is required for a continuous Gamma model.
  u <- stats::pgamma(y, shape = shape, scale = mu / shape)

  # Protect qnorm() from exact machine-boundary values.
  eps <- 1e-12
  u <- pmin(1 - eps, pmax(eps, u))
  r_pit <- stats::qnorm(u)

  # Pearson and deviance residuals as complementary diagnostics.
  r_pearson <- (y - mu) / (mu / sqrt(shape))
  ratio <- y / mu
  dev_component <- 2 * ((y - mu) / mu - log(ratio))
  dev_component <- pmax(dev_component, 0)
  r_deviance <- sign(y - mu) * sqrt(dev_component)

  need(all(is.finite(r_pit)) &&
       all(is.finite(r_pearson)) &&
       all(is.finite(r_deviance)),
       "Residual calculation produced non-finite values.")

  res <- p$frame
  res$fitted_mean <- mu
  res$pit_normal_residual <- r_pit
  res$pearson_residual <- r_pearson
  res$deviance_residual <- r_deviance

  summary_res <- data.frame(
    n = length(r_pit),
    PIT_mean = mean(r_pit),
    PIT_SD = stats::sd(r_pit),
    PIT_q01 = unname(stats::quantile(r_pit, .01)),
    PIT_q025 = unname(stats::quantile(r_pit, .025)),
    PIT_median = stats::median(r_pit),
    PIT_q975 = unname(stats::quantile(r_pit, .975)),
    PIT_q99 = unname(stats::quantile(r_pit, .99)),
    Pearson_mean = mean(r_pearson),
    Pearson_SD = stats::sd(r_pearson),
    Deviance_mean = mean(r_deviance),
    Deviance_SD = stats::sd(r_deviance),
    cor_PIT_with_log_fitted = stats::cor(r_pit, log(mu)),
    stringsAsFactors = FALSE
  )
  csv(summary_res, file.path(out, "residual_distribution_summary.csv"))

  # Fitted-value bins: equal-count deciles.
  breaks <- unique(stats::quantile(
    log(mu), probs = seq(0, 1, .1), na.rm = TRUE
  ))
  need(length(breaks) >= 3L, "Too little fitted-value variation for diagnostics.")

  res$fitted_decile <- cut(
    log(mu), breaks = breaks, include.lowest = TRUE, ordered_result = TRUE
  )

  binned <- do.call(rbind, lapply(split(seq_len(nrow(res)), res$fitted_decile), function(ii) {
    data.frame(
      fitted_decile = as.character(res$fitted_decile[ii[1]]),
      n = length(ii),
      mean_log_fitted = mean(log(res$fitted_mean[ii])),
      mean_PIT_residual = mean(res$pit_normal_residual[ii]),
      SD_PIT_residual = stats::sd(res$pit_normal_residual[ii]),
      mean_Pearson_residual = mean(res$pearson_residual[ii]),
      SD_Pearson_residual = stats::sd(res$pearson_residual[ii]),
      stringsAsFactors = FALSE
    )
  }))
  rownames(binned) <- NULL
  csv(binned, file.path(out, "residuals_by_fitted_decile.csv"))

  # ---------------------------------------------------------------------------
  # 3. Temporal residual screen
  # ---------------------------------------------------------------------------

  previous <- as.integer(p$data$previous)
  gap <- as.integer(p$data$gap)

  need(length(previous) == nrow(res) && length(gap) == nrow(res),
       "AR1 link arrays do not match residual rows.")

  has_prev <- which(previous >= 0L)
  prev_R <- previous[has_prev] + 1L

  need(all(prev_R >= 1L & prev_R <= nrow(res)),
       "Invalid previous-row indices.")

  temporal <- do.call(rbind, lapply(sort(unique(gap[has_prev])), function(gp) {
    ii <- has_prev[gap[has_prev] == gp]
    jj <- previous[ii] + 1L
    if (length(ii) < 20L) return(NULL)

    data.frame(
      calendar_lag_years = gp,
      n_pairs = length(ii),
      pooled_PIT_residual_correlation =
        stats::cor(r_pit[jj], r_pit[ii]),
      pooled_Pearson_residual_correlation =
        stats::cor(r_pearson[jj], r_pearson[ii]),
      stringsAsFactors = FALSE
    )
  }))

  rownames(temporal) <- NULL
  csv(temporal, file.path(out, "temporal_residual_correlation_by_gap.csv"))

  lag1 <- temporal[temporal$calendar_lag_years == 1, , drop = FALSE]
  if (nrow(lag1) == 1L) {
    csv(lag1, file.path(out, "temporal_lag1_summary.csv"))
  }

  # ---------------------------------------------------------------------------
  # 4. Annual residual spatial autocorrelation
  # ---------------------------------------------------------------------------

  need(all(c("SITE_ID", "year_num") %in% names(res)),
       "Frame lacks SITE_ID/year_num for spatial diagnostics.")
  need(all(c("SITE_ID", "x_km", "y_km") %in% names(p$sites)),
       "Prepared input lacks site coordinates.")

  site_year <- stats::aggregate(
    pit_normal_residual ~ SITE_ID + year_num,
    data = res,
    FUN = mean
  )

  sites <- as.data.frame(p$sites)
  sites$SITE_ID <- as.character(sites$SITE_ID)
  site_year$SITE_ID <- as.character(site_year$SITE_ID)

  site_year <- merge(
    site_year,
    sites[c("SITE_ID", "x_km", "y_km")],
    by = "SITE_ID",
    all.x = TRUE,
    sort = FALSE
  )

  need(!anyNA(site_year$x_km) && !anyNA(site_year$y_km),
       "Site-year residual coordinate join failed.")

  years <- sort(unique(site_year$year_num))
  set.seed(seed)

  moran_rows <- vector("list", length(years))

  for (j in seq_along(years)) {
    yr <- years[j]
    z <- site_year[site_year$year_num == yr, , drop = FALSE]

    # One row per site-year should already be guaranteed.
    need(!anyDuplicated(z$SITE_ID),
         paste("Duplicated site-year residual rows in", yr))

    n <- nrow(z)
    if (n <= spatial_k + 2L) {
      moran_rows[[j]] <- data.frame(
        year = yr, n_sites = n, k = NA_integer_,
        Moran_I = NA_real_, permutation_p_positive = NA_real_
      )
      next
    }

    xy <- as.matrix(z[c("x_km", "y_km")])

    kn <- spdep::knearneigh(xy, k = min(spatial_k, n - 1L))
    nb <- spdep::knn2nb(kn)
    lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)

    mc <- spdep::moran.mc(
      z$pit_normal_residual,
      lw,
      nsim = spatial_permutations,
      alternative = "greater",
      zero.policy = TRUE
    )

    moran_rows[[j]] <- data.frame(
      year = yr,
      n_sites = n,
      k = min(spatial_k, n - 1L),
      Moran_I = as.numeric(mc$statistic),
      permutation_p_positive = as.numeric(mc$p.value)
    )

    message(
      "Spatial residual screen ", j, "/", length(years),
      " | year ", yr,
      " | Moran I = ", round(as.numeric(mc$statistic), 4)
    )
  }

  moran <- do.call(rbind, moran_rows)
  moran$q_BH_positive <- stats::p.adjust(
    moran$permutation_p_positive, method = "BH"
  )
  csv(moran, file.path(out, "annual_spatial_moran.csv"))

  valid_moran <- moran[is.finite(moran$Moran_I), , drop = FALSE]

  spatial_summary <- data.frame(
    n_years_tested = nrow(valid_moran),
    median_Moran_I = stats::median(valid_moran$Moran_I),
    mean_Moran_I = mean(valid_moran$Moran_I),
    min_Moran_I = min(valid_moran$Moran_I),
    max_Moran_I = max(valid_moran$Moran_I),
    n_positive_nominal_p_lt_0_05 = sum(
      valid_moran$Moran_I > 0 &
        valid_moran$permutation_p_positive < .05
    ),
    n_positive_BH_q_lt_0_05 = sum(
      valid_moran$Moran_I > 0 &
        valid_moran$q_BH_positive < .05
    ),
    stringsAsFactors = FALSE
  )
  csv(spatial_summary, file.path(out, "spatial_residual_summary.csv"))

  # ---------------------------------------------------------------------------
  # 5. Figures
  # ---------------------------------------------------------------------------

  if (make_figures) {
    font <- if (.Platform$OS.type == "windows") "Garamond" else "serif"
    if (.Platform$OS.type == "windows") {
      try(
        grDevices::windowsFonts(
          Garamond = grDevices::windowsFont("Garamond")
        ),
        silent = TRUE
      )
    }

    # QQ plot: deterministic thinning only for graphics.
    set.seed(seed)
    qn <- min(100000L, length(r_pit))
    qidx <- if (length(r_pit) > qn) sample.int(length(r_pit), qn) else seq_along(r_pit)
    qq <- data.frame(residual = r_pit[qidx])

    qplot <- ggplot2::ggplot(qq, ggplot2::aes(sample = residual)) +
      ggplot2::stat_qq(alpha = .15, size = .5, colour = "black") +
      ggplot2::stat_qq_line(colour = "#2F5AA8", linewidth = 1) +
      ggplot2::theme_classic(base_family = font, base_size = 15) +
      ggplot2::labs(x = "Theoretical normal quantiles",
                    y = "PIT-normal residual quantiles")
    ggplot2::ggsave(
      file.path(out, "residual_QQ.png"),
      qplot, width = 6.5, height = 5.5, dpi = 350, bg = "white"
    )

    bplot <- ggplot2::ggplot(
      binned,
      ggplot2::aes(x = mean_log_fitted, y = SD_PIT_residual)
    ) +
      ggplot2::geom_hline(
        yintercept = 1, linetype = "dashed", colour = "grey45"
      ) +
      ggplot2::geom_point(size = 2, colour = "black") +
      ggplot2::geom_line(linewidth = .8, colour = "#2F5AA8") +
      ggplot2::theme_classic(base_family = font, base_size = 15) +
      ggplot2::labs(
        x = "Mean log fitted abundance",
        y = "SD of PIT-normal residuals"
      )
    ggplot2::ggsave(
      file.path(out, "residual_dispersion_by_fitted.png"),
      bplot, width = 6.5, height = 5.5, dpi = 350, bg = "white"
    )

    tplot <- ggplot2::ggplot(
      temporal,
      ggplot2::aes(
        x = calendar_lag_years,
        y = pooled_PIT_residual_correlation
      )
    ) +
      ggplot2::geom_hline(
        yintercept = 0, linetype = "dashed", colour = "grey45"
      ) +
      ggplot2::geom_point(size = 2, colour = "black") +
      ggplot2::geom_line(linewidth = .8, colour = "#2F5AA8") +
      ggplot2::theme_classic(base_family = font, base_size = 15) +
      ggplot2::labs(
        x = "Calendar lag (years)",
        y = "Residual correlation"
      )
    ggplot2::ggsave(
      file.path(out, "temporal_residual_correlation.png"),
      tplot, width = 6.5, height = 5.5, dpi = 350, bg = "white"
    )

    mplot <- ggplot2::ggplot(
      valid_moran,
      ggplot2::aes(x = year, y = Moran_I)
    ) +
      ggplot2::geom_hline(
        yintercept = 0, linetype = "dashed", colour = "grey45"
      ) +
      ggplot2::geom_point(
        ggplot2::aes(
          shape = q_BH_positive < .05 & Moran_I > 0
        ),
        size = 2, colour = "black"
      ) +
      ggplot2::geom_line(linewidth = .7, colour = "#2F5AA8") +
      ggplot2::scale_shape_manual(
        values = c(`FALSE` = 1, `TRUE` = 16),
        name = "Positive Moran\nBH q < 0.05"
      ) +
      ggplot2::theme_classic(base_family = font, base_size = 15) +
      ggplot2::labs(
        x = "Year",
        y = "Moran's I of site-year residuals"
      )
    ggplot2::ggsave(
      file.path(out, "annual_spatial_moran.png"),
      mplot, width = 8, height = 5.5, dpi = 350, bg = "white"
    )
  }

  # ---------------------------------------------------------------------------
  # 6. Concise final review summary
  # ---------------------------------------------------------------------------

  latitude_file <- file.path(analysis_run, "latitude_effect.csv")
  latitude <- if (file.exists(latitude_file)) {
    utils::read.csv(latitude_file, stringsAsFactors = FALSE)
  } else NULL

  final_summary <- data.frame(
    source_run = analysis_run,
    numerical_fit_ok = all(numerical$passed),
    scaled_gradient = checks$gradient_covariance_norm,
    PIT_residual_mean = summary_res$PIT_mean,
    PIT_residual_SD = summary_res$PIT_SD,
    lag1_residual_correlation =
      if (nrow(lag1) == 1L) lag1$pooled_PIT_residual_correlation else NA_real_,
    median_annual_Moran_I = spatial_summary$median_Moran_I,
    positive_Moran_years_BH_q_lt_0_05 =
      spatial_summary$n_positive_BH_q_lt_0_05,
    latitude_estimate =
      if (!is.null(latitude)) latitude$estimate[1] else NA_real_,
    latitude_low =
      if (!is.null(latitude)) latitude$conf.low[1] else NA_real_,
    latitude_high =
      if (!is.null(latitude)) latitude$conf.high[1] else NA_real_,
    stringsAsFactors = FALSE
  )
  csv(final_summary, file.path(out, "FINAL_REVIEW_SUMMARY.csv"))

  writeLines(c(
    "FINAL RESIDUAL DIAGNOSTICS — ALL-POPULATION GEOGRAPHY MODEL",
    "",
    paste("Source run:", analysis_run),
    "No model was refitted.",
    "No forced numerical refinement was performed because the saved fit already",
    "satisfied the predefined numerical criterion: convergence code 0, positive",
    "Hessian, finite covariance and covariance-scaled gradient < 0.05.",
    "",
    "Residuals are continuous Gamma PIT values transformed to N(0,1), conditional",
    "on the fitted species/site random effects, population AR1 and annual spatial field.",
    "",
    "Spatial diagnostic: Moran's I of mean PIT-normal residual by site within each year,",
    paste0("using ", spatial_k, "-nearest-neighbour row-standardized weights and ",
           spatial_permutations, " permutations; BH correction across years."),
    "",
    "Temporal diagnostic: pooled correlation of conditional residual pairs by true",
    "calendar lag, using the saved AR1 previous-row and gap arrays.",
    "",
    "These diagnostics assess adequacy of the fitted model; they do not alter the",
    "latitude-effect estimate or propagate additional first-stage uncertainty."
  ), file.path(out, "README_final_diagnostics.txt"))

  # Archive only review-size outputs, not the 487k-row residual table.
  review_files <- list.files(
    out,
    pattern = "\\.(csv|txt|png)$",
    full.names = FALSE
  )

  zip_path <- file.path(out, "final_diagnostics_to_review.zip")

  if (requireNamespace("zip", quietly = TRUE)) {
    zip::zipr(
      zip_path,
      files = review_files,
      root = out,
      mode = "mirror",
      include_directories = FALSE
    )
  } else {
    old <- getwd()
    on.exit(setwd(old), add = TRUE)
    setwd(out)
    utils::zip(zipfile = zip_path, files = review_files)
  }

  cat("\nFinal diagnostic summary:\n")
  print(final_summary, row.names = FALSE)
  cat("\nFILE TO SEND:\n", zip_path, "\n", sep = "")

  invisible(list(
    source_run = analysis_run,
    output_directory = out,
    numerical_audit = numerical,
    residual_summary = summary_res,
    temporal = temporal,
    annual_moran = moran,
    spatial_summary = spatial_summary,
    final_summary = final_summary,
    zip_to_send = zip_path
  ))
}
