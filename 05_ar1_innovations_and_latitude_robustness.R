# =============================================================================
# phenoIMPACT | Final lightweight robustness checks
# Version 1.0 | 2026-10-08
#
# NO REFITTING.
#
# This script performs two final checks for the all-population geography model:
#
#   A) Latent AR(1) innovation diagnostic
#      Uses the fitted conditional modes of the population AR(1) process and
#      removes the fitted rho^gap dependence exactly:
#
#        innovation_t =
#          (u_t - rho^gap * u_prev) /
#          (sd_AR * sqrt(1 - rho^(2*gap)))
#
#      The first observed year of each population is standardized by sd_AR.
#
#      IMPORTANT:
#      u_t are CONDITIONAL MODES, not posterior draws. Therefore innovation
#      variance is expected to be shrunken and this is NOT a formal N(0,1) test.
#      The useful diagnostic is whether strong serial correlation remains after
#      the fitted AR(1) transition is removed.
#
#   B) Independent empirical latitude robustness
#      Calculates a simple log-abundance temporal slope for every population
#      directly from the annual abundance indices:
#
#        log(ABUND_INDEX) ~ year_decade
#
#      Then tests the within-species relationship of those empirical slopes with
#      latitude + longitude, using equal total weight per species and a
#      species-cluster bootstrap.
#
#      Repeated for:
#        - populations with >=5 or >=10 observed years
#        - full geographic support
#        - central 98% of within-species latitude
#        - central 90% of within-species latitude
#
#      This second analysis is DESCRIPTIVE / ROBUSTNESS ONLY. It deliberately
#      ignores the Gamma likelihood, AR(1), random effects and Matern field.
#      It therefore does not replace the main TMB result; it asks whether the
#      north-south signal is also visible in a much simpler, independent summary
#      and whether geographic extremes drive it.
#
# RUN
#   source(file.choose(), encoding = "UTF-8")
#   final_robust <- run_allpop_final_light_checks()
#
# SEND BACK
#   final_light_checks_to_review.zip
# =============================================================================


run_allpop_final_light_checks <- function(
  project_root = "E:/phenoIMPACT project/code/phenoIMPACT",
  analysis_run = NULL,
  n_boot = 2000L,
  seed = 20261008L,
  make_figures = TRUE
) {

  need <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)
  csv <- function(x, p) utils::write.csv(x, p, row.names = FALSE, na = "")

  required <- c("data.table", "ggplot2")
  missing <- required[
    !vapply(required, requireNamespace, logical(1), quietly = TRUE)
  ]
  need(
    !length(missing),
    paste(
      "Missing package(s):", paste(missing, collapse = ", "),
      "\nInstall them manually before running this script."
    )
  )

  project_root <- normalizePath(
    project_root, winslash = "/", mustWork = TRUE
  )

  root <- file.path(
    project_root,
    "output", "population_trends",
    "independent_geography_all_populations_15km"
  )
  need(dir.exists(root), paste("Analysis root not found:", root))

  # ---------------------------------------------------------------------------
  # Locate latest completed FIT_OK run
  # ---------------------------------------------------------------------------

  if (is.null(analysis_run)) {

    dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)

    good <- dirs[vapply(dirs, function(d) {

      req <- file.path(
        d,
        c("fit.rds", "prepared_all_populations.rds", "workflow_status.csv")
      )
      if (!all(file.exists(req))) return(FALSE)

      st <- tryCatch(
        utils::read.csv(
          file.path(d, "workflow_status.csv"),
          stringsAsFactors = FALSE
        ),
        error = function(e) NULL
      )

      fit <- tryCatch(
        readRDS(file.path(d, "fit.rds")),
        error = function(e) NULL
      )

      !is.null(st) &&
        nrow(st) == 1L &&
        grepl("^FIT_OK", as.character(st$status[1])) &&
        !is.null(fit) &&
        isTRUE(fit$checks$numerical_checks_ok)

    }, logical(1))]

    need(
      length(good) > 0L,
      "No completed FIT_OK all-population geography run found."
    )

    mt <- file.info(file.path(good, "fit.rds"))$mtime
    analysis_run <- good[which.max(mt)]

  } else {

    analysis_run <- normalizePath(
      analysis_run, winslash = "/", mustWork = TRUE
    )
  }

  analysis_run <- normalizePath(
    analysis_run, winslash = "/", mustWork = TRUE
  )

  message("Source run: ", analysis_run)

  fit_path <- file.path(analysis_run, "fit.rds")
  prep_path <- file.path(
    analysis_run, "prepared_all_populations.rds"
  )

  need(file.exists(fit_path), "fit.rds missing.")
  need(file.exists(prep_path), "prepared_all_populations.rds missing.")

  fit <- readRDS(fit_path)
  p <- readRDS(prep_path)

  need(
    identical(fit$signature, p$signature),
    "Fit and prepared input signatures disagree."
  )
  need(
    isTRUE(fit$checks$numerical_checks_ok),
    "Saved fit does not pass numerical checks."
  )
  need(
    nrow(p$frame) == length(p$data$y),
    "Prepared frame and response have different lengths."
  )

  # ---------------------------------------------------------------------------
  # Output
  # ---------------------------------------------------------------------------

  out <- file.path(analysis_run, "final_light_checks")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  out <- normalizePath(out, winslash = "/", mustWork = TRUE)

  # ===========================================================================
  # A. LATENT AR(1) INNOVATION DIAGNOSTIC
  # ===========================================================================

  message("A) Checking latent AR(1) innovations...")

  need(
    is.list(fit$parameters) &&
      all(c("ar_obs", "log_sd") %in% names(fit$parameters)),
    "Saved latent AR(1) arrays are unavailable."
  )

  u <- as.numeric(fit$parameters$ar_obs)
  previous <- as.integer(p$data$previous)
  gap <- as.integer(p$data$gap)

  n <- nrow(p$frame)

  need(
    length(u) == n &&
      length(previous) == n &&
      length(gap) == n,
    "AR(1) arrays do not match fitted rows."
  )

  rho <- as.numeric(fit$checks$AR1_rho)
  sd_ar <- exp(as.numeric(fit$parameters$log_sd[5L]))

  need(
    length(rho) == 1L &&
      is.finite(rho) &&
      abs(rho) < 1 &&
      length(sd_ar) == 1L &&
      is.finite(sd_ar) &&
      sd_ar > 0,
    "Invalid fitted AR(1) parameters."
  )

  innovation <- rep(NA_real_, n)

  first <- which(previous < 0L)
  later <- which(previous >= 0L)

  innovation[first] <- u[first] / sd_ar

  prev_R <- previous[later] + 1L
  need(
    all(prev_R >= 1L & prev_R <= n),
    "Invalid previous-row indices."
  )

  r_gap <- rho ^ gap[later]
  denom <- sd_ar * sqrt(1 - r_gap^2)

  need(
    all(is.finite(denom)) &&
      all(denom > 0),
    "Invalid AR(1) innovation denominator."
  )

  innovation[later] <-
    (u[later] - r_gap * u[prev_R]) / denom

  need(
    all(is.finite(innovation)),
    "Non-finite latent innovations."
  )

  ar_frame <- p$frame
  ar_frame$latent_AR1_mode <- u
  ar_frame$standardized_AR1_innovation <- innovation
  ar_frame$previous_row <- previous
  ar_frame$calendar_gap <- gap

  # Relation between consecutive innovations along the observed population
  # sequence. This is the key lightweight serial-dependence diagnostic.
  pair_cor <- stats::cor(
    innovation[prev_R],
    innovation[later]
  )

  # Same diagnostic by true calendar gap.
  by_gap <- do.call(
    rbind,
    lapply(sort(unique(gap[later])), function(gp) {

      ii <- later[gap[later] == gp]
      jj <- previous[ii] + 1L

      if (length(ii) < 20L) return(NULL)

      data.frame(
        calendar_gap_years = gp,
        n_pairs = length(ii),
        correlation_consecutive_innovations =
          stats::cor(innovation[jj], innovation[ii]),
        innovation_mean = mean(innovation[ii]),
        innovation_SD = stats::sd(innovation[ii]),
        stringsAsFactors = FALSE
      )
    })
  )

  ar_summary <- data.frame(
    n_observations = n,
    n_populations = length(first),
    fitted_rho = rho,
    fitted_AR1_SD = sd_ar,
    conditional_mode_innovation_mean = mean(innovation),
    conditional_mode_innovation_SD = stats::sd(innovation),
    conditional_mode_innovation_q025 =
      unname(stats::quantile(innovation, .025)),
    conditional_mode_innovation_median =
      stats::median(innovation),
    conditional_mode_innovation_q975 =
      unname(stats::quantile(innovation, .975)),
    correlation_of_consecutive_innovations = pair_cor,
    interpretation = paste0(
      "Conditional-mode latent-process diagnostic only; ",
      "innovation variance is expected to be shrunken. ",
      "Serial correlation is the key quantity."
    ),
    stringsAsFactors = FALSE
  )

  csv(ar_summary, file.path(out, "AR1_innovation_summary.csv"))
  csv(by_gap, file.path(out, "AR1_innovation_by_gap.csv"))

  # Do not export the full 487k-row innovation table unless needed.
  # Export compact population-level summaries instead.
  dt_ar <- data.table::as.data.table(ar_frame)

  pop_ar <- dt_ar[, .(
    n_observations = .N,
    innovation_mean = mean(standardized_AR1_innovation),
    innovation_SD = if (.N > 1L)
      stats::sd(standardized_AR1_innovation) else NA_real_
  ), by = .(pop_id, SPECIES, SITE_ID)]

  csv(
    as.data.frame(pop_ar),
    file.path(out, "AR1_innovation_population_summary.csv")
  )

  # ===========================================================================
  # B. EMPIRICAL POPULATION TREND ROBUSTNESS
  # ===========================================================================

  message("B) Calculating empirical population slopes...")

  d <- data.table::as.data.table(p$frame)

  need(
    all(c(
      "pop_id", "SPECIES", "SITE_ID", "year_decade",
      "ABUND_INDEX",
      "latitude_deg", "longitude_deg",
      "latitude_within_species_deg",
      "longitude_within_species_deg"
    ) %in% names(d)),
    "Prepared frame lacks columns needed for empirical-slope sensitivity."
  )

  need(
    all(is.finite(d$ABUND_INDEX)) &&
      all(d$ABUND_INDEX > 0),
    "Invalid abundance indices in prepared frame."
  )

  d[, log_abundance := log(ABUND_INDEX)]

  # Closed-form OLS slope per population:
  # slope = cov(x,y) / var(x)
  empirical <- d[, {
    nn <- .N
    sx <- sum(year_decade)
    sy <- sum(log_abundance)
    sxx <- sum(year_decade^2)
    sxy <- sum(year_decade * log_abundance)

    den <- sxx - sx^2 / nn

    list(
      n_years = data.table::uniqueN(year_num),
      first_year = min(year_num),
      last_year = max(year_num),
      empirical_log_slope_per_decade =
        if (is.finite(den) && den > 0)
          (sxy - sx * sy / nn) / den
        else NA_real_,
      latitude_deg = latitude_deg[1L],
      longitude_deg = longitude_deg[1L],
      latitude_within_species_deg =
        latitude_within_species_deg[1L],
      longitude_within_species_deg =
        longitude_within_species_deg[1L]
    )
  }, by = .(pop_id, SPECIES, SITE_ID)]

  empirical <- empirical[
    is.finite(empirical_log_slope_per_decade) &
      is.finite(latitude_within_species_deg) &
      is.finite(longitude_within_species_deg)
  ]

  need(
    nrow(empirical) > 1000L,
    "Unexpectedly few empirical population slopes."
  )

  csv(
    as.data.frame(empirical),
    file.path(out, "empirical_population_slopes.csv")
  )

  # ---------------------------------------------------------------------------
  # Within-species weighted projection helpers
  # ---------------------------------------------------------------------------

  centre_subset <- function(x) {

    z <- data.table::copy(x)

    z[, `:=`(
      lat_w = latitude_deg - mean(latitude_deg),
      lon_w = longitude_deg - mean(longitude_deg),
      slope_w =
        empirical_log_slope_per_decade -
        mean(empirical_log_slope_per_decade)
    ), by = SPECIES]

    # Species represented by only one retained population carry no within-species
    # geographic information and are removed.
    z[, n_species_pops := .N, by = SPECIES]
    z <- z[n_species_pops >= 2L]

    z
  }

  crossproducts_by_species <- function(z) {

    z[, {
      X <- cbind(
        latitude = lat_w / 10,
        longitude = lon_w / 10
      )

      # Equal total weight per species.
      w <- rep(1 / .N, .N)

      list(
        A11 = sum(w * X[,1]^2),
        A12 = sum(w * X[,1] * X[,2]),
        A22 = sum(w * X[,2]^2),
        c1 = sum(w * X[,1] * slope_w),
        c2 = sum(w * X[,2] * slope_w),
        n_pops = .N
      )
    }, by = SPECIES]
  }

  solve_parts <- function(parts, counts = NULL) {

    if (is.null(counts)) counts <- rep(1, nrow(parts))

    A <- matrix(
      c(
        sum(counts * parts$A11),
        sum(counts * parts$A12),
        sum(counts * parts$A12),
        sum(counts * parts$A22)
      ),
      nrow = 2, byrow = TRUE
    )

    cc <- c(
      sum(counts * parts$c1),
      sum(counts * parts$c2)
    )

    if (qr(A)$rank < 2L) return(c(latitude = NA_real_, longitude = NA_real_))

    b <- as.numeric(solve(A, cc))
    names(b) <- c("latitude", "longitude")
    b
  }

  bootstrap_species <- function(parts, n_boot, seed) {

    set.seed(seed)

    ns <- nrow(parts)
    need(ns >= 10L, "Too few species for cluster bootstrap.")

    outb <- matrix(
      NA_real_, n_boot, 2L,
      dimnames = list(NULL, c("latitude", "longitude"))
    )

    for (b in seq_len(n_boot)) {

      draw <- sample.int(ns, ns, replace = TRUE)
      counts <- tabulate(draw, nbins = ns)

      outb[b, ] <- solve_parts(parts, counts)
    }

    outb <- outb[stats::complete.cases(outb), , drop = FALSE]

    need(
      nrow(outb) >= .95 * n_boot,
      "More than 5% of bootstrap replicates were rank deficient."
    )

    outb
  }

  # ---------------------------------------------------------------------------
  # Scenarios: minimum years x geographic trimming
  # ---------------------------------------------------------------------------

  scenarios <- list()

  for (min_years in c(5L, 10L)) {

    base <- empirical[n_years >= min_years]

    need(
      nrow(base) > 1000L,
      paste("Too few populations with >=", min_years, "years.")
    )

    # Trimming is based on the ORIGINAL within-species latitude coordinate,
    # before re-centring the retained subset.
    q01 <- stats::quantile(
      base$latitude_within_species_deg, c(.01, .99), na.rm = TRUE
    )
    q05 <- stats::quantile(
      base$latitude_within_species_deg, c(.05, .95), na.rm = TRUE
    )

    defs <- list(
      full_support = rep(TRUE, nrow(base)),
      central_98pct =
        base$latitude_within_species_deg >= q01[1] &
        base$latitude_within_species_deg <= q01[2],
      central_90pct =
        base$latitude_within_species_deg >= q05[1] &
        base$latitude_within_species_deg <= q05[2]
    )

    for (nm in names(defs)) {

      z <- base[defs[[nm]]]
      z <- centre_subset(z)

      parts <- crossproducts_by_species(z)
      est <- solve_parts(parts)

      boot <- bootstrap_species(
        parts,
        n_boot = n_boot,
        seed = seed + min_years * 10L + match(nm, names(defs))
      )

      lat_boot <- boot[, "latitude"]
      lon_boot <- boot[, "longitude"]

      scenarios[[length(scenarios) + 1L]] <- data.frame(
        minimum_years = min_years,
        geographic_support = nm,
        n_populations = nrow(z),
        n_species = data.table::uniqueN(z$SPECIES),
        latitude_estimate_log_slope_per_decade_per_10deg =
          unname(est["latitude"]),
        latitude_lower_95_species_bootstrap =
          unname(stats::quantile(lat_boot, .025)),
        latitude_upper_95_species_bootstrap =
          unname(stats::quantile(lat_boot, .975)),
        latitude_bootstrap_prop_positive =
          mean(lat_boot > 0),
        longitude_estimate_log_slope_per_decade_per_10deg =
          unname(est["longitude"]),
        longitude_lower_95_species_bootstrap =
          unname(stats::quantile(lon_boot, .025)),
        longitude_upper_95_species_bootstrap =
          unname(stats::quantile(lon_boot, .975)),
        longitude_bootstrap_prop_positive =
          mean(lon_boot > 0),
        latitude_growth_factor_ratio_per_10deg =
          exp(unname(est["latitude"])),
        latitude_percent_difference_per_10deg =
          100 * expm1(unname(est["latitude"])),
        stringsAsFactors = FALSE
      )
    }
  }

  scenario_table <- do.call(rbind, scenarios)

  csv(
    scenario_table,
    file.path(out, "empirical_latitude_robustness.csv")
  )

  # ---------------------------------------------------------------------------
  # Leave-one-species-out for primary descriptive scenario
  # >=5 years, full support.
  # ---------------------------------------------------------------------------

  primary <- centre_subset(empirical[n_years >= 5L])
  primary_parts <- crossproducts_by_species(primary)

  loo <- vector("list", nrow(primary_parts))

  for (i in seq_len(nrow(primary_parts))) {

    keep <- seq_len(nrow(primary_parts)) != i
    b <- solve_parts(primary_parts[keep])

    loo[[i]] <- data.frame(
      omitted_species = primary_parts$SPECIES[i],
      latitude_estimate = unname(b["latitude"]),
      longitude_estimate = unname(b["longitude"]),
      stringsAsFactors = FALSE
    )
  }

  loo <- do.call(rbind, loo)

  csv(
    loo,
    file.path(out, "empirical_leave_one_species_out.csv")
  )

  # ---------------------------------------------------------------------------
  # Comparison with main TMB latitude estimate
  # ---------------------------------------------------------------------------

  latitude_file <- file.path(analysis_run, "latitude_effect.csv")
  need(file.exists(latitude_file), "Main latitude_effect.csv missing.")

  main <- utils::read.csv(
    latitude_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  need(nrow(main) == 1L, "Unexpected latitude_effect.csv.")

  comparison <- data.frame(
    analysis = c(
      "Main spatial TMB model",
      paste0(
        "Empirical OLS slopes: ",
        scenario_table$minimum_years,
        "+ yr / ",
        scenario_table$geographic_support
      )
    ),
    estimate = c(
      main$estimate[1],
      scenario_table$latitude_estimate_log_slope_per_decade_per_10deg
    ),
    lower_95 = c(
      main$conf.low[1],
      scenario_table$latitude_lower_95_species_bootstrap
    ),
    upper_95 = c(
      main$conf.high[1],
      scenario_table$latitude_upper_95_species_bootstrap
    ),
    inference = c(
      "Model-based Wald CI",
      rep(
        "Species-cluster bootstrap on empirical population slopes",
        nrow(scenario_table)
      )
    ),
    stringsAsFactors = FALSE
  )

  csv(
    comparison,
    file.path(out, "main_vs_empirical_latitude_comparison.csv")
  )

  # ===========================================================================
  # FIGURES
  # ===========================================================================

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

    # AR1 innovation serial correlation by true calendar gap.
    p_ar <- ggplot2::ggplot(
      by_gap,
      ggplot2::aes(
        x = calendar_gap_years,
        y = correlation_consecutive_innovations
      )
    ) +
      ggplot2::geom_hline(
        yintercept = 0,
        linetype = "dashed",
        colour = "grey45"
      ) +
      ggplot2::geom_point(size = 2, colour = "black") +
      ggplot2::geom_line(linewidth = .8, colour = "#2F5AA8") +
      ggplot2::theme_classic(base_family = font, base_size = 15) +
      ggplot2::labs(
        x = "Calendar gap between observations (years)",
        y = "Correlation of consecutive\nstandardized AR1 innovations"
      )

    ggplot2::ggsave(
      file.path(out, "AR1_innovation_serial_correlation.png"),
      p_ar, width = 7, height = 5.5, dpi = 350, bg = "white"
    )

    # Forest plot comparing main TMB and empirical sensitivities.
    comparison$analysis <- factor(
      comparison$analysis,
      levels = rev(comparison$analysis)
    )

    p_forest <- ggplot2::ggplot(
      comparison,
      ggplot2::aes(x = estimate, y = analysis)
    ) +
      ggplot2::geom_vline(
        xintercept = 0,
        linetype = "dashed",
        colour = "grey45"
      ) +
      ggplot2::geom_errorbarh(
        ggplot2::aes(xmin = lower_95, xmax = upper_95),
        height = .16,
        linewidth = .75
      ) +
      ggplot2::geom_point(size = 2.4) +
      ggplot2::theme_classic(base_family = font, base_size = 14) +
      ggplot2::labs(
        x = "Latitude effect on log abundance trend\nper +10° latitude",
        y = NULL
      )

    ggplot2::ggsave(
      file.path(out, "latitude_extreme_robustness.png"),
      p_forest, width = 9, height = 6.5, dpi = 350, bg = "white"
    )
  }

  # ===========================================================================
  # FINAL SUMMARY
  # ===========================================================================

  primary_emp <- scenario_table[
    scenario_table$minimum_years == 5L &
      scenario_table$geographic_support == "full_support",
    ,
    drop = FALSE
  ]

  trim90 <- scenario_table[
    scenario_table$minimum_years == 5L &
      scenario_table$geographic_support == "central_90pct",
    ,
    drop = FALSE
  ]

  final <- data.frame(
    source_run = analysis_run,
    main_TMB_latitude_estimate = main$estimate[1],
    main_TMB_latitude_low = main$conf.low[1],
    main_TMB_latitude_high = main$conf.high[1],
    fitted_AR1_rho = rho,
    conditional_mode_AR1_innovation_SD =
      stats::sd(innovation),
    consecutive_AR1_innovation_correlation =
      pair_cor,
    empirical_full_latitude_estimate =
      primary_emp$latitude_estimate_log_slope_per_decade_per_10deg,
    empirical_full_latitude_low =
      primary_emp$latitude_lower_95_species_bootstrap,
    empirical_full_latitude_high =
      primary_emp$latitude_upper_95_species_bootstrap,
    empirical_central90_latitude_estimate =
      trim90$latitude_estimate_log_slope_per_decade_per_10deg,
    empirical_central90_latitude_low =
      trim90$latitude_lower_95_species_bootstrap,
    empirical_central90_latitude_high =
      trim90$latitude_upper_95_species_bootstrap,
    leave_one_species_out_latitude_min =
      min(loo$latitude_estimate, na.rm = TRUE),
    leave_one_species_out_latitude_max =
      max(loo$latitude_estimate, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  csv(
    final,
    file.path(out, "FINAL_LIGHT_CHECK_SUMMARY.csv")
  )

  writeLines(
    c(
      "FINAL LIGHTWEIGHT ROBUSTNESS CHECKS — NO REFITTING",
      "",
      paste("Source run:", analysis_run),
      "",
      "A) AR1 innovations",
      "Standardized innovations are calculated from the fitted conditional modes",
      "using the exact rho^gap transition and fitted AR1 SD.",
      "Because these are conditional modes rather than posterior draws, their",
      "variance is expected to be shrunken. The useful quantity is remaining",
      "serial correlation, not whether SD equals exactly 1.",
      "",
      "B) Geographic-extreme robustness",
      "Independent simple population slopes are estimated as",
      "log(ABUND_INDEX) ~ year_decade.",
      "Their within-species latitude + longitude association is then estimated",
      "with equal total weight per species and species-cluster bootstrap intervals.",
      "This is intentionally simpler than the main TMB model and is only a",
      "descriptive robustness check.",
      "Results are repeated after removing the outer 1% or 5% of the original",
      "within-species latitude distribution at each tail.",
      "",
      "The main manuscript inference remains the spatial Gamma/TMB model."
    ),
    file.path(out, "README_final_light_checks.txt")
  )

  # Review archive: compact tables + figures only.
  review_files <- list.files(
    out,
    pattern = "\\.(csv|txt|png)$",
    full.names = FALSE
  )

  zip_path <- file.path(out, "final_light_checks_to_review.zip")

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

  cat("\nAR1 innovation summary:\n")
  print(ar_summary, row.names = FALSE)

  cat("\nEmpirical latitude robustness:\n")
  print(scenario_table, row.names = FALSE)

  cat("\nFILE TO SEND:\n", zip_path, "\n", sep = "")

  invisible(list(
    source_run = analysis_run,
    output_directory = out,
    AR1_summary = ar_summary,
    AR1_by_gap = by_gap,
    empirical_scenarios = scenario_table,
    leave_one_species_out = loo,
    final_summary = final,
    zip_to_send = zip_path
  ))
}
