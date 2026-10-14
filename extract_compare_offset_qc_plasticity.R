# =============================================================================
# phenoIMPACT | Extract QC offset plasticities and compare with pre-QC fits
# Version 1.0 | 2026-10-07
#
# PURPOSE
#   1) Extract contextual population-level OFFSET plasticity from the FINAL
#      QC-filtered 15-km offset models.
#   2) Re-extract the SAME quantity from the immediately previous onset-adjusted
#      15-km offset models (pre-QC).
#   3) Compare OLD vs QC before any abundance-model refit.
#
# NO MODEL IS FITTED.
#
# Sources:
#   OLD:
#     output/phenology_plasticity/offset_with_onset_15km/run_3911d7bc87ac
#   QC:
#     output/phenology_plasticity/offset_with_onset_qc_lastpeak_15km/
#       run_fbf8b3101079
#
# Plasticity definition:
#   population contextual anomaly slope =
#       fixed anomaly coefficient
#     + species random anomaly slope
#     + sum(interaction coefficient * population mean moderator)
#
# Annual ONSET_mean_z is additive, so it does NOT enter the derivative directly.
#
# Biological orientation used additionally for comparison:
#   univoltine: higher = stronger offset advancement  -> -raw slope
#   multivoltine: higher = stronger offset delay      ->  raw slope
#
# The z comparison is standardized separately for OLD and QC over the SAME
# matched-population set within each voltinism group.
#
# Outputs:
#   comparison_summary.csv
#   population_comparison_all.csv
#   top_changed_populations.csv
#   unmatched_populations.csv
#   old_population_plasticity_*.csv
#   qc_population_plasticity_*.csv
#
# Run:
#   source(file.choose(), encoding = "UTF-8")
# =============================================================================

compare_offset_qc_plasticity <- function(
  root = "E:/phenoIMPACT project/code/phenoIMPACT",
  old_run = file.path(
    root, "output", "phenology_plasticity",
    "offset_with_onset_15km", "run_3911d7bc87ac"
  ),
  qc_run = file.path(
    root, "output", "phenology_plasticity",
    "offset_with_onset_qc_lastpeak_15km", "run_fbf8b3101079"
  ),
  out_dir = file.path(
    root, "output", "phenology_plasticity",
    "offset_with_onset_qc_lastpeak_15km", "run_fbf8b3101079",
    "plasticity_comparison_pre_abundance"
  ),
  top_n = 30L
) {

  pkgs <- c("sdmTMB", "lme4", "dplyr")
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) {
    stop("Missing packages: ", paste(miss, collapse = ", "), call. = FALSE)
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  groups <- c("univoltine", "multivoltine")
  analyses <- paste0("offset_", groups)

  need <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)
  same_num <- function(a, b, tol = 1e-7) {
    length(a) == length(b) &&
      all(is.finite(a)) && all(is.finite(b)) &&
      (length(a) == 0L || max(abs(as.numeric(a) - as.numeric(b))) <= tol)
  }
  write_csv <- function(x, path) {
    utils::write.csv(x, path, row.names = FALSE, na = "")
  }

  specification <- function(analysis) {
    need(analysis %in% analyses, "Unknown analysis.")
    group <- sub("^offset_", "", analysis)
    list(
      analysis = analysis,
      group = group,
      response = "OFFSET_mean",
      anomaly = "clim_anomaly_tw90",
      moderators = c(
        "photo_tw90",
        "clim_background_tw90",
        "clim_predictability_tw90",
        "clim_trend_tw90"
      )
    )
  }

  numerical_check <- function(fit, label) {
    need(inherits(fit, "sdmTMB"), paste(label, "is not an sdmTMB fit."))
    need(identical(fit$family$family, "gaussian") &&
           identical(fit$family$link, "identity"),
         paste(label, "has unexpected family/link."))
    need(isTRUE(fit$model$convergence == 0L),
         paste(label, "optimizer convergence code is not 0."))
    need(isTRUE(fit$sd_report$pdHess),
         paste(label, "Hessian is not positive definite."))
    need(length(fit$gradients) > 0L &&
           all(is.finite(fit$gradients)) &&
           max(abs(fit$gradients)) < .001,
         paste(label, "does not pass gradient check."))
    need(
      identical(names(fit$model$par), names(fit$sd_report$par.fixed)) &&
        same_num(fit$model$par, fit$sd_report$par.fixed, 1e-6),
      paste(label, "saved covariance/point estimates are stale or mismatched.")
    )
    invisible(TRUE)
  }

  extract_coefficients <- function(fit, spec) {
    sf <- fit$split_formula[[1]]

    fixed_names <- colnames(stats::model.matrix(
      sf$form_no_bars,
      fit$data[seq_len(min(3L, nrow(fit$data))), , drop = FALSE]
    ))
    ix <- which(names(fit$sd_report$par.fixed) == "b_j")

    need(length(ix) == length(fixed_names) && !anyDuplicated(fixed_names),
         "Fixed-coefficient mapping failed.")

    beta <- stats::setNames(
      as.numeric(fit$sd_report$par.fixed[ix]),
      fixed_names
    )

    V <- as.matrix(fit$sd_report$cov.fixed[ix, ix, drop = FALSE])
    need(all(is.finite(beta)) && all(is.finite(V)) && all(diag(V) > 0),
         "Invalid fixed estimates/covariance.")
    dimnames(V) <- list(fixed_names, fixed_names)

    # Extract one random anomaly slope per species using the mapping stored
    # inside the sdmTMB fit; never assume alphabetical factor ordering.
    blocks <- sf$re_cov_terms$re_b_df
    need(all(c("start", "end", "group_indices", "level") %in% names(blocks)),
         "Random-effect mapping is missing.")

    all_ix <- unlist(
      Map(function(a, b) seq.int(a, b), blocks$start, blocks$end),
      use.names = FALSE
    )

    nm <- names(fit$sd_report$value)
    values <- fit$sd_report$value[
      grepl("^re_b_pars($|\\[)", nm)
    ]
    errors <- fit$sd_report$sd[
      grepl("^re_b_pars($|\\[)", nm)
    ]

    active <- !is.na(fit$tmb_map$re_b_pars)

    need(
      length(values) == length(active) &&
        length(errors) == length(active) &&
        identical(as.integer(all_ix + 1L), which(as.vector(active))),
      "Random-effect report does not match saved parameter mapping."
    )

    group_index <- match("SPECIES_slope", names(sf$re_cov_terms$cnms))
    need(!is.na(group_index), "SPECIES_slope random effect not found.")

    b <- blocks[blocks$group_indices == group_index, , drop = FALSE]
    need(nrow(b) > 0L &&
           all(b$start == b$end) &&
           !anyDuplicated(as.character(b$level)),
         "Expected exactly one random anomaly slope per species.")

    slopes <- data.frame(
      SPECIES = as.character(b$level),
      species_slope_deviation = as.numeric(values[b$start + 1L]),
      species_slope_marginal_SE = as.numeric(errors[b$start + 1L]),
      stringsAsFactors = FALSE
    )

    need(
      setequal(slopes$SPECIES, unique(as.character(fit$data$SPECIES))) &&
        all(is.finite(slopes$species_slope_deviation)),
      "Species random slopes/IDs are invalid."
    )

    list(beta = beta, V = V, species_slopes = slopes)
  }

  validate_formula <- function(fit, spec, label) {
    sf <- fit$split_formula[[1]]
    fixed_terms <- attr(stats::terms(sf$form_no_bars), "term.labels")

    expected <- stats::as.formula(paste(
      "OFFSET_mean ~ ONSET_mean_z + clim_anomaly_tw90 *",
      "(photo_tw90 + clim_background_tw90 +",
      "clim_predictability_tw90 + clim_trend_tw90)"
    ))
    expected_terms <- attr(stats::terms(expected), "term.labels")

    need(setequal(fixed_terms, expected_terms),
         paste(label, "fixed-effect formula is not the expected onset-adjusted model."))

    cn <- sf$re_cov_terms$cnms
    need(
      length(cn) == 4L &&
        setequal(names(cn),
                 c("SITE_ID", "SPECIES", "SPECIES_slope", "site_year_id")) &&
        identical(cn$SPECIES_slope, spec$anomaly),
      paste(label, "random-effects specification is unexpected.")
    )

    invisible(TRUE)
  }

  population_contextual_slopes <- function(fit, coef, spec) {
    d <- as.data.frame(fit$data)

    needed <- c(
      "SPECIES", "SITE_ID", "YEAR", "bms_id",
      spec$anomaly, spec$moderators
    )
    need(all(needed %in% names(d)),
         "Required fitted-data columns are missing.")

    d <- d[, needed, drop = FALSE]
    d$SPECIES <- as.character(d$SPECIES)
    d$SITE_ID <- as.character(d$SITE_ID)
    d$bms_id <- as.character(d$bms_id)
    d$YEAR <- as.integer(d$YEAR)

    need(
      !anyNA(d) &&
        all(is.finite(as.matrix(
          d[, c(spec$anomaly, spec$moderators), drop = FALSE]
        )),
      "Missing/non-finite fitted rows."
    )

    need(
      !anyDuplicated(d[c("SPECIES", "SITE_ID", "YEAR")]),
      "Repeated species-site-year rows found."
    )

    context <- d |>
      dplyr::group_by(.data$SPECIES, .data$SITE_ID, .data$bms_id) |>
      dplyr::summarise(
        n_years = dplyr::n_distinct(.data$YEAR),
        first_year = min(.data$YEAR),
        last_year = max(.data$YEAR),
        anomaly_SD_within_population = stats::sd(.data[[spec$anomaly]]),
        dplyr::across(dplyr::all_of(spec$moderators), mean),
        .groups = "drop"
      ) |>
      as.data.frame()

    context$pop_id <- paste(context$SPECIES, context$SITE_ID, sep = "_")
    need(!anyDuplicated(context$pop_id),
         "Population IDs are not unique.")

    sp_ix <- match(context$SPECIES, coef$species_slopes$SPECIES)
    need(!anyNA(sp_ix),
         "Could not match species random slopes.")

    context$species_slope_deviation <-
      coef$species_slopes$species_slope_deviation[sp_ix]
    context$species_slope_marginal_SE <-
      coef$species_slopes$species_slope_marginal_SE[sp_ix]

    need(spec$anomaly %in% names(coef$beta),
         "Fixed anomaly slope not found.")

    context$fixed_anomaly_slope <- unname(coef$beta[spec$anomaly])
    context$interaction_contribution <- 0

    for (m in spec$moderators) {
      candidates <- c(
        paste(spec$anomaly, m, sep = ":"),
        paste(m, spec$anomaly, sep = ":")
      )
      term <- intersect(candidates, names(coef$beta))
      need(length(term) == 1L,
           paste("Missing/ambiguous interaction:", m))

      contribution <- as.numeric(coef$beta[term]) * context[[m]]
      context[[paste0("contribution__", m)]] <- contribution
      context$interaction_contribution <-
        context$interaction_contribution + contribution
    }

    context$plasticity_raw_days_per_C <-
      context$fixed_anomaly_slope +
      context$species_slope_deviation +
      context$interaction_contribution

    # Biological orientation used by the abundance workflow:
    # higher = stronger advance for uni; stronger delay for multi.
    context$plasticity_bio_days_per_C <-
      if (spec$group == "univoltine") {
        -context$plasticity_raw_days_per_C
      } else {
        context$plasticity_raw_days_per_C
      }

    context$voltinism <- spec$group
    context$analysis <- spec$analysis

    need(all(is.finite(context$plasticity_raw_days_per_C)),
         "Non-finite population plasticity.")

    context[order(context$SPECIES, context$SITE_ID), , drop = FALSE]
  }

  read_source <- function(run, analysis, version = c("old", "qc")) {
    version <- match.arg(version)
    d <- file.path(run, analysis)

    if (version == "old") {
      fit_file <- file.path(d, "fit_with_onset.rds")
      input_file <- file.path(d, "input.rds")
      done_file <- file.path(d, "completed.rds")

      need(all(file.exists(c(fit_file, input_file, done_file))),
           paste("Missing OLD source files for", analysis))

      saved <- readRDS(fit_file)
      input <- readRDS(input_file)
      done <- readRDS(done_file)

      need(identical(done$status, "FIT_OK"),
           paste("OLD fit is not FIT_OK:", analysis))
      need(identical(saved$signature, input$signature) &&
             identical(done$signature, input$signature),
           paste("OLD source signature mismatch:", analysis))

      fit <- saved$fit

    } else {
      fit_file <- file.path(d, "fit_qc.rds")
      input_file <- file.path(d, "input.rds")
      done_file <- file.path(d, "completed.rds")

      need(all(file.exists(c(fit_file, input_file, done_file))),
           paste("Missing QC source files for", analysis))

      saved <- readRDS(fit_file)
      input <- readRDS(input_file)
      done <- readRDS(done_file)

      need(identical(done$status, "FIT_OK"),
           paste("QC fit is not FIT_OK:", analysis))
      need(identical(saved$signature, input$signature) &&
             identical(done$signature, input$signature),
           paste("QC source signature mismatch:", analysis))

      fit <- saved$fit
    }

    list(
      fit = fit,
      input = input,
      fit_file = normalizePath(fit_file, winslash = "/", mustWork = TRUE),
      input_file = normalizePath(input_file, winslash = "/", mustWork = TRUE)
    )
  }

  extract_one <- function(run, analysis, version) {
    spec <- specification(analysis)
    message("Reading ", toupper(version), ": ", analysis)

    src <- read_source(run, analysis, version)
    fit <- src$fit

    numerical_check(fit, paste(version, analysis))
    validate_formula(fit, spec, paste(version, analysis))

    # Verify the fit data match the saved prepared input.
    input_data <- as.data.frame(src$input$data)
    fit_data <- as.data.frame(fit$data)

    need(
      nrow(input_data) == nrow(fit_data) &&
        identical(names(input_data), names(fit_data)),
      paste(version, analysis, "input/fit data dimensions differ.")
    )

    # source_row is the strongest lightweight identity check if present.
    if ("source_row" %in% names(input_data) && "source_row" %in% names(fit_data)) {
      need(
        identical(as.integer(input_data$source_row),
                  as.integer(fit_data$source_row)),
        paste(version, analysis, "source rows differ between fit/input.")
      )
    }

    coef <- extract_coefficients(fit, spec)
    pop <- population_contextual_slopes(fit, coef, spec)

    fixed <- data.frame(
      term = names(coef$beta),
      estimate = as.numeric(coef$beta),
      std.error = sqrt(diag(coef$V)),
      stringsAsFactors = FALSE
    )

    list(
      populations = pop,
      fixed = fixed,
      species_slopes = coef$species_slopes,
      n_obs = nrow(fit$data),
      source_fit = src$fit_file
    )
  }

  compare_group <- function(group) {
    analysis <- paste0("offset_", group)

    old <- extract_one(old_run, analysis, "old")
    # Release the large old fit before loading QC.
    invisible(gc())

    qc <- extract_one(qc_run, analysis, "qc")
    invisible(gc())

    oldp <- old$populations
    qcp <- qc$populations

    write_csv(
      oldp,
      file.path(out_dir, paste0("old_population_plasticity_", group, ".csv"))
    )
    write_csv(
      qcp,
      file.path(out_dir, paste0("qc_population_plasticity_", group, ".csv"))
    )

    old_ids <- oldp$pop_id
    qc_ids <- qcp$pop_id
    common <- intersect(old_ids, qc_ids)

    need(length(common) > 1L,
         paste("Too few common populations in", group))

    unmatched <- rbind(
      data.frame(
        group = group,
        pop_id = setdiff(old_ids, qc_ids),
        status = "present_old_missing_QC",
        stringsAsFactors = FALSE
      ),
      data.frame(
        group = group,
        pop_id = setdiff(qc_ids, old_ids),
        status = "present_QC_missing_old",
        stringsAsFactors = FALSE
      )
    )

    io <- match(common, oldp$pop_id)
    iq <- match(common, qcp$pop_id)

    ctab <- data.frame(
      group = group,
      pop_id = common,
      SPECIES = oldp$SPECIES[io],
      SITE_ID = oldp$SITE_ID[io],
      bms_id = oldp$bms_id[io],

      old_n_years = oldp$n_years[io],
      qc_n_years = qcp$n_years[iq],
      years_removed_by_QC = oldp$n_years[io] - qcp$n_years[iq],

      old_raw_days_per_C = oldp$plasticity_raw_days_per_C[io],
      qc_raw_days_per_C = qcp$plasticity_raw_days_per_C[iq],

      old_bio_days_per_C = oldp$plasticity_bio_days_per_C[io],
      qc_bio_days_per_C = qcp$plasticity_bio_days_per_C[iq],

      old_species_slope_dev = oldp$species_slope_deviation[io],
      qc_species_slope_dev = qcp$species_slope_deviation[iq],

      old_interaction_contribution = oldp$interaction_contribution[io],
      qc_interaction_contribution = qcp$interaction_contribution[iq],

      stringsAsFactors = FALSE
    )

    need(
      identical(as.character(oldp$SPECIES[io]), as.character(qcp$SPECIES[iq])) &&
        identical(as.character(oldp$SITE_ID[io]), as.character(qcp$SITE_ID[iq])),
      paste("Population identities differ within", group)
    )

    # Standardize each version over exactly the same matched population set.
    ctab$old_bio_z <- as.numeric(scale(ctab$old_bio_days_per_C))
    ctab$qc_bio_z <- as.numeric(scale(ctab$qc_bio_days_per_C))

    ctab$raw_change_days_per_C <-
      ctab$qc_raw_days_per_C - ctab$old_raw_days_per_C
    ctab$abs_raw_change_days_per_C <- abs(ctab$raw_change_days_per_C)

    ctab$bio_change_days_per_C <-
      ctab$qc_bio_days_per_C - ctab$old_bio_days_per_C
    ctab$abs_bio_change_days_per_C <- abs(ctab$bio_change_days_per_C)

    ctab$z_change <- ctab$qc_bio_z - ctab$old_bio_z
    ctab$abs_z_change <- abs(ctab$z_change)

    # Rank by standardized change: easiest way to detect populations that moved.
    ctab <- ctab[
      order(-ctab$abs_z_change, -ctab$abs_raw_change_days_per_C),
      ,
      drop = FALSE
    ]

    summary <- data.frame(
      group = group,
      old_n_observations = old$n_obs,
      qc_n_observations = qc$n_obs,
      observations_removed = old$n_obs - qc$n_obs,
      percent_observations_removed =
        100 * (old$n_obs - qc$n_obs) / old$n_obs,

      old_n_populations = nrow(oldp),
      qc_n_populations = nrow(qcp),
      matched_n_populations = length(common),
      populations_lost_after_QC = length(setdiff(old_ids, qc_ids)),
      populations_new_after_QC = length(setdiff(qc_ids, old_ids)),

      cor_raw_old_vs_QC =
        stats::cor(ctab$old_raw_days_per_C, ctab$qc_raw_days_per_C),
      cor_bio_z_old_vs_QC =
        stats::cor(ctab$old_bio_z, ctab$qc_bio_z),

      median_abs_raw_change_days_per_C =
        stats::median(ctab$abs_raw_change_days_per_C),
      p95_abs_raw_change_days_per_C =
        as.numeric(stats::quantile(
          ctab$abs_raw_change_days_per_C, 0.95, names = FALSE
        )),
      max_abs_raw_change_days_per_C =
        max(ctab$abs_raw_change_days_per_C),

      median_abs_z_change =
        stats::median(ctab$abs_z_change),
      p95_abs_z_change =
        as.numeric(stats::quantile(
          ctab$abs_z_change, 0.95, names = FALSE
        )),
      max_abs_z_change =
        max(ctab$abs_z_change),

      n_abs_z_change_gt_0_25 = sum(ctab$abs_z_change > 0.25),
      n_abs_z_change_gt_0_50 = sum(ctab$abs_z_change > 0.50),
      n_abs_z_change_gt_1_00 = sum(ctab$abs_z_change > 1.00),

      max_years_removed_single_population =
        max(ctab$years_removed_by_QC),

      stringsAsFactors = FALSE
    )

    top <- utils::head(ctab, top_n)

    list(
      comparison = ctab,
      summary = summary,
      top = top,
      unmatched = unmatched
    )
  }

  message("==============================================================")
  message("OFFSET PLASTICITY QC COMPARISON")
  message("No models will be fitted.")
  message("==============================================================")

  ans <- lapply(groups, compare_group)
  names(ans) <- groups

  summary <- do.call(rbind, lapply(ans, `[[`, "summary"))
  comparison <- do.call(rbind, lapply(ans, `[[`, "comparison"))
  top <- do.call(rbind, lapply(ans, `[[`, "top"))
  unmatched <- do.call(rbind, lapply(ans, `[[`, "unmatched"))

  write_csv(summary, file.path(out_dir, "comparison_summary.csv"))
  write_csv(comparison, file.path(out_dir, "population_comparison_all.csv"))
  write_csv(top, file.path(out_dir, "top_changed_populations.csv"))
  write_csv(unmatched, file.path(out_dir, "unmatched_populations.csv"))

  # Compact human-readable report.
  report <- c(
    "phenoIMPACT | Offset plasticity OLD vs QC",
    "",
    paste("Created:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste("OLD source:", normalizePath(old_run, winslash = "/", mustWork = TRUE)),
    paste("QC source:", normalizePath(qc_run, winslash = "/", mustWork = TRUE)),
    "",
    "Plasticity is re-extracted independently from both saved fits.",
    "No abundance model is loaded or fitted.",
    "Z scores are standardized separately OLD/QC over the SAME matched populations.",
    "",
    capture.output(print(summary, row.names = FALSE))
  )
  writeLines(report, file.path(out_dir, "README_results.txt"))

  message("\nDONE. No model fitted.")
  message("Output: ", normalizePath(out_dir, winslash = "/", mustWork = TRUE))
  message("\nSUMMARY:")
  print(summary, row.names = FALSE)

  invisible(list(
    summary = summary,
    comparison = comparison,
    top = top,
    unmatched = unmatched,
    output = out_dir
  ))
}

offset_qc_plasticity_comparison <- compare_offset_qc_plasticity()
