# =============================================================================
# phenoIMPACT | FINAL abundance-input preparation after OFFSET LAST_PEAK QC
# Version 1.0 | 2026-10-07
#
# NO MODEL IS FITTED.
#
# This script:
#   1) quantifies how many abundance population-year observations are removed
#      because their populations lost all valid OFFSET estimates after QC;
#   2) uses the FINAL QC offset plasticities;
#   3) standardizes OFFSET plasticity over ALL QC-extracted populations within
#      each voltinism group, BEFORE the abundance join (same principle as before);
#   4) preserves ONSET plasticity z EXACTLY for retained abundance populations;
#   5) drops only abundance populations absent from the FINAL QC offset extraction;
#   6) rebuilds population AR1 previous/gap links after filtering;
#   7) rebuilds species/site indices and subsets their warm-start random effects;
#   8) preserves the same calendar years, 15-km mesh, response, year scale,
#      Gamma model design, and spatial-field warm start;
#   9) writes definitive input_univoltine.rds / input_multivoltine.rds.
#
# It does NOT compile TMB, optimize, estimate uncertainty, or fit any model.
# =============================================================================

prepare_final_qc_abundance_inputs <- function(
  root = "E:/phenoIMPACT project/code/phenoIMPACT",

  plasticity_dir = file.path(
    root, "output", "phenology_plasticity",
    "offset_with_onset_qc_lastpeak_15km", "run_fbf8b3101079",
    "plasticity_comparison_pre_abundance"
  ),

  current_abundance_run = file.path(
    root, "output", "population_trends",
    "spatial_abundance_onset_adjusted_15km",
    "run_20261004_214012_44136"
  ),

  refined_univoltine_fit = file.path(
    root, "output", "population_trends",
    "spatial_abundance_onset_adjusted_15km",
    "run_20261004_214012_44136",
    "refinement_univoltine", "run_20261005_073241_44136",
    "univoltine__full", "fit.rds"
  ),

  out_root = file.path(
    root, "output", "population_trends",
    "spatial_abundance_qc_offset_final_15km"
  )
) {

  need <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)
  same_num <- function(a, b, tol = 1e-8) {
    length(a) == length(b) &&
      all(is.finite(a)) && all(is.finite(b)) &&
      (length(a) == 0L || max(abs(as.numeric(a) - as.numeric(b))) <= tol)
  }
  write_csv <- function(x, path) {
    utils::write.csv(x, path, row.names = FALSE, na = "")
  }
  atomic_rds <- function(x, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    tmp <- tempfile(".writing_", tmpdir = dirname(path), fileext = ".rds")
    on.exit(unlink(tmp), add = TRUE)
    saveRDS(x, tmp, compress = "gzip")
    need(file.rename(tmp, path), paste("Could not finalize:", path))
    invisible(path)
  }

  groups <- c("univoltine", "multivoltine")

  need(requireNamespace("digest", quietly = TRUE),
       "Package 'digest' is required.")

  # Explicit current warm-start sources.
  current_fit_file <- c(
    univoltine = refined_univoltine_fit,
    multivoltine = file.path(
      current_abundance_run, "multivoltine__full", "fit.rds"
    )
  )

  required <- c(
    file.path(current_abundance_run, paste0("input_", groups, ".rds")),
    file.path(current_abundance_run, "mesh.rds"),
    file.path(current_abundance_run, "compiled", "trend_spde.cpp"),
    current_fit_file,
    file.path(plasticity_dir, paste0("qc_population_plasticity_", groups, ".csv"))
  )

  absent <- required[!file.exists(required)]
  need(!length(absent),
       paste("Missing required source file(s):\n", paste(absent, collapse = "\n")))

  # Fingerprint only the relevant preparation inputs.
  manifest <- data.frame(
    path = normalizePath(required, winslash = "/", mustWork = TRUE),
    md5 = unname(tools::md5sum(required)),
    stringsAsFactors = FALSE
  )
  need(!anyNA(manifest$md5), "Could not fingerprint one or more source files.")

  signature <- digest::digest(
    list(
      version = "final_qc_abundance_preparation_1.0",
      manifest = manifest,
      scaling_rule = "offset QC extraction universe within voltinism",
      onset_rule = "retain current onset_plasticity_z exactly",
      population_rule = "drop abundance populations missing from QC offset extraction"
    ),
    algo = "sha256"
  )

  out_dir <- file.path(out_root, paste0("prep_", substr(signature, 1, 12)))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write_csv(manifest, file.path(out_dir, "source_manifest.csv"))

  # Copy read-only structural inputs for the eventual fit workflow.
  need(file.copy(
    file.path(current_abundance_run, "mesh.rds"),
    file.path(out_dir, "mesh.rds"),
    overwrite = TRUE
  ), "Could not copy mesh.rds.")

  dir.create(file.path(out_dir, "compiled"), showWarnings = FALSE)
  need(file.copy(
    file.path(current_abundance_run, "compiled", "trend_spde.cpp"),
    file.path(out_dir, "compiled", "trend_spde.cpp"),
    overwrite = TRUE
  ), "Could not copy trend_spde.cpp.")

  all_summary <- list()
  all_removed <- list()
  all_scaling <- list()
  all_population_plasticity <- list()
  all_checks <- list()

  for (g in groups) {

    message("--------------------------------------------------------------")
    message("Preparing final abundance input: ", g)

    input_file <- file.path(current_abundance_run, paste0("input_", g, ".rds"))
    old <- readRDS(input_file)
    fit <- readRDS(current_fit_file[[g]])
    qc <- utils::read.csv(
      file.path(plasticity_dir, paste0("qc_population_plasticity_", g, ".csv")),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    need(identical(old$group, g), paste("Wrong current abundance input group:", g))
    need(identical(fit$group, g) && identical(fit$variant, "full"),
         paste("Wrong current abundance fit:", g))
    need(isTRUE(fit$checks$numerical_checks_ok),
         paste("Current abundance fit is not numerically accepted:", g))
    need(identical(fit$data_signature, old$data_signature),
         paste("Current fit/input data signature mismatch:", g))
    need(
      identical(
        old$data_signature,
        digest::digest(list(old$frame, old$data$X), algo = "sha256")
      ),
      paste("Current abundance input signature cannot be reproduced:", g)
    )

    f <- old$frame
    need(all(c(
      "SPECIES", "SITE_ID", "pop_id", "YEAR", "year_decade",
      "ABUND_INDEX", "onset_plasticity_z", "offset_plasticity_bio_z",
      "bms_id", "x_km", "y_km"
    ) %in% names(f)), paste("Current frame lacks required columns:", g))

    need(nrow(f) == length(old$data$y),
         paste("Current frame/data row mismatch:", g))
    need(length(fit$parameters$ar_obs) == nrow(f),
         paste("Current AR1 warm start does not align with rows:", g))

    need(all(c(
      "pop_id", "SPECIES", "SITE_ID",
      "plasticity_raw_days_per_C", "plasticity_bio_days_per_C"
    ) %in% names(qc)), paste("QC population table lacks required columns:", g))

    qc$pop_id <- as.character(qc$pop_id)
    qc$SPECIES <- as.character(qc$SPECIES)
    qc$SITE_ID <- as.character(qc$SITE_ID)

    need(!anyDuplicated(qc$pop_id), paste("Duplicate QC population IDs:", g))
    need(all(is.finite(qc$plasticity_bio_days_per_C)),
         paste("Non-finite QC plasticity:", g))

    # FINAL scaling: ALL populations surviving QC extraction, before abundance join.
    offset_mu <- mean(qc$plasticity_bio_days_per_C)
    offset_sd <- stats::sd(qc$plasticity_bio_days_per_C)
    need(is.finite(offset_sd) && offset_sd > 0,
         paste("Invalid QC offset SD:", g))

    qc$offset_plasticity_bio_z_final <-
      (qc$plasticity_bio_days_per_C - offset_mu) / offset_sd

    all_scaling[[g]] <- data.frame(
      group = g,
      n_QC_extracted_populations = nrow(qc),
      offset_bio_mean_days_per_C = offset_mu,
      offset_bio_SD_days_per_C = offset_sd,
      scaling_universe =
        "all QC-extracted offset populations within group; before abundance join",
      stringsAsFactors = FALSE
    )

    old_pop <- unique(as.character(f$pop_id))
    qc_pop <- qc$pop_id
    lost <- setdiff(old_pop, qc_pop)

    removed <- f[f$pop_id %in% lost, , drop = FALSE]
    retained_old_rows <- which(!f$pop_id %in% lost)

    need(length(retained_old_rows) > 0L, paste("No retained abundance rows:", g))
    need(nrow(removed) > 0L, paste(
      "No abundance rows were removed in", g,
      "- this conflicts with the prior overlap check."
    ))

    # Population-level removed detail.
    rem_detail <- if (nrow(removed)) {
      x <- aggregate(
        YEAR ~ pop_id + SPECIES + SITE_ID + bms_id,
        data = removed,
        FUN = function(z) length(unique(z))
      )
      names(x)[names(x) == "YEAR"] <- "n_abundance_years_removed"
      yr_first <- aggregate(
        YEAR ~ pop_id, data = removed, FUN = min
      )
      names(yr_first)[2] <- "first_abundance_year_removed"
      yr_last <- aggregate(
        YEAR ~ pop_id, data = removed, FUN = max
      )
      names(yr_last)[2] <- "last_abundance_year_removed"
      x <- merge(x, yr_first, by = "pop_id", all.x = TRUE, sort = FALSE)
      x <- merge(x, yr_last, by = "pop_id", all.x = TRUE, sort = FALSE)
      x$group <- g
      x
    } else {
      data.frame()
    }
    all_removed[[g]] <- rem_detail

    # Sort retained rows exactly as required for population AR1 reconstruction.
    d <- f[retained_old_rows, , drop = FALSE]
    source_row <- retained_old_rows
    ord <- order(as.character(d$pop_id), as.integer(d$YEAR))
    d <- d[ord, , drop = FALSE]
    source_row <- source_row[ord]

    need(
      identical(
        order(as.character(d$pop_id), as.integer(d$YEAR)),
        seq_len(nrow(d))
      ),
      paste("Could not establish population/year ordering:", g)
    )
    need(!anyDuplicated(d[c("pop_id", "YEAR")]),
         paste("Duplicate retained population-years:", g))

    # Join FINAL QC offset plasticity.
    ix <- match(as.character(d$pop_id), qc$pop_id)
    need(!anyNA(ix), paste("A retained abundance population lacks QC plasticity:", g))
    need(
      identical(as.character(d$SPECIES), qc$SPECIES[ix]) &&
        identical(as.character(d$SITE_ID), qc$SITE_ID[ix]),
      paste("QC plasticity identity mismatch:", g)
    )

    onset_before <- d$onset_plasticity_z
    d$offset_plasticity_bio_z <- qc$offset_plasticity_bio_z_final[ix]

    need(same_num(d$onset_plasticity_z, onset_before, 0),
         paste("ONSET z changed unexpectedly:", g))
    need(all(is.finite(d$offset_plasticity_bio_z)),
         paste("Non-finite final offset z:", g))

    # Build the same fixed-effect design.
    X <- stats::model.matrix(
      ~ year_decade * onset_plasticity_z +
        year_decade * offset_plasticity_bio_z,
      d
    )
    need(setequal(colnames(X), colnames(old$data$X)),
         paste("Fixed-effect terms changed:", g))
    X <- X[, colnames(old$data$X), drop = FALSE]
    need(qr(X)$rank == ncol(X), paste("Rank-deficient final design:", g))

    # Recover OLD species/site parameter ordering from row-level indices.
    sp_map <- unique(data.frame(
      SPECIES = as.character(f$SPECIES),
      old_index = as.integer(old$data$species),
      stringsAsFactors = FALSE
    ))
    site_map <- unique(data.frame(
      SITE_ID = as.character(f$SITE_ID),
      old_index = as.integer(old$data$site),
      stringsAsFactors = FALSE
    ))

    need(!anyDuplicated(sp_map$SPECIES) && !anyDuplicated(sp_map$old_index),
         paste("Invalid old species index map:", g))
    need(!anyDuplicated(site_map$SITE_ID) && !anyDuplicated(site_map$old_index),
         paste("Invalid old site index map:", g))

    used_species <- unique(as.character(d$SPECIES))
    used_sites <- unique(as.character(d$SITE_ID))

    sp_keep <- sp_map[sp_map$SPECIES %in% used_species, , drop = FALSE]
    sp_keep <- sp_keep[order(sp_keep$old_index), , drop = FALSE]

    site_keep <- site_map[site_map$SITE_ID %in% used_sites, , drop = FALSE]
    site_keep <- site_keep[order(site_keep$old_index), , drop = FALSE]

    need(setequal(sp_keep$SPECIES, used_species),
         paste("Species remapping failed:", g))
    need(setequal(site_keep$SITE_ID, used_sites),
         paste("Site remapping failed:", g))

    species_new <- match(as.character(d$SPECIES), sp_keep$SPECIES) - 1L
    site_new <- match(as.character(d$SITE_ID), site_keep$SITE_ID) - 1L
    need(!anyNA(species_new) && !anyNA(site_new),
         paste("New random-effect indices contain NA:", g))

    # Preserve the existing calendar-year universe.
    years <- as.integer(old$years)
    missing_entire_year <- setdiff(years, unique(as.integer(d$YEAR)))
    need(!length(missing_entire_year), paste(
      "QC filtering removed all abundance observations for entire calendar year(s) in",
      g, ":", paste(missing_entire_year, collapse = ", "),
      "- stop for review."
    ))
    year_new <- match(as.integer(d$YEAR), years) - 1L
    need(!anyNA(year_new), paste("Year remapping failed:", g))

    # Rebuild stationary population AR1 links using true calendar gaps.
    n <- nrow(d)
    previous <- rep(-1L, n)
    gap <- rep(1L, n)

    if (n > 1L) {
      same_pop <- as.character(d$pop_id[-1]) ==
        as.character(d$pop_id[-n])
      ids <- which(same_pop) + 1L
      previous[ids] <- ids - 2L  # zero-based index of immediately previous row
      gap[ids] <- as.integer(d$YEAR[ids] - d$YEAR[ids - 1L])
    }
    need(all(gap > 0L), paste("Non-positive AR1 gap:", g))

    # Row-aligned AR1 warm start comes from the reviewed current fit.
    ar_new <- as.numeric(fit$parameters$ar_obs[source_row])

    # Subset warm-start random effects to retained species/sites.
    pars <- fit$parameters
    need(nrow(pars$b_species) == nrow(sp_map),
         paste("Species warm-start dimension mismatch:", g))
    need(nrow(pars$b_site) == nrow(site_map),
         paste("Site warm-start dimension mismatch:", g))

    pars$b_species <- pars$b_species[sp_keep$old_index + 1L, , drop = FALSE]
    pars$b_site <- pars$b_site[site_keep$old_index + 1L, , drop = FALSE]
    pars$ar_obs <- ar_new

    # beta design is structurally identical; use reviewed current beta as warm start.
    need(length(pars$beta) == ncol(X),
         paste("Beta warm-start dimension mismatch:", g))
    pars$beta <- as.numeric(fit$beta[colnames(X)])
    need(all(is.finite(pars$beta)),
         paste("Could not align beta warm start:", g))

    # Annual spatial field is defined on mesh vertices x years and remains unchanged.
    need(
      is.matrix(pars$field) &&
        ncol(pars$field) == length(years),
      paste("Spatial-field warm-start year dimension mismatch:", g)
    )

    # Subset site coordinate table in the SAME retained site order.
    site_ix <- match(site_keep$SITE_ID, as.character(old$sites$SITE_ID))
    need(!anyNA(site_ix), paste("Retained site coordinates missing:", g))
    sites_new <- old$sites[site_ix, , drop = FALSE]
    need(identical(as.character(sites_new$SITE_ID), site_keep$SITE_ID),
         paste("Site coordinate order mismatch:", g))

    newdata <- list(
      y = as.numeric(d$ABUND_INDEX),
      X = X,
      time = as.numeric(d$year_decade),
      species = as.integer(species_new),
      site = as.integer(site_new),
      year = as.integer(year_new),
      previous = as.integer(previous),
      gap = as.integer(gap)
    )

    # Direct invariants against retained current rows.
    oldXsub <- old$data$X[source_row, , drop = FALSE]
    unchanged_cols <- setdiff(
      colnames(X),
      c("offset_plasticity_bio_z",
        "year_decade:offset_plasticity_bio_z")
    )
    need(
      same_num(X[, unchanged_cols, drop = FALSE],
               oldXsub[, unchanged_cols, drop = FALSE], 1e-10),
      paste("A non-offset fixed-effect column changed:", g)
    )
    need(same_num(newdata$y, old$data$y[source_row], 0),
         paste("Abundance response changed:", g))
    need(same_num(newdata$time, old$data$time[source_row], 0),
         paste("Year scale changed:", g))

    input_final <- list(
      signature = signature,
      schema = "final_qc_abundance_preparation_1.0",
      group = g,
      data = newdata,
      frame = d,
      sites = sites_new,
      years = years,
      par = pars,
      data_signature = digest::digest(list(d, X), algo = "sha256"),
      source_data_signature = old$data_signature,
      source_abundance_run = normalizePath(
        current_abundance_run, winslash = "/", mustWork = TRUE
      ),
      source_fit = normalizePath(
        current_fit_file[[g]], winslash = "/", mustWork = TRUE
      ),
      offset_QC = TRUE,
      offset_adjusted_for_onset = TRUE,
      offset_scaling_universe = "all QC-extracted populations before abundance join",
      onset_predictor_preserved_exactly = TRUE,
      removed_populations = lost,
      retained_source_rows = source_row,
      diagnostics_status = "PENDING"
    )

    atomic_rds(
      input_final,
      file.path(out_dir, paste0("input_", g, ".rds"))
    )

    # Full final QC plasticity table for provenance.
    ptab <- qc[, c(
      "pop_id", "SPECIES", "SITE_ID",
      "plasticity_raw_days_per_C", "plasticity_bio_days_per_C"
    )]
    ptab$group <- g
    ptab$offset_plasticity_bio_z_final <-
      qc$offset_plasticity_bio_z_final
    all_population_plasticity[[g]] <- ptab

    species_removed_entirely <- setdiff(sp_map$SPECIES, sp_keep$SPECIES)
    sites_removed_entirely <- setdiff(site_map$SITE_ID, site_keep$SITE_ID)

    all_summary[[g]] <- data.frame(
      group = g,
      old_abundance_observations = nrow(f),
      final_abundance_observations = nrow(d),
      observations_removed = nrow(f) - nrow(d),
      percent_observations_removed =
        100 * (nrow(f) - nrow(d)) / nrow(f),

      old_abundance_populations = length(unique(f$pop_id)),
      final_abundance_populations = length(unique(d$pop_id)),
      populations_removed = length(lost),
      percent_populations_removed =
        100 * length(lost) / length(unique(f$pop_id)),

      species_removed_entirely = length(species_removed_entirely),
      sites_removed_entirely = length(sites_removed_entirely),

      max_abundance_years_removed_single_population =
        if (nrow(rem_detail)) max(rem_detail$n_abundance_years_removed) else 0,

      mean_abundance_years_removed_per_population =
        if (nrow(rem_detail)) mean(rem_detail$n_abundance_years_removed) else 0,

      median_abundance_years_removed_per_population =
        if (nrow(rem_detail)) stats::median(rem_detail$n_abundance_years_removed) else 0,

      n_QC_offset_scaling_populations = nrow(qc),
      offset_QC_mean_days_per_C = offset_mu,
      offset_QC_SD_days_per_C = offset_sd,

      stringsAsFactors = FALSE
    )

    all_checks[[g]] <- data.frame(
      group = g,
      response_identical_for_retained_rows = TRUE,
      year_scale_identical_for_retained_rows = TRUE,
      onset_z_identical_for_retained_rows = TRUE,
      only_offset_design_columns_changed = TRUE,
      AR1_links_rebuilt_after_filter = TRUE,
      true_calendar_gaps_used = TRUE,
      same_calendar_year_universe = TRUE,
      same_15km_mesh_file = TRUE,
      current_reviewed_fit_used_as_warm_start = TRUE,
      final_input_signature = input_final$data_signature,
      stringsAsFactors = FALSE
    )

    message(
      "  removed ", nrow(f) - nrow(d), " / ", nrow(f),
      " abundance rows (",
      round(100 * (nrow(f) - nrow(d)) / nrow(f), 4), "%)"
    )
    message(
      "  removed ", length(lost), " / ", length(unique(f$pop_id)),
      " populations (",
      round(100 * length(lost) / length(unique(f$pop_id)), 4), "%)"
    )
    message(
      "  whole species lost: ", length(species_removed_entirely),
      " | whole sites lost: ", length(sites_removed_entirely)
    )
  }

  summary <- do.call(rbind, all_summary)
  removed <- do.call(rbind, all_removed)
  scaling <- do.call(rbind, all_scaling)
  plasticity <- do.call(rbind, all_population_plasticity)
  checks <- do.call(rbind, all_checks)

  write_csv(summary, file.path(out_dir, "final_abundance_filter_summary.csv"))
  write_csv(removed, file.path(out_dir, "removed_abundance_populations_detail.csv"))
  write_csv(scaling, file.path(out_dir, "final_offset_plasticity_scaling.csv"))
  write_csv(plasticity, file.path(out_dir, "final_qc_offset_population_plasticity.csv"))
  write_csv(checks, file.path(out_dir, "final_input_checks.csv"))

  writeLines(c(
    "phenoIMPACT | FINAL abundance preparation after offset LAST_PEAK QC",
    "",
    paste("Signature:", signature),
    "NO MODEL FITTED.",
    "",
    "Population rule:",
    "  Abundance populations absent from the FINAL QC offset extraction are removed.",
    "",
    "OFFSET scaling:",
    "  Standardized within voltinism using ALL populations surviving the QC offset extraction,",
    "  before matching to abundance. This follows the previous scaling principle.",
    "",
    "ONSET:",
    "  Existing onset_plasticity_z is preserved bit-for-bit for every retained abundance row.",
    "",
    "Abundance model structure prepared:",
    "  Same response, year_decade, fixed formula, population AR1 concept, calendar-year universe,",
    "  and 15-km mesh. AR1 previous/gap links are rebuilt after filtering.",
    "",
    "Warm starts:",
    "  Multivoltine = reviewed current FIT_OK abundance fit.",
    "  Univoltine = reviewed refined FIT_OK abundance fit.",
    "  Species/site random-effect arrays are subset/reindexed if required.",
    "  Annual spatial field remains on the same mesh vertices and years.",
    "",
    "Files input_univoltine.rds and input_multivoltine.rds are definitive PREPARED inputs.",
    "They have NOT yet been passed to TMB or optimized.",
    "",
    capture.output(print(summary, row.names = FALSE))
  ), file.path(out_dir, "README.txt"))

  cat("\n==============================================================\n")
  cat("FINAL QC ABUNDANCE INPUT PREPARATION\n")
  cat("NO MODEL FITTED\n")
  cat("==============================================================\n\n")
  print(summary, row.names = FALSE)
  cat("\nOutput: ", normalizePath(out_dir, winslash = "/", mustWork = TRUE), "\n", sep = "")
  cat("\nPrepared:\n")
  cat("  input_univoltine.rds\n")
  cat("  input_multivoltine.rds\n")
  cat("  mesh.rds\n")
  cat("  compiled/trend_spde.cpp\n")

  invisible(list(
    summary = summary,
    removed = removed,
    scaling = scaling,
    checks = checks,
    output = out_dir
  ))
}

final_qc_abundance_preparation <- prepare_final_qc_abundance_inputs()
