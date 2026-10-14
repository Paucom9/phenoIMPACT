# =============================================================================
# phenoIMPACT | Audit of geographic variation in voltinism classification
# Version 1.0 | 2026-10-07
#
# Fast diagnostic: NO model fitting.
#
# It reproduces the voltinism join used in pheno_plasticity.R and asks:
#   1) Which species are always univoltine, always multivoltine, or change
#      classification among BMS/country units represented in the analysis?
#   2) How many populations/offset observations belong to each class?
#   3) For species that change class, are multivoltine populations generally
#      located farther south than univoltine populations?
#   4) Which species are safe for a strict stable-voltinism sensitivity?
#
# RUN
#   source(file.choose(), encoding = "UTF-8")
#   volt_audit <- run_voltinism_audit()
#
# SEND BACK
#   .../voltinism_audit/run_*/results_to_review.zip
# =============================================================================

run_voltinism_audit <- function(
  project_root = "E:/phenoIMPACT project/code/phenoIMPACT"
) {
  need <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)
  csv <- function(x, p) utils::write.csv(x, p, row.names = FALSE, na = "")

  req <- c("dplyr", "tidyr", "sf", "ggplot2")
  miss <- req[!vapply(req, requireNamespace, logical(1), quietly = TRUE)]
  need(!length(miss), paste("Missing packages:", paste(miss, collapse = ", ")))

  project_root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
  volt_file <- file.path(project_root, "data", "voltinism", "species_country_voltinism.csv")
  coord_file <- file.path(project_root, "data", "ebms_transect_coord.csv")
  pheno_file <- file.path(project_root, "output", "pheno_estimates_allspp.csv")
  need(file.exists(volt_file), paste("Missing:", volt_file))
  need(file.exists(coord_file), paste("Missing:", coord_file))
  need(file.exists(pheno_file), paste("Missing:", pheno_file))

  out_root <- file.path(project_root, "output", "phenology_plasticity", "voltinism_audit")
  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
  out <- file.path(out_root, paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  dir.create(out, recursive = TRUE, showWarnings = FALSE)

  message("Reading voltinism table, coordinates and phenology populations...")
  volt <- utils::read.csv(volt_file, sep = ";", dec = ".", stringsAsFactors = FALSE,
                          check.names = FALSE)
  coords <- utils::read.csv(coord_file, stringsAsFactors = FALSE, check.names = FALSE)
  pheno <- utils::read.csv(pheno_file, stringsAsFactors = FALSE, check.names = FALSE)

  need("SPECIES" %in% names(volt), "Voltinism table lacks SPECIES.")
  need(all(c("transect_id", "bms_id", "transect_lon", "transect_lat") %in% names(coords)),
       "Coordinate table has unexpected columns.")
  need(all(c("SPECIES", "SITE_ID") %in% names(pheno)),
       "Phenology table lacks SPECIES/SITE_ID.")

  # Exact long-format logic inherited from the original phenology workflow.
  vlong <- tidyr::pivot_longer(
    volt,
    cols = -SPECIES,
    names_to = "bms_id",
    values_to = "voltinism_raw"
  ) |>
    dplyr::mutate(
      SPECIES = as.character(SPECIES),
      bms_id = dplyr::recode(
        as.character(bms_id),
        "ES.CTBMS" = "ES-CTBMS",
        "ES.ZEBMS" = "ES-ZEBMS"
      ),
      voltinism_raw = trimws(tolower(as.character(voltinism_raw))),
      voltinism = dplyr::case_when(
        voltinism_raw == "univoltine" ~ "univoltine",
        voltinism_raw == "multivoltine" ~ "multivoltine",
        TRUE ~ NA_character_
      )
    )

  # Coordinates: source CRS is EPSG:3035 metres.
  coord_sf <- coords |>
    dplyr::filter(!is.na(transect_lon), !is.na(transect_lat)) |>
    dplyr::distinct(bms_id, transect_id, transect_lon, transect_lat) |>
    sf::st_as_sf(coords = c("transect_lon", "transect_lat"), crs = 3035) |>
    sf::st_transform(4326)

  xy <- sf::st_coordinates(coord_sf)
  coord_site <- coord_sf |>
    sf::st_drop_geometry() |>
    dplyr::transmute(
      SITE_ID = as.character(transect_id),
      bms_id = as.character(bms_id),
      latitude_deg = xy[, 2],
      longitude_deg = xy[, 1]
    ) |>
    dplyr::distinct()

  need(!anyDuplicated(coord_site$SITE_ID),
       "SITE_ID maps to more than one coordinate/BMS; review coordinate table.")

  # Population universe represented in phenology estimates.
  pop <- pheno |>
    dplyr::transmute(
      SPECIES = as.character(SPECIES),
      SITE_ID = as.character(SITE_ID),
      pheno_type = if ("pheno_type" %in% names(pheno)) as.character(pheno_type) else NA_character_
    ) |>
    dplyr::filter(!is.na(SPECIES), nzchar(SPECIES), !is.na(SITE_ID), nzchar(SITE_ID)) |>
    dplyr::distinct()

  pop <- pop |>
    dplyr::left_join(coord_site, by = "SITE_ID", relationship = "many-to-one") |>
    dplyr::left_join(
      vlong |>
        dplyr::select(SPECIES, bms_id, voltinism_raw, voltinism),
      by = c("SPECIES", "bms_id"),
      relationship = "many-to-one"
    )

  need(!anyNA(pop$bms_id), "Some phenology sites lack BMS identity/coordinates.")

  # Use OFFSET populations as the biologically relevant audit universe when possible.
  if (any(!is.na(pop$pheno_type)) && any(pop$pheno_type == "OFFSET_mean", na.rm = TRUE)) {
    offset_pop <- pop |>
      dplyr::filter(pheno_type == "OFFSET_mean") |>
      dplyr::distinct(SPECIES, SITE_ID, bms_id, latitude_deg, longitude_deg,
                      voltinism_raw, voltinism)
  } else {
    warning("pheno_type/OFFSET_mean unavailable; using all phenology populations.")
    offset_pop <- pop |>
      dplyr::distinct(SPECIES, SITE_ID, bms_id, latitude_deg, longitude_deg,
                      voltinism_raw, voltinism)
  }

  # Species classification among BMS units actually represented by its analysed populations.
  analysed_class <- offset_pop |>
    dplyr::group_by(SPECIES) |>
    dplyr::summarise(
      n_populations = dplyr::n_distinct(SITE_ID),
      n_bms = dplyr::n_distinct(bms_id),
      n_bms_with_known_voltinism = dplyr::n_distinct(bms_id[!is.na(voltinism)]),
      n_uni_bms = dplyr::n_distinct(bms_id[voltinism == "univoltine"]),
      n_multi_bms = dplyr::n_distinct(bms_id[voltinism == "multivoltine"]),
      n_uni_populations = dplyr::n_distinct(SITE_ID[voltinism == "univoltine"]),
      n_multi_populations = dplyr::n_distinct(SITE_ID[voltinism == "multivoltine"]),
      latitude_min = min(latitude_deg, na.rm = TRUE),
      latitude_max = max(latitude_deg, na.rm = TRUE),
      latitude_span = latitude_max - latitude_min,
      analysis_stability_class = dplyr::case_when(
        n_uni_bms > 0 & n_multi_bms > 0 ~ "variable",
        n_uni_bms > 0 & n_multi_bms == 0 ~ "stable_univoltine",
        n_multi_bms > 0 & n_uni_bms == 0 ~ "stable_multivoltine",
        TRUE ~ "unknown"
      ),
      .groups = "drop"
    )

  # Independent classification from the complete source lookup table (all networks in table).
  source_class <- vlong |>
    dplyr::group_by(SPECIES) |>
    dplyr::summarise(
      source_n_uni = sum(voltinism == "univoltine", na.rm = TRUE),
      source_n_multi = sum(voltinism == "multivoltine", na.rm = TRUE),
      source_stability_class = dplyr::case_when(
        source_n_uni > 0 & source_n_multi > 0 ~ "variable",
        source_n_uni > 0 & source_n_multi == 0 ~ "stable_univoltine",
        source_n_multi > 0 & source_n_uni == 0 ~ "stable_multivoltine",
        TRUE ~ "unknown"
      ),
      .groups = "drop"
    )

  species_class <- analysed_class |>
    dplyr::left_join(source_class, by = "SPECIES") |>
    dplyr::mutate(
      stable_across_analysed_networks = analysis_stability_class %in%
        c("stable_univoltine", "stable_multivoltine"),
      stable_across_full_lookup = source_stability_class %in%
        c("stable_univoltine", "stable_multivoltine")
    ) |>
    dplyr::arrange(analysis_stability_class, SPECIES)

  # Detailed BMS/site support for species that switch category in analysed data.
  variable_species <- species_class |>
    dplyr::filter(analysis_stability_class == "variable") |>
    dplyr::pull(SPECIES)

  variable_detail <- offset_pop |>
    dplyr::filter(SPECIES %in% variable_species, !is.na(voltinism)) |>
    dplyr::group_by(SPECIES, bms_id, voltinism) |>
    dplyr::summarise(
      n_populations = dplyr::n_distinct(SITE_ID),
      mean_latitude = mean(latitude_deg),
      min_latitude = min(latitude_deg),
      max_latitude = max(latitude_deg),
      mean_longitude = mean(longitude_deg),
      .groups = "drop"
    ) |>
    dplyr::arrange(SPECIES, mean_latitude)

  # One independent geographic contrast per variable species.
  variable_shift <- offset_pop |>
    dplyr::filter(SPECIES %in% variable_species, !is.na(voltinism)) |>
    dplyr::group_by(SPECIES) |>
    dplyr::summarise(
      n_uni_populations = sum(voltinism == "univoltine"),
      n_multi_populations = sum(voltinism == "multivoltine"),
      mean_latitude_uni = mean(latitude_deg[voltinism == "univoltine"]),
      mean_latitude_multi = mean(latitude_deg[voltinism == "multivoltine"]),
      multi_minus_uni_latitude_deg = mean_latitude_multi - mean_latitude_uni,
      multivoltine_is_south = multi_minus_uni_latitude_deg < 0,
      .groups = "drop"
    ) |>
    dplyr::arrange(multi_minus_uni_latitude_deg)

  # Overall summaries.
  class_summary <- species_class |>
    dplyr::group_by(analysis_stability_class) |>
    dplyr::summarise(
      n_species = dplyr::n(),
      n_populations = sum(n_populations),
      .groups = "drop"
    )

  pop_summary <- offset_pop |>
    dplyr::left_join(
      species_class |>
        dplyr::select(SPECIES, analysis_stability_class),
      by = "SPECIES"
    ) |>
    dplyr::group_by(analysis_stability_class, voltinism) |>
    dplyr::summarise(
      n_species = dplyr::n_distinct(SPECIES),
      n_populations = dplyr::n_distinct(paste(SPECIES, SITE_ID)),
      .groups = "drop"
    )

  shift_summary <- if (nrow(variable_shift)) {
    data.frame(
      n_variable_species = nrow(variable_shift),
      n_multi_south = sum(variable_shift$multivoltine_is_south),
      proportion_multi_south = mean(variable_shift$multivoltine_is_south),
      median_multi_minus_uni_latitude_deg = stats::median(variable_shift$multi_minus_uni_latitude_deg),
      mean_multi_minus_uni_latitude_deg = mean(variable_shift$multi_minus_uni_latitude_deg),
      wilcoxon_V = if (nrow(variable_shift) > 1) unname(stats::wilcox.test(
        variable_shift$multi_minus_uni_latitude_deg, mu = 0, exact = FALSE)$statistic) else NA_real_,
      wilcoxon_p = if (nrow(variable_shift) > 1) stats::wilcox.test(
        variable_shift$multi_minus_uni_latitude_deg, mu = 0, exact = FALSE)$p.value else NA_real_
    )
  } else {
    data.frame(
      n_variable_species = 0, n_multi_south = 0, proportion_multi_south = NA_real_,
      median_multi_minus_uni_latitude_deg = NA_real_, mean_multi_minus_uni_latitude_deg = NA_real_,
      wilcoxon_V = NA_real_, wilcoxon_p = NA_real_
    )
  }

  stable_list <- species_class |>
    dplyr::select(
      SPECIES,
      analysis_stability_class,
      source_stability_class,
      stable_across_analysed_networks,
      stable_across_full_lookup,
      n_populations, n_bms, n_uni_populations, n_multi_populations,
      latitude_min, latitude_max, latitude_span
    )

  # Figure for variable species only: every point is an analysed population.
  if (length(variable_species)) {
    order_species <- variable_shift |>
      dplyr::arrange(multi_minus_uni_latitude_deg) |>
      dplyr::pull(SPECIES)

    figdat <- offset_pop |>
      dplyr::filter(SPECIES %in% variable_species, !is.na(voltinism)) |>
      dplyr::mutate(SPECIES = factor(SPECIES, levels = rev(order_species)))

    p <- ggplot2::ggplot(figdat, ggplot2::aes(latitude_deg, SPECIES, colour = voltinism)) +
      ggplot2::geom_point(alpha = .28, size = 1.2, position = ggplot2::position_jitter(height = .10, width = 0)) +
      ggplot2::scale_colour_manual(values = c(univoltine = "#2F5AA8", multivoltine = "#D55E00")) +
      ggplot2::theme_classic(base_size = 14, base_family = "Garamond") +
      ggplot2::labs(
        x = "Latitude (°)", y = NULL, colour = NULL,
        title = "Species with geographically variable voltinism"
      ) +
      ggplot2::theme(legend.position = "top")

    ggplot2::ggsave(file.path(out, "variable_species_voltinism_latitude.png"),
                    p, width = 10, height = max(5, 0.28 * length(variable_species) + 2),
                    dpi = 350, bg = "white", limitsize = FALSE)
  }

  csv(species_class, file.path(out, "species_voltinism_classification.csv"))
  csv(stable_list, file.path(out, "species_stable_for_sensitivity.csv"))
  csv(offset_pop, file.path(out, "population_voltinism_geography.csv"))
  csv(variable_detail, file.path(out, "variable_species_by_BMS.csv"))
  csv(variable_shift, file.path(out, "variable_species_latitude_shift.csv"))
  csv(class_summary, file.path(out, "summary_species_classes.csv"))
  csv(pop_summary, file.path(out, "summary_population_classes.csv"))
  csv(shift_summary, file.path(out, "summary_variable_species_latitude_shift.csv"))

  writeLines(c(
    "phenoIMPACT voltinism audit",
    paste("Voltinism lookup:", volt_file),
    paste("Phenology populations:", pheno_file),
    paste("Coordinate source:", coord_file),
    "",
    "analysis_stability_class is defined using BMS units represented by each species in the analysed OFFSET population universe.",
    "source_stability_class uses all non-missing uni/multi classifications in species_country_voltinism.csv.",
    "'not present', blank and other non-binary values are treated as unavailable, not as a third voltinism class.",
    "For variable species, multi_minus_uni_latitude_deg < 0 means multivoltine populations occur farther south on average.",
    "The Wilcoxon test uses one latitude contrast per variable species and is descriptive; BMS classifications are coarse regional assignments, not measured realized generations per site-year."
  ), file.path(out, "README.txt"))

  if (requireNamespace("zip", quietly = TRUE)) {
    rel <- list.files(out, recursive = TRUE, full.names = FALSE)
    rel <- rel[grepl("\\.(csv|txt|png)$", rel, ignore.case = TRUE)]
    zip::zipr(file.path(out, "results_to_review.zip"), rel,
              root = out, mode = "mirror", include_directories = FALSE)
  }

  cat("\nSpecies classes:\n")
  print(class_summary, row.names = FALSE)
  cat("\nVariable-species geographic summary:\n")
  print(shift_summary, row.names = FALSE)
  cat("\nOutput:\n", normalizePath(out, winslash = "/"), "\n", sep = "")
  if (file.exists(file.path(out, "results_to_review.zip")))
    cat("FILE TO SEND:\n", normalizePath(file.path(out, "results_to_review.zip"), winslash = "/"), "\n", sep = "")

  invisible(list(
    output = out,
    species = species_class,
    variable_species = variable_shift,
    summary = class_summary,
    zip = file.path(out, "results_to_review.zip")
  ))
}
