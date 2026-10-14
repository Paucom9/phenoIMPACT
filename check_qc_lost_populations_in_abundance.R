# phenoIMPACT | Check QC-lost offset populations against abundance models
# NO MODEL IS FITTED.

check_qc_lost_populations_in_abundance <- function(
  root = "E:/phenoIMPACT project/code/phenoIMPACT",
  plasticity_comparison_dir = file.path(
    root, "output", "phenology_plasticity",
    "offset_with_onset_qc_lastpeak_15km", "run_fbf8b3101079",
    "plasticity_comparison_pre_abundance"
  ),
  abundance_reference = file.path(
    root, "output", "population_trends",
    "spatial_abundance_15km", "run_650c41f15949"
  ),
  current_abundance_run = file.path(
    root, "output", "population_trends",
    "spatial_abundance_onset_adjusted_15km", "run_20261004_214012_44136"
  ),
  out_dir = file.path(
    plasticity_comparison_dir, "abundance_overlap_check"
  )
) {
  groups <- c("univoltine", "multivoltine")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  need <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)

  read_pop <- function(path) {
    need(file.exists(path), paste("Missing file:", path))
    x <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    need("pop_id" %in% names(x), paste("pop_id missing from:", path))
    unique(as.character(x$pop_id))
  }

  abundance_pop_ids <- function(run, group) {
    p <- file.path(run, paste0("input_", group, ".rds"))
    need(file.exists(p), paste("Missing abundance input:", p))
    x <- readRDS(p)
    need(!is.null(x$frame) && "pop_id" %in% names(x$frame),
         paste("Expected frame$pop_id in:", p))
    unique(as.character(x$frame$pop_id))
  }

  summaries <- list()
  details <- list()

  for (g in groups) {
    old_ids <- read_pop(file.path(
      plasticity_comparison_dir,
      paste0("old_population_plasticity_", g, ".csv")
    ))
    qc_ids <- read_pop(file.path(
      plasticity_comparison_dir,
      paste0("qc_population_plasticity_", g, ".csv")
    ))

    lost <- setdiff(old_ids, qc_ids)
    ref_ids <- abundance_pop_ids(abundance_reference, g)
    overlap_ref <- intersect(lost, ref_ids)

    current_file <- file.path(current_abundance_run, paste0("input_", g, ".rds"))
    if (file.exists(current_file)) {
      current_ids <- abundance_pop_ids(current_abundance_run, g)
      same_current_universe <- setequal(ref_ids, current_ids)
      overlap_current <- intersect(lost, current_ids)
    } else {
      current_ids <- character(0)
      same_current_universe <- NA
      overlap_current <- character(0)
    }

    summaries[[g]] <- data.frame(
      group = g,
      n_old_offset_populations = length(old_ids),
      n_qc_offset_populations = length(qc_ids),
      n_populations_lost_by_QC = length(lost),
      n_abundance_populations_reference = length(ref_ids),
      n_lost_present_in_abundance_reference = length(overlap_ref),
      percent_lost_present_in_abundance_reference =
        if (length(lost)) 100 * length(overlap_ref) / length(lost) else 0,
      current_abundance_input_found = file.exists(current_file),
      current_vs_reference_same_population_universe = same_current_universe,
      n_lost_present_in_current_abundance =
        if (file.exists(current_file)) length(overlap_current) else NA_integer_,
      stringsAsFactors = FALSE
    )

    details[[g]] <- data.frame(
      group = rep(g, length(lost)),
      pop_id = lost,
      used_in_reference_abundance = lost %in% ref_ids,
      used_in_current_onset_adjusted_abundance =
        if (file.exists(current_file)) lost %in% current_ids else rep(NA, length(lost)),
      stringsAsFactors = FALSE
    )
  }

  summary <- do.call(rbind, summaries)
  detail <- do.call(rbind, details)
  used <- detail[detail$used_in_reference_abundance %in% TRUE, , drop = FALSE]

  utils::write.csv(summary,
    file.path(out_dir, "qc_lost_population_abundance_overlap_summary.csv"),
    row.names = FALSE)
  utils::write.csv(detail,
    file.path(out_dir, "qc_lost_population_abundance_overlap_detail.csv"),
    row.names = FALSE)
  utils::write.csv(used,
    file.path(out_dir, "qc_lost_populations_used_in_abundance.csv"),
    row.names = FALSE)

  cat("\n==============================================================\n")
  cat("QC-LOST OFFSET POPULATIONS x ABUNDANCE UNIVERSE\n")
  cat("NO MODEL FITTED\n")
  cat("==============================================================\n\n")
  print(summary, row.names = FALSE)
  cat("\nLost QC populations actually used in abundance models: ", nrow(used), "\n", sep = "")
  if (nrow(used)) {
    cat("\nThese populations need a decision before the final abundance rerun:\n")
    print(used, row.names = FALSE)
  } else {
    cat("\nExcellent: none of the QC-lost populations enter the abundance models.\n")
  }
  cat("\nOutput: ", normalizePath(out_dir, winslash = "/", mustWork = TRUE), "\n", sep = "")

  invisible(list(summary = summary, detail = detail, used = used, output = out_dir))
}

qc_abundance_overlap <- check_qc_lost_populations_in_abundance()
