# Shared read-only source and worker checks for the six phylogenetic sensitivities.
# Sourcing defines an API only. No prepared data or fitted models are read here.
pheno_mundlak_phylo_sources <- local({
  analyses <- c("onset_multivoltine", "offset_multivoltine",
    "onset_univoltine", "offset_univoltine")
  need <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
  metadata <- function(path) {
    s <- file.info(path)
    need(nrow(s) == 1L && is.finite(s$size) && s$size > 0 && s$size < 5 * 1024^2,
      paste("Missing, empty or oversized metadata:", path))
    readRDS(path)
  }
  geo <- function(root) {
    e <- new.env(parent = globalenv())
    sys.source(file.path(root, "geography_mundlak_by_voltinism_qc_final.R"), envir = e)
    g <- environment(e$pheno_geo_mundlak$run)
    need(identical(g$VERSION, "geography_mundlak_by_voltinism_qc_final_15km_1.1.0"),
      "The v1/v2 geographic source checker version 1.1.0 is required.")
    g
  }
  manifest <- function(paths) {
    paths <- sort(unique(normalizePath(paths, winslash = "/", mustWork = TRUE)))
    z <- file.info(paths)
    need(!anyNA(z$size) && all(!z$isdir), "Source metadata unavailable.")
    data.frame(path = paths, bytes = z$size, modified = as.numeric(z$mtime),
      stringsAsFactors = FALSE)
  }
  environment_row <- function(root, run, analysis, note = "") {
    need(analysis %in% analyses, "Unknown environmental analysis.")
    directory <- file.path(run, analysis)
    row <- data.frame(analysis = analysis, source_directory = directory,
      source_checkpoint = file.path(directory, "fit_mundlak.rds"),
      prepared_file = file.path(directory, "prepared_input.rds"),
      original_completed_file = file.path(directory, "completed.rds"),
      completed_file = file.path(directory, "completed.rds"),
      fit_checks_file = file.path(directory, "fit_checks.csv"),
      run_config_file = file.path(run, "run_config.rds"),
      environment_run = run, kind = "original", refinement_version = "none",
      version = NA_character_, status = "NOT_COMPLETED", run_signature = NA_character_,
      numerical_checks_ok = FALSE, accepted = FALSE, reason = note,
      validation_scope = "metadata_eligibility_only", stringsAsFactors = FALSE)
    tryCatch({
      src <- geo(root)$source_status(run)
      z <- src$table[src$table$analysis == analysis, , drop = FALSE]
      need(nrow(z) == 1L, "Ambiguous environmental selection.")
      row$source_checkpoint <- z$fit
      row$prepared_file <- z$prepared
      row$completed_file <- z$marker
      row$fit_checks_file <- z$checks
      row$kind <- z$kind
      row$refinement_version <- z$refinement_version
      row$version <- src$config$version
      row$status <- z$status
      row$run_signature <- src$config$signature
      signature_ok <- grepl("^[a-f0-9]{64}$", src$config$signature)
      row$numerical_checks_ok <- z$ready
      paths <- unlist(row[c("source_checkpoint", "prepared_file",
        "original_completed_file", "completed_file", "fit_checks_file", "run_config_file")])
      positive <- all(file.exists(paths)) && all(file.info(paths)$size > 0)
      row$accepted <- isTRUE(z$ready) && signature_ok && positive && !nzchar(note)
      row$reason <- if (row$accepted) "metadata eligible; fitted arrays have not been read" else
        paste(c(note, z$status, if (!signature_ok) "invalid run signature",
          if (!positive) "missing or empty source artifact"), collapse = "; ")
      row
    }, error = function(e) {
      row$reason <- paste(note, conditionMessage(e), sep = "; ")
      row
    })
  }
  environment_files <- function(root, run, analysis) {
    g <- geo(root)
    row <- environment_row(root, run, analysis)
    need(isTRUE(row$accepted), row$reason)
    unique(c(unlist(row[c("source_checkpoint", "prepared_file", "original_completed_file",
      "completed_file", "fit_checks_file", "run_config_file")], use.names = FALSE),
      g$refinement_roots(run)$config_files))
  }
  post_versions <- function() c("post_environment_mundlak_1.1")
  post_environment <- function(root, config) {
    need(config$version %in% post_versions(), "Unsupported post-environment source version.")
    need(is.character(config$environment_run) && length(config$environment_run) == 1L,
      "Post-environment provenance is missing.")
    current <- geo(root)$source_status(config$environment_run, TRUE)$table
    saved <- config$environment
    columns <- c("analysis", "kind", "refinement_version", "ready",
      "fit", "prepared", "checks", "marker", "original_marker")
    need(is.data.frame(saved) && nrow(saved) == 4L && !anyDuplicated(saved$analysis) &&
      setequal(saved$analysis, analyses) && all(columns %in% names(saved)),
      "Post-environment run lacks the four selected environmental sources.")
    saved <- saved[match(current$analysis, saved$analysis), columns, drop = FALSE]
    current <- current[, columns, drop = FALSE]
    rownames(saved) <- rownames(current) <- NULL
    need(identical(saved, current), "Post-environment environmental selection is stale.")
    need(is.data.frame(config$manifest) &&
      all(c("path", "bytes", "modified") %in% names(config$manifest)),
      "Post-environment source manifest is missing.")
    recorded <- config$manifest
    actual <- manifest(recorded$path)
    recorded <- recorded[order(recorded$path), c("path", "bytes", "modified")]
    rownames(recorded) <- rownames(actual) <- NULL
    need(identical(actual, recorded), "Post-environment upstream source metadata changed.")
    invisible(current)
  }
  demography_row <- function(root, run, group, note = "") {
    directory <- file.path(run, paste0(group, "__mundlak"))
    row <- data.frame(group = group, source_directory = directory,
      source_checkpoint = file.path(directory, "fit.rds"),
      prepared_file = file.path(run, paste0("prepared_mundlak_", group, ".rds")),
      run_config_file = file.path(run, "run_config.rds"),
      fit_checks_file = file.path(directory, "fit_checks.csv"),
      result_file = file.path(run, "result_fit.rds"),
      preparation_marker = file.path(run, paste0("prepared_", group, "_completed.rds")),
      source_kind = "post_environment", version = NA_character_, run_signature = NA_character_,
      status = "NOT_COMPLETED", numerical_checks_ok = FALSE, accepted = FALSE,
      expected_nobs = NA_real_, expected_npop = NA_real_, reason = note,
      validation_scope = "metadata_eligibility_only", stringsAsFactors = FALSE)
    tryCatch({
      config <- metadata(row$run_config_file)
      row$version <- config$version
      row$run_signature <- config$signature
      need(config$version %in% post_versions() && identical(config$mode, "fit") &&
        is.character(config$signature) && length(config$signature) == 1L &&
        grepl("^[a-f0-9]{64}$", config$signature), "An accepted post-environment fit run is required.")
      post_environment(root, config)
      result <- metadata(row$result_file)
      prep <- metadata(row$preparation_marker)
      need(identical(result$signature, config$signature) &&
        identical(result$status, "FITTED_RESIDUAL_REVIEW_PENDING") &&
        is.data.frame(result$groups) && !anyDuplicated(result$groups$group) &&
        setequal(result$groups$group, c("multivoltine", "univoltine")) &&
        all(result$groups$status %in% c("FIT_OK", "FIT_OK_CACHED")),
        "The two post-environment demographic fits have not both passed.")
      need(identical(prep$signature, config$signature) && identical(prep$group, group) &&
        length(prep$nobs) == 1L && is.finite(prep$nobs) && prep$nobs > 0 &&
        length(prep$npop) == 1L && is.finite(prep$npop) && prep$npop > 0 &&
        prep$nobs == floor(prep$nobs) && prep$npop == floor(prep$npop) && prep$npop <= prep$nobs,
        "Post-environment preparation sample/identity is invalid.")
      row$expected_nobs <- prep$nobs
      row$expected_npop <- prep$npop
      for (a in analyses) {
        mark <- metadata(file.path(run, "extraction", a, "completed.rds"))
        need(identical(mark$signature, config$signature) && identical(mark$analysis, a),
          "Post-environment extraction provenance mismatch.")
      }
      size <- file.info(row$fit_checks_file)$size
      need(is.finite(size) && size > 0 && size < 2 * 1024^2, "Missing/oversized demographic checks.")
      checks <- utils::read.csv(row$fit_checks_file, stringsAsFactors = FALSE)
      need(nrow(checks) == 1L && identical(checks$group, group) &&
        isTRUE(checks$numerical_checks_ok) && isTRUE(checks$convergence_code == 0) &&
        isTRUE(checks$pdHess) && isTRUE(is.finite(checks$gradient_covariance_norm)) &&
        isTRUE(checks$gradient_covariance_norm < .05) && isTRUE(checks$nobs == prep$nobs),
        "Post-environment numerical checks/sample do not pass.")
      required <- c(row$source_checkpoint, row$prepared_file)
      need(all(file.exists(required)) && all(file.info(required)$size > 0),
        "Post-environment fit/prepared file is missing or empty.")
      row$status <- "FIT_OK"
      row$numerical_checks_ok <- TRUE
      row$accepted <- !nzchar(note)
      row$reason <- if (row$accepted) "post-environment metadata eligible; residual review remains separate" else note
      row
    }, error = function(e) {
      row$reason <- paste(note, conditionMessage(e), sep = "; ")
      row
    })
  }
  demography_files <- function(root, run, group) {
    row <- demography_row(root, run, group)
    need(isTRUE(row$accepted), row$reason)
    config <- metadata(row$run_config_file)
    unique(c(unlist(row[c("source_checkpoint", "prepared_file", "run_config_file",
      "fit_checks_file", "result_file", "preparation_marker")], use.names = FALSE),
      file.path(run, "extraction", analyses, "completed.rds"), config$manifest$path,
      unlist(lapply(analyses, function(a) environment_files(root, config$environment_run, a)),
        use.names = FALSE)))
  }
  idle <- function(root, environment_run = NULL) {
    jobs <- geo(root)$active_jobs(root, environment_run)
    need(!length(jobs), paste("An independent/environmental/refinement worker is active:",
      paste(jobs, collapse = "; ")))
    invisible(TRUE)
  }
  snapshot <- function(root, family, source_directory, analysis) {
    run <- dirname(source_directory)
    if (family == "demography" && identical(metadata(file.path(run, "run_config.rds"))$version,
      "mundlak_geography_plasticity_qc_final_spatial_15km_2.0.2")) {
      return(manifest(c(file.path(source_directory, c("fit.rds", "fit_checks.csv")),
        file.path(run, c("run_config.rds", paste0("prepared_mundlak_", analysis, ".rds"))))))
    }
    paths <- if (family == "environment") environment_files(root, run, analysis) else
      demography_files(root, run, analysis)
    manifest(paths)
  }
  guard <- function(root, family, source_directory, analysis, expected = NULL) {
    run <- dirname(source_directory)
    env_run <- if (family == "environment") run else metadata(file.path(run, "run_config.rds"))$environment_run
    idle(root, env_run)
    current <- snapshot(root, family, source_directory, analysis)
    if (!is.null(expected)) need(identical(current, expected),
      "Selected source changed after audit; restart the phylogenetic launch.")
    invisible(current)
  }
  list(analyses = analyses, environment_row = environment_row, environment_files = environment_files,
    demography_row = demography_row, demography_files = demography_files,
    post_versions = post_versions, post_environment = post_environment,
    metadata = metadata, idle = idle, snapshot = snapshot, guard = guard)
})
