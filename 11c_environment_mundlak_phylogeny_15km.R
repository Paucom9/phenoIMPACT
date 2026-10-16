# phenoIMPACT | Phylogenetic sensitivity of the environmental Mundlak model (11b).
# Version 1.0.0 | 2026-10-15. Loading defines functions; no analysis starts.
# This models phenological dates, not the demographic response in script 10f.
# One selected analysis at a time: onset, offset_univoltine, offset_multivoltine.
# Requires its ALREADY FIT_OK independent environmental-Mundlak fit, version 1.1/1.2/1.3.
# Reuses every row, scale, fixed interaction, IID random effect and the exact 15-km mesh.
# Adds only an anomaly slope with covariance sigma_phylo^2 A alongside sigma_iid^2 I.
# A retains the original-root Brownian covariance, normalized to diagonal one.
# Requires existing pheno_phylogeny_qc_15km_v2.R and pheno_plasticity_phylogeny_pilot_15km.R.
# No package installation or updates are performed.
#
# source("11c_environment_mundlak_phylogeny_15km.R", encoding = "UTF-8")
# job <- run_environment_mundlak_phylogeny(check_only = TRUE)
# watch_environment_mundlak_phylogeny(job)
# job <- run_environment_mundlak_phylogeny(check_only = FALSE)
# watch_environment_mundlak_phylogeny(job)
# Repeat the same run call to resume, including after restarting R.
# If multiple independent sources exist, specify source_run = ".../run_<hash>".
#
# Checkpoints: portable preparation, best finite outer evaluations, completed
# optimizer BEFORE covariance/Hessian reporting, then fitted checkpoint.
# A native Hessian construction/evaluation cannot be checkpointed mid-call;
# that particular call must be rebuilt after a crash. Elapsed time alone does
# not imply that an optimizer checkpoint has yet been reached.
#
# Numerical gate: convergence zero, positive Hessian and covariance-scaled
# gradient norm sqrt(g' V g) < 0.05. Raw gradient < 0.001 remains an advisory
# diagnostic. Both are recorded; old fit/check files are not modified.
# Approximate Wald/delta intervals and ML/AIC are descriptive sensitivity
# diagnostics. No automatic boundary-LRT p value or tree-uncertainty propagation.

pheno_env_mundlak_phylogeny <- local({
# Source/prepared-cache helpers only. Source into the private phylogeny engine.
# Metadata discovery never opens a fitted-model RDS or launches a model.

env_phylo_need <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
}

env_phylo_spec <- function(analysis) {
  env_phylo_need(is.character(analysis) && length(analysis) == 1L &&
    analysis %in% c("onset", "offset_univoltine", "offset_multivoltine"),
    "Select exactly one environmental-Mundlak analysis.")
  window <- if (analysis == "onset") 60L else 90L
  list(analysis = analysis, response = if (analysis == "onset") "ONSET_mean" else "OFFSET_mean",
    window = window, anomaly = paste0("clim_anomaly_tw", window),
    orientation = if (analysis == "offset_multivoltine") 1 else -1,
    orientation_label = if (analysis == "offset_multivoltine") "delay / extension" else "advancement",
    control = if (analysis == "onset") character() else "ONSET_mean_z",
    moderators = c("photoperiod_within_species", "photoperiod_species_mean",
      "background_within_species", "background_species_mean",
      "predictability_within_species", "predictability_species_mean",
      "trend_within_species", "trend_species_mean"))
}

env_phylo_read_metadata <- function(path) {
  info <- file.info(path)
  env_phylo_need(nrow(info) == 1L && is.finite(info$size) && info$size > 0 &&
    info$size <= 2 * 1024^2, paste("Missing or oversized metadata:", path))
  readRDS(path)
}

env_phylo_source_versions <- function() {
  paste0("environment_mundlak_qc_final_spatial_15km_", c("1.1", "1.2", "1.3"))
}

env_phylo_audit_source_row <- function(directory, analysis, discovery_note = "") {
  files <- c(source_checkpoint = file.path(directory, "fit_mundlak.rds"),
    prepared_file = file.path(directory, "prepared_input.rds"),
    completed_file = file.path(directory, "completed.rds"),
    fit_checks_file = file.path(directory, "fit_checks.csv"),
    run_config_file = file.path(dirname(directory), "run_config.rds"))
  info <- file.info(files)
  present <- stats::setNames(file.exists(files) & !is.na(info$isdir) & !info$isdir, names(files))
  sizes <- as.numeric(info$size)
  missing <- basename(files[!present])
  empty <- basename(files[present & (is.na(sizes) | sizes <= 0)])
  errors <- character()
  issues <- if (nzchar(discovery_note)) discovery_note else character()
  scalar_character <- function(x) if (is.character(x) && length(x) == 1L &&
    !is.na(x) && nzchar(x)) x else NA_character_
  read_metadata <- function(name) {
    if (!isTRUE(present[[name]])) return(NULL)
    tryCatch({
      x <- env_phylo_read_metadata(files[[name]])
      if (!is.list(x)) stop("Expected a metadata list.")
      x
    }, error = function(e) {
      errors <<- c(errors, paste0(basename(files[[name]]), ": ", conditionMessage(e)))
      NULL
    })
  }
  config <- read_metadata("run_config_file")
  done <- read_metadata("completed_file")
  # Fitted-model and prepared-data RDS files are deliberately never opened.
  checks <- if (isTRUE(present[["fit_checks_file"]])) tryCatch({
    j <- match("fit_checks_file", names(files))
    if (!is.finite(sizes[[j]]) || sizes[[j]] <= 0 || sizes[[j]] > 2 * 1024^2)
      stop("Empty or oversized checks CSV.")
    utils::read.csv(files[["fit_checks_file"]], stringsAsFactors = FALSE)
  }, error = function(e) {
    errors <<- c(errors, paste0("fit_checks.csv: ", conditionMessage(e)))
    NULL
  }) else NULL
  source_version <- scalar_character(config$version)
  status <- scalar_character(done$status)
  run_signature <- scalar_character(config$signature)
  completed_signature <- scalar_character(done$signature)
  signature_valid <- !is.na(run_signature) && grepl("^[a-f0-9]{64}$", run_signature)
  signature_match <- signature_valid && !is.na(completed_signature) &&
    identical(config$signature, done$signature)
  version_ok <- !is.na(source_version) && source_version %in% env_phylo_source_versions()
  status_ok <- !is.na(status) && identical(status, "FIT_OK")
  checks_rows <- if (is.data.frame(checks)) nrow(checks) else NA_integer_
  numerical <- is.data.frame(checks) && nrow(checks) == 1L &&
    isTRUE(checks$numerical_checks_ok)
  if (length(missing)) issues <- c(issues, paste("missing files:", paste(missing, collapse = ", ")))
  if (length(empty)) issues <- c(issues, paste("empty files:", paste(empty, collapse = ", ")))
  if (length(errors)) issues <- c(issues, paste("metadata read errors:", paste(errors, collapse = " | ")))
  if (!version_ok) issues <- c(issues, paste0("unsupported/missing source version: ",
    if (is.na(source_version)) "[missing]" else source_version,
    "; accepted scientific versions are 1.1, 1.2 and 1.3"))
  if (!status_ok) issues <- c(issues, paste("completion requires FIT_OK; found",
    if (is.na(status)) "[missing]" else status))
  if (!signature_valid) issues <- c(issues, "run signature is missing or is not a 64-character SHA256")
  if (!signature_match) issues <- c(issues, "run/completion signatures do not match")
  if (!numerical) issues <- c(issues, "fit_checks.csv requires exactly one row with numerical_checks_ok=TRUE")
  accepted <- !length(issues)
  data.frame(analysis = analysis,
    source_directory = normalizePath(directory, winslash = "/", mustWork = FALSE),
    source_checkpoint = unname(files[["source_checkpoint"]]),
    prepared_file = unname(files[["prepared_file"]]),
    run_config_file = unname(files[["run_config_file"]]),
    completed_file = unname(files[["completed_file"]]),
    fit_checks_file = unname(files[["fit_checks_file"]]),
    missing_files = paste(missing, collapse = "; "),
    source_checkpoint_bytes = sizes[[1L]], prepared_bytes = sizes[[2L]],
    completed_bytes = sizes[[3L]], fit_checks_bytes = sizes[[4L]], run_config_bytes = sizes[[5L]],
    read_errors = paste(errors, collapse = " | "), version = source_version,
    status = status, run_signature = run_signature, completed_signature = completed_signature,
    run_signature_valid = signature_valid, signature_match = signature_match,
    checks_rows = checks_rows, numerical_checks_ok = numerical,
    accepted = accepted, reason = if (accepted) "metadata eligible; fitted arrays have not been read or validated" else
      paste(unique(issues), collapse = "; "),
    validation_scope = "metadata_eligibility_only", full_fit_validation_run = FALSE,
    stringsAsFactors = FALSE, row.names = NULL)
}

env_phylo_audit_sources <- function(root, analysis = NULL, source_run = NULL) {
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  analyses <- c("onset", "offset_univoltine", "offset_multivoltine")
  if (!is.null(analysis)) {
    env_phylo_need(is.character(analysis) && length(analysis) > 0L &&
      !anyNA(analysis) && all(analysis %in% analyses), "Unknown environmental-Mundlak analysis.")
    analyses <- unique(analysis)
  }
  note <- ""
  explicit_analysis <- explicit_directory <- NULL
  if (is.null(source_run)) {
    parent <- file.path(root, "output/phenology_plasticity/environment_mundlak_qc_final_15km")
    runs <- if (dir.exists(parent)) list.dirs(parent, recursive = FALSE, full.names = TRUE) else character()
    if (!length(runs)) {
      runs <- parent
      note <- if (dir.exists(parent)) "no candidate independent run folders found" else
        "independent environmental-Mundlak output root does not exist"
    }
  } else {
    env_phylo_need(is.character(source_run) && length(source_run) == 1L &&
      !is.na(source_run) && nzchar(source_run), "Invalid explicit independent source folder.")
    path <- normalizePath(source_run, winslash = "/", mustWork = FALSE)
    if (basename(path) %in% c("onset", "offset_univoltine", "offset_multivoltine")) {
      if (is.null(analysis)) analyses <- basename(path)
      explicit_analysis <- basename(path)
      explicit_directory <- path
      runs <- dirname(path)
    } else runs <- path
  }
  rows <- lapply(runs, function(run) lapply(analyses, function(a) {
    row_note <- if (!is.null(explicit_analysis) && !identical(a, explicit_analysis))
      paste0("explicit source folder is for ", explicit_analysis, "; requested ", a) else note
    directory <- if (is.null(explicit_directory)) file.path(run, a) else explicit_directory
    env_phylo_audit_source_row(directory, a, row_note)
  }))
  result <- do.call(rbind, unlist(rows, recursive = FALSE))
  rownames(result) <- NULL
  attr(result, "validation_scope") <- "Metadata eligibility only: no model/prepared RDS, native callbacks, or fits are read/run."
  result
}

env_phylo_resolve_source <- function(root, analysis, source_run = NULL) {
  env_phylo_spec(analysis)
  audit <- env_phylo_audit_sources(root, analysis, source_run)
  accepted <- which(audit$accepted)
  if (length(accepted) != 1L) {
    details <- vapply(seq_len(nrow(audit)), function(i) paste0(
      audit$source_directory[[i]], " | version=", audit$version[[i]],
      " | status=", audit$status[[i]], " | ", audit$reason[[i]]), character(1))
    stop(paste0("Expected exactly one metadata-eligible already FIT_OK independent environmental-Mundlak source for ",
      analysis, "; found ", length(accepted), ".\n", paste(details, collapse = "\n"),
      "\nInspect audit_environment_mundlak_phylogeny_sources() (public audit), or env_phylo_audit_sources(root, analysis, source_run).",
      "\nSupply an explicit compatible source_run if several are eligible. No baseline is fitted automatically."), call. = FALSE)
  }
  row <- audit[accepted, , drop = FALSE]
  list(source_directory = row$source_directory[[1L]], source_checkpoint = row$source_checkpoint[[1L]],
    prepared_file = row$prepared_file[[1L]], completed_file = row$completed_file[[1L]],
    fit_checks_file = row$fit_checks_file[[1L]], run_config_file = row$run_config_file[[1L]],
    analysis = analysis, source_version = row$version[[1L]], source_signature = row$run_signature[[1L]])
}

env_phylo_paths <- function(cfg) {
  resolved <- env_phylo_resolve_source(cfg$project_root, cfg$analysis, cfg$source_directory)
  if (!is.null(cfg$source_checkpoint)) env_phylo_need(identical(
    normalizePath(cfg$source_checkpoint, winslash = "/", mustWork = TRUE),
    resolved$source_checkpoint), "Source directory and source checkpoint disagree.")
  resolved
}

env_phylo_formula_text <- function(formula) {
  paste(deparse(formula, width.cutoff = 500L), collapse = " ")
}

env_phylo_data_equal <- function(x, y, columns) {
  if (!is.data.frame(x) || !is.data.frame(y) || nrow(x) != nrow(y) ||
      !all(columns %in% names(x)) || !all(columns %in% names(y))) return(FALSE)
  all(vapply(columns, function(name) {
    a <- x[[name]]; b <- y[[name]]
    if (is.factor(a) || is.factor(b) || is.character(a) || is.character(b))
      identical(as.character(a), as.character(b))
    else identical(unname(as.vector(a)), unname(as.vector(b)))
  }, logical(1)))
}

env_phylo_prepare_uncached <- function(cfg, packages) {
  resolved <- env_phylo_paths(cfg)
  s <- env_phylo_spec(cfg$analysis)
  env_phylo_need(identical(anomaly, s$anomaly) && setequal(moderators, s$moderators),
    "Engine anomaly/moderator settings do not match the selected analysis.")
  run <- env_phylo_read_metadata(resolved$run_config_file)
  done <- env_phylo_read_metadata(resolved$completed_file)
  prepared <- readRDS(resolved$prepared_file)
  progress(paste("Reading saved independent environmental-Mundlak baseline:", cfg$analysis))
  memory_snapshot(cfg$monitor_directory, "before_independent_Mundlak_read")
  saved <- readRDS(resolved$source_checkpoint)
  # Strip large/unusable native state BEFORE taking a fit alias.
  env_phylo_need(is.list(saved) && is.list(saved$fit) && is.list(saved$fit$sd_report),
    "Expected the environmental-Mundlak saved-fit wrapper.")
  saved$fit$sd_report[["cov"]] <- NA_real_
  saved$fit$sd_report$env <- NULL
  saved$fit$tmb_obj <- NULL
  invisible(gc())
  fit <- saved$fit
  memory_snapshot(cfg$monitor_directory, "after_independent_Mundlak_read")
  env_phylo_need(identical(saved$signature, run$signature) &&
    identical(done$signature, run$signature) && identical(done$status, "FIT_OK") &&
    identical(saved$version, run$version) && identical(prepared$version, run$version) &&
    identical(saved$analysis, cfg$analysis) && identical(prepared$analysis, cfg$analysis),
    "Independent Mundlak fit/prepared/run/completion identity mismatch.")
  env_phylo_need(inherits(fit, "sdmTMB") && as.character(fit$version) == "1.1.0" &&
    isFALSE(fit$reml) && identical(fit$time, "YEAR") &&
    all(fit$spatial == "off") && all(fit$spatiotemporal == "iid") &&
    identical(fit$family$family, "gaussian") && identical(fit$family$link, "identity") &&
    !isTRUE(as.logical(fit$tmb_data$has_smooths)) && is.null(fit$nonlocal_formula) &&
    is.null(fit$time_varying) && is.null(fit$spatial_varying),
    "Unexpected independent environmental-Mundlak model structure/version.")
  env_phylo_need(identical(fit$priors, sdmTMB::sdmTMBpriors()),
    "Unexpected priors; marginal-ML comparison needs review.")
  if (cfg$analysis != "onset") {
    path <- gsub("\\\\", "/", prepared$source_path)
    expected_path <- paste0("/offset_with_onset_qc_lastpeak_15km/run_fbf8b3101079/",
      cfg$analysis, "/fit_qc.rds")
    env_phylo_need(is.character(path) && length(path) == 1L && !is.na(path) &&
      endsWith(path, expected_path), "Independent offset baseline is not derived from the locked LAST_PEAK QC fit.")
  } else {
    env_phylo_need(is.character(prepared$source_path) && length(prepared$source_path) == 1L &&
      basename(prepared$source_path) == "fit_spatial_field.rds" &&
      grepl("primary__onset_", prepared$source_path, fixed = TRUE),
      "Independent onset baseline is not derived from the primary onset fit.")
  }
  f0 <- if (inherits(fit$formula, "formula")) fit$formula else fit$formula[[1L]]
  env_phylo_need(inherits(f0, "formula") && inherits(prepared$formula, "formula") &&
    identical(env_phylo_formula_text(f0), env_phylo_formula_text(prepared$formula)),
    "Saved and prepared Mundlak formulas differ.")
  environment(f0) <- baseenv()
  expected <- stats::as.formula(paste(s$response, "~", s$anomaly, "*(",
    paste(s$moderators, collapse = "+"), ")",
    if (length(s$control)) paste("+", s$control) else ""), env = baseenv())
  sf <- fit$split_formula[[1L]]
  env_phylo_need(identical(all.vars(sf$form_no_bars[[2L]]), s$response) &&
    setequal(attr(stats::terms(sf$form_no_bars), "term.labels"),
      attr(stats::terms(expected), "term.labels")), "Source lacks the exact eight Mundlak moderator interactions.")
  random <- function(f) sort(vapply(reformulas::findbars(f), env_phylo_formula_text, character(1)))
  reference <- stats::as.formula(paste(s$response, "~ (1|SITE_ID)+(1|SPECIES)+(0+",
    s$anomaly, "|SPECIES_slope)+(1|site_year_id)"), env = baseenv())
  env_phylo_need(identical(random(f0), random(reference)), "Unexpected independent random terms.")
  d <- fit$data
  columns <- unique(c(all.vars(f0), "YEAR", "x_km", "y_km",
    if ("source_row" %in% names(d)) "source_row"))
  env_phylo_need(env_phylo_data_equal(d, prepared$data, columns) && !anyNA(d[, columns]) &&
    identical(as.character(d$SPECIES), as.character(d$SPECIES_slope)),
    "Prepared and fitted Mundlak rows/scales/grouping are inconsistent.")
  numeric_columns <- unique(c(s$response, s$anomaly, s$moderators, s$control, "YEAR", "x_km", "y_km"))
  env_phylo_need(all(vapply(d[numeric_columns], function(x) is.numeric(x) &&
    all(is.finite(x)), logical(1))) && all(d$YEAR == as.integer(d$YEAR)), "Invalid numeric Mundlak data.")
  if ("source_row" %in% names(d)) env_phylo_need(!anyDuplicated(d$source_row), "Duplicated original row identifiers.")
  env_phylo_need((is.null(fit$offset) || all(fit$offset == 0)) &&
    (is.null(fit$tmb_data$weights) || all(fit$tmb_data$weights == 1)), "Unexpected offsets or weights.")
  mesh <- fit$spde
  env_phylo_need(isTRUE(all.equal(mesh, prepared$mesh, tolerance = 0, check.attributes = TRUE)) &&
    nrow(mesh$loc_xy) == nrow(d) && same_num(as.vector(as.matrix(mesh$loc_xy)),
      as.vector(as.matrix(d[, c("x_km", "y_km")]))), "The exact fitted 15-km mesh/projection was not preserved.")
  env_phylo_need(same_num(as.numeric(fit$response), d[[s$response]]), "Saved response differs from fitted data.")
  scaling <- prepared$scaling
  env_phylo_need(is.data.frame(scaling) && nrow(scaling) == 4L &&
    all(c("moderator", "within_SD", "between_population_weighted_SD") %in% names(scaling)) &&
    !anyDuplicated(scaling$moderator) && setequal(scaling$moderator,
      c("photoperiod", "background", "predictability", "trend")) &&
    all(is.finite(scaling$within_SD) & scaling$within_SD > 0) &&
    all(is.finite(scaling$between_population_weighted_SD) & scaling$between_population_weighted_SD > 0),
    "Invalid original four-row Mundlak component scaling.")
  verified <- phylo_cache_validate_source(fit, gradient_threshold = cfg$gradient_threshold)
  ck <- checks(fit, cfg$gradient_threshold)
  env_phylo_need(isTRUE(ck$fit_ok), "Independent Mundlak baseline fails the original numerical gate.")
  fixed <- fixed_table(fit)
  bi <- which(names(fit$sd_report$par.fixed) == "b_j")
  V <- fit$sd_report$cov.fixed[bi, bi, drop = FALSE]
  dimnames(V) <- list(fixed$term, fixed$term)
  scaled <- sqrt(sum((chol(fit$sd_report$cov.fixed) %*% as.numeric(fit$gradients))^2))
  env_phylo_need(is.finite(scaled), "Non-finite local baseline numerical distance.")
  original_data_hash <- digest::digest(d, algo = "sha256")
  mesh_hash <- digest::digest(mesh, algo = "sha256")
  species <- levels(droplevels(factor(as.character(d$SPECIES))))
  tree <- ape::read.tree(cfg$tree_file)
  env_phylo_need(inherits(tree, "phylo") && !inherits(tree, "multiPhylo") &&
    ape::is.rooted(tree) && ape::is.ultrametric(tree) &&
    length(tree$edge.length) == nrow(tree$edge) && all(is.finite(tree$edge.length)) &&
    all(tree$edge.length > 0), "Invalid rooted ultrametric tree.")
  sp <- normalize_species(species); tips <- normalize_species(tree$tip.label)
  env_phylo_need(!anyDuplicated(sp) && !anyDuplicated(tips) && !anyNA(match(sp, tips)),
    "Incomplete or ambiguous species/tree matching.")
  tip_names <- tree$tip.label[match(sp, tips)]
  edge <- build_phylo_edge_basis(tree, tip_names)
  A <- edge$A
  dimnames(A) <- list(species, species)
  attr(A, "edge_basis") <- unname(edge$B)
  attr(A, "edge_audit") <- edge$audit
  d$PHYLO_ID <- match(as.character(d$SPECIES), species)
  basis_gib <- 8 * nrow(d) * ncol(edge$B) / 1024^3
  env_phylo_need(basis_gib <= cfg$maximum_single_basis_GiB, "The added ancestral-increment basis exceeds the configured size limit.")
  start <- verified$start
  iid <- exp(fit$sd_report$par.fixed[slope_logsd_index(fit)])
  start$ln_smooth_sigma <- matrix(log(max(.1 * iid, .01)), 1L, 1L)
  saved_checks <- utils::read.csv(resolved$fit_checks_file, stringsAsFactors = FALSE)
  source <- list(fixed = fixed, V = V,
    stats = data.frame(n_observations = nrow(d), n_parameters = length(fit$model$par),
      logLik = -fit$model$objective, AIC = 2 * fit$model$objective + 2 * length(fit$model$par)),
    checks = ck, saved_checks = saved_checks, row_hash = original_data_hash,
    mesh_hash = mesh_hash, run_id = run$signature, baseline_version = run$version,
    checkpoint = resolved$source_checkpoint, formula = f0,
    local_gradient_covariance_norm = scaled,
    fingerprint = digest::digest(list(run$signature, cfg$analysis, original_data_hash,
      mesh_hash, fit$model$par, fit$model$objective, fit$parlist, V, scaling), algo = "sha256"))
  signature <- digest::digest(list(version, cfg$analysis, source$fingerprint, A, packages,
    cfg$gradient_threshold, cfg$scaled_gradient_threshold, cfg$signature_file_hashes), algo = "sha256")
  out <- file.path(cfg$output_root, paste0(cfg$analysis, "_", substr(signature, 1L, 12L)))
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  write_csv(ck, file.path(out, "source_fit_checks.csv"))
  write_csv(saved_checks, file.path(out, "source_saved_checks.csv"))
  write_csv(scaling, file.path(out, "source_component_scaling.csv"))
  write_csv(data.frame(SPECIES = species, tip_label = tip_names, PHYLO_ID = seq_along(species)),
    file.path(out, "species_mapping.csv"))
  write_csv(edge$audit, file.path(out, "tree_increment_basis_checks.csv"))
  write_csv(edge$clades, file.path(out, "tree_increment_clades.csv"))
  write_csv(edge$original_edges, file.path(out, "tree_increment_original_edges.csv"))
  write_csv(data.frame(analysis = cfg$analysis, n_observations = nrow(d), n_species = length(species),
    mesh_reused_exactly = TRUE, single_added_basis_GiB = basis_gib,
    source_native_replay_run = FALSE), file.path(out, "input_summary.csv"))
  writeLines(c("Saved independent environmental-Mundlak baseline; exact fitted rows, scales, random terms and 15-km mesh retained.",
    "Annual IID spatial fields; Gaussian ML; only the species anomaly slope gains a phylogenetic component.",
    "Cov(species slope) = sigma_iid^2 I + sigma_phylo^2 A; original tree root retained and diagonal(A)=1.",
    "All eight within/between moderator interactions and additive annual onset control for offsets are retained.",
    paste("Independent formula:", env_phylo_formula_text(f0)),
    paste("Baseline checkpoint:", resolved$source_checkpoint),
    "Source parameter/report checks are portable only: no baseline native replay or baseline refit."),
    file.path(out, "model_specification.txt"))
  write_csv(data.frame(method = "Saved mapped arrays/report checked; no native source replay",
    max_abs_gradient = verified$max_abs_gradient, gradient_covariance_norm = scaled,
    n_outer_parameters = length(verified$outer), parameter_sha256 = digest::digest(fit$parlist, algo = "sha256")),
    file.path(out, "source_start_validation.csv"))
  result <- list(data = d, formula = add_phylo(f0, A, s$anomaly), A = A,
    start = start, mesh = mesh, source = source, component_scaling = scaling,
    analysis = cfg$analysis, spec = s, out = out, signature = signature,
    priors = fit$priors, extra_time = fit$extra_time)
  phylo_cache_validate_portable_prepared(result)
  rm(fit, saved, prepared)
  invisible(gc())
  memory_snapshot(cfg$monitor_directory, "independent_Mundlak_fit_released")
  result
}

env_phylo_content_hash <- function(prepared) {
  x <- prepared
  x$out <- NULL
  digest::digest(x, algo = "sha256")
}

env_phylo_prepare <- function(cfg, packages) {
  resolved <- env_phylo_paths(cfg)
  cfg <- utils::modifyList(cfg, resolved)
  paths <- unlist(resolved[c("source_checkpoint", "prepared_file", "run_config_file",
    "completed_file", "fit_checks_file")], use.names = FALSE)
  paths <- c(paths, normalizePath(cfg$tree_file, winslash = "/", mustWork = TRUE))
  if (length(cfg$signature_files)) paths <- c(paths,
    normalizePath(cfg$signature_files, winslash = "/", mustWork = TRUE))
  hashes <- unname(tools::md5sum(paths))
  env_phylo_need(!anyNA(hashes), "Cannot fingerprint environmental-Mundlak source files.")
  manifest <- list(schema = "environment_mundlak_phylo_source_1.0", engine_version = version,
    analysis = cfg$analysis, source_signature = resolved$source_signature,
    source_version = resolved$source_version, paths = paths, md5 = hashes, packages = packages,
    raw_gradient_threshold = cfg$gradient_threshold,
    scaled_gradient_threshold = cfg$scaled_gradient_threshold,
    implementation_hashes = cfg$signature_file_hashes)
  key <- digest::digest(manifest, algo = "sha256")
  parent <- file.path(cfg$output_root, "source_cache")
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  cache_file <- file.path(parent, paste0("prepared_", substr(key, 1L, 12L), ".rds"))
  reused <- file.exists(cache_file)
  if (reused) {
    progress("Reading validated portable environmental-Mundlak source cache; no source replay")
    cached <- readRDS(cache_file)
    env_phylo_need(identical(cached$manifest, manifest), "Portable Mundlak source-cache manifest mismatch.")
    p <- cached$prepared
    phylo_cache_validate_portable_prepared(p)
    env_phylo_need(identical(env_phylo_content_hash(p), cached$fingerprint),
      "Portable Mundlak source-cache content mismatch.")
    env_phylo_need(8 * nrow(p$data) * ncol(attr(p$A, "edge_basis")) / 1024^3 <=
      cfg$maximum_single_basis_GiB, "Cached basis exceeds configured size limit.")
    p$out <- file.path(cfg$output_root, basename(p$out))
    dir.create(p$out, recursive = TRUE, showWarnings = FALSE)
  } else {
    p <- env_phylo_prepare_uncached(cfg, packages)
    save_atomic(list(manifest = manifest, fingerprint = env_phylo_content_hash(p), prepared = p), cache_file)
  }
  write_csv(data.frame(cache_file = cache_file, manifest_sha256 = key, reused = reused,
    native_source_replay_run = FALSE), file.path(p$out, "source_cache_used.csv"))
  p
}


env_phylo_numeric_equal <- function(a, b, tolerance = 1e-7) {
  is.numeric(a) && is.numeric(b) && length(a) == length(b) &&
    all(is.finite(a)) && all(is.finite(b)) &&
    (!length(a) || max(abs(a - b)) <= tolerance)
}

env_phylo_covariance <- function(V, n, labels = NULL) {
  env_phylo_need(is.matrix(V) && is.numeric(V) && identical(dim(V), c(n, n)) &&
    all(is.finite(V)), "Fixed covariance has invalid dimensions or values.")
  if (!is.null(labels) && !is.null(dimnames(V))) {
    env_phylo_need(identical(rownames(V), labels) && identical(colnames(V), labels),
      "Fixed covariance parameter ordering differs from its estimates.")
  }
  magnitude <- max(abs(V))
  env_phylo_need(magnitude > 0 && max(abs(V - t(V))) / magnitude <= 1e-8,
    "Fixed covariance is not symmetric within rounding tolerance.")
  V <- (V + t(V)) / 2
  env_phylo_need(all(diag(V) > 0), "Fixed covariance has non-positive diagonal entries.")
  R <- tryCatch(chol(V), error = function(e) NULL)
  env_phylo_need(!is.null(R), "Fixed covariance is not positive definite.")
  list(V = V, R = R)
}

env_phylo_fixed_names <- function(fit) {
  formula <- fit$split_formula[[1L]]$form_no_bars
  env_phylo_need(inherits(formula, "formula") && is.data.frame(fit$data) &&
    nrow(fit$data) > 0L, "Saved fixed design metadata are missing.")
  answer <- colnames(stats::model.matrix(formula,
    fit$data[seq_len(min(3L, nrow(fit$data))), , drop = FALSE]))
  env_phylo_need(length(answer) > 0L && !anyNA(answer) && !anyDuplicated(answer),
    "Saved fixed design column names are invalid.")
  answer
}

env_phylo_table <- function(beta, V) {
  se <- sqrt(diag(V))
  data.frame(term = names(beta), estimate = unname(beta), SE = unname(se),
    lower_95 = unname(beta - 1.96 * se), upper_95 = unname(beta + 1.96 * se),
    p_Wald = unname(2 * stats::pnorm(-abs(beta / se))), stringsAsFactors = FALSE)
}

env_phylo_state <- function(fit) {
  env_phylo_need(isFALSE(fit$reml) && identical(fit$family$family, "gaussian") &&
    identical(fit$family$link, "identity"), "Saved comparison requires Gaussian identity-link ML.")
  theta <- fit$model$par
  report <- fit$sd_report
  env_phylo_need(is.numeric(theta) && length(theta) > 0L && !is.null(names(theta)) &&
    !anyNA(names(theta)) && all(is.finite(theta)), "Saved outer parameters are invalid.")
  env_phylo_need(is.list(report) && identical(names(theta), names(report$par.fixed)) &&
    env_phylo_numeric_equal(theta, report$par.fixed),
    "Optimizer parameters and saved uncertainty report differ.")
  covariance <- env_phylo_covariance(report$cov.fixed, length(theta), names(theta))
  V <- covariance$V
  g <- fit$gradients
  env_phylo_need(is.numeric(g) && length(g) == length(theta) && all(is.finite(g)),
    "Stored outer gradient is missing, non-finite, or misaligned.")
  if (!is.null(names(g))) env_phylo_need(identical(names(g), names(theta)),
    "Stored gradient parameter names are misordered.")
  report_difference <- NA_real_
  if (!is.null(report$gradient.fixed)) {
    env_phylo_need(env_phylo_numeric_equal(g, report$gradient.fixed),
      "Stored gradient and uncertainty-report gradient differ.")
    if (!is.null(names(report$gradient.fixed))) env_phylo_need(
      identical(names(report$gradient.fixed), names(theta)),
      "Uncertainty-report gradient parameter names are misordered.")
    report_difference <- max(abs(g - report$gradient.fixed))
  }
  terms <- env_phylo_fixed_names(fit)
  ix <- which(sub("\\[.*$", "", names(theta)) == "b_j")
  env_phylo_need(length(ix) == length(terms), "Cannot align the saved fixed coefficients.")
  beta <- stats::setNames(unname(theta[ix]), terms)
  Vf <- V[ix, ix, drop = FALSE]
  dimnames(Vf) <- list(terms, terms)
  objective <- fit$model$objective
  env_phylo_need(is.numeric(objective) && length(objective) == 1L && is.finite(objective),
    "Saved marginal likelihood is invalid.")
  se <- sqrt(diag(V))
  decrement <- sqrt(sum(as.vector(covariance$R %*% g)^2))
  env_phylo_need(is.finite(decrement), "Covariance-scaled gradient is not finite.")
  list(theta = theta, V = V, g = g, beta = beta, Vf = Vf,
    fixed = env_phylo_table(beta, Vf), objective = objective,
    raw_gradient = max(abs(g)), scaled_gradient = decrement,
    max_score_times_SE = max(abs(g * se)),
    max_newton_shift_SE = max(abs(as.vector(V %*% g) / se)),
    gradient_report_difference = report_difference,
    theta_report_difference = max(abs(theta - report$par.fixed)))
}

env_phylo_limits <- function() {
  scaled <- if (exists("scaled_gradient_threshold", inherits = TRUE))
    get("scaled_gradient_threshold", inherits = TRUE) else 0.05
  raw <- if (exists("raw_gradient_threshold", inherits = TRUE))
    get("raw_gradient_threshold", inherits = TRUE) else 0.001
  env_phylo_need(length(scaled) == 1L && is.finite(scaled) && scaled > 0 &&
    length(raw) == 1L && is.finite(raw) && raw > 0, "Invalid numerical thresholds.")
  list(scaled = scaled, raw = raw)
}

env_phylo_checks <- function(fit, threshold = NULL) {
  limits <- env_phylo_limits()
  error <- ""
  state <- tryCatch(env_phylo_state(fit), error = function(e) {
    error <<- conditionMessage(e)
    NULL
  })
  code <- fit$model$convergence
  if (length(code) != 1L || !is.numeric(code) || !is.finite(code)) code <- NA_integer_
  pd <- isTRUE(fit$sd_report$pdHess)
  good <- !is.null(state) && isTRUE(code == 0L) && pd &&
    state$scaled_gradient < limits$scaled
  if (!nzchar(error) && !good) error <- paste(c(
    if (!isTRUE(code == 0L)) "Optimizer convergence code is not zero.",
    if (!pd) "Saved uncertainty report does not confirm a positive Hessian.",
    if (!is.null(state) && state$scaled_gradient >= limits$scaled)
      "Covariance-scaled gradient exceeds the prescribed threshold."), collapse = " ")
  number <- function(name) if (is.null(state)) NA_real_ else state[[name]]
  data.frame(convergence_code = code, positive_definite_Hessian = pd, pdHess = pd,
    max_abs_gradient = number("raw_gradient"), raw_gradient_threshold = limits$raw,
    raw_gradient_advisory_passed = !is.null(state) && state$raw_gradient < limits$raw,
    gradient_covariance_norm = number("scaled_gradient"),
    scaled_gradient_threshold = limits$scaled,
    max_abs_score_times_SE = number("max_score_times_SE"),
    max_abs_approx_newton_shift_SE = number("max_newton_shift_SE"),
    gradient_report_max_difference = number("gradient_report_difference"),
    theta_report_max_difference = number("theta_report_difference"),
    finite_fixed_SE = if (is.null(state)) NA else TRUE,
    covariance_positive_definite = if (is.null(state)) NA else TRUE,
    fit_ok = good, failure_reason = error, native_evaluations_run = FALSE,
    check_scope = "Saved-state numerical checks; raw gradient advisory; no global optimum or biological adequacy guarantee",
    stringsAsFactors = FALSE)
}

env_phylo_fixed_table <- function(fit) {
  env_phylo_state(fit)$fixed
}

env_phylo_source_fixed <- function(source) {
  tab <- source$fixed
  env_phylo_need(is.data.frame(tab) && all(c("term", "estimate", "SE") %in% names(tab)) &&
    nrow(tab) > 0L && !anyNA(tab$term) && !anyDuplicated(tab$term) &&
    all(is.finite(tab$estimate)) && all(is.finite(tab$SE) & tab$SE > 0),
    "Independent-source fixed estimates are invalid.")
  env_phylo_need(!is.null(dimnames(source$V)),
    "Independent-source fixed covariance must carry the fixed term names.")
  V <- env_phylo_covariance(source$V, nrow(tab), as.character(tab$term))$V
  env_phylo_need(env_phylo_numeric_equal(tab$SE, sqrt(diag(V))),
    "Independent-source standard errors and fixed covariance differ.")
  beta <- stats::setNames(tab$estimate, as.character(tab$term))
  list(beta = beta, V = V, fixed = env_phylo_table(beta, V))
}

env_phylo_interaction <- function(terms, anomaly, moderator) {
  found <- intersect(c(paste(anomaly, moderator, sep = ":"),
    paste(moderator, anomaly, sep = ":")), terms)
  env_phylo_need(length(found) == 1L, paste("Cannot uniquely identify interaction:", moderator))
  found
}

env_phylo_biological <- function(estimate, lower, upper, orientation) {
  data.frame(estimate_biological = orientation * estimate,
    lower_95_biological = pmin(orientation * lower, orientation * upper),
    upper_95_biological = pmax(orientation * lower, orientation * upper))
}

env_phylo_contrast <- function(beta, V, within, between) {
  w <- stats::setNames(rep(0, length(beta)), names(beta))
  w[within] <- 1
  w[between] <- -1
  variance <- sum(as.vector(chol(V) %*% w)^2)
  env_phylo_need(is.finite(variance) && variance >= 0, "Invalid within-minus-between contrast variance.")
  estimate <- sum(w * beta)
  se <- sqrt(variance)
  data.frame(estimate = estimate, SE = se, lower_95 = estimate - 1.96 * se,
    upper_95 = estimate + 1.96 * se,
    p_Wald = if (se > 0) 2 * stats::pnorm(-abs(estimate / se)) else NA_real_)
}

env_phylo_wide_pair <- function(old, new) {
  names(old) <- paste0(names(old), "_independent")
  names(new) <- paste0(names(new), "_phylogenetic")
  cbind(old, new)
}

env_phylo_payload <- function(saved, p, cfg) {
  fit <- saved$fit
  state <- env_phylo_state(fit)
  original <- env_phylo_source_fixed(p$source)
  env_phylo_need(identical(names(original$beta), names(state$beta)),
    "Independent and phylogenetic fixed coefficient names/order differ.")
  expected <- if (identical(p$analysis, "onset")) 18L else 19L
  env_phylo_need(p$analysis %in% c("onset", "offset_univoltine", "offset_multivoltine") &&
    length(state$beta) == expected, "Unexpected environmental Mundlak fixed design.")
  env_phylo_need(is.list(p$spec) && length(p$spec$orientation) == 1L &&
    p$spec$orientation %in% c(-1, 1) && is.character(p$spec$anomaly) &&
    length(p$spec$anomaly) == 1L, "Environmental analysis specification is invalid.")
  window <- if (p$analysis == "onset") 60L else 90L
  orientation <- if (p$analysis == "offset_multivoltine") 1 else -1
  env_phylo_need(identical(p$spec$anomaly, paste0("clim_anomaly_tw", window)) &&
    length(p$spec$window) == 1L && p$spec$window == window &&
    p$spec$orientation == orientation,
    "Analysis, thermal window, anomaly, or biological orientation are inconsistent.")
  if (!is.null(saved$signature) && !is.null(p$signature)) env_phylo_need(
    identical(saved$signature, p$signature), "Saved phylogenetic fit signature differs from prepared input.")
  source_checks <- p$source$checks
  source_flags <- intersect(c("fit_ok", "numerical_checks_ok"), names(source_checks))
  source_good <- is.data.frame(source_checks) && nrow(source_checks) == 1L &&
    length(source_flags) > 0L && all(vapply(source_flags,
      function(name) isTRUE(source_checks[[name]]), logical(1)))
  env_phylo_need(source_good, "Independent source has not passed its numerical checks.")
  ss <- p$source$stats
  env_phylo_need(is.data.frame(ss) && nrow(ss) == 1L &&
    all(c("n_observations", "n_parameters", "logLik", "AIC") %in% names(ss)) &&
    all(is.finite(unlist(ss[c("n_observations", "n_parameters", "logLik", "AIC")]))),
    "Independent-source likelihood metadata are invalid.")
  env_phylo_need(ss$n_observations == nrow(fit$data) &&
    length(state$theta) == ss$n_parameters + 1L,
    "Compared models must retain the observations and add exactly one outer parameter.")
  env_phylo_need(abs(ss$AIC - (-2 * ss$logLik + 2 * ss$n_parameters)) < 1e-5,
    "Independent-source AIC and likelihood metadata disagree.")
  tolerance <- if (is.null(cfg$likelihood_tolerance)) 1e-4 else cfg$likelihood_tolerance
  env_phylo_need(length(tolerance) == 1L && is.finite(tolerance) && tolerance >= 0,
    "Invalid nested-likelihood comparison tolerance.")
  env_phylo_need(state$objective <= -ss$logLik + tolerance,
    "Phylogenetic marginal likelihood is worse than the nested independent source.")
  fixed <- data.frame(analysis = p$analysis, term = names(state$beta),
    env_phylo_wide_pair(original$fixed[-1L], state$fixed[-1L]),
    change_in_independent_SE = (state$beta - original$beta) / original$fixed$SE,
    SE_ratio = state$fixed$SE / original$fixed$SE, row.names = NULL)
  labels <- c(photoperiod = "Photoperiod", background = "Background temperature",
    predictability = "Temperature predictability", trend = "Warming trend")
  scaling <- p$component_scaling
  env_phylo_need(is.data.frame(scaling) && all(c("moderator", "within_SD",
    "between_population_weighted_SD") %in% names(scaling)) &&
    !anyDuplicated(scaling$moderator) && all(names(labels) %in% scaling$moderator),
    "Environmental component scaling metadata are missing or duplicated.")
  moderator_rows <- list()
  difference_rows <- list()
  for (key in names(labels)) {
    terms <- c(within_species = env_phylo_interaction(names(state$beta), p$spec$anomaly,
      paste0(key, "_within_species")), between_species = env_phylo_interaction(
      names(state$beta), p$spec$anomaly, paste0(key, "_species_mean")))
    sc <- scaling[match(key, scaling$moderator), , drop = FALSE]
    for (scale_name in names(terms)) {
      term <- terms[[scale_name]]
      old <- original$fixed[match(term, original$fixed$term), -1L, drop = FALSE]
      new <- state$fixed[match(term, state$fixed$term), -1L, drop = FALSE]
      sd <- if (scale_name == "within_species") sc$within_SD else sc$between_population_weighted_SD
      env_phylo_need(length(sd) == 1L && is.finite(sd) && sd >= 0,
        "Environmental component SD is invalid.")
      ob <- env_phylo_biological(old$estimate, old$lower_95, old$upper_95, p$spec$orientation)
      nb <- env_phylo_biological(new$estimate, new$lower_95, new$upper_95, p$spec$orientation)
      pair <- env_phylo_wide_pair(cbind(old, ob), cbind(new, nb))
      moderator_rows[[length(moderator_rows) + 1L]] <- data.frame(
        analysis = p$analysis, moderator = labels[[key]], moderator_key = key,
        effect_scale = if (scale_name == "within_species") "within species" else "between species",
        term = term, component_SD = sd, pair,
        biological_effect_per_component_SD_independent = sd * ob$estimate_biological,
        biological_effect_per_component_SD_phylogenetic = sd * nb$estimate_biological,
        change_in_independent_SE = (new$estimate - old$estimate) / old$SE,
        SE_ratio = new$SE / old$SE, row.names = NULL)
    }
    old <- env_phylo_contrast(original$beta, original$V, terms[[1L]], terms[[2L]])
    new <- env_phylo_contrast(state$beta, state$Vf, terms[[1L]], terms[[2L]])
    ob <- env_phylo_biological(old$estimate, old$lower_95, old$upper_95, p$spec$orientation)
    nb <- env_phylo_biological(new$estimate, new$lower_95, new$upper_95, p$spec$orientation)
    difference_rows[[length(difference_rows) + 1L]] <- data.frame(analysis = p$analysis,
      moderator = labels[[key]], contrast = "within minus between",
      env_phylo_wide_pair(cbind(old, ob), cbind(new, nb)), row.names = NULL)
  }
  iid <- slope_logsd_index(fit)
  phylo <- which(sub("\\[.*$", "", names(state$theta)) == "ln_smooth_sigma")
  env_phylo_need(length(iid) == 1L && length(phylo) == 1L && iid != phylo &&
    iid %in% seq_along(state$theta), "Cannot identify the iid and phylogenetic slope SDs.")
  ids <- c(iid, phylo)
  logsd <- unname(state$theta[ids])
  logse <- sqrt(diag(state$V)[ids])
  sd <- exp(logsd)
  env_phylo_need(all(is.finite(sd)) && all(is.finite(sd^2)), "Slope variance estimates overflow.")
  variance <- data.frame(analysis = p$analysis, component = c("independent_species_slope",
    "phylogenetic_species_slope"), log_SD = logsd, log_SD_SE = logse, SD = sd,
    SD_lower_95 = exp(logsd - 1.96 * logse), SD_upper_95 = exp(logsd + 1.96 * logse),
    variance = sd^2, interval_method = "Approximate normal-Wald interval on log SD; variance boundary caveat")
  logit <- 2 * (logsd[2L] - logsd[1L])
  logit_gradient <- rep(0, length(state$theta))
  logit_gradient[phylo] <- 2
  logit_gradient[iid] <- -2
  logit_variance <- sum(as.vector(chol(state$V) %*% logit_gradient)^2)
  env_phylo_need(is.finite(logit_variance) && logit_variance >= 0,
    "Phylogenetic variance-fraction uncertainty is invalid.")
  logit_se <- sqrt(logit_variance)
  fraction <- stats::plogis(logit)
  fraction_table <- data.frame(analysis = p$analysis, phylogenetic_fraction = fraction,
    lower_95 = stats::plogis(logit - 1.96 * logit_se),
    upper_95 = stats::plogis(logit + 1.96 * logit_se),
    near_variance_boundary_descriptive = fraction < 0.01 || fraction > 0.99,
    interpretation = "Conditional between-species anomaly-slope variance fraction; not total phenological variance",
    interval_method = "Joint fixed-covariance delta method on logit fraction; approximate near boundaries")
  ll <- -state$objective
  k <- length(state$theta)
  comparison <- data.frame(analysis = p$analysis, model = c("independent", "independent_plus_phylogeny"),
    n_observations = c(ss$n_observations, nrow(fit$data)), n_parameters = c(ss$n_parameters, k),
    logLik = c(ss$logLik, ll), AIC = c(ss$AIC, -2 * ll + 2 * k),
    delta_AIC_from_independent = c(0, -2 * ll + 2 * k - ss$AIC),
    LR_vs_independent_descriptive = c(0, 2 * (ll - ss$logLik)),
    boundary_LRT_p_computed = FALSE,
    interpretation = "Gaussian ML comparison on unchanged observations; no automatic variance-boundary significance test")
  list(fixed_effect_comparison.csv = fixed,
    moderator_comparison.csv = do.call(rbind, moderator_rows),
    within_minus_between_comparison.csv = do.call(rbind, difference_rows),
    variance_components.csv = variance, phylogenetic_fraction.csv = fraction_table,
    model_comparison.csv = comparison)
}

env_phylo_export_results <- function(saved, p, cfg) {
  out <- p$out
  env_phylo_need(is.character(out) && length(out) == 1L && nzchar(out),
    "Phylogenetic output directory is missing.")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  checks <- env_phylo_checks(saved$fit)
  failure <- checks$failure_reason
  payload <- if (isTRUE(checks$fit_ok)) tryCatch(env_phylo_payload(saved, p, cfg),
    error = function(e) {failure <<- conditionMessage(e); NULL}) else NULL
  good <- !is.null(payload)
  paths <- c("fixed_effect_comparison.csv", "moderator_comparison.csv",
    "within_minus_between_comparison.csv", "variance_components.csv",
    "phylogenetic_fraction.csv", "model_comparison.csv")
  if (good) {
    write_csv(data.frame(analysis = p$analysis, status = "EXPORTING",
      inference_exports_written = FALSE, native_evaluations_run = FALSE, reason = ""),
      file.path(out, "export_status.csv"))
    tryCatch(for (name in names(payload)) write_csv(payload[[name]], file.path(out, name)),
      error = function(e) {
        good <<- FALSE
        failure <<- paste("Comparison CSV writing failed:", conditionMessage(e))
      })
  }
  checks$comparison_checks_ok <- good
  checks$fit_ok <- isTRUE(checks$fit_ok) && good
  checks$failure_reason <- failure
  if (!good) {
    old <- file.path(out, paths)
    old <- old[file.exists(old)]
    if (length(old)) env_phylo_need(all(vapply(old, function(x) unlink(x) == 0L, logical(1))),
      "Cannot remove stale comparison exports after a failed check.")
  }
  write_csv(checks, file.path(out, "fit_checks.csv"))
  write_csv(data.frame(analysis = p$analysis, status = if (good) "EXPORTED" else "NUMERICAL_OR_INPUT_REVIEW_REQUIRED",
    inference_exports_written = good, native_evaluations_run = FALSE, reason = failure,
    note = "Coefficient changes divided by original SE are descriptive, not tests between dependent fits"),
    file.path(out, "export_status.csv"))
  writeLines(c(if (good) "Comparison exports passed the prescribed saved-state numerical and input checks."
      else "Comparison exports were suppressed: numerical or input review is required.",
    paste("Failure reason:", if (nzchar(failure)) failure else "none"),
    paste("Covariance-scaled gradient norm:", checks$gradient_covariance_norm,
      "threshold:", checks$scaled_gradient_threshold),
    paste("Raw maximum gradient:", checks$max_abs_gradient,
      "advisory threshold:", checks$raw_gradient_threshold,
      "advisory passed:", checks$raw_gradient_advisory_passed),
    "The raw gradient advisory and the covariance-scaled passing criterion are reported separately.",
    "These exports use the fixed-parameter covariance; they do not establish availability of ADREPORT/random-effect standard errors.",
    "No native likelihood, gradient, Hessian, optimizer, or package likelihood method was evaluated by this export.",
    "Passing these local numerical checks does not prove a global optimum or biological/spatial model adequacy.",
    "No automatic variance-boundary likelihood-ratio p-value is computed."),
    file.path(out, "inference_status.txt"))
  invisible(good)
}


  VERSION <- "environment_mundlak_phylogeny_15km_1.0.0"
  ROOT <- "E:/phenoIMPACT project/code/phenoIMPACT"
  ANALYSES <- c("onset", "offset_univoltine", "offset_multivoltine")

  need <- function(ok, message) {
    if (!isTRUE(ok)) stop(message, call. = FALSE)
  }

  env_phylo_preflight_mundlak <- function(directory) {
    progress("Checking all three Mundlak formulas against exact Gaussian likelihoods")
    A <- kronecker(diag(3L), matrix(c(1, .65, .65, 1), 2L))
    B <- matrix(0, 6L, 9L)
    for (j in 1:3) B[(2L*j-1L):(2L*j), j] <- sqrt(.65)
    B[cbind(1:6, 4:9)] <- sqrt(.35)
    attr(A, "edge_basis") <- B
    need(max(abs(tcrossprod(B) - A)) < 1e-12, "Invalid preflight increment covariance.")
    result <- list()
    for (response in c("onset", "offset_univoltine", "offset_multivoltine")) {
      d <- expand.grid(YEAR = 2018:2021, population = 1:4, PHYLO_ID = 1:6)
      d$g <- factor(d$PHYLO_ID)
      d$x <- sin(seq_len(nrow(d)) * .73) + .2 * (d$YEAR - 2019)
      components <- character()
      for (k in seq_along(c("photoperiod", "background", "predictability", "trend"))) {
        nm <- c("photoperiod", "background", "predictability", "trend")[[k]]
        context <- sin((d$PHYLO_ID * 5 + d$population) * (.13 + .07*k)) +
          cos(d$PHYLO_ID * (.21 + .11*k))
        mean_context <- ave(context, d$PHYLO_ID, FUN = mean)
        d[[paste0(nm, "_within_species")]] <- context - mean_context
        d[[paste0(nm, "_species_mean")]] <- mean_context
        components <- c(components, paste0(nm, c("_within_species", "_species_mean")))
      }
      d$ONSET_mean_z <- cos(seq_len(nrow(d)) * .37)
      d$y <- sin(seq_len(nrow(d)) * .7)
      fixed <- stats::as.formula(paste0("y ~ x * (", paste(components, collapse = " + "), ")",
        if (response != "onset") " + ONSET_mean_z" else ""), env = baseenv())
      X <- stats::model.matrix(fixed, d)
      need(qr(X)$rank == ncol(X) && ncol(X) == if (response == "onset") 18L else 19L,
           "Mundlak preflight design changed.")
      f <- add_phylo(stats::update(fixed, . ~ . + (0 + x | g)), A, "x")
      fit <- sdmTMB::sdmTMB(f, data = d, family = stats::gaussian(), spatial = "off",
        spatiotemporal = "off", reml = FALSE, do_fit = FALSE, silent = TRUE,
        control = sdmTMB::sdmTMBcontrol(get_joint_precision = FALSE, newton_loops = 0L))
      gaussian_inner(fit$tmb_obj)
      tryCatch({
        for (setting in 1:2) {
          beta <- if (setting == 1L) rep(0, ncol(X)) else seq_len(ncol(X)) * .013 - .1
          sig <- if (setting == 1L) c(1, 1, 1) else c(1.2, .4, .85)
          theta <- fit$tmb_obj$par
          values <- list(b_j = beta, ln_phi = log(sig[1]), re_cov_pars = log(sig[2]),
                         ln_smooth_sigma = log(sig[3]))
          for (nm in names(values)) {
            ids <- par_indices(theta, nm)
            need(length(ids) == length(values[[nm]]), "Mundlak preflight parameter mapping changed.")
            theta[ids] <- values[[nm]]
          }
          actual <- fit$tmb_obj$fn(theta)
          gradient <- fit$tmb_obj$gr(theta)
          id <- d$PHYLO_ID
          covariance <- sig[1]^2 * diag(nrow(d)) + outer(d$x, d$x) *
            (sig[2]^2 * outer(id, id, "==") + sig[3]^2 * A[id, id])
          R <- chol(covariance)
          residual <- as.numeric(d$y - X %*% beta)
          z <- forwardsolve(t(R), residual)
          exact <- .5 * (nrow(d) * log(2*pi) + 2*sum(log(diag(R))) + sum(z*z))
          fixed_gradient <- -as.numeric(crossprod(X, solve(covariance, residual)))
          result[[length(result)+1L]] <- data.frame(analysis = response, setting = setting,
            fixed_coefficients = ncol(X), increment_columns = ncol(B),
            nll_error = abs(actual-exact),
            fixed_gradient_error = max(abs(gradient[par_indices(theta, "b_j")]-fixed_gradient)),
            gradient_finite = all(is.finite(gradient)))
        }
      }, finally = release_tmb(fit$tmb_obj))
    }
    result <- do.call(rbind, result)
    write_csv(result, file.path(directory, "preflight_mundlak_likelihood.csv"))
    # Legacy worker copies this filename into the result directory.
    write_csv(result, file.path(directory, "preflight_onset_likelihood.csv"))
    need(all(result$nll_error < 1e-6 & result$fixed_gradient_error < 1e-6 & result$gradient_finite),
         "Native Mundlak likelihood validation failed; no large model attempted.")
    invisible(result)
  }

  # Reuse the tested v2 engine without executing its hard-coded old-model launch.
  # Cut only at the unique final configuration assignment, after its definitions.
  make_engine <- function(root, analysis, scaled_tolerance, raw_tolerance) {
    private <- new.env(parent = globalenv())
    path <- file.path(root, "pheno_phylogeny_qc_15km_v2.R")
    need(file.exists(path), paste("Missing existing phylogenetic engine:", path))
    sys.source(path, envir = private)
    factory <- private$run_phylogeny_qc_15km_v2
    need(is.function(factory), "Expected the existing v2 phylogenetic engine.")
    statements <- as.list(body(factory))[-1L]
    cutoff <- which(vapply(statements, function(x) {
      is.call(x) && identical(x[[1L]], as.name("<-")) &&
        identical(x[[2L]], as.name("cfg")) &&
        identical(x[[3L]], quote(private$phylogeny_pilot_config))
    }, logical(1)))
    need(length(cutoff) == 1L && cutoff > 1L,
         "The v2 engine structure changed; refusing an unreviewed adapter.")
    body(factory) <- as.call(c(list(as.name("{")),
      statements[seq_len(cutoff - 1L)],
      list(quote(list(engine = engine, defaults = private$phylogeny_pilot_config)))))
    template <- suppressMessages(factory(project_root = root, minimum_free_GiB = 0))
    e <- environment(template$engine$start)
    need(identical(e$version, "phylo_qc_multivoltine_edges_15km_1.1.0"),
         "Unexpected v2 engine version.")
    helpers <- ls(environment(make_engine), pattern = "^env_phylo_")
    for (nm in helpers) {
      f <- get(nm, envir = environment(make_engine))
      if (is.function(f)) {
        environment(f) <- e
        assign(nm, f, envir = e)
      }
    }
    e$version <- VERSION
    e$anomaly <- paste0("clim_anomaly_tw", if (analysis == "onset") "60" else "90")
    e$moderators <- as.vector(rbind(
      paste0(c("photoperiod", "background", "predictability", "trend"), "_within_species"),
      paste0(c("photoperiod", "background", "predictability", "trend"), "_species_mean")))
    e$scaled_gradient_threshold <- scaled_tolerance
    e$raw_gradient_threshold <- raw_tolerance
    e$fixed_table <- e$env_phylo_fixed_table
    e$checks <- e$env_phylo_checks
    e$prepare <- e$env_phylo_prepare
    e$export_results <- e$env_phylo_export_results
    e$preflight_onset <- e$env_phylo_preflight_mundlak

    # The first native construction can fail before an outer iteration exists.
    # Keep a portable phase marker and release only this worker's lock on exit.
    e$worker_base <- e$worker
    e$worker <- function(cfg) {
      on.exit(unlink(cfg$launch_lock, recursive = TRUE), add = TRUE)
      worker_base(cfg)
    }
    environment(e$worker) <- e
    list(engine = template$engine, defaults = template$defaults, environment = e)
  }

  absolute_path <- function(path, root) {
    if (!grepl("^([A-Za-z]:[/\\\\]|/|\\\\\\\\)", path)) path <- file.path(root, path)
    normalizePath(path, winslash = "/", mustWork = TRUE)
  }

  process_matches <- function(pid, created) {
    tryCatch({
      if (length(pid) != 1L || !is.numeric(pid) || !is.finite(pid) || pid <= 0) return(NA)
      if (!pid %in% ps::ps_pids()) return(FALSE)
      if (length(created) != 1L || !is.numeric(created) || !is.finite(created)) return(NA)
      h <- ps::ps_handle(pid)
      isTRUE(abs(as.numeric(ps::ps_create_time(h)) - created) < 0.01) &&
        isTRUE(ps::ps_is_running(h)) && !ps::ps_status(h) %in% c("zombie", "dead")
    }, error = function(e) NA)
  }

  acquire_lock <- function(path) {
    if (dir.exists(path)) {
      owner <- tryCatch(readRDS(file.path(path, "owner.rds")), error = function(e) NULL)
      need(!is.null(owner), paste("A launch lock already exists; inspect:", path))
      active <- process_matches(owner$pid, owner$created)
      need(identical(active, FALSE), paste("Another phylogenetic worker may still be active:", path))
      unlink(path, recursive = TRUE)
    }
    need(dir.create(path), paste("Cannot acquire phylogenetic worker lock:", path))
    owner <- list(pid = Sys.getpid(), created = as.numeric(ps::ps_create_time(ps::ps_handle())))
    saveRDS(owner, file.path(path, "owner.rds"))
    invisible(path)
  }

  run <- function(
    project_root = ROOT, analysis = "offset_multivoltine", source_run = NULL,
    tree_file = "data/EUROPEAN_BUTTERFLIES_MCC_CA.nwk", check_only = FALSE,
    checkpoint_seconds = 120, minimum_free_GiB = 30, maximum_single_basis_GiB = 2,
    extra_optimization_attempts = 0L, scaled_gradient_tolerance = 0.05,
    raw_gradient_tolerance = 0.001
  ) {
    need(length(analysis) == 1L && analysis %in% ANALYSES, "Choose one supported analysis.")
    need(is.logical(check_only) && length(check_only) == 1L && !is.na(check_only), "Invalid check_only.")
    for (x in list(checkpoint_seconds, minimum_free_GiB, maximum_single_basis_GiB,
                   scaled_gradient_tolerance, raw_gradient_tolerance)) {
      need(is.numeric(x) && length(x) == 1L && is.finite(x) && x > 0, "Invalid numeric setting.")
    }
    need(is.numeric(extra_optimization_attempts) && length(extra_optimization_attempts) == 1L &&
      is.finite(extra_optimization_attempts) && extra_optimization_attempts >= 0 &&
      extra_optimization_attempts == as.integer(extra_optimization_attempts), "Invalid refinement count.")
    need(requireNamespace("ps", quietly = TRUE), "The existing package ps is required.")
    root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
    tree <- absolute_path(tree_file, root)
    if (!is.null(source_run)) source_run <- absolute_path(source_run, root)
    resolved <- env_phylo_resolve_source(root, analysis, source_run)
    free <- ps::ps_system_memory()$avail / 1024^3
    message("RAM disponible: ", round(free, 1), " GiB")
    if (!check_only) need(free >= minimum_free_GiB, "Insufficient available RAM for the full fit.")
    for (option in c("pheno.env.mundlak.phylo.job", "pheno.demo.mundlak.phylo.job",
                     "pheno.demography.mundlak.phylo.job", "pheno.phylo.qc.job",
                     "pheno.phylo.pilot.job", "phenoimpact.env.mundlak.job", "phenoimpact.mundlak.qc.job")) {
      previous <- getOption(option)
      if (!is.null(previous)) {
        active <- tryCatch(previous$process$is_alive(), error = function(e) FALSE)
        need(!isTRUE(active), paste("A model worker is already active:", option))
      }
    }
    template <- make_engine(root, analysis, scaled_gradient_tolerance, raw_gradient_tolerance)
    cfg <- utils::modifyList(template$defaults, resolved)
    cfg$project_root <- root
    cfg$analysis <- analysis
    cfg$tree_file <- tree
    cfg$output_root <- file.path(root, "output", "phenology_plasticity", "environment_mundlak_phylogeny_15km")
    dir.create(cfg$output_root, recursive = TRUE, showWarnings = FALSE)
    cfg$output_root <- normalizePath(cfg$output_root, winslash = "/", mustWork = TRUE)
    cfg$check_only <- check_only
    cfg$diagnostic_only <- FALSE
    cfg$fit_seed <- 202610150L + match(analysis, ANALYSES)
    cfg$checkpoint_seconds <- checkpoint_seconds
    cfg$maximum_single_basis_GiB <- maximum_single_basis_GiB
    cfg$extra_optimization_attempts <- as.integer(extra_optimization_attempts)
    cfg$gradient_threshold <- raw_gradient_tolerance
    cfg$scaled_gradient_threshold <- scaled_gradient_tolerance
    cfg$signature_files <- file.path(root, c("pheno_phylogeny_qc_15km_v2.R",
      "pheno_plasticity_phylogeny_pilot_15km.R"))
    # Native engine changes invalidate fitted checkpoints. Reporting-only edits
    # of this adapter invalidate its prepared cache, but preserve the fit.
    cfg$signature_file_hashes <- unname(tools::md5sum(cfg$signature_files))
    adapter_file <- file.path(root, "11c_environment_mundlak_phylogeny_15km.R")
    if (file.exists(adapter_file)) cfg$signature_files <- c(cfg$signature_files, adapter_file)
    cfg$launch_lock <- file.path(cfg$output_root, "worker.lock")
    acquire_lock(cfg$launch_lock)
    launched <- FALSE
    on.exit(if (!launched) unlink(cfg$launch_lock, recursive = TRUE), add = TRUE)
    job <- suppressMessages(template$engine$start(cfg))
    launched <- TRUE
    job$engine <- template$engine
    job$config <- cfg
    owner <- list(pid = job$process$get_pid(), created = NA_real_)
    if (job$process$is_alive()) {
      owner$created <- tryCatch(as.numeric(ps::ps_create_time(ps::ps_handle(owner$pid))),
                               error = function(e) NA_real_)
    }
    if (dir.exists(cfg$launch_lock)) try(saveRDS(owner, file.path(cfg$launch_lock, "owner.rds")), silent = TRUE)
    saveRDS(cfg, file.path(job$monitor_directory, "launch_config.rds"))
    options(pheno.env.mundlak.phylo.job = job)
    message(if (check_only) "Comprovacio previa iniciada; cap ajustament complet." else
      "Un model Mundlak ambiental filogenetic de 15 km iniciat.",
      "\nResposta: ", analysis, "\nModel independent: ", cfg$source_checkpoint,
      "\nSeguiment: watch_environment_mundlak_phylogeny()",
      "\nLog: ", job$log_file)
    invisible(job)
  }

  status <- function(job = getOption("pheno.env.mundlak.phylo.job"), quiet = FALSE) {
    need(!is.null(job), "First run run_environment_mundlak_phylogeny().")
    job$engine$status(job, quiet = quiet)
  }

  watch <- function(job = getOption("pheno.env.mundlak.phylo.job"), interval_seconds = 60) {
    need(!is.null(job), "First run run_environment_mundlak_phylogeny().")
    need(is.numeric(interval_seconds) && length(interval_seconds) == 1L &&
      is.finite(interval_seconds) && interval_seconds >= 1, "Invalid monitoring interval.")
    tryCatch({
      repeat {
        status(job)
        if (!job$process$is_alive()) break
        until <- Sys.time() + interval_seconds
        while (Sys.time() < until && job$process$is_alive()) Sys.sleep(1)
      }
      result <- tryCatch(job$process$get_result(), error = function(e) {
        message("Worker aturat. Checkpoints conservats. Log: ", job$log_file,
                "\n", conditionMessage(e))
        NULL
      })
      invisible(result)
    }, interrupt = function(e) {
      message("Monitor pausat; el worker continua. Torna a executar watch_environment_mundlak_phylogeny().")
      invisible(job)
    })
  }

  stop_job <- function(job = getOption("pheno.env.mundlak.phylo.job")) {
    if (is.null(job)) return(invisible(FALSE))
    if (job$process$is_alive()) job$process$kill_tree()
    if (!job$process$is_alive()) unlink(job$config$launch_lock, recursive = TRUE)
    message("Worker aturat; checkpoints conservats. Repeteix la mateixa crida run per reprendre.")
    invisible(TRUE)
  }

  audit_sources <- function(project_root = ROOT, analysis = NULL, source_run = NULL) {
    root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
    if (!is.null(source_run)) source_run <- absolute_path(source_run, root)
    table <- env_phylo_audit_sources(root, analysis = analysis, source_run = source_run)
    columns <- intersect(c("analysis", "source_directory", "version", "status",
      "numerical_checks_ok", "accepted", "reason"), names(table))
    print(table[columns], row.names = FALSE)
    message("Auditoria de metadades: no s'ha obert cap model ni s'ha iniciat cap ajustament.")
    invisible(table)
  }

  list(run = run, status = status, watch = watch, stop = stop_job, audit_sources = audit_sources)

})

run_environment_mundlak_phylogeny <- function(...) pheno_env_mundlak_phylogeny$run(...)
audit_environment_mundlak_phylogeny_sources <- function(...) pheno_env_mundlak_phylogeny$audit_sources(...)
resume_environment_mundlak_phylogeny <- function(...) pheno_env_mundlak_phylogeny$run(...)
status_environment_mundlak_phylogeny <- function(...) pheno_env_mundlak_phylogeny$status(...)
watch_environment_mundlak_phylogeny <- function(...) pheno_env_mundlak_phylogeny$watch(...)
stop_environment_mundlak_phylogeny <- function(...) pheno_env_mundlak_phylogeny$stop(...)

message("Funcions carregades. Primer: job <- run_environment_mundlak_phylogeny(check_only = TRUE). Seguiment: watch_environment_mundlak_phylogeny(job).")
