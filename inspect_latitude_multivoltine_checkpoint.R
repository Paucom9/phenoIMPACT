# phenoIMPACT | Inspect ONE saved latitude checkpoint, without refitting.
# The original fit_full.rds, its .previous file and all other fits are read-only.
# Sourcing defines the function only. To run:
#   inspection <- inspect_latitude_multivoltine_checkpoint()
# Requires the existing callr and ps installations. No packages are installed.

inspect_latitude_multivoltine_checkpoint <- function(
  root = paste0("E:/phenoIMPACT project/code/phenoIMPACT/",
                "output/phenology_plasticity/latitude_models_15km"),
  reserve_gib = 8,
  poll_seconds = 1
) {
  for (pkg in c("callr", "ps")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package missing: ", pkg, ". No packages were installed.")
    }
  }
  if (length(reserve_gib) != 1L || !is.finite(reserve_gib) || reserve_gib < 4)
    stop("reserve_gib must be a single finite number >= 4.")
  if (length(poll_seconds) != 1L || !is.finite(poll_seconds) || poll_seconds <= 0)
    stop("poll_seconds must be a single positive finite number.")

  path <- file.path(root, "primary__offset_multivoltine_f93f4498e835", "fit_full.rds")
  if (!file.exists(path)) stop("Checkpoint not found: ", path)
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  before <- file.info(path)
  if (isTRUE(before$isdir) || is.na(before$size) || before$size <= 0)
    stop("Checkpoint is not a nonempty file: ", path)

  # These are known project workers, not a comprehensive list of running jobs.
  for (key in c("pheno.latitude.active_job", "pheno.spatial.active_job")) {
    job <- getOption(key)
    alive <- tryCatch(!is.null(job) && isTRUE(job$process$is_alive()),
                      error = function(e) FALSE)
    if (alive) stop("Another known project worker is running: ", key)
  }
  available <- function() {
    a <- as.numeric(ps::ps_system_memory()[["avail"]]) / 1024^3
    if (length(a) != 1L || !is.finite(a)) stop("Cannot read available RAM.")
    a
  }
  free_start <- available()
  # A scheduling precaution, NOT an estimate of the memory needed to read this file.
  if (free_start < 2 * reserve_gib) {
    stop(sprintf("Only %.1f GiB RAM available. This inspection requires at least %.1f GiB free to start; that does NOT guarantee sufficient RAM.",
                 free_start, 2 * reserve_gib))
  }

  out <- tempfile(pattern = paste0("inspection_multivoltine_",
                                   format(Sys.time(), "%Y%m%d_%H%M%S"), "_"),
                  tmpdir = dirname(path))
  if (!dir.create(out)) stop("Cannot create the small-report directory: ", out)
  out <- normalizePath(out, winslash = "/", mustWork = TRUE)
  log <- file.path(out, "inspection.log")
  cat(sprintf("Available RAM: %.1f GiB. Inspecting only the multivoltine original.\n", free_start))
  cat("Reports:", out, "\n")
  cat("No refitting, model-method calls, checkpoint copying or overwriting.\n")

  worker <- callr::r_bg(function(path) {
    phase <- "readRDS"
    read_warnings <- character()
    tryCatch({
      cat("Reading original checkpoint; no model computations.\n")
      x <- withCallingHandlers(readRDS(path), warning = function(w) {
        read_warnings <<- c(read_warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
      phase <- "inspect saved fields"
      if (!is.list(x) || !all(c("fit", "checks") %in% names(x)))
        stop("Readable RDS, but not the expected latitude checkpoint structure.")
      f <- x[["fit"]]
      checks <- x[["checks"]]
      if (!is.list(f) || !inherits(f, "sdmTMB"))
        stop("The saved fit is not an sdmTMB list.")
      if (!is.data.frame(checks) || nrow(checks) > 10L || ncol(checks) > 50L)
        stop("Unexpected stored numerical-check table; no model inference attempted.")
      if (!is.data.frame(f[["data"]])) stop("Unexpected saved data structure.")
      d <- f[["data"]]
      g <- f[["gradients"]]
      fm <- f[["formula"]]
      scalar <- function(z) if (length(z) == 1L && is.atomic(z)) z else NA
      sd <- f[["sd_report"]]
      pf <- sd[["par.fixed"]]
      cv <- sd[["cov.fixed"]]
      # Cache only a small fixed-parameter covariance, never random/joint precision.
      cache_ok <- is.numeric(pf) && length(pf) > 0L && length(pf) <= 200L &&
        is.matrix(cv) && identical(dim(cv), c(length(pf), length(pf)))
      fixed_cache <- if (cache_ok) list(parameters = pf, covariance = cv) else NULL
      report <- list(
        status = "READABLE_EXPECTED_STRUCTURE",
        source = path,
        inspected_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
        r_version = R.version.string,
        label = scalar(x[["label"]]),
        run_id = scalar(x[["run_id"]]),
        fingerprint = scalar(x[["fingerprint"]]),
        fit_class = class(f),
        formula = if (inherits(fm, "formula")) paste(deparse(fm), collapse = " ") else NA_character_,
        n_observations = nrow(d),
        n_sites = if ("SITE_ID" %in% names(d)) length(unique(d[["SITE_ID"]])) else NA_integer_,
        n_species = if ("SPECIES" %in% names(d)) length(unique(d[["SPECIES"]])) else NA_integer_,
        checks_saved = checks,
        extra_optimization_saved = scalar(x[["extra_optimization"]]),
        elapsed_minutes_saved = scalar(x[["elapsed_minutes"]]),
        convergence_cached = scalar(f[["model"]][["convergence"]]),
        optimizer_message_cached = scalar(f[["model"]][["message"]]),
        objective_cached = scalar(f[["model"]][["objective"]]),
        positive_definite_Hessian_cached = scalar(sd[["pdHess"]]),
        max_abs_gradient_cached = if (is.numeric(g) && length(g) && all(is.finite(g))) max(abs(g)) else NA_real_,
        fixed_parameter_cache = fixed_cache,
        warnings_saved = if (is.character(x[["warnings"]])) x[["warnings"]] else character(),
        read_warnings = unique(read_warnings),
        note = "Stored diagnostics only. No reevaluation, refit, full validation or inferential export."
      )
      # Only this small report returns to the parent session, never x or f.
      cat("Checkpoint read; returning stored checks and small cached fields.\n")
      report
    }, error = function(e) {
      list(status = "INSPECTION_ERROR", phase = phase, source = path,
           message = conditionMessage(e), read_warnings = unique(read_warnings))
    })
  }, args = list(path = path), libpath = .libPaths(),
     stdout = log, stderr = "2>&1", supervise = TRUE,
     user_profile = FALSE, system_profile = FALSE)

  # A separate process does NOT isolate system RAM. Monitor it explicitly.
  on.exit({
    if (isTRUE(tryCatch(worker$is_alive(), error = function(e) FALSE))) worker$kill()
  }, add = TRUE)
  minimum_free <- free_start
  last_update <- Sys.time()
  while (worker$is_alive()) {
    free <- available()
    minimum_free <- min(minimum_free, free)
    if (free < reserve_gib) {
      worker$kill()
      msg <- sprintf("Inspection stopped: available RAM fell to %.1f GiB (reserve %.1f GiB). This does not establish file corruption. Original files were not changed.",
                     free, reserve_gib)
      writeLines(msg, file.path(out, "stopped_for_memory.txt"))
      stop(msg, "\nReport directory: ", out)
    }
    if (as.numeric(difftime(Sys.time(), last_update, units = "secs")) >= 30) {
      cat(sprintf("Still reading/inspecting; available RAM %.1f GiB.\n", free))
      last_update <- Sys.time()
    }
    Sys.sleep(poll_seconds)
  }
  report <- tryCatch(worker$get_result(), error = function(e) {
    list(status = "CHILD_PROCESS_ERROR", source = path,
         exit_status = worker$get_exit_status(), message = conditionMessage(e))
  })
  after <- file.info(path)
  report$source_unchanged_by_size_and_mtime <-
    isTRUE(before$size == after$size) &&
    isTRUE(as.numeric(before$mtime) == as.numeric(after$mtime))
  report$source_size_gib <- before$size / 1024^3
  report$minimum_available_ram_gib_observed <- minimum_free
  report$output_directory <- out
  if (!report$source_unchanged_by_size_and_mtime)
    report$status <- "SOURCE_CHANGED_DURING_INSPECTION"

  saveRDS(report, file.path(out, "inspection_summary.rds"))
  if (is.data.frame(report$checks_saved))
    utils::write.csv(report$checks_saved, file.path(out, "stored_checks.csv"), row.names = FALSE)
  brief <- report[setdiff(names(report), "fixed_parameter_cache")]
  text <- capture.output(print(brief))
  writeLines(text, file.path(out, "inspection_summary.txt"))
  cat("\n", paste(text, collapse = "\n"), "\n", sep = "")
  invisible(report)
}
