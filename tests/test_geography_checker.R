# Run from the repository root: Rscript --vanilla tests/test_geography_checker.R
# Small synthetic metadata only. No sdmTMB, model loading or worker launch.
source("geography_mundlak_by_voltinism_qc_final.R")
e <- environment(pheno_geo_mundlak$run)
checks <- 0L
check <- function(ok, msg) {
  if (!isTRUE(ok)) stop(msg, call. = FALSE)
  checks <<- checks + 1L
}
fails <- function(expr, pattern) {
  msg <- tryCatch({ force(expr); NA_character_ }, error = function(x) conditionMessage(x))
  check(!is.na(msg) && grepl(pattern, msg), paste("Expected", pattern, "got", msg))
}
write_checks <- function(path, gradient = .0001, ok = TRUE, convergence = 0L, hessian = TRUE) {
  write.csv(data.frame(numerical_checks_ok = ok, convergence_code = convergence,
    pdHess = hessian, max_abs_gradient = gradient), path, row.names = FALSE)
}
write_refinement <- function(path, analysis, status = "FIT_OK_REFINED") {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  # Intentionally invalid RDS: selection must never deserialize fitted objects.
  writeLines("MODEL MUST NOT BE LOADED", file.path(path, "fit_refined_final.rds"))
  saveRDS(list(summary = data.frame(analysis = analysis, refinement_status = status)),
    file.path(path, "refinement_completed.rds"))
  write_checks(file.path(path, "fit_checks_refined.csv"))
}

root <- tempfile("geo-checker-")
dir.create(root)
root <- normalizePath(root, winslash = "/", mustWork = TRUE)
env_run <- file.path(root, "output", "phenology_plasticity", "environment", "run_test")
dir.create(env_run, recursive = TRUE)
saveRDS(list(version = e$SOURCE_VERSION, signature = "source-signature"),
  file.path(env_run, "run_config.rds"))
for (a in e$ANALYSES) {
  d <- file.path(env_run, a)
  dir.create(d)
  saveRDS(list(status = "FIT_OK", signature = "source-signature"),
    file.path(d, "completed.rds"))
  writeLines("INPUT MUST NOT BE LOADED", file.path(d, "prepared_input.rds"))
  writeLines("MODEL MUST NOT BE LOADED", file.path(d, "fit_mundlak.rds"))
  write_checks(file.path(d, "fit_checks.csv"))
}
selected <- e$source_status(env_run, TRUE)$table
check(all(selected$ready) && all(selected$refinement_version == "none"),
  "Original accepted sources rejected")
original_prepared <- selected$prepared
first <- file.path(env_run, e$ANALYSES[1])
for (bad in list(
    list(gradient = .001), list(gradient = NA_real_), list(ok = FALSE),
    list(convergence = 1L), list(hessian = FALSE))) {
  do.call(write_checks, c(list(path = file.path(first, "fit_checks.csv")), bad))
  check(!e$source_status(env_run)$table$ready[1], "Invalid numerical checks accepted")
}
write_checks(file.path(first, "fit_checks.csv"))
saveRDS(list(status = "FIT_OK", signature = "wrong"), file.path(first, "completed.rds"))
fails(e$source_status(env_run), "signature mismatch")
saveRDS(list(status = "NUMERICAL_REVIEW_REQUIRED", signature = "source-signature"),
  file.path(first, "completed.rds"))
check(!e$source_status(env_run)$table$ready[1], "Unaccepted original passed")

v1 <- file.path(env_run, "postfit_refinement_queue_v1")
ref <- file.path(v1, "offset_multivoltine")
write_refinement(ref, "offset_multivoltine")
selected <- e$source_status(env_run, TRUE)$table
check(selected$ready[1] && selected$refinement_version[1] == "v1",
  "Accepted v1 did not supersede original review status")
check(identical(selected$prepared, original_prepared), "Prepared-input source changed")
write_refinement(ref, "offset_multivoltine", "NUMERICAL_REVIEW_REQUIRED_AFTER_MAX_PASSES")
check(!e$source_status(env_run)$table$ready[1], "Failed v1 accepted")
write_refinement(ref, "offset_multivoltine")
write_refinement(file.path(v1, "onset_multivoltine"), "onset_multivoltine")

v2 <- file.path(env_run, "postfit_refinement_queue_v2_isolated")
dir.create(v2)
selected <- e$source_status(env_run)$table
check(identical(selected$refinement_version, c("v1", "v2", "v2", "v2")),
  "Pending v2 did not supersede v1/original")
check(identical(selected$ready, c(TRUE, FALSE, FALSE, FALSE)), "Unstarted v2 accepted")
check(all(selected$status[-1] == "REFINEMENT_NOT_STARTED"), "Missing pending status")
fails(e$source_status(env_run, TRUE), "selected environmental sources")
for (a in e$V2_ANALYSES) write_refinement(file.path(v2, a), a)
selected <- e$source_status(env_run, TRUE)$table
check(all(selected$ready), "Mixed v1/v2 sources rejected")
check(selected$fit[1] == file.path(ref, "fit_refined_final.rds"), "Wrong carried-forward v1")
bad_v2 <- file.path(v2, "onset_multivoltine")
write_refinement(bad_v2, "onset_multivoltine", "NUMERICAL_REVIEW_REQUIRED_AFTER_MAX_PASSES")
check(!e$source_status(env_run)$table$ready[2], "Failed v2 fell back to accepted v1")
write_refinement(bad_v2, "onset_multivoltine", "FIT_OK_AFTER_REBUILD")
check(e$source_status(env_run)$table$ready[2], "Rebuild acceptance rejected")
for (missing in c("fit_refined_final.rds", "refinement_completed.rds", "fit_checks_refined.csv")) {
  unlink(file.path(bad_v2, missing))
  check(!e$source_status(env_run)$table$ready[2], paste("Accepted v2 without", missing))
  write_refinement(bad_v2, "onset_multivoltine")
}
write_refinement(bad_v2, "onset_univoltine")
fails(e$source_status(env_run), "identity mismatch")
write_refinement(bad_v2, "onset_multivoltine")
write_checks(file.path(bad_v2, "fit_checks_refined.csv"), gradient = .001)
check(!e$source_status(env_run)$table$ready[2], "Refined gradient threshold relaxed")
write_refinement(bad_v2, "onset_multivoltine")
saveRDS(list(status = "OK"), file.path(bad_v2, "isolated_worker_result.rds"))
unlink(file.path(bad_v2, "refinement_completed.rds"))
check(!e$source_status(env_run)$table$ready[2], "Worker success replaced numerical acceptance")
write_refinement(bad_v2, "onset_multivoltine")

# A configured v1 can be outside the scanned project output roots.
configured_v1 <- file.path(root, "external_v1")
write_refinement(file.path(configured_v1, "offset_multivoltine"), "offset_multivoltine")
v2_config <- file.path(v2, "refinement_config.rds")
saveRDS(list(source_run = env_run, old_refinement = configured_v1), v2_config)
selected <- e$source_status(env_run, TRUE)$table
check(selected$fit[1] == file.path(configured_v1, "offset_multivoltine", "fit_refined_final.rds"),
  "Configured v1 path ignored")
saveRDS(list(source_run = v2, old_refinement = configured_v1), v2_config)
fails(e$source_status(env_run), "different environmental run")
saveRDS(list(source_run = env_run, old_refinement = configured_v1), v2_config)
v1_config <- file.path(configured_v1, "refinement_config.rds")
saveRDS(list(source_run = v2), v1_config)
fails(e$source_status(env_run), "different environmental run")
saveRDS(list(source_run = env_run), v1_config)
cat("PASS: original/v1/v2 selection, precedence, provenance, missing artifacts, numerical gates\n")

# Simulate liveness and resources to exercise the public launch gate without ps,
# fitted objects, heavy dependencies, or background processes.
real_pid_alive <- e$pid_alive
real_resources <- e$available_resources
e$pid_alive <- function(pid, create_time) isTRUE(pid == 42L)
e$requireNamespace <- function(package, ...) TRUE
e$available_resources <- function(root)
  list(available_RAM_GiB = 100, available_disk_GiB = 100)
owner_paths <- c(file.path(env_run, "worker_owner.rds"),
  file.path(v1, "worker_owner.rds"), file.path(configured_v1, "worker_owner.rds"),
  file.path(v2, "worker_owner.rds"),
  file.path(v2, "offset_univoltine", "isolated_worker_owner.rds"))
for (path in owner_paths) {
  saveRDS(list(pid = 42L, create_time = 1), path)
  check(path %in% e$active_jobs(root, env_run), paste("Missed active owner", path))
  fails(e$require_environment_idle(env_run), "worker is active")
  saveRDS(list(pid = 43L, create_time = 1), path)
}
check(!length(e$active_jobs(root, env_run)), "Dead ownership blocks execution")
child_owner <- tail(owner_paths, 1)
saveRDS(list(pid = 42L, create_time = 1), child_owner)
invisible(capture.output(pf <- pheno_geo_mundlak$preflight(root, env_run)))
check(pf$sources_ready && !pf$ready && child_owner %in% pf$active_jobs,
  "Public preflight missed orphan isolated worker")
before_files <- sort(list.files(root, recursive = TRUE, all.files = TRUE))
fails(pheno_geo_mundlak$run(root, env_run), "Another known")
check(identical(before_files, sort(list.files(root, recursive = TRUE, all.files = TRUE))),
  "Blocked launch wrote files")
saveRDS(list(pid = 43L, create_time = 1), child_owner)
invisible(capture.output(pf <- pheno_geo_mundlak$preflight(root, env_run)))
check(pf$ready, "Accepted idle sources not ready with sufficient resources")
e$available_resources <- function(root)
  list(available_RAM_GiB = 10, available_disk_GiB = 100)
invisible(capture.output(pf <- pheno_geo_mundlak$preflight(root, env_run)))
check(pf$sources_ready && !pf$resources_ready && !pf$ready, "Low memory reported ready")
e$available_resources <- real_resources

# Workers recheck selection/metadata before touching inputs or loading fit packages.
src <- e$source_status(env_run, TRUE)
cfg <- list(environment_run = env_run, source_selection = src$table, source_files = src$files,
  source_manifest = e$file_manifest(c(src$files, src$table$fit)))
e$require_sources_current(cfg)
saveRDS(list(pid = 42L, create_time = 1), child_owner)
fails(e$fit_one(cfg, "offset_univoltine"), "worker is active")
saveRDS(list(pid = 43L, create_time = 1), child_owner)
write_refinement(bad_v2, "onset_multivoltine", "NUMERICAL_REVIEW_REQUIRED_AFTER_MAX_PASSES")
fails(e$require_sources_current(cfg), "selected environmental sources")
write_refinement(bad_v2, "onset_multivoltine")
src <- e$source_status(env_run, TRUE)
cfg$source_manifest <- e$file_manifest(c(src$files, src$table$fit))
writeLines("CHANGED MODEL MUST STILL NOT BE LOADED", src$table$fit[2])
fails(e$require_sources_current(cfg), "source changed")
writeLines("invalid RDS", child_owner)
fails(e$active_jobs(root, env_run), "unknown input format")
unlink(child_owner)
e$pid_alive <- real_pid_alive
rm("requireNamespace", envir = e)
fails(e$pid_alive(NULL, NULL), "Invalid worker ownership")
fails(e$pid_alive(-1L, 1), "Invalid worker ownership")
fails(e$pid_alive(1L, NA_real_), "Invalid worker ownership")
saveRDS(list(pid = 42L), child_owner)
fails(e$active_owner(child_owner), "Invalid worker ownership")
unlink(child_owner)
cat("PASS: controllers/orphan children, blocked launch, resource readiness, stale-source and child guards\n")

# The serialized geographic worker must contain the new constants and guards.
dump_dir <- file.path(root, "dump")
dir.create(dump_dir)
dump_path <- e$dump_engine(dump_dir)
dumped <- new.env(parent = globalenv())
sys.source(dump_path, envir = dumped)
check(identical(dumped$V2_ANALYSES, e$V2_ANALYSES), "Dumped worker lacks v2 analyses")
check(all(dumped$source_status(env_run, TRUE)$table$ready), "Dumped selector broken")
cat("PASS: serialized worker includes v1/v2 checker\n")

native <- requireNamespace("ps", quietly = TRUE)
if ("--require-native" %in% commandArgs(trailingOnly = TRUE))
  check(native, "Native ps checks required")
if (native) {
  handle <- ps::ps_handle()
  pid <- ps::ps_pid(handle)
  created <- ps::ps_create_time(handle)
  check(e$pid_alive(pid, created), "Current process not recognized")
  check(!e$pid_alive(pid, as.numeric(created) + 1), "Reused PID identity accepted")
  check(!e$pid_alive(.Machine$integer.max, created), "Missing PID accepted")
  # Put this existing process in an isolated-owner marker; never start a worker.
  for (path in owner_paths) if (file.exists(path)) unlink(path)
  saveRDS(list(pid = pid, create_time = created), child_owner)
  invisible(capture.output(pf <- pheno_geo_mundlak$preflight(root, env_run)))
  check(child_owner %in% pf$active_jobs && !pf$ready, "Native orphan owner not blocking")
  fails(pheno_geo_mundlak$run(root, env_run), "Another known")
  cat("PASS: native process identity, missing PID, isolated-owner public gate\n")
} else cat("SKIP: native ps unavailable; rerun with --require-native in the configured R library\n")
unlink(root, recursive = TRUE)
cat("ALL GEOGRAPHIC CHECKER TESTS PASSED:", checks, "assertions\n")
