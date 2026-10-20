# Metadata and tiny in-memory contract tests. No native likelihood, fit or subprocess.
source("mundlak_phylogeny_sources.R")
source("geography_mundlak_by_voltinism_qc_final.R")
source("11c_environment_mundlak_phylogeny_15km.R")
source("10g_mundlak_demography_phylogeny_15km.R")
source("run_mundlak_phylogeny_weekend.R")
h <- pheno_mundlak_phylo_sources
he <- environment(h$guard)
g <- environment(pheno_geo_mundlak$run)
ee <- environment(pheno_env_mundlak_phylogeny$run)
de <- environment(pheno_demography_mundlak_phylogeny$run)
qe <- environment(pheno_mundlak_phylogeny_weekend$start)
nchecks <- 0L
check <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
  nchecks <<- nchecks + 1L
}
fails <- function(expr, pattern) {
  msg <- tryCatch({ force(expr); NA_character_ }, error = function(e) conditionMessage(e))
  check(!is.na(msg) && grepl(pattern, msg), paste("Expected", pattern, "got", msg))
}
root <- tempfile("phylo-contract-")
dir.create(root)
root <- normalizePath(root, winslash = "/", mustWork = TRUE)
scripts <- c("mundlak_phylogeny_sources.R", "geography_mundlak_by_voltinism_qc_final.R",
  "11c_environment_mundlak_phylogeny_15km.R", "10g_mundlak_demography_phylogeny_15km.R",
  "run_mundlak_phylogeny_weekend.R", "pheno_phylogeny_qc_15km_v2.R",
  "pheno_plasticity_phylogeny_pilot_15km.R")
check(all(file.copy(scripts, root)), "Cannot copy fixture scripts")
env_run <- file.path(root, "output/phenology_plasticity/environment_mundlak_by_voltinism_qc_final_15km/run_test")
dir.create(env_run, recursive = TRUE)
signature <- strrep("a", 64)
saveRDS(list(version = g$SOURCE_VERSION, signature = signature), file.path(env_run, "run_config.rds"))
write_checks <- function(path, gradient = .0001) write.csv(data.frame(
  numerical_checks_ok = gradient < .001, convergence_code = 0L,
  pdHess = TRUE, max_abs_gradient = gradient), path, row.names = FALSE)
refinement <- function(path, a, status = "FIT_OK_REFINED") {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  writeLines("DO NOT READ FIT", file.path(path, "fit_refined_final.rds"))
  saveRDS(list(summary = data.frame(analysis = a, refinement_status = status)),
    file.path(path, "refinement_completed.rds"))
  write_checks(file.path(path, "fit_checks_refined.csv"))
}
for (a in h$analyses) {
  d <- file.path(env_run, a)
  dir.create(d)
  writeLines("DO NOT READ INPUT", file.path(d, "prepared_input.rds"))
  writeLines("DO NOT READ FIT", file.path(d, "fit_mundlak.rds"))
  saveRDS(list(status = "FIT_OK", signature = signature), file.path(d, "completed.rds"))
  write_checks(file.path(d, "fit_checks.csv"))
}
audit <- ee$env_phylo_audit_sources(root)
check(nrow(audit) == 4L && all(audit$accepted), "Four original sources not recognized")
fails(ee$env_phylo_spec("onset"), "Select exactly")
for (a in h$analyses) {
  s <- ee$env_phylo_spec(a)
  onset <- startsWith(a, "onset_")
  check(s$window == if (onset) 60L else 90L, "Wrong thermal window")
  check(length(s$control) == if (onset) 0L else 1L, "Wrong onset control")
  check(s$orientation == if (a == "offset_multivoltine") 1 else -1, "Wrong biological sign")
  d <- as.data.frame(matrix(sin(seq_len(4000) * .73), 200, 20))
  names(d)[seq_along(s$moderators)] <- s$moderators
  d[[s$anomaly]] <- cos(seq_len(200) * .39)
  d$ONSET_mean_z <- sin(seq_len(200) * .23)
  d[[s$response]] <- seq_len(200)
  f <- as.formula(paste(s$response, "~", s$anomaly, "*(",
    paste(s$moderators, collapse = "+"), ")",
    if (!onset) "+ ONSET_mean_z" else ""))
  check(ncol(model.matrix(f, d)) == if (onset) 18L else 19L, "Wrong fixed design size")
}
v1 <- file.path(env_run, "postfit_refinement_queue_v1")
v2 <- file.path(env_run, "postfit_refinement_queue_v2_isolated")
refinement(file.path(v1, "offset_multivoltine"), "offset_multivoltine")
refinement(file.path(v1, "onset_multivoltine"), "onset_multivoltine")
dir.create(v2)
audit <- ee$env_phylo_audit_sources(root)
check(sum(audit$accepted) == 1L && all(audit$refinement_version[audit$analysis != "offset_multivoltine"] == "v2"),
  "Pending v2 fell back to accepted originals/v1")
for (a in g$V2_ANALYSES) refinement(file.path(v2, a), a)
saveRDS(list(source_run = env_run, old_refinement = v1), file.path(v2, "refinement_config.rds"))
audit <- ee$env_phylo_audit_sources(root)
check(all(audit$accepted), "Mixed v1/v2 not accepted")
for (a in h$analyses) {
  r <- ee$env_phylo_resolve_source(root, a, env_run)
  check(identical(r$prepared_file, file.path(env_run, a, "prepared_input.rds")), "Refinement changed original input")
  check(grepl("fit_refined_final.rds$", r$source_checkpoint), "Original fit selected over refinement")
  p <- list(version = g$SOURCE_VERSION, analysis = a, spec = ee$env_phylo_spec(a),
    data = data.frame(voltinism_lookup = rep(ee$env_phylo_spec(a)$voltinism, 3)),
    classification_audit = data.frame(group = ee$env_phylo_spec(a)$voltinism, retained_rows = 3L))
  saved <- list(analysis = a, source_signature = signature, source_run = env_run,
    refinement_status = "FIT_OK_REFINED")
  cfg <- list(analysis = a)
  ee$env_phylo_validate_identity(cfg, r, readRDS(r$run_config_file),
    readRDS(r$original_completed_file), p, saved)
  bad <- saved; bad$source_signature <- "wrong"
  fails(ee$env_phylo_validate_identity(cfg, r, readRDS(r$run_config_file),
    readRDS(r$original_completed_file), p, bad), "provenance")
  if (startsWith(a, "onset_")) {
    p$data$voltinism_lookup[1] <- "wrong"
    fails(ee$env_phylo_validate_identity(cfg, r, readRDS(r$run_config_file),
      readRDS(r$original_completed_file), p, saved), "voltinism group")
  }
}
bad_v2 <- file.path(v2, "onset_multivoltine")
for (status in c("NUMERICAL_REVIEW_REQUIRED_AFTER_MAX_PASSES", "REFINING")) {
  refinement(bad_v2, "onset_multivoltine", status)
  check(!h$environment_row(root, env_run, "onset_multivoltine")$accepted, "Bad v2 accepted")
}
refinement(bad_v2, "onset_multivoltine")
write_checks(file.path(bad_v2, "fit_checks_refined.csv"), .001)
check(!h$environment_row(root, env_run, "onset_multivoltine")$accepted, "Threshold relaxed")
refinement(bad_v2, "onset_multivoltine")
saveRDS(list(status = "OK"), file.path(bad_v2, "isolated_worker_result.rds"))
unlink(file.path(bad_v2, "refinement_completed.rds"))
check(!h$environment_row(root, env_run, "onset_multivoltine")$accepted, "Worker status substituted for acceptance")
refinement(bad_v2, "onset_multivoltine")
cat("PASS: four formulas/specifications, v1/v2 precedence, refined wrapper and group identities\n")

# No fallback to old demographics when updated fits do not exist.
audit <- de$demo_phylo_audit_sources(root)
check(nrow(audit) == 2L && !any(audit$accepted) &&
  all(grepl("post_environment_mundlak_15km", audit$source_directory)), "Wrong demographic default")
demo_run <- file.path(root, "output/population_trends/post_environment_mundlak_15km/run_demo")
dir.create(demo_run, recursive = TRUE)
demo_signature <- strrep("b", 64)
src <- g$source_status(env_run, TRUE)
demo_config <- list(version = "post_environment_mundlak_1.1", mode = "fit", signature = demo_signature,
  environment_run = env_run, environment = src$table, manifest = g$file_manifest(c(src$files, src$table$fit)))
saveRDS(demo_config, file.path(demo_run, "run_config.rds"))
for (a in h$analyses) {
  d <- file.path(demo_run, "extraction", a)
  dir.create(d, recursive = TRUE)
  saveRDS(list(signature = demo_signature, analysis = a, n_populations = 2L),
    file.path(d, "completed.rds"))
}
for (group in c("univoltine", "multivoltine")) {
  d <- file.path(demo_run, paste0(group, "__mundlak"))
  dir.create(d)
  writeLines("DO NOT READ FIT", file.path(d, "fit.rds"))
  writeLines("DO NOT READ INPUT", file.path(demo_run, paste0("prepared_mundlak_", group, ".rds")))
  saveRDS(list(signature = demo_signature, group = group, nobs = 6L, npop = 2L),
    file.path(demo_run, paste0("prepared_", group, "_completed.rds")))
  write.csv(data.frame(group = group, numerical_checks_ok = TRUE, convergence_code = 0L,
    pdHess = TRUE, gradient_covariance_norm = .01, nobs = 6L), file.path(d, "fit_checks.csv"), row.names = FALSE)
}
result <- list(signature = demo_signature, status = "FITTED_RESIDUAL_REVIEW_PENDING",
  groups = data.frame(group = c("multivoltine", "univoltine"), status = "FIT_OK"))
saveRDS(result, file.path(demo_run, "result_prepare.rds"))
check(!any(de$demo_phylo_audit_sources(root)$accepted), "Preparation-only run accepted")
saveRDS(result, file.path(demo_run, "result_fit.rds"))
audit <- de$demo_phylo_audit_sources(root)
check(all(audit$accepted) && all(audit$expected_nobs == 6L), "Updated sample rejected")
bad <- demo_config; bad$environment$fit[1] <- "old-source"
saveRDS(bad, file.path(demo_run, "run_config.rds"))
check(!any(de$demo_phylo_audit_sources(root)$accepted), "Stale plasticity provenance accepted")
saveRDS(demo_config, file.path(demo_run, "run_config.rds"))
bad <- result; bad$groups$status[1] <- "NUMERICAL_REVIEW_REQUIRED"
saveRDS(bad, file.path(demo_run, "result_fit.rds"))
check(!any(de$demo_phylo_audit_sources(root)$accepted), "Unaccepted demographics passed")
saveRDS(result, file.path(demo_run, "result_fit.rds"))
marker <- file.path(demo_run, "extraction", h$analyses[1], "completed.rds")
good <- readRDS(marker)
saveRDS(list(signature = "wrong", analysis = h$analyses[1]), marker)
check(!any(de$demo_phylo_audit_sources(root)$accepted), "Mixed extraction signatures accepted")
saveRDS(good, marker)
cat("PASS: updated demographic discovery, completed fit requirement, sample and plasticity provenance\n")

# Mock process queries only; audit source validation itself remains real.
real_geo <- he$geo
real_pid <- g$pid_alive
g$pid_alive <- function(pid, create_time) isTRUE(pid == 42L)
g$requireNamespace <- function(...) TRUE
he$geo <- function(root) g
ee$env_phylo_sources <- function(root) h
de$demo_phylo_sources <- function(root) h
qe$sources_api <- function(root) h
real_load <- qe$load_api
qe$load_api <- function(root, family) {
  api <- new.env(parent = globalenv())
  if (family == "environment") api$audit_environment_mundlak_phylogeny_sources <-
    function(project_root, analysis, source_run) ee$env_phylo_audit_sources(project_root, analysis, source_run)
  else api$audit_demography_mundlak_phylogeny_sources <-
    function(project_root, group, source_run) de$demo_phylo_audit_sources(project_root, group, source_run)
  api
}
audit <- qe$audit_impl(root, env_run, demo_run)
check(nrow(audit) == 6L && all(audit$ready), "Six accepted sources did not make queue ready")
check(length(attr(audit, "source_snapshots")) == 6L, "Source snapshots missing")
check(sum(audit$family == "environment") == 4L && sum(audit$family == "demography") == 2L,
  "Incorrect queue composition")
other_env <- file.path(root, "other_environment_run")
dir.create(other_env)
saveRDS(list(version = g$SOURCE_VERSION, signature = signature), file.path(other_env, "run_config.rds"))
for (a in h$analyses) {
  d <- file.path(other_env, a)
  dir.create(d)
  file.copy(file.path(env_run, a, c("fit_mundlak.rds", "prepared_input.rds", "completed.rds", "fit_checks.csv")), d)
}
mixed <- qe$audit_impl(root, other_env, demo_run)
check(!any(mixed$ready) && all(grepl("provenance differ", mixed$reason)), "Mixed environmental runs accepted")
legacy <- file.path(root, "legacy_demo")
dir.create(legacy)
saveRDS(list(version = de$demo_phylo_version(), run_signature = strrep("c", 64)),
  file.path(legacy, "run_config.rds"))
for (group in c("univoltine", "multivoltine")) {
  d <- file.path(legacy, paste0(group, "__mundlak"))
  dir.create(d)
  writeLines("DO NOT READ FIT", file.path(d, "fit.rds"))
  writeLines("DO NOT READ INPUT", file.path(legacy, paste0("prepared_mundlak_", group, ".rds")))
  write.csv(data.frame(group = group, numerical_checks_ok = TRUE), file.path(d, "fit_checks.csv"), row.names = FALSE)
}
check(all(de$demo_phylo_audit_sources(root, source_run = legacy)$accepted), "Explicit legacy audit broken")
mixed <- qe$audit_impl(root, env_run, legacy)
check(!any(mixed$ready[mixed$family == "demography"]) &&
  all(grepl("updated post-environment", mixed$reason[mixed$family == "demography"])),
  "Six-model queue accepted outdated demographic predictors")
child_owner <- file.path(v2, "offset_univoltine", "isolated_worker_owner.rds")
saveRDS(list(pid = 42L, create_time = 1), child_owner)
saveRDS(list(pid = 43L, create_time = 1), file.path(v2, "worker_owner.rds"))
audit <- qe$audit_impl(root, env_run, demo_run)
check(!any(audit$ready) && all(grepl("worker is active", audit$reason)), "Orphan child did not block queue")
before <- sort(list.files(root, recursive = TRUE, all.files = TRUE))
invisible(capture.output(pf <- start_mundlak_phylogeny_weekend(root, environment_source_run = env_run,
  demography_source_run = demo_run)))
check(nrow(pf) == 6L && !any(pf$ready), "Default start launched instead of auditing")
fails(suppressMessages(invisible(capture.output(start_mundlak_phylogeny_weekend(root, check_only = FALSE,
  environment_source_run = env_run, demography_source_run = demo_run)))), "All six")
check(identical(before, sort(list.files(root, recursive = TRUE, all.files = TRUE))), "Blocked queue wrote files")
fails(de$run(root, source_run = demo_run, check_only = TRUE), "worker is active")
cfg <- list(project_root = root, source_directory = file.path(env_run, "onset_multivoltine"),
  analysis = "onset_multivoltine")
fails(ee$env_phylo_prepare(cfg, NULL), "worker is active")
cfg$source_directory <- file.path(demo_run, "univoltine__mundlak"); cfg$group <- "univoltine"
fails(de$demo_phylo_prepare(cfg, NULL), "worker is active")
unlink(child_owner)
snapshot <- h$snapshot(root, "environment", file.path(env_run, "onset_multivoltine"), "onset_multivoltine")
writeLines("CHANGED FIT PLACEHOLDER", file.path(bad_v2, "fit_refined_final.rds"))
fails(h$guard(root, "environment", file.path(env_run, "onset_multivoltine"), "onset_multivoltine", snapshot),
  "changed after audit")
check(!any(de$demo_phylo_audit_sources(root)$accepted), "Changed environmental artifact did not invalidate demographics")
writeLines("DO NOT READ FIT", file.path(bad_v2, "fit_refined_final.rds"))
demo_config$manifest <- g$file_manifest(c(src$files, src$table$fit))
saveRDS(demo_config, file.path(demo_run, "run_config.rds"))
g$pid_alive <- real_pid
rm("requireNamespace", envir = g)
he$geo <- real_geo
qe$load_api <- real_load
cat("PASS: six-row queue, idle/source guards, orphan children, no writes on blocked launch\n")

native <- requireNamespace("ps", quietly = TRUE) && requireNamespace("digest", quietly = TRUE)
if ("--require-native" %in% commandArgs(trailingOnly = TRUE)) check(native, "Native ps/digest required")
if (native) {
  handle <- ps::ps_handle()
  saveRDS(list(pid = ps::ps_pid(handle), create_time = ps::ps_create_time(handle)), child_owner)
  audit <- qe$audit_impl(root, env_run, demo_run)
  check(!any(audit$ready) && all(grepl("worker is active", audit$reason)), "Native orphan process missed")
  fails(ee$run(root, source_run = env_run, check_only = TRUE), "worker is active")
  unlink(child_owner)
  # Tiny post-environment input identity, with a deliberately changed sample.
  group <- "univoltine"
  cfg <- list(project_root = root, group = group, expected_samples = NULL)
  resolved <- de$demo_phylo_resolve_source(root, group, demo_run)
  p <- list(version = demo_config$version, group = group, predictors_updated = TRUE,
    source_bridge_signature = "old", source_data_signature = "old",
    frame = data.frame(pop_id = rep(c("a", "b"), each = 3)), data = list(y = 1:6))
  p$signature <- digest::digest(list(demo_signature, group, p$frame, p$data), algo = "sha256")
  fit <- p[c("version", "group", "signature", "source_bridge_signature", "source_data_signature")]
  identity <- de$demo_phylo_validate_identity(cfg, resolved, demo_config, p, fit)
  check(identity$post && identity$sample[["nobs"]] == 6L, "Updated input signature/sample not accepted")
  bad <- p; bad$data$y[1] <- 10
  fails(de$demo_phylo_validate_identity(cfg, resolved, demo_config, bad, fit), "saved signature")
  bad <- p; bad$predictors_updated <- FALSE
  fails(de$demo_phylo_validate_identity(cfg, resolved, demo_config, bad, fit), "Updated demographic")
  wrong <- resolved; wrong$expected_nobs <- 7L
  fails(de$demo_phylo_validate_identity(cfg, wrong, demo_config, p, fit), "sample differs")
  # A source change between audit and controller launch must fail before the runner.
  dir <- file.path(root, "controller_test")
  dir.create(dir)
  lock <- file.path(dir, "controller.lock")
  dir.create(lock)
  token <- "synthetic-controller"
  saveRDS(list(token = token), file.path(lock, "owner.rds"))
  handoff <- file.path(dir, "handoff.rds")
  saveRDS(TRUE, handoff)
  task <- qe$tasks()[1, , drop = FALSE]
  task$source_directory <- file.path(env_run, task$analysis)
  snapshot <- h$snapshot(root, "environment", task$source_directory, task$analysis)
  snapshot$bytes[1] <- snapshot$bytes[1] + 1
  controller_cfg <- list(project_root = root, lock = lock, token = token, handoff_file = handoff,
    tasks = task, source_snapshots = setNames(list(snapshot), task$task),
    attempt_directory = dir, output_directory = dir, progress_file = file.path(dir, "progress.rds"),
    stop_on_error = TRUE)
  real_load <- qe$load_api
  qe$load_api <- function(...) stop("RUNNER MUST NOT BE REACHED")
  result <- qe$controller(controller_cfg)
  qe$load_api <- real_load
  check(!result$all_fit_ok && result$results$status == "WORKER_FAILED" &&
    grepl("changed after audit", result$results$reason), "Controller did not reject changed source before launch")
  # Build the R engine API only: do not call start(), native preflight or fit.
  for (a in h$analyses) {
    engine <- ee$make_engine(root, a, .05, .001)$environment
    check(engine$anomaly == ee$env_phylo_spec(a)$anomaly, "Wrong adapter engine anomaly")
    check(is.function(engine$env_phylo_validate_identity) && is.function(engine$env_phylo_sources),
      "Engine omitted new helpers")
  }
  cat("PASS: native process identity, post-environment content signatures, four engine adapters (no fits)\n")
} else cat("SKIP: native ps/digest checks; rerun with --require-native\n")
unlink(root, recursive = TRUE)
cat("ALL PHYLOGENETIC SOURCE TESTS PASSED:", nchecks, "assertions\n")
