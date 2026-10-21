# Contract tests: synthetic metadata/objective callbacks only.
# No real fits, native likelihoods, compilation, optimizers or child processes.
source("env_mundlak_rescue_offset_multivoltine.R")
api <- pheno_env_offset_rescue
e <- environment(api$run)
n <- 0L
check <- function(ok, msg) {
  if (!isTRUE(ok)) stop(msg, call. = FALSE)
  n <<- n + 1L
}
fails <- function(expr, pattern) {
  message <- tryCatch({ force(expr); NA_character_ }, error = function(err) conditionMessage(err))
  check(!is.na(message) && grepl(pattern, message), paste("Expected", pattern, "got", message))
}
ctl <- e$settings()
check(ctl$max_passes == 2 && ctl$stagnation_patience == 1, "Unbounded/late stagnation defaults")
check(ctl$rel.tol == 1e-12 && ctl$x.tol == 1e-10, "Explicit tolerance defaults lost")
fails(e$settings(list(raw_gradient_threshold = .01)), "Unknown")
fails(e$settings(list(max_passes = 1.5)), "Integer")
fails(e$settings(list(maximum_total_seconds = Inf)), "Invalid")
fails(e$settings(list(repeat_evaluations = 1)), "Invalid")
fails(e$settings(list(gradient_progress_rel = 1)), "Invalid")

fake_objective <- function() {
  state <- new.env(parent = emptyenv())
  state$last.par <- c(b_j = 1, b_j = 1)
  state$random <- integer()
  state$parList <- function(par) list(b_j = par)
  list(env = state, parameter_labels = c("intercept", "slope"),
    fn = function(theta) { state$last.par <- theta; sum(theta^2) },
    gr = function(theta) { state$last.par <- theta; 2 * theta })
}
obj <- fake_objective()
theta <- c(b_j = -1, b_j = .2)
bounds <- list(lower = rep(-Inf, 2), upper = rep(Inf, 2))
d <- e$repeat_diagnostics(obj, theta, bounds$lower, bounds$upper, ctl)
check(d$stable && d$max_gradient == 2, "Absolute negative gradient not checked")
check(identical(d$parameters$id, c("b_j", "b_j.1")), "Duplicated parameter IDs")
check(identical(d$parameters$term, c("intercept", "slope")), "Coefficient labels missing")
check(nrow(d$evaluations) == 3 && all(d$parameters$gradient_spread == 0), "Repeatability evidence missing")
fails(e$repeat_diagnostics(obj, theta, c(0, 0), bounds$upper, ctl), "Invalid outer")
bad <- obj
bad$gr <- function(theta) c(NA_real_, 0)
fails(e$repeat_diagnostics(bad, theta, bounds$lower, bounds$upper, ctl), "Non-finite")
index <- 0L
bad$gr <- function(theta) { index <<- index + 1L; c(.01 * index, .001) }
unstable <- e$repeat_diagnostics(bad, theta, bounds$lower, bounds$upper, ctl)
check(!unstable$stable, "Unstable gradient accepted")
bad <- obj
index <- 0L
bad$fn <- function(theta) { index <<- index + 1L; sum(theta^2) + index * .001 }
check(!e$repeat_diagnostics(bad, theta, bounds$lower, bounds$upper, ctl)$stable,
  "Unstable objective accepted")

point <- function(f = 681680.578630163, g = .00299802165837006, stable = TRUE) {
  list(objective = f, max_gradient = g, stable = stable, objective_noise = 0)
}
before <- point()
after <- point()
decision <- e$decision(before, after, 0L, ctl)
check(decision$stagnated && decision$reason == "STAGNATION" && !decision$qualifies,
  "Observed v1 plateau would run again")
check(e$decision(point(g = .024939), after, 0L, ctl)$reason == "PROGRESS",
  "Real gradient progress ignored when objective barely moves")
check(e$decision(before, point(f = before$objective + .01), 0L, ctl)$reason == "OBJECTIVE_WORSE",
  "Worse objective accepted")
check(!e$decision(before, point(g = .01), 0L, ctl)$retain, "Worse gradient without progress retained")
check(e$decision(before, point(g = .0009), 0L, ctl)$qualifies, "Target not detected")
check(!e$decision(before, point(g = .001), 0L, ctl)$qualifies, "Strict gradient threshold relaxed")
check(!e$decision(before, point(g = .0001), 1L, ctl)$qualifies, "Nonzero convergence accepted")
check(!e$decision(before, point(g = .0001, stable = FALSE), 0L, ctl)$retain,
  "Unstable candidate retained")

root <- tempfile("offset-rescue-contract-")
dir.create(root)
root <- normalizePath(root, winslash = "/")
run <- file.path(root, "source")
original <- file.path(run, "offset_multivoltine")
refined <- file.path(run, "postfit_refinement_queue_v1", "offset_multivoltine")
dir.create(original, recursive = TRUE)
dir.create(refined, recursive = TRUE)
signature <- strrep("a", 64)
saveRDS(list(version = e$SOURCE_VERSION, signature = signature), file.path(run, "run_config.rds"))
saveRDS(list(signature = signature, status = "NUMERICAL_REVIEW_REQUIRED"),
  file.path(original, "completed.rds"))
saveRDS(list(source_status = "NUMERICAL_REVIEW_REQUIRED", summary = data.frame(
  analysis = "offset_multivoltine", refinement_status = "NUMERICAL_REVIEW_REQUIRED_AFTER_MAX_PASSES")),
  file.path(refined, "refinement_completed.rds"))
writeLines("DO NOT LOAD THIS FIT", file.path(refined, "fit_refined_final.rds"))
writeLines("DO NOT LOAD THESE DATA", file.path(original, "prepared_input.rds"))
saveRDS(list(b_j = 1), file.path(refined, "parameters_refined.rds"))
write.csv(data.frame(term = "intercept", estimate = 1),
  file.path(original, "fixed_effects.csv"), row.names = FALSE)
write.csv(data.frame(objective = 681680.578630163, max_abs_gradient = .002998,
  nobs = 10, convergence_code = 0L, pdHess = TRUE),
  file.path(refined, "fit_checks_refined.csv"), row.names = FALSE)
real_geo <- e$geo
real_resources <- e$resources
active <- character()
e$geo <- function(root) list(
  refinement_roots = function(run) list(v1 = file.path(run, "postfit_refinement_queue_v1"),
    config_files = character()),
  active_jobs = function(root, environment_run = NULL) active)
e$resources <- function(root) list(packages = c(sdmTMB = "1.1.0"),
  available_ram_GiB = 100, available_disk_GiB = 100)
audit <- api$audit(root, run, quiet = TRUE)
check(audit$ready && identical(audit$analysis, "offset_multivoltine"), "Fixture audit failed")
check(!dir.exists(file.path(root, "output")), "Audit wrote workflow files")
loaded <- FALSE
original_loader <- e$load_source
e$load_source <- function(...) { loaded <<- TRUE; stop("MODEL MUST NOT LOAD") }
check(api$run(root, run)$ready && !loaded, "Default run loaded model")
active <- "orphan/isolated_worker_owner.rds"
fails(api$run(root, run, mode = "diagnose"), "blocked before loading")
fails(api$run(root, run, mode = "rescue"), "blocked before loading")
check(!loaded && !dir.exists(file.path(root, "output")), "Blocked mode caused writes/load")
fails(e$guard(audit), "Active")
active <- character()
e$resources <- function(root) list(available_ram_GiB = 2, available_disk_GiB = 100)
check(!api$audit(root, run, quiet = TRUE)$ready, "RAM guard failed")
e$resources <- function(root) list(available_ram_GiB = 100, available_disk_GiB = 2)
check(!api$audit(root, run, quiet = TRUE)$ready, "Disk guard failed")
e$resources <- function(root) stop("Cannot inspect resources")
check(!api$audit(root, run, quiet = TRUE)$ready, "Resource inspection failure allowed")
e$resources <- function(root) list(available_ram_GiB = 100, available_disk_GiB = 100)
e$guard(audit)
check(TRUE, "Unchanged manifest")
metadata <- file.path(original, "fixed_effects.csv")
old_time <- file.info(metadata)$mtime
writeLines(sub("intercept", "Intercept", readLines(metadata)), metadata)
Sys.setFileTime(metadata, old_time)
fails(e$guard(audit), "Source changed")
write.csv(data.frame(term = "intercept", estimate = 1), metadata, row.names = FALSE)

out <- file.path(root, "synthetic")
dir.create(out)
fake_a <- list(signature = signature, manifest = data.frame(path = "synthetic"))
obj <- fake_objective()
start <- e$repeat_diagnostics(obj, c(b_j = 1, b_j = 1), bounds$lower, bounds$upper, ctl)
calls <- 0L
stalled_optimizer <- function(obj, current, bounds, a, out, pass, ctl, tick) {
  calls <<- calls + 1L
  list(optimizer = list(par = current$theta, convergence = 0L,
    objective = current$objective, message = "relative convergence"), seconds = 0)
}
loop <- e$rescue_loop(obj, start, bounds, fake_a, out, ctl, function() NULL, stalled_optimizer)
check(loop$status == "STAGNATION" && calls == 1L, "Repeated optimizer on plateau")
check(file.exists(file.path(out, "retained_checkpoint.rds")), "Small checkpoint missing")
check(!file.exists(file.path(out, "fit_rescued_final.rds")), "Plateau exported accepted fit")
cp <- readRDS(file.path(out, "retained_checkpoint.rds"))
check(!cp$numerically_accepted && identical(cp$theta, start$theta), "Checkpoint mislabeled")
progress_optimizer <- function(obj, current, bounds, a, out, pass, ctl, tick) {
  calls <<- calls + 1L
  x <- current$theta * if (pass == 1) .5 else 0
  list(optimizer = list(par = x, convergence = 0L, objective = sum(x^2), message = "mock"),
    seconds = 0)
}
calls <- 0L
loop <- e$rescue_loop(obj, start, bounds, fake_a, out, ctl, function() NULL, progress_optimizer)
check(loop$status == "GRADIENT_TARGET_REACHED" && calls == 2L, "Progress/target transition failed")
ctl_one <- e$settings(list(max_passes = 1))
loop_limited <- e$rescue_loop(obj, start, bounds, fake_a, out, ctl_one, function() NULL, progress_optimizer)
check(loop_limited$status == "MAX_PASSES_REACHED", "Pass limit not enforced")
fails(e$checkpoint(obj, c(b_j = 99, b_j = 99), fake_a, "test"), "current state differs")

fake_minimize <- function(start, objective, gradient, lower, upper, control) {
  check(control$rel.tol == 1e-12 && control$iter.max == 100, "Optimizer ignored controls")
  value <- objective(start)
  g <- gradient(start)
  warning("synthetic warning")
  list(par = start, objective = value, convergence = 0L, iterations = 1L,
    evaluations = c(functions = 1L, gradient = 1L), message = "mock")
}
record <- e$optimize_once(obj, start, bounds, fake_a, out, 9L, ctl,
  function() NULL, minimize = fake_minimize)
check(record$callback_counts[["fn"]] == 1L && record$callback_counts[["gr"]] == 1L,
  "Per-pass counts missing")
check(identical(record$warnings, "synthetic warning"), "Warnings not preserved")
check(file.exists(file.path(out, "evaluation_checkpoint.rds")), "Evaluation checkpoint missing")
writeLines("request", file.path(out, "STOP_REQUESTED"))
fails(e$optimize_once(obj, start, bounds, fake_a, out, 10L, ctl,
  function() e$budget(Inf, out), minimize = fake_minimize), "STOP_REQUESTED")
unlink(file.path(out, "STOP_REQUESTED"))
fails(e$budget(-Inf, out), "TIME_BUDGET_REACHED")
old_elapsed <- e$elapsed
fake_clock <- 0
e$elapsed <- function() { fake_clock <<- fake_clock + 1; fake_clock }
fails(e$optimize_once(obj, start, bounds, fake_a, out, 11L,
  e$settings(list(maximum_pass_seconds = .1)), function() NULL,
  minimize = fake_minimize), "PASS_TIME_BUDGET_REACHED")
e$elapsed <- old_elapsed

final <- e$repeat_diagnostics(obj, c(b_j = 0, b_j = 0), bounds$lower, bounds$upper, ctl)
opt <- list(convergence = 0L)
sd <- list(par.fixed = final$theta, gradient.fixed = final$gradient,
  pdHess = TRUE, cov.fixed = diag(2))
check(e$final_gate(final, opt, sd, 2, ctl)$accepted, "Valid final gate failed")
bad_sd <- sd; bad_sd$par.fixed[[1]] <- .01
check(!e$final_gate(final, opt, bad_sd, 2, ctl)$accepted, "Wrong covariance point accepted")
bad_sd <- sd; bad_sd$pdHess <- FALSE
check(!e$final_gate(final, opt, bad_sd, 2, ctl)$accepted, "Bad Hessian accepted")
bad_sd <- sd; bad_sd$cov.fixed[1, 1] <- -1
check(!e$final_gate(final, opt, bad_sd, 2, ctl)$accepted, "Indefinite covariance accepted")
bad_sd <- sd; bad_sd$gradient.fixed[[1]] <- -.003
check(!e$final_gate(final, opt, bad_sd, 2, ctl)$accepted, "Negative excessive sdreport gradient accepted")
check(!e$final_gate(final, list(convergence = 1L), sd, 2, ctl)$accepted, "Failed convergence accepted")
check(!e$final_gate(final, opt, sd, -2, ctl)$accepted, "Worsened final objective accepted")
check(!e$final_gate(utils::modifyList(final, list(stable = FALSE)), opt, sd, 2, ctl)$accepted,
  "Unstable final point accepted")

if ("--require-native" %in% commandArgs(trailingOnly = TRUE)) {
  check(requireNamespace("ps", quietly = TRUE), "Native process inspection required")
  check(file.copy("geography_mundlak_by_voltinism_qc_final.R", root), "Cannot copy worker checker")
  e$geo <- real_geo
  child <- file.path(run, "postfit_refinement_queue_v2_isolated", "offset_univoltine")
  dir.create(child, recursive = TRUE)
  owner <- file.path(child, "isolated_worker_owner.rds")
  saveRDS(list(pid = Sys.getpid(), create_time = ps::ps_create_time(ps::ps_handle())), owner)
  native <- api$audit(root, run, quiet = TRUE)
  check(!native$ready && any(grepl("isolated_worker_owner", native$active_workers)),
    "Orphan isolated v2 worker not detected")
  fails(api$run(root, run, mode = "rescue"), "blocked before loading")
  check(!loaded && !dir.exists(file.path(root, "output")), "Native blocked run loaded/wrote")
  e$resources <- real_resources
  res <- e$resources(root)
  check(res$packages[["sdmTMB"]] == "1.1.0", "Wrong installed sdmTMB")
}
e$geo <- real_geo
e$resources <- real_resources
e$load_source <- original_loader
cat("PASS:", n, "rescue assertions; no real fits, optimizers, native likelihoods or workers started.\n")
