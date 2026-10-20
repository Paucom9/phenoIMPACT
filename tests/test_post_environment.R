# Run from the repository root: Rscript tests/test_post_environment.R
# Synthetic contracts, not a substitute for replaying the real saved TMB fits.
source("post_environment_mundlak.R")
source("10f_mundlak_plasticity_demography_QC_FINAL_spatial_15km_RECOVERY.R")
e <- environment(pheno_post_env$run)
engine <- environment(pheno_mundlak_qc$run)
check <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)
fails <- function(expr, pattern) {
  msg <- tryCatch({ force(expr); NA_character_ }, error = function(x) conditionMessage(x))
  check(!is.na(msg) && grepl(pattern, msg), paste("Expected failure:", pattern, "got", msg))
}
set.seed(72)

# 1. Analytic contextual derivative, species mapping, and biological signs.
for (a in e$ANALYSES) {
  s <- e$spec(a)
  pop <- data.frame(SPECIES = rep(c("spC", "spA", "spB"), each = 3),
    SITE_ID = rep(c("site2", "site3", "site1"), 3), bms_id = "test")
  for (v in s$moderators) pop[[v]] <- rnorm(nrow(pop))
  d <- pop[rep(seq_len(nrow(pop)), each = 3), ]
  d$YEAR <- rep(c(2001L, 2002L, 2004L), nrow(pop))
  d[[s$anomaly]] <- rnorm(nrow(d)); d[[s$response]] <- 100 + rnorm(nrow(d))
  if (length(s$control)) d[[s$control]] <- rnorm(nrow(d))
  formula <- stats::as.formula(paste(s$response, "~", s$anomaly, "*(",
    paste(s$moderators, collapse = "+"), ")",
    if (length(s$control)) paste("+", s$control) else ""))
  terms <- colnames(model.matrix(formula, d))
  beta <- setNames(seq_along(terms) / 7, terms)
  if (length(s$control)) beta[s$control] <- 1000
  sp <- c("spB", "spC", "spA")
  modes <- c(.25, -.75, 1.25)
  cn <- list(SITE_ID = "(Intercept)", SPECIES = "(Intercept)",
             SPECIES_slope = s$anomaly, site_year_id = "(Intercept)")
  fit <- list(data = d, split_formula = list(list(form_no_bars = formula,
    re_cov_terms = list(cnms = cn, re_b_df = data.frame(start = 0:2, end = 0:2,
      group_indices = 3L, level = sp)))),
    model = list(par = setNames(unname(beta), rep("b_j", length(beta)))),
    tmb_map = list(re_b_pars = factor(1:3)), parlist = list(re_b_pars = modes),
    sd_report = list(par.fixed = setNames(unname(beta), rep("b_j", length(beta))),
      cov.fixed = diag(.02, length(beta)), value = setNames(modes, rep("re_b_pars", 3)),
      sd = rep(.1, 3)))
  co <- e$coefficients(fit, s)
  out <- e$population_slopes(d, co, s)
  expected <- rep(beta[s$anomaly], nrow(pop))
  for (v in s$moderators) expected <- expected + pop[[v]] * beta[paste(s$anomaly, v, sep = ":")]
  expected <- expected + modes[match(pop$SPECIES, sp)]
  check(e$same(out$plasticity_raw_days_per_C, expected), paste(a, "raw slope incorrect"))
  check(e$same(out$plasticity_bio_days_per_C, s$orientation * expected), paste(a, "sign incorrect"))
  check(all(out$n_years == 3), "Population year count wrong")
  changed <- d; changed[[s$moderators[1]]][2] <- changed[[s$moderators[1]]][2] + 1
  fails(e$population_slopes(changed, co, s), "varies across years")
  bad <- fit; bad$sd_report$par.fixed[1] <- bad$sd_report$par.fixed[1] + .1
  fails(e$coefficients(bad, s), "different point estimates")
  bad <- fit; bad$parlist$re_b_pars[1] <- 9
  fails(e$coefficients(bad, s), "saved modes disagree")
}
cat("PASS: four environmental derivatives, signs, additive onset, species maps, stale reports\n")

# 2. Old demographic data, deliberately nonalphabetical native index order.
species <- paste0("sp", c(4, 2, 7, 1, 5, 3, 6))
sites <- data.frame(SITE_ID = paste0("site", c(3, 1, 5, 2, 4)),
                    x_km = runif(5, 100, 300), y_km = runif(5, 50, 200),
                    latitude_deg = runif(5, 35, 60), longitude_deg = runif(5, -5, 20))
pop <- expand.grid(SPECIES = species, SITE_ID = sites$SITE_ID, stringsAsFactors = FALSE)
pop$pop_id <- paste(pop$SPECIES, pop$SITE_ID, sep = "_")
pop$onset_plasticity_z <- rnorm(nrow(pop)); pop$offset_plasticity_bio_z <- rnorm(nrow(pop))
f <- pop[rep(seq_len(nrow(pop)), each = 3), ]
f$YEAR <- rep(c(2001L, 2002L, 2004L), nrow(pop))
f$year_decade <- (f$YEAR - 2000) / 10; f$ABUND_INDEX <- exp(rnorm(nrow(f)))
f$bms_id <- "test"
dec <- e$decompose(f, sites); f <- dec$frame
X <- dec$X
schema <- list(y = f$ABUND_INDEX, X = unname(X), time = f$year_decade,
  species = match(f$SPECIES, species) - 1L, site = match(f$SITE_ID, sites$SITE_ID) - 1L,
  year = match(f$YEAR, c(2001L, 2002L, 2004L)) - 1L,
  previous = rep(-1L, nrow(f)), gap = rep(1L, nrow(f)),
  A = matrix(runif(15), 5, 3), M0 = diag(3), M1 = diag(3), M2 = diag(3), use_spatial = 1L)
old <- list(group = "univoltine", signature = "old-sig", frame = f, X = X,
  data = schema, sites = sites, years = c(2001L, 2002L, 2004L))
beta <- setNames(rnorm(ncol(X), 0, .02), colnames(X))
pars <- list(beta = unname(beta), b_species = matrix(1:14, 7, 2),
  b_site = matrix(1:10, 5, 2), ar_obs = seq_len(nrow(f)) / 100,
  field = matrix(1:9, 3, 3), log_sd = rep(0, 5))
fit <- list(group = "univoltine", signature = "old-sig", beta = beta, parameters = pars,
  checks = data.frame(numerical_checks_ok = TRUE, Gamma_shape = 3, AR1_rho = .4,
                     gradient_covariance_norm = .0001),
  V = diag(.01, length(beta)), eta = as.numeric(X %*% beta), warnings = character())
dimnames(fit$V) <- list(names(beta), names(beta))
new <- pop[c("SPECIES", "SITE_ID", "pop_id")]; new$bms_id <- "test"
new$plasticity_raw_days_per_C <- rnorm(nrow(pop), -3, 1)
new$plasticity_bio_days_per_C <- -new$plasticity_raw_days_per_C
on <- new
off <- new; off$plasticity_raw_days_per_C <- rnorm(nrow(pop), -2, 1)
off$plasticity_bio_days_per_C <- -off$plasticity_raw_days_per_C
# Drop an entire species from onset and one site from offset.
on <- on[on$SPECIES != "sp7", ]
off <- off[off$SITE_ID != "site5", ]
# Extra extracted population participates in z scaling, but is not added to abundance.
extra <- on[1, ]; extra$SITE_ID <- "extra"; extra$pop_id <- "sp4_extra"
extra$plasticity_raw_days_per_C <- -9; extra$plasticity_bio_days_per_C <- 9
on <- rbind(on, extra)
p <- e$rebuild_input(old, fit, on, off, "univoltine")
check(nrow(p$frame) == 6 * 4 * 3, "Population filtering wrong")
check(nrow(p$parameters$b_species) == 6 && nrow(p$parameters$b_site) == 4, "Latent dimensions wrong")
keep <- p$retained_source_rows
check(identical(p$data$y, old$data$y[keep]), "Response changed")
check(identical(p$parameters$ar_obs, pars$ar_obs[keep]), "AR1 states incorrectly reordered")
check(e$same(p$parameters$b_species, pars$b_species[species != "sp7", , drop = FALSE]), "Species modes reordered")
check(e$same(p$parameters$b_site, pars$b_site[sites$SITE_ID != "site5", , drop = FALSE]), "Site modes reordered")
check(e$same(p$data$A, old$data$A[sites$SITE_ID != "site5", , drop = FALSE]), "Spatial projection reordered")
check(identical(p$parameters$field, pars$field), "Annual spatial field changed")
check(all(p$data$previous[seq(1, nrow(p$frame), 3)] == -1L), "AR1 chains cross populations")
check(identical(p$data$gap, rep(c(1L, 1L, 2L), nrow(p$populations))), "Calendar gaps lost")
expected_z <- (on$plasticity_bio_days_per_C - mean(on$plasticity_bio_days_per_C)) /
  sd(on$plasticity_bio_days_per_C)
check(e$same(p$frame$onset_plasticity_z, expected_z[match(e$key(p$frame), e$key(on))]),
      "Scaling incorrectly restricted to abundance overlap")
check(p$extracted_outside_abundance$n[1] == 1, "Outside-abundance population not audited")
check(all(abs(tapply(p$populations$onset_plasticity_within_species_z,
                    p$populations$SPECIES, mean)) < 1e-10), "Wrong species weighting")
check(e$same(p$frame$latitude_within_species_10deg,
  (p$frame$latitude_deg - p$frame$latitude_species_mean) / 10), "Geography decomposition wrong")
check(identical(colnames(p$X), colnames(old$X)), "Demographic formula changed")
bad <- on; bad$bms_id[1] <- "different"
fails(e$rebuild_input(old, fit, bad, off, "univoltine"), "BMS identity")
bad <- fit; bad$signature <- "other"
fails(e$rebuild_input(old, bad, on, off, "univoltine"), "identity mismatch")
bad <- old; bad$data$unknown_row_array <- rep(0, nrow(f))
fails(e$rebuild_input(bad, fit, on, off, "univoltine"), "Unknown demographic data schema")
cat("PASS: population join, pre-join scaling, within/between, random/AR1 remapping, year gaps\n")

# 3. Lightweight environmental gates: never fall back from a bad refined fit.
tmp <- tempfile(); dir.create(tmp)
for (a in e$ANALYSES) {
  d <- file.path(tmp, a); dir.create(d)
  saveRDS(list(status = "FIT_OK", signature = "env-sig"), file.path(d, "completed.rds"))
  file.create(file.path(d, "fit_mundlak.rds"), file.path(d, "prepared_input.rds"))
  write.csv(data.frame(numerical_checks_ok = TRUE, convergence_code = 0,
    pdHess = TRUE, max_abs_gradient = .0001), file.path(d, "fit_checks.csv"), row.names = FALSE)
}
check(all(e$select_environment(tmp)$ready), "Original FIT_OK sources not selected")
ref <- file.path(tmp, "postfit_refinement_queue_v1", "offset_multivoltine")
dir.create(ref, recursive = TRUE)
file.create(file.path(ref, "fit_refined_final.rds"))
saveRDS(list(summary = data.frame(analysis = "offset_multivoltine",
  refinement_status = "NUMERICAL_REVIEW_REQUIRED_AFTER_MAX_PASSES")), file.path(ref, "refinement_completed.rds"))
write.csv(data.frame(numerical_checks_ok = FALSE, convergence_code = 0,
  pdHess = TRUE, max_abs_gradient = .002998022), file.path(ref, "fit_checks_refined.csv"), row.names = FALSE)
selected <- e$select_environment(tmp)
check(!selected$ready[2] && selected$kind[2] == "refined", "Silently fell back from bad refinement")
saveRDS(list(summary = data.frame(analysis = "offset_multivoltine", refinement_status = "FIT_OK_REFINED")),
        file.path(ref, "refinement_completed.rds"))
write.csv(data.frame(numerical_checks_ok = TRUE, convergence_code = 0,
  pdHess = TRUE, max_abs_gradient = .0001), file.path(ref, "fit_checks_refined.csv"), row.names = FALSE)
check(all(e$select_environment(tmp)$ready), "Accepted refinement not selected")
unlink(file.path(ref, "refinement_completed.rds"))
check(!e$select_environment(tmp)$ready[2], "Incomplete refinement was accepted")
cat("PASS: original/refined source selection and strict numerical gates\n")

# 3b. Mixed v1/v2 queues, pending v2 work, provenance, and orphan child gates.
write_refinement <- function(path, analysis, accepted = TRUE) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  file.create(file.path(path, "fit_refined_final.rds"))
  saveRDS(list(summary = data.frame(analysis = analysis,
    refinement_status = if (accepted) "FIT_OK_REFINED" else
      "NUMERICAL_REVIEW_REQUIRED_AFTER_MAX_PASSES")),
    file.path(path, "refinement_completed.rds"))
  write.csv(data.frame(numerical_checks_ok = accepted, convergence_code = 0L,
    pdHess = TRUE, max_abs_gradient = if (accepted) .0001 else .002998022),
    file.path(path, "fit_checks_refined.csv"), row.names = FALSE)
}
write_refinement(ref, "offset_multivoltine")
onset_v1 <- file.path(tmp, "postfit_refinement_queue_v1", "onset_multivoltine")
write_refinement(onset_v1, "onset_multivoltine")
v2 <- file.path(tmp, "postfit_refinement_queue_v2_isolated")
dir.create(v2)
selected <- e$select_environment(tmp)
check(identical(selected$refinement_version, c("v2", "v1", "v2", "v2")),
      "Pending v2 queue did not take precedence over v1/original")
check(identical(selected$ready, c(FALSE, TRUE, FALSE, FALSE)),
      "Scheduled but unstarted v2 model was accepted")
check(all(selected$status[c(1, 3, 4)] == "REFINEMENT_NOT_STARTED"),
      "Pending v2 status is not explicit")
for (a in e$V2_ANALYSES) write_refinement(file.path(v2, a), a)
selected <- e$select_environment(tmp)
check(all(selected$ready), "Mixed accepted v1/v2 sources were rejected")
check(selected$fit[2] == file.path(ref, "fit_refined_final.rds"),
      "Carried-forward offset multivoltine lost its v1 source")
check(selected$fit[1] == file.path(v2, "onset_multivoltine", "fit_refined_final.rds"),
      "Old onset v1 took precedence over accepted v2")

bad_v2 <- file.path(v2, "onset_multivoltine")
write_refinement(bad_v2, "onset_multivoltine", accepted = FALSE)
check(!e$select_environment(tmp)$ready[1], "Failed v2 fell back to accepted v1")
write_refinement(bad_v2, "onset_multivoltine")
for (missing_file in c("fit_refined_final.rds", "refinement_completed.rds",
                       "fit_checks_refined.csv")) {
  unlink(file.path(bad_v2, missing_file))
  check(!e$select_environment(tmp)$ready[1], paste("Accepted v2 without", missing_file))
  write_refinement(bad_v2, "onset_multivoltine")
}
write_refinement(bad_v2, "onset_univoltine")
fails(e$select_environment(tmp), "identity mismatch")
write_refinement(bad_v2, "onset_multivoltine")
write.csv(data.frame(numerical_checks_ok = TRUE, convergence_code = 0L,
  pdHess = TRUE, max_abs_gradient = .001), file.path(bad_v2, "fit_checks_refined.csv"),
  row.names = FALSE)
check(!e$select_environment(tmp)$ready[1], "Exact gradient threshold was relaxed")
write_refinement(bad_v2, "onset_multivoltine")

# The isolated worker can report OK while the actual refinement still fails.
saveRDS(list(status = "OK"), file.path(bad_v2, "isolated_worker_result.rds"))
unlink(file.path(bad_v2, "refinement_completed.rds"))
check(!e$select_environment(tmp)$ready[1], "Worker completion replaced numerical acceptance")
write_refinement(bad_v2, "onset_multivoltine")

old_ref <- file.path(tmp, "configured_v1")
write_refinement(file.path(old_ref, "offset_multivoltine"), "offset_multivoltine")
saveRDS(list(source_run = tmp, old_refinement = old_ref),
        file.path(v2, "refinement_config.rds"))
check(e$select_environment(tmp)$fit[2] ==
      file.path(old_ref, "offset_multivoltine", "fit_refined_final.rds"),
      "V2 configured carried-forward v1 path ignored")
saveRDS(list(source_run = v2, old_refinement = old_ref),
        file.path(v2, "refinement_config.rds"))
fails(e$select_environment(tmp), "different environmental run")
saveRDS(list(source_run = tmp, old_refinement = old_ref),
        file.path(v2, "refinement_config.rds"))

# Process liveness is simulated; no workers are launched or terminated.
real_alive <- e$alive
e$alive <- function(owner) isTRUE(owner$test_alive)
owner_paths <- c(file.path(tmp, "worker_owner.rds"),
  file.path(tmp, "postfit_refinement_queue_v1", "worker_owner.rds"),
  file.path(old_ref, "worker_owner.rds"), file.path(v2, "worker_owner.rds"),
  file.path(v2, "offset_univoltine", "isolated_worker_owner.rds"))
for (owner_path in owner_paths) {
  saveRDS(list(pid = 42L, create_time = 1, test_alive = TRUE), owner_path)
  workers <- e$environmental_workers(tmp)
  check(any(workers$active) && workers$owner_file[workers$active] == owner_path,
        paste("Active worker missed:", owner_path))
  fails(e$require_environment_idle(tmp), "worker is active")
  fails(e$child_stage(list(environment_run = tmp), "extract_one", "onset_univoltine"),
        "worker is active")
  saveRDS(list(pid = 42L, create_time = 1, test_alive = FALSE), owner_path)
}
check(!any(e$environmental_workers(tmp)$active), "Dead owners block the workflow")
e$require_environment_idle(tmp)

# Accepted models must remain blocked by a live child with a dead controller.
child_owner <- tail(owner_paths, 1)
saveRDS(list(pid = 42L, create_time = 1, test_alive = TRUE), child_owner)
demo <- file.path(tmp, "demography")
dir.create(demo)
saveRDS(list(source_run = tmp), file.path(demo, "run_config.rds"))
for (g in e$GROUPS) {
  file.create(file.path(demo, paste0("prepared_mundlak_", g, ".rds")))
  dir.create(file.path(demo, paste0(g, "__mundlak")))
  file.create(file.path(demo, paste0(g, "__mundlak"), "fit.rds"))
}
invisible(capture.output(pf <- e$preflight(tmp, tmp, demo)))
check(all(pf$environment$ready) && pf$active && !pf$ready,
      "Preflight did not block accepted sources while an orphan v2 child is active")
real_packages <- e$packages
e$packages <- function(names) invisible(TRUE)
fails(e$run(project_root = tmp, environment_run = tmp, demography_run = demo,
            mode = "prepare"), "Sources are not ready")
e$packages <- real_packages
check(!dir.exists(file.path(tmp, "output")), "Blocked run created output or launched work")

saveRDS(list(pid = 42L, create_time = 1, test_alive = FALSE), child_owner)
selected <- e$select_environment(tmp)
write_refinement(bad_v2, "onset_multivoltine", accepted = FALSE)
fails(e$resources(list(environment_run = tmp, environment = selected)),
      "source selection changed")
write_refinement(bad_v2, "onset_multivoltine")
saveRDS(list(pid = 42L), child_owner)
fails(e$environmental_workers(tmp), "Invalid worker ownership")
writeLines("invalid RDS", child_owner)
fails(e$environmental_workers(tmp), "unknown input format")
unlink(child_owner)
e$alive <- real_alive
cat("PASS: mixed v1/v2, pending/failed/missing v2, provenance, controller/child guards, blocked launch\n")

if (requireNamespace("ps", quietly = TRUE)) {
  owner <- list(pid = Sys.getpid(), create_time = ps::ps_create_time(ps::ps_handle()))
  check(e$alive(owner), "Current process not detected")
  owner$create_time <- owner$create_time - 10
  check(!e$alive(owner), "Reused PID with different creation time accepted")
  check(!e$alive(list(pid = 2147483647L, create_time = 1)), "Missing PID treated as alive")
  cat("PASS: native process identity and missing PID handling\n")
} else cat("SKIP: native ps checks (package not accessible); synthetic worker gates passed\n")

# 4. Export integration and recovery: new predictors must not produce an old/new LRT.
out <- tempfile(); dir.create(out); dir.create(file.path(out, "univoltine__mundlak"))
cfg <- list(out = out, progress = file.path(out, "progress.rds"))
p$predictors_updated <- TRUE; p$signature <- "new-sig"
fit$signature <- p$signature
status <- engine$export_fitted(cfg, "univoltine", p, fit)
check(status == "FIT_OK", "Exporter status wrong")
needed <- c("fixed_effects.csv", "key_within_between_effects.csv", "joint_Wald_tests.csv",
            "onset_offset_balance_within_between.csv", "plasticity_effects_per_component_SD.csv")
check(all(file.exists(file.path(out, "univoltine__mundlak", needed))), "Missing exports")
check(!file.exists(file.path(out, "univoltine__mundlak", "nested_within_between_LRT.csv")),
      "Invalid old/new nested test exported")
fx <- read.csv(file.path(out, "univoltine__mundlak", "fixed_effects.csv"))
unlink(file.path(out, "univoltine__mundlak", "fixed_effects.csv"))
engine$export_fitted(cfg, "univoltine", p, fit)
check(identical(fx, read.csv(file.path(out, "univoltine__mundlak", "fixed_effects.csv"))),
      "Export recovery changed results")
# The original unchanged-predictor export path still performs its valid nested test.
legacy <- p; legacy$predictors_updated <- NULL
legacy$source_fit_beta <- setNames(rep(.01, 6), engine$source_fixed_terms)
legacy$source_fit_V <- diag(.01, 6)
legacy$source_fit_checks <- data.frame(logLik = -110, parameter_df = 20, AIC = 260)
fit$checks$logLik <- -100; fit$checks$parameter_df <- 28; fit$checks$AIC <- 256
engine$export_fitted(cfg, "univoltine", legacy, fit)
check(file.exists(file.path(out, "univoltine__mundlak", "nested_within_between_LRT.csv")),
      "Legacy export regression")
cat("PASS: new/legacy export integration and interrupted-export recovery\n")

# 5. New conditional diagnostic rows use the updated fit and AR1 links.
fit$eta <- rep(0, nrow(p$frame)); fit$parameters <- p$parameters
diag <- e$diagnostic_tables(p, fit)
check(diag$distribution$n == nrow(p$frame) && diag$temporal$n_pairs[2] > 0,
      "Diagnostic row/gap mapping incorrect")
check(all(is.finite(diag$site_year$residual)), "Non-finite residuals")
fails(e$diagnostic_tables(p, within(fit, signature <- "wrong")), "different fits")
cat("PASS: conditional residuals and latent innovation diagnostics\n")
cleanup_paths <- normalizePath(c(tmp, out), winslash = "/", mustWork = TRUE)
cleanup_parent <- paste0(normalizePath(tempdir(), winslash = "/", mustWork = TRUE), "/")
check(all(startsWith(cleanup_paths, cleanup_parent)), "Test cleanup escaped temporary directory")
unlink(cleanup_paths, recursive = TRUE)
cat("ALL POST-ENVIRONMENT CONTRACT TESTS PASSED\n")
