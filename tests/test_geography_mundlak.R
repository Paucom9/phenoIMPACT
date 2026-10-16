# Run from the repository root:
# Rscript tests/test_geography_mundlak.R
# Base-R tests do not launch processes or require the real project data.
# Formula/CRS checks additionally run when lme4/sf are installed.
source("geography_mundlak_by_voltinism_qc_final.R", encoding = "UTF-8")
e <- environment(pheno_geo_mundlak$run)
n_checks <- 0L
ok <- function(x) {
  stopifnot(isTRUE(x))
  n_checks <<- n_checks + 1L
}
same <- function(a, b, tolerance = 1e-10) {
  comparison <- all.equal(a, b, tolerance = tolerance, check.attributes = FALSE)
  if (!isTRUE(comparison)) stop(paste("Comparison", n_checks + 1L,
    deparse(substitute(a)), "vs", deparse(substitute(b)), ":", paste(comparison, collapse = "; ")))
  ok(TRUE)
}
fails <- function(expr, pattern) {
  msg <- tryCatch({ force(expr); "NO ERROR" }, error = function(err) conditionMessage(err))
  ok(grepl(pattern, msg, fixed = TRUE))
}
ok(nrow(pheno_geo_mundlak$plan()) == 4L)
ok(!exists("pheno_env_mundlak", inherits = FALSE))
ok(is.null(getOption("phenoimpact.geo.mundlak.job")))

set.seed(42)
pop <- data.frame(SPECIES = rep(paste0("sp", 1:8), each = 6), SITE_ID = paste0("site", 1:48))
sites <- data.frame(SITE_ID = pop$SITE_ID,
  latitude_deg = 35 + 1.7 * rep(1:8, each = 6) + rep(seq(-1.5, 1.5, length.out = 6), 8),
  longitude_deg = -5 + .4 * rep((1:8)^2, each = 6) + rnorm(48))
annual <- pop[rep(seq_len(48), times = rep(1:8, length.out = 48)), ]
annual$YEAR <- ave(seq_len(nrow(annual)), annual$SITE_ID, FUN = seq_along) + 2000L
annual$row_marker <- seq_len(nrow(annual))
d <- e$decompose_coordinates(annual, sites)
once <- e$decompose_coordinates(pop, sites)
same(d$populations, once$populations)
same(d$data[names(annual)], annual)
ok(!anyDuplicated(d$populations[c("SPECIES", "SITE_ID")]))
for (axis in c("latitude", "longitude")) {
  p <- d$populations
  centre <- d$scaling$centre_degrees[d$scaling$axis == axis]
  w <- p[[paste0(axis, "_within_species_10deg")]]
  b <- p[[paste0(axis, "_species_mean_10deg")]]
  same(10 * (w + b) + centre, p[[paste0(axis, "_deg")]])
  same(as.numeric(tapply(w, p$SPECIES, mean)), rep(0, 8))
  raw <- sites[[paste0(axis, "_deg")]]
  means <- tapply(raw, pop$SPECIES, mean)
  same(centre, mean(means))
}
singleton <- e$decompose_coordinates(pop[1, ], sites)
same(singleton$populations$latitude_within_species_10deg, 0)
same(singleton$populations$longitude_within_species_10deg, 0)
fails(e$decompose_coordinates(annual, sites[-1, ]), "Missing coordinates")
fails(e$decompose_coordinates(annual, rbind(sites, sites[1, ])), "Duplicate site")
bad <- sites; bad$latitude_deg[1] <- NA_real_
fails(e$decompose_coordinates(annual, bad), "Invalid latitude")

have_lme4 <- requireNamespace("lme4", quietly = TRUE)
for (analysis in e$ANALYSES) {
  s <- e$expected_spec(analysis)
  z <- d$data
  z[[s$anomaly]] <- rnorm(nrow(z))
  z$ONSET_mean_z <- rnorm(nrow(z))
  z$SPECIES_slope <- factor(z$SPECIES)
  z$site_year_id <- interaction(z$SITE_ID, z$YEAR, drop = TRUE)
  rhs <- paste0(s$anomaly, " * (", paste(e$GEO_VARS, collapse = "+"), ")")
  if (length(s$control)) rhs <- paste(rhs, "+ ONSET_mean_z")
  fixed <- as.formula(paste("~", rhs))
  X <- model.matrix(fixed, z)
  beta <- setNames(seq_len(ncol(X)) / 20, colnames(X))
  beta["(Intercept)"] <- 150
  z[[s$response]] <- as.numeric(X %*% beta)
  same(unname(lm.fit(X, z[[s$response]])$coefficients), unname(beta), tolerance = 1e-8)
  ok(ncol(X) == if (length(s$control)) 11L else 10L)

  if (have_lme4) {
    source_formula <- as.formula(paste(s$response, "~", s$anomaly,
      "* (photoperiod_within_species + photoperiod_species_mean)",
      if (length(s$control)) "+ ONSET_mean_z" else "",
      "+ (1 | SITE_ID) + (1 | SPECIES) + (0 +", s$anomaly,
      "| SPECIES_slope) + (1 | site_year_id)"))
    f <- e$make_formula(source_formula, s)
    same(e$design_matrix(f, z), X)
    same(e$random_terms_text(f), e$random_terms_text(source_formula))
    ok(!any(grepl("photoperiod", all.vars(f))))
    incomplete <- z; incomplete[[s$anomaly]][1] <- NA_real_
    fails(e$design_matrix(f, incomplete), "missing values")
    singular <- z; singular$longitude_within_species_10deg <- singular$latitude_within_species_10deg
    fails(e$design_matrix(f, singular), "rank deficient")
  }

  V <- diag(.04, length(beta)) + matrix(.003, length(beta), length(beta))
  dimnames(V) <- list(names(beta), names(beta))
  results <- e$geographic_results(beta, V, s, analysis, d$populations)
  terms <- e$interaction_terms(beta, s$anomaly)
  same(results$effects$estimate, unname(beta[terms]))
  same(results$effects$estimate_biological, unname(s$orientation * beta[terms]))
  same(results$effects$effect_scale, c("within species", "between species", "within species", "between species"))
  same(results$joint$df, c(2L, 2L, 4L))
  same(results$differences$estimate[1], unname(beta[terms[1]] - beta[terms[2]]))
  same(results$differences$std.error[1]^2,
       unname(V[terms[1], terms[1]] + V[terms[2], terms[2]] - 2 * V[terms[1], terms[2]]))
  C <- setNames(numeric(length(beta)), names(beta))
  C[s$anomaly] <- 1
  C[terms] <- as.numeric(d$populations[1, e$GEO_VARS])
  same(results$slopes$estimate[1], sum(C * beta))
  same(results$slopes$std.error[1]^2, as.numeric(t(C) %*% V %*% C))
  ok(all(results$effects$conf.low_biological <= results$effects$conf.high_biological))
  # Finite difference in the annual-phenology linear predictor agrees with the
  # exported derivative, including the additive onset covariate for offsets.
  z1 <- z[1, ]; z2 <- z1; z2[[s$anomaly]] <- z1[[s$anomaly]] + 1
  delta <- as.numeric((model.matrix(fixed, z2) - model.matrix(fixed, z1)) %*% beta)
  same(delta, results$slopes$estimate[1])
  Vbad <- V[rev(seq_along(beta)), rev(seq_along(beta))]
  L <- matrix(C, nrow = 1, dimnames = list(NULL, names(beta)))
  fails(e$linear_estimates(L, beta, Vbad), "not aligned")
  Vsingular <- matrix(1, length(beta), length(beta), dimnames = dimnames(V))
  singular_results <- e$geographic_results(beta, Vsingular, s, analysis, d$populations)
  ok(all(is.na(singular_results$joint$p_Wald)))
}
if (!have_lme4) cat("SKIP: lme4 formula/design checks (package unavailable).\n")

# Source-readiness checks only read small metadata. The prepared-input placeholder
# is deliberately NOT an RDS; source_status must never try to load it.
root <- tempfile("geo_test_")
dir.create(root)
env <- file.path(root, "environment")
dir.create(env)
saveRDS(list(version = e$SOURCE_VERSION, signature = "source_sig"), file.path(env, "run_config.rds"))
for (a in e$ANALYSES) {
  dir.create(file.path(env, a))
  writeLines("not a serialized model", file.path(env, a, "prepared_input.rds"))
}
pending <- e$source_status(env)
ok(all(pending$table$status == "NOT_COMPLETED"))
fails(e$source_status(env, TRUE), "Wait until all four")
fails(pheno_geo_mundlak$run(project_root = root, environment_run = env), "Wait until all four")
for (a in e$ANALYSES) saveRDS(list(signature = "source_sig", status = "FIT_OK"),
                              file.path(env, a, "completed.rds"))
ok(all(e$source_status(env, TRUE)$table$status == "FIT_OK"))
saveRDS(list(signature = "wrong", status = "FIT_OK"), file.path(env, e$ANALYSES[1], "completed.rds"))
fails(e$source_status(env), "signature mismatch")

cfg <- list(out = file.path(root, "geography"), signature = "geo_sig")
dir.create(cfg$out)
a <- e$ANALYSES[1]
out <- file.path(cfg$out, a)
dir.create(out)
done <- list(version = e$VERSION, signature = cfg$signature, analysis = a, status = "FIT_OK")
saveRDS(done, file.path(out, "completed.rds"))
fails(e$completed_cache_status(cfg, a), "Incomplete completed cache")
file.create(file.path(out, c("prepared_input.rds", "fit_geography.rds")))
ok(is.null(e$completed_cache_status(cfg, a)))
required <- c("geography_within_between_effects.csv", "within_minus_between_contrasts.csv",
  "joint_geography_tests.csv", "population_fixed_geographic_sensitivity.csv",
  "fixed_effects.csv", "fixed_effect_covariance.csv", "fit_checks.csv")
file.create(file.path(out, required))
same(e$completed_cache_status(cfg, a), "FIT_OK")
done$signature <- "wrong"; saveRDS(done, file.path(out, "completed.rds"))
fails(e$completed_cache_status(cfg, a), "signature/analysis mismatch")

# Serialization preserves the worker's numerical helpers and four-model settings.
engine <- e$dump_engine(root)
other <- new.env(parent = globalenv())
sys.source(engine, other)
same(other$plan(), pheno_geo_mundlak$plan())
same(other$decompose_coordinates(annual, sites)$populations, d$populations)

if (requireNamespace("sf", quietly = TRUE)) {
  ll <- sf::st_as_sf(sites, coords = c("longitude_deg", "latitude_deg"), crs = 4326)
  xy <- sf::st_coordinates(sf::st_transform(ll, 3035))
  projected <- data.frame(SITE_ID = sites$SITE_ID, x_3035 = xy[, 1], y_3035 = xy[, 2],
                           x_km = xy[, 1] / 1000, y_km = xy[, 2] / 1000)
  recovered <- e$recover_coordinates(projected)
  same(recovered$latitude_deg, sites$latitude_deg, tolerance = 1e-7)
  same(recovered$longitude_deg, sites$longitude_deg, tolerance = 1e-7)
  bad <- projected; bad$x_km[1] <- bad$x_km[1] + 1
  fails(e$recover_coordinates(bad), "metre/km coordinates disagree")
} else cat("SKIP: sf coordinate transformation checks (package unavailable).\n")
unlink(root, recursive = TRUE)
cat("PASS:", n_checks, "checks. No spatial fit or real-data process was launched.\n")
