# Compare the original pheno_plasticity.R model with the selected final model.
# Saved estimates only: NO refitting, optimization, predict() or likelihood calls.
# Focus: fixed anomaly slopes at matched observed population contexts.
# source('compare_saved_offset_models.R')
# offset_comparison <- compare_saved_offset_models()
compare_saved_offset_models <- function(
  project_root = "E:/phenoIMPACT project/code/phenoIMPACT",
  environment_run = file.path(project_root, "output", "phenology_plasticity",
    "environment_mundlak_by_voltinism_qc_final_15km", "run_de0890c85d4a"),
  old_bundle = file.path(project_root, "output", "phenology_plasticity", "plasticity_main_models.rds")) {
  need <- function(x, msg) if (!isTRUE(x)) stop(msg, call. = FALSE)
  need(requireNamespace("lme4", quietly = TRUE), "lme4 is required.")
  need(file.exists(old_bundle), paste("Original bundle missing:", old_bundle))
  e <- new.env(parent = globalenv())
  sys.source(file.path(project_root, "post_environment_mundlak.R"), e)
  pf <- e$pheno_post_env$preflight(project_root = project_root, environment_run = environment_run)
  need(!pf$active, "Environmental worker active; retry after completion.")
  src <- pf$environment[pf$environment$analysis == "offset_multivoltine", , drop = FALSE]
  need(nrow(src) == 1 && isTRUE(src$ready), "Final multivoltine offset source not accepted.")
  out <- file.path(environment_run, paste0("saved_offset_comparison_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  need(!dir.exists(out), "Output exists; retry in one second.")
  dir.create(out, recursive = TRUE)
  csv <- function(x, n) write.csv(x, file.path(out, n), row.names = FALSE)
  mods <- c(photoperiod = "photo_tw90", background = "clim_background_tw90",
            predictability = "clim_predictability_tw90", trend = "clim_trend_tw90")
  an <- "clim_anomaly_tw90"
  key <- function(d) paste(d$SPECIES, d$SITE_ID, sep = "\r")
  population_table <- function(d, variables) {
    z <- unique(d[c("SPECIES", "SITE_ID", variables)])
    need(!anyDuplicated(key(z)), "Environmental context varies within population.")
    z
  }
  # Release the historical bundle before loading the much larger final fit.
  bundle <- readRDS(old_bundle)
  old <- bundle$primary$offset_multivoltine
  need(inherits(old$model, "merMod"), "Original offset lmer model not found.")
  need(old$best_window == 90, "Original selected window is not 90 days; comparison needs adaptation.")
  old_form <- lme4::nobars(formula(old$model))
  old_b <- lme4::fixef(old$model); old_V <- as.matrix(vcov(old$model))
  mf <- model.frame(old$model)
  ix <- match(rownames(mf), rownames(old$data))
  need(!anyNA(ix), "Cannot align historical fitted observations with saved data.")
  old_d <- old$data[ix, , drop = FALSE]
  for (v in intersect(names(mf), names(old_d)))
    need(isTRUE(all.equal(mf[[v]], old_d[[v]], check.attributes = FALSE)),
         paste("Historical fitted-data mismatch:", v))
  need(all(unname(mods) %in% names(old_d)), "Historical moderator names differ.")
  old_pop <- population_table(old_d, unname(mods))
  old_n <- nrow(old_d)
  old_keys <- if ("YEAR" %in% names(old_d))
    data.frame(pop_key = key(old_d), YEAR = old_d$YEAR) else NULL
  old_formula_text <- paste(deparse(old_form), collapse = " ")
  rm(bundle, old, mf, old_d); gc()

  p <- readRDS(src$prepared)
  within <- paste0(names(mods), "_within_species")
  between <- paste0(names(mods), "_species_mean")
  final_pop <- population_table(p$data, c(within, between))
  final_n <- nrow(p$data)
  final_keys <- data.frame(pop_key = key(p$data), YEAR = p$data$YEAR)
  final_form <- lme4::nobars(p$formula)
  saved <- readRDS(src$fit)
  fit <- saved$fit
  need(inherits(fit, "sdmTMB"), "Selected final fit is invalid.")
  marker <- readRDS(src$original_marker)
  sig <- if (src$kind == "refined") saved$source_signature else saved$signature
  need(!is.null(sig) && identical(sig, marker$signature), "Final source signature mismatch.")
  need(nrow(fit$data) == nrow(p$data), "Final fitted/prepared rows differ.")
  for (v in all.vars(final_form))
    need(isTRUE(all.equal(fit$data[[v]], p$data[[v]], check.attributes = FALSE)),
         paste("Final fitted/prepared data mismatch:", v))
  need(identical(names(fit$model$par), names(fit$sd_report$par.fixed)) &&
    isTRUE(all.equal(as.numeric(fit$model$par), as.numeric(fit$sd_report$par.fixed), tolerance = 1e-6)),
    "Final covariance is not synchronized with point estimates.")
  terms <- colnames(model.matrix(final_form, p$data[1, , drop = FALSE]))
  ix <- which(names(fit$sd_report$par.fixed) == "b_j")
  need(length(ix) == length(terms), "Cannot align final coefficients.")
  final_b <- setNames(as.numeric(fit$sd_report$par.fixed[ix]), terms)
  final_V <- as.matrix(fit$sd_report$cov.fixed[ix, ix, drop = FALSE])
  dimnames(final_V) <- list(terms, terms)
  rm(saved, fit, p); gc()

  common <- intersect(key(old_pop), key(final_pop))
  need(length(common) > 0, "No common populations.")
  old_context <- old_pop[match(common, key(old_pop)), , drop = FALSE]
  final_context <- final_pop[match(common, key(final_pop)), , drop = FALSE]
  # Validate that the reconstructed original moderator scales are identical.
  differences <- lapply(seq_along(mods), function(j) {
    recovered <- final_context[[within[j]]] + final_context[[between[j]]]
    delta <- recovered - old_context[[mods[j]]]
    data.frame(moderator = names(mods)[j], max_abs_difference = max(abs(delta)))
  })
  differences <- do.call(rbind, differences); csv(differences, "moderator_scale_check.csv")
  need(all(differences$max_abs_difference < 1e-7),
       "Historical and final environmental scales differ; scale-check CSV saved. Do not compare slopes yet.")
  slope <- function(context, form, b, V) {
    tt <- delete.response(terms(form))
    context[[an]] <- 0; context$ONSET_mean_z <- 0
    x0 <- model.matrix(tt, context)
    context[[an]] <- 1; x1 <- model.matrix(tt, context)
    L <- (x1-x0)[, names(b), drop = FALSE]
    need(all(is.finite(V)) && all(is.finite(b)), "Invalid coefficients/covariance.")
    l <- colMeans(L); se <- sqrt(drop(t(l)%*%V%*%l))
    list(values = as.vector(L%*%b), mean = sum(l*b), SE = se)
  }
  a <- slope(old_context, old_form, old_b, old_V)
  b <- slope(final_context, final_form, final_b, final_V)
  summary_row <- function(x, name, n) data.frame(model = name, fitted_rows = n,
    matched_populations = length(common), mean_fixed_slope = x$mean,
    CI_low = x$mean-1.96*x$SE, CI_high = x$mean+1.96*x$SE,
    fraction_positive_fixed = mean(x$values > 0))
  summary <- rbind(summary_row(a, "original", old_n), summary_row(b, "final", final_n))
  csv(summary, "matched_population_summary.csv")
  csv(data.frame(SPECIES = old_context$SPECIES, SITE_ID = old_context$SITE_ID,
    original_fixed_slope = a$values, final_fixed_slope = b$values,
    difference = b$values-a$values), "matched_population_slopes.csv")
  # Slope at the zero used by the ORIGINAL plotting function, for each model.
  csv(data.frame(model = c("original", "final"),
    anomaly_coefficient_at_all_moderators_zero = c(old_b[an], final_b[an]),
    SE = c(sqrt(old_V[an,an]), sqrt(final_V[an,an]))), "zero_reference_slopes.csv")
  csv(src, "selected_final_source.csv")
  overlap <- if (!is.null(old_keys)) nrow(merge(unique(old_keys), unique(final_keys),
                                               by = c("pop_key", "YEAR"))) else NA_integer_
  csv(data.frame(original_populations = nrow(old_pop), final_populations = nrow(final_pop),
                 matched_populations = length(common), shared_population_years = overlap), "sample_overlap.csv")
  paths <- c(old_bundle, src$prepared, src$fit)
  info <- file.info(paths)
  csv(data.frame(source = paths, size_bytes = info$size, modified = as.character(info$mtime)),
      "input_manifest.csv")
  writeLines(c("SAVED-MODEL COMPARISON, NO REFITTING.",
    paste("Original formula:", old_formula_text),
    paste("Final formula:", paste(deparse(final_form), collapse = " ")),
    "Matched populations, identical observed raw environmental contexts, equal weight per population.",
    "Final within/between decomposition retained from its actual fitted sample; not recalculated on overlap.",
    "Slopes are fixed-effect partial derivatives; no species random slopes included.",
    "Final offset slopes condition on onset. This comparison does not isolate individual model changes.",
    "CI for each mean uses its full fixed covariance. No formal CI/test for differences between models.",
    "Units: days per inherited standardized 90-day anomaly, not days per degree C.",
    "The anomaly coefficient at all moderators zero is a reference slope, not an average population response.",
    "Historical/final fits may use different samples and random/spatial structures. Inspect sample_overlap.csv.",
    "No new inference on overall species responses or generation numbers is supplied by this diagnostic."
  ), file.path(out, "README.txt"))
  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
  if (requireNamespace("zip", quietly = TRUE)) zip::zipr(
    file.path(out, "saved_offset_comparison_review.zip"),
    list.files(out, pattern = "\\.(csv|txt)$", full.names = TRUE), root = out)
  print(summary, row.names = FALSE)
  message("Output: ", normalizePath(out, winslash = "/"))
  invisible(list(summary = summary, output_dir = out))
}
