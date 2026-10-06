# Population trends using the extracted 15-km SPATIAL PLASTICITY estimates.
# Adaptation of Pau Colom's population_trend_models.R (2026-05-13).
# The abundance models retain the original Gamma/log + species/site slopes +
# population AR(1). They do NOT include a new spatial field for abundance.
# Run from the phenoIMPACT R project: source(file.choose()). Keep R open.
# 2 full models + 4 reduced fits for LRTs (run_lrt = TRUE, as in the original).
# Models are saved immediately, and compatible checkpoints are reused.
# Escape pauses only the display; the background queue continues.
# To configure before running: options(pheno.trends.functions_only = TRUE).

population_trend_config <- list(
  project_root = NULL,            # NULL = here::here()
  plasticity_file = NULL,         # NULL = latest COMPLETED 15-km extraction
  abundance_file = NULL,          # NULL = output/pheno_abund_estimates_allspp.csv
  output_root = NULL,             # NULL = output/population_trends/spatial_plasticity_15km
  run_lrt = TRUE,                 # FALSE = only the two full models + Wald tests
  prepare_only = FALSE,
  make_figures = TRUE,
  font_family = "Garamond",
  monitor_seconds = 30,
  threads = 1L,                  # One fit at a time; conservative RAM use.
  optimizer_trace = 1L,
  abundance_quantile = 0.999
)

pheno_population_trends <- local({
assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
write_table <- function(x, path) utils::write.csv(x, path, row.names = FALSE, na = "")
formula_text <- function(x) paste(deparse(x, width.cutoff = 500L), collapse = " ")
atomic_save <- function(x, path) {
  tmp <- tempfile("writing_", tmpdir = dirname(path)); on.exit(unlink(tmp), add = TRUE)
  saveRDS(x, tmp)
  backup <- paste0(path, ".previous")
  old <- file.exists(path)
  if (old && !file.copy(path, backup, overwrite = TRUE)) stop("Cannot back up ", path)
  if (old && !file.remove(path)) stop("Cannot replace ", path)
  if (!file.rename(tmp, path)) {
    if (old) file.copy(backup, path, overwrite = TRUE)
    stop("Cannot finalize ", path)
  }
  if (old) unlink(backup)
  invisible(path)
}
set_phase <- function(cfg, phase, task = "", task_started = Sys.time()) {
  atomic_save(list(phase = phase, task = task, task_started = task_started, updated = Sys.time()), cfg$progress_file)
  message("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", task, " ", phase)
}
zscale <- function(x) {
  s <- stats::sd(x)
  assert(length(x) > 1L && is.finite(s) && s > 0, "Plasticity is constant or invalid within a group.")
  as.numeric(scale(x))
}
model_formula <- function(variant = "full") {
  f <- ABUND_INDEX ~ year_decade * onset_plasticity_z + year_decade * offset_plasticity_bio_z +
    (1 + year_decade || SPECIES) + (1 + year_decade || SITE_ID) + ar1(year_fac + 0 | pop_id)
  if (variant == "without_onset_interaction") f <- stats::update(f, . ~ . - year_decade:onset_plasticity_z)
  if (variant == "without_offset_interaction") f <- stats::update(f, . ~ . - year_decade:offset_plasticity_bio_z)
  assert(variant %in% c("full", "without_onset_interaction", "without_offset_interaction"), "Unknown model variant.")
  # Avoid retaining the preparation/worker environment in saved formulas.
  environment(f) <- .GlobalEnv
  f
}
latest_extraction <- function(project) {
  root <- file.path(project, "output", "phenology_plasticity", "spatial_models")
  dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[startsWith(basename(dirs), "population_plasticity_15km_")]
  name <- "phenological_population_plasticity_contextual_spatial_15km.csv"
  dirs <- dirs[vapply(dirs, function(d) {
    summary <- file.path(d, "extraction_summary.csv")
    if (!all(file.exists(c(summary, file.path(d, name))))) return(FALSE)
    s <- tryCatch(utils::read.csv(summary), error = function(e) NULL)
    !is.null(s) && nrow(s) == 3L && all(s$status == "EXTRACTED") &&
      setequal(s$analysis, c("onset", "offset_univoltine", "offset_multivoltine"))
  }, logical(1))]
  assert(length(dirs) > 0L, "No completed 15-km extraction found. Set population_trend_config$plasticity_file explicitly.")
  paths <- file.path(dirs, name)
  paths[order(file.info(paths)$mtime, paths, decreasing = TRUE)][1L]
}

prepare_tables <- function(abundance, plasticity, quantile_probability = 0.999) {
  assert(all(c("SPECIES", "SITE_ID", "YEAR", "ABUND_INDEX") %in% names(abundance)), "Annual abundance columns are missing.")
  required <- c("pop_id", "SPECIES", "SITE_ID", "voltinism", "onset_advancement_plasticity_contextual",
    "offset_termination_plasticity_contextual")
  assert(all(required %in% names(plasticity)), "Spatial plasticity columns are missing.")
  assert(is.numeric(abundance$ABUND_INDEX) && !any(is.infinite(abundance$ABUND_INDEX)), "Abundance must be numeric and finite (NA is allowed before filtering).")
  threshold <- as.numeric(stats::quantile(abundance$ABUND_INDEX, quantile_probability, na.rm = TRUE))
  assert(is.finite(threshold) && threshold > 0, "Invalid abundance quantile threshold.")
  retain <- !is.na(abundance$ABUND_INDEX) & abundance$ABUND_INDEX > 0 & abundance$ABUND_INDEX <= threshold
  filter_summary <- data.frame(n_total = nrow(abundance), n_removed = sum(!retain),
    prop_removed = mean(!retain), abundance_quantile = quantile_probability, abundance_threshold = threshold)
  # Preserve original ordering: compute quantile before distinct()/joins.
  d <- unique(abundance[retain, c("SPECIES", "SITE_ID", "YEAR", "ABUND_INDEX")])
  d$SPECIES <- as.character(d$SPECIES); d$SITE_ID <- as.character(d$SITE_ID)
  yr <- suppressWarnings(as.numeric(as.character(d$YEAR)))
  assert(all(is.na(yr) | (is.finite(yr) & yr == floor(yr))), "YEAR contains non-integer or infinite values.")
  d$year_num <- as.integer(yr)
  d <- d[!is.na(d$SPECIES) & nzchar(d$SPECIES) & !is.na(d$SITE_ID) & nzchar(d$SITE_ID) & !is.na(d$year_num), ]
  assert(nrow(d) > 0L, "No valid annual abundance rows.")
  assert(!anyDuplicated(d[c("SPECIES", "SITE_ID", "year_num")]),
    "Conflicting abundance values for the same SPECIES/SITE_ID/YEAR. Resolve them before fitting; values were not averaged.")
  d$pop_id <- paste(d$SPECIES, d$SITE_ID, sep = "_")
  identities <- unique(d[c("SPECIES", "SITE_ID", "pop_id")])
  assert(!anyDuplicated(identities$pop_id), "Legacy population-ID collision in abundance data.")
  year_center <- mean(d$year_num) # Same pre-join reference as the original script.
  d$year_decade <- (d$year_num - year_center) / 10
  d$log_abund_index <- log(d$ABUND_INDEX)

  p <- as.data.frame(plasticity[, required])
  for (key in c("SPECIES", "SITE_ID", "pop_id", "voltinism")) p[[key]] <- as.character(p[[key]])
  assert(!anyNA(p[c("SPECIES", "SITE_ID", "pop_id", "voltinism")]) &&
    !anyDuplicated(p$pop_id) && !anyDuplicated(p[c("SPECIES", "SITE_ID")]), "Missing or duplicated plasticity keys.")
  assert(identical(p$pop_id, paste(p$SPECIES, p$SITE_ID, sep = "_")), "Plasticity IDs disagree with species/site keys.")
  assert(all(p$voltinism %in% c("univoltine", "multivoltine")) &&
    all(is.finite(p$onset_advancement_plasticity_contextual)) &&
    all(is.finite(p$offset_termination_plasticity_contextual)), "Invalid group or incomplete extracted plasticity.")
  p$offset_plasticity_bio <- ifelse(p$voltinism == "univoltine",
    -p$offset_termination_plasticity_contextual, p$offset_termination_plasticity_contextual)
  # Retain the original scaling population: all extracted populations in each
  # group, BEFORE matching/filtering abundance. Mean/SD are exported explicitly.
  scaling <- p |>
    dplyr::group_by(voltinism) |>
    dplyr::summarise(n_populations = dplyr::n(),
      onset_mean = mean(onset_advancement_plasticity_contextual), onset_SD = stats::sd(onset_advancement_plasticity_contextual),
      offset_bio_mean = mean(offset_plasticity_bio), offset_bio_SD = stats::sd(offset_plasticity_bio), .groups = "drop")
  p <- p |>
    dplyr::group_by(voltinism) |>
    dplyr::mutate(onset_plasticity_z = zscale(onset_advancement_plasticity_contextual),
      offset_plasticity_bio_z = zscale(offset_plasticity_bio)) |>
    dplyr::ungroup() |> as.data.frame()
  correlations <- p |>
    dplyr::group_by(voltinism) |>
    dplyr::summarise(n_populations = dplyr::n(),
      pearson_raw = stats::cor(onset_advancement_plasticity_contextual, offset_termination_plasticity_contextual),
      spearman_raw = stats::cor(onset_advancement_plasticity_contextual, offset_termination_plasticity_contextual, method = "spearman"),
      pearson_bio = stats::cor(onset_plasticity_z, offset_plasticity_bio_z),
      spearman_bio = stats::cor(onset_plasticity_z, offset_plasticity_bio_z, method = "spearman"), .groups = "drop")
  joined <- dplyr::left_join(d, p[, setdiff(names(p), "pop_id"), drop = FALSE], by = c("SPECIES", "SITE_ID"))
  assert(nrow(joined) == nrow(d), "Plasticity join multiplied abundance rows.")
  unmatched <- unique(joined[is.na(joined$voltinism), c("pop_id", "SPECIES", "SITE_ID")])
  joined <- joined[!is.na(joined$voltinism), , drop = FALSE]
  assert(nrow(joined) > 0L, "No abundance rows matched spatial plasticity.")
  unused_plasticity <- p[!p$pop_id %in% joined$pop_id, c("pop_id", "SPECIES", "SITE_ID", "voltinism"), drop = FALSE]
  groups <- list(); summaries <- list(); years <- list()
  for (group in c("univoltine", "multivoltine")) {
    g <- joined[joined$voltinism == group, , drop = FALSE]
    assert(nrow(g) > 0L, paste("No matched abundance data for", group))
    for (key in c("SPECIES", "SITE_ID", "pop_id")) g[[key]] <- factor(g[[key]])
    g$YEAR <- factor(g$year_num)
    calendar <- seq.int(min(g$year_num), max(g$year_num))
    missing_years <- setdiff(calendar, g$year_num)
    # glmmTMB can drop factor levels absent from the ENTIRE data subset, even
    # when explicitly supplied. A year missing in just one population is fine.
    # Stop before fitting if a whole-group gap would compress calendar distance.
    assert(!length(missing_years), paste("AR(1) calendar gap across the entire", group,
      "dataset:", paste(missing_years, collapse = ", "),
      "Review the time structure before fitting; no calendar years were compressed."))
    g$year_fac <- factor(g$year_num, levels = calendar)
    assert(nlevels(g$SPECIES) > 1L && nlevels(g$SITE_ID) > 1L && length(unique(g$year_num)) > 1L,
      paste("Insufficient species, sites or years for", group))
    g <- g[order(g$pop_id, g$year_num), , drop = FALSE]
    X <- stats::model.matrix(~ year_decade * onset_plasticity_z + year_decade * offset_plasticity_bio_z, g)
    assert(qr(X)$rank == ncol(X), paste("Rank-deficient fixed design for", group))
    groups[[group]] <- g
    years[[group]] <- g |>
      dplyr::group_by(pop_id, SPECIES, SITE_ID) |>
      dplyr::summarise(n_years = dplyr::n_distinct(year_num), first_year = min(year_num),
        last_year = max(year_num), .groups = "drop") |>
      dplyr::mutate(voltinism = group)
    summaries[[group]] <- data.frame(voltinism = group, n_observations = nrow(g), n_populations = nlevels(g$pop_id),
      n_species = nlevels(g$SPECIES), n_sites = nlevels(g$SITE_ID), min_year = min(g$year_num), max_year = max(g$year_num),
      n_populations_one_year = sum(years[[group]]$n_years == 1L), minimum_year_filter = "none_as_in_original")
  }
  list(groups = groups, plasticity = p, scaling = scaling, correlations = correlations,
    filter_summary = filter_summary, summary = dplyr::bind_rows(summaries), population_years = dplyr::bind_rows(years),
    unmatched_abundance = unmatched, unused_plasticity = unused_plasticity, year_center = year_center)
}

prepare_worker <- function(cfg) {
  set_phase(cfg, "Reading annual abundance and spatial plasticity")
  a <- utils::read.csv(cfg$abundance_file, stringsAsFactors = FALSE)
  p <- utils::read.csv(cfg$plasticity_file, stringsAsFactors = FALSE)
  result <- prepare_tables(a, p, cfg$abundance_quantile)
  set_phase(cfg, "Saving prepared inputs and join summaries")
  for (group in names(result$groups)) atomic_save(list(data = result$groups[[group]],
    year_center = result$year_center), file.path(cfg$out, paste0("prepared_", group, ".rds")))
  write_table(dplyr::bind_rows(result$groups), file.path(cfg$out, "df_population_trend_combined.csv"))
  for (name in c("plasticity", "scaling", "correlations", "filter_summary", "summary", "population_years", "unmatched_abundance", "unused_plasticity"))
    write_table(result[[name]], file.path(cfg$out, paste0("preparation_", name, ".csv")))
  write_table(data.frame(year_center = result$year_center), file.path(cfg$out, "year_center.csv"))
  if (cfg$make_figures) {
    tryCatch({
      plot <- ggplot2::ggplot(result$plasticity, ggplot2::aes(onset_plasticity_z, offset_plasticity_bio_z)) +
        ggplot2::geom_point(alpha = .15, size = .7) + ggplot2::geom_smooth(method = "lm", formula = y ~ x) +
        ggplot2::facet_wrap(~ voltinism) + ggplot2::theme_classic(base_size = 15, base_family = cfg$font_family) +
        ggplot2::labs(x = "Onset advancement (SD)", y = "Offset response (SD)")
      ggplot2::ggsave(file.path(cfg$out, "figures", "plot_plasticity_correlation.png"), plot, width = 7, height = 5, dpi = 300)
    }, error = function(e) writeLines(conditionMessage(e), file.path(cfg$out, "correlation_figure_error.txt")))
  }
  invisible(result$summary)
}

model_checks <- function(model) {
  coefs <- summary(model)$coefficients$cond
  ll <- stats::logLik(model)
  finite <- all(is.finite(coefs[, "Estimate"])) && all(is.finite(coefs[, "Std. Error"])) && all(coefs[, "Std. Error"] > 0)
  good <- isTRUE(model$fit$convergence == 0L) && isTRUE(model$sdr$pdHess) && finite && is.finite(as.numeric(ll))
  data.frame(convergence_code = model$fit$convergence, pdHess = isTRUE(model$sdr$pdHess), finite_fixed_SE = finite,
    logLik = as.numeric(ll), parameter_df = attr(ll, "df"), AIC = stats::AIC(model), BIC = stats::BIC(model),
    nobs = stats::nobs(model), numerical_checks_ok = good)
}
fixed_table <- function(model) {
  m <- summary(model)$coefficients$cond
  out <- data.frame(term = rownames(m), estimate = m[, 1], std.error = m[, 2], z = m[, 3], p_Wald = m[, 4], row.names = NULL)
  out$conf.low <- out$estimate - stats::qnorm(.975) * out$std.error
  out$conf.high <- out$estimate + stats::qnorm(.975) * out$std.error
  out$percent_change_per_decade_at_mean_plasticity <- ifelse(out$term == "year_decade", 100 * expm1(out$estimate), NA_real_)
  interaction <- grepl("year_decade:", out$term) | grepl(":year_decade", out$term)
  # exp(interaction) is a ratio of decadal abundance-change FACTORS per SD.
  # It is not itself a population's percent change per decade.
  out$ratio_of_decadal_factors_per_SD <- ifelse(interaction, exp(out$estimate), NA_real_)
  out
}
interaction_name <- function(beta, variable) {
  term <- intersect(c(paste0("year_decade:", variable), paste0(variable, ":year_decade")), names(beta))
  assert(length(term) == 1L, paste("Interaction missing:", variable)); term
}
effect_predictions <- function(beta, covariance, data, variable) {
  term <- interaction_name(beta, variable)
  pred <- expand.grid(year_num = seq.int(min(data$year_num), max(data$year_num)), plasticity_z = c(1, 0, -1))
  pred$delta_decades <- (pred$year_num - min(data$year_num)) / 10
  pred$slope <- unname(beta["year_decade"]) + unname(beta[term]) * pred$plasticity_z
  variance <- covariance["year_decade", "year_decade"] + pred$plasticity_z^2 * covariance[term, term] +
    2 * pred$plasticity_z * covariance["year_decade", term]
  assert(all(is.finite(variance)) && all(variance >= -1e-10), "Invalid prediction variance.")
  se <- pred$delta_decades * sqrt(pmax(variance, 0))
  log_rel <- pred$delta_decades * pred$slope
  pred$percent_change <- 100 * expm1(log_rel)
  pred$low_percent <- 100 * expm1(log_rel - stats::qnorm(.975) * se)
  pred$high_percent <- 100 * expm1(log_rel + stats::qnorm(.975) * se)
  pred$group <- factor(pred$plasticity_z, levels = c(1, 0, -1), labels = c("+1 SD", "Mean", "-1 SD"))
  pred$predictor <- variable
  pred
}
population_linear_trends <- function(model, data, group) {
  beta <- glmmTMB::fixef(model)$cond
  re <- glmmTMB::ranef(model, condVar = FALSE)$cond
  sp <- re$SPECIES; site <- re$SITE_ID
  assert("year_decade" %in% colnames(sp) && "year_decade" %in% colnames(site), "Species/site temporal slopes missing.")
  p <- unique(data[c("pop_id", "SPECIES", "SITE_ID", "onset_plasticity_z", "offset_plasticity_bio_z")])
  p$species_year_re <- sp[match(as.character(p$SPECIES), rownames(sp)), "year_decade"]
  p$site_year_re <- site[match(as.character(p$SITE_ID), rownames(site)), "year_decade"]
  assert(!anyNA(p), "Population/species/site slope join failed. Missing effects were not replaced by zero.")
  p$beta_log_per_decade <- unname(beta["year_decade"]) +
    unname(beta[interaction_name(beta, "onset_plasticity_z")]) * p$onset_plasticity_z +
    unname(beta[interaction_name(beta, "offset_plasticity_bio_z")]) * p$offset_plasticity_bio_z +
    p$species_year_re + p$site_year_re
  p$trend_percent_decade <- 100 * expm1(p$beta_log_per_decade)
  p$voltinism <- group
  p
}
export_full <- function(fit, data, group, out, cfg) {
  writeLines(capture.output(summary(fit)), file.path(out, "model_summary.txt"))
  write_table(fixed_table(fit), file.path(out, "fixed_effects.csv"))
  V <- as.matrix(stats::vcov(fit)$cond)
  write_table(data.frame(term = rownames(V), V, check.names = FALSE), file.path(out, "fixed_effect_covariance.csv"))
  write_table(population_linear_trends(fit, data, group), file.path(out, "model_based_linear_population_trends.csv"))
  status <- list()
  if (requireNamespace("performance", quietly = TRUE)) {
    status$collinearity <- tryCatch({
      col <- performance::check_collinearity(fit)
      saveRDS(col, file.path(out, "collinearity.rds"))
      write_table(as.data.frame(col), file.path(out, "collinearity.csv")); "OK"
    }, error = function(e) paste("ERROR:", conditionMessage(e)))
  } else status$collinearity <- "NOT_RUN: optional package performance is unavailable"
  beta <- glmmTMB::fixef(fit)$cond
  for (v in c("onset_plasticity_z", "offset_plasticity_bio_z")) {
    pred <- effect_predictions(beta, V, data, v)
    write_table(pred, file.path(out, paste0("predictions_", v, ".csv")))
    if (cfg$make_figures) status[[paste0("figure_", v)]] <- tryCatch({
      legend <- if (v == "onset_plasticity_z") "Onset advancement" else if (group == "univoltine") "Offset advancement" else "Offset response\n(towards delay)"
      plot <- ggplot2::ggplot(pred, ggplot2::aes(year_num, percent_change, colour = group, fill = group, group = group)) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = low_percent, ymax = high_percent), alpha = .10, colour = NA) +
        ggplot2::geom_line(linewidth = 1.4) +
        ggplot2::scale_colour_manual(values = c("+1 SD" = "#009E73", "Mean" = "#0072B2", "-1 SD" = "#D55E00")) +
        ggplot2::scale_fill_manual(values = c("+1 SD" = "#009E73", "Mean" = "#0072B2", "-1 SD" = "#D55E00")) +
        ggplot2::labs(x = "Year", y = "Relative change in abundance (%)", colour = legend, fill = legend) +
        ggplot2::theme_classic(base_size = 15, base_family = cfg$font_family)
      ggplot2::ggsave(file.path(cfg$out, "figures", paste0("plot_trend_", group, "_", v, ".png")), plot, width = 7, height = 5, dpi = 300)
      "OK"
    }, error = function(e) paste("ERROR:", conditionMessage(e)))
  }
  write_table(data.frame(stage = names(status), status = unlist(status), row.names = NULL), file.path(out, "export_status.csv"))
}

fit_worker <- function(cfg, group, variant) {
  # Register S3 methods even when this fresh worker only reads a checkpoint.
  assert(requireNamespace("glmmTMB", quietly = TRUE), "Package glmmTMB is required.")
  task <- paste(group, variant, sep = "__"); started <- Sys.time()
  out <- file.path(cfg$out, task); dir.create(out, showWarnings = FALSE)
  log_path <- file.path(out, "fit.rds")
  set_phase(cfg, "Reading prepared data", task, started)
  prepared <- readRDS(file.path(cfg$out, paste0("prepared_", group, ".rds")))
  data <- prepared$data
  signature <- digest::digest(list(cfg$signature, group, variant, data, formula_text(model_formula(variant))), algo = "sha256")
  cached <- file.exists(log_path)
  if (cached) {
    set_phase(cfg, "Reading saved model", task, started)
    saved <- readRDS(log_path)
    assert(identical(saved$signature, signature), "Model checkpoint does not match inputs or specification.")
    fit <- saved$fit
  } else {
    set_phase(cfg, "Fitting Gamma model", task, started)
    warnings <- character()
    fit <- withCallingHandlers(glmmTMB::glmmTMB(model_formula(variant), data = data,
      family = stats::Gamma(link = "log"), REML = FALSE, na.action = stats::na.fail,
      control = glmmTMB::glmmTMBControl(parallel = cfg$threads, optCtrl = list(trace = cfg$optimizer_trace))),
      warning = function(w) { warnings <<- unique(c(warnings, conditionMessage(w))); message("WARNING: ", conditionMessage(w)); invokeRestart("muffleWarning") })
    # Save before any export/plot. A failed diagnostic or font cannot lose the fit.
    saved <- list(fit = fit, signature = signature, warnings = warnings,
      elapsed_minutes = as.numeric(difftime(Sys.time(), started, units = "mins")))
    atomic_save(saved, log_path)
  }
  set_phase(cfg, "Checking saved fit", task, started)
  checks <- model_checks(fit)
  write_table(checks, file.path(out, "fit_checks.csv"))
  writeLines(saved$warnings, file.path(out, "warnings.txt"))
  metrics <- list(checks = checks, signature = signature, data_signature = digest::digest(data, algo = "sha256"),
    formula = formula_text(model_formula(variant)), elapsed_minutes = saved$elapsed_minutes)
  atomic_save(metrics, file.path(out, "metrics.rds"))
  status <- if (checks$numerical_checks_ok) "OK" else "NUMERICAL_REVIEW_REQUIRED"
  if (checks$numerical_checks_ok && variant == "full") {
    set_phase(cfg, "Exporting model results and figures", task, started)
    export_full(fit, data, group, out, cfg)
  }
  data.frame(task = task, group = group, variant = variant, status = status, cached = cached,
    elapsed_minutes = as.numeric(difftime(Sys.time(), started, units = "mins")), note = "")
}

lrt_from_metrics <- function(full, reduced, group, variable) {
  assert(isTRUE(full$checks$numerical_checks_ok) && isTRUE(reduced$checks$numerical_checks_ok), "LRT requires two numerically acceptable fits.")
  assert(identical(full$data_signature, reduced$data_signature), "LRT fits use different data.")
  statistic <- 2 * (full$checks$logLik - reduced$checks$logLik)
  df <- full$checks$parameter_df - reduced$checks$parameter_df
  assert(is.finite(statistic) && statistic >= -1e-6 && df == 1, "Invalid nested likelihood comparison; review optimization.")
  data.frame(voltinism = group, interaction = paste0("year_decade:", variable), Chisq = max(0, statistic), LRT_df = df,
    p_LRT = stats::pchisq(max(0, statistic), df = df, lower.tail = FALSE),
    AIC_full = full$checks$AIC, AIC_reduced = reduced$checks$AIC,
    delta_AIC = reduced$checks$AIC - full$checks$AIC,
    BIC_full = full$checks$BIC, BIC_reduced = reduced$checks$BIC,
    delta_BIC = reduced$checks$BIC - full$checks$BIC)
}
combine_exports <- function(cfg, statuses) {
  terms <- list(); trends <- list(); lrts <- list(); issues <- character()
  for (group in c("univoltine", "multivoltine")) {
    full_task <- paste(group, "full", sep = "__")
    if (!any(statuses$task == full_task & statuses$status == "OK")) next
    folder <- file.path(cfg$out, full_task)
    t <- utils::read.csv(file.path(folder, "fixed_effects.csv")); t$model <- group; terms[[group]] <- t
    trends[[group]] <- utils::read.csv(file.path(folder, "model_based_linear_population_trends.csv"))
    full <- readRDS(file.path(folder, "metrics.rds"))
    if (cfg$run_lrt) for (effect in c("onset", "offset")) {
      task <- paste(group, paste0("without_", effect, "_interaction"), sep = "__")
      if (!any(statuses$task == task & statuses$status == "OK")) next
      reduced <- readRDS(file.path(cfg$out, task, "metrics.rds"))
      variable <- if (effect == "onset") "onset_plasticity_z" else "offset_plasticity_bio_z"
      key <- paste(group, effect)
      lrts[[key]] <- tryCatch(lrt_from_metrics(full, reduced, group, variable), error = function(e) {
        issues <<- c(issues, paste(key, conditionMessage(e))); NULL
      })
    }
  }
  writeLines(issues, file.path(cfg$out, "comparison_errors.txt"))
  if (length(terms)) write_table(dplyr::bind_rows(terms), file.path(cfg$out, "summary_population_trend_model_terms.csv"))
  if (length(lrts)) write_table(dplyr::bind_rows(lrts), file.path(cfg$out, "summary_population_trend_LRTs.csv"))
  if (length(trends)) {
    pop <- dplyr::bind_rows(trends)
    write_table(pop, file.path(cfg$out, "model_based_linear_population_trends.csv"))
    summaries <- pop |>
      dplyr::group_by(voltinism) |>
      dplyr::summarise(n_populations = dplyr::n(), mean_trend = mean(trend_percent_decade),
        median_trend = stats::median(trend_percent_decade), q25 = stats::quantile(trend_percent_decade, .25),
        q75 = stats::quantile(trend_percent_decade, .75), percent_positive = 100 * mean(trend_percent_decade > 0), .groups = "drop")
    write_table(summaries, file.path(cfg$out, "summary_linear_population_trends.csv"))
    if (cfg$make_figures && length(trends) == 2L) tryCatch({
      plot <- ggplot2::ggplot(pop, ggplot2::aes(trend_percent_decade, fill = voltinism, colour = voltinism)) +
        ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), binwidth = 5, position = "identity", alpha = .35, colour = "white") +
        ggplot2::geom_density(linewidth = 1.1, alpha = 0) +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey30") +
        ggplot2::geom_vline(data = summaries, ggplot2::aes(xintercept = median_trend, colour = voltinism), linetype = "dashed") +
        ggplot2::coord_cartesian(xlim = stats::quantile(pop$trend_percent_decade, c(.01, .99))) +
        ggplot2::scale_colour_manual(values = c(univoltine = "#2F5D62", multivoltine = "#B85C38")) +
        ggplot2::scale_fill_manual(values = c(univoltine = "#2F5D62", multivoltine = "#B85C38")) +
        ggplot2::theme_classic(base_size = 15, base_family = cfg$font_family) +
        ggplot2::labs(x = "Linear abundance change (% per decade)", y = "Density", fill = NULL, colour = NULL) +
        ggplot2::theme(legend.position = "top")
      ggplot2::ggsave(file.path(cfg$out, "figures", "plot_density_pop_trends_univol_multivol.png"), plot, width = 7, height = 5, dpi = 300)
    }, error = function(e) writeLines(conditionMessage(e), file.path(cfg$out, "distribution_figure_error.txt")))
  }
  if (requireNamespace("writexl", quietly = TRUE)) tryCatch({
    sheets <- list(workflow_status = statuses, fixed_effects = dplyr::bind_rows(terms),
      likelihood_tests = dplyr::bind_rows(lrts),
      data_summary = utils::read.csv(file.path(cfg$out, "preparation_summary.csv")),
      plasticity_scaling = utils::read.csv(file.path(cfg$out, "preparation_scaling.csv")),
      correlations = utils::read.csv(file.path(cfg$out, "preparation_correlations.csv")))
    sheets <- sheets[vapply(sheets, ncol, integer(1)) > 0L]
    writexl::write_xlsx(sheets, file.path(cfg$out, "population_trend_results_spatial_plasticity_15km.xlsx"))
  }, error = function(e) writeLines(conditionMessage(e), file.path(cfg$out, "xlsx_error.txt")))
  length(issues)
}

launch_child <- function(cfg, fun, args, log) {
  child <- callr::r_bg(function(engine, fun, args) {
    e <- new.env(parent = globalenv()); sys.source(engine, envir = e)
    do.call(e[[fun]], args)
  }, args = list(engine = cfg$engine, fun = fun, args = args), stdout = log, stderr = "2>&1", supervise = TRUE, wd = cfg$project_root)
  while (child$is_alive()) Sys.sleep(1)
  child$get_result()
}
workflow <- function(cfg) {
  assert(identical(unname(tools::md5sum(c(cfg$abundance_file, cfg$plasticity_file))), cfg$input_md5),
    "Inputs changed after the run was configured.")
  set_phase(cfg, "Preparing input data")
  launch_child(cfg, "prepare_worker", list(cfg = cfg), file.path(cfg$out, "preparation.log"))
  if (cfg$prepare_only) { set_phase(cfg, "Preparation complete; no models fitted"); return(cfg$out) }
  manifest <- expand.grid(group = c("univoltine", "multivoltine"),
    variant = if (cfg$run_lrt) c("full", "without_onset_interaction", "without_offset_interaction") else "full", stringsAsFactors = FALSE)
  manifest$task <- paste(manifest$group, manifest$variant, sep = "__")
  write_table(manifest, file.path(cfg$out, "task_manifest.csv"))
  statuses <- list()
  for (i in seq_len(nrow(manifest))) {
    task <- manifest[i, ]; group <- task$group; variant <- task$variant
    if (variant != "full" && !identical(statuses[[paste(group, "full", sep = "__")]]$status, "OK")) {
      result <- data.frame(task = task$task, group = group, variant = variant, status = "SKIPPED_FULL_MODEL_FAILED",
        cached = FALSE, elapsed_minutes = 0, note = "The full model requires review.")
    } else {
      dir.create(file.path(cfg$out, task$task), showWarnings = FALSE)
      set_phase(cfg, paste("Starting task", i, "of", nrow(manifest)), task$task)
      result <- tryCatch(launch_child(cfg, "fit_worker", list(cfg = cfg, group = group, variant = variant),
        file.path(cfg$out, task$task, "worker.log")), error = function(e) data.frame(task = task$task,
          group = group, variant = variant, status = "ERROR", cached = FALSE, elapsed_minutes = NA_real_, note = conditionMessage(e)))
    }
    statuses[[task$task]] <- result
    write_table(dplyr::bind_rows(statuses), file.path(cfg$out, "workflow_status.csv"))
  }
  set_phase(cfg, "Combining exports")
  status <- dplyr::bind_rows(statuses)
  issues <- combine_exports(cfg, status)
  done <- all(status$status == "OK") && issues == 0L
  set_phase(cfg, if (done) "Finished; review model and export diagnostics" else "Finished with tasks requiring review; inspect workflow_status.csv")
  list(output_root = cfg$out, status = status, comparisons_ok = issues == 0L)
}

watch <- function(job = getOption("pheno.trends.active_job")) {
  assert(!is.null(job$process), "No population-trend job is registered.")
  tryCatch({
    repeat {
      progress <- if (file.exists(job$cfg$progress_file)) tryCatch(suppressWarnings(readRDS(job$cfg$progress_file)), error = function(e) NULL) else NULL
      alive <- job$process$is_alive()
      sfile <- file.path(job$cfg$out, "workflow_status.csv")
      status <- if (file.exists(sfile)) tryCatch(utils::read.csv(sfile), error = function(e) NULL) else NULL
      elapsed <- as.numeric(difftime(Sys.time(), job$started, units = "mins"))
      message("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ",
        if (is.null(progress)) "Starting" else paste(progress$task, progress$phase),
        "\nTasks resolved: ", if (is.null(status)) 0L else nrow(status), "/", if (job$cfg$prepare_only) 0L else if (job$cfg$run_lrt) 6L else 2L,
        " | cached: ", if (is.null(status)) 0L else sum(status$cached),
        " | failed/skipped: ", if (is.null(status)) 0L else sum(status$status != "OK"),
        " | elapsed: ", round(elapsed, 1), " min | worker running: ", alive)
      if (!is.null(progress) && nzchar(progress$task)) message("Current task elapsed: ",
        round(as.numeric(difftime(Sys.time(), progress$task_started, units = "mins")), 1), " min. No reliable remaining-time estimate yet.")
      if (!alive) break
      Sys.sleep(job$cfg$monitor_seconds)
    }
    result <- job$process$get_result()
    message("Outputs: ", job$cfg$out)
    invisible(result)
  }, interrupt = function(e) {
    message("Display paused; THE QUEUE CONTINUES. Keep R open. Resume with pheno_population_trends$watch().")
    invisible(job)
  })
}
wait_active <- function(job) {
  if (is.null(job$process)) return(invisible(NULL))
  while (isTRUE(tryCatch(job$process$is_alive(), error = function(e) FALSE))) {
    message("Waiting for the existing worker to finish before starting population-trend fits.")
    Sys.sleep(30)
  }
}
run <- function(cfg) {
  pkgs <- c("callr", "dplyr", "glmmTMB", "digest", "here")
  if (cfg$make_figures) pkgs <- c(pkgs, "ggplot2")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  assert(!length(missing), paste("Missing packages:", paste(missing, collapse = ", ")))
  if (is.null(cfg$project_root)) cfg$project_root <- here::here()
  cfg$project_root <- normalizePath(cfg$project_root, winslash = "/", mustWork = TRUE)
  if (is.null(cfg$abundance_file)) cfg$abundance_file <- file.path(cfg$project_root, "output", "pheno_abund_estimates_allspp.csv")
  if (is.null(cfg$plasticity_file)) cfg$plasticity_file <- latest_extraction(cfg$project_root)
  if (is.null(cfg$output_root)) cfg$output_root <- file.path(cfg$project_root, "output", "population_trends", "spatial_plasticity_15km")
  cfg$abundance_file <- normalizePath(cfg$abundance_file, winslash = "/", mustWork = TRUE)
  cfg$plasticity_file <- normalizePath(cfg$plasticity_file, winslash = "/", mustWork = TRUE)
  assert(isTRUE(cfg$abundance_quantile == .999), "This adaptation preserves the original 99.9% abundance threshold.")
  assert(is.finite(cfg$monitor_seconds) && cfg$monitor_seconds >= 1 && cfg$monitor_seconds <= 60, "Monitor interval must be 1-60 seconds.")
  assert(is.finite(cfg$threads) && cfg$threads >= 1, "threads must be positive.")
  for (option in c("pheno.spatial.active_job", "pheno.offset.recovery.active_job", "pheno.review.active_job", "pheno.extract.active_job", "pheno.trends.active_job")) wait_active(getOption(option))
  if (exists("spatial_job", envir = .GlobalEnv, inherits = FALSE)) wait_active(get("spatial_job", envir = .GlobalEnv))
  message("Plasticity input: ", cfg$plasticity_file)
  versions <- stats::setNames(vapply(c("glmmTMB", "TMB", "Matrix", "lme4", "dplyr"),
    function(p) as.character(utils::packageVersion(p)), character(1)), c("glmmTMB", "TMB", "Matrix", "lme4", "dplyr"))
  cfg$input_md5 <- unname(tools::md5sum(c(cfg$abundance_file, cfg$plasticity_file)))
  cfg$signature <- digest::digest(list(version = "spatial_plasticity_trends_1.0", inputs = cfg$input_md5,
    packages = versions, threads = cfg$threads, quantile = cfg$abundance_quantile, formula = formula_text(model_formula())), algo = "sha256")
  cfg$out <- file.path(cfg$output_root, paste0("run_", substr(cfg$signature, 1, 12)))
  dir.create(file.path(cfg$out, "figures"), recursive = TRUE, showWarnings = FALSE)
  cfg$out <- normalizePath(cfg$out, winslash = "/", mustWork = TRUE)
  cfg$progress_file <- file.path(cfg$out, "progress.rds")
  cfg$engine <- file.path(cfg$out, "workflow_engine.R")
  functions <- ls(envir = environment(run), all.names = TRUE)
  functions <- functions[vapply(functions, function(n) is.function(get(n, envir = environment(run))), logical(1))]
  dump(functions, file = cfg$engine, envir = environment(run), control = "all")
  atomic_save(cfg, file.path(cfg$out, "run_config.rds"))
  write_table(data.frame(package = names(versions), version = unname(versions)), file.path(cfg$out, "package_versions.csv"))
  write_table(data.frame(input = c("annual_abundance", "spatial_plasticity"),
    path = c(cfg$abundance_file, cfg$plasticity_file), md5 = cfg$input_md5), file.path(cfg$out, "input_sources.csv"))
  writeLines(c(
    "Population trends using contextual plasticity extracted from the 15-km spatial phenology fits.",
    "Abundance model: Gamma/log, ML, year interactions with onset and offset, independent species/site intercept and temporal slope, population AR(1).",
    "This abundance model has NO explicit spatial field. Spatial correction of predictor estimates does not establish spatial independence of abundance residuals.",
    "Original filter retained: positive abundance <= global 99.9th percentile, computed before deduplication. No minimum-years filter added.",
    "Conflicting abundance values per species/site/year stop the preparation; identical repeated records are removed.",
    "Plasticity z-scores use all extracted populations within each voltinism group, before the abundance join, as in the original script.",
    "Year is centred over cleaned/deduplicated abundance rows before matching plasticity, then expressed per decade.",
    "AR(1) retains calendar spacing for missing observations within populations. Preparation stops if an internal year is absent from the entire voltinism group, since glmmTMB may drop that factor level.",
    "Positive onset = stronger advancement. Univoltine offset is sign-inverted: positive = advancement. Multivoltine offset retains raw sign: positive = delay.",
    "For multivoltines, a higher offset value is a response towards delay, not necessarily greater absolute plasticity; many estimated slopes are negative.",
    "Default six fits: two full models and four models each dropping one year-by-plasticity interaction. Set run_lrt=FALSE for two fits and Wald tests only.",
    "An exponentiated interaction is a ratio of decadal abundance-change factors per one SD of plasticity, not percent population change itself.",
    "Effect curves hold the other plasticity predictor at zero and species/site/AR(1) deviations at zero, relative to the first observed year.",
    "Model-based population trends contain the linear fixed + species + site slope; AR(1) annual deviations are NOT part of that linear trend.",
    "Intervals/tests are conditional on estimated plasticity. First-stage plasticity uncertainty and GAM abundance-index uncertainty are not propagated.",
    "Numerical convergence and collinearity are checked. A residual spatial/temporal adequacy analysis is still needed; DHARMa simulations are not run automatically.",
    "Previously extracted univoltine plasticity retains the known spatial-range warning; this workflow does not resolve it.",
    "Model fits are saved before exports. A failed fit is retained for review, never silently refitted or used for inference.",
    "Rerun the same script to reuse matching fit checkpoints. Every run creates a timestamped controller log; each task has worker.log.",
    "Review workflow_status.csv, each fit_checks.csv, warnings.txt, and export_status.csv. Missing optional collinearity package is recorded.",
    "https://glmmtmb.github.io/glmmTMB/articles/covstruct.html",
    "https://glmmtmb.github.io/glmmTMB/articles/troubleshooting.html"
  ), file.path(cfg$out, "README.txt"))
  # Status files describe this invocation; fit checkpoints remain untouched.
  # Remove only reproducible aggregate exports from previous invocations so
  # a prepare-only/partial run cannot be mistaken for current complete results.
  for (f in c("workflow_status.csv", "progress.rds", "task_manifest.csv",
              "summary_population_trend_model_terms.csv", "summary_population_trend_LRTs.csv",
              "model_based_linear_population_trends.csv", "summary_linear_population_trends.csv",
              "population_trend_results_spatial_plasticity_15km.xlsx", "comparison_errors.txt"))
    if (file.exists(file.path(cfg$out, f))) unlink(file.path(cfg$out, f))
  log <- tempfile(paste0("controller_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_"), tmpdir = cfg$out, fileext = ".log")
  process <- callr::r_bg(function(engine, cfg) {
    e <- new.env(parent = globalenv()); sys.source(engine, envir = e); e$workflow(cfg)
  }, args = list(engine = cfg$engine, cfg = cfg), stdout = log, stderr = "2>&1", supervise = TRUE, wd = cfg$project_root)
  job <- list(process = process, cfg = cfg, started = Sys.time(), log_file = log)
  options(pheno.trends.active_job = job)
  message("Started population-trend workflow. Keep R open. Log: ", log)
  invisible(job)
}
stop_job <- function(job = getOption("pheno.trends.active_job")) {
  assert(!is.null(job$process), "No registered job.")
  job$process$kill_tree()
  message("Workflow stopped. Saved fits remain; an unsaved fit will restart next time.")
  invisible(job)
}
list(run = run, watch = watch, stop = stop_job)
})

if (!isTRUE(getOption("pheno.trends.functions_only", FALSE))) {
  population_trend_job <- pheno_population_trends$run(population_trend_config)
  pheno_population_trends$watch(population_trend_job)
}
