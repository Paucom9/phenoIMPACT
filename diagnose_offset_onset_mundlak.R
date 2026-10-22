# phenoIMPACT: computationally simplified 2 x 2 offset diagnostic.
# Same observations and random effects in all four models; NO spatial field.
# Reads prepared data only. Does not load, overwrite or refit accepted sdmTMB models.
# source('diagnose_offset_onset_mundlak.R')
# offset_diagnostic <- diagnose_offset_onset_mundlak()
# Optional smaller diagnostic: max_populations_per_species = 10

diagnose_offset_onset_mundlak <- function(
  project_root = "E:/phenoIMPACT project/code/phenoIMPACT",
  environment_run = file.path(project_root, "output", "phenology_plasticity",
    "environment_mundlak_by_voltinism_qc_final_15km", "run_de0890c85d4a"),
  groups = c("multivoltine", "univoltine"),
  max_populations_per_species = Inf, seed = 20261022L) {
  need <- function(x, msg) if (!isTRUE(x)) stop(msg, call. = FALSE)
  need(requireNamespace("lme4", quietly = TRUE), "Package lme4 is required.")
  need(all(groups %in% c("multivoltine", "univoltine")) && !anyDuplicated(groups), "Invalid groups.")
  need(length(max_populations_per_species) == 1 && max_populations_per_species >= 2,
       "Use Inf or at least two populations per species.")
  out <- file.path(environment_run, paste0("diagnostic_offset_onset_mundlak_",
                                          format(Sys.time(), "%Y%m%d_%H%M%S")))
  need(!dir.exists(out), "Output already exists; retry in a second.")
  dir.create(out, recursive = TRUE)
  csv <- function(x, name) write.csv(x, file.path(out, name), row.names = FALSE)
  mods <- c("photoperiod", "background", "predictability", "trend")
  anomaly <- "clim_anomaly_tw90"
  spec <- data.frame(model = c("A", "B", "C", "D"),
    mundlak = c(FALSE, FALSE, TRUE, TRUE), onset = c(FALSE, TRUE, FALSE, TRUE))
  summaries <- checks <- interactions <- population_slopes <- manifests <- list()
  set.seed(seed)

  for (g in groups) {
    path <- file.path(environment_run, paste0("offset_", g), "prepared_input.rds")
    need(file.exists(path), paste("Missing prepared input:", path))
    p <- readRDS(path); d <- p$data; rm(p); gc()
    wc <- paste0(mods, "_within_species")
    bc <- paste0(mods, "_species_mean")
    required <- c("OFFSET_mean", "ONSET_mean_z", anomaly, "SPECIES", "SITE_ID", "YEAR", wc, bc)
    need(all(required %in% names(d)), "Prepared data lack required variables.")
    d <- d[, required, drop = FALSE]
    need(all(complete.cases(d)), "Prepared source contains missing values; review before fitting.")
    need(all(vapply(d[setdiff(required, c("SPECIES", "SITE_ID"))],
                    function(x) all(is.finite(x)), logical(1))), "Non-finite values.")
    d$SPECIES <- as.character(d$SPECIES); d$SITE_ID <- as.character(d$SITE_ID)
    d$pop_id <- paste(d$SPECIES, d$SITE_ID, sep = "\r")
    pop <- unique(d[c("pop_id", "SPECIES", "SITE_ID", wc, bc)])
    need(!anyDuplicated(pop$pop_id), "Moderator components vary within population.")
    # Optional balanced cap, selecting complete populations with all available years.
    if (is.finite(max_populations_per_species)) {
      selected <- unlist(lapply(split(pop$pop_id, pop$SPECIES), function(ids) {
        ids <- sort(ids)
        if (length(ids) <= max_populations_per_species) ids else
          sample(ids, size = as.integer(max_populations_per_species))
      }), use.names = FALSE)
      d <- d[d$pop_id %in% selected, , drop = FALSE]
      pop <- pop[pop$pop_id %in% selected, , drop = FALSE]
    }
    # Recover original standardized moderators, then recompute equal-population
    # species means on this SAME diagnostic dataset (essential for subset runs).
    raw <- paste0(mods, "_original")
    for (j in seq_along(mods)) {
      pop[[raw[j]]] <- pop[[wc[j]]] + pop[[bc[j]]]
      sp <- tapply(pop[[raw[j]]], pop$SPECIES, mean)
      pop[[bc[j]]] <- unname(sp[pop$SPECIES])
      pop[[wc[j]]] <- pop[[raw[j]]] - pop[[bc[j]]]
    }
    ix <- match(d$pop_id, pop$pop_id)
    for (v in c(raw, wc, bc)) d[[v]] <- pop[[v]][ix]
    d$SPECIES <- factor(d$SPECIES); d$SPECIES_slope <- d$SPECIES
    d$SITE_ID <- factor(d$SITE_ID)
    d$site_year_id <- interaction(d$SITE_ID, d$YEAR, drop = TRUE)
    manifests[[g]] <- data.frame(group = g, prepared = path,
      prepared_md5 = unname(tools::md5sum(path)), rows = nrow(d),
      populations = nrow(pop), species = nlevels(d$SPECIES),
      max_populations_per_species = max_populations_per_species, seed = seed)
    csv(do.call(rbind, manifests), "input_manifest.csv")
    csv(pop[c("SPECIES", "SITE_ID")], paste0(g, "_selected_populations.csv"))
    all_prediction_rows <- list()

    for (i in seq_len(nrow(spec))) {
      s <- spec[i, ]; label <- paste(g, s$model, sep = "_")
      predictors <- if (s$mundlak) c(wc, bc) else raw
      rhs <- paste0(anomaly, " * (", paste(predictors, collapse = " + "), ")",
        if (s$onset) " + ONSET_mean_z" else "",
        " + (1 | SITE_ID) + (1 | SPECIES) + (0 + ", anomaly,
        " | SPECIES_slope) + (1 | site_year_id)")
      f <- as.formula(paste("OFFSET_mean ~", rhs))
      message("Fitting ", label, " | ", nrow(d), " rows | no spatial field")
      warnings <- character(); t0 <- Sys.time()
      fit <- withCallingHandlers(lme4::lmer(f, data = d, REML = FALSE,
        control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))),
        warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
      need(nobs(fit) == nrow(d), "Unexpected change in fitted observations.")
      b <- lme4::fixef(fit); V <- as.matrix(vcov(fit))
      # Fixed slopes at every observed population context: numerical derivative
      # of the LINEAR design, exactly unit difference in anomaly (not a refit).
      nd <- pop; nd[[anomaly]] <- 0; nd$ONSET_mean_z <- 0
      ff <- lme4::nobars(f); tt <- delete.response(terms(ff))
      X0 <- model.matrix(tt, nd); nd[[anomaly]] <- 1
      X1 <- model.matrix(tt, nd)
      L <- (X1-X0)[, names(b), drop = FALSE]
      fixed <- as.vector(L %*% b)
      meanL <- colMeans(L)
      average <- sum(meanL*b); average_se <- sqrt(drop(t(meanL)%*%V%*%meanL))
      re <- lme4::ranef(fit)$SPECIES_slope
      need(anomaly %in% names(re), "Cannot find species random anomaly slopes.")
      rs <- re[match(pop$SPECIES, rownames(re)), anomaly]
      need(all(is.finite(rs)), "Cannot align species random slopes.")
      conditional <- fixed + rs
      summaries[[label]] <- data.frame(group = g, model = s$model,
        mundlak = s$mundlak, controls_onset = s$onset,
        fixed_mean_slope = average, SE = average_se,
        CI_low = average - 1.96*average_se, CI_high = average + 1.96*average_se,
        conditional_mean_slope = mean(conditional),
        fraction_positive_fixed = mean(fixed > 0),
        fraction_positive_conditional = mean(conditional > 0),
        nobs = nobs(fit), AIC = AIC(fit), elapsed_minutes =
          as.numeric(difftime(Sys.time(), t0, units = "mins")))
      population_slopes[[label]] <- data.frame(group = g, model = s$model,
        SPECIES = pop$SPECIES, SITE_ID = pop$SITE_ID,
        fixed_slope = fixed, conditional_slope = conditional)
      messages <- fit@optinfo$conv$lme4$messages
      checks[[label]] <- data.frame(group = g, model = s$model,
        optimizer_code = paste(fit@optinfo$conv$opt, collapse = ";"),
        singular = lme4::isSingular(fit),
        convergence_messages = paste(messages, collapse = " | "),
        warnings = paste(unique(warnings), collapse = " | "),
        residual_SD = sigma(fit))
      ints <- grep(":", names(b), value = TRUE)
      interactions[[label]] <- data.frame(group = g, model = s$model, term = ints,
        estimate = b[ints], SE = sqrt(diag(V))[ints])
      csv(do.call(rbind, summaries), "model_comparison.csv")
      csv(do.call(rbind, checks), "fit_checks.csv")
      csv(do.call(rbind, interactions), "interaction_coefficients.csv")
      csv(do.call(rbind, population_slopes), "population_slopes.csv")
      # Compact saved coefficients only; no large lmer model objects retained.
      saveRDS(list(formula = f, beta = b, covariance = V,
                   checks = checks[[label]]), file.path(out, paste0(label, "_coefficients.rds")))
      rm(fit, V, X0, X1, L); gc()
    }
    rm(d, pop); gc()
  }
  ans <- do.call(rbind, summaries)
  writeLines(c(
    "DIAGNOSTIC ONLY: Gaussian ML models without spatial fields; not substitutes for accepted models.",
    "A=no Mundlak/no onset; B=no Mundlak/with onset; C=Mundlak/no onset; D=Mundlak/with onset.",
    "All four models use identical observations and random structure within each voltinism group.",
    "Rows are from the current prepared FINAL-QC inputs. This does not isolate earlier QC or spatial changes.",
    "Species means use equal population weights, recomputed for any capped diagnostic subset.",
    "Slope units: days per inherited SD of the 90-day temperature anomaly, not days per degree C.",
    "Mean fixed slope averages observed population contexts with equal population weights.",
    "Its CI uses full fixed covariance but omits spatial dependence and is diagnostic only.",
    "Conditional slopes add species random slopes; fractions positive are descriptive point estimates, not significance tests.",
    "Positive slope=offset delay; negative=advancement. Onset-controlled slopes are partial associations.",
    "Check convergence messages and singularity before interpreting contrasts.",
    "AIC comparisons are confined to these simplified models within the same group.",
    "No original model/input/output is changed. No demographic process is stopped.",
    "Historical saved-model comparison remains a separate check; these models isolate onset/Mundlak under a common simplified structure."
  ), file.path(out, "README.txt"))
  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
  review_files <- list.files(out, pattern = "\\.(csv|txt)$", full.names = TRUE)
  if (requireNamespace("zip", quietly = TRUE))
    zip::zipr(file.path(out, "offset_diagnostic_review.zip"), review_files, root = out)
  print(ans, row.names = FALSE)
  message("Diagnostic output: ", normalizePath(out, winslash = "/"))
  invisible(list(summary = ans, checks = do.call(rbind, checks), output_dir = out))
}
