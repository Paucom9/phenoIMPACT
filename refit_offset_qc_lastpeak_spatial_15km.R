# =============================================================================
# phenoIMPACT | QC-filtered spatial OFFSET refits (15-km mesh)
# Version 1.0 | 2026-10-05
#
# Scientific change only:
#   exclude offset rows where LAST_PEAK is finite and OFFSET_mean < LAST_PEAK.
#
# OFFSET_mean is NOT recalculated. The source models are the current FIT_OK
# onset-adjusted offset models from run_3911d7bc87ac. The 90-day climate window,
# all environmental moderators, annual ONSET_mean covariate, random effects,
# Gaussian likelihood and original 15-km triangulation are retained.
#
# All parameters are re-estimated on the filtered rows. When the random-effect
# grouping universe is unchanged, the current FIT_OK model is used as a warm
# start. Otherwise the filtered model is fitted without a full latent warm start.
#
# This script does NOT refit onset/first peak, recalculate phenology, extract final
# population plasticity, or refit abundance models. Review the two QC fits first.
#
# Usage:
#   source(file.choose(), encoding = "UTF-8")
#   offset_qc_job <- pheno_offset_qc$run(monitor = FALSE)
#   pheno_offset_qc$status()
# =============================================================================

pheno_offset_qc <- local({

  VERSION <- "offset_onset_qc_lastpeak_spatial15_v1.0"
  ROOT <- "E:/phenoIMPACT project/code/phenoIMPACT"
  SOURCE_RUN <- file.path(
    ROOT, "output", "phenology_plasticity", "offset_with_onset_15km",
    "run_3911d7bc87ac"
  )
  ANALYSES <- c("offset_univoltine", "offset_multivoltine")
  CORE <- c("sdmTMB", "TMB", "Matrix", "fmesher", "lme4")
  START_KEYS <- c("b_j", "ln_kappa", "ln_tau_E", "ln_phi",
                  "re_cov_pars", "re_b_pars", "epsilon_st")
  state <- new.env(parent = emptyenv())

  need <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)
  ftxt <- function(x) paste(deparse(x, width.cutoff = 500L), collapse = " ")
  csv <- function(x, path) utils::write.csv(x, path, row.names = FALSE, na = "")
  md5 <- function(x) unname(tools::md5sum(x))
  stamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  same_num <- function(a, b, tol = 1e-8) {
    length(a) == length(b) && all(is.finite(a)) && all(is.finite(b)) &&
      (length(a) == 0L || max(abs(as.numeric(a) - as.numeric(b))) <= tol)
  }

  atomic_rds <- function(x, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    tmp <- tempfile(".writing_", tmpdir = dirname(path), fileext = ".rds")
    on.exit(unlink(tmp), add = TRUE)
    saveRDS(x, tmp, compress = "gzip")
    bak <- paste0(path, ".previous")
    had <- file.exists(path)
    if (had) {
      need(file.copy(path, bak, overwrite = TRUE), paste("Cannot back up", path))
      need(file.remove(path), paste("Cannot replace", path))
    }
    if (!file.rename(tmp, path)) {
      if (had) file.copy(bak, path, overwrite = TRUE)
      stop("Cannot finalize checkpoint: ", path, call. = FALSE)
    }
    if (had) unlink(bak)
    invisible(path)
  }

  phase <- function(cfg, stage, analysis = "", detail = "") {
    rec <- list(time = Sys.time(), stage = stage, analysis = analysis, detail = detail)
    atomic_rds(rec, file.path(cfg$run, "progress.rds"))
    message("[", stamp(), "] ", analysis, " | ", stage, " | ", detail)
  }

  check_packages <- function() {
    packages <- c(CORE, "data.table", "digest", "callr", "ps")
    missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
    need(!length(missing), paste("Missing packages:", paste(missing, collapse = ", ")))
    need(as.character(utils::packageVersion("sdmTMB")) == "1.1.0",
         "Use the original sdmTMB 1.1.0 library.")
    invisible(TRUE)
  }

  versions <- function() data.frame(
    package = CORE,
    version = vapply(CORE, function(x) as.character(utils::packageVersion(x)), character(1))
  )

  disk_check <- function(root, min_gib = 50) {
    drive <- if (.Platform$OS.type == "windows") substr(root, 1L, 3L) else "/"
    d <- ps::ps_disk_usage(drive)
    free <- as.numeric(d$available) / 1024^3
    need(is.finite(free) && free >= min_gib,
         paste0("Only ", round(free, 1), " GiB free on ", drive,
                "; require at least ", min_gib, " GiB."))
    message("Disk preflight: ", round(free, 1), " GiB available on ", drive)
    invisible(free)
  }

  source_paths <- function(analysis) {
    d <- file.path(SOURCE_RUN, analysis)
    list(
      dir = d,
      input = file.path(d, "input.rds"),
      fit = file.path(d, "fit_with_onset.rds"),
      complete = file.path(d, "completed.rds")
    )
  }

  fixed_table <- function(fit, formula, data) {
    Xnames <- colnames(stats::model.matrix(
      lme4::nobars(formula), data[seq_len(min(3L, nrow(data))), , drop = FALSE]
    ))
    ix <- which(names(fit$sd_report$par.fixed) == "b_j")
    need(length(ix) == length(Xnames), "Cannot align fixed coefficients.")
    need(identical(names(fit$model$par), names(fit$sd_report$par.fixed)) &&
           same_num(fit$model$par, fit$sd_report$par.fixed),
         "Point estimates and covariance refer to different parameters.")
    b <- as.numeric(fit$sd_report$par.fixed[ix])
    V <- as.matrix(fit$sd_report$cov.fixed[ix, ix, drop = FALSE])
    need(all(is.finite(b)) && all(is.finite(V)) && all(diag(V) > 0),
         "Invalid fixed-effect estimate/covariance.")
    dimnames(V) <- list(Xnames, Xnames)
    se <- sqrt(diag(V)); z <- stats::qnorm(.975)
    data.frame(
      term = Xnames, estimate = b, std.error = se,
      conf.low = b - z * se, conf.high = b + z * se,
      p_Wald = 2 * stats::pnorm(-abs(b / se)), row.names = NULL
    )
  }

  numerical_checks <- function(fit, fixed = NULL) {
    g <- fit$gradients
    maxg <- if (length(g) && all(is.finite(g))) max(abs(g)) else NA_real_
    code <- if (length(fit$model$convergence) == 1L) fit$model$convergence else NA_integer_
    pd <- isTRUE(fit$sd_report$pdHess)
    se_ok <- !is.null(fixed) && nrow(fixed) > 0L &&
      all(is.finite(fixed$estimate)) && all(is.finite(fixed$std.error) & fixed$std.error > 0)
    data.frame(
      convergence_code = code,
      positive_definite_Hessian = pd,
      max_abs_gradient = maxg,
      finite_fixed_SE = se_ok,
      numerical_ok = isTRUE(code == 0L) && pd && is.finite(maxg) && maxg < .001 && se_ok
    )
  }

  verify_source <- function(analysis) {
    s <- source_paths(analysis)
    need(all(file.exists(unlist(s[c("input", "fit", "complete")]))),
         paste("Missing onset-adjusted source files for", analysis))
    inp <- readRDS(s$input)
    saved <- readRDS(s$fit)
    done <- readRDS(s$complete)
    fit <- saved$fit
    need(inherits(fit, "sdmTMB"), "Source fit is not sdmTMB.")
    need(identical(done$status, "FIT_OK"), paste("Source fit not FIT_OK:", analysis))
    need(identical(saved$signature, inp$signature), "Source fit/input signature mismatch.")
    need(identical(fit$family$family, "gaussian") &&
           identical(fit$family$link, "identity") && !isTRUE(fit$reml),
         "Unexpected source family/link/REML.")
    fx <- fixed_table(fit, inp$formula, inp$data)
    need(numerical_checks(fit, fx)$numerical_ok,
         paste("Source fit fails numerical checks:", analysis))
    need("ONSET_mean_z" %in% fx$term, "Source fit lacks annual onset control.")
    list(paths = s, input = inp, fit = fit, fixed = fx)
  }

  read_start <- function(fit) {
    p <- fit$tmb_obj$env$parList(par = fit$tmb_obj$env$last.par.best)
    need(all(START_KEYS %in% names(p)), "Source start list incomplete.")
    p <- p[START_KEYS]
    need(all(vapply(p, function(x) is.numeric(x) && all(is.finite(x)), logical(1))),
         "Source start list contains non-finite values.")
    p
  }

  join_qc <- function(d, phenology_file, out) {
    keys <- c("SPECIES", "SITE_ID", "YEAR")
    tab <- data.table::fread(
      phenology_file,
      select = c("SPECIES", "SITE_ID", "YEAR", "ONSET_mean",
                 "OFFSET_mean", "LAST_PEAK", "N_PEAKS"),
      colClasses = list(character = c("SPECIES", "SITE_ID")),
      data.table = FALSE, nThread = 1L, showProgress = FALSE
    )
    tab$YEAR <- as.integer(tab$YEAR)
    tab <- unique(tab)
    dup <- duplicated(tab[keys]) | duplicated(tab[keys], fromLast = TRUE)
    if (any(dup)) {
      csv(tab[dup, , drop = FALSE], file.path(out, "conflicting_phenology_records.csv"))
      stop("Conflicting phenology rows; see audit.", call. = FALSE)
    }

    req <- data.frame(
      SPECIES = as.character(d$SPECIES), SITE_ID = as.character(d$SITE_ID),
      YEAR = as.integer(d$YEAR), original_row = seq_len(nrow(d)),
      stringsAsFactors = FALSE
    )
    j <- merge(req, tab, by = keys, all.x = TRUE, sort = FALSE)
    j <- j[order(j$original_row), , drop = FALSE]
    need(nrow(j) == nrow(d) && identical(j$original_row, seq_len(nrow(d))),
         "Phenology join changed source row count/order.")
    need(all(is.finite(j$OFFSET_mean)), "Some source rows lack OFFSET_mean.")
    need(max(abs(j$OFFSET_mean - d$OFFSET_mean)) < 1e-7,
         "Phenology CSV OFFSET_mean differs from source fitted data.")
    need(all(is.finite(j$ONSET_mean)), "Some source offset rows lack ONSET_mean.")

    invalid <- is.finite(j$LAST_PEAK) & j$OFFSET_mean < j$LAST_PEAK
    j$qc_last_peak_valid <- !invalid
    j$gap_last_peak_minus_offset <- j$LAST_PEAK - j$OFFSET_mean
    csv(j[invalid, , drop = FALSE], file.path(out, "excluded_offset_before_last_peak.csv"))

    sm <- data.frame(
      n_source_rows = nrow(d),
      n_excluded = sum(invalid),
      percent_excluded = 100 * mean(invalid),
      n_retained = sum(!invalid),
      n_with_LAST_PEAK = sum(is.finite(j$LAST_PEAK)),
      n_excluded_multipeak = sum(invalid & is.finite(j$N_PEAKS) & j$N_PEAKS >= 2),
      percent_excluded_multipeak = if (sum(invalid))
        100 * mean(j$N_PEAKS[invalid] >= 2, na.rm = TRUE) else NA_real_,
      max_gap_days = if (sum(invalid))
        max(j$LAST_PEAK[invalid] - j$OFFSET_mean[invalid]) else NA_real_
    )
    csv(sm, file.path(out, "qc_filter_summary.csv"))
    list(keep = !invalid, joined = j, summary = sm)
  }

  group_signature <- function(d) list(
    species = sort(unique(as.character(d$SPECIES))),
    sites = sort(unique(as.character(d$SITE_ID))),
    years = sort(unique(as.integer(d$YEAR))),
    site_year = sort(unique(as.character(d$site_year_id)))
  )

  same_groups <- function(a, b) {
    identical(a$species, b$species) && identical(a$sites, b$sites) &&
      identical(a$years, b$years) && identical(a$site_year, b$site_year)
  }

  rescale_onset <- function(source, d, warm_ok) {
    sc <- source$input$onset_scaling
    old_c <- as.numeric(sc$centre[1]); old_sd <- as.numeric(sc$SD[1])
    new_c <- mean(d$ONSET_mean); new_sd <- stats::sd(d$ONSET_mean)
    need(all(is.finite(c(old_c, old_sd, new_c, new_sd))) && old_sd > 0 && new_sd > 0,
         "Invalid onset scaling.")
    d$ONSET_mean_z <- (d$ONSET_mean - new_c) / new_sd

    start <- NULL
    if (warm_ok) {
      start <- read_start(source$fit)
      b <- stats::setNames(source$fixed$estimate, source$fixed$term)
      need(all(c("(Intercept)", "ONSET_mean_z") %in% names(b)),
           "Cannot transform onset warm start.")
      b_on <- b[["ONSET_mean_z"]]
      b[["ONSET_mean_z"]] <- b_on * new_sd / old_sd
      b[["(Intercept)"]] <- b[["(Intercept)"]] + b_on * (new_c - old_c) / old_sd
      Xn <- colnames(stats::model.matrix(lme4::nobars(source$input$formula),
        d[seq_len(min(5L, nrow(d))), , drop = FALSE]))
      need(setequal(names(b), Xn), "Filtered fixed-effect design changed unexpectedly.")
      start$b_j <- unname(b[Xn])
    }

    list(
      data = d, start = start,
      scaling = data.frame(old_centre = old_c, old_SD = old_sd,
                           new_centre = new_c, new_SD = new_sd,
                           n = nrow(d))
    )
  }

  project_mesh <- function(d, source_fit) {
    need(inherits(source_fit$spde, "sdmTMBmesh"), "Source mesh object invalid.")
    m <- sdmTMB::make_mesh(d, xy_cols = c("x_km", "y_km"), mesh = source_fit$spde$mesh)
    need(nrow(m$loc_xy) == nrow(d), "Filtered mesh projection row mismatch.")
    need(same_num(as.numeric(as.matrix(m$loc_xy)),
                  as.numeric(as.matrix(d[, c("x_km", "y_km")]))),
         "Filtered mesh projection coordinates do not align.")
    m
  }

  prepare_worker <- function(cfg, analysis) {
    check_packages()
    out <- file.path(cfg$run, analysis)
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    done <- file.path(out, "preparation_complete.rds")
    if (file.exists(done)) {
      x <- readRDS(done)
      need(identical(x$signature, cfg$signature), "Preparation cache mismatch.")
      return("PREPARED_CACHED")
    }

    phase(cfg, "PREPARING", analysis, "Reading current onset-adjusted FIT_OK source")
    source <- verify_source(analysis)
    d0 <- source$input$data

    phase(cfg, "FILTERING_QC", analysis,
          "Exclude only finite LAST_PEAK > OFFSET_mean")
    qc <- join_qc(d0, cfg$phenology_file, out)
    d <- d0[qc$keep, , drop = FALSE]
    need(nrow(d) > 0L && nrow(d) < nrow(d0),
         paste("QC removed zero/all rows for", analysis, "; review filter."))

    # Rebuild grouping factors on retained rows.
    d$SITE_ID <- factor(as.character(d$SITE_ID))
    d$SPECIES <- factor(as.character(d$SPECIES))
    d$SPECIES_slope <- d$SPECIES
    d$YEAR <- as.integer(d$YEAR)
    d$site_year_id <- interaction(d$SITE_ID, d$YEAR, drop = TRUE, lex.order = TRUE)

    warm_ok <- same_groups(group_signature(d0), group_signature(d))
    rs <- rescale_onset(source, d, warm_ok)
    d <- rs$data

    X <- stats::model.matrix(lme4::nobars(source$input$formula), d)
    need(nrow(X) == nrow(d) && all(is.finite(X)) && qr(X)$rank == ncol(X),
         "Filtered fixed-effect design invalid/rank deficient.")

    mesh <- project_mesh(d, source$fit)

    input <- list(
      signature = cfg$signature,
      analysis = analysis,
      data = d,
      formula = source$input$formula,
      mesh = mesh,
      start = rs$start,
      use_full_warm_start = warm_ok,
      onset_scaling = rs$scaling,
      source_fixed = source$fixed,
      source_fit_path = normalizePath(source$paths$fit, winslash = "/", mustWork = TRUE),
      source_input_path = normalizePath(source$paths$input, winslash = "/", mustWork = TRUE),
      source_n = nrow(d0),
      qc = qc$summary
    )

    csv(rs$scaling, file.path(out, "onset_scaling_after_qc.csv"))
    csv(data.frame(full_warm_start_used = warm_ok),
        file.path(out, "warm_start_check.csv"))
    writeLines(c(
      "QC rule: retain LAST_PEAK missing OR OFFSET_mean >= LAST_PEAK.",
      "OFFSET_mean is not recalculated.",
      paste("Source rows:", nrow(d0)),
      paste("Retained rows:", nrow(d)),
      paste("Excluded rows:", nrow(d0) - nrow(d)),
      paste("Full warm start used:", warm_ok),
      paste("Formula:", ftxt(input$formula)),
      "Same 15-km triangulation; filtered observations reprojected.",
      "ONSET_mean remains additive and is re-centred/scaled on retained rows.",
      "All model parameters are re-estimated."
    ), file.path(out, "model_specification.txt"))

    atomic_rds(input, file.path(out, "input.rds"))
    atomic_rds(list(signature = cfg$signature, time = Sys.time()), done)
    paste0("PREPARED_", nrow(d), "_ROWS")
  }

  export_fit <- function(fit, input, out) {
    fx <- fixed_table(fit, input$formula, input$data)
    ch <- numerical_checks(fit, fx)
    csv(ch, file.path(out, "fit_checks.csv"))
    csv(fx, file.path(out, "fixed_effects.csv"))

    comp <- merge(input$source_fixed, fx, by = "term", all = TRUE,
                  suffixes = c("_before_QC", "_after_QC"), sort = FALSE)
    comp$estimate_change <- comp$estimate_after_QC - comp$estimate_before_QC
    comp$SE_ratio <- comp$std.error_after_QC / comp$std.error_before_QC
    csv(comp, file.path(out, "coefficient_comparison.csv"))

    focal <- fx[grepl("clim_anomaly_tw90", fx$term, fixed = TRUE) |
                  fx$term == "ONSET_mean_z", , drop = FALSE]
    csv(focal, file.path(out, "focal_coefficients.csv"))

    tryCatch(
      csv(as.data.frame(sdmTMB::tidy(fit, effects = "ran_pars", conf.int = TRUE)),
          file.path(out, "random_parameters.csv")),
      error = function(e) writeLines(conditionMessage(e),
        file.path(out, "random_parameter_export_error.txt"))
    )
    writeLines(tryCatch(capture.output(sdmTMB::sanity(fit, gradient_thresh = .001)),
                        error = function(e) paste("Sanity export failed:", conditionMessage(e))),
               file.path(out, "sanity.txt"))
    ch
  }

  fit_worker <- function(cfg, analysis) {
    check_packages()
    options(sdmTMB.cores = 1L)
    TMB::openmp(n = 1L, DLL = "sdmTMB")
    out <- file.path(cfg$run, analysis)
    input <- readRDS(file.path(out, "input.rds"))
    need(identical(input$signature, cfg$signature), "Prepared input signature mismatch.")

    final <- file.path(out, "fit_qc.rds")
    complete <- file.path(out, "completed.rds")
    if (file.exists(complete) && file.exists(final)) {
      x <- readRDS(complete)
      need(identical(x$signature, cfg$signature), "Completed cache mismatch.")
      return(x$status)
    }
    if (file.exists(final)) {
      saved <- readRDS(final)
      need(identical(saved$signature, cfg$signature), "Saved fit signature mismatch.")
      ch <- export_fit(saved$fit, input, out)
      if (isTRUE(ch$numerical_ok)) {
        atomic_rds(list(signature = cfg$signature, status = "FIT_OK", time = Sys.time()), complete)
        return("FIT_OK")
      }
      return("NUMERICAL_REVIEW_REQUIRED")
    }

    phase(cfg, "FITTING", analysis,
          paste(nrow(input$data), "retained rows after LAST_PEAK QC"))
    warnings <- character()
    catcher <- function(w) {
      warnings <<- unique(c(warnings, conditionMessage(w)))
      message("WARNING: ", conditionMessage(w)); invokeRestart("muffleWarning")
    }

    args <- list(multiphase = FALSE, nlminb_loops = 1L, newton_loops = 0L,
                 parallel = 1L, get_joint_precision = FALSE,
                 collapse_spatial_variance = FALSE)
    if (isTRUE(input$use_full_warm_start)) args$start <- input$start
    control <- do.call(sdmTMB::sdmTMBcontrol, args)

    set.seed(cfg$fit_seed)
    t0 <- Sys.time()
    fit <- withCallingHandlers(
      sdmTMB::sdmTMB(
        formula = input$formula, data = input$data, mesh = input$mesh,
        time = "YEAR", family = stats::gaussian(link = "identity"),
        spatial = "off", spatiotemporal = "iid", reml = FALSE,
        silent = FALSE, control = control
      ), warning = catcher
    )

    record <- function(f, attempt) list(
      signature = cfg$signature, analysis = analysis, fit = f,
      attempt = attempt, warnings = unique(warnings), qc = input$qc,
      use_full_warm_start = input$use_full_warm_start,
      elapsed_minutes = as.numeric(difftime(Sys.time(), t0, units = "mins"))
    )
    atomic_rds(record(fit, 0L), file.path(out, "fit_initial.rds"))

    fx <- tryCatch(fixed_table(fit, input$formula, input$data), error = function(e) NULL)
    ch <- numerical_checks(fit, fx)
    hist <- cbind(attempt = 0L, objective = fit$model$objective, ch)
    csv(hist, file.path(out, "optimization_history.csv"))

    for (attempt in seq_len(cfg$max_extra_newton)) {
      if (isTRUE(ch$numerical_ok)) break
      need(!any(is.finite(fit$lower)) && !any(is.finite(fit$upper)),
           "Finite optimizer bounds: stopping before Newton refinement.")
      phase(cfg, "EXTRA_OPTIMIZATION", analysis,
            paste("Newton step", attempt, "of at most", cfg$max_extra_newton))
      nf <- tryCatch(
        withCallingHandlers(sdmTMB::run_extra_optimization(
          fit, nlminb_loops = 0L, newton_loops = 1L), warning = catcher),
        error = function(e) { warnings <<- c(warnings, conditionMessage(e)); NULL }
      )
      if (is.null(nf)) break
      atomic_rds(record(nf, attempt),
                 file.path(out, sprintf("fit_newton_%02d.rds", attempt)))
      if (!is.finite(nf$model$objective) ||
          nf$model$objective > fit$model$objective + .001) break
      fit <- nf
      fx <- tryCatch(fixed_table(fit, input$formula, input$data), error = function(e) NULL)
      ch <- numerical_checks(fit, fx)
      hist <- rbind(hist, cbind(attempt = attempt, objective = fit$model$objective, ch))
      csv(hist, file.path(out, "optimization_history.csv"))
    }

    fit$tmb_obj$fn(fit$model$par)
    fit$parlist <- fit$tmb_obj$env$parList(par = fit$tmb_obj$env$last.par.best)
    fit$last.par.best <- fit$tmb_obj$env$last.par.best
    atomic_rds(record(fit, max(hist$attempt)), final)
    writeLines(unique(warnings), file.path(out, "warnings.txt"))

    phase(cfg, "EXPORTING", analysis, "Comparing QC fit with current onset-adjusted source")
    ch <- export_fit(fit, input, out)
    status <- if (isTRUE(ch$numerical_ok)) "FIT_OK" else "NUMERICAL_REVIEW_REQUIRED"
    if (status == "FIT_OK") atomic_rds(
      list(signature = cfg$signature, status = status, time = Sys.time()), complete
    )
    status
  }

  collect <- function(cfg) {
    for (nm in c("qc_filter_summary", "fit_checks", "coefficient_comparison", "focal_coefficients")) {
      rows <- lapply(ANALYSES, function(a) {
        p <- file.path(cfg$run, a, paste0(nm, ".csv"))
        if (!file.exists(p)) return(NULL)
        cbind(analysis = a,
              utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE))
      })
      rows <- Filter(Negate(is.null), rows)
      if (length(rows)) csv(do.call(rbind, rows), file.path(cfg$run, paste0("all_", nm, ".csv")))
    }
    if (requireNamespace("zip", quietly = TRUE)) tryCatch({
      rel <- list.files(cfg$run, recursive = TRUE, full.names = FALSE)
      keep <- grepl("\\.(csv|txt|log|R)$", rel, ignore.case = TRUE)
      zip::zipr(file.path(cfg$run, "results_to_review.zip"), rel[keep],
                root = cfg$run, include_directories = FALSE)
    }, error = function(e) message("Optional review ZIP failed: ", conditionMessage(e)))
    invisible(NULL)
  }

  invoke_worker <- function(cfg, fun, analysis) {
    log <- file.path(cfg$run, analysis,
                     paste0(fun, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
    dir.create(dirname(log), recursive = TRUE, showWarnings = FALSE)
    child <- callr::r_bg(function(engine, cfg, fun, analysis) {
      e <- new.env(parent = globalenv()); sys.source(engine, envir = e)
      do.call(e[[fun]], list(cfg = cfg, analysis = analysis))
    }, args = list(cfg$engine, cfg, fun, analysis), libpath = cfg$libpath,
    stdout = log, stderr = "2>&1", supervise = TRUE,
    user_profile = FALSE, system_profile = FALSE,
    env = c(callr::rcmd_safe_env(), OMP_NUM_THREADS = "1",
            OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1"))
    on.exit(if (child$is_alive()) child$kill(), add = TRUE)
    child$wait(); child$get_result()
  }

  controller <- function(cfg) {
    on.exit(unlink(cfg$lock, recursive = TRUE), add = TRUE)
    st <- data.frame(analysis = ANALYSES, preparation = "PENDING",
                     fit = "PENDING", message = "", stringsAsFactors = FALSE)
    save_status <- function() csv(st, file.path(cfg$run, "workflow_status.csv"))
    save_status()
    tryCatch({
      check_packages(); phase(cfg, "PREFLIGHT", detail = "QC specification validated")
      for (i in seq_along(ANALYSES)) {
        st$preparation[i] <- "RUNNING"; save_status()
        ans <- tryCatch(invoke_worker(cfg, "prepare_worker", ANALYSES[i]),
          error = function(e) {
            st$preparation[i] <<- "PREPARATION_FAILED"
            st$message[i] <<- conditionMessage(e)
            st$fit[st$fit == "PENDING"] <<- "NOT_STARTED"
            save_status(); stop(conditionMessage(e), call. = FALSE)
          })
        st$preparation[i] <- ans; save_status()
      }
      if (cfg$prepare_only) {
        phase(cfg, "PREPARED_ONLY", detail = "QC counts exported; no fits run")
        return(invisible(st))
      }
      for (i in seq_along(ANALYSES)) {
        st$fit[i] <- "RUNNING"; save_status()
        ans <- tryCatch(invoke_worker(cfg, "fit_worker", ANALYSES[i]),
          error = function(e) { st$message[i] <<- conditionMessage(e); "REVIEW_OR_ERROR" })
        st$fit[i] <- ans; save_status(); collect(cfg)
      }
      ok <- all(st$fit == "FIT_OK")
      phase(cfg, if (ok) "FITS_COMPLETE" else "FINISHED_WITH_ISSUES",
            detail = if (ok) "Review QC effect sizes before downstream extraction"
                     else "Inspect logs/checkpoints before downstream work")
      invisible(st)
    }, error = function(e) {
      writeLines(conditionMessage(e), file.path(cfg$run, "controller_error.txt"))
      phase(cfg, "STOPPED", detail = conditionMessage(e)); invisible(st)
    }, finally = { save_status(); collect(cfg) })
  }

  engine_file <- function(run) {
    env <- environment(engine_file)
    n <- ls(env, all.names = TRUE)
    fun <- n[vapply(n, function(x) is.function(get(x, envir = env)), logical(1))]
    path <- file.path(run, "workflow_engine.R")
    dump(c("VERSION", "ROOT", "SOURCE_RUN", "ANALYSES", "CORE", "START_KEYS", fun),
         file = path, envir = env, control = "all")
    normalizePath(path, winslash = "/", mustWork = TRUE)
  }

  status <- function(job = state$job) {
    need(!is.null(job), "No QC offset job in this R session.")
    message("Elapsed: ", round(as.numeric(difftime(Sys.time(), job$started, units = "hours")), 2),
            " h | controller running: ", job$process$is_alive())
    p <- file.path(job$run, "progress.rds")
    if (file.exists(p)) {
      x <- tryCatch(readRDS(p), error = function(e) NULL)
      if (!is.null(x)) message(x$stage, " | ", x$analysis, " | ", x$detail)
    }
    w <- file.path(job$run, "workflow_status.csv")
    if (file.exists(w)) try(print(utils::read.csv(w, stringsAsFactors = FALSE), row.names = FALSE), silent = TRUE)
    message("Output: ", job$run)
    invisible(job)
  }

  watch <- function(job = state$job, every = 60) {
    need(!is.null(job), "No QC offset job in this R session.")
    tryCatch({ repeat { status(job); if (!job$process$is_alive()) break; Sys.sleep(every) } },
      interrupt = function(e) message("Monitoring paused; queue continues. Use pheno_offset_qc$watch()."))
    invisible(job)
  }

  stop_job <- function(job = state$job) {
    need(!is.null(job), "No QC offset job in this R session.")
    if (job$process$is_alive()) job$process$kill()
    message("QC queue stopped. Completed checkpoints retained.")
    invisible(job)
  }

  run <- function(root = ROOT, phenology_file = NULL, prepare_only = FALSE,
                  max_extra_newton = 3L, monitor = TRUE) {
    check_packages(); need(dir.exists(root), paste("Project root missing:", root))
    need(dir.exists(SOURCE_RUN), paste("Source onset-adjusted run missing:", SOURCE_RUN))
    disk_check(root, 50)
    if (is.null(phenology_file)) phenology_file <- file.path(root, "output", "pheno_estimates_allspp.csv")
    need(file.exists(phenology_file), paste("Phenology file missing:", phenology_file))
    phenology_file <- normalizePath(phenology_file, winslash = "/", mustWork = TRUE)

    source_files <- unlist(lapply(ANALYSES, function(a) {
      s <- source_paths(a); c(s$input, s$fit, s$complete)
    }))
    inputs <- unique(c(phenology_file, source_files))
    hashes <- md5(inputs); need(!anyNA(hashes), "Cannot fingerprint inputs.")
    sig <- digest::digest(list(version = VERSION, files = inputs, md5 = hashes,
      packages = versions(), qc_rule = "is.na(LAST_PEAK) | OFFSET_mean >= LAST_PEAK",
      max_extra_newton = max_extra_newton), algo = "sha256")

    outroot <- file.path(root, "output", "phenology_plasticity",
                         "offset_with_onset_qc_lastpeak_15km")
    dir.create(outroot, recursive = TRUE, showWarnings = FALSE)
    lock <- file.path(outroot, "queue.lock")
    if (dir.exists(lock)) {
      ownerfile <- file.path(lock, "owner.rds")
      owner <- if (file.exists(ownerfile)) tryCatch(readRDS(ownerfile), error = function(e) NULL) else NULL
      alive <- if (is.null(owner)) FALSE else tryCatch(
        ps::ps_is_running(ps::ps_handle(owner$pid, time = owner$create_time)), error = function(e) FALSE)
      need(!alive, paste("A QC offset queue may already be running:", lock))
      unlink(lock, recursive = TRUE)
    }
    need(dir.create(lock), "Cannot acquire QC queue lock.")
    launched <- FALSE; on.exit(if (!launched) unlink(lock, recursive = TRUE), add = TRUE)

    run_dir <- file.path(outroot, paste0("run_", substr(sig, 1, 12)))
    dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
    cfg <- list(signature = sig, root = root, run = run_dir, lock = lock,
      phenology_file = phenology_file, prepare_only = prepare_only,
      max_extra_newton = as.integer(max_extra_newton), fit_seed = 20261005L,
      libpath = .libPaths(), input_paths = inputs, input_md5 = hashes)
    cfg$engine <- engine_file(run_dir)
    atomic_rds(cfg, file.path(run_dir, "run_config.rds"))
    csv(data.frame(file = inputs, md5 = hashes), file.path(run_dir, "source_manifest.csv"))
    writeLines(c(
      "QC-FILTERED OFFSET REFITS WITH ANNUAL ONSET CONTROL",
      "QC rule: exclude finite LAST_PEAK > OFFSET_mean; LAST_PEAK NA retained.",
      "OFFSET_mean is not recalculated.",
      "Source = current FIT_OK onset-adjusted 15-km offset models.",
      "Same 90-day climate window, fixed/random effects and 15-km triangulation.",
      "ONSET_mean remains additive; re-centred/scaled on retained rows.",
      "All parameters are re-estimated.",
      "No final plasticity extraction or abundance refit is done here."
    ), file.path(run_dir, "README.txt"))

    log <- file.path(run_dir, paste0("controller_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
    proc <- callr::r_bg(function(engine, cfg) {
      e <- new.env(parent = globalenv()); sys.source(engine, envir = e)
      e$need(identical(e$md5(cfg$input_paths), cfg$input_md5), "Inputs changed after configuration.")
      e$controller(cfg)
    }, args = list(cfg$engine, cfg), libpath = cfg$libpath,
    stdout = log, stderr = "2>&1", supervise = TRUE,
    user_profile = FALSE, system_profile = FALSE)

    h <- ps::ps_handle(proc$get_pid())
    atomic_rds(list(pid = proc$get_pid(), create_time = ps::ps_create_time(h)), file.path(lock, "owner.rds"))
    job <- list(process = proc, run = run_dir, started = Sys.time(), log = log)
    state$job <- job; launched <- TRUE
    message("QC-filtered OFFSET queue launched.")
    message("Output: ", run_dir)
    message("No onset/peak/abundance model is changed. Keep R open.")
    if (monitor) watch(job)
    invisible(job)
  }

  list(run = run, status = status, watch = watch, stop = stop_job)
})

message("Functions loaded. Start with:")
message("  offset_qc_job <- pheno_offset_qc$run(monitor = FALSE)")
message("  pheno_offset_qc$status()")
