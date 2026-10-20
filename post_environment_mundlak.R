# phenoIMPACT | Accepted environmental Mundlak fits -> updated demography
# 2026-10-20. Sourcing this file NEVER starts a job.
# See README_post_environment_CA.md for the scientific contract and commands.
# One large R process at a time. Original fits/inputs are read-only.

pheno_post_env <- local({
  VERSION <- "post_environment_mundlak_1.1"
  ROOT <- "E:/phenoIMPACT project/code/phenoIMPACT"
  GROUPS <- c("multivoltine", "univoltine")
  ANALYSES <- c("onset_multivoltine", "offset_multivoltine",
                "onset_univoltine", "offset_univoltine")
  V2_ANALYSES <- c("offset_univoltine", "onset_multivoltine", "onset_univoltine")
  ENGINE <- "10f_mundlak_plasticity_demography_QC_FINAL_spatial_15km_RECOVERY.R"
  state <- new.env(parent = emptyenv())
  need <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
  csv <- function(x, p) utils::write.csv(x, p, row.names = FALSE, na = "")
  same <- function(x, y, tol = 1e-10) {
    length(x) == length(y) && all(is.finite(x)) && all(is.finite(y)) &&
      (length(x) == 0L || max(abs(x - y)) <= tol)
  }
  key <- function(d) {
    need(all(c("SPECIES", "SITE_ID") %in% names(d)), "Missing population keys.")
    need(!anyNA(d[c("SPECIES", "SITE_ID")]), "Missing species/site identity.")
    need(!any(grepl("\r", as.character(d$SPECIES), fixed = TRUE)) &&
           !any(grepl("\r", as.character(d$SITE_ID), fixed = TRUE)), "Invalid key separator.")
    paste(as.character(d$SPECIES), as.character(d$SITE_ID), sep = "\r")
  }
  atomic <- function(x, p) {
    dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
    tmp <- tempfile(".writing_", tmpdir = dirname(p), fileext = ".rds")
    on.exit(unlink(tmp), add = TRUE)
    saveRDS(x, tmp, compress = "gzip")
    bak <- paste0(p, ".previous")
    if (file.exists(bak)) unlink(bak)
    had <- file.exists(p)
    if (had) need(file.rename(p, bak), paste("Cannot back up", p))
    if (!file.rename(tmp, p)) {
      if (had) file.rename(bak, p)
      stop("Cannot save ", p, call. = FALSE)
    }
    if (file.exists(bak)) unlink(bak)
    invisible(p)
  }
  phase <- function(cfg, stage, detail = "", group = "", terminal = FALSE) {
    atomic(list(stage = stage, detail = detail, group = group,
                time = Sys.time(), terminal = terminal), cfg$progress)
    message(format(Sys.time(), "%F %T"), " | ", stage, " | ", group, " | ", detail)
  }
  packages <- function(names) {
    missing <- names[!vapply(names, requireNamespace, logical(1), quietly = TRUE)]
    need(!length(missing), paste("Missing packages:", paste(missing, collapse = ", ")))
  }
  alive <- function(owner) {
    if (is.null(owner$pid) || is.null(owner$create_time)) return(FALSE)
    packages("ps")
    tryCatch({
      h <- ps::ps_handle(owner$pid)
      ps::ps_is_running(h) && isTRUE(abs(as.numeric(ps::ps_create_time(h)) -
        as.numeric(owner$create_time)) < 0.001)
    }, error = function(e) {
      if (inherits(e, "no_such_process")) return(FALSE)
      stop("Cannot verify worker process: ", conditionMessage(e), call. = FALSE)
    })
  }
  active_owner <- function(path) {
    if (!file.exists(path)) return(FALSE)
    owner <- readRDS(path)
    need(is.list(owner) && length(owner$pid) == 1L && is.finite(owner$pid) &&
      length(owner$create_time) == 1L && is.finite(as.numeric(owner$create_time)),
      paste("Invalid worker ownership metadata:", path))
    alive(owner)
  }
  active_run <- function(run) {
    # A corrupt ownership file is not evidence that a worker has stopped.
    active_owner(file.path(run, "worker_owner.rds"))
  }
  refinement_roots <- function(run) {
    v1 <- file.path(run, "postfit_refinement_queue_v1")
    v2 <- file.path(run, "postfit_refinement_queue_v2_isolated")
    cfg_path <- file.path(v2, "refinement_config.rds")
    if (file.exists(cfg_path)) {
      cfg <- readRDS(cfg_path)
      need(is.character(cfg$source_run) && length(cfg$source_run) == 1L &&
        identical(normalizePath(cfg$source_run, winslash = "/", mustWork = TRUE),
                  normalizePath(run, winslash = "/", mustWork = TRUE)),
        "V2 refinement configuration belongs to a different environmental run.")
      need(is.character(cfg$old_refinement) && length(cfg$old_refinement) == 1L &&
        !is.na(cfg$old_refinement) && nzchar(cfg$old_refinement),
        "V2 refinement configuration lacks its v1 source.")
      v1 <- cfg$old_refinement
    }
    list(v1 = v1, v2 = v2)
  }
  environmental_workers <- function(run) {
    roots <- refinement_roots(run)
    paths <- unique(c(
      file.path(run, "worker_owner.rds"),
      file.path(run, "postfit_refinement_queue_v1", "worker_owner.rds"),
      file.path(roots$v1, "worker_owner.rds"),
      file.path(roots$v2, "worker_owner.rds"),
      file.path(roots$v2, ANALYSES, "isolated_worker_owner.rds")))
    paths <- paths[file.exists(paths)]
    data.frame(owner_file = paths,
      active = vapply(paths, active_owner, logical(1), USE.NAMES = FALSE),
      stringsAsFactors = FALSE)
  }
  require_environment_idle <- function(run) {
    need(!any(environmental_workers(run)$active),
         "An environmental/refinement worker is active (original, v1 or v2).")
  }
  file_manifest <- function(paths) {
    paths <- unique(unname(normalizePath(paths, winslash = "/", mustWork = TRUE)))
    z <- file.info(paths)
    need(!anyNA(z$size) && all(!z$isdir), "Source file metadata unavailable.")
    data.frame(path = paths, bytes = z$size, modified = as.numeric(z$mtime),
               stringsAsFactors = FALSE)
  }
  check_manifest <- function(cfg) {
    need(identical(file_manifest(cfg$manifest$path), cfg$manifest),
         "A source changed after preflight. Restart run() to create a new version.")
  }
  accepted_checks <- function(z) {
    is.data.frame(z) && nrow(z) == 1L &&
      isTRUE(z$numerical_checks_ok) && isTRUE(z$convergence_code == 0L) &&
      isTRUE(z$pdHess) && isTRUE(is.finite(z$max_abs_gradient)) &&
      isTRUE(z$max_abs_gradient < 0.001)
  }
  default_env <- function(root) file.path(root, "output", "phenology_plasticity",
    "environment_mundlak_by_voltinism_qc_final_15km", "run_de0890c85d4a")

  # Only small completion/check files are read here; never load a fitted model.
  select_environment <- function(run) {
    roots <- refinement_roots(run)
    v2_exists <- dir.exists(roots$v2)
    rows <- lapply(ANALYSES, function(a) {
      src <- file.path(run, a)
      ref1 <- file.path(roots$v1, a)
      ref2 <- file.path(roots$v2, a)
      # A scheduled v2 model supersedes v1 even before its child directory exists.
      # Offset multivoltine is carried forward from v1 by the isolated v2 queue.
      refinement_version <- if (dir.exists(ref2) || (v2_exists && a %in% V2_ANALYSES)) {
        "v2"
      } else if (dir.exists(ref1) || (v2_exists && a == "offset_multivoltine")) {
        "v1"
      } else "none"
      ref <- if (refinement_version == "v2") ref2 else ref1
      done <- file.path(src, "completed.rds")
      prep <- file.path(src, "prepared_input.rds")
      original <- file.path(src, "fit_mundlak.rds")
      fit <- original
      check <- file.path(src, "fit_checks.csv")
      marker <- done
      kind <- "original"
      status <- "NOT_COMPLETED"
      ready <- FALSE
      original_signature <- NA_character_
      if (file.exists(done)) {
        z <- readRDS(done)
        status <- as.character(z$status)
        original_signature <- z$signature
      }
      # Any refinement artifact takes precedence. Never silently fall back.
      if (refinement_version != "none") {
        kind <- "refined"
        fit <- file.path(ref, "fit_refined_final.rds")
        check <- file.path(ref, "fit_checks_refined.csv")
        marker <- file.path(ref, "refinement_completed.rds")
        status <- if (dir.exists(ref)) "REFINEMENT_NOT_COMPLETED" else "REFINEMENT_NOT_STARTED"
        if (file.exists(marker)) {
          z <- readRDS(marker)
          need(is.data.frame(z$summary) && nrow(z$summary) == 1L &&
                 identical(z$summary$analysis, a),
               paste("Refinement marker identity mismatch:", a))
          status <- as.character(z$summary$refinement_status)
        }
      }
      all_files <- all(file.exists(c(done, prep, fit, check, marker)))
      if (all_files) {
        z <- utils::read.csv(check, stringsAsFactors = FALSE)
        ready <- status %in% c("FIT_OK", "FIT_OK_REFINED", "FIT_OK_AFTER_REBUILD") &&
          accepted_checks(z)
      }
      data.frame(analysis = a, kind = kind, refinement_version = refinement_version,
                 status = status, ready = ready,
                 fit = fit, prepared = prep, checks = check, marker = marker,
                 original_marker = done, original_signature = original_signature,
                 stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
  }
  resolve_demography <- function(root, run = NULL) {
    required <- c("run_config.rds", paste0("prepared_mundlak_", GROUPS, ".rds"),
                  file.path(paste0(GROUPS, "__mundlak"), "fit.rds"))
    if (is.null(run)) {
      base <- file.path(root, "output", "population_trends",
                        "mundlak_geography_plasticity_qc_final_spatial_15km")
      candidates <- list.dirs(base, recursive = FALSE, full.names = TRUE)
      candidates <- candidates[vapply(candidates, function(d)
        all(file.exists(file.path(d, required))), logical(1))]
      need(length(candidates) == 1L, paste(
        "Specify demography_run: expected exactly one complete previous Mundlak run; found",
        length(candidates), paste(candidates, collapse = "\n")))
      run <- candidates
    }
    need(all(file.exists(file.path(run, required))), "Previous demographic run is incomplete.")
    run <- normalizePath(run, winslash = "/", mustWork = TRUE)
    cfg <- readRDS(file.path(run, "run_config.rds"))
    need(is.character(cfg$source_run) && length(cfg$source_run) == 1L,
         "Previous demographic config lacks its likelihood source.")
    list(run = run, likelihood_run = cfg$source_run, paths = file.path(run, required))
  }
  preflight <- function(project_root = ROOT, environment_run = default_env(project_root),
                        demography_run = NULL) {
    env <- select_environment(environment_run)
    workers <- environmental_workers(environment_run)
    active <- any(workers$active)
    demo_error <- NULL
    demo <- tryCatch(resolve_demography(project_root, demography_run), error = function(e) {
      demo_error <<- conditionMessage(e); NULL
    })
    print(env[c("analysis", "kind", "refinement_version", "status", "ready")], row.names = FALSE)
    message("Environmental/refinement worker active: ", active)
    if (active) print(workers[workers$active, , drop = FALSE], row.names = FALSE)
    if (!is.null(demo_error)) message(demo_error)
    if (!is.null(demo)) message("Previous demographic run: ", demo$run)
    message("Metadata only: no model loaded and no process started.")
    invisible(list(ready = all(env$ready) && !active && !is.null(demo),
                   environment = env, active = active, workers = workers, demography = demo))
  }

  spec <- function(analysis) {
    need(analysis %in% ANALYSES, "Unknown environmental analysis.")
    onset <- startsWith(analysis, "onset_")
    group <- sub("^(onset|offset)_", "", analysis)
    list(analysis = analysis, group = group,
      response = if (onset) "ONSET_mean" else "OFFSET_mean",
      anomaly = if (onset) "clim_anomaly_tw60" else "clim_anomaly_tw90",
      control = if (onset) character() else "ONSET_mean_z",
      orientation = if (onset || group == "univoltine") -1 else 1,
      moderators = as.vector(outer(c("photoperiod", "background", "predictability", "trend"),
                                    c("within_species", "species_mean"), paste, sep = "_")))
  }
  # Saved coefficients/reports only: no sdmTMB reconstruction or new optimization.
  coefficients <- function(fit, s) {
    sf <- fit$split_formula[[1]]
    expected <- stats::as.formula(paste(s$response, "~", s$anomaly, "*(",
      paste(s$moderators, collapse = "+"), ")",
      if (length(s$control)) paste("+", s$control) else ""))
    need(setequal(attr(stats::terms(sf$form_no_bars), "term.labels"),
                  attr(stats::terms(expected), "term.labels")), "Unexpected environmental formula.")
    cn <- sf$re_cov_terms$cnms
    need(length(cn) == 4L && setequal(names(cn),
      c("SITE_ID", "SPECIES", "SPECIES_slope", "site_year_id")) &&
      identical(cn$SPECIES_slope, s$anomaly), "Unexpected random-effects structure.")
    need(same(fit$model$par, fit$sd_report$par.fixed, 1e-6) &&
           identical(names(fit$model$par), names(fit$sd_report$par.fixed)),
         "Saved covariance belongs to different point estimates.")
    fixed_names <- colnames(stats::model.matrix(sf$form_no_bars,
                                               fit$data[1L, , drop = FALSE]))
    ix <- which(names(fit$sd_report$par.fixed) == "b_j")
    need(length(ix) == length(fixed_names) && !anyDuplicated(fixed_names),
         "Cannot map fixed coefficients.")
    beta <- stats::setNames(as.numeric(fit$sd_report$par.fixed[ix]), fixed_names)
    V <- as.matrix(fit$sd_report$cov.fixed[ix, ix, drop = FALSE])
    need(all(is.finite(beta)) && all(is.finite(V)) && all(diag(V) > 0),
         "Invalid fixed estimates/covariance.")
    dimnames(V) <- list(fixed_names, fixed_names)
    blocks <- sf$re_cov_terms$re_b_df
    need(all(c("start", "end", "group_indices", "level") %in% names(blocks)),
         "Missing random-effect mapping.")
    ids <- unlist(Map(seq.int, blocks$start, blocks$end), use.names = FALSE)
    report_ids <- grepl("^re_b_pars($|\\[)", names(fit$sd_report$value))
    values <- as.numeric(fit$sd_report$value[report_ids])
    errors <- as.numeric(fit$sd_report$sd[report_ids])
    active <- !is.na(fit$tmb_map$re_b_pars)
    need(length(values) == length(active) && length(errors) == length(active) &&
           identical(as.integer(ids + 1L), which(as.vector(active))),
         "Random-effect report/mapping mismatch.")
    need(same(as.numeric(fit$parlist$re_b_pars)[active], values[active], 1e-6),
         "Random-effect report and saved modes disagree.")
    b <- blocks[blocks$group_indices == match("SPECIES_slope", names(cn)), , drop = FALSE]
    need(nrow(b) > 0 && all(b$start == b$end) && !anyDuplicated(as.character(b$level)),
         "Expected one anomaly random slope per species.")
    slopes <- data.frame(SPECIES = as.character(b$level),
      species_slope_deviation = values[b$start + 1L],
      species_slope_marginal_SE = errors[b$start + 1L], stringsAsFactors = FALSE)
    need(setequal(slopes$SPECIES, as.character(fit$data$SPECIES)) &&
           all(is.finite(slopes$species_slope_deviation)), "Species slope identities invalid.")
    list(beta = beta, V = V, slopes = slopes, formula = sf$form_no_bars)
  }
  population_slopes <- function(d, coef, s) {
    needed <- c("SPECIES", "SITE_ID", "bms_id", "YEAR", s$moderators, s$anomaly)
    need(all(needed %in% names(d)) && !anyNA(d[needed]), "Missing environmental data.")
    need(!anyDuplicated(d[c("SPECIES", "SITE_ID", "YEAR")]), "Duplicated population-years.")
    k <- key(d)
    first <- !duplicated(k)
    pop <- as.data.frame(d[first, c("SPECIES", "SITE_ID", "bms_id", s$moderators)])
    for (v in c("SPECIES", "SITE_ID", "bms_id")) pop[[v]] <- as.character(pop[[v]])
    ix <- match(k, key(pop))
    need(identical(as.character(d$bms_id), pop$bms_id[ix]), "BMS changes within a population.")
    for (v in s$moderators) need(same(d[[v]], pop[[v]][ix]),
      paste("Environmental moderator varies across years:", v))
    pop$n_years <- tabulate(ix, nbins = nrow(pop))
    pop$pop_id <- paste(pop$SPECIES, pop$SITE_ID, sep = "_")
    need(!anyDuplicated(pop$pop_id), "Ambiguous legacy population IDs.")
    C <- matrix(0, nrow(pop), length(coef$beta), dimnames = list(NULL, names(coef$beta)))
    need(s$anomaly %in% names(coef$beta), "Missing main anomaly effect.")
    C[, s$anomaly] <- 1
    for (v in s$moderators) {
      term <- intersect(c(paste(s$anomaly, v, sep = ":"), paste(v, s$anomaly, sep = ":")),
                        names(coef$beta))
      need(length(term) == 1L, paste("Missing/ambiguous anomaly interaction:", v))
      C[, term] <- pop[[v]]
    }
    # Independent design-matrix derivative check, holding annual onset constant.
    zero <- one <- pop
    zero[[s$anomaly]] <- 0; one[[s$anomaly]] <- 1
    for (v in s$control) { zero[[v]] <- 0; one[[v]] <- 0 }
    rhs <- stats::delete.response(stats::terms(coef$formula))
    derivative <- stats::model.matrix(rhs, one) - stats::model.matrix(rhs, zero)
    need(identical(colnames(derivative), colnames(C)) && same(derivative, C, 1e-9),
         "Analytical slope differs from the design-matrix derivative.")
    sp <- match(pop$SPECIES, coef$slopes$SPECIES)
    need(!anyNA(sp), "Species random slopes did not join.")
    pop$fixed_contextual_slope <- as.numeric(C %*% coef$beta)
    pop$fixed_component_SE <- sqrt(pmax(0, rowSums((C %*% coef$V) * C)))
    pop$species_slope_deviation <- coef$slopes$species_slope_deviation[sp]
    pop$species_slope_marginal_SE <- coef$slopes$species_slope_marginal_SE[sp]
    pop$plasticity_raw_days_per_C <- pop$fixed_contextual_slope + pop$species_slope_deviation
    pop$plasticity_bio_days_per_C <- s$orientation * pop$plasticity_raw_days_per_C
    pop$analysis <- s$analysis; pop$voltinism <- s$group
    pop$biological_orientation <- if (s$orientation == -1) "advancement" else "delay"
    need(all(is.finite(pop$plasticity_raw_days_per_C)), "Non-finite contextual plasticity.")
    # Separate SEs are not a total-slope SE: fixed/random covariance is not available.
    pop
  }
  extract_one <- function(cfg, analysis) {
    check_manifest(cfg)
    row <- cfg$environment[cfg$environment$analysis == analysis, , drop = FALSE]
    out <- file.path(cfg$out, "extraction", analysis)
    done <- file.path(out, "completed.rds")
    if (file.exists(done) && file.exists(file.path(out, "extraction.rds")) &&
        file.exists(file.path(out, "population_plasticity.csv"))) {
      need(identical(readRDS(done)$signature, cfg$signature), "Extraction cache mismatch.")
      return(invisible(out))
    }
    packages("sdmTMB")
    need(as.character(utils::packageVersion("sdmTMB")) == "1.1.0",
         "Use the source sdmTMB 1.1.0 library for these saved fits.")
    phase(cfg, "EXTRACTING", analysis)
    p <- readRDS(row$prepared); saved <- readRDS(row$fit); fit <- saved$fit
    need(inherits(fit, "sdmTMB") && identical(saved$analysis, analysis), "Wrong environmental fit.")
    source_signature <- if (row$kind == "refined") saved$source_signature else saved$signature
    need(identical(source_signature, row$original_signature), "Environmental provenance mismatch.")
    checks <- if (row$kind == "refined") saved$checks else utils::read.csv(row$checks)
    need(accepted_checks(checks) && isTRUE(fit$model$convergence == 0L) &&
           isTRUE(fit$sd_report$pdHess), "Environmental fit has not passed numerical checks.")
    need(identical(fit$family$family, "gaussian") && identical(fit$family$link, "identity"),
         "Unexpected environmental likelihood.")
    s <- spec(analysis)
    need(identical(p$analysis, analysis) && identical(p$spec$anomaly, s$anomaly) &&
           identical(p$spec$orientation, s$orientation), "Environmental specification changed.")
    cols <- c("SPECIES", "SITE_ID", "YEAR", "bms_id", s$response, s$anomaly,
              s$moderators, s$control)
    need(all(cols %in% names(p$data)) && all(cols %in% names(fit$data)), "Fit/input columns missing.")
    need(nrow(p$data) == nrow(fit$data), "Fit/input row count mismatch.")
    for (v in cols) {
      ok <- if (is.numeric(p$data[[v]])) same(p$data[[v]], fit$data[[v]], 0) else
        identical(as.character(p$data[[v]]), as.character(fit$data[[v]]))
      need(ok, paste("Fit/input mismatch in", v))
    }
    co <- coefficients(fit, s)
    pop <- population_slopes(as.data.frame(fit$data), co, s)
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    csv(pop, file.path(out, "population_plasticity.csv"))
    csv(co$slopes, file.path(out, "species_slopes.csv"))
    csv(data.frame(term = names(co$beta), estimate = as.numeric(co$beta),
                   SE = sqrt(diag(co$V))), file.path(out, "fixed_coefficients.csv"))
    atomic(list(signature = cfg$signature, populations = pop, beta = co$beta, V = co$V,
                species_slopes = co$slopes, source = row,
                uncertainty_propagated_to_demography = FALSE), file.path(out, "extraction.rds"))
    atomic(list(signature = cfg$signature, analysis = analysis, n_populations = nrow(pop),
                time = Sys.time()), done)
    invisible(out)
  }

  scale_plasticity <- function(pop) {
    need(!anyDuplicated(key(pop)), "Duplicate extracted populations.")
    value <- pop$plasticity_bio_days_per_C
    need(length(value) > 1 && all(is.finite(value)), "Invalid extracted plasticity.")
    mu <- mean(value); sd <- stats::sd(value)
    need(is.finite(sd) && sd > 0, "Zero or invalid plasticity SD.")
    pop$z <- (value - mu) / sd
    list(populations = pop, scaling = data.frame(n_populations = nrow(pop), mean = mu,
      SD = sd, universe = "all extracted populations of this analysis, before abundance join"))
  }
  decompose <- function(d, sites) {
    pop <- unique(d[c("pop_id", "SPECIES", "SITE_ID", "onset_plasticity_z",
                      "offset_plasticity_bio_z")])
    need(!anyDuplicated(key(pop)) && !anyDuplicated(pop$pop_id), "Plasticity varies within population.")
    need(!anyDuplicated(as.character(sites$SITE_ID)), "Duplicate coordinate identities.")
    m <- match(as.character(pop$SITE_ID), as.character(sites$SITE_ID))
    need(!anyNA(m), "Missing site coordinates.")
    for (v in c("latitude_deg", "longitude_deg")) {
      pop[[v]] <- sites[[v]][m]
      need(length(pop[[v]]) == nrow(pop) && all(is.finite(pop[[v]])), "Invalid geographic coordinates.")
      axis <- sub("_deg$", "", v)
      pop[[paste0(axis, "_species_mean")]] <- unname(tapply(pop[[v]], pop$SPECIES, mean)[pop$SPECIES])
      pop[[paste0(axis, "_within_species_deg")]] <- pop[[v]] - pop[[paste0(axis, "_species_mean")]]
      pop[[paste0(axis, "_within_species_10deg")]] <- pop[[paste0(axis, "_within_species_deg")]] / 10
    }
    scales <- list()
    for (metric in c("onset", "offset")) {
      v <- if (metric == "onset") "onset_plasticity_z" else "offset_plasticity_bio_z"
      b <- paste0(metric, "_plasticity_species_mean_z")
      w <- paste0(metric, "_plasticity_within_species_z")
      pop[[b]] <- unname(tapply(pop[[v]], pop$SPECIES, mean)[pop$SPECIES])
      pop[[w]] <- pop[[v]] - pop[[b]]
      need(same(pop[[v]], pop[[w]] + pop[[b]]) &&
             max(abs(tapply(pop[[w]], pop$SPECIES, mean))) < 1e-10,
           "Within/between reconstruction failed.")
      scales[[metric]] <- data.frame(
        component = paste0(metric, c("_within_species", "_between_species")),
        fitted_variable = c(w, b),
        SD_for_descriptive_rescaling = c(stats::sd(pop[[w]]), stats::sd(pop[[b]][!duplicated(pop$SPECIES)])),
        SD_weighting = c("equal population", "equal species"),
        fitted_unit = "one updated group-wide plasticity SD; no component re-standardization")
    }
    m <- match(key(d), key(pop))
    for (v in setdiff(names(pop), c("pop_id", "SPECIES", "SITE_ID"))) d[[v]] <- pop[[v]][m]
    X <- stats::model.matrix(~ year_decade * (
      latitude_within_species_10deg + longitude_within_species_10deg +
      onset_plasticity_within_species_z + onset_plasticity_species_mean_z +
      offset_plasticity_within_species_z + offset_plasticity_species_mean_z), d)
    need(all(is.finite(X)) && qr(X)$rank == ncol(X), "Updated demographic design is rank deficient.")
    list(frame = d, populations = pop, X = X, scales = do.call(rbind, scales))
  }
  # Explicitly remap the established Gamma/log TMB schema and latent states.
  # This accepts whole-population loss and keeps all annual rows of retained populations.
  rebuild_input <- function(old, fit, onset, offset, group) {
    need(identical(old$group, group) && identical(fit$group, group) &&
           identical(fit$signature, old$signature), "Old demographic fit/input identity mismatch.")
    need(isTRUE(fit$checks$numerical_checks_ok), "Old demographic warm start is not numerically accepted.")
    schema <- c("y", "X", "time", "species", "site", "year", "previous", "gap",
                "A", "M0", "M1", "M2", "use_spatial")
    need(setequal(names(old$data), schema) && isTRUE(old$data$use_spatial == 1L),
         "Unknown demographic data schema: do not silently subset new data fields.")
    f <- as.data.frame(old$frame)
    for (v in c("SPECIES", "SITE_ID", "pop_id", "bms_id")) f[[v]] <- as.character(f[[v]])
    need(!anyNA(f[c("SPECIES", "SITE_ID", "YEAR", "pop_id", "bms_id")]) &&
           !anyDuplicated(f[c("pop_id", "YEAR")]), "Invalid source abundance identities.")
    need(same(f$ABUND_INDEX, old$data$y, 0) && same(f$year_decade, old$data$time, 0),
         "Source response/time and frame disagree.")
    need(identical(colnames(old$X), names(fit$beta)) &&
           same(fit$parameters$beta, unname(fit$beta), 0), "Old fixed coefficient mapping mismatch.")
    need(length(fit$parameters$ar_obs) == nrow(f), "Old AR1 states/rows disagree.")
    on <- scale_plasticity(onset); off <- scale_plasticity(offset)
    a <- match(key(f), key(on$populations)); b <- match(key(f), key(off$populations))
    keep <- which(!is.na(a) & !is.na(b))
    need(length(keep) > 0, "No abundance populations have both updated plasticities.")
    keep <- keep[order(f$pop_id[keep], f$YEAR[keep])]
    d <- f[keep, , drop = FALSE]
    need(identical(d$bms_id, as.character(on$populations$bms_id[a[keep]])) &&
           identical(d$bms_id, as.character(off$populations$bms_id[b[keep]])), "BMS identity mismatch.")
    d$onset_plasticity_z <- on$populations$z[a[keep]]
    d$offset_plasticity_bio_z <- off$populations$z[b[keep]]
    d$onset_plasticity_raw_days_per_C <- on$populations$plasticity_raw_days_per_C[a[keep]]
    d$offset_plasticity_raw_days_per_C <- off$populations$plasticity_raw_days_per_C[b[keep]]
    # Recover saved index order, never infer random-effect order alphabetically.
    index_map <- function(label, index, n) {
      m <- unique(data.frame(label = label, index = as.integer(index)))
      need(nrow(m) == n && !anyDuplicated(m$label) && !anyDuplicated(m$index) &&
             identical(sort(m$index), seq_len(n) - 1L), "Invalid saved latent index mapping.")
      m[order(m$index), , drop = FALSE]
    }
    sp <- index_map(f$SPECIES, old$data$species, nrow(fit$parameters$b_species))
    si <- index_map(f$SITE_ID, old$data$site, nrow(fit$parameters$b_site))
    sp <- sp[sp$label %in% d$SPECIES, , drop = FALSE]
    si <- si[si$label %in% d$SITE_ID, , drop = FALSE]
    need(identical(as.character(old$sites$SITE_ID)[old$data$site + 1L], f$SITE_ID),
         "Old spatial matrix/site table order disagrees with site indices.")
    sites <- old$sites[si$index + 1L, , drop = FALSE]
    # The old annual frame/population table already stores the reviewed coordinates.
    if (!all(c("latitude_deg", "longitude_deg") %in% names(sites))) {
      coords <- unique(f[c("SITE_ID", "latitude_deg", "longitude_deg")])
      need(!anyDuplicated(coords$SITE_ID), "Coordinates vary within site.")
      m <- match(as.character(sites$SITE_ID), coords$SITE_ID)
      sites$latitude_deg <- coords$latitude_deg[m]
      sites$longitude_deg <- coords$longitude_deg[m]
    }
    dec <- decompose(d, sites); d <- dec$frame; X <- dec$X
    need(identical(colnames(X), names(fit$beta)), "Updated fixed terms changed unexpectedly.")
    years <- as.integer(old$years)
    need(!length(setdiff(years, d$YEAR)), "An entire source year was lost; review the sample before fitting.")
    need(identical(as.integer(f$YEAR), years[old$data$year + 1L]), "Old field/year mapping mismatch.")
    n <- nrow(d); previous <- rep(-1L, n); gap <- rep(1L, n)
    if (n > 1L) {
      ids <- which(d$pop_id[-1L] == d$pop_id[-n]) + 1L
      previous[ids] <- ids - 2L
      gap[ids] <- as.integer(d$YEAR[ids] - d$YEAR[ids - 1L])
    }
    need(all(gap > 0L), "Non-positive calendar gap.")
    data <- old$data
    data$X <- unname(X); data$y <- old$data$y[keep]; data$time <- old$data$time[keep]
    data$species <- as.integer(match(d$SPECIES, sp$label) - 1L)
    data$site <- as.integer(match(d$SITE_ID, si$label) - 1L)
    data$year <- as.integer(match(d$YEAR, years) - 1L)
    data$previous <- previous; data$gap <- gap
    data$A <- old$data$A[si$index + 1L, , drop = FALSE]
    pars <- fit$parameters
    pars$b_species <- pars$b_species[sp$index + 1L, , drop = FALSE]
    pars$b_site <- pars$b_site[si$index + 1L, , drop = FALSE]
    pars$ar_obs <- pars$ar_obs[keep]
    # Projection preserves the old fixed predictor as closely as the new design allows.
    target_eta <- as.numeric(old$X[keep, , drop = FALSE] %*% unname(fit$beta))
    pars$beta <- as.numeric(qr.coef(qr(X), target_eta))
    need(all(is.finite(pars$beta)), "Cannot construct the updated fixed-effect warm start.")
    need(nrow(pars$field) == nrow(data$M0) && ncol(pars$field) == length(years),
         "Mesh/year dimensions changed.")
    audit <- unique(f[c("pop_id", "SPECIES", "SITE_ID", "bms_id")])
    aa <- match(key(audit), key(on$populations)); bb <- match(key(audit), key(off$populations))
    audit$has_new_onset <- !is.na(aa); audit$has_new_offset <- !is.na(bb)
    audit$retained <- audit$has_new_onset & audit$has_new_offset
    audit$n_abundance_years <- as.integer(table(f$pop_id)[audit$pop_id])
    scaling <- rbind(cbind(analysis = paste0("onset_", group), on$scaling),
                     cbind(analysis = paste0("offset_", group), off$scaling))
    list(frame = d, populations = dec$populations, X = X, data = data, parameters = pars,
         sites = sites, years = years, component_scales = cbind(group = group, dec$scales),
         scaling = scaling, sample_audit = audit, retained_source_rows = keep,
         warm_projection_RMSE = sqrt(mean((as.numeric(X %*% pars$beta) - target_eta)^2)),
         extracted_outside_abundance = data.frame(
           analysis = c(paste0("onset_", group), paste0("offset_", group)),
           n = c(sum(!key(on$populations) %in% key(audit)),
                 sum(!key(off$populations) %in% key(audit)))))
  }
  load_engine <- function(cfg) {
    e <- new.env(parent = globalenv())
    sys.source(cfg$engine, envir = e)
    engine <- environment(e$pheno_mundlak_qc$run)
    engine$VERSION <- VERSION
    engine
  }
  prepare_group <- function(cfg, group) {
    check_manifest(cfg)
    path <- file.path(cfg$out, paste0("prepared_mundlak_", group, ".rds"))
    done <- file.path(cfg$out, paste0("prepared_", group, "_completed.rds"))
    if (file.exists(done) && file.exists(path)) {
      need(identical(readRDS(done)$signature, cfg$signature), "Preparation cache mismatch.")
      return(invisible(path))
    }
    phase(cfg, "PREPARING_DEMOGRAPHY", "Updating plasticities, scaling and latent indices", group)
    old_file <- file.path(cfg$demography_run, paste0("prepared_mundlak_", group, ".rds"))
    fit_file <- file.path(cfg$demography_run, paste0(group, "__mundlak"), "fit.rds")
    old <- readRDS(old_file); fit <- readRDS(fit_file)
    read_pop <- function(metric) {
      x <- readRDS(file.path(cfg$out, "extraction", paste0(metric, "_", group), "extraction.rds"))
      need(identical(x$signature, cfg$signature), "Extraction/preparation signature mismatch.")
      x$populations
    }
    p <- rebuild_input(old, fit, read_pop("onset"), read_pop("offset"), group)
    p$version <- VERSION; p$group <- group
    p$signature <- digest::digest(list(cfg$signature, group, p$frame, p$data), algo = "sha256")
    p$source_bridge_signature <- old$signature
    p$source_data_signature <- old$signature
    p$source_fit_objective <- fit$optimizer$objective
    p$source_prepared_file <- old_file; p$source_fit_file <- fit_file
    p$predictors_updated <- TRUE
    p$plasticity_uncertainty_propagated <- FALSE
    p$residual_diagnostics <- "PENDING"
    csv(p$sample_audit, file.path(cfg$out, paste0("population_overlap_", group, ".csv")))
    csv(p$scaling, file.path(cfg$out, paste0("plasticity_scaling_", group, ".csv")))
    csv(p$component_scales, file.path(cfg$out, paste0("component_scales_", group, ".csv")))
    csv(p$extracted_outside_abundance, file.path(cfg$out, paste0("extracted_outside_abundance_", group, ".csv")))
    vars <- colnames(p$X)[3:8]
    eng <- load_engine(cfg)
    csv(eng$predictor_vif(p$populations, vars), file.path(cfg$out, paste0("predictor_VIF_", group, ".csv")))
    cor <- stats::cor(p$populations[vars])
    csv(data.frame(predictor = rownames(cor), cor, check.names = FALSE),
        file.path(cfg$out, paste0("predictor_correlations_", group, ".csv")))
    atomic(p, path)
    atomic(list(signature = cfg$signature, group = group, nobs = nrow(p$frame),
                npop = nrow(p$populations), removed = sum(!p$sample_audit$retained),
                warm_projection_RMSE = p$warm_projection_RMSE), done)
    invisible(path)
  }

  # Replay the OLD model on its OLD data to verify native-code compatibility.
  # The new fit is NOT required to reproduce the old likelihood on changed X.
  verify_old_likelihood <- function(cfg, group) {
    file <- file.path(cfg$out, paste0("old_likelihood_replay_", group, ".rds"))
    if (file.exists(file)) {
      z <- readRDS(file)
      need(identical(z$signature, cfg$signature) && isTRUE(z$ok), "Invalid old-likelihood replay cache.")
      return(invisible(TRUE))
    }
    eng <- load_engine(cfg); eng$load_dll(cfg$source_run)
    phase(cfg, "VERIFYING_OLD_LIKELIHOOD", "One evaluation on original data; no refit", group)
    p <- readRDS(file.path(cfg$demography_run, paste0("prepared_mundlak_", group, ".rds")))
    f <- readRDS(file.path(cfg$demography_run, paste0(group, "__mundlak"), "fit.rds"))
    need(identical(p$signature, f$signature), "Old fit/preparation signature mismatch.")
    obj <- eng$make_objective(p$data, f$parameters)
    on.exit(try(TMB::FreeADFun(obj), silent = TRUE), add = TRUE)
    value <- obj$fn(obj$par)
    report <- obj$report(obj$env$last.par)
    error <- abs(value - f$optimizer$objective)
    eta_error <- max(abs(as.numeric(report$eta) - f$eta))
    ok <- is.finite(error) && error < 0.001 && is.finite(eta_error) && eta_error < 1e-5
    z <- data.frame(group = group, objective_error = error, eta_error = eta_error, ok = ok)
    csv(z, sub("rds$", "csv", file))
    need(ok, "Old likelihood replay failed. Check original R/TMB libraries and compiled engine.")
    atomic(list(signature = cfg$signature, ok = ok, check = z), file)
    invisible(TRUE)
  }
  fit_group <- function(cfg, group) {
    check_manifest(cfg)
    eng <- load_engine(cfg)
    eng$fit_group(cfg, group)
  }
  diagnostic_tables <- function(p, fit) {
    need(identical(fit$signature, p$signature), "Residuals and input refer to different fits.")
    d <- p$frame
    mu <- exp(fit$eta); shape <- fit$checks$Gamma_shape
    need(length(mu) == nrow(d) && all(is.finite(mu) & mu > 0) &&
           length(shape) == 1L && is.finite(shape) && shape > 0, "Invalid Gamma fitted values.")
    pit <- stats::pgamma(d$ABUND_INDEX, shape = shape, scale = mu / shape)
    r <- stats::qnorm(pmin(1 - 1e-10, pmax(1e-10, pit)))
    pearson <- (d$ABUND_INDEX - mu) / (mu / sqrt(shape))
    distribution <- data.frame(n = length(r), mean = mean(r), SD = stats::sd(r),
      q025 = unname(stats::quantile(r, .025)), median = stats::median(r),
      q975 = unname(stats::quantile(r, .975)), mean_Pearson = mean(pearson),
      SD_Pearson = stats::sd(pearson), conditioning = "fitted random effects, AR1 and spatial field")
    bin <- pmin(10L, ceiling(rank(mu, ties.method = "first") * 10 / length(mu)))
    binned <- do.call(rbind, lapply(sort(unique(bin)), function(b) {
      i <- which(bin == b)
      data.frame(decile = b, n = length(i), mean_fitted = mean(mu[i]),
                 mean_PIT = mean(r[i]), SD_PIT = stats::sd(r[i]))
    }))
    cor_safe <- function(x, y) {
      if (length(x) < 3 || !isTRUE(stats::sd(x) > 0) || !isTRUE(stats::sd(y) > 0)) return(NA_real_)
      stats::cor(x, y)
    }
    rowkey <- paste(key(d), d$YEAR, sep = "\r")
    temporal <- do.call(rbind, lapply(1:5, function(lag) {
      j <- match(paste(key(d), d$YEAR - lag, sep = "\r"), rowkey)
      i <- which(!is.na(j))
      data.frame(lag_years = lag, n_pairs = length(i),
                 correlation = cor_safe(r[i], r[j[i]]))
    }))
    rho <- fit$checks$AR1_rho
    latent <- fit$parameters$ar_obs
    latent_sd <- exp(fit$parameters$log_sd[5L])
    need(length(latent) == nrow(d) && is.finite(rho) && abs(rho) < 1 &&
           is.finite(latent_sd) && latent_sd > 0, "Invalid AR1 parameters.")
    innov <- latent / latent_sd
    i <- which(p$data$previous >= 0L); prev <- p$data$previous[i] + 1L
    rg <- rho ^ p$data$gap[i]
    innov[i] <- (latent[i] - rg * latent[prev]) / (latent_sd * sqrt(1 - rg^2))
    innovations <- data.frame(n = length(innov), mean = mean(innov), SD = stats::sd(innov),
      successive_correlation = cor_safe(innov[i], innov[prev]),
      interpretation = "descriptive latent-mode innovations; not independent observed residuals")
    site_year <- aggregate(r, d[c("SITE_ID", "YEAR")], mean)
    names(site_year)[3] <- "residual"
    list(distribution = distribution, binned = binned, temporal = temporal,
         innovations = innovations, site_year = site_year)
  }
  diagnostics <- function(cfg, group) {
    out <- file.path(cfg$out, paste0(group, "__mundlak"), "residual_diagnostics")
    marker <- file.path(out, "status.rds")
    outputs <- file.path(out, c("distribution.csv", "binned.csv", "temporal.csv",
                                "innovations.csv", "annual_spatial_Moran.csv"))
    if (file.exists(marker) && all(file.exists(outputs))) {
      need(identical(readRDS(marker)$signature, cfg$signature), "Diagnostic cache mismatch.")
      return(invisible(TRUE))
    }
    p <- readRDS(file.path(cfg$out, paste0("prepared_mundlak_", group, ".rds")))
    f <- readRDS(file.path(cfg$out, paste0(group, "__mundlak"), "fit.rds"))
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    phase(cfg, "RESIDUAL_SCREENS", "Conditional residuals and two-sided spatial screen", group)
    z <- diagnostic_tables(p, f)
    for (nm in c("distribution", "binned", "temporal", "innovations"))
      csv(z[[nm]], file.path(out, paste0(nm, ".csv")))
    packages("spdep")
    set.seed(20261020L)
    moran <- lapply(sort(unique(z$site_year$YEAR)), function(year) {
      d <- z$site_year[z$site_year$YEAR == year, , drop = FALSE]
      m <- match(as.character(d$SITE_ID), as.character(p$sites$SITE_ID))
      need(!anyNA(m), "Missing spatial diagnostic coordinates.")
      d$x <- p$sites$x_km[m]; d$y <- p$sites$y_km[m]
      # Aggregate co-located transects so nearest-neighbour distances are unambiguous.
      d <- aggregate(residual ~ x + y, d, mean)
      if (nrow(d) < 4 || !isTRUE(stats::sd(d$residual) > 0)) return(data.frame(
        YEAR = year, n_locations = nrow(d), Moran_I = NA_real_, p_two_sided = NA_real_))
      nb <- spdep::knn2nb(spdep::knearneigh(as.matrix(d[c("x", "y")]),
                                          k = min(8L, nrow(d) - 1L)), sym = TRUE)
      mc <- spdep::moran.mc(d$residual, spdep::nb2listw(nb, style = "W"),
                            nsim = 499L, alternative = "two.sided")
      data.frame(YEAR = year, n_locations = nrow(d),
                 Moran_I = as.numeric(mc$statistic), p_two_sided = mc$p.value)
    })
    csv(do.call(rbind, moran), file.path(out, "annual_spatial_Moran.csv"))
    writeLines(c("Conditional descriptive screens; not calibrated goodness-of-fit tests.",
      "No IID uniformity/dispersion p-values are claimed.",
      "Latent AR1 innovations are based on estimated modes.",
      "New residual review remains PENDING. Numerical acceptance is separate.",
      "Plasticity-estimation uncertainty has not been propagated into demographic SEs."),
      file.path(out, "interpretation.txt"))
    atomic(list(signature = cfg$signature, numerical_checks_ok = f$checks$numerical_checks_ok,
                residual_review = "PENDING", status = "DESCRIPTIVE_SCREENS_EXPORTED"),
           file.path(out, "status.rds"))
    invisible(TRUE)
  }

  resources <- function(cfg) {
    require_environment_idle(cfg$environment_run)
    selected <- select_environment(cfg$environment_run)
    need(all(selected$ready) && identical(selected, cfg$environment),
         "Environmental source selection changed after preflight. Restart run().")
    packages("ps")
    free <- as.numeric(ps::ps_system_memory()$avail) / 1024^3
    disk <- as.numeric(ps::ps_disk_usage(cfg$out)$available) / 1024^3
    need(is.finite(free) && free >= cfg$minimum_free_ram_GiB,
         paste("Insufficient available RAM:", round(free, 1), "GiB. Resume when memory is available."))
    need(is.finite(disk) && disk >= cfg$min_disk,
         paste("Insufficient disk:", round(disk, 1), "GiB."))
  }
  child_stage <- function(cfg, stage, value) {
    resources(cfg); check_manifest(cfg)
    # Separate processes release the complete native model/tape between stages.
    callr::r(function(script, cfg, stage, value) {
      e <- new.env(parent = globalenv()); sys.source(script, envir = e)
      env <- environment(e$pheno_post_env$run)
      get(stage, envir = env, inherits = FALSE)(cfg, value)
    }, args = list(script = cfg$script, cfg = cfg, stage = stage, value = value),
    libpath = cfg$libpaths, user_profile = FALSE, system_profile = FALSE, show = TRUE)
  }
  controller <- function(cfg) {
    handle <- ps::ps_handle()
    owner <- list(pid = Sys.getpid(), create_time = ps::ps_create_time(handle), token = cfg$token)
    atomic(owner, file.path(cfg$lock, "owner.rds"))
    atomic(c(owner, list(out = cfg$out, log = cfg$log)), file.path(cfg$out, "worker_owner.rds"))
    on.exit({
      path <- file.path(cfg$lock, "owner.rds")
      if (file.exists(path) && identical(readRDS(path)$token, cfg$token)) unlink(cfg$lock, recursive = TRUE)
    }, add = TRUE)
    tryCatch({
      check_manifest(cfg)
      for (a in ANALYSES) child_stage(cfg, "extract_one", a)
      statuses <- list()
      for (g in GROUPS) {
        child_stage(cfg, "prepare_group", g)
        status <- "PREPARED"
        if (cfg$mode == "fit") {
          child_stage(cfg, "verify_old_likelihood", g)
          status <- child_stage(cfg, "fit_group", g)
          child_stage(cfg, "diagnostics", g)
        }
        statuses[[g]] <- data.frame(group = g, status = status, residual_review = "PENDING")
        csv(do.call(rbind, statuses), file.path(cfg$out, "workflow_status.csv"))
      }
      tab <- do.call(rbind, statuses)
      eng <- load_engine(cfg); eng$combine_outputs(cfg, tab)
      final <- if (cfg$mode == "prepare") "PREPARED" else if (
        all(tab$status %in% c("FIT_OK", "FIT_OK_CACHED"))) "FITTED_RESIDUAL_REVIEW_PENDING" else
          "NUMERICAL_REVIEW_REQUIRED"
      atomic(list(signature = cfg$signature, status = final, groups = tab, time = Sys.time()),
             file.path(cfg$out, paste0("result_", cfg$mode, ".rds")))
      phase(cfg, final, "Finished; inspect current diagnostics before final inference", terminal = TRUE)
    }, error = function(e) {
      msg <- conditionMessage(e)
      atomic(list(signature = cfg$signature, error = msg, time = Sys.time()),
             file.path(cfg$out, "last_error.rds"))
      phase(cfg, "STOPPED_WITH_ERROR", msg, terminal = TRUE)
      stop(e)
    })
  }

  run <- function(project_root = ROOT, environment_run = default_env(project_root),
                  demography_run = NULL, mode = c("fit", "prepare"),
                  minimum_free_ram_GiB = 20, minimum_free_disk_GiB = 40,
                  optimizer_iterations = 1200L, optimizer_evaluations = 1800L,
                  background = TRUE) {
    mode <- match.arg(mode)
    packages(c("callr", "ps", "digest", "sdmTMB", "TMB", "Matrix"))
    if (mode == "fit") packages("spdep")
    pf <- preflight(project_root, environment_run, demography_run)
    need(pf$ready, "Sources are not ready. No extraction or demographic fit started.")
    project_root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
    script <- file.path(project_root, "post_environment_mundlak.R")
    engine <- file.path(project_root, ENGINE)
    native <- TMB::dynlib(file.path(pf$demography$likelihood_run, "compiled", "trend_spde"))
    paths <- c(unlist(pf$environment[c("fit", "prepared", "checks", "marker", "original_marker")]),
               pf$demography$paths, native)
    manifest <- file_manifest(paths)
    need(all(file.exists(c(script, engine))), "Required pipeline scripts are missing from project_root.")
    code_md5 <- unname(tools::md5sum(c(script, engine)))
    libraries <- vapply(c("sdmTMB", "TMB", "Matrix"), function(p) as.character(utils::packageVersion(p)), "")
    signature <- digest::digest(list(VERSION, manifest, code_md5, libraries,
                                     optimizer_iterations, optimizer_evaluations), algo = "sha256")
    base <- file.path(project_root, "output", "population_trends", "post_environment_mundlak_15km")
    out <- file.path(base, paste0("run_", substr(signature, 1, 12)))
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    out <- normalizePath(out, winslash = "/", mustWork = TRUE)
    completed <- file.path(out, paste0("result_", mode, ".rds"))
    if (file.exists(completed)) {
      z <- readRDS(completed)
      need(identical(z$signature, signature), "Completed run signature mismatch.")
      message("Already completed: ", z$status, "\n", out)
      job <- list(process = NULL, out = out, log = NULL)
      state$job <- job; return(invisible(job))
    }
    lock <- file.path(base, "pipeline.lock")
    if (dir.exists(lock)) {
      owner_path <- file.path(lock, "owner.rds")
      need(file.exists(owner_path), paste("Lock has no owner; inspect before removing:", lock))
      need(!alive(readRDS(owner_path)), "A post-environment workflow is already active.")
      unlink(lock, recursive = TRUE)
    }
    need(dir.create(lock), "Cannot acquire post-environment workflow lock.")
    token <- paste(Sys.getpid(), format(Sys.time(), "%Y%m%d%H%M%OS6"), sep = "_")
    atomic(list(pid = Sys.getpid(), create_time = ps::ps_create_time(ps::ps_handle()), token = token),
           file.path(lock, "owner.rds"))
    launched <- FALSE
    on.exit(if (!launched) unlink(lock, recursive = TRUE), add = TRUE)
    snapshot <- file.path(out, "code")
    dir.create(snapshot, showWarnings = FALSE)
    need(all(file.copy(c(script, engine), snapshot, overwrite = TRUE)), "Cannot snapshot pipeline code.")
    cfg <- list(version = VERSION, signature = signature, project_root = project_root,
      environment_run = normalizePath(environment_run, winslash = "/", mustWork = TRUE),
      environment = pf$environment, demography_run = pf$demography$run,
      source_run = pf$demography$likelihood_run, manifest = manifest, out = out,
      progress = file.path(out, "progress.rds"), lock = lock, token = token, mode = mode,
      script = file.path(snapshot, basename(script)), engine = file.path(snapshot, ENGINE),
      minimum_free_ram_GiB = minimum_free_ram_GiB, min_disk = minimum_free_disk_GiB,
      iter = as.integer(optimizer_iterations), eval = as.integer(optimizer_evaluations),
      libpaths = .libPaths(), libraries = libraries,
      log = file.path(out, paste0("controller_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log")))
    resources(cfg)
    atomic(cfg, file.path(out, "run_config.rds"))
    csv(manifest, file.path(out, "source_manifest.csv"))
    csv(pf$environment, file.path(out, "environmental_sources.csv"))
    writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
    phase(cfg, "QUEUED", paste("Mode:", mode))
    atomic(list(out = out, log = cfg$log), file.path(base, "latest_job.rds"))
    if (!background) {
      controller(cfg)
      job <- list(process = NULL, out = out, log = cfg$log)
      state$job <- job
      return(invisible(job))
    }
    child <- callr::r_bg(function(script, cfg) {
      e <- new.env(parent = globalenv()); sys.source(script, envir = e)
      environment(e$pheno_post_env$run)$controller(cfg)
    }, args = list(script = cfg$script, cfg = cfg), libpath = cfg$libpaths,
    stdout = cfg$log, stderr = "2>&1", supervise = FALSE,
    user_profile = FALSE, system_profile = FALSE)
    launched <- TRUE
    job <- list(process = child, out = out, log = cfg$log)
    state$job <- job
    message("Post-environment pipeline started. Keep the computer on.\n",
            "Monitor: pheno_post_env$status()\n", out)
    invisible(job)
  }
  status <- function(job = NULL, project_root = ROOT) {
    if (is.null(job)) job <- state$job
    if (is.null(job)) {
      f <- file.path(project_root, "output", "population_trends",
                     "post_environment_mundlak_15km", "latest_job.rds")
      need(file.exists(f), "No post-environment job has been started.")
      job <- readRDS(f)
    }
    p <- file.path(job$out, "progress.rds")
    if (file.exists(p)) print(readRDS(p))
    message("Alive: ", active_run(job$out))
    s <- file.path(job$out, "workflow_status.csv")
    if (file.exists(s)) print(utils::read.csv(s), row.names = FALSE)
    invisible(job)
  }
  list(preflight = preflight, run = run, status = status)
})
