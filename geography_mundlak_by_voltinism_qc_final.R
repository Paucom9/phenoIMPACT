# phenoIMPACT | Direct geographic patterns of phenological plasticity
# 2026-10-20 | 1.1.0
#
# source("geography_mundlak_by_voltinism_qc_final.R", encoding = "UTF-8")
# pheno_geo_mundlak$preflight()  # small metadata only; safe during other fits
# geo_job <- pheno_geo_mundlak$run()  # AFTER all four selected sources pass QC and workers stop
# pheno_geo_mundlak$status()
#
# Four Gaussian ML models: onset/offset x uni/multivoltine.
# phenology ~ anomaly * (latitude_within + latitude_between +
#                       longitude_within + longitude_between)
#             [+ annual ONSET_mean_z for offsets] + original random terms
# Annual IID Matern field, unchanged 15-km mesh, spatial = "off".
# The response is ANNUAL OBSERVED PHENOLOGY, not extracted plasticity.
# Exact rows, anomaly units/windows and mesh come from the accepted environmental
# prepared inputs. No environmental coefficient or predicted slope is used.
# Coordinates: EPSG:3035 metres -> EPSG:4326 degrees, effects per 10 degrees.
# Species means give equal weight to each species x site population, within group.
# All four main geographic effects accompany their anomaly interactions.
#
# Sourcing only defines this new API. It never starts a process or modifies the
# environmental workflow. One detached controller runs one fit process at a time.
# Completed fits are reused; a saved fit can be re-exported without refitting.
# An interrupted optimization without a saved fit restarts from its initial values.
# Numerical FIT_OK is not a substitute for residual/model diagnostics.
# Preflight recognizes v1/v2 refinements and isolated workers using small metadata.
# See README_geography_checker_CA.md. Loading this script never starts a fit.

pheno_geo_mundlak <- local({
  VERSION <- "geography_mundlak_by_voltinism_qc_final_15km_1.1.0"
  SOURCE_VERSION <- "environment_mundlak_by_voltinism_qc_final_15km_2.1_resilient"
  ROOT <- "E:/phenoIMPACT project/code/phenoIMPACT"
  ENV_RELATIVE <- file.path("output", "phenology_plasticity",
    "environment_mundlak_by_voltinism_qc_final_15km", "run_de0890c85d4a")
  OUTPUT_NAME <- "geography_mundlak_by_voltinism_qc_final_15km"
  ANALYSES <- c("offset_multivoltine", "onset_multivoltine",
                "offset_univoltine", "onset_univoltine")
  V2_ANALYSES <- c("offset_univoltine", "onset_multivoltine", "onset_univoltine")
  GEO_VARS <- c("latitude_within_species_10deg", "latitude_species_mean_10deg",
                "longitude_within_species_10deg", "longitude_species_mean_10deg")
  state <- new.env(parent = emptyenv())

  need <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)
  csv <- function(x, path) utils::write.csv(x, path, row.names = FALSE, na = "")
  text_formula <- function(x) paste(deparse(x, width.cutoff = 500L), collapse = " ")

  atomic_save <- function(x, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    tmp <- tempfile(".writing_", tmpdir = dirname(path), fileext = ".rds")
    on.exit(unlink(tmp), add = TRUE)
    saveRDS(x, tmp, compress = "gzip")
    backup <- paste0(path, ".previous")
    if (file.exists(backup)) unlink(backup)
    had <- file.exists(path)
    if (had) need(file.rename(path, backup), paste("Cannot retain", path))
    if (!file.rename(tmp, path)) {
      if (had) file.rename(backup, path)
      stop("Cannot finalize: ", path, call. = FALSE)
    }
    if (file.exists(backup)) unlink(backup)
    invisible(path)
  }

  read_small <- function(path) {
    need(file.exists(path), paste("Missing metadata:", path))
    need(file.info(path)$size < 5 * 1024^2, paste("Unexpectedly large metadata:", path))
    readRDS(path)
  }

  phase <- function(cfg, stage, analysis = "", detail = "", terminal = FALSE) {
    p <- list(stage = stage, analysis = analysis, detail = detail,
              time = Sys.time(), terminal = terminal)
    atomic_save(p, file.path(cfg$out, "progress.rds"))
    message(format(p$time, "%Y-%m-%d %H:%M:%S"), " | ", stage,
            " | ", analysis, " | ", detail)
    invisible(p)
  }

  packages <- function() c("sdmTMB", "TMB", "Matrix", "fmesher", "lme4",
                            "sf", "digest", "callr", "ps")
  check_packages <- function() {
    missing <- packages()[!vapply(packages(), requireNamespace, logical(1), quietly = TRUE)]
    need(!length(missing), paste("Missing packages:", paste(missing, collapse = ", ")))
    need(as.character(utils::packageVersion("sdmTMB")) == "1.1.0",
         "Use the same sdmTMB 1.1.0 library as the accepted sources.")
  }

  expected_spec <- function(analysis) {
    need(analysis %in% ANALYSES, "Unknown analysis.")
    onset <- startsWith(analysis, "onset_")
    list(response = if (onset) "ONSET_mean" else "OFFSET_mean",
         anomaly = if (onset) "clim_anomaly_tw60" else "clim_anomaly_tw90",
         window = if (onset) 60L else 90L,
         control = if (onset) character() else "ONSET_mean_z",
         voltinism = sub("^(onset|offset)_", "", analysis),
         orientation = if (analysis == "offset_multivoltine") 1 else -1,
         orientation_label = if (analysis == "offset_multivoltine")
           "delay conditional on annual onset" else "advancement")
  }

  plan <- function() {
    do.call(rbind, lapply(ANALYSES, function(a) {
      s <- expected_spec(a)
      data.frame(analysis = a, response = s$response, anomaly = s$anomaly,
        window_days = s$window, additive_control = paste(s$control, collapse = "+"),
        n_geographic_interactions = 4L, stringsAsFactors = FALSE)
    }))
  }

  pid_alive <- function(pid, create_time) {
    need(is.numeric(pid) && length(pid) == 1L && is.finite(pid) &&
      pid > 0 && pid == floor(pid) && length(create_time) == 1L &&
      is.finite(as.numeric(create_time)), "Invalid worker ownership metadata.")
    tryCatch({
      h <- ps::ps_handle(pid)
      isTRUE(ps::ps_is_running(h)) && isTRUE(abs(
        as.numeric(ps::ps_create_time(h)) - as.numeric(create_time)) < 0.001)
    }, error = function(e) {
      if (inherits(e, "no_such_process")) return(FALSE)
      stop("Cannot verify worker process: ", conditionMessage(e), call. = FALSE)
    })
  }

  active_owner <- function(path) {
    owner <- read_small(path)
    need(is.list(owner), paste("Invalid worker ownership metadata:", path))
    pid_alive(owner$pid, owner$create_time)
  }

  refinement_roots <- function(environment_run) {
    v1 <- file.path(environment_run, "postfit_refinement_queue_v1")
    v2 <- file.path(environment_run, "postfit_refinement_queue_v2_isolated")
    config_files <- character()
    check_source <- function(path) {
      cfg <- read_small(path)
      need(is.character(cfg$source_run) && length(cfg$source_run) == 1L &&
        !is.na(cfg$source_run) && identical(
          normalizePath(cfg$source_run, winslash = "/", mustWork = TRUE),
          normalizePath(environment_run, winslash = "/", mustWork = TRUE)),
        "Refinement configuration belongs to a different environmental run.")
      cfg
    }
    v2_config <- file.path(v2, "refinement_config.rds")
    if (file.exists(v2_config)) {
      cfg <- check_source(v2_config)
      need(is.character(cfg$old_refinement) && length(cfg$old_refinement) == 1L &&
        !is.na(cfg$old_refinement) && nzchar(cfg$old_refinement),
        "V2 refinement configuration lacks its v1 source.")
      v1 <- cfg$old_refinement
      config_files <- c(config_files, v2_config)
    }
    v1_config <- file.path(v1, "refinement_config.rds")
    if (file.exists(v1_config)) {
      check_source(v1_config)
      config_files <- c(config_files, v1_config)
    }
    list(v1 = v1, v2 = v2, config_files = config_files)
  }

  owner_files <- function(roots) {
    unique(unlist(lapply(roots[dir.exists(roots)], function(root)
      list.files(root, pattern = "^(worker|fit|isolated_worker)_owner\\.rds$",
        recursive = TRUE, full.names = TRUE)), use.names = FALSE))
  }

  environmental_owner_files <- function(environment_run) {
    roots <- refinement_roots(environment_run)
    owner_files(c(environment_run, roots$v1, roots$v2))
  }

  require_environment_idle <- function(environment_run) {
    paths <- environmental_owner_files(environment_run)
    need(!any(vapply(paths, active_owner, logical(1))),
      "An environmental/refinement worker is active (original, v1 or v2); no new fit launched.")
  }

  active_jobs <- function(project_root, environment_run = NULL) {
    need(requireNamespace("ps", quietly = TRUE), "Install ps for process checks.")
    found <- character()
    for (key in grep("pheno.*(job|active)", names(options()), value = TRUE)) {
      job <- getOption(key)
      alive <- tryCatch(!is.null(job$process) && job$process$is_alive(),
        error = function(e) stop("Cannot verify workflow option ", key, ": ",
          conditionMessage(e), call. = FALSE))
      if (isTRUE(alive)) found <- c(found, paste("R option", key))
    }
    files <- owner_files(file.path(project_root, "output",
      c("phenology_plasticity", "population_trends")))
    if (!is.null(environment_run))
      files <- unique(c(files, environmental_owner_files(environment_run)))
    for (f in files) {
      # Corrupt/inaccessible ownership is not evidence that a worker stopped.
      if (active_owner(f)) found <- c(found, f)
    }
    unique(found)
  }

  accepted_checks <- function(path) {
    need(file.info(path)$size < 5 * 1024^2, paste("Unexpectedly large checks:", path))
    z <- utils::read.csv(path, stringsAsFactors = FALSE)
    is.data.frame(z) && nrow(z) == 1L &&
      isTRUE(z$numerical_checks_ok) && isTRUE(z$convergence_code == 0L) &&
      isTRUE(z$pdHess) && isTRUE(is.finite(z$max_abs_gradient)) &&
      isTRUE(z$max_abs_gradient < 0.001)
  }

  source_status <- function(environment_run, require_complete = FALSE) {
    cfg <- read_small(file.path(environment_run, "run_config.rds"))
    need(identical(cfg$version, SOURCE_VERSION), "Unsupported environmental source version.")
    need(is.character(cfg$signature) && length(cfg$signature) == 1L &&
           !is.na(cfg$signature) && nzchar(cfg$signature),
         "Environmental source has no valid signature.")
    roots <- refinement_roots(environment_run)
    v2_exists <- dir.exists(roots$v2)
    tab <- do.call(rbind, lapply(ANALYSES, function(a) {
      prep <- file.path(environment_run, a, "prepared_input.rds")
      completion <- file.path(environment_run, a, "completed.rds")
      done <- if (file.exists(completion)) read_small(completion) else NULL
      if (!is.null(done)) {
        need(identical(done$signature, cfg$signature), paste("Source completion signature mismatch:", a))
        need(is.character(done$status) && length(done$status) == 1L && !is.na(done$status),
             paste("Invalid source completion status:", a))
      }
      ref1 <- file.path(roots$v1, a)
      ref2 <- file.path(roots$v2, a)
      # Scheduled v2 models supersede v1 before the per-model directory exists.
      # Offset multivoltine is carried forward from the v1 path in the v2 config.
      refinement_version <- if (dir.exists(ref2) || (v2_exists && a %in% V2_ANALYSES)) {
        "v2"
      } else if (dir.exists(ref1) || (v2_exists && a == "offset_multivoltine")) {
        "v1"
      } else "none"
      status <- if (is.null(done)) "NOT_COMPLETED" else done$status
      kind <- "original"
      fit <- file.path(environment_run, a, "fit_mundlak.rds")
      checks <- file.path(environment_run, a, "fit_checks.csv")
      marker <- completion
      if (refinement_version != "none") {
        kind <- "refined"
        ref <- if (refinement_version == "v2") ref2 else ref1
        fit <- file.path(ref, "fit_refined_final.rds")
        checks <- file.path(ref, "fit_checks_refined.csv")
        marker <- file.path(ref, "refinement_completed.rds")
        status <- if (dir.exists(ref)) "REFINEMENT_NOT_COMPLETED" else "REFINEMENT_NOT_STARTED"
        if (file.exists(marker)) {
          z <- read_small(marker)
          need(is.data.frame(z$summary) && nrow(z$summary) == 1L &&
            identical(z$summary$analysis, a), paste("Refinement marker identity mismatch:", a))
          status <- z$summary$refinement_status
          need(is.character(status) && length(status) == 1L && !is.na(status),
            paste("Invalid refinement status:", a))
        }
      }
      all_files <- all(file.exists(c(completion, prep, fit, checks, marker)))
      ready <- all_files &&
        status %in% c("FIT_OK", "FIT_OK_REFINED", "FIT_OK_AFTER_REBUILD") &&
        accepted_checks(checks)
      data.frame(analysis = a, kind = kind, refinement_version = refinement_version,
        status = status, ready = ready, prepared_available = file.exists(prep),
        prepared_MiB = if (file.exists(prep)) round(file.info(prep)$size / 1024^2, 1) else NA_real_,
        prepared = prep, fit = fit, checks = checks, marker = marker,
        original_marker = completion,
        stringsAsFactors = FALSE)
    }))
    if (require_complete) need(all(tab$ready),
      "Wait until all four selected environmental sources pass numerical checks with prepared inputs; no fit was launched.")
    files <- unique(c(file.path(environment_run, "run_config.rds"),
      roots$config_files, tab$prepared, tab$original_marker, tab$marker, tab$checks))
    list(config = cfg, table = tab, files = files)
  }

  file_manifest <- function(paths) {
    paths <- unique(unname(normalizePath(paths, winslash = "/", mustWork = TRUE)))
    z <- file.info(paths)
    need(!anyNA(z$size) && all(!z$isdir), "Source file metadata unavailable.")
    data.frame(path = paths, bytes = z$size, modified = as.numeric(z$mtime),
      stringsAsFactors = FALSE)
  }

  require_sources_current <- function(cfg) {
    require_environment_idle(cfg$environment_run)
    src <- source_status(cfg$environment_run, require_complete = TRUE)
    need(identical(src$table, cfg$source_selection) && identical(src$files, cfg$source_files),
      "Environmental source selection changed after preflight; restart run().")
    need(identical(file_manifest(c(src$files, src$table$fit)), cfg$source_manifest),
      "An environmental source changed after preflight; restart run().")
    invisible(TRUE)
  }

  available_resources <- function(project_root) {
    mem <- ps::ps_system_memory()
    disk <- ps::ps_disk_usage(project_root)
    list(available_RAM_GiB = as.numeric(mem$avail) / 1024^3,
      available_disk_GiB = as.numeric(disk$available[1]) / 1024^3)
  }

  preflight <- function(project_root = ROOT,
                        environment_run = file.path(project_root, ENV_RELATIVE)) {
    project_root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
    environment_run <- normalizePath(environment_run, winslash = "/", mustWork = TRUE)
    src <- source_status(environment_run)
    jobs <- active_jobs(project_root, environment_run)
    out <- c(list(source = src$table, active_jobs = jobs,
      sources_ready = all(src$table$ready)), available_resources(project_root))
    out$resources_ready <- isTRUE(out$available_RAM_GiB >= 35) &&
      isTRUE(out$available_disk_GiB >= 50)
    out$ready <- out$sources_ready && !length(jobs) && out$resources_ready
    print(src$table[c("analysis", "kind", "refinement_version", "status", "ready",
      "prepared_available", "prepared_MiB")], row.names = FALSE)
    message("Available RAM: ", round(out$available_RAM_GiB, 1), " GiB; disk: ",
            round(out$available_disk_GiB, 1), " GiB.")
    if (length(jobs)) message("Active workflow(s):\n", paste(jobs, collapse = "\n"))
    message("Sources ready: ", out$sources_ready, "; ready with default resource limits: ", out$ready)
    message("Metadata only: no prepared data or fitted models loaded; no process started.")
    invisible(out)
  }

  resource_check <- function(path, minimum_ram_GiB, minimum_disk_GiB) {
    invisible(gc())
    mem <- ps::ps_system_memory()
    free <- as.numeric(mem$avail) / 1024^3
    disk <- as.numeric(ps::ps_disk_usage(path)$available[1]) / 1024^3
    message("Available RAM: ", round(free, 1), " GiB; disk: ", round(disk, 1), " GiB.")
    need(is.finite(free) && free >= minimum_ram_GiB,
         paste("Insufficient RAM; require", minimum_ram_GiB, "GiB available before preparing a fit."))
    need(is.finite(disk) && disk >= minimum_disk_GiB,
         paste("Insufficient disk; require", minimum_disk_GiB, "GiB available."))
    invisible(TRUE)
  }

  random_terms_text <- function(formula) {
    bars <- lme4::findbars(formula)
    need(length(bars) > 0L, "Source formula has no random effects.")
    paste(vapply(bars, function(x) paste0("(", text_formula(x), ")"), character(1)),
          collapse = " + ")
  }

  make_formula <- function(source_formula, s) {
    fixed <- paste0(s$anomaly, " * (", paste(GEO_VARS, collapse = " + "), ")")
    if (length(s$control)) fixed <- paste(fixed, "+", paste(s$control, collapse = " + "))
    stats::as.formula(paste(s$response, "~", fixed, "+", random_terms_text(source_formula)),
                      env = baseenv())
  }

  # Separate pure decomposition from the CRS transformation for exact unit tests.
  decompose_coordinates <- function(data, sites) {
    required <- c("SPECIES", "SITE_ID")
    need(all(required %in% names(data)), "Data lack SPECIES/SITE_ID.")
    need(all(c("SITE_ID", "latitude_deg", "longitude_deg") %in% names(sites)),
         "Site table lacks geographic coordinates.")
    need(!anyNA(data[required]) &&
           all(nzchar(as.character(data$SPECIES))) && all(nzchar(as.character(data$SITE_ID))),
         "Missing/empty population identifiers.")
    need(!anyDuplicated(as.character(sites$SITE_ID)), "Duplicate site coordinates.")
    need(all(is.finite(sites$latitude_deg)) && all(abs(sites$latitude_deg) <= 90) &&
           all(is.finite(sites$longitude_deg)) && all(abs(sites$longitude_deg) <= 180),
         "Invalid latitude/longitude.")
    pop <- unique(data[required])
    pop$SPECIES <- as.character(pop$SPECIES)
    pop$SITE_ID <- as.character(pop$SITE_ID)
    m <- match(pop$SITE_ID, as.character(sites$SITE_ID))
    need(!anyNA(m), "Missing coordinates: no rows may be dropped.")
    pop$latitude_deg <- sites$latitude_deg[m]
    pop$longitude_deg <- sites$longitude_deg[m]
    scaling <- list()
    for (axis in c("latitude", "longitude")) {
      raw <- pop[[paste0(axis, "_deg")]]
      mean_species <- tapply(raw, pop$SPECIES, mean)
      # Centre between-species means at the equal-species grand mean.
      # This only sets the intercept reference, not either geographic slope.
      centre <- mean(mean_species)
      mean_pop <- as.numeric(mean_species[pop$SPECIES])
      within <- paste0(axis, "_within_species_10deg")
      between <- paste0(axis, "_species_mean_10deg")
      pop[[paste0(axis, "_species_mean_deg")]] <- mean_pop
      pop[[within]] <- (raw - mean_pop) / 10
      pop[[between]] <- (mean_pop - centre) / 10
      need(max(abs(10 * (pop[[within]] + pop[[between]]) + centre - raw)) < 1e-10,
           paste("Coordinate reconstruction failed:", axis))
      need(max(abs(tapply(pop[[within]], pop$SPECIES, mean))) < 1e-10,
           paste("Species centering failed:", axis))
      scaling[[axis]] <- data.frame(axis = axis, centre_degrees = centre,
        divisor_degrees = 10, within_SD = stats::sd(pop[[within]]),
        between_population_weighted_SD = stats::sd(pop[[between]]),
        species_mean_weighting = "equal population within each voltinism group",
        between_centre_weighting = "equal species", stringsAsFactors = FALSE)
    }
    key <- function(d) paste(as.character(d$SPECIES), as.character(d$SITE_ID), sep = "\r")
    need(!any(grepl("\r", c(pop$SPECIES, pop$SITE_ID), fixed = TRUE)), "Invalid key delimiter in IDs.")
    m <- match(key(data), key(pop))
    need(!anyNA(m) && !anyDuplicated(key(pop)), "Population join failed.")
    d <- data
    for (v in c("latitude_deg", "longitude_deg", GEO_VARS)) d[[v]] <- pop[[v]][m]
    list(data = d, populations = pop, scaling = do.call(rbind, scaling))
  }

  recover_coordinates <- function(data) {
    vars <- c("SITE_ID", "x_3035", "y_3035", "x_km", "y_km")
    need(all(vars %in% names(data)), "Source data must include verified EPSG:3035 metre and km coordinates.")
    need(all(vapply(data[vars[-1]], function(x) is.numeric(x) && all(is.finite(x)), logical(1))),
         "Non-numeric or missing projected coordinates.")
    need(max(abs(data$x_3035 / 1000 - data$x_km)) < 1e-8 &&
           max(abs(data$y_3035 / 1000 - data$y_km)) < 1e-8,
         "Projected metre/km coordinates disagree.")
    sites <- unique(data[vars])
    need(!anyDuplicated(as.character(sites$SITE_ID)), "A site has conflicting projected coordinates.")
    if ("bms_id" %in% names(data)) {
      ids <- unique(data[c("SITE_ID", "bms_id")])
      need(!anyNA(ids) && !anyDuplicated(as.character(ids$SITE_ID)),
           "SITE_ID is not unique across monitoring networks.")
    }
    need(max(abs(sites$x_3035)) > 180 && max(abs(sites$y_3035)) > 90,
         "Coordinates resemble degrees, not EPSG:3035 metres.")
    sf_sites <- sf::st_as_sf(sites, coords = c("x_3035", "y_3035"), crs = 3035)
    ll_sf <- sf::st_transform(sf_sites, 4326)
    ll <- sf::st_coordinates(ll_sf)
    back <- sf::st_coordinates(sf::st_transform(ll_sf, 3035))
    need(max(abs(back - as.matrix(sites[c("x_3035", "y_3035")]))) < .01,
         "Coordinate transformation round-trip failed.")
    sites$longitude_deg <- ll[, 1]
    sites$latitude_deg <- ll[, 2]
    sites
  }

  check_source_input <- function(p, analysis) {
    s <- expected_spec(analysis)
    need(identical(p$version, SOURCE_VERSION) && identical(p$analysis, analysis),
         "Prepared source version/analysis mismatch.")
    for (nm in c("response", "anomaly", "window", "control", "voltinism"))
      need(isTRUE(all.equal(p$spec[[nm]], s[[nm]])), paste("Source specification mismatch:", nm))
    d <- p$data
    need(is.data.frame(d) && nrow(d) > 0, "Empty source data.")
    required <- unique(c("SPECIES", "SITE_ID", "YEAR", s$response, s$anomaly,
                          s$control, all.vars(p$formula)))
    need(all(required %in% names(d)), "Missing source model variables.")
    need(!anyNA(d[required]), "Missing values in accepted model rows; no filtering allowed.")
    need(!anyDuplicated(d[c("SPECIES", "SITE_ID", "YEAR")]), "Duplicated population-year rows.")
    need(is.numeric(d$YEAR) && all(is.finite(d$YEAR)) && all(d$YEAR == floor(d$YEAR)),
         "YEAR must retain integer calendar years.")
    for (v in c(s$response, s$anomaly, s$control))
      need(is.numeric(d[[v]]) && all(is.finite(d[[v]])), paste("Invalid source variable:", v))
    need(identical(as.character(p$formula[[2]]), s$response), "Source formula response mismatch.")
    need(!is.null(p$mesh$loc_xy) && all(c("x_km", "y_km") %in% names(d)) &&
      isTRUE(all.equal(unname(as.matrix(p$mesh$loc_xy)),
        unname(as.matrix(d[c("x_km", "y_km")])), check.attributes = FALSE, tolerance = 1e-10)),
      "Source mesh projection does not match the exact row order.")
    invisible(s)
  }

  design_matrix <- function(formula, data) {
    fixed <- lme4::nobars(formula)
    mf <- stats::model.frame(fixed, data, na.action = stats::na.fail)
    X <- stats::model.matrix(fixed, mf)
    need(nrow(X) == nrow(data) && all(is.finite(X)), "Design changed rows or contains non-finite values.")
    need(qr(X)$rank == ncol(X), "Geographic fixed-effect design is rank deficient; review support.")
    X
  }

  starting_parameters <- function(source, formula, data, beta) {
    TMB::openmp(n = 1L, DLL = "sdmTMB")
    skeleton <- sdmTMB::sdmTMB(formula, data = data, mesh = source$mesh, time = "YEAR",
      family = stats::gaussian(), spatial = "off", spatiotemporal = "iid", reml = FALSE,
      do_fit = FALSE, silent = TRUE, control = sdmTMB::sdmTMBcontrol(multiphase = FALSE,
        newton_loops = 0L, parallel = 1L, get_joint_precision = FALSE,
        collapse_spatial_variance = FALSE))
    base <- if (is.list(skeleton$tmb_params)) skeleton$tmb_params else
      skeleton$tmb_obj$env$parList(par = skeleton$tmb_obj$par)
    need(length(base$b_j) == length(beta), "Incompatible fixed-effect start dimensions.")
    # Only non-fixed starting parameters are reused. They are starting values,
    # not fixed estimates; the geographic likelihood re-estimates every parameter.
    for (nm in c("ln_phi", "ln_kappa", "ln_tau_E", "re_cov_pars", "re_b_pars", "epsilon_st")) {
      old <- source$start[[nm]]
      if (nm %in% names(base) && is.numeric(old) &&
          identical(dim(old), dim(base[[nm]])) && length(old) == length(base[[nm]]) &&
          all(is.finite(old))) base[[nm]] <- old
    }
    base$b_j <- unname(beta)
    rm(skeleton)
    invisible(gc())
    base
  }

  support_tables <- function(pop, out) {
    R <- stats::cor(pop[GEO_VARS])
    csv(data.frame(predictor = rownames(R), R, check.names = FALSE),
        file.path(out, "geography_component_correlations.csv"))
    vif <- do.call(rbind, lapply(GEO_VARS, function(v) {
      r2 <- summary(stats::lm(stats::reformulate(setdiff(GEO_VARS, v), response = v), data = pop))$r.squared
      data.frame(predictor = v, VIF = 1 / (1 - r2))
    }))
    csv(vif, file.path(out, "geography_component_VIF.csv"))
    by_species <- do.call(rbind, lapply(split(pop, pop$SPECIES), function(z) {
      data.frame(SPECIES = z$SPECIES[1], n_populations = nrow(z),
        latitude_min = min(z$latitude_deg), latitude_max = max(z$latitude_deg),
        latitude_span = diff(range(z$latitude_deg)), longitude_min = min(z$longitude_deg),
        longitude_max = max(z$longitude_deg), longitude_span = diff(range(z$longitude_deg)))
    }))
    csv(by_species, file.path(out, "species_geographic_support.csv"))
    invisible(vif)
  }

  prepare_one <- function(cfg, analysis) {
    out <- file.path(cfg$out, analysis)
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    phase(cfg, "PREPARING", analysis, "Reading exact environmental prepared input (no fitted model)")
    source_file <- file.path(cfg$environment_run, analysis, "prepared_input.rds")
    source <- readRDS(source_file)
    s <- check_source_input(source, analysis)
    sites <- recover_coordinates(source$data)
    dec <- decompose_coordinates(source$data, sites)
    need(identical(dec$data[names(source$data)], source$data), "Source columns/row order changed.")
    formula <- make_formula(source$formula, s)
    X <- design_matrix(formula, dec$data)
    need(ncol(X) == if (length(s$control)) 11L else 10L, "Unexpected geographic design size.")
    need(identical(random_terms_text(formula), random_terms_text(source$formula)),
         "Random-effects structure changed.")
    beta <- stats::lm.fit(X, dec$data[[s$response]])$coefficients
    need(all(is.finite(beta)), "Non-finite OLS initialization.")
    start <- starting_parameters(source, formula, dec$data, beta)
    p <- list(version = VERSION, signature = cfg$signature, analysis = analysis, spec = s,
      data = dec$data, populations = dec$populations, sites = sites, scaling = dec$scaling,
      formula = formula, source_formula = source$formula, mesh = source$mesh,
      source_file = source_file, start = start,
      sample_digest = digest::digest(source$data, algo = "sha256"))
    csv(dec$populations, file.path(out, "population_geography_components.csv"))
    csv(dec$scaling, file.path(out, "geography_component_scaling.csv"))
    csv(sites, file.path(out, "site_coordinates.csv"))
    support_tables(dec$populations, out)
    if (!is.null(source$classification_audit))
      csv(source$classification_audit, file.path(out, "voltinism_subset_audit.csv"))
    if (!is.null(source$missing_voltinism_populations) && nrow(source$missing_voltinism_populations))
      csv(source$missing_voltinism_populations, file.path(out, "source_unclassified_populations.csv"))
    csv(data.frame(analysis = analysis, n_rows = nrow(p$data), n_populations = nrow(p$populations),
      n_species = length(unique(p$data$SPECIES)), fixed_columns = ncol(X), fixed_rank = qr(X)$rank,
      same_rows_and_original_columns = TRUE, same_random_terms = TRUE,
      same_mesh_projection = TRUE, source_data_digest = p$sample_digest),
      file.path(out, "preparation_checks.csv"))
    writeLines(c(paste("Source:", source_file), paste("Source formula:", text_formula(source$formula)),
      paste("Geographic formula:", text_formula(formula)),
      "Gaussian identity ML; annual IID Matern field; unchanged 15-km mesh.",
      "Geographic gradients are per 10 degrees; temperature anomaly units are inherited unchanged.",
      "OLS coefficients initialize optimization only; every parameter is re-estimated.",
      "No environmental fixed coefficients or predicted plasticities enter this model.",
      "No AIC/LRT comparison with the environmental fit is exported: the fixed designs differ."),
      file.path(out, "model_specification.txt"))
    atomic_save(p, file.path(out, "prepared_input.rds"))
    p
  }

  fit_checks <- function(fit) {
    grad <- tryCatch(as.numeric(fit$tmb_obj$gr(fit$model$par)), error = function(e) NA_real_)
    max_grad <- if (length(grad) && all(is.finite(grad))) max(abs(grad)) else Inf
    conv <- if (length(fit$model$convergence) == 1L) as.integer(fit$model$convergence) else NA_integer_
    pd <- isTRUE(fit$sd_report$pdHess)
    ok <- identical(conv, 0L) && pd && is.finite(max_grad) && max_grad < .001
    data.frame(convergence_code = conv, pdHess = pd, max_abs_gradient = max_grad,
      numerical_checks_ok = ok, logLik = as.numeric(stats::logLik(fit)),
      AIC = stats::AIC(fit), nobs = nrow(fit$data))
  }

  fixed_estimates <- function(fit, p) {
    terms <- colnames(design_matrix(p$formula, p$data))
    ix <- which(names(fit$sd_report$par.fixed) == "b_j")
    need(length(ix) == length(terms), "Cannot align fitted fixed effects.")
    beta <- stats::setNames(as.numeric(fit$sd_report$par.fixed[ix]), terms)
    V <- as.matrix(fit$sd_report$cov.fixed[ix, ix, drop = FALSE])
    dimnames(V) <- list(terms, terms)
    need(all(is.finite(beta)) && all(is.finite(V)) && all(diag(V) > 0),
         "Invalid fitted fixed-effect covariance; numerical review required.")
    list(beta = beta, V = V)
  }

  linear_estimates <- function(L, beta, V) {
    need(is.matrix(L) && identical(colnames(L), names(beta)) &&
           identical(rownames(V), names(beta)) && identical(colnames(V), names(beta)),
         "Linear contrast/covariance names are not aligned.")
    est <- as.numeric(L %*% beta)
    variance <- rowSums((L %*% V) * L)
    need(all(is.finite(variance)) && all(variance >= -1e-10), "Invalid contrast variances.")
    se <- sqrt(pmax(variance, 0))
    data.frame(estimate = est, std.error = se, conf.low = est - 1.96 * se,
      conf.high = est + 1.96 * se, p_Wald = 2 * stats::pnorm(-abs(est / se)))
  }

  interaction_terms <- function(beta, anomaly) {
    stats::setNames(vapply(GEO_VARS, function(v) {
      term <- intersect(c(paste(anomaly, v, sep = ":"), paste(v, anomaly, sep = ":")), names(beta))
      need(length(term) == 1L, paste("Missing geographic interaction:", v))
      term
    }, character(1)), GEO_VARS)
  }

  orient <- function(tab, s) {
    tab$estimate_biological <- s$orientation * tab$estimate
    tab$conf.low_biological <- pmin(s$orientation * tab$conf.low, s$orientation * tab$conf.high)
    tab$conf.high_biological <- pmax(s$orientation * tab$conf.low, s$orientation * tab$conf.high)
    tab$biological_response <- s$orientation_label
    tab
  }

  geographic_results <- function(beta, V, s, analysis, populations) {
    identity <- diag(length(beta))
    dimnames(identity) <- list(names(beta), names(beta))
    terms <- interaction_terms(beta, s$anomaly)
    key <- linear_estimates(identity[unname(terms), , drop = FALSE], beta, V)
    key$analysis <- analysis
    key$axis <- c("latitude", "latitude", "longitude", "longitude")
    key$effect_scale <- c("within species", "between species", "within species", "between species")
    key$term <- unname(terms)
    key$units <- "change in days per source-anomaly unit per 10 degrees"
    key <- orient(key, s)
    differences <- lapply(c("latitude", "longitude"), function(axis) {
      L <- identity[terms[paste0(axis, "_within_species_10deg")], , drop = FALSE] -
           identity[terms[paste0(axis, "_species_mean_10deg")], , drop = FALSE]
      z <- orient(linear_estimates(L, beta, V), s)
      z$analysis <- analysis
      z$contrast <- paste0(axis, "_within_minus_between")
      z
    })
    groups <- list(within_geography = terms[c(1, 3)], between_geography = terms[c(2, 4)],
                   all_geography = terms)
    joint <- do.call(rbind, lapply(names(groups), function(nm) {
      idx <- unname(groups[[nm]])
      cov <- V[idx, idx, drop = FALSE]
      ev <- eigen(cov, symmetric = TRUE, only.values = TRUE)$values
      valid <- all(is.finite(ev)) && min(ev) > max(ev) * 1e-10
      stat <- if (valid) as.numeric(crossprod(beta[idx], solve(cov, beta[idx]))) else NA_real_
      data.frame(analysis = analysis, test = nm, df = length(idx), Wald_chisq = stat,
        p_Wald = if (valid) stats::pchisq(stat, df = length(idx), lower.tail = FALSE) else NA_real_,
        covariance_invertible = valid)
    }))
    # Fixed systematic sensitivity evaluated ONLY at observed population contexts.
    # These are not independent observations and do not include random species slopes.
    L <- matrix(0, nrow(populations), length(beta), dimnames = list(NULL, names(beta)))
    L[, s$anomaly] <- 1
    L[, unname(terms)] <- as.matrix(populations[GEO_VARS])
    slopes <- orient(linear_estimates(L, beta, V), s)
    slopes$p_Wald <- NULL
    slopes <- cbind(populations[c("SPECIES", "SITE_ID", "latitude_deg", "longitude_deg")], slopes)
    slopes$analysis <- analysis
    slopes$units <- "days per source-anomaly unit"
    slopes$estimand <- "fixed geographic component; excludes species random slope"
    list(effects = key, differences = do.call(rbind, differences), joint = joint, slopes = slopes)
  }

  export_results <- function(fit, p, out) {
    checks <- fit_checks(fit)
    csv(checks, file.path(out, "fit_checks.csv"))
    fx <- fixed_estimates(fit, p)
    I <- diag(length(fx$beta))
    colnames(I) <- names(fx$beta)
    fixed <- cbind(term = names(fx$beta), linear_estimates(I, fx$beta, fx$V))
    csv(fixed, file.path(out, "fixed_effects.csv"))
    csv(data.frame(term = rownames(fx$V), fx$V, check.names = FALSE),
        file.path(out, "fixed_effect_covariance.csv"))
    z <- geographic_results(fx$beta, fx$V, p$spec, p$analysis, p$populations)
    files <- c(effects = "geography_within_between_effects.csv",
      differences = "within_minus_between_contrasts.csv", joint = "joint_geography_tests.csv",
      slopes = "population_fixed_geographic_sensitivity.csv")
    for (nm in names(files)) {
      z[[nm]]$numerical_checks_ok <- checks$numerical_checks_ok
      csv(z[[nm]], file.path(out, files[[nm]]))
    }
    writeLines(tryCatch(capture.output(sdmTMB::sanity(fit, gradient_thresh = .001)),
      error = function(e) paste("Sanity export failed:", conditionMessage(e))), file.path(out, "sanity.txt"))
    tryCatch(csv(as.data.frame(sdmTMB::tidy(fit, effects = "ran_pars", conf.int = TRUE)),
      file.path(out, "random_parameters.csv")), error = function(e)
        writeLines(conditionMessage(e), file.path(out, "random_parameters_error.txt")))
    writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
    if (isTRUE(checks$numerical_checks_ok)) "FIT_OK" else "NUMERICAL_REVIEW_REQUIRED"
  }

  validate_cache <- function(x, cfg, analysis) {
    need(identical(x$signature, cfg$signature) && identical(x$version, VERSION) &&
           identical(x$analysis, analysis), "Geographic cache version/signature/analysis mismatch.")
    invisible(TRUE)
  }

  completed_cache_status <- function(cfg, analysis) {
    out <- file.path(cfg$out, analysis)
    path <- file.path(out, "completed.rds")
    if (!file.exists(path)) return(NULL)
    done <- read_small(path)
    validate_cache(done, cfg, analysis)
    need(done$status %in% c("FIT_OK", "NUMERICAL_REVIEW_REQUIRED"), "Unknown cached completion status.")
    need(all(file.exists(file.path(out, c("prepared_input.rds", "fit_geography.rds")))),
         "Incomplete completed cache.")
    required <- c("geography_within_between_effects.csv", "within_minus_between_contrasts.csv",
      "joint_geography_tests.csv", "population_fixed_geographic_sensitivity.csv",
      "fixed_effects.csv", "fixed_effect_covariance.csv", "fit_checks.csv")
    # Missing exports are regenerated from the saved fit; never silently accepted.
    if (!all(file.exists(file.path(out, required)))) return(NULL)
    done$status
  }

  fit_one <- function(cfg, analysis) {
    require_sources_current(cfg)
    check_packages()
    resource_check(cfg$out, cfg$minimum_ram_GiB, cfg$minimum_disk_GiB)
    out <- file.path(cfg$out, analysis)
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    done_file <- file.path(out, "completed.rds")
    fit_file <- file.path(out, "fit_geography.rds")
    prep_file <- file.path(out, "prepared_input.rds")
    cached <- completed_cache_status(cfg, analysis)
    if (!is.null(cached)) {
      phase(cfg, "REUSING", analysis, cached)
      return(cached)
    }
    p <- if (file.exists(prep_file)) readRDS(prep_file) else prepare_one(cfg, analysis)
    validate_cache(p, cfg, analysis)
    if (file.exists(fit_file)) {
      phase(cfg, "EXPORTING_CACHED", analysis, "Saved fit reused; no optimization")
      saved <- readRDS(fit_file)
      validate_cache(saved, cfg, analysis)
      fit <- saved$fit
    } else {
      phase(cfg, "FITTING", analysis, paste(nrow(p$data), "rows | direct geographic Mundlak model"))
      options(sdmTMB.cores = 1L)
      TMB::openmp(n = 1L, DLL = "sdmTMB")
      set.seed(20261016L)
      warnings <- character()
      fit <- withCallingHandlers(sdmTMB::sdmTMB(formula = p$formula, data = p$data,
        mesh = p$mesh, time = "YEAR", family = stats::gaussian(link = "identity"),
        spatial = "off", spatiotemporal = "iid", reml = FALSE, silent = FALSE,
        control = sdmTMB::sdmTMBcontrol(start = p$start, multiphase = FALSE,
          nlminb_loops = 1L, newton_loops = 0L, parallel = 1L,
          get_joint_precision = FALSE, collapse_spatial_variance = FALSE)),
        warning = function(w) {
          warnings <<- unique(c(warnings, conditionMessage(w)))
          message("WARNING: ", conditionMessage(w))
          invokeRestart("muffleWarning")
        })
      fit$tmb_obj$fn(fit$model$par)
      fit$parlist <- fit$tmb_obj$env$parList(par = fit$tmb_obj$env$last.par.best)
      fit$last.par.best <- fit$tmb_obj$env$last.par.best
      atomic_save(list(version = VERSION, signature = cfg$signature, analysis = analysis,
        fit = fit, warnings = warnings), fit_file)
      writeLines(warnings, file.path(out, "warnings.txt"))
    }
    phase(cfg, "EXPORTING", analysis, "Geographic gradients, joint tests and covariance")
    status <- export_results(fit, p, out)
    atomic_save(list(version = VERSION, signature = cfg$signature, analysis = analysis,
      status = status, time = Sys.time()), done_file)
    phase(cfg, status, analysis)
    status
  }

  combine_outputs <- function(cfg, statuses) {
    files <- c("geography_within_between_effects.csv", "within_minus_between_contrasts.csv",
                "joint_geography_tests.csv", "fit_checks.csv", "preparation_checks.csv")
    for (f in files) {
      parts <- lapply(ANALYSES, function(a) {
        path <- file.path(cfg$out, a, f)
        if (!file.exists(path)) return(NULL)
        x <- utils::read.csv(path, stringsAsFactors = FALSE)
        x$analysis <- a
        x
      })
      parts <- Filter(Negate(is.null), parts)
      if (length(parts)) csv(do.call(rbind, parts), file.path(cfg$out, paste0("ALL_", f)))
    }
    csv(statuses, file.path(cfg$out, "workflow_status.csv"))
    writeLines(c("DIRECT GEOGRAPHIC PATTERNS OF PHENOLOGICAL PLASTICITY",
      "Annual phenology response; no environmental moderators in the fixed design.",
      "Same environmental-source rows, anomaly units/windows, random structure and 15-km mesh.",
      "Within/between latitude and longitude interactions, per 10 degrees.",
      "Primary result: ALL_geography_within_between_effects.csv.",
      "Joint tests: within (2 df), between (2 df), all geography (4 df), Wald tests.",
      "Within-minus-between contrasts explicitly test differences between scales.",
      "Population sensitivity predictions are fixed geographic components, not independent data.",
      "They exclude random species slopes; CIs propagate the full fixed covariance.",
      "No new inference comparing significance codes or treating predictions as observations.",
      "FIT_OK denotes numerical convergence only; inspect residual/spatial diagnostics separately."),
      file.path(cfg$out, "README_results.txt"))
    files <- list.files(cfg$out, pattern = "\\.(csv|txt)$", recursive = TRUE, full.names = FALSE)
    zipfile <- file.path(cfg$out, "geography_mundlak_results_to_review.zip")
    if (requireNamespace("zip", quietly = TRUE)) {
      zip::zipr(zipfile, files = files, root = cfg$out, mode = "mirror", include_directories = FALSE)
    } else {
      message("Optional package zip unavailable; review tables remain in ", cfg$out)
    }
    invisible(zipfile)
  }

  dump_engine <- function(out) {
    e <- environment(dump_engine)
    nms <- ls(e, all.names = TRUE)
    funs <- nms[vapply(nms, function(nm) is.function(get(nm, e)), logical(1))]
    path <- file.path(out, "geography_mundlak_engine.R")
    dump(c("VERSION", "SOURCE_VERSION", "ROOT", "ENV_RELATIVE", "OUTPUT_NAME", "ANALYSES", "V2_ANALYSES",
           "GEO_VARS", funs), file = path, envir = e, control = "all")
    normalizePath(path, winslash = "/", mustWork = TRUE)
  }

  worker <- function(cfg) {
    # The controller remains small; each fit runs in a fresh serial child process.
    tryCatch({
      # Parent transfers the launch lock before the controller may complete.
      # This also prevents a fast cached run from deleting/recreating that lock.
      transferred <- FALSE
      for (attempt in seq_len(100L)) {
        owner <- tryCatch(read_small(file.path(cfg$lock_dir, "owner.rds")), error = function(e) NULL)
        if (!is.null(owner) && identical(as.integer(owner$pid), as.integer(Sys.getpid()))) {
          transferred <- TRUE
          break
        }
        Sys.sleep(.1)
      }
      need(transferred, "Controller did not receive the launch lock.")
      require_sources_current(cfg)
      need(identical(unname(tools::md5sum(cfg$source_files)), cfg$source_md5),
           "Source files changed after launch configuration.")
      statuses <- data.frame(analysis = character(), status = character())
      for (analysis in ANALYSES) {
        require_sources_current(cfg)
        cached <- completed_cache_status(cfg, analysis)
        if (!is.null(cached)) {
          result <- cached
          phase(cfg, "REUSING", analysis, result)
        } else {
          resource_check(cfg$out, cfg$minimum_ram_GiB, cfg$minimum_disk_GiB)
          log_file <- file.path(cfg$out, paste0(analysis, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
          child <- callr::r_bg(function(engine, cfg, analysis) {
            e <- new.env(parent = globalenv())
            sys.source(engine, envir = e)
            e$fit_one(cfg, analysis)
          }, args = list(cfg$engine, cfg, analysis), libpath = cfg$libpath,
            stdout = log_file, stderr = "2>&1", supervise = TRUE,
            user_profile = FALSE, system_profile = FALSE, wd = cfg$root,
            env = c(callr::rcmd_safe_env(), OMP_NUM_THREADS = "1",
                    OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1"))
          h <- ps::ps_handle(child$get_pid())
          atomic_save(list(pid = child$get_pid(), create_time = ps::ps_create_time(h),
            analysis = analysis, log = log_file), file.path(cfg$out, "fit_owner.rds"))
          child$wait()
          result <- child$get_result()
        }
        statuses <- rbind(statuses, data.frame(analysis = analysis, status = result))
        csv(statuses, file.path(cfg$out, "workflow_status_partial.csv"))
      }
      zip <- combine_outputs(cfg, statuses)
      final <- if (all(statuses$status == "FIT_OK")) "COMPLETED" else "COMPLETED_WITH_NUMERICAL_REVIEW"
      phase(cfg, final, detail = zip, terminal = TRUE)
      list(status = final, out = cfg$out)
    }, error = function(e) {
      msg <- conditionMessage(e)
      try(writeLines(msg, file.path(cfg$out, "fatal_error.txt")), silent = TRUE)
      try(phase(cfg, "STOPPED_WITH_ERROR", detail = msg, terminal = TRUE), silent = TRUE)
      list(status = "STOPPED_WITH_ERROR", error = msg, out = cfg$out)
    }, finally = {
      # Do not leave a child running if the controller encounters an R error.
      if (exists("child", inherits = FALSE)) {
        try(if (child$is_alive()) child$kill_tree(), silent = TRUE)
      }
      unlink(cfg$lock_dir, recursive = TRUE)
    })
  }

  acquire_launch_lock <- function(root) {
    lock <- file.path(root, "launch_lock")
    if (dir.exists(lock)) {
      owner <- tryCatch(read_small(file.path(lock, "owner.rds")), error = function(e) NULL)
      need(!is.null(owner), paste("Launch lock has no readable owner; inspect before retrying:", lock))
      need(!pid_alive(owner$pid, owner$create_time), "A launch/controller already owns the geographic lock.")
      unlink(lock, recursive = TRUE)
    }
    need(dir.create(lock, showWarnings = FALSE), "Could not acquire geographic launch lock.")
    h <- ps::ps_handle()
    atomic_save(list(pid = ps::ps_pid(h), create_time = ps::ps_create_time(h)), file.path(lock, "owner.rds"))
    lock
  }

  run <- function(project_root = ROOT, environment_run = file.path(project_root, ENV_RELATIVE),
                   minimum_ram_GiB = 35, minimum_disk_GiB = 50) {
    need(is.numeric(minimum_ram_GiB) && length(minimum_ram_GiB) == 1L &&
           is.finite(minimum_ram_GiB) && minimum_ram_GiB >= 35, "RAM threshold must be at least 35 GiB.")
    need(is.numeric(minimum_disk_GiB) && length(minimum_disk_GiB) == 1L &&
           is.finite(minimum_disk_GiB) && minimum_disk_GiB >= 20, "Disk threshold must be at least 20 GiB.")
    project_root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
    environment_run <- normalizePath(environment_run, winslash = "/", mustWork = TRUE)
    src <- source_status(environment_run, require_complete = TRUE)
    jobs <- active_jobs(project_root, environment_run)
    need(!length(jobs), paste("Another known phenoIMPACT workflow is active; no new fit launched:",
                             paste(jobs, collapse = "\n")))
    check_packages()
    resource_check(project_root, minimum_ram_GiB, minimum_disk_GiB)
    out_root <- file.path(project_root, "output", "phenology_plasticity", OUTPUT_NAME)
    dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
    lock <- acquire_launch_lock(out_root)
    launched <- FALSE
    on.exit(if (!launched) unlink(lock, recursive = TRUE), add = TRUE)
    source_files <- src$files
    source_manifest <- file_manifest(c(source_files, src$table$fit))
    hashes <- unname(tools::md5sum(source_files))
    need(!anyNA(hashes), "Could not fingerprint source files.")
    e <- environment(run)
    funs <- sort(ls(e, all.names = TRUE))
    funs <- funs[vapply(funs, function(nm) is.function(get(nm, e)), logical(1))]
    code <- vapply(funs, function(nm) paste(deparse(get(nm, e)), collapse = "\n"), character(1))
    versions <- vapply(packages(), function(nm) as.character(utils::packageVersion(nm)), character(1))
    signature <- digest::digest(list(version = VERSION, source_signature = src$config$signature,
      source_files = source_files, md5 = hashes, selected_sources = src$table,
      analyses = ANALYSES, variables = GEO_VARS,
      code = code, packages = versions), algo = "sha256")
    out <- file.path(out_root, paste0("run_", substr(signature, 1, 12)))
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    cfg <- list(version = VERSION, signature = signature, root = project_root,
      environment_run = environment_run, source_files = source_files, source_md5 = hashes,
      source_selection = src$table, source_manifest = source_manifest,
      source_signature = src$config$signature, out = out, lock_dir = lock,
      minimum_ram_GiB = minimum_ram_GiB, minimum_disk_GiB = minimum_disk_GiB,
      libpath = .libPaths(), package_versions = versions)
    require_sources_current(cfg)
    cfg$engine <- dump_engine(out)
    atomic_save(cfg, file.path(out, "run_config.rds"))
    csv(data.frame(file = source_files, md5 = hashes), file.path(out, "source_manifest.csv"))
    csv(src$table, file.path(out, "environmental_source_selection.csv"))
    if (file.exists(file.path(out, "fatal_error.txt"))) unlink(file.path(out, "fatal_error.txt"))
    phase(cfg, "QUEUED", detail = paste(ANALYSES, collapse = " -> "))
    log <- file.path(out, paste0("controller_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
    child <- callr::r_bg(function(engine, cfg) {
      e <- new.env(parent = globalenv())
      sys.source(engine, envir = e)
      e$worker(cfg)
    }, args = list(cfg$engine, cfg), libpath = cfg$libpath, stdout = log, stderr = "2>&1",
      supervise = FALSE, cleanup = FALSE, cleanup_tree = FALSE, poll_connection = FALSE,
      user_profile = FALSE, system_profile = FALSE, wd = project_root,
      env = c(callr::rcmd_safe_env(), OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1"))
    # Until ownership metadata is saved, a failed launch must kill its controller.
    on.exit(if (!launched) try(child$kill_tree(), silent = TRUE), add = TRUE)
    h <- ps::ps_handle(child$get_pid())
    owner <- list(pid = child$get_pid(), create_time = ps::ps_create_time(h),
                   started = Sys.time(), out = out, log = log)
    atomic_save(owner, file.path(out, "worker_owner.rds"))
    atomic_save(owner, file.path(lock, "owner.rds"))
    job <- c(owner, list(process = child))
    state$job <- job
    options(phenoimpact.geo.mundlak.job = job)
    launched <- TRUE
    message("Four geographic fits queued sequentially. Output: ", out,
            "\nMonitor with pheno_geo_mundlak$status().")
    invisible(job)
  }

  status <- function(job = NULL, project_root = ROOT) {
    if (is.null(job)) job <- state$job
    if (is.null(job)) job <- getOption("phenoimpact.geo.mundlak.job")
    if (is.null(job)) {
      root <- file.path(project_root, "output", "phenology_plasticity", OUTPUT_NAME)
      runs <- list.dirs(root, recursive = FALSE, full.names = TRUE)
      runs <- runs[startsWith(basename(runs), "run_")]
      need(length(runs) > 0, "No geographic run found.")
      files <- file.path(runs, "worker_owner.rds")
      files <- files[file.exists(files)]
      need(length(files) > 0, "No saved geographic worker metadata.")
      job <- read_small(files[which.max(file.info(files)$mtime)])
    }
    progress <- file.path(job$out, "progress.rds")
    if (file.exists(progress)) print(read_small(progress))
    alive <- if (!is.null(job$process)) tryCatch(job$process$is_alive(), error = function(e) FALSE) else
      pid_alive(job$pid, job$create_time)
    cat("\nController alive:", alive, "\nOutput:", job$out, "\n")
    fit_owner <- file.path(job$out, "fit_owner.rds")
    if (file.exists(fit_owner)) {
      owner <- read_small(fit_owner)
      cat("Fit process alive:", pid_alive(owner$pid, owner$create_time),
          "\nAnalysis:", owner$analysis, "\nFit log:", owner$log, "\n")
    }
    partial <- file.path(job$out, "workflow_status_partial.csv")
    if (file.exists(partial)) print(utils::read.csv(partial), row.names = FALSE)
    err <- file.path(job$out, "fatal_error.txt")
    if (!alive && file.exists(err)) cat(paste(readLines(err, warn = FALSE), collapse = "\n"), "\n")
    zip <- file.path(job$out, "geography_mundlak_results_to_review.zip")
    if (!alive && file.exists(zip)) cat("\nResults to share:\n", zip, "\n", sep = "")
    invisible(job)
  }

  list(plan = plan, preflight = preflight, run = run, status = status)
})
