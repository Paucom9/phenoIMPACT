# =============================================================================
# phenoIMPACT | Direct FLIGHT-LENGTH plasticity with annual spatial control
# Version 1.1 recovery patch | 2026-10-12
# IMPORTANT: internal model VERSION/signature intentionally remains 1.0 so
# the completed run_bc5b3c03313a univoltine checkpoint is reused exactly.
#
# PURPOSE
# Finalize direct flight-length plasticity before testing its demographic effect.
#
# Input:
#   output/phenology_plasticity/flight_length_direct_pilot/
#     flight_length_model_data.csv
#
# Model, fitted separately for univoltines and multivoltines:
#
# FLIGHT_LENGTH_mean ~ season_temp_anomaly_C *
#   (season_photoperiod_z + season_background_z +
#    season_predictability_z + season_trend_z) +
#   (1 | SITE_ID) + (1 | SPECIES) +
#   (0 + season_temp_anomaly_C | SPECIES_slope) +
#   (1 | site_year_id)
#
# Gaussian ML + annual IID Matern field:
#   spatial = "off"
#   spatiotemporal = "iid"
#   time = "YEAR"
#   mesh cutoff = 15 km
#
# Positive contextual population plasticity = LONGER flight period in warmer
# years (days per degree C).
#
# The contextual slope is:
#   fixed anomaly slope
#   + species anomaly random-slope deviation
#   + sum(anomaly x moderator beta * population moderator value)
#
# This exactly follows the population-plasticity extraction logic used for the
# reviewed onset/offset spatial models: no population random slope is invented,
# and the additive spatial field is not part of the temperature derivative.
#
# This script DOES NOT fit abundance trends yet. It produces the spatially
# corrected population flight-length plasticities needed for the next model:
#   abundance trend ~ onset plasticity + flight-length plasticity
#
# Run:
#   source(file.choose(), encoding = "UTF-8")
#   flight_length_spatial_job <- pheno_flight_length_spatial$run(monitor = FALSE)
#   pheno_flight_length_spatial$status()
#
# Keep R open. Do not run another large TMB/sdmTMB job concurrently.
# =============================================================================

pheno_flight_length_spatial <- local({

  VERSION <- "flight_length_spatial_15km_1.0"
  ROOT <- "E:/phenoIMPACT project/code/phenoIMPACT"
  INPUT <- file.path(
    ROOT, "output", "phenology_plasticity",
    "flight_length_direct_pilot",
    "flight_length_model_data.csv"
  )
  OUTPUT_ROOT <- file.path(
    ROOT, "output", "phenology_plasticity",
    "flight_length_direct_spatial_15km"
  )
  GROUPS <- c("univoltine", "multivoltine")
  GRADIENT_THRESHOLD <- 0.001
  EXTRA_OPTIMIZATION_ATTEMPTS <- 2L
  MESH_CUTOFF_KM <- 15
  FIT_SEED <- 20261008L

  state <- new.env(parent = emptyenv())

  need <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)
  csv <- function(x, path) {
    utils::write.csv(x, path, row.names = FALSE, na = "")
  }
  ftext <- function(x) paste(deparse(x, width.cutoff = 500L), collapse = " ")

  atomic_rds <- function(x, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    tmp <- tempfile(".writing_", tmpdir = dirname(path), fileext = ".rds")
    on.exit(unlink(tmp), add = TRUE)
    saveRDS(x, tmp, compress = "gzip")
    bak <- paste0(path, ".previous")
    if (file.exists(bak)) unlink(bak)
    had <- file.exists(path)
    if (had) need(file.rename(path, bak), paste("Cannot back up:", path))
    if (!file.rename(tmp, path)) {
      if (had) file.rename(bak, path)
      stop("Cannot finalize: ", path, call. = FALSE)
    }
    if (file.exists(bak)) unlink(bak)
    invisible(path)
  }

  check_packages <- function() {
    pkgs <- c(
      "sdmTMB", "TMB", "Matrix", "fmesher",
      "lme4", "dplyr", "digest", "callr"
    )
    miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
    need(!length(miss),
         paste("Missing packages:", paste(miss, collapse = ", ")))
    # The static random-effect extractor below is the mapping already validated
    # in this project for sdmTMB 1.1.0.
    need(
      identical(as.character(utils::packageVersion("sdmTMB")), "1.1.0"),
      paste0(
        "This extractor is validated for sdmTMB 1.1.0; found ",
        as.character(utils::packageVersion("sdmTMB")),
        ". Stop rather than silently remap random effects."
      )
    )
    invisible(TRUE)
  }

  package_versions <- function() {
    pkgs <- c("sdmTMB", "TMB", "Matrix", "fmesher", "lme4")
    data.frame(
      package = pkgs,
      version = vapply(
        pkgs,
        function(p) as.character(utils::packageVersion(p)),
        character(1)
      ),
      stringsAsFactors = FALSE
    )
  }

  progress <- function(path, stage, group = "", detail = "", terminal = FALSE) {
    atomic_rds(
      list(
        time = Sys.time(),
        stage = stage,
        group = group,
        detail = detail,
        terminal = terminal
      ),
      path
    )
    message(
      "[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ",
      if (nzchar(group)) paste0(group, " | ") else "",
      stage,
      if (nzchar(detail)) paste0(" | ", detail) else ""
    )
  }

  normal_ci <- function(x) {
    need(all(c("term", "estimate", "std.error") %in% names(x)),
         "Unexpected fixed-effect table.")
    z <- stats::qnorm(.975)
    x$conf.low <- x$estimate - z * x$std.error
    x$conf.high <- x$estimate + z * x$std.error
    x$p_Wald <- 2 * stats::pnorm(-abs(x$estimate / x$std.error))
    x
  }

  formula_model <- function() {
    stats::as.formula(
      paste0(
        "FLIGHT_LENGTH_mean ~ season_temp_anomaly_C * (",
        "season_photoperiod_z + season_background_z + ",
        "season_predictability_z + season_trend_z) + ",
        "(1 | SITE_ID) + (1 | SPECIES) + ",
        "(0 + season_temp_anomaly_C | SPECIES_slope) + ",
        "(1 | site_year_id)"
      ),
      env = baseenv()
    )
  }

  fixed_formula <- function() {
    stats::as.formula(
      paste0(
        "~ season_temp_anomaly_C * (",
        "season_photoperiod_z + season_background_z + ",
        "season_predictability_z + season_trend_z)"
      ),
      env = baseenv()
    )
  }

  validate_data <- function(d, group) {
    req <- c(
      "SPECIES", "SITE_ID", "bms_id", "voltinism", "YEAR",
      "FLIGHT_LENGTH_mean", "season_temp_anomaly_C",
      "season_photoperiod_z", "season_background_z",
      "season_predictability_z", "season_trend_z",
      "x_km", "y_km", "source_row"
    )
    need(all(req %in% names(d)),
         paste("Missing input columns for", group))
    need(nrow(d) > 0L, paste("No rows for", group))
    need(!anyDuplicated(d[c("SPECIES", "SITE_ID", "YEAR")]),
         paste("Duplicate population-year rows:", group))

    nums <- c(
      "YEAR", "FLIGHT_LENGTH_mean", "season_temp_anomaly_C",
      "season_photoperiod_z", "season_background_z",
      "season_predictability_z", "season_trend_z",
      "x_km", "y_km", "source_row"
    )
    need(
      all(vapply(d[nums], function(x) all(is.finite(as.numeric(x))), logical(1))),
      paste("Non-finite model values:", group)
    )

    X <- stats::model.matrix(fixed_formula(), d)
    need(
      nrow(X) == nrow(d) &&
        all(is.finite(X)) &&
        qr(X)$rank == ncol(X),
      paste("Fixed design is invalid/rank-deficient:", group)
    )
    invisible(TRUE)
  }

  prepare_group <- function(all_data, group) {
    d <- all_data[all_data$voltinism == group, , drop = FALSE]
    d$SPECIES <- factor(as.character(d$SPECIES))
    d$SITE_ID <- factor(as.character(d$SITE_ID))
    d$bms_id <- as.character(d$bms_id)
    d$YEAR <- as.integer(d$YEAR)
    d$SPECIES_slope <- d$SPECIES
    d$site_year_id <- interaction(
      d$SITE_ID, d$YEAR, drop = TRUE, lex.order = TRUE
    )
    d$source_row <- as.integer(d$source_row)
    validate_data(d, group)

    cn <- c("SPECIES", "SITE_ID", "YEAR")
    d <- d[order(as.character(d$SPECIES),
                 as.character(d$SITE_ID),
                 d$YEAR), , drop = FALSE]
    need(!anyDuplicated(d[cn]), "Unexpected duplicates after ordering.")
    rownames(d) <- NULL
    d
  }

  build_mesh <- function(d) {
    xy <- unique(d[c("x_km", "y_km")])
    set.seed(FIT_SEED)
    tri <- sdmTMB::make_mesh(
      xy,
      xy_cols = c("x_km", "y_km"),
      cutoff = MESH_CUTOFF_KM,
      type = "cutoff"
    )
    mesh <- sdmTMB::make_mesh(
      d,
      xy_cols = c("x_km", "y_km"),
      mesh = tri$mesh
    )
    need(
      nrow(mesh$loc_xy) == nrow(d) &&
        isTRUE(all.equal(
          unname(as.matrix(mesh$loc_xy)),
          unname(as.matrix(d[c("x_km", "y_km")])),
          check.attributes = FALSE
        )),
      "Mesh/data projection is misaligned."
    )
    list(triangulation = tri, projected = mesh)
  }

  fit_checks <- function(fit) {
    g <- fit$gradients
    grad <- if (length(g) && all(is.finite(g))) max(abs(g)) else NA_real_
    code <- fit$model$convergence
    if (length(code) != 1L) code <- NA_integer_
    pd <- isTRUE(fit$sd_report$pdHess)

    fixed <- tryCatch(
      sdmTMB::tidy(fit, effects = "fixed", conf.int = FALSE),
      error = function(e) NULL
    )
    se_ok <- !is.null(fixed) &&
      all(c("estimate", "std.error") %in% names(fixed)) &&
      nrow(fixed) > 0L &&
      all(is.finite(fixed$estimate)) &&
      all(is.finite(fixed$std.error) & fixed$std.error > 0)

    sanity <- tryCatch(
      sdmTMB::sanity(
        fit,
        gradient_thresh = GRADIENT_THRESHOLD,
        silent = TRUE
      ),
      error = function(e) NULL
    )
    sanity_ok <- !is.null(sanity) && isTRUE(all(unlist(sanity)))

    data.frame(
      convergence_code = code,
      positive_definite_Hessian = pd,
      max_abs_gradient = grad,
      finite_fixed_SE = se_ok,
      all_sanity_checks_pass = sanity_ok,
      fit_ok =
        isTRUE(code == 0L) &&
        pd &&
        is.finite(grad) &&
        grad < GRADIENT_THRESHOLD &&
        se_ok,
      stringsAsFactors = FALSE
    )
  }

  extract_coefficients <- function(fit) {
    fixed_names <- colnames(stats::model.matrix(
      fit$split_formula[[1]]$form_no_bars,
      fit$data[seq_len(min(2L, nrow(fit$data))), , drop = FALSE]
    ))
    ix <- which(names(fit$sd_report$par.fixed) == "b_j")
    need(
      length(ix) == length(fixed_names) &&
        !anyDuplicated(fixed_names),
      "Fixed-coefficient mapping failed."
    )

    beta <- stats::setNames(
      as.numeric(fit$sd_report$par.fixed[ix]),
      fixed_names
    )
    V <- fit$sd_report$cov.fixed[ix, ix, drop = FALSE]
    need(
      all(is.finite(beta)) &&
        all(is.finite(V)) &&
        all(diag(V) > 0),
      "Invalid fixed estimates/covariance."
    )
    dimnames(V) <- list(fixed_names, fixed_names)

    sf <- fit$split_formula[[1]]
    blocks <- sf$re_cov_terms$re_b_df
    need(
      all(c("start", "end", "group_indices", "level") %in% names(blocks)),
      "Saved random-effect mapping is missing."
    )

    all_ix <- unlist(
      Map(function(a, b) seq.int(a, b), blocks$start, blocks$end),
      use.names = FALSE
    )
    values <- fit$sd_report$value[
      grepl("^re_b_pars($|\\[)", names(fit$sd_report$value))
    ]
    errors <- fit$sd_report$sd[
      grepl("^re_b_pars($|\\[)", names(fit$sd_report$value))
    ]
    active <- !is.na(fit$tmb_map$re_b_pars)

    need(
      length(values) == length(active) &&
        length(errors) == length(active) &&
        identical(
          as.integer(all_ix + 1L),
          which(as.vector(active))
        ),
      "Random-effect parameter blocks do not match the saved report."
    )

    group_index <- match(
      "SPECIES_slope",
      names(sf$re_cov_terms$cnms)
    )
    need(!is.na(group_index), "SPECIES_slope random effect not found.")

    b <- blocks[
      blocks$group_indices == group_index,
      ,
      drop = FALSE
    ]
    need(
      nrow(b) > 0L &&
        all(b$start == b$end) &&
        !anyDuplicated(as.character(b$level)),
      "Expected exactly one anomaly random slope per species."
    )

    slopes <- data.frame(
      SPECIES = as.character(b$level),
      sp_slope_dev = as.numeric(values[b$start + 1L]),
      sp_slope_dev_marginal_SE =
        as.numeric(errors[b$start + 1L]),
      stringsAsFactors = FALSE
    )

    need(
      setequal(
        slopes$SPECIES,
        unique(as.character(fit$data$SPECIES))
      ) &&
        all(is.finite(slopes$sp_slope_dev)) &&
        all(is.finite(slopes$sp_slope_dev_marginal_SE)),
      "Species random-slope extraction failed."
    )

    list(beta = beta, covariance = V, slopes = slopes)
  }

  population_plasticity <- function(d, coef, group) {
    mods <- c(
      "season_photoperiod_z",
      "season_background_z",
      "season_predictability_z",
      "season_trend_z"
    )
    anomaly <- "season_temp_anomaly_C"

    context <- d |>
      dplyr::group_by(SPECIES, SITE_ID, bms_id) |>
      dplyr::summarise(
        n_years = dplyr::n_distinct(YEAR),
        first_year = min(YEAR),
        last_year = max(YEAR),
        anomaly_SD_within_population =
          stats::sd(season_temp_anomaly_C),
        dplyr::across(dplyr::all_of(mods), mean),
        .groups = "drop"
      ) |>
      as.data.frame()

    context$SPECIES <- as.character(context$SPECIES)
    context$SITE_ID <- as.character(context$SITE_ID)
    context$bms_id <- as.character(context$bms_id)
    context$pop_id <- paste(
      context$SPECIES, context$SITE_ID, sep = "_"
    )
    need(!anyDuplicated(context$pop_id),
         "Population-ID collision.")

    idx <- match(context$SPECIES, coef$slopes$SPECIES)
    need(!anyNA(idx), "Species slope join failed.")
    context$species_slope_dev <-
      coef$slopes$sp_slope_dev[idx]
    context$species_slope_dev_marginal_SE <-
      coef$slopes$sp_slope_dev_marginal_SE[idx]

    need(anomaly %in% names(coef$beta),
         "Fixed anomaly slope missing.")
    context$fixed_anomaly_slope <-
      unname(coef$beta[[anomaly]])
    context$interaction_contribution <- 0

    for (mod in mods) {
      term <- intersect(
        c(
          paste(anomaly, mod, sep = ":"),
          paste(mod, anomaly, sep = ":")
        ),
        names(coef$beta)
      )
      need(length(term) == 1L,
           paste("Missing/ambiguous interaction:", mod))
      val <- as.numeric(coef$beta[[term]]) * context[[mod]]
      context[[paste0("contribution__", mod)]] <- val
      context$interaction_contribution <-
        context$interaction_contribution + val
    }

    context$flight_length_plasticity_days_per_C <-
      context$fixed_anomaly_slope +
      context$species_slope_dev +
      context$interaction_contribution

    need(
      all(is.finite(context$flight_length_plasticity_days_per_C)),
      "Non-finite contextual flight-length plasticity."
    )

    mu <- mean(context$flight_length_plasticity_days_per_C)
    ss <- stats::sd(context$flight_length_plasticity_days_per_C)
    need(is.finite(ss) && ss > 0,
         "Invalid flight-length plasticity SD.")

    context$flight_length_plasticity_z <-
      (context$flight_length_plasticity_days_per_C - mu) / ss
    context$voltinism <- group

    scaling <- data.frame(
      group = group,
      n_populations = nrow(context),
      mean_days_per_C = mu,
      SD_days_per_C = ss,
      scaling_universe =
        "all populations in the accepted spatial flight-length fit; before abundance join",
      stringsAsFactors = FALSE
    )

    list(populations = context, scaling = scaling)
  }

  compare_existing <- function(pop, group, out_dir) {
    pilot_file <- file.path(
      ROOT, "output", "phenology_plasticity",
      "flight_length_direct_pilot",
      "pilot_population_flight_length_plasticity.csv"
    )
    axes_file <- file.path(
      ROOT, "output", "population_trends",
      "spatial_abundance_qc_offset_final_15km",
      "run_eee755f6aac4",
      "global_onset_offset_balance",
      "population_global_balance_axes.csv"
    )

    rows <- list()

    if (file.exists(pilot_file)) {
      old <- utils::read.csv(
        pilot_file,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      old <- old[old$voltinism == group, , drop = FALSE]
      keep <- c(
        "pop_id",
        "flight_length_plasticity_days_per_C",
        "flight_length_plasticity_z"
      )
      if (all(keep %in% names(old))) {
        old_cmp <- old[, keep, drop = FALSE]
        names(old_cmp)[
          names(old_cmp) == "flight_length_plasticity_days_per_C"
        ] <- "pilot_days_per_C"
        names(old_cmp)[
          names(old_cmp) == "flight_length_plasticity_z"
        ] <- "pilot_z"

        # Fail explicitly if parsing ever changes; never pass text to cor().
        need(
          is.numeric(old_cmp$pilot_days_per_C) &&
            is.numeric(old_cmp$pilot_z),
          "Pilot plasticity columns are not numeric after CSV import."
        )

        x <- merge(
          pop[, c(
            "pop_id",
            "flight_length_plasticity_days_per_C",
            "flight_length_plasticity_z"
          )],
          old_cmp[, c("pop_id", "pilot_days_per_C", "pilot_z")],
          by = "pop_id",
          all = FALSE
        )
        rows[["pilot"]] <- data.frame(
          group = group,
          comparison = "spatial_vs_nonspatial_pilot",
          n_populations = nrow(x),
          correlation = stats::cor(
            x$flight_length_plasticity_z,
            x$pilot_z
          ),
          stringsAsFactors = FALSE
        )
        csv(
          x,
          file.path(out_dir, "comparison_nonspatial_pilot.csv")
        )
      }
    }

    if (file.exists(axes_file)) {
      ax <- utils::read.csv(
        axes_file,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      ax <- unique(ax[
        ax$group == group,
        c(
          "pop_id",
          "season_extension",
          "coordinated_advancement"
        )
      ])
      x <- merge(
        pop[, c("pop_id", "flight_length_plasticity_z")],
        ax,
        by = "pop_id",
        all = FALSE
      )
      if (nrow(x)) {
        need(
          is.numeric(x$flight_length_plasticity_z) &&
            is.numeric(x$season_extension) &&
            is.numeric(x$coordinated_advancement),
          "Previous balance-axis columns are not numeric after CSV import."
        )
        rows[["extension"]] <- data.frame(
          group = group,
          comparison =
            "spatial_flight_length_vs_relative_season_extension",
          n_populations = nrow(x),
          correlation = stats::cor(
            x$flight_length_plasticity_z,
            x$season_extension
          ),
          stringsAsFactors = FALSE
        )
        rows[["advancement"]] <- data.frame(
          group = group,
          comparison =
            "spatial_flight_length_vs_coordinated_advancement",
          n_populations = nrow(x),
          correlation = stats::cor(
            x$flight_length_plasticity_z,
            x$coordinated_advancement
          ),
          stringsAsFactors = FALSE
        )
        csv(
          x,
          file.path(out_dir, "comparison_previous_balance_axes.csv")
        )
      }
    }

    if (length(rows)) {
      out <- do.call(rbind, rows)
      csv(out, file.path(out_dir, "plasticity_comparisons.csv"))
      out
    } else {
      data.frame()
    }
  }

  fit_one <- function(d, group, group_dir, progress_file, run_id) {
    dir.create(group_dir, recursive = TRUE, showWarnings = FALSE)
    fit_file <- file.path(group_dir, "fit_spatial_field.rds")
    meta_file <- file.path(group_dir, "metadata.rds")
    input_file <- file.path(group_dir, "prepared_input.rds")

    form <- formula_model()

    signature <- digest::digest(
      list(
        version = VERSION,
        group = group,
        data = d,
        formula = ftext(form),
        mesh_cutoff_km = MESH_CUTOFF_KM,
        spatial = "off",
        spatiotemporal = "iid",
        time = "YEAR",
        reml = FALSE,
        seed = FIT_SEED,
        gradient_threshold = GRADIENT_THRESHOLD,
        package_versions = package_versions()
      ),
      algo = "sha256"
    )

    if (file.exists(fit_file) &&
        file.exists(meta_file) &&
        file.exists(input_file)) {
      meta <- readRDS(meta_file)
      saved <- readRDS(fit_file)
      old_input <- readRDS(input_file)
      need(
        identical(meta$signature, signature) &&
          identical(saved$signature, signature) &&
          isTRUE(all.equal(
            old_input,
            d,
            check.attributes = TRUE
          )),
        paste("Checkpoint mismatch:", group)
      )
      checks <- fit_checks(saved$fit)
      saved$checks <- checks
      need(isTRUE(checks$fit_ok),
           paste("Cached spatial fit is not numerically accepted:", group))
      message("Reusing accepted spatial flight-length fit: ", group)
      return(saved)
    }

    progress(
      progress_file,
      "BUILDING_MESH",
      group,
      paste(nrow(d), "observations")
    )
    bm <- build_mesh(d)
    tri <- bm$triangulation
    mesh <- bm$projected

    validation <- data.frame(
      group = group,
      n_observations = nrow(d),
      n_populations = length(unique(
        paste(as.character(d$SPECIES),
              as.character(d$SITE_ID), sep = "_")
      )),
      n_species = nlevels(d$SPECIES),
      n_sites = nlevels(d$SITE_ID),
      n_site_years = nlevels(d$site_year_id),
      n_mesh_vertices = tri$mesh$n,
      cutoff_km = MESH_CUTOFF_KM,
      stringsAsFactors = FALSE
    )
    csv(validation, file.path(group_dir, "validation.csv"))
    atomic_rds(d, input_file)
    atomic_rds(
      list(
        version = VERSION,
        group = group,
        run_id = run_id,
        signature = signature,
        formula = form,
        triangulation = tri,
        package_versions = package_versions(),
        validation = validation
      ),
      meta_file
    )

    progress(
      progress_file,
      "FITTING",
      group,
      "Gaussian ML + annual IID Matern field"
    )

    warns <- character()
    capture_warning <- function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }

    started <- proc.time()[["elapsed"]]
    set.seed(FIT_SEED)

    fit <- withCallingHandlers(
      sdmTMB::sdmTMB(
        formula = form,
        data = d,
        mesh = mesh,
        time = "YEAR",
        family = stats::gaussian(link = "identity"),
        spatial = "off",
        spatiotemporal = "iid",
        reml = FALSE,
        silent = FALSE
      ),
      warning = capture_warning
    )

    checks <- fit_checks(fit)
    elapsed <- (proc.time()[["elapsed"]] - started) / 60

    saved <- list(
      version = VERSION,
      group = group,
      run_id = run_id,
      signature = signature,
      fit = fit,
      checks = checks,
      warnings = unique(warns),
      extra_optimization = 0L,
      elapsed_minutes = elapsed
    )
    atomic_rds(saved, fit_file)

    for (attempt in seq_len(EXTRA_OPTIMIZATION_ATTEMPTS)) {
      if (isTRUE(checks$fit_ok)) break
      progress(
        progress_file,
        "EXTRA_OPTIMIZATION",
        group,
        paste("attempt", attempt, "of", EXTRA_OPTIMIZATION_ATTEMPTS)
      )

      improved <- tryCatch(
        withCallingHandlers(
          sdmTMB::run_extra_optimization(
            fit,
            nlminb_loops = 1,
            newton_loops = 0
          ),
          warning = capture_warning
        ),
        error = function(e) {
          warns <<- c(
            warns,
            paste("Extra optimization failed:", conditionMessage(e))
          )
          NULL
        }
      )
      if (is.null(improved)) break

      fit <- improved
      checks <- fit_checks(fit)
      elapsed <- (proc.time()[["elapsed"]] - started) / 60

      saved$fit <- fit
      saved$checks <- checks
      saved$warnings <- unique(warns)
      saved$extra_optimization <- attempt
      saved$elapsed_minutes <- elapsed
      atomic_rds(saved, fit_file)
    }

    csv(
      cbind(
        group = group,
        checks,
        elapsed_minutes = saved$elapsed_minutes,
        extra_optimization = saved$extra_optimization
      ),
      file.path(group_dir, "fit_checks.csv")
    )
    writeLines(
      saved$warnings,
      file.path(group_dir, "warnings.txt")
    )

    need(isTRUE(saved$checks$fit_ok),
         paste(
           "Spatial flight-length fit did not pass numerical checks:",
           group,
           "- checkpoint retained."
         ))

    progress(
      progress_file,
      "EXTRACTING_PLASTICITY",
      group,
      "accepted spatial fit"
    )

    fixed <- normal_ci(
      sdmTMB::tidy(
        saved$fit,
        effects = "fixed",
        conf.int = FALSE
      )
    )
    fixed$group <- group
    csv(fixed, file.path(group_dir, "fixed_effects.csv"))

    ran <- sdmTMB::tidy(
      saved$fit,
      effects = "ran_pars",
      conf.int = TRUE
    )
    ran$group <- group
    csv(ran, file.path(group_dir, "random_parameters.csv"))

    coef <- extract_coefficients(saved$fit)
    csv(
      coef$slopes,
      file.path(group_dir, "species_random_slopes.csv")
    )
    csv(
      data.frame(
        term = rownames(coef$covariance),
        coef$covariance,
        check.names = FALSE
      ),
      file.path(group_dir, "fixed_effect_covariance.csv")
    )

    pp <- population_plasticity(d, coef, group)
    csv(
      pp$populations,
      file.path(
        group_dir,
        "population_flight_length_plasticity_spatial.csv"
      )
    )
    csv(
      pp$scaling,
      file.path(group_dir, "plasticity_scaling.csv")
    )

    comparison <- compare_existing(
      pp$populations,
      group,
      group_dir
    )

    focal_terms <- c(
      "season_temp_anomaly_C",
      paste0(
        "season_temp_anomaly_C:",
        c(
          "season_photoperiod_z",
          "season_background_z",
          "season_predictability_z",
          "season_trend_z"
        )
      )
    )
    focal <- fixed[fixed$term %in% focal_terms, , drop = FALSE]
    csv(focal, file.path(group_dir, "focal_temperature_effects.csv"))

    writeLines(
      c(
        "phenoIMPACT | Direct flight-length plasticity | spatial 15 km",
        paste("Group:", group),
        "",
        "Positive population plasticity = longer flight period in warmer years.",
        "Units = days per degree C because season_temp_anomaly_C is unstandardized.",
        "",
        "Model: Gaussian ML; site/species intercepts; independent species anomaly slope;",
        "site-year intercept; annual IID Matern field; spatial='off'; time=YEAR.",
        "",
        "Contextual population slope = fixed anomaly effect + species anomaly slope",
        "+ all four anomaly x context interactions evaluated at population context.",
        "No population random slope is added. The additive spatial field is not part",
        "of the derivative with respect to temperature anomaly.",
        "",
        capture.output(print(
          cbind(
            group = group,
            saved$checks,
            elapsed_minutes = saved$elapsed_minutes,
            extra_optimization = saved$extra_optimization
          ),
          row.names = FALSE
        )),
        "",
        "FOCAL FIXED EFFECTS:",
        capture.output(print(focal, row.names = FALSE)),
        "",
        "PLASTICITY SCALING:",
        capture.output(print(pp$scaling, row.names = FALSE)),
        "",
        "COMPARISONS:",
        if (nrow(comparison)) {
          capture.output(print(comparison, row.names = FALSE))
        } else {
          "No comparison files were available."
        }
      ),
      file.path(group_dir, "README_results.txt")
    )

    saved
  }

  worker <- function(cfg) {
    check_packages()
    dir.create(cfg$out, recursive = TRUE, showWarnings = FALSE)
    progress_file <- file.path(cfg$out, "progress.rds")

    progress(progress_file, "READING_INPUT")

    need(file.exists(cfg$input),
         paste("Missing pilot model data:", cfg$input))
    d0 <- utils::read.csv(
      cfg$input,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    need(
      all(GROUPS %in% unique(d0$voltinism)),
      "Pilot model data do not contain both voltinism groups."
    )

    input_manifest <- data.frame(
      path = normalizePath(cfg$input, winslash = "/", mustWork = TRUE),
      md5 = unname(tools::md5sum(cfg$input)),
      stringsAsFactors = FALSE
    )
    csv(input_manifest, file.path(cfg$out, "input_manifest.csv"))
    csv(package_versions(), file.path(cfg$out, "package_versions.csv"))

    status_rows <- list()
    populations <- list()
    scalings <- list()
    comparisons <- list()
    fixed_tables <- list()

    for (g in GROUPS) {
      result <- tryCatch({
        progress(progress_file, "PREPARING", g)
        d <- prepare_group(d0, g)

        group_hash <- digest::digest(
          list(
            version = VERSION,
            input_md5 = input_manifest$md5,
            group = g,
            rows = d$source_row,
            cutoff = MESH_CUTOFF_KM,
            seed = FIT_SEED
          ),
          algo = "sha256"
        )
        group_dir <- file.path(
          cfg$out,
          paste0(g, "_", substr(group_hash, 1, 12))
        )
        dir.create(group_dir, recursive = TRUE, showWarnings = FALSE)

        saved <- fit_one(
          d = d,
          group = g,
          group_dir = group_dir,
          progress_file = progress_file,
          run_id = cfg$run_id
        )

        pop_file <- file.path(
          group_dir,
          "population_flight_length_plasticity_spatial.csv"
        )
        populations[[g]] <- utils::read.csv(
          pop_file,
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
        scalings[[g]] <- utils::read.csv(
          file.path(group_dir, "plasticity_scaling.csv"),
          stringsAsFactors = FALSE
        )
        fixed_tables[[g]] <- utils::read.csv(
          file.path(group_dir, "fixed_effects.csv"),
          stringsAsFactors = FALSE
        )
        # Recompute lightweight comparison outputs even when the expensive
        # spatial fit is reused from checkpoint. This is the recovery path for
        # the v1.0 post-processing failure after the univoltine fit completed.
        comparisons[[g]] <- compare_existing(
          populations[[g]],
          g,
          group_dir
        )

        # Also regenerate the focal fixed-effect table from the saved fixed CSV.
        focal_terms <- c(
          "season_temp_anomaly_C",
          paste0(
            "season_temp_anomaly_C:",
            c(
              "season_photoperiod_z",
              "season_background_z",
              "season_predictability_z",
              "season_trend_z"
            )
          )
        )
        focal_cached <- fixed_tables[[g]][
          fixed_tables[[g]]$term %in% focal_terms,
          ,
          drop = FALSE
        ]
        csv(
          focal_cached,
          file.path(group_dir, "focal_temperature_effects.csv")
        )

        data.frame(
          group = g,
          status = "FIT_OK",
          n_observations = nrow(d),
          n_populations = nrow(populations[[g]]),
          output_directory = group_dir,
          message = "",
          stringsAsFactors = FALSE
        )
      }, error = function(e) {
        data.frame(
          group = g,
          status = "ERROR",
          n_observations = NA_integer_,
          n_populations = NA_integer_,
          output_directory = "",
          message = conditionMessage(e),
          stringsAsFactors = FALSE
        )
      })

      status_rows[[g]] <- result
      csv(
        do.call(rbind, status_rows),
        file.path(cfg$out, "workflow_status.csv")
      )

      if (result$status != "FIT_OK") {
        progress(
          progress_file,
          "ERROR",
          g,
          result$message,
          terminal = TRUE
        )
        stop(result$message, call. = FALSE)
      }

      rm(list = intersect(c("d", "saved"), ls()))
      invisible(gc())
    }

    all_pop <- do.call(rbind, populations)
    all_scale <- do.call(rbind, scalings)
    all_fixed <- do.call(rbind, fixed_tables)
    nonempty_comp <- comparisons[
      vapply(comparisons, nrow, integer(1)) > 0L
    ]
    all_comp <- if (length(nonempty_comp)) {
      do.call(rbind, nonempty_comp)
    } else data.frame()

    csv(
      all_pop,
      file.path(
        cfg$out,
        "all_population_flight_length_plasticity_spatial.csv"
      )
    )
    csv(
      all_scale,
      file.path(cfg$out, "all_plasticity_scaling.csv")
    )
    csv(
      all_fixed,
      file.path(cfg$out, "all_fixed_effects.csv")
    )
    if (nrow(all_comp)) {
      csv(
        all_comp,
        file.path(cfg$out, "all_plasticity_comparisons.csv")
      )
    }

    progress(
      progress_file,
      "FITS_COMPLETE",
      "",
      paste(
        "Both spatial flight-length fits accepted;",
        "population plasticities ready for abundance-trend model"
      ),
      terminal = TRUE
    )

    list(
      output = cfg$out,
      status = do.call(rbind, status_rows)
    )
  }

  make_cfg <- function() {
    check_packages()
    need(file.exists(INPUT), paste("Missing:", INPUT))
    dir.create(OUTPUT_ROOT, recursive = TRUE, showWarnings = FALSE)

    sig <- digest::digest(
      list(
        version = VERSION,
        input = normalizePath(INPUT, winslash = "/", mustWork = TRUE),
        md5 = unname(tools::md5sum(INPUT)),
        groups = GROUPS,
        gradient = GRADIENT_THRESHOLD,
        extra = EXTRA_OPTIMIZATION_ATTEMPTS,
        cutoff = MESH_CUTOFF_KM,
        seed = FIT_SEED,
        packages = package_versions()
      ),
      algo = "sha256"
    )

    out <- file.path(
      OUTPUT_ROOT,
      paste0("run_", substr(sig, 1, 12))
    )
    dir.create(out, recursive = TRUE, showWarnings = FALSE)

    list(
      version = VERSION,
      run_id = sig,
      input = INPUT,
      out = out
    )
  }

  run <- function(monitor = FALSE) {
    check_packages()

    old <- getOption("pheno.flight.length.spatial.active_job")
    if (!is.null(old) &&
        isTRUE(tryCatch(old$process$is_alive(), error = function(e) FALSE))) {
      stop(
        "A flight-length spatial job is already running. Use pheno_flight_length_spatial$status().",
        call. = FALSE
      )
    }

    cfg <- make_cfg()
    p <- callr::r_bg(
      function(worker_fun, cfg) worker_fun(cfg),
      args = list(worker_fun = worker, cfg = cfg),
      stdout = file.path(cfg$out, "worker.log"),
      stderr = "2>&1",
      supervise = TRUE,
      wd = ROOT
    )

    job <- list(
      process = p,
      output = cfg$out,
      progress_file = file.path(cfg$out, "progress.rds"),
      started = Sys.time()
    )
    state$job <- job
    options(pheno.flight.length.spatial.active_job = job)

    message("Spatial flight-length job started.")
    message("Output: ", cfg$out)

    if (isTRUE(monitor)) {
      repeat {
        status()
        if (!p$is_alive()) break
        Sys.sleep(60)
      }
    }

    invisible(job)
  }

  status <- function() {
    job <- if (exists("job", envir = state, inherits = FALSE)) {
      state$job
    } else {
      getOption("pheno.flight.length.spatial.active_job")
    }
    need(!is.null(job), "No flight-length spatial job is available.")

    alive <- tryCatch(job$process$is_alive(), error = function(e) FALSE)
    elapsed <- as.numeric(
      difftime(Sys.time(), job$started, units = "hours")
    )

    prog <- NULL
    if (file.exists(job$progress_file)) {
      prog <- tryCatch(readRDS(job$progress_file), error = function(e) NULL)
    }

    cat(
      "Elapsed:", round(elapsed, 2), "h | controller running:", alive, "\n"
    )
    if (!is.null(prog)) {
      cat(
        prog$stage, " | ", prog$group, " | ", prog$detail, "\n",
        sep = ""
      )
    }

    st_file <- file.path(job$output, "workflow_status.csv")
    if (file.exists(st_file)) {
      print(utils::read.csv(st_file), row.names = FALSE)
    }
    cat("Output:", job$output, "\n")

    if (!alive) {
      exit <- tryCatch(
        job$process$get_exit_status(),
        error = function(e) NA_integer_
      )
      if (!is.na(exit) && exit != 0L) {
        cat("Worker exit status:", exit, "\n")
        log_file <- file.path(job$output, "worker.log")
        if (file.exists(log_file)) {
          z <- readLines(log_file, warn = FALSE)
          cat("\nLast worker log lines:\n")
          cat(tail(z, 25), sep = "\n")
          cat("\n")
        }
      }
    }

    invisible(list(
      alive = alive,
      progress = prog,
      output = job$output
    ))
  }

  stop <- function() {
    job <- if (exists("job", envir = state, inherits = FALSE)) {
      state$job
    } else {
      getOption("pheno.flight.length.spatial.active_job")
    }
    if (is.null(job)) return(invisible(FALSE))
    if (isTRUE(tryCatch(job$process$is_alive(), error = function(e) FALSE))) {
      job$process$kill_tree()
    }
    message("Worker stopped. Completed checkpoints are retained.")
    invisible(TRUE)
  }

  list(
    run = run,
    status = status,
    stop = stop
  )
})
