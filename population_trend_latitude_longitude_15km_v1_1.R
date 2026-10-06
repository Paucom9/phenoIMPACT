# =============================================================================
# phenoIMPACT | Geographic patterns in model-based population abundance trends
# Version 1.1 | 2026-10-05
#
# PURPOSE
#   Test whether model-based population abundance trends show within-species
#   geographic structure:
#
#       relative population trend ~ within-species latitude + longitude
#
#   separately for univoltine and multivoltine butterflies.
#
# IMPORTANT
#   - NO abundance model is refitted.
#   - NO phenology/plasticity model is refitted.
#   - Reads the CURRENT onset-adjusted 15-km spatial abundance models.
#   - Univoltine MUST come from the accepted refined fit.
#   - Multivoltine comes from the reviewed source run.
#   - The annual spatial field and population AR1 are NOT interpreted as a
#     linear population trend. The extracted linear trend is:
#
#       fixed year slope
#       + year x onset-plasticity
#       + year x offset-plasticity
#       + species temporal random slope
#       + site temporal random slope
#
#     This matches the interpretation used previously for model-based linear
#     population trends. Annual AR1 and annual spatial deviations describe
#     year-specific departures from that long-term linear component.
#
# INFERENCE
#   The abundance trends are model-derived estimates. Their full joint
#   uncertainty (especially conditional species/site random-slope covariance)
#   was not stored by the source workflow. Therefore this script DOES NOT
#   pretend that a second-stage ordinary regression has exact model-based SEs.
#
#   Instead it:
#     1) estimates the geographic projection on the fitted population trends;
#     2) gives a SPECIES-CLUSTER BOOTSTRAP interval, preserving all populations
#        of a resampled species together;
#     3) repeats the estimate with equal-population and equal-species weighting;
#     4) checks leave-one-species-out and leave-one-network-out sensitivity;
#     5) screens residual spatial autocorrelation (Moran's I) if spdep exists.
#
#   Bootstrap intervals quantify robustness across the sampled species, NOT
#   full first-stage model uncertainty.
#
# OUTPUT
#   output/population_trends/latitude_population_trends_15km/run_*/
#       population_trends_with_geography.csv
#       geographic_summary.csv
#       species_support.csv
#       leave_one_species_out.csv
#       leave_one_network_out.csv
#       residual_spatial_screen.csv
#       population_trend_latitude.png / .pdf
#       README.txt
#       results_to_review.zip
#
# RUN
#   source(file.choose(), encoding = "UTF-8")
#   pop_geo <- run_population_trend_geography()
#
# =============================================================================

run_population_trend_geography <- function(
  project_root = "E:/phenoIMPACT project/code/phenoIMPACT",
  source_run = file.path(
    project_root, "output", "population_trends",
    "spatial_abundance_onset_adjusted_15km",
    "run_20261004_214012_44136"
  ),
  n_boot = 2000L,
  seed = 20261005L,
  make_figure = TRUE
) {

  stopf <- function(...) stop(..., call. = FALSE)
  need <- function(x, ...) if (!isTRUE(x)) stopf(...)
  same_num <- function(x, y, tol = 1e-8) {
    is.numeric(x) && is.numeric(y) && length(x) == length(y) &&
      all(is.finite(x)) && all(is.finite(y)) &&
      (length(x) == 0L || max(abs(as.numeric(x) - as.numeric(y))) <= tol)
  }
  csv <- function(x, path) utils::write.csv(x, path, row.names = FALSE, na = "")

  need(requireNamespace("sf", quietly = TRUE), "Package 'sf' is required.")
  if (make_figure) {
    need(requireNamespace("ggplot2", quietly = TRUE), "Package 'ggplot2' is required for the figure.")
  }

  project_root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
  source_run <- normalizePath(source_run, winslash = "/", mustWork = TRUE)

  need(dir.exists(source_run), "Source abundance run not found: ", source_run)

  # ---------------------------------------------------------------------------
  # Locate CURRENT reviewed fits
  # ---------------------------------------------------------------------------

  multi_fit_path <- file.path(source_run, "multivoltine__full", "fit.rds")
  multi_input_path <- file.path(source_run, "input_multivoltine.rds")
  uni_input_path <- file.path(source_run, "input_univoltine.rds")

  need(file.exists(multi_fit_path), "Current multivoltine fit missing: ", multi_fit_path)
  need(file.exists(multi_input_path), "Multivoltine input missing.")
  need(file.exists(uni_input_path), "Univoltine input missing.")

  refinement_root <- file.path(source_run, "refinement_univoltine")
  need(dir.exists(refinement_root),
    "No univoltine refinement directory found. This workflow refuses to fall back to the pre-refinement fit.")

  refinement_dirs <- list.dirs(refinement_root, recursive = FALSE, full.names = TRUE)
  refinement_dirs <- refinement_dirs[vapply(refinement_dirs, function(d) {
    result <- file.path(d, "result.rds")
    fit <- file.path(d, "univoltine__full", "fit.rds")
    if (!all(file.exists(c(result, fit)))) return(FALSE)
    z <- tryCatch(readRDS(result), error = function(e) NULL)
    !is.null(z) && identical(z$status, "FIT_OK") &&
      identical(normalizePath(z$fit_file, winslash = "/", mustWork = TRUE),
                normalizePath(fit, winslash = "/", mustWork = TRUE))
  }, logical(1))]

  need(length(refinement_dirs) == 1L,
    "Expected exactly one accepted FIT_OK univoltine refinement; found ",
    length(refinement_dirs), ". Review: ", refinement_root)

  uni_refinement <- normalizePath(refinement_dirs, winslash = "/", mustWork = TRUE)
  uni_fit_path <- file.path(uni_refinement, "univoltine__full", "fit.rds")

  message("Source abundance run: ", source_run)
  message("Univoltine fit: ", uni_fit_path)
  message("Multivoltine fit: ", multi_fit_path)

  # ---------------------------------------------------------------------------
  # New output directory
  # ---------------------------------------------------------------------------

  out_root <- file.path(
    project_root, "output", "population_trends",
    "latitude_population_trends_15km"
  )
  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

  out <- file.path(
    out_root,
    paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid())
  )
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  out <- normalizePath(out, winslash = "/", mustWork = TRUE)

  # ---------------------------------------------------------------------------
  # Helpers
  # ---------------------------------------------------------------------------

  interaction_name <- function(beta, variable) {
    z <- intersect(
      c(paste0("year_decade:", variable),
        paste0(variable, ":year_decade")),
      names(beta)
    )
    need(length(z) == 1L, "Interaction term missing/ambiguous for ", variable)
    z
  }

  site_coordinates <- function(input) {
    s <- as.data.frame(input$sites)
    need(all(c("SITE_ID", "x_km", "y_km") %in% names(s)),
      "Expected SITE_ID, x_km and y_km in input$sites.")
    s$SITE_ID <- as.character(s$SITE_ID)
    need(!anyDuplicated(s$SITE_ID), "Duplicated SITE_ID in input$sites.")
    need(all(is.finite(s$x_km)) && all(is.finite(s$y_km)),
      "Non-finite abundance-model coordinates.")

    sm <- s
    sm$x_m <- sm$x_km * 1000
    sm$y_m <- sm$y_km * 1000

    sfobj <- sf::st_as_sf(sm, coords = c("x_m", "y_m"), crs = 3035)
    ll <- sf::st_coordinates(sf::st_transform(sfobj, 4326))
    s$longitude_deg <- ll[, 1]
    s$latitude_deg <- ll[, 2]

    need(all(is.finite(s$longitude_deg)) && all(is.finite(s$latitude_deg)) &&
      all(abs(s$longitude_deg) <= 180) && all(abs(s$latitude_deg) <= 90),
      "Invalid transformed longitude/latitude.")
    s
  }

  extract_population_trends <- function(group, fit_path, input_path, refined = FALSE) {
    message("Reading ", group, " saved arrays; no refitting.")
    fit <- readRDS(fit_path)
    input <- readRDS(input_path)

    need(is.list(fit) && is.list(input), "Unexpected saved object structure for ", group)
    need(identical(fit$group, group) && identical(fit$variant, "full"),
      "Wrong fit group/variant for ", group)
    need(isTRUE(fit$offset_adjusted_for_onset),
      "Fit is not the onset-adjusted abundance model: ", group)
    need(isTRUE(fit$checks$numerical_checks_ok),
      "Current ", group, " fit does not pass numerical checks.")
    need(identical(fit$data_signature, input$data_signature),
      "Fit/input data signature mismatch for ", group)
    need(identical(input$group, group), "Input group mismatch for ", group)

    beta <- fit$beta
    need(is.numeric(beta) && !is.null(names(beta)) && all(is.finite(beta)),
      "Invalid fixed coefficients for ", group)

    pars <- fit$parameters
    need(is.list(pars) &&
      all(c("b_species", "b_site") %in% names(pars)),
      "Saved random-effect arrays missing for ", group)

    b_species <- as.matrix(pars$b_species)
    b_site <- as.matrix(pars$b_site)
    need(ncol(b_species) == 2L && ncol(b_site) == 2L,
      "Expected intercept + temporal slope for species/site in ", group)
    need(all(is.finite(b_species)) && all(is.finite(b_site)),
      "Non-finite random effects for ", group)

    f <- as.data.frame(input$frame)
    need(all(c(
      "pop_id", "SPECIES", "SITE_ID",
      "onset_plasticity_z", "offset_plasticity_bio_z"
    ) %in% names(f)), "Abundance input frame lacks required columns for ", group)

    need(length(input$data$species) == nrow(f) &&
      length(input$data$site) == nrow(f),
      "Row-to-random-effect indices do not match the input frame.")

    f$SPECIES <- as.character(f$SPECIES)
    f$SITE_ID <- as.character(f$SITE_ID)
    f$pop_id <- as.character(f$pop_id)
    f$species_index <- as.integer(input$data$species) + 1L
    f$site_index <- as.integer(input$data$site) + 1L

    need(all(f$species_index >= 1L & f$species_index <= nrow(b_species)),
      "Species index outside random-effect matrix.")
    need(all(f$site_index >= 1L & f$site_index <= nrow(b_site)),
      "Site index outside random-effect matrix.")

    # Verify the integer random-effect indices identify one and only one level.
    sp_map <- unique(f[c("SPECIES", "species_index")])
    si_map <- unique(f[c("SITE_ID", "site_index")])
    need(!anyDuplicated(sp_map$SPECIES) && !anyDuplicated(sp_map$species_index),
      "Species-to-random-effect mapping is not one-to-one.")
    need(!anyDuplicated(si_map$SITE_ID) && !anyDuplicated(si_map$site_index),
      "Site-to-random-effect mapping is not one-to-one.")

    # Verify site-index ordering against input$sites before using coordinates.
    sites <- site_coordinates(input)
    need(nrow(sites) == nrow(b_site),
      "input$sites and b_site have different row counts.")
    site_by_index <- data.frame(
      site_index = seq_len(nrow(sites)),
      SITE_ID_from_sites = as.character(sites$SITE_ID),
      stringsAsFactors = FALSE
    )
    chk <- merge(si_map, site_by_index, by = "site_index", sort = FALSE)
    need(nrow(chk) == nrow(si_map) &&
      all(chk$SITE_ID == chk$SITE_ID_from_sites),
      "Site random-effect order does not match input$sites.")

    # One row per population.
    key <- c(
      "pop_id", "SPECIES", "SITE_ID",
      "onset_plasticity_z", "offset_plasticity_bio_z",
      "species_index", "site_index"
    )
    p <- unique(f[key])
    need(!anyDuplicated(p$pop_id), "Population has multiple predictor/index combinations: ", group)

    int_on <- interaction_name(beta, "onset_plasticity_z")
    int_off <- interaction_name(beta, "offset_plasticity_bio_z")
    need("year_decade" %in% names(beta), "Fixed year slope missing.")

    p$species_year_slope <- b_species[p$species_index, 2L]
    p$site_year_slope <- b_site[p$site_index, 2L]

    # Model-based long-term linear trend on the fitted log-abundance scale.
    p$trend_log_per_decade <-
      unname(beta["year_decade"]) +
      unname(beta[int_on]) * p$onset_plasticity_z +
      unname(beta[int_off]) * p$offset_plasticity_bio_z +
      p$species_year_slope +
      p$site_year_slope

    p$trend_percent_decade <- 100 * expm1(p$trend_log_per_decade)
    p$voltinism <- group
    p$fit_source <- if (refined) "refined_current" else "current"
    p$fit_file <- normalizePath(fit_path, winslash = "/", mustWork = TRUE)

    # Coordinates/network.
    add <- sites[c("SITE_ID", intersect(c(
      "bms_id", "x_km", "y_km", "longitude_deg", "latitude_deg"
    ), names(sites)))]
    p <- merge(p, add, by = "SITE_ID", all.x = TRUE, sort = FALSE)
    need(nrow(p) > 0L && !anyNA(p$latitude_deg) && !anyNA(p$longitude_deg),
      "Coordinate join failed for ", group)

    p
  }

  centre_within_species <- function(d) {
    d <- d[is.finite(d$trend_log_per_decade) &
      is.finite(d$latitude_deg) & is.finite(d$longitude_deg), , drop = FALSE]
    need(nrow(d) > 2L, "Too few rows after geography checks.")

    # Preallocate full-length columns before assigning species blocks.
    # Without this, assigning the first non-contiguous index block to a new
    # data-frame column can create a shorter vector and fail immediately.
    n <- nrow(d)
    d$latitude_species_mean <- rep(NA_real_, n)
    d$longitude_species_mean <- rep(NA_real_, n)
    d$trend_species_mean <- rep(NA_real_, n)
    d$trend_percent_species_mean <- rep(NA_real_, n)

    sp <- split(seq_len(n), d$SPECIES, drop = TRUE)
    for (idx in sp) {
      d$latitude_species_mean[idx] <- mean(d$latitude_deg[idx])
      d$longitude_species_mean[idx] <- mean(d$longitude_deg[idx])
      d$trend_species_mean[idx] <- mean(d$trend_log_per_decade[idx])
      d$trend_percent_species_mean[idx] <- mean(d$trend_percent_decade[idx])
    }

    need(
      all(is.finite(d$latitude_species_mean)) &&
      all(is.finite(d$longitude_species_mean)) &&
      all(is.finite(d$trend_species_mean)) &&
      all(is.finite(d$trend_percent_species_mean)),
      "Within-species centring failed for at least one population."
    )

    d$latitude_within_species_deg <- d$latitude_deg - d$latitude_species_mean
    d$longitude_within_species_deg <- d$longitude_deg - d$longitude_species_mean
    d$trend_within_species_log_decade <- d$trend_log_per_decade - d$trend_species_mean
    d$trend_within_species_percent_decade <-
      d$trend_percent_decade - d$trend_percent_species_mean

    # Numerical identity checks: centred variables must sum to ~0 within species.
    check_cols <- c(
      "latitude_within_species_deg",
      "longitude_within_species_deg",
      "trend_within_species_log_decade",
      "trend_within_species_percent_decade"
    )
    for (nm in check_cols) {
      sums <- rowsum(d[[nm]], group = d$SPECIES, reorder = FALSE)
      need(max(abs(sums)) < 1e-8 * max(1, n),
        "Within-species centring identity failed for ", nm)
    }

    d
  }

  wls_projection <- function(d, specification = c("lat_lon", "lat_only"),
                             weighting = c("equal_population", "equal_species")) {
    specification <- match.arg(specification)
    weighting <- match.arg(weighting)

    xlat <- d$latitude_within_species_deg / 10
    xlon <- d$longitude_within_species_deg / 10
    y <- d$trend_within_species_log_decade

    X <- if (specification == "lat_lon") cbind(latitude = xlat, longitude = xlon) else
      cbind(latitude = xlat)

    nsp <- table(d$SPECIES)
    w <- if (weighting == "equal_population") rep(1, nrow(d)) else
      1 / as.numeric(nsp[d$SPECIES])

    XtWX <- crossprod(X, w * X)
    XtWy <- crossprod(X, w * y)
    need(qr(XtWX)$rank == ncol(X), "Second-stage geographic design is rank deficient.")

    b <- as.numeric(solve(XtWX, XtWy))
    names(b) <- colnames(X)
    fitted <- as.numeric(X %*% b)
    list(
      beta = b,
      fitted = fitted,
      residual = y - fitted,
      weights = w,
      X = X
    )
  }

  species_bootstrap <- function(d, specification, weighting, n_boot, seed) {
    set.seed(seed)
    species <- unique(d$SPECIES)
    ns <- length(species)
    need(ns >= 10L, "Too few species for species-cluster bootstrap.")

    # Precompute each species contribution to X'WX and X'Wy.
    # We do not fit each species separately because some individual species
    # may not span enough 2-D geographic space for a full latitude+longitude fit.
    parts <- lapply(species, function(s) {
      z <- d[d$SPECIES == s, , drop = FALSE]
      xlat <- z$latitude_within_species_deg / 10
      xlon <- z$longitude_within_species_deg / 10
      X <- if (specification == "lat_lon") cbind(latitude = xlat, longitude = xlon) else
        cbind(latitude = xlat)
      y <- z$trend_within_species_log_decade
      w <- if (weighting == "equal_population") rep(1, nrow(z)) else rep(1 / nrow(z), nrow(z))
      list(A = crossprod(X, w * X), c = crossprod(X, w * y))
    })

    k <- if (specification == "lat_lon") 2L else 1L
    boot <- matrix(NA_real_, n_boot, k,
      dimnames = list(NULL, if (k == 2L) c("latitude", "longitude") else "latitude"))

    for (b in seq_len(n_boot)) {
      draw <- sample.int(ns, ns, replace = TRUE)
      A <- Reduce(`+`, lapply(draw, function(i) parts[[i]]$A))
      cc <- Reduce(`+`, lapply(draw, function(i) parts[[i]]$c))
      if (qr(A)$rank == k) boot[b, ] <- as.numeric(solve(A, cc))
    }

    boot <- boot[stats::complete.cases(boot), , drop = FALSE]
    need(nrow(boot) >= 0.9 * n_boot,
      "More than 10% of species-bootstrap replicates were rank deficient.")
    boot
  }

  summarize_fit <- function(d, group, specification, weighting, n_boot, seed) {
    fit <- wls_projection(d, specification, weighting)
    boot <- species_bootstrap(d, specification, weighting, n_boot, seed)

    rows <- lapply(names(fit$beta), function(term) {
      vals <- boot[, term]
      est <- unname(fit$beta[term])
      data.frame(
        voltinism = group,
        specification = if (specification == "lat_lon")
          "latitude_plus_longitude_within_species" else "latitude_only_within_species",
        weighting = weighting,
        term = term,
        n_populations = nrow(d),
        n_species = length(unique(d$SPECIES)),
        estimate_log_trend_per_decade_per_10deg = est,
        lower_95_species_bootstrap = unname(stats::quantile(vals, .025)),
        upper_95_species_bootstrap = unname(stats::quantile(vals, .975)),
        bootstrap_median = stats::median(vals),
        bootstrap_prop_positive = mean(vals > 0),
        bootstrap_prop_negative = mean(vals < 0),
        growth_factor_ratio_per_10deg = exp(est),
        growth_factor_percent_difference_per_10deg = 100 * expm1(est),
        inference_scope = paste0(
          "Species-cluster bootstrap across fitted population trends; ",
          "does not propagate full first-stage random-effect uncertainty"
        ),
        stringsAsFactors = FALSE
      )
    })
    list(summary = do.call(rbind, rows), fit = fit, boot = boot)
  }

  # ---------------------------------------------------------------------------
  # Read/extract current model-based trends
  # ---------------------------------------------------------------------------

  uni <- extract_population_trends(
    "univoltine", uni_fit_path, uni_input_path, refined = TRUE
  )
  multi <- extract_population_trends(
    "multivoltine", multi_fit_path, multi_input_path, refined = FALSE
  )

  allp <- rbind(uni, multi)
  allp <- do.call(rbind, lapply(split(allp, allp$voltinism), centre_within_species))
  rownames(allp) <- NULL

  csv(allp, file.path(out, "population_trends_with_geography.csv"))

  # ---------------------------------------------------------------------------
  # Species support
  # ---------------------------------------------------------------------------

  support <- do.call(rbind, lapply(split(allp, list(allp$voltinism, allp$SPECIES), drop = TRUE),
    function(z) data.frame(
      voltinism = z$voltinism[1],
      SPECIES = z$SPECIES[1],
      n_populations = nrow(z),
      n_sites = length(unique(z$SITE_ID)),
      latitude_min = min(z$latitude_deg),
      latitude_max = max(z$latitude_deg),
      latitude_span_deg = diff(range(z$latitude_deg)),
      longitude_span_deg = diff(range(z$longitude_deg)),
      trend_mean_log_per_decade = mean(z$trend_log_per_decade),
      trend_mean_percent_per_decade = mean(z$trend_percent_decade),
      stringsAsFactors = FALSE
    )
  ))
  rownames(support) <- NULL
  csv(support, file.path(out, "species_support.csv"))

  # ---------------------------------------------------------------------------
  # Main + sensitivity projections
  # ---------------------------------------------------------------------------

  fits <- list()
  summaries <- list()

  for (g in c("univoltine", "multivoltine")) {
    d <- allp[allp$voltinism == g, , drop = FALSE]
    for (spec in c("lat_lon", "lat_only")) {
      for (wt in c("equal_population", "equal_species")) {
        key <- paste(g, spec, wt, sep = "__")
        z <- summarize_fit(
          d, g, spec, wt, n_boot,
          seed + match(g, c("univoltine", "multivoltine")) * 100L +
            match(spec, c("lat_lon", "lat_only")) * 10L +
            match(wt, c("equal_population", "equal_species"))
        )
        fits[[key]] <- z
        summaries[[key]] <- z$summary
      }
    }
  }

  summary_tab <- do.call(rbind, summaries)
  rownames(summary_tab) <- NULL
  csv(summary_tab, file.path(out, "geographic_summary.csv"))

  # ---------------------------------------------------------------------------
  # Leave-one-species-out, using the preferred lat+lon equal-species projection
  # ---------------------------------------------------------------------------

  loo_species <- list()
  for (g in c("univoltine", "multivoltine")) {
    d <- allp[allp$voltinism == g, , drop = FALSE]
    spp <- unique(d$SPECIES)
    for (s in spp) {
      z <- d[d$SPECIES != s, , drop = FALSE]
      if (length(unique(z$SPECIES)) < 5L) next
      z <- centre_within_species(z)
      f <- tryCatch(
        wls_projection(z, "lat_lon", "equal_species"),
        error = function(e) NULL
      )
      if (is.null(f)) next
      loo_species[[length(loo_species) + 1L]] <- data.frame(
        voltinism = g,
        omitted_species = s,
        n_populations = nrow(z),
        n_species = length(unique(z$SPECIES)),
        latitude_estimate = unname(f$beta["latitude"]),
        longitude_estimate = unname(f$beta["longitude"]),
        stringsAsFactors = FALSE
      )
    }
  }
  loo_species <- do.call(rbind, loo_species)
  csv(loo_species, file.path(out, "leave_one_species_out.csv"))

  # ---------------------------------------------------------------------------
  # Leave-one-network-out
  # Re-centre species AFTER removing the network.
  # ---------------------------------------------------------------------------

  loo_network <- list()
  if ("bms_id" %in% names(allp) && any(nzchar(as.character(allp$bms_id)))) {
    networks <- sort(unique(as.character(allp$bms_id)))
    networks <- networks[!is.na(networks) & nzchar(networks)]
    for (g in c("univoltine", "multivoltine")) {
      base <- allp[allp$voltinism == g, , drop = FALSE]
      for (net in networks) {
        z <- base[as.character(base$bms_id) != net, , drop = FALSE]
        if (nrow(z) < 100L || length(unique(z$SPECIES)) < 5L) next
        z <- centre_within_species(z)
        f <- tryCatch(
          wls_projection(z, "lat_lon", "equal_species"),
          error = function(e) NULL
        )
        if (is.null(f)) next
        loo_network[[length(loo_network) + 1L]] <- data.frame(
          voltinism = g,
          omitted_network = net,
          n_populations = nrow(z),
          n_species = length(unique(z$SPECIES)),
          latitude_estimate = unname(f$beta["latitude"]),
          longitude_estimate = unname(f$beta["longitude"]),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (length(loo_network)) {
    loo_network <- do.call(rbind, loo_network)
    csv(loo_network, file.path(out, "leave_one_network_out.csv"))
  } else {
    csv(data.frame(note = "No network sensitivity available"),
      file.path(out, "leave_one_network_out.csv"))
  }

  # ---------------------------------------------------------------------------
  # Residual spatial screen
  # Preferred fit = latitude + longitude, equal species.
  # Moran is on SITE-MEAN residuals, so duplicate species at the same site do
  # not enter as coincident spatial points.
  # ---------------------------------------------------------------------------

  moran <- list()
  if (requireNamespace("spdep", quietly = TRUE)) {
    set.seed(seed + 9000L)
    for (g in c("univoltine", "multivoltine")) {
      d <- allp[allp$voltinism == g, , drop = FALSE]
      key <- paste(g, "lat_lon", "equal_species", sep = "__")
      f <- fits[[key]]$fit
      d$residual_geo <- f$residual

      site_res <- aggregate(
        residual_geo ~ SITE_ID + x_km + y_km,
        data = d,
        FUN = mean
      )
      site_res <- site_res[stats::complete.cases(site_res), , drop = FALSE]
      n <- nrow(site_res)

      if (n >= 20L) {
        k <- min(8L, n - 1L)
        kn <- spdep::knearneigh(as.matrix(site_res[c("x_km", "y_km")]), k = k)
        nb <- spdep::knn2nb(kn)
        lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
        mc <- spdep::moran.mc(
          site_res$residual_geo, lw,
          nsim = 999L,
          zero.policy = TRUE
        )
        moran[[g]] <- data.frame(
          voltinism = g,
          n_sites = n,
          k_neighbours = k,
          moran_I = as.numeric(mc$statistic),
          permutation_p = mc$p.value,
          note = paste0(
            "Moran screen of site-mean second-stage residuals; ",
            "diagnostic only, not first-stage uncertainty propagation"
          ),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(moran)) {
    moran <- do.call(rbind, moran)
  } else {
    moran <- data.frame(
      note = "Moran screen not run (package spdep unavailable or insufficient sites)."
    )
  }
  csv(moran, file.path(out, "residual_spatial_screen.csv"))

  # ---------------------------------------------------------------------------
  # Figure: preferred lat+lon equal-species projection
  # Plot Frisch-Waugh-Lovell partial latitude relationship on LOG trend scale.
  # ---------------------------------------------------------------------------

  if (make_figure) {
    font <- if (.Platform$OS.type == "windows") "Garamond" else "serif"
    if (.Platform$OS.type == "windows") {
      try(grDevices::windowsFonts(Garamond = grDevices::windowsFont("Garamond")), silent = TRUE)
    }

    pp <- list()
    for (g in c("univoltine", "multivoltine")) {
      d <- allp[allp$voltinism == g, , drop = FALSE]

      nsp <- table(d$SPECIES)
      w <- 1 / as.numeric(nsp[d$SPECIES])
      x <- d$latitude_within_species_deg / 10
      z <- d$longitude_within_species_deg / 10
      y <- d$trend_within_species_log_decade

      # Weighted FWL residualisation against longitude.
      denz <- sum(w * z^2)
      need(is.finite(denz) && denz > 0, "No longitude variation for partial plot.")
      xres <- x - z * sum(w * z * x) / denz
      yres <- y - z * sum(w * z * y) / denz

      key <- paste(g, "lat_lon", "equal_species", sep = "__")
      lat_row <- summary_tab[
        summary_tab$voltinism == g &
          summary_tab$specification == "latitude_plus_longitude_within_species" &
          summary_tab$weighting == "equal_species" &
          summary_tab$term == "latitude", , drop = FALSE
      ]
      need(nrow(lat_row) == 1L, "Preferred latitude result missing for figure.")

      b <- lat_row$estimate_log_trend_per_decade_per_10deg
      lo <- lat_row$lower_95_species_bootstrap
      hi <- lat_row$upper_95_species_bootstrap

      plotdat <- data.frame(
        latitude_partial_deg = xres * 10,
        trend_partial_log_decade = yres
      )

      xx <- seq(
        stats::quantile(plotdat$latitude_partial_deg, .005),
        stats::quantile(plotdat$latitude_partial_deg, .995),
        length.out = 250
      )
      line <- data.frame(
        x = xx,
        fit = b * xx / 10,
        low = lo * xx / 10,
        high = hi * xx / 10
      )
      line$ymin <- pmin(line$low, line$high)
      line$ymax <- pmax(line$low, line$high)

      label <- sprintf(
        "%.3f [%.3f, %.3f] per 10°",
        b, lo, hi
      )

      pp[[g]] <- ggplot2::ggplot(
        plotdat,
        ggplot2::aes(latitude_partial_deg, trend_partial_log_decade)
      ) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey45") +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey45") +
        ggplot2::geom_point(size = .55, alpha = .15, colour = "black") +
        ggplot2::geom_ribbon(
          data = line,
          ggplot2::aes(x = x, ymin = ymin, ymax = ymax),
          inherit.aes = FALSE,
          fill = "#4C78A8", alpha = .22
        ) +
        ggplot2::geom_line(
          data = line,
          ggplot2::aes(x = x, y = fit),
          inherit.aes = FALSE,
          colour = "#2F5AA8", linewidth = 1.1
        ) +
        ggplot2::annotate(
          "label",
          x = -Inf, y = Inf,
          label = label,
          hjust = -.05, vjust = 1.15,
          family = font, size = 4
        ) +
        ggplot2::labs(
          title = if (g == "univoltine") "Univoltine" else "Multivoltine",
          x = "Latitude within species (°)",
          y = "Relative abundance trend\n(log change per decade)"
        ) +
        ggplot2::theme_classic(base_family = font, base_size = 15) +
        ggplot2::theme(
          panel.border = ggplot2::element_rect(
            colour = "black", fill = NA, linewidth = .5
          ),
          plot.title = ggplot2::element_text(
            hjust = .5, face = "bold", size = 17
          )
        )
    }

    # Avoid requiring patchwork: save each and a combined horizontal figure
    # only when patchwork is already installed.
    for (g in names(pp)) {
      ggplot2::ggsave(
        file.path(out, paste0("population_trend_latitude_", g, ".png")),
        pp[[g]], width = 7, height = 5, dpi = 350, bg = "white"
      )
    }

    if (requireNamespace("patchwork", quietly = TRUE)) {
      combined <- pp$univoltine + pp$multivoltine +
        patchwork::plot_annotation(tag_levels = "a")
      ggplot2::ggsave(
        file.path(out, "population_trend_latitude.png"),
        combined, width = 13, height = 5.4, dpi = 350, bg = "white"
      )
      if (capabilities("cairo")) {
        ggplot2::ggsave(
          file.path(out, "population_trend_latitude.pdf"),
          combined, width = 13, height = 5.4,
          device = grDevices::cairo_pdf, bg = "white"
        )
      }
    }
  }

  # ---------------------------------------------------------------------------
  # README + archive
  # ---------------------------------------------------------------------------

  readme <- c(
    "phenoIMPACT | Geographic patterns in model-based abundance trends",
    "",
    paste("Source abundance run:", source_run),
    paste("Univoltine refined fit:", uni_fit_path),
    paste("Multivoltine current fit:", multi_fit_path),
    "",
    "No model was refitted.",
    "Population linear trend = fixed year slope + year:plasticity terms + species year slope + site year slope.",
    "Population AR1 and annual spatial-field realizations are year-specific deviations and are not added to the linear trend.",
    "Primary geographic specification: species-centred latitude + species-centred longitude.",
    "Both response trend and coordinates are centred within species before projection.",
    "Preferred weighting: equal total weight per species; equal-population weighting is exported as sensitivity.",
    "Intervals are species-cluster bootstrap intervals across the fitted population trends.",
    "They do NOT propagate the full joint first-stage uncertainty of random slopes/plasticity/abundance indices.",
    "A residual Moran screen is included when spdep is available. If residual spatial structure remains, consider a dedicated spatial second-stage sensitivity rather than assuming this projection is spatially independent.",
    "Trend scale for the primary projection is log abundance change per decade, matching the Gamma/log abundance model.",
    "growth_factor_percent_difference_per_10deg = 100*(exp(latitude coefficient)-1).",
    "",
    paste("Bootstrap replicates:", n_boot),
    paste("Seed:", seed)
  )
  writeLines(readme, file.path(out, "README.txt"))

  # Input provenance.
  manifest <- data.frame(
    role = c("univoltine_fit", "univoltine_input", "multivoltine_fit", "multivoltine_input"),
    path = c(uni_fit_path, uni_input_path, multi_fit_path, multi_input_path),
    bytes = file.info(c(uni_fit_path, uni_input_path, multi_fit_path, multi_input_path))$size,
    md5 = unname(tools::md5sum(c(uni_fit_path, uni_input_path, multi_fit_path, multi_input_path))),
    stringsAsFactors = FALSE
  )
  csv(manifest, file.path(out, "input_manifest.csv"))

  zip_path <- file.path(out, "results_to_review.zip")
  if (requireNamespace("zip", quietly = TRUE)) {
    rel <- list.files(out, recursive = TRUE, full.names = FALSE)
    rel <- rel[grepl("\\.(csv|txt|png|pdf)$", rel, ignore.case = TRUE)]
    rel <- setdiff(rel, basename(zip_path))
    zip::zipr(zip_path, files = rel, root = out, mode = "mirror",
      include_directories = FALSE)
  } else {
    # Base R fallback.
    old <- getwd(); on.exit(setwd(old), add = TRUE)
    setwd(out)
    rel <- list.files(".", recursive = TRUE)
    rel <- rel[grepl("\\.(csv|txt|png|pdf)$", rel, ignore.case = TRUE)]
    utils::zip(zipfile = zip_path, files = rel)
  }

  cat("\nOutput directory:\n", out, "\n", sep = "")
  cat("\nMain geographic results:\n")
  print(summary_tab[
    summary_tab$specification == "latitude_plus_longitude_within_species" &
      summary_tab$weighting == "equal_species", ], row.names = FALSE)
  cat("\nFILE TO SEND:\n", zip_path, "\n", sep = "")

  invisible(list(
    output_directory = out,
    trends = allp,
    summary = summary_tab,
    moran = moran,
    zip_to_send = zip_path
  ))
}
