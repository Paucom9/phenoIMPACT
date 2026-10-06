# phenoIMPACT | Independent geographic test of abundance trends
# 2026-10-06
#
# Direct model, separately for multivoltines and univoltines:
#   abundance ~ year * latitude_within_species + year * longitude_within_species
#   + species/site temporal slopes + population AR(1) + annual IID Matern field
#
# IMPORTANT: phenological plasticity is NOT included in the new fitted design.
# Existing reviewed abundance fits are used only as warm starts; all parameters
# are re-estimated.
#
# Run:
#   source(file.choose(), encoding = "UTF-8")
#   geo_abund <- run_independent_abundance_geography()
#
# The script runs multivoltines first, then univoltines, and saves each fit
# immediately after it is completed. Keep R open.

run_independent_abundance_geography <- function(
  project_root = "E:/phenoIMPACT project/code/phenoIMPACT",
  source_run = file.path(
    project_root, "output", "population_trends",
    "spatial_abundance_onset_adjusted_15km",
    "run_20261004_214012_44136"
  ),
  groups = c("multivoltine", "univoltine"),
  iter_max = 1000L,
  eval_max = 1500L
) {

  need <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)
  csv <- function(x, p) utils::write.csv(x, p, row.names = FALSE, na = "")
  pkgs <- c("TMB", "Matrix", "fmesher", "sf", "digest")
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  need(!length(miss), paste("Missing packages:", paste(miss, collapse = ", ")))

  project_root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
  source_run <- normalizePath(source_run, winslash = "/", mustWork = TRUE)
  need(all(groups %in% c("multivoltine", "univoltine")), "Unknown group.")

  # -------------------------------------------------------------------------
  # Current reviewed fits
  # -------------------------------------------------------------------------
  multi_fit <- file.path(source_run, "multivoltine__full", "fit.rds")
  need(file.exists(multi_fit), "Current multivoltine fit missing.")

  ref_root <- file.path(source_run, "refinement_univoltine")
  need(dir.exists(ref_root), "Univoltine refinement directory missing.")
  ref_dirs <- list.dirs(ref_root, recursive = FALSE, full.names = TRUE)
  ref_dirs <- ref_dirs[vapply(ref_dirs, function(d) {
    r <- file.path(d, "result.rds")
    f <- file.path(d, "univoltine__full", "fit.rds")
    if (!all(file.exists(c(r, f)))) return(FALSE)
    x <- tryCatch(readRDS(r), error = function(e) NULL)
    !is.null(x) && identical(x$status, "FIT_OK")
  }, logical(1))]
  need(length(ref_dirs) == 1L,
       paste("Expected one FIT_OK univoltine refinement; found", length(ref_dirs)))
  uni_fit <- file.path(ref_dirs, "univoltine__full", "fit.rds")

  fit_paths <- list(multivoltine = normalizePath(multi_fit, winslash = "/"),
                    univoltine = normalizePath(uni_fit, winslash = "/"))

  message("Source run: ", source_run)
  message("Multivoltine warm start: ", fit_paths$multivoltine)
  message("Univoltine warm start: ", fit_paths$univoltine)
  message("Plasticity is NOT included in the new fixed-effect design.")

  # -------------------------------------------------------------------------
  # Output
  # -------------------------------------------------------------------------
  out_root <- file.path(project_root, "output", "population_trends",
                        "independent_geography_spatial_15km")
  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
  out <- file.path(out_root, paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  out <- normalizePath(out, winslash = "/", mustWork = TRUE)

  # -------------------------------------------------------------------------
  # TMB engine from the reviewed abundance workflow
  # -------------------------------------------------------------------------
  dll <- TMB::dynlib(file.path(source_run, "compiled", "trend_spde"))
  need(file.exists(dll), "Reviewed trend_spde DLL missing.")
  if (!"trend_spde" %in% names(getLoadedDLLs())) dyn.load(dll)
  TMB::openmp(n = 1L, DLL = "trend_spde")
  mesh <- readRDS(file.path(source_run, "mesh.rds"))

  sparse_TMB <- function(x) methods::as(
    methods::as(methods::as(x, "dMatrix"), "generalMatrix"), "TsparseMatrix")

  add_mesh <- function(input) {
    d <- input$data
    d$use_spatial <- 1L
    d$A <- methods::as(
      fmesher::fm_basis(mesh$mesh,
                        loc = as.matrix(input$sites[c("x_km", "y_km")])),
      "CsparseMatrix"
    )
    need(all(abs(Matrix::rowSums(d$A) - 1) < 1e-8),
         "Some sites lie outside the saved mesh.")
    d$M0 <- mesh$M0; d$M1 <- mesh$M1; d$M2 <- mesh$M2
    for (nm in c("A", "M0", "M1", "M2")) d[[nm]] <- sparse_TMB(d[[nm]])
    d
  }

  make_obj <- function(data, parameters) TMB::MakeADFun(
    data = data, parameters = parameters,
    random = c("b_species", "b_site", "ar_obs", "field"),
    DLL = "trend_spde", silent = TRUE,
    inner.control = list(maxit = 1000, trace = FALSE)
  )

  # -------------------------------------------------------------------------
  # Population-level geography, centred within species
  # -------------------------------------------------------------------------
  make_geography <- function(input) {
    f <- as.data.frame(input$frame)
    for (v in c("SPECIES", "SITE_ID", "pop_id")) f[[v]] <- as.character(f[[v]])

    sites <- as.data.frame(input$sites)
    sites$SITE_ID <- as.character(sites$SITE_ID)
    need(all(c("SITE_ID", "x_km", "y_km") %in% names(sites)),
         "Site coordinates missing.")

    # Verify site ordering against TMB zero-based index.
    smap <- unique(data.frame(SITE_ID = f$SITE_ID,
                              idx = as.integer(input$data$site) + 1L))
    chk <- merge(smap,
      data.frame(idx = seq_len(nrow(sites)), SITE_ID2 = sites$SITE_ID),
      by = "idx", sort = FALSE)
    need(nrow(chk) == nrow(smap) && all(chk$SITE_ID == chk$SITE_ID2),
         "Site ordering does not match TMB indices.")

    ss <- sites
    ss$x_m <- ss$x_km * 1000; ss$y_m <- ss$y_km * 1000
    sfobj <- sf::st_as_sf(ss, coords = c("x_m", "y_m"), crs = 3035)
    ll <- sf::st_coordinates(sf::st_transform(sfobj, 4326))
    sites$longitude_deg <- ll[,1]; sites$latitude_deg <- ll[,2]

    pop <- unique(f[c("pop_id", "SPECIES", "SITE_ID")])
    pop <- merge(pop, sites[c("SITE_ID", "x_km", "y_km",
                              "longitude_deg", "latitude_deg")],
                 by = "SITE_ID", all.x = TRUE, sort = FALSE)
    need(!anyDuplicated(pop$pop_id) && !anyNA(pop$latitude_deg),
         "Population-coordinate join failed.")

    latm <- tapply(pop$latitude_deg, pop$SPECIES, mean)
    lonm <- tapply(pop$longitude_deg, pop$SPECIES, mean)
    pop$latitude_within_species_deg <- pop$latitude_deg - unname(latm[pop$SPECIES])
    pop$longitude_within_species_deg <- pop$longitude_deg - unname(lonm[pop$SPECIES])
    pop$latitude_within_species_10deg <- pop$latitude_within_species_deg / 10
    pop$longitude_within_species_10deg <- pop$longitude_within_species_deg / 10

    need(max(abs(rowsum(pop$latitude_within_species_deg, pop$SPECIES,
                        reorder = FALSE))) < 1e-8, "Latitude centring failed.")
    need(max(abs(rowsum(pop$longitude_within_species_deg, pop$SPECIES,
                        reorder = FALSE))) < 1e-8, "Longitude centring failed.")

    m <- match(f$pop_id, pop$pop_id)
    need(!anyNA(m), "Annual rows could not be matched to population geography.")
    for (v in c("latitude_deg", "longitude_deg",
                "latitude_within_species_deg", "longitude_within_species_deg",
                "latitude_within_species_10deg", "longitude_within_species_10deg"))
      f[[v]] <- pop[[v]][m]

    X <- stats::model.matrix(
      ~ year_decade * (latitude_within_species_10deg +
                       longitude_within_species_10deg), data = f)
    expected <- c("(Intercept)", "year_decade",
      "latitude_within_species_10deg", "longitude_within_species_10deg",
      "year_decade:latitude_within_species_10deg",
      "year_decade:longitude_within_species_10deg")
    need(identical(colnames(X), expected), "Unexpected geography design matrix.")
    need(qr(X)$rank == ncol(X), "Geography design is rank deficient.")

    r <- cor(pop$latitude_within_species_deg, pop$longitude_within_species_deg)
    diag <- data.frame(n_observations = nrow(f), n_populations = nrow(pop),
      n_species = length(unique(pop$SPECIES)), n_sites = length(unique(pop$SITE_ID)),
      lat_lon_correlation = r, lat_lon_VIF = 1/(1-r^2))
    list(frame = f, pop = pop, X = X, diagnostics = diag)
  }

  fixed_table <- function(beta, V) {
    se <- sqrt(diag(V)); z <- beta/se
    data.frame(term = names(beta), estimate = as.numeric(beta),
      std.error = as.numeric(se), conf.low = as.numeric(beta-1.96*se),
      conf.high = as.numeric(beta+1.96*se), z = as.numeric(z),
      p_Wald = as.numeric(2*pnorm(-abs(z))))
  }

  all_effects <- list(); all_checks <- list(); all_pops <- list()

  # -------------------------------------------------------------------------
  # Fit each group
  # -------------------------------------------------------------------------
  for (group in groups) {
    message("\n==================== ", toupper(group), " ====================")
    input <- readRDS(file.path(source_run, paste0("input_", group, ".rds")))
    oldfit <- readRDS(fit_paths[[group]])
    need(identical(input$group, group) && identical(oldfit$group, group),
         "Wrong saved group.")
    need(identical(input$data_signature, oldfit$data_signature),
         "Current fit/input signature mismatch.")
    need(isTRUE(oldfit$checks$numerical_checks_ok),
         "Current warm-start fit requires numerical review.")

    g <- make_geography(input)
    gd <- g$diagnostics; gd$group <- group
    csv(gd, file.path(out, paste0("geography_diagnostics_", group, ".csv")))
    csv(g$pop, file.path(out, paste0("population_geography_", group, ".csv")))

    # Sanity: reproduce the CURRENT reviewed source likelihood before changing X.
    source_data <- add_mesh(input)
    source_obj <- make_obj(source_data, oldfit$parameters)
    source_val <- source_obj$fn(source_obj$par)
    tol <- max(.001, abs(oldfit$optimizer$objective)*1e-8)
    replay <- data.frame(group = group,
      saved_negative_logLik = oldfit$optimizer$objective,
      replay_negative_logLik = source_val,
      absolute_difference = abs(source_val-oldfit$optimizer$objective),
      tolerance = tol,
      equivalent = is.finite(source_val) &&
        abs(source_val-oldfit$optimizer$objective) <= tol)
    csv(replay, file.path(out, paste0("source_replay_", group, ".csv")))
    TMB::FreeADFun(source_obj); rm(source_obj); gc()
    need(replay$equivalent, paste("Source likelihood replay failed for", group))

    # New independent geography model: overwrite the entire fixed-effect design.
    data <- add_mesh(input); data$X <- unname(g$X)
    pars <- oldfit$parameters
    beta0 <- rep(0, ncol(g$X)); names(beta0) <- colnames(g$X)
    beta0["(Intercept)"] <- oldfit$beta["(Intercept)"]
    beta0["year_decade"] <- oldfit$beta["year_decade"]
    pars$beta <- unname(beta0)

    obj <- make_obj(data, pars)
    message("Fitting geography model: ", nrow(g$frame), " annual rows; ",
            nrow(g$pop), " populations; ", length(unique(g$pop$SPECIES)), " species")
    started <- Sys.time()
    opt <- stats::nlminb(obj$par, obj$fn, obj$gr,
      control = list(iter.max = iter_max, eval.max = eval_max, rel.tol = 1e-10))
    if (opt$convergence != 0L) {
      message("nlminb code ", opt$convergence, "; trying one BFGS refinement.")
      bfgs <- stats::optim(opt$par, obj$fn, obj$gr, method = "BFGS",
                          control = list(maxit = 300L, reltol = 1e-10))
      if (is.finite(bfgs$value) && bfgs$value <= opt$objective + 1e-6)
        opt <- list(par = bfgs$par, objective = bfgs$value,
                    convergence = bfgs$convergence, message = "BFGS refinement")
    }
    need(is.finite(obj$fn(opt$par)), "Final geography objective is non-finite.")

    message("Computing Hessian/covariance...")
    sd <- TMB::sdreport(obj, par.fixed = opt$par,
      getJointPrecision = FALSE, getReportCovariance = FALSE)
    grad <- as.numeric(obj$gr(opt$par)); Vfull <- sd$cov.fixed
    ids <- which(names(opt$par) == "beta")
    beta <- setNames(as.numeric(opt$par[ids]), colnames(g$X))
    V <- Vfull[ids,ids,drop=FALSE]; dimnames(V) <- list(names(beta),names(beta))
    scaled <- if (all(is.finite(Vfull)))
      sqrt(max(0, as.numeric(crossprod(grad, Vfull %*% grad)))) else Inf

    obj$fn(opt$par)
    parameters <- obj$env$parList(par = obj$env$last.par)
    report <- obj$report(obj$env$last.par)
    extent <- max(diff(range(input$sites$x_km)), diff(range(input$sites$y_km)))
    good <- opt$convergence == 0L && isTRUE(sd$pdHess) && all(is.finite(V)) &&
      all(diag(V)>0) && is.finite(scaled) && scaled < .05

    checks <- data.frame(group = group, convergence_code = opt$convergence,
      pdHess = isTRUE(sd$pdHess), max_abs_gradient = max(abs(grad)),
      gradient_covariance_norm = scaled, numerical_checks_ok = good,
      logLik = -opt$objective, AIC = 2*opt$objective + 2*length(opt$par),
      nobs = nrow(g$frame), spatial_range_km = as.numeric(report$range),
      spatial_sd = as.numeric(report$spatial_sd),
      range_greater_than_1_5_extent = as.numeric(report$range)>1.5*extent,
      AR1_rho = as.numeric(report$rho), Gamma_shape = as.numeric(report$shape),
      elapsed_minutes = as.numeric(difftime(Sys.time(), started, units = "mins")))
    csv(checks, file.path(out, paste0("fit_checks_", group, ".csv")))

    fx <- fixed_table(beta, V); fx$group <- group
    csv(fx, file.path(out, paste0("fixed_effects_", group, ".csv")))

    targets <- c(latitude = "year_decade:latitude_within_species_10deg",
                 longitude = "year_decade:longitude_within_species_10deg")
    effects <- do.call(rbind, lapply(names(targets), function(direction) {
      z <- fx[fx$term == targets[[direction]],,drop=FALSE]
      data.frame(group = group, direction = direction, term = z$term,
        estimate = z$estimate, std.error = z$std.error,
        conf.low = z$conf.low, conf.high = z$conf.high, p_Wald = z$p_Wald,
        growth_factor_ratio_per_10deg = exp(z$estimate),
        lower_ratio = exp(z$conf.low), upper_ratio = exp(z$conf.high),
        percent_difference_in_decadal_growth_factor_per_10deg = 100*expm1(z$estimate),
        lower_percent = 100*expm1(z$conf.low), upper_percent = 100*expm1(z$conf.high))
    }))
    csv(effects, file.path(out, paste0("geographic_effects_", group, ".csv")))

    # Model-based population linear trends for later visualization.
    spi <- as.integer(input$data$species)+1L; sii <- as.integer(input$data$site)+1L
    pops <- unique(data.frame(pop_id = as.character(g$frame$pop_id),
      SPECIES = as.character(g$frame$SPECIES), SITE_ID = as.character(g$frame$SITE_ID),
      latitude_deg = g$frame$latitude_deg, longitude_deg = g$frame$longitude_deg,
      latitude_within_species_deg = g$frame$latitude_within_species_deg,
      longitude_within_species_deg = g$frame$longitude_within_species_deg,
      lat10 = g$frame$latitude_within_species_10deg,
      lon10 = g$frame$longitude_within_species_10deg,
      species_index = spi, site_index = sii))
    need(!anyDuplicated(pops$pop_id), "Population trend mapping not unique.")
    bsp <- as.matrix(parameters$b_species); bsi <- as.matrix(parameters$b_site)
    pops$trend_log_per_decade <- beta["year_decade"] +
      beta[targets["latitude"]]*pops$lat10 + beta[targets["longitude"]]*pops$lon10 +
      bsp[pops$species_index,2] + bsi[pops$site_index,2]
    pops$trend_percent_per_decade <- 100*expm1(pops$trend_log_per_decade)
    pops$group <- group
    csv(pops, file.path(out, paste0("population_trends_", group, ".csv")))

    saveRDS(list(beta=beta,V=V,parameters=parameters,optimizer=opt,
                 checks=checks,group=group,design=colnames(g$X),
                 plasticity_in_design=FALSE),
            file.path(out, paste0("fit_geography_", group, ".rds")))

    all_effects[[group]] <- effects; all_checks[[group]] <- checks; all_pops[[group]] <- pops
    TMB::FreeADFun(obj); rm(obj, sd, Vfull, parameters); gc()
    message(group, " finished. Latitude effect = ", signif(effects$estimate[1],4),
            " [", signif(effects$conf.low[1],4), ", ", signif(effects$conf.high[1],4), "]")
    if (!good) message("WARNING: numerical checks require review before inference.")
  }

  if (length(all_effects)) csv(do.call(rbind,all_effects), file.path(out,"geographic_effects_all.csv"))
  if (length(all_checks)) csv(do.call(rbind,all_checks), file.path(out,"fit_checks_all.csv"))
  if (length(all_pops)) csv(do.call(rbind,all_pops), file.path(out,"population_trends_all.csv"))

  writeLines(c(
    "INDEPENDENT ABUNDANCE-GEOGRAPHY TEST",
    paste("Source run:", source_run), "",
    "Phenological plasticity is NOT included in the fitted design matrix.",
    "The same population universe is retained only for direct comparability.",
    "Existing abundance fits are warm starts; all parameters are re-estimated.",
    "Primary test: year_decade:latitude_within_species_10deg, controlling longitude.",
    "Positive = northern populations have more positive/less negative trends within species.",
    "100*(exp(beta)-1) is the percent difference in the decadal growth factor per +10 degrees latitude.",
    "Same Gamma/log likelihood, species/site time slopes, population AR1 and annual IID Matern field.",
    "15 km is the mesh construction cutoff, not the fitted spatial range."
  ), file.path(out,"README.txt"))

  if (requireNamespace("zip", quietly = TRUE)) {
    rel <- list.files(out, recursive = TRUE, full.names = FALSE)
    rel <- rel[grepl("\\.(csv|txt)$", rel, ignore.case = TRUE)]
    zip::zipr(file.path(out,"results_to_review.zip"), rel, root=out,
              mode="mirror", include_directories=FALSE)
  }

  cat("\nOutput directory:\n", out, "\n", sep="")
  if (file.exists(file.path(out,"geographic_effects_all.csv"))) {
    cat("\nMain results:\n")
    print(read.csv(file.path(out,"geographic_effects_all.csv")), row.names=FALSE)
  }
  if (file.exists(file.path(out,"results_to_review.zip")))
    cat("\nFILE TO SEND:\n", file.path(out,"results_to_review.zip"), "\n", sep="")

  invisible(list(output=out,
    effects=if(length(all_effects))do.call(rbind,all_effects) else NULL,
    checks=if(length(all_checks))do.call(rbind,all_checks) else NULL))
}
