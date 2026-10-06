# ============================================================================
# phenoIMPACT | Post-hoc diagnostics of FOUR CURRENT spatial fits
# Version 1.0 -- 2026-10-05
#
# source(file.choose(), encoding = "UTF-8")
# diagnostic_job <- pheno_current_diagnostics$run(monitor = FALSE)
# pheno_current_diagnostics$status()
#
# No models are fitted, updated or optimized. No native TMB objective is called.
# Old fit files are READ ONLY. No model engine, DLL or optimizer is sourced.
# Only small diagnostic caches, tables, PNGs and logs are written to a NEW run.
#
# EXACT sources (not the latest file found):
# OFFSET: offset_with_onset_15km/run_3911d7bc87ac (both groups)
# ABUNDANCE: spatial_abundance_onset_adjusted_15km/run_20261004_214012_44136
# UNIVOLTINE abundance is the REFINED fit:
#   refinement_univoltine/run_20261005_073241_44136/univoltine__full/fit.rds
#
# Joint conditional random-effect draws:
#   H = Q_prior + Z' W Z, sparse Cholesky; no dense inverse.
#   Gaussian identity: W = 1 / sigma^2; conditional Gaussian is exact at fixed
#     estimated parameters for the verified linear-Gaussian specification.
#   Gamma log: W = shape * y / exp(eta_mode); conditional Laplace approximation.
# All fixed coefficients, dispersion/variance/range/correlation parameters stay
# fixed at their saved estimates. A solve using the conditional score CHECKS the
# stored mode but NEVER updates it.
#
# Draw 1 is predeclared primary. Draws 2:5 assess Monte Carlo sensitivity, not
# independent replicated evidence. Do not select a passing draw, average draws
# or pool their p-values. EB-mode residuals are a secondary comparator only.
#
# Spatial: site/year means -> coordinate/network/year means; standardize by
# sqrt(sum of squared observation weights). Primary Moran graph: symmetrized
# up-to-8-neighbours within 100 km, row standardized; yearly all/within networks;
# 999 permutations; bilateral equal-tail p-values; BH across years per scope
# and draw. These are exploratory residual-permutation screens, not a calibrated
# parametric bootstrap. Heterogeneous/nonexchangeable residuals affect p-values.
# Distance-band products, QQ, binned patterns, species/network/year dispersion
# and exact calendar-lag correlations are DESCRIPTIVE, without extra p-values.
#
# This does NOT propagate measurement uncertainty in abundance, onset/offset
# dates or plasticity predictors; it does not test latent-effect normality or
# give external predictive validation. No model passes all assumptions merely
# because a diagnostic p-value is above 0.05.
#
# Method reference:
# https://sdmtmb.github.io/sdmTMB/articles/residual-checking.html
# The Gaussian layout is checked against saved sdmTMB 1.1.0 data and mappings.
# Gamma system reused from check_one_sample_spatial_residuals_fixed.R (v1.0.1).
#
# Dependencies already used in this project: Matrix, fmesher, lme4, sdmTMB 1.1.0,
# callr, ps, digest; zip is optional. Nothing is installed or upgraded.
# Definitions only on source(); explicit run() launches a background controller.
# Jobs run in fresh, sequential R subprocesses. Keep the parent R session open.
# ============================================================================

pheno_current_diagnostics <- local({
  VERSION <- "current_four_diagnostics_1.0"
  ROOT <- "E:/phenoIMPACT project/code/phenoIMPACT"
  GAMMA_CPP_MD5 <- "88e15b8506c2e85a2e1e0d3839101530"
  state <- new.env(parent = emptyenv())
  check <- function(ok, text) if (!isTRUE(ok)) stop(text, call. = FALSE)
  same <- function(a, b, tol = 1e-8) {
    length(a) == length(b) && all(is.finite(a)) && all(is.finite(b)) &&
      (length(a) == 0L || max(abs(as.numeric(a) - as.numeric(b))) <= tol)
  }
  text_formula <- function(x) paste(deparse(x, width.cutoff = 500), collapse = " ")
  csv <- function(x, path) utils::write.csv(x, path, row.names = FALSE, na = "")
  atomic <- function(x, path) {
    # Only derived diagnostics/configuration, never an input or a fitted model.
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    tmp <- tempfile(".diagnostic_", tmpdir = dirname(path), fileext = ".rds")
    on.exit(unlink(tmp), add = TRUE)
    saveRDS(x, tmp, compress = "gzip")
    bak <- paste0(path, ".previous")
    if (file.exists(path)) {
      if (file.exists(bak)) check(unlink(bak)==0L, paste("Cannot remove old diagnostic backup:",bak))
      check(file.rename(path, bak), paste("Cannot preserve previous diagnostic:", path))
    }
    if (!file.rename(tmp, path)) {
      if (file.exists(bak)) file.rename(bak, path)
      stop("Cannot finalize diagnostic cache: ", path, call. = FALSE)
    }
    if (file.exists(bak)) unlink(bak)
    invisible(path)
  }
  phase <- function(cfg, stage, job = "", detail = "") {
    p <- list(stage = stage, job = job, detail = detail, time = Sys.time())
    atomic(p, file.path(cfg$out, "progress.rds"))
    message("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", job, " | ", stage, " | ", detail)
  }
  resource_check <- function(path, disk_GiB = 15, ram_GiB = 10) {
    disk <- ps::ps_disk_usage(path)$available[1] / 1024^3
    ram <- ps::ps_system_memory()$avail / 1024^3
    check(is.finite(disk) && disk >= disk_GiB,
      paste("Only", round(disk,1), "GiB disk free; need", disk_GiB, "before diagnostics."))
    check(is.finite(ram) && ram >= ram_GiB,
      paste("Only", round(ram,1), "GiB RAM available; close other large jobs first."))
    # This is a guard, NOT an estimate of peak sparse-Cholesky memory.
    invisible(c(disk_GiB = disk, RAM_GiB = ram))
  }
  packages <- function() {
    required <- c("Matrix", "fmesher", "lme4", "sdmTMB", "callr", "ps", "digest")
    missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
    check(!length(missing), paste("Missing packages:", paste(missing, collapse = ", "),
      "-- use the original project library. No packages were installed."))
    check(as.character(utils::packageVersion("sdmTMB")) == "1.1.0",
      "Use the original sdmTMB 1.1.0 library, not the separate error-sensitivity build.")
    data.frame(package = required, version = vapply(required,
      function(x) as.character(utils::packageVersion(x)), character(1)))
  }
  verify_frame <- function(frame) {
    keys <- c("SPECIES", "SITE_ID", "YEAR", "bms_id", "x_km", "y_km")
    check(all(keys %in% names(frame)) && nrow(frame)>0 && !anyNA(frame[keys]),
      "Missing observation identifiers, years or coordinates.")
    check(!anyDuplicated(frame[c("SPECIES","SITE_ID","YEAR")]),
      "Duplicate population/year observations; nothing was averaged or dropped.")
    check(is.numeric(frame$YEAR) && all(is.finite(frame$YEAR)) &&
      all(frame$YEAR == floor(frame$YEAR)), "YEAR must be integer calendar years.")
    check(all(is.finite(frame$x_km)) && all(is.finite(frame$y_km)), "Non-finite coordinates.")
    for (key in c("SPECIES","SITE_ID","bms_id")) frame[[key]] <- as.character(frame[[key]])
    frame
  }
  precision_transform <- function(ch, e) {
    a <- Matrix::solve(ch, e, system = "Lt")
    Matrix::solve(ch, a, system = "Pt")
  }

  sparse_triplets <- function(x) {
    check(methods::is(x, "sparseMatrix"), "Projection must be a sparse Matrix object.")
    # General form expands symmetric/triangular storage if ever encountered;
    # dMatrix ensures that an explicit numeric x slot exists.
    tr <- methods::as(methods::as(methods::as(x, "generalMatrix"),
                                  "dMatrix"), "TsparseMatrix")
    ir <- methods::slot(tr, "i") + 1L
    jc <- methods::slot(tr, "j") + 1L
    val <- methods::slot(tr, "x")
    check(length(ir) == length(jc) && length(ir) == length(val) &&
            all(ir >= 1L & ir <= nrow(x)) &&
            all(jc >= 1L & jc <= ncol(x)) &&
            is.numeric(val) && all(is.finite(val)),
          "Invalid sparse projection indices or values.")
    list(i = ir, j = jc, x = val)
  }

  make_aggregation <- function(frame) {
    cols <- c("SITE_ID", "YEAR", "bms_id", "x_km", "y_km")
    check(all(cols %in% names(frame)) && !anyNA(frame[cols]),
      "Missing identifiers/coordinates in input frame.")
    d <- frame[cols]
    d$.row_index <- seq_len(nrow(d))
    a1 <- stats::aggregate(.row_index ~ SITE_ID + YEAR + bms_id + x_km + y_km,
                           d, FUN = min)
    a1$.sy <- seq_len(nrow(a1))
    joined <- merge(d, a1[c(cols, ".sy")], by = cols, sort = FALSE)
    joined <- joined[order(joined$.row_index), , drop = FALSE]
    check(identical(as.integer(joined$.row_index), seq_len(nrow(d))),
      "Observation-to-site/year aggregation is not one-to-one.")
    sy <- as.integer(joined$.sy)
    count_sy <- tabulate(sy, nbins = nrow(a1))
    loc_cols <- c("YEAR", "bms_id", "x_km", "y_km")
    a2 <- stats::aggregate(.sy ~ YEAR + bms_id + x_km + y_km, a1, FUN = min)
    a2$.loc <- seq_len(nrow(a2))
    site_map <- merge(a1[c(loc_cols, ".sy")], a2[c(loc_cols, ".loc")],
                      by = loc_cols, sort = FALSE)
    site_map <- site_map[order(site_map$.sy), , drop = FALSE]
    check(identical(as.integer(site_map$.sy), seq_len(nrow(a1))),
      "Site-to-location aggregation is not one-to-one.")
    sites_per_loc <- tabulate(site_map$.loc, nbins = nrow(a2))
    loc <- as.integer(site_map$.loc[sy])
    weights <- 1 / count_sy[sy] / sites_per_loc[loc]
    B <- Matrix::sparseMatrix(i = loc, j = seq_len(nrow(d)), x = weights,
                              dims = c(nrow(a2), nrow(d)))
    check(max(abs(Matrix::rowSums(B) - 1)) < 1e-10, "Aggregation weights do not sum to one.")
    variance <- as.numeric(Matrix::rowSums(B * B))
    check(all(is.finite(variance) & variance > 0), "Invalid aggregation variance.")
    meta <- a2[loc_cols]
    meta$n_observations <- tabulate(loc, nbins = nrow(a2))
    meta$n_sites <- sites_per_loc
    meta$null_sd_of_mean <- sqrt(variance)
    list(B = B, meta = meta)
  }

  gamma_system_raw <- function(input, pars, spatial, mesh = NULL, A = NULL) {
    d <- input$data
    y <- as.numeric(d$y); n <- length(y)
    ns <- nrow(pars$b_species); nt <- nrow(pars$b_site)
    check(n == nrow(input$frame) && all(is.finite(y) & y > 0),
      "The saved model requires finite, positive Gamma observations.")
    check(is.matrix(pars$b_species) && ncol(pars$b_species) == 2L &&
            is.matrix(pars$b_site) && ncol(pars$b_site) == 2L,
      "Expected two independent random coefficients for species and sites.")
    check(length(pars$log_sd) == 5L && length(pars$ar_obs) == n &&
            ncol(d$X) == length(pars$beta), "Unexpected parameter dimensions.")
    for (nm in c("time", "species", "site", "year", "previous", "gap"))
      check(length(d[[nm]]) == n && all(is.finite(d[[nm]])), paste("Invalid", nm))
    check(all(d$species >= 0 & d$species < ns & d$species == floor(d$species)) &&
            all(d$site >= 0 & d$site < nt & d$site == floor(d$site)),
      "Invalid zero-based species/site indices.")
    check(all(d$gap >= 1 & d$gap == floor(d$gap)), "AR(1) gaps must be positive integers.")
    sd_re <- exp(as.numeric(pars$log_sd))
    rho <- as.numeric(pars$rho_raw / sqrt(1 + pars$rho_raw^2))
    shape <- as.numeric(exp(pars$log_shape))
    check(all(is.finite(sd_re) & sd_re > 0) && is.finite(rho) && abs(rho) < 1 &&
            length(shape) == 1L && is.finite(shape) && shape > 0, "Invalid variance/shape parameters.")
    ii <- seq_len(n); prev <- as.integer(d$previous) + 1L
    nonfirst <- which(prev > 0L)
    check(all(d$previous >= -1 & d$previous == floor(d$previous)) &&
            all(prev[nonfirst] < nonfirst), "Invalid AR(1) previous-observation indices.")
    if (length(nonfirst)) {
      check(all(input$frame$pop_id[nonfirst] == input$frame$pop_id[prev[nonfirst]]) &&
        all(input$frame$YEAR[nonfirst] - input$frame$YEAR[prev[nonfirst]] == d$gap[nonfirst]),
        "AR(1) links do not match the saved population/year frame.")
    }
    r <- rho^d$gap[nonfirst]
    invv <- rep(1 / sd_re[5]^2, n)
    invv[nonfirst] <- 1 / (sd_re[5]^2 * (1 - r^2))
    check(all(is.finite(invv) & invv > 0), "Invalid AR(1) innovation precision.")
    Q_ar <- Matrix::sparseMatrix(
      i = c(ii, prev[nonfirst], nonfirst, prev[nonfirst]),
      j = c(ii, prev[nonfirst], prev[nonfirst], nonfirst),
      x = c(invv, r^2 * invv[nonfirst], -r * invv[nonfirst], -r * invv[nonfirst]),
      dims = c(n, n))
    Q_coef <- Matrix::Diagonal(x = c(rep(1 / sd_re[1]^2, ns), rep(1 / sd_re[2]^2, ns),
                                     rep(1 / sd_re[3]^2, nt), rep(1 / sd_re[4]^2, nt)))
    ar_start <- 2L * ns + 2L * nt
    field_start <- ar_start + n
    u <- c(as.vector(pars$b_species), as.vector(pars$b_site), as.vector(pars$ar_obs))
    zi <- rep(ii, 5L)
    zj <- c(d$species + 1L, ns + d$species + 1L,
            2L * ns + d$site + 1L, 2L * ns + nt + d$site + 1L, ar_start + ii)
    zx <- c(rep(1, n), d$time, rep(1, n), d$time, rep(1, n))
    prior_parts <- list(Q_coef, Q_ar)
    if (spatial) {
      check(!is.null(mesh) && !is.null(A), "Spatial model requires saved mesh and projection.")
      nv <- nrow(mesh$M0); ny <- length(input$years)
      check(identical(dim(pars$field), as.integer(c(nv, ny))) &&
              nrow(A) == nt && ncol(A) == nv, "Saved field or projection dimensions disagree.")
      check(all(d$year >= 0 & d$year < ny & d$year == floor(d$year)) &&
              all(input$years[d$year + 1L] == input$frame$YEAR), "Incorrect field year mapping.")
      check(max(abs(Matrix::rowSums(A) - 1)) < 1e-8, "Some site coordinates lie outside the mesh.")
      kappa <- sqrt(8) / exp(pars$log_range)
      tau <- 1 / (sqrt(4 * pi) * kappa * exp(pars$log_spatial_sd))
      Q_one_field <- tau^2 * (kappa^4 * mesh$M0 + 2 * kappa^2 * mesh$M1 + mesh$M2)
      Q_field <- kronecker(Matrix::Diagonal(ny), Q_one_field)
      prior_parts[[3L]] <- Q_field
      u <- c(u, as.vector(pars$field))
      aa <- sparse_triplets(A[d$site + 1L, , drop = FALSE])
      zi <- c(zi, aa$i)
      zj <- c(zj, field_start + aa$j + nv * d$year[aa$i])
      zx <- c(zx, aa$x)
    }
    check(all(is.finite(u)) && all(is.finite(zx)), "Non-finite latent values/design.")
    Z <- Matrix::drop0(Matrix::sparseMatrix(i = zi, j = zj, x = zx, dims = c(n, length(u))))
    Q <- Matrix::forceSymmetric(do.call(Matrix::bdiag, prior_parts), uplo = "L")
    fixed_eta <- as.numeric(d$X %*% pars$beta)
    eta <- fixed_eta + as.numeric(Z %*% u)
    check(all(is.finite(eta)), "Non-finite reconstructed predictions.")
    w <- shape * exp(log(y) - eta)
    check(all(is.finite(w) & w > 0), "Non-finite Gamma curvature.")
    WZ <- Matrix::Diagonal(x = sqrt(w)) %*% Z
    H <- Matrix::forceSymmetric(Q + Matrix::crossprod(WZ), uplo = "L")
    gradient <- as.numeric(Q %*% u + Matrix::crossprod(Z, shape - w))
    list(Z = Z, H = H, gradient = gradient, eta = eta, u = u,
         shape = shape, rho = rho, n_random = length(u), fixed_eta = fixed_eta)
  }

  gamma_residuals <- function(y, eta, shape) {
    # Stable tails; no clipping of extreme residuals and no dropped observations.
    scale <- exp(eta) / shape
    check(all(is.finite(scale) & scale > 0), "Overflow/underflow in Gamma means.")
    lo <- stats::pgamma(y, shape = shape, scale = scale, log.p = TRUE)
    hi <- stats::pgamma(y, shape = shape, scale = scale, lower.tail = FALSE, log.p = TRUE)
    z <- numeric(length(y)); left <- lo < log(0.5)
    z[left] <- stats::qnorm(lo[left], log.p = TRUE)
    z[!left] <- stats::qnorm(hi[!left], lower.tail = FALSE, log.p = TRUE)
    check(all(is.finite(z)), "Non-finite PIT residuals. Nothing was silently clipped or omitted.")
    z
  }
  # Gaussian model: use the stored model matrices and parameter-block mappings.
  # No calls to predict(), residuals(), tidy(), sdreport(), or native TMB methods.
  gaussian_system_raw <- function(fit, input, pars) {
    d <- verify_frame(input$data); n <- nrow(d); td <- fit$tmb_data
    check(inherits(fit, "sdmTMB") && identical(fit$family$family,"gaussian") &&
      identical(fit$family$link,"identity") && !isTRUE(fit$reml),
      "Offset diagnostic requires the saved Gaussian identity-link ML model.")
    check(length(fit$split_formula)==1L, "Only one Gaussian component is supported.")
    sf <- fit$split_formula[[1]]
    want <- OFFSET_mean ~ ONSET_mean_z + clim_anomaly_tw90 *
      (photo_tw90 + clim_background_tw90 + clim_predictability_tw90 + clim_trend_tw90)
    check(setequal(attr(stats::terms(sf$form_no_bars),"term.labels"),
      attr(stats::terms(want),"term.labels")), "Unexpected onset-adjusted offset formula.")
    check(isTRUE(all.equal(fit$data[,names(input$data),drop=FALSE],input$data)),
      "Offset fit/data disagree; no rows were removed or reordered.")
    check(same(d$ONSET_mean_z,(d$ONSET_mean-input$onset_scaling$centre)/input$onset_scaling$SD),
      "Onset covariate/scaling does not match the saved data.")
    check(isTRUE(input$settings$mesh_cutoff_km==15) &&
      isTRUE(input$settings$best_window==90), "Unexpected offset mesh/window.")
    check(isTRUE(fit$model$convergence==0L) && isTRUE(fit$sd_report$pdHess) &&
      length(fit$gradients)>0 && all(is.finite(fit$gradients)) &&
      max(abs(fit$gradients))<0.001, "Offset fit fails saved numerical checks.")
    check(same(fit$model$par,fit$sd_report$par.fixed), "Stale offset covariance report.")
    cn <- sf$re_cov_terms$cnms
    check(length(cn)==4L && setequal(names(cn),
      c("SITE_ID","SPECIES","SPECIES_slope","site_year_id")) &&
      identical(cn$SPECIES_slope,"clim_anomaly_tw90") &&
      all(vapply(cn[setdiff(names(cn),"SPECIES_slope")],
        function(x) identical(x,"(Intercept)"), logical(1))),
      "Unexpected Gaussian random coefficients; refusing to guess their covariance.")
    required_flags <- c("anisotropy","barrier","spatial_only","ar1_fields","rw_fields",
      "spatial_covariate","has_smooths","random_walk","ar1_time","est_epsilon_model",
      "include_spatial","no_spatial")
    check(all(required_flags %in% names(td)), "Missing sdmTMB structure flags.")
    check(all(vapply(td[required_flags], function(x) length(x)>0 && all(x==0),logical(1))),
      "Expected only isotropic IID annual fields, no static field or other latent terms.")
    check(length(td$weights_i)==n && length(td$offset_i)==n &&
      all(td$weights_i==1) && all(td$offset_i==0),
      "Non-unit likelihood weights or an offset require a different diagnostic.")
    y <- as.numeric(td$y_i)
    check(length(y)==n && same(y,d$OFFSET_mean), "Gaussian response mismatch.")
    check(length(pars$ln_phi)==1L && is.finite(pars$ln_phi), "Invalid Gaussian dispersion.")
    sigma <- exp(as.numeric(pars$ln_phi))
    X <- stats::model.matrix(sf$form_no_bars,d)
    check(length(td$X_ij)==1L && isTRUE(all.equal(unname(as.matrix(td$X_ij[[1]])),
      unname(X),check.attributes=FALSE,tolerance=1e-10)), "Saved Gaussian fixed design mismatch.")
    ix <- which(names(fit$sd_report$par.fixed)=="b_j")
    check(length(ix)==ncol(X) && same(pars$b_j,fit$sd_report$par.fixed[ix],1e-6),
      "Portable Gaussian beta does not match the saved fit report.")
    for (nm in c("b_j","re_cov_pars","re_b_pars","epsilon_st","ln_kappa","ln_tau_E","ln_phi")) {
      check(!is.null(pars[[nm]]) && same(pars[[nm]],fit$parlist[[nm]],1e-6),
        paste("Portable parameter file differs from synchronized fit$parlist:",nm))
    }
    b <- sf$re_cov_terms$re_b_df
    check(all(c("start","end","group_indices","level") %in% names(b)) &&
      all(b$start==b$end), "Expected independent scalar random-effect blocks.")
    nre <- length(pars$re_b_pars)
    active <- !is.na(as.vector(fit$tmb_map$re_b_pars))
    check(length(active)==nre && all(active),
      "Unexpected inactive Gaussian random-effect slots.")
    check(identical(sort(as.integer(b$start+1L)),seq_len(nre)),
      "Random-effect blocks do not cover the portable parameter vector exactly.")
    check(length(td$Zt_list)==1L &&
      identical(dim(td$Zt_list[[1]]),as.integer(c(nre,n))), "Saved Zt dimensions disagree.")
    Zre <- Matrix::t(td$Zt_list[[1]])
    # Independently rebuild its scalar design using labels, not assumed sorting.
    zi <- zj <- integer(); zx <- numeric()
    qre <- rep(NA_real_,nre)
    map <- as.matrix(td$re_cov_df_map)
    check(nrow(map)==4L && ncol(map)>=4L && all(map[,2]==1L) &&
      all(map[,3]==map[,4]), "Unexpected random-effect variance parameter map.")
    for (g in seq_along(cn)) {
      name <- names(cn)[g]; blocks <- b[b$group_indices==g,,drop=FALSE]
      check(nrow(blocks)>0 && !anyDuplicated(as.character(blocks$level)), "Invalid group labels.")
      index <- match(as.character(d[[name]]),as.character(blocks$level))
      check(!anyNA(index),paste("Unmatched saved random-effect levels:",name))
      cols <- as.integer(blocks$start[index]+1L)
      vals <- if (name=="SPECIES_slope") d$clim_anomaly_tw90 else rep(1,n)
      zi <- c(zi,seq_len(n)); zj <- c(zj,cols); zx <- c(zx,vals)
      sd_index <- as.integer(map[g,3]+1L)
      check(sd_index>=1L && sd_index<=length(pars$re_cov_pars), "SD index is outside parameter vector.")
      qre[blocks$start+1L] <- exp(-2*as.numeric(pars$re_cov_pars[sd_index]))
    }
    check(all(is.finite(qre)&qre>0), "Missing or invalid Gaussian random-effect precision.")
    Zmanual <- Matrix::sparseMatrix(i=zi,j=zj,x=zx,dims=c(n,nre))
    design_error <- max(abs(Zmanual-Zre))
    check(is.finite(design_error) && design_error<1e-10,
      "Saved Gaussian random design differs from label-based reconstruction.")
    rm(Zmanual,zi,zj,zx)
    ny <- as.integer(td$n_t); year <- as.integer(td$year_i)
    eps <- pars$epsilon_st; dims <- dim(eps)
    check(length(dims)==3L && dims[2]==ny && dims[3]==1L,
      "Unexpected epsilon_st dimensions (vertices x years x one component required).")
    nv <- dims[1]; yy <- sort(unique(d$YEAR))
    check(length(yy)==ny && length(year)==n && all(year>=0 & year<ny) &&
      all(yy[year+1L]==d$YEAR), "Gaussian field-year mapping does not match calendar years.")
    A <- td$A_st
    si <- as.integer(td$A_spatial_index)+1L
    check(ncol(A)==nv && length(si)==n && all(si>=1L & si<=nrow(A)),
      "Invalid saved Gaussian observation projection.")
    AO <- A[si,,drop=FALSE]
    check(max(abs(Matrix::rowSums(AO)-1))<1e-8, "Gaussian projection rows do not sum to one.")
    entries <- sparse_triplets(AO)
    Zf <- Matrix::sparseMatrix(i=entries$i,j=entries$j+nv*year[entries$i],
      x=entries$x,dims=c(n,nv*ny))
    check(length(pars$ln_kappa)==2L && length(pars$ln_tau_E)==1L,
      "Unexpected Matern parameter dimensions.")
    kappa <- exp(as.numeric(pars$ln_kappa)[2]); tau <- exp(as.numeric(pars$ln_tau_E))
    M <- td$spde
    check(all(c("M0","M1","M2") %in% names(M)) &&
      all(vapply(M[c("M0","M1","M2")], function(x) identical(dim(x),c(nv,nv)),logical(1))),
      "Saved Gaussian SPDE matrices are missing or incompatible.")
    Qf <- tau^2*(kappa^4*M$M0+2*kappa^2*M$M1+M$M2)
    Q <- Matrix::forceSymmetric(Matrix::bdiag(Matrix::Diagonal(x=qre),
      kronecker(Matrix::Diagonal(ny),Qf)),uplo="L")
    Z <- methods::as(cbind(Zre,Zf),"CsparseMatrix")
    u <- c(as.numeric(pars$re_b_pars),as.numeric(eps))
    eta_fixed <- as.numeric(X%*%pars$b_j)
    eta <- eta_fixed+as.numeric(Z%*%u)
    # Independent projector multiplication checks the vectorized field indexing.
    field_direct <- as.matrix(A %*% matrix(as.numeric(eps),nv,ny))
    eta_direct <- eta_fixed+as.numeric(Zre%*%as.numeric(pars$re_b_pars))+
      field_direct[cbind(si,year+1L)]
    eta_error <- max(abs(eta-eta_direct))
    check(all(is.finite(eta)) && eta_error<1e-8, "Gaussian predictor reconstruction failed.")
    W <- rep(1/sigma^2,n)
    H <- Matrix::forceSymmetric(Q+Matrix::crossprod(Z)/sigma^2,uplo="L")
    gradient <- as.numeric(Q%*%u+Matrix::crossprod(Z,(eta-y)/sigma^2))
    list(Z=Z,H=H,gradient=gradient,eta=eta,u=u,fixed_eta=eta_fixed,
      y=y,frame=d,sigma=sigma,shape=NA_real_,rho=NA_real_,n_random=length(u),
      family="gaussian",prediction_error=eta_error,design_error=design_error,
      range_km=sqrt(8)/kappa,field_sd=1/(sqrt(4*pi)*kappa*tau))
  }

  load_system <- function(job, cfg) {
    phase(cfg,"READING_SAVED_FIT",job$key,"Read-only; no fitting or native TMB evaluation")
    input <- readRDS(job$input)
    saved <- readRDS(job$fit)
    if (job$kind=="offset") {
      done <- readRDS(job$complete)
      check(identical(input$signature,saved$signature) &&
        identical(done$signature,saved$signature) && identical(done$status,"FIT_OK") &&
        identical(input$analysis,job$key) && identical(saved$analysis,job$key),
        "Offset input/fit/completion signatures disagree.")
      pars <- readRDS(job$parameters)
      phase(cfg,"RECONSTRUCTING",job$key,"Gaussian conditional precision from stored matrices")
      ans <- gaussian_system_raw(saved$fit,input,pars)
    } else {
      check(identical(input$signature,saved$signature) &&
        identical(input$data_signature,saved$data_signature) &&
        identical(saved$group,job$group) && identical(saved$variant,"full") &&
        isTRUE(saved$checks$numerical_checks_ok) && isTRUE(saved$offset_adjusted_for_onset),
        "Abundance fit is not the numerically accepted onset-adjusted model for these inputs.")
      check(identical(input$data_signature,digest::digest(list(input$frame,input$data$X),algo="sha256")),
        "Abundance data/X signature has changed.")
      if (job$group=="univoltine") check(!is.null(saved$refinement_signature) &&
        identical(saved$refinement_reason,"SCALED_GRADIENT_CRITERION_MET"),
        "The univoltine abundance source must be the accepted REFINED fit.")
      check(identical(unname(tools::md5sum(job$cpp)),GAMMA_CPP_MD5),
        "Abundance C++ template differs from the validated analytic diagnostic model.")
      mesh <- readRDS(job$mesh)
      A <- methods::as(fmesher::fm_basis(mesh$mesh,
        loc=as.matrix(input$sites[c("x_km","y_km")])),"CsparseMatrix")
      check(all(as.character(input$sites$SITE_ID[input$data$site+1L])==
        as.character(input$frame$SITE_ID)),"Abundance site projection ordering changed.")
      phase(cfg,"RECONSTRUCTING",job$key,"Gamma conditional precision; same saved 15-km mesh")
      ans <- gamma_system_raw(input,saved$parameters,spatial=TRUE,mesh=mesh,A=A)
      ans$prediction_error <- max(abs(ans$eta-saved$eta))
      ans$design_error <- 0
      check(length(saved$eta)==length(ans$eta) && ans$prediction_error<1e-6,
        "Gamma linear predictions do not reproduce the saved fit.")
      check(same(saved$beta,saved$parameters$beta,1e-8),"Gamma fixed coefficients are inconsistent.")
      ans$frame <- verify_frame(input$frame)
      ans$y <- as.numeric(input$data$y)
      ans$family <- "Gamma"; ans$sigma <- NA_real_
      ans$range_km <- exp(saved$parameters$log_range)
      ans$field_sd <- exp(saved$parameters$log_spatial_sd)
    }
    # Return no fitted-model environments, optimizer closures or dense covariance.
    rm(input,saved); invisible(gc())
    ans
  }

  self_tests <- function() {
    Q <- Matrix::Matrix(matrix(c(5,1,0,1,1,4,1,0,0,1,3,1,1,0,1,4),4,4),sparse=TRUE)
    ch <- Matrix::Cholesky(Matrix::forceSymmetric(Q),LDL=FALSE,super=TRUE)
    root <- as.matrix(precision_transform(ch,diag(4)))
    e <- max(abs(tcrossprod(root)-solve(as.matrix(Q))))
    check(is.finite(e) && e<1e-10,"Precision-root self-test failed.")
    small <- Matrix::sparseMatrix(i=c(1L,1L,2L,3L),j=c(1L,3L,2L,4L),
      x=c(.4,.6,1,1),dims=c(3L,4L))
    for (a in list(small,small[1,,drop=FALSE],small[,1,drop=FALSE],small*0)) {
      t <- sparse_triplets(a)
      b <- Matrix::sparseMatrix(i=t$i,j=t$j,x=t$x,dims=dim(a))
      check(max(abs(b-a))<1e-12,"Sparse triplet self-test failed.")
    }
    # Directional derivative of Gamma score versus analytic Hessian.
    Z <- matrix(c(1,0,1,1,0,1,1,-1,1,2,0,1),4,3,byrow=TRUE)
    prior <- diag(c(.5,2,3)); u <- c(.2,-.1,.3); b <- c(.5,.7,1,.2)
    y <- c(1.2,.8,2.1,1.5); shape <- 2.3; v <- c(.6,-.3,.2)
    score <- function(u) as.vector(prior%*%u+t(Z)%*%(shape-shape*y/exp(b+Z%*%u)))
    eta <- as.vector(b+Z%*%u); w <- shape*y/exp(eta)
    H <- prior+t(Z)%*%diag(w)%*%Z
    diff <- (score(u+1e-5*v)-score(u-1e-5*v))/(2e-5)
    eh <- max(abs(diff-H%*%v))
    check(eh<1e-7,"Gamma analytic Hessian self-test failed.")
    # Gaussian conditional moments and score at a directly solved toy mean.
    sigma <- 1.4; Hn <- prior+crossprod(Z)/sigma^2
    mu <- solve(Hn,crossprod(Z,y-b)/sigma^2)
    check(max(abs(prior%*%mu+crossprod(Z,b+Z%*%mu-y)/sigma^2))<1e-10,
      "Gaussian conditional mean/precision self-test failed.")
    a <- data.frame(SITE_ID=c("a","a","b","c"),YEAR=2000L,bms_id="n",
      x_km=c(1,1,1,2),y_km=0)
    ag <- make_aggregation(a)
    means <- as.numeric(ag$B%*%c(1,3,5,7))
    check(same(sort(means),c(3.5,7)) &&
      same(sort(ag$meta$null_sd_of_mean^2),c(.375,1)),"Aggregation-weight self-test failed.")
    # Explicit singleton final permutation batch regression test.
    P <- matrix(seq_len(20),nrow=20,ncol=1)
    WP <- as.matrix(Matrix::Diagonal(20)%*%P); dim(WP)<-c(20,1)
    check(length(colSums(P*WP))==1,"Singleton permutation-batch self-test failed.")
    data.frame(test=c("precision_square_root","gamma_Hessian_derivative","Gaussian_score"),
      maximum_error=c(e,eh,max(abs(prior%*%mu+crossprod(Z,b+Z%*%mu-y)/sigma^2))))
  }

  make_graphs <- function(meta) {
    graphs <- list()
    for (yr in sort(unique(meta$YEAR))) {
      rows <- which(meta$YEAR == yr)
      a <- meta[rows, , drop = FALSE]; n <- nrow(a)
      D <- as.matrix(stats::dist(as.matrix(a[c("x_km", "y_km")])))
      diag(D) <- Inf
      for (scope in c("all_networks", "within_networks")) {
        DD <- D
        if (scope == "within_networks") DD[outer(a$bms_id, a$bms_id, "!=")] <- Inf
        edges <- lapply(seq_len(n), function(i) {
          ids <- which(DD[i, ] > 0 & DD[i, ] <= 100)
          head(ids[order(DD[i, ids])], 8L)
        })
        jj <- as.integer(unlist(edges, use.names = FALSE))
        W <- Matrix::sparseMatrix(i = rep(seq_len(n), lengths(edges)), j = jj,
                                 x = rep(1, length(jj)), dims = c(n, n))
        W <- methods::as((W + Matrix::t(W)) > 0, "dMatrix")
        keep <- Matrix::rowSums(W) > 0
        W <- W[keep, keep, drop = FALSE]
        if (any(keep)) W <- Matrix::Diagonal(x = 1 / Matrix::rowSums(W)) %*% W
        bms <- a$bms_id[keep]
        blocks <- if (scope == "within_networks") {
          split(seq_along(bms), bms)
        } else {
          list(seq_along(bms))
        }
        graphs[[paste(yr, scope, sep = "__")]] <- list(
          YEAR = yr, scope = scope, rows = rows[keep], bms = bms,
          W = W, blocks = blocks, n_excluded = sum(!keep))
      }
    }
    graphs
  }


  moran_tests <- function(tab, graphs, job, draw, cfg) {
    group <- job$group; model <- job$key
    seed <- cfg$permutation_seed; nperm <- cfg$nperm
    ans <- list()
    for (aggregation in c("variance_standardized", "original_mean")) {
      # Reset identically for every model/draw. Draw generation has its own seeds.
      set.seed(seed + 800000L)
      zall <- if (aggregation == "variance_standardized") tab$z_standardized else tab$z_mean
      for (g in graphs) {
        v <- zall[g$rows]; nn <- length(v)
        row <- data.frame(group = group, model = model,
          residual_type = if (draw == 0L) "EB_mode" else "one_sample", draw = draw,
          aggregation = aggregation, YEAR = g$YEAR, scope = g$scope,
          n_connected = nn, n_excluded = g$n_excluded,
          I = NA_real_, null_mean = NA_real_, null_sd = NA_real_,
          I_minus_null_mean = NA_real_, p_positive = NA_real_, p_negative = NA_real_,
          p_two_sided = NA_real_, nperm = nperm, status = "too_few_connected",
          stringsAsFactors = FALSE)
        if (nn >= 15L) {
          if (g$scope == "within_networks") v <- v - ave(v, g$bms, FUN = mean)
          v <- v - mean(v); den <- sum(v^2)
          row$status <- "zero_variance"
          if (is.finite(den) && den > 1e-20) {
            I <- as.numeric(crossprod(v, g$W %*% v) / den)
            sims <- numeric(nperm)
            for (start in seq.int(1L, nperm, by = 100L)) {
              m <- min(100L, nperm - start + 1L)
              P <- matrix(v, nrow = nn, ncol = m)
              for (j in seq_len(m)) for (ind in g$blocks)
                P[ind, j] <- v[ind[sample.int(length(ind))]]
              WP <- as.matrix(g$W %*% P)
              dim(WP) <- c(nn, m)
              sims[start + seq_len(m) - 1L] <- colSums(P * WP) / den
            }
            p_hi <- (1 + sum(sims >= I)) / (nperm + 1)
            p_lo <- (1 + sum(sims <= I)) / (nperm + 1)
            row$I <- I; row$null_mean <- mean(sims); row$null_sd <- stats::sd(sims)
            row$I_minus_null_mean <- I - row$null_mean
            row$p_positive <- p_hi; row$p_negative <- p_lo
            row$p_two_sided <- min(1, 2 * min(p_hi, p_lo)); row$status <- "ok"
          }
        }
        ans[[length(ans) + 1L]] <- row
      }
    }
    out <- do.call(rbind, ans)
    # Separate BH family across years for each group/model/draw/aggregation/scope.
    family <- interaction(out$aggregation, out$scope, drop = TRUE)
    for (tail in c("positive", "negative", "two_sided"))
      out[[paste0("q_BH_", tail)]] <- ave(out[[paste0("p_", tail)]], family,
        FUN = function(p) stats::p.adjust(p, method = "BH"))
    out
  }


  summarize_tests <- function(tests) {
    med <- function(x) if (any(is.finite(x))) stats::median(x[is.finite(x)]) else NA_real_
    keys <- c("group", "model", "residual_type", "draw", "aggregation", "scope")
    pieces <- split(seq_len(nrow(tests)), interaction(tests[keys], drop = TRUE, lex.order = TRUE))
    rows <- lapply(pieces, function(ind) {
      a <- tests[ind, , drop = FALSE]
      cbind(a[1L, keys, drop = FALSE], data.frame(
        years_tested = sum(a$status == "ok"), years_skipped = sum(a$status != "ok"),
        median_I = med(a$I), median_null_mean = med(a$null_mean),
        median_I_minus_null_mean = med(a$I_minus_null_mean),
        n_sig_BH_positive = sum(a$q_BH_positive < 0.05, na.rm = TRUE),
        n_sig_BH_negative = sum(a$q_BH_negative < 0.05, na.rm = TRUE),
        n_sig_BH_two_sided = sum(a$q_BH_two_sided < 0.05, na.rm = TRUE)))
    })
    out <- do.call(rbind, rows); rownames(out) <- NULL; out
  }


  distribution <- function(z) {
    m <- mean(z); s <- stats::sd(z); c2 <- mean((z-m)^2)
    data.frame(n=length(z),mean=m,sd=s,
      q025=unname(stats::quantile(z,.025)),median=stats::median(z),
      q975=unname(stats::quantile(z,.975)),minimum=min(z),maximum=max(z),
      fraction_abs_z_gt_1_96=mean(abs(z)>1.96),fraction_abs_z_gt_3=mean(abs(z)>3),
      skewness=if(c2>0) mean((z-m)^3)/c2^1.5 else NA_real_,
      excess_kurtosis=if(c2>0) mean((z-m)^4)/c2^2-3 else NA_real_)
  }
  grouped_distribution <- function(z, labels, key) {
    pieces <- split(seq_along(z),as.character(labels),drop=TRUE)
    ans <- lapply(names(pieces),function(k) {
      a <- distribution(z[pieces[[k]]]); a[[key]]<-k; a
    })
    do.call(rbind,ans)
  }
  binned_patterns <- function(frame,z,eta,fixed_eta,kind) {
    vars <- if(kind=="offset") c("ONSET_mean","clim_anomaly_tw90","photo_tw90",
      "clim_background_tw90","clim_predictability_tw90","clim_trend_tw90") else
      c("onset_plasticity_z","offset_plasticity_bio_z","year_decade")
    xlist <- c(list(fixed_linear_predictor=fixed_eta,conditional_linear_predictor=eta,
      calendar_year=frame$YEAR),frame[intersect(vars,names(frame))])
    rows <- lapply(names(xlist),function(nm) {
      x <- as.numeric(xlist[[nm]])
      check(length(x)==length(z) && all(is.finite(x)),paste("Invalid diagnostic predictor",nm))
      if(nm=="calendar_year") g<-as.character(x) else {
        cuts<-unique(stats::quantile(x,seq(0,1,length.out=21),names=FALSE))
        g<-if(length(cuts)<2L) rep("constant",length(x)) else
          cut(x,cuts,include.lowest=TRUE,labels=FALSE)
      }
      by<-split(seq_along(z),g,drop=TRUE)
      do.call(rbind,lapply(names(by),function(k) {
        i<-by[[k]]
        cbind(predictor=nm,bin=k,x_min=min(x[i]),x_median=stats::median(x[i]),
          x_max=max(x[i]),distribution(z[i]))
      }))
    })
    do.call(rbind,rows)
  }
  temporal_patterns <- function(frame,z) {
    keys<-c("SPECIES","SITE_ID","YEAR")
    base<-frame[keys]; base$idx<-seq_len(nrow(base))
    cor_ok<-function(x,y) if(length(x)>=3L && stats::sd(x)>0 && stats::sd(y)>0)
      stats::cor(x,y) else NA_real_
    out<-list()
    for(lag in 1:3) {
      earlier<-base; earlier$YEAR<-earlier$YEAR+lag
      pair<-merge(base,earlier,by=keys,suffixes=c("_current","_earlier"),sort=FALSE)
      all_groups<-c(list(all_populations=seq_len(nrow(pair))),
        split(seq_len(nrow(pair)),paste0("species:",pair$SPECIES),drop=TRUE))
      for(k in names(all_groups)) {
        i<-all_groups[[k]]; x<-z[pair$idx_earlier[i]]; y<-z[pair$idx_current[i]]
        out[[length(out)+1L]]<-data.frame(lag_calendar_years=lag,stratum=k,n_pairs=length(i),
          Pearson_r=cor_ok(x,y),mean_product=if(length(i)) mean(x*y) else NA_real_)
      }
    }
    do.call(rbind,out)
  }
  distance_products <- function(tab) {
    bands<-data.frame(lo=c(0,15,30,60,100,200),hi=c(15,30,60,100,200,500))
    out<-list()
    for(yr in sort(unique(tab$YEAR))) {
      a<-tab[tab$YEAR==yr,,drop=FALSE]; n<-nrow(a)
      D<-as.matrix(stats::dist(as.matrix(a[c("x_km","y_km")])))
      for(scope in c("all_networks","within_networks")) {
        v<-a$z_standardized
        if(scope=="within_networks") v<-v-ave(v,a$bms_id,FUN=mean)
        v<-v-mean(v); den<-mean(v^2)
        allowed<-upper.tri(D)
        if(scope=="within_networks") allowed<-allowed & outer(a$bms_id,a$bms_id,"==")
        for(j in seq_len(nrow(bands))) {
          ind<-which(allowed & D>bands$lo[j] & D<=bands$hi[j],arr.ind=TRUE)
          np<-nrow(ind)
          out[[length(out)+1L]]<-data.frame(YEAR=yr,scope=scope,d_min_km=bands$lo[j],
            d_max_km=bands$hi[j],n_pairs=np,
            standardized_product=if(np>=100L && is.finite(den) && den>0)
              mean(v[ind[,1]]*v[ind[,2]])/den else NA_real_,
            status=if(np<100L) "fewer_than_100_pairs" else if(den<=0) "zero_variance" else "descriptive")
        }
      }
    }
    do.call(rbind,out)
  }
  png_plot <- function(file, draw_fun) {
    grDevices::png(file,width=1600,height=1100,res=180)
    on.exit(grDevices::dev.off(),add=TRUE)
    graphics::par(mar=c(5,5,3,1)+.1)
    draw_fun()
  }
  render_primary <- function(out,key,qq,bins,sp,temporal,moran) {
    png_plot(file.path(out,"QQ_primary.png"),function() {
      graphics::plot(qq$normal_quantile,qq$residual_quantile,pch=16,cex=.5,
        xlab="Standard normal quantile",ylab="One-sample residual quantile",main=key)
      graphics::abline(0,1,lty=2)
    })
    for(nm in unique(bins$predictor)) {
      d<-bins[bins$predictor==nm,,drop=FALSE]; d<-d[order(d$x_median),,drop=FALSE]
      safe<-gsub("[^A-Za-z0-9_]","_",nm)
      png_plot(file.path(out,paste0("mean_vs_",safe,".png")),function() {
        graphics::plot(d$x_median,d$mean,type="b",pch=16,
          xlab=nm,ylab="Mean one-sample residual (binned)",main=key)
        graphics::abline(h=0,lty=2)
      })
      png_plot(file.path(out,paste0("dispersion_vs_",safe,".png")),function() {
        graphics::plot(d$x_median,d$sd,type="b",pch=16,
          xlab=nm,ylab="Residual SD (binned)",main=key)
        graphics::abline(h=1,lty=2)
      })
    }
    png_plot(file.path(out,"species_dispersion.png"),function() {
      graphics::plot(sp$n,sp$sd,log="x",pch=16,
        xlab="Observations per species (log scale)",ylab="One-sample residual SD",main=key)
      graphics::abline(h=1,lty=2)
    })
    for(scope in c("all_networks","within_networks")) {
      d<-moran[moran$scope==scope & moran$aggregation=="variance_standardized",,drop=FALSE]
      d<-d[order(d$YEAR),,drop=FALSE]
      png_plot(file.path(out,paste0("Moran_",scope,".png")),function() {
        graphics::plot(d$YEAR,d$I,type="b",pch=16,xlab="Year",ylab="Moran's I",main=paste(key,scope))
        graphics::lines(d$YEAR,d$null_mean,lty=2)
      })
    }
    t<-temporal[temporal$stratum=="all_populations",,drop=FALSE]
    if(any(is.finite(t$Pearson_r))) png_plot(file.path(out,"calendar_lag_correlation.png"),function() {
      graphics::plot(t$lag_calendar_years,t$Pearson_r,type="b",pch=16,xlab="Calendar lag (years)",
        ylab="Pooled residual correlation (descriptive)",main=key)
      graphics::abline(h=0,lty=2)
    })
    invisible(NULL)
  }
  source_jobs <- function(root) {
    off <- file.path(root,"output","phenology_plasticity","offset_with_onset_15km","run_3911d7bc87ac")
    ab <- file.path(root,"output","population_trends","spatial_abundance_onset_adjusted_15km",
      "run_20261004_214012_44136")
    refined <- file.path(ab,"refinement_univoltine","run_20261005_073241_44136")
    jobs <- list()
    # The two smaller, previously tested Gamma diagnostic systems run first.
    for(g in c("univoltine","multivoltine")) {
      key <- paste0("abundance_",g)
      fit <- if(g=="univoltine") file.path(refined,"univoltine__full","fit.rds") else
        file.path(ab,"multivoltine__full","fit.rds")
      jobs[[key]] <- list(key=key,kind="abundance",group=g,fit=fit,
        input=file.path(ab,paste0("input_",g,".rds")),mesh=file.path(ab,"mesh.rds"),
        cpp=file.path(ab,"compiled","trend_spde.cpp"))
    }
    for(g in c("univoltine","multivoltine")) {
      key <- paste0("offset_",g); src <- file.path(off,key)
      jobs[[key]] <- list(key=key,kind="offset",group=g,fit=file.path(src,"fit_with_onset.rds"),
        input=file.path(src,"input.rds"),parameters=file.path(src,"parameters.rds"),
        complete=file.path(src,"completed.rds"))
    }
    jobs
  }
  manifest <- function(cfg) {
    roles <- c("input","fit","parameters","complete","mesh","cpp")
    tab <- do.call(rbind,lapply(cfg$jobs,function(j) {
      r<-intersect(roles,names(j))
      data.frame(job=j$key,role=r,path=unname(unlist(j[r])),stringsAsFactors=FALSE)
    }))
    check(all(file.exists(tab$path)),paste("Required source file(s) missing:",
      paste(tab$path[!file.exists(tab$path)],collapse="\n")))
    u<-unique(tab$path); hashes<-character(length(u))
    for(i in seq_along(u)) {
      phase(cfg,"FINGERPRINTING",detail=paste(i,"/",length(u),basename(u[i])))
      hashes[i]<-unname(tools::md5sum(u[i]))
      check(!is.na(hashes[i]),paste("Cannot read source:",u[i]))
    }
    info<-file.info(tab$path)
    tab$bytes<-info$size; tab$mtime<-format(info$mtime,"%Y-%m-%d %H:%M:%OS6",tz="UTC")
    tab$md5<-hashes[match(tab$path,u)]
    tab
  }
  source_unchanged <- function(man) {
    inf<-file.info(man$path)
    check(all(file.exists(man$path)) && all(inf$size==man$bytes) &&
      identical(format(inf$mtime,"%Y-%m-%d %H:%M:%OS6",tz="UTC"),as.character(man$mtime)),
      "An input changed during diagnostics. Review sources before using the outputs.")
  }
  residual_worker <- function(cfg,key) {
    packages(); job<-cfg$jobs[[key]]
    out<-file.path(cfg$out,key); dir.create(out,recursive=TRUE,showWarnings=FALSE)
    resource_check(cfg$out,cfg$minimum_disk_GiB,cfg$minimum_RAM_GiB)
    src<-cfg$manifest[cfg$manifest$job==key,,drop=FALSE]; source_unchanged(src)
    read_cache<-function(path) {
      if(!file.exists(path)) return(NULL)
      a<-readRDS(path)
      check(identical(a$signature,cfg$signature),paste("Incompatible diagnostic cache:",path))
      a
    }
    done<-read_cache(file.path(out,"complete.rds"))
    if(!is.null(done)) return("DIAGNOSTICS_CACHED")
    small_file<-file.path(out,"observation_metadata.rds")
    meta<-read_cache(small_file)
    paths<-file.path(out,sprintf("draw_%02d.rds",0:cfg$n_draws))
    missing<-which(!file.exists(paths))-1L
    if(is.null(meta)||length(missing)) {
      s<-load_system(job,cfg)
      frame<-s$frame
      phase(cfg,"FACTORING_CONDITIONAL_PRECISION",key,
        paste(s$n_random,"random effects; analytical sparse matrix; no optimization"))
      ch<-Matrix::Cholesky(s$H,perm=TRUE,LDL=FALSE,super=TRUE)
      correction<-as.numeric(Matrix::solve(ch,s$gradient,system="A"))
      decrement<-sqrt(max(0,sum(s$gradient*correction)))
      eta_error<-max(abs(as.numeric(s$Z%*%correction)))
      scale<-if(s$family=="gaussian") s$sigma else 1
      checks<-data.frame(model=key,family=s$family,n_observations=nrow(frame),
        n_random=s$n_random,max_eta_reconstruction_error=s$prediction_error,
        max_random_design_error=s$design_error,
        max_conditional_gradient=max(abs(s$gradient)),conditional_score_norm=decrement,
        max_mode_check_step=max(abs(correction)),max_eta_mode_check_step=eta_error,
        precision_MiB=as.numeric(utils::object.size(s$H))/1024^2,
        Cholesky_MiB=as.numeric(utils::object.size(ch))/1024^2,
        fitted_range_km=s$range_km,fitted_field_SD=s$field_sd,Gaussian_SD=s$sigma,
        Gamma_shape=s$shape,AR1_rho=s$rho,mode_was_updated=FALSE)
      csv(checks,file.path(out,"reconstruction_checks.csv"))
      check(all(is.finite(correction)) && is.finite(decrement) && decrement<0.01 &&
        is.finite(eta_error) && eta_error/scale<1e-4,
        "Analytical system does not reproduce a stationary saved mode; inspect reconstruction_checks.csv. No parameter was updated.")
      phase(cfg,"RECONSTRUCTION_VERIFIED",key,paste("Conditional score norm",signif(decrement,4),
        "| saved effects retained unchanged"))
      keys<-c("SPECIES","SITE_ID","YEAR","bms_id","x_km","y_km","source_row",
        "pop_id","ONSET_mean","ONSET_mean_z","clim_anomaly_tw90","photo_tw90",
        "clim_background_tw90","clim_predictability_tw90","clim_trend_tw90",
        "year_decade","onset_plasticity_z","offset_plasticity_bio_z")
      frame<-frame[intersect(keys,names(frame))]
      meta_new<-list(signature=cfg$signature,frame=frame,y=s$y,eta=s$eta,
        fixed_eta=s$fixed_eta,family=s$family,sigma=s$sigma,shape=s$shape)
      if(is.null(meta)) {
        meta<-meta_new; atomic(meta,small_file)
      } else check(same(meta$eta,s$eta,1e-6),"Cached predictions changed on resume.")
      s$H<-NULL;s$gradient<-NULL;s$u<-NULL;s$fixed_eta<-NULL
      rm(correction,meta_new);invisible(gc())
      for(draw in missing) {
        draw_seed<-as.integer(cfg$seed+match(key,names(cfg$jobs))*1000L+draw)
        phase(cfg,"GENERATING_RESIDUALS",key,if(draw==0L) "EB-mode comparator" else
          paste("Joint conditional draw",draw,"/",cfg$n_draws))
        if(draw==0L) eta<-s$eta else {
          set.seed(draw_seed)
          delta<-as.numeric(precision_transform(ch,stats::rnorm(s$n_random)))
          eta<-s$eta+as.numeric(s$Z%*%delta)
        }
        z<-if(s$family=="gaussian") (s$y-eta)/s$sigma else gamma_residuals(s$y,eta,s$shape)
        check(length(z)==nrow(frame) && all(is.finite(z)),
          "Non-finite/misaligned residuals; none were dropped or clipped.")
        rec<-list(signature=cfg$signature,draw=draw,seed=if(draw==0L) NA_integer_ else draw_seed,
          z=z,distribution=distribution(z))
        atomic(rec,paths[draw+1L])
      }
      rm(s,ch,z,eta,rec,frame);invisible(gc())
    }
    phase(cfg,"SPATIAL_GRAPHS",key,"Same 8-neighbour / 100-km definition; site/year variance standardization")
    agg<-make_aggregation(meta$frame)
    graphs<-make_graphs(agg$meta)
    csv(agg$meta,file.path(out,"aggregation_metadata.csv"))
    all_moran<-dist_rows<-species_rows<-temporal_rows<-list()
    for(draw in 0:cfg$n_draws) {
      rec<-read_cache(paths[draw+1L]);check(!is.null(rec),"Missing residual cache.")
      z<-rec$z
      tab<-agg$meta;tab$z_mean<-as.numeric(agg$B%*%z)
      tab$z_standardized<-tab$z_mean/tab$null_sd_of_mean
      dist_rows[[draw+1L]]<-cbind(model=key,draw=draw,seed=rec$seed,rec$distribution)
      testfile<-file.path(out,sprintf("draw_%02d_moran.rds",draw))
      saved<-read_cache(testfile)
      if(is.null(saved)) {
        phase(cfg,"MORAN_BILATERAL",key,paste("draw",draw,"|",cfg$nperm,
          "permutations per year/scope; positive AND negative tails"))
        mo<-moran_tests(tab,graphs,job,draw,cfg)
        atomic(list(signature=cfg$signature,tests=mo),testfile)
      } else mo<-saved$tests
      all_moran[[draw+1L]]<-mo
      sp<-grouped_distribution(z,meta$frame$SPECIES,"SPECIES")
      species_rows[[draw+1L]]<-cbind(model=key,draw=draw,sp)
      temporal<-temporal_patterns(meta$frame,z)
      temporal_rows[[draw+1L]]<-cbind(model=key,draw=draw,temporal)
      if(draw==1L) {
        phase(cfg,"PRIMARY_PATTERNS",key,"QQ, dispersion, predictors, temporal lags and descriptive distance bands")
        p<-sort(unique(c(.0001,.0005,seq(.001,.999,length.out=999),.9995,.9999)))
        qq<-data.frame(probability=p,normal_quantile=stats::qnorm(p),
          residual_quantile=as.numeric(stats::quantile(z,p)))
        csv(qq,file.path(out,"QQ_primary.csv"))
        bins<-binned_patterns(meta$frame,z,meta$eta,meta$fixed_eta,job$kind)
        csv(bins,file.path(out,"binned_patterns_primary.csv"))
        csv(grouped_distribution(z,meta$frame$bms_id,"bms_id"),file.path(out,"network_residuals_primary.csv"))
        csv(grouped_distribution(z,meta$frame$YEAR,"YEAR"),file.path(out,"year_residuals_primary.csv"))
        csv(distance_products(tab),file.path(out,"distance_products_primary.csv"))
        csv(tab,file.path(out,"site_year_residuals_primary.csv"))
        ix<-head(order(abs(z),decreasing=TRUE),100L)
        csv(cbind(meta$frame[ix,,drop=FALSE],observation=meta$y[ix],
          fitted_at_mode=if(meta$family=="Gamma") exp(meta$eta[ix]) else meta$eta[ix],
          one_sample_z=z[ix]),file.path(out,"largest_residuals_for_audit.csv"))
        if(cfg$make_plots) tryCatch(render_primary(out,key,qq,bins,sp,temporal,mo),
          error=function(e)writeLines(conditionMessage(e),file.path(out,"plot_warning.txt")))
      }
    }
    mo<-do.call(rbind,all_moran);summ<-summarize_tests(mo)
    csv(mo,file.path(out,"moran_all_draws.csv"))
    csv(summ,file.path(out,"moran_summary_all_draws.csv"))
    csv(summ[summ$draw==1L & summ$aggregation=="variance_standardized",,drop=FALSE],
      file.path(out,"moran_summary_primary.csv"))
    csv(do.call(rbind,dist_rows),file.path(out,"residual_distribution_all_draws.csv"))
    csv(do.call(rbind,species_rows),file.path(out,"species_residuals_all_draws.csv"))
    csv(do.call(rbind,temporal_rows),file.path(out,"temporal_patterns_all_draws.csv"))
    source_unchanged(src)
    atomic(list(signature=cfg$signature,finished_at=Sys.time(),status="DIAGNOSTICS_EXPORTED"),
      file.path(out,"complete.rds"))
    "DIAGNOSTICS_EXPORTED"
  }
  write_readme <- function(cfg) {
    writeLines(c(
      "CURRENT FOUR-FIT POST-HOC DIAGNOSTICS -- phenoIMPACT",
      paste("Version:",VERSION),"No fitted model or input is modified; no optimization is performed.",
      "Sources are fixed explicitly in source_manifest.csv; univoltine abundance is the refined fit.",
      "Gaussian offsets: exact conditional multivariate normal given estimated parameters.",
      "Gamma abundance: multivariate normal / Laplace approximation at the saved random-effect mode.",
      "All fixed, variance, range, correlation and dispersion parameters are held at their estimates.",
      "H = Q_prior + Z' W Z is reconstructed analytically; no saved native TMB pointer is called.",
      "The conditional score solve verifies the stored mode; it is NEVER applied as an update.",
      "Draw 0 = mode-residual comparator only. Draw 1 = predeclared primary. Draws 2:5 = sensitivity.",
      paste("Configured draws:",cfg$n_draws,"; permutations:",cfg$nperm),
      "Primary aggregation: mean per site/year, then equal means of colocated sites in each network/year,",
      "divided by sqrt(sum of squared observation weights). This variance standardization assumes iid N(0,1)",
      "observation residuals under the ideal reference; heterogeneity or dependence can invalidate it.",
      "Original unstandardized means are also retained as a secondary sensitivity, not a replacement.",
      "Moran: up to 8 nearest neighbours within 100 km; symmetrized; row-standardized; isolates excluded.",
      "Global tests centre/permute within each year. Within-network tests centre and permute within network/year.",
      "Directional p = (1 + permutation-tail count) / (nperm + 1), including ties.",
      "Bilateral p = min(1, 2*min(p_positive,p_negative)); BH across years per model/draw/aggregation/scope.",
      "No p-values are pooled across draws, and no passing draw is selected. No latent covariance is exported.",
      "Distance products: six bands 0-15/15-30/30-60/60-100/100-200/200-500 km, primary standardized scores,",
      "each pair counted once; mean(v_i*v_j)/mean(v^2), at least 100 pairs. DESCRIPTIVE, not p-values/CIs.",
      "Distance and Moran radii are NOT the model's 15-km mesh-construction cutoff or its estimated range.",
      "QQ and dispersion summaries use quantile residuals: N(0,1) is the reference even for Gamma models.",
      "Conditional-fitted values on binned plots use the EB-mode predictor; patterns may be affected by fitting.",
      "Fixed-linear-predictor plots exclude fitted random effects. No smoothing/regression is fitted to the residuals.",
      "Calendar-lag correlations use only exact same-population pairs separated by 1/2/3 years. Gaps are not compressed.",
      "Temporal correlations are descriptive and pooled across available pairs; species-level summaries are supplied.",
      "No formal temporal-null calibration, residual bootstrap, model refit, R2 or causal inference is performed.",
      "Large residuals are exported for inspection only: no observation is removed, winsorized or imputed.",
      "Limitations: no propagation of measurement/plasticity estimation error; no validation of the latent distributions;",
      "not an external validation; aggregation can conceal species-specific spatial dependence. Moran p-values rely",
      "on an approximate exchangeability reference, and non-significance does not prove absence of dependence.",
      "DIAGNOSTICS_EXPORTED is a workflow status, NOT a declaration that assumptions are satisfied.",
      "Heavy Gaussian loading/factorization can require substantial RAM. Run serially without other large jobs.",
      "Derived residual caches are in this diagnostic run only and allow resuming without refactorizing completed jobs.",
      "Source hashes are recalculated when resuming; incompatible input/settings are not silently reused.",
      "No packages or fonts are installed. Figures are base-R diagnostic PNGs, not manuscript figures.",
      "Method: https://sdmtmb.github.io/sdmTMB/articles/residual-checking.html",
      "Implementation reference: https://raw.githubusercontent.com/sdmTMB/sdmTMB/v1.1.0/src/sdmTMB.cpp"
    ),file.path(cfg$out,"README.txt"))
  }
  collect <- function(cfg) {
    read_all<-function(filename) {
      rows<-lapply(names(cfg$jobs),function(key) {
        p<-file.path(cfg$out,key,filename)
        if(!file.exists(p)) return(NULL)
        d<-utils::read.csv(p,check.names=FALSE,stringsAsFactors=FALSE)
        if(!"model" %in% names(d)) d$model<-key
        d
      })
      rows<-Filter(Negate(is.null),rows)
      if(!length(rows)) return(NULL)
      do.call(rbind,rows)
    }
    for(nm in c("reconstruction_checks.csv","moran_summary_primary.csv",
      "moran_summary_all_draws.csv","residual_distribution_all_draws.csv",
      "species_residuals_all_draws.csv","temporal_patterns_all_draws.csv",
      "distance_products_primary.csv")) {
      d<-read_all(nm)
      if(!is.null(d)) csv(d,file.path(cfg$out,paste0("all_",nm)))
    }
    all<-read_all("moran_summary_all_draws.csv")
    if(!is.null(all)) {
      a<-all[all$draw>0,,drop=FALSE]
      pieces<-split(seq_len(nrow(a)),interaction(a$model,a$aggregation,a$scope,drop=TRUE))
      sens<-do.call(rbind,lapply(pieces,function(i) {
        d<-a[i,,drop=FALSE]
        data.frame(model=d$model[1],aggregation=d$aggregation[1],scope=d$scope[1],n_draws=nrow(d),
          min_median_I=min(d$median_I),max_median_I=max(d$median_I),
          min_n_sig_BH_two_sided=min(d$n_sig_BH_two_sided),
          median_n_sig_BH_two_sided=stats::median(d$n_sig_BH_two_sided),
          max_n_sig_BH_two_sided=max(d$n_sig_BH_two_sided))
      }))
      csv(sens,file.path(cfg$out,"moran_draw_sensitivity.csv"))
    }
    logs<-list.files(cfg$out,pattern="\\.log$",full.names=TRUE,recursive=TRUE)
    for(log in logs) tryCatch(writeLines(tail(readLines(log,warn=FALSE),500),
      paste0(log,"_snapshot.txt")),error=function(e)NULL)
    if(requireNamespace("zip",quietly=TRUE)) tryCatch({
      files<-list.files(cfg$out,recursive=TRUE,full.names=FALSE)
      files<-files[grepl("\\.(csv|txt|png|R)$",files,ignore.case=TRUE)]
      zip::zipr(file.path(cfg$out,"results_to_review.zip"),files=files,root=cfg$out,
        mode="mirror",include_directories=FALSE)
    },error=function(e)message("Optional archive failed; tables remain available: ",conditionMessage(e)))
    invisible(NULL)
  }
  controller <- function(cfg) {
    on.exit(unlink(cfg$lock,recursive=TRUE),add=TRUE)
    status<-data.frame(task=names(cfg$jobs),status="PENDING",message="",stringsAsFactors=FALSE)
    save_status<-function() csv(status,file.path(cfg$out,"workflow_status.csv"))
    save_status()
    end_status<-"DIAGNOSTICS_INCOMPLETE"
    tryCatch({
      cfg$packages<-packages();resource_check(cfg$out,cfg$minimum_disk_GiB,cfg$minimum_RAM_GiB)
      csv(self_tests(),file.path(cfg$out,"implementation_self_tests.csv"))
      csv(cfg$packages,file.path(cfg$out,"package_versions.csv"))
      cfg$manifest<-manifest(cfg)
      sig<-digest::digest(list(version=VERSION,source=cfg$manifest[c("job","role","path","md5")],
        packages=cfg$packages,n_draws=cfg$n_draws,nperm=cfg$nperm,seed=cfg$seed,
        permutation_seed=cfg$permutation_seed),algo="sha256")
      if(!is.null(cfg$signature)) check(identical(sig,cfg$signature),
        "Input/settings/package versions changed; start a separate diagnostic run.")
      cfg$signature<-sig
      atomic(cfg,file.path(cfg$out,"diagnostic_config.rds"))
      csv(cfg$manifest,file.path(cfg$out,"source_manifest.csv"))
      write_readme(cfg)
      for(key in names(cfg$jobs)) {
        i<-match(key,status$task);status$status[i]<-"RUNNING";save_status()
        log<-file.path(cfg$out,paste0(key,"_",format(Sys.time(),"%Y%m%d_%H%M%S"),".log"))
        value<-tryCatch({
          child<-callr::r_bg(function(engine,cfg,key) {
            e<-new.env(parent=globalenv());sys.source(engine,envir=e)
            e$residual_worker(cfg,key)
          },args=list(cfg$engine,cfg,key),libpath=cfg$libpath,stdout=log,stderr="2>&1",
            supervise=TRUE,user_profile=FALSE,system_profile=FALSE,
            env=c(callr::rcmd_safe_env(),OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1"))
          atomic(list(pid=child$get_pid(),task=key,log=log),file.path(cfg$out,"active_worker.rds"))
          child$wait();child$get_result()
        },error=function(e) {
          status$message[i]<<-conditionMessage(e);"DIAGNOSTIC_ERROR"
        },finally={if(exists("child",inherits=FALSE) && child$is_alive()) child$kill()})
        status$status[i]<-value;save_status()
      }
      end_status<-if(all(status$status %in% c("DIAGNOSTICS_EXPORTED","DIAGNOSTICS_CACHED")))
        "DIAGNOSTICS_COMPLETE" else "FINISHED_WITH_DIAGNOSTIC_ERRORS"
    },error=function(e) {
      writeLines(conditionMessage(e),file.path(cfg$out,"controller_error.txt"))
      message(conditionMessage(e))
    })
    finished<-Sys.time()
    save_status()
    writeLines(capture.output(utils::sessionInfo()),file.path(cfg$out,"sessionInfo.txt"))
    atomic(list(status=end_status,started_at=cfg$started_at,finished_at=finished,
      output=cfg$out,workflow=status),file.path(cfg$out,"result.rds"))
    phase(cfg,end_status,detail="No model was refitted. Diagnostic tables still require scientific review.")
    collect(cfg)
    invisible(end_status)
  }
  engine_file <- function(out) {
    e<-environment(engine_file);ns<-ls(e,all.names=TRUE)
    ns<-ns[vapply(ns,function(n)is.function(get(n,e)),logical(1))]
    path<-file.path(out,"diagnostic_engine.R")
    dump(c("VERSION","ROOT","GAMMA_CPP_MD5",ns),file=path,envir=e,control="all")
    normalizePath(path,winslash="/",mustWork=TRUE)
  }
  get_job <- function() {
    j<-state$job
    if(is.null(j)) j<-getOption("phenoimpact.current_diagnostics.job")
    j
  }
  status <- function(job=get_job()) {
    check(!is.null(job),"No diagnostic job in this session.")
    alive<-job$process$is_alive()
    resultpath<-file.path(job$out,"result.rds")
    r<-if(!alive && file.exists(resultpath)) tryCatch(readRDS(resultpath),error=function(e)NULL) else NULL
    duration<-if(!is.null(r)) difftime(r$finished_at,r$started_at,units="hours") else
      difftime(Sys.time(),job$started_at,units="hours")
    message(if(is.null(r)) "Elapsed since launch: " else "Recorded duration: ",
      round(as.numeric(duration),2)," h | controller running: ",alive)
    pp<-file.path(job$out,"progress.rds")
    p<-if(file.exists(pp))tryCatch(readRDS(pp),error=function(e)NULL)else NULL
    if(!is.null(p))message(p$stage," | ",p$job," | ",p$detail)
    q<-file.path(job$out,"workflow_status.csv")
    if(file.exists(q)) print(utils::read.csv(q,stringsAsFactors=FALSE),row.names=FALSE)
    message("Output: ",job$out)
    if(!alive && !identical(job$process$get_exit_status(),0L))message("Controller exit code: ",
      job$process$get_exit_status(),". Read log: ",job$log)
    invisible(p)
  }
  watch <- function(job=get_job(),every=60) {
    check(!is.null(job),"No diagnostic job in this session.")
    tryCatch(repeat {status(job);if(!job$process$is_alive())break;Sys.sleep(every)},
      interrupt=function(e)message("Display paused; diagnostics continue. Use pheno_current_diagnostics$status()."))
    invisible(job)
  }
  stop_job <- function(job=get_job()) {
    check(!is.null(job),"No diagnostic job in this session.")
    if(job$process$is_alive())job$process$kill()
    message("Diagnostic controller stopped. Saved models are untouched; derived caches retained.")
    invisible(job)
  }
  launch <- function(cfg,monitor=FALSE) {
    lock<-cfg$lock
    if(dir.exists(lock)) {
      p<-file.path(lock,"owner.rds")
      owner<-if(file.exists(p))tryCatch(readRDS(p),error=function(e)NULL)else NULL
      check(!is.null(owner),paste("Unidentified diagnostic lock; inspect:",lock))
      alive<-tryCatch(ps::ps_is_running(ps::ps_handle(owner$pid,time=owner$create_time)),error=function(e)FALSE)
      check(!alive,"A diagnostic controller is already running; do not start a second copy.")
      unlink(lock,recursive=TRUE)
    }
    check(dir.create(lock),"Cannot acquire diagnostic lock.")
    launched<-FALSE
    on.exit(if(!launched)unlink(lock,recursive=TRUE),add=TRUE)
    cfg$started_at<-Sys.time()
    cfg$engine<-engine_file(cfg$out)
    atomic(cfg,file.path(cfg$out,"diagnostic_config.rds"))
    # Remove only this derived status so a resumed run does not display an old ending.
    if(file.exists(file.path(cfg$out,"result.rds")))unlink(file.path(cfg$out,"result.rds"))
    phase(cfg,"QUEUED",detail="Read-only diagnostics; fingerprinting and matrix work run in the background")
    log<-file.path(cfg$out,paste0("controller_",format(Sys.time(),"%Y%m%d_%H%M%S"),".log"))
    ready<-file.path(cfg$out,"launch_ready.rds")
    if(file.exists(ready))unlink(ready)
    proc<-callr::r_bg(function(engine,cfg,ready) {
      # Allow the parent to register the lock owner before expensive work.
      for(i in 1:200) {if(file.exists(ready))break;Sys.sleep(.05)}
      if(!file.exists(ready))stop("Launch registration timed out.")
      e<-new.env(parent=globalenv());sys.source(engine,envir=e);e$controller(cfg)
    },args=list(cfg$engine,cfg,ready),libpath=cfg$libpath,stdout=log,stderr="2>&1",
      supervise=TRUE,user_profile=FALSE,system_profile=FALSE)
    h<-ps::ps_handle(proc$get_pid())
    atomic(list(pid=proc$get_pid(),create_time=ps::ps_create_time(h)),file.path(lock,"owner.rds"))
    atomic(TRUE,ready)
    job<-list(process=proc,out=cfg$out,started_at=cfg$started_at,log=log)
    state$job<-job;options(phenoimpact.current_diagnostics.job=job);launched<-TRUE
    message("Diagnostics queued: two abundance fits (refined univoltine) + two onset-adjusted offsets.")
    message("No model fitting or optimization. Output: ",cfg$out)
    message("Keep R open; avoid other large models. Inspect pheno_current_diagnostics$status().")
    if(monitor)watch(job)
    invisible(job)
  }
  run <- function(root=ROOT,n_draws=5L,nperm=999L,make_plots=TRUE,monitor=FALSE,
                  minimum_disk_GiB=15,minimum_RAM_GiB=10) {
    packages()
    j<-get_job()
    if(!is.null(j) && j$process$is_alive()) {
      message("Diagnostics already running.");if(monitor)watch(j);return(invisible(j))
    }
    check(dir.exists(root),paste("Project root not found:",root))
    check(length(n_draws)==1L && is.finite(n_draws) && n_draws==floor(n_draws) &&
      n_draws>=1 && n_draws<=20,"n_draws must be an integer from 1 to 20.")
    check(length(nperm)==1L && is.finite(nperm) && nperm==floor(nperm) &&
      nperm>=99 && nperm<=99999,"nperm must be an integer from 99 to 99999.")
    check(is.logical(monitor)&&length(monitor)==1L&&!is.na(monitor) &&
      is.logical(make_plots)&&length(make_plots)==1L&&!is.na(make_plots),"Invalid logical settings.")
    check(length(minimum_disk_GiB)==1L && is.finite(minimum_disk_GiB) && minimum_disk_GiB>=1 &&
      length(minimum_RAM_GiB)==1L && is.finite(minimum_RAM_GiB) && minimum_RAM_GiB>=1,
      "Invalid resource guards.")
    root<-normalizePath(root,winslash="/",mustWork=TRUE)
    resource_check(root,minimum_disk_GiB,minimum_RAM_GiB)
    parent<-file.path(root,"output","diagnostics","onset_adjusted_spatial_15km")
    dir.create(parent,recursive=TRUE,showWarnings=FALSE)
    out<-file.path(parent,paste0("run_",format(Sys.time(),"%Y%m%d_%H%M%S"),"_",Sys.getpid()))
    check(!dir.exists(out)&&dir.create(out),"Could not create a unique diagnostic output folder.")
    cfg<-list(version=VERSION,root=root,out=out,jobs=source_jobs(root),signature=NULL,
      n_draws=as.integer(n_draws),nperm=as.integer(nperm),make_plots=make_plots,
      seed=20261005L,permutation_seed=20261901L,minimum_disk_GiB=minimum_disk_GiB,
      minimum_RAM_GiB=minimum_RAM_GiB,libpath=.libPaths(),lock=file.path(parent,"queue.lock"))
    launch(cfg,monitor)
  }
  resume <- function(out=NULL,monitor=FALSE) {
    packages();j<-get_job()
    if(!is.null(j)&&j$process$is_alive()) {message("Diagnostics already running.");return(invisible(j))}
    if(is.null(out)) {check(!is.null(j),"Provide the existing diagnostic run directory in out=.");out<-j$out}
    cfg<-readRDS(file.path(out,"diagnostic_config.rds"))
    check(identical(cfg$version,VERSION),"Diagnostic engine version differs; do not mix caches.")
    resource_check(out,cfg$minimum_disk_GiB,cfg$minimum_RAM_GiB)
    launch(cfg,monitor)
  }
  list(run=run,resume=resume,status=status,watch=watch,stop=stop_job,self_test=self_tests)
})
message("Functions loaded. Start with pheno_current_diagnostics$run(monitor = FALSE).")
