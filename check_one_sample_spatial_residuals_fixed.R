# phenoIMPACT | Population trends | Joint one-sample residual diagnostics
# Version 1.0.1 -- 2026-09-29
# Fix: extract sparse (i,j,x) directly from Matrix slots, not summary().
# Compatible with version 1.0 checkpoints: model, draws, seeds and Moran tests
# are unchanged. Completed baseline diagnostics are reused automatically.
#
# NO model fitting, fixed-parameter optimisation, package installation, DLL loading,
# or changes to existing fits. Only Matrix and (for the saved mesh) fmesher are used.
#
# The conditional Laplace precision is computed analytically for the EXACT
# Gamma/log + Gaussian random-effect model in this run's trend_spde.cpp:
#   H = Q_prior + Z' diag(shape * y / exp(eta_mode)) Z.
# A JOINT draw u* = u_mode + H^(-1/2) e includes species/site intercepts and
# slopes, population AR(1) effects, and (for the spatial model) ALL annual fields.
# All fixed coefficients AND variance/correlation/shape parameters stay fixed.
# No dense covariance matrix is constructed. The sparse Cholesky can still be
# memory-intensive; jobs run sequentially and completed draws are cached.
#
# Method: https://sdmtmb.github.io/sdmTMB/reference/residuals.sdmTMB.html
#         https://sdmtmb.github.io/sdmTMB/articles/residual-checking.html
# This is a custom implementation for the saved Gamma model, NOT a call to sdmTMB.
# Laplace-normal approximation, not MCMC or a refit-based bootstrap.
#
# Draw 1 is predeclared as primary. Draws 2:5 show sensitivity to the random draw.
# NEVER average residuals across draws, pool their p-values, or pick a passing draw.
#
# To run: source(file.choose(), encoding = "UTF-8")
# To define functions only: options(phenoimpact.one_sample_skip_autorun = TRUE)

pheno_one_sample_diagnostics <- function(
    run_dir = "E:/phenoIMPACT project/code/phenoIMPACT/output/population_trends/spatial_abundance_15km/run_650c41f15949",
    groups = c("univoltine", "multivoltine"),
    models = c("baseline", "spatial"),
    n_draws = 5L,
    nperm = 999L,
    seed = 20260929L,
    out_dir = file.path(run_dir, "posthoc_residual_moran", "one_sample")) {

  # Numerical/cache schema is unchanged by this sparse-extraction bug fix.
  # Keep 1.0 so valid baseline checkpoints from the interrupted run are reused.
  engine_version <- "1.0"
  implementation_version <- "1.0.1"
  fail <- function(...) stop(..., call. = FALSE)
  check <- function(ok, text) if (!isTRUE(ok)) fail(text)
  scalar_integer <- function(x, lower, upper = .Machine$integer.max) {
    is.numeric(x) && length(x) == 1L && is.finite(x) &&
      x == floor(x) && x >= lower && x <= upper
  }
  check(dir.exists(run_dir), paste("Run directory not found:", run_dir))
  check(scalar_integer(n_draws, 1, 100), "n_draws must be an integer from 1 to 100.")
  check(scalar_integer(nperm, 99, 1000000), "nperm must be an integer from 99 to 1000000.")
  check(scalar_integer(seed, 1, .Machine$integer.max - 1000000), "Invalid seed.")
  check(length(groups) > 0L && !anyDuplicated(groups) &&
          all(groups %in% c("univoltine", "multivoltine")), "Invalid groups.")
  check(length(models) > 0L && !anyDuplicated(models) &&
          all(models %in% c("baseline", "spatial")), "Invalid models.")
  for (pkg in c("Matrix", if ("spatial" %in% models) "fmesher")) {
    if (!requireNamespace(pkg, quietly = TRUE)) fail(
      "Required package is missing: ", pkg, ". Nothing has been fitted or overwritten.")
  }
  n_draws <- as.integer(n_draws); nperm <- as.integer(nperm); seed <- as.integer(seed)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  check(dir.exists(out_dir), "Could not create the output directory.")
  log_file <- file.path(out_dir, "diagnostic_log.txt")
  note <- function(...) {
    txt <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste0(..., collapse = ""))
    message(txt); cat(txt, "\n", file = log_file, append = TRUE)
  }
  write_csv <- function(x, file) utils::write.csv(x, file, row.names = FALSE, na = "")
  save_cache <- function(x, path) {
    tmp <- tempfile(pattern = ".checkpoint_", tmpdir = dirname(path))
    on.exit(unlink(tmp), add = TRUE)
    saveRDS(x, tmp, compress = FALSE)
    # Only this diagnostic's own derived cache is replaced, never an input/fit.
    check(file.copy(tmp, path, overwrite = TRUE), paste("Cannot save cache:", path))
  }
  had_rng <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_rng) old_rng <- get(".Random.seed", envir = .GlobalEnv)
  old_kind <- RNGkind()
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_rng) {
      assign(".Random.seed", old_rng, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(list = ".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")

  # Read sparse entries without depending on summary() S3/S4 dispatch or on
  # which packages/functions happen to be attached in the user's session.
  # Triplet Matrix slots use zero-based row/column indices; return ONE-based
  # indices to match sparseMatrix() below. Never densify the real projection.
  # Matrix reference:
  # https://stat.ethz.ch/R-manual/R-devel/library/Matrix/html/TsparseMatrix-class.html
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

  # Small deterministic regression test, including repeated site rows and
  # the zero-based year offsets used by the annual-field design. No RNG used.
  check_sparse_triplets <- function() {
    a <- Matrix::sparseMatrix(
      i = c(1L, 1L, 2L, 2L, 3L, 3L),
      j = c(1L, 3L, 2L, 4L, 1L, 4L),
      x = c(0.25, 0.75, 0.4, 0.6, 0.8, 0.2), dims = c(3L, 4L))
    site <- c(3L, 1L, 3L, 2L, 1L)
    year <- c(0L, 1L, 2L, 1L, 0L)
    repeated <- a[site, , drop = FALSE]
    tr <- sparse_triplets(repeated)
    restored <- Matrix::sparseMatrix(i = tr$i, j = tr$j, x = tr$x,
                                     dims = dim(repeated))
    err <- max(abs(as.matrix(restored) - as.matrix(repeated)))
    field <- matrix(seq_len(12L) / 7, nrow = 4L, ncol = 3L)
    zf <- Matrix::sparseMatrix(i = tr$i,
      j = tr$j + nrow(field) * year[tr$i], x = tr$x,
      dims = c(length(site), length(field)))
    got <- as.numeric(zf %*% as.vector(field))
    expected <- vapply(seq_along(site), function(i) {
      as.numeric(a[site[i], , drop = FALSE] %*% field[, year[i] + 1L])
    }, numeric(1))
    err <- max(err, max(abs(got - expected)))
    # Preserve singleton and empty sparse cases as matrices, not vectors.
    for (b in list(a[1L, , drop = FALSE], a[, 1L, drop = FALSE], a * 0)) {
      tb <- sparse_triplets(b)
      rb <- Matrix::sparseMatrix(i = tb$i, j = tb$j, x = tb$x, dims = dim(b))
      err <- max(err, max(abs(as.matrix(rb) - as.matrix(b))))
    }
    check(is.finite(err) && err < 1e-12,
      "Sparse projection extraction self-test failed. No model was processed.")
    err
  }
  projection_error <- check_sparse_triplets()
  note("Sparse projection extraction self-test passed (v", implementation_version, ").")

  # Sparse square-root solve. If P H P' = L L', delta = P' L'^(-1) e.
  precision_transform <- function(ch, e) {
    a <- Matrix::solve(ch, e, system = "Lt")
    Matrix::solve(ch, a, system = "Pt")
  }
  # Deterministic check of factor orientation and fill-reducing permutation.
  check_precision_transform <- function() {
    small <- Matrix::Matrix(matrix(c(
      5, 1, 0, 1, 1, 4, 1, 0, 0, 1, 3, 1, 1, 0, 1, 4), 4, 4), sparse = TRUE)
    ch <- Matrix::Cholesky(Matrix::forceSymmetric(small), LDL = FALSE, super = TRUE)
    root <- as.matrix(precision_transform(ch, diag(4)))
    target <- solve(as.matrix(small))
    err <- max(abs(tcrossprod(root) - target))
    check(is.finite(err) && err < 1e-10,
      "Sparse precision sampler self-test failed. No model was processed.")
    err
  }
  sampler_error <- check_precision_transform()
  note("Precision sampler self-test passed. No model will be refitted.")

  # Build observation aggregation ONCE. Match the previous script exactly:
  # observation means -> SITE_ID/year means -> means across colocated sites.
  # The primary score divides each final mean by sqrt(sum of squared observation
  # weights). Under iid N(0,1) observation residuals it has variance one, even
  # when sites contain different numbers of species. Original means are also
  # tested as a secondary screen for direct comparability with earlier output.
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

  # Exact prior precision and design for the template in this specific run.
  make_system <- function(input, pars, spatial, mesh = NULL, A = NULL) {
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
  residual_record <- function(z, agg, group, model, draw, draw_seed, stamp) {
    zmean <- as.numeric(agg$B %*% z)
    tab <- agg$meta
    tab$z_mean <- zmean
    tab$z_standardized <- zmean / tab$null_sd_of_mean
    dist <- data.frame(group = group, model = model,
      residual_type = if (draw == 0L) "EB_mode" else "one_sample",
      draw = draw, seed = draw_seed, n = length(z), mean = mean(z), sd = stats::sd(z),
      q025 = unname(stats::quantile(z, 0.025)), median = stats::median(z),
      q975 = unname(stats::quantile(z, 0.975)), fraction_abs_z_gt_1_96 = mean(abs(z) > 1.96),
      min = min(z), max = max(z), stringsAsFactors = FALSE)
    list(stamp = stamp, tab = tab, distribution = dist)
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

  moran_tests <- function(tab, graphs, group, model, draw) {
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

  all_tests <- distributions <- reconstruction_checks <- list()
  mesh_file <- file.path(run_dir, "mesh.rds")
  if ("spatial" %in% models) check(file.exists(mesh_file), "Saved mesh.rds is missing.")
  note("Settings: ", n_draws, " draws (draw 1 primary); ", nperm, " permutations; sequential jobs.")
  writeLines(c(
    "JOINT ONE-SAMPLE RESIDUAL DIAGNOSTICS -- phenoIMPACT",
    paste("Implementation version:", implementation_version),
    paste("Numerical/cache version:", engine_version),
    "Version 1.0.1 fixes sparse extraction only; version 1.0 checkpoints remain compatible.",
    "INPUT FILES ARE READ ONLY. No fixed-parameter optimisation and no model refits.",
    "The sparse conditional Hessian is rebuilt analytically for Gamma(log).",
    "H = Q_prior + Z' diag(shape * y / exp(eta_mode)) Z.",
    "Prior: independent species and site intercepts/slopes; stationary AR(1) per population",
    "with irregular integer gaps; independent annual Matern fields with the SAVED SPDE mesh.",
    "The annual-field precision is tau^2 * (kappa^4*M0 + 2*kappa^2*M1 + M2).",
    "Each draw includes all latent components jointly; fixed coefficients, shape, SDs,",
    "range and correlations are held at their saved point estimates.",
    "This uses a conditional multivariate-normal/Laplace approximation, NOT MCMC.",
    "Draw 1 is primary. All additional draws are retained as random-draw sensitivity checks.",
    "Do not select the best draw, average residuals across draws, or pool p-values.",
    "PRIMARY AGGREGATION: same two-stage site/year-location means as before, divided by",
    "sqrt(sum of squared observation weights). This accounts for unequal species counts",
    "and equally weighted colocated sites under an iid standard-normal residual reference.",
    "SECONDARY AGGREGATION: original unstandardized means, for comparison with earlier screens.",
    "The secondary raw-mean permutation test can be affected by unequal residual variances.",
    "Graph unchanged: up to 8 nearest neighbours within 100 km; symmetrized; row-standardized.",
    "The graph radius is NOT the mesh cutoff (15 km).",
    "Within-network tests centre and permute within network/year; global tests within year.",
    "Both tails include ties and +1; two-sided p = min(1, 2*min(p_positive,p_negative)).",
    "BH is across years separately for each group/model/draw/aggregation/scope and each tail.",
    "Focus on two-sided results. Directional tests are provided for interpretation.",
    "LIMITATIONS: estimated fixed parameters, Laplace approximation and repeated years mean",
    "these remain diagnostic screens, not a fully calibrated bootstrap or proof of validity.",
    "No significance at this graph/aggregation does not exclude dependence at other scales",
    "or temporal dependence. These diagnostics do not validate the latent-effect distributions",
    "or the abundance/plasticity measurement process.",
    "Cache contains derived residual summaries only. Re-sourcing with the same settings resumes.",
    "Changing data, seeds, draw count or package versions requires a different out_dir.",
    "References:",
    "https://sdmtmb.github.io/sdmTMB/reference/residuals.sdmTMB.html",
    "https://sdmtmb.github.io/sdmTMB/articles/residual-checking.html"
  ), file.path(out_dir, "README.txt"))

  for (group in groups) {
    input_file <- file.path(run_dir, paste0("input_", group, ".rds"))
    check(file.exists(input_file), paste("Missing input:", input_file))
    note("=== ", group, " === Reading saved input; constructing aggregation and graphs.")
    input <- readRDS(input_file)
    check(identical(as.character(input$group), group), "Input group does not match file name.")
    check(all(c("pop_id", "YEAR") %in% names(input$frame)), "Population/year identifiers missing.")
    check(!anyDuplicated(input$frame[c("pop_id", "YEAR")]), "Duplicate population/year observations.")
    agg <- make_aggregation(input$frame)
    graphs <- make_graphs(agg$meta)
    write_csv(agg$meta, file.path(out_dir, paste0(group, "__aggregation_metadata.csv")))
    mesh <- A <- NULL
    if ("spatial" %in% models) {
      mesh <- readRDS(mesh_file)
      A <- methods::as(fmesher::fm_basis(mesh$mesh,
        loc = as.matrix(input$sites[c("x_km", "y_km")])), "CsparseMatrix")
      check(all(input$sites$SITE_ID[input$data$site + 1L] == input$frame$SITE_ID),
        "Site ordering does not match the saved frame.")
    }
    for (model in models) {
      job <- paste(group, model, sep = "__")
      job_dir <- file.path(out_dir, job)
      dir.create(job_dir, showWarnings = FALSE)
      spatial <- model == "spatial"
      fit_file <- if (spatial) file.path(run_dir, paste0(group, "__full"), "fit.rds") else NULL
      check(!spatial || file.exists(fit_file), paste("Missing saved fit:", fit_file))
      stamp <- list(engine_version = engine_version,
        input_md5 = unname(tools::md5sum(input_file)),
        fit_md5 = if (spatial) unname(tools::md5sum(fit_file)) else NULL,
        mesh_md5 = if (spatial) unname(tools::md5sum(mesh_file)) else NULL,
        group = group, model = model, n_draws = n_draws, nperm = nperm, seed = seed,
        Matrix = as.character(utils::packageVersion("Matrix")),
        fmesher = if (spatial) as.character(utils::packageVersion("fmesher")) else NULL)
      cache_file <- function(draw) file.path(job_dir, sprintf("draw_%02d.rds", draw))
      read_cache <- function(draw) {
        path <- cache_file(draw)
        if (!file.exists(path)) return(NULL)
        a <- readRDS(path)
        check(identical(a$stamp, stamp), paste(
          "Cache settings/input mismatch in", path, "-- choose a new out_dir; old results are untouched."))
        a
      }
      cached <- lapply(0:n_draws, read_cache)
      missing <- which(vapply(cached, is.null, logical(1))) - 1L
      fit <- NULL
      if (length(missing)) {
        note(job, ": reconstructing the saved linear predictor and joint conditional precision.")
        if (spatial) {
          fit <- readRDS(fit_file)
          check(identical(fit$data_signature, input$data_signature), "Fit/input signatures disagree.")
          check(isTRUE(fit$checks$numerical_checks_ok), "Saved spatial fit was flagged for numerical review.")
          check(identical(fit$group, group) && identical(fit$variant, "full"), "Wrong saved fit variant/group.")
          pars <- fit$parameters; eta_saved <- fit$eta
        } else {
          pars <- input$par; eta_saved <- input$baseline_eta
        }
        sys <- make_system(input, pars, spatial, mesh, A)
        check(length(eta_saved) == length(sys$eta) && all(is.finite(eta_saved)),
          "Saved eta is missing or invalid.")
        eta_error <- max(abs(sys$eta - eta_saved))
        check(eta_error < 1e-5,
          paste("Prediction reconstruction mismatch (max |delta eta| =", signif(eta_error, 5),
                "). Diagnostics stopped; no parameters were modified."))
        note(job, ": predictions match (max difference = ", signif(eta_error, 3),
             "). Sparse Cholesky of ", sys$n_random, " latent effects.")
        # No dense inverse, no parameter search, no TMB objective is evaluated.
        ch <- Matrix::Cholesky(sys$H, perm = TRUE, LDL = FALSE, super = TRUE)
        correction <- as.numeric(Matrix::solve(ch, sys$gradient, system = "A"))
        max_step <- max(abs(correction))
        decrement <- sqrt(max(0, sum(sys$gradient * correction)))
        check(all(is.finite(correction)) && max_step < 0.01,
          paste("Saved latent modes need review: max conditional Newton step =", signif(max_step, 5),
                ". No optimisation/refit was attempted."))
        if (max_step > 0.001) warning(job, ": saved modes have a non-negligible Newton step; inspect checks.")
        checks <- data.frame(group = group, model = model, n_observations = length(sys$eta),
          n_random = sys$n_random, max_abs_eta_reconstruction_error = eta_error,
          max_abs_conditional_gradient = max(abs(sys$gradient)),
          max_abs_conditional_Newton_step = max_step, conditional_Newton_decrement = decrement,
          sparse_H_size_MB = as.numeric(utils::object.size(sys$H)) / 1024^2,
          Cholesky_size_MB = as.numeric(utils::object.size(ch)) / 1024^2,
          shape = sys$shape, AR1_rho = sys$rho, sampler_selftest_error = sampler_error,
          stringsAsFactors = FALSE)
        write_csv(checks, file.path(job_dir, "reconstruction_checks.csv"))
        # Retain only Z and the factor needed for draws; release the precision.
        sys$H <- NULL; sys$gradient <- NULL; sys$u <- NULL; sys$fixed_eta <- NULL
        correction <- NULL; fit <- NULL; invisible(gc())
        for (draw in missing) {
          draw_seed <- if (draw == 0L) NA_integer_ else as.integer(
            seed + match(group, c("univoltine", "multivoltine")) * 10000L +
              match(model, c("baseline", "spatial")) * 1000L + draw)
          note(job, ": residuals for ", if (draw == 0L) "EB-mode comparator" else paste("joint draw", draw, "/", n_draws))
          if (draw == 0L) eta <- sys$eta else {
            set.seed(draw_seed)
            delta <- as.numeric(precision_transform(ch, stats::rnorm(sys$n_random)))
            check(length(delta) == sys$n_random && all(is.finite(delta)), "Invalid joint Gaussian draw.")
            eta <- sys$eta + as.numeric(sys$Z %*% delta)
          }
          z <- gamma_residuals(input$data$y, eta, sys$shape)
          rec <- residual_record(z, agg, group, model, draw, draw_seed, stamp)
          save_cache(rec, cache_file(draw))
          cached[[draw + 1L]] <- rec
          write_csv(rec$tab, file.path(job_dir, sprintf("draw_%02d__site_year_residuals.csv", draw)))
          rm(z, eta, rec); invisible(gc())
        }
        rm(sys, ch, pars); invisible(gc())
      } else note(job, ": reusing all cached draws; no precision reconstruction needed.")
      checks_path <- file.path(job_dir, "reconstruction_checks.csv")
      if (file.exists(checks_path)) reconstruction_checks[[job]] <- utils::read.csv(checks_path)
      for (draw in 0:n_draws) {
        rec <- cached[[draw + 1L]]
        distributions[[paste(job, draw)]] <- rec$distribution
        tests_path <- file.path(job_dir, sprintf("draw_%02d__moran.csv", draw))
        complete_path <- file.path(job_dir, sprintf("draw_%02d__moran_complete.rds", draw))
        if (file.exists(complete_path)) {
          complete <- readRDS(complete_path)
          check(identical(complete$stamp, stamp), "Moran cache settings mismatch; choose a new out_dir.")
          tests <- complete$tests
        } else {
          # Never trust a CSV without its completed checkpoint.
          note(job, ": bilateral Moran screens, ", if (draw == 0L) "EB-mode" else paste("draw", draw))
          tests <- moran_tests(rec$tab, graphs, group, model, draw)
          save_cache(list(stamp = stamp, tests = tests), complete_path)
        }
        write_csv(tests, tests_path)
        all_tests[[paste(job, draw)]] <- tests
      }
      rm(cached, fit); invisible(gc())
    }
    rm(input, agg, graphs, mesh, A); invisible(gc())
  }

  tests <- do.call(rbind, all_tests); rownames(tests) <- NULL
  summary <- summarize_tests(tests)
  dist <- do.call(rbind, distributions); rownames(dist) <- NULL
  checks <- do.call(rbind, reconstruction_checks); rownames(checks) <- NULL
  primary <- summary[summary$draw == 1L & summary$aggregation == "variance_standardized", , drop = FALSE]
  eb_vs_primary <- summary[summary$draw %in% c(0L, 1L), , drop = FALSE]
  # Every draw is retained. Quantiles below describe Monte Carlo sensitivity,
  # not confidence intervals, new replicated tests, or posterior model probabilities.
  sens <- summary[summary$residual_type == "one_sample", , drop = FALSE]
  sens_keys <- c("group", "model", "aggregation", "scope")
  parts <- split(seq_len(nrow(sens)), interaction(sens[sens_keys], drop = TRUE, lex.order = TRUE))
  sensitivity <- do.call(rbind, lapply(parts, function(ind) {
    a <- sens[ind, , drop = FALSE]
    cbind(a[1, sens_keys, drop = FALSE], data.frame(n_draws = nrow(a),
      median_I_across_draws = stats::median(a$median_I),
      min_median_I = min(a$median_I), max_median_I = max(a$median_I),
      min_n_sig_two_sided = min(a$n_sig_BH_two_sided),
      median_n_sig_two_sided = stats::median(a$n_sig_BH_two_sided),
      max_n_sig_two_sided = max(a$n_sig_BH_two_sided)))
  }))
  rownames(sensitivity) <- NULL
  write_csv(tests, file.path(out_dir, "moran_all_draws.csv"))
  write_csv(summary, file.path(out_dir, "moran_summary_all_draws.csv"))
  write_csv(primary, file.path(out_dir, "moran_summary_primary.csv"))
  write_csv(eb_vs_primary, file.path(out_dir, "moran_EB_vs_primary.csv"))
  write_csv(sensitivity, file.path(out_dir, "draw_sensitivity.csv"))
  write_csv(dist, file.path(out_dir, "residual_distribution_all_draws.csv"))
  write_csv(checks, file.path(out_dir, "model_reconstruction_checks.csv"))
  writeLines(capture.output(utils::sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
  note("DONE. No models were refitted. All completed draws and their tests have been retained.")
  message("\nPRIMARY: joint draw 1, variance-standardized site/year scores, bilateral tests + BH")
  print(primary[c("group", "model", "scope", "years_tested", "median_I",
                  "n_sig_BH_positive", "n_sig_BH_negative", "n_sig_BH_two_sided")], row.names = FALSE)
  message("\nRANDOM-DRAW SENSITIVITY (not pooled p-values)")
  print(sensitivity[sensitivity$aggregation == "variance_standardized", , drop = FALSE], row.names = FALSE)
  message("\nResults: ", normalizePath(out_dir, winslash = "/", mustWork = TRUE))
  invisible(list(primary = primary, sensitivity = sensitivity, output_dir = out_dir))
}

if (!isTRUE(getOption("phenoimpact.one_sample_skip_autorun", FALSE))) {
  pheno_one_sample_diagnostics()
}
