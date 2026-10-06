# phenoIMPACT: two-tailed post-hoc residual Moran diagnostics.
# Reads the four CSVs produced by check_residual_spatial_autocorrelation.R.
# DOES NOT load fitted models, refit models, or overwrite the earlier diagnostics.
# Preserves the original graph, centring, permutation blocks, seed and nperm.
#
# Alternatives: upper tail, lower tail, and equal-tail two-sided test.
# https://r-spatial.github.io/spdep/reference/moran.mc.html
# These are residual-permutation SCREENS, not a model-based bootstrap.

pheno_moran_two_sided <- function(
    run_dir = "E:/phenoIMPACT project/code/phenoIMPACT/output/population_trends/spatial_abundance_15km/run_650c41f15949",
    nperm = 999L,
    seed = 20260928L) {

  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required. No models have been touched.")
  }
  if (length(nperm) != 1L || !is.finite(nperm) || nperm < 1 ||
      nperm != floor(nperm) || nperm > .Machine$integer.max) {
    stop("nperm must be one positive integer.")
  }
  if (length(seed) != 1L || !is.finite(seed)) stop("seed must be a finite number.")
  nperm <- as.integer(nperm)
  input_dir <- file.path(run_dir, "posthoc_residual_moran")
  if (!dir.exists(input_dir)) stop("Residual CSV directory not found: ", input_dir)
  out_dir <- file.path(input_dir, "two_sided")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(out_dir)) stop("Cannot create output directory: ", out_dir)

  # Leave the user's random-number state unchanged after this function returns.
  had_rng <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_rng) old_rng <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_rng) {
      assign(".Random.seed", old_rng, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(list = ".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  jobs <- expand.grid(
    model = c("baseline", "spatial"),
    group = c("univoltine", "multivoltine"),
    stringsAsFactors = FALSE
  )
  jobs$file <- file.path(input_dir, paste0(
    jobs$group, "__", jobs$model, "__site_year_residuals.csv"))
  if (any(!file.exists(jobs$file))) {
    stop("Missing residual CSV(s):\n", paste(jobs$file[!file.exists(jobs$file)], collapse = "\n"))
  }

  screen <- function(tab) {
    required <- c("YEAR", "bms_id", "x_km", "y_km", "z")
    if (!all(required %in% names(tab))) stop("Missing columns: ",
      paste(setdiff(required, names(tab)), collapse = ", "))
    if (!nrow(tab) || anyNA(tab[required])) stop("Empty data or missing values in residual CSV.")
    if (!all(vapply(tab[c("x_km", "y_km", "z")], function(x) {
      is.numeric(x) && all(is.finite(x))
    }, logical(1)))) stop("Coordinates and z must be finite numeric values.")

    # The earlier script resets this seed independently for each group/model.
    set.seed(seed)
    answers <- list()
    null_draws <- list()

    for (yr in sort(unique(tab$YEAR))) {
      a <- tab[tab$YEAR == yr, , drop = FALSE]
      n <- nrow(a)
      D <- as.matrix(stats::dist(as.matrix(a[c("x_km", "y_km")])))
      diag(D) <- Inf

      for (scope in c("all_networks", "within_networks")) {
        DD <- D
        if (scope == "within_networks") {
          DD[outer(a$bms_id, a$bms_id, "!=")] <- Inf
        }
        # Same as the earlier screen: up to 8 nearest neighbours within 100 km,
        # then symmetrise the graph and row-standardise. Symmetrising can give
        # a site more than 8 neighbours. The model mesh cutoff remains 15 km.
        edges <- lapply(seq_len(n), function(i) {
          ids <- which(DD[i, ] > 0 & DD[i, ] <= 100)
          head(ids[order(DD[i, ids])], 8L)
        })
        jj <- as.integer(unlist(edges, use.names = FALSE))
        W <- Matrix::sparseMatrix(
          i = rep(seq_len(n), lengths(edges)), j = jj,
          x = rep(1, length(jj)), dims = c(n, n))
        W <- methods::as((W + Matrix::t(W)) > 0, "dMatrix")
        keep <- Matrix::rowSums(W) > 0
        W <- W[keep, keep, drop = FALSE]
        v <- a$z[keep]
        bms <- a$bms_id[keep]
        nn <- length(v)

        ans <- data.frame(
          YEAR = yr, scope = scope, n_connected = nn, n_excluded = n - nn,
          I = NA_real_, null_mean = NA_real_, null_sd = NA_real_,
          null_q025 = NA_real_, null_q975 = NA_real_,
          I_minus_null_mean = NA_real_, p_positive = NA_real_,
          p_negative = NA_real_, p_two_sided = NA_real_,
          nperm = nperm, status = "too_few_connected",
          stringsAsFactors = FALSE)

        if (nn >= 15L && stats::sd(v) > 1e-12) {
          if (scope == "within_networks") v <- v - ave(v, bms, FUN = mean)
          v <- v - mean(v)
          den <- sum(v^2)
          ans$status <- "zero_variance_after_centring"
          if (den > 1e-20) {
            W <- Matrix::Diagonal(x = 1 / Matrix::rowSums(W)) %*% W
            I <- as.numeric(crossprod(v, W %*% v) / den)
            blocks <- if (scope == "within_networks") {
              split(seq_len(nn), bms)
            } else list(seq_len(nn))
            sims <- numeric(nperm)

            for (start in seq.int(1L, nperm, by = 100L)) {
              m <- min(100L, nperm - start + 1L)
              P <- matrix(v, nrow = nn, ncol = m)
              for (j in seq_len(m)) {
                for (ind in blocks) P[ind, j] <- v[ind[sample.int(length(ind))]]
              }
              # Preserve two dimensions, including a final batch of size one.
              WP <- as.matrix(W %*% P)
              dim(WP) <- c(nn, m)
              sims[start + seq_len(m) - 1L] <- colSums(P * WP) / den
            }
            if (any(!is.finite(sims))) stop("Non-finite permutation statistics.")

            # Include the observed statistic (+1) and include ties in each tail.
            # Do not infer the negative-tail p as simply 1 - p_positive.
            p_hi <- (1 + sum(sims >= I)) / (nperm + 1)
            p_lo <- (1 + sum(sims <= I)) / (nperm + 1)
            ans$I <- I
            ans$null_mean <- mean(sims)
            ans$null_sd <- stats::sd(sims)
            ans$null_q025 <- unname(stats::quantile(sims, 0.025))
            ans$null_q975 <- unname(stats::quantile(sims, 0.975))
            ans$I_minus_null_mean <- I - ans$null_mean
            ans$p_positive <- p_hi
            ans$p_negative <- p_lo
            ans$p_two_sided <- min(1, 2 * min(p_hi, p_lo))
            ans$status <- "ok"
            null_draws[[paste(yr, scope, sep = "__")]] <- sims
          }
        } else if (nn >= 15L) {
          ans$status <- "zero_variance"
        }
        answers[[length(answers) + 1L]] <- ans
      }
    }
    out <- do.call(rbind, answers)
    # Separate BH family for each group x model x scope, across years.
    # Separate correction for each of the three alternative hypotheses.
    for (tail in c("positive", "negative", "two_sided")) {
      out[[paste0("q_BH_", tail)]] <- ave(out[[paste0("p_", tail)]], out$scope,
        FUN = function(p) stats::p.adjust(p, method = "BH"))
    }
    list(table = out, null_draws = null_draws)
  }

  results <- summaries <- list()
  med <- function(x) if (any(is.finite(x))) stats::median(x[is.finite(x)]) else NA_real_

  for (i in seq_len(nrow(jobs))) {
    group <- jobs$group[i]
    model <- jobs$model[i]
    key <- paste(group, model, sep = "__")
    message("\n=== ", group, " | ", model, " ===")
    tab <- utils::read.csv(jobs$file[i], stringsAsFactors = FALSE)
    res <- screen(tab)
    mo <- res$table
    mo$group <- group
    mo$model <- model
    utils::write.csv(mo, file.path(out_dir, paste0(key, "__moran_two_sided.csv")), row.names = FALSE)
    saveRDS(res$null_draws, file.path(out_dir, paste0(key, "__permutation_reference.rds")))
    results[[key]] <- mo

    for (scope in c("all_networks", "within_networks")) {
      d <- mo[mo$scope == scope, , drop = FALSE]
      summaries[[paste(key, scope, sep = "__")]] <- data.frame(
        group = group, model = model, scope = scope,
        years_tested = sum(d$status == "ok"),
        years_skipped = sum(d$status != "ok"),
        median_I = med(d$I), median_null_mean = med(d$null_mean),
        median_I_minus_null_mean = med(d$I_minus_null_mean),
        n_sig_BH_positive = sum(d$q_BH_positive < 0.05, na.rm = TRUE),
        n_sig_BH_negative = sum(d$q_BH_negative < 0.05, na.rm = TRUE),
        n_sig_BH_two_sided = sum(d$q_BH_two_sided < 0.05, na.rm = TRUE))
    }
    message("Saved ", key, "; ", sum(mo$status == "ok"), " year/scope tests.")
  }

  all_tests <- do.call(rbind, results)
  summary <- do.call(rbind, summaries)
  rownames(all_tests) <- rownames(summary) <- NULL
  utils::write.csv(all_tests, file.path(out_dir, "moran_all_models_two_sided.csv"), row.names = FALSE)
  utils::write.csv(summary, file.path(out_dir, "moran_summary_two_sided.csv"), row.names = FALSE)

  # Check that the original I and positive-tail results are reproduced.
  previous_file <- file.path(input_dir, "moran_all_models.csv")
  if (file.exists(previous_file)) {
    previous <- utils::read.csv(previous_file, stringsAsFactors = FALSE)
    keys <- c("group", "model", "YEAR", "scope")
    cols <- c(keys, "I", "p_positive")
    if (all(cols %in% names(previous))) {
      cmp <- merge(previous[cols], all_tests[cols], by = keys,
                   suffixes = c("_previous", "_new"), all = TRUE)
      cmp$delta_I <- cmp$I_new - cmp$I_previous
      cmp$delta_p_positive <- cmp$p_positive_new - cmp$p_positive_previous
      utils::write.csv(cmp, file.path(out_dir, "agreement_with_previous_screen.csv"), row.names = FALSE)
      changed <- any(abs(cmp$delta_I) > 1e-9, na.rm = TRUE) ||
                 any(abs(cmp$delta_p_positive) > 1e-9, na.rm = TRUE)
      if (changed) warning("Original results differ: inspect agreement_with_previous_screen.csv. ",
                           "Changes to nperm, seed or RNG settings can change permutation p-values.")
    }
  }
  utils::write.csv(data.frame(
    residual_file = jobs$file,
    md5 = unname(tools::md5sum(jobs$file)),
    nperm = nperm, seed = seed),
    file.path(out_dir, "input_manifest.csv"), row.names = FALSE)
  writeLines(capture.output(utils::sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
  writeLines(c(
    "POST-HOC TWO-TAILED MORAN RESIDUAL SCREEN",
    "No fitted models were loaded or refitted. Earlier diagnostics are unchanged.",
    "Input: previously saved site/year residual averages, not raw observations.",
    "Graph: up to 8 nearest neighbours within 100 km; symmetrise; row-standardise.",
    "All-networks: global centring and unrestricted permutations within each year.",
    "Within-networks: network-specific centring and permutations within each network/year.",
    "Isolated locations are excluded; at least 15 connected locations are required.",
    "p_positive = (1 + count(I_perm >= I_obs)) / (nperm + 1).",
    "p_negative = (1 + count(I_perm <= I_obs)) / (nperm + 1).",
    "p_two_sided = min(1, 2 * min(p_positive, p_negative)); ties included.",
    "BH: across years separately for each group/model/scope and each alternative.",
    "The null quantiles are permutation reference limits, not confidence intervals for I.",
    "The permutation null mean need not be zero, especially with network blocking.",
    "LIMITATIONS: fitted residuals need not be exchangeable under the fitted model.",
    "This extends the earlier diagnostic screen; it is not a model-based calibration.",
    "Negative I, even significant in this screen, is not by itself proof of overfitting.",
    "A null result here does not rule out dependence at other scales or aggregation levels.",
    "Reference: https://r-spatial.github.io/spdep/reference/moran.mc.html"
  ), file.path(out_dir, "README.txt"))

  message("\nDONE. No models were loaded or refitted.")
  message("Results written to: ", normalizePath(out_dir, winslash = "/", mustWork = TRUE))
  print(summary, row.names = FALSE)
  invisible(list(summary = summary, tests = all_tests, output_dir = out_dir))
}

# source(file.choose(), encoding = "UTF-8") runs with the configured run directory.
# To define the function without running it, first set:
# options(phenoimpact.moran_skip_autorun = TRUE)
if (!isTRUE(getOption("phenoimpact.moran_skip_autorun", FALSE))) {
  pheno_moran_two_sided()
}
