# Post-hoc residual spatial autocorrelation diagnostics
# phenoIMPACT population trends: spatial_abundance_15km
# Reuses saved fits. DOES NOT refit any model.

rm(list = ls())
gc()

# ---- EDIT ONLY IF THE RUN IS IN A DIFFERENT LOCATION -------------------------
run_dir <- "E:/phenoIMPACT project/code/phenoIMPACT/output/population_trends/spatial_abundance_15km/run_650c41f15949"
# -----------------------------------------------------------------------------

nperm <- 999L
seed  <- 20260928L
groups <- c("univoltine", "multivoltine")

stopifnot(dir.exists(run_dir))
if (!requireNamespace("Matrix", quietly = TRUE)) stop("Package 'Matrix' is required.")

out_dir <- file.path(run_dir, "posthoc_residual_moran")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Same residual definition as the fitted workflow.
residual_table <- function(input, eta, shape) {
  mu <- exp(eta)
  y  <- input$data$y
  lo <- pgamma(y, shape = shape, scale = mu / shape, log.p = TRUE)
  hi <- pgamma(y, shape = shape, scale = mu / shape,
               lower.tail = FALSE, log.p = TRUE)
  z <- ifelse(lo < log(.5),
              qnorm(lo, log.p = TRUE),
              qnorm(hi, lower.tail = FALSE, log.p = TRUE))

  if (!all(is.finite(z))) stop("Non-finite Gamma PIT residuals.")

  d <- input$frame[c("SITE_ID", "YEAR", "bms_id", "x_km", "y_km")]
  d$z <- z
  sy <- aggregate(z ~ SITE_ID + YEAR + bms_id + x_km + y_km, d, mean)
  xy <- aggregate(z ~ YEAR + bms_id + x_km + y_km, sy, mean)

  list(
    site_year = xy,
    summary = data.frame(
      n = length(z),
      mean = mean(z),
      sd = sd(z),
      fraction_abs_z_gt_1_96 = mean(abs(z) > 1.96)
    )
  )
}

# Corrected version of moran_screen() from workflow_engine.R.
# The original failed when the final permutation batch had one column because
# W %*% P was simplified; drop = FALSE / as.matrix preserves matrix dimensions.
moran_screen_fixed <- function(tab, nperm = 999L, seed = 20260928L) {
  set.seed(seed)
  answer <- list()

  for (yr in sort(unique(tab$YEAR))) {
    a  <- tab[tab$YEAR == yr, , drop = FALSE]
    xy <- as.matrix(a[c("x_km", "y_km")])
    n  <- nrow(a)
    D  <- as.matrix(dist(xy))
    diag(D) <- Inf

    for (scope in c("all_networks", "within_networks")) {
      DD <- D
      if (scope == "within_networks") {
        DD[outer(a$bms_id, a$bms_id, "!=")] <- Inf
      }

      edges <- lapply(seq_len(n), function(i) {
        ids <- which(DD[i, ] > 0 & DD[i, ] <= 100)
        head(ids[order(DD[i, ids])], 8L)
      })

      W <- Matrix::sparseMatrix(
        i = rep(seq_len(n), lengths(edges)),
        j = unlist(edges), x = 1, dims = c(n, n)
      )
      W <- methods::as((W + Matrix::t(W)) > 0, "dMatrix")

      keep <- Matrix::rowSums(W) > 0
      W <- W[keep, keep, drop = FALSE]
      v <- a$z[keep]
      bms <- a$bms_id[keep]
      nn <- length(v)
      I <- p <- NA_real_

      if (nn >= 15 && sd(v) > 1e-12) {
        if (scope == "within_networks") v <- v - ave(v, bms, FUN = mean)
        v <- v - mean(v)
        den <- sum(v^2)

        if (den > 1e-20) {
          W <- Matrix::Diagonal(x = 1 / Matrix::rowSums(W)) %*% W
          I <- as.numeric(crossprod(v, W %*% v) / den)
          blocks <- if (scope == "within_networks") split(seq_len(nn), bms) else list(seq_len(nn))

          more <- 0L
          for (start in seq.int(1L, nperm, by = 100L)) {
            m <- min(100L, nperm - start + 1L)
            P <- matrix(v, nrow = nn, ncol = m)
            for (j in seq_len(m)) {
              for (ind in blocks) P[ind, j] <- v[ind[sample.int(length(ind))]]
            }

            # BUG FIX: explicitly preserve a 2-D matrix for m = 1.
            WP <- as.matrix(W %*% P)
            dim(WP) <- c(nn, m)
            ip <- colSums(P * WP) / den
            more <- more + sum(ip >= I)
          }
          p <- (1 + more) / (1 + nperm)
        }
      }

      answer[[length(answer) + 1L]] <- data.frame(
        YEAR = yr, scope = scope, n_connected = nn,
        n_excluded = n - nn, I = I, p_positive = p
      )
    }
  }

  out <- do.call(rbind, answer)
  out$q_BH <- ave(out$p_positive, out$scope,
                  FUN = function(p) p.adjust(p, "BH"))
  out
}

all_moran <- list()
all_summary <- list()

for (group in groups) {
  message("\n=== ", group, " ===")

  input_file <- file.path(run_dir, paste0("input_", group, ".rds"))
  fit_file   <- file.path(run_dir, paste0(group, "__full"), "fit.rds")
  stopifnot(file.exists(input_file), file.exists(fit_file))

  input <- readRDS(input_file)
  fit   <- readRDS(fit_file)

  model_specs <- list(
    baseline = list(
      eta = input$baseline_eta,
      shape = exp(input$par$log_shape)
    ),
    spatial = list(
      eta = fit$eta,
      shape = fit$checks$Gamma_shape
    )
  )

  for (model in names(model_specs)) {
    message("  Residuals + Moran: ", model)
    rr <- residual_table(input, model_specs[[model]]$eta,
                         model_specs[[model]]$shape)

    sy_file <- file.path(out_dir, paste0(group, "__", model, "__site_year_residuals.csv"))
    ds_file <- file.path(out_dir, paste0(group, "__", model, "__residual_distribution.csv"))
    mo_file <- file.path(out_dir, paste0(group, "__", model, "__moran.csv"))

    write.csv(rr$site_year, sy_file, row.names = FALSE)
    write.csv(rr$summary, ds_file, row.names = FALSE)

    mo <- moran_screen_fixed(rr$site_year, nperm = nperm, seed = seed)
    mo$group <- group
    mo$model <- model
    write.csv(mo, mo_file, row.names = FALSE)

    all_moran[[paste(group, model, sep = "__")]] <- mo

    all_summary[[paste(group, model, sep = "__")]] <- data.frame(
      group = group,
      model = model,
      years_tested = length(unique(mo$YEAR)),
      median_I_all_networks = median(mo$I[mo$scope == "all_networks"], na.rm = TRUE),
      median_I_within_networks = median(mo$I[mo$scope == "within_networks"], na.rm = TRUE),
      n_sig_BH_all_networks = sum(mo$q_BH[mo$scope == "all_networks"] < 0.05, na.rm = TRUE),
      n_sig_BH_within_networks = sum(mo$q_BH[mo$scope == "within_networks"] < 0.05, na.rm = TRUE),
      n_nominal_p_all_networks = sum(mo$p_positive[mo$scope == "all_networks"] < 0.05, na.rm = TRUE),
      n_nominal_p_within_networks = sum(mo$p_positive[mo$scope == "within_networks"] < 0.05, na.rm = TRUE)
    )
  }

  rm(input, fit)
  gc()
}

moran_all <- do.call(rbind, all_moran)
summary_all <- do.call(rbind, all_summary)
rownames(moran_all) <- NULL
rownames(summary_all) <- NULL

write.csv(moran_all, file.path(out_dir, "moran_all_models.csv"), row.names = FALSE)
write.csv(summary_all, file.path(out_dir, "moran_summary.csv"), row.names = FALSE)

# Direct baseline-vs-spatial change in Moran's I by group/year/scope.
b <- moran_all[moran_all$model == "baseline",
               c("group", "YEAR", "scope", "I", "p_positive", "q_BH")]
s <- moran_all[moran_all$model == "spatial",
               c("group", "YEAR", "scope", "I", "p_positive", "q_BH")]
names(b)[4:6] <- paste0(names(b)[4:6], "_baseline")
names(s)[4:6] <- paste0(names(s)[4:6], "_spatial")
comparison <- merge(b, s, by = c("group", "YEAR", "scope"), all = TRUE)
comparison$delta_I_spatial_minus_baseline <- comparison$I_spatial - comparison$I_baseline
write.csv(comparison, file.path(out_dir, "moran_baseline_vs_spatial.csv"), row.names = FALSE)

message("\nDONE. No models were refitted.")
message("Results written to: ", normalizePath(out_dir, winslash = "/", mustWork = TRUE))
print(summary_all)
