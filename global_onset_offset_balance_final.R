# =============================================================================
# phenoIMPACT | Global onset-offset balance from FINAL abundance models
# Version 1.0 | 2026-10-08
#
# NO MODEL IS REFITTED.
#
# This is an EXACT LINEAR REPARAMETERIZATION of the final abundance model:
#
#   original:
#     year * onset_advancement_z
#     year * offset_plasticity_bio_z
#
# First orient OFFSET identically in both groups:
#   offset_advancement_z =
#       offset_plasticity_bio_z       (univoltine)
#      -offset_plasticity_bio_z       (multivoltine)
#
# Then define two complementary phenological axes:
#
#   coordinated_advancement =
#       (onset_advancement_z + offset_advancement_z) / 2
#
#   season_extension =
#       (onset_advancement_z - offset_advancement_z) / 2
#
# Therefore:
#   +1 coordinated_advancement =
#       +1 SD stronger onset advancement AND
#       +1 SD stronger offset advancement
#
#   +1 season_extension =
#       +1 SD stronger onset advancement AND
#       -1 SD offset advancement
#       (= relatively later termination / greater seasonal extension)
#
# Both axes MUST be included simultaneously. Dropping one would change the
# biological question and confound season displacement with season duration.
#
# Because this transformation is invertible, the transformed model has exactly:
#   - the same fitted values
#   - the same likelihood
#   - the same random effects / AR1 / spatial field
#   - the same information as the final model
#
# The script transforms beta and its full covariance analytically and verifies
# the fixed-effect linear predictor row-by-row.
#
# IMPORTANT INTERPRETATION:
# OFFSET plasticity is the onset-conditioned termination response from the
# final two-stage analysis. Thus "coordinated advancement" means coordinated
# onset advancement + additional/conditional offset advancement in that
# framework; it is not a newly estimated unconditional offset thermal slope.
# =============================================================================

analyse_global_phenology_balance <- function(
  root = "E:/phenoIMPACT project/code/phenoIMPACT",
  final_run = file.path(
    root, "output", "population_trends",
    "spatial_abundance_qc_offset_final_15km",
    "run_eee755f6aac4"
  ),
  out_dir = file.path(
    root, "output", "population_trends",
    "spatial_abundance_qc_offset_final_15km",
    "run_eee755f6aac4",
    "global_onset_offset_balance"
  )
) {

  groups <- c("univoltine", "multivoltine")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  need <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)
  csv <- function(x, path) utils::write.csv(x, path, row.names = FALSE, na = "")

  interaction_name <- function(names_beta, variable) {
    x <- intersect(
      c(
        paste0("year_decade:", variable),
        paste0(variable, ":year_decade")
      ),
      names_beta
    )
    need(length(x) == 1L, paste("Could not identify interaction:", variable))
    x
  }

  fixed_table <- function(beta, V) {
    se <- sqrt(diag(V))
    z <- beta / se
    out <- data.frame(
      term = names(beta),
      estimate = as.numeric(beta),
      std.error = as.numeric(se),
      conf.low = as.numeric(beta - 1.96 * se),
      conf.high = as.numeric(beta + 1.96 * se),
      p_Wald = as.numeric(2 * stats::pnorm(-abs(z))),
      stringsAsFactors = FALSE
    )
    out$ratio_of_decadal_factors <- ifelse(
      grepl("year_decade:", out$term, fixed = TRUE),
      exp(out$estimate),
      NA_real_
    )
    out$ratio_conf.low <- ifelse(
      grepl("year_decade:", out$term, fixed = TRUE),
      exp(out$conf.low),
      NA_real_
    )
    out$ratio_conf.high <- ifelse(
      grepl("year_decade:", out$term, fixed = TRUE),
      exp(out$conf.high),
      NA_real_
    )
    out
  }

  read_final <- function(group) {
    fit_file <- file.path(final_run, paste0(group, "__full"), "fit.rds")
    done_file <- file.path(final_run, paste0(group, "__full"), "completed.rds")
    config_file <- file.path(final_run, "run_config.rds")

    need(all(file.exists(c(fit_file, done_file, config_file))),
         paste("Missing final model files for", group))

    cfg <- readRDS(config_file)
    fit <- readRDS(fit_file)
    done <- readRDS(done_file)

    need(identical(done$status, "FIT_OK"),
         paste("Final model is not FIT_OK:", group))
    need(isTRUE(fit$checks$numerical_checks_ok),
         paste("Final numerical checks failed:", group))
    need(identical(fit$group, group) && identical(fit$variant, "full"),
         paste("Unexpected final fit:", group))

    prep <- cfg$prep
    input_file <- file.path(prep, paste0("input_", group, ".rds"))
    need(file.exists(input_file),
         paste("Prepared final input missing:", input_file))
    input <- readRDS(input_file)

    need(identical(fit$data_signature, input$data_signature),
         paste("Final fit/input signature mismatch:", group))

    list(fit = fit, input = input)
  }

  transform_group <- function(group) {

    src <- read_final(group)
    fit <- src$fit
    input <- src$input
    frame <- input$frame

    beta_old <- fit$beta
    V_old <- fit$V

    need(
      identical(names(beta_old), rownames(V_old)) &&
        identical(names(beta_old), colnames(V_old)),
      paste("Fixed covariance names do not match beta:", group)
    )

    required <- c(
      "(Intercept)",
      "year_decade",
      "onset_plasticity_z",
      "offset_plasticity_bio_z"
    )
    need(all(required %in% names(beta_old)),
         paste("Missing original fixed term:", group))

    ion <- interaction_name(names(beta_old), "onset_plasticity_z")
    ioff <- interaction_name(names(beta_old), "offset_plasticity_bio_z")

    # s converts current offset predictor to a common ADVANCEMENT orientation:
    # uni: current positive = advance, so s=+1
    # multi: current positive = delay, so s=-1
    s <- if (group == "univoltine") 1 else -1

    old_names <- names(beta_old)
    new_names <- c(
      "(Intercept)",
      "year_decade",
      "coordinated_advancement",
      "season_extension",
      "year_decade:coordinated_advancement",
      "year_decade:season_extension"
    )

    need(length(old_names) == 6L,
         paste("Expected six final fixed coefficients in", group))

    # Transformation beta_new = L %*% beta_old
    L <- matrix(
      0,
      nrow = 6L,
      ncol = 6L,
      dimnames = list(new_names, old_names)
    )

    L["(Intercept)", "(Intercept)"] <- 1
    L["year_decade", "year_decade"] <- 1

    L["coordinated_advancement", "onset_plasticity_z"] <- 1
    L["coordinated_advancement", "offset_plasticity_bio_z"] <- s

    L["season_extension", "onset_plasticity_z"] <- 1
    L["season_extension", "offset_plasticity_bio_z"] <- -s

    L["year_decade:coordinated_advancement", ion] <- 1
    L["year_decade:coordinated_advancement", ioff] <- s

    L["year_decade:season_extension", ion] <- 1
    L["year_decade:season_extension", ioff] <- -s

    beta_new <- as.numeric(L %*% beta_old)
    names(beta_new) <- rownames(L)

    V_new <- L %*% V_old %*% t(L)
    dimnames(V_new) <- list(new_names, new_names)

    need(all(is.finite(beta_new)) &&
           all(is.finite(V_new)) &&
           all(diag(V_new) > 0),
         paste("Invalid transformed estimates/covariance:", group))

    # Build population axes.
    x <- frame$onset_plasticity_z
    offset_adv <- s * frame$offset_plasticity_bio_z

    G <- (x + offset_adv) / 2
    E <- (x - offset_adv) / 2

    axes_rows <- data.frame(
      group = group,
      pop_id = as.character(frame$pop_id),
      SPECIES = as.character(frame$SPECIES),
      SITE_ID = as.character(frame$SITE_ID),
      YEAR = as.integer(frame$YEAR),
      year_decade = frame$year_decade,
      onset_advancement_z = x,
      offset_advancement_z = offset_adv,
      coordinated_advancement = G,
      season_extension = E,
      stringsAsFactors = FALSE
    )

    # Verify the exact fixed-effect linear predictor.
    X_old <- input$data$X
    need(
      identical(colnames(X_old), names(beta_old)),
      paste("Original X/beta columns differ:", group)
    )

    X_new <- cbind(
      "(Intercept)" = 1,
      "year_decade" = frame$year_decade,
      "coordinated_advancement" = G,
      "season_extension" = E,
      "year_decade:coordinated_advancement" =
        frame$year_decade * G,
      "year_decade:season_extension" =
        frame$year_decade * E
    )

    X_new <- X_new[, names(beta_new), drop = FALSE]

    eta_old <- as.numeric(X_old %*% beta_old)
    eta_new <- as.numeric(X_new %*% beta_new)
    max_error <- max(abs(eta_old - eta_new))

    need(is.finite(max_error) && max_error < 1e-10,
         paste("Reparameterized fixed predictor is not exact:", group,
               "max error =", signif(max_error, 6)))

    tab <- fixed_table(beta_new, V_new)
    tab$group <- group
    tab <- tab[, c("group", setdiff(names(tab), "group")), drop = FALSE]

    focal <- tab[
      tab$term %in% c(
        "year_decade:coordinated_advancement",
        "year_decade:season_extension"
      ),
      ,
      drop = FALSE
    ]

    focal$interpretation <- c(
      "Effect on population trend of +1 SD stronger onset advancement AND +1 SD stronger offset advancement, holding the extension axis constant.",
      "Effect on population trend of +1 SD stronger onset advancement AND 1 SD relatively later termination, holding coordinated advancement constant."
    )[match(
      focal$term,
      c(
        "year_decade:coordinated_advancement",
        "year_decade:season_extension"
      )
    )]

    summary <- data.frame(
      group = group,
      n_observations = nrow(frame),
      n_populations = length(unique(frame$pop_id)),
      onset_offset_advancement_cor =
        stats::cor(x, offset_adv),
      coordinated_extension_cor =
        stats::cor(G, E),
      coordinated_advancement_SD = stats::sd(G),
      season_extension_SD = stats::sd(E),
      fixed_predictor_max_abs_error = max_error,
      exact_reparameterization = TRUE,
      stringsAsFactors = FALSE
    )

    list(
      table = tab,
      focal = focal,
      axes = axes_rows,
      summary = summary,
      covariance = V_new,
      L = L
    )
  }

  ans <- lapply(groups, transform_group)
  names(ans) <- groups

  all_fixed <- do.call(rbind, lapply(ans, `[[`, "table"))
  focal <- do.call(rbind, lapply(ans, `[[`, "focal"))
  axes <- do.call(rbind, lapply(ans, `[[`, "axes"))
  summary <- do.call(rbind, lapply(ans, `[[`, "summary"))

  csv(all_fixed, file.path(out_dir, "global_balance_fixed_effects.csv"))
  csv(focal, file.path(out_dir, "global_balance_focal_interactions.csv"))
  csv(axes, file.path(out_dir, "population_global_balance_axes.csv"))
  csv(summary, file.path(out_dir, "global_balance_summary.csv"))

  for (g in groups) {
    V <- ans[[g]]$covariance
    csv(
      data.frame(term = rownames(V), V, check.names = FALSE),
      file.path(out_dir, paste0(g, "_global_balance_covariance.csv"))
    )
    csv(
      data.frame(
        new_term = rownames(ans[[g]]$L),
        ans[[g]]$L,
        check.names = FALSE
      ),
      file.path(out_dir, paste0(g, "_transformation_matrix.csv"))
    )
  }

  writeLines(c(
    "phenoIMPACT | Global onset-offset balance",
    "",
    "This is an exact linear reparameterization of the FINAL QC abundance models.",
    "No model was fitted and no likelihood was re-optimized.",
    "",
    "OFFSET orientation:",
    "  univoltine: positive = stronger advancement",
    "  multivoltine: original positive = delay, therefore sign is reversed so positive = stronger advancement",
    "",
    "Axes:",
    "  coordinated_advancement = (onset_advancement_z + offset_advancement_z) / 2",
    "  season_extension = (onset_advancement_z - offset_advancement_z) / 2",
    "",
    "Interpretation of +1 coordinated_advancement:",
    "  onset advancement +1 SD AND offset advancement +1 SD.",
    "",
    "Interpretation of +1 season_extension:",
    "  onset advancement +1 SD AND offset advancement -1 SD;",
    "  i.e. relatively later termination / longer season.",
    "",
    "Both axes are simultaneous predictors and together contain exactly the same",
    "fixed-effect information as the original onset + offset parameterization.",
    "",
    "OFFSET remains the onset-conditioned termination plasticity from the two-stage analysis.",
    "",
    capture.output(print(focal, row.names = FALSE))
  ), file.path(out_dir, "README_results.txt"))

  cat("\n==============================================================\n")
  cat("GLOBAL ONSET-OFFSET BALANCE\n")
  cat("EXACT REPARAMETERIZATION — NO MODEL REFITTED\n")
  cat("==============================================================\n\n")

  print(
    focal[, c(
      "group", "term", "estimate", "std.error",
      "conf.low", "conf.high", "p_Wald",
      "ratio_of_decadal_factors",
      "ratio_conf.low", "ratio_conf.high"
    )],
    row.names = FALSE
  )

  cat("\nExact-predictor checks:\n")
  print(summary, row.names = FALSE)

  cat("\nOutput: ",
      normalizePath(out_dir, winslash = "/", mustWork = TRUE),
      "\n", sep = "")

  invisible(list(
    focal = focal,
    fixed = all_fixed,
    axes = axes,
    summary = summary,
    output = out_dir
  ))
}

global_phenology_balance <- analyse_global_phenology_balance()
