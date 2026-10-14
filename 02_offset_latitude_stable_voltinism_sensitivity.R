# =============================================================================
# phenoIMPACT | Offset latitude sensitivity restricted to stable-voltinism spp.
# Version 1.0 | 2026-10-07
#
# Fast sensitivity: NO phenology or abundance model is refitted.
#
# Uses the already extracted population-level contextual plasticities and asks
# whether the within-species latitude + longitude pattern in OFFSET remains when
# species that switch voltinism category geographically are removed.
#
# Two restrictions are reported:
#   1) stable_across_analysed_networks  [primary sensitivity]
#   2) stable_across_full_lookup       [strict sensitivity]
#
# CIs are species-cluster bootstrap intervals across the fitted population
# plasticities. They are robustness intervals, not full first-stage model-based
# uncertainty intervals.
#
# RUN AFTER THE AUDIT:
#   source(file.choose(), encoding = "UTF-8")
#   stable_offset <- run_stable_voltinism_offset_sensitivity()
# =============================================================================

run_stable_voltinism_offset_sensitivity <- function(
  project_root = "E:/phenoIMPACT project/code/phenoIMPACT",
  latitude_run = file.path(
    project_root, "output", "phenology_plasticity",
    "latitude_from_extracted_plasticity_15km",
    "run_20261005_124834_165c6dd85b83"
  ),
  longitude_results = file.path(project_root, "data", "longitude_sensitivity_results.csv"),
  n_boot = 2000L,
  seed = 20261007L
) {
  need <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)
  csv <- function(x, p) utils::write.csv(x, p, row.names = FALSE, na = "")

  req <- c("dplyr", "ggplot2")
  miss <- req[!vapply(req, requireNamespace, logical(1), quietly = TRUE)]
  need(!length(miss), paste("Missing packages:", paste(miss, collapse = ", ")))

  project_root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
  latitude_run <- normalizePath(latitude_run, winslash = "/", mustWork = TRUE)
  need(file.exists(longitude_results), paste("Missing:", longitude_results))

  # Locate latest completed voltinism audit.
  audit_root <- file.path(project_root, "output", "phenology_plasticity", "voltinism_audit")
  need(dir.exists(audit_root), "Run 01_audit_voltinism_geography.R first.")
  audit_dirs <- list.dirs(audit_root, recursive = FALSE, full.names = TRUE)
  audit_dirs <- audit_dirs[file.exists(file.path(audit_dirs, "species_stable_for_sensitivity.csv"))]
  need(length(audit_dirs) > 0, "No completed voltinism audit found.")
  audit_dir <- audit_dirs[order(file.info(audit_dirs)$mtime, decreasing = TRUE)][1]
  stable_file <- file.path(audit_dir, "species_stable_for_sensitivity.csv")
  stable <- utils::read.csv(stable_file, stringsAsFactors = FALSE, check.names = FALSE)

  needed_stable <- c(
    "SPECIES", "analysis_stability_class", "source_stability_class",
    "stable_across_analysed_networks", "stable_across_full_lookup"
  )
  need(all(needed_stable %in% names(stable)), "Unexpected stable-species table schema.")

  coord_files <- list.files(file.path(latitude_run, "inputs", "coordinates"),
                            pattern = "\\.csv$", full.names = TRUE)
  need(length(coord_files) > 0, "No coordinate files in latitude run.")
  coords <- dplyr::bind_rows(lapply(coord_files, function(f) {
    x <- utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
    x[, c("SITE_ID", "bms_id", "longitude_deg", "latitude_deg")]
  })) |>
    dplyr::mutate(SITE_ID = as.character(SITE_ID), bms_id = as.character(bms_id)) |>
    dplyr::distinct()

  dup <- coords |>
    dplyr::count(SITE_ID, bms_id) |>
    dplyr::filter(n > 1)
  need(nrow(dup) == 0, "Coordinate table has duplicate site/BMS rows.")

  main <- utils::read.csv(longitude_results, stringsAsFactors = FALSE, check.names = FALSE)
  need(all(c("analysis", "term", "weighting", "specification",
             "estimate_raw_per_10deg", "lower_95_raw", "upper_95_raw") %in% names(main)),
       "Unexpected longitude-results schema.")

  analyses <- c("offset_univoltine", "offset_multivoltine")
  expected_class <- c(
    offset_univoltine = "stable_univoltine",
    offset_multivoltine = "stable_multivoltine"
  )
  orient <- c(offset_univoltine = -1, offset_multivoltine = 1)

  fit_projection <- function(d, weighting = c("equal_species", "equal_population")) {
    weighting <- match.arg(weighting)
    xlat <- d$latitude_within_species_deg / 10
    xlon <- d$longitude_within_species_deg / 10
    y <- d$response_bio_centered
    X <- cbind(latitude = xlat, longitude = xlon)
    nsp <- table(d$SPECIES)
    w <- if (weighting == "equal_population") rep(1, nrow(d)) else
      1 / as.numeric(nsp[d$SPECIES])
    A <- crossprod(X, w * X)
    b <- solve(A, crossprod(X, w * y))
    setNames(as.numeric(b), colnames(X))
  }

  bootstrap_species <- function(d, weighting, n_boot, seed) {
    set.seed(seed)
    spp <- unique(d$SPECIES)
    parts <- lapply(spp, function(s) {
      z <- d[d$SPECIES == s, , drop = FALSE]
      X <- cbind(
        latitude = z$latitude_within_species_deg / 10,
        longitude = z$longitude_within_species_deg / 10
      )
      y <- z$response_bio_centered
      w <- if (weighting == "equal_population") rep(1, nrow(z)) else rep(1/nrow(z), nrow(z))
      list(A = crossprod(X, w * X), c = crossprod(X, w * y))
    })
    ans <- matrix(NA_real_, n_boot, 2, dimnames = list(NULL, c("latitude", "longitude")))
    for (i in seq_len(n_boot)) {
      draw <- sample.int(length(spp), length(spp), replace = TRUE)
      A <- Reduce(`+`, lapply(draw, function(j) parts[[j]]$A))
      cc <- Reduce(`+`, lapply(draw, function(j) parts[[j]]$c))
      if (qr(A)$rank == 2) ans[i, ] <- solve(A, cc)
    }
    ans[stats::complete.cases(ans), , drop = FALSE]
  }

  out_root <- file.path(latitude_run, "stable_voltinism_offset_sensitivity")
  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
  out <- file.path(out_root, paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  dir.create(out, recursive = TRUE, showWarnings = FALSE)

  all_results <- list()
  all_data <- list()

  restrictions <- c(
    analysed_stable = "stable_across_analysed_networks",
    full_lookup_stable = "stable_across_full_lookup"
  )

  for (a in analyses) {
    pop_file <- file.path(latitude_run, a, "population_plasticity_with_latitude.csv")
    need(file.exists(pop_file), paste("Missing:", pop_file))
    p0 <- utils::read.csv(pop_file, stringsAsFactors = FALSE, check.names = FALSE)
    p0$SPECIES <- as.character(p0$SPECIES)
    p0$SITE_ID <- as.character(p0$SITE_ID)
    p0$bms_id <- as.character(p0$bms_id)

    need(all(c("SPECIES", "SITE_ID", "bms_id", "plasticity_within_species_raw") %in% names(p0)),
         paste("Unexpected population plasticity schema:", a))

    p0 <- dplyr::left_join(
      p0,
      coords,
      by = c("SITE_ID", "bms_id"),
      relationship = "many-to-one",
      suffix = c("", "_coord")
    )
    need(!anyNA(p0$longitude_deg), paste("Missing longitude:", a))

    for (rname in names(restrictions)) {
      flag <- restrictions[[rname]]
      class_col <- if (rname == "analysed_stable") "analysis_stability_class" else "source_stability_class"

      keep <- stable[[flag]] & stable[[class_col]] == expected_class[[a]]
      spp_keep <- stable$SPECIES[keep]
      d <- p0[p0$SPECIES %in% spp_keep, , drop = FALSE]
      need(nrow(d) > 0, paste("No stable species retained for", a, rname))
      need(length(unique(d$SPECIES)) >= 5,
           paste("Too few stable species for", a, rname))

      # Re-centre geography and response within species after restriction.
      d <- d |>
        dplyr::group_by(SPECIES) |>
        dplyr::mutate(
          latitude_within_species_deg = latitude_deg - mean(latitude_deg),
          longitude_within_species_deg = longitude_deg - mean(longitude_deg),
          response_bio = orient[[a]] * plasticity_within_species_raw,
          response_bio_centered = response_bio - mean(response_bio)
        ) |>
        dplyr::ungroup()

      for (wt in c("equal_species", "equal_population")) {
        b <- fit_projection(d, wt)
        boot <- bootstrap_species(d, wt, n_boot,
                                  seed + match(a, analyses)*100 + match(rname,names(restrictions))*10 + match(wt,c("equal_species","equal_population")))
        need(nrow(boot) >= .9*n_boot, "Too many rank-deficient bootstrap replicates.")
        for (term in c("latitude", "longitude")) {
          vals <- boot[, term]
          all_results[[length(all_results)+1L]] <- data.frame(
            analysis = a,
            restriction = rname,
            weighting = wt,
            term = term,
            n_species = length(unique(d$SPECIES)),
            n_populations = nrow(d),
            estimate_biological_per_10deg = unname(b[term]),
            lower_95_species_bootstrap = unname(stats::quantile(vals, .025)),
            upper_95_species_bootstrap = unname(stats::quantile(vals, .975)),
            bootstrap_prop_positive = mean(vals > 0),
            bootstrap_prop_negative = mean(vals < 0),
            stringsAsFactors = FALSE
          )
        }
      }
      d$analysis <- a
      d$restriction <- rname
      all_data[[length(all_data)+1L]] <- d
    }
  }

  results <- dplyr::bind_rows(all_results)
  dat <- dplyr::bind_rows(all_data)

  # Main reviewed full-set reference, biologically oriented.
  ref <- main |>
    dplyr::filter(
      analysis %in% analyses,
      specification == "latitude_plus_longitude_within_species",
      weighting == "equal_population",
      term %in% c("latitude", "longitude")
    ) |>
    dplyr::mutate(
      orientation = orient[analysis],
      full_set_estimate_biological = orientation * estimate_raw_per_10deg,
      full_set_lower_biological = pmin(orientation * lower_95_raw, orientation * upper_95_raw),
      full_set_upper_biological = pmax(orientation * lower_95_raw, orientation * upper_95_raw)
    ) |>
    dplyr::select(analysis, term, full_set_estimate_biological,
                  full_set_lower_biological, full_set_upper_biological)

  comparison <- results |>
    dplyr::left_join(ref, by = c("analysis", "term")) |>
    dplyr::mutate(
      estimate_difference_from_full = estimate_biological_per_10deg - full_set_estimate_biological
    )

  csv(results, file.path(out, "stable_voltinism_geographic_results.csv"))
  csv(comparison, file.path(out, "comparison_with_full_set.csv"))
  csv(dat, file.path(out, "retained_population_values.csv"))
  csv(stable, file.path(out, "species_stability_used.csv"))

  # Figure for primary analysed-stable / equal-species result.
  figdata <- dat |>
    dplyr::filter(restriction == "analysed_stable")
  panels <- list()
  for (a in analyses) {
    z <- figdata[figdata$analysis == a, , drop = FALSE]
    rr <- results |>
      dplyr::filter(analysis == a, restriction == "analysed_stable",
                    weighting == "equal_species", term == "latitude")
    need(nrow(rr) == 1, "Primary stable result missing.")
    slope <- rr$estimate_biological_per_10deg
    lo <- rr$lower_95_species_bootstrap
    hi <- rr$upper_95_species_bootstrap
    xx <- seq(stats::quantile(z$latitude_within_species_deg,.005),
              stats::quantile(z$latitude_within_species_deg,.995),length.out=250)
    line <- data.frame(
      x=xx,fit=slope*xx/10,
      low=lo*xx/10,high=hi*xx/10
    )
    line$ymin <- pmin(line$low,line$high); line$ymax <- pmax(line$low,line$high)
    title <- if (a == "offset_univoltine") "Offset - stable univoltine species" else
      "Offset - stable multivoltine species"
    ylab <- if (a == "offset_univoltine") "Relative offset advancement" else "Relative offset delay"
    panels[[a]] <- ggplot2::ggplot(z, ggplot2::aes(latitude_within_species_deg,response_bio_centered)) +
      ggplot2::geom_hline(yintercept=0,linetype="dashed",colour="grey45") +
      ggplot2::geom_vline(xintercept=0,linetype="dashed",colour="grey45") +
      ggplot2::geom_point(colour="black",alpha=.12,size=.55) +
      ggplot2::geom_ribbon(data=line,ggplot2::aes(x=x,ymin=ymin,ymax=ymax),
                           inherit.aes=FALSE,fill="#4C78A8",alpha=.25) +
      ggplot2::geom_line(data=line,ggplot2::aes(x=x,y=fit),inherit.aes=FALSE,
                         colour="#2F5AA8",linewidth=1.1) +
      ggplot2::annotate("label",x=-Inf,y=Inf,
                        label=sprintf("%.2f [%.2f, %.2f] per 10°",slope,lo,hi),
                        hjust=-.05,vjust=1.15,family="Garamond",size=4) +
      ggplot2::theme_classic(base_size=15,base_family="Garamond") +
      ggplot2::theme(panel.border=ggplot2::element_rect(fill=NA,colour="black",linewidth=.5),
                     plot.title=ggplot2::element_text(hjust=.5,face="bold")) +
      ggplot2::labs(title=title,x="Latitude within species (°)",y=ylab)
  }
  for (a in names(panels)) {
    ggplot2::ggsave(file.path(out,paste0(a,"_stable_voltinism.png")),panels[[a]],
                    width=7,height=5,dpi=350,bg="white")
  }

  writeLines(c(
    "Stable-voltinism sensitivity of offset geographic plasticity patterns.",
    paste("Voltinism audit:", audit_dir),
    paste("Latitude extraction run:", latitude_run),
    "No phenology or abundance model was refitted.",
    "Species classified as geographically variable in voltinism are excluded.",
    "Primary sensitivity: stable across BMS units represented in analysed offset populations.",
    "Strict sensitivity: stable across the complete species_country_voltinism lookup.",
    "Biological orientation: positive = more offset advancement (univoltine) or more offset delay (multivoltine).",
    "Latitude and longitude are re-centred within each retained species.",
    "Intervals are species-cluster bootstrap robustness intervals and do not propagate full first-stage fixed-effect covariance."
  ), file.path(out,"README.txt"))

  if(requireNamespace("zip",quietly=TRUE)) {
    rel <- list.files(out,recursive=TRUE,full.names=FALSE)
    rel <- rel[grepl("\\.(csv|txt|png)$",rel,ignore.case=TRUE)]
    zip::zipr(file.path(out,"results_to_review.zip"),rel,root=out,mode="mirror",include_directories=FALSE)
  }

  cat("\nPrimary stable-voltinism sensitivity (analysed networks, equal species):\n")
  print(results |>
          dplyr::filter(restriction=="analysed_stable",weighting=="equal_species",term=="latitude"),
        row.names=FALSE)
  cat("\nOutput:\n",normalizePath(out,winslash="/"),"\n",sep="")
  if(file.exists(file.path(out,"results_to_review.zip")))
    cat("FILE TO SEND:\n",normalizePath(file.path(out,"results_to_review.zip"),winslash="/"),"\n",sep="")

  invisible(list(output=out,results=results,comparison=comparison,
                 zip=file.path(out,"results_to_review.zip")))
}
