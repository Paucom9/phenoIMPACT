# =============================================================================
# phenoIMPACT | Audit ONLY two anomalous phenology curves
# Uses the same GAM + PELT settings as pheno_gams.R
#
# Cases:
#   Polygonia c-album | UKBMS.2537 | 2008
#   Pontia daplidice   | ES-CTBMS.9 | 2013
#
# This script does NOT rerun the full phenology workflow and does NOT change data.
# It reads the raw eBMS count/visit files, reconstructs only these two curves,
# exports two PNGs and a small CSV summary.
# =============================================================================

audit_two_pheno_curves <- function(
  project_root = "E:/phenoIMPACT project/code/phenoIMPACT",
  out_dir = file.path(project_root, "output", "diagnostics", "two_offset_curve_audit")
) {

  pkgs <- c("data.table", "dplyr", "lubridate", "mgcv", "changepoint", "splus2R", "ggplot2")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Missing packages: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  count_file <- file.path(project_root, "data", "ebms_count.csv")
  visit_file <- file.path(project_root, "data", "ebms_visit.csv")

  if (!file.exists(count_file)) stop("Missing: ", count_file, call. = FALSE)
  if (!file.exists(visit_file)) stop("Missing: ", visit_file, call. = FALSE)

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  cases <- data.frame(
    SPECIES = c("Polygonia c-album", "Pontia daplidice"),
    SITE_ID = c("UKBMS.2537", "ES-CTBMS.9"),
    YEAR = c(2008L, 2013L),
    stringsAsFactors = FALSE
  )

  # Read only required columns where possible.
  counts <- data.table::fread(
    count_file,
    select = c("transect_id", "visit_date", "species_name", "count", "year"),
    data.table = FALSE,
    showProgress = FALSE
  )
  visits <- data.table::fread(
    visit_file,
    select = c("transect_id", "visit_date", "visit_id", "year"),
    data.table = FALSE,
    showProgress = FALSE
  )

  names(counts)[match(c("transect_id","visit_date","species_name","count","year"),
                      names(counts))] <-
    c("SITE_ID","DATE","SPECIES","COUNT","YEAR")
  names(visits)[match(c("transect_id","visit_date","visit_id","year"),
                      names(visits))] <-
    c("SITE_ID","DATE","visit_id","YEAR")

  counts$SITE_ID <- as.character(counts$SITE_ID)
  counts$SPECIES <- as.character(counts$SPECIES)
  counts$DATE <- as.Date(counts$DATE)
  counts$YEAR <- as.integer(as.character(counts$YEAR))

  visits$SITE_ID <- as.character(visits$SITE_ID)
  visits$DATE <- as.Date(visits$DATE)
  visits$YEAR <- as.integer(as.character(visits$YEAR))

  # Keep only the two target site-years/species.
  keep_counts <- Reduce(`|`, lapply(seq_len(nrow(cases)), function(i) {
    counts$SPECIES == cases$SPECIES[i] &
      counts$SITE_ID == cases$SITE_ID[i] &
      counts$YEAR == cases$YEAR[i]
  }))
  keep_visits <- Reduce(`|`, lapply(seq_len(nrow(cases)), function(i) {
    visits$SITE_ID == cases$SITE_ID[i] &
      visits$YEAR == cases$YEAR[i]
  }))

  counts <- counts[keep_counts, , drop = FALSE]
  visits <- visits[keep_visits, , drop = FALSE]

  if (!nrow(counts)) stop("No target count rows found.", call. = FALSE)
  if (!nrow(visits)) stop("No target visit rows found.", call. = FALSE)

  find_peaks_same_as_original <- function(x, ignore_threshold = 0.2,
                                          span = 11, strict = TRUE) {
    x[!is.finite(x)] <- min(x, na.rm = TRUE)
    pks <- splus2R::peaks(x = x, span = span, strict = strict)
    max_x <- max(x, na.rm = TRUE)
    which(ifelse(x > ignore_threshold * max_x, pks, FALSE))
  }

  audit_one <- function(species_pick, site_pick, year_pick) {

    sub_count <- counts[
      counts$SPECIES == species_pick &
        counts$SITE_ID == site_pick &
        counts$YEAR == year_pick, , drop = FALSE
    ]
    sub_visit <- visits[
      visits$SITE_ID == site_pick &
        visits$YEAR == year_pick, , drop = FALSE
    ]

    if (!nrow(sub_count)) stop("No counts for ", species_pick, " / ", site_pick,
                               " / ", year_pick, call. = FALSE)
    if (!nrow(sub_visit)) stop("No visits for ", site_pick, " / ", year_pick,
                               call. = FALSE)

    # Exact logic used in pheno_gams.R:
    # zero-fill visits at which the species had no count record.
    missing_dates <- sub_visit[!sub_visit$DATE %in% sub_count$DATE, , drop = FALSE]

    all_counts <- dplyr::bind_rows(
      data.frame(
        julian_day = lubridate::yday(sub_count$DATE),
        COUNT = as.numeric(sub_count$COUNT)
      ),
      data.frame(
        julian_day = lubridate::yday(missing_dates$DATE),
        COUNT = 0
      )
    )

    # Structural zeros outside the main monitoring season.
    all_counts <- dplyr::bind_rows(
      all_counts,
      data.frame(julian_day = c(1:30, 335:365), COUNT = 0)
    )

    gam_model <- mgcv::gam(
      COUNT ~ s(julian_day),
      data = all_counts,
      family = mgcv::nb()
    )

    julian_days <- 1:365
    pred <- as.numeric(stats::predict(
      gam_model,
      newdata = data.frame(julian_day = julian_days),
      type = "response"
    ))

    cp_mean <- changepoint::cpt.mean(
      pred, method = "PELT", penalty = "Manual", pen.value = 2
    )
    cp_var <- changepoint::cpt.var(
      pred, method = "PELT", penalty = "Manual", pen.value = 0.05
    )

    cps_mean <- changepoint::cpts(cp_mean)
    cps_var <- changepoint::cpts(cp_var)

    peaks <- find_peaks_same_as_original(pred)
    first_peak <- if (length(peaks)) min(peaks) else NA_integer_
    last_peak <- if (length(peaks)) max(peaks) else NA_integer_
    onset_mean <- if (length(cps_mean)) min(cps_mean) else NA_integer_
    offset_mean <- if (length(cps_mean)) max(cps_mean) else NA_integer_
    onset_var <- if (length(cps_var)) min(cps_var) else NA_integer_
    offset_var <- if (length(cps_var)) max(cps_var) else NA_integer_

    positive <- sub_count[is.finite(sub_count$COUNT) & sub_count$COUNT > 0, , drop = FALSE]
    first_obs <- if (nrow(positive)) min(lubridate::yday(positive$DATE)) else NA_integer_
    last_obs <- if (nrow(positive)) max(lubridate::yday(positive$DATE)) else NA_integer_

    curve <- data.frame(Julian_Day = julian_days, Predicted_Count = pred)

    p <- ggplot2::ggplot() +
      ggplot2::geom_line(
        data = curve,
        ggplot2::aes(Julian_Day, Predicted_Count),
        linewidth = 1.25
      ) +
      ggplot2::geom_point(
        data = all_counts[all_counts$julian_day > 30 & all_counts$julian_day < 335, ],
        ggplot2::aes(julian_day, COUNT),
        alpha = 0.55, size = 1.8
      ) +
      ggplot2::geom_vline(
        xintercept = c(onset_mean, offset_mean),
        linetype = "dashed", linewidth = 0.9
      ) +
      ggplot2::geom_vline(
        xintercept = c(onset_var, offset_var),
        linetype = "dotted", linewidth = 0.8
      ) +
      ggplot2::geom_vline(
        xintercept = last_obs,
        linetype = "dotdash", linewidth = 0.8
      ) +
      ggplot2::geom_point(
        data = curve[curve$Julian_Day %in% peaks, , drop = FALSE],
        ggplot2::aes(Julian_Day, Predicted_Count),
        size = 2.8
      ) +
      ggplot2::annotate(
        "text", x = offset_mean, y = Inf,
        label = paste0("mean offset = ", offset_mean),
        angle = 90, vjust = 1.3, hjust = 1.05, size = 3.5
      ) +
      ggplot2::annotate(
        "text", x = last_peak, y = Inf,
        label = paste0("last peak = ", last_peak),
        angle = 90, vjust = -0.3, hjust = 1.05, size = 3.5
      ) +
      ggplot2::labs(
        title = species_pick,
        subtitle = paste0(site_pick, " | ", year_pick,
                          " | mean offset ", offset_mean,
                          " | variance offset ", offset_var,
                          " | last peak ", last_peak,
                          " | last positive obs. ", last_obs),
        x = "Day of year",
        y = "Predicted / observed abundance",
        caption = paste(
          "Dashed = PELT mean onset/offset;",
          "dotted = PELT variance onset/offset;",
          "dot-dash = last positive observation;",
          "points on curve = detected peaks."
        )
      ) +
      ggplot2::theme_classic(base_size = 13)

    safe_sp <- gsub("[^A-Za-z0-9]+", "_", species_pick)
    safe_site <- gsub("[^A-Za-z0-9]+", "_", site_pick)

    fig_path <- file.path(
      out_dir,
      paste0("audit_", safe_sp, "_", safe_site, "_", year_pick, ".png")
    )
    ggplot2::ggsave(fig_path, p, width = 9, height = 5.5, dpi = 300, bg = "white")

    curve_path <- file.path(
      out_dir,
      paste0("curve_", safe_sp, "_", safe_site, "_", year_pick, ".csv")
    )
    utils::write.csv(curve, curve_path, row.names = FALSE)

    data.frame(
      SPECIES = species_pick,
      SITE_ID = site_pick,
      YEAR = year_pick,
      n_visits = length(unique(sub_visit$visit_id)),
      n_positive_dates = length(unique(positive$DATE)),
      first_positive_obs = first_obs,
      last_positive_obs = last_obs,
      n_detected_peaks = length(peaks),
      first_peak = first_peak,
      last_peak = last_peak,
      onset_mean = onset_mean,
      offset_mean = offset_mean,
      onset_var = onset_var,
      offset_var = offset_var,
      last_peak_minus_offset_mean = last_peak - offset_mean,
      last_obs_minus_offset_mean = last_obs - offset_mean,
      figure = fig_path,
      stringsAsFactors = FALSE
    )
  }

  ans <- do.call(rbind, lapply(seq_len(nrow(cases)), function(i) {
    message("Auditing: ", cases$SPECIES[i], " | ", cases$SITE_ID[i], " | ", cases$YEAR[i])
    audit_one(cases$SPECIES[i], cases$SITE_ID[i], cases$YEAR[i])
  }))

  utils::write.csv(ans, file.path(out_dir, "audit_summary.csv"), row.names = FALSE)

  message("\nDONE. Only two curves were fitted.")
  message("Output: ", normalizePath(out_dir, winslash = "/", mustWork = TRUE))
  print(ans, row.names = FALSE)

  invisible(ans)
}

two_curve_audit <- audit_two_pheno_curves()
