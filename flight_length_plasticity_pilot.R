# =============================================================================
# phenoIMPACT | Direct FLIGHT-LENGTH plasticity pilot
# Version 1.0 | 2026-10-08
#
# PURPOSE
# Test directly whether warmer-than-usual years lengthen or shorten the annual
# flight period, and extract population-level contextual flight-length
# plasticity before any spatial refit or abundance-trend model.
#
# IMPORTANT:
# FLIGHT_LENGTH_mean is a duration, not a day-of-year. Therefore the existing
# climate_variables.R phenophase code must NOT be applied directly to
# FLIGHT_LENGTH_mean (it would interpret e.g. "120 days long" as DOY 120).
#
# Thermal exposure used here:
#   - For each SPECIES x SITE_ID population, define a FIXED "typical flight
#     season" from its mean QC-valid ONSET_mean to mean QC-valid OFFSET_mean.
#   - For each calendar year, calculate mean temperature over this fixed window.
#   - Annual anomaly = that annual seasonal mean minus the population's
#     long-term mean seasonal temperature.
#
# Fixed population-specific windows avoid annual circularity: the thermal window
# does not move with that year's observed flight length.
#
# Context moderators:
#   background temperature = mean annual fixed-season temperature
#   warming trend          = linear trend in fixed-season temp, deg C / decade
#   predictability         = - mean squared daily departure from site x DOY
#                            climatology within the fixed typical flight season
#   photoperiod            = mean daylength within the fixed typical season
#
# QC:
#   - finite onset, offset and flight length
#   - verify FLIGHT_LENGTH_mean == OFFSET_mean - ONSET_mean + 1
#   - exclude only known offset failure: finite LAST_PEAK & OFFSET_mean < LAST_PEAK
#   - no new phenology outlier filter
#
# PILOT MODEL (separate univoltine / multivoltine):
# FLIGHT_LENGTH_mean ~ season_temp_anomaly_C *
#   (season_photoperiod_z + season_background_z +
#    season_predictability_z + season_trend_z) +
#   (1 | SITE_ID) + (1 | SPECIES) +
#   (0 + season_temp_anomaly_C | SPECIES_slope) +
#   (1 | site_year_id)
#
# Gaussian ML. NO spatial field yet. This pilot is to validate the new direct
# metric and inspect the signal before starting expensive 15-km spatial fits.
#
# Positive population plasticity = flight season becomes LONGER in warmer years.
#
# Run:
#   source(file.choose(), encoding = "UTF-8")
# =============================================================================

run_flight_length_plasticity_pilot <- function(
  root = "E:/phenoIMPACT project/code/phenoIMPACT",
  min_climate_coverage = 0.90,
  out_dir = file.path(
    root, "output", "phenology_plasticity",
    "flight_length_direct_pilot"
  )
) {

  pkgs <- c("data.table", "dplyr", "lme4", "sf", "digest")
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) {
    stop("Missing packages: ", paste(miss, collapse = ", "), call. = FALSE)
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  need <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)
  write_csv <- function(x, path) {
    utils::write.csv(x, path, row.names = FALSE, na = "")
  }
  zscale <- function(x) {
    s <- stats::sd(x, na.rm = TRUE)
    need(is.finite(s) && s > 0, "Cannot standardize constant/invalid variable.")
    (x - mean(x, na.rm = TRUE)) / s
  }

  pheno_file <- file.path(root, "output", "pheno_estimates_allspp.csv")
  temp_file <- file.path(root, "output", "climate", "daily_temperature.csv")
  coord_file <- file.path(root, "data", "ebms_transect_coord.csv")
  volt_file <- file.path(root, "data", "voltinism", "species_country_voltinism.csv")

  required_files <- c(pheno_file, temp_file, coord_file, volt_file)
  need(all(file.exists(required_files)),
       paste("Missing input file(s):\n",
             paste(required_files[!file.exists(required_files)], collapse = "\n")))

  message("Reading phenology...")
  ph <- data.table::fread(pheno_file, showProgress = FALSE)

  required_pheno <- c(
    "SPECIES", "SITE_ID", "YEAR",
    "ONSET_mean", "OFFSET_mean", "FLIGHT_LENGTH_mean",
    "LAST_PEAK"
  )
  need(all(required_pheno %in% names(ph)),
       "pheno_estimates_allspp.csv lacks required columns.")

  ph[, SPECIES := as.character(SPECIES)]
  ph[, SITE_ID := as.character(SITE_ID)]
  ph[, YEAR := as.integer(YEAR)]
  ph[, bms_id := sub("\\..*$", "", SITE_ID)]

  finite_core <- is.finite(ph$ONSET_mean) &
    is.finite(ph$OFFSET_mean) &
    is.finite(ph$FLIGHT_LENGTH_mean)

  ph_core <- ph[finite_core]
  need(nrow(ph_core) > 0L, "No finite flight-length phenology rows.")

  identity_error <- max(abs(
    ph_core$FLIGHT_LENGTH_mean -
      (ph_core$OFFSET_mean - ph_core$ONSET_mean + 1)
  ))
  need(is.finite(identity_error) && identity_error < 1e-8,
       paste("FLIGHT_LENGTH_mean identity failed; max error =", identity_error))

  invalid_offset <- is.finite(ph_core$LAST_PEAK) &
    ph_core$OFFSET_mean < ph_core$LAST_PEAK

  qc_summary <- data.frame(
    n_finite_core = nrow(ph_core),
    n_removed_offset_before_last_peak = sum(invalid_offset),
    percent_removed = 100 * mean(invalid_offset),
    max_flight_length_identity_error = identity_error
  )

  write_csv(
    as.data.frame(ph_core[invalid_offset]),
    file.path(out_dir, "excluded_offset_before_last_peak.csv")
  )
  write_csv(qc_summary, file.path(out_dir, "phenology_qc_summary.csv"))

  phq <- ph_core[!invalid_offset]
  need(all(phq$OFFSET_mean >= phq$ONSET_mean),
       "Some retained offsets precede onset.")

  message("Preparing voltinism lookup...")
  vol <- data.table::fread(volt_file, sep = ";", showProgress = FALSE)
  need("SPECIES" %in% names(vol), "Voltinism file lacks SPECIES.")

  vol_long <- data.table::melt(
    vol,
    id.vars = "SPECIES",
    variable.name = "bms_id",
    value.name = "voltinism",
    variable.factor = FALSE
  )
  vol_long[, bms_id := data.table::fcase(
    bms_id == "ES.CTBMS", "ES-CTBMS",
    bms_id == "ES.ZEBMS", "ES-ZEBMS",
    default = bms_id
  )]
  vol_long[, SPECIES := as.character(SPECIES)]
  vol_long[, voltinism := as.character(voltinism)]
  vol_long <- unique(vol_long[, .(SPECIES, bms_id, voltinism)])

  dupv <- vol_long[, .N, by = .(SPECIES, bms_id)][N > 1L]
  need(nrow(dupv) == 0L, "Conflicting voltinism entries.")

  phq <- merge(
    phq, vol_long,
    by = c("SPECIES", "bms_id"),
    all.x = TRUE, sort = FALSE
  )

  phq <- phq[voltinism %in% c("univoltine", "multivoltine")]
  need(nrow(phq) > 0L, "No uni/multivoltine rows after voltinism join.")

  message("Defining fixed typical flight seasons...")
  windows <- phq[
    ,
    .(
      n_valid_pheno_years = data.table::uniqueN(YEAR),
      typical_onset = mean(ONSET_mean),
      typical_offset = mean(OFFSET_mean),
      first_pheno_year = min(YEAR),
      last_pheno_year = max(YEAR)
    ),
    by = .(SPECIES, SITE_ID, bms_id, voltinism)
  ]

  windows[, start_day := pmax(1L, as.integer(ceiling(typical_onset)))]
  windows[, end_day := pmin(365L, as.integer(floor(typical_offset)))]
  windows[, expected_days := end_day - start_day + 1L]

  need(all(windows$expected_days > 0L),
       "At least one typical flight window is empty.")

  message("Preparing coordinates and photoperiod...")
  coords <- data.table::fread(coord_file, showProgress = FALSE)

  need(all(c(
    "transect_id", "transect_lon", "transect_lat"
  ) %in% names(coords)), "Coordinate file lacks required columns.")

  coords <- unique(coords[
    is.finite(transect_lon) & is.finite(transect_lat),
    .(
      SITE_ID = as.character(transect_id),
      x_3035 = as.numeric(transect_lon),
      y_3035 = as.numeric(transect_lat)
    )
  ])

  need(!anyDuplicated(coords$SITE_ID), "Duplicate/conflicting site coordinates.")

  coords_sf <- sf::st_as_sf(
    as.data.frame(coords),
    coords = c("x_3035", "y_3035"),
    crs = 3035,
    remove = FALSE
  )
  coords_wgs <- sf::st_transform(coords_sf, 4326)
  xy <- sf::st_coordinates(coords_wgs)

  coords$latitude <- xy[, 2]
  coords$x_km <- coords$x_3035 / 1000
  coords$y_km <- coords$y_3035 / 1000

  windows <- merge(
    windows,
    coords[, .(SITE_ID, latitude, x_km, y_km)],
    by = "SITE_ID",
    all.x = TRUE,
    sort = FALSE
  )

  need(all(is.finite(windows$latitude)) &&
         all(is.finite(windows$x_km)) &&
         all(is.finite(windows$y_km)),
       "Some flight-length populations lack valid coordinates.")

  daylength <- function(lat, doy) {
    lat_rad <- lat * pi / 180
    decl <- 23.44 * pi / 180 *
      sin(2 * pi * (doy - 80) / 365)
    xx <- -tan(lat_rad) * tan(decl)
    xx <- pmin(pmax(xx, -1), 1)
    24 * acos(xx) / pi
  }

  windows[, season_photoperiod_h := vapply(
    seq_len(.N),
    function(i) {
      mean(daylength(
        latitude[i],
        seq.int(start_day[i], end_day[i])
      ))
    },
    numeric(1)
  )]

  message("Reading daily temperature and building cumulative climate tables...")
  temp <- data.table::fread(
    temp_file,
    select = c("transect_id", "year", "julian_day", "temp"),
    showProgress = FALSE
  )

  data.table::setnames(
    temp,
    c("transect_id", "year"),
    c("SITE_ID", "YEAR")
  )
  temp[, SITE_ID := as.character(SITE_ID)]
  temp[, YEAR := as.integer(YEAR)]
  temp[, julian_day := as.integer(julian_day)]

  valid_sites <- unique(windows$SITE_ID)
  temp <- temp[
    SITE_ID %in% valid_sites &
      julian_day >= 1L & julian_day <= 365L
  ]
  need(nrow(temp) > 0L, "No daily temperature rows for retained sites.")

  temp[, clim_doy := mean(temp, na.rm = TRUE),
       by = .(SITE_ID, julian_day)]
  temp[, daily_resid := temp - clim_doy]

  temp[, temp_ok := as.integer(is.finite(temp))]
  temp[, temp_value := data.table::fifelse(is.finite(temp), temp, 0)]
  temp[, resid_ok := as.integer(is.finite(daily_resid))]
  temp[, resid2_value := data.table::fifelse(
    is.finite(daily_resid), daily_resid^2, 0
  )]

  data.table::setorder(temp, SITE_ID, YEAR, julian_day)

  temp[
    ,
    `:=`(
      cs_temp = cumsum(temp_value),
      cs_n = cumsum(temp_ok),
      cs_resid2 = cumsum(resid2_value),
      cs_resid_n = cumsum(resid_ok)
    ),
    by = .(SITE_ID, YEAR)
  ]

  cum <- temp[, .(
    SITE_ID, YEAR, julian_day,
    cs_temp, cs_n, cs_resid2, cs_resid_n
  )]

  zero <- unique(cum[, .(SITE_ID, YEAR)])
  zero[, `:=`(
    julian_day = 0L,
    cs_temp = 0,
    cs_n = 0L,
    cs_resid2 = 0,
    cs_resid_n = 0L
  )]
  cum <- data.table::rbindlist(list(zero, cum), use.names = TRUE)
  data.table::setkey(cum, SITE_ID, YEAR, julian_day)

  site_years <- unique(cum[YEAR > 0, .(SITE_ID, YEAR)])

  q <- merge(
    windows,
    site_years,
    by = "SITE_ID",
    allow.cartesian = TRUE,
    sort = FALSE
  )
  q[, row_id := .I]

  end_q <- q[, .(
    row_id,
    SITE_ID,
    YEAR,
    julian_day = end_day
  )]
  start_q <- q[, .(
    row_id,
    SITE_ID,
    YEAR,
    julian_day = start_day - 1L
  )]

  end_v <- cum[
    end_q,
    on = .(SITE_ID, YEAR, julian_day),
    roll = Inf
  ][, .(
    row_id,
    end_temp = cs_temp,
    end_n = cs_n,
    end_resid2 = cs_resid2,
    end_resid_n = cs_resid_n
  )]

  start_v <- cum[
    start_q,
    on = .(SITE_ID, YEAR, julian_day),
    roll = Inf
  ][, .(
    row_id,
    start_temp = cs_temp,
    start_n = cs_n,
    start_resid2 = cs_resid2,
    start_resid_n = cs_resid_n
  )]

  annual <- merge(q, end_v, by = "row_id", sort = FALSE)
  annual <- merge(annual, start_v, by = "row_id", sort = FALSE)

  annual[, n_temp_days := end_n - start_n]
  annual[, n_resid_days := end_resid_n - start_resid_n]
  annual[, temp_sum := end_temp - start_temp]
  annual[, resid2_sum := end_resid2 - start_resid2]
  annual[, climate_coverage := n_temp_days / expected_days]
  annual[, season_temp_C := data.table::fifelse(
    n_temp_days > 0,
    temp_sum / n_temp_days,
    NA_real_
  )]

  climate_qc <- annual[
    ,
    .(
      n_population_years = .N,
      prop_coverage_ge_threshold =
        mean(climate_coverage >= min_climate_coverage, na.rm = TRUE),
      min_coverage = min(climate_coverage, na.rm = TRUE),
      median_coverage = stats::median(climate_coverage, na.rm = TRUE)
    )
  ]
  write_csv(
    as.data.frame(climate_qc),
    file.path(out_dir, "climate_coverage_summary.csv")
  )

  annual <- annual[
    is.finite(season_temp_C) &
      climate_coverage >= min_climate_coverage
  ]

  message("Calculating population climate context...")
  context <- annual[
    ,
    {
      yy <- YEAR
      tt <- season_temp_C
      trend <- if (sum(is.finite(tt)) >= 6L &&
                   length(unique(yy[is.finite(tt)])) >= 2L) {
        unname(stats::coef(stats::lm(tt ~ yy))[["yy"]]) * 10
      } else {
        NA_real_
      }

      pred <- if (sum(n_resid_days) > 1L) {
        -sum(resid2_sum, na.rm = TRUE) /
          sum(n_resid_days, na.rm = TRUE)
      } else {
        NA_real_
      }

      list(
        n_climate_years = .N,
        season_background_C = mean(tt),
        season_trend_C_decade = trend,
        season_predictability = pred
      )
    },
    by = .(
      SPECIES, SITE_ID, bms_id, voltinism,
      typical_onset, typical_offset,
      start_day, end_day, expected_days,
      season_photoperiod_h,
      latitude, x_km, y_km,
      n_valid_pheno_years
    )
  ]

  need(all(is.finite(context$season_background_C)),
       "Non-finite background temperature.")
  need(all(is.finite(context$season_photoperiod_h)),
       "Non-finite photoperiod.")

  complete_context <- is.finite(context$season_background_C) &
    is.finite(context$season_trend_C_decade) &
    is.finite(context$season_predictability) &
    is.finite(context$season_photoperiod_h)

  context_excluded <- context[!complete_context]
  write_csv(
    as.data.frame(context_excluded),
    file.path(out_dir, "populations_missing_climate_context.csv")
  )
  context <- context[complete_context]

  scaling <- data.frame(
    variable = c(
      "season_photoperiod_h",
      "season_background_C",
      "season_predictability",
      "season_trend_C_decade"
    ),
    mean = c(
      mean(context$season_photoperiod_h),
      mean(context$season_background_C),
      mean(context$season_predictability),
      mean(context$season_trend_C_decade)
    ),
    SD = c(
      stats::sd(context$season_photoperiod_h),
      stats::sd(context$season_background_C),
      stats::sd(context$season_predictability),
      stats::sd(context$season_trend_C_decade)
    )
  )
  need(all(is.finite(scaling$SD) & scaling$SD > 0),
       "Invalid moderator scaling.")
  write_csv(scaling, file.path(out_dir, "context_scaling.csv"))

  context[, season_photoperiod_z := zscale(season_photoperiod_h)]
  context[, season_background_z := zscale(season_background_C)]
  context[, season_predictability_z := zscale(season_predictability)]
  context[, season_trend_z := zscale(season_trend_C_decade)]

  annual_model <- merge(
    annual,
    context[, .(
      SPECIES, SITE_ID, bms_id, voltinism,
      season_background_C,
      season_photoperiod_z,
      season_background_z,
      season_predictability_z,
      season_trend_z,
      x_km, y_km
    )],
    by = c("SPECIES", "SITE_ID", "bms_id", "voltinism"),
    all = FALSE,
    sort = FALSE
  )
  annual_model[
    ,
    season_temp_anomaly_C :=
      season_temp_C - season_background_C
  ]

  response <- phq[, .(
    SPECIES, SITE_ID, bms_id, voltinism,
    YEAR,
    FLIGHT_LENGTH_mean,
    ONSET_mean,
    OFFSET_mean,
    LAST_PEAK
  )]

  model_data <- merge(
    response,
    annual_model[, .(
      SPECIES, SITE_ID, bms_id, voltinism,
      YEAR,
      season_temp_C,
      season_temp_anomaly_C,
      season_photoperiod_z,
      season_background_z,
      season_predictability_z,
      season_trend_z,
      x_km, y_km
    )],
    by = c("SPECIES", "SITE_ID", "bms_id", "voltinism", "YEAR"),
    all = FALSE,
    sort = FALSE
  )

  model_data[, source_row := .I]
  model_data[, SITE_ID := factor(SITE_ID)]
  model_data[, SPECIES := factor(SPECIES)]
  model_data[, SPECIES_slope := SPECIES]
  model_data[, site_year_id := interaction(
    SITE_ID, YEAR, drop = TRUE, lex.order = TRUE
  )]

  numeric_needed <- c(
    "FLIGHT_LENGTH_mean", "season_temp_anomaly_C",
    "season_photoperiod_z", "season_background_z",
    "season_predictability_z", "season_trend_z",
    "x_km", "y_km"
  )
  need(all(vapply(
    model_data[, ..numeric_needed],
    function(x) all(is.finite(x)),
    logical(1)
  )), "Non-finite values remain in final pilot data.")

  need(!anyDuplicated(
    model_data[, .(SPECIES, SITE_ID, YEAR)]
  ), "Duplicate population-year rows in final pilot data.")

  write_csv(
    as.data.frame(model_data),
    file.path(out_dir, "flight_length_model_data.csv")
  )
  write_csv(
    as.data.frame(context),
    file.path(out_dir, "flight_length_population_context.csv")
  )

  message("Fitting non-spatial pilot models...")
  formula <- stats::as.formula(
    paste(
      "FLIGHT_LENGTH_mean ~ season_temp_anomaly_C *",
      "(season_photoperiod_z + season_background_z +",
      "season_predictability_z + season_trend_z) +",
      "(1 | SITE_ID) + (1 | SPECIES) +",
      "(0 + season_temp_anomaly_C | SPECIES_slope) +",
      "(1 | site_year_id)"
    )
  )

  fit_one <- function(group) {
    d <- as.data.frame(model_data[voltinism == group])
    need(nrow(d) > 0L, paste("No model rows for", group))

    m <- lme4::lmer(
      formula,
      data = d,
      REML = FALSE,
      control = lme4::lmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = 2e6)
      )
    )

    beta <- lme4::fixef(m)
    V <- as.matrix(stats::vcov(m))
    se <- sqrt(diag(V))

    fixed <- data.frame(
      group = group,
      term = names(beta),
      estimate = as.numeric(beta),
      std.error = as.numeric(se),
      conf.low = as.numeric(beta - 1.96 * se),
      conf.high = as.numeric(beta + 1.96 * se),
      p_Wald = as.numeric(
        2 * stats::pnorm(-abs(beta / se))
      ),
      stringsAsFactors = FALSE
    )

    pop <- unique(d[, c(
      "SPECIES", "SITE_ID", "bms_id",
      "season_photoperiod_z",
      "season_background_z",
      "season_predictability_z",
      "season_trend_z"
    )])
    pop$SPECIES <- as.character(pop$SPECIES)
    pop$SITE_ID <- as.character(pop$SITE_ID)
    pop$pop_id <- paste(pop$SPECIES, pop$SITE_ID, sep = "_")

    re <- lme4::ranef(m)$SPECIES_slope
    need("season_temp_anomaly_C" %in% names(re),
         "Species anomaly random slope not found.")

    pop$species_slope_dev <- re[
      match(pop$SPECIES, rownames(re)),
      "season_temp_anomaly_C"
    ]
    need(all(is.finite(pop$species_slope_dev)),
         "Species random-slope join failed.")

    anomaly <- "season_temp_anomaly_C"
    mods <- c(
      "season_photoperiod_z",
      "season_background_z",
      "season_predictability_z",
      "season_trend_z"
    )

    pop$fixed_anomaly_slope <- unname(beta[[anomaly]])
    pop$interaction_contribution <- 0

    for (mod in mods) {
      term <- intersect(
        c(
          paste(anomaly, mod, sep = ":"),
          paste(mod, anomaly, sep = ":")
        ),
        names(beta)
      )
      need(length(term) == 1L,
           paste("Missing interaction:", mod))

      contribution <- beta[[term]] * pop[[mod]]
      pop[[paste0("contribution__", mod)]] <- contribution
      pop$interaction_contribution <-
        pop$interaction_contribution + contribution
    }

    pop$flight_length_plasticity_days_per_C <-
      pop$fixed_anomaly_slope +
      pop$species_slope_dev +
      pop$interaction_contribution

    need(all(is.finite(pop$flight_length_plasticity_days_per_C)),
         "Non-finite contextual flight-length plasticity.")

    pop$flight_length_plasticity_z <- as.numeric(scale(
      pop$flight_length_plasticity_days_per_C
    ))
    pop$voltinism <- group

    convergence <- data.frame(
      group = group,
      n_observations = stats::nobs(m),
      n_populations = length(unique(
        paste(d$SPECIES, d$SITE_ID, sep = "_")
      )),
      n_species = length(unique(d$SPECIES)),
      n_sites = length(unique(d$SITE_ID)),
      singular = lme4::isSingular(m, tol = 1e-4),
      opt_code = if (!is.null(m@optinfo$conv$opt)) {
        as.integer(m@optinfo$conv$opt)
      } else {
        NA_integer_
      },
      opt_message = paste(
        unlist(m@optinfo$conv$lme4$messages),
        collapse = " | "
      ),
      logLik = as.numeric(stats::logLik(m)),
      AIC = stats::AIC(m),
      stringsAsFactors = FALSE
    )

    list(
      model = m,
      fixed = fixed,
      population = pop,
      convergence = convergence
    )
  }

  fits <- lapply(c("univoltine", "multivoltine"), fit_one)
  names(fits) <- c("univoltine", "multivoltine")

  fixed_all <- dplyr::bind_rows(lapply(fits, `[[`, "fixed"))
  pop_all <- dplyr::bind_rows(lapply(fits, `[[`, "population"))
  convergence_all <- dplyr::bind_rows(
    lapply(fits, `[[`, "convergence")
  )

  write_csv(
    fixed_all,
    file.path(out_dir, "pilot_fixed_effects.csv")
  )
  write_csv(
    pop_all,
    file.path(out_dir, "pilot_population_flight_length_plasticity.csv")
  )
  write_csv(
    convergence_all,
    file.path(out_dir, "pilot_fit_checks.csv")
  )

  saveRDS(
    lapply(fits, `[[`, "model"),
    file.path(out_dir, "pilot_models_nonspatial.rds")
  )

  axis_file <- file.path(
    root, "output", "population_trends",
    "spatial_abundance_qc_offset_final_15km",
    "run_eee755f6aac4",
    "global_onset_offset_balance",
    "population_global_balance_axes.csv"
  )

  comparison <- data.frame()
  if (file.exists(axis_file)) {
    ax <- utils::read.csv(
      axis_file,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    need(all(c(
      "group", "pop_id", "season_extension",
      "coordinated_advancement"
    ) %in% names(ax)),
    "Global-balance axis file has unexpected columns.")

    ax <- unique(ax[, c(
      "group", "pop_id",
      "season_extension",
      "coordinated_advancement"
    )])

    tmp <- merge(
      pop_all,
      ax,
      by.x = c("voltinism", "pop_id"),
      by.y = c("group", "pop_id"),
      all = FALSE
    )

    comparison <- tmp |>
      dplyr::group_by(voltinism) |>
      dplyr::summarise(
        n_populations = dplyr::n(),
        cor_direct_flightlength_vs_relative_extension =
          stats::cor(
            flight_length_plasticity_z,
            season_extension
          ),
        cor_direct_flightlength_vs_coordinated_advancement =
          stats::cor(
            flight_length_plasticity_z,
            coordinated_advancement
          ),
        .groups = "drop"
      )

    write_csv(
      tmp,
      file.path(out_dir, "pilot_vs_previous_balance_population_data.csv")
    )
    write_csv(
      comparison,
      file.path(out_dir, "pilot_vs_previous_balance_correlations.csv")
    )
  }

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
  focal <- fixed_all[fixed_all$term %in% focal_terms, ]
  write_csv(
    focal,
    file.path(out_dir, "pilot_focal_temperature_effects.csv")
  )

  writeLines(c(
    "phenoIMPACT | Direct flight-length plasticity PILOT",
    "",
    "This is intentionally non-spatial. It validates the direct metric and",
    "new fixed-season thermal exposure before expensive 15-km spatial refits.",
    "",
    "Positive plasticity = longer flight season in warmer-than-usual years.",
    "",
    "Thermal exposure is mean annual temperature across a FIXED population-",
    "specific typical flight season (mean QC onset to mean QC offset).",
    "The window does not depend on that year's observed flight length.",
    "",
    "Known OFFSET_mean < LAST_PEAK failures are excluded; no other new",
    "phenology outlier filter is imposed.",
    "",
    "If this pilot is coherent and not singular/pathological, the next step is",
    "to fit the same two models with the reviewed annual IID Matern field and",
    "15-km mesh, then use the spatial contextual plasticities in abundance trends.",
    "",
    capture.output(print(convergence_all, row.names = FALSE)),
    "",
    "FOCAL FIXED EFFECTS:",
    capture.output(print(focal, row.names = FALSE)),
    "",
    "COMPARISON WITH PREVIOUS ONSET-OFFSET BALANCE:",
    if (nrow(comparison)) {
      capture.output(print(comparison, row.names = FALSE))
    } else {
      "Previous global-balance axis file was not found; comparison skipped."
    }
  ), file.path(out_dir, "README_results.txt"))

  cat("\n==============================================================\n")
  cat("DIRECT FLIGHT-LENGTH PLASTICITY PILOT\n")
  cat("NON-SPATIAL — NO ABUNDANCE MODEL FITTED\n")
  cat("==============================================================\n\n")

  cat("Phenology QC:\n")
  print(qc_summary, row.names = FALSE)

  cat("\nPilot fit checks:\n")
  print(convergence_all, row.names = FALSE)

  cat("\nFocal temperature effects:\n")
  print(focal, row.names = FALSE)

  if (nrow(comparison)) {
    cat("\nCorrelation with previous onset-offset balance axes:\n")
    print(comparison, row.names = FALSE)
  }

  cat("\nOutput: ",
      normalizePath(out_dir, winslash = "/", mustWork = TRUE),
      "\n", sep = "")

  invisible(list(
    qc = qc_summary,
    fit_checks = convergence_all,
    focal = focal,
    population_plasticity = pop_all,
    comparison = comparison,
    output = out_dir
  ))
}

flight_length_pilot <- run_flight_length_plasticity_pilot()
