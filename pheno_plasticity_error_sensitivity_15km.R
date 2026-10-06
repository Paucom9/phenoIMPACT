# phenoIMPACT: Gaussian / species dispersion / Student-t sensitivity.
# Source this file, then:
#   pheno_error_sensitivity$setup()
#   pheno_error_sensitivity$run()
# Default: offset_multivoltine pilot, three sequential fits in fresh R workers.
# Existing data, models and the ordinary R library are never overwritten.

pheno_error_sensitivity <- local({
  CODE_VERSION <- "error_sensitivity_1.0.1"
  PIN_SHA <- "929257861c690137a9f8661058097ec3c163993c"
  PIN_VERSION <- "1.1.0.9015"
  MIN_DISP_N <- 100L
  PIN_WIN45 <- "https://r2.ropensci.org/5eb9dd8254a375276e1a99e7cb5f42ab658da7390ad0a161727e6a25d8a9c285"
  DEFAULT_ROOT <- "E:/phenoIMPACT project/code/phenoIMPACT"
  VARIANTS <- c("gaussian_constant", "gaussian_species", "student_species")
  state <- new.env(parent = emptyenv())

  assert <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)
  txt <- function(x) paste(deparse(x, width.cutoff = 500L), collapse = " ")
  read_csv <- function(path) utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  write_csv <- function(x, path) utils::write.csv(x, path, row.names = FALSE, na = "")
  pid_alive <- function(pid) tryCatch(ps::ps_is_running(ps::ps_handle(pid)),error=function(e)FALSE)
  atomic_rds <- function(x, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    tmp <- tempfile(tmpdir = dirname(path), fileext = ".rds")
    on.exit(unlink(tmp), add = TRUE)
    saveRDS(x, tmp, compress = "gzip")
    backup <- paste0(path, ".previous")
    if (file.exists(path)) {
      assert(file.copy(path, backup, overwrite = TRUE), paste("Cannot back up", path))
      assert(file.remove(path), paste("Cannot replace", path))
    }
    if (!file.rename(tmp, path)) {
      if (file.exists(backup)) file.copy(backup, path, overwrite = TRUE)
      stop("Cannot save ", path)
    }
    if (file.exists(backup)) unlink(backup)
    invisible(path)
  }
  paths <- function(root = DEFAULT_ROOT) {
    root <- normalizePath(root, winslash = "/", mustWork = TRUE)
    list(root = root,
      source = file.path(root, "output", "phenology_plasticity", "spatial_models"),
      output = file.path(root, "output", "phenology_plasticity", "error_sensitivity_15km"),
      lib = file.path(root, "r_libraries", paste0("error_sensitivity_", substr(PIN_SHA, 1, 12))))
  }
  engine_file <- function(dir) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    env <- environment(engine_file)
    n <- ls(env, all.names = TRUE)
    n <- n[vapply(n, function(k) is.function(get(k, env)), logical(1))]
    n <- c("CODE_VERSION", "PIN_SHA", "PIN_VERSION", "MIN_DISP_N", "PIN_WIN45", "DEFAULT_ROOT", "VARIANTS", n)
    path <- file.path(dir, "workflow_engine.R")
    dump(n, file = path, envir = env, control = "all")
    normalizePath(path, winslash = "/", mustWork = TRUE)
  }
  package_record <- function() {
    pk <- c("sdmTMB", "TMB", "RTMB", "Matrix", "fmesher", "reformulas", "lme4")
    data.frame(package = pk, version = vapply(pk, function(x)
      as.character(utils::packageVersion(x)), character(1)))
  }
  verify_pin <- function(lib) {
    .libPaths(c(lib, .libPaths()))
    assert(requireNamespace("sdmTMB", quietly = TRUE), "sdmTMB is missing from the sensitivity library; run setup().")
    d <- utils::packageDescription("sdmTMB")
    assert(identical(as.character(d$Version), PIN_VERSION) && identical(as.character(d$RemoteSha), PIN_SHA),
      "Unexpected sdmTMB build. This workflow requires its pinned commit; run setup().")
    assert("dispformula" %in% names(formals(sdmTMB::sdmTMB)), "This build does not support dispformula.")
    assert(normalizePath(find.package("sdmTMB"), winslash = "/") ==
      normalizePath(file.path(lib, "sdmTMB"), winslash = "/"), "sdmTMB was loaded from the ordinary library.")
    options(sdmTMB.cores = 1L, sdmTMB.backend = "tmb")
    invisible(TRUE)
  }
  missing_local_packages <- function(required, lib) {
    installed <- utils::installed.packages(lib.loc = lib, noCache = TRUE)
    setdiff(required, installed[, "Package"])
  }
  install_pinned_binary <- function(lib, download = utils::download.file,
                                    install = utils::install.packages) {
    stage <- tempfile("sdmTMB_binary_")
    assert(dir.create(stage), "Cannot create the temporary binary download directory.")
    on.exit(unlink(stage, recursive = TRUE), add = TRUE)
    # Windows infers the package directory from the ZIP basename, even with repos=NULL.
    # Keep the package name; a random tempfile basename leads to a missing DESCRIPTION.
    dest <- file.path(stage, paste0("sdmTMB_", PIN_VERSION, ".zip"))
    download(PIN_WIN45, dest, mode = "wb", quiet = FALSE)
    assert("sdmTMB/DESCRIPTION" %in% utils::unzip(dest, list = TRUE)$Name,
      "The downloaded ZIP does not contain sdmTMB/DESCRIPTION.")
    install(dest, repos = NULL, type = "win.binary", lib = lib)
    invisible(NULL)
  }
  setup_worker <- function(lib) {
    dir.create(lib, recursive = TRUE, showWarnings = FALSE)
    .libPaths(c(lib, .libPaths()))
    options(timeout = max(600, getOption("timeout")), Ncpus = 1L)
    existing <- file.path(lib, "sdmTMB", "DESCRIPTION")
    good <- FALSE
    if (file.exists(existing)) {
      dd <- read.dcf(existing)
      good <- "RemoteSha" %in% colnames(dd) && identical(unname(dd[1, "RemoteSha"]), PIN_SHA)
    }
    if (!good) {
      # All explicitly installed packages go into the isolated library.
      # Reuse successful installations when retrying an interrupted setup.
      req <- c("RTMB", "TMB", "Matrix", "RcppEigen", "Rcpp", "abind", "fishMod", "igraph",
        "assertthat", "cli", "fmesher", "generics", "lifecycle", "mvtnorm",
        "reformulas", "rlang", "lme4", "digest", "callr", "processx", "zip")
      missing <- missing_local_packages(req, lib)
      if (length(missing)) {
        utils::install.packages(missing, lib = lib, repos = "https://cloud.r-project.org",
          dependencies = NA)
      } else {
        message("Reusing dependencies already installed in the sensitivity library.")
      }
      if (.Platform$OS.type == "windows" && startsWith(as.character(getRversion()), "4.5.")) {
        install_pinned_binary(lib)
      } else {
        if (!requireNamespace("remotes", quietly = TRUE))
          utils::install.packages("remotes", lib = lib, repos = "https://cloud.r-project.org")
        remotes::install_github(paste0("sdmTMB/sdmTMB@", PIN_SHA), lib = lib,
          dependencies = NA, upgrade = "never", build_vignettes = FALSE)
      }
    }
    verify_pin(lib)
    versions <- package_record()
    write_csv(versions, file.path(lib, "sensitivity_package_versions.csv"))
    versions
  }
  setup <- function(root = DEFAULT_ROOT) {
    assert(requireNamespace("callr", quietly = TRUE), "Install callr in your ordinary R library first.")
    p <- paths(root)
    e <- engine_file(file.path(p$output, "setup"))
    message("Preparing separate sensitivity library: ", p$lib)
    ans <- callr::r(function(engine, lib) {
      x <- new.env(parent = globalenv()); sys.source(engine, x); x$setup_worker(lib)
    }, args = list(e, p$lib), show = TRUE, libpath = .libPaths(),user_profile=FALSE,system_profile=FALSE)
    message("Pinned sensitivity environment ready. The ordinary R library was not changed.")
    invisible(ans)
  }

  locate <- function(source_root, analysis) {
    dirs <- list.dirs(source_root, recursive = FALSE, full.names = TRUE)
    dirs <- dirs[startsWith(basename(dirs), paste0("primary__", analysis, "_"))]
    need <- c("metadata.rds", "prepared_input.rds", "fit_spatial_field.rds")
    dirs <- dirs[vapply(dirs, function(d) all(file.exists(file.path(d, need))), logical(1))]
    assert(length(dirs) == 1L, paste("Expected exactly one primary input for", analysis, "; found", length(dirs)))
    cp <- if (analysis == "offset_multivoltine")
      file.path(source_root, "offset_recovery_15km", basename(dirs), "fit_spatial_field.rds") else
      file.path(dirs, "fit_spatial_field.rds")
    assert(file.exists(cp), paste("Required checkpoint missing:", cp))
    list(source = normalizePath(dirs, winslash = "/"), checkpoint = normalizePath(cp, winslash = "/"))
  }
  dispersion_groups <- function(data, minimum = MIN_DISP_N) {
    # Keep EVERY observation. Low-n species share one scale parameter only;
    # their species identities in the mean/random effects remain unchanged.
    n <- table(data$SPECIES)
    common <- names(n)[n >= minimum]
    lab <- ifelse(as.character(data$SPECIES) %in% common,
      paste0("species:",as.character(data$SPECIES)),"pooled_low_n")
    data$disp_group <- factor(lab)
    assert(nlevels(data$disp_group)>1L,"Not enough supported dispersion groups for this sensitivity.")
    data
  }
  extract_input <- function(job, out) {
    # Fresh process in the ORIGINAL library. parList is the R-level parameter
    # extractor. Never call a serialized native fn/gr/report pointer or optimizer.
    requireNamespace("sdmTMB", quietly = TRUE)
    p <- readRDS(file.path(job$source, "prepared_input.rds"))
    meta <- readRDS(file.path(job$source, "metadata.rds"))
    s <- readRDS(job$checkpoint); f <- s$fit
    assert(inherits(f, "sdmTMB") && identical(s$run_id, meta$run_id) &&
      identical(s$label, "spatial_field"), "Original checkpoint identity mismatch.")
    assert(isTRUE(meta$settings$mesh_cutoff_km == 15), "Expected a 15-km primary fit.")
    assert(identical(txt(p$formula), txt(meta$formula)), "Original formula mismatch.")
    assert(isTRUE(all.equal(f$data[, names(p$data), drop = FALSE], p$data)), "Original data differ from prepared input.")
    assert(identical(f$family$family, "gaussian") && identical(f$family$link, "identity"), "Expected Gaussian identity baseline.")
    assert(!isTRUE(f$reml) && !any(f$offset != 0), "Unexpected REML or offset specification.")
    assert(isTRUE(f$model$convergence == 0) && isTRUE(f$sd_report$pdHess) &&
      length(f$gradients) > 0 && all(is.finite(f$gradients)) && max(abs(f$gradients)) < .001,
      "Original checkpoint fails numerical checks.")
    assert(!anyNA(p$data) && !anyDuplicated(p$data$source_row), "Missing data or duplicate source_row.")
    X <- stats::model.matrix(lme4::nobars(p$formula), p$data)
    # run_extra_optimization() can leave $parlist stale. Extract the saved best
    # parameter vector instead, and verify it against the saved covariance.
    old_parameters <- f$tmb_obj$env$parList(par = f$tmb_obj$env$last.par.best)
    b <- as.numeric(old_parameters$b_j)
    idx <- which(names(f$sd_report$par.fixed) == "b_j")
    assert(length(b) == ncol(X) && length(idx) == length(b), "Cannot align saved fixed coefficients.")
    assert(max(abs(b - f$sd_report$par.fixed[idx])) < 1e-6,
      "Saved point estimates and covariance refer to different parameter values.")
    V <- f$sd_report$cov.fixed[idx, idx, drop = FALSE]
    se <- sqrt(diag(V)); names(b) <- colnames(X)
    assert(all(is.finite(b)) && all(is.finite(se) & se > 0), "Invalid baseline coefficient SE.")
    fixed <- data.frame(term = names(b), estimate = unname(b), std.error = se,
      conf.low = b - qnorm(.975) * se, conf.high = b + qnorm(.975) * se, row.names = NULL)
    signature <- digest::digest(list(data = p$data, formula = txt(p$formula), mesh = f$spde,
      checkpoint = unname(tools::md5sum(job$checkpoint))), algo = "sha256")
    d <- dispersion_groups(p$data)
    mapping <- unique(data.frame(SPECIES=as.character(d$SPECIES),disp_group=as.character(d$disp_group)))
    mapping$n <- as.integer(table(d$SPECIES)[mapping$SPECIES])
    mapping$minimum_for_separate_scale <- MIN_DISP_N
    write_csv(mapping,file.path(out,"dispersion_group_mapping.csv"))
    neutral <- list(data = d, formula = p$formula, response = p$response,
      mesh = f$spde, old_start = old_parameters, old_fixed = fixed,
      old_objective = f$model$objective, old_run_id = meta$run_id,
      old_versions = meta$package_versions, signature = signature,
      checkpoint = job$checkpoint, source = job$source)
    atomic_rds(neutral, file.path(out, "input.rds"))
    write_csv(fixed, file.path(out, "original_fixed_effects.csv"))
    write_csv(meta$package_versions, file.path(out, "original_package_versions.csv"))
    write_csv(data.frame(n = nrow(p$data), species = nlevels(p$data$SPECIES),
      source_signature = signature, checkpoint = job$checkpoint, formula = txt(p$formula)),
      file.path(out, "input_manifest.csv"))
    signature
  }
  fit_checks <- function(fit) {
    fixed <- as.data.frame(sdmTMB::tidy(fit, effects = "fixed"))
    g <- fit$gradients
    grad <- if (length(g) && all(is.finite(g))) max(abs(g)) else NA_real_
    code <- fit$model$convergence
    if (length(code) != 1L) code <- NA_integer_
    pd <- isTRUE(fit$sd_report$pdHess)
    se_ok <- all(is.finite(fixed$estimate)) && all(is.finite(fixed$std.error) & fixed$std.error > 0)
    list(table = data.frame(convergence_code = code, positive_definite_Hessian = pd,
      max_abs_gradient = grad, finite_fixed_SE = se_ok,
      fit_ok = isTRUE(code == 0) && pd && is.finite(grad) && grad < .001 && se_ok), fixed = fixed)
  }
  fixed_comparison <- function(base, other) {
    cols <- c("term", "estimate", "std.error", "conf.low", "conf.high")
    assert(all(cols %in% names(base)) && all(cols %in% names(other)), "Missing fixed-effect columns.")
    assert(!anyDuplicated(base$term) && !anyDuplicated(other$term) && setequal(base$term, other$term),
      "Fixed-effect terms differ between comparable models.")
    z <- merge(base[, cols], other[, cols], by = "term", suffixes = c("_reference", "_alternative"), sort = FALSE)
    z$delta_estimate <- z$estimate_alternative - z$estimate_reference
    z$delta_in_reference_SE <- z$delta_estimate / z$std.error_reference
    z$SE_ratio <- z$std.error_alternative / z$std.error_reference
    z$CI_width_ratio <- (z$conf.high_alternative - z$conf.low_alternative) /
      (z$conf.high_reference - z$conf.low_reference)
    z$sign_changed <- sign(z$estimate_reference) != sign(z$estimate_alternative)
    z$reference_excludes_zero <- z$conf.low_reference > 0 | z$conf.high_reference < 0
    z$alternative_excludes_zero <- z$conf.low_alternative > 0 | z$conf.high_alternative < 0
    z$zero_exclusion_changed <- z$reference_excludes_zero != z$alternative_excludes_zero
    z$is_anomaly_term <- grepl("clim_anomaly", z$term, fixed = TRUE)
    z$is_interaction <- grepl(":", z$term, fixed = TRUE)
    z
  }
  make_start <- function(previous, variant, data) {
    # Only shared, explicitly supported parameters; no previous_fit across
    # families, versions or dispersion dimensions.
    keys <- c("b_j", "ln_kappa", "ln_tau_E", "re_cov_pars", "re_b_pars", "epsilon_st")
    assert(all(keys %in% names(previous)), "Starting parameter record is incomplete.")
    start <- previous[keys]
    if (variant == "gaussian_constant") start$ln_phi <- previous$ln_phi else {
      k <- nlevels(data$disp_group)
      if (!is.null(previous$b_disp_k) && length(previous$b_disp_k) == k) {
        start$b_disp_k <- previous$b_disp_k
      } else {
        assert(length(previous$ln_phi) == 1, "Cannot initialize dispersion.")
        start$b_disp_k <- rep(as.numeric(previous$ln_phi), k)
      }
    }
    if (variant == "student_species") {
      # df starts at 8, but is estimated (>1), not fixed. Scale is NOT Student-t SD.
      start$ln_student_df <- log(8 - 1)
      start$b_disp_k <- start$b_disp_k + log(sqrt((8 - 2) / 8))
    }
    start
  }
  normal_quantile <- function(y, mu, scale, family, df = NA_real_) {
    assert(length(y) == length(mu) && length(scale) %in% c(1L, length(y)), "Residual dimensions differ.")
    assert(all(is.finite(scale) & scale > 0), "Invalid observation scale.")
    u <- (y - mu) / scale
    if (family == "gaussian") return(u)
    assert(family == "student" && length(df) == 1 && is.finite(df) && df > 1, "Invalid Student-t df.")
    # Tail-stable probability integral transform: avoids qnorm(1) rounding.
    z <- numeric(length(u)); neg <- u <= 0
    z[neg] <- qnorm(pt(u[neg], df = df, log.p = TRUE), log.p = TRUE)
    z[!neg] <- qnorm(pt(u[!neg], df = df, lower.tail = FALSE, log.p = TRUE),
      lower.tail = FALSE, log.p = TRUE)
    z
  }
  residual_draw <- function(fit, seed) {
    # Public simulate API gives one approximate posterior random-effect draw.
    # Read eta/scale from its TMB report, then apply the exact Gaussian/Student CDF.
    # residuals(type='mle-mvn') is deliberately NOT called with dispformula.
    sim <- stats::simulate(fit, nsim = 1L, seed = seed, type = "mle-mvn",
      mle_mvn_samples = "single", observation_error = FALSE,
      return_tmb_report = TRUE, silent = TRUE)[[1L]]
    n <- nrow(fit$data)
    mu <- as.numeric(sim$eta_i)
    scale <- if (isTRUE(fit$has_dispformula)) as.numeric(sim$phi_i) else rep(as.numeric(sim$phi), n)
    df <- if (fit$family$family == "student") as.numeric(sim$student_df) else NA_real_
    assert(length(mu) == n && length(scale) == n, "Pinned simulation report has unexpected dimensions.")
    z <- normal_quantile(as.numeric(fit$response), mu, scale, fit$family$family, df)
    assert(all(is.finite(z)), "Non-finite posterior-draw quantile residuals.")
    list(z = z, mu_draw = mu, scale = scale, df = df)
  }
  residual_summary <- function(z) {
    q <- quantile(z, c(.001, .01, .5, .99, .999), names = FALSE)
    m <- mean(z); s <- sd(z)
    data.frame(n = length(z), mean = m, sd = s,
      skew = mean((z-m)^3) / mean((z-m)^2)^1.5,
      excess_kurtosis = mean((z-m)^4) / mean((z-m)^2)^2 - 3,
      fraction_abs_gt_3 = mean(abs(z) > 3), q001=q[1], q01=q[2], median=q[3], q99=q[4], q999=q[5])
  }
  bin_residuals <- function(x, z, groups = 20L) {
    # Fixed equal-count ranks; ties kept deterministic in input row order.
    ord <- order(x, seq_along(x)); g <- integer(length(x))
    g[ord] <- pmin(groups, ceiling(seq_along(x) * groups / length(x)))
    do.call(rbind, lapply(split(seq_along(x), g), function(i)
      data.frame(x = mean(x[i]), n = length(i), mean = mean(z[i]), sd = sd(z[i]))))
  }
  diagnostics <- function(fit, input, out, seeds) {
    pred <- stats::predict(fit)
    assert(nrow(pred) == nrow(input$data) && identical(pred$source_row, input$data$source_row), "Prediction row order differs.")
    summaries <- list(); qq <- list()
    for (k in seq_along(seeds)) {
      r <- residual_draw(fit, seeds[k]); z <- r$z
      summaries[[k]] <- cbind(seed=seeds[k], student_df=r$df, residual_summary(z))
      prob <- sort(unique(c(seq(.001,.999,length.out=999),.0001,.0005,.9995,.9999)))
      qq[[k]] <- data.frame(seed=seeds[k], theoretical=qnorm(prob), observed=as.numeric(quantile(z,prob)))
      if (k == 1L) {
        keep <- intersect(c("source_row","SITE_ID","SPECIES","YEAR","x_km","y_km"), names(input$data))
        rows <- input$data[,keep,drop=FALSE]
        rows$observed <- input$data[[input$response]]; rows$fitted <- pred$est
        rows$response_residual <- rows$observed - rows$fitted
        rows$observation_scale <- r$scale; rows$mvn_quantile_residual <- z
        con <- gzfile(file.path(out,"residuals.csv.gz"),"wt")
        tryCatch(write_csv(rows,con),finally=close(con))
        write_csv(bin_residuals(pred$est,z),file.path(out,"residuals_by_fitted_bin.csv"))
        sp <- do.call(rbind,lapply(split(seq_along(z),input$data$SPECIES),function(i)
          data.frame(SPECIES=as.character(input$data$SPECIES[i[1]]),n=length(i),mean=mean(z[i]),sd=sd(z[i]),
            scale=mean(r$scale[i]))))
        write_csv(sp,file.path(out,"residuals_by_species.csv"))
        write_csv(data.frame(student_df=r$df,finite_residual_variance=is.na(r$df)||r$df>2),
          file.path(out,"tail_parameters.csv"))
      }
      rm(r); invisible(gc())
    }
    write_csv(do.call(rbind,summaries),file.path(out,"residual_distribution.csv"))
    write_csv(do.call(rbind,qq),file.path(out,"residual_qq.csv"))
    invisible(TRUE)
  }
  smoke_test <- function(lib, out) {
    verify_pin(lib)
    set.seed(7429)
    d <- data.frame(SPECIES=factor(rep(letters[1:3],each=300)),x=rnorm(900),
      SITE_ID=factor(rep(1:30,30)),YEAR=rep(1:3,each=300),source_row=1:900)
    site <- data.frame(x_km=runif(30,0,100),y_km=runif(30,0,100))
    d$x_km <- site$x_km[as.integer(d$SITE_ID)];d$y_km <- site$y_km[as.integer(d$SITE_ID)]
    d <- dispersion_groups(d)
    d$y <- 3 + 1.4*d$x + rnorm(30)[as.integer(d$SITE_ID)] +
      rt(900,df=5)*sqrt(3/5)*c(.5,1.5,3)[as.integer(d$SPECIES)]
    mesh <- sdmTMB::make_mesh(d,xy_cols=c("x_km","y_km"),cutoff=20)
    fs <- list()
    for (v in VARIANTS) {
      fit <- sdmTMB::sdmTMB(y~x+(1|SITE_ID),data=d,mesh=mesh,spatial="off",
        family=if(v=="student_species")sdmTMB::student(df=NULL) else gaussian(),
        dispformula=if(v=="gaussian_constant")~1 else ~0+disp_group,
        control=sdmTMB::sdmTMBcontrol(parallel=1L,backend="tmb"),silent=TRUE)
      cc <- fit_checks(fit)
      assert(cc$table$fit_ok,paste("Small numerical test failed:",v))
      dr <- residual_draw(fit,901)
      if(v=="gaussian_constant") {
        set.seed(901); native <- as.numeric(residuals(fit,type="mle-mvn"))
        assert(max(abs(native-dr$z))<1e-5,"Analytical MVN residuals do not match native Gaussian residuals.")
      }
      if(v=="gaussian_species") {
        ss <- tapply(dr$scale,d$SPECIES,mean)
        assert(ss[3]/ss[1]>3,"Small test did not recover the simulated dispersion ordering.")
      }
      fs[[v]] <- cbind(variant=v,cc$table,student_df=dr$df)
    }
    write_csv(do.call(rbind,fs),file.path(out,"small_model_checks.csv"))
    atomic_rds(list(pin=PIN_SHA,versions=package_record(),passed=TRUE),file.path(out,"preflight_passed.rds"))
    TRUE
  }
  worker <- function(cfg, analysis, variant) {
    verify_pin(cfg$lib)
    adir <- file.path(cfg$run,analysis); out <- file.path(adir,variant)
    dir.create(out,recursive=TRUE,showWarnings=FALSE)
    phase <- function(s) {atomic_rds(list(phase=s,time=Sys.time()),file.path(out,"progress.rds"));message(s)}
    phase("Loading verified input")
    input <- readRDS(file.path(adir,"input.rds"))
    fingerprint <- digest::digest(list(version=CODE_VERSION,input=input$signature,variant=variant,minimum=MIN_DISP_N,
      pin=PIN_SHA,packages=package_record(),seeds=cfg$seeds),algo="sha256")
    checkpoint <- file.path(out,"fit.rds")
    if (file.exists(file.path(out,"completed.rds"))) {
      done <- readRDS(file.path(out,"completed.rds"))
      assert(identical(done$fingerprint,fingerprint),"Completed job signature changed; start a new run.")
      return(done)
    }
    assert(!file.exists(checkpoint),paste("An incomplete checkpoint is retained at",checkpoint,
      "Run a new sensitivity queue; no completed optimization is silently overwritten."))
    prev <- if(variant=="gaussian_constant")input$old_start else {
      predecessor <- if(variant=="gaussian_species")"gaussian_constant" else "gaussian_species"
      readRDS(file.path(adir,predecessor,"parameters.rds"))
    }
    start <- make_start(prev,variant,input$data)
    rm(prev)
    fml <- if(variant=="gaussian_constant")~1 else ~0+disp_group
    fam <- if(variant=="student_species")sdmTMB::student(df=NULL) else gaussian()
    warnings <- character()
    phase(paste("Fitting",variant,"with",nrow(input$data),"rows"))
    t0 <- Sys.time()
    fit <- withCallingHandlers(sdmTMB::sdmTMB(formula=input$formula,data=input$data,mesh=input$mesh,
      time="YEAR",family=fam,spatial="off",spatiotemporal="iid",reml=FALSE,dispformula=fml,
      control=sdmTMB::sdmTMBcontrol(start=start,multiphase=FALSE,parallel=1L,
        newton_loops=0L,backend="tmb",collapse_spatial_variance=FALSE),silent=FALSE),
      warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart("muffleWarning")})
    # Save immediately before any further optimization or diagnostics.
    atomic_rds(list(fit=fit,fingerprint=fingerprint,variant=variant),checkpoint)
    check <- fit_checks(fit)
    if(!check$table$fit_ok) {
      phase("One additional optimization pass; initial checkpoint retained")
      fit <- sdmTMB::run_extra_optimization(fit,nlminb_loops=1L,newton_loops=0L)
      fit$parlist <- fit$tmb_obj$env$parList(par=fit$tmb_obj$env$last.par.best)
      atomic_rds(list(fit=fit,fingerprint=fingerprint,variant=variant),checkpoint)
      check <- fit_checks(fit)
    }
    assert(identical(fit$data$source_row,input$data$source_row),"Fitted rows changed.")
    check$table$elapsed_minutes <- as.numeric(difftime(Sys.time(),t0,units="mins"))
    check$table$logLik <- as.numeric(logLik(fit));check$table$AIC <- AIC(fit)
    check$table$n_parameters <- attr(logLik(fit),"df")
    write_csv(check$table,file.path(out,"fit_checks.csv"))
    writeLines(unique(warnings),file.path(out,"warnings.txt"))
    write_csv(check$fixed,file.path(out,"fixed_effects.csv"))
    writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
    detail <- sdmTMB::sanity(fit,gradient_thresh=.001,silent=TRUE)
    write_csv(data.frame(check=names(detail),passed=as.logical(unlist(detail))),file.path(out,"sanity_checks.csv"))
    write_csv(as.data.frame(sdmTMB::tidy(fit,effects="ran_pars")),file.path(out,"random_parameters.csv"))
    if(variant!="gaussian_constant") write_csv(as.data.frame(sdmTMB::tidy(fit,effects="dispersion")),
      file.path(out,"dispersion_coefficients.csv"))
    atomic_rds(fit$parlist,file.path(out,"parameters.rds"))
    assert(check$table$fit_ok,"Numerical checks failed; the queue will stop and retain this checkpoint.")
    if(variant=="gaussian_constant") {
      bridge <- fixed_comparison(input$old_fixed,check$fixed)
      write_csv(bridge,file.path(adir,"version_bridge_coefficients.csv"))
      objective_change <- fit$model$objective-input$old_objective
      # Explicit computational comparison gate, not a hypothesis test.
      gate <- data.frame(max_abs_change_in_old_SE=max(abs(bridge$delta_in_reference_SE)),
        max_relative_SE_change=max(abs(bridge$SE_ratio-1)), objective_change=objective_change)
      gate$passed <- gate$max_abs_change_in_old_SE<.05 && gate$max_relative_SE_change<.05 && abs(objective_change)<.1
      write_csv(gate,file.path(adir,"version_bridge_checks.csv"))
      assert(gate$passed,"The new-version Gaussian baseline differs from the original. Queue stopped for review.")
    }
    phase("Exporting posterior-draw quantile residuals")
    diagnostic_error <- NULL
    tryCatch(diagnostics(fit,input,out,cfg$seeds),error=function(e)diagnostic_error<<-conditionMessage(e))
    if(!is.null(diagnostic_error))writeLines(diagnostic_error,file.path(out,"diagnostic_error.txt"))
    done <- list(fingerprint=fingerprint,status=if(is.null(diagnostic_error))"OK" else "FIT_OK_DIAGNOSTICS_FAILED",
      variant=variant,time=Sys.time())
    atomic_rds(done,file.path(out,"completed.rds"));phase(done$status)
    done
  }
  collect <- function(run = state$run) {
    assert(length(run)==1 && dir.exists(run),"Supply an existing run directory.")
    cfg <- readRDS(file.path(run,"run_config.rds")); all <- list(); fitrows <- list(); residualrows <- list()
    for(a in cfg$analyses) {
      adir <- file.path(run,a)
      refs <- list(original=file.path(adir,"original_fixed_effects.csv"),
        gaussian_constant=file.path(adir,"gaussian_constant","fixed_effects.csv"),
        gaussian_species=file.path(adir,"gaussian_species","fixed_effects.csv"))
      for(v in VARIANTS) {
        od <- file.path(adir,v); chk <- file.path(od,"fit_checks.csv")
        if(file.exists(chk))fitrows[[paste(a,v)]] <- cbind(analysis=a,variant=v,read_csv(chk))
        if(file.exists(file.path(od,"residual_distribution.csv")))
          residualrows[[paste(a,v)]] <- cbind(analysis=a,variant=v,read_csv(file.path(od,"residual_distribution.csv")))
        if(!file.exists(chk)||!isTRUE(read_csv(chk)$fit_ok))next
        target <- file.path(od,"fixed_effects.csv")
        for(ref in names(refs)) {
          if(ref==v || !file.exists(refs[[ref]]))next
          if(ref=="gaussian_species" && v!="student_species")next
          z <- fixed_comparison(read_csv(refs[[ref]]),read_csv(target))
          all[[paste(a,v,ref)]] <- cbind(analysis=a,reference=ref,alternative=v,z)
        }
      }
    }
    if(length(all)) {
      tab <- do.call(rbind,all);rownames(tab)<-NULL
      write_csv(tab,file.path(run,"coefficient_comparison.csv"))
      write_csv(tab[tab$is_anomaly_term,],file.path(run,"plasticity_comparison.csv"))
      # Absolute changes and SE ratios are the main outputs. No binary 'robust' flag.
    }
    if(length(fitrows))write_csv(do.call(rbind,fitrows),file.path(run,"model_comparison.csv"))
    if(length(residualrows))write_csv(do.call(rbind,residualrows),file.path(run,"residual_comparison.csv"))
    make_figures(run)
    invisible(run)
  }
  make_figures <- function(run) {
    path <- file.path(run,"plasticity_comparison.csv")
    if(!file.exists(path))return(invisible(NULL))
    z <- read_csv(path);z <- z[z$reference=="gaussian_constant",]
    for(a in unique(z$analysis)) {
      q <- z[z$analysis==a,]; terms <- unique(q$term)
      grDevices::png(file.path(run,paste0(a,"_coefficient_comparison.png")),width=1800,height=1100,res=160)
      tryCatch({
        par(mar=c(4,18,1,1)); lim <- range(c(q$conf.low_reference,q$conf.high_reference,q$conf.low_alternative,q$conf.high_alternative))
        plot(NA,xlim=lim,ylim=c(.5,length(terms)+.5),yaxt="n",ylab="",xlab="Coefficient and 95% Wald interval")
        axis(2,at=seq_along(terms),labels=terms,las=1,cex.axis=.7);abline(v=0,lty=2,col="grey60")
        for(k in seq_along(VARIANTS)) {
          v <- VARIANTS[k]; col <- c("#777777","#237e76","#b66b39")[k]
          if(k==1) {
            w <- q[!duplicated(q$term),]; est<-w$estimate_reference;lo<-w$conf.low_reference;hi<-w$conf.high_reference
          } else {w<-q[q$alternative==v,];est<-w$estimate_alternative;lo<-w$conf.low_alternative;hi<-w$conf.high_alternative}
          yy <- match(w$term,terms)+(k-2)*.16
          segments(lo,yy,hi,yy,col=col);points(est,yy,pch=16,col=col,cex=.8)
        }
        legend("topright",legend=VARIANTS,col=c("#777777","#237e76","#b66b39"),pch=16,bty="n",cex=.7)
      },finally=dev.off())
    }
    invisible(NULL)
  }
  supervisor <- function(cfg) {
    # Supervisor holds no fitted models. Every fit is a fresh, sequential process.
    lock <- file.path(cfg$output,"queue.lock")
    on.exit(unlink(lock,recursive=TRUE),add=TRUE)
    q <- expand.grid(variant=VARIANTS,analysis=cfg$analyses,stringsAsFactors=FALSE)
    q$status <- "PENDING"; q$message <- ""
    save_status <- function()write_csv(q,file.path(cfg$run,"queue_status.csv"))
    save_status()
    tryCatch({
      phase <- "Small-model compatibility and residual-method checks"
      message(phase)
      callr::r(function(engine,lib,out){x<-new.env();sys.source(engine,x);x$smoke_test(lib,out)},
        args=list(cfg$engine,cfg$lib,cfg$run),libpath=c(cfg$lib,cfg$original_libpaths),show=TRUE,
        user_profile=FALSE,system_profile=FALSE)
      for(a in cfg$analyses) {
        adir <- file.path(cfg$run,a);dir.create(adir,recursive=TRUE,showWarnings=FALSE)
        callr::r(function(engine,job,out){x<-new.env();sys.source(engine,x);x$extract_input(job,out)},
          args=list(cfg$engine,cfg$jobs[[a]],adir),libpath=cfg$original_libpaths,show=TRUE,
          user_profile=FALSE,system_profile=FALSE)
        for(v in VARIANTS) {
          i<-which(q$analysis==a&q$variant==v);q$status[i]<-"RUNNING";save_status()
          od<-file.path(adir,v);dir.create(od,recursive=TRUE,showWarnings=FALSE)
          log<-file.path(od,"worker.log")
          proc<-callr::r_bg(function(engine,cfg,a,v){x<-new.env();sys.source(engine,x);x$worker(cfg,a,v)},
            args=list(cfg$engine,cfg,a,v),libpath=c(cfg$lib,cfg$original_libpaths),
            stdout=log,stderr="2>&1",supervise=TRUE,user_profile=FALSE,system_profile=FALSE,
            env=c(callr::rcmd_safe_env(),OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1"))
          atomic_rds(list(pid=proc$get_pid(),analysis=a,variant=v,log=log),file.path(cfg$run,"active_worker.rds"))
          proc$wait()
          ans<-tryCatch(proc$get_result(),error=function(e){q$status[i]<<-"FAILED";q$message[i]<<-conditionMessage(e);save_status();stop(e)})
          q$status[i]<-ans$status;save_status();collect(cfg$run)
          message(format(Sys.time())," ",a," / ",v,": ",ans$status)
          assert(ans$status=="OK","Fit saved, but residual export failed. Queue stopped; inspect diagnostic_error.txt.")
        }
      }
      writeLines("QUEUE_FINISHED",file.path(cfg$run,"queue_finished.txt"))
    },error=function(e){writeLines(conditionMessage(e),file.path(cfg$run,"queue_error.txt"));message("QUEUE STOPPED: ",conditionMessage(e))})
    collect(cfg$run)
    pack_results(cfg$run)
    invisible(q)
  }
  pack_results <- function(run) {
    cfg <- readRDS(file.path(run,"run_config.rds"))
    if(!requireNamespace("zip",lib.loc=cfg$lib,quietly=TRUE)) {
      writeLines("Automatic ZIP creation unavailable; send the exported CSV/PNG files.",
        file.path(run,"zip_unavailable.txt"));return(invisible(NULL))
    }
    files<-list.files(run,recursive=TRUE,full.names=TRUE)
    keep<-grepl("\\.(csv|csv\\.gz|png|txt|log|R)$",files,ignore.case=TRUE)
    files<-files[keep & !grepl("/results_to_review/",files,fixed=TRUE)]
    # Excludes all large model checkpoints and prepared input RDS.
    relative <- substring(files,nchar(run)+2L)
    zip::zipr(file.path(run,"results_to_review.zip"),files=relative,root=run,include_directories=FALSE)
    invisible(NULL)
  }
  run <- function(root = DEFAULT_ROOT, analyses = "offset_multivoltine", seeds = c(20260922L,20260923L,20260924L)) {
    assert(requireNamespace("callr",quietly=TRUE)&&requireNamespace("processx",quietly=TRUE),"callr and processx are required.")
    assert(length(analyses)>0 && !anyDuplicated(analyses) && all(analyses %in%
      c("onset","first_peak","offset_univoltine","offset_multivoltine")),"Invalid analysis names.")
    assert(length(seeds)>0 && all(is.finite(seeds)),"Invalid residual seeds.")
    p<-paths(root);assert(dir.exists(p$lib),"Run pheno_error_sensitivity$setup() first.")
    dir.create(p$output,recursive=TRUE,showWarnings=FALSE)
    lock<-file.path(p$output,"queue.lock")
    if(dir.exists(lock)) {
      lp<-file.path(lock,"owner.rds")
      assert(file.exists(lp),paste("Queue lock needs review:",lock))
      owner<-readRDS(lp)
      assert(!pid_alive(owner$pid),paste("A sensitivity queue is already running; PID",owner$pid))
      unlink(lock,recursive=TRUE)
    }
    jobs<-setNames(lapply(analyses,function(a)locate(p$source,a)),analyses)
    assert(dir.create(lock),"Cannot acquire the sensitivity queue lock.")
    keep_lock<-FALSE;on.exit(if(!keep_lock)unlink(lock,recursive=TRUE),add=TRUE)
    out<-file.path(p$output,paste0("run_",format(Sys.time(),"%Y%m%d_%H%M%S"),"_",Sys.getpid()))
    dir.create(out,recursive=TRUE,showWarnings=FALSE)
    engine<-engine_file(out)
    cfg<-c(p,list(run=normalizePath(out,winslash="/"),engine=engine,analyses=analyses,
      seeds=as.integer(seeds),jobs=jobs,original_libpaths=.libPaths(),code_version=CODE_VERSION,pin=PIN_SHA))
    atomic_rds(cfg,file.path(out,"run_config.rds"))
    log<-file.path(out,"queue.log")
    proc<-callr::r_bg(function(engine,cfg){x<-new.env();sys.source(engine,x);x$supervisor(cfg)},
      args=list(engine,cfg),libpath=.libPaths(),stdout=log,stderr="2>&1",supervise=FALSE,
      user_profile=FALSE,system_profile=FALSE)
    atomic_rds(list(pid=proc$get_pid(),run=out),file.path(lock,"owner.rds"))
    atomic_rds(list(run=out,pid=proc$get_pid()),file.path(p$output,"latest_run.rds"))
    state$run<-out;state$process<-proc;keep_lock<-TRUE
    message("Queue started. ",length(analyses)*3," sequential fits; one R thread per fit.\nResults: ",out,
      "\nUse pheno_error_sensitivity$status() or $watch().")
    invisible(out)
  }
  current_run <- function(root=DEFAULT_ROOT) {
    if(!is.null(state$run))return(state$run)
    p<-paths(root);f<-file.path(p$output,"latest_run.rds")
    assert(file.exists(f),"No sensitivity run found.")
    readRDS(f)$run
  }
  status <- function(root=DEFAULT_ROOT,run=current_run(root)) {
    message("Results: ",run)
    f<-file.path(run,"queue_status.csv")
    if(file.exists(f))print(read_csv(f),row.names=FALSE)
    er<-file.path(run,"queue_error.txt")
    if(file.exists(er))message(paste(readLines(er,warn=FALSE),collapse="\n"))
    af<-file.path(run,"active_worker.rds")
    if(file.exists(af)) {
      w<-readRDS(af);pf<-file.path(run,w$analysis,w$variant,"progress.rds")
      if(file.exists(pf))print(readRDS(pf))
    }
    invisible(run)
  }
  watch <- function(root=DEFAULT_ROOT,run=current_run(root),seconds=30) {
    repeat {
      status(root,run)
      if(file.exists(file.path(run,"queue_finished.txt"))||file.exists(file.path(run,"queue_error.txt")))break
      Sys.sleep(seconds)
    }
    invisible(run)
  }
  list(setup=setup,run=run,status=status,watch=watch,collect=collect,
    internal=list(fixed_comparison=fixed_comparison,normal_quantile=normal_quantile,
      make_start=make_start,engine_file=engine_file,smoke_test=smoke_test,dispersion_groups=dispersion_groups,
      residual_summary=residual_summary,bin_residuals=bin_residuals,paths=paths))
})
message("Loaded pheno_error_sensitivity. Run $setup() once, then $run() for the multivoltine pilot.")
