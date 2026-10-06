# phenoIMPACT: joint phylogenetic slope pilot, offset_multivoltine only.
# Version 1.0.6, 2026-09-29. Requires sdmTMB 1.1.0.
#
# QUICK START (a fresh R session is recommended):
# source("pheno_plasticity_phylogeny_pilot_15km.R")
# phylo_job <- run_phylogeny_pilot()
# watch_phylogeny_pilot()
#
# Optional: run_phylogeny_pilot(check_only = TRUE) validates without a full fit.
# After a failed initialization: phylo_job <- diagnose_phylogeny_pilot()
# then watch_phylogeny_pilot(). This never starts an outer optimizer.
# Send diagnostic_report.txt from the printed Results directory.
# Ctrl+C in watch pauses the display only. Keep R open and the computer awake.
# status_phylogeny_pilot(); stop_phylogeny_pilot()
#
# ONE NEW FIT: reuse the recovered 15-km multivoltine checkpoint as the null.
# Preserve every fitted row, predictor scale, moderator, grouping term, and mesh.
# Add anomaly * u_phylo[species], Cov(u_phylo)=sigma_phylo^2 * A.
# Existing independent species slopes remain: Cov(u_iid)=sigma_iid^2 * I.
# A is the Brownian covariance of the supplied ultrametric tree, diagonal 1.
# Species intercepts remain independent, as in the source model.
# This is a conditional model of slope variation after the existing moderators.
# The first peak is intentionally excluded.
#
# IMPLEMENTATION: a documented mgcv custom basis B=t(chol(A)), full-rank
# identity penalty, no centering and no penalty rescaling. sdmTMB integrates
# its Gaussian coefficients through the standard smooth2random interface.
# This does NOT use the unsupported bs="re" smoother or patch package code.
# The covariance is embedded as serialized bytes in the formula, so sdmTMB
# can reconstruct it exactly without hidden variables in another environment.
# A hexadecimal string avoids decimal rounding and deeply nested byte lists
# when R updates/deparses long formulas.
# The formula uses mgcv::s explicitly so a clean worker needs no attached mgcv.
#
# Before any large fit, mandatory checks compare the transformed basis with A
# and compare TMB's small Gaussian marginal likelihood with an exact dense
# multivariate-normal likelihood, including independent slopes and nonzero means.
# A failed check stops the worker. Serialized native TMB pointers are never used.
# The source likelihood is replayed with a fresh TMB object at its saved optimum
# (no outer refit), verifying both likelihood agreement and a finite gradient.
# The full joint model is built before optimization; its starting likelihood and
# gradient must be finite. One alternative zero-latent start is tried if needed.
# Gaussian identity-link latent effects have a quadratic conditional objective;
# their inner Newton solver uses smartsearch=FALSE, maxit=20.
# Version 1.0.5 enables TMBad sparse-Hessian compression and serial tape creation
# in the isolated worker. This trades speed for memory; the model is unchanged.
# It also releases temporary native TMB objects explicitly and forces Hessian
# construction before evaluating the marginal density, exposing allocation errors
# that TMB's internal Newton wrapper can otherwise turn into a generic NaN.
# Version 1.0.6 fixes backend detection: getFramework() returns an attributed
# string (including openmp), so compare its character value without attributes.
#
# OUTPUT: output/phenology_plasticity/phylogeny_pilot_15km/offset_multivoltine_<hash>/
#   phylogenetic_fraction.csv: sigma_phylo^2/(sigma_phylo^2+sigma_iid^2),
#       with an approximate logit delta-method CI from the joint covariance.
#   variance_components.csv; fixed_effect_comparison.csv; moderator_comparison.csv
#   model_comparison.csv: ML likelihood/AIC comparison, NO automatic LRT p-value.
#   fit_checks.csv; source_fit_checks.csv; preflight_*.csv; model_specification.txt
#   source_replay_*.csv; joint_initial_*.csv: numerical initialization diagnostics.
#   warnings.txt; optimizer_failure.txt/.rds: retained immediately on failure.
#   fit_phylogenetic.rds: checkpoint, saved before any additional optimization.
# Near-zero variance components can invalidate Wald/delta intervals. This pilot
# does not implement a parametric-bootstrap significance test or propagate tree
# uncertainty. A small estimated component is not proof of absence.
# The phylogenetic fraction is not the fraction of total phenological variance.
#
# MEMORY: one worker, one TMB thread, no joint precision matrix requested.
# The original full model object is released before the new fit starts.
# The added basis is dense (n_observations x n_species); memory use is logged.
# Loading the original checkpoint and fitting still require substantial RAM.
# Sources:
# https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/smooth.construct.html
# https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/smooth2random.html
# https://github.com/sdmTMB/sdmTMB/blob/v1.1.0/R/smoothers.R
# https://github.com/sdmTMB/sdmTMB/blob/v1.1.0/src/sdmTMB.cpp
# https://github.com/kaskr/adcomp/blob/master/TMB/R/TMB.R
# Validation: parsed and sourced in R 4.6.0 (WebAssembly); mgcv 1.9-4 plus
# sdmTMB 1.1.0's smoother parser reconstructed the expected covariance exactly
# and reproduced the prediction basis. Native sdmTMB fitting was unavailable.
# Version 1.0.1 also reproduces/fixes the missing-s error with 73 species,
# without attaching mgcv, including formula serialization and smooth removal.
# Version 1.0.2's finite/non-finite/error gates and diagnostic snapshots were
# tested in R; the full spatial TMB model must still be validated on user data.
# Version 1.0.3 adds diagnosis only, not a claimed fix for the nonfinite density.
# Its report and branching logic were tested with controlled R objects; native
# TMB evaluation still requires the actual source checkpoint on the user's PC.
# TMB integration is checked at runtime by the mandatory likelihood test.

phylogeny_pilot_config <- list(
  project_root = "E:/phenoIMPACT project/code/phenoIMPACT",
  source_checkpoint = NULL,
  tree_file = "data/EUROPEAN_BUTTERFLIES_MCC_CA.nwk",
  output_root = NULL,
  expected_run_id = "fd831801913eb3cf70dc64fca70c2116c5826207823d15602f5a843076e7eecd",
  expected_fingerprint = "4089eaf13bde7ef6a911a5103ebdd33b767f519a33aba3c716987429c4fbabae",
  expected_species = 73L,
  gradient_threshold = 0.001,
  extra_optimization_attempts = 1L,
  maximum_single_basis_GiB = 2,
  check_only = FALSE,
  diagnostic_only = FALSE,
  fit_seed = 20260929L
)

pheno_phylogeny_pilot <- local({
  version <- "phylo_slope_pilot_1.0.6"
  anomaly <- "clim_anomaly_tw90"
  moderators <- c("photo_tw90", "clim_background_tw90",
                  "clim_predictability_tw90", "clim_trend_tw90")
  fixed_formula <- stats::as.formula(paste0("OFFSET_mean ~ ", anomaly,
    " * (", paste(moderators, collapse = " + "), ")"), env = baseenv())
  base_formula <- stats::as.formula(paste0(paste(deparse(fixed_formula), collapse = " "),
    " + (1 | SITE_ID) + (1 | SPECIES) + (0 + ", anomaly,
    " | SPECIES_slope) + (1 | site_year_id)"), env = baseenv())
  monitor <- new.env(parent = emptyenv())
  monitor$path <- NULL; monitor$state <- NULL

  save_atomic <- function(x, path) {
    tmp <- tempfile("tmp_", dirname(path), fileext = ".rds")
    on.exit(unlink(tmp), add = TRUE)
    saveRDS(x, tmp, compress = "gzip")
    backup <- paste0(path, ".previous")
    if (file.exists(backup)) unlink(backup)
    if (file.exists(path) && !file.rename(path, backup)) stop("Cannot back up: ", path)
    if (!file.rename(tmp, path)) {
      if (file.exists(backup)) file.rename(backup, path)
      stop("Cannot write checkpoint: ", path)
    }
    if (file.exists(backup)) unlink(backup)
    invisible(path)
  }
  write_csv <- function(x, path) utils::write.csv(x, path, row.names = FALSE, na = "")
  progress <- function(phase) {
    message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", phase)
    monitor$state$phase <- phase
    monitor$state$updated_at <- as.numeric(Sys.time())
    save_atomic(monitor$state, monitor$path)
  }
  normalize_species <- function(x) gsub("[[:space:]_]+", "_", tolower(trimws(x)))
  par_indices <- function(x, name) which(sub("\\[.*$", "", names(x)) == name)

  # Documented extension points; B and the penalty define the covariance exactly.
  construct_phylo <- function(object, data, knots) {
    A <- object$xt$A
    if (!is.matrix(A) || nrow(A) != ncol(A) || any(!is.finite(A))) stop("Invalid A in phylogenetic basis")
    if (max(abs(A - t(A))) > 1e-10) stop("A is not symmetric")
    id <- data[[object$term]]
    if (!is.numeric(id) || anyNA(id) || any(id != as.integer(id)) ||
        any(id < 1 | id > nrow(A))) stop("Invalid PHYLO_ID")
    B <- t(chol(A))
    object$X <- B[as.integer(id), , drop = FALSE]
    object$S <- list(diag(ncol(B)))
    object$rank <- ncol(B); object$null.space.dim <- 0L
    object$df <- ncol(B); object$bs.dim <- ncol(B)
    object$C <- matrix(0, 0L, ncol(B))
    object$no.rescale <- TRUE; object$side.constrain <- FALSE
    object$plot.me <- FALSE; object$te.ok <- 0L
    object$B <- B
    class(object) <- "phenoBM.smooth"
    object
  }
  predict_phylo <- function(object, data) {
    id <- data[[object$term]]
    if (!is.numeric(id) || anyNA(id) || any(id != as.integer(id)) ||
        any(id < 1 | id > nrow(object$B))) stop("Unknown PHYLO_ID in prediction")
    object$B[as.integer(id), , drop = FALSE]
  }
  register_basis <- function() {
    base::registerS3method("smooth.construct", "phenoBM.smooth.spec", construct_phylo,
                          envir = asNamespace("mgcv"))
    base::registerS3method("Predict.matrix", "phenoBM.smooth", predict_phylo,
                          envir = asNamespace("mgcv"))
  }
  add_phylo <- function(f, A, by) {
    bytes <- serialize(unname(A), NULL, version = 2L)
    hex <- paste(sprintf("%02x", as.integer(bytes)), collapse = "")
    last <- paste0(nchar(hex), "L")
    literal <- paste0("unserialize(as.raw(strtoi(substring('", hex,
      "',seq.int(1L,", last, ",2L),seq.int(2L,", last, ",2L)),base=16L)))")
    stats::as.formula(paste0(paste(deparse(f, width.cutoff = 500L), collapse = " "),
      " + mgcv::s(PHYLO_ID, by = ", by, ", bs = 'phenoBM', xt = list(A = ", literal, "))"),
      env = baseenv())
  }
  dependency_check <- function() {
    pkgs <- c("sdmTMB", "TMB", "mgcv", "ape", "Matrix", "digest", "callr", "reformulas")
    missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing)) stop("Install required packages: ", paste(missing, collapse = ", "))
    if (as.character(utils::packageVersion("sdmTMB")) != "1.1.0") {
      stop("This pilot targets sdmTMB 1.1.0. Review compatibility before using another version.")
    }
    register_basis()
    TMB::openmp(n = 1L, DLL = "sdmTMB")
    data.frame(package = pkgs, version = vapply(pkgs, function(p)
      as.character(utils::packageVersion(p)), character(1)))
  }
  fixed_table <- function(fit) {
    terms <- colnames(stats::model.matrix(fixed_formula, fit$data))
    b <- fit$sd_report$par.fixed; V <- fit$sd_report$cov.fixed
    i <- par_indices(b, "b_j")
    if (length(i) != length(terms) || !is.matrix(V) || nrow(V) != length(b)) stop("Cannot align fixed effects")
    se <- sqrt(diag(V)[i]); est <- unname(b[i])
    data.frame(term = terms, estimate = est, SE = se,
      lower_95 = est - 1.96 * se, upper_95 = est + 1.96 * se,
      p_Wald = 2 * stats::pnorm(-abs(est / se)))
  }
  checks <- function(fit, threshold) {
    z <- fixed_table(fit)
    grad <- if (length(fit$gradients) && all(is.finite(fit$gradients))) max(abs(fit$gradients)) else NA_real_
    code <- fit$model$convergence
    if (length(code) != 1L) code <- NA_integer_
    goodSE <- all(is.finite(z$estimate)) && all(is.finite(z$SE) & z$SE > 0)
    pd <- isTRUE(fit$sd_report$pdHess)
    data.frame(convergence_code = code, positive_definite_Hessian = pd,
      max_abs_gradient = grad, finite_fixed_SE = goodSE,
      fit_ok = isTRUE(code == 0L) && pd && is.finite(grad) && grad < threshold && goodSE,
      check_scope = "numerical checks; not a test of spatial or biological adequacy")
  }
  fit_stats <- function(fit) {
    k <- length(fit$model$par); objective <- fit$model$objective
    if (length(objective) != 1L || !is.finite(objective)) stop("Non-finite likelihood")
    data.frame(n_observations = nrow(fit$data), n_parameters = k,
      logLik = -objective, AIC = 2 * objective + 2 * k)
  }
  slope_logsd_index <- function(fit) {
    rc <- fit$split_formula[[1]]$re_cov_terms
    g <- match("SPECIES_slope", names(rc$cnms))
    if (is.na(g) || !identical(unname(rc$cnms[[g]]), anomaly)) stop("Independent species slope mapping differs")
    r <- rc$re_df[rc$re_df$group_indices == g - 1L & rc$re_df$is_sd == 1L, , drop = FALSE]
    j <- par_indices(fit$sd_report$par.fixed, "re_cov_pars")
    if (nrow(r) != 1L || length(j) != nrow(rc$re_df)) stop("Unexpected covariance parameter mapping")
    j[r$par_index + 1L]
  }
  evaluate_start <- function(obj, theta, directory, label, reference_nll = NA_real_) {
    value <- NA_real_; gradient <- rep(NA_real_, length(theta)); error <- ""
    likelihood_evaluated <- FALSE; gradient_evaluated <- FALSE
    tryCatch({
      value <- obj$fn(theta)
      likelihood_evaluated <- TRUE
      if (length(value) != 1L || !is.finite(value)) stop("Initial likelihood is not finite")
      gradient_evaluated <- TRUE
      gradient <- as.numeric(obj$gr(theta))
      if (length(gradient) != length(theta) || any(!is.finite(gradient))) stop("Initial gradient is not finite")
    }, error = function(e) error <<- conditionMessage(e))
    if (length(value) != 1L) value <- NA_real_
    if (length(gradient) != length(theta)) gradient <- rep(NA_real_, length(theta))
    good <- is.finite(value) && all(is.finite(gradient)) && !nzchar(error)
    record <- data.frame(stage = label, nll = value, reference_nll = reference_nll,
      nll_difference = value - reference_nll,
      max_abs_gradient = if (all(is.finite(gradient))) max(abs(gradient)) else NA_real_,
      likelihood_evaluated = likelihood_evaluated, gradient_evaluated = gradient_evaluated,
      finite = good, error = error)
    write_csv(record, file.path(directory, paste0(label, "_evaluation.csv")))
    write_csv(data.frame(index = seq_along(theta), parameter = names(theta),
      value = unname(theta), gradient = gradient, finite_gradient = is.finite(gradient)),
      file.path(directory, paste0(label, "_parameters.csv")))
    list(value = value, gradient = gradient, finite = good, record = record)
  }
  gaussian_inner <- function(obj) {
    TMB::newtonOption(obj, smartsearch = FALSE, maxit = 20L)
    invisible(obj)
  }
  configure_memory <- function(directory) {
    # The registered query is also used by TMB itself. Compression only applies
    # to TMBad; refuse an unsupported backend rather than claiming it is active.
    # The native result carries an 'openmp' attribute even when it is FALSE.
    # Strip attributes before identical(); they do not identify the AD backend.
    framework <- as.character(.Call("getFramework", PACKAGE = "sdmTMB"))
    if (!identical(framework, "TMBad")) stop("This memory configuration requires an sdmTMB binary built with TMBad; detected: ", framework)
    before <- TMB::config(DLL = "sdmTMB")
    settings <- list(tmbad.sparse_hessian_compress = 1L,
      tape.parallel = 0L, optimize.parallel = 0L, autopar = 0L, nthreads = 1L)
    if (!all(names(settings) %in% names(before))) stop("Installed TMB configuration lacks required memory options")
    do.call(TMB::config, c(settings, list(DLL = "sdmTMB")))
    after <- TMB::config(DLL = "sdmTMB")
    verified <- vapply(names(settings), function(nm)
      identical(as.integer(after[[nm]]), settings[[nm]]), logical(1))
    write_csv(data.frame(framework = framework, setting = names(settings),
      before = vapply(names(settings), function(nm) as.integer(before[[nm]]), integer(1)),
      after = vapply(names(settings), function(nm) as.integer(after[[nm]]), integer(1)),
      applied = verified), file.path(directory, "tmb_memory_settings.csv"))
    if (!all(verified)) stop("TMB memory configuration was not applied")
    message("TMBad Hessian compression enabled; serial tapes; one native thread.")
    invisible(after)
  }
  release_tmb <- function(obj) {
    # Only use for live objects created in this worker, never serialized pointers.
    invisible(try(TMB::FreeADFun(obj), silent = TRUE))
  }
  memory_snapshot <- function(directory, label) {
    tryCatch({
      x <- list(stage = label, time = as.character(Sys.time()), R_gc = gc())
      if (requireNamespace("ps", quietly = TRUE)) {
        x$system_bytes <- ps::ps_system_memory()
        x$worker_bytes <- ps::ps_memory_info(ps::ps_handle())
      }
      cat(paste(capture.output(dput(x)), collapse = "\n"), "\n\n",
        file = file.path(directory, "memory_checkpoints.txt"), append = TRUE)
    }, error = function(e) message("Memory snapshot unavailable: ", conditionMessage(e)))
    invisible(NULL)
  }
  prime_joint_hessian <- function(obj, directory) {
    progress("Building compressed random-effect Hessian before marginal evaluation")
    memory_snapshot(directory, "before_Hessian_construction")
    problem <- NULL
    tryCatch({
      # Access forces TMB's delayed sparseHessianFun promise. Do not materialize
      # another Hessian matrix or factorization merely to test construction.
      hessian_function <- obj$env$spHess
      if (!is.function(hessian_function)) stop("TMB sparse Hessian function unavailable")
    }, error = function(e) problem <<- conditionMessage(e))
    memory_snapshot(directory, "after_Hessian_construction")
    if (!is.null(problem)) {
      writeLines(problem, file.path(directory, "hessian_construction_error.txt"))
      # FreeADFun itself accesses spHess. Replace only this failed promise on
      # the discarded object so cleanup cannot retry the failed allocation.
      obj$env$spHess <- local({ ADHess <- NULL; function(...) NULL })
      release_tmb(obj)
      stop("Random-effect Hessian construction failed before optimization: ", problem,
        ". No alternative start was attempted. See hessian_construction_error.txt and memory_checkpoints.txt in: ",
        directory, call. = FALSE)
    }
    invisible(TRUE)
  }
  # Replay the exact saved model at its saved outer optimum. This also recovers
  # current latent modes after extra optimization, without trusting stale parlist.
  recover_start <- function(fit, directory) {
    progress("Replaying the source likelihood at its saved optimum (no refit)")
    parameters <- fit$parlist
    if (is.null(parameters)) parameters <- fit$tmb_params
    if (!is.list(parameters) || !all(names(fit$tmb_params) %in% names(parameters))) stop("Source parameter template unavailable")
    parameters <- parameters[names(fit$tmb_params)]
    for (attempt in 1:2) {
      if (attempt == 2L) {
        progress("Retrying source replay with zero latent starting values")
        for (nm in fit$tmb_random) parameters[[nm]][] <- 0
      }
      obj <- TMB::MakeADFun(data = fit$tmb_data, parameters = parameters,
        map = fit$tmb_map, random = fit$tmb_random, profile = fit$control$profile,
        DLL = "sdmTMB", silent = TRUE)
      gaussian_inner(obj)
      if (!identical(names(obj$par), names(fit$model$par))) stop("Rebuilt source parameter order differs")
      diagnostic <- evaluate_start(obj, fit$model$par, directory,
        paste0("source_replay_", attempt), fit$model$objective)
      if (isTRUE(diagnostic$finite) && abs(diagnostic$value - fit$model$objective) < .01) {
        obj$fn(fit$model$par)
        parameters <- obj$env$parList(par = obj$env$last.par.best)
        keep <- c("b_j", "ln_phi", "ln_tau_E", "ln_kappa", "re_cov_pars", "re_b_pars", "epsilon_st")
        if (!all(keep %in% names(parameters)) ||
            any(!vapply(parameters[keep], function(x) all(is.finite(x)), logical(1)))) stop("Replayed starting values invalid")
        out <- parameters[keep]
        release_tmb(obj); rm(obj); invisible(gc())
        return(out)
      }
      release_tmb(obj); rm(obj); invisible(gc())
    }
    stop("Cannot reproduce source likelihood with a finite gradient. Inspect source_replay_*.csv in: ", directory)
  }
  build_joint <- function(p, start = p$start) {
    control <- sdmTMB::sdmTMBcontrol(start = start, get_joint_precision = FALSE,
      multiphase = FALSE, nlminb_loops = 1L, newton_loops = 0L, trace = 1L)
    fit <- suppressMessages(sdmTMB::sdmTMB(formula = p$formula, data = p$data,
      mesh = p$mesh, time = "YEAR", family = stats::gaussian(link = "identity"),
      spatial = "off", spatiotemporal = "iid", reml = FALSE, priors = p$priors,
      extra_time = p$extra_time, control = control, silent = TRUE, do_fit = FALSE))
    gaussian_inner(fit$tmb_obj)
    fit
  }
  initialize_joint <- function(p) {
    for (attempt in 1:2) {
      start <- p$start
      if (attempt == 2L) {
        start$re_b_pars[] <- 0
        start$epsilon_st[] <- 0
        start$ln_smooth_sigma[] <- start$ln_smooth_sigma + log(3)
      }
      progress(paste("Building and checking full-model start", attempt, "of at most 2"))
      fit <- build_joint(p, start)
      prime_joint_hessian(fit$tmb_obj, p$out)
      diagnostic <- evaluate_start(fit$tmb_obj, fit$tmb_obj$par, p$out,
        paste0("joint_initial_", attempt))
      if (isTRUE(diagnostic$finite)) {
        # run_extra_optimization accepts this live model object. No estimates or
        # successful-convergence claim are exported until its real fit completes.
        fit$model <- list(par = fit$tmb_obj$par, objective = diagnostic$value,
          iterations = 0L, evaluations = c("function" = 0L, "gradient" = 0L))
        fit$sd_report <- list()
        return(fit)
      }
      release_tmb(fit$tmb_obj); rm(fit); invisible(gc())
    }
    stop("Both full-model starting evaluations failed. No outer fit was attempted. Inspect joint_initial_*.csv in: ", p$out)
  }
  # These diagnostics evaluate the existing parameter values. They never call
  # nlminb, optim, or run_extra_optimization, and never create a fitted checkpoint.
  # A marginal evaluation still solves for conditional latent modes internally.
  compare_inputs <- function(old, new, expected = character(), output_file = NULL) {
    fields <- union(names(old), names(new))
    if (!length(fields)) return(data.frame(field = character(), expected_change = logical(),
      same = logical(), source_shape = character(), joint_shape = character(),
      comparison_error = logical(), detail = character()))
    shape <- function(x) paste0(paste(class(x), collapse = "/"), ":",
      if (is.null(dim(x))) length(x) else paste(dim(x), collapse = "x"))
    records <- list()
    for (nm in fields) {
      present <- nm %in% names(old) && nm %in% names(new)
      comparison_error <- FALSE
      comparison <- if (!present) "Field absent" else
        if (!identical(dim(old[[nm]]), dim(new[[nm]]))) "Dimensions differ; numerical subtraction not attempted" else
        tryCatch(all.equal(old[[nm]], new[[nm]], tolerance = 1e-12), error = function(e) {
          comparison_error <<- TRUE
          paste("Comparison could not be evaluated:", conditionMessage(e))
        })
      records[[length(records) + 1L]] <- data.frame(field = nm, expected_change = nm %in% expected,
        same = present && isTRUE(comparison), source_shape = shape(old[[nm]]),
        joint_shape = shape(new[[nm]]), comparison_error = comparison_error,
        detail = if (isTRUE(comparison)) "" else
          substr(paste(comparison, collapse = "; "), 1L, 400L))
      if (!is.null(output_file)) write_csv(do.call(rbind, records), output_file)
    }
    do.call(rbind, records)
  }
  parameter_summary <- function(parameters, gradient = NULL) {
    groups <- unique(names(parameters))
    do.call(rbind, lapply(groups, function(nm) {
      i <- which(names(parameters) == nm); x <- parameters[i]; finite <- is.finite(x)
      data.frame(parameter = nm, n = length(x), nonfinite = sum(!finite),
        minimum = if (any(finite)) min(x[finite]) else NA_real_,
        maximum = if (any(finite)) max(x[finite]) else NA_real_,
        nonfinite_gradient = if (is.null(gradient)) NA_integer_ else sum(!is.finite(gradient[i])))
    }))
  }
  conditional_evaluation <- function(data, parameters, map, directory, label) {
    # random=NULL evaluates the joint density directly; no latent integration.
    obj <- TMB::MakeADFun(data = data, parameters = parameters, map = map,
      random = NULL, profile = NULL, DLL = "sdmTMB", silent = TRUE)
    on.exit({release_tmb(obj); rm(obj); invisible(gc())}, add = TRUE)
    value <- NA_real_; gradient <- NULL; error <- ""
    tryCatch({
      value <- obj$fn(obj$par)
      if (length(value) != 1L || !is.finite(value)) stop("Conditional joint density is not finite")
      gradient <- obj$gr(obj$par)
      if (length(gradient) != length(obj$par) || any(!is.finite(gradient))) stop("Conditional gradient is not finite")
    }, error = function(e) error <<- conditionMessage(e))
    if (length(value) != 1L) value <- NA_real_
    if (!is.null(gradient) && length(gradient) != length(obj$par)) gradient <- NULL
    out <- data.frame(stage = label, nll = value, finite = is.finite(value),
      gradient_finite = !is.null(gradient) && all(is.finite(gradient)), error = error)
    write_csv(out, file.path(directory, paste0(label, ".csv")))
    write_csv(parameter_summary(obj$par, gradient),
      file.path(directory, paste0(label, "_parameter_groups.csv")))
    out
  }
  diagnose_joint <- function(p) {
    report <- file.path(p$out, "diagnostic_report.txt")
    writeLines(c(paste("Phylogeny initialization diagnosis", version),
      "No outer optimization. Conditional densities and one marginal evaluation only.",
      "Blank gradients in earlier CSVs meant not evaluated after a nonfinite likelihood.",
      paste("Results:", p$out)), report)
    note <- function(...) cat(paste(..., collapse = " "), "\n", file = report, append = TRUE)
    stage <- "diagnostic_start"
    mark <- function(label) {
      stage <<- label
      note("STAGE:", label)
      progress(label)
    }
    on.exit({
      for (path in list.files(p$out, pattern = "[.]csv$", full.names = TRUE)) {
        note("\nFILE:", basename(path))
        cat(paste(readLines(path, warn = FALSE), collapse = "\n"), "\n", file = report, append = TRUE)
      }
      for (name in c("memory_checkpoints.txt", "hessian_construction_error.txt")) {
        path <- file.path(p$out, name)
        if (file.exists(path)) {
          note("\nFILE:", name)
          cat(paste(readLines(path, warn = FALSE), collapse = "\n"), "\n", file = report, append = TRUE)
        }
      }
    }, add = TRUE)
    capture <- function(w) {
      note("WARNING:", conditionMessage(w)); invokeRestart("muffleWarning")
    }
    result <- tryCatch(withCallingHandlers({
      mark("Diagnosing reconstruction: building once, without fitting")
      fit <- build_joint(p)
      on.exit(release_tmb(fit$tmb_obj), add = TRUE)
      note("Model construction completed.")
      base <- p$source_template
      mark("Comparing source and expanded TMB data, recording each field")
      td <- compare_inputs(base$data, fit$tmb_data,
        c("has_smooths", "Xs", "Zs", "b_smooth_start"),
        output_file = file.path(p$out, "diagnostic_data_comparison.csv"))
      mark("Comparing source and expanded parameter maps, recording each field")
      tm <- compare_inputs(base$map, fit$tmb_map, c("bs", "b_smooth", "ln_smooth_sigma"),
        output_file = file.path(p$out, "diagnostic_map_comparison.csv"))
      write_csv(td, file.path(p$out, "diagnostic_data_comparison.csv"))
      write_csv(tm, file.path(p$out, "diagnostic_map_comparison.csv"))
      note("Unexpected data differences:", paste(td$field[!td$same & !td$expected_change], collapse = ", "))
      note("Unexpected parameter-map differences:", paste(tm$field[!tm$same & !tm$expected_change], collapse = ", "))
      note("Fields whose comparison raised an error:",
        paste(c(td$field[td$comparison_error], tm$field[tm$comparison_error]), collapse = ", "))
      note("Source random blocks:", paste(base$random, collapse = ", "))
      note("Joint random blocks:", paste(fit$tmb_random, collapse = ", "))
      note("Source profile:", paste(base$profile, collapse = ", "))
      note("Joint profile:", paste(fit$control$profile, collapse = ", "))
      mark("Checking dimensions and finite values of the added basis")
      basis <- fit$tmb_data$Zs[[1L]]
      write_csv(data.frame(rows = nrow(basis), columns = ncol(basis),
        nonfinite = sum(!is.finite(basis)), max_abs = max(abs(basis)),
        zero_columns = sum(colSums(abs(basis)) == 0),
        unpenalized_columns = ncol(fit$tmb_data$Xs)),
        file.path(p$out, "diagnostic_basis.csv"))
      rm(basis)
      for (nm in intersect(names(p$start), names(base$parameters))) base$parameters[[nm]] <- p$start[[nm]]
      mark("Diagnosing source joint density without latent integration")
      original <- conditional_evaluation(base$data, base$parameters, base$map,
        p$out, "diagnostic_source_conditional")
      mark("Diagnosing expanded joint density without latent integration")
      expanded <- conditional_evaluation(fit$tmb_data, fit$tmb_params, fit$tmb_map,
        p$out, "diagnostic_joint_conditional")
      # New coefficients start at zero: only their normal prior should be added.
      stopifnot(all(fit$tmb_params$b_smooth == 0))
      expected <- length(fit$tmb_params$b_smooth) *
        (.5 * log(2 * pi) + as.numeric(fit$tmb_params$ln_smooth_sigma))
      difference <- expanded$nll - original$nll - expected
      write_csv(data.frame(source_joint_nll = original$nll, expanded_joint_nll = expanded$nll,
        expected_added_prior_nll = expected, unexplained_difference = difference,
        tolerance = 1e-5, agreement = is.finite(difference) && abs(difference) < 1e-5),
        file.path(p$out, "diagnostic_conditional_identity.csv"))
      if (!isTRUE(expanded$finite) || !isTRUE(expanded$gradient_finite)) {
        # Turning off the added component is a diagnostic identity check only.
        # This temporary object is never optimized or used for inference.
        mark("Isolating the added smooth term; no optimization")
        disabled_data <- fit$tmb_data; disabled_data$has_smooths <- 0L
        conditional_evaluation(disabled_data, fit$tmb_params, fit$tmb_map,
          p$out, "diagnostic_added_term_disabled")
        note("STOP: failure occurs in the conditional joint model, before latent integration.")
        list(stage = "conditional_failure", report = report)
      } else if (!isTRUE(original$finite) || !is.finite(difference) || abs(difference) >= 1e-5) {
        note("STOP: source-to-expanded conditional identity failed. Inspect data and map differences.")
        list(stage = "reconstruction_difference", report = report)
      } else {
        mark("Diagnosing latent integration at one fixed parameter vector; no outer fit")
        prime_joint_hessian(fit$tmb_obj, p$out)
        initial_full <- fit$tmb_obj$env$last.par
        z <- evaluate_start(fit$tmb_obj, fit$tmb_obj$par, p$out, "diagnostic_marginal")
        write_csv(parameter_summary(fit$tmb_obj$env$last.par),
          file.path(p$out, "diagnostic_latent_after_evaluation.csv"))
        if (!isTRUE(z$finite)) {
          mark("Checking the conditional random-effect Hessian at the finite initial state")
          tryCatch({
            H <- fit$tmb_obj$env$spHess(par = initial_full, random = TRUE)
            hv <- methods::slot(H, "x"); hd <- Matrix::diag(H)
            factor_error <- ""; pd <- FALSE
            if (all(is.finite(hv))) {
              pd <- tryCatch({
                factor <- Matrix::Cholesky(Matrix::forceSymmetric(H), LDL = FALSE, perm = TRUE)
                rm(factor); TRUE
              }, error = function(e) { factor_error <<- conditionMessage(e); FALSE })
            }
            write_csv(data.frame(n = nrow(H), stored_entries = length(hv),
              nonfinite_entries = sum(!is.finite(hv)),
              diagonal_min = min(hd), diagonal_max = max(hd),
              nonpositive_diagonal = sum(hd <= 0, na.rm = TRUE),
              Cholesky_passed = pd, factor_error = factor_error),
              file.path(p$out, "diagnostic_random_hessian.csv"))
            rm(H, hv, hd); invisible(gc())
          }, error = function(e) note("HESSIAN CHECK ERROR:", conditionMessage(e)))
        }
        note(if (isTRUE(z$finite)) "Starting marginal likelihood and gradient passed; no fit was attempted." else
          "Conditional identity passed; failure occurs during marginal evaluation/latent integration.")
        list(stage = if (isTRUE(z$finite)) "initialization_passed" else "marginal_failure", report = report)
      }
    }, warning = capture, error = function(e) {
      note("ERROR CALL:", substr(paste(deparse(conditionCall(e)), collapse = " "), 1L, 600L))
      calls <- tail(sys.calls(), 12L)
      note("CALL STACK (function names):", paste(vapply(calls, function(cl)
        substr(paste(deparse(cl[[1L]]), collapse = " "), 1L, 120L), character(1)), collapse = " -> "))
    }), error = function(e) {
      note("ERROR STAGE:", stage)
      note("DIAGNOSTIC ERROR:", conditionMessage(e))
      list(stage = "diagnostic_error", report = report)
    })
    note("Final diagnostic stage:", result$stage)
    result
  }
  save_optimizer_failure <- function(fit, error, directory) {
    # Pure numeric state only: an unfinished model is never labelled a valid fit.
    e <- fit$tmb_obj$env
    writeLines(conditionMessage(error), file.path(directory, "optimizer_failure.txt"))
    save_atomic(list(error = conditionMessage(error), timestamp = Sys.time(),
      initial_outer = fit$model$par, last_full_parameters = e$last.par,
      best_full_parameters = e$last.par.best, random_indices = e$random),
      file.path(directory, "optimizer_failure.rds"))
  }
  preflight <- function(A, directory) {
    # Check by multiplication, absence of centering, and exact covariance scale.
    n <- nrow(A)
    d <- data.frame(PHYLO_ID = rep(seq_len(n), 2L), x = rep(c(1, -0.7), each = n))
    f <- add_phylo(stats::as.formula("y ~ 1 + x", env = baseenv()), A, "x")
    # sdmTMB removes smooths by updating the formula before model.matrix().
    # Test this with the full species covariance, not only a tiny example.
    d$y <- 0
    fixed_only <- getFromNamespace("remove_s_and_t2", "sdmTMB")(f)
    if (!isTRUE(all.equal(stats::model.matrix(fixed_only, d),
        stats::model.matrix(~ 1 + x, d)))) stop("Phylogenetic smoother was not removed cleanly from the fixed formula")
    parse_sm <- getFromNamespace("parse_smoothers", "sdmTMB")
    sm <- parse_sm(f, d)
    if (length(sm$Zs) != 1L || ncol(sm$Xs) != 0L || ncol(sm$Zs[[1]]) != n) stop("Phylogenetic basis was constrained or split")
    expected <- A[d$PHYLO_ID, d$PHYLO_ID] * outer(d$x, d$x)
    error <- max(abs(tcrossprod(sm$Zs[[1]]) - expected))
    if (!is.finite(error) || error > 1e-8) stop("Phylogenetic covariance reconstruction failed: ", error)
    write_csv(data.frame(n_species = n, maximum_covariance_error = error,
      unpenalized_columns = ncol(sm$Xs), formula_removal_passed = TRUE,
      passed = TRUE), file.path(directory, "preflight_basis.csv"))
    # Exact Gaussian likelihood check on a small known tree, including IID slopes.
    At <- kronecker(diag(3), matrix(c(1, .65, .65, 1), 2, 2))
    tiny <- data.frame(PHYLO_ID = rep(1:6, each = 8), x = rep(seq(-1, 1, length.out = 8), 6))
    tiny$g <- factor(tiny$PHYLO_ID)
    tiny$y <- sin(seq_len(nrow(tiny)) * .7) + .3 * tiny$x
    form <- add_phylo(stats::as.formula("y ~ 1 + x + (0 + x | g)", env = baseenv()), At, "x")
    obj <- sdmTMB::sdmTMB(form, data = tiny, family = stats::gaussian(), spatial = "off",
      spatiotemporal = "off", reml = FALSE, do_fit = FALSE, silent = TRUE,
      control = sdmTMB::sdmTMBcontrol(get_joint_precision = FALSE, newton_loops = 0L))
    gaussian_inner(obj$tmb_obj)
    output <- list()
    for (j in 1:2) {
      b <- if (j == 1L) c(0, 0) else c(.7, -.3)
      sig <- if (j == 1L) c(1, 1, 1) else c(1.2, .5, .9) # residual, IID, phylo
      theta <- obj$tmb_obj$par
      specs <- list(b_j = b, ln_phi = log(sig[1]), re_cov_pars = log(sig[2]), ln_smooth_sigma = log(sig[3]))
      for (nm in names(specs)) {
        i <- par_indices(theta, nm)
        if (length(i) != length(specs[[nm]])) stop("Small likelihood-test parameter mismatch: ", nm)
        theta[i] <- specs[[nm]]
      }
      actual <- obj$tmb_obj$fn(theta)
      gradient_finite <- all(is.finite(obj$tmb_obj$gr(theta)))
      id <- tiny$PHYLO_ID
      V <- sig[1]^2 * diag(nrow(tiny)) + outer(tiny$x, tiny$x) *
        (sig[3]^2 * At[id, id] + sig[2]^2 * outer(id, id, "=="))
      R <- chol(V); residual <- tiny$y - b[1] - b[2] * tiny$x
      whitened <- forwardsolve(t(R), residual)
      exact <- .5 * (nrow(tiny) * log(2 * pi) + 2 * sum(log(diag(R))) + sum(whitened^2))
      difference <- abs(actual - exact)
      output[[j]] <- data.frame(setting = j, TMB_nll = actual, exact_nll = exact,
        absolute_error = difference, gradient_finite = gradient_finite,
        passed = is.finite(difference) && difference < 1e-6 && gradient_finite)
    }
    output <- do.call(rbind, output)
    write_csv(output, file.path(directory, "preflight_likelihood.csv"))
    if (!all(output$passed)) stop("Exact small-model likelihood check failed; no large fit attempted")
    release_tmb(obj$tmb_obj); rm(obj, sm); invisible(gc())
  }
  prepare <- function(cfg, packages) {
    progress("Reading the recovered multivoltine checkpoint")
    saved <- readRDS(cfg$source_checkpoint)
    if (!identical(saved$run_id, cfg$expected_run_id) ||
        !identical(saved$fingerprint, cfg$expected_fingerprint) || !inherits(saved$fit, "sdmTMB")) {
      stop("Unexpected source checkpoint; the recovered 15-km multivoltine fit is required")
    }
    fit <- saved$fit
    if (as.character(fit$version) != "1.1.0" || isTRUE(fit$reml) ||
        !identical(fit$time, "YEAR") || !all(fit$spatial == "off") ||
        !all(fit$spatiotemporal == "iid") || isTRUE(fit$family$delta) ||
        !all(fit$family$family == "gaussian") || !all(fit$family$link == "identity") ||
        !is.null(fit$time_varying) || !is.null(fit$spatial_varying) ||
        isTRUE(as.logical(fit$tmb_data$has_smooths)) || !is.null(fit$nonlocal_formula)) {
      stop("Source model structure/version differs from the intended pilot")
    }
    f0 <- if (inherits(fit$formula, "formula")) fit$formula else fit$formula[[1]]
    if (!identical(all.vars(f0), all.vars(base_formula)) ||
        !setequal(attr(stats::terms(reformulas::nobars(f0)), "term.labels"),
                  attr(stats::terms(fixed_formula), "term.labels"))) stop("Source formula differs")
    bar_text <- function(f) sort(vapply(reformulas::findbars(f), function(z) paste(deparse(z), collapse = ""), character(1)))
    if (!identical(bar_text(f0), bar_text(base_formula))) stop("Source random-effect terms differ")
    d <- fit$data
    cols <- unique(c(all.vars(base_formula), "source_row", "x_km", "y_km"))
    if (!all(cols %in% names(d)) || anyNA(d[, cols]) || anyDuplicated(d$source_row)) stop("Invalid source rows")
    if (any(!vapply(d[, c("OFFSET_mean", anomaly, moderators, "x_km", "y_km")],
                    function(x) is.numeric(x) && all(is.finite(x)), logical(1)))) stop("Non-finite source values")
    if (!identical(as.character(d$SPECIES), as.character(d$SPECIES_slope))) stop("Species grouping differs")
    if (!is.null(fit$offset) && any(fit$offset != 0)) stop("Unexpected model offset")
    if (!is.null(fit$tmb_data$weights) && any(fit$tmb_data$weights != 1)) stop("Unexpected observation weights")
    ck <- checks(fit, cfg$gradient_threshold)
    if (!isTRUE(ck$fit_ok)) stop("Source fit fails numerical checks")
    species <- levels(droplevels(factor(d$SPECIES)))
    if (length(species) != cfg$expected_species) stop("Unexpected number of species: ", length(species))
    tree <- ape::read.tree(cfg$tree_file)
    if (!inherits(tree, "phylo") || !ape::is.rooted(tree) || !ape::is.ultrametric(tree) ||
        is.null(tree$edge.length) || any(!is.finite(tree$edge.length)) || any(tree$edge.length <= 0)) stop("Tree is not a valid positive-branch chronogram")
    sp <- normalize_species(species); tips <- normalize_species(tree$tip.label)
    if (anyDuplicated(sp) || anyDuplicated(tips) || anyNA(match(sp, tips))) stop("Species/tree matching failed; no fuzzy matching is performed")
    tip_names <- tree$tip.label[match(sp, tips)]
    tr <- ape::keep.tip(tree, tip_names)
    C <- ape::vcv.phylo(tr)[tip_names, tip_names, drop = FALSE]
    if (diff(range(diag(C))) / mean(diag(C)) > 1e-7) stop("Pruned tree is not ultrametric")
    A <- C / mean(diag(C)); dimnames(A) <- list(species, species)
    chol(A)
    d$PHYLO_ID <- match(as.character(d$SPECIES), species)
    X <- stats::model.matrix(fixed_formula, d)
    if (qr(X)$rank != ncol(X)) stop("Fixed effects are rank deficient")
    if (!identical(unname(fit$response), unname(d$OFFSET_mean))) {
      if (!isTRUE(all.equal(as.numeric(fit$response), d$OFFSET_mean, tolerance = 1e-12))) stop("Source response differs")
    }
    mesh <- fit$spde
    if (nrow(mesh$loc_xy) != nrow(d) || !isTRUE(all.equal(unname(as.matrix(mesh$loc_xy)),
        unname(as.matrix(d[, c("x_km", "y_km")])), tolerance = 1e-10))) stop("Saved mesh projection does not match source rows")
    basis_gib <- 8 * nrow(d) * length(species) / 1024^3
    if (basis_gib > cfg$maximum_single_basis_GiB) stop("Added basis exceeds configured memory guard: ", round(basis_gib, 2), " GiB")
    start <- recover_start(fit, cfg$monitor_directory)
    iid <- exp(fit$sd_report$par.fixed[slope_logsd_index(fit)])
    start$ln_smooth_sigma <- matrix(log(max(.1 * iid, .01)), 1L, 1L)
    source <- list(fixed = fixed_table(fit), stats = fit_stats(fit), checks = ck,
      saved_checks = saved$checks, row_hash = digest::digest(d$source_row, algo = "sha256"),
      run_id = saved$run_id, fingerprint = saved$fingerprint,
      checkpoint = cfg$source_checkpoint, formula = f0)
    signature <- digest::digest(list(version, source$fingerprint, data = d[, cols],
      tree_covariance = A, mesh = mesh$mesh, packages = packages,
      threshold = cfg$gradient_threshold, extra = cfg$extra_optimization_attempts,
      seed = cfg$fit_seed), algo = "sha256")
    out <- file.path(cfg$output_root, paste0("offset_multivoltine_", substr(signature, 1L, 12L)))
    if (isTRUE(cfg$diagnostic_only)) out <- tempfile("diagnosis_", tmpdir = out)
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    replay_files <- list.files(cfg$monitor_directory,
      pattern = "^(source_replay_.*|tmb_memory_settings)[.]csv$", full.names = TRUE)
    if (length(replay_files)) file.copy(replay_files, out, overwrite = TRUE)
    write_csv(source$checks, file.path(out, "source_fit_checks.csv"))
    if (!is.null(source$saved_checks)) write_csv(source$saved_checks, file.path(out, "source_saved_checks.csv"))
    write_csv(packages, file.path(out, "package_versions.csv"))
    write_csv(data.frame(SPECIES = species, tip_label = tip_names, PHYLO_ID = seq_along(species)), file.path(out, "species_mapping.csv"))
    write_csv(data.frame(n_observations = nrow(d), n_species = length(species),
      n_sites = length(unique(d$SITE_ID)), n_years = length(unique(d$YEAR)),
      mesh_reused_exactly = TRUE, single_added_basis_GiB = basis_gib,
      note = "Several basis copies plus original spatial model memory may be needed; this is not total RAM"), file.path(out, "input_summary.csv"))
    save_atomic(list(A = A, tree = tr, source = source, signature = signature), file.path(out, "provenance.rds"))
    f <- add_phylo(base_formula, A, anomaly)
    writeLines(c("Cov(species slope) = sigma_phylo^2 * A + sigma_iid^2 * I; diagonal(A)=1.",
      "Only the species slope gains a phylogenetic component; species intercepts remain IID.",
      "All environmental interactions and the annual IID spatial field are preserved.",
      "One source MCC-CA tree. Predictor units unchanged. Gaussian ML. No peak analysis.",
      "Likelihood/AIC comparison only; no automatic chi-square variance-component p-value.",
      "Intervals are approximate and unreliable at variance boundaries or failed convergence.",
      paste("Source checkpoint:", cfg$source_checkpoint), paste("Code:", version),
      paste("One added dense basis:", round(basis_gib, 3), "GiB; total memory is larger.")),
      file.path(out, "model_specification.txt"))
    priors <- fit$priors; extra_time <- fit$extra_time
    source_template <- if (isTRUE(cfg$diagnostic_only)) list(data = fit$tmb_data,
      parameters = fit$tmb_params, map = fit$tmb_map, random = fit$tmb_random,
      profile = fit$control$profile) else NULL
    rm(fit, saved, X, C, tree); invisible(gc())
    list(data = d, A = A, formula = f, mesh = mesh, start = start, source = source,
         out = out, signature = signature, priors = priors, extra_time = extra_time,
         source_template = source_template)
  }
  export_results <- function(saved, p, cfg) {
    fit <- saved$fit; out <- p$out
    ck <- checks(fit, cfg$gradient_threshold)
    write_csv(ck, file.path(out, "fit_checks.csv"))
    if (!isTRUE(ck$fit_ok)) {
      writeLines("Numerical checks failed. Checkpoint retained; inference is suppressed.", file.path(out, "inference_status.txt"))
      return(invisible(FALSE))
    }
    b <- fit$sd_report$par.fixed; V <- fit$sd_report$cov.fixed
    ip <- par_indices(b, "ln_smooth_sigma"); ii <- slope_logsd_index(fit)
    if (length(ip) != 1L || any(!is.finite(V[c(ip, ii), c(ip, ii)]))) stop("Variance uncertainty unavailable")
    z <- b[c(ii, ip)]; se <- sqrt(diag(V)[c(ii, ip)])
    if (any(!is.finite(se))) stop("Invalid variance-component SE")
    variances <- data.frame(component = c("independent_species_slope", "phylogenetic_species_slope"),
      SD = exp(z), variance = exp(2 * z), log_SD_SE = se,
      SD_lower_95 = exp(z - 1.96 * se), SD_upper_95 = exp(z + 1.96 * se))
    write_csv(variances, file.path(out, "variance_components.csv"))
    eta <- 2 * (b[ip] - b[ii])
    v_eta <- 4 * (V[ip, ip] + V[ii, ii] - 2 * V[ip, ii])
    if (!is.finite(v_eta) || v_eta < -1e-10) stop("Invalid variance-fraction uncertainty")
    se_eta <- sqrt(max(0, v_eta))
    near_boundary <- stats::plogis(eta) < .001 || stats::plogis(eta) > .999
    write_csv(data.frame(phylogenetic_fraction_species_slope = stats::plogis(eta),
      lower_95_approx = stats::plogis(eta - 1.96 * se_eta),
      upper_95_approx = stats::plogis(eta + 1.96 * se_eta),
      logit_SE = se_eta, near_variance_boundary = near_boundary,
      interval = "Joint-covariance delta method on logit scale; unreliable near boundaries",
      interpretation = "Fraction of conditional between-species slope variance; not total phenological variance"),
      file.path(out, "phylogenetic_fraction.csv"))
    fx <- fixed_table(fit)
    comparison <- merge(p$source$fixed, fx, by = "term", suffixes = c("_original", "_phylogenetic"), sort = FALSE)
    comparison$estimate_change <- comparison$estimate_phylogenetic - comparison$estimate_original
    comparison$change_in_original_SE <- comparison$estimate_change / comparison$SE_original
    comparison$sign_changed <- sign(comparison$estimate_original) != sign(comparison$estimate_phylogenetic)
    comparison$SE_ratio <- comparison$SE_phylogenetic / comparison$SE_original
    write_csv(comparison, file.path(out, "fixed_effect_comparison.csv"))
    terms <- unique(c(paste(anomaly, moderators, sep = ":"), paste(moderators, anomaly, sep = ":")))
    interaction <- comparison[comparison$term %in% terms, , drop = FALSE]
    if (nrow(interaction) != 4L) stop("Expected four moderator interactions")
    write_csv(interaction, file.path(out, "moderator_comparison.csv"))
    new <- fit_stats(fit); old <- p$source$stats
    LR <- 2 * (new$logLik - old$logLik)
    if (new$n_parameters - old$n_parameters != 1L) stop("Models differ by an unexpected parameter count")
    valid <- is.finite(LR) && LR >= -1e-3
    write_csv(data.frame(logLik_original = old$logLik, logLik_phylogenetic = new$logLik,
      AIC_original = old$AIC, AIC_phylogenetic = new$AIC,
      delta_AIC_original_minus_phylogenetic = old$AIC - new$AIC,
      likelihood_ratio = if (valid) max(0, LR) else NA_real_,
      likelihood_comparison_valid = valid, p_LRT = NA_real_,
      reason_no_p = "Zero variance is a boundary; calibrated significance needs a suitable profile/bootstrap analysis"),
      file.path(out, "model_comparison.csv"))
    writeLines(if (valid) "Pilot estimates available. Check variance boundaries, uncertainty and source spatial adequacy before interpretation." else
      "New likelihood is below the nested source likelihood: optimization comparison failed; interpret neither variance nor coefficient differences yet.",
      file.path(out, "inference_status.txt"))
    invisible(valid)
  }
  worker <- function(cfg) {
    monitor$path <- cfg$progress_file; monitor$state <- readRDS(monitor$path)
    complete <- FALSE
    on.exit({
      monitor$state$finished <- TRUE; monitor$state$failed <- !complete
      monitor$state$ended_at <- as.numeric(Sys.time())
      if (!complete) monitor$state$phase <- "Stopped with error; inspect worker.log"
      save_atomic(monitor$state, monitor$path)
    }, add = TRUE)
    packages <- dependency_check()
    configure_memory(cfg$monitor_directory)
    writeLines(capture.output(sessionInfo()), file.path(cfg$monitor_directory, "sessionInfo.txt"))
    p <- prepare(cfg, packages)
    monitor$state$result_directory <- p$out
    if (isTRUE(cfg$diagnostic_only)) {
      result <- diagnose_joint(p)
      progress(paste("Diagnosis finished:", result$stage, "- send diagnostic_report.txt"))
      complete <- TRUE
      return(list(output_dir = p$out, diagnostic_only = TRUE,
        stage = result$stage, report = result$report))
    }
    progress("Checking exact phylogenetic covariance and small-model likelihood")
    preflight(p$A, p$out)
    if (cfg$check_only) {
      progress("Preflight passed; full model not fitted")
      complete <- TRUE
      return(list(output_dir = p$out, check_only = TRUE))
    }
    path <- file.path(p$out, "fit_phylogenetic.rds")
    if (file.exists(path)) {
      progress("Reading matching pilot checkpoint")
      cached <- readRDS(path)
      if (!identical(cached$signature, p$signature)) stop("Pilot checkpoint fingerprint differs")
      cached_checks <- checks(cached$fit, cfg$gradient_threshold)
      write_csv(cached_checks, file.path(p$out, "fit_checks.csv"))
      if (!isTRUE(cached_checks$fit_ok)) {
        stop("Pilot checkpoint failed numerical checks. It is retained; inspect fit_checks.csv before scheduling a retry.")
      }
      okay <- export_results(cached, p, cfg)
      progress(if (isTRUE(okay)) "Completed from cached pilot" else "Cached fit has diagnostic limitations; inspect inference_status.txt")
      complete <- TRUE
      return(list(output_dir = p$out, cached = TRUE, fit_ok = cached_checks$fit_ok, inference_ready = okay))
    }
    set.seed(cfg$fit_seed)
    warnings <- character()
    collect <- function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      writeLines(unique(warnings), file.path(p$out, "warnings.txt"))
      message("Model warning: ", conditionMessage(w))
      invokeRestart("muffleWarning")
    }
    started <- proc.time()[["elapsed"]]
    fit <- withCallingHandlers(initialize_joint(p), warning = collect)
    progress("Fitting joint phylogenetic + independent species slopes; initial gradient verified")
    monitor$state$fit_started_at <- as.numeric(Sys.time())
    save_atomic(monitor$state, monitor$path)
    fit <- tryCatch(withCallingHandlers(sdmTMB::run_extra_optimization(fit,
      nlminb_loops = 1L, newton_loops = 0L), warning = collect), error = function(e) {
        save_optimizer_failure(fit, e, p$out)
        stop("Optimization failed: ", conditionMessage(e),
             ". Diagnostics retained in: ", p$out, call. = FALSE)
      })
    fit$parlist <- fit$tmb_obj$env$parList(par = fit$tmb_obj$env$last.par.best)
    fit$last.par.best <- fit$tmb_obj$env$last.par.best
    saved <- list(signature = p$signature, fit = fit, checks = checks(fit, cfg$gradient_threshold),
      warnings = unique(warnings), elapsed_minutes = (proc.time()[["elapsed"]] - started) / 60,
      source = p$source, A = p$A, code_version = version)
    progress("Saving the initial fit before further optimization")
    save_atomic(saved, path)
    for (attempt in seq_len(cfg$extra_optimization_attempts)) {
      if (isTRUE(saved$checks$fit_ok)) break
      progress(paste("Additional optimization attempt", attempt))
      improved <- tryCatch(withCallingHandlers(sdmTMB::run_extra_optimization(fit,
        nlminb_loops = 1L, newton_loops = 0L), warning = collect), error = function(e) {
          warnings <<- c(warnings, conditionMessage(e)); NULL
        })
      if (is.null(improved)) break
      fit <- improved; saved$fit <- fit; saved$checks <- checks(fit, cfg$gradient_threshold)
      fit$parlist <- fit$tmb_obj$env$parList(par = fit$tmb_obj$env$last.par.best)
      fit$last.par.best <- fit$tmb_obj$env$last.par.best
      saved$fit <- fit
      saved$warnings <- unique(warnings)
      saved$elapsed_minutes <- (proc.time()[["elapsed"]] - started) / 60
      save_atomic(saved, path)
    }
    if (!identical(fit$data$source_row, p$data$source_row)) stop("Fitted row order changed")
    writeLines(unique(warnings), file.path(p$out, "warnings.txt"))
    progress("Exporting variance partition and moderator comparison")
    okay <- export_results(saved, p, cfg)
    progress(if (isTRUE(okay)) "Completed; review results and uncertainty" else "Completed fit with diagnostic limitations; inspect inference_status.txt")
    complete <- TRUE
    list(output_dir = p$out, fit_ok = saved$checks$fit_ok, inference_ready = okay)
  }
  start <- function(cfg) {
    if (!requireNamespace("callr", quietly = TRUE)) stop("Install callr first")
    if (is.null(cfg$diagnostic_only)) cfg$diagnostic_only <- FALSE
    if (!is.logical(cfg$diagnostic_only) || length(cfg$diagnostic_only) != 1L ||
        is.na(cfg$diagnostic_only)) stop("Invalid diagnostic_only")
    if (!is.logical(cfg$check_only) || length(cfg$check_only) != 1L || is.na(cfg$check_only)) stop("Invalid check_only")
    if (!is.numeric(cfg$extra_optimization_attempts) || length(cfg$extra_optimization_attempts) != 1L ||
        !is.finite(cfg$extra_optimization_attempts) || cfg$extra_optimization_attempts < 0 ||
        cfg$extra_optimization_attempts %% 1 != 0) stop("Invalid extra_optimization_attempts")
    if (!is.finite(cfg$gradient_threshold) || cfg$gradient_threshold <= 0 ||
        !is.finite(cfg$maximum_single_basis_GiB) || cfg$maximum_single_basis_GiB <= 0) stop("Invalid numerical/memory threshold")
    cfg$project_root <- normalizePath(cfg$project_root, winslash = "/", mustWork = TRUE)
    absolute <- function(x) grepl("^([A-Za-z]:[/\\\\]|/|\\\\\\\\)", x)
    if (is.null(cfg$source_checkpoint)) cfg$source_checkpoint <- file.path(cfg$project_root,
      "output/phenology_plasticity/spatial_models/offset_recovery_15km",
      "primary__offset_multivoltine_fd831801913e/fit_spatial_field.rds")
    if (!absolute(cfg$source_checkpoint)) cfg$source_checkpoint <- file.path(cfg$project_root, cfg$source_checkpoint)
    if (!absolute(cfg$tree_file)) cfg$tree_file <- file.path(cfg$project_root, cfg$tree_file)
    cfg$source_checkpoint <- normalizePath(cfg$source_checkpoint, winslash = "/", mustWork = TRUE)
    cfg$tree_file <- normalizePath(cfg$tree_file, winslash = "/", mustWork = TRUE)
    if (is.null(cfg$output_root)) cfg$output_root <- file.path(cfg$project_root, "output/phenology_plasticity/phylogeny_pilot_15km")
    if (!absolute(cfg$output_root)) cfg$output_root <- file.path(cfg$project_root, cfg$output_root)
    dir.create(cfg$output_root, recursive = TRUE, showWarnings = FALSE)
    cfg$output_root <- normalizePath(cfg$output_root, winslash = "/", mustWork = TRUE)
    for (option in c("pheno.phylo.pilot.job", "pheno.latitude.active_job", "pheno.spatial.active_job")) {
      old <- getOption(option)
      if (!is.null(old) && isTRUE(tryCatch(old$process$is_alive(), error = function(e) FALSE))) {
        stop("A model worker is already running in this session: ", option)
      }
    }
    cfg$monitor_directory <- tempfile(paste0("monitor_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_"), cfg$output_root)
    dir.create(cfg$monitor_directory)
    cfg$progress_file <- file.path(cfg$monitor_directory, "progress.rds")
    save_atomic(list(started_at = as.numeric(Sys.time()), updated_at = as.numeric(Sys.time()),
      ended_at = NA_real_, fit_started_at = NA_real_, phase = "Starting worker",
      finished = FALSE, failed = FALSE, result_directory = NA_character_), cfg$progress_file)
    log <- file.path(cfg$monitor_directory, "worker.log")
    process <- callr::r_bg(worker, args = list(cfg = cfg), libpath = .libPaths(),
      stdout = log, stderr = log, supervise = TRUE, user_profile = FALSE,
      system_profile = FALSE, package = TRUE, wd = cfg$project_root)
    job <- list(process = process, progress_file = cfg$progress_file,
      log_file = log, monitor_directory = cfg$monitor_directory)
    options(pheno.phylo.pilot.job = job)
    message(if (cfg$diagnostic_only) "Started initialization diagnosis: no outer optimization or model refit." else
      if (cfg$check_only) "Started preflight worker." else "Started one phylogeny pilot worker: mandatory preflight, then one new fit.",
      "\nUse watch_phylogeny_pilot(). Keep R open and the computer awake.\nLog: ", log)
    invisible(job)
  }
  status <- function(job, quiet = FALSE) {
    if (is.null(job)) stop("Run run_phylogeny_pilot() first")
    s <- tryCatch(readRDS(job$progress_file), error = function(e) NULL)
    if (is.null(s)) return(invisible(NULL))
    alive <- job$process$is_alive()
    if (!alive && !s$finished) { s$phase <- "Worker stopped unexpectedly; inspect worker.log"; s$failed <- TRUE }
    now <- if (isTRUE(s$finished) && is.finite(s$ended_at)) s$ended_at else as.numeric(Sys.time())
    if (!quiet) {
      cat(sprintf("Elapsed: %.2f h | worker running: %s\n", (now - s$started_at) / 3600, alive))
      cat("Current:", s$phase, "\n")
      if (!is.na(s$result_directory)) cat("Results:", s$result_directory, "\n")
      if (!alive && isTRUE(s$failed)) cat("Log:", job$log_file, "\n")
    }
    invisible(s)
  }
  watch <- function(job, interval) {
    if (is.null(job)) stop("Run run_phylogeny_pilot() first")
    if (length(interval) != 1L || !is.finite(interval) || interval < 1) stop("Invalid interval")
    tryCatch({
      repeat {
        status(job)
        if (!job$process$is_alive()) break
        until <- Sys.time() + interval
        while (Sys.time() < until && job$process$is_alive()) Sys.sleep(1)
      }
      invisible(job$process$get_result())
    }, interrupt = function(e) {
      message("Display paused; fitting continues. Use watch_phylogeny_pilot() to resume.")
      invisible(job)
    })
  }
  stop_job <- function(job) {
    if (is.null(job)) return(invisible(FALSE))
    if (job$process$is_alive()) job$process$kill_tree()
    message("Worker stopped. Completed checkpoints retained; an unfinished fit must restart.")
    invisible(TRUE)
  }
  read_fit <- function(path) {
    dependency_check()
    saved <- readRDS(path)
    if (!identical(saved$code_version, version)) stop("Unexpected pilot checkpoint")
    saved
  }
  list(start = start, status = status, watch = watch, stop = stop_job, read_fit = read_fit,
       register = register_basis, preflight = preflight)
})

run_phylogeny_pilot <- function(cfg = phylogeny_pilot_config, check_only = cfg$check_only) {
  cfg$check_only <- check_only
  pheno_phylogeny_pilot$start(cfg)
}
diagnose_phylogeny_pilot <- function(cfg = phylogeny_pilot_config) {
  cfg$diagnostic_only <- TRUE
  cfg$check_only <- FALSE
  pheno_phylogeny_pilot$start(cfg)
}
status_phylogeny_pilot <- function(job = getOption("pheno.phylo.pilot.job"), quiet = FALSE) {
  pheno_phylogeny_pilot$status(job, quiet)
}
watch_phylogeny_pilot <- function(job = getOption("pheno.phylo.pilot.job"), interval_seconds = 60) {
  pheno_phylogeny_pilot$watch(job, interval_seconds)
}
stop_phylogeny_pilot <- function(job = getOption("pheno.phylo.pilot.job")) pheno_phylogeny_pilot$stop(job)
read_phylogeny_pilot_fit <- function(path) pheno_phylogeny_pilot$read_fit(path)

message("Phylogeny pilot v1.0.6 functions loaded. Fit with compressed Hessian: run_phylogeny_pilot(). Diagnosis only: diagnose_phylogeny_pilot(). Then watch_phylogeny_pilot().")
