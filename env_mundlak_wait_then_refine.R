# phenoIMPACT | WAIT -> REFINE environmental Mundlak fits
# 2026-10-17
#
# Start this SECOND queue while the current environmental Mundlak queue is running.
# It waits for the source worker to stop, then sequentially refines ONLY fits marked
# NUMERICAL_REVIEW_REQUIRED. Scientific model/data/mesh are unchanged.
# Refinement uses nlminb only: newton_loops = 0. Original fit_mundlak.rds files
# and original result tables are NEVER overwritten.
#
# Run:
#   source(file.choose(), encoding = "UTF-8")
#   env_refine_job <- pheno_env_refine$run(monitor = FALSE)
# Status:
#   pheno_env_refine$status()
#
# The computer must remain ON and not sleep/hibernate.

pheno_env_refine <- local({
  VERSION <- "environment_mundlak_wait_refine_1.0"
  ROOT <- "E:/phenoIMPACT project/code/phenoIMPACT"
  DEFAULT_SOURCE_RUN <- file.path(
    ROOT, "output", "phenology_plasticity",
    "environment_mundlak_by_voltinism_qc_final_15km",
    "run_de0890c85d4a"
  )
  ANALYSES <- c(
    "offset_multivoltine", "onset_multivoltine",
    "offset_univoltine", "onset_univoltine"
  )
  state <- new.env(parent = emptyenv())

  need <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)
  csv <- function(x, p) utils::write.csv(x, p, row.names = FALSE, na = "")

  atomic_rds <- function(x, p) {
    dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
    tmp <- tempfile(".writing_", tmpdir = dirname(p), fileext = ".rds")
    on.exit(unlink(tmp), add = TRUE)
    saveRDS(x, tmp, compress = "gzip")
    bak <- paste0(p, ".previous")
    if (file.exists(bak)) unlink(bak)
    had <- file.exists(p)
    if (had) need(file.rename(p, bak), paste("Cannot back up", p))
    if (!file.rename(tmp, p)) {
      if (had && file.exists(bak)) file.rename(bak, p)
      stop("Cannot finalize ", p, call. = FALSE)
    }
    if (file.exists(bak)) unlink(bak)
    invisible(p)
  }

  phase <- function(cfg, stage, analysis = "", detail = "", terminal = FALSE) {
    z <- list(stage = stage, analysis = analysis, detail = detail,
              time = Sys.time(), terminal = terminal)
    atomic_rds(z, file.path(cfg$out, "progress.rds"))
    message("[", format(z$time, "%Y-%m-%d %H:%M:%S"), "] ", stage,
            if (nzchar(analysis)) paste0(" | ", analysis) else "",
            if (nzchar(detail)) paste0(" | ", detail) else "")
  }

  check_packages <- function() {
    pkgs <- c("sdmTMB", "TMB", "Matrix", "fmesher", "callr", "ps")
    miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
    need(!length(miss), paste("Missing package(s):", paste(miss, collapse = ", ")))
    need(as.character(utils::packageVersion("sdmTMB")) == "1.1.0",
         "Use the original sdmTMB 1.1.0 library.")
    invisible(TRUE)
  }

  pid_alive <- function(pid, create_time) {
    if (is.null(pid) || is.null(create_time)) return(FALSE)
    tryCatch({
      h <- ps::ps_handle(pid)
      isTRUE(ps::ps_is_running(h)) &&
        isTRUE(all.equal(as.numeric(ps::ps_create_time(h)),
                         as.numeric(create_time), tolerance = 1e-6))
    }, error = function(e) FALSE)
  }

  source_alive <- function(run) {
    f <- file.path(run, "worker_owner.rds")
    if (!file.exists(f)) return(FALSE)
    x <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(x)) return(FALSE)
    pid_alive(x$pid, x$create_time)
  }

  ram_free <- function() as.numeric(ps::ps_system_memory()$avail) / 1024^3
  disk_free <- function(path) as.numeric(ps::ps_disk_usage(path)$available) / 1024^3

  wait_source <- function(cfg) {
    repeat {
      alive <- source_alive(cfg$source_run)
      p <- file.path(cfg$source_run, "progress.rds")
      s <- if (file.exists(p)) tryCatch(readRDS(p), error = function(e) NULL) else NULL
      if (!alive) break
      phase(cfg, "WAITING_FOR_SOURCE_QUEUE",
            if (!is.null(s)) s$analysis else "",
            paste0("source stage = ", if (!is.null(s)) s$stage else "unknown",
                   "; next check in ", cfg$wait_seconds, " s"))
      Sys.sleep(cfg$wait_seconds)
    }
    phase(cfg, "SOURCE_QUEUE_FINISHED", "", "Source worker is no longer running")
    Sys.sleep(60)
    invisible(TRUE)
  }

  wait_ram <- function(cfg, analysis) {
    repeat {
      a <- ram_free()
      if (is.finite(a) && a >= cfg$minimum_free_ram_GiB) break
      phase(cfg, "WAITING_FOR_RAM", analysis,
            paste(round(a, 1), "GiB free; require", cfg$minimum_free_ram_GiB, "GiB"))
      Sys.sleep(cfg$wait_seconds)
    }
    phase(cfg, "RAM_READY", analysis, paste(round(ram_free(), 1), "GiB available"))
    invisible(TRUE)
  }

  sync_fit <- function(fit) {
    fit$tmb_obj$fn(fit$model$par)
    fit$parlist <- fit$tmb_obj$env$parList(par = fit$tmb_obj$env$last.par.best)
    fit$last.par.best <- fit$tmb_obj$env$last.par.best
    fit
  }

  fit_checks <- function(fit) {
    g <- tryCatch(as.numeric(fit$tmb_obj$gr(fit$model$par)), error = function(e) NA_real_)
    mg <- if (length(g) && all(is.finite(g))) max(abs(g)) else Inf
    conv <- if (length(fit$model$convergence) == 1L) fit$model$convergence else NA_integer_
    pd <- isTRUE(fit$sd_report$pdHess)
    ok <- isTRUE(conv == 0L) && pd && is.finite(mg) && mg < .001
    data.frame(
      convergence_code = conv, pdHess = pd, max_abs_gradient = mg,
      numerical_checks_ok = ok, objective = fit$model$objective,
      logLik = tryCatch(as.numeric(stats::logLik(fit)), error = function(e) NA_real_),
      AIC = tryCatch(stats::AIC(fit), error = function(e) NA_real_),
      nobs = nrow(fit$data), stringsAsFactors = FALSE
    )
  }

  rebuild_fit <- function(fit, p, out) {
    rebuilt <- sdmTMB::sdmTMB(
      formula = p$formula, data = p$data, mesh = p$mesh, time = "YEAR",
      family = stats::gaussian(link = "identity"), spatial = "off",
      spatiotemporal = "iid", reml = FALSE, previous_fit = fit,
      do_fit = FALSE, silent = TRUE,
      control = sdmTMB::sdmTMBcontrol(
        multiphase = FALSE, nlminb_loops = 1L, newton_loops = 0L,
        parallel = 1L, get_joint_precision = FALSE,
        collapse_spatial_variance = FALSE
      )
    )
    need(identical(names(rebuilt$tmb_obj$par), names(fit$model$par)) &&
           length(rebuilt$tmb_obj$par) == length(fit$model$par),
         "Rebuilt outer parameter vector differs from saved fit.")
    val <- rebuilt$tmb_obj$fn(fit$model$par)
    need(is.finite(val) && is.finite(fit$model$objective) &&
           abs(val - fit$model$objective) < .05,
         paste("Saved/rebuilt likelihood mismatch:", fit$model$objective, val))
    csv(data.frame(saved_objective = fit$model$objective,
                   rebuilt_objective = val,
                   absolute_difference = abs(val - fit$model$objective)),
        file.path(out, "source_rebuild_check.csv"))
    fit$tmb_obj <- rebuilt$tmb_obj
    fit$lower <- rebuilt$lower
    fit$upper <- rebuilt$upper
    fit$model$objective <- val
    rm(rebuilt); invisible(gc())
    sync_fit(fit)
  }

  refine_one <- function(cfg, analysis) {
    src <- file.path(cfg$source_run, analysis)
    donef <- file.path(src, "completed.rds")
    fitf <- file.path(src, "fit_mundlak.rds")
    prepf <- file.path(src, "prepared_input.rds")

    if (!all(file.exists(c(donef, fitf, prepf)))) {
      return(data.frame(analysis = analysis, source_status = "NOT_COMPLETED",
                        refinement_status = "SKIPPED_NOT_AVAILABLE",
                        initial_gradient = NA, final_gradient = NA,
                        passes_used = 0L, stringsAsFactors = FALSE))
    }
    done <- readRDS(donef)
    source_status <- as.character(done$status)
    if (identical(source_status, "FIT_OK")) {
      return(data.frame(analysis = analysis, source_status = source_status,
                        refinement_status = "SKIPPED_ALREADY_FIT_OK",
                        initial_gradient = NA, final_gradient = NA,
                        passes_used = 0L, stringsAsFactors = FALSE))
    }
    if (!identical(source_status, "NUMERICAL_REVIEW_REQUIRED")) {
      return(data.frame(analysis = analysis, source_status = source_status,
                        refinement_status = "SKIPPED_UNEXPECTED_STATUS",
                        initial_gradient = NA, final_gradient = NA,
                        passes_used = 0L, stringsAsFactors = FALSE))
    }

    wait_ram(cfg, analysis)
    need(disk_free(cfg$out) >= cfg$minimum_free_disk_GiB,
         paste("Insufficient disk before refining", analysis))

    out <- file.path(cfg$out, analysis)
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    completed_refine <- file.path(out, "refinement_completed.rds")
    if (file.exists(completed_refine)) {
      old <- readRDS(completed_refine)
      if (identical(old$source_status, source_status)) return(old$summary)
    }

    phase(cfg, "LOADING_FOR_REFINEMENT", analysis, "Loading fit_mundlak.rds")
    p <- readRDS(prepf)
    saved <- readRDS(fitf)
    need(is.list(saved) && inherits(saved$fit, "sdmTMB"), "Invalid saved Mundlak fit.")
    fit <- rebuild_fit(saved$fit, p, out)
    ch <- fit_checks(fit)
    initial_gradient <- ch$max_abs_gradient

    hist <- data.frame(
      pass = 0L, objective = ch$objective, max_abs_gradient = ch$max_abs_gradient,
      convergence_code = ch$convergence_code, pdHess = ch$pdHess,
      numerical_checks_ok = ch$numerical_checks_ok,
      accepted_as_current = TRUE, timestamp = as.character(Sys.time()),
      stringsAsFactors = FALSE
    )
    csv(hist, file.path(out, "refinement_history.csv"))
    passes_used <- 0L
    final_status <- if (isTRUE(ch$numerical_checks_ok)) "FIT_OK_AFTER_REBUILD" else "NUMERICAL_REVIEW_REQUIRED"

    if (!isTRUE(ch$numerical_checks_ok)) for (pass in seq_len(cfg$max_passes)) {
      phase(cfg, "REFINING", analysis,
            paste("nlminb pass", pass, "of", cfg$max_passes,
                  "| gradient", signif(ch$max_abs_gradient, 5),
                  "| objective", signif(fit$model$objective, 12)))
      warnings <- character()
      catcher <- function(w) {
        warnings <<- unique(c(warnings, conditionMessage(w)))
        message("WARNING: ", conditionMessage(w)); invokeRestart("muffleWarning")
      }
      cand <- tryCatch(
        withCallingHandlers(
          sdmTMB::run_extra_optimization(fit, nlminb_loops = 1L, newton_loops = 0L),
          warning = catcher
        ),
        error = function(e) NULL
      )
      if (is.null(cand)) {
        final_status <- "REFINEMENT_ERROR"; break
      }
      cand <- sync_fit(cand)
      cc <- fit_checks(cand)
      atomic_rds(list(version = VERSION, analysis = analysis, pass = pass,
                      warnings = warnings, checks = cc, fit = cand),
                 file.path(out, sprintf("fit_refinement_%02d.rds", pass)))
      accept <- is.finite(cc$objective) &&
        cc$objective <= fit$model$objective + cfg$objective_tolerance
      hist <- rbind(hist, data.frame(
        pass = pass, objective = cc$objective, max_abs_gradient = cc$max_abs_gradient,
        convergence_code = cc$convergence_code, pdHess = cc$pdHess,
        numerical_checks_ok = cc$numerical_checks_ok,
        accepted_as_current = accept, timestamp = as.character(Sys.time()),
        stringsAsFactors = FALSE
      ))
      csv(hist, file.path(out, "refinement_history.csv"))
      if (!accept) { final_status <- "STOPPED_OBJECTIVE_WORSE"; rm(cand); gc(); break }
      rm(fit); fit <- cand; ch <- cc; passes_used <- pass; gc()
      if (isTRUE(ch$numerical_checks_ok)) { final_status <- "FIT_OK_REFINED"; break }
      final_status <- "NUMERICAL_REVIEW_REQUIRED_AFTER_MAX_PASSES"
    }

    fit <- sync_fit(fit)
    final_checks <- fit_checks(fit)
    atomic_rds(list(version = VERSION, source_run = cfg$source_run,
                    source_signature = saved$signature, analysis = analysis,
                    passes_used = passes_used, refinement_status = final_status,
                    checks = final_checks, fit = fit),
               file.path(out, "fit_refined_final.rds"))
    csv(final_checks, file.path(out, "fit_checks_refined.csv"))
    atomic_rds(fit$parlist, file.path(out, "parameters_refined.rds"))

    summary <- data.frame(
      analysis = analysis, source_status = source_status,
      refinement_status = final_status,
      initial_gradient = initial_gradient,
      final_gradient = final_checks$max_abs_gradient,
      passes_used = passes_used,
      stringsAsFactors = FALSE
    )
    atomic_rds(list(source_status = source_status, summary = summary, time = Sys.time()),
               completed_refine)
    rm(fit, p, saved); invisible(gc())
    summary
  }

  controller <- function(cfg) {
    check_packages(); options(sdmTMB.cores = 1L); TMB::openmp(n = 1L, DLL = "sdmTMB")
    tryCatch({
      phase(cfg, "WAITING_FOR_SOURCE_QUEUE", "", "Second queue active; no fit loaded")
      wait_source(cfg)
      rows <- list()
      for (a in ANALYSES) {
        phase(cfg, "CHECKING", a, "Inspecting source status")
        rows[[a]] <- tryCatch(
          refine_one(cfg, a),
          error = function(e) {
            msg <- conditionMessage(e)
            writeLines(
              c(paste0("ERROR: ", msg), paste0("TIME: ", Sys.time())),
              file.path(cfg$out, paste0(a, "_refinement_error.txt"))
            )
            data.frame(
              analysis = a, source_status = "UNKNOWN_AFTER_ERROR",
              refinement_status = paste0("REFINEMENT_ERROR: ", msg),
              initial_gradient = NA, final_gradient = NA,
              passes_used = 0L, stringsAsFactors = FALSE
            )
          }
        )
        tab <- do.call(rbind, rows)
        csv(tab, file.path(cfg$out, "refinement_summary_partial.csv"))
        invisible(gc())
      }
      tab <- do.call(rbind, rows)
      csv(tab, file.path(cfg$out, "refinement_summary.csv"))
      final <- if (all(tab$refinement_status %in% c(
        "SKIPPED_ALREADY_FIT_OK", "SKIPPED_NOT_AVAILABLE",
        "FIT_OK_AFTER_REBUILD", "FIT_OK_REFINED"
      ))) "COMPLETED" else "COMPLETED_WITH_NUMERICAL_REVIEW"
      atomic_rds(list(status = final, summary = tab, source_run = cfg$source_run,
                      output = cfg$out, time = Sys.time()),
                 file.path(cfg$out, "result.rds"))
      phase(cfg, final, "", "Refinement queue finished", terminal = TRUE)
    }, error = function(e) {
      msg <- conditionMessage(e)
      try(writeLines(c(paste0("ERROR: ", msg), paste0("TIME: ", Sys.time())),
                     file.path(cfg$out, "fatal_error.txt")), silent = TRUE)
      try(phase(cfg, "STOPPED_WITH_ERROR", "", msg, terminal = TRUE), silent = TRUE)
      atomic_rds(list(status = "STOPPED_WITH_ERROR", error = msg,
                      source_run = cfg$source_run, output = cfg$out, time = Sys.time()),
                 file.path(cfg$out, "result.rds"))
    })
  }

  dump_engine <- function(out) {
    e <- environment(dump_engine)
    n <- ls(e, all.names = TRUE)
    f <- n[vapply(n, function(k) is.function(get(k, e)), logical(1))]
    p <- file.path(out, "refinement_engine.R")
    dump(c("VERSION", "ROOT", "DEFAULT_SOURCE_RUN", "ANALYSES", f),
         file = p, envir = e, control = "all")
    normalizePath(p, winslash = "/", mustWork = TRUE)
  }

  recover <- function(source_run = DEFAULT_SOURCE_RUN) {
    source_run <- normalizePath(source_run, winslash = "/", mustWork = TRUE)
    out <- file.path(source_run, "postfit_refinement_queue_v1")
    need(dir.exists(out), "No refinement queue directory exists yet.")
    ownf <- file.path(out, "worker_owner.rds")
    own <- if (file.exists(ownf)) tryCatch(readRDS(ownf), error = function(e) NULL) else NULL
    logs <- list.files(out, pattern = "^refinement_controller_.*\\.log$", full.names = TRUE)
    log <- if (length(logs)) logs[which.max(file.info(logs)$mtime)] else NA_character_
    list(process = NULL,
         pid = if (!is.null(own)) own$pid else NULL,
         create_time = if (!is.null(own)) own$create_time else NULL,
         out = out, log = log,
         started = if (!is.null(own)) own$started else NA)
  }

  run <- function(source_run = DEFAULT_SOURCE_RUN, max_passes = 4L,
                  minimum_free_ram_GiB = 35, minimum_free_disk_GiB = 40,
                  wait_seconds = 120, objective_tolerance = .001,
                  monitor = FALSE) {
    check_packages()
    source_run <- normalizePath(source_run, winslash = "/", mustWork = TRUE)
    need(max_passes >= 1L && max_passes <= 10L, "max_passes must be 1..10.")
    out <- file.path(source_run, "postfit_refinement_queue_v1")
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    ownf <- file.path(out, "worker_owner.rds")
    if (file.exists(ownf)) {
      old <- tryCatch(readRDS(ownf), error = function(e) NULL)
      if (!is.null(old) && pid_alive(old$pid, old$create_time)) {
        message("Refinement queue already running.")
        job <- list(process = NULL, pid = old$pid, create_time = old$create_time,
                    out = out, log = old$log, started = old$started)
        state$job <- job; return(invisible(job))
      }
    }
    engine <- dump_engine(out)
    cfg <- list(source_run = source_run, out = normalizePath(out, winslash = "/"),
                max_passes = as.integer(max_passes),
                minimum_free_ram_GiB = minimum_free_ram_GiB,
                minimum_free_disk_GiB = minimum_free_disk_GiB,
                wait_seconds = wait_seconds, objective_tolerance = objective_tolerance)
    atomic_rds(cfg, file.path(out, "refinement_config.rds"))
    log <- file.path(out, paste0("refinement_controller_",
                                 format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
    child <- callr::r_bg(function(engine, cfg) {
      e <- new.env(parent = globalenv()); sys.source(engine, envir = e); e$controller(cfg)
    }, args = list(engine = engine, cfg = cfg), stdout = log, stderr = "2>&1",
    supervise = FALSE, user_profile = FALSE, system_profile = FALSE)
    h <- ps::ps_handle(child$get_pid())
    own <- list(pid = child$get_pid(), create_time = ps::ps_create_time(h),
                started = Sys.time(), out = out, log = log)
    atomic_rds(own, ownf)
    job <- list(process = child, pid = own$pid, create_time = own$create_time,
                out = out, log = log, started = own$started)
    state$job <- job
    message("WAIT -> REFINE queue launched.")
    message("It will sleep until the current source worker stops.")
    message("Monitor with: pheno_env_refine$status()")
    if (monitor) repeat {
      status(job); if (!pid_alive(job$pid, job$create_time)) break; Sys.sleep(60)
    }
    invisible(job)
  }

  status <- function(job = NULL) {
    if (is.null(job)) job <- state$job
    if (is.null(job)) { job <- recover(); message("Recovered refinement queue from disk: ", job$out) }
    pf <- file.path(job$out, "progress.rds")
    if (file.exists(pf)) print(tryCatch(readRDS(pf), error = function(e) NULL))
    alive <- if (!is.null(job$process)) tryCatch(job$process$is_alive(), error = function(e) FALSE)
             else pid_alive(job$pid, job$create_time)
    cat("\nAlive:", alive, "\n")
    sf <- if (file.exists(file.path(job$out, "refinement_summary.csv")))
      file.path(job$out, "refinement_summary.csv") else
      file.path(job$out, "refinement_summary_partial.csv")
    if (file.exists(sf)) {
      cat("\nRefinement summary:\n")
      print(utils::read.csv(sf, stringsAsFactors = FALSE), row.names = FALSE)
    }
    if (!alive && file.exists(file.path(job$out, "fatal_error.txt"))) {
      cat("\nError:\n"); cat(readLines(file.path(job$out, "fatal_error.txt"), warn = FALSE), sep = "\n")
    }
    invisible(job)
  }

  list(run = run, status = status, recover = recover)
})

run_environment_mundlak_refinement_queue <- function(...) pheno_env_refine$run(...)
status_environment_mundlak_refinement_queue <- function(...) pheno_env_refine$status(...)
