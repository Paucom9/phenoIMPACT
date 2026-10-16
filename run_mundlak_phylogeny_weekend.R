# phenoIMPACT | Sequential queue for the five Mundlak phylogenetic sensitivities.
# Loading defines functions only. No package installation or model starts.
# First: audit_mundlak_phylogeny_weekend(project_root = "...")
# Then: job <- start_mundlak_phylogeny_weekend(check_only = FALSE)
# Keep R open and Windows awake. A paused watch does not stop the queue.
# Repeating start resumes via each model's validated checkpoints; failed native
# Hessian calls must be rebuilt. One model worker runs at a time.

pheno_mundlak_phylogeny_weekend <- local({
  VERSION <- "mundlak_phylogeny_weekend_1.0.0"
  ROOT <- "E:/phenoIMPACT project/code/phenoIMPACT"
  ENV_SCRIPT <- "11c_environment_mundlak_phylogeny_15km.R"
  DEM_SCRIPT <- "10g_mundlak_demography_phylogeny_15km.R"
  ACTIVE_OPTIONS <- c("phenoimpact.env.mundlak.job", "phenoimpact.mundlak.qc.job",
    "pheno.env.mundlak.phylo.job", "pheno.demography.mundlak.phylo.job",
    "pheno.demo.mundlak.phylo.job",
    "pheno.phylo.qc.job", "pheno.phylo.pilot.job", "pheno.mundlak.phylogeny.weekend.job")
  need <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
  tasks <- function() data.frame(
    task = c("environment_onset", "environment_offset_univoltine",
      "environment_offset_multivoltine", "demography_univoltine", "demography_multivoltine"),
    family = c(rep("environment", 3L), rep("demography", 2L)),
    analysis = c("onset", "offset_univoltine", "offset_multivoltine", "univoltine", "multivoltine"),
    stringsAsFactors = FALSE)
  absolute <- function(path, root, mustWork = TRUE) {
    if (!grepl("^([A-Za-z]:[/\\\\]|/|\\\\\\\\)", path)) path <- file.path(root, path)
    normalizePath(path, winslash = "/", mustWork = mustWork)
  }
  load_api <- function(root, family) {
    env <- new.env(parent = globalenv())
    path <- file.path(root, if (family == "environment") ENV_SCRIPT else DEM_SCRIPT)
    need(file.exists(path), paste("Missing queue dependency:", path))
    suppressMessages(sys.source(path, envir = env))
    env
  }
  save_atomic <- function(value, path) {
    tmp <- tempfile("queue_", dirname(path))
    on.exit(unlink(tmp), add = TRUE)
    saveRDS(value, tmp)
    # Windows cannot atomically replace an existing file with file.rename.
    if (file.exists(path)) unlink(path)
    need(file.rename(tmp, path), paste("Cannot write queue metadata:", path))
    invisible(path)
  }
  process_matches <- function(pid, created) {
    if (!is.numeric(pid) || length(pid) != 1L || !is.finite(pid)) return(NA)
    tryCatch({
      if (!pid %in% ps::ps_pids()) return(FALSE)
      handle <- ps::ps_handle(pid)
      if (ps::ps_status(handle) %in% c("zombie", "dead")) return(FALSE)
      if (!is.numeric(created) || length(created) != 1L || !is.finite(created)) return(NA)
      isTRUE(abs(as.numeric(ps::ps_create_time(handle)) - created) < 0.01) &&
        isTRUE(ps::ps_is_running(handle))
    }, error = function(e) NA)
  }
  owner_self <- function(token) list(token = token, pid = Sys.getpid(),
    created = as.numeric(ps::ps_create_time(ps::ps_handle())), child_pid = NULL,
    child_created = NULL)
  lock_owned <- function(lock, token) {
    owner <- tryCatch(readRDS(file.path(lock, "owner.rds")), error = function(e) NULL)
    !is.null(owner) && identical(owner$token, token)
  }
  release_lock <- function(lock, token) {
    if (lock_owned(lock, token)) unlink(lock, recursive = TRUE)
    invisible(NULL)
  }
  acquire_lock <- function(lock, token) {
    if (dir.exists(lock)) {
      old <- tryCatch(readRDS(file.path(lock, "owner.rds")), error = function(e) NULL)
      need(!is.null(old), paste("Cannot identify existing queue lock. Inspect:", lock))
      alive <- process_matches(old$pid, old$created)
      child_alive <- if (is.null(old$child_pid)) FALSE else
        process_matches(old$child_pid, old$child_created)
      need(identical(alive, FALSE) && identical(child_alive, FALSE),
        paste("A queue controller/model may still be running. Inspect:", lock))
      unlink(lock, recursive = TRUE)
    }
    need(dir.create(lock), paste("Cannot acquire queue lock:", lock))
    save_atomic(owner_self(token), file.path(lock, "owner.rds"))
  }
  known_workers_idle <- function() {
    for (option in ACTIVE_OPTIONS) {
      old <- getOption(option)
      if (!is.null(old)) {
        alive <- tryCatch(old$process$is_alive(), error = function(e) NA)
        need(identical(alive, FALSE), paste("Another registered model/controller may be active:", option))
      }
    }
    invisible(TRUE)
  }
  audit_impl <- function(project_root, environment_source_run, demography_source_run) {
    root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
    sources <- list(environment = environment_source_run, demography = demography_source_run)
    t <- tasks()
    rows <- candidates <- vector("list", nrow(t))
    apis <- list()
    for (family in unique(t$family)) apis[[family]] <- tryCatch(load_api(root, family), error = identity)
    for (i in seq_len(nrow(t))) {
      family <- t$family[[i]]
      api <- apis[[family]]
      tab <- if (inherits(api, "error")) api else tryCatch({
        fun <- get(if (family == "environment") "audit_environment_mundlak_phylogeny_sources" else
          "audit_demography_mundlak_phylogeny_sources", envir = api, inherits = FALSE)
        args <- list(project_root = root, source_run = sources[[family]])
        args[[if (family == "environment") "analysis" else "group"]] <- t$analysis[[i]]
        value <- NULL
        invisible(capture.output(value <- suppressMessages(do.call(fun, args))))
        need(is.data.frame(value) && all(c("source_directory", "accepted", "reason") %in% names(value)),
          "Unexpected source-audit schema; no model was started.")
        value
      }, error = identity)
      okay <- if (inherits(tab, "error")) integer() else which(!is.na(tab$accepted) & tab$accepted)
      reason <- if (inherits(tab, "error")) conditionMessage(tab) else if (length(okay) == 1L) "eligible FIT_OK source" else
        if (length(okay) > 1L) "multiple eligible sources; specify the independent source_run" else
          paste(unique(tab$reason), collapse = " | ")
      rows[[i]] <- data.frame(task = t$task[[i]], family = family, analysis = t$analysis[[i]],
        accepted_sources = length(okay), ready = length(okay) == 1L,
        source_directory = if (length(okay) == 1L) as.character(tab$source_directory[[okay]]) else NA_character_,
        reason = reason, stringsAsFactors = FALSE)
      if (!inherits(tab, "error")) {
        tab$queue_task <- t$task[[i]]
        candidates[[i]] <- tab
      }
    }
    result <- do.call(rbind, rows)
    attr(result, "source_audits") <- candidates
    attr(result, "all_ready") <- all(result$ready)
    attr(result, "validation_scope") <- "Metadata eligibility only; no fitted model is opened and no native call or compilation runs."
    result
  }
  audit <- function(project_root = ROOT, environment_source_run = NULL, demography_source_run = NULL) {
    result <- audit_impl(project_root, environment_source_run, demography_source_run)
    print(result, row.names = FALSE)
    message("Auditoria de les cinc fonts independents; cap model iniciat.")
    invisible(result)
  }
  record_table <- function(records) {
    if (!length(records)) return(data.frame())
    do.call(rbind, lapply(records, function(x) data.frame(task = x$task, family = x$family,
      analysis = x$analysis, status = x$status, started_at = x$started_at,
      ended_at = x$ended_at, output_dir = x$output_dir, log_file = x$log_file,
      reason = x$reason, stringsAsFactors = FALSE)))
  }
  controller <- function(cfg) {
    need(lock_owned(cfg$lock, cfg$token), "Queue launch lock has changed.")
    # The parent records the controller PID before allowing any model launch.
    # This also makes a failed R subprocess startup safely resumable.
    deadline <- Sys.time() + 15
    while (!file.exists(cfg$handoff_file)) {
      need(Sys.time() < deadline, "Queue controller launch handoff was not completed.")
      Sys.sleep(0.1)
    }
    owner <- owner_self(cfg$token)
    save_atomic(owner, file.path(cfg$lock, "owner.rds"))
    active_job <- NULL
    on.exit({
      if (!is.null(active_job) && isTRUE(tryCatch(active_job$process$is_alive(), error = function(e) FALSE)))
        try(active_job$process$kill_tree(), silent = TRUE)
      release_lock(cfg$lock, cfg$token)
    }, add = TRUE)
    records <- list()
    state <- list(stage = "STARTING", current_task = NA_character_, updated_at = as.character(Sys.time()),
      terminal = FALSE, current_log = NA_character_, current_output = NA_character_,
      attempt_directory = cfg$attempt_directory, controller_pid = Sys.getpid())
    update <- function(...) {
      state <<- utils::modifyList(state, list(...))
      state$updated_at <<- as.character(Sys.time())
      save_atomic(state, cfg$progress_file)
    }
    on.exit(if (!isTRUE(state$terminal)) {
      state$stage <- "CONTROLLER_STOPPED"; state$terminal <- TRUE
      state$updated_at <- as.character(Sys.time())
      try(save_atomic(state, cfg$progress_file), silent = TRUE)
    }, add = TRUE)
    for (i in seq_len(nrow(cfg$tasks))) {
      task <- cfg$tasks[i, , drop = FALSE]
      record <- list(task = task$task[[1L]], family = task$family[[1L]], analysis = task$analysis[[1L]],
        status = "STARTING", started_at = as.character(Sys.time()), ended_at = NA_character_,
        output_dir = NA_character_, log_file = NA_character_, reason = "", result = NULL)
      update(stage = "STARTING_MODEL", current_task = record$task, current_log = NA_character_, current_output = NA_character_, model_progress = NULL)
      api <- NULL
      error <- tryCatch({
        api <- load_api(cfg$project_root, record$family)
        fun <- get(if (record$family == "environment") "run_environment_mundlak_phylogeny" else
          "run_demography_mundlak_phylogeny", envir = api, inherits = FALSE)
        args <- list(project_root = cfg$project_root, source_run = task$source_directory[[1L]], check_only = FALSE)
        args[[if (record$family == "environment") "analysis" else "group"]] <- record$analysis
        args <- c(args, cfg$model_arguments[[record$family]])
        active_job <- do.call(fun, args)
        need(!is.null(active_job$process), "Model runner did not return a worker process.")
        record$log_file <- active_job$log_file
        owner$child_pid <- active_job$process$get_pid()
        owner$child_created <- if (active_job$process$is_alive())
          as.numeric(ps::ps_create_time(ps::ps_handle(owner$child_pid))) else NA_real_
        save_atomic(owner, file.path(cfg$lock, "owner.rds"))
        save_atomic(list(task = record$task, monitor_directory = active_job$monitor_directory,
          log_file = active_job$log_file, progress_file = active_job$progress_file),
          file.path(cfg$attempt_directory, paste0(record$task, "_worker.rds")))
        update(stage = "RUNNING_MODEL", current_log = record$log_file)
        while (active_job$process$is_alive()) {
          p <- tryCatch(readRDS(active_job$progress_file), error = function(e) NULL)
          if (!is.null(p)) {
            directory <- p$result_directory
            if (is.character(directory) && length(directory) == 1L && !is.na(directory)) record$output_dir <- directory
            update(model_progress = p, current_output = record$output_dir)
          }
          Sys.sleep(cfg$poll_seconds)
        }
        result <- active_job$process$get_result()
        need(is.list(result), "Unexpected model result; inspect the worker log.")
        record$result <- result
        if (is.character(result$output_dir) && length(result$output_dir) == 1L) record$output_dir <- result$output_dir
        passed <- isTRUE(result$inference_ready) &&
          (isTRUE(result$fit_ok) || (is.null(result$fit_ok) && isTRUE(result$cached))) &&
          !isTRUE(result$check_only)
        record$status <- if (passed) "FIT_OK" else "INFERENCE_NOT_READY"
        if (!passed) record$reason <- "Fit/checkpoints retained; the model did not report inference_ready=TRUE. Inspect its numerical checks and log."
        NULL
      }, error = identity)
      if (inherits(error, "error")) {
        record$status <- "WORKER_FAILED"
        record$reason <- conditionMessage(error)
      }
      # Never advance while a previous worker could still consume model RAM.
      if (!is.null(active_job) && isTRUE(tryCatch(active_job$process$is_alive(), error = function(e) FALSE))) {
        try(active_job$process$kill_tree(), silent = TRUE)
        active_job$process$wait(timeout = 10000)
        need(!active_job$process$is_alive(), "Cannot stop previous model worker; refusing to overlap models.")
      }
      if (!is.null(active_job)) {
        p <- tryCatch(readRDS(active_job$progress_file), error = function(e) NULL)
        if (is.na(record$output_dir) && !is.null(p) && is.character(p$result_directory) &&
            length(p$result_directory) == 1L) record$output_dir <- p$result_directory
      }
      record$ended_at <- as.character(Sys.time())
      records[[length(records) + 1L]] <- record
      save_atomic(record, file.path(cfg$attempt_directory, paste0(record$task, "_result.rds")))
      utils::write.csv(record_table(records), file.path(cfg$attempt_directory, "queue_results.csv"), row.names = FALSE)
      utils::write.csv(record_table(records), file.path(cfg$output_directory, "queue_results_latest.csv"), row.names = FALSE)
      owner$child_pid <- owner$child_created <- NULL
      save_atomic(owner, file.path(cfg$lock, "owner.rds"))
      for (option in ACTIVE_OPTIONS) options(structure(list(NULL), names = option))
      active_job <- api <- NULL
      invisible(gc())
      update(stage = "MODEL_FINISHED", model_status = record$status, model_reason = record$reason)
      if (cfg$stop_on_error && record$status != "FIT_OK") break
    }
    okay <- length(records) == nrow(cfg$tasks) && all(vapply(records, function(x) identical(x$status, "FIT_OK"), logical(1)))
    update(stage = if (okay) "ALL_FIT_OK" else "COMPLETED_WITH_FAILURES", terminal = TRUE, current_task = NA_character_)
    result <- list(all_fit_ok = okay, output_directory = cfg$output_directory,
      attempt_directory = cfg$attempt_directory, results = record_table(records))
    save_atomic(result, file.path(cfg$attempt_directory, "queue_result.rds"))
    result
  }
  start <- function(project_root = ROOT, check_only = TRUE,
    environment_source_run = NULL, demography_source_run = NULL,
    stop_on_error = FALSE, minimum_free_GiB = NULL,
    environment_arguments = list(), demography_arguments = list(), poll_seconds = 5) {
    need(is.logical(check_only) && length(check_only) == 1L && !is.na(check_only), "Invalid check_only.")
    need(is.logical(stop_on_error) && length(stop_on_error) == 1L && !is.na(stop_on_error), "Invalid stop_on_error.")
    need(is.numeric(poll_seconds) && length(poll_seconds) == 1L && is.finite(poll_seconds) && poll_seconds >= 1,
      "poll_seconds must be at least one second.")
    root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
    audit_table <- audit(root, environment_source_run, demography_source_run)
    if (check_only) return(invisible(audit_table))
    need(all(audit_table$ready), "All five independent sources must be uniquely FIT_OK before the queue can start. No model was started.")
    for (package in c("callr", "ps", "digest")) need(requireNamespace(package, quietly = TRUE), paste("Existing package required:", package))
    required <- file.path(root, c(ENV_SCRIPT, DEM_SCRIPT, "pheno_phylogeny_qc_15km_v2.R", "pheno_plasticity_phylogeny_pilot_15km.R"))
    need(all(file.exists(required)), paste("Missing queue scripts:", paste(required[!file.exists(required)], collapse = ", ")))
    known_workers_idle()
    need(is.list(environment_arguments) && is.list(demography_arguments), "Model arguments must be lists.")
    arguments <- list(environment = environment_arguments, demography = demography_arguments)
    for (family in names(arguments)) {
      x <- arguments[[family]]
      need(!length(x) || (!is.null(names(x)) && all(nzchar(names(x))) && !anyDuplicated(names(x))), "Name each model argument once.")
      need(!any(names(x) %in% c("project_root", "check_only", "source_run", "analysis", "group")), "Queue-controlled model arguments cannot be overridden.")
      if (!is.null(minimum_free_GiB)) {
        need(is.numeric(minimum_free_GiB) && length(minimum_free_GiB) == 1L && is.finite(minimum_free_GiB) && minimum_free_GiB > 0,
          "Invalid minimum_free_GiB.")
        need(!"minimum_free_GiB" %in% names(x), "Set minimum_free_GiB once, globally or in each family's arguments.")
        x$minimum_free_GiB <- minimum_free_GiB
      }
      arguments[[family]] <- x
    }
    signature <- digest::digest(list(version = VERSION, tasks = audit_table[c("task", "family", "analysis", "source_directory")],
      scripts = unname(tools::md5sum(required)), arguments = arguments), algo = "xxhash64")
    parent <- file.path(root, "output", "mundlak_phylogeny_weekend")
    dir.create(parent, recursive = TRUE, showWarnings = FALSE)
    lock <- file.path(parent, "controller.lock")
    token <- paste(Sys.getpid(), format(Sys.time(), "%Y%m%d%H%M%OS6"), basename(tempfile()), sep = "_")
    acquire_lock(lock, token)
    launched <- FALSE
    process <- NULL
    on.exit(if (!launched) {
      if (!is.null(process) && isTRUE(tryCatch(process$is_alive(), error = function(e) FALSE))) try(process$kill_tree(), silent = TRUE)
      release_lock(lock, token)
    }, add = TRUE)
    output <- file.path(parent, paste0("queue_", signature))
    dir.create(output, recursive = TRUE, showWarnings = FALSE)
    attempt <- tempfile(paste0("attempt_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_"), output)
    need(dir.create(attempt), "Cannot create queue attempt directory.")
    cfg <- list(version = VERSION, project_root = root, tasks = audit_table,
      output_directory = output, attempt_directory = attempt, lock = lock, token = token,
      progress_file = file.path(attempt, "queue_progress.rds"), handoff_file = file.path(attempt, "controller_handoff.rds"), stop_on_error = stop_on_error,
      model_arguments = arguments, poll_seconds = poll_seconds)
    save_atomic(cfg, file.path(attempt, "queue_config.rds"))
    save_atomic(list(stage = "STARTING_CONTROLLER", terminal = FALSE,
      updated_at = as.character(Sys.time()), current_task = NA_character_), cfg$progress_file)
    log <- file.path(attempt, "controller.log")
    process <- callr::r_bg(controller, args = list(cfg = cfg), libpath = .libPaths(),
      stdout = log, stderr = "2>&1", supervise = TRUE, user_profile = FALSE,
      system_profile = FALSE, package = TRUE, wd = root)
    owner <- list(token = token, pid = process$get_pid(), created = tryCatch(
      as.numeric(ps::ps_create_time(ps::ps_handle(process$get_pid()))), error = function(e) NA_real_),
      child_pid = NULL, child_created = NULL)
    save_atomic(owner, file.path(lock, "owner.rds"))
    save_atomic(TRUE, cfg$handoff_file)
    launched <- TRUE
    job <- list(process = process, config = cfg, progress_file = cfg$progress_file,
      log_file = log, output_directory = output, attempt_directory = attempt)
    options(pheno.mundlak.phylogeny.weekend.job = job)
    message("Cua de cinc models iniciada; un worker de model alhora. Manteniu R obert i l'ordinador despert.",
      "\nSeguiment: watch_mundlak_phylogeny_weekend(job)", "\nRegistre: ", log)
    invisible(job)
  }
  status <- function(job = getOption("pheno.mundlak.phylogeny.weekend.job"), quiet = FALSE) {
    need(!is.null(job), "First start_mundlak_phylogeny_weekend(check_only = FALSE).")
    p <- tryCatch(readRDS(job$progress_file), error = function(e) NULL)
    alive <- tryCatch(job$process$is_alive(), error = function(e) FALSE)
    if (!quiet) {
      if (!is.null(p)) print(p)
      cat("Controller alive:", alive, "\nLog:", job$log_file, "\n")
    }
    invisible(list(progress = p, alive = alive, log_file = job$log_file))
  }
  watch <- function(job = getOption("pheno.mundlak.phylogeny.weekend.job"), interval_seconds = 60) {
    need(!is.null(job), "First start_mundlak_phylogeny_weekend(check_only = FALSE).")
    need(is.numeric(interval_seconds) && length(interval_seconds) == 1L && is.finite(interval_seconds) && interval_seconds >= 1,
      "Invalid monitoring interval.")
    tryCatch({
      repeat {
        status(job)
        if (!job$process$is_alive()) break
        until <- Sys.time() + interval_seconds
        while (Sys.time() < until && job$process$is_alive()) Sys.sleep(1)
      }
      tryCatch(invisible(job$process$get_result()), error = function(e) {
        message("Controller aturat; checkpoints conservats. Consulteu ", job$log_file, "\n", conditionMessage(e))
        invisible(NULL)
      })
    }, interrupt = function(e) {
      message("Monitor pausat; la cua continua. Torneu a executar watch_mundlak_phylogeny_weekend(job).")
      invisible(job)
    })
  }
  stop_queue <- function(job = getOption("pheno.mundlak.phylogeny.weekend.job")) {
    if (is.null(job)) return(invisible(FALSE))
    if (job$process$is_alive()) {
      job$process$kill_tree()
      job$process$wait(timeout = 10000)
    }
    need(!job$process$is_alive(), "The queue controller could not be stopped.")
    old <- tryCatch(readRDS(file.path(job$config$lock, "owner.rds")), error = function(e) NULL)
    child <- if (is.null(old) || is.null(old$child_pid)) FALSE else process_matches(old$child_pid, old$child_created)
    until <- Sys.time() + 5
    while (!identical(child, FALSE) && Sys.time() < until) {
      Sys.sleep(0.1)
      child <- process_matches(old$child_pid, old$child_created)
    }
    if (identical(child, FALSE)) release_lock(job$config$lock, job$config$token)
    p <- tryCatch(readRDS(job$progress_file), error = function(e) list())
    p$stage <- "STOPPED_BY_USER"; p$terminal <- TRUE; p$updated_at <- as.character(Sys.time())
    save_atomic(p, job$progress_file)
    message("Cua aturada; checkpoints conservats. Repetiu start amb les mateixes opcions per reprendre.")
    invisible(TRUE)
  }
  list(audit = audit, start = start, status = status, watch = watch, stop = stop_queue)
})

audit_mundlak_phylogeny_weekend <- function(...) pheno_mundlak_phylogeny_weekend$audit(...)
start_mundlak_phylogeny_weekend <- function(...) pheno_mundlak_phylogeny_weekend$start(...)
resume_mundlak_phylogeny_weekend <- function(...) pheno_mundlak_phylogeny_weekend$start(...)
status_mundlak_phylogeny_weekend <- function(...) pheno_mundlak_phylogeny_weekend$status(...)
watch_mundlak_phylogeny_weekend <- function(...) pheno_mundlak_phylogeny_weekend$watch(...)
stop_mundlak_phylogeny_weekend <- function(...) pheno_mundlak_phylogeny_weekend$stop(...)

message("Cua carregada; cap model iniciat. Primer: audit_mundlak_phylogeny_weekend(). La cua completa requereix check_only = FALSE.")
