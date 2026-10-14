resume_phylogeny_qc_15km <-
function(
    output_dir = paste0(
      "E:/phenoIMPACT project/code/phenoIMPACT/output/",
      "phenology_plasticity/phylogeny_qc_15km_v2/",
      "offset_multivoltine_b0dfeb3c6754"),
    project_root = "E:/phenoIMPACT project/code/phenoIMPACT",
    qc_run = "run_fbf8b3101079", checkpoint_seconds = 300,
    minimum_free_GiB = 30) {
  
  need <- function(ok, msg)
    if (!isTRUE(ok)) stop(msg, call. = FALSE)
  
  need(exists("run_phylogeny_qc_15km_v2", mode = "function"),
       "Carrega primer el script v2 amb source().")
  
  for (option in c("pheno.phylo.qc.job", "pheno.phylo.pilot.job")) {
    old <- getOption(option)
    need(is.null(old) || !isTRUE(tryCatch(old$process$is_alive(),
                                          error = function(e) FALSE)), "Ja hi ha un worker en marxa.")
  }
  
  root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
  out <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
  expected_parent <- normalizePath(file.path(root,
                                             "output/phenology_plasticity/phylogeny_qc_15km_v2"),
                                   winslash = "/", mustWork = TRUE)
  
  need(identical(dirname(out), expected_parent),
       "Carpeta de resultats inesperada.")
  
  state_path <- file.path(out, "portable_optimizer_state.rds")
  need(file.exists(state_path), "Falta el checkpoint portable.")
  state <- readRDS(state_path)
  
  need(
    is.character(state$signature) && length(state$signature) == 1L &&
      !is.na(state$signature) &&
      grepl("^[a-f0-9]{64}$", state$signature) &&
      identical(basename(out),
                paste0("offset_multivoltine_", substr(state$signature, 1L, 12L))) &&
      is.numeric(state$theta) && length(state$theta) > 0L &&
      !is.null(names(state$theta)) && all(is.finite(state$theta)) &&
      length(state$nll) == 1L && is.finite(state$nll) &&
      length(state$attempt) == 1L && is.finite(state$attempt) &&
      state$attempt >= 0 && state$attempt %% 1 == 0 &&
      length(state$stage) == 1L &&
      state$stage %in% c("NUMERICAL_REFINEMENT_REQUIRED",
                         "FINITE_OPTIMIZER_EVALUATION", "OPTIMIZED_BEFORE_SDREPORT"),
    "Checkpoint inconsistent o etapa no reconeguda.")
  
  ck <- utils::read.csv(file.path(out, "fit_checks.csv"))
  need(nrow(ck) == 1L && isTRUE(ck$convergence_code == 0L) &&
         isTRUE(ck$positive_definite_Hessian) && isTRUE(ck$finite_fixed_SE) &&
         is.finite(ck$max_abs_gradient) && ck$max_abs_gradient >= .001 &&
         !isTRUE(ck$fit_ok),
       "Aquesta represa es per un ajust que nomes falla el gradient.")
  
  used <- utils::read.csv(file.path(out, "source_cache_used.csv"))
  need(nrow(used) == 1L && file.exists(used$cache_file),
       "Falta la cache portable del model font.")
  
  runner <- get("run_phylogeny_qc_15km_v2", mode = "function")
  top <- as.list(body(runner))
  idx <- which(vapply(top, function(x) is.call(x) &&
                        identical(x[[1L]], as.name("for")), logical(1)))
  need(length(idx) == 1L, "Estructura del runner v2 no reconeguda.")
  
  backup <- tempfile("before_gradient_refinement_", tmpdir = out)
  dir.create(backup)
  files <- file.path(out, c(
    "fit_phylogenetic.rds", "portable_optimizer_state.rds",
    "fit_checks.csv", "inference_status.txt",
    "warnings.txt", "model_comparison.csv"))
  files <- files[file.exists(files)]
  need(all(file.copy(files, backup, overwrite = FALSE)),
       "No es pot preservar l'ajust anterior.")
  
  patch <- substitute({
    joint <- e$build_joint
    parts <- as.list(body(joint))
    parts <- append(parts, list(quote({
      fit$nlminb_control <- utils::modifyList(fit$nlminb_control,
                                              list(rel.tol = 1e-12, x.tol = 1e-10,
                                                   iter.max = 1000L, eval.max = 3000L, trace = 0L))
      fit$control <- utils::modifyList(fit$control,
                                       list(rel.tol = 1e-12, x.tol = 1e-10,
                                            iter.max = 1000L, eval.max = 3000L, trace = 0L))
    })), after = length(parts) - 1L)
    body(joint) <- as.call(parts)
    e$build_joint <- joint
    
    w <- e$worker
    parts <- as.list(body(w))
    i <- which(vapply(parts, function(x) is.call(x) &&
                        identical(x[[1L]], as.name("<-")) &&
                        identical(x[[2L]], as.name("p")), logical(1)))
    if (length(i) != 1L) stop("Worker v2 no reconegut.")
    
    guard <- substitute({
      need(identical(
        normalizePath(p$out, winslash = "/", mustWork = TRUE), TARGET) &&
          identical(p$signature, SIGNATURE),
        "Les entrades o versions han canviat; no es repren un altre model.")
      write_csv(data.frame(
        rel_tol = 1e-12, x_tol = 1e-10,
        iter_max = 1000L, eval_max = 3000L,
        gradient_threshold = cfg$gradient_threshold,
        maximum_outer_attempts_this_launch = 1L,
        original_backup = BACKUP),
        file.path(p$out, "gradient_refinement_controls.csv"))
    }, list(TARGET = EXPECTED, SIGNATURE = SIG, BACKUP = SAVED))
    
    parts <- append(parts, list(guard), after = i)
    body(w) <- as.call(parts)
    e$worker <- w
    
    sync <- e$synchronize_fit
    parts <- as.list(body(sync))
    audit <- substitute({
      theta <- fit$model$par
      g <- as.numeric(fit$gradients)
      V <- fit$sd_report$cov.fixed
      d <- diag(V)
      s <- sqrt(ifelse(is.finite(d) & d >= 0, d, NA_real_))
      j <- ave(seq_along(theta), names(theta), FUN = seq_along)
      
      write_csv(data.frame(
        parameter = paste0(names(theta), "[", j, "]"),
        gradient = g, absolute_gradient = abs(g),
        SE_times_gradient = s * g),
        file.path(TARGET, "gradient_diagnostics.csv"))
      
      write_csv(data.frame(
        convergence = fit$model$convergence,
        message = paste(fit$model$message, collapse = "; "),
        iterations = fit$model$iterations,
        objective = fit$model$objective),
        file.path(TARGET, "refinement_optimizer_summary.csv"))
    }, list(TARGET = EXPECTED))
    
    parts <- append(parts, list(audit), after = length(parts) - 1L)
    body(sync) <- as.call(parts)
    e$synchronize_fit <- sync
  }, list(EXPECTED = out, SIG = state$signature, SAVED = backup))
  
  top <- append(top, list(patch), after = idx)
  body(runner) <- as.call(top)
  message("Ajust anterior preservat: ", backup)
  
  job <- runner(project_root = root, qc_run = qc_run,
                checkpoint_seconds = checkpoint_seconds,
                extra_optimization_attempts = 0L,
                minimum_free_GiB = minimum_free_GiB)
  
  job$previous_fit_backup <- backup
  invisible(job)
}
