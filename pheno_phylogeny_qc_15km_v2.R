run_phylogeny_qc_15km_v2 <-
function (project_root = "E:/phenoIMPACT project/code/phenoIMPACT", 
    qc_run = "run_fbf8b3101079", check_only = FALSE, checkpoint_seconds = 300, 
    maximum_single_basis_GiB = 2, extra_optimization_attempts = 1L, 
    minimum_free_GiB = 30) 
{
    {
        if (!requireNamespace("ps", quietly = TRUE)) 
            stop("Cal el paquet ps.")
        free <- ps::ps_system_memory()$avail/1024^3
        message("RAM disponible abans d'iniciar: ", round(free, 
            1), " GiB")
        if (free < minimum_free_GiB) 
            stop(paste("RAM disponible insuficient per aquesta prova.", 
                "Usa una sessio R neta abans d'iniciar-la."))
    }
    root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
    stopifnot(length(qc_run) == 1L, grepl("^run_[a-f0-9]{12}$", 
        qc_run), length(checkpoint_seconds) == 1L, is.finite(checkpoint_seconds), 
        checkpoint_seconds >= 0)
    private <- new.env(parent = globalenv())
    suppressMessages(sys.source(file.path(root, "pheno_plasticity_phylogeny_pilot_15km.R"), 
        envir = private))
    engine <- private$pheno_phylogeny_pilot
    e <- environment(engine$start)
    if (!identical(e$version, "phylo_slope_pilot_1.0.6")) {
        stop("Cal el pilot existent v1.0.6.")
    }
    e$version <- "phylo_qc_multivoltine_edges_15km_1.1.0"
    e$need <- function(ok, msg) {
        if (!isTRUE(ok)) 
            stop(msg, call. = FALSE)
    }
    e$same_num <- function(a, b, tol = 9.9999999999999995e-07) {
        length(a) == length(b) && all(is.finite(a)) && all(is.finite(b)) && 
            (length(a) == 0L || max(abs(as.numeric(a) - as.numeric(b))) <= 
                tol)
    }
    fixed_table <- function(fit) {
        terms <- colnames(stats::model.matrix(fit$split_formula[[1]]$form_no_bars, 
            fit$data[seq_len(min(3L, nrow(fit$data))), , drop = FALSE]))
        b <- fit$sd_report$par.fixed
        V <- fit$sd_report$cov.fixed
        i <- par_indices(b, "b_j")
        need(length(i) == length(terms) && !anyDuplicated(terms) && 
            is.matrix(V) && nrow(V) == length(b) && identical(names(fit$model$par), 
            names(b)) && same_num(fit$model$par, b) && all(is.finite(b[i])), 
            "Coeficients fixos o report inconsistents.")
        v <- diag(V)[i]
        se <- sqrt(ifelse(is.finite(v) & v >= 0, v, NA_real_))
        est <- unname(b[i])
        data.frame(term = terms, estimate = est, SE = se, lower_95 = est - 
            1.96 * se, upper_95 = est + 1.96 * se, p_Wald = 2 * 
            stats::pnorm(-abs(est/se)))
    }
    prepare_uncached_source <- function(cfg, packages) {
        progress("Reading completed post-QC multivoltine fit")
        memory_snapshot(cfg$monitor_directory, "before_source_read")
        src <- dirname(cfg$source_checkpoint)
        run <- dirname(src)
        rc <- readRDS(file.path(run, "run_config.rds"))
        done <- readRDS(file.path(src, "completed.rds"))
        saved <- readRDS(cfg$source_checkpoint)
        saved <- phylo_cache_strip_unused_source(saved)
        invisible(gc())
        fit <- saved$fit
        memory_snapshot(cfg$monitor_directory, "after_source_read")
        need(identical(done$status, "FIT_OK") && identical(done$signature, 
            rc$signature) && identical(saved$signature, rc$signature) && 
            identical(saved$analysis, "offset_multivoltine") && 
            identical(basename(run), paste0("run_", substr(rc$signature, 
                1L, 12L))), "Identitat/complecio QC inconsistent.")
        need(inherits(fit, "sdmTMB") && as.character(fit$version) == 
            "1.1.0" && isFALSE(fit$reml) && identical(fit$time, 
            "YEAR") && all(fit$spatial == "off") && all(fit$spatiotemporal == 
            "iid") && identical(fit$family$family, "gaussian") && 
            identical(fit$family$link, "identity") && !isTRUE(as.logical(fit$tmb_data$has_smooths)) && 
            is.null(fit$nonlocal_formula) && is.null(fit$time_varying) && 
            is.null(fit$spatial_varying), "Estructura/version del model QC inesperada.")
        need(identical(fit$priors, sdmTMB::sdmTMBpriors()), "Priors no estandard: cal revisar la comparacio ML/AIC.")
        f0 <- if (inherits(fit$formula, "formula")) {
            fit$formula
        }
        else fit$formula[[1L]]
        environment(f0) <- baseenv()
        sf <- fit$split_formula[[1L]]
        expected <- OFFSET_mean ~ ONSET_mean_z + clim_anomaly_tw90 * 
            (photo_tw90 + clim_background_tw90 + clim_predictability_tw90 + 
                clim_trend_tw90)
        need(identical(all.vars(sf$form_no_bars[[2L]]), "OFFSET_mean") && 
            setequal(attr(stats::terms(sf$form_no_bars), "term.labels"), 
                attr(stats::terms(expected), "term.labels")), 
            "Formula QC inesperada o falta ONSET_mean_z.")
        bars <- sort(vapply(reformulas::findbars(f0), function(z) paste(deparse(z), 
            collapse = ""), character(1)))
        reference <- OFFSET_mean ~ (1 | SITE_ID) + (1 | SPECIES) + 
            (0 + clim_anomaly_tw90 | SPECIES_slope) + (1 | site_year_id)
        need(identical(bars, sort(vapply(reformulas::findbars(reference), 
            function(z) paste(deparse(z), collapse = ""), character(1)))), 
            "Termes aleatoris QC inesperats.")
        d <- fit$data
        cols <- unique(c(all.vars(f0), "source_row", "YEAR", 
            "ONSET_mean", "x_km", "y_km"))
        need(all(cols %in% names(d)) && !anyNA(d[, cols]) && 
            !anyDuplicated(d$source_row) && all(vapply(d[, c("OFFSET_mean", 
            "ONSET_mean", "ONSET_mean_z", "YEAR", anomaly, moderators, 
            "x_km", "y_km")], function(x) is.numeric(x) && all(is.finite(x)), 
            logical(1))) && all(d$YEAR == as.integer(d$YEAR)) && 
            identical(as.character(d$SPECIES), as.character(d$SPECIES_slope)), 
            "Dades/grouping QC invalids.")
        q <- saved$qc
        qcols <- c("n_source_rows", "n_excluded", "n_retained")
        qcsv <- utils::read.csv(file.path(src, "qc_filter_summary.csv"))
        need(nrow(q) == 1L && nrow(qcsv) == 1L && all(qcols %in% 
            names(q)) && same_num(unlist(q[qcols]), unlist(qcsv[qcols])) && 
            q$n_retained == nrow(d) && q$n_excluded > 0 && q$n_source_rows == 
            q$n_retained + q$n_excluded, "Comptatges QC inconsistents.")
        sc <- utils::read.csv(file.path(src, "onset_scaling_after_qc.csv"))
        need(nrow(sc) == 1L && sc$new_SD > 0 && same_num(d$ONSET_mean_z, 
            (d$ONSET_mean - sc$new_centre)/sc$new_SD), "Escala ONSET_mean_z inconsistent.")
        need(any(grepl("Same 15-km triangulation", readLines(file.path(src, 
            "model_specification.txt")), fixed = TRUE)), "No es confirma la malla de 15 km.")
        need((is.null(fit$offset) || all(fit$offset == 0)) && 
            (is.null(fit$tmb_data$weights) || all(fit$tmb_data$weights == 
                1)), "Offsets/weights inesperats.")
        ck <- checks(fit, cfg$gradient_threshold)
        need(isTRUE(ck$fit_ok), "El model QC no supera els controls numerics.")
        ix <- par_indices(fit$sd_report$par.fixed, "b_j")
        values <- fit$sd_report$value[grepl("^re_b_pars($|\\[)", 
            names(fit$sd_report$value))]
        need(same_num(fit$parlist$b_j, fit$sd_report$par.fixed[ix]) && 
            same_num(as.vector(fit$parlist$re_b_pars), values), 
            "Parametres portables QC/report desactualitzats.")
        species <- levels(droplevels(factor(d$SPECIES)))
        tree <- ape::read.tree(cfg$tree_file)
        need(inherits(tree, "phylo") && !inherits(tree, "multiPhylo") && 
            ape::is.rooted(tree) && ape::is.ultrametric(tree) && 
            length(tree$edge.length) == nrow(tree$edge) && all(is.finite(tree$edge.length) & 
            tree$edge.length > 0), "Arbre invalid.")
        sp <- normalize_species(species)
        tips <- normalize_species(tree$tip.label)
        need(!anyDuplicated(sp) && !anyDuplicated(tips) && !anyNA(match(sp, 
            tips)), "Matching especie/arbre incomplet.")
        tip_names <- tree$tip.label[match(sp, tips)]
        C <- ape::vcv.phylo(tree)[tip_names, tip_names, drop = FALSE]
        need(diff(range(diag(C)))/mean(diag(C)) < 9.9999999999999995e-08, 
            "Covariancia no ultrametrica.")
        A <- C/mean(diag(C))
        dimnames(A) <- list(species, species)
        chol(A)
        {
            edge <- build_phylo_edge_basis(tree, tip_names)
            need(max(abs(unname(edge$A) - unname(A))) < 1e-10, 
                "Covariancia d'increments incorrecta.")
            attr(A, "edge_basis") <- unname(edge$B)
            attr(A, "edge_audit") <- edge$audit
        }
        d$PHYLO_ID <- match(as.character(d$SPECIES), species)
        mesh <- fit$spde
        need(nrow(mesh$loc_xy) == nrow(d) && same_num(as.vector(as.matrix(mesh$loc_xy)), 
            as.vector(as.matrix(d[, c("x_km", "y_km")]))), "Projeccio de la malla/ordre de files inconsistent.")
        need(same_num(as.numeric(fit$response), d$OFFSET_mean), 
            "Resposta QC inconsistent.")
        basis_gib <- 8 * nrow(d) * ncol(attr(A, "edge_basis"))/1024^3
        need(basis_gib <= cfg$maximum_single_basis_GiB, paste("La base afegida excedeix el limit:", 
            basis_gib, "GiB."))
        start <- portable_source_start(fit, cfg$monitor_directory)
        replay_files <- list.files(cfg$monitor_directory, pattern = "^source_replay_[12]_evaluation[.]csv$", 
            full.names = TRUE)
        iid <- exp(fit$sd_report$par.fixed[slope_logsd_index(fit)])
        start$ln_smooth_sigma <- matrix(log(max(0.10000000000000001 * 
            iid, 0.01)), 1L, 1L)
        source <- list(fixed = fixed_table(fit), stats = fit_stats(fit), 
            checks = ck, saved_checks = utils::read.csv(file.path(src, 
                "fit_checks.csv")), row_hash = digest::digest(d$source_row, 
                algo = "sha256"), run_id = saved$signature, fingerprint = digest::digest(list(saved$signature, 
                fit$model$par, fit$model$objective, q), algo = "sha256"), 
            checkpoint = cfg$source_checkpoint, formula = f0)
        signature <- digest::digest(list(version, source$fingerprint, 
            d[, cols], A, mesh$mesh, packages, cfg$gradient_threshold, 
            cfg$fit_seed), algo = "sha256")
        out <- file.path(cfg$output_root, paste0("offset_multivoltine_", 
            substr(signature, 1L, 12L)))
        dir.create(out, recursive = TRUE, showWarnings = FALSE)
        {
            write_csv(attr(A, "edge_audit"), file.path(out, "tree_increment_basis_checks.csv"))
            writeLines(c("Original-root Brownian covariance represented by ancestral increments.", 
                "The scientific model and original 15-km mesh are preserved.", 
                "One additional outer variance parameter for Gaussian ML/AIC comparison.", 
                "Dense input storage remains; sufficient RAM is not guaranteed.", 
                "Fixed-parameter covariance retained; ADREPORT SEs skipped.", 
                "Saved-parameter/report checks replace native source replay."), 
                file.path(out, "RECOVERY_README.txt"))
        }
        write_csv(packages, file.path(out, "package_versions.csv"))
        write_csv(source$checks, file.path(out, "source_fit_checks.csv"))
        write_csv(source$saved_checks, file.path(out, "source_saved_checks.csv"))
        write_csv(q, file.path(out, "source_qc_filter_summary.csv"))
        write_csv(data.frame(SPECIES = species, tip_label = tip_names, 
            PHYLO_ID = seq_along(species)), file.path(out, "species_mapping.csv"))
        write_csv(data.frame(n_observations = nrow(d), n_species = length(species), 
            n_sites = length(unique(d$SITE_ID)), n_years = length(unique(d$YEAR)), 
            mesh_reused_exactly = TRUE, single_added_basis_GiB = basis_gib, 
            note = "One dense basis only; not a total RAM estimate"), 
            file.path(out, "input_summary.csv"))
        save_atomic(list(A = A, tree = tree, source = source, 
            signature = signature), file.path(out, "provenance.rds"))
        f <- add_phylo(f0, A, anomaly)
        writeLines(c("Post-QC multivoltines; same rows, scales, onset, random terms and exact 15-km mesh.", 
            "Annual spatial fields preserved: spatial=off, spatiotemporal=iid.", 
            "Cov(species slope) = sigma_iid^2 I + sigma_phylo^2 A.", 
            "Original tree root; diagonal(A)=1. Species intercepts remain independent.", 
            "Gaussian ML. Approximate variance intervals; no automatic significance test.", 
            "Periodic optimizer snapshots use portable outer parameters.", 
            "Stopping preserves completed snapshots; restart creates fresh native objects.", 
            paste("Source:", cfg$source_checkpoint), paste("Formula:", 
                paste(deparse(f0), collapse = " ")), paste("One added basis GiB:", 
                basis_gib, "(not total RAM; several copies and spatial Hessians need more).")), 
            file.path(out, "model_specification.txt"))
        priors <- fit$priors
        extra_time <- fit$extra_time
        rm(fit, saved, C)
        invisible(gc())
        memory_snapshot(cfg$monitor_directory, "source_fit_released")
        list(data = d, A = A, formula = f, mesh = mesh, start = start, 
            source = source, out = out, signature = signature, 
            priors = priors, extra_time = extra_time)
    }
    synchronize_fit <- function(fit) {
        theta <- fit$model$par
        obj <- fit$tmb_obj
        nll <- obj$fn(theta)
        need(is.finite(nll), "Versemblanca final no finita.")
        fit$model$objective <- nll
        fit$sd_report <- TMB::sdreport(obj, par.fixed = theta, 
            getJointPrecision = FALSE, getReportCovariance = FALSE, 
            skip.delta.method = TRUE)
        fit$gradients <- as.numeric(obj$gr(theta))
        obj$fn(theta)
        fit$parlist <- obj$env$parList(par = obj$env$last.par)
        fit$last.par.best <- obj$env$last.par
        need(identical(names(theta), names(fit$sd_report$par.fixed)) && 
            same_num(theta, fit$sd_report$par.fixed), "Report final no correspon als parametres de l'optimitzador.")
        values <- fit$sd_report$value[grepl("^re_b_pars($|\\[)", 
            names(fit$sd_report$value))]
        need(same_num(as.vector(fit$parlist$re_b_pars), values), 
            "Modes aleatoris i report final inconsistents.")
        fit
    }
    preflight_onset <- function(directory) {
        A <- kronecker(diag(3L), matrix(c(1, 0.65000000000000002, 
            0.65000000000000002, 1), 2L))
        d <- data.frame(PHYLO_ID = rep(1:6, each = 8), x = rep(seq(-1, 
            1, length.out = 8), 6))
        d$g <- factor(d$PHYLO_ID)
        d$ONSET_mean_z <- cos(seq_len(nrow(d)) * 0.37)
        d$y <- sin(seq_len(nrow(d)) * 0.69999999999999996)
        f <- add_phylo(y ~ 1 + x + ONSET_mean_z + (0 + x | g), 
            A, "x")
        fit <- sdmTMB::sdmTMB(f, data = d, family = stats::gaussian(), 
            spatial = "off", spatiotemporal = "off", reml = FALSE, 
            do_fit = FALSE, silent = TRUE, control = sdmTMB::sdmTMBcontrol(get_joint_precision = FALSE, 
                newton_loops = 0L))
        on.exit(release_tmb(fit$tmb_obj), add = TRUE)
        gaussian_inner(fit$tmb_obj)
        X <- stats::model.matrix(~1 + x + ONSET_mean_z, d)
        results <- lapply(1:2, function(k) {
            b <- if (k == 1L) 
                c(0, 0, 0)
            else c(0.69999999999999996, -0.29999999999999999, 
                0.45000000000000001)
            s <- if (k == 1L) 
                c(1, 1, 1)
            else c(1.2, 0.5, 0.90000000000000002)
            theta <- fit$tmb_obj$par
            pars <- list(b_j = b, ln_phi = log(s[1]), re_cov_pars = log(s[2]), 
                ln_smooth_sigma = log(s[3]))
            for (nm in names(pars)) {
                j <- par_indices(theta, nm)
                need(length(j) == length(pars[[nm]]), "Mapping del test onset inesperat.")
                theta[j] <- pars[[nm]]
            }
            actual <- fit$tmb_obj$fn(theta)
            g <- fit$tmb_obj$gr(theta)
            id <- d$PHYLO_ID
            V <- s[1]^2 * diag(nrow(d)) + outer(d$x, d$x) * (s[3]^2 * 
                A[id, id] + s[2]^2 * outer(id, id, "=="))
            R <- chol(V)
            r <- as.numeric(d$y - X %*% b)
            z <- forwardsolve(t(R), r)
            exact <- 0.5 * (nrow(d) * log(2 * pi) + 2 * sum(log(diag(R))) + 
                sum(z^2))
            ge <- -as.numeric(crossprod(X, solve(V, r)))
            data.frame(setting = k, nll_error = abs(actual - 
                exact), fixed_gradient_error = max(abs(g[par_indices(theta, 
                "b_j")] - ge)), gradient_finite = all(is.finite(g)))
        })
        results <- do.call(rbind, results)
        write_csv(results, file.path(directory, "preflight_onset_likelihood.csv"))
        need(all(results$nll_error < 9.9999999999999995e-07 & 
            results$fixed_gradient_error < 9.9999999999999995e-07 & 
            results$gradient_finite), "Test natiu amb onset incorrecte; no s'ha carregat el model QC.")
    }
    worker <- function(cfg) {
        monitor$path <- cfg$progress_file
        monitor$state <- readRDS(monitor$path)
        finished <- FALSE
        on.exit({
            monitor$state$finished <- TRUE
            monitor$state$failed <- !finished
            monitor$state$ended_at <- as.numeric(Sys.time())
            if (!finished) {
                monitor$state$phase <- "Stopped with error; inspect worker.log and optimizer snapshot"
            }
            save_atomic(monitor$state, monitor$path)
        }, add = TRUE)
        packages <- dependency_check()
        configure_memory(cfg$monitor_directory)
        progress("Checking native phylogenetic likelihood with onset")
        preflight_onset(cfg$monitor_directory)
        writeLines(capture.output(utils::sessionInfo()), file.path(cfg$monitor_directory, 
            "sessionInfo.txt"))
        p <- prepare(cfg, packages)
        monitor$state$result_directory <- p$out
        file.copy(file.path(cfg$monitor_directory, "preflight_onset_likelihood.csv"), 
            p$out, overwrite = TRUE)
        audit_files <- list.files(cfg$monitor_directory, pattern = "^(source_replay_.*[.]csv|tmb_memory_settings[.]csv|memory_checkpoints[.]txt)$", 
            full.names = TRUE)
        if (length(audit_files)) {
            file.copy(audit_files, p$out, overwrite = TRUE)
        }
        progress("Checking full species covariance and exact Gaussian likelihood")
        preflight(p$A, p$out)
        if (cfg$check_only) {
            progress("Checks passed; no full fit attempted")
            finished <- TRUE
            return(list(output_dir = p$out, check_only = TRUE))
        }
        path <- file.path(p$out, "fit_phylogenetic.rds")
        cached <- if (file.exists(path)) 
            readRDS(path)
        if (!is.null(cached)) {
            need(identical(cached$signature, p$signature), "Checkpoint del model inconsistent.")
            if (isTRUE(checks(cached$fit, cfg$gradient_threshold)$fit_ok) && 
                cached$fit$model$objective <= -p$source$stats$logLik + 
                  0.00050000000000000001) {
                okay <- export_results(cached, p, cfg)
                finished <- TRUE
                progress("Completed from matching checkpoint")
                return(list(output_dir = p$out, cached = TRUE, 
                  inference_ready = okay))
            }
        }
        rm(cached)
        invisible(gc())
        set.seed(cfg$fit_seed)
        warnings <- character()
        collect <- function(w) {
            warnings <<- unique(c(warnings, conditionMessage(w)))
            writeLines(warnings, file.path(p$out, "warnings.txt"))
            message("Warning: ", conditionMessage(w))
            invokeRestart("muffleWarning")
        }
        fit <- withCallingHandlers(initialize_joint(p), warning = collect)
        memory_snapshot(p$out, "joint_initialized")
        initial_names <- colnames(stats::model.matrix(fit$split_formula[[1]]$form_no_bars, 
            fit$data[seq_len(min(3L, nrow(fit$data))), , drop = FALSE]))
        need(identical(initial_names, p$source$fixed$term), "L'ordre dels coeficients fixos ha canviat.")
        state_path <- file.path(p$out, "portable_optimizer_state.rds")
        state <- if (file.exists(state_path)) 
            readRDS(state_path)
        attempt <- 0L
        optimized <- FALSE
        if (!is.null(state)) {
            need(identical(state$signature, p$signature) && identical(names(state$theta), 
                names(fit$tmb_obj$par)) && length(state$theta) == 
                length(fit$tmb_obj$par) && all(is.finite(state$theta)), 
                "Snapshot de l'optimitzador inconsistent.")
            replay <- evaluate_start(fit$tmb_obj, state$theta, 
                p$out, "joint_restart", state$nll)
            need(isTRUE(replay$finite) && abs(replay$value - 
                state$nll) < 0.01, "No es pot reproduir el snapshot guardat.")
            fit$model$par <- state$theta
            fit$model$objective <- replay$value
            attempt <- state$attempt
            if (identical(state$stage, "OPTIMIZED_BEFORE_SDREPORT")) {
                fit$model <- state$optimizer
                optimized <- TRUE
            }
            progress("Restarting from saved outer parameters; fresh native model")
        }
        started <- proc.time()[["elapsed"]]
        monitor$state$fit_started_at <- as.numeric(Sys.time())
        for (offset in seq.int(0L, cfg$extra_optimization_attempts)) {
            k <- attempt + offset
            obj <- fit$tmb_obj
            original_fn <- obj$fn
            tracker <- new.env(parent = emptyenv())
            tracker$theta <- fit$model$par
            tracker$nll <- original_fn(tracker$theta)
            tracker$last_write <- -Inf
            tracker$evaluations <- 0L
            persist <- function(stage = "FINITE_OPTIMIZER_EVALUATION", 
                optimizer = NULL) {
                save_atomic(list(signature = p$signature, theta = tracker$theta, 
                  nll = tracker$nll, stage = stage, attempt = k, 
                  optimizer = optimizer, time = Sys.time(), evaluations_this_attempt = tracker$evaluations), 
                  state_path)
                tracker$last_write <- as.numeric(Sys.time())
            }
            obj$fn <- function(x, ...) {
                if (missing(x)) {
                  value <- original_fn(...)
                  x <- obj$env$last.par[obj$env$lfixed()]
                }
                else value <- original_fn(x, ...)
                tracker$evaluations <- tracker$evaluations + 
                  1L
                if (length(value) == 1L && is.finite(value) && 
                  all(is.finite(x))) {
                  if (value < tracker$nll) {
                    tracker$theta <- x
                    tracker$nll <- value
                  }
                  if (as.numeric(Sys.time()) - tracker$last_write >= 
                    cfg$checkpoint_seconds) {
                    persist()
                    progress(paste("Optimizer attempt", k, "evaluations", 
                      tracker$evaluations, "snapshot saved"))
                  }
                }
                value
            }
            if (!optimized) {
                progress(paste("Fitting 15-km joint model; optimizer attempt", 
                  k))
                fit$model <- withCallingHandlers(stats::nlminb(start = fit$model$par, 
                  objective = obj$fn, gradient = obj$gr, control = fit$nlminb_control, 
                  lower = fit$lower, upper = fit$upper), warning = collect)
                tracker$theta <- fit$model$par
                tracker$nll <- original_fn(fit$model$par)
                fit$model$objective <- tracker$nll
                persist("OPTIMIZED_BEFORE_SDREPORT", fit$model)
            }
            obj$fn <- original_fn
            fit$tmb_obj <- obj
            progress("Optimizer saved; computing uncertainty at final parameters")
            memory_snapshot(p$out, "before_sdreport")
            fit <- withCallingHandlers(synchronize_fit(fit), 
                warning = collect)
            memory_snapshot(p$out, "after_sdreport")
            saved <- list(signature = p$signature, fit = fit, 
                checks = checks(fit, cfg$gradient_threshold), 
                source = p$source, A = p$A, code_version = version, 
                warnings = warnings, elapsed_minutes_this_worker = (proc.time()[["elapsed"]] - 
                  started)/60)
            save_atomic(saved, path)
            write_csv(saved$checks, file.path(p$out, "fit_checks.csv"))
            if (isTRUE(saved$checks$fit_ok) && fit$model$objective <= 
                -p$source$stats$logLik + 0.00050000000000000001) 
                break
            optimized <- FALSE
            save_atomic(list(signature = p$signature, theta = fit$model$par, 
                nll = fit$model$objective, stage = "NUMERICAL_REFINEMENT_REQUIRED", 
                attempt = k + 1L), state_path)
        }
        need(identical(fit$data$source_row, p$data$source_row), 
            "L'ordre de les files ha canviat.")
        okay <- export_results(saved, p, cfg)
        progress(if (okay) {
            "Completed; compare coefficients and variance uncertainty"
        }
        else {
            "Fit saved; inference not ready, inspect fit_checks/inference_status"
        })
        finished <- TRUE
        list(output_dir = p$out, fit_ok = saved$checks$fit_ok, 
            inference_ready = okay)
    }
    replace_phylo_basis <- function(expr) {
        walk <- function(x) {
            if (is.call(x) && identical(x[[1L]], as.name("<-")) && 
                identical(x[[2L]], as.name("B"))) 
                return(quote(B <- if (is.null(attr(A, "edge_basis"))) t(chol(A)) else attr(A, 
                  "edge_basis")))
            if (is.call(x)) 
                for (i in seq_along(x)) if (!identical(x[[i]], 
                  quote(expr = ))) 
                  x[[i]] <- walk(x[[i]])
            x
        }
        walk(expr)
    }
    replace_phylo_preflight <- function(expr) {
        walk <- function(x) {
            if (identical(x, quote(ncol(sm$Zs[[1]]) != n))) 
                return(quote(ncol(sm$Zs[[1]]) != if (is.null(attr(A, 
                  "edge_basis"))) n else ncol(attr(A, "edge_basis"))))
            if (is.call(x) && identical(x[[1L]], as.name("<-")) && 
                identical(x[[2L]], as.name("At"))) 
                return(quote(At <- local({
                  ids <- seq_len(min(6L, nrow(A)))
                  out <- A[ids, ids, drop = FALSE]
                  attr(out, "edge_basis") <- attr(A, "edge_basis")[ids, 
                    , drop = FALSE]
                  out
                })))
            if (is.call(x)) 
                for (i in seq_along(x)) if (!identical(x[[i]], 
                  quote(expr = ))) 
                  x[[i]] <- walk(x[[i]])
            x
        }
        walk(expr)
    }
    replace_phylo_call_argument <- function(expr, target, name, 
        value) {
        walk <- function(x) {
            if (is.call(x) && identical(x[[1L]], target)) {
                args <- as.list(x)
                args[[name]] <- value
                x <- as.call(args)
            }
            if (is.call(x)) 
                for (i in seq_along(x)) if (!identical(x[[i]], 
                  quote(expr = ))) 
                  x[[i]] <- walk(x[[i]])
            x
        }
        walk(expr)
    }
    build_phylo_edge_basis <- function(tree, tip_names, tolerance = 1e-10) {
        need <- function(ok, message) if (!isTRUE(ok)) 
            stop(message, call. = FALSE)
        need(requireNamespace("ape", quietly = TRUE), "Cal el paquet ape.")
        need(inherits(tree, "phylo") && !inherits(tree, "multiPhylo") && 
            ape::is.rooted(tree) && ape::is.ultrametric(tree), 
            "Cal un arbre arrelat i ultrametric.")
        need(is.character(tree$tip.label) && !anyNA(tree$tip.label) && 
            all(nzchar(tree$tip.label)) && !anyDuplicated(tree$tip.label), 
            "Noms de l'arbre invalids.")
        need(is.character(tip_names) && length(tip_names) >= 
            2L && !anyNA(tip_names) && all(nzchar(tip_names)) && 
            !anyDuplicated(tip_names) && all(tip_names %in% tree$tip.label), 
            "Matching d'especies invalid.")
        E <- tree$edge
        lengths <- tree$edge.length
        need(is.matrix(E) && ncol(E) == 2L && nrow(E) > 0L && 
            length(lengths) == nrow(E) && all(is.finite(lengths)) && 
            all(lengths >= 0) && !anyDuplicated(E[, 2L]), "Branques invalides.")
        need(length(tolerance) == 1L && is.finite(tolerance) && 
            tolerance > 0, "Tolerancia invalida.")
        roots <- setdiff(E[, 1L], E[, 2L])
        need(length(roots) == 1L, "Cal una sola arrel original.")
        root <- roots[[1L]]
        child_edge <- stats::setNames(seq_len(nrow(E)), as.character(E[, 
            2L]))
        descendants <- matrix(FALSE, nrow(E), length(tip_names))
        for (i in seq_along(tip_names)) {
            node <- match(tip_names[[i]], tree$tip.label)
            steps <- 0L
            while (node != root) {
                j <- unname(child_edge[as.character(node)])
                need(length(j) == 1L && !is.na(j), "El cami ancestral no arriba a l'arrel.")
                descendants[j, i] <- TRUE
                node <- E[j, 1L]
                steps <- steps + 1L
                need(steps <= nrow(E), "Cicle en el cami ancestral.")
            }
        }
        active <- which(rowSums(descendants) > 0L & lengths > 
            0)
        need(length(active) > 0L, "No hi ha increments ancestrals.")
        keys <- vapply(active, function(j) paste(which(descendants[j, 
            ]), collapse = ","), character(1))
        key_order <- sort(unique(keys))
        groups <- lapply(key_order, function(k) active[keys == 
            k])
        full_C <- ape::vcv.phylo(tree, corr = FALSE)
        C <- full_C[tip_names, tip_names, drop = FALSE]
        H <- mean(diag(C))
        need(all(is.finite(C)) && is.finite(H) && H > 0 && diff(range(diag(C)))/H < 
            9.9999999999999995e-08, "Profunditats de l'arbre inconsistents.")
        A <- C/H
        clade_ids <- sprintf("clade_%04d", seq_along(groups))
        B <- matrix(0, length(tip_names), length(groups), dimnames = list(tip_names, 
            clade_ids))
        clades <- vector("list", length(groups))
        edge_clade <- rep(NA_character_, nrow(E))
        for (j in seq_along(groups)) {
            rows <- groups[[j]]
            mask <- descendants[rows[[1L]], ]
            weight <- sum(lengths[rows])
            B[mask, j] <- sqrt(weight/H)
            edge_clade[rows] <- clade_ids[[j]]
            clades[[j]] <- data.frame(clade_id = clade_ids[[j]], 
                descendant_key = key_order[[j]], n_matched_descendants = sum(mask), 
                summed_branch_length = weight, normalized_increment_variance = weight/H, 
                n_original_edges = length(rows), original_edge_rows = paste(rows, 
                  collapse = ","), stringsAsFactors = FALSE)
        }
        error <- max(abs(tcrossprod(B) - A))
        need(is.finite(error) && error <= tolerance, paste("Covariancia d'increments incorrecta:", 
            error))
        need(!inherits(try(chol(A), silent = TRUE), "try-error"), 
            "La covariancia no es definida positiva.")
        dense <- t(chol(A))
        list(B = B, A = A, H = H, tip_order = tip_names, clades = do.call(rbind, 
            clades), original_edges = data.frame(edge_row = seq_len(nrow(E)), 
            parent_node = E[, 1L], child_node = E[, 2L], branch_length = lengths, 
            n_matched_descendants = rowSums(descendants), clade_id = edge_clade), 
            audit = data.frame(n_original_tips = length(tree$tip.label), 
                n_selected_tips = nrow(B), n_original_edges = nrow(E), 
                n_active_original_edges = length(active), n_increment_columns = ncol(B), 
                n_collapsed_edges = length(active) - ncol(B), 
                original_root = root, common_original_tip_depth = H, 
                maximum_covariance_error = error, basis_rank = qr(B)$rank, 
                edge_nonzeros = sum(B != 0), edge_density = mean(B != 
                  0), dense_Cholesky_nonzeros = sum(dense != 
                  0), dense_Cholesky_density = mean(dense != 
                  0), mean_edge_nonzeros_per_tip = mean(rowSums(B != 
                  0)), mean_Cholesky_nonzeros_per_tip = mean(rowSums(dense != 
                  0))))
    }
    phylo_cache_need <- function(ok, message) {
        if (!isTRUE(ok)) 
            stop(message, call. = FALSE)
    }
    phylo_cache_same_numeric <- function(x, y, tolerance = 9.9999999999999995e-07) {
        is.numeric(x) && is.numeric(y) && length(x) == length(y) && 
            all(is.finite(x)) && all(is.finite(y)) && (length(x) == 
            0L || max(abs(as.numeric(x) - as.numeric(y))) <= 
            tolerance)
    }
    phylo_cache_parameter_shapes <- function(parameters) {
        phylo_cache_need(is.list(parameters) && length(parameters) > 
            0L && !is.null(names(parameters)) && !anyNA(names(parameters)) && 
            all(nzchar(names(parameters))) && !anyDuplicated(names(parameters)), 
            "Blocs de parametres invalids.")
        data.frame(block = names(parameters), storage_type = vapply(parameters, 
            typeof, character(1)), length = vapply(parameters, 
            length, integer(1)), dimensions = vapply(parameters, 
            function(x) {
                if (is.null(dim(x))) 
                  "vector"
                else paste(dim(x), collapse = "x")
            }, character(1)), stringsAsFactors = FALSE)
    }
    phylo_cache_mapped_parameters <- function(parameters, template, 
        map, tolerance = 9.9999999999999995e-07) {
        phylo_cache_parameter_shapes(parameters)
        phylo_cache_parameter_shapes(template)
        phylo_cache_need(identical(names(parameters), names(template)), 
            "Blocs o ordre dels parametres inconsistents.")
        phylo_cache_need(is.list(map) && !is.null(names(map)) && 
            !anyDuplicated(names(map)) && all(names(map) %in% 
            names(template)), "Mapa de parametres invalid.")
        pieces <- lapply(names(template), function(name) {
            current <- parameters[[name]]
            expected <- template[[name]]
            phylo_cache_need(is.numeric(current) && is.numeric(expected) && 
                typeof(current) %in% c("integer", "double") && 
                typeof(expected) %in% c("integer", "double") && 
                identical(typeof(current), typeof(expected)) && 
                length(current) == length(expected) && identical(dim(current), 
                dim(expected)) && all(is.finite(current)), paste("Parametre o dimensions invalides:", 
                name))
            values <- as.numeric(current)
            mapping <- map[[name]]
            if (!is.null(mapping)) {
                phylo_cache_need(is.factor(mapping) && length(mapping) == 
                  length(values) && (is.null(dim(mapping)) || 
                  identical(dim(mapping), dim(expected))), paste("Mapa invalid:", 
                  name))
                codes <- as.integer(mapping)
                levels <- seq_along(levels(mapping))
                values <- vapply(levels, function(level) {
                  group <- values[which(codes == level)]
                  phylo_cache_need(length(group) > 0L && all(is.finite(group)) && 
                    max(abs(group - group[1L])) <= tolerance, 
                    paste("Valors compartits inconsistents:", 
                      name))
                  group[1L]
                }, numeric(1))
            }
            stats::setNames(values, rep(name, length(values)))
        })
        unlist(pieces, use.names = TRUE)
    }
    phylo_cache_validate_source <- function(fit, tolerance = 9.9999999999999995e-07, 
        gradient_threshold = 0.001) {
        phylo_cache_need(inherits(fit, "sdmTMB") && is.list(fit$model) && 
            is.list(fit$sd_report), "Cal un model sdmTMB guardat amb el seu report.")
        phylo_cache_need(isFALSE(fit$reml) && identical(fit$family$family, 
            "gaussian") && identical(fit$family$link, "identity") && 
            (is.null(fit$control$profile) || isFALSE(fit$control$profile)), 
            "Cal Gaussian ML sense profiling.")
        full <- phylo_cache_mapped_parameters(fit$parlist, fit$tmb_params, 
            fit$tmb_map, tolerance)
        stored <- fit$last.par.best
        phylo_cache_need(identical(names(full), names(stored)) && 
            phylo_cache_same_numeric(full, stored, tolerance), 
            "Els modes guardats no concorden amb last.par.best.")
        phylo_cache_need(is.character(fit$tmb_random) && !anyNA(fit$tmb_random) && 
            !anyDuplicated(fit$tmb_random) && all(fit$tmb_random %in% 
            names(fit$tmb_params)), "Blocs aleatoris inesperats.")
        outer <- full[!names(full) %in% fit$tmb_random]
        theta <- fit$model$par
        report <- fit$sd_report
        phylo_cache_need(identical(names(outer), names(theta)) && 
            phylo_cache_same_numeric(outer, theta, tolerance), 
            "L'estat latent correspon a altres parametres externs.")
        phylo_cache_need(identical(names(theta), names(report$par.fixed)) && 
            phylo_cache_same_numeric(theta, report$par.fixed, 
                tolerance), "Estimacions i report corresponen a punts diferents.")
        phylo_cache_need(is.matrix(report$cov.fixed) && nrow(report$cov.fixed) == 
            length(theta) && ncol(report$cov.fixed) == length(theta) && 
            all(is.finite(report$cov.fixed)) && all(diag(report$cov.fixed) > 
            0), "Covariancia fixa invalida.")
        fixed <- which(names(theta) == "b_j")
        phylo_cache_need(length(fixed) > 0L && phylo_cache_same_numeric(fit$parlist$b_j, 
            report$par.fixed[fixed], tolerance), "Coeficients fixos inconsistents.")
        labels <- names(report$value)
        phylo_cache_need(is.numeric(report$value) && !is.null(labels) && 
            !anyNA(labels) && length(report$sd) == length(report$value), 
            "Arrays ADREPORT invalids.")
        random <- grepl("^re_b_pars($|\\[)", labels)
        phylo_cache_need(any(random) && phylo_cache_same_numeric(as.vector(fit$parlist$re_b_pars), 
            report$value[random], tolerance) && all(is.finite(report$sd[random])) && 
            all(report$sd[random] >= 0), "Modes aleatoris o errors marginals inconsistents.")
        phylo_cache_need(isTRUE(fit$model$convergence == 0L) && 
            isTRUE(report$pdHess) && length(fit$gradients) > 
            0L && all(is.finite(fit$gradients)) && max(abs(fit$gradients)) < 
            gradient_threshold && length(fit$model$objective) == 
            1L && is.finite(fit$model$objective), "El model font no supera els controls numerics.")
        if (!is.null(report$gradient.fixed)) 
            phylo_cache_need(phylo_cache_same_numeric(as.numeric(fit$gradients), 
                as.numeric(report$gradient.fixed), tolerance), 
                "Gradients guardats inconsistents.")
        keep <- c("b_j", "ln_phi", "ln_tau_E", "ln_kappa", "re_cov_pars", 
            "re_b_pars", "epsilon_st")
        phylo_cache_need(all(keep %in% names(fit$parlist)), "Falten blocs de parametres inicials.")
        list(start = fit$parlist[keep], outer = outer, parameter_shapes = phylo_cache_parameter_shapes(fit$parlist), 
            random_blocks = fit$tmb_random, max_abs_gradient = max(abs(fit$gradients)), 
            source_objective = fit$model$objective, checks = "Saved arrays/report verified; no native source replay.")
    }
    phylo_cache_validate_portable_prepared <- function(prepared) {
        required <- c("data", "A", "formula", "mesh", "start", 
            "source", "signature", "priors", "extra_time")
        phylo_cache_need(is.list(prepared) && all(required %in% 
            names(prepared)) && is.data.frame(prepared$data) && 
            is.matrix(prepared$A) && is.list(prepared$start) && 
            is.list(prepared$source), "Estructura de cache incompleta.")
        inspect <- function(x, path) {
            if (inherits(x, "formula")) {
                phylo_cache_need(identical(environment(x), baseenv()), 
                  paste("Entorn de formula no portable:", path))
                entries <- as.list(x)
                attributes(entries) <- NULL
                inspect(entries, path)
                return(invisible())
            }
            phylo_cache_need(!typeof(x) %in% c("externalptr", 
                "environment", "closure", "builtin", "special", 
                "promise"), paste("Punter o entorn no portable:", 
                path))
            if (isS4(x)) {
                for (name in methods::slotNames(x)) inspect(methods::slot(x, 
                  name), paste0(path, "@", name))
            }
            else if (is.list(x) || is.pairlist(x) || is.call(x)) {
                entries <- as.list(x)
                labels <- names(entries)
                if (is.null(labels)) 
                  labels <- seq_along(entries)
                for (i in seq_along(entries)) inspect(entries[[i]], 
                  paste0(path, "$", labels[i]))
            }
            attrs <- attributes(x)
            if (length(attrs)) 
                for (name in names(attrs)) inspect(attrs[[name]], 
                  paste0(path, " attr:", name))
            invisible()
        }
        inspect(prepared, "prepared")
        invisible(TRUE)
    }
    phylo_cache_strip_unused_source <- function(saved) {
        phylo_cache_need(is.list(saved) && is.list(saved$fit), 
            "Registre del model font invalid.")
        saved$fit$sd_report[["cov"]] <- NA_real_
        saved$fit$sd_report$env <- NULL
        saved$fit$tmb_obj <- NULL
        saved
    }
    portable_source_start <- function(fit, directory) {
        verified <- phylo_cache_validate_source(fit)
        write_csv(data.frame(method = "Saved mapped state and report verified without native replay", 
            native_source_replay_run = FALSE, max_abs_gradient_saved = verified$max_abs_gradient, 
            n_outer_parameters = length(verified$outer), parameter_sha256 = digest::digest(fit$parlist, 
                algo = "sha256")), file.path(directory, "source_start_validation.csv"))
        verified$start
    }
    prepare <- function(cfg, packages) {
        src <- dirname(cfg$source_checkpoint)
        inputs <- c(cfg$source_checkpoint, cfg$tree_file, file.path(dirname(src), 
            "run_config.rds"), file.path(src, c("completed.rds", 
            "fit_checks.csv", "qc_filter_summary.csv", "onset_scaling_after_qc.csv", 
            "model_specification.txt")))
        need(all(file.exists(inputs)), "Falten entrades QC per validar la cache.")
        hashes <- unname(tools::md5sum(inputs))
        need(!anyNA(hashes), "No es poden calcular els hashes.")
        manifest <- list(schema = "portable_phylo_source_edges_1.0", 
            engine_version = version, paths = inputs, md5 = hashes, 
            packages = packages, seed = cfg$fit_seed, gradient_threshold = cfg$gradient_threshold)
        key <- digest::digest(manifest, algo = "sha256")
        parent <- file.path(cfg$output_root, "source_cache")
        dir.create(parent, recursive = TRUE, showWarnings = FALSE)
        cache_file <- file.path(parent, paste0("prepared_", substr(key, 
            1L, 12L), ".rds"))
        if (file.exists(cache_file)) {
            progress("Reading validated portable source cache; no source replay")
            cached <- readRDS(cache_file)
            need(identical(cached$manifest, manifest), "Manifest de cache inconsistent.")
            p <- cached$prepared
            phylo_cache_validate_portable_prepared(p)
            fingerprint <- digest::digest(list(p$data, p$A, p$formula, 
                p$mesh, p$start, p$source, p$signature, p$priors, 
                p$extra_time), algo = "sha256")
            need(identical(fingerprint, cached$fingerprint), 
                "Contingut de cache inconsistent.")
            need(8 * nrow(p$data) * ncol(attr(p$A, "edge_basis"))/1024^3 <= 
                cfg$maximum_single_basis_GiB, "La base excedeix el limit configurat.")
            p$out <- file.path(cfg$output_root, basename(p$out))
            dir.create(p$out, recursive = TRUE, showWarnings = FALSE)
            write_csv(data.frame(cache_file = cache_file, manifest_sha256 = key, 
                reused = TRUE, native_source_replay_run = FALSE), 
                file.path(p$out, "source_cache_used.csv"))
            return(p)
        }
        progress(paste("Preparing portable source once;", "checking saved arrays instead of replaying native likelihood"))
        p <- prepare_uncached_source(cfg, packages)
        phylo_cache_validate_portable_prepared(p)
        fingerprint <- digest::digest(list(p$data, p$A, p$formula, 
            p$mesh, p$start, p$source, p$signature, p$priors, 
            p$extra_time), algo = "sha256")
        save_atomic(list(manifest = manifest, fingerprint = fingerprint, 
            prepared = p), cache_file)
        write_csv(data.frame(cache_file = cache_file, manifest_sha256 = key, 
            reused = FALSE, native_source_replay_run = FALSE), 
            file.path(p$out, "source_cache_used.csv"))
        progress("Portable source cache saved before full joint Hessian construction")
        p
    }
    {
        construct_phylo <- e$construct_phylo
        body(construct_phylo) <- replace_phylo_basis(body(construct_phylo))
        preflight <- e$preflight
        body(preflight) <- replace_phylo_preflight(body(preflight))
        build_joint <- e$build_joint
        body(build_joint) <- replace_phylo_call_argument(body(build_joint), 
            quote(sdmTMB::sdmTMBcontrol), "trace", 0L)
        starter <- engine$start
        body(starter) <- replace_phylo_call_argument(body(starter), 
            quote(callr::r_bg), "stderr", "2>&1")
        engine$start <- starter
    }
    for (name in c("fixed_table", "prepare", "prepare_uncached_source", 
        "synchronize_fit", "preflight_onset", "worker", "construct_phylo", 
        "preflight", "build_joint", "build_phylo_edge_basis", 
        "phylo_cache_need", "phylo_cache_same_numeric", "phylo_cache_parameter_shapes", 
        "phylo_cache_mapped_parameters", "phylo_cache_validate_source", 
        "phylo_cache_validate_portable_prepared", "phylo_cache_strip_unused_source", 
        "portable_source_start")) {
        f <- get(name)
        environment(f) <- e
        assign(name, f, envir = e)
    }
    cfg <- private$phylogeny_pilot_config
    cfg$project_root <- root
    cfg$source_checkpoint <- file.path(root, "output/phenology_plasticity/offset_with_onset_qc_lastpeak_15km", 
        qc_run, "offset_multivoltine/fit_qc.rds")
    cfg$output_root <- file.path(root, "output/phenology_plasticity/phylogeny_qc_15km_v2")
    cfg$check_only <- check_only
    cfg$diagnostic_only <- FALSE
    cfg$fit_seed <- 20261007L
    cfg$checkpoint_seconds <- checkpoint_seconds
    cfg$maximum_single_basis_GiB <- maximum_single_basis_GiB
    cfg$extra_optimization_attempts <- extra_optimization_attempts
    job <- suppressMessages(engine$start(cfg))
    job$engine <- engine
    options(pheno.phylo.qc.job = job)
    message("Model filogenetic post-QC de 15 km iniciat en un worker.\n", 
        "Seguiment: watch_phylogeny_qc_15km(). Mante R obert.\n", 
        "Log: ", job$log_file)
    invisible(job)
}
status_phylogeny_qc_15km <-
function(
        job = getOption("pheno.phylo.qc.job"), quiet = FALSE) {
    if (is.null(job)) stop("Executa run_phylogeny_qc_15km() primer.")
    job$engine$status(job, quiet)
}
watch_phylogeny_qc_15km <-
function(
        job = getOption("pheno.phylo.qc.job"), interval_seconds = 60) {
    if (is.null(job)) stop("Executa run_phylogeny_qc_15km() primer.")
    job$engine$watch(job, interval_seconds)
}
stop_phylogeny_qc_15km <-
function(job = getOption("pheno.phylo.qc.job")) {
    if (is.null(job)) return(invisible(FALSE))
    job$engine$stop(job)
}
