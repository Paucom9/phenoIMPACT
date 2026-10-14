# phenoIMPACT: diagnostic filogenetic lleuger.
# Llegeix nomes Newick i CSV petits; no carrega models ni fa refits.
# Dependencia: ape.
#
# Lambda i la permutacio descriuen estimacions parcialment agrupades,
# no plasticitats latents observades. Les SE marginals no proporcionen
# la covariancia conjunta entre pendents.

pheno_phylogenetic_screen <- local({
  
  ANALYSES <- c("onset", "offset_univoltine", "offset_multivoltine")
  
  need <- function(ok, msg) {
    if (!isTRUE(ok)) stop(msg, call. = FALSE)
  }
  
  csv <- function(x, path) {
    utils::write.csv(
      x, path, row.names = FALSE, na = "",
      fileEncoding = "UTF-8"
    )
  }
  
  read_csv <- function(path) {
    utils::read.csv(
      path, stringsAsFactors = FALSE,
      check.names = FALSE, fileEncoding = "UTF-8-BOM"
    )
  }
  
  norm <- function(x) {
    tolower(trimws(gsub(
      "[[:space:]]+", " ",
      gsub("_", " ", as.character(x), fixed = TRUE)
    )))
  }
  
  valid_names <- function(x) {
    length(x) > 0 &&
      !anyNA(x) &&
      all(nzchar(norm(x))) &&
      !anyDuplicated(norm(x))
  }
  
  absolute <- function(x) {
    grepl("^(/|[A-Za-z]:[/\\\\]|\\\\\\\\)", x)
  }
  
  resolve <- function(x, root) {
    if (!absolute(x)) x <- file.path(root, x)
    normalizePath(x, winslash = "/", mustWork = TRUE)
  }
  
  provenance <- function(folder, analysis = NULL) {
    
    p <- read_csv(file.path(folder, "provenance.csv"))
    
    need(
      nrow(p) == 1 &&
        all(c(
          "analysis", "formula", "offset_adjusted_for_onset"
        ) %in% names(p)),
      "Provenance incompleta: cal analysis, formula i offset_adjusted_for_onset."
    )
    
    need(
      p$analysis %in% ANALYSES &&
        (is.null(analysis) || p$analysis == analysis),
      "La resposta de provenance.csv no coincideix amb la carpeta seleccionada."
    )
    
    need(
      !is.na(p$formula) && nzchar(trimws(p$formula)),
      "Falta la formula original a provenance.csv."
    )
    
    adjusted <- as.logical(p$offset_adjusted_for_onset)
    need(!is.na(adjusted), "offset_adjusted_for_onset no es TRUE/FALSE.")
    
    old <- all(c(
      "cutoff_km", "no_refit", "checkpoint_fingerprint"
    ) %in% names(p))
    
    updated <- all(c(
      "source_run", "source_checkpoint", "source_signature",
      "new_phenology_fit", "max_abs_gradient"
    ) %in% names(p))
    
    if (old) {
      
      need(
        isTRUE(p$cutoff_km == 15) &&
          isTRUE(as.logical(p$no_refit)) &&
          !is.na(p$checkpoint_fingerprint) &&
          nzchar(trimws(p$checkpoint_fingerprint)),
        "Provenance original invalida o malla diferent de 15 km."
      )
      
      need(
        file.exists(file.path(folder, "saved_fit_checks.csv")),
        "Falta saved_fit_checks.csv de l'extraccio original."
      )
      
      stage <- if (adjusted) {
        "original_export_onset_adjusted"
      } else {
        "original_export_not_onset_adjusted"
      }
      
    } else if (updated) {
      
      need(
        p$analysis != "onset" &&
          isTRUE(adjusted) &&
          isFALSE(as.logical(p$new_phenology_fit)) &&
          !is.na(p$source_signature) &&
          nzchar(trimws(p$source_signature)) &&
          !is.na(p$source_checkpoint) &&
          nzchar(trimws(p$source_checkpoint)) &&
          !is.na(p$source_run) &&
          nzchar(trimws(p$source_run)) &&
          is.finite(p$max_abs_gradient) &&
          p$max_abs_gradient >= 0 &&
          p$max_abs_gradient < 0.001,
        "Provenance invalida dels offsets ajustats per onset."
      )
      
      stage <- "updated_onset_adjusted_export_QC_stage_unverified"
      
    } else {
      stop(
        "Format de provenance no reconegut; no es pressuposa la procedencia.",
        call. = FALSE
      )
    }
    
    list(table = p, stage = stage, adjusted = adjusted)
  }
  
  inventory <- function(root) {
    
    files <- list.files(
      file.path(root, "output"),
      pattern = "^species_random_slopes[.]csv$",
      recursive = TRUE, full.names = TRUE
    )
    
    rows <- lapply(files, function(f) {
      
      folder <- dirname(f)
      
      info <- tryCatch(
        provenance(folder),
        error = function(e) e
      )
      
      data.frame(
        analysis = if (inherits(info, "error")) {
          basename(folder)
        } else {
          info$table$analysis
        },
        source_dir = normalizePath(
          folder, winslash = "/", mustWork = TRUE
        ),
        usable = !inherits(info, "error"),
        model_stage = if (inherits(info, "error")) {
          NA_character_
        } else {
          info$stage
        },
        note = if (inherits(info, "error")) {
          conditionMessage(info)
        } else {
          ""
        },
        stringsAsFactors = FALSE
      )
    })
    
    if (!length(rows)) {
      return(data.frame(
        analysis = character(),
        source_dir = character(),
        usable = logical(),
        model_stage = character(),
        note = character()
      ))
    }
    
    do.call(rbind, rows)
  }
  
  # C_lambda = lambda*C + (1-lambda)*diag(diag(C)).
  # Es conserva l'arrel original i nomes s'aplica una escala comuna.
  # La descomposicio espectral evita inversions a cada permutacio.
  
  basis <- function(C) {
    
    need(
      is.matrix(C) &&
        nrow(C) == ncol(C) &&
        nrow(C) >= 4 &&
        all(is.finite(C)) &&
        isTRUE(all.equal(C, t(C), tolerance = 1e-10)) &&
        all(diag(C) > 0),
      "Covariancia browniana invalida."
    )
    
    C <- C / mean(diag(C))
    d <- sqrt(diag(C))
    R <- C / outer(d, d)
    e <- eigen(R, symmetric = TRUE)
    
    need(
      min(e$values) > 0,
      "Covariancia singular: no es reparen branques ni matrius automaticament."
    )
    
    list(
      U = e$vectors,
      values = e$values,
      d = d,
      z1 = as.numeric(crossprod(e$vectors, 1 / d)),
      logD = sum(log(d * d))
    )
  }
  
  fit_lambda <- function(y, b, profile = FALSE) {
    
    need(
      length(y) == length(b$d) && all(is.finite(y)),
      "Pendents invalides."
    )
    
    scale <- stats::sd(y)
    
    need(
      is.finite(scale) && scale > 0,
      "Pendents constants: lambda no es identificable."
    )
    
    y <- (y - mean(y)) / scale
    z <- as.numeric(crossprod(b$U, y / b$d))
    n <- length(z)
    
    ll <- function(lambda) {
      
      w <- 1 + lambda * (b$values - 1)
      if (any(w <= 0)) return(-Inf)
      
      mu <- sum(b$z1 * z / w) / sum(b$z1 * b$z1 / w)
      q <- sum((z - mu * b$z1)^2 / w)
      
      if (!is.finite(q) || q <= 0) return(-Inf)
      
      -0.5 * (
        n * (log(2 * pi) + 1 + log(q / n)) +
          b$logD + sum(log(w))
      )
    }
    
    grid <- seq(0, 1, length.out = 21)
    vals <- vapply(grid, ll, numeric(1))
    
    need(
      all(is.finite(vals)),
      "Versemblanca no finita; revisar arbre i pendents."
    )
    
    maxima <- which(
      vals[2:20] >= vals[1:19] &
        vals[2:20] >= vals[3:21]
    ) + 1L
    
    brackets <- unique(c(maxima, which.max(vals)))
    candidates <- c(0, 1, grid[which.max(vals)])
    
    for (j in brackets) {
      
      opt <- stats::optimize(
        ll,
        c(grid[max(1, j - 1)], grid[min(21, j + 1)]),
        maximum = TRUE, tol = 1e-7
      )
      
      candidates <- c(candidates, opt$maximum)
    }
    
    scores <- vapply(candidates, ll, numeric(1))
    best <- max(scores)
    
    # Si el perfil es numericament pla, s'escull el lambda menor.
    chosen <- which(scores >= best - 1e-9)
    lambda <- min(candidates[chosen])
    
    list(
      lambda = lambda,
      LR = max(0, 2 * (best - ll(0))),
      at_boundary = lambda < 1e-5 || lambda > 1 - 1e-5,
      profile_logLik_span = diff(range(vals)),
      profile = if (profile) {
        data.frame(
          lambda = grid,
          relative_logLik = vals - best
        )
      } else {
        NULL
      }
    )
  }
  
  tree_checks <- function(tree) {
    
    need(
      inherits(tree, "phylo") &&
        !inherits(tree, "multiPhylo") &&
        valid_names(tree$tip.label),
      "Cal un sol arbre amb noms d'especie unics."
    )
    
    need(
      ape::is.rooted(tree),
      "L'arbre no esta arrelat; no s'arrela automaticament."
    )
    
    need(
      length(tree$edge.length) == nrow(tree$edge) &&
        all(is.finite(tree$edge.length)) &&
        all(tree$edge.length >= 0),
      "Calen longituds de branca finites i no negatives."
    )
    
    data.frame(
      n_tips = length(tree$tip.label),
      rooted = TRUE,
      ultrametric = ape::is.ultrametric(tree),
      zero_length_edges = sum(tree$edge.length == 0)
    )
  }
  
  read_map <- function(path, tree) {
    
    if (is.null(path) || !file.exists(path)) return(NULL)
    
    m <- read_csv(path)
    
    need(
      all(c("SPECIES", "tip_label") %in% names(m)) &&
        valid_names(m$SPECIES) &&
        valid_names(m$tip_label) &&
        all(norm(m$tip_label) %in% norm(tree$tip.label)),
      "Mapa de noms invalid: calen parelles SPECIES,tip_label verificades i sense duplicats."
    )
    
    m
  }
  
  plot_result <- function(tree, x, fit, permutations, out) {
    
    tr <- ape::ladderize(ape::keep.tip(tree, x$tip_label))
    x <- x[match(tr$tip.label, x$tip_label), , drop = FALSE]
    
    grDevices::pdf(
      file.path(out, "diagnostic.pdf"),
      width = 12,
      height = max(7, 0.16 * nrow(x) + 2),
      useDingbats = FALSE
    )
    
    on.exit(grDevices::dev.off(), add = TRUE)
    
    graphics::layout(matrix(1:2, 1), widths = c(0.6, 0.4))
    graphics::par(mar = c(4, 1, 2, 1))
    
    ape::plot.phylo(
      tr, cex = 0.65, font = 3, no.margin = FALSE
    )
    
    graphics::title("Species tree and conditional slope deviations")
    
    pp <- get(
      "last_plot.phylo",
      envir = get(".PlotPhyloEnv", envir = asNamespace("ape"))
    )
    
    yy <- pp$yy[seq_len(nrow(x))]
    lo <- x$sp_slope_dev - 1.96 * x$sp_slope_dev_marginal_SE
    hi <- x$sp_slope_dev + 1.96 * x$sp_slope_dev_marginal_SE
    
    graphics::par(mar = c(4, 1, 2, 1))
    
    graphics::plot(
      range(c(lo, hi, 0)),
      range(yy) + c(-0.5, 0.5),
      type = "n", yaxt = "n",
      xlab = "Days / fitted anomaly unit",
      ylab = "",
      main = "Deviation +/- 1.96 marginal SE"
    )
    
    graphics::abline(v = 0, lty = 2, col = "grey60")
    graphics::segments(lo, yy, hi, yy, col = "grey60")
    
    graphics::points(
      x$sp_slope_dev, yy, pch = 16, col = "#246B91"
    )
    
    graphics::layout(matrix(1:2, 1))
    graphics::par(mar = c(4, 4, 3, 1))
    
    graphics::plot(
      fit$profile$lambda,
      fit$profile$relative_logLik,
      type = "l",
      xlab = "Lambda (restricted to 0-1)",
      ylab = "Relative ML log-likelihood"
    )
    
    graphics::abline(v = fit$lambda, col = "#AE6A1B", lty = 2)
    graphics::title("Profile of exported point estimates")
    
    graphics::hist(
      permutations$LR,
      breaks = "FD",
      main = "Conditional label permutations",
      xlab = "Likelihood-ratio statistic",
      col = "grey85", border = "white"
    )
    
    graphics::abline(v = fit$LR, col = "#AE6A1B", lwd = 2)
  }
  
  one_analysis <- function(
    a, folder, tree, C, map, out,
    B, seed, loo, allow_unmatched
  ) {
    
    dir.create(out)
    
    paths <- file.path(
      folder, c("species_random_slopes.csv", "provenance.csv")
    )
    
    info <- provenance(folder, a)
    x <- read_csv(paths[1])
    
    needed <- c(
      "SPECIES", "sp_slope_dev", "sp_slope_dev_marginal_SE"
    )
    
    need(
      all(needed %in% names(x)),
      "Falten columnes de species_random_slopes.csv."
    )
    
    x <- x[, needed, drop = FALSE]
    
    need(
      valid_names(x$SPECIES) &&
        is.numeric(x$sp_slope_dev) &&
        is.numeric(x$sp_slope_dev_marginal_SE) &&
        all(is.finite(x$sp_slope_dev)) &&
        all(is.finite(x$sp_slope_dev_marginal_SE)) &&
        all(x$sp_slope_dev_marginal_SE >= 0),
      "Identificadors, pendents o SE invalids/duplicats."
    )
    
    if ("n_species" %in% names(info$table)) {
      need(
        isTRUE(info$table$n_species == nrow(x)),
        "n_species de provenance.csv no coincideix amb la taula de pendents."
      )
    }
    
    key <- norm(x$SPECIES)
    
    if (!is.null(map)) {
      
      j <- match(key, norm(map$SPECIES))
      
      key[!is.na(j)] <- norm(
        map$tip_label[j[!is.na(j)]]
      )
    }
    
    x$tip_label <- tree$tip.label[
      match(key, norm(tree$tip.label))
    ]
    
    x$matched <- !is.na(x$tip_label)
    
    need(
      !anyDuplicated(x$tip_label[x$matched]),
      "Diverses especies coincideixen amb una sola punta de l'arbre."
    )
    
    csv(x, file.path(out, "species_matching.csv"))
    
    csv(
      x[!x$matched, , drop = FALSE],
      file.path(out, "unmatched_species.csv")
    )
    
    csv(info$table, file.path(out, "source_provenance.csv"))
    
    checks_file <- file.path(folder, "saved_fit_checks.csv")
    source_pass <- NA
    
    if (file.exists(checks_file)) {
      
      checks <- read_csv(checks_file)
      csv(checks, file.path(out, "source_saved_fit_checks.csv"))
      
      if (
        "all_sanity_checks_pass" %in% names(checks) &&
        nrow(checks) > 0
      ) {
        source_pass <- all(
          !is.na(checks$all_sanity_checks_pass) &
            as.logical(checks$all_sanity_checks_pass)
        )
      }
      
      paths <- c(paths, checks_file)
    }
    
    csv(
      data.frame(
        path = normalizePath(
          paths, winslash = "/", mustWork = TRUE
        ),
        md5 = unname(tools::md5sum(paths))
      ),
      file.path(out, "input_files.csv")
    )
    
    need(
      all(x$matched) || allow_unmatched,
      paste(
        "Hi ha especies sense correspondencia.",
        "Revisa unmatched_species.csv i el mapa de noms;",
        "no s'eliminen automaticament."
      )
    )
    
    total <- nrow(x)
    x <- x[x$matched, , drop = FALSE]
    x <- x[order(x$tip_label), , drop = FALSE]
    
    need(
      nrow(x) >= 5,
      "Calen almenys cinc especies coincidents; amb poques especies el diagnostic es feble."
    )
    
    # Submatriu de la covariancia de l'arbre ORIGINAL.
    cov <- C[x$tip_label, x$tip_label, drop = FALSE]
    b <- basis(cov)
    
    observed <- fit_lambda(
      x$sp_slope_dev, b, profile = TRUE
    )
    
    permutations <- data.frame(
      iteration = seq_len(B),
      lambda = NA_real_,
      LR = NA_real_
    )
    
    set.seed(seed)
    
    for (i in seq_len(B)) {
      
      f <- fit_lambda(
        sample(x$sp_slope_dev, replace = FALSE), b
      )
      
      permutations$lambda[i] <- f$lambda
      permutations$LR[i] <- f$LR
      
      if (i %% 250 == 0) {
        message(a, ": ", i, "/", B, " permutacions")
      }
    }
    
    p <- (
      1 + sum(permutations$LR >= observed$LR - 1e-9)
    ) / (B + 1)
    
    loo_table <- data.frame(
      omitted_species = character(),
      lambda = numeric()
    )
    
    if (loo) {
      
      estimates <- vapply(
        seq_len(nrow(x)),
        function(i) {
          
          y <- x$sp_slope_dev[-i]
          if (stats::sd(y) == 0) return(NA_real_)
          
          fit_lambda(
            y,
            basis(cov[-i, -i, drop = FALSE])
          )$lambda
        },
        numeric(1)
      )
      
      loo_table <- data.frame(
        omitted_species = x$SPECIES,
        lambda = estimates
      )
    }
    
    x$deviation_lower_approx95 <-
      x$sp_slope_dev - 1.96 * x$sp_slope_dev_marginal_SE
    
    x$deviation_upper_approx95 <-
      x$sp_slope_dev + 1.96 * x$sp_slope_dev_marginal_SE
    
    x$marginal_SE_over_observed_slope_SD_descriptive <-
      x$sp_slope_dev_marginal_SE / stats::sd(x$sp_slope_dev)
    
    csv(x, file.path(out, "matched_species_slopes.csv"))
    csv(observed$profile, file.path(out, "lambda_profile.csv"))
    csv(permutations, file.path(out, "label_permutations.csv"))
    csv(loo_table, file.path(out, "leave_one_species_out.csv"))
    
    plot_result(tree, x, observed, permutations, out)
    
    bounds <- if (any(is.finite(loo_table$lambda))) {
      range(loo_table$lambda, na.rm = TRUE)
    } else {
      c(NA, NA)
    }
    
    data.frame(
      analysis = a,
      status = "EXPLORATORY_SCREEN_COMPLETED",
      model_stage = info$stage,
      offset_adjusted_for_onset = info$adjusted,
      n_species_exported = total,
      n_species_matched = nrow(x),
      lambda_point_estimate_exploratory = observed$lambda,
      lambda_at_boundary = observed$at_boundary,
      profile_logLik_span = observed$profile_logLik_span,
      LR_point_estimates = observed$LR,
      p_label_permutation_exploratory = p,
      n_permutations = B,
      permutation_MC_SE_approx = sqrt(p * (1 - p) / (B + 1)),
      lambda_LOO_min = bounds[1],
      lambda_LOO_max = bounds[2],
      median_marginal_SE_over_slope_SD_descriptive =
        stats::median(
          x$marginal_SE_over_observed_slope_SD_descriptive
        ),
      source_fit_checks_all_pass = source_pass,
      joint_slope_uncertainty_available = FALSE,
      formal_latent_signal_test_run = FALSE,
      note = if (nrow(x) < 20) {
        "Poques especies: perfil i potencia limitats."
      } else {
        ""
      },
      stringsAsFactors = FALSE
    )
  }
  
  run <- function(
    project_root = "E:/phenoIMPACT project/code/phenoIMPACT",
    source_dirs = NULL,
    analyses = ANALYSES,
    tree_file = "data/EUROPEAN_BUTTERFLIES_MCC_CA.nwk",
    name_map_file = "data/phylogeny_species_map.csv",
    n_permutations = 999L,
    seed = 20261006L,
    leave_one_out = TRUE,
    allow_unmatched = FALSE,
    choose_sources = interactive()
  ) {
    
    need(
      requireNamespace("ape", quietly = TRUE),
      "Instal.la ape: install.packages('ape')"
    )
    
    need(
      length(analyses) > 0 &&
        !anyDuplicated(analyses) &&
        all(analyses %in% ANALYSES),
      "Analyses invalides."
    )
    
    need(
      length(n_permutations) == 1 &&
        is.finite(n_permutations) &&
        n_permutations >= 19 &&
        n_permutations == as.integer(n_permutations),
      "n_permutations ha de ser un enter >=19; recomanat 999."
    )
    
    need(
      length(seed) == 1 &&
        is.finite(seed) &&
        seed >= 0 &&
        seed <= .Machine$integer.max - length(ANALYSES),
      "Seed invalid."
    )
    
    need(
      is.logical(leave_one_out) &&
        length(leave_one_out) == 1 &&
        !is.na(leave_one_out) &&
        is.logical(allow_unmatched) &&
        length(allow_unmatched) == 1 &&
        !is.na(allow_unmatched),
      "Opcions logiques invalides."
    )
    
    need(
      is.logical(choose_sources) &&
        length(choose_sources) == 1 &&
        !is.na(choose_sources),
      "choose_sources ha de ser TRUE/FALSE."
    )
    
    root <- normalizePath(
      project_root, winslash = "/", mustWork = TRUE
    )
    
    had_seed <- exists(
      ".Random.seed", envir = .GlobalEnv, inherits = FALSE
    )
    
    if (had_seed) {
      previous_seed <- get(".Random.seed", envir = .GlobalEnv)
    }
    
    on.exit({
      if (had_seed) {
        assign(".Random.seed", previous_seed, envir = .GlobalEnv)
      } else if (exists(
        ".Random.seed", envir = .GlobalEnv, inherits = FALSE
      )) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    
    parent <- file.path(
      root, "output", "phenology_plasticity",
      "phylogenetic_screen"
    )
    
    dir.create(parent, recursive = TRUE, showWarnings = FALSE)
    
    out <- tempfile(
      paste0(
        "screen_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_"
      ),
      tmpdir = parent
    )
    
    need(
      dir.create(out),
      "No es pot crear una carpeta nova de resultats."
    )
    
    out <- normalizePath(out, winslash = "/", mustWork = TRUE)
    
    candidates <- inventory(root)
    csv(candidates, file.path(out, "source_candidates.csv"))
    
    qc_root <- file.path(
      root, "output", "phenology_plasticity",
      "offset_with_onset_qc_lastpeak_15km"
    )
    
    qc_dirs <- if (dir.exists(qc_root)) {
      list.dirs(qc_root, recursive = FALSE, full.names = TRUE)
    } else {
      character()
    }
    
    csv(
      data.frame(QC_run_directory = qc_dirs),
      file.path(out, "QC_run_directories_present.csv")
    )
    
    if (length(qc_dirs)) {
      message(
        "Hi ha carpetes de refits QC. ",
        "La seva presencia no confirma finalitzacio; ",
        "els CSV disponibles poden ser pre-QC."
      )
    }
    
    if (is.null(source_dirs)) {
      
      selected <- lapply(analyses, function(a) {
        candidates$source_dir[
          candidates$usable & candidates$analysis == a
        ]
      })
      
      if (
        choose_sources &&
        interactive() &&
        all(lengths(selected) >= 1L)
      ) {
        
        for (i in seq_along(analyses)) {
          
          if (length(selected[[i]]) > 1L) {
            
            rows <- candidates[
              candidates$usable &
                candidates$analysis == analyses[i],
              , drop = FALSE
            ]
            
            choices <- paste(
              rows$model_stage, rows$source_dir, sep = " | "
            )
            
            choice <- utils::select.list(
              choices,
              multiple = FALSE,
              graphics = FALSE,
              title = paste(
                "Tria la font de", analyses[i],
                "(0 = cancel.lar)"
              )
            )
            
            if (!nzchar(choice)) break
            
            selected[[i]] <- rows$source_dir[
              match(choice, choices)
            ]
          }
        }
      }
      
      if (!all(lengths(selected) == 1L)) {
        
        print(candidates, row.names = FALSE)
        
        message(
          "No s'ha calculat cap diagnostic. ",
          "Tria source_dirs amb una carpeta per resposta.\n",
          "Inventari: ", out
        )
        
        return(invisible(list(
          status = "SOURCE_SELECTION_REQUIRED",
          output_dir = out,
          candidates = candidates
        )))
      }
      
      source_dirs <- stats::setNames(
        vapply(selected, `[`, character(1), 1),
        analyses
      )
    }
    
    need(
      is.character(source_dirs) &&
        !is.null(names(source_dirs)) &&
        !anyDuplicated(names(source_dirs)) &&
        setequal(names(source_dirs), analyses),
      paste(
        "source_dirs ha de ser un vector de carpetes",
        "amb els noms de les analyses seleccionades."
      )
    )
    
    source_dirs <- vapply(
      source_dirs[analyses], resolve,
      character(1), root = root
    )
    
    csv(
      data.frame(
        analysis = analyses,
        source_dir = unname(source_dirs)
      ),
      file.path(out, "selected_sources.csv")
    )
    
    tree_path <- resolve(tree_file, root)
    tree <- ape::read.tree(tree_path)
    audit <- tree_checks(tree)
    
    csv(audit, file.path(out, "tree_checks.csv"))
    
    if (!audit$ultrametric) {
      message(
        "Arbre no ultrametric: es conserva C browniana original, ",
        "amb variancies diagonals desiguals; ",
        "permutacio nomes exploratoria."
      )
    }
    
    C <- ape::vcv.phylo(tree, corr = FALSE)
    need(all(is.finite(C)), "Covariancia browniana no finita.")
    
    csv(
      data.frame(
        path = tree_path,
        md5 = unname(tools::md5sum(tree_path))
      ),
      file.path(out, "tree_source.csv")
    )
    
    map_path <- if (is.null(name_map_file)) {
      NULL
    } else if (absolute(name_map_file)) {
      name_map_file
    } else {
      file.path(root, name_map_file)
    }
    
    map <- read_map(map_path, tree)
    
    if (!is.null(map)) {
      
      csv(map, file.path(out, "name_map_used.csv"))
      
      csv(
        data.frame(
          path = normalizePath(
            map_path, winslash = "/", mustWork = TRUE
          ),
          md5 = unname(tools::md5sum(map_path))
        ),
        file.path(out, "name_map_source.csv")
      )
    }
    
    writeLines(c(
      "DIAGNOSTIC EXPLORATORI DE FILOGENIA EN PENDENTS DE PLASTICITAT",
      "Llegeix nomes arbre i CSV petits. No carrega RDS ni reajusta models.",
      "Resposta: sp_slope_dev, desviacio per especie de la pendent d'anomalia termica, condicionada als moderadors.",
      "Zero no vol dir plasticitat nul.la. Sumar el mateix coeficient fix o invertir el signe no canvia lambda.",
      "Unitats: dies per unitat FITTED d'anomalia; no es pressuposa dies/grau C.",
      "",
      "LAMBDA: ML amb mitjana GLS i variancia estimades; restringit a [0,1].",
      "C_lambda = lambda*C + (1-lambda)*diag(diag(C)), C browniana amb l'arrel original.",
      "Lambda=0 elimina covariancia entre especies; lambda=1 conserva el patro brownia. No es percentatge de variancia.",
      "Als arbres no ultrametrics, lambda=0 encara conserva variancies diagonals desiguals.",
      "Un limit lambda=0/1 o un perfil pla requereix prudencia; no es dona un IC formal.",
      "",
      "INFERENCIA: pendents parcialment agrupades, no trets observats independents.",
      "p_label_permutation_exploratory compara l'associacio entre arbre i estimacions puntuals amb etiquetes permutades.",
      "NO es un test calibrat d'absencia de senyal en les plasticitats latents: precisio, shrinkage i covariancia desiguals trenquen intercanviabilitat.",
      "No s'aplica una decisio automatica per P<0.05 ni cap llindar de lambda.",
      "SE/SD i barres +/-1.96 SE son descriptives per a la desviacio; no intervals de plasticitat total.",
      "No es simulen errors independents: els CSV no contenen la incertesa conjunta de les pendents.",
      "Leave-one-species-out mesura influencia d'una especie, no propaga aquesta incertesa.",
      "Un sol arbre no cobreix incertesa filogenetica; resultats no significatius no proven absencia de senyal.",
      "Referencia: Houslay & Wilson 2017, doi:10.1093/beheco/arx023.",
      "",
      "PROCEDENCIA: revisar selected_sources.csv i cada source_provenance.csv abans d'interpretar.",
      "Original export: offsets no ajustats per onset. Updated export: offsets ajustats per onset, revisio residual pot ser PENDING.",
      "El refit QC LAST_PEAK no exporta aquests CSV automaticament. No es confirma que una exportacio sigui la definitiva.",
      "Es copien els avisos originals; aquest diagnostic no valida els ajusts font.",
      "Sense parelles taxonomiques verificades no es fa matching fuzzy ni s'eliminen especies per defecte.",
      paste("Permutacions:", n_permutations, "Seed base:", seed),
      "Totes les sortides es creen en una carpeta nova. Arbre, models i taules originals queden intactes."
    ), file.path(out, "LLEGEIX_PRIMER.txt"), useBytes = TRUE)
    
    summaries <- list()
    
    for (a in analyses) {
      
      message("Diagnostic: ", a)
      
      result <- tryCatch(
        one_analysis(
          a, source_dirs[[a]], tree, C, map,
          file.path(out, a),
          as.integer(n_permutations),
          as.integer(seed + match(a, ANALYSES)),
          leave_one_out, allow_unmatched
        ),
        error = function(e) {
          data.frame(
            analysis = a,
            status = "ERROR",
            note = conditionMessage(e),
            stringsAsFactors = FALSE
          )
        }
      )
      
      summaries[[a]] <- result
      keys <- unique(unlist(lapply(summaries, names)))
      
      aligned <- lapply(summaries, function(x) {
        
        for (k in setdiff(keys, names(x))) x[[k]] <- NA
        
        x[, keys, drop = FALSE]
      })
      
      summary <- do.call(rbind, aligned)
      csv(summary, file.path(out, "screen_summary.csv"))
    }
    
    writeLines(
      capture.output(utils::sessionInfo()),
      file.path(out, "sessionInfo.txt")
    )
    
    print(
      summary[, intersect(c(
        "analysis",
        "status",
        "lambda_point_estimate_exploratory",
        "p_label_permutation_exploratory",
        "lambda_LOO_min",
        "lambda_LOO_max",
        "note"
      ), names(summary)), drop = FALSE],
      row.names = FALSE
    )
    
    message(
      "Resultats: ", out,
      "\nLlegeix LLEGEIX_PRIMER.txt i la procedencia. ",
      "No s'ha fet un test formal sobre plasticitats latents."
    )
    
    invisible(list(
      status = if (all(
        summary$status == "EXPLORATORY_SCREEN_COMPLETED"
      )) {
        "COMPLETED"
      } else {
        "INCOMPLETE"
      },
      output_dir = out,
      summary = summary,
      candidates = candidates
    ))
  }
  
  list(run = run)
})

run_phylogenetic_screen <- function(...) {
  pheno_phylogenetic_screen$run(...)
}

message(
  "Funcions carregades. Executa result <- run_phylogenetic_screen(). ",
  "No s'ha iniciat cap analisi."
)
