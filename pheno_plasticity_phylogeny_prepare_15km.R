# Lightweight phylogenetic preparation for phenoIMPACT, 2026-09-28.
# Run from the project directory:
#   source("pheno_plasticity_phylogeny_prepare_15km.R")
#   phylogeny_result <- run_phylogeny()
# Requires ape. Reads ONLY the Newick tree and small coefficient/metadata CSVs.
# Does not load, fit or modify any phenology model or running job.
# This is a descriptive audit, NOT a formal test of phylogenetic signal.
# The exported species effects are partially pooled conditional predictions.
# Their marginal SEs are not independent sampling errors of unpooled slopes.
# References: Houslay & Wilson (2017), doi:10.1093/beheco/arx023;
# ape documentation: https://cran.r-project.org/package=ape

pheno_phylogeny_prepare <- local({
  analyses <- c("onset", "offset_univoltine", "offset_multivoltine")
  assert <- function(ok, msg) {
    if (!isTRUE(ok)) stop(msg, call. = FALSE)
  }
  write_csv <- function(x, path) {
    utils::write.csv(x, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  }
  read_csv <- function(path) {
    utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
                   fileEncoding = "UTF-8-BOM")
  }
  norm_name <- function(x) {
    tolower(trimws(gsub("[[:space:]]+", " ", gsub("_", " ", as.character(x), fixed = TRUE))))
  }
  usable_names <- function(x) {
    length(x) > 0L && !anyNA(x) && all(nzchar(norm_name(x)))
  }
  path_clean <- function(x) normalizePath(x, winslash = "/", mustWork = TRUE)

  inventory <- function(root) {
    dirs <- if (dir.exists(root)) list.dirs(root, recursive = FALSE, full.names = TRUE) else character()
    dirs <- dirs[startsWith(basename(dirs), "population_plasticity_15km_")]
    required <- c("species_random_slopes.csv", "fixed_effects.csv", "provenance.csv", "saved_fit_checks.csv")
    ready <- vapply(dirs, function(d) {
      all(vapply(analyses, function(a) all(file.exists(file.path(d, a, required))), logical(1)))
    }, logical(1))
    data.frame(extraction_directory = dirs, all_three_CSV_sets_present = ready,
               stringsAsFactors = FALSE)
  }

  read_map <- function(path, tree) {
    if (is.null(path) || !file.exists(path)) return(NULL)
    x <- read_csv(path)
    assert(all(c("SPECIES", "tip_label") %in% names(x)),
           "Name map must contain SPECIES and tip_label columns.")
    assert(usable_names(x$SPECIES) && usable_names(x$tip_label), "Name map contains empty names.")
    assert(!anyDuplicated(norm_name(x$SPECIES)), "Duplicated SPECIES in the name map.")
    assert(!anyDuplicated(norm_name(x$tip_label)), "Name map merges species into one tree tip; review it.")
    assert(all(norm_name(x$tip_label) %in% norm_name(tree$tip.label)),
           "Some name-map targets do not occur in the tree.")
    x
  }

  align_species <- function(x, tree, name_map) {
    target <- norm_name(x$SPECIES)
    method <- rep("spaces_underscores_case_only", nrow(x))
    if (!is.null(name_map)) {
      j <- match(target, norm_name(name_map$SPECIES))
      mapped <- !is.na(j)
      target[mapped] <- norm_name(name_map$tip_label[j[mapped]])
      method[mapped] <- "explicit_name_map"
    }
    j <- match(target, norm_name(tree$tip.label))
    found <- !is.na(j)
    labels <- tree$tip.label[j]
    assert(!anyDuplicated(labels[found]), "Multiple model species match one tip; no merging is allowed.")
    data.frame(SPECIES = x$SPECIES, tip_label = labels, matched = found,
               match_method = ifelse(found, method, "unmatched"), stringsAsFactors = FALSE)
  }

  tree_audit <- function(tree) {
    lengths_ok <- !is.null(tree$edge.length) && length(tree$edge.length) == nrow(tree$edge) &&
      all(is.finite(tree$edge.length)) && all(tree$edge.length >= 0)
    depths <- if (lengths_ok) ape::node.depth.edgelength(tree)[seq_along(tree$tip.label)] else NA_real_
    rooted <- ape::is.rooted(tree)
    data.frame(n_tips = length(tree$tip.label), rooted = rooted,
      branch_lengths_finite_nonnegative = lengths_ok,
      zero_length_edges = if (lengths_ok) sum(tree$edge.length == 0) else NA_integer_,
      ultrametric = if (lengths_ok) ape::is.ultrametric(tree) else NA,
      min_root_to_tip = if (lengths_ok) min(depths) else NA_real_,
      max_root_to_tip = if (lengths_ok) max(depths) else NA_real_,
      branch_length_units = "as supplied; not inferred from filename",
      stringsAsFactors = FALSE)
  }

  plot_slopes <- function(tree, x, path, use_lengths) {
    # x has exactly the same row order as tree$tip.label.
    n <- length(tree$tip.label)
    yhat <- x$sp_slope_dev
    se <- x$sp_slope_dev_marginal_SE
    lo <- yhat - stats::qnorm(0.975) * se
    hi <- yhat + stats::qnorm(0.975) * se
    colours <- ifelse(yhat < 0, "#246B91", "#AE6A1B")
    limits <- range(c(lo, hi, 0))
    if (diff(limits) == 0) limits <- limits + c(-1, 1)
    limits <- limits + c(-1, 1) * diff(limits) * 0.06
    grDevices::pdf(path, width = 15, height = max(7, 0.18 * n + 1.7),
                   useDingbats = FALSE, pointsize = 10)
    on.exit(grDevices::dev.off(), add = TRUE)
    graphics::layout(matrix(1:2, nrow = 1), widths = c(0.60, 0.40))
    graphics::par(mar = c(4.8, 0.7, 0.7, 0.7), yaxs = "i")
    ape::plot.phylo(tree, type = "phylogram", direction = "rightwards",
      use.edge.length = use_lengths, show.tip.label = TRUE, font = 3,
      cex = 0.68, y.lim = c(0.5, n + 0.5), no.margin = FALSE,
      tip.color = colours)
    pp <- get("last_plot.phylo", envir = get(".PlotPhyloEnv", envir = asNamespace("ape")))
    y <- pp$yy[seq_len(n)]
    if (use_lengths) {
      ape::axisPhylo(cex.axis = 0.75)
      graphics::mtext("Branch length (source units)", side = 1, line = 2.4, cex = 0.8)
    } else {
      graphics::mtext("Topology only", side = 1, line = 2.4, cex = 0.8)
    }
    graphics::par(mar = c(4.8, 0.7, 0.7, 1), yaxs = "i", xaxs = "i")
    graphics::plot(NA_real_, NA_real_, xlim = limits, ylim = c(0.5, n + 0.5),
                   axes = FALSE, xlab = "", ylab = "")
    graphics::abline(v = 0, col = "grey65", lty = 2)
    graphics::segments(lo, y, hi, y, col = grDevices::adjustcolor(colours, alpha.f = 0.65))
    graphics::points(yhat, y, pch = 16, col = colours, cex = 0.62)
    graphics::axis(1, cex.axis = 0.8)
    graphics::mtext("Species slope deviation", side = 1, line = 2, cex = 0.9)
    graphics::mtext("Days / fitted anomaly unit; bars: +/- 1.96 marginal SE", side = 1, line = 3.3, cex = 0.68)
    invisible(NULL)
  }

  one_analysis <- function(a, source, tree, name_map, out_root) {
    out <- file.path(out_root, a)
    dir.create(out)
    source <- file.path(source, a)
    paths <- file.path(source, c("species_random_slopes.csv", "fixed_effects.csv", "provenance.csv", "saved_fit_checks.csv"))
    x <- read_csv(paths[1])
    fixed <- read_csv(paths[2])
    provenance <- read_csv(paths[3])
    checks <- read_csv(paths[4])
    assert(nrow(provenance) == 1L && all(c("analysis", "cutoff_km", "no_refit") %in% names(provenance)),
           paste(a, "has missing or invalid provenance."))
    assert(identical(as.character(provenance$analysis), a) &&
             isTRUE(provenance$cutoff_km == 15) && isTRUE(as.logical(provenance$no_refit)),
           paste(a, "is not the expected extraction of existing 15-km models."))
    needed <- c("SPECIES", "sp_slope_dev", "sp_slope_dev_marginal_SE")
    assert(all(needed %in% names(x)), paste(a, "is missing species slope columns."))
    x <- x[, needed, drop = FALSE]
    x$SPECIES <- as.character(x$SPECIES)
    assert(usable_names(x$SPECIES) && !anyDuplicated(norm_name(x$SPECIES)),
           paste(a, "has missing or duplicated species identifiers."))
    assert(is.numeric(x$sp_slope_dev) && is.numeric(x$sp_slope_dev_marginal_SE) &&
      all(is.finite(x$sp_slope_dev)) && all(is.finite(x$sp_slope_dev_marginal_SE)) &&
      all(x$sp_slope_dev_marginal_SE >= 0), paste(a, "has invalid estimates or SEs."))
    audit <- align_species(x, tree, name_map)
    write_csv(audit, file.path(out, "species_matching.csv"))
    write_csv(audit[!audit$matched, , drop = FALSE], file.path(out, "unmatched_species.csv"))
    write_csv(provenance, file.path(out, "source_provenance.csv"))
    write_csv(checks, file.path(out, "source_saved_fit_checks.csv"))
    write_csv(data.frame(path = vapply(paths, path_clean, character(1)),
                         md5 = unname(tools::md5sum(paths))), file.path(out, "input_files.csv"))
    x$tip_label <- audit$tip_label
    x$matched <- audit$matched
    x$deviation_lower_approx95 <- x$sp_slope_dev - stats::qnorm(0.975) * x$sp_slope_dev_marginal_SE
    x$deviation_upper_approx95 <- x$sp_slope_dev + stats::qnorm(0.975) * x$sp_slope_dev_marginal_SE
    anomaly <- paste0("clim_anomaly_tw", if (a == "onset") 60L else 90L)
    assert(all(c("term", "estimate") %in% names(fixed)), "Invalid fixed-effects table.")
    b <- fixed$estimate[fixed$term == anomaly]
    assert(length(b) == 1L && is.finite(b), "Expected one fixed anomaly coefficient.")
    x$slope_at_zero_moderators_point_only <- as.numeric(b) + x$sp_slope_dev
    x$full_slope_SE_available <- FALSE
    write_csv(x, file.path(out, "species_slopes_audited.csv"))
    n_match <- sum(x$matched)
    cov_ready <- FALSE
    rooted <- ultrametric <- NA
    if (n_match >= 2L) {
      pruned <- ape::keep.tip(tree, x$tip_label[x$matched])
      pruned <- ape::ladderize(pruned)
      matched <- x[match(pruned$tip.label, x$tip_label), , drop = FALSE]
      assert(identical(pruned$tip.label, matched$tip_label), "Tree/table order mismatch.")
      ape::write.tree(pruned, file = file.path(out, "matched_tree.nwk"))
      write_csv(matched, file.path(out, "matched_species_slopes.csv"))
      tc <- tree_audit(pruned)
      rooted <- tc$rooted
      ultrametric <- tc$ultrametric
      write_csv(tc, file.path(out, "matched_tree_checks.csv"))
      plot_slopes(pruned, matched, file.path(out, "tree_and_species_slope_deviations.pdf"),
                  tc$branch_lengths_finite_nonnegative)
      # Small matrices for a later, uncertainty-aware comparative model.
      # No rerooting, ultrametric transformation, polytomy resolution, or PD repair.
      if (rooted && tc$branch_lengths_finite_nonnegative) {
        v <- ape::vcv.phylo(pruned, corr = FALSE)
        v <- v[pruned$tip.label, pruned$tip.label, drop = FALSE]
        positive <- all(is.finite(v)) && all(diag(v) > 0)
        if (positive) {
          correlation <- stats::cov2cor(v)
          cov_ready <- !inherits(try(chol(correlation), silent = TRUE), "try-error")
          if (cov_ready) {
            saveRDS(list(tip_order = pruned$tip.label, brownian_covariance = v,
                         brownian_correlation = correlation, tree = pruned,
                         note = "Prepared inputs only; no signal test was performed."),
                    file.path(out, "phylogenetic_covariance_inputs.rds"))
          }
        }
      }
    }
    fit_ok <- if ("all_sanity_checks_pass" %in% names(checks))
      all(!is.na(checks$all_sanity_checks_pass) & as.logical(checks$all_sanity_checks_pass)) else NA
    data.frame(analysis = a, status = if (n_match >= 2L) "PREPARED" else "TOO_FEW_MATCHES",
      n_species = nrow(x), n_matched = n_match, n_unmatched = nrow(x) - n_match,
      matched_percent = round(100 * n_match / nrow(x), 1),
      saved_fit_checks_all_pass = fit_ok, rooted = rooted, ultrametric = ultrametric,
      phylogenetic_covariance_ready = cov_ready, formal_signal_test_run = FALSE,
      note = if (n_match < nrow(x)) "Matched subset only; review unmatched_species.csv." else "",
      stringsAsFactors = FALSE)
  }

  bind_summary <- function(items) {
    keys <- unique(unlist(lapply(items, names), use.names = FALSE))
    items <- lapply(items, function(x) {
      for (k in setdiff(keys, names(x))) x[[k]] <- NA
      x[, keys, drop = FALSE]
    })
    do.call(rbind, items)
  }

  run <- function(project_root = "E:/phenoIMPACT project/code/phenoIMPACT",
                  tree_file = file.path(project_root, "data", "EUROPEAN_BUTTERFLIES_FULLMCC_DROPTIPED.nwk"),
                  extraction_dir = NULL,
                  name_map_file = file.path(project_root, "data", "phylogeny_species_map.csv")) {
    assert(requireNamespace("ape", quietly = TRUE), "Install ape once: install.packages('ape')")
    project_root <- path_clean(project_root)
    tree_file <- path_clean(tree_file)
    tree <- ape::read.tree(tree_file)
    assert(inherits(tree, "phylo") && !inherits(tree, "multiPhylo"), "Expected exactly one Newick tree.")
    assert(usable_names(tree$tip.label) && length(tree$tip.label) >= 2L, "Tree has invalid tip labels.")
    parent <- file.path(project_root, "output", "phenology_plasticity", "phylogeny_15km")
    dir.create(parent, recursive = TRUE, showWarnings = FALSE)
    out <- tempfile(paste0("preparation_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_"), tmpdir = parent)
    assert(dir.create(out), "Cannot create output directory.")
    out <- path_clean(out)
    message("Output: ", out)
    write_csv(data.frame(tip_label = tree$tip.label, matching_key = norm_name(tree$tip.label)),
              file.path(out, "all_tree_tips.csv"))
    assert(!anyDuplicated(norm_name(tree$tip.label)),
           paste("Tree names collide after space/underscore/case normalization. Review", out))
    write_csv(tree_audit(tree), file.path(out, "full_tree_checks.csv"))
    write_csv(data.frame(path = tree_file, md5 = unname(tools::md5sum(tree_file))),
              file.path(out, "tree_source.csv"))
    name_map <- read_map(name_map_file, tree)
    if (!is.null(name_map)) {
      write_csv(name_map, file.path(out, "name_map_used.csv"))
      write_csv(data.frame(path = path_clean(name_map_file), md5 = unname(tools::md5sum(name_map_file))),
                file.path(out, "name_map_source.csv"))
    }
    writeLines(c(
      "PHYLOGENETIC PREPARATION: onset, univoltine offset, multivoltine offset.",
      "Uses only the supplied Newick and small CSV exports from existing 15-km fits.",
      "No refit, large-model read, worker launch, or changes to existing model files.",
      "First peak is not included: the source extraction does not export it.",
      "",
      "ESTIMAND AND PLOT:",
      "The plots show species deviations in the temperature-anomaly slope, conditional on all fitted moderators.",
      "Zero is the model's species-average deviation, NOT zero total thermal plasticity.",
      "Negative deviation means a more negative slope than the fixed slope at the same moderator values.",
      "Bars are the exported deviation +/- 1.96 marginal SE; they are approximate prediction intervals for the species effect.",
      "They are NOT full intervals for total slopes, and do not represent tree or model-specification uncertainty.",
      "The point-only slope at zero moderators equals the fixed anomaly coefficient + species deviation.",
      "Moderator zero is used on the original fitted scale; it is not assumed to equal a biological mean.",
      "Adding the same fixed slope to every species does not alter their relative pattern within an analysis.",
      "Units are days per one FITTED anomaly unit. Days/degree C requires confirming the original scaling.",
      "",
      "INFERENCE:",
      "These are partially pooled model predictions (BLUP-type effects), not independent observed species traits.",
      "Marginal SEs do not supply the joint uncertainty or undo shrinkage.",
      "No K, lambda, phylogenetic heritability, or P value is calculated here.",
      "A formal analysis must model/propagate slope estimation and preferably estimate phylogenetic slope covariance jointly.",
      "This conditional comparison concerns differences remaining after the environmental moderators, not total unadjusted plasticity.",
      "A dated single MCC tree does not represent uncertainty across possible trees.",
      "Reference: Houslay & Wilson (2017), Avoiding the misuse of BLUP in behavioural ecology, doi:10.1093/beheco/arx023.",
      "",
      "MATCHING AND SOURCE DIAGNOSTICS:",
      "Only spaces, underscores and case are normalized automatically. No fuzzy taxonomy or genus-based grafting.",
      "Unmatched species are listed explicitly; plots/matrices use the matched subset.",
      "To resolve synonyms, create data/phylogeny_species_map.csv with SPECIES and tip_label columns and verified pairs only.",
      "Branches are preserved; trees are not automatically rerooted, made ultrametric, or made positive definite.",
      "Full and pruned tree diagnostics are retained; missing lengths allow a topology-only exploratory plot.",
      "Original fit warnings remain in source_saved_fit_checks.csv (including any known univoltine range warning).",
      "A successful preparation does not validate the original model fits or establish phylogenetic signal."
    ), file.path(out, "README.txt"))
    candidates <- inventory(file.path(project_root, "output", "phenology_plasticity", "spatial_models"))
    write_csv(candidates, file.path(out, "extraction_candidates.csv"))
    if (is.null(extraction_dir)) {
      complete <- candidates$extraction_directory[candidates$all_three_CSV_sets_present]
      if (length(complete) != 1L) {
        state <- if (length(complete) == 0L) "CSV_EXPORTS_NOT_FOUND" else "CHOOSE_EXTRACTION_DIRECTORY"
        result <- data.frame(status = state, tree_tips = length(tree$tip.label),
                             complete_extractions_found = length(complete), formal_signal_test_run = FALSE)
        write_csv(result, file.path(out, "preparation_summary.csv"))
        print(result, row.names = FALSE)
        if (length(complete) > 1L) {
          print(complete)
          message("Rerun with run_phylogeny(extraction_dir = 'the/chosen/directory').")
        } else {
          message("Tree audit saved. Existing species_random_slopes.csv exports are required; no large extraction was started.")
          message("If stored elsewhere, pass their parent extraction directory via extraction_dir.")
        }
        return(invisible(list(output_dir = out, summary = result, extraction_candidates = candidates)))
      }
      extraction_dir <- complete[[1]]
    }
    extraction_dir <- path_clean(extraction_dir)
    message("Using existing CSV exports: ", extraction_dir)
    writeLines(extraction_dir, file.path(out, "selected_extraction.txt"))
    results <- list()
    for (a in analyses) {
      message("Preparing ", a, " ...")
      results[[a]] <- tryCatch(one_analysis(a, extraction_dir, tree, name_map, out),
        error = function(e) data.frame(analysis = a, status = "ERROR", formal_signal_test_run = FALSE,
                                       note = conditionMessage(e), stringsAsFactors = FALSE))
      write_csv(bind_summary(results), file.path(out, "preparation_summary.csv"))
    }
    result <- bind_summary(results)
    writeLines(capture.output(utils::sessionInfo()), file.path(out, "sessionInfo.txt"))
    print(result, row.names = FALSE)
    message("Preparation finished. Output: ", out)
    message("Review the summary, unmatched species and source fit warnings before interpreting the plots.")
    invisible(list(output_dir = out, summary = result, extraction_dir = extraction_dir))
  }
  list(run = run)
})

run_phylogeny <- function(...) pheno_phylogeny_prepare$run(...)
message("Phylogeny preparation loaded. Run phylogeny_result <- run_phylogeny().")
