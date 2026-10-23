# Supplementary environmental interactions: within vs between species.
# Saved final sources only, one model in memory. No refitting or likelihood calls.
# source("supplement_environment_mundlak.R")
# env_supp <- supplement_environment_mundlak()
supplement_environment_mundlak <- function(
  project_root = "E:/phenoIMPACT project/code/phenoIMPACT",
  environment_run = file.path(project_root, "output", "phenology_plasticity",
    "environment_mundlak_by_voltinism_qc_final_15km", "run_de0890c85d4a"),
  font_family = "Arial", width_mm = 180, height_mm = 140, dpi = 300,
  show_plot = interactive()) {
  need <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)
  pkgs <- c("ggplot2", "patchwork", "lme4", "writexl")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  need(!length(missing), paste("Required packages:", paste(missing, collapse = ", ")))
  need(capabilities("cairo"), "Cairo is required for the vector PDF.")
  selector <- file.path(project_root, "post_environment_mundlak.R")
  need(file.exists(selector), "Local post_environment_mundlak.R not found.")
  e <- new.env(parent = globalenv())
  sys.source(selector, envir = e)
  pf <- e$pheno_post_env$preflight(project_root = project_root,
                                  environment_run = environment_run)
  need(!isTRUE(pf$active), "An environmental worker is active; export after it finishes.")
  sources <- pf$environment
  analyses <- c("onset_univoltine", "onset_multivoltine",
                "offset_univoltine", "offset_multivoltine")
  need(all(analyses %in% sources$analysis), "Incomplete source selection.")
  sources <- sources[match(analyses, sources$analysis), , drop = FALSE]
  need(all(sources$ready), "Some selected environmental sources are not accepted.")
  print(sources[c("analysis", "kind", "refinement_version", "status", "fit")], row.names = FALSE)
  mods <- c("photoperiod", "background", "predictability", "trend")
  labels <- c("Photoperiod", "Background temperature", "Temperature predictability", "Warming trend")
  colours <- setNames(c("#1b9e77", "#E6AB02", "#7570b3", "#d95f02"), labels)
  rows <- all_fixed <- contrasts <- list()

  for (i in seq_along(analyses)) {
    a <- analyses[i]
    message("Reading selected source: ", a, " (one model in memory)")
    p <- readRDS(sources$prepared[i])
    saved <- readRDS(sources$fit[i])
    fit <- if (inherits(saved, "sdmTMB")) saved else saved$fit
    need(inherits(fit, "sdmTMB"), paste("Invalid fit:", a))
    if (!is.null(saved$analysis)) need(identical(saved$analysis, a), "Analysis identity mismatch.")
    original <- readRDS(sources$original_marker[i])
    signature <- if (sources$kind[i] == "refined") saved$source_signature else saved$signature
    need(!is.null(signature) && identical(signature, original$signature),
         paste("Source provenance mismatch:", a))
    form <- p$formula
    need(nrow(fit$data) == nrow(p$data), "Fit/prepared row counts differ.")
    vars <- unique(all.vars(lme4::nobars(form)))
    for (v in vars) need(isTRUE(all.equal(fit$data[[v]], p$data[[v]],
                                         check.attributes = FALSE)),
                         paste("Fit/prepared data mismatch:", v))
    terms <- colnames(model.matrix(lme4::nobars(form), p$data[1, , drop = FALSE]))
    ix <- which(names(fit$sd_report$par.fixed) == "b_j")
    need(length(ix) == length(terms), "Cannot align fixed coefficients.")
    need(identical(names(fit$model$par), names(fit$sd_report$par.fixed)) &&
         isTRUE(all.equal(as.numeric(fit$model$par),
                          as.numeric(fit$sd_report$par.fixed), tolerance = 1e-6)),
         "Covariance and point estimates are not synchronized.")
    beta <- setNames(as.numeric(fit$sd_report$par.fixed[ix]), terms)
    V <- as.matrix(fit$sd_report$cov.fixed[ix, ix, drop = FALSE])
    dimnames(V) <- list(terms, terms)
    need(all(is.finite(beta)) && all(is.finite(V)) && all(diag(V) > 0),
         "Invalid fixed-effect covariance.")
    need(min(eigen((V + t(V))/2, symmetric = TRUE, only.values = TRUE)$values) >
           -1e-8 * max(diag(V)), "Fixed covariance is not positive semidefinite.")
    anomaly <- if (startsWith(a, "onset_")) "clim_anomaly_tw60" else "clim_anomaly_tw90"
    interaction <- function(v) {
      hit <- intersect(c(paste(anomaly, v, sep = ":"), paste(v, anomaly, sep = ":")), terms)
      need(length(hit) == 1, paste("Missing/ambiguous interaction:", v))
      hit
    }
    envvars <- as.vector(outer(mods, c("within_species", "species_mean"), paste, sep = "_"))
    need(all(envvars %in% names(p$data)), "Missing environmental components.")
    pop <- unique(p$data[c("SPECIES", "SITE_ID", envvars)])
    need(!anyDuplicated(pop[c("SPECIES", "SITE_ID")]),
         "Environmental components vary within a population.")

    all_fixed[[a]] <- data.frame(analysis = a, term = terms,
      estimate = as.numeric(beta), SE = sqrt(diag(V)),
      CI_low = beta - 1.96*sqrt(diag(V)), CI_high = beta + 1.96*sqrt(diag(V)),
      p_Wald = 2*pnorm(-abs(beta/sqrt(diag(V)))))
    for (j in seq_along(mods)) {
      for (component in c("within_species", "species_mean")) {
        v <- paste(mods[j], component, sep = "_"); term <- interaction(v)
        sc <- sd(pop[[v]])
        need(is.finite(sc) && sc > 0, paste("Invalid component SD:", a, v))
        b <- unname(beta[term]); se <- sqrt(V[term,term])
        rows[[length(rows)+1L]] <- data.frame(
          analysis = a, phenology = if (startsWith(a,"onset")) "Onset" else "Offset",
          voltinism = if (grepl("_univoltine$",a)) "Univoltine" else "Multivoltine",
          moderator = labels[j],
          component = if (component == "within_species") "Within species" else "Between species",
          term = term, estimate_original_scale = b, SE_original_scale = se,
          CI_low_original_scale = b-1.96*se, CI_high_original_scale = b+1.96*se,
          component_SD = sc, beta_per_component_SD = b*sc,
          SE_per_component_SD = se*sc, CI_low = (b-1.96*se)*sc,
          CI_high = (b+1.96*se)*sc, z_Wald = b/se,
          p_Wald = 2*pnorm(-abs(b/se)), n_observations = nrow(p$data),
          n_populations = nrow(pop), n_species = length(unique(pop$SPECIES)),
          source_kind = sources$kind[i], source_status = sources$status[i])
      }
      # Difference tested on the COMMON original moderator scale.
      w <- interaction(paste0(mods[j],"_within_species"))
      bt <- interaction(paste0(mods[j],"_species_mean"))
      diff <- beta[w]-beta[bt]
      se_diff <- sqrt(max(0,V[w,w]+V[bt,bt]-2*V[w,bt]))
      contrasts[[length(contrasts)+1L]] <- data.frame(analysis=a,moderator=labels[j],
        contrast="Within minus between, original common moderator scale",
        estimate=unname(diff),SE=se_diff,CI_low=unname(diff)-1.96*se_diff,
        CI_high=unname(diff)+1.96*se_diff,
        p_Wald=if(se_diff>0) 2*pnorm(-abs(diff/se_diff)) else NA_real_)
    }
    rm(p,saved,fit,V,pop); invisible(gc())
  }
  tab <- do.call(rbind,rows)
  need(nrow(tab)==32L,"Expected 32 environmental interactions.")
  tab$moderator <- factor(tab$moderator,levels=labels)
  tab$component <- factor(tab$component,levels=c("Within species","Between species"))
  # Keep all fixed effects and contrasts so the table does not hide model terms.
  fixed <- do.call(rbind,all_fixed); contrasts <- do.call(rbind,contrasts)
  caption <- paste(
    "Supplementary Fig. X | Within- and between-species environmental associations with phenological plasticity.",
    "Points and horizontal bars show interaction coefficients and pointwise 95% normal Wald confidence intervals",
    "for flight onset (a,b) and offset (c,d), in univoltine (a,c) and multivoltine (b,d) populations.",
    "Filled circles show within-species associations; open circles show between-species associations.",
    "Coefficients are expressed as the change in the temperature-anomaly response (days per inherited SD of anomaly)",
    "per SD of the corresponding environmental component, calculated across populations with equal population weight.",
    "Between-species component SDs are therefore population-weighted, not equal-species-weighted.",
    "Positive interaction coefficients shift responses towards later dates; negative coefficients towards earlier dates.",
    "They do not by themselves establish that the conditional response is a delay or an advancement.",
    "Offset models control for annual onset. Dashed vertical lines indicate zero interaction.",
    "Different component SDs mean equal plotted coefficients need not represent equal environmental changes.",
    "Tests of within-minus-between differences use the common original moderator scale (Supplementary Table X).")
  plots <- list()
  xrange <- range(c(tab$CI_low,tab$CI_high,0))
  xrange <- xrange + c(-1,1)*diff(xrange)*.05
  for(i in seq_along(analyses)) {
    z <- tab[tab$analysis==analyses[i],]
    z$y <- match(as.character(z$moderator),rev(labels)) +
      ifelse(z$component=="Within species", .12, -.12)
    plots[[i]] <- ggplot2::ggplot(z,ggplot2::aes(x=beta_per_component_SD,y=y,
      colour=moderator,shape=component)) +
      ggplot2::geom_vline(xintercept=0,linetype="dashed",colour="grey60",linewidth=.36) +
      ggplot2::geom_segment(ggplot2::aes(x=CI_low,xend=CI_high,yend=y),
        linewidth=.45,show.legend=FALSE) +
      ggplot2::geom_point(size=2.1,stroke=.7) +
      ggplot2::scale_colour_manual(values=colours,guide="none") +
      ggplot2::scale_shape_manual(values=c("Within species"=16,"Between species"=1),name=NULL) +
      ggplot2::scale_y_continuous(breaks=1:4,labels=rev(labels),limits=c(.5,4.5)) +
      ggplot2::coord_cartesian(xlim=xrange) +
      ggplot2::labs(tag=letters[i],
        title=paste(if(i<=2) "Onset" else "Offset",
          if(i%%2==1) "· Univoltine" else "· Multivoltine"),
        x=if(i>2) "Interaction coefficient (per component SD)" else NULL,y=NULL) +
      ggplot2::theme_classic(base_family=font_family,base_size=8) +
      ggplot2::theme(axis.line=ggplot2::element_line(linewidth=.36),
        axis.ticks=ggplot2::element_line(linewidth=.36),
        plot.title=ggplot2::element_text(face="bold",size=9),
        plot.tag=ggplot2::element_text(face="bold",size=8),
        legend.position="bottom",
        axis.text.y=if(i%%2==1) ggplot2::element_text() else ggplot2::element_blank(),
        axis.ticks.y=ggplot2::element_blank())
  }
  figure <- patchwork::wrap_plots(plots,ncol=2,guides="collect") &
    ggplot2::theme(legend.position="bottom")
  out <- file.path(environment_run,paste0("supplement_environment_",format(Sys.time(),"%Y%m%d_%H%M%S")))
  need(!dir.exists(out),"Output already exists; retry in one second.")
  dir.create(out,recursive=TRUE)
  csv <- function(x,n) write.csv(x,file.path(out,n),row.names=FALSE)
  csv(tab,"environmental_interactions.csv")
  csv(fixed,"all_fixed_effects_original_scale.csv")
  csv(contrasts,"within_between_contrasts.csv")
  csv(sources,"selected_sources.csv")
  guide <- data.frame(notes=c(
    "Primary table: 32 anomaly-by-environment interactions; includes within and between species.",
    "beta_per_component_SD and CI_low/CI_high correspond to plotted effects.",
    "Original-scale coefficients and SEs are retained for reproducibility.",
    "Units: days per inherited SD of thermal anomaly per SD of the environmental component.",
    "Component SD is calculated across unique populations with equal population weights.",
    "Anomaly standardization is inherited; do not assume a common degree-C unit across analyses.",
    "Wald P values are two-sided, unadjusted for multiplicity; CIs are pointwise.",
    "Within-minus-between contrasts use original common scale, not separately SD-rescaled coefficients.",
    "No significance of one component and non-significance of the other establishes a difference.",
    "All_fixed_effects contains all terms on original model scales, not only environmental interactions.",
    "Local source selection and its documented numerical-review warnings are retained.",
    "Collinearity and sensitivity limitations remain; supplementary reporting does not resolve them.",
    "No fitting, optimization, likelihood evaluation or demographic uncertainty propagation performed."))
  writexl::write_xlsx(list(Interactions=tab,All_fixed_effects=fixed,
    Within_between_contrasts=contrasts,Sources=sources,Notes=guide),
    file.path(out,"Supplementary_Table_environment.xlsx"))
  ggplot2::ggsave(file.path(out,"Supplementary_Figure_environment.pdf"),figure,
    width=width_mm,height=height_mm,units="mm",device=grDevices::cairo_pdf,bg="white")
  ggplot2::ggsave(file.path(out,"Supplementary_Figure_environment.png"),figure,
    width=width_mm,height=height_mm,units="mm",dpi=dpi,bg="white")
  writeLines(caption,file.path(out,"figure_caption.txt"))
  writeLines(guide$notes,file.path(out,"README.txt"))
  writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
  if(requireNamespace("zip",quietly=TRUE)) zip::zipr(file.path(out,"supplement_environment_review.zip"),
    list.files(out,full.names=TRUE),root=out)
  if(show_plot) print(figure)
  message("Supplementary figure and table: ",normalizePath(out,winslash="/"))
  invisible(list(figure=figure,table=tab,contrasts=contrasts,output_dir=out))
}
