# ============================================================================
# phenoIMPACT | Focused review of the TWO onset-adjusted OFFSET plasticity models
# Version 1.0 | 2026-10-05
#
# INPUT: observation_metadata.rds and draw_01.rds ... draw_05.rds produced by
# diagnose_current_spatial_fits.R. These are small observation/residual caches,
# NOT fitted model objects. The original fit files are NEVER opened or changed.
#
# Questions:
# 1. Does unequal dispersion remain after separating population mean residuals?
# 2. Are extreme residuals concentrated in species, networks, periods or records
#    with weaker documented survey support?
# 3. Does calendar-lag dependence remain after removing each population's mean?
#
# This is NOT a new model, a variance correction, a refit, an influence analysis
# of coefficients, or a significance-based model-selection procedure.
# No abundance models, Moran tests, TMB objects or large matrices are loaded.
#
# Temporal screen: exact 1-, 2-, 3-year pairs, never compressed gaps. Primary
# eligibility: >=5 observed years and >=3 pairs per population at that lag.
# >=10 years is a descriptive sensitivity; neither changes a fitted data set.
# Raw and population-centred correlations use the SAME eligible pairs.
# Pair-weighted and equal-population-weighted results are both retained.
# 199 within-population permutations for draw 1 provide REFERENCE limits, not
# confidence intervals. They retain population means, dispersion, sample size
# and observed years, and reproduce the effects of centring short series.
# The reference assumes exchangeable years within populations; it does NOT
# calibrate all model-estimation, cross-population or nonstationarity effects.
# No p-values, FDR, pooled-draw tests or automatic model-change decisions.
#
# Draw 1 is primary; all five existing one-sample draws are reviewed separately.
# Species re-scaling is descriptive, estimated on the same data, never used to
# update coefficients, uncertainty, plasticity estimates or model likelihoods.
#
# RStudio:
#   source(file.choose(), encoding = "UTF-8")
#   offset_review_job <- pheno_offset_review$run(monitor = FALSE)
#   pheno_offset_review$status()
#
# Requires callr only for background execution. Optional data.table reads the
# survey-support CSV efficiently; optional zip packages the small review files.
# No packages/fonts are installed. Synthetic self-tests run before the review.
# ============================================================================

pheno_offset_review <- local({
  VERSION <- "offset_residual_followup_1.0"
  ROOT <- "E:/phenoIMPACT project/code/phenoIMPACT"
  DIAGNOSTIC_RUN <- "run_20261005_103726_44136"
  OFFSET_RUN <- "run_3911d7bc87ac"
  MODELS <- c("offset_univoltine", "offset_multivoltine")
  state <- new.env(parent = emptyenv())

  need <- function(ok, text) if (!isTRUE(ok)) stop(text, call. = FALSE)
  write_csv <- function(x, path) utils::write.csv(x, path, row.names = FALSE, na = "")
  atomic <- function(x, path) {
    # Only derived files in this NEW review folder are ever written.
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    tmp <- tempfile(".writing_", tmpdir = dirname(path), fileext = ".rds")
    on.exit(unlink(tmp), add = TRUE)
    saveRDS(x, tmp, compress = "gzip")
    bak <- paste0(path, ".previous")
    had <- file.exists(path)
    if (had) {
      need(file.copy(path, bak, overwrite = TRUE), paste("Cannot back up", path))
      need(file.remove(path), paste("Cannot replace derived file", path))
    }
    if (!file.rename(tmp, path)) {
      if (had) file.copy(bak, path, overwrite = TRUE)
      stop("Cannot finalize derived file: ", path, call. = FALSE)
    }
    if (had) unlink(bak)
    invisible(path)
  }
  phase <- function(cfg, stage, model = "", detail = "", finished = FALSE) {
    rec <- list(stage = stage, model = model, detail = detail,
                started = cfg$launched, time = Sys.time(), finished = finished)
    atomic(rec, file.path(cfg$out, "progress.rds"))
    message("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ",
            stage, " | ", model, " | ", detail)
  }
  finite_num <- function(x) is.numeric(x) && all(is.finite(x))
  key_id <- function(frame) {
    cols <- c("SPECIES", "SITE_ID")
    need(all(cols %in% names(frame)), "Missing population identifiers.")
    vals <- lapply(frame[cols], as.character)
    need(!anyNA(vals[[1]]) && !anyNA(vals[[2]]) &&
           !any(grepl("\r", vals[[1]], fixed = TRUE)) &&
           !any(grepl("\r", vals[[2]], fixed = TRUE)), "Invalid population identifiers.")
    paste(vals[[1]], vals[[2]], sep = "\r")
  }
  # Return sums for EVERY integer group 1..ng, including absent groups as zero.
  group_sums <- function(x, g, ng) {
    x <- as.matrix(x)
    if (!length(g)) return(matrix(0, ng, ncol(x), dimnames = list(NULL, colnames(x))))
    z <- rowsum(x, g, reorder = TRUE)
    ans <- matrix(0, ng, ncol(x), dimnames = list(NULL, colnames(x)))
    ans[as.integer(rownames(z)), ] <- z
    ans
  }
  moments <- function(x) {
    n <- length(x); mu <- mean(x); v <- mean((x-mu)^2)
    data.frame(n = n, mean = mu, SD = if(n>1) stats::sd(x) else NA_real_,
      MAD_normal_scaled = stats::mad(x), q025 = unname(stats::quantile(x,.025)),
      median = stats::median(x), q975 = unname(stats::quantile(x,.975)),
      percent_abs_z_gt_3 = 100*mean(abs(x)>3),
      percent_z_below_minus3 = 100*mean(x < -3), percent_z_above_3 = 100*mean(x > 3),
      skewness = if(v>0) mean((x-mu)^3)/v^1.5 else NA_real_,
      excess_kurtosis = if(v>0) mean((x-mu)^4)/v^2-3 else NA_real_)
  }
  concentration <- function(z, labels, label_name) {
    idx <- split(seq_along(z), as.character(labels), drop = TRUE)
    nt <- sum(abs(z)>3)
    out <- do.call(rbind, lapply(names(idx), function(k) {
      i <- idx[[k]]; a <- moments(z[i]); ne <- sum(abs(z[i])>3)
      cbind(data.frame(level = k), a, n_extreme = ne,
        observation_share = length(i)/length(z),
        extreme_share = if(nt>0) ne/nt else NA_real_,
        tail_enrichment = if(nt>0) (ne/length(i))/(nt/length(z)) else NA_real_)
    }))
    names(out)[1] <- label_name; rownames(out) <- NULL
    out
  }
  population_stats <- function(z, pop, np) {
    n <- tabulate(pop, np)
    sm <- as.numeric(group_sums(z,pop,np)); mu <- sm/n
    w <- z-mu[pop]
    ss <- as.numeric(group_sums(w*w,pop,np))
    list(n=n, mean=mu, within=w, SS=ss,
         SD=ifelse(n>1,sqrt(ss/pmax(n-1,1)),NA_real_))
  }
  species_stats <- function(z, frame, pop, ps, sigma) {
    idx <- split(seq_along(z), as.character(frame$SPECIES), drop=TRUE)
    do.call(rbind, lapply(names(idx), function(sp) {
      i <- idx[[sp]]; pp <- unique(pop[i]); s <- moments(z[i])
      dfw <- sum(ps$n[pp]-1L); df5 <- sum(ps$n[pp][ps$n[pp]>=5L]-1L)
      ss_total <- sum((z[i]-mean(z[i]))^2)
      ssw <- sum(ps$SS[pp]); between <- sum(ps$n[pp]*(ps$mean[pp]-mean(z[i]))^2)
      need(abs(ss_total-ssw-between) < 1e-7*max(1,ss_total),
           "Population sum-of-squares decomposition failed.")
      sdc <- if(dfw>0) sqrt(ssw/dfw) else NA_real_
      sd5 <- if(df5>0) sqrt(sum(ps$SS[pp][ps$n[pp]>=5L])/df5) else NA_real_
      sdz <- stats::sd(z[i]); rescaled <- if(length(i)>1L && is.finite(sdz) && sdz>1e-12)
        (z[i]-mean(z[i]))/sdz else rep(NA_real_,length(i))
      cbind(data.frame(SPECIES=sp, n_populations=length(pp),
        n_populations_ge5=sum(ps$n[pp]>=5L), n_networks=length(unique(frame$bms_id[i])),
        n_calendar_years=length(unique(frame$YEAR[i]))), s,
        within_population_SD=sdc, within_population_SD_ge5=sd5,
        residual_SD_days=s$SD*sigma, within_population_SD_days=sdc*sigma,
        # Descriptive sums of squares; NOT variance explained by the fitted model.
        between_population_SS_share=if(ss_total>0) between/ss_total else NA_real_,
        species_rescaled_percent_abs_gt3=if(all(is.finite(rescaled)))
          100*mean(abs(rescaled)>3) else NA_real_,
        stable_dispersion_support=length(i)>=100L && length(pp)>=5L)
    }))
  }
  # Weighted Pearson correlations for all pairs AND for each species. No model.
  correlations <- function(x,y,w,sp,ns) {
    if (!length(x)) return(rep(NA_real_,ns+1L))
    a <- cbind(W=w,X=w*x,Y=w*y,XX=w*x*x,YY=w*y*y,XY=w*x*y,N=rep(1,length(x)))
    s <- rbind(colSums(a),group_sums(a,sp,ns))
    xx <- s[,"XX"]-s[,"X"]^2/pmax(s[,"W"],.Machine$double.eps)
    yy <- s[,"YY"]-s[,"Y"]^2/pmax(s[,"W"],.Machine$double.eps)
    xy <- s[,"XY"]-s[,"X"]*s[,"Y"]/pmax(s[,"W"],.Machine$double.eps)
    r <- rep(NA_real_,nrow(s)); ok <- s[,"N"]>=3 & xx>1e-14 & yy>1e-14
    r[ok] <- pmax(-1,pmin(1,xy[ok]/sqrt(xx[ok]*yy[ok])))
    r
  }
  pair_design <- function(frame,pop,pop_n,species,thresholds,lags,min_pairs) {
    np <- length(pop_n); ns <- max(species)
    keys <- paste(pop,frame$YEAR,sep=":")
    need(!anyDuplicated(keys), "More than one observation in a population/year.")
    designs <- list()
    for(lag in lags) {
      earlier <- match(paste(pop,frame$YEAR-lag,sep=":"),keys)
      current <- which(!is.na(earlier)); previous <- earlier[current]
      need(all(frame$YEAR[current]-frame$YEAR[previous]==lag) &&
             all(pop[current]==pop[previous]), "Calendar pair alignment failed.")
      counts <- tabulate(pop[current],np)
      for(th in thresholds) {
        eligible <- pop_n>=th & counts>=min_pairs
        ii <- which(eligible[pop[current]])
        cur <- current[ii]; pre <- previous[ii]; sp <- species[cur]
        pop_pairs <- tabulate(pop[cur],np)
        distinct <- if(length(cur)) !duplicated(pop[cur]) else logical()
        nsppop <- if(length(cur)) tabulate(sp[distinct],ns) else integer(ns)
        supp <- data.frame(stratum=c("ALL_SPECIES",levels(factor(as.character(frame$SPECIES)))),
          lag_years=lag,min_years=th,min_pairs_per_population=min_pairs,
          n_pairs=c(length(cur),tabulate(sp,ns)),
          n_populations=c(sum(pop_pairs>0),nsppop),
          weighting=NA_character_)
        # The factor levels below must be in the same order as species integers.
        need(nrow(supp)==ns+1L,"Unexpected species support dimensions.")
        weights <- list(pair_weighted=rep(1,length(cur)),
                        equal_population=if(length(cur)) 1/pop_pairs[pop[cur]] else numeric())
        for(wname in names(weights)) {
          ss <- supp; ss$weighting <- wname
          name <- paste(lag,th,wname,sep="__")
          designs[[name]] <- list(current=cur,previous=pre,species=sp,
              weight=weights[[wname]],support=ss,ns=ns)
        }
      }
    }
    designs
  }
  temporal_estimates <- function(z,within,designs) {
    do.call(rbind,lapply(designs,function(a) {
      cbind(a$support,
        r_raw=correlations(z[a$previous],z[a$current],a$weight,a$species,a$ns),
        r_within_population=correlations(within[a$previous],within[a$current],
                                         a$weight,a$species,a$ns))
    }))
  }
  permutation_reference <- function(within,pop,designs,cfg,model) {
    # Input rows are sorted by population. Random ordering WITHIN each population
    # preserves the complete residual distribution and the calendar-year pattern.
    # Do NOT compare centred small-sample correlations to an assumed zero null.
    set.seed(cfg$seed + match(model,MODELS)*1000L)
    mats <- lapply(designs,function(a) matrix(NA_real_,a$ns+1L,cfg$nperm))
    for(b in seq_len(cfg$nperm)) {
      ord <- order(pop,stats::runif(length(pop)),method="radix")
      if(b==1L) need(identical(pop[ord],pop), "Permutation crossed population blocks.")
      perm <- within[ord]
      for(k in seq_along(designs)) {
        a<-designs[[k]]
        mats[[k]][,b]<-correlations(perm[a$previous],perm[a$current],
                                   a$weight,a$species,a$ns)
      }
      if(b==1L || b %% 25L==0L || b==cfg$nperm)
        phase(cfg,"TEMPORAL_REFERENCE",model,paste(b,"/",cfg$nperm,"within-population shuffles; no model fitting"))
    }
    tab <- do.call(rbind,lapply(seq_along(designs),function(k) {
      a<-designs[[k]]; m<-mats[[k]]
      limits<-t(apply(m,1L,function(x) {
        x<-x[is.finite(x)]
        if(!length(x)) return(c(n=0,mean=NA,lo=NA,hi=NA))
        c(n=length(x),mean=mean(x),lo=unname(stats::quantile(x,.025)),
          hi=unname(stats::quantile(x,.975)))
      }))
      cbind(a$support,n_reference_draws=limits[,1],reference_mean=limits[,2],
        reference_q025=limits[,3],reference_q975=limits[,4])
    }))
    rownames(tab)<-NULL
    list(table=tab,draws=mats)
  }
  read_support <- function(cfg) {
    p<-cfg$support_file
    if(is.null(p) || !file.exists(p)) return(list(status="NOT_AVAILABLE",data=NULL,
      message="Survey-support file absent/not requested. No detectability or observer adjustment is claimed."))
    phase(cfg,"READING_SURVEY_SUPPORT",detail="Optional join; observation and residual caches are unchanged")
    header<-names(utils::read.csv(p,nrows=0L,check.names=FALSE))
    keys<-c("SPECIES","SITE_ID","YEAR")
    if(!all(keys %in% header)) return(list(status="NOT_JOINED_MISSING_KEYS",data=NULL,
       message="Survey-support file lacks SPECIES, SITE_ID or YEAR."))
    fields<-intersect(c("n_visits_site_year","n_zero_visits_after_last_obs",
                         "last_visit_doy","last_obs_doy","first_visit_doy"),header)
    if(!length(fields)) return(list(status="NOT_JOINED_NO_SUPPORT_COLUMNS",data=NULL,
       message="No recognized original survey-support fields; no replacements were inferred."))
    take<-c(keys,fields)
    if(requireNamespace("data.table",quietly=TRUE)) {
      tab<-data.table::fread(p,select=take,colClasses=list(character=c("SPECIES","SITE_ID")),
                            data.table=FALSE,nThread=1L,showProgress=FALSE)
    } else {
      types<-rep("NULL",length(header));names(types)<-header
      types[take]<-NA_character_;types[c("SPECIES","SITE_ID")]<-"character"
      tab<-utils::read.csv(p,colClasses=types,stringsAsFactors=FALSE,check.names=FALSE)
    }
    tab<-unique(as.data.frame(tab))
    for(k in c("SPECIES","SITE_ID")) tab[[k]]<-as.character(tab[[k]])
    yr<-suppressWarnings(as.numeric(as.character(tab$YEAR)))
    if(anyNA(tab[keys]) || any(!is.finite(yr) | yr!=floor(yr)))
      return(list(status="NOT_JOINED_INVALID_KEYS",data=NULL,message="Missing/non-integer survey identifiers; review source file."))
    tab$YEAR<-as.integer(yr)
    # Conflicting key matches are NEVER averaged or multiplied into the data.
    kk<-paste(key_id(tab),tab$YEAR,sep="\r"); dup<-duplicated(kk)|duplicated(kk,fromLast=TRUE)
    if(any(dup)) {
      write_csv(tab[dup,,drop=FALSE],file.path(cfg$out,"conflicting_survey_support_keys.csv"))
      return(list(status="NOT_JOINED_AMBIGUOUS_KEYS",data=NULL,message="Conflicting species/site/year support rows exported; optional join skipped."))
    }
    list(status="AVAILABLE",data=tab,keys=kk,fields=fields,message="Original support fields available; observer and detection biases are not estimated.")
  }
  attach_support <- function(frame,y,support) {
    out<-data.frame(quality_class=rep("not_available",nrow(frame)),stringsAsFactors=FALSE)
    if(is.null(support$data)) return(out)
    hit<-match(paste(key_id(frame),frame$YEAR,sep="\r"),support$keys)
    a<-support$data[hit,support$fields,drop=FALSE];rownames(a)<-NULL
    out<-cbind(out,support_match=!is.na(hit),a)
    need(nrow(out)==nrow(frame),"Optional survey join changed row count.")
    if(all(c("n_zero_visits_after_last_obs","last_visit_doy") %in% names(a)) &&
       is.numeric(a$n_zero_visits_after_last_obs) && is.numeric(a$last_visit_doy)) {
      ok<-is.finite(a$n_zero_visits_after_last_obs)&is.finite(a$last_visit_doy)
      out$quality_class<-"unknown"
      out$quality_class[ok]<-ifelse(a$n_zero_visits_after_last_obs[ok]>=1 & y[ok]<=a$last_visit_doy[ok],
        "supported_1zero_within_visit_boundary","not_supported_1zero_or_outside_boundary")
    } else out$quality_class<-ifelse(out$support_match,"matched_but_flag_unavailable","unmatched")
    out
  }
  make_plots <- function(out,species,temporal,primary) {
    png<-function(name,expr) {
      grDevices::png(file.path(out,name),width=1800,height=1300,res=180)
      on.exit(grDevices::dev.off(),add=TRUE)
      graphics::par(mar=c(5,5,2,1)); force(expr)
    }
    s<-species[species$draw==1L & species$stable_dispersion_support,,drop=FALSE]
    if(nrow(s)) png("dispersion_before_vs_within_population.png",{
      lim<-range(c(s$SD,s$within_population_SD),finite=TRUE)
      graphics::plot(s$SD,s$within_population_SD,pch=16,cex=.7,xlim=lim,ylim=lim,
        xlab="Residual SD within species (before population centring)",
        ylab="Pooled within-population SD")
      graphics::abline(a=0,b=1,lty=2)
    })
    t<-temporal[temporal$draw==1L & temporal$lag_years==1L & temporal$min_years==5L &
        temporal$weighting=="pair_weighted" & temporal$stratum!="ALL_SPECIES" &
        temporal$n_pairs>=100L & temporal$n_populations>=5L,,drop=FALSE]
    if(nrow(t)) png("temporal_before_vs_population_centred.png",{
      lim<-range(c(t$r_raw,t$r_within_population),finite=TRUE)
      graphics::plot(t$r_raw,t$r_within_population,pch=16,cex=.7,xlim=lim,ylim=lim,
        xlab="Lag-1 pooled residual correlation (same eligible pairs)",
        ylab="Lag-1 correlation after population centring")
      graphics::abline(a=0,b=1,lty=2)
    })
    a<-primary[primary$lag_years==1 & primary$min_years==5 &
        primary$weighting=="pair_weighted" & primary$stratum!="ALL_SPECIES" &
        primary$n_pairs>=100 & primary$n_populations>=5 &
        is.finite(primary$r_minus_reference_mean),,drop=FALSE]
    if(nrow(a)) png("temporal_lag1_with_shuffle_reference.png",{
      a<-head(a[order(abs(a$r_minus_reference_mean),decreasing=TRUE),],12)
      a<-a[rev(seq_len(nrow(a))),];yy<-seq_len(nrow(a))
      graphics::par(mar=c(5,11,2,1))
      lim<-range(c(a$r_within_population,a$reference_q025,a$reference_q975),finite=TRUE)
      graphics::plot(a$r_within_population,yy,pch=16,yaxt="n",ylab="",xlim=lim,
        xlab="Centred lag-1 r: point observed; line = permutation reference, NOT CI")
      graphics::segments(a$reference_q025,yy,a$reference_q975,yy,lty=2)
      graphics::points(a$reference_mean,yy,pch=3)
      graphics::axis(2,at=yy,labels=a$stratum,las=1,cex.axis=.65)
    })
    invisible(NULL)
  }
  model_review <- function(model,cfg,support) {
    src<-file.path(cfg$diagnostic_dir,model);out<-file.path(cfg$out,model)
    dir.create(out,showWarnings=FALSE)
    input_names<-c("complete.rds","observation_metadata.rds",sprintf("draw_%02d.rds",1:5))
    paths<-file.path(src,input_names)
    need(all(file.exists(paths)),paste("Missing diagnostic cache(s), not model fits:\n",
      paste(paths[!file.exists(paths)],collapse="\n"),
      "\nUse the original diagnostic folder on E:, not the review ZIP. Nothing will be refitted."))
    phase(cfg,"CHECKING_RESIDUAL_CACHES",model,"Hashing observation metadata and five saved residual vectors, not model files")
    md5<-unname(tools::md5sum(paths));need(!anyNA(md5),"Cannot hash diagnostic caches.")
    stamp<-list(version=VERSION,paths=paths,md5=md5,nperm=cfg$nperm,seed=cfg$seed,
      min_years=cfg$thresholds,min_pairs=cfg$min_pairs,lags=cfg$lags,support_md5=cfg$support_md5)
    donepath<-file.path(out,"complete.rds")
    if(file.exists(donepath)) {
      done<-readRDS(donepath);need(identical(done$stamp,stamp),"Input/settings changed; use a new review output.")
      return("REVIEW_CACHED")
    }
    meta<-readRDS(paths[2]);old_done<-readRDS(paths[1])
    need(identical(meta$signature,old_done$signature) && identical(meta$family,"gaussian"),
         "Caches do not describe the completed Gaussian OFFSET diagnostics.")
    needed<-c("SPECIES","SITE_ID","YEAR","bms_id","ONSET_mean","ONSET_mean_z")
    frame<-as.data.frame(meta$frame)
    need(all(needed %in% names(frame)) && !anyNA(frame[needed]),"Missing annual onset or identifiers in metadata.")
    n<-nrow(frame)
    need(length(meta$y)==n && length(meta$eta)==n && finite_num(meta$y) &&
      finite_num(meta$eta) && finite_num(meta$sigma) && length(meta$sigma)==1L && meta$sigma>0,
      "Invalid Gaussian observation metadata.")
    need(finite_num(frame$YEAR) && all(frame$YEAR==floor(frame$YEAR)),"Non-integer calendar years.")
    pk<-key_id(frame)
    ord<-order(pk,frame$YEAR,method="radix"); frame<-frame[ord,,drop=FALSE]
    y<-meta$y[ord];eta<-meta$eta[ord];pk<-pk[ord]
    pop<-match(pk,unique(pk));np<-max(pop)
    sf<-factor(as.character(frame$SPECIES));species<-as.integer(sf);ns<-nlevels(sf)
    need(!anyDuplicated(paste(pop,frame$YEAR,sep=":")),"Duplicate population/calendar-year rows.")
    pop_n<-tabulate(pop,np);designs<-pair_design(frame,pop,pop_n,species,cfg$thresholds,cfg$lags,cfg$min_pairs)
    opt_support<-attach_support(frame,y,support)
    write_csv(data.frame(file=paths,md5=md5),file.path(out,"residual_input_manifest.csv"))
    write_csv(data.frame(model=model,n_observations=n,n_populations=np,n_species=ns,
      n_populations_lt5=sum(pop_n<5),n_populations_ge5=sum(pop_n>=5),
      n_populations_ge10=sum(pop_n>=10),primary_draw=1,all_draws_retained=5,
      Gaussian_SD_days=meta$sigma),file.path(out,"input_checks.csv"))
    spall<-tmall<-global<-quality<-list(); primary_z<-primary_w<-NULL
    for(draw in 1:5) {
      phase(cfg,"REVIEWING_RESIDUAL_DRAW",model,paste(draw,"/ 5; dispersion, tails and exact calendar lags"))
      rec<-readRDS(paths[draw+2L])
      need(identical(rec$signature,meta$signature) && isTRUE(rec$draw==draw) &&
        length(rec$z)==n && finite_num(rec$z),"Residual draw and metadata are incompatible.")
      z<-rec$z[ord];ps<-population_stats(z,pop,np);within<-ps$within
      s<-species_stats(z,frame,pop,ps,meta$sigma)
      s<-cbind(model=model,draw=draw,s);spall[[draw]]<-s
      tm<-temporal_estimates(z,within,designs)
      tmall[[draw]]<-cbind(model=model,draw=draw,tm)
      mu_s<-as.numeric(group_sums(z,species,ns))/tabulate(species,ns)
      ss_s<-as.numeric(group_sums((z-mu_s[species])^2,species,ns))
      sd_s<-sqrt(ss_s/pmax(tabulate(species,ns)-1,1))
      ok_s<-tabulate(species,ns)>=100L & sd_s>1e-12
      rz<-(z-mu_s[species])/sd_s[species]
      gr<-cbind(model=model,draw=draw,scale="original_one_sample_z",moments(z))
      # Same supported-species rows on BOTH sides of the re-scaling comparison.
      ii<-which(ok_s[species]);
      if(length(ii)) gr<-rbind(gr,
        cbind(model=model,draw=draw,scale="original_z_species_n_ge100",moments(z[ii])),
        cbind(model=model,draw=draw,scale="empirically_species_rescaled_DESCRIPTIVE",moments(rz[ii])))
      global[[draw]]<-gr
      quality[[draw]]<-cbind(model=model,draw=draw,
          concentration(z,opt_support$quality_class,"quality_class"))
      if(draw==1L) {
        primary_z<-z;primary_w<-within
        first<-match(seq_len(np),pop)
        pd<-data.frame(SPECIES=as.character(frame$SPECIES[first]),SITE_ID=as.character(frame$SITE_ID[first]),
          n_years=ps$n,residual_mean=ps$mean,residual_SD=ps$SD,
          first_year=as.integer(tapply(frame$YEAR,pop,min)),last_year=as.integer(tapply(frame$YEAR,pop,max)))
        write_csv(pd,file.path(out,"population_summary_primary.csv"))
        for(name in c("SPECIES","bms_id")) write_csv(concentration(z,frame[[name]],name),
          file.path(out,paste0("tail_concentration_",name,".csv")))
        decade<-paste0(floor(frame$YEAR/10)*10,"s")
        write_csv(concentration(z,decade,"decade"),file.path(out,"tail_concentration_decade.csv"))
        comb<-paste(frame$SPECIES,frame$bms_id,sep=" | ")
        write_csv(concentration(z,comb,"species_network"),file.path(out,"tail_concentration_species_network.csv"))
        audit<-order(abs(z),decreasing=TRUE)[seq_len(min(200L,n))]
        write_csv(cbind(frame[audit,,drop=FALSE],OFFSET_mean=y[audit],
          fitted_OFFSET_at_mode=eta[audit],one_sample_z=z[audit],
          residual_days=z[audit]*meta$sigma,opt_support[audit,,drop=FALSE]),
          file.path(out,"largest_residuals_for_inspection.csv"))
      }
      rm(rec,z,ps,within);invisible(gc())
    }
    sp<-do.call(rbind,spall);tm<-do.call(rbind,tmall);gs<-do.call(rbind,global)
    write_csv(sp,file.path(out,"species_dispersion_all_draws.csv"))
    write_csv(tm,file.path(out,"temporal_all_draws.csv"))
    write_csv(gs,file.path(out,"tail_rescaling_all_draws.csv"))
    write_csv(do.call(rbind,quality),file.path(out,"sampling_support_tails_all_draws.csv"))
    refpath<-file.path(out,"temporal_reference_cache.rds")
    if(file.exists(refpath)) {
      rr<-readRDS(refpath);need(identical(rr$stamp,stamp),"Reference cache mismatch.")
      ref<-rr$result
    } else {
      ref<-permutation_reference(primary_w,pop,designs,cfg,model)
      atomic(list(stamp=stamp,result=ref),refpath)
    }
    keys<-c("stratum","lag_years","min_years","min_pairs_per_population","weighting","n_pairs","n_populations")
    primary<-merge(tm[tm$draw==1L,,drop=FALSE],ref$table,by=keys,all.x=TRUE,sort=FALSE)
    need(nrow(primary)==sum(tm$draw==1L),"Temporal reference join changed result count.")
    primary$r_minus_reference_mean<-primary$r_within_population-primary$reference_mean
    primary$reference_label<-"Within-population shuffle reference; NOT a model CI or formal calibrated test"
    write_csv(primary,file.path(out,"temporal_primary_with_reference.csv"))
    summary_keys<-interaction(tm$stratum,tm$lag_years,tm$min_years,tm$weighting,drop=TRUE)
    sens<-do.call(rbind,lapply(split(seq_len(nrow(tm)),summary_keys),function(ii) {
      a<-tm[ii,,drop=FALSE];v<-a$r_within_population;v<-v[is.finite(v)]
      cbind(a[1L,c("model",keys),drop=FALSE],n_draws=nrow(a),
        minimum_r=if(length(v)) min(v) else NA_real_,median_r=if(length(v)) stats::median(v) else NA_real_,
        maximum_r=if(length(v)) max(v) else NA_real_)
    }))
    write_csv(sens,file.path(out,"temporal_draw_sensitivity.csv"))
    ss<-sp[sp$draw==1L & sp$stable_dispersion_support,,drop=FALSE]
    range_summary<-function(v) {
      v<-v[is.finite(v)&v>0]
      if(!length(v)) return(c(p10=NA,p90=NA,variance_ratio_p90_p10=NA))
      qq<-stats::quantile(v,c(.1,.9),names=FALSE)
      c(p10=qq[1],p90=qq[2],variance_ratio_p90_p10=(qq[2]/qq[1])^2)
    }
    disp<-rbind(data.frame(scale="within_species_SD",t(range_summary(ss$SD))),
                data.frame(scale="pooled_within_population_SD",t(range_summary(ss$within_population_SD))))
    write_csv(cbind(model=model,n_supported_species=nrow(ss),disp),file.path(out,"dispersion_review_summary.csv"))
    if(cfg$make_plots) tryCatch(make_plots(out,sp,tm,primary),error=function(e)
      writeLines(conditionMessage(e),file.path(out,"plot_warning.txt")))
    need(identical(unname(tools::md5sum(paths)),md5),"Diagnostic input caches changed during review.")
    atomic(list(stamp=stamp,finished=Sys.time(),status="REVIEW_EXPORTED"),donepath)
    "REVIEW_EXPORTED"
  }
  self_test <- function() {
    # Exact sums-of-squares and pair alignment; no synthetic model is fitted.
    z<-c(1,3,5,10,12,14);p<-rep(1:2,each=3)
    a<-population_stats(z,p,2)
    need(max(abs(a$mean-c(3,12)))<1e-12 && max(abs(a$within-rep(c(-2,0,2),2)))<1e-12,
         "Population-centering self-test failed.")
    need(abs(sum((z-mean(z))^2)-sum(a$SS)-sum(a$n*(a$mean-mean(z))^2))<1e-10,
         "Sum-of-squares self-test failed.")
    x<-c(1,2,4,2,6,8);y<-c(2,5,7,8,3,1);g<-rep(1:2,each=3)
    r<-correlations(x,y,rep(1,6),g,2)
    want<-c(stats::cor(x,y),stats::cor(x[1:3],y[1:3]),stats::cor(x[4:6],y[4:6]))
    need(max(abs(r-want))<1e-12,"Grouped correlation self-test failed.")
    f<-data.frame(SPECIES=rep(c("a","b"),each=3),SITE_ID=rep(c("x","y"),each=3),
                   YEAR=c(2000,2001,2003,2000,2002,2003))
    des<-pair_design(f,p,c(3,3),g,thresholds=3,lags=1:3,min_pairs=1)
    pairs<-unique(des[["1__3__pair_weighted"]]$current)
    need(identical(pairs,c(2L,6L)),"Calendar gaps were compressed in pair self-test.")
    # Centring iid observations induces negative expected distinct-pair products.
    v<-c(-1,0,1);perms<-rbind(c(1,2,3),c(1,3,2),c(2,1,3),c(2,3,1),c(3,1,2),c(3,2,1))
    ex<-mean(apply(perms,1,function(o) v[o[1]]*v[o[2]]))
    need(abs(ex + mean(v^2)/(length(v)-1))<1e-12,"Centring reference self-test failed.")
    invisible(TRUE)
  }
  readme <- function(cfg) {
    writeLines(c(
      "FOCUSED REVIEW: ONSET-ADJUSTED OFFSET PLASTICITY ONLY",VERSION,
      "No original fit, input, residual or abundance file is modified. No model is fitted.",
      "Reads only existing Gaussian observation metadata and five one-sample residual caches.",
      "Draw 1 remains primary; all five draws are retained. Draws are not independent replicates.",
      "Dispersion: within-species SD and robust MAD; pooled within-population SD after mean removal.",
      "Within-population SD divides pooled squared deviations by sum(n_population - 1).",
      "Between-population SS share is descriptive; it is NOT variance explained by plasticity or a model R2.",
      "Species re-scaling uses the same observations to estimate mean and SD. Its QQ/tail improvement",
      "would be partly mechanical; it does NOT validate a species-dispersion model or correct any CI.",
      "Stable species-dispersion summaries require >=100 rows and >=5 populations for display only.",
      "Tail concentration reports both the fraction of all observations and the fraction of extremes.",
      "Network/period associations may reflect species composition; species-network summaries are supplied.",
      "Optional survey support uses only existing recognized source fields. Conflicting joins are skipped",
      "and exported, never averaged. This does NOT estimate observer effects or detectability.",
      "The 1-zero flag requires >=1 zero after last observation and OFFSET_mean <= last_visit_doy.",
      "No observation is deleted, winsorized or imputed. Large residuals are audit candidates only.",
      "TEMPORAL: exact calendar lags 1/2/3 years. Pairs never cross populations; gaps are not compressed.",
      "Primary eligibility: >=5 observed years AND >=3 pairs at a given lag. >=10 years is sensitivity.",
      "Raw and centred correlations use exactly the same pairs. Centre over all observed years per population.",
      "Pair-weighted estimates emphasize long/well-sampled series; equal-population estimates give",
      "each eligible population a total weight of 1, but residual amplitudes can still affect correlation.",
      paste("Reference:",cfg$nperm,"within-population shuffles of draw-1 residuals, at fixed observed years."),
      "This preserves population sample sizes, means and residual distributions. It reproduces the",
      "negative short-series centring reference rather than assuming that centred r should equal zero.",
      "Reference quantiles are NOT confidence intervals, and are not corrected for scanning many species.",
      "They assume exchangeability over years within a population; nonstationarity and shared network/year",
      "effects can invalidate this approximation. No formally calibrated temporal test, p-values or FDR.",
      "Differences before/after centring support diagnostic interpretation, not proof of a specific cause.",
      "No automatic species exclusions, new AR1 model, distribution change or numerical optimization.",
      "This review does NOT establish robustness of coefficient estimates to an alternative error model.",
      "Next scientific decision: whether a targeted variance sensitivity is justified; not whether every",
      "diagnostic is perfect. A strong residual pattern alone does not dictate a model change.",
      "SOURCE/APIs:","https://sdmtmb.github.io/sdmTMB/articles/residual-checking.html",
      "https://stat.ethz.ch/R-manual/R-devel/library/base/html/order.html",
      "https://callr.r-lib.org/reference/r_bg.html"),file.path(cfg$out,"README.txt"))
  }
  collect <- function(cfg) {
    names<-c("input_checks","dispersion_review_summary","temporal_primary_with_reference",
             "temporal_draw_sensitivity","tail_rescaling_all_draws","sampling_support_tails_all_draws")
    for(nm in names) {
      rows<-lapply(MODELS,function(k) {
        f<-file.path(cfg$out,k,paste0(nm,".csv"))
        if(file.exists(f)) utils::read.csv(f,stringsAsFactors=FALSE,check.names=FALSE) else NULL
      });rows<-Filter(Negate(is.null),rows)
      if(length(rows)) write_csv(do.call(rbind,rows),file.path(cfg$out,paste0("all_",nm,".csv")))
    }
    if(requireNamespace("zip",quietly=TRUE)) tryCatch({
      f<-list.files(cfg$out,recursive=TRUE,full.names=FALSE)
      f<-f[grepl("\\.(csv|txt|png|R)$",f)]
      zip::zipr(file.path(cfg$out,"results_to_review.zip"),files=f,root=cfg$out,
                mode="mirror",include_directories=FALSE)
    },error=function(e)message("Optional ZIP failed; all tables retained: ",conditionMessage(e)))
    invisible(NULL)
  }
  controller <- function(cfg) {
    status<-data.frame(model=MODELS,status="PENDING",message="",stringsAsFactors=FALSE)
    flush<-function() write_csv(status,file.path(cfg$out,"workflow_status.csv"))
    flush();readme(cfg)
    tryCatch({
      self_test();phase(cfg,"SELF_TESTS_PASSED",detail="No model fitting, no covariance reconstruction")
      manfile<-file.path(cfg$diagnostic_dir,"source_manifest.csv")
      need(file.exists(manfile),"Original diagnostic source_manifest.csv is missing.")
      man<-utils::read.csv(manfile,stringsAsFactors=FALSE)
      need(all(c("job","role","path") %in% names(man)),"Invalid original source manifest.")
      mf<-man[man$job %in% MODELS & man$role=="fit",,drop=FALSE]
      need(nrow(mf)==2L && setequal(mf$job,MODELS) &&
        all(grepl(paste0("/offset_with_onset_15km/",OFFSET_RUN,"/"),gsub("\\\\","/",mf$path),fixed=TRUE)),
        "Diagnostic sources are not the two expected onset-adjusted offset fits.")
      write_csv(mf,file.path(cfg$out,"original_fit_provenance_NOT_LOADED.csv"))
      cfg$support_md5<-if(!is.null(cfg$support_file) && file.exists(cfg$support_file))
        unname(tools::md5sum(cfg$support_file)) else NA_character_
      if(!is.null(cfg$support_file) && file.exists(cfg$support_file))
        need(!is.na(cfg$support_md5), "Cannot fingerprint the optional survey-support file.")
      support<-read_support(cfg)
      write_csv(data.frame(file=if(is.null(cfg$support_file)) "" else cfg$support_file,
        md5=cfg$support_md5,status=support$status,message=support$message),file.path(cfg$out,"survey_support_status.csv"))
      for(i in seq_along(MODELS)) {
        status$status[i]<-"RUNNING";flush()
        ans<-tryCatch(model_review(MODELS[i],cfg,support),error=function(e) {
          status$message[i]<<-conditionMessage(e);"REVIEW_FAILED"
        })
        status$status[i]<-ans;flush();invisible(gc())
      }
      if(!is.na(cfg$support_md5)) need(identical(unname(tools::md5sum(cfg$support_file)),cfg$support_md5),
        "Survey-support file changed during review.")
      collect(cfg)
      ok<-all(status$status %in% c("REVIEW_EXPORTED","REVIEW_CACHED"))
      phase(cfg,if(ok) "REVIEW_COMPLETE" else "REVIEW_WITH_ERRORS",detail=
        if(ok) "Tables require scientific interpretation. No model was fitted or changed." else
          "Inspect workflow_status.csv; completed outputs retained.",finished=TRUE)
      atomic(list(finished=Sys.time(),ok=ok),file.path(cfg$out,"finished.rds"))
      invisible(status)
    },error=function(e) {
      writeLines(conditionMessage(e),file.path(cfg$out,"error.txt"))
      phase(cfg,"STOPPED",detail=conditionMessage(e),finished=TRUE)
      collect(cfg);stop(conditionMessage(e),call.=FALSE)
    })
  }
  engine_file <- function(out) {
    e<-environment(engine_file);n<-ls(e,all.names=TRUE)
    ff<-n[vapply(n,function(k) is.function(get(k,envir=e)),logical(1))]
    p<-file.path(out,"review_engine.R")
    dump(c("VERSION","ROOT","DIAGNOSTIC_RUN","OFFSET_RUN","MODELS",ff),file=p,envir=e,control="all")
    p
  }
  launch <- function(cfg) {
    cfg$launched<-Sys.time();atomic(cfg,file.path(cfg$out,"run_config.rds"))
    phase(cfg,"QUEUED",detail="Only residual caches; no fits or matrices will be loaded")
    log<-file.path(cfg$out,paste0("worker_",format(Sys.time(),"%Y%m%d_%H%M%S"),".log"))
    p<-callr::r_bg(function(engine,cfg) {
      e<-new.env(parent=globalenv());sys.source(engine,envir=e);e$controller(cfg)
    },args=list(cfg$engine,cfg),stdout=log,stderr="2>&1",supervise=TRUE,
      user_profile=FALSE,system_profile=FALSE,libpath=.libPaths())
    job<-list(process=p,out=cfg$out,log=log,started=cfg$launched)
    state$job<-job
    message("Offset residual review launched. Output: ",cfg$out)
    message("Keep R open. This review never refits a model.")
    invisible(job)
  }
  status <- function(job=state$job) {
    need(!is.null(job),"No offset-review job in this R session.")
    alive<-job$process$is_alive();f<-file.path(job$out,"progress.rds")
    a<-tryCatch(readRDS(f),error=function(e)NULL)
    end<-if(!alive && !is.null(a) && isTRUE(a$finished)) a$time else Sys.time()
    message(if(alive) "Elapsed: " else "Recorded elapsed: ",
      round(as.numeric(difftime(end,job$started,units="mins")),1)," min | worker running: ",alive)
    if(!is.null(a)) message(a$stage," | ",a$model," | ",a$detail)
    w<-file.path(job$out,"workflow_status.csv")
    if(file.exists(w)) try(print(utils::read.csv(w,stringsAsFactors=FALSE),row.names=FALSE),silent=TRUE)
    if(!alive && !identical(job$process$get_exit_status(),0L)) message("Worker error; see ",job$log)
    message("Output: ",job$out);invisible(a)
  }
  watch <- function(job=state$job,every=30) {
    need(!is.null(job),"No active review.")
    tryCatch(repeat {status(job);if(!job$process$is_alive())break;Sys.sleep(every)},
      interrupt=function(e)message("Display paused; the review continues."))
    invisible(job)
  }
  run <- function(root=ROOT,diagnostic_dir=file.path(root,"output","diagnostics",
      "onset_adjusted_spatial_15km",DIAGNOSTIC_RUN),
      sampling_support_file=file.path(root,"output","phenology_sampling_support_allspp.csv"),
      nperm=199L,make_plots=TRUE,monitor=FALSE) {
    need(requireNamespace("callr",quietly=TRUE),"Package callr is required; no packages were installed.")
    if(!is.null(state$job) && state$job$process$is_alive()) {
      message("This review is already running.");return(invisible(state$job))
    }
    need(dir.exists(diagnostic_dir),paste("Diagnostic cache folder not found:",diagnostic_dir))
    need(is.numeric(nperm)&&length(nperm)==1L&&is.finite(nperm)&&nperm==floor(nperm)&&nperm>=99&&nperm<=9999,
      "nperm must be an integer from 99 to 9999.")
    need(is.logical(make_plots)&&length(make_plots)==1L&&!is.na(make_plots)&&
      is.logical(monitor)&&length(monitor)==1L&&!is.na(monitor),"Invalid logical setting.")
    d<-normalizePath(diagnostic_dir,winslash="/",mustWork=TRUE)
    # A new review subfolder; the diagnostic caches remain immutable inputs.
    out<-file.path(d,"offset_followup",paste0("run_",format(Sys.time(),"%Y%m%d_%H%M%S"),"_",Sys.getpid()))
    need(!dir.exists(out),"Review folder exists; resume it or wait one second before a fresh launch.")
    dir.create(out,recursive=TRUE,showWarnings=FALSE);need(dir.exists(out),"Cannot create review folder.")
    cfg<-list(out=out,diagnostic_dir=d,support_file=sampling_support_file,nperm=as.integer(nperm),
      make_plots=make_plots,thresholds=c(5L,10L),lags=1:3,min_pairs=3L,seed=20261005L)
    cfg$engine<-engine_file(out);job<-launch(cfg)
    if(monitor)watch(job);invisible(job)
  }
  resume <- function(out=if(!is.null(state$job)) state$job$out else NULL,monitor=FALSE) {
    need(!is.null(out),"Specify the existing review output folder to resume.")
    if(!is.null(state$job)&&state$job$process$is_alive()) return(invisible(state$job))
    need(requireNamespace("callr",quietly=TRUE),"Package callr is required.")
    cfg<-readRDS(file.path(out,"run_config.rds"));need(file.exists(cfg$engine),"Review engine missing.")
    job<-launch(cfg);if(monitor)watch(job);invisible(job)
  }
  stop_job <- function(job=state$job) {
    need(!is.null(job),"No active review.")
    if(job$process$is_alive()) job$process$kill()
    message("Review stopped; original residuals and models are untouched. Completed reviews can be resumed.")
    invisible(job)
  }
  list(run=run,status=status,watch=watch,resume=resume,stop=stop_job,self_test=self_test)
})
message("Functions loaded. Start with pheno_offset_review$run(monitor = FALSE).")
