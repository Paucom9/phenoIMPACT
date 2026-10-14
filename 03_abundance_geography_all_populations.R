# =============================================================================
# phenoIMPACT | Independent abundance-geography test using ALL abundance pops
# Version 1.0 | 2026-10-07
#
# This is the strict sensitivity proposed after the voltinism discussion.
# It does NOT split populations by voltinism and does NOT require a plasticity
# estimate. The model uses every eligible abundance population in
# output/pheno_abund_estimates_allspp.csv.
#
# MODEL
#   ABUND_INDEX ~ year_decade * latitude_within_species_10deg
#               + year_decade * longitude_within_species_10deg
#               + species intercept + year slope
#               + site intercept + year slope
#               + population AR(1)
#               + annual IID Matern spatial field
#
# Primary coefficient:
#   year_decade:latitude_within_species_10deg
#
# Positive = populations 10 degrees farther north within the same species have
# a more positive / less negative decadal abundance trajectory, controlling
# longitude.
#
# RUN
#   source(file.choose(), encoding = "UTF-8")
#   all_geo_job <- pheno_allpop_geo$run(monitor = FALSE)
#
# MONITOR
#   pheno_allpop_geo$status()
#
# Keep R open. This is ONE large spatial TMB fit.
# =============================================================================

pheno_allpop_geo <- local({

  VERSION <- "all_population_abundance_geography_spatial_15km_1.0"
  ROOT <- "E:/phenoIMPACT project/code/phenoIMPACT"
  SOURCE_RUN <- paste0(
    ROOT,
    "/output/population_trends/spatial_abundance_onset_adjusted_15km/",
    "run_20261004_214012_44136"
  )
  GRAD_TOL <- 0.05
  state <- new.env(parent = emptyenv())

  need <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)
  csv <- function(x, p) utils::write.csv(x, p, row.names = FALSE, na = "")

  atomic_save <- function(x, p) {
    dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
    tmp <- tempfile(".writing_", tmpdir = dirname(p), fileext = ".rds")
    on.exit(unlink(tmp), add = TRUE)
    saveRDS(x, tmp, compress = "gzip")
    bak <- paste0(p, ".previous")
    if (file.exists(bak)) unlink(bak)
    had <- file.exists(p)
    if (had) need(file.rename(p, bak), paste("Cannot retain previous:", p))
    if (!file.rename(tmp, p)) {
      if (had && file.exists(bak)) file.rename(bak, p)
      stop("Cannot finalize: ", p, call. = FALSE)
    }
    if (file.exists(bak)) unlink(bak)
    invisible(p)
  }

  phase <- function(cfg, stage, detail = "", terminal = FALSE) {
    z <- list(stage=stage, detail=detail, time=Sys.time(), terminal=terminal)
    atomic_save(z, cfg$progress)
    message("[", format(z$time,"%Y-%m-%d %H:%M:%S"), "] ", stage,
            if(nzchar(detail)) paste0(" | ",detail) else "")
  }

  check_packages <- function() {
    p <- c("TMB","Matrix","fmesher","sf","callr","ps","digest")
    miss <- p[!vapply(p, requireNamespace, logical(1), quietly=TRUE)]
    need(!length(miss), paste("Missing packages:",paste(miss,collapse=", ")))
  }

  disk_check <- function(path, min_gib=25, verbose=TRUE) {
    d <- ps::ps_disk_usage(path)
    need(nrow(d)==1L && is.finite(d$available),"Cannot check free disk space.")
    free <- as.numeric(d$available)/1024^3
    if(verbose) message("Disk preflight: ",round(free,1)," GiB free.")
    need(free>=min_gib,paste("Only",round(free,1),"GiB free; require",min_gib,"GiB."))
    invisible(free)
  }

  current_fits <- function(source_run) {
    multi <- file.path(source_run,"multivoltine__full","fit.rds")
    need(file.exists(multi),"Current multivoltine fit missing.")
    rr <- file.path(source_run,"refinement_univoltine")
    need(dir.exists(rr),"Univoltine refinement directory missing.")
    dirs <- list.dirs(rr,recursive=FALSE,full.names=TRUE)
    good <- dirs[vapply(dirs,function(d) {
      r <- file.path(d,"result.rds"); f <- file.path(d,"univoltine__full","fit.rds")
      if(!all(file.exists(c(r,f))))return(FALSE)
      z <- tryCatch(readRDS(r),error=function(e)NULL)
      !is.null(z) && identical(z$status,"FIT_OK")
    },logical(1))]
    need(length(good)==1L,paste("Expected one FIT_OK univoltine refinement; found",length(good)))
    list(multivoltine=normalizePath(multi,winslash="/",mustWork=TRUE),
         univoltine=normalizePath(file.path(good,"univoltine__full","fit.rds"),winslash="/",mustWork=TRUE))
  }

  make_mesh <- function(loc, cutoff=15) {
    xy <- unique(as.matrix(loc[c("x_km","y_km")]))
    mesh <- fmesher::fm_rcdt_2d_inla(loc=xy,cutoff=cutoff)
    fem <- fmesher::fm_fem(mesh,order=2)
    list(mesh=mesh,
         M0=methods::as(fem$c0,"CsparseMatrix"),
         M1=methods::as(fem$g1,"CsparseMatrix"),
         M2=methods::as(fem$g2,"CsparseMatrix"),
         cutoff_km=cutoff)
  }

  sparse_TMB <- function(x) methods::as(methods::as(methods::as(x,"dMatrix"),"generalMatrix"),"TsparseMatrix")

  make_objective <- function(data, parameters) {
    for(nm in c("A","M0","M1","M2"))data[[nm]] <- sparse_TMB(data[[nm]])
    TMB::MakeADFun(data=data,parameters=parameters,
      random=c("b_species","b_site","ar_obs","field"),
      DLL="trend_spde",silent=TRUE,
      inner.control=list(maxit=1000,trace=FALSE))
  }

  load_dll <- function(source_run) {
    dll <- TMB::dynlib(file.path(source_run,"compiled","trend_spde"))
    need(file.exists(dll),paste("Missing source DLL:",dll))
    if(!"trend_spde" %in% names(getLoadedDLLs()))dyn.load(dll)
    TMB::openmp(n=1L,DLL="trend_spde")
    invisible(dll)
  }

  prepare_all <- function(cfg) {
    phase(cfg,"PREPARING","Reading ALL abundance populations and coordinates")

    abundance_file <- file.path(cfg$project_root,"output","pheno_abund_estimates_allspp.csv")
    coord_file <- file.path(cfg$project_root,"data","ebms_transect_coord.csv")
    need(file.exists(abundance_file),paste("Missing:",abundance_file))
    need(file.exists(coord_file),paste("Missing:",coord_file))

    a <- utils::read.csv(abundance_file,stringsAsFactors=FALSE,check.names=FALSE)
    need(all(c("SPECIES","SITE_ID","YEAR","ABUND_INDEX") %in% names(a)),
         "Unexpected abundance columns.")
    need(is.numeric(a$ABUND_INDEX),"ABUND_INDEX must be numeric.")

    # Preserve the existing abundance workflow's 99.9% global threshold.
    threshold <- as.numeric(stats::quantile(a$ABUND_INDEX,.999,na.rm=TRUE))
    keep <- !is.na(a$ABUND_INDEX) & is.finite(a$ABUND_INDEX) &
      a$ABUND_INDEX>0 & a$ABUND_INDEX<=threshold

    d <- unique(a[keep,c("SPECIES","SITE_ID","YEAR","ABUND_INDEX")])
    d$SPECIES <- as.character(d$SPECIES)
    d$SITE_ID <- as.character(d$SITE_ID)
    d$year_num <- suppressWarnings(as.integer(as.character(d$YEAR)))
    d <- d[!is.na(d$SPECIES)&nzchar(d$SPECIES)&!is.na(d$SITE_ID)&nzchar(d$SITE_ID)&!is.na(d$year_num),]
    need(nrow(d)>0,"No abundance observations retained.")
    need(!anyDuplicated(d[c("SPECIES","SITE_ID","year_num")]),
         "Conflicting abundance values for the same species/site/year.")
    d$pop_id <- paste(d$SPECIES,d$SITE_ID,sep="_")

    coords <- utils::read.csv(coord_file,stringsAsFactors=FALSE,check.names=FALSE)
    need(all(c("transect_id","bms_id","transect_lon","transect_lat")%in%names(coords)),
         "Unexpected coordinate columns.")
    sites0 <- unique(data.frame(
      SITE_ID=as.character(coords$transect_id),
      bms_id=as.character(coords$bms_id),
      x_km=as.numeric(coords$transect_lon)/1000,
      y_km=as.numeric(coords$transect_lat)/1000,
      stringsAsFactors=FALSE
    ))
    sites0 <- sites0[is.finite(sites0$x_km)&is.finite(sites0$y_km),]
    need(!anyDuplicated(sites0$SITE_ID),"SITE_ID has conflicting coordinates/BMS.")

    # Keep abundance rows with a valid coordinate.
    mm <- match(d$SITE_ID,sites0$SITE_ID)
    missing_sites <- unique(d$SITE_ID[is.na(mm)])
    need(length(missing_sites)==0,
         paste("Abundance sites lack coordinates:",paste(head(missing_sites,20),collapse=", ")))

    # Geographic degrees for within-species centring.
    sf0 <- sites0
    sf0$x_m <- sf0$x_km*1000; sf0$y_m <- sf0$y_km*1000
    sfo <- sf::st_as_sf(sf0,coords=c("x_m","y_m"),crs=3035)
    ll <- sf::st_coordinates(sf::st_transform(sfo,4326))
    sites0$longitude_deg <- ll[,1];sites0$latitude_deg <- ll[,2]

    # Population coordinates, equal population weight when centring within species.
    pop <- unique(d[c("pop_id","SPECIES","SITE_ID")])
    pop <- merge(pop,sites0,by="SITE_ID",all.x=TRUE,sort=FALSE)
    need(!anyDuplicated(pop$pop_id),"Population coordinate join duplicated rows.")
    latm <- tapply(pop$latitude_deg,pop$SPECIES,mean)
    lonm <- tapply(pop$longitude_deg,pop$SPECIES,mean)
    pop$latitude_species_mean <- unname(latm[pop$SPECIES])
    pop$longitude_species_mean <- unname(lonm[pop$SPECIES])
    pop$latitude_within_species_deg <- pop$latitude_deg-pop$latitude_species_mean
    pop$longitude_within_species_deg <- pop$longitude_deg-pop$longitude_species_mean
    pop$latitude_within_species_10deg <- pop$latitude_within_species_deg/10
    pop$longitude_within_species_10deg <- pop$longitude_within_species_deg/10
    need(max(abs(rowsum(pop$latitude_within_species_deg,pop$SPECIES,reorder=FALSE)))<1e-8,
         "Latitude centring failed.")
    need(max(abs(rowsum(pop$longitude_within_species_deg,pop$SPECIES,reorder=FALSE)))<1e-8,
         "Longitude centring failed.")

    m <- match(d$pop_id,pop$pop_id)
    for(v in c("latitude_deg","longitude_deg","latitude_species_mean","longitude_species_mean",
               "latitude_within_species_deg","longitude_within_species_deg",
               "latitude_within_species_10deg","longitude_within_species_10deg"))d[[v]] <- pop[[v]][m]

    # Centre time on this complete abundance universe.
    year_center <- mean(d$year_num)
    d$year_decade <- (d$year_num-year_center)/10

    # Stable deterministic ordering for all latent structures.
    d <- d[order(d$pop_id,d$year_num),]
    rownames(d) <- NULL
    species_levels <- sort(unique(d$SPECIES))
    site_levels <- sort(unique(d$SITE_ID))
    years <- sort(unique(d$year_num))
    d$species_index <- match(d$SPECIES,species_levels)-1L
    d$site_index <- match(d$SITE_ID,site_levels)-1L
    d$year_index <- match(d$year_num,years)-1L

    # AR1 previous observed row and true calendar gap for every population.
    previous <- rep(-1L,nrow(d)); gap <- rep(1L,nrow(d))
    spl <- split(seq_len(nrow(d)),d$pop_id)
    for(idx in spl) {
      if(length(idx)>1L) {
        previous[idx[-1L]] <- idx[-length(idx)]-1L
        gap[idx[-1L]] <- diff(d$year_num[idx])
      }
    }
    need(all(gap>=1L),"Invalid AR1 calendar gaps.")

    X <- stats::model.matrix(
      ~ year_decade*(latitude_within_species_10deg+longitude_within_species_10deg),data=d)
    expected <- c("(Intercept)","year_decade","latitude_within_species_10deg",
                  "longitude_within_species_10deg",
                  "year_decade:latitude_within_species_10deg",
                  "year_decade:longitude_within_species_10deg")
    need(identical(colnames(X),expected),"Unexpected all-population fixed design.")
    need(qr(X)$rank==ncol(X),"All-population fixed design is rank deficient.")

    sites <- sites0[match(site_levels,sites0$SITE_ID),]
    need(all(sites$SITE_ID==site_levels),"Site ordering failed.")

    phase(cfg,"MESH","Building 15-km mesh for all abundance sites")
    mesh <- make_mesh(sites,15)
    atomic_save(mesh,file.path(cfg$out,"mesh.rds"))

    A <- methods::as(fmesher::fm_basis(mesh$mesh,loc=as.matrix(sites[c("x_km","y_km")])),"CsparseMatrix")
    need(all(abs(Matrix::rowSums(A)-1)<1e-8),"Some all-population sites lie outside new mesh.")

    data <- list(
      y=as.numeric(d$ABUND_INDEX),X=unname(X),time=as.numeric(d$year_decade),
      species=as.integer(d$species_index),site=as.integer(d$site_index),
      year=as.integer(d$year_index),previous=as.integer(previous),gap=as.integer(gap),
      A=A,M0=mesh$M0,M1=mesh$M1,M2=mesh$M2,use_spatial=1L
    )

    # Hyperparameter warm starts from the two reviewed current models.
    fits <- lapply(cfg$current_fits,readRDS)
    nobs_src <- vapply(fits,function(f)f$checks$nobs,numeric(1))
    w <- nobs_src/sum(nobs_src)
    avg <- function(name) Reduce(`+`,Map(function(f,ww)as.numeric(f$parameters[[name]])*ww,fits,w))
    beta0 <- rep(0,ncol(X));names(beta0)<-colnames(X)
    beta0["(Intercept)"] <- sum(vapply(fits,function(f)unname(f$beta["(Intercept)"]),numeric(1))*w)
    beta0["year_decade"] <- sum(vapply(fits,function(f)unname(f$beta["year_decade"]),numeric(1))*w)

    parameters <- list(
      beta=unname(beta0),
      log_sd=avg("log_sd"),rho_raw=avg("rho_raw"),log_shape=avg("log_shape"),
      b_species=matrix(0,length(species_levels),2),
      b_site=matrix(0,length(site_levels),2),
      ar_obs=rep(0,nrow(d)),
      log_range=avg("log_range"),log_spatial_sd=avg("log_spatial_sd"),
      field=matrix(0,nrow(mesh$M0),length(years))
    )

    prepared <- list(
      version=VERSION,frame=d,populations=pop,sites=sites,years=years,
      species_levels=species_levels,site_levels=site_levels,X=X,data=data,parameters=parameters,
      abundance_threshold=threshold,year_center=year_center,
      signature=digest::digest(list(d[c("SPECIES","SITE_ID","year_num","ABUND_INDEX",
                                          "latitude_within_species_10deg","longitude_within_species_10deg")],X),algo="sha256")
    )
    atomic_save(prepared,file.path(cfg$out,"prepared_all_populations.rds"))
    csv(pop,file.path(cfg$out,"population_geography_all.csv"))
    csv(data.frame(
      n_observations=nrow(d),n_populations=length(unique(d$pop_id)),
      n_species=length(species_levels),n_sites=length(site_levels),n_years=length(years),
      min_year=min(years),max_year=max(years),abundance_99_9_threshold=threshold,
      mesh_vertices=nrow(mesh$M0),mesh_cutoff_km=15,
      within_lat_lon_correlation=cor(pop$latitude_within_species_deg,pop$longitude_within_species_deg),
      within_lat_lon_VIF=1/(1-cor(pop$latitude_within_species_deg,pop$longitude_within_species_deg)^2)
    ),file.path(cfg$out,"data_summary.csv"))
    invisible(TRUE)
  }

  fixed_table <- function(beta,V) {
    se <- sqrt(diag(V));z <- beta/se
    data.frame(term=names(beta),estimate=as.numeric(beta),std.error=as.numeric(se),
      conf.low=as.numeric(beta-1.96*se),conf.high=as.numeric(beta+1.96*se),
      z=as.numeric(z),p_Wald=as.numeric(2*pnorm(-abs(z))))
  }

  fit_worker <- function(cfg) {
    check_packages();disk_check(cfg$out,10,FALSE);load_dll(cfg$source_run)
    p <- readRDS(file.path(cfg$out,"prepared_all_populations.rds"))
    mesh <- readRDS(file.path(cfg$out,"mesh.rds"))
    fitfile <- file.path(cfg$out,"fit.rds")
    if(file.exists(fitfile)) {
      f <- readRDS(fitfile)
      need(identical(f$version,VERSION)&&identical(f$signature,p$signature),"Cached fit mismatch.")
      return(if(isTRUE(f$checks$numerical_checks_ok))"FIT_OK_CACHED" else "NUMERICAL_REVIEW_REQUIRED")
    }

    data <- p$data
    pars <- p$parameters
    warmfile <- file.path(cfg$out,"last_evaluation.rds")
    optfile <- file.path(cfg$out,"optimized_parameters.rds")
    warm <- if(file.exists(warmfile)) {z<-readRDS(warmfile);if(identical(z$version,VERSION)&&identical(z$signature,p$signature))z else NULL} else NULL
    optimized <- if(file.exists(optfile)) {z<-readRDS(optfile);if(identical(z$version,VERSION)&&identical(z$signature,p$signature))z else NULL} else NULL
    if(!is.null(warm))pars <- warm$parameters
    if(!is.null(optimized))pars <- optimized$parameters

    o <- make_objective(data,pars)
    on.exit(try(TMB::FreeADFun(o),silent=TRUE),add=TRUE)
    if(!is.null(warm)&&is.null(optimized))o$par <- warm$par
    if(!is.null(optimized))o$par <- optimized$optimizer$par

    warns <- character();capwarn <- function(w){warns<<-unique(c(warns,conditionMessage(w)));message("WARNING: ",conditionMessage(w));invokeRestart("muffleWarning")}
    current <- function(par){o$fn(par);x<-o$env$parList(par=o$env$last.par);need(all(vapply(x,function(v)is.numeric(v)&&all(is.finite(v)),logical(1))),"Non-finite parameter arrays.");x}

    if(is.null(optimized)) {
      nfn<-0L;ngr<-0L;last<-Sys.time()-3600
      fn <- function(par) {
        val<-o$fn(par);nfn<<-nfn+1L
        if(is.finite(val)&&(nfn==1L||as.numeric(difftime(Sys.time(),last,units="secs"))>=90)) {
          atomic_save(list(version=VERSION,signature=p$signature,par=par,parameters=current(par),objective=val,time=Sys.time()),warmfile)
          phase(cfg,"FITTING",paste("fn",nfn,"| gr",ngr,"| objective",signif(val,10),"| checkpoint saved"))
          last<<-Sys.time()
        }
        val
      }
      gr <- function(par){ngr<<-ngr+1L;o$gr(par)}
      phase(cfg,"FITTING",paste(nrow(p$frame),"annual rows;",length(unique(p$frame$pop_id)),"populations;",length(p$species_levels),"species"))
      set.seed(20261007L)
      opt <- withCallingHandlers(nlminb(o$par,fn,gr,control=list(iter.max=cfg$iter,eval.max=cfg$eval,rel.tol=1e-10)),warning=capwarn)
      if(opt$convergence!=0L) {
        phase(cfg,"REFINING",paste("nlminb code",opt$convergence,"-> BFGS"))
        second <- withCallingHandlers(optim(opt$par,fn,gr,method="BFGS",control=list(maxit=300L,reltol=1e-10)),warning=capwarn)
        if(is.finite(second$value)&&second$value<=opt$objective+1e-6)opt<-list(par=second$par,objective=second$value,convergence=second$convergence,message="BFGS refinement",first_optimizer=opt,counts=second$counts)
      }
      val<-o$fn(opt$par);need(is.finite(val)&&abs(val-opt$objective)<.01,"Final objective mismatch.")
      atomic_save(list(version=VERSION,signature=p$signature,parameters=current(opt$par),optimizer=opt),optfile)
    } else {
      phase(cfg,"RECOVERING_UNCERTAINTY","Optimizer already completed; no repeated optimization")
      opt <- optimized$optimizer
      val<-o$fn(opt$par);need(is.finite(val)&&abs(val-opt$objective)<.01,"Saved optimum replay failed.")
    }

    phase(cfg,"UNCERTAINTY","Computing Hessian/covariance")
    sd <- withCallingHandlers(TMB::sdreport(o,par.fixed=opt$par,getJointPrecision=FALSE,getReportCovariance=FALSE),warning=capwarn)
    g <- as.numeric(o$gr(opt$par));Vfull<-sd$cov.fixed
    ids <- which(names(opt$par)=="beta");need(length(ids)==ncol(p$X),"Beta mapping changed.")
    beta <- setNames(as.numeric(opt$par[ids]),colnames(p$X));V<-Vfull[ids,ids,drop=FALSE];dimnames(V)<-list(names(beta),names(beta))
    scaled <- if(all(is.finite(g))&&all(is.finite(Vfull)))sqrt(max(0,as.numeric(crossprod(g,Vfull%*%g)))) else Inf
    o$fn(opt$par);parameters<-current(opt$par);report<-o$report(o$env$last.par)
    extent<-max(diff(range(p$sites$x_km)),diff(range(p$sites$y_km)))
    good<-opt$convergence==0L&&isTRUE(sd$pdHess)&&all(is.finite(V))&&all(diag(V)>0)&&is.finite(scaled)&&scaled<GRAD_TOL&&is.finite(opt$objective)
    checks<-data.frame(convergence_code=opt$convergence,pdHess=isTRUE(sd$pdHess),max_abs_gradient=max(abs(g)),gradient_covariance_norm=scaled,numerical_checks_ok=good,logLik=-opt$objective,parameter_df=length(opt$par),AIC=2*opt$objective+2*length(opt$par),nobs=nrow(p$frame),spatial_range_km=as.numeric(report$range),spatial_sd=as.numeric(report$spatial_sd),range_greater_than_1_5_extent=as.numeric(report$range)>1.5*extent,spatial_sd_below_0_001=as.numeric(report$spatial_sd)<.001,AR1_rho=as.numeric(report$rho),Gamma_shape=as.numeric(report$shape))
    fit<-list(version=VERSION,signature=p$signature,parameters=parameters,optimizer=opt,beta=beta,V=V,cov_fixed=Vfull,checks=checks,eta=as.numeric(report$eta),field_at_sites=report$spatial,warnings=unique(warns),plasticity_in_design=FALSE,voltinism_in_design=FALSE)
    atomic_save(fit,fitfile);csv(checks,file.path(cfg$out,"fit_checks.csv"));writeLines(unique(warns),file.path(cfg$out,"warnings.txt"))
    fx<-fixed_table(beta,V);csv(fx,file.path(cfg$out,"fixed_effects.csv"));csv(data.frame(term=rownames(V),V,check.names=FALSE),file.path(cfg$out,"fixed_effect_covariance.csv"))

    term <- "year_decade:latitude_within_species_10deg";z<-fx[fx$term==term,,drop=FALSE];need(nrow(z)==1L,"Latitude interaction missing.")
    result<-data.frame(term=term,estimate=z$estimate,std.error=z$std.error,conf.low=z$conf.low,conf.high=z$conf.high,p_Wald=z$p_Wald,growth_factor_ratio_per_10deg=exp(z$estimate),lower_ratio=exp(z$conf.low),upper_ratio=exp(z$conf.high),percent_difference_in_decadal_growth_factor_per_10deg=100*expm1(z$estimate),lower_percent=100*expm1(z$conf.low),upper_percent=100*expm1(z$conf.high))
    csv(result,file.path(cfg$out,"latitude_effect.csv"))
    lonterm<-"year_decade:longitude_within_species_10deg";lz<-fx[fx$term==lonterm,,drop=FALSE];csv(data.frame(term=lonterm,estimate=lz$estimate,std.error=lz$std.error,conf.low=lz$conf.low,conf.high=lz$conf.high,p_Wald=lz$p_Wald),file.path(cfg$out,"longitude_effect.csv"))

    # Population model-based linear slopes (annual AR1/spatial deviations excluded).
    d<-p$frame;bsp<-as.matrix(parameters$b_species);bsi<-as.matrix(parameters$b_site)
    pops<-unique(d[c("pop_id","SPECIES","SITE_ID","latitude_deg","longitude_deg","latitude_within_species_deg","longitude_within_species_deg","latitude_within_species_10deg","longitude_within_species_10deg","species_index","site_index")])
    pops$species_year_slope<-bsp[pops$species_index+1L,2];pops$site_year_slope<-bsi[pops$site_index+1L,2]
    pops$trend_log_per_decade<-beta["year_decade"]+beta[term]*pops$latitude_within_species_10deg+beta[lonterm]*pops$longitude_within_species_10deg+pops$species_year_slope+pops$site_year_slope
    pops$trend_percent_per_decade<-100*expm1(pops$trend_log_per_decade)
    csv(pops,file.path(cfg$out,"model_based_population_linear_trends.csv"))

    phase(cfg,if(good)"FIT_OK" else "NUMERICAL_REVIEW_REQUIRED",paste("latitude beta",signif(z$estimate,5),"[",signif(z$conf.low,5),",",signif(z$conf.high,5),"] | scaled gradient",signif(scaled,5)))
    if(good)"FIT_OK" else "NUMERICAL_REVIEW_REQUIRED"
  }

  worker <- function(cfg) {
    check_packages();disk_check(cfg$out,cfg$min_disk);load_dll(cfg$source_run)
    tryCatch({
      prepare_all(cfg)
      status<-fit_worker(cfg)
      csv(data.frame(status=status),file.path(cfg$out,"workflow_status.csv"))
      writeLines(c(
        "ALL-POPULATION INDEPENDENT ABUNDANCE-GEOGRAPHY MODEL",
        paste("Source likelihood engine:",cfg$source_run),
        "Abundance input: output/pheno_abund_estimates_allspp.csv",
        "No phenological plasticity and no voltinism term/grouping are used.",
        "All eligible abundance populations are retained after the same >0 and <=99.9th-percentile filter.",
        "Latitude/longitude are centred across all retained populations within each species.",
        "Model retains species/site temporal slopes, population AR1, and annual IID Matern spatial field.",
        "A fresh 15-km mesh is built because the all-population site universe can exceed the previous plasticity-selected universe.",
        "15 km is the mesh construction cutoff, not the fitted correlation range.",
        "Primary test: year_decade:latitude_within_species_10deg."
      ),file.path(cfg$out,"README.txt"))
      if(requireNamespace("zip",quietly=TRUE)) {
        rel<-list.files(cfg$out,recursive=TRUE,full.names=FALSE);rel<-rel[grepl("\\.(csv|txt)$",rel,ignore.case=TRUE)]
        zip::zipr(file.path(cfg$out,"results_to_review.zip"),rel,root=cfg$out,mode="mirror",include_directories=FALSE)
      }
      phase(cfg,if(grepl("^FIT_OK",status))"COMPLETED" else "COMPLETED_WITH_NUMERICAL_REVIEW",file.path(cfg$out,"results_to_review.zip"),TRUE)
      list(status=status,output=cfg$out)
    },error=function(e){writeLines(conditionMessage(e),file.path(cfg$out,"fatal_error.txt"));try(phase(cfg,"STOPPED_WITH_ERROR",conditionMessage(e),TRUE),silent=TRUE);list(status="STOPPED_WITH_ERROR",error=conditionMessage(e),output=cfg$out)})
  }

  dump_engine <- function(out) {
    e<-environment(dump_engine);n<-ls(e,all.names=TRUE);fun<-n[vapply(n,function(k)is.function(get(k,envir=e)),logical(1))]
    p<-file.path(out,"workflow_engine.R");dump(c("VERSION","ROOT","SOURCE_RUN","GRAD_TOL",fun),file=p,envir=e,control="all");normalizePath(p,winslash="/",mustWork=TRUE)
  }

  run <- function(project_root=ROOT,source_run=SOURCE_RUN,monitor=FALSE,minimum_free_GiB=25,optimizer_iterations=1200L,optimizer_evaluations=1800L) {
    check_packages();cur<-state$job;if(!is.null(cur)&&cur$process$is_alive()){message("Workflow already running.");return(invisible(cur))}
    project_root<-normalizePath(project_root,winslash="/",mustWork=TRUE);source_run<-normalizePath(source_run,winslash="/",mustWork=TRUE);fits<-current_fits(source_run)
    outroot<-file.path(project_root,"output","population_trends","independent_geography_all_populations_15km");dir.create(outroot,recursive=TRUE,showWarnings=FALSE)
    out<-file.path(outroot,paste0("run_",format(Sys.time(),"%Y%m%d_%H%M%S"),"_",Sys.getpid()));dir.create(out,recursive=TRUE,showWarnings=FALSE);out<-normalizePath(out,winslash="/",mustWork=TRUE)
    cfg<-list(version=VERSION,project_root=project_root,source_run=source_run,current_fits=fits,out=out,progress=file.path(out,"progress.rds"),min_disk=minimum_free_GiB,iter=as.integer(optimizer_iterations),eval=as.integer(optimizer_evaluations));cfg$engine<-dump_engine(out);atomic_save(cfg,file.path(out,"run_config.rds"));phase(cfg,"QUEUED","One all-population geography fit; no voltinism/plasticity")
    log<-file.path(out,paste0("controller_",format(Sys.time(),"%Y%m%d_%H%M%S"),".log"))
    child<-callr::r_bg(function(engine,cfg){e<-new.env(parent=globalenv());sys.source(engine,envir=e);e$worker(cfg)},args=list(engine=cfg$engine,cfg=cfg),libpath=.libPaths(),stdout=log,stderr="2>&1",supervise=TRUE,user_profile=FALSE,system_profile=FALSE,wd=project_root,env=c(callr::rcmd_safe_env(),OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1"))
    job<-list(process=child,out=out,log=log,started=Sys.time());state$job<-job;options(phenoimpact.allpop.geo.job=job)
    message("All-population geography model launched. Keep R open.");message("Output: ",out);message("Monitor: pheno_allpop_geo$status()")
    if(monitor)status(job);invisible(job)
  }

  status <- function(job=NULL) {
    if(is.null(job))job<-state$job;if(is.null(job))job<-getOption("phenoimpact.allpop.geo.job");need(!is.null(job),"No registered job.")
    p<-file.path(job$out,"progress.rds");if(file.exists(p)){z<-tryCatch(readRDS(p),error=function(e)NULL);if(!is.null(z))print(z)}
    cat("\nAlive:",job$process$is_alive(),"\n")
    if(!job$process$is_alive()) {
      cat("Exit status:",job$process$get_exit_status(),"\n")
      f<-file.path(job$out,"fit_checks.csv");if(file.exists(f)){cat("\nFit checks:\n");print(read.csv(f),row.names=FALSE)}
      l<-file.path(job$out,"latitude_effect.csv");if(file.exists(l)){cat("\nLatitude effect:\n");print(read.csv(l),row.names=FALSE)}
      z<-file.path(job$out,"results_to_review.zip");if(file.exists(z))cat("\nFILE TO SEND:\n",z,"\n",sep="")
      e<-file.path(job$out,"fatal_error.txt");if(file.exists(e))cat("\nError:\n",paste(readLines(e,warn=FALSE),collapse="\n"),"\n",sep="")
    }
    invisible(job)
  }

  list(run=run,status=status)
})

run_all_population_abundance_geography <- function(...) pheno_allpop_geo$run(...)
status_all_population_abundance_geography <- function(...) pheno_allpop_geo$status(...)
