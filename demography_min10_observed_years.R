# phenoIMPACT: >=10 distinct observed abundance years per species x site.
# Sourcing defines the API only. The default run mode is preparation, not fitting.
# Reuses the baseline's saved code/likelihood, not mutable root scripts.
pheno_demo_min10 <- local({
  VERSION <- 'demography_min10_observed_years_1.0.0'
  GROUPS <- c('multivoltine', 'univoltine')
  ROOT <- 'E:/phenoIMPACT project/code/phenoIMPACT'
  REL <- 'output/population_trends/post_environment_mundlak_15km/run_d984406cdf33'
  NAME <- 'demography_min10_observed_years_15km'
  SCRIPT <- 'demography_min10_observed_years.R'
  state <- new.env(parent=emptyenv())
  need <- function(ok,msg) if(!isTRUE(ok)) stop(msg,call.=FALSE)
  packages <- function(x) {
    absent <- x[!vapply(x,requireNamespace,logical(1),quietly=TRUE)]
    need(!length(absent),paste('Missing packages:',paste(absent,collapse=', ')))
  }
  csv <- function(x,p) utils::write.csv(x,p,row.names=FALSE,na='')
  atomic <- function(x,p) {
    dir.create(dirname(p),recursive=TRUE,showWarnings=FALSE)
    tmp <- tempfile('.writing_',tmpdir=dirname(p),fileext='.rds')
    on.exit(unlink(tmp),add=TRUE); saveRDS(x,tmp,compress='gzip')
    backup <- paste0(p,'.previous')
    if(file.exists(backup)) unlink(backup)
    had <- file.exists(p)
    if(had) need(file.rename(p,backup),paste('Cannot back up',p))
    if(!file.rename(tmp,p)) {
      if(had) file.rename(backup,p)
      stop('Cannot save: ',p,call.=FALSE)
    }
    if(file.exists(backup)) unlink(backup)
    invisible(p)
  }
  phase <- function(cfg,stage,detail='',group='',terminal=FALSE) {
    atomic(list(stage=stage,detail=detail,group=group,time=Sys.time(),terminal=terminal),cfg$progress)
    message(format(Sys.time(),'%F %T'),' | ',stage,' | ',group,' | ',detail)
  }
  alive <- function(x) {
    packages('ps')
    need(is.list(x)&&length(x$pid)==1L&&is.finite(x$pid)&&
      length(x$create_time)==1L&&is.finite(as.numeric(x$create_time)), 'Invalid process ownership metadata')
    tryCatch({h<-ps::ps_handle(x$pid); ps::ps_is_running(h)&&
      abs(as.numeric(ps::ps_create_time(h))-as.numeric(x$create_time))<.001},
      error=function(e) {if(inherits(e,'no_such_process')) FALSE else stop(e)})
  }
  owners <- function(root,own=NULL) {
    paths <- list.files(file.path(root,'output'),
      pattern='^(worker_owner|fit_owner|isolated_worker_owner)\\.rds$',recursive=TRUE,full.names=TRUE)
    if(!is.null(own)) paths <- paths[!startsWith(normalizePath(paths,winslash='/'),paste0(own,'/'))]
    paths[vapply(paths,function(p) alive(readRDS(p)),logical(1))]
  }
  manifest <- function(paths) {
    need(all(file.exists(paths)),paste('Missing source files:',paste(paths[!file.exists(paths)],collapse=', ')))
    data.frame(path=normalizePath(paths,winslash='/',mustWork=TRUE),
      bytes=unname(file.info(paths)$size),md5=unname(tools::md5sum(paths)),stringsAsFactors=FALSE)
  }
  check_sources <- function(cfg) {
    need(identical(manifest(cfg$manifest$path),cfg$manifest),'Source files changed; prepare a new signed run')
    current <- vapply(names(cfg$libraries),function(n) as.character(utils::packageVersion(n)),'')
    need(identical(current,cfg$libraries),'R package versions differ from the prepared run')
  }
  helpers <- function(script) {
    e<-new.env(parent=globalenv()); sys.source(script,envir=e)
    need(exists('pheno_post_env',e,inherits=FALSE),'Missing baseline post-environment API')
    environment(e$pheno_post_env$run)
  }
  engine <- function(cfg) {
    e<-new.env(parent=globalenv()); sys.source(cfg$engine,envir=e)
    x<-environment(e$pheno_mundlak_qc$run); x$VERSION<-VERSION; x
  }
  counts <- function(d) {
    need(!anyNA(d[c('SPECIES','SITE_ID','pop_id','YEAR')])&&
      !anyDuplicated(d[c('SPECIES','SITE_ID','YEAR')]),'Missing identities or duplicate population-year')
    # The threshold is on observations, not last year minus first year.
    n<-vapply(split(d$YEAR,d$pop_id),function(y) length(unique(y)),integer(1))
    data.frame(pop_id=names(n),observed_years=unname(n),stringsAsFactors=FALSE)
  }
  sample_row <- function(d,group,scenario) data.frame(group=group,scenario=scenario,
    observations=nrow(d),populations=length(unique(d$pop_id)),species=length(unique(d$SPECIES)),
    sites=length(unique(d$SITE_ID)),networks=length(unique(d$bms_id)))
  subset_input <- function(old,fit,minimum_years,helper,run_signature) {
    packages(c('Matrix','digest')) # Register S4 sparse subsetting in fresh workers.
    need(identical(old$signature,fit$signature)&&identical(old$group,fit$group),'Baseline fit/input identity mismatch')
    need(isTRUE(fit$checks$numerical_checks_ok),'Baseline demographic fit requires numerical review')
    d<-as.data.frame(old$frame); audit<-counts(d)
    need(all(is.finite(d$ABUND_INDEX)&d$ABUND_INDEX>0),'Invalid Gamma response')
    need(identical(colnames(old$X),names(fit$beta)),'Baseline fixed coefficient mapping mismatch')
    need(isTRUE(all.equal(unname(d$ABUND_INDEX),old$data$y,tolerance=0))&&
      isTRUE(all.equal(unname(d$year_decade),old$data$time,tolerance=0)), 'Baseline response/time mismatch')
    audit$retained<-audit$observed_years>=minimum_years
    keep<-which(d$pop_id %in% audit$pop_id[audit$retained])
    need(length(keep)>0,'No population satisfies the observed-year threshold')
    keep<-keep[order(d$pop_id[keep],d$YEAR[keep])]
    f<-d[keep,,drop=FALSE]
    index_map<-function(labels,index,n) {
      z<-unique(data.frame(label=as.character(labels),index=as.integer(index)))
      need(nrow(z)==n&&!anyDuplicated(z$label)&&!anyDuplicated(z$index)&&
        identical(sort(z$index),seq_len(n)-1L),'Invalid saved latent index mapping')
      z[order(z$index),,drop=FALSE]
    }
    sp<-index_map(d$SPECIES,old$data$species,nrow(fit$parameters$b_species))
    si<-index_map(d$SITE_ID,old$data$site,nrow(fit$parameters$b_site))
    sp<-sp[sp$label %in% f$SPECIES,,drop=FALSE]; si<-si[si$label %in% f$SITE_ID,,drop=FALSE]
    need(identical(as.character(old$sites$SITE_ID)[old$data$site+1L],as.character(d$SITE_ID)),
      'Saved site/spatial projection order mismatch')
    sites<-old$sites[si$index+1L,,drop=FALSE]
    # Keep population plasticity z scores (original group-wide scales), but
    # recalculate species means and deviations over retained populations.
    dec<-helper$decompose(f,sites); f<-dec$frame; X<-dec$X
    need(identical(colnames(X),names(fit$beta))&&qr(X)$rank==ncol(X),
      'Filtered fixed-effect design changed or is rank deficient')
    dat<-old$data; dat$X<-unname(X); dat$y<-old$data$y[keep]; dat$time<-old$data$time[keep]
    dat$species<-as.integer(match(f$SPECIES,sp$label)-1L)
    dat$site<-as.integer(match(f$SITE_ID,si$label)-1L)
    dat$year<-as.integer(match(f$YEAR,old$years)-1L)
    dat$previous<-rep(-1L,nrow(f)); dat$gap<-rep(1L,nrow(f))
    if(nrow(f)>1L) {
      i<-which(f$pop_id[-1L]==f$pop_id[-nrow(f)])+1L
      dat$previous[i]<-i-2L; dat$gap[i]<-as.integer(f$YEAR[i]-f$YEAR[i-1L])
    }
    need(all(dat$gap>0L)&&!anyNA(dat$year),'Invalid calendar gaps or year mapping')
    dat$A<-old$data$A[si$index+1L,,drop=FALSE]
    pars<-fit$parameters
    pars$b_species<-pars$b_species[sp$index+1L,,drop=FALSE]
    pars$b_site<-pars$b_site[si$index+1L,,drop=FALSE]; pars$ar_obs<-pars$ar_obs[keep]
    target<-drop(old$X[keep,,drop=FALSE] %*% fit$beta)
    pars$beta<-unname(qr.coef(qr(X),target))
    need(all(is.finite(pars$beta))&&ncol(pars$field)==length(old$years)&&
      nrow(pars$field)==nrow(dat$M0),'Invalid warm start or spatial mesh dimensions')
    p<-old; p$frame<-f; p$populations<-dec$populations; p$X<-X
    p$data<-dat; p$parameters<-pars; p$sites<-sites; p$sample_audit<-audit
    p$component_scales<-cbind(group=old$group,dec$scales)
    p$version<-VERSION; p$minimum_observed_years<-minimum_years
    p$signature<-digest::digest(list(VERSION,run_signature,old$group,f,dat,pars),algo='sha256')
    p$source_bridge_signature<-old$signature; p$source_data_signature<-old$signature
    p$source_fit_objective<-fit$optimizer$objective
    p$predictors_updated<-TRUE; p$sample_filtered<-TRUE
    p$plasticity_uncertainty_propagated<-FALSE; p$residual_diagnostics<-'PENDING'
    p$retained_source_rows<-keep
    p$warm_projection_RMSE<-sqrt(mean((drop(X %*% pars$beta)-target)^2))
    p
  }
  baseline_files <- function(root,baseline_run) {
    need(file.exists(file.path(baseline_run,'run_config.rds')),'Baseline run_config.rds missing')
    old<-readRDS(file.path(baseline_run,'run_config.rds'))
    need(identical(old$version,'post_environment_mundlak_1.2.1'),'Unsupported baseline workflow version')
    script<-file.path(baseline_run,'code',basename(old$script))
    eng<-file.path(baseline_run,'code',basename(old$engine))
    dll<-file.path(old$source_run,'compiled',paste0('trend_spde',.Platform$dynlib.ext))
    paths<-c(file.path(root,SCRIPT),file.path(baseline_run,'run_config.rds'),script,eng,
      file.path(old$source_run,'compiled','trend_spde.cpp'),dll,
      unlist(lapply(GROUPS,function(g)c(file.path(baseline_run,paste0('prepared_mundlak_',g,'.rds')),
        file.path(baseline_run,paste0(g,'__mundlak'),'fit.rds')))))
    list(old=old,script=script,engine=eng,manifest=manifest(paths))
  }
  preflight <- function(project_root=ROOT,baseline_run=file.path(project_root,REL),minimum_years=10L,
                        minimum_ram_GiB=35,minimum_disk_GiB=40) {
    need(length(minimum_years)==1L&&is.finite(minimum_years)&&minimum_years>=10&&minimum_years==as.integer(minimum_years),
      'minimum_years must be an integer >=10')
    need(is.finite(minimum_ram_GiB)&&minimum_ram_GiB>=20&&is.finite(minimum_disk_GiB)&&minimum_disk_GiB>=20,
      'Resource admission thresholds must be >=20 GiB')
    packages(c('Matrix','digest','ps','callr','TMB','sdmTMB','spdep'))
    root<-normalizePath(project_root,winslash='/',mustWork=TRUE)
    baseline<-normalizePath(baseline_run,winslash='/',mustWork=TRUE)
    src<-baseline_files(root,baseline)
    versions<-vapply(names(src$old$libraries),function(n) as.character(utils::packageVersion(n)),'')
    need(identical(versions,src$old$libraries),'Use the same R package versions as the baseline fit')
    rows<-list()
    for(g in GROUPS) {
      p<-readRDS(file.path(baseline,paste0('prepared_mundlak_',g,'.rds')))
      audit<-counts(p$frame); keep<-p$frame$pop_id %in% audit$pop_id[audit$observed_years>=minimum_years]
      rows[[g]]<-rbind(sample_row(p$frame,g,'baseline'),sample_row(p$frame[keep,],g,paste0('min',minimum_years)))
      rm(p); gc()
    }
    active<-owners(root)
    ram<-as.numeric(ps::ps_system_memory()$avail)/1024^3
    disk<-as.numeric(ps::ps_disk_usage(root)$available)/1024^3
    result<-list(ready_to_prepare=TRUE,ready_to_fit=!length(active)&&ram>=minimum_ram_GiB&&disk>=minimum_disk_GiB,
      sample=do.call(rbind,rows),blocking_jobs=active,available_RAM_GiB=ram,available_disk_GiB=disk,
      minimum_years=as.integer(minimum_years),root=root,baseline=baseline,sources=src,libraries=versions,
      minimum_ram_GiB=minimum_ram_GiB,minimum_disk_GiB=minimum_disk_GiB)
    print(result$sample,row.names=FALSE)
    message('Ready to prepare: TRUE; ready to fit: ',result$ready_to_fit,
      '. RAM ',round(ram,1),' GiB; disk ',round(disk,1),' GiB. No process launched.')
    invisible(result)
  }
  prepare_group <- function(cfg,g) {
    check_sources(cfg)
    path<-file.path(cfg$out,paste0('prepared_mundlak_',g,'.rds'))
    marker<-file.path(cfg$out,paste0('prepared_',g,'_completed.rds'))
    if(file.exists(path)&&file.exists(marker)) {
      m<-readRDS(marker); p<-readRDS(path)
      need(identical(m$run_signature,cfg$signature)&&identical(m$signature,p$signature), 'Preparation cache mismatch')
      return(invisible(path))
    }
    old<-readRDS(file.path(cfg$baseline,paste0('prepared_mundlak_',g,'.rds')))
    fit<-readRDS(file.path(cfg$baseline,paste0(g,'__mundlak'),'fit.rds'))
    p<-subset_input(old,fit,cfg$minimum_years,helpers(cfg$helper),cfg$signature)
    p$source_prepared_file<-file.path(cfg$baseline,paste0('prepared_mundlak_',g,'.rds'))
    p$source_fit_file<-file.path(cfg$baseline,paste0(g,'__mundlak'),'fit.rds')
    csv(p$sample_audit,file.path(cfg$out,paste0('population_selection_',g,'.csv')))
    csv(p$scaling,file.path(cfg$out,paste0('plasticity_scaling_',g,'.csv')))
    csv(p$component_scales,file.path(cfg$out,paste0('component_scales_',g,'.csv')))
    csv(engine(cfg)$predictor_vif(p$populations,colnames(p$X)[3:8]),
      file.path(cfg$out,paste0('predictor_VIF_',g,'.csv')))
    atomic(p,path); atomic(list(run_signature=cfg$signature,signature=p$signature),marker)
    phase(cfg,'PREPARED',paste(nrow(p$frame),'rows;',nrow(p$populations),'populations'),g)
    invisible(path)
  }
  resource_check <- function(cfg) {
    check_sources(cfg)
    active<-owners(cfg$root,own=cfg$out)
    need(!length(active),paste('Another workflow is active:',paste(active,collapse='; ')))
    ram<-as.numeric(ps::ps_system_memory()$avail)/1024^3
    disk<-as.numeric(ps::ps_disk_usage(cfg$out)$available)/1024^3
    need(ram>=cfg$minimum_ram_GiB&&disk>=cfg$minimum_disk_GiB,
      paste('Insufficient admission resources. RAM:',round(ram,1),'GiB; disk:',round(disk,1),'GiB'))
  }
  stage <- function(cfg,g,action) {
    check_sources(cfg)
    if(action=='prepare') return(prepare_group(cfg,g))
    if(action=='diagnostics') return(helpers(cfg$helper)$diagnostics(cfg,g))
    need(action=='fit','Unknown stage')
    resource_check(cfg)
    atomic(list(pid=Sys.getpid(),create_time=ps::ps_create_time(ps::ps_handle()),
      group=g,out=cfg$out),file.path(cfg$out,'fit_owner.rds'))
    eng<-engine(cfg)
    # Verify DLL compatibility on the exact original sample before fitting.
    replay<-file.path(cfg$out,paste0('baseline_likelihood_replay_',g,'.rds'))
    if(!file.exists(replay)) {
      eng$load_dll(cfg$source_run)
      old<-readRDS(file.path(cfg$baseline,paste0('prepared_mundlak_',g,'.rds')))
      fit<-readRDS(file.path(cfg$baseline,paste0(g,'__mundlak'),'fit.rds'))
      obj<-eng$make_objective(old$data,fit$parameters)
      on.exit(try(TMB::FreeADFun(obj),silent=TRUE),add=TRUE)
      val<-obj$fn(obj$par); eta<-as.numeric(obj$report(obj$env$last.par)$eta)
      error<-abs(val-fit$optimizer$objective); eta_error<-max(abs(eta-fit$eta))
      need(error<.001&&eta_error<1e-5,'Baseline DLL likelihood replay failed')
      atomic(list(run_signature=cfg$signature,objective_error=error,eta_error=eta_error),replay)
      TMB::FreeADFun(obj); obj<-NULL; rm(old,fit); gc()
    } else need(identical(readRDS(replay)$run_signature,cfg$signature),'Replay cache mismatch')
    eng$fit_group(cfg,g)
  }
  review <- function(job=NULL) {
    if(is.null(job)) job<-state$job
    need(!is.null(job),'Pass a job returned by run()')
    cfg<-readRDS(file.path(job$out,'run_config.rds')); check_sources(cfg)
    rows<-list(); checks<-list(); residuals<-list(); balances<-list()
    for(g in GROUPS) {
      path<-file.path(cfg$out,paste0(g,'__mundlak'),'fit.rds')
      if(!file.exists(path)) next
      p<-readRDS(file.path(cfg$out,paste0('prepared_mundlak_',g,'.rds'))); f<-readRDS(path)
      need(identical(p$signature,f$signature),'Review fit/input mismatch')
      old<-readRDS(file.path(cfg$baseline,paste0(g,'__mundlak'),'fit.rds'))
      terms<-grep('^year_decade:',names(f$beta),value=TRUE)
      se<-sqrt(diag(f$V))[terms]; oldse<-sqrt(diag(old$V))[terms]
      rows[[g]]<-data.frame(group=g,term=terms,baseline_beta=unname(old$beta[terms]),
        min_years_beta=unname(f$beta[terms]),baseline_low=unname(old$beta[terms]-1.96*oldse),
        baseline_high=unname(old$beta[terms]+1.96*oldse),min_years_low=unname(f$beta[terms]-1.96*se),
        min_years_high=unname(f$beta[terms]+1.96*se),
        same_sign=unname(sign(old$beta[terms])==sign(f$beta[terms])),
        baseline_percent=unname(100*expm1(old$beta[terms])),min_years_percent=unname(100*expm1(f$beta[terms])),
        SE_ratio=unname(se/oldse),new_nobs=nrow(p$frame),new_npop=nrow(p$populations),
        numerical_checks_ok=isTRUE(f$checks$numerical_checks_ok))
      checks[[g]]<-f$checks
      for(scenario in c('baseline','min_years')) {
        folder<-if(scenario=='baseline') cfg$baseline else cfg$out
        rd<-file.path(folder,paste0(g,'__mundlak'),'residual_diagnostics')
        if(all(file.exists(file.path(rd,c('distribution.csv','temporal.csv','annual_spatial_Moran.csv'))))) {
          dist<-utils::read.csv(file.path(rd,'distribution.csv'))
          tm<-utils::read.csv(file.path(rd,'temporal.csv'))
          mo<-utils::read.csv(file.path(rd,'annual_spatial_Moran.csv'))
          residuals[[paste(g,scenario)]]<-data.frame(group=g,scenario=scenario,
            residual_mean=dist$mean,residual_SD=dist$SD,
            lag1=tm$correlation[tm$lag_years==1],mean_Moran_I=mean(mo$Moran_I,na.rm=TRUE),
            years_spatial_BH_below_05=sum(stats::p.adjust(mo$p_two_sided,'BH')<.05,na.rm=TRUE))
        }
        balance_path<-file.path(folder,paste0(g,'__mundlak'),'onset_offset_balance_within_between.csv')
        if(file.exists(balance_path)) {
          z<-utils::read.csv(balance_path); z$scenario<-scenario
          balances[[paste(g,scenario)]]<-z
        }
      }
    }
    if(length(rows)) csv(do.call(rbind,rows),file.path(cfg$out,'comparison_with_baseline.csv'))
    if(length(checks)) csv(do.call(rbind,checks),file.path(cfg$out,'review_fit_checks.csv'))
    if(length(residuals)) csv(do.call(rbind,residuals),file.path(cfg$out,'comparison_residual_screens.csv'))
    if(length(balances)) csv(do.call(rbind,balances),file.path(cfg$out,'comparison_onset_offset_balance.csv'))
    writeLines(c('Descriptive comparison on different samples; no between-fit LRT or AIC ranking.',
      'Same group-wide plasticity SD; within/between components recomputed on retained populations.',
      'Changes reflect sample composition as well as observed-year support.',
      'Numerical convergence is separate from residual review. Inspect residual_diagnostics for both groups.',
      'Plasticity uncertainty propagation is excluded by user decision.'),file.path(cfg$out,'review_interpretation.txt'))
    invisible(rows)
  }
  release_lock <- function(cfg,token) {
    lock<-cfg$lock; owner<-file.path(lock,'owner.rds')
    if(!file.exists(owner)||!identical(readRDS(owner)$token,token)) return(invisible(FALSE))
    need(identical(normalizePath(lock,winslash='/',mustWork=TRUE),paste0(cfg$out,'/pipeline.lock')),
      'Refusing lock removal outside the signed output directory')
    unlink(lock,recursive=TRUE); invisible(TRUE)
  }
  controller <- function(cfg) {
    owner<-list(pid=Sys.getpid(),create_time=ps::ps_create_time(ps::ps_handle()),token=cfg$token,
      out=cfg$out,log=cfg$log)
    atomic(owner,file.path(cfg$lock,'owner.rds')); atomic(owner,file.path(cfg$out,'worker_owner.rds'))
    on.exit({
      release_lock(cfg,cfg$token)
      path<-file.path(cfg$out,'worker_owner.rds')
      if(file.exists(path)&&identical(readRDS(path)$pid,Sys.getpid())) unlink(path)
    },add=TRUE)
    tryCatch({
      statuses<-list()
      for(g in GROUPS) {
        for(action in c('prepare',if(cfg$mode=='fit') c('fit','diagnostics'))) {
          if(action=='fit') resource_check(cfg)
          result<-callr::r(function(script,cfg,g,action) {
            e<-new.env(parent=globalenv()); sys.source(script,envir=e)
            environment(e$pheno_demo_min10$run)$stage(cfg,g,action)
          },args=list(script=cfg$script,cfg=cfg,g=g,action=action),libpath=cfg$libpaths,
            user_profile=FALSE,system_profile=FALSE,show=TRUE)
          if(action=='fit') fit_status<-result
        }
        statuses[[g]]<-data.frame(group=g,status=if(cfg$mode=='fit') fit_status else 'PREPARED',
          residual_review='PENDING')
        csv(do.call(rbind,statuses),file.path(cfg$out,'workflow_status.csv'))
      }
      engine(cfg)$combine_outputs(cfg,do.call(rbind,statuses))
      if(cfg$mode=='fit') review(list(out=cfg$out))
      final<-if(cfg$mode=='prepare') 'PREPARED' else if(all(vapply(statuses,
        function(z) z$status %in% c('FIT_OK','FIT_OK_CACHED'),logical(1)))) 'FITTED_RESIDUAL_REVIEW_PENDING' else 'NUMERICAL_REVIEW_REQUIRED'
      phase(cfg,final,'Finished; review exported comparisons and residual diagnostics',terminal=TRUE)
    },error=function(e) {
      atomic(list(error=conditionMessage(e),time=Sys.time()),file.path(cfg$out,'last_error.rds'))
      phase(cfg,'STOPPED_WITH_ERROR',conditionMessage(e),terminal=TRUE); stop(e)
    })
    invisible(list(out=cfg$out,log=cfg$log))
  }
  run <- function(project_root=ROOT,baseline_run=file.path(project_root,REL),minimum_years=10L,
                  mode=c('prepare','fit'),background=TRUE,minimum_ram_GiB=35,minimum_disk_GiB=40) {
    mode<-match.arg(mode)
    pf<-preflight(project_root,baseline_run,minimum_years,minimum_ram_GiB,minimum_disk_GiB)
    if(mode=='fit') need(pf$ready_to_fit,'Fit admission blocked: inspect preflight resources and active jobs')
    signature<-digest::digest(list(VERSION,pf$sources$manifest,pf$libraries,pf$minimum_years),algo='sha256')
    out<-file.path(pf$root,'output','population_trends',NAME,paste0('run_',substr(signature,1,12)))
    dir.create(out,recursive=TRUE,showWarnings=FALSE); out<-normalizePath(out,winslash='/',mustWork=TRUE)
    lock<-file.path(out,'pipeline.lock')
    cfg<-list(version=VERSION,signature=signature,run_signature=signature,root=pf$root,baseline=pf$baseline,
      out=out,lock=lock,progress=file.path(out,'progress.rds'),minimum_years=pf$minimum_years,
      minimum_ram_GiB=minimum_ram_GiB,minimum_disk_GiB=minimum_disk_GiB,mode=mode,
      source_run=pf$sources$old$source_run,manifest=pf$sources$manifest,libraries=pf$libraries,
      libpaths=.libPaths(),iter=pf$sources$old$iter,eval=pf$sources$old$eval,
      token=paste(Sys.getpid(),format(Sys.time(),'%Y%m%d%H%M%OS6'),sep='_'),
      log=file.path(out,paste0('controller_',format(Sys.time(),'%Y%m%d_%H%M%S'),'.log')))
    if(dir.exists(lock)) {
      owner<-file.path(lock,'owner.rds'); need(file.exists(owner),'Lock has no owner; inspect manually')
      previous<-readRDS(owner); need(!alive(previous),'This sensitivity workflow is already active')
      release_lock(cfg,previous$token)
    }
    need(dir.create(lock,showWarnings=FALSE),'Cannot acquire sensitivity lock')
    atomic(list(pid=Sys.getpid(),create_time=ps::ps_create_time(ps::ps_handle()),token=cfg$token),file.path(lock,'owner.rds'))
    launched<-FALSE; on.exit(if(!launched) release_lock(cfg,cfg$token),add=TRUE)
    code<-file.path(out,'code'); dir.create(code,showWarnings=FALSE)
    files<-c(file.path(pf$root,SCRIPT),pf$sources$script,pf$sources$engine)
    dest<-file.path(code,basename(files))
    for(i in seq_along(files)) {
      if(file.exists(dest[i])) need(identical(unname(tools::md5sum(dest[i])),unname(tools::md5sum(files[i]))),
        'Saved code differs from the signed sources') else need(file.copy(files[i],dest[i]),'Cannot snapshot code')
    }
    cfg$script<-dest[1]; cfg$helper<-dest[2]; cfg$engine<-dest[3]
    atomic(cfg,file.path(out,'run_config.rds')); csv(pf$sample,file.path(out,'sample_before_after.csv'))
    csv(cfg$manifest,file.path(out,'source_manifest.csv'))
    if(mode=='prepare'||!background) {
      job<-controller(cfg); state$job<-job; return(invisible(job))
    }
    phase(cfg,'QUEUED','Two sequential demographic fits; >=10 distinct observed years per population')
    child<-callr::r_bg(function(script,cfg) {
      e<-new.env(parent=globalenv()); sys.source(script,envir=e)
      environment(e$pheno_demo_min10$run)$controller(cfg)
    },args=list(script=cfg$script,cfg=cfg),libpath=cfg$libpaths,stdout=cfg$log,stderr='2>&1',
      supervise=FALSE,user_profile=FALSE,system_profile=FALSE,
      env=c(callr::rcmd_safe_env(),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1'))
    on.exit(if(!launched) try(child$kill_tree(),silent=TRUE),add=TRUE)
    owner<-list(pid=child$get_pid(),create_time=ps::ps_create_time(ps::ps_handle(child$get_pid())),
      token=cfg$token,out=out,log=cfg$log)
    atomic(owner,file.path(out,'worker_owner.rds'))
    launched<-TRUE; job<-c(owner,list(process=child)); state$job<-job
    message('Started sequential sensitivity queue. Monitor with pheno_demo_min10$status(job).')
    invisible(job)
  }
  status <- function(job=NULL) {
    if(is.null(job)) job<-state$job
    need(!is.null(job),'Pass the job returned by run(), or list(out="absolute run directory")')
    p<-file.path(job$out,'progress.rds'); if(file.exists(p)) print(readRDS(p))
    owner<-file.path(job$out,'worker_owner.rds')
    cat('Controller alive:',if(file.exists(owner)) alive(readRDS(owner)) else FALSE,'\nOutput:',job$out,'\n')
    s<-file.path(job$out,'workflow_status.csv'); if(file.exists(s)) print(utils::read.csv(s),row.names=FALSE)
    invisible(job)
  }
  list(preflight=preflight,run=run,status=status,review=review)
})
