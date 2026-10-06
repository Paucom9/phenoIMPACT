# phenoIMPACT: compare existing population trends with annual spatial models.
# Source this file from the phenoIMPACT R project. The original fits are preserved.
# Requires R >= 4.1, TMB/glmmTMB/Matrix/fmesher/callr/digest/here; ggplot2 for figures.
# Windows R 4.5 requires Rtools45 to compile this engine once per configuration.
# Default: TWO new full fits, multivoltines then univoltines. No original refits.
# Same Gamma/log ML, species/site independent intercepts and time slopes,
# population AR(1), plus an IID annual Matern field (nu=1, common range/variance).
# Mesh cutoff=15 km controls construction, not estimated correlation range.
# Baseline likelihood equivalence is checked before any spatial optimization.
# Original data/predictor scaling are read directly from the saved model frames.
# Missing population-years are integrated out analytically under the SAME AR(1).
# Embedded coordinates: the sites in the two reviewed baseline fits, EPSG:3035 km.
# Keep R open. Escape pauses the display; the queue continues.
# Functions: watch_spatial_trends(), spatial_trends_status(), stop_spatial_trends().
# To configure first: options(pheno.trends.spatial.functions_only=TRUE); source(...)
# Then: spatial_trend_job <- run_spatial_trends(spatial_trend_config)
# run_lrt=TRUE adds four reduced fits; matching full fits are reused.
# Spatial variance is NOT tested with an ordinary chi-square likelihood test.
# Inference remains conditional on estimated plasticity and abundance indices.

spatial_trend_config <- list(
  project_root = NULL,
  baseline_root = NULL,           # Default: output/population_trends/spatial_plasticity_15km
  baseline_files = NULL,          # Optional named paths: c(multivoltine='.../fit.rds', univoltine='.../fit.rds')
  coordinates_file = NULL,        # NULL: embedded reviewed coordinates; no manual join required
  output_root = NULL,             # Default: output/population_trends/spatial_abundance_15km
  groups = c('multivoltine','univoltine'),
  mesh_cutoff_km = 15,
  run_lrt = FALSE,                # First compare full-model estimates and Wald intervals
  prepare_only = FALSE,
  diagnostics = TRUE,
  permutations = 999L,
  seed = 20260928L,
  monitor_seconds = 30,
  optimizer_iterations = 1000L,
  optimizer_evaluations = 1500L,
  make_figures = TRUE,
  font_family = 'Garamond'
)

pheno_spatial_trends <- local({
assert <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call.=FALSE)
write_table <- function(x, path) utils::write.csv(x, path, row.names=FALSE, na='')
atomic_save <- function(x, path) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  tmp<-tempfile('writing_',tmpdir=dirname(path));on.exit(unlink(tmp),add=TRUE)
  saveRDS(x,tmp,compress=FALSE)
  backup<-paste0(path,'.previous')
  if(file.exists(path)) {
    assert(file.copy(path,backup,overwrite=TRUE),paste('Cannot back up',path))
    assert(file.remove(path),paste('Cannot replace',path))
  }
  if(!file.rename(tmp,path)) {
    if(file.exists(backup))file.copy(backup,path,overwrite=TRUE)
    stop('Cannot save ',path)
  }
  if(file.exists(backup))unlink(backup)
}
phase <- function(cfg, text, task='', details=list()) {
  atomic_save(c(list(phase=text,task=task,updated=Sys.time()),details),cfg$progress_file)
  message('[',format(Sys.time(),'%Y-%m-%d %H:%M:%S'),'] ',task,' ',text)
}
fixed_table <- function(beta,V) {
  se<-sqrt(diag(V));z<-beta/se
  data.frame(term=names(beta),estimate=as.numeric(beta),std.error=as.numeric(se),
    conf.low=as.numeric(beta-1.96*se),conf.high=as.numeric(beta+1.96*se),
    p_Wald=as.numeric(2*pnorm(-abs(z))))
}
# Read serialized matrices only; never invoke a saved external TMB pointer.
read_baseline <- function(path, coords, group) {
  saved<-readRDS(path);fit<-if(inherits(saved,'glmmTMB'))saved else saved$fit
  assert(inherits(fit,'glmmTMB'),'Expected a saved glmmTMB fit or list containing $fit.')
  assert(fit$fit$convergence==0 && isTRUE(fit$sdr$pdHess),'The baseline requires numerical review.')
  assert(fit$modelInfo$family$family=='Gamma' && fit$modelInfo$family$link=='log' &&
    !isTRUE(fit$modelInfo$REML),'Expected Gamma/log ML baseline.')
  native<-fit$obj$env$data;d<-fit$frame;pars<-fit$fit$parfull
  assert(all(native$weights==1) && all(native$offset==0) && ncol(native$Xzi)==0 &&
    ncol(if(is.null(native$Xdisp))native$Xd else native$Xdisp)==1,'Unexpected weights, offset, zero inflation or dispersion model.')
  expected<-c('(Intercept)','year_decade','onset_plasticity_z','offset_plasticity_bio_z',
    'year_decade:onset_plasticity_z','year_decade:offset_plasticity_bio_z')
  assert(setequal(colnames(native$X),expected),'Unexpected baseline fixed-effects design.')
  cnms<-fit$modelInfo$reTrms$cond$cnms;struc<-fit$modelInfo$reStruc$condReStruc
  flist<-fit$modelInfo$reTrms$cond$flist
  assert(identical(names(cnms),c('SPECIES','SITE_ID','pop_id')) &&
    all(vapply(cnms[1:2],identical,logical(1),c('(Intercept)','year_decade'))),
    'Expected independent species/site intercepts and slopes plus population AR(1).')
  theta<-unname(fit$fit$par[names(fit$fit$par)=='theta'])
  assert(length(theta)==6,'Unexpected covariance parameterization.')
  assert(as.integer(struc[[1]]$blockCode)==0L && as.integer(struc[[2]]$blockCode)==0L &&
    as.integer(struc[[3]]$blockCode)==3L,'Baseline covariance must be diag, diag, ar1.')
  beta<-setNames(as.numeric(pars[names(pars)=='beta']),colnames(native$X))
  ds<-intersect(c('betad','betadisp'),names(pars))
  assert(length(ds)==1 && sum(names(pars)==ds)==1,'Unexpected Gamma shape parameter.')
  b<-unname(pars[names(pars)=='b'])
  sizes<-vapply(struc,function(s)s$blockSize*s$blockReps,numeric(1));ends<-cumsum(sizes)
  assert(sum(sizes)==length(b),'Random effect dimensions disagree.')
  bsp<-matrix(b[seq_len(sizes[1])],ncol=2,byrow=TRUE)
  bsi<-matrix(b[seq.int(ends[1]+1,ends[2])],ncol=2,byrow=TRUE)
  ar<-as.numeric(native$Z[,seq.int(ends[2]+1,ends[3]),drop=FALSE]%*%b[seq.int(ends[2]+1,ends[3])])
  d$year_num<-as.integer(as.character(d$year_fac))
  ord<-order(as.character(d$pop_id),d$year_num);d<-d[ord,,drop=FALSE]
  assert(!anyDuplicated(d[c('pop_id','year_num')]),'Duplicated population-years.')
  years<-seq.int(min(d$year_num),max(d$year_num))
  assert(identical(as.integer(sub('^year_fac','',cnms[[3]])),years),'Saved AR(1) levels do not match calendar years.')
  sites<-levels(flist$SITE_ID);sp<-levels(flist$SPECIES)
  assert(!anyDuplicated(coords$SITE_ID),'Duplicate coordinate IDs.')
  loc<-coords[match(sites,coords$SITE_ID),,drop=FALSE]
  assert(!anyNA(loc) && all(is.finite(as.matrix(loc[c('x_km','y_km')]))),
    'Missing site coordinates; no observations were silently dropped.')
  n<-nrow(d);prev<-c(-1L,seq_len(n-1L)-1L)
  first<-c(TRUE,as.character(d$pop_id[-1])!=as.character(d$pop_id[-n]));prev[first]<- -1L
  gaps<-c(1L,diff(d$year_num));gaps[first]<-1L
  assert(all(gaps>0),'Non-increasing calendar years within a population.')
  X<-as.matrix(native$X[ord,,drop=FALSE]);ar<-ar[ord]
  isp<-match(as.character(d$SPECIES),sp);isi<-match(as.character(d$SITE_ID),sites)
  eta<-as.numeric(X%*%beta)+bsp[isp,1]+d$year_decade*bsp[isp,2]+bsi[isi,1]+d$year_decade*bsi[isi,2]+ar
  original_eta<-as.numeric(native$X%*%beta+native$Z%*%b)[ord]
  assert(max(abs(eta-original_eta))<1e-8,'Prediction reconstruction failed.')
  pidx<-which(names(fit$sdr$par.fixed)=='beta');V<-fit$sdr$cov.fixed[pidx,pidx,drop=FALSE]
  dimnames(V)<-list(names(beta),names(beta))
  frame<-data.frame(SPECIES=as.character(d$SPECIES),SITE_ID=as.character(d$SITE_ID),pop_id=as.character(d$pop_id),
    YEAR=d$year_num,year_decade=d$year_decade,ABUND_INDEX=d$ABUND_INDEX,
    onset_plasticity_z=d$onset_plasticity_z,offset_plasticity_bio_z=d$offset_plasticity_bio_z,
    bms_id=loc$bms_id[isi],x_km=loc$x_km[isi],y_km=loc$y_km[isi],stringsAsFactors=FALSE)
  data<-list(y=as.numeric(d$ABUND_INDEX),X=X,time=as.numeric(d$year_decade),
    species=as.integer(isp-1),site=as.integer(isi-1),year=as.integer(match(d$year_num,years)-1),
    previous=as.integer(prev),gap=as.integer(gaps))
  par<-list(beta=unname(beta),log_sd=theta[1:5],rho_raw=theta[6],log_shape=unname(pars[names(pars)==ds]),
    b_species=bsp,b_site=bsi,ar_obs=ar,log_range=log(100),log_spatial_sd=log(.15),field=matrix(numeric(),0,0))
  list(data=data,par=par,frame=frame,sites=loc,years=years,beta=beta,V=V,baseline_eta=eta,
    baseline_logLik= -fit$fit$objective,baseline_df=length(fit$fit$par),
    baseline_AIC=2*fit$fit$objective+2*length(fit$fit$par),
    source_version=if(length(fit$modelInfo$packageVersion))as.character(fit$modelInfo$packageVersion)else 'not_recorded',group=group,
    data_signature=digest::digest(list(frame,X),algo='sha256'))
}
with_mesh <- function(input, mesh=NULL, spatial=TRUE) {
  data<-input$data;par<-input$par
  zero<-Matrix::sparseMatrix(i=integer(),j=integer(),dims=c(0L,0L))
  data$use_spatial<-as.integer(spatial)
  if(spatial) {
    assert(!is.null(mesh),'A mesh is required.')
    data$A<-methods::as(fmesher::fm_basis(mesh$mesh,loc=as.matrix(input$sites[c('x_km','y_km')])), 'CsparseMatrix')
    assert(all(abs(Matrix::rowSums(data$A)-1)<1e-8),'Some coordinates lie outside the mesh.')
    data$M0<-mesh$M0;data$M1<-mesh$M1;data$M2<-mesh$M2
    par$field<-matrix(0,nrow(mesh$M0),length(input$years))
  } else {
    data$A<-Matrix::sparseMatrix(i=integer(),j=integer(),dims=c(nrow(input$sites),0L))
    data$M0<-zero;data$M1<-zero;data$M2<-zero
    par$field<-matrix(numeric(),0,0)
  }
  list(data=data,par=par)
}
make_objective <- function(input, spatial=TRUE, variant='full') {
  data<-input$data;par<-input$par
  for(nm in c('A','M0','M1','M2')) data[[nm]]<-methods::as(methods::as(methods::as(data[[nm]],'dMatrix'),'generalMatrix'),'TsparseMatrix')
  if(variant!='full') {
    variable<-switch(variant,without_onset_interaction='onset_plasticity_z',
      without_offset_interaction='offset_plasticity_bio_z',stop('Unknown variant.'))
    drop<-match(paste0('year_decade:',variable),colnames(data$X))
    assert(!is.na(drop),'Missing interaction.')
    data$X<-data$X[,-drop,drop=FALSE];par$beta<-par$beta[-drop]
  }
  random<-c('b_species','b_site','ar_obs')
  map<-NULL
  if(spatial) random<-c(random,'field') else map<-list(log_range=factor(NA),log_spatial_sd=factor(NA))
  obj<-TMB::MakeADFun(data=data,parameters=par,random=random,map=map,DLL='trend_spde',silent=TRUE,
    inner.control=list(maxit=1000,trace=FALSE))
  list(obj=obj,data=data,par=par)
}
load_dll <- function(cfg) {
  path<-TMB::dynlib(file.path(cfg$out,'compiled','trend_spde'))
  if(!'trend_spde'%in%names(getLoadedDLLs()))dyn.load(path)
}
compile_engine <- function(cfg) {
  folder<-file.path(cfg$out,'compiled');dir.create(folder,recursive=TRUE,showWarnings=FALSE)
  old<-setwd(folder);on.exit(setwd(old),add=TRUE)
  writeLines(cpp_source(),'trend_spde.cpp')
  if(!file.exists(TMB::dynlib('trend_spde'))) {
    phase(cfg,'Compiling model engine')
    code<-tryCatch(TMB::compile('trend_spde.cpp',flags='-O1 -g0'),error=function(e)stop(
      'Compilation failed. On Windows with R 4.5, install Rtools45 from https://cran.r-project.org/bin/windows/Rtools/ and restart R. Details: ',conditionMessage(e)))
    assert(code==0 && file.exists(TMB::dynlib('trend_spde')),'Model engine did not compile; inspect the controller log.')
  }
  invisible(NULL)
}
make_mesh <- function(loc,cutoff) {
  xy<-unique(as.matrix(loc[c('x_km','y_km')]))
  mesh<-fmesher::fm_rcdt_2d_inla(loc=xy,cutoff=cutoff)
  fem<-fmesher::fm_fem(mesh,order=2)
  list(mesh=mesh,M0=methods::as(fem$c0,'CsparseMatrix'),M1=methods::as(fem$g1,'CsparseMatrix'),
    M2=methods::as(fem$g2,'CsparseMatrix'),cutoff_km=cutoff)
}
prepare_worker <- function(cfg) {
  coords<-if(is.null(cfg$coordinates_file))default_coordinates() else utils::read.csv(cfg$coordinates_file,stringsAsFactors=FALSE)
  if(all(c('transect_id','transect_lon','transect_lat')%in%names(coords)))coords<-data.frame(
    SITE_ID=as.character(coords$transect_id),bms_id=coords$bms_id,x_km=coords$transect_lon/1000,y_km=coords$transect_lat/1000)
  assert(all(c('SITE_ID','bms_id','x_km','y_km')%in%names(coords)),'Coordinate columns are missing.')
  locations<-list();rows<-list()
  for(group in cfg$groups) {
    phase(cfg,'Reading saved baseline and preserving fitted rows',group)
    input<-read_baseline(cfg$baseline_files[[group]],coords,group)
    atomic_save(input,file.path(cfg$out,paste0('input_',group,'.rds')))
    write_table(fixed_table(input$beta,input$V),file.path(cfg$out,paste0('baseline_',group,'_fixed_effects.csv')))
    locations[[group]]<-input$sites
    rows[[group]]<-data.frame(group=group,n=nrow(input$frame),sites=nrow(input$sites),
      species=nrow(input$par$b_species),populations=length(unique(input$frame$pop_id)),
      latent_AR1_observed=length(input$par$ar_obs),years=length(input$years),
      original_glmmTMB_version=input$source_version,data_signature=input$data_signature)
    rm(input);gc()
  }
  loc<-unique(do.call(rbind,locations));loc<-loc[order(loc$SITE_ID),]
  assert(!anyDuplicated(loc$SITE_ID),'Coordinate mismatch between groups.')
  phase(cfg,'Building the common 15-km mesh')
  mesh<-make_mesh(loc,cfg$mesh_cutoff_km)
  atomic_save(mesh,file.path(cfg$out,'mesh.rds'))
  write_table(do.call(rbind,rows),file.path(cfg$out,'input_summary.csv'))
  write_table(loc,file.path(cfg$out,'coordinates_used.csv'))
  write_table(data.frame(cutoff_km=cfg$mesh_cutoff_km,vertices=nrow(mesh$M0),
    triangles=nrow(mesh$mesh$graph$tv)),file.path(cfg$out,'mesh_summary.csv'))
  grDevices::png(file.path(cfg$out,'mesh.png'),width=1400,height=1100,res=140)
  tryCatch({plot(mesh$mesh,main='',asp=1);points(loc$x_km,loc$y_km,pch=16,cex=.25,col='#31666b')},finally=grDevices::dev.off())
  invisible(TRUE)
}
validate_baseline <- function(input,cfg,task) {
  phase(cfg,'Verifying baseline likelihood in the spatial engine (no refit)',task)
  a<-make_objective(with_mesh(input,spatial=FALSE),spatial=FALSE)
  calculated<- -a$obj$fn(a$obj$par)
  tolerance<-max(.001,abs(input$baseline_logLik)*1e-8)
  ok<-is.finite(calculated) && abs(calculated-input$baseline_logLik)<=tolerance
  tab<-data.frame(group=input$group,saved_logLik=input$baseline_logLik,
    reconstructed_logLik=calculated,absolute_difference=abs(calculated-input$baseline_logLik),
    tolerance=tolerance,equivalent=ok)
  write_table(tab,file.path(cfg$out,paste0('baseline_equivalence_',input$group,'.csv')))
  TMB::FreeADFun(a$obj);rm(a);gc()
  assert(ok,'The new engine did not reproduce the baseline likelihood. No spatial fit was launched.')
  invisible(tab)
}
# Exact continuous Gamma PIT conditional on all fitted effects.
residual_table <- function(input,eta,shape) {
  mu<-exp(eta);y<-input$data$y
  lo<-pgamma(y,shape=shape,scale=mu/shape,log.p=TRUE)
  hi<-pgamma(y,shape=shape,scale=mu/shape,lower.tail=FALSE,log.p=TRUE)
  z<-ifelse(lo<log(.5),qnorm(lo,log.p=TRUE),qnorm(hi,lower.tail=FALSE,log.p=TRUE))
  d<-input$frame[c('SITE_ID','YEAR','bms_id','x_km','y_km')];d$z<-z
  assert(all(is.finite(z)),'Non-finite Gamma PIT residuals.')
  sy<-aggregate(z~SITE_ID+YEAR+bms_id+x_km+y_km,d,mean)
  xy<-aggregate(z~YEAR+bms_id+x_km+y_km,sy,mean)
  list(site_year=xy,summary=data.frame(n=length(z),mean=mean(z),sd=sd(z),
    fraction_abs_z_gt_1_96=mean(abs(z)>1.96)))
}
# Symmetric union of eight nearest neighbours <=100 km. Blocks restrict
# both the graph and permutation; no edge crosses a BMS in within-network tests.
moran_screen <- function(tab,nperm=999L,seed=20260928L) {
  set.seed(seed);answer<-list()
  for(yr in sort(unique(tab$YEAR))) {
    a<-tab[tab$YEAR==yr,];xy<-as.matrix(a[c('x_km','y_km')]);n<-nrow(a)
    D<-as.matrix(dist(xy));diag(D)<-Inf
    for(scope in c('all_networks','within_networks')) {
      DD<-D
      if(scope=='within_networks')DD[outer(a$bms_id,a$bms_id,'!=')]<-Inf
      edges<-lapply(seq_len(n),function(i) {
        ids<-which(DD[i,]>0 & DD[i,]<=100)
        head(ids[order(DD[i,ids])],8L)
      })
      W<-Matrix::sparseMatrix(i=rep(seq_len(n),lengths(edges)),j=unlist(edges),x=1,dims=c(n,n))
      W<-methods::as((W+Matrix::t(W))>0,'dMatrix')
      keep<-Matrix::rowSums(W)>0;W<-W[keep,keep,drop=FALSE];v<-a$z[keep];bms<-a$bms_id[keep]
      nn<-length(v);I<-p<-NA_real_
      if(nn>=15 && sd(v)>1e-12) {
        if(scope=='within_networks')v<-v-ave(v,bms,FUN=mean)
        v<-v-mean(v);den<-sum(v^2)
        if(den>1e-20) {
          W<-Matrix::Diagonal(x=1/Matrix::rowSums(W))%*%W
          I<-as.numeric(crossprod(v,W%*%v)/den)
          blocks<-if(scope=='within_networks')split(seq_len(nn),bms) else list(seq_len(nn))
          # Vectorized sparse multiplication in bounded batches.
          more<-0L
          for(start in seq.int(1L,nperm,by=100L)) {
            m<-min(100L,nperm-start+1L);P<-matrix(v,nrow=nn,ncol=m)
            for(j in seq_len(m))for(ind in blocks)P[ind,j]<-v[ind[sample.int(length(ind))]]
            ip<-colSums(P*(W%*%P))/den;more<-more+sum(ip>=I)
          }
          p<-(1+more)/(1+nperm)
        }
      }
      answer[[length(answer)+1L]]<-data.frame(YEAR=yr,scope=scope,n_connected=nn,n_excluded=n-nn,I=I,p_positive=p)
    }
  }
  out<-do.call(rbind,answer);out$q_BH<-ave(out$p_positive,out$scope,FUN=function(p)p.adjust(p,'BH'));out
}
export_predictions <- function(beta,V,years) {
  rows<-list()
  for(v in c('onset_plasticity_z','offset_plasticity_bio_z'))for(z in c(-1,0,1)) {
    it<-paste0('year_decade:',v);t<-(years-min(years))/10
    slope<-beta['year_decade']+beta[it]*z
    variance<-V['year_decade','year_decade']+z^2*V[it,it]+2*z*V['year_decade',it]
    assert(variance>=-1e-10,'Invalid fixed-effect prediction variance.')
    se<-sqrt(max(variance,0))*t;eta<-as.numeric(slope)*t
    rows[[length(rows)+1L]]<-data.frame(variable=v,plasticity_z=z,YEAR=years,
      percent_change=100*expm1(eta),low=100*expm1(eta-1.96*se),high=100*expm1(eta+1.96*se))
  }
  do.call(rbind,rows)
}
fit_worker <- function(cfg,group,variant='full') {
  task<-paste(group,variant,sep='__');out<-file.path(cfg$out,task);dir.create(out,showWarnings=FALSE)
  load_dll(cfg)
  phase(cfg,'Reading prepared input',task)
  input<-readRDS(file.path(cfg$out,paste0('input_',group,'.rds')))
  cache<-file.path(out,'fit.rds');cached<-file.exists(cache)
  if(cached) {
    fit<-readRDS(cache)
    assert(identical(fit$signature,cfg$signature) && identical(fit$data_signature,input$data_signature),
      'Saved fit does not match the current inputs/settings.')
    phase(cfg,'Reusing saved spatial estimates',task)
  } else {
    if(variant=='full')validate_baseline(input,cfg,task)
    a<-with_mesh(input,readRDS(file.path(cfg$out,'mesh.rds')))
    if(variant!='full') {
      full<-readRDS(file.path(cfg$out,paste0(group,'__full'),'fit.rds'))
      assert(isTRUE(full$checks$numerical_checks_ok),'The full spatial model requires review.')
      a$par<-full$parameters
    }
    built<-make_objective(a,variant=variant);obj<-built$obj
    checkpoint<-file.path(out,'last_evaluation.rds')
    # A recovery file stores outer parameters at successful evaluations.
    # This is a warm start, not an exact continuation of an interrupted optimizer.
    if(file.exists(checkpoint)) {
      warm<-readRDS(checkpoint)
      if(identical(warm$signature,cfg$signature) && length(warm$par)==length(obj$par))obj$par<-warm$par
    }
    count<-0L;grad_count<-0L;last<-Sys.time()-60;started<-Sys.time()
    objective<-function(par) {
      value<-obj$fn(par);count<<-count+1L
      if(is.finite(value) && (count==1L || as.numeric(difftime(Sys.time(),last,units='secs'))>=30)) {
        atomic_save(list(signature=cfg$signature,par=par,objective=value,time=Sys.time()),checkpoint)
        phase(cfg,'Fitting spatial Gamma model',task,list(fit_started=started,fn_evaluations=count,
          gradient_evaluations=grad_count,objective=value));last<<-Sys.time()
      }
      value
    }
    gradient<-function(par){grad_count<<-grad_count+1L;obj$gr(par)}
    phase(cfg,'Starting spatial optimization',task,list(fit_started=started))
    warnings<-character()
    opt<-withCallingHandlers(nlminb(obj$par,objective,gradient,
      control=list(iter.max=cfg$optimizer_iterations,eval.max=cfg$optimizer_evaluations,rel.tol=1e-10)),
      warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
    if(opt$convergence!=0L) {
      phase(cfg,'Refining a non-converged nlminb solution with BFGS',task,list(fit_started=started))
      second<-optim(opt$par,objective,gradient,method='BFGS',control=list(maxit=300,reltol=1e-10))
      if(is.finite(second$value) && second$value<=opt$objective+1e-6)opt<-list(par=second$par,
        objective=second$value,convergence=second$convergence,message='BFGS refinement',
        first_optimizer=opt,counts=second$counts)
    }
    atomic_save(list(signature=cfg$signature,par=opt$par,objective=opt$objective,time=Sys.time()),checkpoint)
    obj$fn(opt$par)
    parameters<-obj$env$parList()
    assert(max(abs(parameters$beta-opt$par[names(opt$par)=='beta']))<1e-10,'Saved fixed parameters do not match the optimizer.')
    # Save fitted parameter values before Hessian/exports, even if those fail.
    atomic_save(list(signature=cfg$signature,parameters=parameters,optimizer=opt,
      data_signature=input$data_signature),file.path(out,'optimized_parameters.rds'))
    phase(cfg,'Computing Hessian and fixed-effect uncertainty',task,list(fit_started=started))
    sdr<-withCallingHandlers(TMB::sdreport(obj,par.fixed=opt$par,getJointPrecision=FALSE,getReportCovariance=FALSE),
      warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
    g<-as.numeric(obj$gr(opt$par));cov<-sdr$cov.fixed
    ids<-which(names(opt$par)=='beta');beta<-setNames(opt$par[ids],colnames(built$data$X))
    V<-cov[ids,ids,drop=FALSE];dimnames(V)<-list(names(beta),names(beta))
    scaled<-if(all(is.finite(cov)))sqrt(max(0,as.numeric(crossprod(g,cov%*%g)))) else Inf
    good<-opt$convergence==0 && isTRUE(sdr$pdHess) && all(is.finite(V)) &&
      all(diag(V)>0) && is.finite(scaled) && scaled<.05 && is.finite(opt$objective)
    obj$fn(opt$par)
    report<-obj$report()
    extent<-max(diff(range(input$sites$x_km)),diff(range(input$sites$y_km)))
    checks<-data.frame(convergence_code=opt$convergence,pdHess=isTRUE(sdr$pdHess),
      max_abs_gradient=max(abs(g)),gradient_covariance_norm=scaled,numerical_checks_ok=good,
      logLik= -opt$objective,parameter_df=length(opt$par),AIC=2*opt$objective+2*length(opt$par),
      nobs=nrow(input$frame),spatial_range_km=report$range,spatial_sd=report$spatial_sd,
      range_greater_than_1_5_extent=report$range>1.5*extent,
      spatial_sd_below_0_001=report$spatial_sd<.001,
      AR1_rho=report$rho,Gamma_shape=report$shape)
    fit<-list(signature=cfg$signature,data_signature=input$data_signature,parameters=parameters,
      optimizer=opt,beta=beta,V=V,cov_fixed=cov,checks=checks,eta=as.numeric(report$eta),
      field_at_sites=report$spatial,group=group,variant=variant,
      fit_minutes=as.numeric(difftime(Sys.time(),started,units='mins')),warnings=unique(warnings))
    atomic_save(fit,cache)
    TMB::FreeADFun(obj);rm(obj,built,sdr,a);gc()
  }
  write_table(fit$checks,file.path(out,'fit_checks.csv'));writeLines(fit$warnings,file.path(out,'warnings.txt'))
  # Keep estimates from problematic fits for review but suppress comparative inference.
  write_table(fixed_table(fit$beta,fit$V),file.path(out,'fixed_effects.csv'))
  if(!isTRUE(fit$checks$numerical_checks_ok))return(data.frame(task=task,group=group,variant=variant,
    status='NUMERICAL_REVIEW_REQUIRED',cached=cached,fit_minutes=fit$fit_minutes,note='Inspect fit_checks.csv; excluded from comparisons.'))
  if(variant=='full') {
    phase(cfg,'Exporting effect comparisons',task)
    write_table(export_predictions(fit$beta,fit$V,input$years),file.path(out,'predictions_spatial.csv'))
    write_table(export_predictions(input$beta,input$V,input$years),file.path(out,'predictions_baseline.csv'))
    fields<-data.frame(SITE_ID=rep(input$sites$SITE_ID,times=length(input$years)),
      YEAR=rep(input$years,each=nrow(input$sites)),field=as.vector(fit$field_at_sites))
    write_table(fields,file.path(out,'annual_spatial_field.csv'))
    if(cfg$diagnostics) {
      phase(cfg,'Checking residual spatial dependence in both models',task)
      for(model in c('baseline','spatial')) {
        dest<-file.path(out,paste0('residual_',model,'.csv'))
        if(!file.exists(dest)) {
          rr<-residual_table(input,if(model=='baseline')input$baseline_eta else fit$eta,
            if(model=='baseline')exp(input$par$log_shape) else fit$checks$Gamma_shape)
          write_table(rr$site_year,file.path(out,paste0('site_year_',model,'.csv')))
          write_table(rr$summary,file.path(out,paste0('distribution_',model,'.csv')))
          sc<-moran_screen(rr$site_year,cfg$permutations,cfg$seed)
          write_table(sc,dest)
        }
      }
    }
  }
  data.frame(task=task,group=group,variant=variant,status='OK',cached=cached,
    fit_minutes=fit$fit_minutes,note='Spatial range/variance flags still require scientific review.')
}
combine_results <- function(cfg,status) {
  coefficients<-list();comparisons<-list();residuals<-list();lrts<-list();errors<-character()
  for(group in cfg$groups) {
    task<-paste0(group,'__full')
    if(!any(status$task==task & status$status=='OK'))next
    input<-readRDS(file.path(cfg$out,paste0('input_',group,'.rds')))
    out<-file.path(cfg$out,task);fit<-readRDS(file.path(out,'fit.rds'))
    b<-fixed_table(input$beta,input$V);s<-fixed_table(fit$beta,fit$V)
    b$model<-'Baseline';s$model<-'Spatial annual';b$group<-s$group<-group
    coefficients[[group]]<-rbind(b,s)
    comparison<-merge(b,s,by=c('group','term'),suffixes=c('_baseline','_spatial'))
    comparison$change_in_estimate<-comparison$estimate_spatial-comparison$estimate_baseline
    comparison$SE_ratio<-comparison$std.error_spatial/comparison$std.error_baseline
    comparison$CI_excludes_zero_baseline<-comparison$conf.low_baseline>0 | comparison$conf.high_baseline<0
    comparison$CI_excludes_zero_spatial<-comparison$conf.low_spatial>0 | comparison$conf.high_spatial<0
    write_table(comparison,file.path(out,'coefficient_comparison.csv'))
    comparisons[[group]]<-data.frame(group=group,n=nrow(input$frame),
      baseline_AIC=input$baseline_AIC,spatial_AIC=fit$checks$AIC,
      delta_AIC_spatial_minus_baseline=fit$checks$AIC-input$baseline_AIC,
      spatial_range_km=fit$checks$spatial_range_km,spatial_sd=fit$checks$spatial_sd,
      range_requires_review=fit$checks$range_greater_than_1_5_extent,
      spatial_sd_near_zero=fit$checks$spatial_sd_below_0_001,
      baseline_AR1_rho=input$par$rho_raw/sqrt(1+input$par$rho_raw^2),spatial_AR1_rho=fit$checks$AR1_rho)
    if(cfg$diagnostics)for(model in c('baseline','spatial')) {
      x<-utils::read.csv(file.path(out,paste0('residual_',model,'.csv')))
      for(scope in unique(x$scope)) {
        z<-x[x$scope==scope & is.finite(x$I),]
        residuals[[length(residuals)+1L]]<-data.frame(group=group,model=model,scope=scope,
          years=nrow(z),median_I=median(z$I),n_positive_q05=sum(z$I>0 & z$q_BH<.05,na.rm=TRUE))
      }
    }
    if(cfg$run_lrt)for(v in c('onset','offset')) {
      reduced_task<-paste0(group,'__without_',v,'_interaction')
      if(!any(status$task==reduced_task & status$status=='OK'))next
      red<-readRDS(file.path(cfg$out,reduced_task,'fit.rds'))
      stat<-2*(fit$checks$logLik-red$checks$logLik)
      df<-fit$checks$parameter_df-red$checks$parameter_df
      if(!identical(fit$data_signature,red$data_signature) || !is.finite(stat) || stat< -1e-5 || df!=1) {
        errors<-c(errors,paste('Invalid interaction comparison:',reduced_task));next
      }
      lrts[[length(lrts)+1L]]<-data.frame(group=group,interaction=v,Chisq=max(0,stat),df=df,
        p_LRT=pchisq(max(0,stat),df,lower.tail=FALSE))
    }
    rm(input,fit);gc()
  }
  if(length(coefficients)) {
    tab<-do.call(rbind,coefficients);write_table(tab,file.path(cfg$out,'all_fixed_effects.csv'))
    write_table(do.call(rbind,comparisons),file.path(cfg$out,'model_comparison.csv'))
    if(cfg$make_figures)tryCatch({
      p<-tab[grepl('year_decade:',tab$term),]
      p$effect<-ifelse(grepl('onset',p$term),'Onset advancement','Offset response')
      if(.Platform$OS.type=='windows')do.call(grDevices::windowsFonts,setNames(list(grDevices::windowsFont(cfg$font_family)),cfg$font_family))
      plot<-ggplot2::ggplot(p,ggplot2::aes(estimate,effect,colour=model))+
        ggplot2::geom_vline(xintercept=0,linetype='dotted',colour='grey45')+
        ggplot2::geom_errorbar(ggplot2::aes(xmin=conf.low,xmax=conf.high),width=.18,
          position=ggplot2::position_dodge(width=.5),orientation='y')+
        ggplot2::geom_point(position=ggplot2::position_dodge(width=.5),size=2)+
        ggplot2::facet_wrap(~group)+ggplot2::scale_colour_manual(values=c('Baseline'='#777777','Spatial annual'='#31666b'))+
        ggplot2::theme_classic(base_size=14,base_family=cfg$font_family)+
        ggplot2::theme(legend.position='top')+ggplot2::labs(x='Year × plasticity coefficient (per decade and SD)',y=NULL,colour=NULL)
      ggplot2::ggsave(file.path(cfg$out,'interaction_comparison.png'),plot,width=9,height=4,dpi=300)
    },error=function(e)errors<<-c(errors,paste('Figure:',conditionMessage(e))))
  }
  if(length(residuals))write_table(do.call(rbind,residuals),file.path(cfg$out,'residual_spatial_comparison.csv'))
  if(length(lrts))write_table(do.call(rbind,lrts),file.path(cfg$out,'interaction_LRTs.csv'))
  writeLines(errors,file.path(cfg$out,'export_errors.txt'))
  invisible(errors)
}
launch_child <- function(cfg,fun,args,log) {
  child<-callr::r_bg(function(engine,fun,args){
    e<-new.env(parent=globalenv());sys.source(engine,envir=e);do.call(e[[fun]],args)
  },args=list(engine=cfg$engine,fun=fun,args=args),stdout=log,stderr='2>&1',supervise=TRUE,wd=cfg$project_root)
  on.exit(if(child$is_alive())child$kill_tree(),add=TRUE)
  while(child$is_alive())Sys.sleep(1)
  child$get_result()
}
workflow <- function(cfg) {
  assert(identical(unname(tools::md5sum(cfg$baseline_files)),cfg$input_md5),'Baseline files changed after launch.')
  compile_engine(cfg)
  phase(cfg,'Preparing identical model inputs')
  if(!file.exists(file.path(cfg$out,'preparation_complete.rds'))) {
    launch_child(cfg,'prepare_worker',list(cfg=cfg),file.path(cfg$out,'preparation.log'))
    atomic_save(list(signature=cfg$signature),file.path(cfg$out,'preparation_complete.rds'))
  }
  if(cfg$prepare_only) {phase(cfg,'Preparation complete; no fits launched');return(cfg$out)}
  manifest<-expand.grid(group=cfg$groups,variant=if(cfg$run_lrt)c('full','without_onset_interaction','without_offset_interaction')else 'full',stringsAsFactors=FALSE)
  manifest$task<-paste(manifest$group,manifest$variant,sep='__')
  write_table(manifest,file.path(cfg$out,'task_manifest.csv'))
  results<-list()
  for(i in seq_len(nrow(manifest))) {
    task<-manifest[i,];out<-file.path(cfg$out,task$task);dir.create(out,showWarnings=FALSE)
    if(task$variant!='full' && !identical(results[[paste0(task$group,'__full')]]$status,'OK')) {
      s<-data.frame(task=task$task,group=task$group,variant=task$variant,status='SKIPPED_FULL_MODEL_FAILED',cached=FALSE,fit_minutes=NA_real_,note='Full model requires review.')
    } else {
      phase(cfg,paste('Starting task',i,'of',nrow(manifest)),task$task)
      log<-tempfile(paste0('worker_',format(Sys.time(),'%Y%m%d_%H%M%S'),'_'),tmpdir=out,fileext='.log')
      atomic_save(log,file.path(out,'latest_log.rds'))
      s<-tryCatch(launch_child(cfg,'fit_worker',list(cfg=cfg,group=task$group,variant=task$variant),log),
        error=function(e)data.frame(task=task$task,group=task$group,variant=task$variant,status='ERROR',cached=FALSE,
          fit_minutes=NA_real_,note=conditionMessage(e)))
    }
    results[[task$task]]<-s
    write_table(do.call(rbind,results),file.path(cfg$out,'workflow_status.csv'))
    # Partial comparison is updated after each full model.
    phase(cfg,'Updating model comparisons',task$task)
    errors<-tryCatch(combine_results(cfg,do.call(rbind,results)),error=function(e){
      writeLines(conditionMessage(e),file.path(cfg$out,'comparison_error.txt'));conditionMessage(e)})
  }
  ok<-all(vapply(results,function(s)s$status=='OK',logical(1))) && !length(errors)
  phase(cfg,if(ok)'Finished; review comparisons and diagnostics' else 'Finished with errors or checks requiring review')
  list(output_root=cfg$out,status=do.call(rbind,results),complete=ok)
}
status <- function(job=getOption('pheno.trends.spatial.active_job')) {
  assert(!is.null(job$process),'No spatial abundance job is registered.')
  read_safe<-function(path,fun)if(file.exists(path))tryCatch(suppressWarnings(fun(path)),error=function(e)NULL)else NULL
  p<-read_safe(job$cfg$progress_file,readRDS)
  s<-read_safe(file.path(job$cfg$out,'workflow_status.csv'),utils::read.csv)
  elapsed<-as.numeric(difftime(Sys.time(),job$started,units='mins'));alive<-job$process$is_alive()
  total<-if(job$cfg$prepare_only)0L else length(job$cfg$groups)*if(job$cfg$run_lrt)3L else 1L
  n<-if(is.null(s))0L else nrow(s);cached<-if(is.null(s))0L else sum(s$cached)
  failed<-if(is.null(s))0L else sum(s$status!='OK')
  message('[',format(Sys.time(),'%Y-%m-%d %H:%M:%S'),'] ',if(is.null(p))'Starting' else paste(p$task,p$phase),
    '\nTasks resolved: ',n,'/',total,' | cached: ',cached,' | failed/review: ',failed,
    ' | elapsed: ',round(elapsed,1),' min | worker running: ',alive)
  if(!is.null(p$fit_started))message('Current fit: ',round(as.numeric(difftime(Sys.time(),p$fit_started,units='mins')),1),
    ' min | completed likelihood evaluations: ',if(is.null(p$fn_evaluations))0L else p$fn_evaluations,
    ' | no reliable remaining-time estimate.')
  if(!is.null(p$updated))message('Last internal update: ',format(p$updated,'%H:%M:%S'),
    '. A likelihood/Hessian calculation can take several minutes between updates.')
  if(!alive)message('Process exit status: ',job$process$get_exit_status(),'. Results: ',job$cfg$out)
  invisible(list(progress=p,status=s,alive=alive))
}
watch <- function(job=getOption('pheno.trends.spatial.active_job')) {
  assert(!is.null(job$process),'No spatial abundance job is registered.')
  tryCatch({
    repeat {state<-status(job);if(!state$alive)break;Sys.sleep(job$cfg$monitor_seconds)}
    result<-job$process$get_result();message('Output: ',job$cfg$out);invisible(result)
  },interrupt=function(e){message('Display paused; FITTING CONTINUES. Keep R open. Resume with watch_spatial_trends().');invisible(job)})
}
stop_job <- function(job=getOption('pheno.trends.spatial.active_job')) {
  assert(!is.null(job$process),'No job is registered.');job$process$kill_tree()
  message('Worker stopped. Finished fits are retained; incomplete fits can restart from a saved parameter evaluation.')
  invisible(job)
}
resolve_baselines <- function(cfg) {
  if(!is.null(cfg$baseline_files)) {
    assert(all(cfg$groups%in%names(cfg$baseline_files)),'baseline_files must be named by voltinism group.')
    return(vapply(cfg$baseline_files[cfg$groups],normalizePath,character(1),winslash='/',mustWork=TRUE))
  }
  root<-cfg$baseline_root
  if(is.null(root))root<-file.path(cfg$project_root,'output','population_trends','spatial_plasticity_15km')
  assert(dir.exists(root),paste('Baseline folder missing:',root))
  candidates<-unique(c(root,list.dirs(root,recursive=FALSE,full.names=TRUE)))
  candidates<-candidates[vapply(candidates,function(d)all(file.exists(file.path(d,paste0(cfg$groups,'__full'),'fit.rds'))),logical(1))]
  assert(length(candidates)>0,'No complete saved baseline set found. Set baseline_root or named baseline_files.')
  modified<-vapply(candidates,function(d)max(as.numeric(file.info(file.path(d,paste0(cfg$groups,'__full'),'fit.rds'))$mtime)),numeric(1))
  chosen<-candidates[order(modified,candidates,decreasing=TRUE)][1]
  message('Baseline run: ',chosen)
  setNames(normalizePath(file.path(chosen,paste0(cfg$groups,'__full'),'fit.rds'),winslash='/',mustWork=TRUE),cfg$groups)
}
wait_existing <- function() {
  opts<-c('pheno.spatial.active_job','pheno.trends.active_job','pheno.trends.spatial.active_job',
    'pheno.offset.recovery.active_job','pheno.review.active_job','pheno.extract.active_job')
  jobs<-lapply(opts,getOption)
  for(nm in c('spatial_job','population_trend_job','spatial_trend_job'))
    if(exists(nm,envir=.GlobalEnv,inherits=FALSE))jobs[[length(jobs)+1L]]<-get(nm,envir=.GlobalEnv)
  for(job in jobs) {
    if(!is.list(job)||is.null(job$process))next
    while(isTRUE(tryCatch(job$process$is_alive(),error=function(e)FALSE))) {
      message('Waiting for the existing fitting worker to finish; no concurrent fits are launched.')
      Sys.sleep(30)
    }
  }
}
run <- function(cfg) {
  pkgs<-c('TMB','glmmTMB','Matrix','fmesher','callr','digest','here')
  if(cfg$make_figures)pkgs<-c(pkgs,'ggplot2')
  missing<-pkgs[!vapply(pkgs,requireNamespace,logical(1),quietly=TRUE)]
  assert(!length(missing),paste('Missing R packages:',paste(missing,collapse=', ')))
  assert(length(cfg$groups)>0 && all(cfg$groups%in%c('multivoltine','univoltine')) && !anyDuplicated(cfg$groups),'Invalid groups.')
  assert(cfg$mesh_cutoff_km>0 && cfg$monitor_seconds>=1 && cfg$monitor_seconds<=60 && cfg$permutations>=99,'Invalid mesh/monitor/permutation settings.')
  if(is.null(cfg$project_root))cfg$project_root<-here::here()
  cfg$project_root<-normalizePath(cfg$project_root,winslash='/',mustWork=TRUE)
  cfg$baseline_files<-resolve_baselines(cfg)
  if(!is.null(cfg$coordinates_file))cfg$coordinates_file<-normalizePath(cfg$coordinates_file,winslash='/',mustWork=TRUE)
  wait_existing()
  cfg$input_md5<-unname(tools::md5sum(cfg$baseline_files))
  if(is.null(cfg$output_root))cfg$output_root<-file.path(cfg$project_root,'output','population_trends','spatial_abundance_15km')
  env<-environment(run);nms<-ls(env,all.names=TRUE);nms<-nms[vapply(nms,function(n)is.function(get(n,env)),logical(1))]
  engine_text<-capture.output(dump(nms,file='',envir=env,control='all'))
  versions<-setNames(vapply(c('TMB','glmmTMB','Matrix','fmesher'),function(p)as.character(utils::packageVersion(p)),character(1)),c('TMB','glmmTMB','Matrix','fmesher'))
  cfg$signature<-digest::digest(list(code=engine_text,cpp=cpp_source(),baseline=cfg$input_md5,
    coords=if(is.null(cfg$coordinates_file))default_coordinates()else unname(tools::md5sum(cfg$coordinates_file)),
    groups=cfg$groups,cutoff=cfg$mesh_cutoff_km,versions=versions,R=as.character(getRversion()),
    optimizer=c(cfg$optimizer_iterations,cfg$optimizer_evaluations),permutations=cfg$permutations,seed=cfg$seed),algo='sha256')
  cfg$out<-file.path(cfg$output_root,paste0('run_',substr(cfg$signature,1,12)))
  dir.create(cfg$out,recursive=TRUE,showWarnings=FALSE);cfg$out<-normalizePath(cfg$out,winslash='/',mustWork=TRUE)
  cfg$engine<-file.path(cfg$out,'workflow_engine.R');writeLines(engine_text,cfg$engine)
  cfg$progress_file<-file.path(cfg$out,'progress.rds');atomic_save(cfg,file.path(cfg$out,'run_config.rds'))
  write_table(data.frame(group=cfg$groups,baseline=cfg$baseline_files,md5=cfg$input_md5),file.path(cfg$out,'baseline_sources.csv'))
  write_table(data.frame(package=names(versions),version=versions),file.path(cfg$out,'package_versions.csv'))
  # Replace invocation summaries only; never remove saved model fits or inputs.
  for(f in c('workflow_status.csv','progress.rds','task_manifest.csv','comparison_error.txt',
    'all_fixed_effects.csv','model_comparison.csv','residual_spatial_comparison.csv','interaction_LRTs.csv'))
    if(file.exists(file.path(cfg$out,f)))unlink(file.path(cfg$out,f))
  log<-tempfile(paste0('controller_',format(Sys.time(),'%Y%m%d_%H%M%S'),'_'),tmpdir=cfg$out,fileext='.log')
  process<-callr::r_bg(function(engine,cfg){e<-new.env(parent=globalenv());sys.source(engine,envir=e);e$workflow(cfg)},
    args=list(engine=cfg$engine,cfg=cfg),stdout=log,stderr='2>&1',supervise=TRUE,wd=cfg$project_root)
  job<-list(process=process,cfg=cfg,started=Sys.time(),log_file=log)
  options(pheno.trends.spatial.active_job=job)
  message('Started spatial abundance comparison. Keep R open. Detailed log: ',log)
  invisible(job)
}

cpp_source <- function() "// Gamma/log abundance model with the original independent species/site\n// intercepts and time slopes, population AR(1), and optional annual IID\n// Matern SPDE fields (nu=1). Full normalizing constants are retained.\n#include <TMB.hpp>\ntemplate<class Type>\nType objective_function<Type>::operator() () {\n  using namespace density;\n  DATA_VECTOR(y); DATA_MATRIX(X); DATA_VECTOR(time);\n  DATA_IVECTOR(species); DATA_IVECTOR(site); DATA_IVECTOR(year);\n  DATA_IVECTOR(previous); DATA_IVECTOR(gap);\n  DATA_SPARSE_MATRIX(A); DATA_SPARSE_MATRIX(M0);\n  DATA_SPARSE_MATRIX(M1); DATA_SPARSE_MATRIX(M2);\n  DATA_INTEGER(use_spatial);\n  PARAMETER_VECTOR(beta); PARAMETER_VECTOR(log_sd);\n  PARAMETER(rho_raw); PARAMETER(log_shape);\n  PARAMETER_MATRIX(b_species); PARAMETER_MATRIX(b_site);\n  PARAMETER_VECTOR(ar_obs);\n  PARAMETER(log_range); PARAMETER(log_spatial_sd);\n  PARAMETER_MATRIX(field);\n  Type nll=0;\n  Type rho=rho_raw/sqrt(Type(1)+rho_raw*rho_raw);\n  Type shape=exp(log_shape);\n  vector<Type> sd=exp(log_sd);\n  for(int j=0;j<b_species.rows();j++)\n    for(int k=0;k<2;k++) nll-=dnorm(b_species(j,k),Type(0),sd(k),true);\n  for(int j=0;j<b_site.rows();j++)\n    for(int k=0;k<2;k++) nll-=dnorm(b_site(j,k),Type(0),sd(k+2),true);\n  // Missing latent years integrate out analytically. Integer powers also\n  // handle negative rho correctly. The first observed year is stationary.\n  for(int i=0;i<y.size();i++) {\n    if(previous(i)<0) nll-=dnorm(ar_obs(i),Type(0),sd(4),true);\n    else {\n      Type r=1;\n      for(int k=0;k<gap(i);k++) r*=rho;\n      nll-=dnorm(ar_obs(i),r*ar_obs(previous(i)),sd(4)*sqrt(Type(1)-r*r),true);\n    }\n  }\n  matrix<Type> spatial(A.rows(),field.cols()); spatial.setZero();\n  Type range=exp(log_range), spatial_sd=exp(log_spatial_sd);\n  if(use_spatial) {\n    Type kappa=sqrt(Type(8))/range;\n    Type k2=kappa*kappa;\n    Type tau=Type(1)/(sqrt(Type(4)*Type(M_PI))*kappa*spatial_sd);\n    Eigen::SparseMatrix<Type> Q=k2*k2*M0+Type(2)*k2*M1+M2;\n    auto gmrf=GMRF(Q);\n    for(int t=0;t<field.cols();t++) {\n      vector<Type> w=field.col(t);\n      nll+=SCALE(gmrf,Type(1)/tau)(w);\n    }\n    spatial=A*field;\n  }\n  vector<Type> eta=X*beta;\n  for(int i=0;i<y.size();i++) {\n    eta(i)+=b_species(species(i),0)+time(i)*b_species(species(i),1)\n      +b_site(site(i),0)+time(i)*b_site(site(i),1)+ar_obs(i);\n    if(use_spatial) eta(i)+=spatial(site(i),year(i));\n    nll-=dgamma(y(i),shape,exp(eta(i))/shape,true);\n  }\n  REPORT(eta); REPORT(spatial); REPORT(rho); REPORT(shape);\n  REPORT(range); REPORT(spatial_sd);\n  return nll;\n}\n"
default_coordinates <- function() utils::read.csv(text="SITE_ID,bms_id,x_km,y_km\nDEBMS.100035,DEBMS,4278.445,3186.072\nDEBMS.100179,DEBMS,4115.373,3101.734\nDEBMS.100186,DEBMS,4277.869,2751.117\nDEBMS.100359,DEBMS,4207.141,2757.367\nDEBMS.100412,DEBMS,4116.708,3106.872\nDEBMS.101572,DEBMS,4411.492,2963.08\nDEBMS.104523,DEBMS,4426.267,3337.758\nDEBMS.104526,DEBMS,4234.448,3038.864\nDEBMS.104528,DEBMS,4310.455,3091.337\nDEBMS.104601,DEBMS,4428.192,2771.666\nDEBMS.104690,DEBMS,4452.616,3153.529\nDEBMS.104861,DEBMS,4489.985,3141.918\nDEBMS.105236,DEBMS,4094.27,2912.617\nDEBMS.108034,DEBMS,4244.073,2991.094\nDEBMS.108226,DEBMS,4459.269,3157.561\nDEBMS.108227,DEBMS,4459.659,3158.948\nDEBMS.108427,DEBMS,4560.268,3268.068\nDEBMS.108619,DEBMS,4456.843,3148.791\nDEBMS.108919,DEBMS,4427.836,3390.759\nDEBMS.109739,DEBMS,4224.654,2929.491\nDEBMS.119471,DEBMS,4094.84,3077.527\nDEBMS.14762,DEBMS,4347.613,3005.391\nDEBMS.14763,DEBMS,4208.447,3004.556\nDEBMS.14764,DEBMS,4233.263,3005.573\nDEBMS.14775,DEBMS,4267.981,2890.704\nDEBMS.14777,DEBMS,4219.146,2872.779\nDEBMS.14782,DEBMS,4123.227,3050.66\nDEBMS.14788,DEBMS,4444.417,3020.538\nDEBMS.14789,DEBMS,4384.943,2778.013\nDEBMS.14794,DEBMS,4207.359,3007.353\nDEBMS.14798,DEBMS,4222.058,2963.191\nDEBMS.14813,DEBMS,4434.41,3222.803\nDEBMS.14821,DEBMS,4257.563,2916.937\nDEBMS.14834,DEBMS,4393.0,2953.844\nDEBMS.14891,DEBMS,4456.839,3443.904\nDEBMS.14892,DEBMS,4078.016,2964.806\nDEBMS.14902,DEBMS,4158.948,2943.873\nDEBMS.14903,DEBMS,4226.062,2857.673\nDEBMS.14944,DEBMS,4208.961,2811.286\nDEBMS.14958,DEBMS,4341.573,3423.122\nDEBMS.14979,DEBMS,4482.432,3142.189\nDEBMS.14980,DEBMS,4490.269,3144.085\nDEBMS.14981,DEBMS,4488.592,3141.602\nDEBMS.14983,DEBMS,4489.056,3142.094\nDEBMS.14988,DEBMS,4573.478,3118.513\nDEBMS.14994,DEBMS,4515.393,3086.82\nDEBMS.14995,DEBMS,4532.343,3084.131\nDEBMS.14997,DEBMS,4540.571,3055.168\nDEBMS.15000,DEBMS,4395.107,3191.675\nDEBMS.15001,DEBMS,4396.347,3190.938\nDEBMS.15002,DEBMS,4406.972,3195.362\nDEBMS.15004,DEBMS,4394.214,3184.439\nDEBMS.15007,DEBMS,4455.081,3155.122\nDEBMS.15009,DEBMS,4462.211,3161.67\nDEBMS.15011,DEBMS,4436.699,3140.636\nDEBMS.15014,DEBMS,4476.391,3104.688\nDEBMS.15038,DEBMS,4257.504,2871.127\nDEBMS.15040,DEBMS,4403.06,2955.544\nDEBMS.15043,DEBMS,4305.142,2972.124\nDEBMS.15089,DEBMS,4221.507,2865.172\nDEBMS.15100,DEBMS,4123.507,3052.583\nDEBMS.15103,DEBMS,4182.161,2978.776\nDEBMS.15105,DEBMS,4187.78,2931.988\nDEBMS.15110,DEBMS,4283.656,2837.434\nDEBMS.15112,DEBMS,4267.802,2801.769\nDEBMS.15114,DEBMS,4277.356,2787.85\nDEBMS.15116,DEBMS,4312.464,2791.776\nDEBMS.15123,DEBMS,4139.066,2750.373\nDEBMS.15127,DEBMS,4440.266,3031.924\nDEBMS.15128,DEBMS,4441.004,3025.527\nDEBMS.15129,DEBMS,4444.878,3019.396\nDEBMS.15143,DEBMS,4305.005,2985.145\nDEBMS.15144,DEBMS,4361.593,2987.581\nDEBMS.15161,DEBMS,4390.287,2803.265\nDEBMS.15164,DEBMS,4417.289,2771.854\nDEBMS.15165,DEBMS,4439.123,2772.1\nDEBMS.15167,DEBMS,4518.905,2759.282\nDEBMS.15171,DEBMS,4525.972,2753.446\nDEBMS.15173,DEBMS,4232.519,3092.209\nDEBMS.15181,DEBMS,4221.927,2960.578\nDEBMS.15182,DEBMS,4221.717,2961.917\nDEBMS.15187,DEBMS,4361.859,3314.416\nDEBMS.15188,DEBMS,4238.931,3307.768\nDEBMS.15189,DEBMS,4237.939,3308.797\nDEBMS.15195,DEBMS,4332.028,3226.455\nDEBMS.15205,DEBMS,4399.622,2957.114\nDEBMS.15206,DEBMS,4349.045,2937.601\nDEBMS.15207,DEBMS,4372.392,2937.044\nDEBMS.15209,DEBMS,4461.758,2898.652\nDEBMS.15210,DEBMS,4450.598,2883.736\nDEBMS.15211,DEBMS,4476.113,2887.66\nDEBMS.17404,DEBMS,4337.411,3142.509\nDEBMS.17408,DEBMS,4547.462,3298.995\nDEBMS.17409,DEBMS,4554.187,3281.724\nDEBMS.17410,DEBMS,4535.862,3262.512\nDEBMS.17411,DEBMS,4541.923,3289.148\nDEBMS.17412,DEBMS,4541.923,3289.148\nDEBMS.17525,DEBMS,4377.667,3148.135\nDEBMS.17526,DEBMS,4262.477,3278.875\nDEBMS.17532,DEBMS,4337.72,3142.883\nDEBMS.17565,DEBMS,4392.987,3096.713\nDEBMS.17603,DEBMS,4138.772,2751.197\nDEBMS.17604,DEBMS,4534.738,3253.413\nDEBMS.17608,DEBMS,4311.481,2806.493\nDEBMS.17613,DEBMS,4262.648,2854.283\nDEBMS.17614,DEBMS,4594.014,3115.723\nDEBMS.17620,DEBMS,4265.768,3027.437\nDEBMS.18802,DEBMS,4538.388,3256.68\nDEBMS.20471,DEBMS,4233.187,3025.323\nDEBMS.20495,DEBMS,4314.55,3268.016\nDEBMS.20496,DEBMS,4626.01,3113.166\nDEBMS.20575,DEBMS,4562.703,3438.658\nDEBMS.20886,DEBMS,4342.917,3134.372\nDEBMS.20981,DEBMS,4140.302,2774.422\nDEBMS.21051,DEBMS,4457.621,3440.973\nDEBMS.21058,DEBMS,4227.056,2855.434\nDEBMS.21059,DEBMS,4226.331,2857.608\nDEBMS.21325,DEBMS,4548.124,3289.859\nDEBMS.21344,DEBMS,4196.501,2874.044\nDEBMS.21346,DEBMS,4538.911,3299.076\nDEBMS.21347,DEBMS,4218.171,2970.98\nDEBMS.21348,DEBMS,4205.3,2972.574\nDEBMS.21516,DEBMS,4458.05,3136.715\nDEBMS.21636,DEBMS,4455.69,3135.631\nDEBMS.21836,DEBMS,4438.159,3191.63\nDEBMS.21838,DEBMS,4438.17,3193.64\nDEBMS.21839,DEBMS,4437.924,3193.514\nDEBMS.22535,DEBMS,4524.276,3466.891\nDEBMS.22918,DEBMS,4438.999,2794.847\nDEBMS.23434,DEBMS,4534.993,3240.563\nDEBMS.23491,DEBMS,4167.809,2761.725\nDEBMS.23493,DEBMS,4122.72,3053.237\nDEBMS.23495,DEBMS,4118.08,3051.173\nDEBMS.23497,DEBMS,4535.445,3055.478\nDEBMS.29548,DEBMS,4442.199,3011.664\nDEBMS.29554,DEBMS,4547.487,3289.777\nDEBMS.29560,DEBMS,4367.915,2737.963\nDEBMS.29561,DEBMS,4368.321,2737.761\nDEBMS.29564,DEBMS,4528.838,3098.855\nDEBMS.29565,DEBMS,4084.417,2909.334\nDEBMS.29841,DEBMS,4284.481,2768.709\nDEBMS.30015,DEBMS,4403.605,3111.315\nDEBMS.30025,DEBMS,4159.0,2730.6\nDEBMS.30145,DEBMS,4185.077,2987.649\nDEBMS.31269,DEBMS,4457.269,2886.84\nDEBMS.31634,DEBMS,4215.169,2980.524\nDEBMS.31703,DEBMS,4477.28,3140.289\nDEBMS.31704,DEBMS,4497.216,3111.101\nDEBMS.33175,DEBMS,4207.5,3048.066\nDEBMS.33176,DEBMS,4214.036,2958.617\nDEBMS.33447,DEBMS,4266.075,3279.362\nDEBMS.33830,DEBMS,4421.989,2776.916\nDEBMS.33860,DEBMS,4202.514,2896.303\nDEBMS.33863,DEBMS,4513.275,2777.437\nDEBMS.35718,DEBMS,4417.111,2841.967\nDEBMS.35743,DEBMS,4636.093,3118.352\nDEBMS.35960,DEBMS,4538.936,3264.592\nDEBMS.35961,DEBMS,4539.877,3267.259\nDEBMS.35962,DEBMS,4542.177,3265.921\nDEBMS.36189,DEBMS,4203.236,2886.765\nDEBMS.38701,DEBMS,4595.141,3113.494\nDEBMS.66897,DEBMS,4462.238,3195.295\nDEBMS.66898,DEBMS,4439.59,3201.973\nDEBMS.67535,DEBMS,4574.043,3311.167\nDEBMS.67702,DEBMS,4251.29,2823.49\nDEBMS.67709,DEBMS,4578.262,3302.508\nDEBMS.67710,DEBMS,4578.798,3303.862\nDEBMS.68144,DEBMS,4573.297,3276.265\nDEBMS.68338,DEBMS,4224.8,2816.019\nDEBMS.68956,DEBMS,4115.37,3054.272\nDEBMS.69155,DEBMS,4153.454,2979.782\nDEBMS.69343,DEBMS,4599.908,3133.25\nDEBMS.69474,DEBMS,4624.262,3246.743\nDEBMS.69915,DEBMS,4168.506,3215.241\nDEBMS.69930,DEBMS,4094.322,3197.479\nDEBMS.69931,DEBMS,4092.467,3195.727\nDEBMS.69940,DEBMS,4154.45,3184.336\nDEBMS.69946,DEBMS,4066.714,3163.483\nDEBMS.69953,DEBMS,4131.14,3160.304\nDEBMS.69956,DEBMS,4091.363,3151.991\nDEBMS.69958,DEBMS,4137.003,3153.574\nDEBMS.69961,DEBMS,4059.745,3141.594\nDEBMS.69964,DEBMS,4119.595,3142.097\nDEBMS.69987,DEBMS,4113.147,3105.151\nDEBMS.69990,DEBMS,4182.198,3091.103\nDEBMS.70019,DEBMS,4121.143,3063.202\nDEBMS.70093,DEBMS,4456.622,3138.134\nDEBMS.70112,DEBMS,4197.094,3015.109\nDEBMS.70198,DEBMS,4304.098,3261.014\nDEBMS.70200,DEBMS,4516.058,3468.94\nDEBMS.70232,DEBMS,4311.4,3407.971\nDEBMS.70315,DEBMS,4115.602,2933.241\nDEBMS.70522,DEBMS,4327.353,3396.023\nDEBMS.70922,DEBMS,4214.036,2958.617\nDEBMS.70951,DEBMS,4436.676,2774.601\nDEBMS.70952,DEBMS,4437.086,2774.322\nDEBMS.71013,DEBMS,4587.063,3210.337\nDEBMS.71889,DEBMS,4339.583,2988.101\nDEBMS.72160,DEBMS,4206.809,2996.755\nDEBMS.72570,DEBMS,4206.951,2996.059\nDEBMS.72845,DEBMS,4361.067,2988.656\nDEBMS.75890,DEBMS,4090.045,2908.011\nDEBMS.77853,DEBMS,4203.854,2986.96\nDEBMS.77854,DEBMS,4206.367,2985.674\nDEBMS.78510,DEBMS,4205.005,2872.001\nDEBMS.78832,DEBMS,4389.428,3194.87\nDEBMS.79035,DEBMS,4441.211,3115.81\nDEBMS.79036,DEBMS,4442.932,3116.949\nDEBMS.79225,DEBMS,4378.228,3178.093\nDEBMS.79327,DEBMS,4323.921,2952.081\nDEBMS.80936,DEBMS,4494.725,3146.403\nDEBMS.80937,DEBMS,4497.551,3146.333\nDEBMS.80938,DEBMS,4492.663,3144.989\nDEBMS.80940,DEBMS,4497.254,3145.33\nDEBMS.80941,DEBMS,4495.594,3143.318\nDEBMS.80942,DEBMS,4497.549,3141.225\nDEBMS.84264,DEBMS,4554.136,3298.24\nDEBMS.85421,DEBMS,4554.11,3298.904\nDEBMS.86641,DEBMS,4461.214,3051.746\nDEBMS.87053,DEBMS,4455.14,3155.093\nDEBMS.87054,DEBMS,4455.198,3155.067\nDEBMS.87055,DEBMS,4455.038,3155.078\nDEBMS.87792,DEBMS,4332.74,3218.806\nDEBMS.87902,DEBMS,4370.186,2748.557\nDEBMS.87958,DEBMS,4125.195,3044.953\nDEBMS.88257,DEBMS,4118.506,2911.645\nDEBMS.88446,DEBMS,4573.744,3118.756\nDEBMS.88674,DEBMS,4248.072,3337.246\nDEBMS.88675,DEBMS,4248.828,3336.187\nDEBMS.88809,DEBMS,4328.723,3255.323\nDEBMS.88956,DEBMS,4216.121,2973.18\nDEBMS.89229,DEBMS,4105.762,2919.987\nDEBMS.92880,DEBMS,4540.048,3294.457\nDEBMS.93284,DEBMS,4329.851,3383.131\nDEBMS.93285,DEBMS,4519.664,3445.686\nDEBMS.93602,DEBMS,4381.556,2745.985\nDEBMS.93821,DEBMS,4557.396,3316.378\nDEBMS.93863,DEBMS,4226.615,3044.802\nDEBMS.93936,DEBMS,4332.74,3218.806\nDEBMS.94230,DEBMS,4376.571,3201.676\nDEBMS.96119,DEBMS,4224.966,2989.317\nDEBMS.96611,DEBMS,4160.161,2764.373\nDEBMS.96612,DEBMS,4160.401,2764.412\nDEBMS.96666,DEBMS,4155.999,2785.461\nDEBMS.97906,DEBMS,4563.748,3279.001\nDEBMS.98165,DEBMS,4520.543,2905.36\nES-CTBMS.1,ES-CTBMS,3749.693,2150.047\nES-CTBMS.10,ES-CTBMS,3685.802,2097.293\nES-CTBMS.101,ES-CTBMS,3824.065,1893.924\nES-CTBMS.106,ES-CTBMS,3683.752,2084.802\nES-CTBMS.107,ES-CTBMS,3653.739,2103.375\nES-CTBMS.108,ES-CTBMS,3699.601,2094.998\nES-CTBMS.109,ES-CTBMS,3646.285,2148.952\nES-CTBMS.11,ES-CTBMS,3686.899,2102.771\nES-CTBMS.110,ES-CTBMS,3683.54,2124.826\nES-CTBMS.112,ES-CTBMS,3679.635,2111.55\nES-CTBMS.113,ES-CTBMS,3715.83,2106.63\nES-CTBMS.114,ES-CTBMS,3646.523,2089.669\nES-CTBMS.115,ES-CTBMS,3701.149,2140.487\nES-CTBMS.116,ES-CTBMS,3609.799,2129.007\nES-CTBMS.117,ES-CTBMS,3588.661,2211.51\nES-CTBMS.12,ES-CTBMS,3690.233,2109.105\nES-CTBMS.120,ES-CTBMS,3650.18,2060.298\nES-CTBMS.122,ES-CTBMS,3674.663,2100.472\nES-CTBMS.123,ES-CTBMS,3662.136,2114.302\nES-CTBMS.125,ES-CTBMS,3663.678,2069.419\nES-CTBMS.126,ES-CTBMS,3631.686,2202.335\nES-CTBMS.127,ES-CTBMS,3648.068,2137.508\nES-CTBMS.128,ES-CTBMS,3685.938,2087.37\nES-CTBMS.129,ES-CTBMS,3675.657,2089.078\nES-CTBMS.13,ES-CTBMS,3698.106,2095.539\nES-CTBMS.130,ES-CTBMS,3745.42,2169.498\nES-CTBMS.147,ES-CTBMS,3666.018,2075.78\nES-CTBMS.155,ES-CTBMS,3689.863,1873.623\nES-CTBMS.168,ES-CTBMS,3616.971,2179.87\nES-CTBMS.18,ES-CTBMS,3539.917,2097.938\nES-CTBMS.19,ES-CTBMS,3695.831,2103.731\nES-CTBMS.20,ES-CTBMS,3690.799,2106.803\nES-CTBMS.21,ES-CTBMS,3658.256,2072.122\nES-CTBMS.24,ES-CTBMS,3652.597,2098.084\nES-CTBMS.26,ES-CTBMS,3639.148,2057.376\nES-CTBMS.28,ES-CTBMS,3682.023,2107.432\nES-CTBMS.29,ES-CTBMS,3680.836,2094.141\nES-CTBMS.33,ES-CTBMS,3691.484,2088.863\nES-CTBMS.34,ES-CTBMS,3671.631,2075.039\nES-CTBMS.36,ES-CTBMS,3636.817,2062.123\nES-CTBMS.38,ES-CTBMS,3591.906,2045.127\nES-CTBMS.40,ES-CTBMS,3649.189,2117.348\nES-CTBMS.41,ES-CTBMS,3546.268,2090.904\nES-CTBMS.42,ES-CTBMS,3648.149,2137.662\nES-CTBMS.43,ES-CTBMS,3605.904,2117.4\nES-CTBMS.45,ES-CTBMS,3635.684,2057.541\nES-CTBMS.47,ES-CTBMS,3516.138,2084.86\nES-CTBMS.48,ES-CTBMS,3525.618,2065.385\nES-CTBMS.5,ES-CTBMS,3729.096,2168.814\nES-CTBMS.51,ES-CTBMS,3576.525,2067.464\nES-CTBMS.52,ES-CTBMS,3594.583,2045.424\nES-CTBMS.53,ES-CTBMS,3681.801,2103.007\nES-CTBMS.55,ES-CTBMS,3640.122,2149.074\nES-CTBMS.58,ES-CTBMS,3610.573,2067.816\nES-CTBMS.59,ES-CTBMS,3749.764,2147.274\nES-CTBMS.60,ES-CTBMS,3803.356,1896.425\nES-CTBMS.61,ES-CTBMS,3828.584,1891.603\nES-CTBMS.66,ES-CTBMS,3550.734,2157.551\nES-CTBMS.67,ES-CTBMS,3538.296,1997.916\nES-CTBMS.68,ES-CTBMS,3672.45,2076.335\nES-CTBMS.69,ES-CTBMS,3679.291,2077.969\nES-CTBMS.70,ES-CTBMS,3714.145,2156.124\nES-CTBMS.72,ES-CTBMS,3649.982,2155.696\nES-CTBMS.75,ES-CTBMS,3661.717,2087.82\nES-CTBMS.76,ES-CTBMS,3659.884,2079.704\nES-CTBMS.77,ES-CTBMS,3727.2,2127.438\nES-CTBMS.79,ES-CTBMS,3646.183,2103.54\nES-CTBMS.8,ES-CTBMS,3656.035,2069.317\nES-CTBMS.80,ES-CTBMS,3656.234,2099.086\nES-CTBMS.82,ES-CTBMS,3502.045,2030.021\nES-CTBMS.85,ES-CTBMS,3583.767,2177.651\nES-CTBMS.86,ES-CTBMS,3584.788,2208.067\nES-CTBMS.87,ES-CTBMS,3565.34,2163.883\nES-CTBMS.88,ES-CTBMS,3695.738,2094.034\nES-CTBMS.89,ES-CTBMS,3709.993,2090.487\nES-CTBMS.9,ES-CTBMS,3700.585,2146.389\nES-CTBMS.90,ES-CTBMS,3651.995,2179.052\nES-CTBMS.91,ES-CTBMS,3628.306,2207.269\nES-CTBMS.92,ES-CTBMS,3620.957,2194.382\nES-CTBMS.93,ES-CTBMS,3620.455,2203.956\nES-CTBMS.94,ES-CTBMS,3548.088,2068.804\nES-CTBMS.95,ES-CTBMS,3659.746,2075.311\nES-CTBMS.96,ES-CTBMS,3619.867,2189.843\nES-CTBMS.98,ES-CTBMS,3622.993,2195.038\nES-ZEBMS.PV_16,ES-ZEBMS,3262.005,2340.303\nES-ZEBMS.PV_22,ES-ZEBMS,3238.142,2336.062\nES-ZEBMS.PV_23,ES-ZEBMS,3287.795,2335.597\nES-ZEBMS.PV_24,ES-ZEBMS,3289.44,2339.192\nES-ZEBMS.PV_6,ES-ZEBMS,3304.753,2261.068\nESBMS.223573,ESBMS,2804.609,2359.08\nESBMS.223863,ESBMS,3114.086,2346.821\nESBMS.223865,ESBMS,3115.044,2351.427\nESBMS.223875,ESBMS,3163.617,2043.373\nFIBMS.107,FIBMS,5117.917,4217.79\nFIBMS.232,FIBMS,5051.24,4234.058\nFIBMS.246,FIBMS,5149.249,4250.343\nFIBMS.249,FIBMS,5016.026,4218.577\nFIBMS.255,FIBMS,5172.638,4315.045\nFIBMS.258,FIBMS,5000.635,4292.71\nFIBMS.261,FIBMS,5316.258,4361.813\nFIBMS.262,FIBMS,5295.702,4517.957\nFIBMS.266,FIBMS,5038.28,4288.663\nFIBMS.267,FIBMS,5062.021,4261.064\nFIBMS.269,FIBMS,5001.71,4290.083\nFIBMS.271,FIBMS,5077.58,4211.704\nFIBMS.280,FIBMS,5044.626,4210.147\nFIBMS.284,FIBMS,5072.046,4354.261\nFIBMS.285,FIBMS,4997.702,4210.522\nFIBMS.286,FIBMS,5182.804,4235.967\nFIBMS.288,FIBMS,5052.066,4292.802\nFIBMS.290,FIBMS,5063.212,4702.835\nFIBMS.292,FIBMS,5242.778,4286.64\nFIBMS.293,FIBMS,5307.113,4510.849\nFIBMS.295,FIBMS,5355.008,4452.273\nFIBMS.298,FIBMS,5086.74,4725.775\nFIBMS.33,FIBMS,4913.142,4489.754\nFIBMS.51,FIBMS,4911.374,4493.12\nFIBMS.55,FIBMS,4991.101,4193.95\nFIBMS.57,FIBMS,5235.552,4273.376\nFIBMS.58,FIBMS,5305.286,4376.492\nFIBMS.60,FIBMS,5328.441,4499.112\nFIBMS.61,FIBMS,5366.783,4485.773\nFIBMS.63,FIBMS,5057.257,4404.896\nFIBMS.68,FIBMS,5201.769,4444.363\nFIBMS.72,FIBMS,5204.904,4357.909\nFIBMS.77,FIBMS,5202.238,4350.06\nFIBMS.80,FIBMS,5038.303,4288.692\nFIBMS.89,FIBMS,5200.794,4362.319\nIEBMS.C03,IEBMS,3069.718,3342.762\nIEBMS.C16,IEBMS,3013.52,3359.506\nIEBMS.C23,IEBMS,2986.411,3366.102\nIEBMS.C29,IEBMS,2990.289,3342.964\nIEBMS.C38,IEBMS,3016.378,3375.238\nIEBMS.CE12,IEBMS,3060.02,3487.955\nIEBMS.D01,IEBMS,3258.193,3494.292\nIEBMS.D08,IEBMS,3245.637,3476.439\nIEBMS.D13,IEBMS,3245.998,3472.97\nIEBMS.D19,IEBMS,3255.465,3467.745\nIEBMS.DL02,IEBMS,3160.52,3669.993\nIEBMS.DL03,IEBMS,3201.029,3663.775\nIEBMS.DL05,IEBMS,3158.476,3686.569\nIEBMS.DL06,IEBMS,3144.422,3657.098\nIEBMS.DL07,IEBMS,3194.535,3702.071\nIEBMS.G01,IEBMS,3060.071,3533.048\nIEBMS.KE01,IEBMS,3200.633,3484.029\nIEBMS.KE02,IEBMS,3201.526,3484.294\nIEBMS.KE04,IEBMS,3233.303,3487.294\nIEBMS.KE06,IEBMS,3229.509,3488.928\nIEBMS.LD01,IEBMS,3153.216,3527.803\nIEBMS.LM01,IEBMS,3150.868,3616.008\nIEBMS.S01,IEBMS,3127.543,3626.174\nIEBMS.T16,IEBMS,3124.119,3463.981\nIEBMS.W03,IEBMS,3161.678,3372.53\nIEBMS.WW04,IEBMS,3234.353,3444.99\nIEBMS.WW07,IEBMS,3232.594,3428.949\nIEBMS.WX09,IEBMS,3194.802,3368.413\nIEBMS.WX10,IEBMS,3207.889,3374.437\nLUBMS.EBMS:Luxembourg:70_S2010,LUBMS,4024.424,2945.903\nLUBMS.EBMS:Luxembourg:70_S2011,LUBMS,4024.196,2945.994\nNLBMS.100,NLBMS,4035.744,3224.773\nNLBMS.1002,NLBMS,3948.694,3169.695\nNLBMS.1004,NLBMS,3961.918,3158.485\nNLBMS.1007,NLBMS,4068.837,3229.813\nNLBMS.1036,NLBMS,4084.022,3219.0\nNLBMS.1047,NLBMS,4016.773,3087.525\nNLBMS.105,NLBMS,3994.5,3212.946\nNLBMS.1058,NLBMS,3999.177,3157.845\nNLBMS.1078,NLBMS,3957.514,3293.44\nNLBMS.1086,NLBMS,4036.392,3224.756\nNLBMS.1087,NLBMS,4018.757,3219.756\nNLBMS.1091,NLBMS,4070.353,3370.654\nNLBMS.1095,NLBMS,4070.413,3370.514\nNLBMS.1102,NLBMS,3890.912,3169.122\nNLBMS.1112,NLBMS,3965.264,3157.905\nNLBMS.1116,NLBMS,4114.896,3257.37\nNLBMS.1119,NLBMS,4052.739,3312.797\nNLBMS.1122,NLBMS,4054.384,3228.688\nNLBMS.1124,NLBMS,3957.715,3240.034\nNLBMS.1125,NLBMS,3966.652,3337.199\nNLBMS.1127,NLBMS,4052.958,3312.995\nNLBMS.1156,NLBMS,4033.301,3089.974\nNLBMS.1171,NLBMS,3973.703,3341.457\nNLBMS.1222,NLBMS,4116.302,3287.973\nNLBMS.1224,NLBMS,4116.511,3288.351\nNLBMS.1228,NLBMS,4116.288,3288.023\nNLBMS.1231,NLBMS,4118.879,3287.603\nNLBMS.1232,NLBMS,4118.764,3287.645\nNLBMS.126,NLBMS,4071.755,3211.721\nNLBMS.1283,NLBMS,4044.14,3183.591\nNLBMS.1293,NLBMS,4023.854,3182.997\nNLBMS.1303,NLBMS,3988.712,3249.194\nNLBMS.1321,NLBMS,4023.265,3090.794\nNLBMS.1339,NLBMS,3969.995,3327.368\nNLBMS.1387,NLBMS,3948.337,3263.741\nNLBMS.1391,NLBMS,3947.668,3262.667\nNLBMS.1395,NLBMS,3947.288,3261.847\nNLBMS.1398,NLBMS,3949.941,3261.472\nNLBMS.140,NLBMS,4094.836,3308.954\nNLBMS.1401,NLBMS,3951.203,3261.39\nNLBMS.1411,NLBMS,3888.156,3181.689\nNLBMS.1426,NLBMS,3951.645,3266.936\nNLBMS.1430,NLBMS,3950.157,3267.805\nNLBMS.1434,NLBMS,4029.862,3249.008\nNLBMS.1440,NLBMS,4016.869,3091.368\nNLBMS.1443,NLBMS,4022.906,3217.094\nNLBMS.1462,NLBMS,4016.99,3087.222\nNLBMS.1484,NLBMS,3941.311,3232.657\nNLBMS.1487,NLBMS,3957.372,3287.096\nNLBMS.1490,NLBMS,3959.519,3294.429\nNLBMS.150,NLBMS,3987.785,3156.499\nNLBMS.152,NLBMS,3983.718,3145.628\nNLBMS.1528,NLBMS,3991.882,3216.18\nNLBMS.1557,NLBMS,3958.1,3295.56\nNLBMS.1580,NLBMS,3941.258,3232.018\nNLBMS.1596,NLBMS,3987.219,3236.603\nNLBMS.160,NLBMS,3994.964,3215.143\nNLBMS.1604,NLBMS,3968.195,3255.822\nNLBMS.1609,NLBMS,3967.94,3255.709\nNLBMS.1613,NLBMS,3966.275,3326.29\nNLBMS.173,NLBMS,4037.68,3096.552\nNLBMS.176,NLBMS,4042.55,3191.213\nNLBMS.1778,NLBMS,4034.929,3218.81\nNLBMS.1784,NLBMS,4068.113,3344.787\nNLBMS.1785,NLBMS,4068.062,3344.942\nNLBMS.1795,NLBMS,4014.438,3166.723\nNLBMS.1798,NLBMS,4014.771,3164.887\nNLBMS.1837,NLBMS,4070.345,3220.504\nNLBMS.186,NLBMS,3961.479,3291.161\nNLBMS.1912,NLBMS,4059.729,3285.936\nNLBMS.1922,NLBMS,3964.352,3314.499\nNLBMS.1932,NLBMS,3977.138,3181.975\nNLBMS.1963,NLBMS,4079.661,3370.308\nNLBMS.1968,NLBMS,4080.398,3369.41\nNLBMS.1972,NLBMS,4026.221,3212.082\nNLBMS.1976,NLBMS,4081.123,3365.131\nNLBMS.1979,NLBMS,4078.796,3366.514\nNLBMS.1991,NLBMS,3964.263,3180.559\nNLBMS.20,NLBMS,3986.596,3245.243\nNLBMS.201,NLBMS,4083.771,3338.506\nNLBMS.209,NLBMS,3992.12,3185.797\nNLBMS.212,NLBMS,4030.955,3138.313\nNLBMS.220,NLBMS,3978.759,3204.928\nNLBMS.2244,NLBMS,3994.507,3212.911\nNLBMS.2251,NLBMS,4089.606,3288.235\nNLBMS.2254,NLBMS,4003.593,3161.385\nNLBMS.2258,NLBMS,4003.555,3161.33\nNLBMS.2278,NLBMS,4087.087,3255.629\nNLBMS.2298,NLBMS,3987.889,3235.846\nNLBMS.242,NLBMS,3937.231,3204.772\nNLBMS.2470,NLBMS,4023.095,3165.273\nNLBMS.248,NLBMS,4025.526,3255.715\nNLBMS.2496,NLBMS,4027.736,3216.832\nNLBMS.2510,NLBMS,4035.433,3194.807\nNLBMS.253,NLBMS,4068.511,3287.967\nNLBMS.2531,NLBMS,4062.192,3345.193\nNLBMS.2535,NLBMS,4023.573,3217.535\nNLBMS.2546,NLBMS,4041.114,3167.122\nNLBMS.2549,NLBMS,3998.178,3187.941\nNLBMS.2552,NLBMS,4001.265,3186.579\nNLBMS.257,NLBMS,4044.604,3355.937\nNLBMS.2601,NLBMS,3982.338,3172.214\nNLBMS.2604,NLBMS,3976.451,3170.368\nNLBMS.2612,NLBMS,3979.812,3174.667\nNLBMS.2614,NLBMS,3982.003,3171.091\nNLBMS.262,NLBMS,4095.759,3248.605\nNLBMS.2648,NLBMS,4069.025,3286.953\nNLBMS.265,NLBMS,4020.1,3141.766\nNLBMS.2658,NLBMS,4032.466,3228.831\nNLBMS.2663,NLBMS,4036.971,3224.827\nNLBMS.2665,NLBMS,4036.79,3225.126\nNLBMS.2668,NLBMS,4022.082,3215.776\nNLBMS.269,NLBMS,4018.628,3137.603\nNLBMS.273,NLBMS,4018.554,3137.386\nNLBMS.2789,NLBMS,4113.277,3253.102\nNLBMS.279,NLBMS,4018.21,3137.54\nNLBMS.2818,NLBMS,4017.232,3093.554\nNLBMS.2824,NLBMS,3941.745,3206.418\nNLBMS.2832,NLBMS,4038.842,3098.063\nNLBMS.2855,NLBMS,4049.634,3241.45\nNLBMS.288,NLBMS,4120.978,3252.97\nNLBMS.2880,NLBMS,4038.447,3229.545\nNLBMS.2911,NLBMS,4130.674,3333.474\nNLBMS.292,NLBMS,4090.328,3277.925\nNLBMS.2936,NLBMS,4018.488,3320.678\nNLBMS.2986,NLBMS,4022.058,3167.153\nNLBMS.2994,NLBMS,4026.502,3162.917\nNLBMS.3001,NLBMS,4115.854,3323.986\nNLBMS.3005,NLBMS,3979.158,3264.622\nNLBMS.3009,NLBMS,4025.265,3163.762\nNLBMS.3012,NLBMS,4024.466,3164.389\nNLBMS.3017,NLBMS,4033.18,3165.89\nNLBMS.3020,NLBMS,3937.196,3222.113\nNLBMS.3023,NLBMS,3937.788,3222.293\nNLBMS.3037,NLBMS,4032.11,3228.824\nNLBMS.3093,NLBMS,4006.275,3155.159\nNLBMS.3095,NLBMS,3936.555,3235.998\nNLBMS.3099,NLBMS,3937.033,3236.386\nNLBMS.3102,NLBMS,3937.302,3236.632\nNLBMS.3118,NLBMS,3968.3,3339.948\nNLBMS.3129,NLBMS,3911.461,3214.765\nNLBMS.3131,NLBMS,3939.782,3237.763\nNLBMS.3149,NLBMS,3939.684,3209.844\nNLBMS.3162,NLBMS,3986.512,3241.053\nNLBMS.317,NLBMS,4049.531,3172.891\nNLBMS.3211,NLBMS,3939.554,3210.637\nNLBMS.3213,NLBMS,3932.231,3208.725\nNLBMS.3221,NLBMS,3941.107,3218.527\nNLBMS.3245,NLBMS,4022.849,3165.51\nNLBMS.326,NLBMS,3968.609,3327.935\nNLBMS.3290,NLBMS,3889.64,3141.22\nNLBMS.3292,NLBMS,4018.312,3089.271\nNLBMS.332,NLBMS,3999.015,3210.135\nNLBMS.3328,NLBMS,4042.691,3216.461\nNLBMS.3332,NLBMS,4084.613,3217.799\nNLBMS.3338,NLBMS,4068.71,3287.328\nNLBMS.3340,NLBMS,4100.145,3231.485\nNLBMS.335,NLBMS,4072.658,3264.374\nNLBMS.3360,NLBMS,3958.715,3290.031\nNLBMS.3394,NLBMS,3973.647,3340.002\nNLBMS.3434,NLBMS,3947.126,3220.695\nNLBMS.346,NLBMS,4021.301,3140.688\nNLBMS.3468,NLBMS,4019.721,3130.082\nNLBMS.353,NLBMS,4021.622,3252.383\nNLBMS.3573,NLBMS,4023.475,3263.076\nNLBMS.3574,NLBMS,4022.556,3263.143\nNLBMS.36,NLBMS,4018.759,3220.007\nNLBMS.3649,NLBMS,3959.331,3295.326\nNLBMS.3660,NLBMS,4049.708,3291.717\nNLBMS.3677,NLBMS,3973.049,3338.888\nNLBMS.3718,NLBMS,3940.511,3219.883\nNLBMS.3741,NLBMS,3968.22,3260.431\nNLBMS.3747,NLBMS,4050.633,3343.318\nNLBMS.3748,NLBMS,4054.36,3272.877\nNLBMS.3750,NLBMS,4009.209,3245.177\nNLBMS.3757,NLBMS,4036.663,3226.168\nNLBMS.3760,NLBMS,4023.933,3089.149\nNLBMS.3767,NLBMS,3992.366,3216.156\nNLBMS.3788,NLBMS,3948.138,3220.653\nNLBMS.3790,NLBMS,4133.829,3326.016\nNLBMS.3792,NLBMS,4004.863,3151.628\nNLBMS.3804,NLBMS,4092.345,3213.394\nNLBMS.383,NLBMS,3961.101,3298.478\nNLBMS.3847,NLBMS,4061.214,3325.403\nNLBMS.385,NLBMS,3957.472,3281.692\nNLBMS.3855,NLBMS,3945.599,3221.109\nNLBMS.389,NLBMS,3958.326,3289.981\nNLBMS.391,NLBMS,4063.38,3334.764\nNLBMS.393,NLBMS,3933.974,3241.405\nNLBMS.3945,NLBMS,4020.514,3141.721\nNLBMS.3947,NLBMS,4018.656,3137.601\nNLBMS.3950,NLBMS,4018.597,3137.362\nNLBMS.3951,NLBMS,4018.447,3233.019\nNLBMS.3953,NLBMS,3968.98,3340.711\nNLBMS.396,NLBMS,3956.848,3289.443\nNLBMS.3975,NLBMS,4045.345,3165.111\nNLBMS.399,NLBMS,4019.184,3212.896\nNLBMS.4018,NLBMS,3994.037,3246.902\nNLBMS.4020,NLBMS,3989.576,3251.91\nNLBMS.4022,NLBMS,3992.325,3250.561\nNLBMS.4023,NLBMS,3991.096,3249.501\nNLBMS.407,NLBMS,4036.413,3224.809\nNLBMS.4089,NLBMS,3946.305,3257.075\nNLBMS.4126,NLBMS,3948.504,3169.39\nNLBMS.4131,NLBMS,3992.352,3248.753\nNLBMS.4142,NLBMS,3948.753,3259.148\nNLBMS.4153,NLBMS,4068.487,3287.934\nNLBMS.4159,NLBMS,3942.285,3170.317\nNLBMS.4161,NLBMS,3942.384,3170.109\nNLBMS.4191,NLBMS,3996.942,3252.263\nNLBMS.4199,NLBMS,3964.468,3262.292\nNLBMS.420,NLBMS,4035.7,3224.712\nNLBMS.4205,NLBMS,3968.314,3261.865\nNLBMS.423,NLBMS,4033.709,3223.409\nNLBMS.4240,NLBMS,4035.436,3199.582\nNLBMS.4245,NLBMS,4102.795,3292.821\nNLBMS.4254,NLBMS,4013.396,3153.503\nNLBMS.4257,NLBMS,4013.793,3152.695\nNLBMS.426,NLBMS,3966.942,3327.234\nNLBMS.4260,NLBMS,4010.54,3146.544\nNLBMS.4266,NLBMS,4010.438,3146.586\nNLBMS.4278,NLBMS,3957.122,3286.735\nNLBMS.4284,NLBMS,4018.413,3156.208\nNLBMS.4287,NLBMS,4018.352,3156.516\nNLBMS.432,NLBMS,3958.534,3286.663\nNLBMS.4325,NLBMS,4024.095,3262.383\nNLBMS.4326,NLBMS,4126.61,3346.426\nNLBMS.4333,NLBMS,4016.931,3274.495\nNLBMS.435,NLBMS,3959.6,3289.488\nNLBMS.4352,NLBMS,4049.052,3124.113\nNLBMS.4374,NLBMS,4075.077,3359.155\nNLBMS.4377,NLBMS,4074.749,3359.201\nNLBMS.4382,NLBMS,4018.606,3144.249\nNLBMS.4384,NLBMS,4017.741,3144.813\nNLBMS.4396,NLBMS,3966.28,3327.026\nNLBMS.4398,NLBMS,4070.557,3210.344\nNLBMS.4401,NLBMS,4106.941,3231.384\nNLBMS.4416,NLBMS,3965.525,3321.941\nNLBMS.4469,NLBMS,3919.758,3207.632\nNLBMS.4470,NLBMS,3987.645,3250.847\nNLBMS.4476,NLBMS,3937.663,3221.831\nNLBMS.4480,NLBMS,3963.367,3288.329\nNLBMS.4494,NLBMS,4060.814,3249.91\nNLBMS.4496,NLBMS,4056.544,3251.501\nNLBMS.4506,NLBMS,3917.353,3225.036\nNLBMS.4510,NLBMS,4107.541,3347.893\nNLBMS.4516,NLBMS,4016.596,3089.387\nNLBMS.4525,NLBMS,4070.043,3257.04\nNLBMS.4527,NLBMS,4070.17,3257.345\nNLBMS.4541,NLBMS,3953.042,3177.925\nNLBMS.4545,NLBMS,3983.349,3185.767\nNLBMS.4546,NLBMS,3983.061,3185.245\nNLBMS.4548,NLBMS,3982.71,3185.263\nNLBMS.4553,NLBMS,4037.269,3192.88\nNLBMS.4554,NLBMS,4037.547,3193.345\nNLBMS.4557,NLBMS,3986.441,3235.682\nNLBMS.4559,NLBMS,4023.026,3214.691\nNLBMS.4585,NLBMS,3993.207,3170.332\nNLBMS.4586,NLBMS,4110.237,3245.085\nNLBMS.4587,NLBMS,4110.076,3244.561\nNLBMS.4588,NLBMS,3933.016,3239.259\nNLBMS.4589,NLBMS,3993.257,3170.886\nNLBMS.459,NLBMS,4047.326,3178.671\nNLBMS.4609,NLBMS,3986.031,3185.309\nNLBMS.4610,NLBMS,3985.231,3183.432\nNLBMS.4612,NLBMS,3988.023,3171.884\nNLBMS.4616,NLBMS,3953.203,3267.131\nNLBMS.4619,NLBMS,3962.492,3176.393\nNLBMS.4621,NLBMS,3965.451,3186.023\nNLBMS.4622,NLBMS,3962.204,3174.461\nNLBMS.4624,NLBMS,3959.879,3178.424\nNLBMS.4627,NLBMS,3939.475,3169.773\nNLBMS.463,NLBMS,3979.672,3236.802\nNLBMS.4641,NLBMS,3944.266,3209.58\nNLBMS.469,NLBMS,3966.18,3178.198\nNLBMS.4701,NLBMS,4128.988,3333.503\nNLBMS.4708,NLBMS,3985.208,3240.76\nNLBMS.4709,NLBMS,4002.913,3239.314\nNLBMS.4712,NLBMS,4039.048,3197.806\nNLBMS.4714,NLBMS,4002.842,3239.735\nNLBMS.4722,NLBMS,3953.285,3260.079\nNLBMS.4726,NLBMS,4004.12,3238.196\nNLBMS.4727,NLBMS,3937.863,3243.819\nNLBMS.4728,NLBMS,4036.528,3133.52\nNLBMS.4729,NLBMS,4030.376,3120.814\nNLBMS.4735,NLBMS,4020.919,3167.239\nNLBMS.4736,NLBMS,4039.981,3186.608\nNLBMS.4739,NLBMS,3991.207,3249.198\nNLBMS.4740,NLBMS,3955.834,3178.006\nNLBMS.4742,NLBMS,3933.814,3180.422\nNLBMS.4745,NLBMS,4061.099,3270.042\nNLBMS.4753,NLBMS,4097.287,3331.848\nNLBMS.4755,NLBMS,4037.009,3132.241\nNLBMS.4763,NLBMS,4034.109,3087.36\nNLBMS.4764,NLBMS,4033.575,3087.375\nNLBMS.4767,NLBMS,3960.869,3179.822\nNLBMS.477,NLBMS,4053.611,3231.913\nNLBMS.4771,NLBMS,4090.539,3213.05\nNLBMS.4772,NLBMS,4033.847,3087.117\nNLBMS.4775,NLBMS,4113.119,3245.105\nNLBMS.4778,NLBMS,3983.4,3185.727\nNLBMS.4782,NLBMS,3896.55,3165.918\nNLBMS.4783,NLBMS,4105.813,3318.399\nNLBMS.4789,NLBMS,4052.298,3158.888\nNLBMS.4792,NLBMS,4027.633,3216.807\nNLBMS.494,NLBMS,4015.97,3181.894\nNLBMS.497,NLBMS,3940.453,3230.801\nNLBMS.5089,NLBMS,3993.529,3215.04\nNLBMS.517,NLBMS,3972.601,3339.013\nNLBMS.519,NLBMS,4071.623,3211.919\nNLBMS.526,NLBMS,3980.074,3235.925\nNLBMS.533,NLBMS,4017.066,3253.066\nNLBMS.534,NLBMS,3970.493,3257.525\nNLBMS.537,NLBMS,4035.905,3091.688\nNLBMS.538,NLBMS,4033.381,3087.605\nNLBMS.545,NLBMS,4021.25,3259.888\nNLBMS.556,NLBMS,4089.602,3275.539\nNLBMS.563,NLBMS,4092.738,3281.178\nNLBMS.565,NLBMS,4049.961,3259.239\nNLBMS.570,NLBMS,3970.588,3275.012\nNLBMS.577,NLBMS,3940.474,3201.308\nNLBMS.583,NLBMS,3982.506,3236.506\nNLBMS.592,NLBMS,3932.34,3228.073\nNLBMS.595,NLBMS,4024.278,3257.405\nNLBMS.597,NLBMS,4023.667,3258.098\nNLBMS.605,NLBMS,4036.677,3198.047\nNLBMS.620,NLBMS,4022.926,3216.976\nNLBMS.638,NLBMS,3935.227,3240.756\nNLBMS.640,NLBMS,3933.571,3241.033\nNLBMS.643,NLBMS,3956.379,3240.557\nNLBMS.646,NLBMS,3958.017,3239.345\nNLBMS.664,NLBMS,3975.986,3172.316\nNLBMS.667,NLBMS,4032.488,3200.338\nNLBMS.681,NLBMS,4091.848,3213.805\nNLBMS.684,NLBMS,4080.341,3205.367\nNLBMS.686,NLBMS,4108.95,3244.596\nNLBMS.694,NLBMS,4025.457,3131.839\nNLBMS.699,NLBMS,3996.12,3227.463\nNLBMS.705,NLBMS,3974.181,3297.569\nNLBMS.710,NLBMS,3955.373,3275.091\nNLBMS.718,NLBMS,3955.779,3272.706\nNLBMS.738,NLBMS,3954.397,3269.864\nNLBMS.742,NLBMS,4047.959,3251.499\nNLBMS.748,NLBMS,3952.906,3267.104\nNLBMS.753,NLBMS,3953.315,3267.145\nNLBMS.760,NLBMS,3953.362,3265.885\nNLBMS.765,NLBMS,3952.288,3265.133\nNLBMS.793,NLBMS,4022.378,3318.267\nNLBMS.797,NLBMS,4022.111,3316.72\nNLBMS.826,NLBMS,3957.691,3291.48\nNLBMS.838,NLBMS,3950.725,3264.045\nNLBMS.841,NLBMS,3949.128,3263.742\nNLBMS.844,NLBMS,3949.695,3262.41\nNLBMS.846,NLBMS,3950.157,3262.019\nNLBMS.849,NLBMS,3949.456,3261.502\nNLBMS.857,NLBMS,3951.281,3262.164\nNLBMS.858,NLBMS,3948.759,3259.144\nNLBMS.862,NLBMS,3947.799,3256.856\nNLBMS.865,NLBMS,3951.499,3261.274\nNLBMS.869,NLBMS,3948.878,3261.58\nNLBMS.895,NLBMS,4024.539,3224.262\nNLBMS.896,NLBMS,4019.949,3316.144\nNLBMS.902,NLBMS,4052.906,3232.661\nNLBMS.910,NLBMS,3990.609,3229.255\nNLBMS.914,NLBMS,3956.544,3284.556\nNLBMS.917,NLBMS,3993.109,3171.211\nNLBMS.929,NLBMS,3934.349,3240.798\nNLBMS.960,NLBMS,3913.367,3209.995\nNLBMS.964,NLBMS,4020.611,3312.115\nNLBMS.979,NLBMS,4021.883,3315.478\nNLBMS.98,NLBMS,4036.065,3224.59\nNLBMS.997,NLBMS,3954.863,3195.891\nSEBMS.128,SEBMS,4428.869,3901.041\nSEBMS.16,SEBMS,4524.024,3598.771\nSEBMS.161,SEBMS,4682.288,4215.817\nSEBMS.162,SEBMS,4684.652,4216.17\nSEBMS.192,SEBMS,4625.196,4086.862\nSEBMS.198,SEBMS,4592.725,3894.413\nSEBMS.20,SEBMS,4519.478,3693.324\nSEBMS.235,SEBMS,4764.11,4096.353\nSEBMS.24,SEBMS,4547.717,3670.292\nSEBMS.260,SEBMS,4776.071,4051.626\nUKBMS.1,UKBMS,3628.631,3305.507\nUKBMS.100,UKBMS,3381.618,3361.322\nUKBMS.1001,UKBMS,3522.539,3173.874\nUKBMS.1004,UKBMS,3524.888,3140.163\nUKBMS.1005,UKBMS,3542.34,3158.294\nUKBMS.1006,UKBMS,3564.461,3175.242\nUKBMS.1009,UKBMS,3536.046,3154.717\nUKBMS.101,UKBMS,3589.596,3662.344\nUKBMS.1010,UKBMS,3553.251,3173.04\nUKBMS.1013,UKBMS,3526.998,3154.238\nUKBMS.1014,UKBMS,3517.363,3154.276\nUKBMS.1016,UKBMS,3533.743,3166.324\nUKBMS.1017,UKBMS,3534.728,3166.161\nUKBMS.1018,UKBMS,3485.901,3165.227\nUKBMS.1019,UKBMS,3486.773,3163.761\nUKBMS.102,UKBMS,3475.789,3201.57\nUKBMS.1021,UKBMS,3546.904,3157.74\nUKBMS.1022,UKBMS,3557.448,3165.233\nUKBMS.1024,UKBMS,3549.775,3195.755\nUKBMS.1029,UKBMS,3549.642,3195.357\nUKBMS.1031,UKBMS,3530.916,3147.289\nUKBMS.1032,UKBMS,3516.941,3153.999\nUKBMS.1034,UKBMS,3531.389,3115.012\nUKBMS.1035,UKBMS,3520.113,3129.215\nUKBMS.1036,UKBMS,3523.625,3178.133\nUKBMS.1037,UKBMS,3529.098,3144.848\nUKBMS.1039,UKBMS,3532.28,3167.888\nUKBMS.104,UKBMS,3618.702,3330.316\nUKBMS.1040,UKBMS,3528.28,3164.488\nUKBMS.1044,UKBMS,3534.561,3146.685\nUKBMS.1048,UKBMS,3530.631,3141.242\nUKBMS.1049,UKBMS,3549.751,3148.657\nUKBMS.1050,UKBMS,3553.007,3186.184\nUKBMS.1051,UKBMS,3527.192,3167.717\nUKBMS.1052,UKBMS,3530.353,3157.27\nUKBMS.1054,UKBMS,3532.429,3127.535\nUKBMS.1055,UKBMS,3536.555,3142.189\nUKBMS.1057,UKBMS,3522.289,3176.658\nUKBMS.1058,UKBMS,3522.666,3174.706\nUKBMS.106,UKBMS,3697.156,3170.557\nUKBMS.1060,UKBMS,3510.929,3144.51\nUKBMS.1063,UKBMS,3528.402,3130.743\nUKBMS.1065,UKBMS,3559.66,3168.218\nUKBMS.1067,UKBMS,3527.323,3128.484\nUKBMS.107,UKBMS,3486.452,3258.115\nUKBMS.1070,UKBMS,3528.163,3154.552\nUKBMS.1073,UKBMS,3547.869,3139.602\nUKBMS.1074,UKBMS,3553.609,3178.381\nUKBMS.1075,UKBMS,3510.529,3141.786\nUKBMS.1076,UKBMS,3510.911,3140.478\nUKBMS.1081,UKBMS,3559.72,3186.633\nUKBMS.1083,UKBMS,3510.411,3140.126\nUKBMS.1084,UKBMS,3532.005,3164.987\nUKBMS.1089,UKBMS,3527.918,3143.825\nUKBMS.1090,UKBMS,3527.163,3154.007\nUKBMS.1093,UKBMS,3513.224,3173.514\nUKBMS.1094,UKBMS,3550.843,3147.844\nUKBMS.1095,UKBMS,3550.336,3146.607\nUKBMS.1098,UKBMS,3552.249,3151.47\nUKBMS.1099,UKBMS,3550.054,3148.584\nUKBMS.1103,UKBMS,3493.745,3539.833\nUKBMS.1106,UKBMS,3498.348,3526.897\nUKBMS.1108,UKBMS,3487.317,3533.552\nUKBMS.111,UKBMS,3706.981,3227.308\nUKBMS.1110,UKBMS,3490.258,3525.84\nUKBMS.1111,UKBMS,3491.077,3522.362\nUKBMS.1113,UKBMS,3490.417,3522.604\nUKBMS.1114,UKBMS,3487.697,3538.164\nUKBMS.1115,UKBMS,3489.276,3538.679\nUKBMS.1116,UKBMS,3487.918,3537.361\nUKBMS.1117,UKBMS,3489.069,3537.581\nUKBMS.1118,UKBMS,3489.248,3535.972\nUKBMS.1119,UKBMS,3488.272,3537.596\nUKBMS.112,UKBMS,3316.417,3144.877\nUKBMS.1123,UKBMS,3486.416,3541.642\nUKBMS.1124,UKBMS,3492.798,3537.799\nUKBMS.1125,UKBMS,3490.171,3529.921\nUKBMS.1127,UKBMS,3493.031,3537.859\nUKBMS.113,UKBMS,3716.621,3239.724\nUKBMS.1134,UKBMS,3523.226,3553.545\nUKBMS.1136,UKBMS,3474.178,3466.218\nUKBMS.1142,UKBMS,3496.686,3528.218\nUKBMS.1143,UKBMS,3491.214,3526.101\nUKBMS.1144,UKBMS,3496.38,3527.364\nUKBMS.1145,UKBMS,3490.496,3527.834\nUKBMS.1146,UKBMS,3496.88,3526.363\nUKBMS.1147,UKBMS,3490.534,3527.265\nUKBMS.1149,UKBMS,3490.196,3529.133\nUKBMS.115,UKBMS,3558.012,3245.834\nUKBMS.1150,UKBMS,3491.288,3529.217\nUKBMS.1151,UKBMS,3497.055,3526.814\nUKBMS.1157,UKBMS,3487.583,3537.553\nUKBMS.1159,UKBMS,3487.474,3532.318\nUKBMS.116,UKBMS,3353.501,3760.915\nUKBMS.1163,UKBMS,3490.076,3525.862\nUKBMS.1168,UKBMS,3460.337,3464.571\nUKBMS.117,UKBMS,3557.287,3249.017\nUKBMS.1176,UKBMS,3486.926,3534.227\nUKBMS.118,UKBMS,3682.412,3146.729\nUKBMS.12,UKBMS,3673.29,3200.039\nUKBMS.1208,UKBMS,3597.857,3248.212\nUKBMS.1210,UKBMS,3609.667,3275.686\nUKBMS.1213,UKBMS,3596.279,3250.697\nUKBMS.122,UKBMS,3662.79,3145.999\nUKBMS.1220,UKBMS,3605.515,3256.546\nUKBMS.1221,UKBMS,3598.233,3246.005\nUKBMS.1224,UKBMS,3612.125,3269.894\nUKBMS.1231,UKBMS,3615.306,3323.346\nUKBMS.124,UKBMS,3476.767,3534.65\nUKBMS.126,UKBMS,3710.955,3153.192\nUKBMS.128,UKBMS,3519.115,3552.318\nUKBMS.130,UKBMS,3514.204,3149.047\nUKBMS.1301,UKBMS,3568.039,3253.909\nUKBMS.1305,UKBMS,3589.926,3242.919\nUKBMS.1306,UKBMS,3584.853,3240.716\nUKBMS.1307,UKBMS,3569.924,3253.995\nUKBMS.1309,UKBMS,3580.323,3237.204\nUKBMS.131,UKBMS,3500.241,3156.444\nUKBMS.1310,UKBMS,3549.681,3227.849\nUKBMS.1311,UKBMS,3540.047,3201.938\nUKBMS.1312,UKBMS,3546.719,3219.648\nUKBMS.1314,UKBMS,3554.504,3224.638\nUKBMS.1315,UKBMS,3559.167,3244.398\nUKBMS.1317,UKBMS,3586.828,3214.251\nUKBMS.1318,UKBMS,3586.151,3238.469\nUKBMS.1319,UKBMS,3553.988,3208.43\nUKBMS.132,UKBMS,3718.829,3150.127\nUKBMS.1320,UKBMS,3553.617,3222.183\nUKBMS.1321,UKBMS,3550.262,3216.032\nUKBMS.1322,UKBMS,3560.57,3225.128\nUKBMS.1323,UKBMS,3585.758,3239.195\nUKBMS.1325,UKBMS,3579.13,3210.106\nUKBMS.1328,UKBMS,3577.006,3234.597\nUKBMS.1329,UKBMS,3589.942,3243.018\nUKBMS.133,UKBMS,3457.227,3136.779\nUKBMS.1335,UKBMS,3554.533,3204.012\nUKBMS.1336,UKBMS,3552.938,3262.641\nUKBMS.1339,UKBMS,3552.953,3214.215\nUKBMS.134,UKBMS,3356.011,3145.637\nUKBMS.1340,UKBMS,3525.122,3222.314\nUKBMS.1343,UKBMS,3550.171,3198.439\nUKBMS.1344,UKBMS,3572.677,3213.557\nUKBMS.1345,UKBMS,3568.865,3207.683\nUKBMS.1346,UKBMS,3572.55,3194.212\nUKBMS.1349,UKBMS,3561.704,3247.923\nUKBMS.135,UKBMS,3705.356,3173.507\nUKBMS.1351,UKBMS,3580.446,3229.992\nUKBMS.1354,UKBMS,3584.318,3220.085\nUKBMS.1355,UKBMS,3567.002,3258.404\nUKBMS.1356,UKBMS,3571.21,3232.612\nUKBMS.1357,UKBMS,3573.451,3218.254\nUKBMS.1358,UKBMS,3570.146,3231.686\nUKBMS.136,UKBMS,3496.209,3230.967\nUKBMS.1361,UKBMS,3564.281,3221.031\nUKBMS.1362,UKBMS,3524.452,3216.532\nUKBMS.1368,UKBMS,3553.188,3213.937\nUKBMS.1369,UKBMS,3554.169,3208.085\nUKBMS.137,UKBMS,3602.18,3445.473\nUKBMS.1370,UKBMS,3563.861,3249.1\nUKBMS.1371,UKBMS,3551.891,3244.904\nUKBMS.1372,UKBMS,3538.563,3201.689\nUKBMS.1373,UKBMS,3576.623,3224.351\nUKBMS.1375,UKBMS,3557.219,3225.995\nUKBMS.1377,UKBMS,3585.527,3239.385\nUKBMS.1378,UKBMS,3584.869,3275.741\nUKBMS.1379,UKBMS,3536.31,3270.544\nUKBMS.1380,UKBMS,3532.731,3264.84\nUKBMS.1382,UKBMS,3547.331,3199.165\nUKBMS.1383,UKBMS,3545.675,3240.077\nUKBMS.139,UKBMS,3697.509,3171.318\nUKBMS.1391,UKBMS,3553.979,3223.913\nUKBMS.1396,UKBMS,3581.47,3237.578\nUKBMS.1397,UKBMS,3577.873,3235.784\nUKBMS.1399,UKBMS,3578.045,3236.049\nUKBMS.14,UKBMS,3582.553,3221.89\nUKBMS.140,UKBMS,3564.289,3899.21\nUKBMS.1401,UKBMS,3592.858,3134.988\nUKBMS.1402,UKBMS,3631.325,3122.92\nUKBMS.1404,UKBMS,3568.351,3142.604\nUKBMS.141,UKBMS,3427.671,3985.422\nUKBMS.1410,UKBMS,3560.724,3150.386\nUKBMS.1411,UKBMS,3612.26,3130.607\nUKBMS.1412,UKBMS,3611.826,3130.425\nUKBMS.1413,UKBMS,3612.654,3130.262\nUKBMS.1414,UKBMS,3621.753,3133.129\nUKBMS.1415,UKBMS,3626.173,3151.483\nUKBMS.1419,UKBMS,3615.716,3153.807\nUKBMS.142,UKBMS,3608.157,3240.011\nUKBMS.1420,UKBMS,3615.783,3153.762\nUKBMS.1422,UKBMS,3631.477,3128.084\nUKBMS.1423,UKBMS,3614.086,3126.193\nUKBMS.1427,UKBMS,3612.524,3144.654\nUKBMS.143,UKBMS,3496.39,3952.204\nUKBMS.1430,UKBMS,3623.551,3154.377\nUKBMS.1431,UKBMS,3629.12,3153.676\nUKBMS.1432,UKBMS,3606.693,3161.873\nUKBMS.1433,UKBMS,3582.378,3160.817\nUKBMS.1434,UKBMS,3581.287,3159.903\nUKBMS.1438,UKBMS,3578.287,3147.228\nUKBMS.1439,UKBMS,3601.217,3137.879\nUKBMS.144,UKBMS,3518.443,3676.41\nUKBMS.1445,UKBMS,3628.146,3149.729\nUKBMS.1446,UKBMS,3629.625,3150.086\nUKBMS.1447,UKBMS,3628.961,3150.559\nUKBMS.145,UKBMS,3367.167,3432.773\nUKBMS.1451,UKBMS,3616.947,3135.68\nUKBMS.1453,UKBMS,3599.195,3135.36\nUKBMS.1454,UKBMS,3653.99,3130.041\nUKBMS.1459,UKBMS,3638.73,3156.593\nUKBMS.146,UKBMS,3565.682,3806.173\nUKBMS.1460,UKBMS,3637.915,3156.748\nUKBMS.1469,UKBMS,3569.493,3159.364\nUKBMS.147,UKBMS,3467.433,3387.332\nUKBMS.1471,UKBMS,3564.612,3160.489\nUKBMS.1474,UKBMS,3598.21,3160.414\nUKBMS.1476,UKBMS,3585.809,3129.658\nUKBMS.148,UKBMS,3365.196,3320.966\nUKBMS.1480,UKBMS,3610.829,3131.382\nUKBMS.1482,UKBMS,3599.718,3133.154\nUKBMS.1484,UKBMS,3619.906,3140.625\nUKBMS.1487,UKBMS,3578.109,3138.854\nUKBMS.149,UKBMS,3373.856,3245.508\nUKBMS.15,UKBMS,3523.757,3260.757\nUKBMS.150,UKBMS,3318.2,3575.551\nUKBMS.1503,UKBMS,3680.587,3158.818\nUKBMS.151,UKBMS,3347.493,3320.056\nUKBMS.1511,UKBMS,3680.982,3147.693\nUKBMS.1519,UKBMS,3720.993,3168.15\nUKBMS.1520,UKBMS,3676.566,3170.191\nUKBMS.1521,UKBMS,3677.569,3172.232\nUKBMS.1522,UKBMS,3677.271,3170.328\nUKBMS.1523,UKBMS,3675.267,3169.607\nUKBMS.1524,UKBMS,3680.584,3148.258\nUKBMS.1525,UKBMS,3646.662,3163.023\nUKBMS.1526,UKBMS,3645.431,3162.313\nUKBMS.1529,UKBMS,3697.635,3165.835\nUKBMS.153,UKBMS,3373.632,3710.394\nUKBMS.1536,UKBMS,3631.909,3189.699\nUKBMS.1537,UKBMS,3634.631,3191.336\nUKBMS.1538,UKBMS,3721.163,3156.085\nUKBMS.1539,UKBMS,3630.827,3182.4\nUKBMS.1541,UKBMS,3630.17,3182.712\nUKBMS.1542,UKBMS,3646.662,3163.023\nUKBMS.1544,UKBMS,3720.17,3169.909\nUKBMS.1545,UKBMS,3719.761,3169.267\nUKBMS.1546,UKBMS,3701.65,3171.764\nUKBMS.155,UKBMS,3557.762,3245.209\nUKBMS.1552,UKBMS,3702.324,3173.467\nUKBMS.1556,UKBMS,3627.588,3187.303\nUKBMS.156,UKBMS,3695.127,3161.985\nUKBMS.1561,UKBMS,3721.18,3163.757\nUKBMS.1564,UKBMS,3666.357,3181.171\nUKBMS.1565,UKBMS,3665.47,3181.927\nUKBMS.157,UKBMS,3467.499,3388.524\nUKBMS.158,UKBMS,3470.132,3164.493\nUKBMS.160,UKBMS,3451.665,3628.441\nUKBMS.1603,UKBMS,3637.418,3397.94\nUKBMS.1605,UKBMS,3617.15,3346.653\nUKBMS.1606,UKBMS,3618.07,3346.703\nUKBMS.1607,UKBMS,3633.729,3395.884\nUKBMS.1608,UKBMS,3609.511,3350.006\nUKBMS.161,UKBMS,3583.598,3567.891\nUKBMS.1610,UKBMS,3608.616,3351.01\nUKBMS.1613,UKBMS,3614.275,3348.827\nUKBMS.1614,UKBMS,3650.264,3402.223\nUKBMS.1615,UKBMS,3634.607,3392.08\nUKBMS.163,UKBMS,3500.899,3721.814\nUKBMS.165,UKBMS,3602.43,3862.882\nUKBMS.167,UKBMS,3444.195,3850.892\nUKBMS.17,UKBMS,3575.433,3227.391\nUKBMS.170,UKBMS,3661.152,3223.906\nUKBMS.1701,UKBMS,3352.483,3156.508\nUKBMS.1703,UKBMS,3346.927,3137.896\nUKBMS.1705,UKBMS,3394.581,3175.456\nUKBMS.1708,UKBMS,3380.614,3161.587\nUKBMS.1709,UKBMS,3347.092,3166.149\nUKBMS.171,UKBMS,3662.352,3223.869\nUKBMS.1710,UKBMS,3327.931,3154.308\nUKBMS.1717,UKBMS,3325.05,3141.062\nUKBMS.172,UKBMS,3636.795,3285.088\nUKBMS.1720,UKBMS,3362.893,3132.025\nUKBMS.1721,UKBMS,3359.505,3141.833\nUKBMS.1728,UKBMS,3394.323,3175.47\nUKBMS.173,UKBMS,3636.649,3284.178\nUKBMS.1733,UKBMS,3409.429,3152.323\nUKBMS.1735,UKBMS,3356.461,3144.714\nUKBMS.1736,UKBMS,3362.323,3123.855\nUKBMS.1737,UKBMS,3362.343,3124.675\nUKBMS.1743,UKBMS,3326.472,3125.27\nUKBMS.1744,UKBMS,3365.393,3177.733\nUKBMS.1745,UKBMS,3351.49,3136.284\nUKBMS.1748,UKBMS,3351.782,3138.681\nUKBMS.175,UKBMS,3451.186,3628.13\nUKBMS.176,UKBMS,3747.588,3262.35\nUKBMS.177,UKBMS,3466.712,3736.266\nUKBMS.178,UKBMS,3702.249,3147.478\nUKBMS.1780,UKBMS,3351.791,3138.691\nUKBMS.179,UKBMS,3635.939,3285.218\nUKBMS.1811,UKBMS,3275.006,3159.142\nUKBMS.1813,UKBMS,3225.088,3124.432\nUKBMS.182,UKBMS,3390.249,3816.285\nUKBMS.1821,UKBMS,3225.073,3125.573\nUKBMS.1832,UKBMS,3207.319,3112.063\nUKBMS.1835,UKBMS,3261.961,3122.066\nUKBMS.1836,UKBMS,3231.717,3117.432\nUKBMS.184,UKBMS,3680.765,3286.876\nUKBMS.1840,UKBMS,3240.61,3102.036\nUKBMS.1841,UKBMS,3261.509,3118.685\nUKBMS.1842,UKBMS,3264.213,3117.929\nUKBMS.188,UKBMS,3502.729,3140.791\nUKBMS.190,UKBMS,3534.894,3774.343\nUKBMS.1901,UKBMS,3612.369,3240.299\nUKBMS.1903,UKBMS,3618.891,3227.739\nUKBMS.1904,UKBMS,3591.453,3241.751\nUKBMS.1905,UKBMS,3611.607,3212.55\nUKBMS.1906,UKBMS,3619.511,3241.243\nUKBMS.1908,UKBMS,3586.658,3244.58\nUKBMS.1911,UKBMS,3641.922,3246.345\nUKBMS.192,UKBMS,3531.055,3773.955\nUKBMS.1920,UKBMS,3635.471,3262.346\nUKBMS.1922,UKBMS,3611.419,3212.879\nUKBMS.1925,UKBMS,3634.929,3262.132\nUKBMS.1926,UKBMS,3634.978,3261.819\nUKBMS.1927,UKBMS,3633.369,3261.884\nUKBMS.1928,UKBMS,3633.352,3261.785\nUKBMS.193,UKBMS,3495.226,3133.299\nUKBMS.1932,UKBMS,3600.462,3211.205\nUKBMS.1934,UKBMS,3619.248,3245.146\nUKBMS.1935,UKBMS,3622.137,3247.305\nUKBMS.1937,UKBMS,3615.032,3230.107\nUKBMS.1940,UKBMS,3630.595,3236.656\nUKBMS.1941,UKBMS,3629.692,3233.049\nUKBMS.1942,UKBMS,3607.119,3226.956\nUKBMS.1943,UKBMS,3600.401,3204.716\nUKBMS.1945,UKBMS,3623.348,3211.258\nUKBMS.195,UKBMS,3534.912,3773.625\nUKBMS.1954,UKBMS,3602.799,3232.651\nUKBMS.1965,UKBMS,3607.44,3210.165\nUKBMS.1966,UKBMS,3607.073,3210.169\nUKBMS.1968,UKBMS,3598.793,3234.74\nUKBMS.1970,UKBMS,3627.129,3204.741\nUKBMS.1971,UKBMS,3627.82,3238.235\nUKBMS.1977,UKBMS,3587.316,3241.221\nUKBMS.1981,UKBMS,3609.55,3198.828\nUKBMS.1982,UKBMS,3618.018,3236.516\nUKBMS.1984,UKBMS,3624.518,3250.056\nUKBMS.1989,UKBMS,3643.827,3239.53\nUKBMS.1991,UKBMS,3607.048,3209.47\nUKBMS.1992,UKBMS,3613.436,3256.472\nUKBMS.1993,UKBMS,3621.514,3242.331\nUKBMS.1998,UKBMS,3619.561,3238.493\nUKBMS.1999,UKBMS,3628.134,3229.145\nUKBMS.2,UKBMS,3627.181,3303.335\nUKBMS.2000,UKBMS,3605.258,3185.224\nUKBMS.2001,UKBMS,3613.167,3185.97\nUKBMS.2003,UKBMS,3604.286,3176.754\nUKBMS.2004,UKBMS,3603.794,3176.836\nUKBMS.2006,UKBMS,3603.661,3177.873\nUKBMS.2007,UKBMS,3581.632,3190.755\nUKBMS.2008,UKBMS,3619.077,3180.553\nUKBMS.2010,UKBMS,3600.354,3176.544\nUKBMS.2013,UKBMS,3616.832,3181.298\nUKBMS.2014,UKBMS,3594.813,3175.775\nUKBMS.2015,UKBMS,3605.73,3178.88\nUKBMS.2016,UKBMS,3605.268,3179.637\nUKBMS.2017,UKBMS,3625.542,3183.903\nUKBMS.2018,UKBMS,3608.933,3182.154\nUKBMS.2019,UKBMS,3620.087,3181.766\nUKBMS.2020,UKBMS,3617.043,3191.725\nUKBMS.2021,UKBMS,3617.589,3191.354\nUKBMS.2024,UKBMS,3582.273,3162.942\nUKBMS.2025,UKBMS,3623.622,3176.768\nUKBMS.2026,UKBMS,3624.649,3176.768\nUKBMS.2027,UKBMS,3613.563,3182.829\nUKBMS.2028,UKBMS,3620.113,3183.231\nUKBMS.2031,UKBMS,3594.801,3179.141\nUKBMS.2032,UKBMS,3623.578,3191.218\nUKBMS.2034,UKBMS,3595.832,3175.822\nUKBMS.2035,UKBMS,3612.279,3176.832\nUKBMS.2037,UKBMS,3624.21,3179.904\nUKBMS.2038,UKBMS,3594.075,3176.215\nUKBMS.2039,UKBMS,3625.408,3183.909\nUKBMS.2040,UKBMS,3616.512,3197.468\nUKBMS.2041,UKBMS,3612.766,3201.747\nUKBMS.2042,UKBMS,3620.73,3182.553\nUKBMS.2043,UKBMS,3585.536,3192.308\nUKBMS.2044,UKBMS,3617.373,3180.564\nUKBMS.2045,UKBMS,3611.384,3186.643\nUKBMS.2046,UKBMS,3623.081,3166.731\nUKBMS.2047,UKBMS,3574.487,3163.15\nUKBMS.2049,UKBMS,3605.145,3185.027\nUKBMS.2051,UKBMS,3570.103,3179.383\nUKBMS.2052,UKBMS,3586.159,3194.489\nUKBMS.2053,UKBMS,3599.136,3176.401\nUKBMS.2054,UKBMS,3617.197,3181.346\nUKBMS.2055,UKBMS,3608.913,3199.34\nUKBMS.2056,UKBMS,3586.567,3177.609\nUKBMS.2057,UKBMS,3612.181,3197.477\nUKBMS.2058,UKBMS,3577.165,3170.492\nUKBMS.2061,UKBMS,3600.005,3188.229\nUKBMS.2062,UKBMS,3618.444,3186.486\nUKBMS.2065,UKBMS,3567.058,3172.663\nUKBMS.2066,UKBMS,3617.811,3177.234\nUKBMS.2067,UKBMS,3606.18,3186.391\nUKBMS.2068,UKBMS,3605.891,3180.346\nUKBMS.2071,UKBMS,3625.634,3187.12\nUKBMS.2072,UKBMS,3623.438,3195.496\nUKBMS.2073,UKBMS,3606.277,3187.593\nUKBMS.2075,UKBMS,3602.729,3171.809\nUKBMS.2077,UKBMS,3605.567,3176.542\nUKBMS.2079,UKBMS,3575.088,3170.938\nUKBMS.2080,UKBMS,3596.112,3175.674\nUKBMS.2081,UKBMS,3623.064,3166.632\nUKBMS.2082,UKBMS,3608.804,3188.798\nUKBMS.2083,UKBMS,3610.857,3187.746\nUKBMS.2085,UKBMS,3601.705,3180.127\nUKBMS.2086,UKBMS,3601.393,3180.078\nUKBMS.2087,UKBMS,3610.399,3186.198\nUKBMS.2089,UKBMS,3618.317,3182.04\nUKBMS.21,UKBMS,3622.506,3326.958\nUKBMS.2104,UKBMS,3545.351,3309.782\nUKBMS.2105,UKBMS,3545.14,3310.427\nUKBMS.2106,UKBMS,3544.607,3310.025\nUKBMS.2122,UKBMS,3531.314,3296.374\nUKBMS.2127,UKBMS,3525.178,3300.954\nUKBMS.2131,UKBMS,3544.627,3310.242\nUKBMS.2136,UKBMS,3543.318,3297.44\nUKBMS.2200,UKBMS,3492.27,3147.06\nUKBMS.2201,UKBMS,3491.723,3145.506\nUKBMS.2202,UKBMS,3476.15,3149.209\nUKBMS.2204,UKBMS,3432.898,3160.639\nUKBMS.2206,UKBMS,3446.943,3153.093\nUKBMS.2207,UKBMS,3470.764,3166.15\nUKBMS.2209,UKBMS,3471.935,3128.852\nUKBMS.2210,UKBMS,3456.662,3160.728\nUKBMS.2211,UKBMS,3465.507,3170.559\nUKBMS.2214,UKBMS,3460.354,3170.691\nUKBMS.2215,UKBMS,3482.613,3164.509\nUKBMS.2219,UKBMS,3463.413,3136.111\nUKBMS.2221,UKBMS,3466.749,3158.535\nUKBMS.2225,UKBMS,3469.55,3163.929\nUKBMS.2226,UKBMS,3435.207,3151.834\nUKBMS.2227,UKBMS,3434.659,3152.029\nUKBMS.2228,UKBMS,3434.577,3151.125\nUKBMS.2229,UKBMS,3485.291,3140.716\nUKBMS.2230,UKBMS,3444.361,3133.714\nUKBMS.2231,UKBMS,3462.788,3131.739\nUKBMS.2233,UKBMS,3473.546,3166.591\nUKBMS.2234,UKBMS,3468.918,3155.058\nUKBMS.2235,UKBMS,3490.503,3150.136\nUKBMS.2238,UKBMS,3478.735,3146.32\nUKBMS.2239,UKBMS,3461.327,3161.716\nUKBMS.2244,UKBMS,3433.541,3150.336\nUKBMS.2245,UKBMS,3433.816,3150.637\nUKBMS.2246,UKBMS,3433.056,3150.043\nUKBMS.2247,UKBMS,3433.174,3150.41\nUKBMS.2248,UKBMS,3491.818,3140.895\nUKBMS.2249,UKBMS,3480.23,3156.455\nUKBMS.2250,UKBMS,3487.383,3140.549\nUKBMS.2251,UKBMS,3471.887,3162.82\nUKBMS.2255,UKBMS,3476.699,3140.781\nUKBMS.2256,UKBMS,3477.199,3139.977\nUKBMS.2257,UKBMS,3488.848,3148.683\nUKBMS.2258,UKBMS,3485.448,3156.564\nUKBMS.2266,UKBMS,3480.401,3160.006\nUKBMS.2267,UKBMS,3441.282,3155.276\nUKBMS.2268,UKBMS,3469.548,3134.591\nUKBMS.2270,UKBMS,3442.104,3148.659\nUKBMS.2271,UKBMS,3475.128,3125.495\nUKBMS.2272,UKBMS,3496.068,3135.823\nUKBMS.2273,UKBMS,3455.283,3163.234\nUKBMS.2274,UKBMS,3417.458,3149.069\nUKBMS.2278,UKBMS,3454.965,3160.306\nUKBMS.2279,UKBMS,3471.937,3127.877\nUKBMS.2280,UKBMS,3469.299,3132.521\nUKBMS.2282,UKBMS,3487.1,3139.713\nUKBMS.2284,UKBMS,3494.562,3135.263\nUKBMS.2287,UKBMS,3458.891,3128.998\nUKBMS.2293,UKBMS,3475.573,3127.441\nUKBMS.2295,UKBMS,3415.555,3160.737\nUKBMS.2297,UKBMS,3487.085,3146.165\nUKBMS.2299,UKBMS,3446.707,3151.639\nUKBMS.23,UKBMS,3440.49,3201.441\nUKBMS.2305,UKBMS,3474.907,3303.548\nUKBMS.2307,UKBMS,3484.073,3288.5\nUKBMS.2309,UKBMS,3480.416,3286.998\nUKBMS.2310,UKBMS,3498.802,3318.734\nUKBMS.2311,UKBMS,3477.676,3283.634\nUKBMS.2315,UKBMS,3456.776,3269.368\nUKBMS.2316,UKBMS,3491.323,3261.973\nUKBMS.2318,UKBMS,3440.045,3283.886\nUKBMS.2319,UKBMS,3501.976,3300.238\nUKBMS.2322,UKBMS,3461.235,3286.913\nUKBMS.2325,UKBMS,3479.148,3299.886\nUKBMS.2328,UKBMS,3497.231,3263.337\nUKBMS.2331,UKBMS,3478.23,3286.27\nUKBMS.2332,UKBMS,3491.736,3253.776\nUKBMS.2333,UKBMS,3486.136,3307.672\nUKBMS.2338,UKBMS,3491.35,3300.805\nUKBMS.2339,UKBMS,3498.788,3268.871\nUKBMS.2340,UKBMS,3481.062,3249.872\nUKBMS.2341,UKBMS,3545.245,3310.206\nUKBMS.2345,UKBMS,3498.068,3303.953\nUKBMS.2347,UKBMS,3543.217,3299.184\nUKBMS.2348,UKBMS,3543.272,3299.031\nUKBMS.2354,UKBMS,3447.203,3323.08\nUKBMS.2355,UKBMS,3510.28,3290.099\nUKBMS.2358,UKBMS,3482.999,3324.153\nUKBMS.2366,UKBMS,3510.577,3287.536\nUKBMS.2370,UKBMS,3481.067,3289.06\nUKBMS.2372,UKBMS,3459.077,3375.927\nUKBMS.2380,UKBMS,3497.631,3367.938\nUKBMS.2388,UKBMS,3543.839,3295.707\nUKBMS.2389,UKBMS,3539.299,3317.171\nUKBMS.24,UKBMS,3373.856,3360.049\nUKBMS.2401,UKBMS,3508.072,3170.553\nUKBMS.2402,UKBMS,3509.586,3170.452\nUKBMS.2403,UKBMS,3508.292,3170.424\nUKBMS.2404,UKBMS,3508.554,3169.565\nUKBMS.2405,UKBMS,3476.373,3203.8\nUKBMS.2406,UKBMS,3507.131,3170.642\nUKBMS.2407,UKBMS,3503.39,3183.558\nUKBMS.2408,UKBMS,3504.112,3183.031\nUKBMS.2409,UKBMS,3478.141,3198.323\nUKBMS.2410,UKBMS,3505.592,3167.443\nUKBMS.2411,UKBMS,3477.035,3204.638\nUKBMS.2413,UKBMS,3471.999,3207.576\nUKBMS.2415,UKBMS,3505.592,3167.443\nUKBMS.2417,UKBMS,3480.78,3183.589\nUKBMS.2418,UKBMS,3509.156,3180.017\nUKBMS.2424,UKBMS,3474.428,3195.689\nUKBMS.2426,UKBMS,3497.699,3229.708\nUKBMS.2427,UKBMS,3477.148,3223.184\nUKBMS.2428,UKBMS,3465.169,3183.712\nUKBMS.2429,UKBMS,3510.978,3169.335\nUKBMS.2430,UKBMS,3509.255,3172.971\nUKBMS.2431,UKBMS,3506.981,3169.575\nUKBMS.2432,UKBMS,3492.628,3211.711\nUKBMS.2434,UKBMS,3474.271,3215.329\nUKBMS.2436,UKBMS,3506.38,3218.624\nUKBMS.2438,UKBMS,3503.212,3176.068\nUKBMS.2441,UKBMS,3477.666,3175.436\nUKBMS.2443,UKBMS,3501.635,3234.792\nUKBMS.2448,UKBMS,3495.961,3232.48\nUKBMS.2450,UKBMS,3486.049,3165.507\nUKBMS.25,UKBMS,3669.297,3286.17\nUKBMS.2500,UKBMS,3469.931,3216.864\nUKBMS.2502,UKBMS,3437.457,3204.525\nUKBMS.2505,UKBMS,3384.889,3190.886\nUKBMS.2506,UKBMS,3430.515,3210.289\nUKBMS.2510,UKBMS,3402.102,3176.654\nUKBMS.2513,UKBMS,3410.986,3178.279\nUKBMS.2514,UKBMS,3408.735,3177.505\nUKBMS.2525,UKBMS,3432.173,3194.169\nUKBMS.2527,UKBMS,3394.446,3183.489\nUKBMS.2528,UKBMS,3436.993,3204.253\nUKBMS.2529,UKBMS,3434.784,3184.353\nUKBMS.2530,UKBMS,3419.01,3180.367\nUKBMS.2537,UKBMS,3407.891,3174.42\nUKBMS.2540,UKBMS,3446.349,3182.471\nUKBMS.2541,UKBMS,3428.485,3210.558\nUKBMS.2542,UKBMS,3429.376,3195.407\nUKBMS.2543,UKBMS,3383.689,3190.727\nUKBMS.2546,UKBMS,3412.468,3186.785\nUKBMS.2549,UKBMS,3433.751,3186.023\nUKBMS.2556,UKBMS,3443.828,3216.421\nUKBMS.2557,UKBMS,3435.304,3187.47\nUKBMS.2562,UKBMS,3425.87,3222.038\nUKBMS.2566,UKBMS,3438.371,3204.024\nUKBMS.2570,UKBMS,3434.763,3227.922\nUKBMS.2571,UKBMS,3436.16,3228.189\nUKBMS.2572,UKBMS,3446.911,3183.114\nUKBMS.2573,UKBMS,3422.669,3206.942\nUKBMS.2574,UKBMS,3423.594,3221.992\nUKBMS.2575,UKBMS,3437.011,3205.08\nUKBMS.2578,UKBMS,3425.831,3192.388\nUKBMS.2579,UKBMS,3434.941,3187.519\nUKBMS.2580,UKBMS,3435.697,3186.58\nUKBMS.2582,UKBMS,3441.179,3200.812\nUKBMS.2583,UKBMS,3443.867,3197.892\nUKBMS.2586,UKBMS,3453.554,3214.243\nUKBMS.2587,UKBMS,3432.271,3189.073\nUKBMS.2588,UKBMS,3448.684,3226.336\nUKBMS.2589,UKBMS,3466.614,3242.978\nUKBMS.2592,UKBMS,3437.77,3228.85\nUKBMS.2593,UKBMS,3432.769,3186.381\nUKBMS.2594,UKBMS,3448.413,3226.595\nUKBMS.2596,UKBMS,3417.817,3166.318\nUKBMS.2597,UKBMS,3439.251,3202.226\nUKBMS.2599,UKBMS,3435.249,3187.017\nUKBMS.26,UKBMS,3500.956,3430.87\nUKBMS.2605,UKBMS,3459.618,3699.832\nUKBMS.2613,UKBMS,3539.297,3786.898\nUKBMS.2615,UKBMS,3476.685,3824.124\nUKBMS.2616,UKBMS,3508.406,3722.77\nUKBMS.2618,UKBMS,3435.712,3720.509\nUKBMS.2619,UKBMS,3359.989,3828.568\nUKBMS.2620,UKBMS,3506.807,3764.27\nUKBMS.2627,UKBMS,3432.842,3740.647\nUKBMS.2628,UKBMS,3445.303,3727.016\nUKBMS.2630,UKBMS,3448.773,3715.099\nUKBMS.2640,UKBMS,3405.313,3708.571\nUKBMS.2641,UKBMS,3404.01,3825.892\nUKBMS.2643,UKBMS,3431.485,3640.295\nUKBMS.2648,UKBMS,3446.942,3715.559\nUKBMS.2652,UKBMS,3503.376,3724.583\nUKBMS.2656,UKBMS,3504.152,3724.093\nUKBMS.2657,UKBMS,3502.823,3724.369\nUKBMS.2665,UKBMS,3444.013,3703.988\nUKBMS.2668,UKBMS,3452.153,3893.097\nUKBMS.2671,UKBMS,3529.608,3732.425\nUKBMS.2677,UKBMS,3504.167,3722.129\nUKBMS.2678,UKBMS,3505.556,3726.062\nUKBMS.2679,UKBMS,3461.725,3756.97\nUKBMS.2690,UKBMS,3547.709,3720.782\nUKBMS.2691,UKBMS,3535.353,3723.172\nUKBMS.2697,UKBMS,3469.759,3738.566\nUKBMS.2698,UKBMS,3469.262,3755.16\nUKBMS.2699,UKBMS,3426.652,3723.836\nUKBMS.27,UKBMS,3628.754,3313.074\nUKBMS.2702,UKBMS,3536.378,3397.457\nUKBMS.2705,UKBMS,3510.406,3407.355\nUKBMS.2708,UKBMS,3503.983,3413.435\nUKBMS.2709,UKBMS,3492.411,3436.327\nUKBMS.2713,UKBMS,3502.924,3411.763\nUKBMS.2718,UKBMS,3518.766,3415.82\nUKBMS.28,UKBMS,3751.8,3277.008\nUKBMS.2800,UKBMS,3556.53,3376.402\nUKBMS.2801,UKBMS,3552.027,3373.701\nUKBMS.2802,UKBMS,3562.608,3378.231\nUKBMS.2805,UKBMS,3545.215,3391.608\nUKBMS.2806,UKBMS,3557.859,3366.831\nUKBMS.2817,UKBMS,3538.86,3406.188\nUKBMS.2818,UKBMS,3543.143,3390.126\nUKBMS.2819,UKBMS,3561.842,3380.899\nUKBMS.2820,UKBMS,3554.074,3379.785\nUKBMS.2821,UKBMS,3553.484,3380.57\nUKBMS.2822,UKBMS,3572.931,3332.915\nUKBMS.2824,UKBMS,3555.373,3389.4\nUKBMS.2825,UKBMS,3537.957,3413.454\nUKBMS.2826,UKBMS,3536.842,3414.047\nUKBMS.2827,UKBMS,3563.382,3379.219\nUKBMS.2828,UKBMS,3562.315,3402.972\nUKBMS.2830,UKBMS,3571.927,3371.625\nUKBMS.2831,UKBMS,3590.238,3384.782\nUKBMS.2833,UKBMS,3560.098,3387.39\nUKBMS.2836,UKBMS,3537.671,3359.943\nUKBMS.2838,UKBMS,3588.831,3369.778\nUKBMS.2840,UKBMS,3576.726,3407.976\nUKBMS.2841,UKBMS,3574.431,3409.377\nUKBMS.2842,UKBMS,3576.588,3411.962\nUKBMS.2843,UKBMS,3582.362,3370.454\nUKBMS.2844,UKBMS,3528.347,3424.312\nUKBMS.2845,UKBMS,3559.028,3401.49\nUKBMS.2846,UKBMS,3565.532,3369.714\nUKBMS.2849,UKBMS,3547.978,3368.281\nUKBMS.2851,UKBMS,3539.095,3356.962\nUKBMS.2852,UKBMS,3548.145,3359.615\nUKBMS.2853,UKBMS,3544.051,3355.727\nUKBMS.2854,UKBMS,3568.746,3408.094\nUKBMS.2855,UKBMS,3537.142,3401.394\nUKBMS.2856,UKBMS,3549.689,3379.071\nUKBMS.2858,UKBMS,3547.581,3400.257\nUKBMS.2859,UKBMS,3554.526,3384.359\nUKBMS.2860,UKBMS,3553.988,3386.583\nUKBMS.2863,UKBMS,3571.763,3398.951\nUKBMS.2864,UKBMS,3553.216,3377.77\nUKBMS.2865,UKBMS,3569.16,3402.131\nUKBMS.2866,UKBMS,3563.936,3393.352\nUKBMS.2867,UKBMS,3537.068,3397.342\nUKBMS.2870,UKBMS,3569.085,3405.293\nUKBMS.2871,UKBMS,3566.152,3367.375\nUKBMS.2873,UKBMS,3568.085,3398.754\nUKBMS.2874,UKBMS,3549.377,3403.106\nUKBMS.2876,UKBMS,3544.03,3352.58\nUKBMS.2877,UKBMS,3545.299,3355.924\nUKBMS.2878,UKBMS,3553.282,3347.985\nUKBMS.2879,UKBMS,3547.019,3353.503\nUKBMS.2882,UKBMS,3567.265,3399.298\nUKBMS.29,UKBMS,3566.884,3229.637\nUKBMS.2901,UKBMS,3597.225,3558.542\nUKBMS.2905,UKBMS,3582.021,3568.046\nUKBMS.2909,UKBMS,3582.538,3548.964\nUKBMS.2911,UKBMS,3584.466,3570.627\nUKBMS.2915,UKBMS,3548.812,3560.8\nUKBMS.2924,UKBMS,3586.51,3559.263\nUKBMS.2928,UKBMS,3573.574,3591.597\nUKBMS.2930,UKBMS,3591.961,3569.972\nUKBMS.2931,UKBMS,3559.456,3584.15\nUKBMS.2932,UKBMS,3588.462,3571.941\nUKBMS.2934,UKBMS,3564.248,3578.569\nUKBMS.2935,UKBMS,3567.035,3579.879\nUKBMS.2936,UKBMS,3580.947,3580.973\nUKBMS.2937,UKBMS,3579.217,3608.456\nUKBMS.2938,UKBMS,3572.211,3600.694\nUKBMS.2939,UKBMS,3571.243,3587.95\nUKBMS.2940,UKBMS,3570.561,3692.94\nUKBMS.2941,UKBMS,3579.101,3585.915\nUKBMS.2944,UKBMS,3574.611,3553.136\nUKBMS.2945,UKBMS,3590.467,3660.846\nUKBMS.2947,UKBMS,3572.434,3596.083\nUKBMS.2952,UKBMS,3589.442,3561.478\nUKBMS.2954,UKBMS,3584.247,3576.975\nUKBMS.2961,UKBMS,3550.942,3581.206\nUKBMS.2966,UKBMS,3550.597,3610.534\nUKBMS.2969,UKBMS,3582.752,3606.03\nUKBMS.2970,UKBMS,3583.998,3626.855\nUKBMS.3,UKBMS,3704.155,3357.351\nUKBMS.30,UKBMS,3566.109,3228.927\nUKBMS.3001,UKBMS,3366.986,3271.933\nUKBMS.3002,UKBMS,3445.006,3432.71\nUKBMS.3004,UKBMS,3403.95,3394.607\nUKBMS.3005,UKBMS,3302.159,3288.064\nUKBMS.3008,UKBMS,3396.166,3437.802\nUKBMS.3009,UKBMS,3292.063,3279.04\nUKBMS.3012,UKBMS,3382.207,3246.874\nUKBMS.3021,UKBMS,3439.721,3432.946\nUKBMS.3030,UKBMS,3326.622,3317.155\nUKBMS.3034,UKBMS,3427.278,3239.111\nUKBMS.3039,UKBMS,3413.452,3285.028\nUKBMS.3040,UKBMS,3445.092,3421.836\nUKBMS.3042,UKBMS,3416.436,3280.464\nUKBMS.3046,UKBMS,3296.092,3272.266\nUKBMS.3055,UKBMS,3458.027,3409.448\nUKBMS.3061,UKBMS,3415.676,3438.28\nUKBMS.3071,UKBMS,3353.207,3299.144\nUKBMS.3073,UKBMS,3416.679,3264.259\nUKBMS.31,UKBMS,3556.749,3244.081\nUKBMS.3102,UKBMS,3688.137,3278.565\nUKBMS.3104,UKBMS,3695.235,3269.718\nUKBMS.3106,UKBMS,3686.584,3294.385\nUKBMS.3109,UKBMS,3747.678,3267.389\nUKBMS.3110,UKBMS,3705.127,3254.291\nUKBMS.3111,UKBMS,3694.9,3267.566\nUKBMS.3114,UKBMS,3712.914,3245.018\nUKBMS.3115,UKBMS,3706.268,3253.391\nUKBMS.3116,UKBMS,3714.397,3263.3\nUKBMS.3119,UKBMS,3739.221,3261.596\nUKBMS.3120,UKBMS,3553.586,3222.15\nUKBMS.3124,UKBMS,3689.512,3274.837\nUKBMS.3125,UKBMS,3746.945,3278.773\nUKBMS.32,UKBMS,3691.103,3156.376\nUKBMS.3211,UKBMS,3643.98,3275.753\nUKBMS.3212,UKBMS,3640.329,3276.26\nUKBMS.3214,UKBMS,3644.242,3278.798\nUKBMS.3215,UKBMS,3665.047,3279.448\nUKBMS.3220,UKBMS,3691.834,3226.475\nUKBMS.3226,UKBMS,3649.83,3238.726\nUKBMS.3231,UKBMS,3623.471,3294.47\nUKBMS.3232,UKBMS,3634.911,3281.784\nUKBMS.3234,UKBMS,3638.303,3279.081\nUKBMS.3235,UKBMS,3672.734,3201.994\nUKBMS.3237,UKBMS,3656.113,3274.44\nUKBMS.3238,UKBMS,3690.804,3226.738\nUKBMS.3239,UKBMS,3653.829,3268.018\nUKBMS.3243,UKBMS,3654.519,3270.949\nUKBMS.3244,UKBMS,3650.263,3239.833\nUKBMS.3250,UKBMS,3648.889,3209.726\nUKBMS.3252,UKBMS,3688.057,3218.717\nUKBMS.3253,UKBMS,3632.914,3217.513\nUKBMS.3254,UKBMS,3632.23,3218.999\nUKBMS.3255,UKBMS,3632.203,3215.495\nUKBMS.3256,UKBMS,3627.467,3209.231\nUKBMS.3257,UKBMS,3660.866,3289.131\nUKBMS.3258,UKBMS,3668.795,3201.357\nUKBMS.3260,UKBMS,3675.032,3206.398\nUKBMS.3261,UKBMS,3631.967,3223.465\nUKBMS.3262,UKBMS,3631.523,3221.007\nUKBMS.3266,UKBMS,3646.286,3233.249\nUKBMS.3268,UKBMS,3628.317,3305.074\nUKBMS.3270,UKBMS,3642.879,3301.323\nUKBMS.3271,UKBMS,3673.709,3202.38\nUKBMS.3272,UKBMS,3626.698,3210.092\nUKBMS.3273,UKBMS,3627.308,3295.391\nUKBMS.3300,UKBMS,3703.571,3307.567\nUKBMS.3301,UKBMS,3734.631,3316.814\nUKBMS.3303,UKBMS,3716.563,3356.108\nUKBMS.3305,UKBMS,3760.828,3323.02\nUKBMS.3306,UKBMS,3743.591,3311.395\nUKBMS.3309,UKBMS,3681.043,3338.083\nUKBMS.3313,UKBMS,3715.662,3341.101\nUKBMS.3314,UKBMS,3730.911,3312.467\nUKBMS.3316,UKBMS,3697.071,3359.742\nUKBMS.3318,UKBMS,3745.645,3312.776\nUKBMS.3319,UKBMS,3715.791,3341.418\nUKBMS.3320,UKBMS,3717.523,3341.331\nUKBMS.3323,UKBMS,3749.339,3327.425\nUKBMS.3335,UKBMS,3734.942,3316.449\nUKBMS.3337,UKBMS,3682.884,3347.333\nUKBMS.3340,UKBMS,3692.58,3294.192\nUKBMS.3341,UKBMS,3685.229,3355.266\nUKBMS.3343,UKBMS,3704.059,3341.961\nUKBMS.34,UKBMS,3627.137,3273.678\nUKBMS.3400,UKBMS,3566.93,3454.431\nUKBMS.3403,UKBMS,3551.731,3485.945\nUKBMS.3404,UKBMS,3587.921,3442.052\nUKBMS.3405,UKBMS,3521.286,3514.01\nUKBMS.3408,UKBMS,3639.475,3510.515\nUKBMS.3411,UKBMS,3637.604,3455.976\nUKBMS.3412,UKBMS,3539.026,3505.627\nUKBMS.3413,UKBMS,3542.407,3531.259\nUKBMS.3414,UKBMS,3522.388,3510.746\nUKBMS.3415,UKBMS,3517.158,3523.149\nUKBMS.3416,UKBMS,3631.437,3463.818\nUKBMS.3420,UKBMS,3632.362,3514.557\nUKBMS.3421,UKBMS,3597.098,3522.83\nUKBMS.3424,UKBMS,3527.513,3542.305\nUKBMS.3425,UKBMS,3571.02,3499.214\nUKBMS.3428,UKBMS,3626.977,3512.585\nUKBMS.3429,UKBMS,3570.62,3502.556\nUKBMS.3430,UKBMS,3638.812,3517.331\nUKBMS.3434,UKBMS,3575.317,3498.933\nUKBMS.3435,UKBMS,3574.682,3498.693\nUKBMS.3500,UKBMS,3270.139,3606.533\nUKBMS.3501,UKBMS,3287.069,3646.148\nUKBMS.3502,UKBMS,3299.489,3608.456\nUKBMS.3504,UKBMS,3291.519,3611.612\nUKBMS.3511,UKBMS,3244.159,3658.29\nUKBMS.3523,UKBMS,3297.467,3606.769\nUKBMS.3539,UKBMS,3193.271,3588.33\nUKBMS.3541,UKBMS,3269.163,3603.158\nUKBMS.3546,UKBMS,3270.139,3606.533\nUKBMS.3553,UKBMS,3276.027,3676.733\nUKBMS.3554,UKBMS,3251.802,3691.985\nUKBMS.3555,UKBMS,3294.829,3675.416\nUKBMS.3556,UKBMS,3171.94,3603.834\nUKBMS.3558,UKBMS,3300.38,3617.21\nUKBMS.36,UKBMS,3744.954,3321.158\nUKBMS.3602,UKBMS,3496.653,3263.332\nUKBMS.3603,UKBMS,3480.136,3293.629\nUKBMS.3604,UKBMS,3471.365,3242.644\nUKBMS.3606,UKBMS,3486.456,3295.726\nUKBMS.3607,UKBMS,3456.464,3260.17\nUKBMS.3609,UKBMS,3498.919,3269.052\nUKBMS.3610,UKBMS,3461.278,3287.931\nUKBMS.3612,UKBMS,3434.589,3300.209\nUKBMS.3613,UKBMS,3442.573,3338.298\nUKBMS.3619,UKBMS,3474.667,3387.07\nUKBMS.3620,UKBMS,3498.046,3303.768\nUKBMS.3621,UKBMS,3495.98,3290.834\nUKBMS.3622,UKBMS,3504.414,3336.027\nUKBMS.3623,UKBMS,3484.368,3288.197\nUKBMS.3627,UKBMS,3487.246,3249.329\nUKBMS.3628,UKBMS,3457.915,3344.809\nUKBMS.3629,UKBMS,3499.017,3262.938\nUKBMS.3632,UKBMS,3484.296,3281.8\nUKBMS.3633,UKBMS,3480.65,3294.632\nUKBMS.3636,UKBMS,3490.286,3318.539\nUKBMS.3640,UKBMS,3498.348,3291.26\nUKBMS.3641,UKBMS,3489.588,3317.836\nUKBMS.3643,UKBMS,3504.685,3398.318\nUKBMS.3644,UKBMS,3507.575,3333.019\nUKBMS.3645,UKBMS,3484.101,3326.51\nUKBMS.3646,UKBMS,3455.298,3303.668\nUKBMS.3647,UKBMS,3502.674,3300.745\nUKBMS.3648,UKBMS,3479.985,3286.031\nUKBMS.3650,UKBMS,3507.065,3328.877\nUKBMS.3651,UKBMS,3510.067,3326.851\nUKBMS.3652,UKBMS,3494.528,3290.721\nUKBMS.37,UKBMS,3538.543,3571.198\nUKBMS.3800,UKBMS,3519.816,3140.361\nUKBMS.3801,UKBMS,3542.202,3142.065\nUKBMS.3802,UKBMS,3525.929,3155.736\nUKBMS.3803,UKBMS,3516.028,3155.22\nUKBMS.3804,UKBMS,3515.576,3155.75\nUKBMS.3805,UKBMS,3534.039,3166.275\nUKBMS.3806,UKBMS,3524.142,3178.077\nUKBMS.3814,UKBMS,3528.143,3198.037\nUKBMS.3815,UKBMS,3569.12,3164.819\nUKBMS.3816,UKBMS,3550.452,3146.08\nUKBMS.3817,UKBMS,3531.496,3164.97\nUKBMS.3818,UKBMS,3538.376,3152.654\nUKBMS.3819,UKBMS,3519.231,3146.384\nUKBMS.3820,UKBMS,3518.93,3140.237\nUKBMS.3824,UKBMS,3522.159,3152.501\nUKBMS.3825,UKBMS,3551.092,3146.279\nUKBMS.3826,UKBMS,3517.397,3177.674\nUKBMS.3827,UKBMS,3525.183,3160.431\nUKBMS.3828,UKBMS,3520.797,3162.683\nUKBMS.3830,UKBMS,3544.777,3197.541\nUKBMS.3831,UKBMS,3525.744,3170.598\nUKBMS.3833,UKBMS,3510.361,3141.049\nUKBMS.3834,UKBMS,3515.679,3155.741\nUKBMS.3836,UKBMS,3514.078,3167.354\nUKBMS.3837,UKBMS,3547.135,3184.95\nUKBMS.3840,UKBMS,3526.845,3147.152\nUKBMS.3841,UKBMS,3523.923,3136.907\nUKBMS.3842,UKBMS,3521.506,3140.434\nUKBMS.3843,UKBMS,3498.789,3188.793\nUKBMS.3846,UKBMS,3514.139,3141.699\nUKBMS.3847,UKBMS,3514.263,3142.013\nUKBMS.3848,UKBMS,3511.035,3147.93\nUKBMS.3849,UKBMS,3511.353,3144.525\nUKBMS.3852,UKBMS,3502.477,3134.941\nUKBMS.3854,UKBMS,3499.302,3158.734\nUKBMS.3855,UKBMS,3511.324,3143.835\nUKBMS.3856,UKBMS,3532.971,3166.554\nUKBMS.3858,UKBMS,3522.344,3186.199\nUKBMS.3861,UKBMS,3545.999,3184.199\nUKBMS.3862,UKBMS,3523.069,3183.844\nUKBMS.3863,UKBMS,3521.46,3183.908\nUKBMS.3865,UKBMS,3519.093,3160.934\nUKBMS.3866,UKBMS,3523.839,3159.029\nUKBMS.3870,UKBMS,3511.018,3134.472\nUKBMS.3874,UKBMS,3570.175,3188.411\nUKBMS.3876,UKBMS,3503.96,3140.892\nUKBMS.3878,UKBMS,3528.835,3145.095\nUKBMS.3879,UKBMS,3497.738,3151.678\nUKBMS.3880,UKBMS,3497.738,3151.678\nUKBMS.39,UKBMS,3672.123,3409.386\nUKBMS.3900,UKBMS,3544.095,3309.382\nUKBMS.3901,UKBMS,3543.701,3309.448\nUKBMS.3902,UKBMS,3543.381,3297.816\nUKBMS.3906,UKBMS,3544.571,3308.591\nUKBMS.3911,UKBMS,3519.637,3310.009\nUKBMS.3912,UKBMS,3538.52,3304.254\nUKBMS.3913,UKBMS,3552.891,3305.526\nUKBMS.3914,UKBMS,3543.189,3299.27\nUKBMS.3915,UKBMS,3544.21,3309.871\nUKBMS.3918,UKBMS,3545.752,3313.17\nUKBMS.3919,UKBMS,3546.937,3313.978\nUKBMS.3920,UKBMS,3547.824,3316.558\nUKBMS.3921,UKBMS,3547.11,3316.678\nUKBMS.3922,UKBMS,3548.014,3317.176\nUKBMS.3923,UKBMS,3551.718,3306.636\nUKBMS.3924,UKBMS,3530.574,3299.591\nUKBMS.3925,UKBMS,3548.948,3316.599\nUKBMS.3927,UKBMS,3547.405,3301.243\nUKBMS.3936,UKBMS,3531.153,3338.701\nUKBMS.3937,UKBMS,3531.854,3339.473\nUKBMS.3938,UKBMS,3541.845,3308.03\nUKBMS.3939,UKBMS,3547.168,3289.089\nUKBMS.3940,UKBMS,3558.057,3313.577\nUKBMS.3941,UKBMS,3549.408,3301.789\nUKBMS.4,UKBMS,3626.793,3304.255\nUKBMS.4002,UKBMS,3627.738,3237.741\nUKBMS.4009,UKBMS,3611.91,3236.314\nUKBMS.4012,UKBMS,3588.515,3239.294\nUKBMS.4019,UKBMS,3609.759,3239.211\nUKBMS.4020,UKBMS,3597.818,3221.123\nUKBMS.4021,UKBMS,3589.773,3173.701\nUKBMS.4022,UKBMS,3618.132,3241.473\nUKBMS.4023,UKBMS,3619.824,3203.417\nUKBMS.4027,UKBMS,3629.447,3230.348\nUKBMS.4200,UKBMS,3411.408,3725.96\nUKBMS.4206,UKBMS,3546.485,4016.949\nUKBMS.4207,UKBMS,3526.594,3720.356\nUKBMS.4210,UKBMS,3426.469,3725.321\nUKBMS.4211,UKBMS,3438.88,3729.928\nUKBMS.4214,UKBMS,3528.185,3731.229\nUKBMS.4221,UKBMS,3436.618,3713.362\nUKBMS.4222,UKBMS,3497.233,3879.092\nUKBMS.4224,UKBMS,3490.298,3958.495\nUKBMS.4226,UKBMS,3436.147,3708.867\nUKBMS.43,UKBMS,3490.454,3525.074\nUKBMS.4301,UKBMS,3442.143,2972.147\nUKBMS.4306,UKBMS,3439.359,2974.784\nUKBMS.4314,UKBMS,3442.604,2974.203\nUKBMS.4315,UKBMS,3442.674,2972.867\nUKBMS.4334,UKBMS,3441.775,2976.999\nUKBMS.44,UKBMS,3443.858,3131.234\nUKBMS.4402,UKBMS,3455.418,3160.739\nUKBMS.4404,UKBMS,3432.564,3160.166\nUKBMS.4405,UKBMS,3492.084,3146.825\nUKBMS.4506,UKBMS,3431.533,3209.605\nUKBMS.4511,UKBMS,3400.459,3188.173\nUKBMS.4512,UKBMS,3404.253,3183.577\nUKBMS.4514,UKBMS,3408.738,3177.647\nUKBMS.4516,UKBMS,3429.339,3195.261\nUKBMS.4517,UKBMS,3428.453,3195.408\nUKBMS.4518,UKBMS,3433.954,3209.74\nUKBMS.4521,UKBMS,3402.702,3177.228\nUKBMS.48,UKBMS,3618.513,3210.308\nUKBMS.51,UKBMS,3498.384,3871.181\nUKBMS.54,UKBMS,3457.999,3260.054\nUKBMS.59,UKBMS,3683.066,3302.337\nUKBMS.6,UKBMS,3356.022,3144.669\nUKBMS.60,UKBMS,3579.391,3681.301\nUKBMS.61,UKBMS,3479.997,3858.99\nUKBMS.62,UKBMS,3667.215,3192.09\nUKBMS.65,UKBMS,3630.608,3307.279\nUKBMS.69,UKBMS,3688.096,3359.524\nUKBMS.71,UKBMS,3615.372,3130.099\nUKBMS.76,UKBMS,3660.599,3288.608\nUKBMS.8,UKBMS,3546.089,3155.738\nUKBMS.83,UKBMS,3500.957,3206.509\nUKBMS.84,UKBMS,3523.508,3162.428\nUKBMS.85,UKBMS,3541.258,3405.5\nUKBMS.86,UKBMS,3541.964,3169.35\nUKBMS.89,UKBMS,3565.128,3147.201\nUKBMS.9,UKBMS,3561.886,3142.632\nUKBMS.91,UKBMS,3601.468,3139.253\nUKBMS.92,UKBMS,3297.858,3572.485\nUKBMS.93,UKBMS,3487.501,3162.805\nUKBMS.95,UKBMS,3484.354,3324.355\nUKBMS.96,UKBMS,3521.614,3395.781\nUKBMS.989,UKBMS,3557.448,3823.701\nUKBMS.99,UKBMS,3643.051,3386.238\nUKBMS.991,UKBMS,3509.654,3177.648\nUKBMS.992,UKBMS,3552.05,3663.026\nUKBMS.993,UKBMS,3526.528,3578.199\nUKBMS.994,UKBMS,3542.961,3246.482\nUKBMS.995,UKBMS,3346.638,3166.148\nUKBMS.996,UKBMS,3286.982,3598.065\nUKBMS.997,UKBMS,3559.319,3821.353\nUKBMS.998,UKBMS,3520.616,3296.549\nUKBMS.999,UKBMS,3564.66,3172.273\n", stringsAsFactors=FALSE, colClasses=c('character','character','numeric','numeric'))

list(run=run,watch=watch,status=status,stop=stop_job)
})
run_spatial_trends <- pheno_spatial_trends$run
watch_spatial_trends <- pheno_spatial_trends$watch
spatial_trends_status <- pheno_spatial_trends$status
stop_spatial_trends <- pheno_spatial_trends$stop
if(!isTRUE(getOption('pheno.trends.spatial.functions_only',FALSE))) {
  spatial_trend_job <- run_spatial_trends(spatial_trend_config)
  watch_spatial_trends(spatial_trend_job)
}
