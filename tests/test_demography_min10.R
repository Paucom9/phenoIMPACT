# Rscript --vanilla tests/test_demography_min10.R [--real]
# No optimization, DLL evaluation, process launch, or writes to baseline outputs.
source('demography_min10_observed_years.R')
source('post_environment_mundlak.R')
e<-environment(pheno_demo_min10$run); h<-environment(pheno_post_env$run)
check<-function(x,msg) if(!isTRUE(x)) stop(msg,call.=FALSE)
fails<-function(expr,pattern) {
  msg<-tryCatch({force(expr); NA_character_},error=function(e) conditionMessage(e))
  check(!is.na(msg)&&grepl(pattern,msg),paste('Expected:',pattern,'received:',msg))
}
check(identical(names(pheno_demo_min10),c('preflight','run','status','review')),'API mismatch')
check(identical(formals(e$run)$mode,quote(c('prepare','fit'))),'Preparation must be default')
run_text<-paste(deparse(e$run),collapse=' ')
check(all(vapply(c('cleanup = FALSE','cleanup_tree = FALSE','poll_connection = FALSE'),
  function(x) grepl(x,run_text,fixed=TRUE),logical(1))),
  'Detached controller must survive the launching R session')
set.seed(20261023)
species<-paste0('sp',c(4,2,7,1,5,3,6))
sites<-data.frame(SITE_ID=paste0('site',c(3,1,5,2,4)),x_km=1:5,y_km=11:15,
  latitude_deg=runif(5,40,55),longitude_deg=runif(5,-5,15))
pop<-expand.grid(SPECIES=species,SITE_ID=sites$SITE_ID,stringsAsFactors=FALSE)
pop$pop_id<-paste(pop$SPECIES,pop$SITE_ID,sep='_'); pop$bms_id<-'fixture'
pop$onset_plasticity_z<-rnorm(nrow(pop)); pop$offset_plasticity_bio_z<-rnorm(nrow(pop))
d<-do.call(rbind,lapply(seq_len(nrow(pop)),function(i) {
  # Exclude a whole species and a whole site; leave native maps nonalphabetical.
  y<-if(pop$SPECIES[i]==species[3]||pop$SITE_ID[i]==sites$SITE_ID[2]) 1990:1998 else 1990:1999
  if(i==1) y<-seq.int(1990L,2017L,by=3L) # ten nonconsecutive years retained
  z<-pop[rep(i,length(y)),]; z$YEAR<-y; z$year_decade<-(y-2000)/10
  z$ABUND_INDEX<-exp(rnorm(length(y))); z
}))
dec<-h$decompose(d,sites); d<-dec$frame; X<-dec$X
dat<-list(y=d$ABUND_INDEX,X=unname(X),time=d$year_decade,
  species=match(d$SPECIES,species)-1L,site=match(d$SITE_ID,sites$SITE_ID)-1L,
  year=match(d$YEAR,1990:2023)-1L,previous=rep(-1L,nrow(d)),gap=rep(1L,nrow(d)),
  A=matrix(seq_len(15),5,3),M0=diag(3),M1=diag(3),M2=diag(3),use_spatial=1L)
old<-list(signature='fixture',group='univoltine',frame=d,X=X,data=dat,sites=sites,years=1990:2023,
  scaling=data.frame(analysis='fixture',mean=2,SD=3),parameters=list())
beta<-setNames(rnorm(ncol(X),sd=.03),colnames(X))
pars<-list(beta=unname(beta),b_species=matrix(seq_len(14),7,2),b_site=matrix(seq_len(10),5,2),
  ar_obs=seq_len(nrow(d))/100,field=matrix(seq_len(102),3,34),log_sd=rep(0,5))
fit<-list(group='univoltine',signature='fixture',beta=beta,parameters=pars,
  checks=data.frame(numerical_checks_ok=TRUE),optimizer=list(objective=100))
p<-e$subset_input(old,fit,10L,h,'run-fixture')
audit<-e$counts(p$frame)
check(all(audit$observed_years>=10L),'Short populations retained')
check(nrow(p$populations)==24L,'Expected 6 species x 4 sites')
check(pop$pop_id[1] %in% p$frame$pop_id,'Ten gapped observations incorrectly excluded')
check(any(vapply(split(p$frame$YEAR,p$frame$pop_id),function(y) diff(range(y))==9,logical(1))),
  'Ten consecutive observations must be eligible despite span=9')
check(identical(p$scaling,old$scaling),'Original plasticity scaling changed')
idx<-p$retained_source_rows
check(identical(p$frame$ABUND_INDEX,old$frame$ABUND_INDEX[idx]),'Response changed')
check(identical(p$frame$year_decade,old$frame$year_decade[idx]),'Time scale changed')
check(identical(p$frame$onset_plasticity_z,old$frame$onset_plasticity_z[idx]),'Plasticity z changed')
check(identical(p$parameters$ar_obs,fit$parameters$ar_obs[idx]),'AR1 states mismapped')
check(identical(p$parameters$field,fit$parameters$field),'Spatial mesh/year field changed')
sp_ids<-unique(data.frame(label=p$frame$SPECIES,index=p$data$species))
si_ids<-unique(data.frame(label=p$frame$SITE_ID,index=p$data$site))
check(identical(sort(sp_ids$index),0:5),'Species indices not contiguous')
check(identical(sort(si_ids$index),0:3),'Site indices not contiguous')
for(i in seq_len(nrow(sp_ids))) check(identical(unname(p$parameters$b_species[sp_ids$index[i]+1L,]),
  unname(fit$parameters$b_species[match(sp_ids$label[i],species),])),'Species warm state remapping failed')
for(i in seq_len(nrow(si_ids))) {
  oldi<-match(si_ids$label[i],sites$SITE_ID); newi<-si_ids$index[i]+1L
  check(identical(unname(p$parameters$b_site[newi,]),unname(fit$parameters$b_site[oldi,])),
    'Site warm state remapping failed')
  check(identical(unname(p$data$A[newi,]),unname(old$data$A[oldi,])),'Spatial projection remapping failed')
}
edge<-which(p$data$previous>=0L); prev<-p$data$previous[edge]+1L
check(all(p$frame$pop_id[edge]==p$frame$pop_id[prev]),'AR1 crosses populations')
check(all(p$data$gap[edge]==p$frame$YEAR[edge]-p$frame$YEAR[prev]),'Calendar gaps incorrect')
check(any(p$data$gap==3L),'Gapped calendar years lost')
for(metric in c('onset','offset')) {
  w<-p$populations[[paste0(metric,'_plasticity_within_species_z')]]
  check(max(abs(tapply(w,p$populations$SPECIES,mean)))<1e-12,'Within component not re-centred')
}
check(identical(p$plasticity_uncertainty_propagated,FALSE),'Unrequested uncertainty propagation')
check(!identical(p$signature,old$signature),'Filtered signature reuses baseline signature')
bad<-old; bad$frame<-rbind(bad$frame,bad$frame[1,]); fails(e$counts(bad$frame),'duplicate')
badfit<-fit; badfit$signature<-'wrong'; fails(e$subset_input(old,badfit,10,h,'x'),'identity')
badfit<-fit; badfit$checks$numerical_checks_ok<-FALSE
fails(e$subset_input(old,badfit,10,h,'x'),'numerical review')
fails(e$subset_input(old,fit,40,h,'x'),'No population')
cat('PASS: observed-year eligibility, signs/scales, recentering, AR1 gaps, latent indices, mesh, identity guards\n')
if('--real' %in% commandArgs(trailingOnly=TRUE)) {
  b<-file.path(e$ROOT,e$REL); helper<-e$helpers(file.path(b,'code/post_environment_mundlak.R'))
  baseline_marker_before<-unname(tools::md5sum(file.path(b,'progress.rds')))
  expected<-list(multivoltine=c(117796,7475,64,1694),univoltine=c(251015,15146,129,1678))
  for(g in e$GROUPS) {
    old<-readRDS(file.path(b,paste0('prepared_mundlak_',g,'.rds')))
    fit<-readRDS(file.path(b,paste0(g,'__mundlak/fit.rds')))
    p<-e$subset_input(old,fit,10L,helper,'real-integration-test')
    row<-e$sample_row(p$frame,g,'min10')
    check(identical(as.numeric(row[1,c('observations','populations','species','sites')]),expected[[g]]),
      paste(g,'real sample does not match independently audited counts'))
    check(all(e$counts(p$frame)$observed_years>=10),'Real filtered sample retains short populations')
    cat('PASS: real preparation',g,nrow(p$frame),'rows;',nrow(p$populations),'populations\n')
  }
  # Exercise exports/review in an explicitly temporary test directory using
  # stand-in estimates, never saving these as actual sensitivity fit results.
  cfg<-readRDS(file.path(b,'run_config.rds'))
  cfg$baseline<-b; cfg$signature<-'export-smoke-test'
  cfg$out<-tempfile('min10_EXPORT_TEST_ONLY_'); dir.create(cfg$out)
  cfg$progress<-file.path(cfg$out,'progress.rds')
  cfg$manifest<-e$manifest(c(file.path(b,'run_config.rds'),cfg$engine,cfg$script))
  eng<-e$engine(cfg)
  for(g in e$GROUPS) {
    old<-readRDS(file.path(b,paste0('prepared_mundlak_',g,'.rds')))
    f<-readRDS(file.path(b,paste0(g,'__mundlak/fit.rds')))
    p<-e$subset_input(old,f,10L,helper,'export-test')
    f$signature<-p$signature; f$version<-e$VERSION
    f$parameters<-p$parameters; f$eta<-f$eta[p$retained_source_rows]
    f$checks$nobs<-nrow(p$frame)
    saveRDS(p,file.path(cfg$out,paste0('prepared_mundlak_',g,'.rds')))
    dir.create(file.path(cfg$out,paste0(g,'__mundlak')))
    saveRDS(f,file.path(cfg$out,paste0(g,'__mundlak/fit.rds')))
    check(identical(eng$export_fitted(cfg,g,p,f),'FIT_OK'),'Existing fit exporter incompatible')
    diag<-helper$diagnostic_tables(p,f)
    rd<-file.path(cfg$out,paste0(g,'__mundlak/residual_diagnostics')); dir.create(rd)
    for(n in c('distribution','temporal')) e$csv(diag[[n]],file.path(rd,paste0(n,'.csv')))
    e$csv(data.frame(YEAR=1990L,Moran_I=.01,p_two_sided=.5),file.path(rd,'annual_spatial_Moran.csv'))
  }
  saveRDS(cfg,file.path(cfg$out,'run_config.rds'))
  e$review(list(out=cfg$out))
  check(all(file.exists(file.path(cfg$out,c('comparison_with_baseline.csv',
    'comparison_residual_screens.csv','comparison_onset_offset_balance.csv','review_fit_checks.csv')))),
    'Review exports missing')
  comp<-read.csv(file.path(cfg$out,'comparison_with_baseline.csv'))
  check(nrow(comp)==12L&&all(comp$same_sign),'Stand-in estimates not matched by term/group')
  check(identical(unname(tools::md5sum(file.path(b,'progress.rds'))),baseline_marker_before),
    'Export tests modified the baseline progress marker')
  cat('PASS: existing effect/contrast exporters and automatic numerical/residual/balance comparisons\n')
}
