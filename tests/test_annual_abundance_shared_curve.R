source('annual_abundance_shared_curve.R')
api <- shared_curve_abundance
shape <- rep(1/365,365)
stopifnot(isTRUE(all.equal(api$annual_index(c(2,3),c(100,110),shape),912.5)))
stopifnot(api$annual_index(c(0,0),c(100,110),shape)==0)
stopifnot(!api$coverage(100:110,shape)$coverage_ok)
stopifnot(api$coverage(seq(1,365,by=7),shape)$coverage_ok)
stopifnot(!api$coverage(c(1:10,355:365),shape)$coverage_ok)
# End-to-end fixture: one species, ten adequately visited years, sparse and zero years.
root <- tempfile('shared_abundance_test_'); dir.create(root); dir.create(file.path(root,'data')); dir.create(file.path(root,'output'))
vis <- expand.grid(year=2000:2009,doy=seq(40,320,by=7))
vis$visit_date <- as.character(as.Date(paste0(vis$year,'-01-01'))+vis$doy-1)
vis$transect_id <- 'TEST.1'; vis$visit_id <- seq_len(nrow(vis))
cnt <- vis; cnt$species_name <- 'Test species'
cnt$count <- round(20*exp(-.5*((cnt$doy-180)/25)^2))
cnt$count[cnt$year==2008 & cnt$doy!=180] <- 0
cnt$count[cnt$year==2009] <- 0
cnt <- cnt[cnt$count>0,]
write.csv(vis,file.path(root,'data','ebms_visit.csv'),row.names=FALSE)
write.csv(cnt,file.path(root,'data','ebms_count.csv'),row.names=FALSE)
write.csv(data.frame(transect_id='TEST.1'),file.path(root,'data','ebms_transect_climate.csv'),row.names=FALSE)
pf <- file.path(root,'pop.csv'); write.csv(data.frame(SPECIES='Test species',SITE_ID='TEST.1'),pf,row.names=FALSE)
out <- api$run(root,population_file=pf)
a <- read.csv(file.path(out$output_dir,'annual_abundance.csv'))
stopifnot(out$successful==1L,nrow(a)==10L,a$n_positive_dates[a$YEAR==2008]==1L,
          a$ABUND_INDEX[a$YEAR==2008]>0,a$ABUND_INDEX[a$YEAR==2009]==0,
          all(a$coverage_ok))
cat('Shared-curve abundance checks passed.\n')
