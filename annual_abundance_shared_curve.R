# Separate, opt-in sensitivity workflow. source() only defines functions.
# Annual indices retain sparse years; phenological estimates are never modified.
shared_curve_abundance <- local({
  need <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)
  annual_index <- function(count, doy, shape) {
    exposure <- sum(shape[doy])
    if (!is.finite(exposure) || exposure <= 0) return(NA_real_)
    sum(count) / exposure
  }
  coverage <- function(doy, shape, min_visits=10L, min_mass=.8, max_gap=28L) {
    days <- sort(unique(doy)); cumulative <- cumsum(shape)/sum(shape)
    lo <- which(cumulative >= .05)[1]; hi <- which(cumulative >= .95)[1]
    mass <- sum(shape[seq.int(min(days), max(days))]) / sum(shape)
    core <- sort(unique(c(lo, days[days >= lo & days <= hi], hi)))
    gap <- if(length(core)>1L) max(diff(core)) else Inf
    data.frame(coverage_mass=mass, core_max_gap_days=gap,
      coverage_ok=length(days)>=min_visits && mass>=min_mass && gap<=max_gap)
  }
  run <- function(root='.', population_file=NULL, output_dir=NULL,
                  minimum_index_years=10L, min_visits=10L,
                  min_mass=.8, max_gap=28L, k=20L, max_populations=Inf) {
    need(requireNamespace('data.table',quietly=TRUE) && requireNamespace('mgcv',quietly=TRUE),
         'Install data.table and mgcv.')
    root <- normalizePath(root, winslash='/', mustWork=TRUE)
    if(is.null(output_dir)) output_dir <- file.path(root,'output',paste0('annual_abundance_shared_curve_',format(Sys.time(),'%Y%m%d_%H%M%S')))
    need(!dir.exists(output_dir),'Output directory already exists; choose a new directory.')
    need(min_mass>0 && min_mass<=1 && max_gap>0 && min_visits>=3 && minimum_index_years>=1,
         'Invalid sampling thresholds.')
    count <- data.table::fread(file.path(root,'data','ebms_count.csv'))
    visits <- data.table::fread(file.path(root,'data','ebms_visit.csv'))
    climate <- data.table::fread(file.path(root,'data','ebms_transect_climate.csv'))
    for(d in list(count,visits)) {
      d[, DATE := as.Date(visit_date)]
      d[, YEAR := as.integer(year)]
      d[, SITE_ID := as.character(transect_id)]
      d[, doy := as.integer(format(DATE,'%j'))]
      # Match original 365-day GAM grid: leap-day-year day 366 is mapped to 365.
      d[doy==366L,doy:=365L]
    }
    need(!anyNA(visits[,.(DATE,YEAR,SITE_ID,doy)]) && !anyNA(count[,.(DATE,YEAR,SITE_ID,doy)]),
         'Missing visit/count identities or dates.')
    need(all(is.finite(count$count) & count$count>=0),'Invalid raw counts.')
    need(all(as.integer(format(visits$DATE,'%Y'))==visits$YEAR) &&
         all(as.integer(format(count$DATE,'%Y'))==count$YEAR),'Date/year mismatch.')
    visits <- visits[SITE_ID %in% climate$transect_id]
    effort <- visits[,.(n_visits=data.table::uniqueN(visit_id)),by=.(SITE_ID,YEAR)]
    sites <- effort[n_visits>=min_visits, .N,by=SITE_ID][N>=10L,SITE_ID]
    effort <- effort[SITE_ID %in% sites & n_visits>=min_visits]
    # One count unit per site-date; duplicate visits on a date require clarification.
    visit_ids <- unique(visits[,.(SITE_ID,YEAR,DATE,visit_id)])
    duplicate_dates <- visit_ids[, .N,by=.(SITE_ID,YEAR,DATE)][N>1L]
    duplicate_dates <- merge(duplicate_dates,effort[,.(SITE_ID,YEAR)],by=c('SITE_ID','YEAR'))
    visits <- unique(visits[,.(SITE_ID,YEAR,DATE,doy)])
    visits <- merge(visits,effort,by=c('SITE_ID','YEAR'))
    count[, SPECIES := as.character(species_name)]
    count <- count[,.(COUNT=sum(count)),by=.(SPECIES,SITE_ID,YEAR,DATE,doy)]
    if(is.null(population_file)) {
      old <- data.table::fread(file.path(root,'output','pheno_abund_estimates_allspp.csv'))
      pops <- unique(old[is.finite(ABUND_INDEX) & ABUND_INDEX>0,.(SPECIES,SITE_ID,YEAR)])[
        ,.N,by=.(SPECIES,SITE_ID)][N>=minimum_index_years,.(SPECIES,SITE_ID)]
    } else pops <- unique(data.table::fread(population_file)[,.(SPECIES,SITE_ID)])
    pops <- pops[SITE_ID %in% sites]
    if(is.finite(max_populations)) pops <- head(pops,max_populations)
    need(nrow(pops)>0,'No eligible populations.')
    dir.create(output_dir,recursive=TRUE)
    saveRDS(list(minimum_index_years=minimum_index_years,min_visits=min_visits,
      min_mass=min_mass,max_gap=max_gap,k=k,population_file=population_file,
      method='Pooled NB GAM shape; annual sum-count / sampled-shape exposure',
      uncertainty='Point estimates only; shape and abundance uncertainty not propagated'),
      file.path(output_dir,'config.rds'))
    results <- list(); failures <- list()
    for(i in seq_len(nrow(pops))) {
      sp <- pops$SPECIES[i]; site <- pops$SITE_ID[i]
      message(i,'/',nrow(pops),' ',sp,' / ',site)
      result <- tryCatch({
        need(!site %in% duplicate_dates$SITE_ID,
             'Multiple distinct visit IDs per eligible site-date: resolve sampling units before fitting.')
        vv <- visits[SITE_ID==site]
        cc <- count[SITE_ID==site & SPECIES==sp & YEAR %in% vv$YEAR]
        need(!nrow(cc[!vv,on=.(SITE_ID,YEAR,DATE)]), 'Count dates lack eligible visits; inspect sampling coverage.')
        d <- merge(vv,cc[,.(YEAR,DATE,COUNT)],by=c('YEAR','DATE'),all.x=TRUE)
        d[is.na(COUNT),COUNT:=0]
        totals <- d[,.(total=sum(COUNT)),by=YEAR]
        positive_years <- totals[total>0,YEAR]
        need(length(positive_years)>=2L,'Fewer than two years with detections.')
        train <- as.data.frame(d[YEAR %in% positive_years,.(YEAR,doy,COUNT)])
        # Preserve original winter-anchor assumption, now with a shared seasonal shape.
        anchors <- expand.grid(YEAR=positive_years,doy=c(1:30,335:365)); anchors$COUNT<-0
        train <- rbind(train,anchors); train$year_factor <- factor(train$YEAR)
        fit <- mgcv::gam(COUNT ~ year_factor + s(doy,k=k),data=train,
                         family=mgcv::nb(),method='REML')
        need(isTRUE(fit$converged),'GAM did not converge.')
        grid <- data.frame(doy=1:365,year_factor=factor(rep(levels(train$year_factor)[1],365),levels=levels(train$year_factor)))
        shape <- as.numeric(predict(fit,newdata=grid,type='response'))
        need(all(is.finite(shape)) && sum(shape)>0,'Invalid shared flight shape.')
        shape <- shape/sum(shape)
        years <- lapply(sort(unique(d$YEAR)),function(y) {
          z <- d[YEAR==y]; cv <- coverage(z$doy,shape,min_visits,min_mass,max_gap)
          estimate <- annual_index(z$COUNT,z$doy,shape)
          data.frame(SPECIES=sp,SITE_ID=site,YEAR=y,n_visit_dates=nrow(z),
            n_positive_dates=sum(z$COUNT>0),total_observed=sum(z$COUNT),
            index_raw=estimate,ABUND_INDEX=if(cv$coverage_ok) estimate else NA_real_,cv)
        })
        list(years=data.table::rbindlist(years),shape=data.frame(doy=1:365,relative_density=shape))
      },error=function(e)e)
      if(inherits(result,'error')) failures[[length(failures)+1L]] <- data.frame(SPECIES=sp,SITE_ID=site,error=conditionMessage(result))
      else {
        results[[length(results)+1L]] <- result$years
        # Per-population shape files permit independent coverage review.
        data.table::fwrite(result$shape,file.path(output_dir,paste0('shape_',i,'.csv')))
      }
      data.table::fwrite(pops[i],file.path(output_dir,paste0('population_',i,'.csv')))
      if(length(results)) data.table::fwrite(data.table::rbindlist(results),file.path(output_dir,'annual_abundance.csv'))
      if(length(failures)) data.table::fwrite(data.table::rbindlist(failures),file.path(output_dir,'failures.csv'))
    }
    writeLines(capture.output(sessionInfo()),file.path(output_dir,'sessionInfo.txt'))
    invisible(list(output_dir=output_dir,successful=length(results),failed=length(failures)))
  }
  list(run=run,annual_index=annual_index,coverage=coverage)
})
