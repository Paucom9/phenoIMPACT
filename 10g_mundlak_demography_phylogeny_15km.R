# phenoIMPACT | Phylogenetic sensitivity of demographic Mundlak models (10e/10f).
# Version 1.0.0 | 2026-10-15. Sourcing defines functions and starts no analysis.
# Models annual abundance (Gamma/log): multivoltine and univoltine.
# Reuses the completed FIT_OK source, exact 14 fixed terms, population AR1,
# species/site IID intercept and time slopes, and original annual 15-km SPDE.
# Adds only a phylogenetic species TIME slope alongside the existing IID slope.
# Stage-one plasticity predictors remain fixed, isolating covariance sensitivity.
# Original-root Brownian covariance is normalized to unit diagonal.
# No package updates. First native compilation on Windows requires matching Rtools.
#
# source("10g_mundlak_demography_phylogeny_15km.R", encoding = "UTF-8")
# audit_demography_mundlak_phylogeny_sources()
# job <- run_demography_mundlak_phylogeny(group = "multivoltine", check_only = TRUE)
# watch_demography_mundlak_phylogeny(job)
# job <- run_demography_mundlak_phylogeny(group = "multivoltine")
# Repeat the same run call to resume after restarting R.
# check_only runs tiny native preflight and validates source arrays, but no full fit.
#
# Portable checkpoints retain best finite evaluations and complete optimizer
# output before covariance calculation. A native call cannot resume mid-call.
# Numerical gate: optimizer convergence zero, positive definite covariance,
# and sqrt(g' V g) < 0.05. Raw gradient is also reported as an advisory check.
# Boundary variance tests are not assigned an automatic chi-square p value.

pheno_demography_mundlak_phylogeny <- local({
build_phylo_edge_basis <- function(tree, tip_names, tolerance = 1e-10) {
    need <- function(ok, message) if (!isTRUE(ok)) 
        stop(message, call. = FALSE)
    need(requireNamespace("ape", quietly = TRUE), "Cal el paquet ape.")
    need(inherits(tree, "phylo") && !inherits(tree, "multiPhylo") && ape::is.rooted(tree) && ape::is.ultrametric(tree), 
        "Cal un arbre arrelat i ultrametric.")
    need(is.character(tree$tip.label) && !anyNA(tree$tip.label) && all(nzchar(tree$tip.label)) && !anyDuplicated(tree$tip.label), 
        "Noms de l'arbre invalids.")
    need(is.character(tip_names) && length(tip_names) >= 2L && !anyNA(tip_names) && all(nzchar(tip_names)) && 
        !anyDuplicated(tip_names) && all(tip_names %in% tree$tip.label), "Matching d'especies invalid.")
    E <- tree$edge
    lengths <- tree$edge.length
    need(is.matrix(E) && ncol(E) == 2L && nrow(E) > 0L && length(lengths) == nrow(E) && all(is.finite(lengths)) && 
        all(lengths >= 0) && !anyDuplicated(E[, 2L]), "Branques invalides.")
    need(length(tolerance) == 1L && is.finite(tolerance) && tolerance > 0, "Tolerancia invalida.")
    roots <- setdiff(E[, 1L], E[, 2L])
    need(length(roots) == 1L, "Cal una sola arrel original.")
    root <- roots[[1L]]
    child_edge <- stats::setNames(seq_len(nrow(E)), as.character(E[, 2L]))
    descendants <- matrix(FALSE, nrow(E), length(tip_names))
    for (i in seq_along(tip_names)) {
        node <- match(tip_names[[i]], tree$tip.label)
        steps <- 0L
        while (node != root) {
            j <- unname(child_edge[as.character(node)])
            need(length(j) == 1L && !is.na(j), "El cami ancestral no arriba a l'arrel.")
            descendants[j, i] <- TRUE
            node <- E[j, 1L]
            steps <- steps + 1L
            need(steps <= nrow(E), "Cicle en el cami ancestral.")
        }
    }
    active <- which(rowSums(descendants) > 0L & lengths > 0)
    need(length(active) > 0L, "No hi ha increments ancestrals.")
    keys <- vapply(active, function(j) paste(which(descendants[j, ]), collapse = ","), character(1))
    key_order <- sort(unique(keys))
    groups <- lapply(key_order, function(k) active[keys == k])
    full_C <- ape::vcv.phylo(tree, corr = FALSE)
    C <- full_C[tip_names, tip_names, drop = FALSE]
    H <- mean(diag(C))
    need(all(is.finite(C)) && is.finite(H) && H > 0 && diff(range(diag(C)))/H < 1e-07, "Profunditats de l'arbre inconsistents.")
    A <- C/H
    clade_ids <- sprintf("clade_%04d", seq_along(groups))
    B <- matrix(0, length(tip_names), length(groups), dimnames = list(tip_names, clade_ids))
    clades <- vector("list", length(groups))
    edge_clade <- rep(NA_character_, nrow(E))
    for (j in seq_along(groups)) {
        rows <- groups[[j]]
        mask <- descendants[rows[[1L]], ]
        weight <- sum(lengths[rows])
        B[mask, j] <- sqrt(weight/H)
        edge_clade[rows] <- clade_ids[[j]]
        clades[[j]] <- data.frame(clade_id = clade_ids[[j]], descendant_key = key_order[[j]], n_matched_descendants = sum(mask), 
            summed_branch_length = weight, normalized_increment_variance = weight/H, n_original_edges = length(rows), 
            original_edge_rows = paste(rows, collapse = ","), stringsAsFactors = FALSE)
    }
    error <- max(abs(tcrossprod(B) - A))
    need(is.finite(error) && error <= tolerance, paste("Covariancia d'increments incorrecta:", error))
    need(!inherits(try(chol(A), silent = TRUE), "try-error"), "La covariancia no es definida positiva.")
    dense <- t(chol(A))
    list(B = B, A = A, H = H, tip_order = tip_names, clades = do.call(rbind, clades), original_edges = data.frame(edge_row = seq_len(nrow(E)), 
        parent_node = E[, 1L], child_node = E[, 2L], branch_length = lengths, n_matched_descendants = rowSums(descendants), 
        clade_id = edge_clade), audit = data.frame(n_original_tips = length(tree$tip.label), n_selected_tips = nrow(B), 
        n_original_edges = nrow(E), n_active_original_edges = length(active), n_increment_columns = ncol(B), 
        n_collapsed_edges = length(active) - ncol(B), original_root = root, common_original_tip_depth = H, 
        maximum_covariance_error = error, basis_rank = qr(B)$rank, edge_nonzeros = sum(B != 0), edge_density = mean(B != 
            0), dense_Cholesky_nonzeros = sum(dense != 0), dense_Cholesky_density = mean(dense != 0), 
        mean_edge_nonzeros_per_tip = mean(rowSums(B != 0)), mean_Cholesky_nonzeros_per_tip = mean(rowSums(dense != 
            0))))
}


# Native Gamma/log demographic-Mundlak phylogenetic sensitivity backend.
# Definitions only: sourcing never compiles, reads fitted data, or fits a model.

demo_phylo_cpp_source <- function() paste(c(
  '// Exact trend_spde likelihood plus an independent phylogenetic TIME slope.',
  '// B is a root-preserving, row-normalized ancestral-increment basis.',
  '// Its edge coefficients have a proper independent standard-normal prior.',
  '#include <TMB.hpp>',
  'template<class Type>',
  'Type objective_function<Type>::operator() () {',
  '  using namespace density;',
  '  DATA_VECTOR(y); DATA_MATRIX(X); DATA_VECTOR(time);',
  '  DATA_IVECTOR(species); DATA_IVECTOR(site); DATA_IVECTOR(year);',
  '  DATA_IVECTOR(previous); DATA_IVECTOR(gap);',
  '  DATA_SPARSE_MATRIX(A); DATA_SPARSE_MATRIX(M0);',
  '  DATA_SPARSE_MATRIX(M1); DATA_SPARSE_MATRIX(M2);',
  '  DATA_INTEGER(use_spatial); DATA_INTEGER(use_phylo);',
  '  DATA_SPARSE_MATRIX(B);',
  '  PARAMETER_VECTOR(beta); PARAMETER_VECTOR(log_sd);',
  '  PARAMETER(rho_raw); PARAMETER(log_shape);',
  '  PARAMETER_MATRIX(b_species); PARAMETER_MATRIX(b_site);',
  '  PARAMETER_VECTOR(ar_obs);',
  '  PARAMETER(log_range); PARAMETER(log_spatial_sd);',
  '  PARAMETER_MATRIX(field);',
  '  PARAMETER(log_sd_phylo); PARAMETER_VECTOR(z_phylo);',
  '  Type nll=0;',
  '  Type rho=rho_raw/sqrt(Type(1)+rho_raw*rho_raw);',
  '  Type shape=exp(log_shape);',
  '  vector<Type> sd=exp(log_sd);',
  '  for(int j=0;j<b_species.rows();j++)',
  '    for(int k=0;k<2;k++) nll-=dnorm(b_species(j,k),Type(0),sd(k),true);',
  '  for(int j=0;j<b_site.rows();j++)',
  '    for(int k=0;k<2;k++) nll-=dnorm(b_site(j,k),Type(0),sd(k+2),true);',
  '  // Integer-gap transitions analytically integrate out unobserved years.',
  '  // The first observation of each population has the stationary marginal.',
  '  for(int i=0;i<y.size();i++) {',
  '    if(previous(i)<0) nll-=dnorm(ar_obs(i),Type(0),sd(4),true);',
  '    else {',
  '      Type r=1;',
  '      for(int k=0;k<gap(i);k++) r*=rho;',
  '      nll-=dnorm(ar_obs(i),r*ar_obs(previous(i)),sd(4)*sqrt(Type(1)-r*r),true);',
  '    }',
  '  }',
  '  matrix<Type> spatial(A.rows(),field.cols()); spatial.setZero();',
  '  Type range=exp(log_range), spatial_sd=exp(log_spatial_sd);',
  '  if(use_spatial) {',
  '    Type kappa=sqrt(Type(8))/range;',
  '    Type k2=kappa*kappa;',
  '    Type tau=Type(1)/(sqrt(Type(4)*Type(M_PI))*kappa*spatial_sd);',
  '    Eigen::SparseMatrix<Type> Q=k2*k2*M0+Type(2)*k2*M1+M2;',
  '    auto gmrf=GMRF(Q);',
  '    for(int t=0;t<field.cols();t++) {',
  '      vector<Type> w=field.col(t);',
  '      nll+=SCALE(gmrf,Type(1)/tau)(w);',
  '    }',
  '    spatial=A*field;',
  '  }',
  '  vector<Type> phylo_species_slope(b_species.rows());',
  '  phylo_species_slope.setZero();',
  '  Type sd_phylo=0;',
  '  if(use_phylo) {',
  '    sd_phylo=exp(log_sd_phylo);',
  '    for(int j=0;j<z_phylo.size();j++)',
  '      nll-=dnorm(z_phylo(j),Type(0),Type(1),true);',
  '    phylo_species_slope=sd_phylo*(B*z_phylo);',
  '  }',
  '  vector<Type> eta=X*beta;',
  '  for(int i=0;i<y.size();i++) {',
  '    eta(i)+=b_species(species(i),0)+time(i)*b_species(species(i),1)',
  '      +b_site(site(i),0)+time(i)*b_site(site(i),1)+ar_obs(i);',
  '    if(use_spatial) eta(i)+=spatial(site(i),year(i));',
  '    if(use_phylo) eta(i)+=time(i)*phylo_species_slope(species(i));',
  '    nll-=dgamma(y(i),shape,exp(eta(i))/shape,true);',
  '  }',
  '  REPORT(eta); REPORT(spatial); REPORT(rho); REPORT(shape);',
  '  REPORT(range); REPORT(spatial_sd);',
  '  REPORT(phylo_species_slope); REPORT(sd_phylo);',
  '  return nll;',
  '}'
), collapse = '\n')

demo_phylo_compile <- function(directory) {
  if (!requireNamespace('TMB', quietly = TRUE)) stop('TMB is required.')
  if (!requireNamespace('digest', quietly = TRUE)) stop('digest is required.')
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  directory <- normalizePath(directory, winslash = '/', mustWork = TRUE)
  dll_name <- 'mundlak_demography_phylo_v1'
  source <- demo_phylo_cpp_source()
  cpp <- file.path(directory, paste0(dll_name, '.cpp'))
  lib <- TMB::dynlib(file.path(directory, dll_name))
  signature <- digest::digest(source, algo = 'sha256', serialize = FALSE)
  compilation_abi <- list(
    R_version = as.character(getRversion()), platform = R.version$platform,
    arch = R.version$arch, TMB_version = as.character(utils::packageVersion('TMB')),
    Matrix_version = as.character(utils::packageVersion('Matrix')),
    RcppEigen_version = as.character(utils::packageVersion('RcppEigen')),
    framework = 'CppAD', flags = '-O1 -g0', safebounds = TRUE,
    safeunload = TRUE, libtmb = TRUE, libinit = TRUE,
    supernodal = FALSE, longint = FALSE, eigen_disable_warnings = TRUE)
  compile_signature <- digest::digest(list(source_sha256 = signature,
                                           compilation_abi = compilation_abi), algo = 'sha256')
  signature_file <- file.path(directory, paste0(dll_name, '.sha256'))
  build_file <- file.path(directory, paste0(dll_name, '.build.rds'))
  existing <- if (file.exists(signature_file)) readLines(signature_file, warn = FALSE) else character()
  build <- if (file.exists(build_file)) tryCatch(readRDS(build_file), error = function(e) NULL) else NULL
  cache_matches <- identical(existing, signature) && is.list(build) &&
    identical(build$source_sha256, signature) &&
    identical(build$compile_signature, compile_signature) &&
    identical(build$compilation_abi, compilation_abi)
  loaded <- getLoadedDLLs()[[dll_name]]
  if (!is.null(loaded)) {
    same_path <- identical(normalizePath(loaded[['path']], winslash = '/', mustWork = TRUE),
                           normalizePath(lib, winslash = '/', mustWork = TRUE))
    if (!same_path || !cache_matches)
      stop('A different demographic phylogeny DLL is already loaded. Use a fresh R worker.')
  }
  if (!file.exists(lib) || !cache_matches) {
    writeLines(source, cpp, useBytes = TRUE)
    old_dir <- setwd(directory)
    on.exit(setwd(old_dir), add = TRUE)
    env_names <- c('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                   'VECLIB_MAXIMUM_THREADS', 'MAKEFLAGS')
    old_env <- Sys.getenv(env_names, unset = NA_character_)
    on.exit({
      missing_env <- is.na(old_env)
      if (any(missing_env)) Sys.unsetenv(env_names[missing_env])
      if (any(!missing_env)) do.call(Sys.setenv, as.list(setNames(old_env[!missing_env], env_names[!missing_env])))
    }, add = TRUE)
    do.call(Sys.setenv, as.list(setNames(c('1','1','1','1','-j1'), env_names)))
    result <- tryCatch(TMB::compile(basename(cpp), flags = compilation_abi$flags,
                                   framework = compilation_abi$framework), error = function(e)
      stop('Demographic phylogeny compilation failed. On Windows R 4.5, install Rtools45 ',
           '(https://cran.r-project.org/bin/windows/Rtools/) and restart R. Details: ',
           conditionMessage(e)))
    if (!identical(as.integer(result), 0L) || !file.exists(lib))
      stop('Demographic phylogeny compilation failed; inspect the worker log.')
    writeLines(signature, signature_file)
    saveRDS(list(source_sha256 = signature, compile_signature = compile_signature,
                 compilation_abi = compilation_abi), build_file)
  }
  if (is.null(getLoadedDLLs()[[dll_name]])) dyn.load(lib)
  TMB::openmp(1L, DLL = dll_name)
  invisible(list(DLL = dll_name, path = lib, cpp = cpp, source_sha256 = signature,
                 compile_signature = compile_signature, compilation_abi = compilation_abi))
}

demo_phylo_make_objective <- function(data, parameters, phylo = TRUE,
                                     DLL = 'mundlak_demography_phylo_v1',
                                     inner.control = list(maxit = 1000, trace = FALSE),
                                     random = TRUE) {
  need <- function(ok, message) if (!isTRUE(ok)) stop(message)
  need(requireNamespace('TMB', quietly = TRUE) && requireNamespace('Matrix', quietly = TRUE),
       'TMB and Matrix are required.')
  need(is.logical(phylo) && length(phylo) == 1L && !is.na(phylo), 'phylo must be TRUE or FALSE.')
  required_data <- c('y','X','time','species','site','year','previous','gap',
                     'A','M0','M1','M2','use_spatial')
  required_par <- c('beta','log_sd','rho_raw','log_shape','b_species','b_site',
                    'ar_obs','log_range','log_spatial_sd','field')
  need(all(required_data %in% names(data)), 'Incomplete original demographic data.')
  need(all(required_par %in% names(parameters)), 'Incomplete original demographic parameters.')
  n <- length(data$y)
  need(n > 0L && all(is.finite(data$y)) && all(data$y > 0), 'Gamma response must be finite and positive.')
  need(is.matrix(data$X) && nrow(data$X) == n && all(is.finite(data$X)) &&
         length(parameters$beta) == ncol(data$X), 'Invalid fixed design dimensions or values.')
  need(all(vapply(data[c('time','species','site','year','previous','gap')], length, integer(1)) == n),
       'Observation-index vectors must have response length.')
  need(all(is.finite(data$time)), 'Nonfinite time.')
  need(is.matrix(parameters$b_species) && ncol(parameters$b_species) == 2L &&
         is.matrix(parameters$b_site) && ncol(parameters$b_site) == 2L,
       'Expected separate IID species/site intercepts and time slopes.')
  need(length(parameters$log_sd) == 5L && length(parameters$ar_obs) == n,
       'Expected five IID SDs and one AR1 state per observed population-year.')
  need(all(vapply(parameters[c('rho_raw','log_shape','log_range','log_spatial_sd')],
                  length, integer(1)) == 1L), 'Expected scalar AR1, shape, range and spatial-SD parameters.')
  valid_integer <- function(x) is.numeric(x) && all(is.finite(x)) && all(x == as.integer(x))
  need(all(vapply(data[c('species','site','year','previous','gap')], valid_integer, logical(1))),
       'Observation indices and calendar gaps must be finite integers.')
  need(all(data$species >= 0L & data$species < nrow(parameters$b_species)) &&
         all(data$site >= 0L & data$site < nrow(parameters$b_site)), 'Species/site indices out of bounds.')
  need(all(data$previous >= -1L & data$previous < seq_len(n) - 1L) && all(data$gap > 0L),
       'AR1 transitions require preceding observations or -1 and positive calendar gaps.')
  previous_rows <- which(data$previous >= 0L)
  if (length(previous_rows)) {
    predecessors <- data$previous[previous_rows] + 1L
    need(all(data$site[previous_rows] == data$site[predecessors]) &&
           all(data$species[previous_rows] == data$species[predecessors]),
         'AR1 transition crosses population (species/site).')
    need(all(data$year[previous_rows] - data$year[predecessors] == data$gap[previous_rows]),
         'AR1 gap does not match the saved consecutive calendar-year indices.')
  }
  need(length(data$use_spatial) == 1L && data$use_spatial %in% 0:1, 'Invalid spatial flag.')
  need(is.matrix(parameters$field), 'Field must be a matrix.')
  if (data$use_spatial == 1L) {
    need(nrow(data$A) == nrow(parameters$b_site) && ncol(data$A) == nrow(parameters$field) &&
           all(data$year >= 0L & data$year < ncol(parameters$field)) &&
           all(vapply(data[c('M0','M1','M2')], function(x)
             identical(as.integer(dim(x)), rep(as.integer(nrow(parameters$field)), 2L)), logical(1))),
         'Spatial/annual field dimensions do not match saved data.')
  } else need(all(dim(parameters$field) == 0L), 'Nonspatial field must have zero dimensions.')
  need(all(vapply(parameters[required_par], function(x) all(is.finite(x)), logical(1))),
       'Nonfinite original parameters.')
  if (phylo) {
    need('B' %in% names(data) && inherits(data$B, 'sparseMatrix') &&
           nrow(data$B) == nrow(parameters$b_species) && ncol(data$B) > 0L,
         'B must have one row per saved species and at least one ancestral-increment column.')
    need(all(is.finite(data$B@x)), 'Nonfinite phylogenetic basis.')
    row_variances <- as.numeric(Matrix::rowSums(data$B ^ 2))
    need(all(abs(row_variances - 1) < 1e-8), 'Phylogenetic basis rows must have unit marginal variance.')
    if (is.null(parameters$log_sd_phylo)) parameters$log_sd_phylo <- log(0.1)
    if (is.null(parameters$z_phylo)) parameters$z_phylo <- numeric(ncol(data$B))
    need(length(parameters$log_sd_phylo) == 1L && is.finite(parameters$log_sd_phylo) &&
           length(parameters$z_phylo) == ncol(data$B) && all(is.finite(parameters$z_phylo)),
         'Invalid phylogenetic slope parameters.')
  } else {
    data$B <- Matrix::sparseMatrix(i = integer(), j = integer(),
                                  dims = c(nrow(parameters$b_species), 0L))
    parameters$log_sd_phylo <- 0
    parameters$z_phylo <- numeric()
  }
  data$use_phylo <- as.integer(phylo)
  for (nm in c('A','M0','M1','M2','B'))
    data[[nm]] <- methods::as(methods::as(methods::as(data[[nm]], 'dMatrix'),
                                       'generalMatrix'), 'TsparseMatrix')
  random_names <- c('b_species','b_site','ar_obs')
  map <- list()
  if (data$use_spatial == 1L) random_names <- c(random_names, 'field') else
    map <- c(map, list(log_range = factor(NA), log_spatial_sd = factor(NA)))
  if (phylo) random_names <- c(random_names, 'z_phylo') else
    map <- c(map, list(log_sd_phylo = factor(NA)))
  if (!isTRUE(random)) random_names <- NULL
  TMB::MakeADFun(data = data, parameters = parameters, random = random_names,
                map = if (length(map)) map else NULL, DLL = DLL, silent = TRUE,
                inner.control = inner.control)
}

# Bounded native compatibility check: six synthetic rows, two species/sites,
# four calendar years and two mesh vertices. No fitted user data are read and
# no model is optimized. Compilation is cached by exact source hash.
demo_phylo_native_preflight <- function(directory, write_report = TRUE) {
  engine <- demo_phylo_compile(directory)
  diagonal <- Matrix::Diagonal(2L)
  data <- list(y = c(.8,1.1,1.3,1.2,.9,1.05),
    X = cbind(1, rep(c(-.5,0,1),2)), time = rep(c(-.5,0,1),2),
    species = rep(0:1, each=3), site = rep(0:1, each=3),
    year = rep(c(0L,1L,3L),2), previous = c(-1L,0L,1L,-1L,3L,4L),
    gap = rep(c(1L,1L,2L),2), A = diagonal, M0 = diagonal,
    M1 = diagonal, M2 = diagonal, use_spatial = 1L,
    B = Matrix::Matrix(rbind(c(sqrt(.4),sqrt(.6),0),
                            c(sqrt(.4),0,sqrt(.6))), sparse=TRUE))
  pars <- list(beta = c(.1,-.1), log_sd = log(c(.3,.15,.2,.1,.17)),
    rho_raw = -.4/sqrt(1-.4^2), log_shape = log(4),
    b_species = matrix(0,2,2), b_site = matrix(0,2,2),
    ar_obs = numeric(6), log_range = log(3), log_spatial_sd = log(.15),
    field = matrix(0,2,4))
  objects <- list()
  on.exit(for (obj in objects) try(TMB::FreeADFun(obj),silent=TRUE),add=TRUE)
  independent <- demo_phylo_make_objective(data,pars,phylo=FALSE,DLL=engine$DLL)
  objects[[length(objects)+1L]] <- independent
  null_pars <- pars
  null_pars$log_sd_phylo <- log(1e-8)
  null_pars$z_phylo <- numeric(ncol(data$B))
  near_null <- demo_phylo_make_objective(data,null_pars,phylo=TRUE,DLL=engine$DLL)
  objects[[length(objects)+1L]] <- near_null
  v0 <- independent$fn(independent$par); v1 <- near_null$fn(near_null$par)
  g0 <- as.numeric(independent$gr(independent$par))
  g1 <- as.numeric(near_null$gr(near_null$par))
  pars$log_sd_phylo <- log(.2)
  pars$z_phylo <- c(.3,-.2,.1)
  joint <- demo_phylo_make_objective(data,pars,phylo=TRUE,DLL=engine$DLL,random=FALSE)
  objects[[length(objects)+1L]] <- joint
  joint$fn(joint$par)
  report <- joint$report(joint$env$last.par)
  bs <- as.numeric(.2*(data$B %*% pars$z_phylo))
  expected_eta <- as.numeric(data$X %*% pars$beta) + data$time*bs[data$species+1L]
  value <- c(abs(v0-v1), max(abs(g0-g1[seq_along(g0)])),
             abs(tail(g1,1)), abs(length(near_null$par)-length(independent$par)-1L),
             max(abs(as.numeric(report$eta)-expected_eta)))
  tolerance <- c(1e-7,1e-6,1e-7,0,1e-10)
  table <- data.frame(check=c('near_zero_phylo_Laplace_NLL',
    'near_zero_phylo_original_gradient','near_zero_phylo_logSD_score',
    'exactly_one_added_outer_parameter','nonzero_phylo_temporal_predictor'),
    absolute_error=value,tolerance=tolerance,
    passed=is.finite(value)&value<=tolerance,
    scope='six synthetic observations; no user data; no optimization',
    stringsAsFactors=FALSE)
  if (write_report) utils::write.csv(table,file.path(directory,'native_preflight.csv'),row.names=FALSE)
  if (!all(table$passed)) stop('Bounded demographic phylogeny native preflight failed; inspect native_preflight.csv.')
  invisible(table)
}


# Pure-R discovery and portable source preparation for the accepted 10e/10f model.
# No saved native objective, callback, DLL, or baseline optimizer is invoked here.

demo_phylo_need <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
}

demo_phylo_version <- function() "mundlak_geography_plasticity_qc_final_spatial_15km_2.0.2"

demo_phylo_groups <- function() c("multivoltine", "univoltine")

demo_phylo_expected_samples <- function() list(
  multivoltine = c(nobs = 152971L, npop = 15643L),
  univoltine = c(nobs = 304056L, npop = 27026L))

demo_phylo_X_names <- function() {
  main <- c("latitude_within_species_10deg", "longitude_within_species_10deg",
    "onset_plasticity_within_species_z", "onset_plasticity_species_mean_z",
    "offset_plasticity_within_species_z", "offset_plasticity_species_mean_z")
  c("(Intercept)", "year_decade", main, paste0("year_decade:", main))
}

demo_phylo_same <- function(x, y, tolerance = 1e-9) {
  is.numeric(x) && is.numeric(y) && identical(dim(x), dim(y)) && length(x) == length(y) &&
    all(is.finite(x)) && all(is.finite(y)) &&
    all(abs(as.numeric(x) - as.numeric(y)) <= tolerance * (1 + abs(as.numeric(y))))
}

demo_phylo_metadata <- function(file) {
  s <- file.info(file)
  demo_phylo_need(nrow(s) == 1L && is.finite(s$size) && s$size > 0 && s$size <= 2 * 1024^2,
    paste("Missing, empty, or oversized metadata:", file))
  x <- readRDS(file)
  demo_phylo_need(is.list(x), paste("Expected metadata list:", file))
  x
}

demo_phylo_audit_row <- function(run, group, note = "") {
  directory <- file.path(run, paste0(group, "__mundlak"))
  files <- c(source_checkpoint = file.path(directory, "fit.rds"),
    prepared_file = file.path(run, paste0("prepared_mundlak_", group, ".rds")),
    run_config_file = file.path(run, "run_config.rds"),
    fit_checks_file = file.path(directory, "fit_checks.csv"),
    result_file = file.path(run, "result.rds"))
  info <- file.info(files)
  exists <- stats::setNames(file.exists(files) & !is.na(info$isdir) & !info$isdir, names(files))
  issues <- if (nzchar(note)) note else character()
  errors <- character()
  read_meta <- function(name) {
    if (!exists[[name]]) return(NULL)
    tryCatch(demo_phylo_metadata(files[[name]]), error = function(e) {
      errors <<- c(errors, paste0(basename(files[[name]]), ": ", conditionMessage(e))); NULL
    })
  }
  config <- read_meta("run_config_file")
  result <- read_meta("result_file")
  checks <- if (exists[["fit_checks_file"]]) tryCatch({
    s <- file.info(files[["fit_checks_file"]])$size
    demo_phylo_need(is.finite(s) && s > 0 && s <= 2 * 1024^2, "Empty or oversized fit_checks.csv.")
    utils::read.csv(files[["fit_checks_file"]], stringsAsFactors = FALSE)
  }, error = function(e) {
    errors <<- c(errors, paste0("fit_checks.csv: ", conditionMessage(e))); NULL
  }) else NULL
  required <- c("source_checkpoint", "prepared_file", "run_config_file", "fit_checks_file")
  missing <- basename(files[required][!exists[required]])
  empty <- basename(files[required][exists[required] & (is.na(info$size[match(required, names(files))]) |
    info$size[match(required, names(files))] <= 0)])
  version <- if (is.character(config$version) && length(config$version) == 1L) config$version else NA_character_
  signature <- config$run_signature
  valid_signature <- is.character(signature) && length(signature) == 1L && !is.na(signature) &&
    grepl("^[a-f0-9]{64}$", signature)
  checks_ok <- is.data.frame(checks) && nrow(checks) == 1L &&
    identical(as.character(checks$group), group) && isTRUE(checks$numerical_checks_ok)
  if (length(missing)) issues <- c(issues, paste("missing files:", paste(missing, collapse = ", ")))
  if (length(empty)) issues <- c(issues, paste("empty files:", paste(empty, collapse = ", ")))
  if (length(errors)) issues <- c(issues, paste("read errors:", paste(errors, collapse = " | ")))
  if (!identical(version, demo_phylo_version())) issues <- c(issues, "not the accepted 2.0.2 demographic Mundlak baseline")
  if (!valid_signature) issues <- c(issues, "run_signature must be a 64-character SHA256")
  if (!is.null(result) && (!identical(result$signature, signature) ||
      !identical(result$version, config$version))) issues <- c(issues, "run/result signature or version mismatch")
  if (!checks_ok) issues <- c(issues, "one matching group row with numerical_checks_ok=TRUE is required")
  if (!exists[["fit_checks_file"]]) issues <- c(issues,
    "fit_checks.csv is not ready; saved-array metadata recovery can be reviewed separately without refitting")
  data.frame(group = group, source_directory = normalizePath(directory, winslash = "/", mustWork = FALSE),
    source_checkpoint = files[["source_checkpoint"]], prepared_file = files[["prepared_file"]],
    run_config_file = files[["run_config_file"]], fit_checks_file = files[["fit_checks_file"]],
    result_file = files[["result_file"]], version = version,
    run_signature = if (valid_signature) signature else NA_character_,
    numerical_checks_ok = checks_ok, saved_fit_read_for_checks = FALSE,
    source_checkpoint_bytes = as.numeric(info$size[[1L]]), prepared_bytes = as.numeric(info$size[[2L]]),
    missing_files = paste(missing, collapse = "; "), read_errors = paste(errors, collapse = " | "),
    accepted = !length(issues), reason = if (length(issues)) paste(unique(issues), collapse = "; ") else
      "metadata eligible; exact design, identities and portable arrays await preparation",
    validation_scope = "eligibility_only_not_full_array_validation", stringsAsFactors = FALSE)
}

demo_phylo_audit_sources <- function(root, group = NULL, source_run = NULL) {
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  groups <- if (is.null(group)) demo_phylo_groups() else group
  demo_phylo_need(is.character(groups) && length(groups) > 0 && !anyNA(groups) &&
    all(groups %in% demo_phylo_groups()), "Unknown demographic Mundlak group.")
  groups <- unique(groups)
  note <- ""
  if (is.null(source_run)) {
    parent <- file.path(root, "output/population_trends/mundlak_geography_plasticity_qc_final_spatial_15km")
    runs <- if (dir.exists(parent)) list.dirs(parent, recursive = FALSE, full.names = TRUE) else character()
    if (!length(runs)) { runs <- parent; note <- "no independent demographic Mundlak run folders found" }
  } else {
    demo_phylo_need(is.character(source_run) && length(source_run) == 1L && !is.na(source_run) &&
      nzchar(source_run), "Invalid explicit demographic Mundlak source folder.")
    path <- normalizePath(source_run, winslash = "/", mustWork = FALSE)
    supplied <- sub("__mundlak$", "", basename(path))
    if (basename(path) %in% paste0(demo_phylo_groups(), "__mundlak")) {
      if (is.null(group)) groups <- supplied
      if (!all(groups == supplied)) note <- paste("explicit source folder is for", supplied,
        "; requested", paste(groups, collapse = ","))
      runs <- dirname(path)
    } else runs <- path
  }
  rows <- lapply(runs, function(run) lapply(groups, function(g) demo_phylo_audit_row(run, g, note)))
  out <- do.call(rbind, unlist(rows, recursive = FALSE)); rownames(out) <- NULL
  attr(out, "validation_scope") <- "Metadata eligibility only; no fitted-model or prepared-data RDS is read."
  out
}

demo_phylo_resolve_source <- function(root, group, source_run = NULL) {
  demo_phylo_need(is.character(group) && length(group) == 1L && group %in% demo_phylo_groups(),
    "Select exactly one demographic Mundlak group.")
  audit <- demo_phylo_audit_sources(root, group, source_run)
  good <- which(audit$accepted)
  if (length(good) != 1L) stop(paste0("Expected exactly one accepted independent 14-term demographic Mundlak source for ",
    group, "; found ", length(good), ".\n", paste(paste(audit$source_directory, audit$reason, sep = " | "), collapse = "\n"),
    "\nInspect demo_phylo_audit_sources(root, group, source_run), and select an explicit source_run if needed. No baseline is fitted automatically."), call. = FALSE)
  as.list(audit[good, c("group", "source_directory", "source_checkpoint", "prepared_file", "run_config_file",
    "fit_checks_file", "result_file", "run_signature"), drop = FALSE])
}

demo_phylo_plain <- function(x) {
  if (is.environment(x) || is.function(x) || typeof(x) %in% c("externalptr", "weakref")) return(FALSE)
  if (isS4(x)) return(all(vapply(methods::slotNames(x), function(n) demo_phylo_plain(methods::slot(x, n)), logical(1))))
  if (is.list(x)) return(all(vapply(x, demo_phylo_plain, logical(1))))
  TRUE
}

demo_phylo_save <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile("portable-", tmpdir = dirname(file))
  saveRDS(x, tmp, compress = "gzip")
  if (file.exists(file)) unlink(file)
  demo_phylo_need(file.rename(tmp, file), paste("Cannot finalize:", file))
  invisible(file)
}

demo_phylo_mapping <- function(labels, index, count, what) {
  labels <- as.character(labels)
  demo_phylo_need(length(index) == length(labels) && is.numeric(index) && !anyNA(index) &&
    all(index == as.integer(index)) && all(index >= 0) && all(index < count), paste("Invalid", what, "indices."))
  pairs <- unique(data.frame(label = labels, index = as.integer(index), stringsAsFactors = FALSE))
  demo_phylo_need(nrow(pairs) == count && !anyDuplicated(pairs$label) && !anyDuplicated(pairs$index) &&
    identical(sort(pairs$index), seq_len(count) - 1L), paste("Ambiguous", what, "mapping."))
  pairs$label[order(pairs$index)]
}

demo_phylo_prepare_uncached <- function(cfg, packages = NULL, resolved = NULL) {
  demo_phylo_need(requireNamespace("Matrix", quietly = TRUE), "Matrix is required to validate sparse source arrays.")
  if (is.null(resolved)) resolved <- demo_phylo_resolve_source(cfg$project_root, cfg$group, cfg$source_directory)
  run <- demo_phylo_metadata(resolved$run_config_file)
  p <- readRDS(resolved$prepared_file)
  fit <- readRDS(resolved$source_checkpoint)
  version <- demo_phylo_version()
  demo_phylo_need(is.list(p) && is.list(fit) && identical(run$version, version) &&
    identical(p$version, version) && identical(fit$version, version) &&
    identical(p$group, cfg$group) && identical(fit$group, cfg$group) &&
    identical(p$signature, fit$signature) && identical(p$source_bridge_signature, fit$source_bridge_signature) &&
    identical(p$source_data_signature, fit$source_data_signature), "Demographic baseline fit/prepared/run identities disagree.")
  demo_phylo_need(is.character(p$signature) && length(p$signature) == 1L &&
    grepl("^[a-f0-9]{64}$", p$signature), "Invalid demographic prepared signature.")
  expected_signature <- digest::digest(list(version = version, group = cfg$group,
    source_signature = p$source_bridge_signature, X = p$X, response = p$data$y), algo = "sha256")
  demo_phylo_need(identical(p$signature, expected_signature), "Demographic prepared content does not reproduce its saved signature.")
  demo_phylo_need(is.data.frame(p$frame) && is.data.frame(p$populations) && is.data.frame(p$sites) &&
    is.list(p$data) && is.list(fit$parameters), "Missing demographic prepared frame/data/parameters.")
  d <- p$frame; data <- p$data; pars <- fit$parameters; n <- nrow(d)
  fields <- c("SPECIES", "SITE_ID", "pop_id", "YEAR", "year_decade", "ABUND_INDEX",
    "onset_plasticity_z", "offset_plasticity_bio_z", demo_phylo_X_names()[3:8])
  demo_phylo_need(all(fields %in% names(d)) && !anyNA(d[fields]) && n > 0L,
    "Demographic annual frame lacks required finite observations.")
  sample <- if (is.null(cfg$expected_samples)) demo_phylo_expected_samples()[[cfg$group]] else cfg$expected_samples[[cfg$group]]
  demo_phylo_need(is.numeric(sample) && all(c("nobs", "npop") %in% names(sample)) &&
    n == sample[["nobs"]] && length(unique(d$pop_id)) == sample[["npop"]],
    "Demographic source does not contain the reviewed FINAL-QC common sample.")
  X <- stats::model.matrix(~ year_decade * (latitude_within_species_10deg + longitude_within_species_10deg +
    onset_plasticity_within_species_z + onset_plasticity_species_mean_z +
    offset_plasticity_within_species_z + offset_plasticity_species_mean_z), data = d)
  demo_phylo_need(identical(colnames(X), demo_phylo_X_names()) && identical(colnames(p$X), colnames(X)) &&
    demo_phylo_same(X, p$X) && demo_phylo_same(unname(X), data$X) && qr(X)$rank == 14L,
    "Demographic source is not the exact full-rank 14-term Mundlak design.")
  demo_phylo_need(demo_phylo_same(data$y, as.numeric(d$ABUND_INDEX), 0) && all(data$y > 0) &&
    demo_phylo_same(data$time, as.numeric(d$year_decade), 0) &&
    length(p$years) > 0 && !anyDuplicated(p$years), "Demographic response/year scales changed.")
  required_parameters <- c("beta", "log_sd", "rho_raw", "log_shape", "b_species", "b_site", "ar_obs", "log_range", "log_spatial_sd", "field")
  demo_phylo_need(setequal(names(pars), required_parameters) && !anyDuplicated(names(pars)) &&
    all(vapply(pars, function(x) is.numeric(x) && all(is.finite(x)), logical(1))), "Unexpected or non-finite source parameter arrays.")
  demo_phylo_need(is.matrix(pars$b_species) && ncol(pars$b_species) == 2L &&
    is.matrix(pars$b_site) && ncol(pars$b_site) == 2L && length(pars$ar_obs) == n &&
    is.matrix(pars$field) && length(pars$beta) == 14L && length(pars$log_sd) == 5L &&
    all(vapply(pars[c("rho_raw", "log_shape", "log_range", "log_spatial_sd")], length, integer(1)) == 1L),
    "Demographic latent/outer parameter dimensions changed.")
  species <- demo_phylo_mapping(d$SPECIES, data$species, nrow(pars$b_species), "species")
  sites <- demo_phylo_mapping(d$SITE_ID, data$site, nrow(pars$b_site), "site")
  demo_phylo_need(is.numeric(data$year) && length(data$year) == n && !anyNA(data$year) &&
    all(data$year == as.integer(data$year)) && all(data$year >= 0) &&
    all(data$year < ncol(pars$field)) && length(p$years) == ncol(pars$field), "Invalid demographic year indices.")
  demo_phylo_need(identical(sites, as.character(p$sites$SITE_ID)) &&
    identical(as.character(d$YEAR), as.character(p$years[data$year + 1L])),
    "Prepared site/year order differs from fitted random-effect indices.")
  demo_phylo_need(all(c("A", "M0", "M1", "M2", "use_spatial", "previous", "gap") %in% names(data)) &&
    identical(as.integer(data$use_spatial), 1L) && nrow(data$A) == length(sites) &&
    ncol(data$A) == nrow(pars$field) && all(vapply(data[c("M0", "M1", "M2")], function(m)
      identical(dim(m), rep(nrow(pars$field), 2L)), logical(1))) &&
    all(abs(Matrix::rowSums(data$A) - 1) < 1e-8), "Exact 15-km spatial mesh/projection arrays are inconsistent.")
  previous <- rep(-1L, n); gap <- rep(1L, n)
  if (n > 1L) {
    rows <- which(as.character(d$pop_id[-1L]) == as.character(d$pop_id[-n])) + 1L
    previous[rows] <- rows - 2L
    gap[rows] <- as.integer(d$YEAR[rows] - d$YEAR[rows - 1L])
  }
  demo_phylo_need(identical(as.integer(data$previous), previous) && identical(as.integer(data$gap), gap) &&
    all(gap > 0L), "Population AR1 links/calendar gaps changed.")
  popvars <- c("pop_id", "SPECIES", "SITE_ID", "onset_plasticity_z", "offset_plasticity_bio_z", demo_phylo_X_names()[3:8])
  pop <- unique(d[popvars])
  demo_phylo_need(!anyDuplicated(pop$pop_id) && nrow(pop) == sample[["npop"]], "Plasticity/context varies within population.")
  for (key in c("onset", "offset")) {
    total <- pop[[if (key == "onset") "onset_plasticity_z" else "offset_plasticity_bio_z"]]
    within <- pop[[paste0(key, "_plasticity_within_species_z")]]
    between <- pop[[paste0(key, "_plasticity_species_mean_z")]]
    mu <- tapply(total, as.character(pop$SPECIES), mean)
    demo_phylo_need(demo_phylo_same(total, within + between) &&
      demo_phylo_same(between, as.numeric(mu[as.character(pop$SPECIES)])), "Plasticity Mundlak decomposition/weighting changed.")
  }
  ck <- fit$checks; opt <- fit$optimizer
  demo_phylo_need(is.data.frame(ck) && nrow(ck) == 1L && identical(as.character(ck$group), cfg$group) &&
    isTRUE(ck$numerical_checks_ok) && identical(as.integer(opt$convergence), 0L) && isTRUE(ck$pdHess) &&
    is.finite(ck$gradient_covariance_norm) && ck$gradient_covariance_norm < min(.05, cfg$scaled_gradient_threshold) &&
    is.finite(ck$max_abs_gradient) && ck$max_abs_gradient >= 0 && is.finite(opt$objective) &&
    is.finite(p$source_fit_objective) && opt$objective <= p$source_fit_objective + .001 && ck$nobs == n,
    "Independent demographic baseline fails its original numerical/likelihood acceptance gate.")
  outer <- unlist(pars[c("beta", "log_sd", "rho_raw", "log_shape", "log_range", "log_spatial_sd")], use.names = FALSE)
  outer_names <- c(rep("beta", 14L), rep("log_sd", 5L), "rho_raw", "log_shape", "log_range", "log_spatial_sd")
  demo_phylo_need(identical(names(opt$par), outer_names) && demo_phylo_same(as.numeric(opt$par), outer) &&
    identical(names(fit$beta), colnames(X)) && demo_phylo_same(as.numeric(fit$beta), pars$beta),
    "Saved demographic optimizer/parameter/beta mappings disagree.")
  Vfull <- fit$cov_fixed; V <- fit$V
  demo_phylo_need(is.matrix(Vfull) && identical(dim(Vfull), c(23L, 23L)) && all(is.finite(Vfull)) &&
    all(diag(Vfull) > 0) && isTRUE(all.equal(Vfull, t(Vfull), check.attributes = FALSE)) &&
    !inherits(try(chol(Vfull), silent = TRUE), "try-error") &&
    is.matrix(V) && identical(dimnames(V), list(colnames(X), colnames(X))) &&
    demo_phylo_same(V, Vfull[1:14, 1:14, drop = FALSE]), "Saved demographic fixed covariance is inconsistent/non-positive.")
  demo_phylo_need(demo_phylo_same(as.numeric(ck$logLik), -as.numeric(opt$objective)) &&
    ck$parameter_df == length(outer) && demo_phylo_same(as.numeric(ck$AIC), 2 * opt$objective + 2 * length(outer)),
    "Demographic likelihood/AIC summaries do not match the saved optimum.")
  spatial <- as.matrix(data$A %*% pars$field)
  eta <- as.numeric(X %*% fit$beta) + pars$b_species[data$species + 1L, 1L] +
    data$time * pars$b_species[data$species + 1L, 2L] + pars$b_site[data$site + 1L, 1L] +
    data$time * pars$b_site[data$site + 1L, 2L] + pars$ar_obs + spatial[cbind(data$site + 1L, data$year + 1L)]
  demo_phylo_need(demo_phylo_same(spatial, fit$field_at_sites, 1e-7) &&
    demo_phylo_same(eta, fit$eta, 1e-7), "Saved demographic latent arrays do not reproduce the fitted predictor.")
  scales <- p$component_scales
  demo_phylo_need(is.data.frame(scales) && nrow(scales) == 4L && all(c("component", "fitted_variable",
    "SD_for_descriptive_rescaling", "SD_weighting", "fitted_unit") %in% names(scales)) &&
    identical(as.character(scales$fitted_variable), demo_phylo_X_names()[5:8]) &&
    all(is.finite(scales$SD_for_descriptive_rescaling) & scales$SD_for_descriptive_rescaling > 0),
    "Missing or inconsistent four-component demographic scaling table.")
  species_one <- pop[!duplicated(pop$SPECIES), , drop = FALSE]
  expected_sd <- c(stats::sd(pop$onset_plasticity_within_species_z), stats::sd(species_one$onset_plasticity_species_mean_z),
    stats::sd(pop$offset_plasticity_within_species_z), stats::sd(species_one$offset_plasticity_species_mean_z))
  demo_phylo_need(demo_phylo_same(scales$SD_for_descriptive_rescaling, expected_sd), "Demographic component descriptive SDs changed.")
  tree <- ape::read.tree(cfg$tree_file)
  demo_phylo_need(inherits(tree, "phylo") && !inherits(tree, "multiPhylo") && ape::is.rooted(tree) &&
    ape::is.ultrametric(tree) && length(tree$edge.length) == nrow(tree$edge) &&
    all(is.finite(tree$edge.length)) && all(tree$edge.length > 0), "Invalid rooted ultrametric demographic tree.")
  sp <- normalize_species(species); tips <- normalize_species(tree$tip.label)
  demo_phylo_need(!anyDuplicated(sp) && !anyDuplicated(tips) && !anyNA(match(sp, tips)), "Incomplete/ambiguous demographic species-tree match.")
  tip_names <- tree$tip.label[match(sp, tips)]
  edge <- build_phylo_edge_basis(tree, tip_names)
  B <- methods::as(Matrix::Matrix(edge$B, sparse = TRUE), "CsparseMatrix")
  demo_phylo_need(nrow(B) == length(species) && all(is.finite(B@x)) &&
    8 * nrow(B) * ncol(B) / 1024^3 <= cfg$maximum_single_basis_GiB,
    "Demographic tree basis exceeds its configured size/dimension limit.")
  original_data_hash <- digest::digest(data, algo = "sha256")
  mesh_hash <- digest::digest(data[c("A", "M0", "M1", "M2")], algo = "sha256")
  source_hashes <- tools::md5sum(unlist(resolved[c("source_checkpoint", "prepared_file", "run_config_file")], use.names = FALSE))
  source <- list(group = cfg$group, beta = fit$beta, V = V,
    parameters = pars[c("beta", "log_sd", "rho_raw", "log_shape", "log_range", "log_spatial_sd")],
    optimizer = opt, checks = ck,
    stats = data.frame(nobs = n, k = length(outer), logLik = -opt$objective, AIC = ck$AIC),
    component_scales = scales, row_hash = digest::digest(d, algo = "sha256"), data_hash = original_data_hash,
    mesh_hash = mesh_hash, checkpoint = resolved$source_checkpoint, signature = p$signature,
    run_signature = run$run_signature, baseline_version = version,
    scaled_gradient_recomputed = FALSE, numerical_gate = "original saved covariance-scaled gradient norm < 0.05")
  fingerprint <- digest::digest(list(source, source_hashes, tree, B, packages, cfg$signature_file_hashes,
    cfg$scaled_gradient_threshold), algo = "sha256")
  engine_version <- if (exists("VERSION", inherits = TRUE)) get("VERSION", inherits = TRUE) else "demography_phylogeny_source_1.0"
  signature <- digest::digest(list(engine_version, cfg$group, fingerprint), algo = "sha256")
  out <- file.path(cfg$output_root, paste0(cfg$group, "_", substr(signature, 1L, 12L)))
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  data$B <- B; data$use_phylo <- 1L
  pars$log_sd_phylo <- log(max(.1 * exp(pars$log_sd[[2L]]), .001))
  pars$z_phylo <- rep(0, ncol(B))
  result <- list(version = engine_version, group = cfg$group, data = data, parameters = pars,
    frame = d, sites = p$sites, years = p$years, X_names = colnames(X), species = species, B = B,
    A = edge$A, source = source, component_scales = scales, signature = signature, out = out,
    source_file_hashes = source_hashes, source_prepared_file = resolved$prepared_file,
    source_checkpoint = resolved$source_checkpoint, packages = packages)
  demo_phylo_need(demo_phylo_plain(result), "Portable demographic preparation contains a native pointer/environment/function.")
  utils::write.csv(ck, file.path(out, "source_fit_checks.csv"), row.names = FALSE)
  utils::write.csv(scales, file.path(out, "source_component_scales.csv"), row.names = FALSE)
  utils::write.csv(data.frame(SPECIES = species, source_species_index = seq_along(species) - 1L, tip_label = tip_names),
    file.path(out, "species_mapping.csv"), row.names = FALSE)
  utils::write.csv(edge$audit, file.path(out, "tree_increment_basis_checks.csv"), row.names = FALSE)
  if (is.data.frame(edge$clades)) utils::write.csv(edge$clades, file.path(out, "tree_increment_clades.csv"), row.names = FALSE)
  result
}

demo_phylo_prepare <- function(cfg, packages = NULL) {
  if (is.null(cfg$scaled_gradient_threshold)) cfg$scaled_gradient_threshold <- .05
  if (is.null(cfg$maximum_single_basis_GiB)) cfg$maximum_single_basis_GiB <- 2
  cfg$project_root <- normalizePath(cfg$project_root, winslash = "/", mustWork = TRUE)
  if (!file.exists(cfg$tree_file)) cfg$tree_file <- file.path(cfg$project_root, cfg$tree_file)
  cfg$tree_file <- normalizePath(cfg$tree_file, winslash = "/", mustWork = TRUE)
  resolved <- demo_phylo_resolve_source(cfg$project_root, cfg$group, cfg$source_directory)
  inputs <- unlist(resolved[c("source_checkpoint", "prepared_file", "run_config_file")], use.names = FALSE)
  optional <- unlist(resolved[c("fit_checks_file", "result_file")], use.names = FALSE)
  inputs <- c(inputs, optional[file.exists(optional)], cfg$tree_file)
  engine_version <- if (exists("VERSION", inherits = TRUE)) get("VERSION", inherits = TRUE) else "demography_phylogeny_source_1.0"
  manifest <- list(version = engine_version, group = cfg$group, files = tools::md5sum(inputs),
    packages = packages, signature_file_hashes = cfg$signature_file_hashes,
    scaled_gradient_threshold = cfg$scaled_gradient_threshold, expected_samples = cfg$expected_samples)
  key <- digest::digest(manifest, algo = "sha256")
  cache <- file.path(cfg$output_root, "source_cache", paste0("prepared_", key, ".rds"))
  reused <- FALSE
  if (file.exists(cache)) {
    x <- readRDS(cache)
    demo_phylo_need(is.list(x) && identical(x$manifest, manifest) && demo_phylo_plain(x$prepared) &&
      identical(x$prepared$group, cfg$group) &&
      digest::digest(x$prepared, algo = "sha256") == x$content_hash,
      "Demographic portable source cache failed its identity/content checks.")
    prepared <- x$prepared
    demo_phylo_need(8 * nrow(prepared$B) * ncol(prepared$B) / 1024^3 <= cfg$maximum_single_basis_GiB,
      "Cached demographic basis exceeds its configured size limit.")
    reused <- TRUE
  } else {
    prepared <- demo_phylo_prepare_uncached(cfg, packages, resolved)
    demo_phylo_save(list(manifest = manifest, prepared = prepared,
      content_hash = digest::digest(prepared, algo = "sha256")), cache)
  }
  dir.create(prepared$out, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(data.frame(cache_file = cache, reused = reused, native_source_replay_run = FALSE,
    scaled_gradient_recomputed = FALSE), file.path(prepared$out, "source_cache_used.csv"), row.names = FALSE)
  prepared
}


demo_phylo_need <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
}

demo_phylo_numeric_equal <- function(a, b, tolerance = 1e-7) {
  is.numeric(a) && is.numeric(b) && length(a) == length(b) &&
    all(is.finite(a)) && all(is.finite(b)) &&
    (!length(a) || max(abs(a - b)) <= tolerance)
}

demo_phylo_covariance <- function(V, n, labels = NULL) {
  demo_phylo_need(is.matrix(V) && is.numeric(V) && identical(dim(V), c(n, n)) &&
    all(is.finite(V)), "Fixed covariance has invalid dimensions or values.")
  if (!is.null(labels) && !is.null(dimnames(V))) {
    demo_phylo_need(identical(rownames(V), labels) && identical(colnames(V), labels),
      "Fixed covariance parameter ordering differs from its estimates.")
  }
  magnitude <- max(abs(V))
  demo_phylo_need(magnitude > 0 && max(abs(V - t(V))) / magnitude <= 1e-8,
    "Fixed covariance is not symmetric within rounding tolerance.")
  V <- (V + t(V)) / 2
  demo_phylo_need(all(diag(V) > 0), "Fixed covariance has non-positive diagonal entries.")
  R <- tryCatch(chol(V), error = function(e) NULL)
  demo_phylo_need(!is.null(R), "Fixed covariance is not positive definite.")
  list(V = V, R = R)
}

demo_phylo_table <- function(beta, V) {
  se <- sqrt(diag(V))
  data.frame(term = names(beta), estimate = unname(beta), SE = unname(se),
    lower_95 = unname(beta - 1.96 * se), upper_95 = unname(beta + 1.96 * se),
    p_Wald = unname(2 * stats::pnorm(-abs(beta / se))), stringsAsFactors = FALSE)
}

demo_phylo_limits <- function(cfg = list()) {
  scaled <- if (!is.null(cfg$scaled_gradient_threshold)) cfg$scaled_gradient_threshold else
    if (exists("scaled_gradient_threshold", inherits = TRUE))
    get("scaled_gradient_threshold", inherits = TRUE) else 0.05
  raw <- if (!is.null(cfg$raw_gradient_threshold)) cfg$raw_gradient_threshold else
    if (exists("raw_gradient_threshold", inherits = TRUE))
    get("raw_gradient_threshold", inherits = TRUE) else 0.001
  demo_phylo_need(length(scaled) == 1L && is.finite(scaled) && scaled > 0 &&
    length(raw) == 1L && is.finite(raw) && raw > 0, "Invalid numerical thresholds.")
  list(scaled = scaled, raw = raw)
}

demo_phylo_contrast <- function(beta, V, within, between) {
  w <- stats::setNames(rep(0, length(beta)), names(beta))
  w[within] <- 1
  w[between] <- -1
  variance <- sum(as.vector(chol(V) %*% w)^2)
  demo_phylo_need(is.finite(variance) && variance >= 0, "Invalid within-minus-between contrast variance.")
  estimate <- sum(w * beta)
  se <- sqrt(variance)
  data.frame(estimate = estimate, SE = se, lower_95 = estimate - 1.96 * se,
    upper_95 = estimate + 1.96 * se,
    p_Wald = if (se > 0) 2 * stats::pnorm(-abs(estimate / se)) else NA_real_)
}

demo_phylo_wide_pair <- function(old, new) {
  names(old) <- paste0(names(old), "_independent")
  names(new) <- paste0(names(new), "_phylogenetic")
  cbind(old, new)
}

demo_phylo_expected_names <- function() c(
  "(Intercept)", "year_decade", "latitude_within_species_10deg",
  "longitude_within_species_10deg", "onset_plasticity_within_species_z",
  "onset_plasticity_species_mean_z", "offset_plasticity_within_species_z",
  "offset_plasticity_species_mean_z", "year_decade:latitude_within_species_10deg",
  "year_decade:longitude_within_species_10deg", "year_decade:onset_plasticity_within_species_z",
  "year_decade:onset_plasticity_species_mean_z", "year_decade:offset_plasticity_within_species_z",
  "year_decade:offset_plasticity_species_mean_z")

demo_phylo_state <- function(fit, p = NULL) {
  theta <- fit$optimizer$par
  raw <- names(theta)
  labels <- fit$outer_names
  demo_phylo_need(is.numeric(theta) && length(theta) > 0L && !is.null(raw) &&
    !anyNA(raw) && all(nzchar(raw)) && all(is.finite(theta)), "Invalid saved outer parameters.")
  demo_phylo_need(is.character(labels) && length(labels) == length(theta) &&
    !anyNA(labels) && all(nzchar(labels)) && !anyDuplicated(labels) &&
    identical(fit$outer_raw_names, raw), "Invalid or inconsistent saved outer parameter labels.")
  demo_phylo_need(identical(names(fit$outer_parameters), raw) &&
    demo_phylo_numeric_equal(theta, fit$outer_parameters),
    "Optimizer parameters and saved uncertainty-report parameters differ.")
  blocks <- sub("\\[.*$", "", raw)
  expected_blocks <- c("beta", "log_sd", "rho_raw", "log_shape", "log_range",
    "log_spatial_sd", "log_sd_phylo")
  demo_phylo_need(setequal(unique(blocks), expected_blocks) && is.list(fit$parameters),
    "Unexpected Gamma/log demographic outer parameter blocks.")
  for (block in unique(blocks)) demo_phylo_need(
    demo_phylo_numeric_equal(as.numeric(theta[blocks == block]), as.vector(fit$parameters[[block]])),
    paste("Saved parameter block differs from optimizer:", block))
  bi <- which(blocks == "beta")
  sd_ids <- which(blocks == "log_sd")
  phylo <- which(blocks == "log_sd_phylo")
  demo_phylo_need(length(bi) == 14L && length(sd_ids) == 5L && length(phylo) == 1L,
    "Unexpected fixed-effect or variance-parameter dimensions.")
  beta_names <- demo_phylo_expected_names()
  demo_phylo_need(is.numeric(fit$beta) && identical(names(fit$beta), beta_names) &&
    demo_phylo_numeric_equal(unname(fit$beta), unname(theta[bi])),
    "Saved fixed coefficients do not match the canonical 14-column design and optimizer.")
  if (!is.null(p)) {
    design_names <- if (!is.null(p$X_names)) p$X_names else colnames(p[["X", exact = TRUE]])
    demo_phylo_need(identical(design_names, beta_names) && identical(fit$group, p$group),
      "Prepared design or group differs from saved fit.")
    if (!is.null(fit$signature) && !is.null(p$signature)) demo_phylo_need(
      identical(fit$signature, p$signature), "Saved fit signature differs from prepared input.")
  }
  demo_phylo_need(fit$group %in% c("univoltine", "multivoltine") && length(fit$group) == 1L,
    "Invalid demographic group.")
  if (!is.null(dimnames(fit$cov_fixed))) demo_phylo_need(
    (identical(rownames(fit$cov_fixed), raw) && identical(colnames(fit$cov_fixed), raw)) ||
    (identical(rownames(fit$cov_fixed), labels) && identical(colnames(fit$cov_fixed), labels)),
    "Full covariance parameter labels/order are inconsistent.")
  covariance <- demo_phylo_covariance(fit$cov_fixed, length(theta))
  V <- covariance$V
  g <- if (!is.null(fit[["gradient", exact = TRUE]])) fit[["gradient", exact = TRUE]] else
    fit[["gradients", exact = TRUE]]
  demo_phylo_need(is.numeric(g) && length(g) == length(theta) && all(is.finite(g)),
    "Stored full outer gradient is missing, non-finite, or misaligned.")
  if (!is.null(names(g))) demo_phylo_need(identical(names(g), raw) || identical(names(g), labels),
    "Stored gradient labels/order are inconsistent.")
  if (!is.null(fit$gradient_report)) demo_phylo_need(demo_phylo_numeric_equal(g, fit$gradient_report),
    "Stored and uncertainty-report gradients differ.")
  Vf <- demo_phylo_covariance(fit$V, 14L, beta_names)$V
  demo_phylo_need(demo_phylo_numeric_equal(as.vector(Vf), as.vector(V[bi, bi, drop = FALSE])),
    "Saved fixed covariance differs from the full covariance beta block.")
  beta <- stats::setNames(unname(theta[bi]), beta_names)
  objective <- fit$optimizer$objective
  demo_phylo_need(is.numeric(objective) && length(objective) == 1L && is.finite(objective),
    "Saved marginal likelihood is invalid.")
  norm <- sqrt(sum(as.vector(covariance$R %*% g)^2))
  demo_phylo_need(is.finite(norm), "Covariance-scaled gradient is not finite.")
  se <- sqrt(diag(V))
  list(theta = theta, raw = raw, labels = labels, beta = beta, V = V, Vf = Vf,
    g = g, iid = sd_ids[2L], phylo = phylo, objective = objective,
    fixed = demo_phylo_table(beta, Vf), raw_gradient = max(abs(g)),
    scaled_gradient = norm, score_SE = max(abs(g * se)),
    newton_SE = max(abs(as.vector(V %*% g) / se)),
    theta_report_difference = max(abs(theta - fit$outer_parameters)))
}

demo_phylo_checks <- function(fit, p = NULL, cfg = list()) {
  limits <- demo_phylo_limits(cfg)
  failure <- ""
  state <- tryCatch(demo_phylo_state(fit, p), error = function(e) {
    failure <<- conditionMessage(e)
    NULL
  })
  code <- fit$optimizer$convergence
  if (length(code) != 1L || !is.numeric(code) || !is.finite(code)) code <- NA_integer_
  pd <- if (!is.null(fit$pdHess)) isTRUE(fit$pdHess) else isTRUE(fit$checks$pdHess)
  good <- !is.null(state) && isTRUE(code == 0L) && pd && state$scaled_gradient < limits$scaled
  if (!nzchar(failure) && !good) failure <- paste(c(
    if (!isTRUE(code == 0L)) "Optimizer convergence code is not zero.",
    if (!pd) "Saved uncertainty report does not confirm a positive Hessian.",
    if (!is.null(state) && state$scaled_gradient >= limits$scaled)
      "Covariance-scaled gradient exceeds the prescribed threshold."), collapse = " ")
  value <- function(name) if (is.null(state)) NA_real_ else state[[name]]
  nobs <- fit$nobs
  if (is.null(nobs) && !is.null(p$data$y)) nobs <- length(p$data$y)
  if (is.null(nobs) && !is.null(p$frame)) nobs <- nrow(p$frame)
  if (is.null(nobs) && !is.null(p[["X", exact = TRUE]])) nobs <- nrow(p[["X", exact = TRUE]])
  if (length(nobs) != 1L || !is.finite(nobs)) nobs <- NA_real_
  data.frame(group = if (length(fit$group) == 1L) fit$group else NA_character_,
    convergence_code = code, pdHess = pd, max_abs_gradient = value("raw_gradient"),
    raw_gradient_threshold = limits$raw,
    raw_gradient_advisory_passed = !is.null(state) && state$raw_gradient < limits$raw,
    gradient_covariance_norm = value("scaled_gradient"), scaled_gradient_threshold = limits$scaled,
    max_abs_score_times_SE = value("score_SE"), max_abs_approx_newton_shift_SE = value("newton_SE"),
    theta_report_max_difference = value("theta_report_difference"),
    finite_fixed_SE = if (is.null(state)) NA else TRUE,
    covariance_positive_definite = if (is.null(state)) NA else TRUE,
    nobs = nobs, parameter_df = if (is.null(state)) NA_integer_ else length(state$theta),
    logLik = if (is.null(state)) NA_real_ else -state$objective,
    AIC = if (is.null(state)) NA_real_ else 2 * state$objective + 2 * length(state$theta),
    numerical_checks_ok = good, fit_ok = good, failure_reason = failure,
    native_evaluations_run = FALSE,
    scope = "Saved Gamma/log ML fit; local covariance-scaled numerical checks, not global optimum or model adequacy proof",
    stringsAsFactors = FALSE)
}

demo_phylo_fixed_table <- function(fit) demo_phylo_state(fit)$fixed

demo_phylo_source <- function(p) {
  source <- p$source
  beta <- source$beta
  demo_phylo_need(is.numeric(beta) && identical(names(beta), demo_phylo_expected_names()) &&
    all(is.finite(beta)) && identical(source$group, p$group), "Invalid independent demographic fixed estimates/group.")
  V <- demo_phylo_covariance(source$V, 14L, names(beta))$V
  demo_phylo_need(demo_phylo_numeric_equal(unname(beta), as.vector(source$parameters$beta)),
    "Independent source parameter and fixed coefficient arrays disagree.")
  checks <- source$checks
  flags <- intersect(c("fit_ok", "numerical_checks_ok"), names(checks))
  limits <- demo_phylo_limits()
  demo_phylo_need(is.data.frame(checks) && nrow(checks) == 1L && length(flags) > 0L &&
    all(vapply(flags, function(nm) isTRUE(checks[[nm]]), logical(1))) &&
    isTRUE(checks$convergence_code == 0L) && isTRUE(checks$pdHess) &&
    length(checks$gradient_covariance_norm) == 1L && is.finite(checks$gradient_covariance_norm) &&
    checks$gradient_covariance_norm < limits$scaled,
    "Independent source does not pass its stored original numerical checks.")
  stats <- source$stats
  pick <- function(names) {
    found <- names[names %in% names(stats)]
    demo_phylo_need(length(found) >= 1L, paste("Missing source likelihood metadata:", names[1L]))
    answer <- stats[[found[1L]]]
    demo_phylo_need(is.numeric(answer) && length(answer) == 1L && is.finite(answer),
      "Invalid source likelihood metadata.")
    answer
  }
  nobs <- pick(c("nobs", "n_observations"))
  k <- pick(c("k", "n_parameters", "parameter_df"))
  ll <- pick("logLik")
  aic <- pick("AIC")
  demo_phylo_need(nobs > 0 && nobs == as.integer(nobs) && k > 14 && k == as.integer(k) &&
    abs(aic - (-2 * ll + 2 * k)) < 1e-5, "Inconsistent independent likelihood/AIC metadata.")
  list(beta = beta, V = V, fixed = demo_phylo_table(beta, V), nobs = nobs,
    k = k, logLik = ll, AIC = aic,
    baseline_QC_scope = "Original saved covariance-scaled gradient check; baseline full gradient was not stored and was not recomputed")
}

demo_phylo_joint_wald <- function(beta, V, C) {
  estimate <- as.vector(C %*% beta)
  covariance <- C %*% V %*% t(C)
  R <- tryCatch(chol(covariance), error = function(e) NULL)
  demo_phylo_need(!is.null(R), "Joint contrast covariance is not positive definite.")
  statistic <- sum(as.vector(backsolve(R, estimate, transpose = TRUE))^2)
  data.frame(Wald_chisq = statistic, df = nrow(C),
    p_Wald = stats::pchisq(statistic, df = nrow(C), lower.tail = FALSE))
}

demo_phylo_payload <- function(fit, p, cfg) {
  state <- demo_phylo_state(fit, p)
  source <- demo_phylo_source(p)
  nobs <- if (!is.null(fit$nobs)) fit$nobs else fit$checks$nobs
  if (is.null(nobs) && !is.null(p[["X", exact = TRUE]])) nobs <- nrow(p[["X", exact = TRUE]])
  demo_phylo_need(length(nobs) == 1L && is.finite(nobs) && nobs == source$nobs &&
    length(state$theta) == source$k + 1L,
    "Compared demographic models must retain observations and add exactly one outer parameter.")
  if (!is.null(p[["X", exact = TRUE]])) demo_phylo_need(nrow(p[["X", exact = TRUE]]) == nobs,
    "Prepared design row count differs from saved model.")
  tolerance <- if (is.null(cfg$likelihood_tolerance)) 0.001 else cfg$likelihood_tolerance
  demo_phylo_need(length(tolerance) == 1L && is.finite(tolerance) && tolerance >= 0 &&
    state$objective <= -source$logLik + tolerance,
    "Phylogenetic marginal likelihood is worse than its nested independent source.")
  fixed <- data.frame(group = p$group, term = names(state$beta),
    demo_phylo_wide_pair(source$fixed[-1L], state$fixed[-1L]),
    change_in_independent_SE = (state$beta - source$beta) / source$fixed$SE,
    SE_ratio = state$fixed$SE / source$fixed$SE, row.names = NULL)
  variables <- c("onset_plasticity_within_species_z", "onset_plasticity_species_mean_z",
    "offset_plasticity_within_species_z", "offset_plasticity_species_mean_z")
  terms <- paste0("year_decade:", variables)
  scaling <- p$component_scales
  demo_phylo_need(is.data.frame(scaling) && nrow(scaling) == 4L &&
    all(c("fitted_variable", "SD_for_descriptive_rescaling", "SD_weighting") %in% names(scaling)) &&
    !anyDuplicated(scaling$fitted_variable) && setequal(scaling$fitted_variable, variables),
    "Invalid four-row plasticity component scaling metadata.")
  scaling <- scaling[match(variables, scaling$fitted_variable), , drop = FALSE]
  demo_phylo_need(all(is.finite(scaling$SD_for_descriptive_rescaling) &
    scaling$SD_for_descriptive_rescaling > 0) &&
    identical(as.character(scaling$SD_weighting), rep(c("equal population", "equal species"), 2L)),
    "Invalid component SDs or component-weighting metadata.")
  if ("group" %in% names(scaling)) demo_phylo_need(all(scaling$group == p$group),
    "Component scaling belongs to another demographic group.")
  rows <- lapply(seq_along(terms), function(i) {
    term <- terms[i]
    old <- source$fixed[match(term, source$fixed$term), -1L, drop = FALSE]
    new <- state$fixed[match(term, state$fixed$term), -1L, drop = FALSE]
    sd <- scaling$SD_for_descriptive_rescaling[i]
    description <- if (i <= 2L) "greater onset advancement" else if (p$group == "univoltine")
      "greater offset advancement" else "greater offset delay"
    ratio <- function(z) data.frame(ratio_of_decadal_growth_factors = exp(z$estimate),
      ratio_lower_95 = exp(z$lower_95), ratio_upper_95 = exp(z$upper_95),
      estimate_per_component_SD = z$estimate * sd,
      lower_95_per_component_SD = z$lower_95 * sd, upper_95_per_component_SD = z$upper_95 * sd)
    data.frame(group = p$group, plasticity = if (i <= 2L) "onset" else "offset",
      effect_scale = if (i %% 2L == 1L) "within species" else "between species",
      term = term, predictor_direction = description,
      component_SD = sd, component_SD_weighting = scaling$SD_weighting[i],
      demo_phylo_wide_pair(cbind(old, ratio(old)), cbind(new, ratio(new))),
      change_in_independent_SE = (new$estimate - old$estimate) / old$SE,
      SE_ratio = new$SE / old$SE,
      interpretation = "Positive coefficient: larger decadal abundance-change factor with higher already-oriented plasticity; not a causal effect",
      row.names = NULL)
  })
  differences <- lapply(c(1L, 3L), function(i) {
    old <- demo_phylo_contrast(source$beta, source$V, terms[i], terms[i + 1L])
    new <- demo_phylo_contrast(state$beta, state$Vf, terms[i], terms[i + 1L])
    data.frame(group = p$group, plasticity = if (i == 1L) "onset" else "offset",
      contrast = "within minus between", demo_phylo_wide_pair(old, new), row.names = NULL)
  })
  C4 <- matrix(0, 4L, 14L, dimnames = list(NULL, names(state$beta)))
  C4[cbind(seq_len(4L), match(terms, colnames(C4)))] <- 1
  contrasts <- list(within_species_effects_zero = C4[c(1L, 3L), , drop = FALSE],
    between_species_effects_zero = C4[c(2L, 4L), , drop = FALSE],
    within_equals_between = rbind(C4[1L, ] - C4[2L, ], C4[3L, ] - C4[4L, ]),
    all_four_plasticity_effects_zero = C4)
  joint <- lapply(names(contrasts), function(name) data.frame(group = p$group, test = name,
    demo_phylo_wide_pair(demo_phylo_joint_wald(source$beta, source$V, contrasts[[name]]),
      demo_phylo_joint_wald(state$beta, state$Vf, contrasts[[name]])),
    interpretation = "Approximate fixed-effect Wald test within each model; not a test between models or of the phylogenetic variance"))
  ids <- c(state$iid, state$phylo)
  logsd <- unname(state$theta[ids])
  se <- sqrt(diag(state$V)[ids])
  demo_phylo_need(all(is.finite(exp(2 * logsd))), "Temporal slope variance estimate overflows.")
  old_logsd <- p$source$parameters$log_sd
  demo_phylo_need(is.numeric(old_logsd) && length(old_logsd) == 5L && all(is.finite(old_logsd)) &&
    is.finite(exp(2 * old_logsd[2L])), "Independent temporal species slope SD is invalid.")
  variance <- data.frame(group = p$group, model = c("independent", rep("independent_plus_phylogeny", 2L)),
    component = c("iid_species_temporal_slope", "iid_species_temporal_slope", "phylogenetic_species_temporal_slope"),
    log_SD = c(old_logsd[2L], logsd), log_SD_SE = c(NA_real_, se),
    SD = exp(c(old_logsd[2L], logsd)),
    SD_lower_95 = c(NA_real_, exp(logsd - 1.96 * se)),
    SD_upper_95 = c(NA_real_, exp(logsd + 1.96 * se)),
    variance = exp(2 * c(old_logsd[2L], logsd)),
    interval_method = c("Baseline point estimate; full variance covariance not supplied",
      rep("Approximate normal-Wald log-SD interval; variance-boundary caveat", 2L)))
  w <- rep(0, length(state$theta)); w[state$phylo] <- 2; w[state$iid] <- -2
  logit <- 2 * (logsd[2L] - logsd[1L])
  logit_se <- sqrt(sum(as.vector(chol(state$V) %*% w)^2))
  fraction <- stats::plogis(logit)
  fraction_table <- data.frame(group = p$group, phylogenetic_fraction = fraction,
    lower_95 = stats::plogis(logit - 1.96 * logit_se), upper_95 = stats::plogis(logit + 1.96 * logit_se),
    near_variance_boundary_descriptive = fraction < 0.01 || fraction > 0.99,
    interpretation = "Conditional between-species temporal trend variance fraction, not total abundance variance",
    interval_method = "Joint outer-covariance delta method on logit fraction; approximate near variance boundaries")
  ll <- -state$objective; k <- length(state$theta)
  comparison <- data.frame(group = p$group, model = c("independent", "independent_plus_phylogeny"),
    nobs = c(source$nobs, nobs), k = c(source$k, k), logLik = c(source$logLik, ll),
    AIC = c(source$AIC, -2 * ll + 2 * k), delta_AIC_from_independent = c(0, -2 * ll + 2 * k - source$AIC),
    LR_vs_independent_descriptive = c(0, 2 * (ll - source$logLik)), boundary_LRT_p_computed = FALSE,
    likelihood = "Gamma/log marginal ML with the same Laplace approximation, data, predictors and spatial/temporal structure",
    baseline_QC_scope = source$baseline_QC_scope)
  list(fixed_effect_comparison.csv = fixed, focal_effect_comparison.csv = do.call(rbind, rows),
    within_minus_between_comparison.csv = do.call(rbind, differences),
    joint_Wald_comparison.csv = do.call(rbind, joint), variance_components.csv = variance,
    phylogenetic_fraction.csv = fraction_table, model_comparison.csv = comparison)
}

demo_phylo_export_results <- function(fit, p, cfg = list()) {
  out <- p$out
  demo_phylo_need(is.character(out) && length(out) == 1L && nzchar(out), "Missing demographic output directory.")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  checks <- demo_phylo_checks(fit, p, cfg)
  failure <- checks$failure_reason
  payload <- if (isTRUE(checks$fit_ok)) tryCatch(demo_phylo_payload(fit, p, cfg),
    error = function(e) {failure <<- conditionMessage(e); NULL}) else NULL
  good <- !is.null(payload)
  paths <- c("fixed_effect_comparison.csv", "focal_effect_comparison.csv",
    "within_minus_between_comparison.csv", "joint_Wald_comparison.csv", "variance_components.csv",
    "phylogenetic_fraction.csv", "model_comparison.csv")
  writer <- function(x, path) utils::write.csv(x, path, row.names = FALSE, na = "")
  if (good) {
    writer(data.frame(group = p$group, status = "EXPORTING", inference_exports_written = FALSE),
      file.path(out, "export_status.csv"))
    tryCatch(for (name in names(payload)) writer(payload[[name]], file.path(out, name)),
      error = function(e) {good <<- FALSE; failure <<- paste("Comparison CSV writing failed:", conditionMessage(e))})
  }
  if (!good) {
    old <- file.path(out, paths); old <- old[file.exists(old)]
    if (length(old)) demo_phylo_need(all(vapply(old, function(x) unlink(x) == 0L, logical(1))),
      "Cannot remove stale comparison exports after a failed check.")
  }
  checks$comparison_checks_ok <- good
  checks$numerical_checks_ok <- checks$fit_ok <- isTRUE(checks$fit_ok) && good
  checks$failure_reason <- failure
  writer(checks, file.path(out, "fit_checks.csv"))
  writer(data.frame(group = p$group, status = if (good) "EXPORTED" else "NUMERICAL_OR_INPUT_REVIEW_REQUIRED",
    inference_exports_written = good, native_evaluations_run = FALSE, reason = failure),
    file.path(out, "export_status.csv"))
  writeLines(c(if (good) "Demographic comparison exports passed saved-state numerical and comparability checks."
      else "Demographic comparison exports suppressed; retain the fitted checkpoint for review.",
    paste("Failure reason:", if (nzchar(failure)) failure else "none"),
    paste("Full covariance-scaled gradient:", checks$gradient_covariance_norm,
      "threshold:", checks$scaled_gradient_threshold),
    paste("Raw gradient advisory passed:", checks$raw_gradient_advisory_passed,
      "raw maximum:", checks$max_abs_gradient),
    "The baseline full gradient was not stored by 10f; its original saved numerical checks are retained without native reevaluation.",
    "Plasticity predictors are held fixed. This isolates phylogenetic covariance of demographic temporal slopes.",
    "Onset and offset predictors are already biologically oriented; demographic coefficients are not sign-flipped.",
    "No native likelihood, gradient, Hessian or optimizer callback was evaluated by this export.",
    "Wald/delta intervals and fixed-effect tests are approximate; no global optimum or causal/model-adequacy claim is established.",
    "Coefficient changes divided by original SE are descriptive, not a test between dependent fits.",
    "No automatic variance-boundary likelihood-ratio p-value is computed."),
    file.path(out, "inference_status.txt"))
  invisible(good)
}



  VERSION <- "mundlak_demography_phylogeny_15km_1.0.0"
  ROOT <- "E:/phenoIMPACT project/code/phenoIMPACT"
  GROUPS <- c("multivoltine", "univoltine")

  need <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
  normalize_species <- function(x) gsub("[[:space:]_]+", "_", tolower(trimws(x)))

  save_atomic <- function(x, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    temporary <- tempfile(".writing_", tmpdir = dirname(path), fileext = ".rds")
    on.exit(unlink(temporary), add = TRUE)
    saveRDS(x, temporary, compress = FALSE)
    backup <- paste0(path, ".previous")
    if (file.exists(backup)) unlink(backup)
    previous <- file.exists(path)
    if (previous) need(file.rename(path, backup), paste("Cannot back up:", path))
    if (!file.rename(temporary, path)) {
      if (previous) file.rename(backup, path)
      stop("Cannot write checkpoint: ", path, call. = FALSE)
    }
    if (file.exists(backup)) unlink(backup)
    invisible(path)
  }

  write_csv <- function(x, path) utils::write.csv(x, path, row.names = FALSE, na = "")
  phase <- function(cfg, stage, detail = "", terminal = FALSE, result_directory = NA_character_) {
    value <- list(stage = stage, group = cfg$group, detail = detail,
      time = Sys.time(), terminal = terminal, result_directory = result_directory,
      started_at = cfg$started_at)
    save_atomic(value, cfg$progress_file)
    message(format(value$time, "%Y-%m-%d %H:%M:%S"), " | ", stage, " | ", detail)
    invisible(value)
  }

  memory_snapshot <- function(directory, stage) {
    tryCatch({
      value <- list(stage = stage, time = Sys.time(),
        system = ps::ps_system_memory(), worker = ps::ps_memory_info(ps::ps_handle()))
      cat(paste(capture.output(dput(value)), collapse = "\n"), "\n\n",
        file = file.path(directory, "memory_checkpoints.txt"), append = TRUE)
    }, error = function(e) message("Memory snapshot unavailable: ", conditionMessage(e)))
  }

  absolute_path <- function(path, root) {
    if (!grepl("^([A-Za-z]:[/\\\\]|/|\\\\\\\\)", path)) path <- file.path(root, path)
    normalizePath(path, winslash = "/", mustWork = TRUE)
  }

  dependencies <- function() {
    names <- c("TMB", "Matrix", "ape", "digest", "callr", "ps")
    missing <- names[!vapply(names, requireNamespace, logical(1), quietly = TRUE)]
    need(!length(missing), paste("Required installed packages:", paste(missing, collapse = ", ")))
    data.frame(package = names, version = vapply(names, function(x)
      as.character(utils::packageVersion(x)), character(1)), stringsAsFactors = FALSE)
  }

  lock_matches <- function(path, token) {
    x <- tryCatch(readRDS(file.path(path, "owner.rds")), error = function(e) NULL)
    !is.null(x) && identical(x$token, token)
  }

  acquire_lock <- function(path, token) {
    if (dir.exists(path)) {
      owner <- tryCatch(readRDS(file.path(path, "owner.rds")), error = function(e) NULL)
      need(!is.null(owner), paste("Inspect existing worker lock:", path))
      alive <- tryCatch({
        if (length(owner$pid) != 1L || !is.numeric(owner$pid) || !is.finite(owner$pid) || owner$pid <= 0) {
          NA
        } else if (!owner$pid %in% ps::ps_pids()) FALSE else if (
          length(owner$created) != 1L || !is.numeric(owner$created) || !is.finite(owner$created)) NA else {
          h <- ps::ps_handle(owner$pid)
          isTRUE(abs(as.numeric(ps::ps_create_time(h)) - owner$created) < .01) &&
            ps::ps_is_running(h) && !ps::ps_status(h) %in% c("zombie", "dead")
        }
      }, error = function(e) NA)
      need(identical(alive, FALSE), paste("A demographic phylogeny worker may still be active:", path))
      unlink(path, recursive = TRUE)
    }
    need(dir.create(path), paste("Cannot acquire lock:", path))
    save_atomic(list(pid = Sys.getpid(), created = as.numeric(ps::ps_create_time(ps::ps_handle())),
      token = token), file.path(path, "owner.rds"))
  }

  release_lock <- function(cfg) {
    if (lock_matches(cfg$launch_lock, cfg$lock_token)) unlink(cfg$launch_lock, recursive = TRUE)
  }

  outer_labels <- function(theta, terms) {
    raw <- names(theta)
    answer <- raw
    answer[raw == "beta"] <- paste0("fixed:", terms)
    answer[raw == "log_sd"] <- paste0("log_SD:", c("species_intercept", "species_time",
      "site_intercept", "site_time", "population_AR1"))
    make.unique(answer)
  }

  worker <- function(cfg) {
    on.exit(release_lock(cfg), add = TRUE)
    tryCatch({
      packages <- dependencies()
      writeLines(capture.output(utils::sessionInfo()), file.path(cfg$monitor_directory, "sessionInfo.txt"))
      need(ps::ps_disk_usage(cfg$output_root)$available / 1024^3 >= cfg$minimum_disk_GiB,
           "Insufficient free disk space for checkpoints.")
      phase(cfg, "NATIVE_PREFLIGHT", "Compile/cache the extended likelihood; test tiny cases before reading the source")
      engine <- demo_phylo_compile(cfg$compiled_directory)
      cfg$DLL <- engine$DLL
      demo_phylo_native_preflight(cfg$compiled_directory)
      phase(cfg, "PREPARING", "Read or reuse the completed independent demographic Mundlak fit")
      p <- demo_phylo_prepare(cfg, packages)
      cfg$result_directory <- p$out
      phase(cfg, "PREPARED", "Portable data, original 15-km mesh and tree covariance verified", result_directory = p$out)
      file.copy(file.path(cfg$compiled_directory, "native_preflight.csv"),
        file.path(p$out, "native_preflight.csv"), overwrite = TRUE)
      if (cfg$check_only) {
        phase(cfg, "CHECKS_PASSED", "No full model constructed or optimized", TRUE, p$out)
        return(list(output_dir = p$out, check_only = TRUE))
      }
      final_file <- file.path(p$out, "fit_phylogenetic.rds")
      if (file.exists(final_file) && !cfg$refine) {
        phase(cfg, "REUSING", "Validate saved arrays and regenerate exports; no native model construction", result_directory = p$out)
        fit <- readRDS(final_file)
        need(identical(fit$signature, p$signature), "Saved demographic phylogenetic signature mismatch.")
        okay <- demo_phylo_export_results(fit, p, cfg)
        checks <- demo_phylo_checks(fit, p, cfg)
        phase(cfg, "CACHED", if (okay) "Saved fit and inference reused" else "Saved fit retained; inspect fit_checks and inference_status", TRUE, p$out)
        return(list(output_dir = p$out, cached = TRUE, fit_ok = isTRUE(checks$fit_ok), inference_ready = isTRUE(okay)))
      }
      checkpoint_file <- file.path(p$out, "portable_optimizer_state.rds")
      checkpoint <- if (file.exists(checkpoint_file)) readRDS(checkpoint_file) else NULL
      parameters <- p$parameters
      if (!is.null(checkpoint)) {
        need(identical(checkpoint$signature, p$signature) && is.list(checkpoint$parameters),
             "Portable optimizer checkpoint mismatch.")
        parameters <- checkpoint$parameters
      }
      phase(cfg, "BUILDING_NATIVE_MODEL", "Fresh native state; retained random modes are warm starts", result_directory = p$out)
      memory_snapshot(p$out, "before_native_construction")
      obj <- demo_phylo_make_objective(p$data, parameters, phylo = TRUE, DLL = cfg$DLL)
      phase(cfg, "BUILDING_RANDOM_HESSIAN", "This native construction cannot be resumed mid-call", result_directory = p$out)
      hessian_function <- obj$env$spHess
      need(is.function(hessian_function), "Sparse random-effect Hessian unavailable.")
      memory_snapshot(p$out, "after_native_construction")
      theta <- if (is.null(checkpoint)) obj$par else checkpoint$theta
      need(identical(names(theta), names(obj$par)) && length(theta) == length(obj$par) &&
        all(is.finite(theta)), "Restart outer parameter mapping changed.")
      phase(cfg, "INITIAL_EVALUATION", "Marginal likelihood and gradient at the warm start", result_directory = p$out)
      initial <- obj$fn(theta)
      initial_gradient <- obj$gr(theta)
      need(is.finite(initial) && all(is.finite(initial_gradient)), "Initial marginal evaluation is not finite.")
      if (!is.null(checkpoint)) need(abs(initial - checkpoint$nll) < .01,
        "Cannot reproduce the saved marginal objective at the checkpoint.")
      tracker <- new.env(parent = emptyenv())
      tracker$theta <- theta
      tracker$nll <- initial
      tracker$parameters <- obj$env$parList(par = obj$env$last.par)
      tracker$last_write <- -Inf
      tracker$evaluations <- 0L
      persist <- function(stage = "FINITE_OPTIMIZER_EVALUATION", optimizer = NULL) {
        save_atomic(list(version = VERSION, signature = p$signature, stage = stage,
          theta = tracker$theta, nll = tracker$nll, parameters = tracker$parameters,
          optimizer = optimizer, evaluations = tracker$evaluations, time = Sys.time()), checkpoint_file)
        tracker$last_write <- as.numeric(Sys.time())
      }
      optimized <- !is.null(checkpoint) && !cfg$refine &&
        identical(checkpoint$stage, "OPTIMIZED_BEFORE_SDREPORT")
      if (optimized) {
        opt <- checkpoint$optimizer
        need(is.list(opt) && identical(opt$par, theta), "Completed optimizer checkpoint is inconsistent.")
        phase(cfg, "RECOVERING_UNCERTAINTY", "Optimizer already completed; no repeated outer optimization", result_directory = p$out)
      } else {
        persist()
        objective <- function(par) {
          value <- obj$fn(par)
          tracker$evaluations <- tracker$evaluations + 1L
          if (is.finite(value) && all(is.finite(par))) {
            if (value < tracker$nll) {
              tracker$nll <- value
              tracker$theta <- par
              tracker$parameters <- obj$env$parList(par = obj$env$last.par)
            }
            if (as.numeric(Sys.time()) - tracker$last_write >= cfg$checkpoint_seconds) {
              persist()
              phase(cfg, "FITTING", paste("Evaluations", tracker$evaluations,
                "| best objective", signif(tracker$nll, 11), "| checkpoint saved"), result_directory = p$out)
            }
          }
          value
        }
        phase(cfg, "FITTING", paste(length(p$data$y), "annual rows; original 14-term Mundlak design"), result_directory = p$out)
        opt <- stats::nlminb(theta, objective, obj$gr,
          control = list(iter.max = cfg$optimizer_iterations, eval.max = cfg$optimizer_evaluations,
                         rel.tol = 1e-10, trace = 0L))
        tracker$nll <- obj$fn(opt$par)
        need(is.finite(tracker$nll), "Optimizer returned a non-finite marginal objective.")
        opt$objective <- tracker$nll
        tracker$theta <- opt$par
        tracker$parameters <- obj$env$parList(par = obj$env$last.par)
        persist("OPTIMIZED_BEFORE_SDREPORT", opt)
      }
      phase(cfg, "UNCERTAINTY", "Optimizer checkpoint saved; computing covariance", result_directory = p$out)
      memory_snapshot(p$out, "before_sdreport")
      sd <- TMB::sdreport(obj, par.fixed = opt$par, getJointPrecision = FALSE,
        getReportCovariance = FALSE, skip.delta.method = TRUE)
      gradient <- as.numeric(obj$gr(opt$par))
      value <- obj$fn(opt$par)
      parameters <- obj$env$parList(par = obj$env$last.par)
      opt$objective <- value
      raw <- names(opt$par)
      ix <- which(raw == "beta")
      need(length(ix) == length(p$X_names) && length(ix) == 14L, "Fixed coefficient mapping changed.")
      beta <- stats::setNames(unname(opt$par[ix]), p$X_names)
      covariance <- sd$cov.fixed
      V <- covariance[ix, ix, drop = FALSE]
      dimnames(V) <- list(p$X_names, p$X_names)
      fit <- list(version = VERSION, signature = p$signature, group = cfg$group,
        nobs = length(p$data$y), parameters = parameters, optimizer = opt,
        beta = beta, V = V, cov_fixed = covariance, gradient = gradient,
        outer_parameters = sd$par.fixed, outer_raw_names = raw,
        outer_names = outer_labels(opt$par, p$X_names), pdHess = isTRUE(sd$pdHess),
        source = p$source[c("signature", "checkpoint", "stats")],
        numerical_threshold = cfg$scaled_gradient_threshold)
      fit$checks <- demo_phylo_checks(fit, p, cfg)
      save_atomic(fit, final_file)
      memory_snapshot(p$out, "after_fit_saved")
      okay <- demo_phylo_export_results(fit, p, cfg)
      phase(cfg, if (isTRUE(okay)) "COMPLETED" else "NUMERICAL_REVIEW_REQUIRED",
        if (isTRUE(okay)) "Fit saved and sensitivity exports ready" else "Checkpoint retained; inspect fit_checks and inference_status",
        TRUE, p$out)
      list(output_dir = p$out, fit_ok = isTRUE(fit$checks$fit_ok), inference_ready = isTRUE(okay))
    }, error = function(error) {
      writeLines(conditionMessage(error), file.path(cfg$monitor_directory, "worker_error.txt"))
      phase(cfg, "STOPPED_WITH_ERROR", conditionMessage(error), TRUE,
        if (is.null(cfg$result_directory)) NA_character_ else cfg$result_directory)
      list(output_dir = cfg$result_directory, fit_ok = FALSE, inference_ready = FALSE,
           error = conditionMessage(error))
    })
  }

  audit_sources <- function(project_root = ROOT, group = NULL, source_run = NULL) {
    root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
    if (!is.null(source_run)) source_run <- absolute_path(source_run, root)
    value <- demo_phylo_audit_sources(root, group = group, source_run = source_run)
    print(value[intersect(c("group", "source_directory", "version", "status",
      "numerical_checks_ok", "accepted", "reason"), names(value))], row.names = FALSE)
    invisible(value)
  }

  run <- function(project_root = ROOT, group = "multivoltine", source_run = NULL,
    tree_file = "data/EUROPEAN_BUTTERFLIES_MCC_CA.nwk", check_only = FALSE,
    checkpoint_seconds = 120, minimum_free_GiB = 25, minimum_disk_GiB = 25,
    maximum_single_basis_GiB = 2, optimizer_iterations = 1200L,
    optimizer_evaluations = 1800L, scaled_gradient_threshold = 0.05,
    refine = FALSE, expected_samples = NULL) {
    packages <- dependencies()
    need(length(group) == 1L && group %in% GROUPS, "Choose multivoltine or univoltine.")
    for (x in list(check_only, refine)) need(is.logical(x) && length(x) == 1L && !is.na(x), "Invalid logical option.")
    for (x in list(checkpoint_seconds, minimum_free_GiB, minimum_disk_GiB, maximum_single_basis_GiB,
                   optimizer_iterations, optimizer_evaluations, scaled_gradient_threshold))
      need(is.numeric(x) && length(x) == 1L && is.finite(x) && x > 0, "Invalid numeric setting.")
    need(optimizer_iterations == as.integer(optimizer_iterations) &&
      optimizer_evaluations == as.integer(optimizer_evaluations), "Optimizer limits must be integers.")
    root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
    if (!is.null(source_run)) source_run <- absolute_path(source_run, root)
    resolved <- demo_phylo_resolve_source(root, group, source_run)
    free <- ps::ps_system_memory()$avail / 1024^3
    message("RAM disponible: ", round(free, 1), " GiB")
    if (!check_only) need(free >= minimum_free_GiB, "Insufficient available RAM.")
    for (key in c("pheno.demo.mundlak.phylo.job", "pheno.env.mundlak.phylo.job",
                  "phenoimpact.env.mundlak.job", "phenoimpact.mundlak.qc.job")) {
      old <- getOption(key)
      alive <- if (is.null(old)) FALSE else tryCatch(old$process$is_alive(), error = function(e) FALSE)
      need(!isTRUE(alive), paste("A model worker is already active:", key))
    }
    cfg <- resolved
    cfg$project_root <- root
    cfg$group <- group
    cfg$tree_file <- absolute_path(tree_file, root)
    cfg$output_root <- file.path(root, "output/population_trends/mundlak_demography_phylogeny_15km")
    dir.create(cfg$output_root, recursive = TRUE, showWarnings = FALSE)
    cfg$output_root <- normalizePath(cfg$output_root, winslash = "/", mustWork = TRUE)
    cfg$compiled_directory <- file.path(cfg$output_root, "compiled")
    cfg$signature_file_hashes <- digest::digest(demo_phylo_cpp_source(), algo = "sha256", serialize = FALSE)
    cfg$scaled_gradient_threshold <- scaled_gradient_threshold
    cfg$maximum_single_basis_GiB <- maximum_single_basis_GiB
    cfg$expected_samples <- expected_samples
    cfg$check_only <- check_only
    cfg$refine <- refine
    cfg$checkpoint_seconds <- checkpoint_seconds
    cfg$minimum_disk_GiB <- minimum_disk_GiB
    cfg$optimizer_iterations <- as.integer(optimizer_iterations)
    cfg$optimizer_evaluations <- as.integer(optimizer_evaluations)
    cfg$monitor_directory <- tempfile(paste0("monitor_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_"), cfg$output_root)
    dir.create(cfg$monitor_directory)
    cfg$progress_file <- file.path(cfg$monitor_directory, "progress.rds")
    cfg$started_at <- Sys.time()
    cfg$launch_lock <- file.path(cfg$output_root, "worker.lock")
    cfg$lock_token <- paste(Sys.getpid(), format(Sys.time(), "%Y%m%d%H%M%OS6"), sep = "_")
    acquire_lock(cfg$launch_lock, cfg$lock_token)
    launched <- FALSE
    on.exit(if (!launched) release_lock(cfg), add = TRUE)
    phase(cfg, "QUEUED", if (check_only) "Preflight only" else "One demographic phylogenetic fit")
    log <- file.path(cfg$monitor_directory, "worker.log")
    process <- callr::r_bg(worker, args = list(cfg = cfg), libpath = .libPaths(),
      stdout = log, stderr = "2>&1", supervise = TRUE, user_profile = FALSE,
      system_profile = FALSE, package = TRUE, wd = root,
      env = c(callr::rcmd_safe_env(), OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1"))
    launched <- TRUE
    owner <- list(pid = process$get_pid(), created = tryCatch(
      as.numeric(ps::ps_create_time(ps::ps_handle(process$get_pid()))), error = function(e) NA_real_), token = cfg$lock_token)
    if (lock_matches(cfg$launch_lock, cfg$lock_token)) try(save_atomic(owner,
      file.path(cfg$launch_lock, "owner.rds")), silent = TRUE)
    save_atomic(cfg, file.path(cfg$monitor_directory, "launch_config.rds"))
    job <- list(process = process, config = cfg, progress_file = cfg$progress_file,
      log_file = log, monitor_directory = cfg$monitor_directory)
    options(pheno.demo.mundlak.phylo.job = job)
    message(if (check_only) "Comprovacio previa demografica iniciada." else "Un model demografic Mundlak filogenetic iniciat.",
      "\nGrup: ", group, "\nSeguiment: watch_demography_mundlak_phylogeny()\nLog: ", log)
    invisible(job)
  }

  status <- function(job = getOption("pheno.demo.mundlak.phylo.job"), quiet = FALSE) {
    need(!is.null(job), "First run run_demography_mundlak_phylogeny().")
    value <- tryCatch(readRDS(job$progress_file), error = function(e) NULL)
    if (is.null(value)) return(invisible(NULL))
    alive <- job$process$is_alive()
    if (!alive && !isTRUE(value$terminal)) value$stage <- "WORKER_STOPPED_UNEXPECTEDLY"
    if (!quiet) {
      cat("Worker actiu:", alive, "|", value$stage, "|", value$group, "\n", value$detail, "\n")
      if (!is.na(value$result_directory)) cat("Resultats:", value$result_directory, "\n")
      if (!alive && value$stage == "WORKER_STOPPED_UNEXPECTEDLY") cat("Log:", job$log_file, "\n")
    }
    invisible(value)
  }

  watch <- function(job = getOption("pheno.demo.mundlak.phylo.job"), interval_seconds = 60) {
    need(!is.null(job), "First run run_demography_mundlak_phylogeny().")
    need(length(interval_seconds) == 1L && is.finite(interval_seconds) && interval_seconds >= 1, "Invalid interval.")
    tryCatch({
      repeat {
        status(job)
        if (!job$process$is_alive()) break
        until <- Sys.time() + interval_seconds
        while (Sys.time() < until && job$process$is_alive()) Sys.sleep(1)
      }
      invisible(tryCatch(job$process$get_result(), error = function(e) {
        message("Worker aturat; checkpoints conservats. ", conditionMessage(e), "\nLog: ", job$log_file)
        NULL
      }))
    }, interrupt = function(e) {
      message("Monitor pausat; l'ajustament continua. Repeteix watch_demography_mundlak_phylogeny().")
      invisible(job)
    })
  }

  stop_job <- function(job = getOption("pheno.demo.mundlak.phylo.job")) {
    if (is.null(job)) return(invisible(FALSE))
    if (job$process$is_alive()) job$process$kill_tree()
    if (!job$process$is_alive()) release_lock(job$config)
    message("Worker aturat; checkpoints conservats.")
    invisible(TRUE)
  }

  list(run = run, status = status, watch = watch, stop = stop_job, audit_sources = audit_sources)

})

run_demography_mundlak_phylogeny <- function(...) pheno_demography_mundlak_phylogeny$run(...)
audit_demography_mundlak_phylogeny_sources <- function(...) pheno_demography_mundlak_phylogeny$audit_sources(...)
resume_demography_mundlak_phylogeny <- function(...) pheno_demography_mundlak_phylogeny$run(...)
status_demography_mundlak_phylogeny <- function(...) pheno_demography_mundlak_phylogeny$status(...)
watch_demography_mundlak_phylogeny <- function(...) pheno_demography_mundlak_phylogeny$watch(...)
stop_demography_mundlak_phylogeny <- function(...) pheno_demography_mundlak_phylogeny$stop(...)

message("Funcions carregades. Primer: job <- run_demography_mundlak_phylogeny(check_only = TRUE). Seguiment: watch_demography_mundlak_phylogeny(job).")
