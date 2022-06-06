library(gcimputeR)

impute_nominal_cont_gc = function(X, cat_index,
                                  mu = NULL, corr = NULL,
                                  n_MI = 0, fast_MI = FALSE,
                                  maxit=50, eps=0.01, nlevels = 20, runiter = 0, verbose=FALSE,
                                  seed = NULL,
                                  init_trunc_sampling = FALSE,
                                  trunc_method = 'Iterative', n_sample=5000, n_update=1,
                                  Z = NULL,
                                  old = FALSE,
                                  ...){
  if (!is.null(seed)) set.seed(seed)
  n = nrow(X)
  p = ncol(X)
  X = as.numeric(as.matrix(X))
  dim(X) = c(n,p)

  # check cat_index
  if (is.logical(cat_index)){
    if (length(cat_index) != p) stop('invalid cat_index')
    cat_index_int = which(cat_index)
  }else{
    cat_index = as.integer(cat_index)
    if (min(cat_index)<1 | max(cat_index)>p) stop('invalid cat_index')
    cat_index_int = cat_index
  }

  cat_freq = get_cat_index_freq(X[,cat_index_int,drop=FALSE])
  min_nl = min(cat_freq$nlevel)
  if (min_nl <= 1){
    print(paste0('Number of categories: ', cat_freq$nlevel))
    stop('some categorical var has no more than 1 level')
  }else if (min_nl == 2){
    idx = which(cat_freq$nlevel > 2)
    cat_index_int = cat_index_int[idx]
    pcat = length(cat_freq$nlevel)
    print(paste0(pcat-length(idx), ' of ', pcat, ' categoricals have only two categories: treated as binary'))
  }

  # relabel
  cat_labels = list()
  for (j in cat_index_int){
    relabel = cat_to_integers(X[,j])
    cat_labels[[as.character(j)]] = relabel$xlevels
    X[,j] = relabel$x
  }
  cat_index = index_int_to_logi(cat_index_int, p)

  X_cat = X[,cat_index,drop=FALSE]
  if (any(cat_index)) stopifnot(apply(X_cat, 2, check_cat_label))
  X_noncat = X[,!cat_index,drop=FALSE]

  # Do not allow empty row
  if (any(apply(X, 1, function(x){sum(!is.na(x))}) == 0)) stop("remove empty row")
  # Do not allow column with only one level
  if (any(apply(X, 2, function(x){length(unique(x[!is.na(x)]))}) <= 1)){
    stop('remove column with only 0 or 1 unique value')
  }

  cat_level = numeric(p) * NA
  cat_freq = get_cat_index_freq(X_cat)
  if (old) cat_level[cat_index] = cat_freq$nlevel - 1 else cat_level[cat_index] = cat_freq$nlevel #***!
  cat_index_list = create_cat_index_list(cat_level)


  cat_index_all = purrr::reduce(cat_index_list, c)
  d_cat = length(cat_index_all)
  d = ncol(X_noncat)+d_cat

  # TEST POINT 1 !!!
  dcat_index = logical(d)
  dcat_index[cat_index_all] = TRUE
  dord_index = logical(d)
  ord_in_noncat = apply(X_noncat, 2, function(x) {length(unique(x)) <= nlevels})
  dord_index[-cat_index_all] = ord_in_noncat
  d_index = dcat_index | dord_index
  dord = sum(d_index)
  cat_in_d = index_int_to_logi(match(which(dcat_index), which(d_index)), l = dord)

  c_index = !d_index
  cat_input = list(X_cat = X_cat, cat_index_list = cat_index_list, cat_index_all = cat_index_all, old = old)

  # categorical dimensions
  if (any(dcat_index)){
    if (is.null(mu)){
      mu_est = get_cat_mu(cat_freq$freq, old = old, verbose = verbose) #*?
      mu = mu_est$mu
    }
  }else{
    mu = NULL
  }


  Lower = matrix(NA, n, dord)
  Upper = Lower
  bounds = get_cat_bounds(X_cat, mu, cat_index_list,check = TRUE, old = old) #***
  Lower[,cat_in_d] = bounds$lower
  Upper[,cat_in_d] = bounds$upper

  bounds = range_transform(X_noncat[,ord_in_noncat,drop=FALSE], type= 'ordinal')
  Lower[,!cat_in_d] = bounds$Lower
  Upper[,!cat_in_d] = bounds$Upper

  if (is.null(Z)){
    Z = initZ(Lower, Upper, X,
              cat_index, ord_in_noncat, cat_in_d, c_index, dord_index, dcat_index,
              cat_index_list, Z_cont=NULL, m=1, method = 'univariate_mean', old = old)
  }

  if (is.null(corr)){
    fit_em = em_mixedgc(Z, Lower, Upper,
                        d_index=d_index, dcat_index=dcat_index,
                        cat_input = cat_input,
                        maxit = maxit, eps = eps, runiter=runiter, verbose=verbose,
                        trunc_method = trunc_method, n_sample=n_sample, n_update=n_update, ...) #***!
    Zimp = fit_em$Zimp
    corr = fit_em$corr
    loglik = fit_em$loglik
    Z = fit_em$Z
  }else{
    out = latent_operation('fillup',
                           Z, Lower, Upper,
                           d_index=d_index, dcat_index=dcat_index,
                           cat_input = cat_input,
                           corr = corr,
                           n_update = n_update, n_sample = n_sample, trunc_method = trunc_method) #***!
    Zimp = out$Zimp
    loglik = NULL
  }

  # Impute X using Imputed Z
  Ximp = latent_to_observed(Zimp, X, mu, cat_labels,
                            cat_index = cat_index, ord_in_noncat = ord_in_noncat,
                            cat_index_all = cat_index_all, cat_index_list = cat_index_list, old = old)

  if (n_MI > 0){
    call_sample <- function(Z, trunc_method, n_MI){
      out = latent_operation('sample',
                             Z, Lower, Upper,
                             d_index=d_index, dcat_index=dcat_index,
                             cat_input = cat_input,
                             corr = corr,
                             n_update = n_update, n_sample = n_sample, trunc_method = trunc_method,
                             n_MI = n_MI)
      Zfill = drop(out$Zimp_sample)
      Zfill
    }

    if (fast_MI){
      Z_cont = Z[,c_index,drop=FALSE]
      Zinits = initZ(Lower, Upper, X,
                     cat_index, ord_in_noncat, cat_in_d, c_index, dord_index, dcat_index,
                     cat_index_list,
                     Z_cont=Z_cont, m=n_MI, method = 'sampling', old = old)
      Zimps = map(Zinits, ~ call_sample(.x, 'Iterative', 1)) # no sampling
    }else{
      Zimps = call_sample(Z, 'Sampling', n_MI)
      Zimps = map(1:n_MI, ~ Zimps[,,.x])
    }

    call_ZtoX <- function(Zimp){
      ximp = latent_to_observed(Zimp, X, mu, cat_labels,
                                cat_index = cat_index, ord_in_noncat = ord_in_noncat,
                                cat_index_all = cat_index_all, cat_index_list = cat_index_list, old = old)
      ximp
    }
    Ximp_MI = map(Zimps, call_ZtoX)
  }else Ximp_MI = NULL

  noncat_index = setdiff(1:p, cat_index_int)
  ord_index = noncat_index[ord_in_noncat]
  var_types = list(continuous = noncat_index[!ord_in_noncat],
                   ordinal = ord_index,
                   categorical = cat_index_int)
  return(list(Ximp = Ximp, Ximp_MI = Ximp_MI,
              corr=corr, mu = mu,
              loglik=loglik,
              var_types = var_types,
              cat_index_list = cat_index_list,
              Z = Z))
}

cat_to_integers <- function(x){
  x = droplevels(as.factor(x))
  xlevels = levels(x)
  nlevel = nlevels(x)
  levels(x) = 1:nlevel
  list(x=as.numeric(x), xlevels = xlevels)
}

Ximp_transform_cat <- function(Z_cat, X_cat, cat_index_list, old = FALSE){
  if (ncol(Z_cat) != sum(purrr::map_int(cat_index_list, length))) stop('something wrong')
  cat_index_list = adjust_index_list(cat_index_list)
  Ximp_cat = X_cat
  for (j in seq_along(cat_index_list)){
    index_m = is.na(X_cat[,j])
    index_cat = cat_index_list[[j]]
    zmis = Z_cat[index_m,index_cat,drop=FALSE]
    Ximp_cat[index_m,j] = apply(zmis, 1, nominal_z_to_x_col, old = old)
  }
  Ximp_cat
}

nominal_z_to_x_col <- function(z, old = FALSE){
  argmax = which.max(z)
  if (old){
    if (z[argmax]<0) argmax = 1 else argmax = argmax + 1
  }
  argmax
}


latent_to_observed <- function(Zimp, X, mu, cat_labels, ord_in_noncat,
                               cat_index, cat_index_all, cat_index_list, old=FALSE){
  d_cat = length(cat_index_all)
  n = nrow(Zimp)
  Z_cat = Zimp[,cat_index_all] + matrix(mu, n, d_cat, byrow = TRUE)
  Ximp = X
  X_cat = X[,cat_index,drop=FALSE]
  X_noncat = X[,!cat_index,drop=FALSE]
  Ximp[,cat_index] = Ximp_transform_cat(Z_cat = Z_cat, X_cat = X_cat,
                                        cat_index_list = cat_index_list, old = old) #***!
  Ximp[,!cat_index] = Ximp_transform(Z = Zimp[,-cat_index_all], X = X_noncat, d_index = ord_in_noncat)

  cat_index_int = which(cat_index)
  for (j in cat_index_int) Ximp[,j] = relabel(Ximp[,j], cat_labels[[as.character(j)]])

  Ximp
}

index_int_to_logi <- function(index, l){
  out = logical(l)
  out[index] = TRUE
  out
}

cat_index_from_list <- function(l) as.integer(names(l))

check_cat_label <- function(x){
  x = x[!is.na(x)]
  xmax = max(x)
  xmin = min(x)
  nlevel = nlevels(as.factor(x))
  xmin == 1 & xmax == nlevel
}

get_cat_bounds <- function(X_cat, mu, cat_index_list, check=FALSE, old = FALSE){
  # TODO check values of X_cat
  if (old) return(get_cat_bounds_old(X_cat, mu, cat_index_list, check))
  d_cat = sum(purrr::map_int(cat_index_list, length))
  if (d_cat != length(mu)) stop('invalid input')
  p_cat = length(cat_index_list)
  if (p_cat != ncol(X_cat)) stop('invalid input')
  n = nrow(X_cat)
  lower = matrix(NA, n, d_cat)
  upper = matrix(NA, n, d_cat)
  incat_index_list = adjust_index_list(cat_index_list)
  for (j in 1:p_cat){
    index_o = !is.na(X_cat[,j])
    x_cat_obs = X_cat[index_o,j]
    n_obs = length(x_cat_obs)
    # initialize to (0, Inf): at argmax, we want (-inf, inf), at other loc, we want (mu_j - mu_argmax, inf)
    # index in 1,...,d_cat
    index_cat = incat_index_list[[j]]
    dj_cat = length(index_cat)
    l_o = matrix(0, n_obs, dj_cat)
    u_o = l_o + Inf

    # adjust for mean
    # length d_cat
    mu_j = mu[index_cat]
    # z_argmax - z_{-argmax} + mu_j[argmax] - mu_j[-argmax] >=0
    # thus z_argmax - z_{-argmax} >= mu_j[-argmax] - mu_j[argmax] (RHS computed below)
    mu_diff = matrix(mu_j, n_obs, dj_cat, byrow = TRUE)  - matrix(mu_j[x_cat_obs], n_obs, dj_cat, byrow = FALSE)
    l_o = l_o + mu_diff
    # no constraints at argmax, thus -Inf lower
    argmax_coor = matrix(c(1:n_obs, x_cat_obs), nrow = n_obs)
    l_o[argmax_coor] = -Inf

    lower[index_o, index_cat] = l_o
    upper[index_o, index_cat] = u_o
  }

  if (check){
    for (i in 1:n){
      ind1 = is.na(lower[i,])
      ind2 = get_cat_slicing_index(X_cat[i,], incat_index_list, keep = 'missing', d_cat=d_cat)$cat
      if (!all(ind1==ind2)) stop('something wrong!')
    }
  }
  list(lower = lower, upper = upper)
}

relabel <- function(x, label){
  stopifnot(check_cat_label(x))
  x = as.factor(x)
  levels(x) = label
  as.numeric(levels(x))[x]
}


adjust_index_list <- function(index_list){
  start = 1
  for (i in seq_along(index_list)){
    vals = index_list[[i]]
    index_list[[i]] = vals - vals[1] + start
    start = start + length(vals)
  }
  index_list
}

get_cat_slicing_index <- function(x_cat, cat_index_list, keep = 'observed', d_cat=NULL){
  if (is.character(keep)){
    if (keep == 'observed') index_incat = !is.na(x_cat)
    else if (keep == 'missing') index_incat = is.na(x_cat)
    else stop('invalid char keep')
  }else if (is.logical(keep) & length(keep) == length(x_cat)){
    index_incat = keep
  } else stop('invalid keep')

  if (is.null(d_cat)) d_cat= sum(purrr::map_int(cat_index_list, length))
  index_cat = logical(d_cat)
  if (any(index_incat)){
    intindex_cat = purrr::reduce(cat_index_list[index_incat], c)
    index_cat[intindex_cat] = TRUE
  }

  list('incat'=index_incat, 'cat'=index_cat)
}

initZ <- function(Lower, Upper, X,
                  cat_index, ord_in_noncat, cat_in_d, c_index, dord_index, dcat_index,
                  cat_index_list, Z_cont=NULL, m=1, method = 'univariate_mean', old = FALSE){
  X_cat = X[,cat_index,drop=FALSE]
  X_noncat = X[,!cat_index,drop=FALSE]

  if (any(c_index) & is.null(Z_cont)){
    Z_cont = range_transform(X_noncat[,!ord_in_noncat,drop=FALSE], type = 'continuous')$Z
  }else if (is.null(Z_cont)) stopifnot(all(!c_index))
  else stopifnot(ncol(Z_cont) == sum(c_index))

  n = nrow(X)
  d = length(dcat_index)

  call_initZ_noncat <- function(){
    Zinit = initZ_noncat(Lower, Upper, X_cat, cat_in_d, cat_index_list, old = old, method = method)
    Zinit
  }
  setZ <- function(Zinit){
    Z = matrix(NA, n, d)
    Z[,dord_index] = Zinit$Zord
    Z[,dcat_index] = Zinit$Zcat
    Z[,c_index] = Z_cont
    Z
  }

  if (m == 1){
    Zinit = call_initZ_noncat()
    out = setZ(Zinit)
  }else{
    Zinits = purrr::map(1:m, ~ call_initZ_noncat())
    out = purrr::map(Zinits, setZ)
  }

  out
}

initZ_noncat = function(Lower, Upper, X_cat,
                        cat_in_d, cat_index_list, method = 'univariate_mean', old = FALSE){
  Z_init = initZ_interval_truncated(Lower, Upper, method = method)
  Zord = Z_init[,!cat_in_d,drop=FALSE]
  Zcat = Z_init[,cat_in_d,drop=FALSE]
  Zcat = Z_to_original_trunc(Zcat, X_cat, cat_index_list, old=old)
  list(Zord = Zord, Zcat = Zcat)
}


x_to_A <- function(x, cat_index_list, d_cat=NULL, adjust=TRUE, test=TRUE, old = FALSE){
  if (adjust) cat_index_list = adjust_index_list(cat_index_list)
  if (test){
    if (!is.null(d_cat)){
      if (d_cat != sum(purrr::map_int(cat_index_list, length))) stop('invalid input')
    }
  }
  if (is.null(d_cat)) d_cat = sum(purrr::map_int(cat_index_list, length))
  if (any(is.na(x))) stop('invalid x')

  if (old){
    index_notbase = get_cat_slicing_index(x, cat_index_list, keep = x!=1, d_cat = d_cat)
    if (any(index_notbase$incat)){
      A = diag(nrow = d_cat, ncol = d_cat)
      # for each xi != 1
      for (i in which(index_notbase$incat)){
        index = cat_index_list[[i]]
        x_index = x[i]-1
        Ai = -diag(length(index))
        Ai[,x_index] = 1
        A[index,index] = Ai
      }
    }else{
      A = NULL
    }
  }else{
    A = diag(nrow = d_cat, ncol = d_cat)
    p = length(x)
    for (i in 1:p){
      index = cat_index_list[[i]]
      Ai = -diag(length(index))
      Ai[,x[i]] = 1
      A[index,index] = Ai
    }
  }
  A
}

Z_to_original_trunc <- function(Z, X_cat, cat_index_list, old = FALSE){
  n = nrow(Z)
  for (i in 1:n){
    z = Z[i,]
    x_cat = X_cat[i,]
    obs_indices = !is.na(z)
    cat_obs = !is.na(x_cat)
    if (any(cat_obs)){
      A = x_to_A(x = x_cat[cat_obs], cat_index_list = cat_index_list[cat_obs], old = old)
      if (!is.null(A)) z[obs_indices] = A %*% z[obs_indices]
      Z[i,] = z
    }
  }
  Z
}


create_cat_index_list <- function(cat_index_level){
  cat_index_list = list()
  start = 1
  for (i in seq_along(cat_index_level)){
    l = cat_index_level[i]
    if (!is.na(l)){
      cat_index_list[[as.character(i)]] = start:(start+l-1)
      start = start+l
    }else start = start+1
  }
  cat_index_list
}
