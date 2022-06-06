trunc_to_integer <- function(X, Ximp, n_level = 20, cat_index = NULL){
  if (is.null(Ximp)) return(NULL)
  p = ncol(Ximp)
  noncont_cols = get_cat_cols(X, n_level)
  trunc_cols = setdiff(noncont_cols, cat_index)
  for (j in trunc_cols){
    x = X[,j]
    xu = unique(x[!is.na(x)])
    Ximp[,j] = trunc_rating(Ximp[,j], xmin = min(xu), xmax = max(xu))
  }
  Ximp
}

trunc_rating <- function(x, xmin, xmax){
  xnew = round(x)
  xnew[x<=xmin] = xmin
  xnew[x>=xmax] = xmax
  xnew
}

to_df_with_factor <- function(X, n_level=20, cat_index=NULL){
  X = as.data.frame(X)
  noncont_cols = get_cat_cols(X, n_level)
  for (i in noncont_cols) X[,i] = as.factor(X[,i])
  r = list()
  p = ncol(X)
  if (is.null(cat_index)) r[['ord_cols']] = noncont_cols
  else{
    stopifnot(all(cat_index %in% noncont_cols))
    r[['cat_cols']] = cat_index
    r[['ord_cols']] = setdiff(noncont_cols, cat_index)
  }
  for (i in r[['ord_cols']]) X[,i] = ordered(X[,i])
  r[['X']] = X
  r[['nfactor']] = length(noncont_cols)
  r
}

get_cat_cols <- function(X, n_level=10){
  X = as.data.frame(X)
  which(purrr::map_lgl(X, ~ length(unique(.x[!is.na(.x)]))<=n_level))
}


map_ordinal_FAMD = function(x){
  x = as.character(x)
  l = length(x)
  y = unlist(strsplit(x, split = '_', fixed = TRUE))
  as.numeric(y[(1:l)*2])
}

to_original_factor_FAMD <- function(X, nfactor, n_level = 10){
  if (nfactor == 1) return(X)
  cat_cols = which(purrr::map_lgl(X, ~ length(unique(.x[!is.na(.x)]))<=n_level))
  for (i in cat_cols) X[,i] = map_ordinal_FAMD(X[,i])
  X
}

post_processing_FAMD <- function(res, nfactor, n_level){
  if (is.null(res$error)){
    est = res$result
    Ximp = to_original_factor_FAMD(est$completeObs, nfactor = nfactor, n_level = n_level)
    to_numeric_matrix(Ximp)
  }else{
    NULL
  }
}

majority_vote = function(x){
  xt = table(x)
  as.integer(names(xt)[which.max(xt)])
}

MI_to_single <- function(Xlist, Xobs, n_level=20){
  l = length(Xlist)
  Ximp = to_numeric_matrix(Xobs)
  p = ncol(Xobs)
  noncont_cols = get_cat_cols(Xobs, n_level)
  for (name in c('cat', 'cont')){
    if (name == 'cat'){
      index = noncont_cols
      fpool = majority_vote
    }
    else{
      index = setdiff(1:p, noncont_cols)
      fpool = mean
    }
    X = Xobs[,index,drop=FALSE]
    Xlist_ = map(Xlist, ~ .x[,index,drop=FALSE])
    mis = is.na(X)
    xmis_all = as.data.frame(map(Xlist_, ~ .x[mis]))
    xmis_all   = do.call(cbind, xmis_all)
    ximp = apply(xmis_all, 1, fpool)
    #ximp = apply(map_dfr(Xlist_, ~ .x[mis]), 1, fpool)
    Ximp[,index][mis] = ximp
  }
  Ximp
}


is_mask_of <- function(Xmask, X){
  o_xmask = which(!is.na(Xmask))
  crit1 = all(o_xmask %in% which(!is.na(X)))
  if (crit1){
    r = all(Xmask[o_xmask] == X[o_xmask])
  }else{
    r = FALSE
  }
  r
}

test_at_best_para <- function(allres_list, eval_name){
  summary_list = map(allres_list, colMeans)
  best_para = names(which.min(summary_list[[eval_name]]))
  names_use = names(summary_list)
  test_name = map_lgl(names_use, ~ grepl('test', .x, fixed = TRUE))
  keep_name = c(names_use[test_name], 'time')
  #test_name = c('test_all', 'test_cat', 'time')
  e = map(allres_list[keep_name], ~ unname(.x[,best_para]))
  #e[['time']] = unname(drop(allres_list[['time']]))
  list(result = e, best_para = best_para)
}

test_at_best <- function(allres_list){
  r = list(bestall = test_at_best_para(allres_list, 'val_all')$result,
           bestcat = test_at_best_para(allres_list, 'val_cat')$result)
  r
}

one_hot_matrix <- function(X, cat_index){
  Xcat = X[,cat_index,drop=FALSE]
  p = ncol(X)
  n = nrow(X)
  cat_level = numeric(p) * NA
  cat_freq = get_cat_index_freq(Xcat)
  cat_level[cat_index] = cat_freq$nlevel
  cat_index_list = create_cat_index_list(cat_level)
  cat_index_all = purrr::reduce(cat_index_list, c)
  d_cat = length(cat_index_all)
  p_cat = length(cat_index_list)
  d = p-p_cat+d_cat

  Xencoded = matrix(0, n, d)
  Xencoded[,-cat_index_all] = X[,-cat_index]
  cat_labels = list()
  for (j in names(cat_index_list)){
    jint = as.integer(j)
    relabel = cat_to_integers(X[,jint])
    cat_labels[[as.character(j)]] = relabel$xlevels
    Xencoded[,cat_index_list[[j]]] = one_hot_vec(relabel$x, K=length(relabel$xlevels))
  }
  list(X=Xencoded, cat_index_list=cat_index_list, cat_labels=cat_labels)
}

inverse_one_hot_matrix <- function(Xencoded, cat_index_list, cat_labels){
  cat_index = cat_index_from_list(cat_index_list)
  d = ncol(Xencoded)
  cat_index_all = purrr::reduce(cat_index_list, c)
  d_cat = length(cat_index_all)
  p_cat = length(cat_index_list)
  p = d - d_cat + p_cat
  n = nrow(Xencoded)
  X = matrix(0, n, p)
  X[,-cat_index] = Xencoded[,-cat_index_all]
  for (j in names(cat_index_list)){
    jint = as.integer(j)
    Xj = Xencoded[,cat_index_list[[j]]]
    Xj_argmax = apply(Xj, 1, which.max)
    X[,jint] = relabel(Xj_argmax, cat_labels[[j]])
  }
  X
}

one_hot_vec <- function(x, K){
  n = length(x)
  mis = is.na(x)
  xencoded = matrix(0, n, K)
  xencoded[mis,]=NA
  xencoded[matrix(ncol = 2, c(1:n, x))] = 1
  xencoded
}

cal_misclass = function(xhat, xobs, xtrue){
  xobs = as.numeric(as.matrix(xobs))
  xhat = as.numeric(as.matrix(xhat))
  xtrue = as.numeric(as.matrix(xtrue))
  if (is.null(xobs)) loc = !is.na(xtrue) else loc = is.na(xobs) & (!is.na(xtrue))
  1 - mean(xhat[loc] == xtrue[loc])
}

cal_misclass_scaled = function(xhat, xobs, xtrue, base_from_true = FALSE, reduce =TRUE){
  xobs = to_numeric_matrix(xobs)
  xhat = to_numeric_matrix(xhat)
  xtrue = to_numeric_matrix(xtrue)

  if (base_from_true) xbase = xtrue else xbase = xobs
  base = apply(xbase, 2, function(x) names(which.max(table(x))))
  base = as.integer(base)
  p = ncol(xhat)
  err = numeric(p)
  n = nrow(xobs)

  for (j in 1:p){
    err_base = cal_misclass(xhat = rep(base[j],n), xobs = xobs[,j], xtrue = xtrue[,j])
    err[j] = cal_misclass(xhat = xhat[,j], xobs = xobs[,j], xtrue = xtrue[,j])/err_base
  }

  remove = is.infinite(err)
  if (any(remove)){
    #stop('Perfect median imputation appears')
    warning(paste0('Perfect majority vote imputation appears in ', sum(remove), ' vars'))
  }
  if (reduce) err = mean(err[!remove])
  err
}
