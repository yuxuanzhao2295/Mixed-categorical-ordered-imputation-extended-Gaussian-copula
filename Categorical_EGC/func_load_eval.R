eval_cont_ord_cat <- function(est, xtrue, xobs,  prefix='', average = FALSE){
  ximp = est$Ximp
  cat_index = est$cat_index
  if (is.null(cat_index)){
    print('must provide cat_index for cont-cat evaluation')
    stop('stopped')
  }
  if (is.logical(cat_index)) cat_index = which(cat_index)

  noncat_cols = get_cat_cols(xobs, n_level = 20)
  stopifnot(all(cat_index %in% noncat_cols))

  err = list()
  err[[paste0(prefix, 'cat')]] = gcimputeR::cal_misclass(ximp[,cat_index],
                                                         xobs = xobs[,cat_index],
                                                         xtrue = xtrue[,cat_index])
  if (length(cat_index) < length(noncat_cols)){
    ord_index = setdiff(noncat_cols, cat_index)
    err[[paste0(prefix, 'ord')]] = gcimputeR::cal_mae(ximp[,ord_index],
                                                           xobs = xobs[,ord_index],
                                                           xtrue = xtrue[,ord_index])
  }

  if (length(noncat_cols) < ncol(xtrue)){
    err[[paste0(prefix, 'cont')]] = gcimputeR::cal_mae(ximp[,-noncat_cols],
                                                           xobs = xobs[,-noncat_cols],
                                                           xtrue = xtrue[,-noncat_cols])
  }

  err = unlist(err)
  if (average) err[paste0(prefix, 'all')] = mean(err)

  err
}


eval_cat_noncat <- function(est, xtrue, xobs,  noncont_levels = 20,
                            round = FALSE, reduce = TRUE,
                            base_from_true = FALSE, verbose = FALSE,
                            prefix=''){
  ximp = est$Ximp
  cat_index = est$cat_index
  if (is.null(cat_index)){
    print('must provide cat_index for cont-cat evaluation')
    stop('stopped')
  }
  if (is.logical(cat_index)) cat_index = which(cat_index)

  noncat_cols = get_cat_cols(xobs, n_level = noncont_levels)
  stopifnot(all(cat_index %in% noncat_cols))

  ximp = to_numeric_matrix(ximp)
  xtrue = to_numeric_matrix(xtrue)
  xobs = to_numeric_matrix(xobs)

  p = ncol(ximp)
  e = numeric(p)
  e[-cat_index] = cal_mae_scaled(xhat = ximp[,-cat_index,drop=FALSE],
                                 xobs = xobs[,-cat_index,drop=FALSE],
                                 xtrue = xtrue[,-cat_index,drop=FALSE],
                                 base_from_true = base_from_true,
                                 round = round,
                                 reduce = FALSE, verbose = verbose)
  e[cat_index] = cal_misclass_scaled(xhat = ximp[,cat_index,drop=FALSE],
                                     xobs = xobs[,cat_index,drop=FALSE],
                                     xtrue = xtrue[,cat_index,drop=FALSE],
                                     base_from_true = base_from_true,
                                     reduce = FALSE)
  if (reduce){
    e[!is.finite(e)] = NA
    r = unlist(list(cat = mean(e[cat_index], na.rm = TRUE),
                noncat = mean(e[-cat_index], na.rm = TRUE),
                all = mean(e, na.rm = TRUE)
                ))
    names(r) = paste0(prefix, names(r))
  }else{
    r = e
  }

  r
}
