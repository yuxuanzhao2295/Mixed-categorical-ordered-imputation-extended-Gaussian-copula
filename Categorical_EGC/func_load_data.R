library(gcimputeR)

data_masker_general <- function(data, cat_index, seed = 1,
                                mask_fraction = 0.2, mask_type = 'MCAR', val_ratio = 0, label=NULL){
  d_masked = gcimputeR::mask_MCAR(data, mask_fraction = mask_fraction, seed = seed)

  if (mask_type == 'MCAR'){
    d_masked = mask_MCAR(data, mask_fraction = mask_fraction, seed = seed)
  }else if (mask_type == 'MNAR'){
    maskby = add_maskby(data, cat_index)
    d_masked = mask_MNAR_three_piece(X = data, BY = maskby)
  }else if (mask_type == 'MAR'){
    maskby = add_maskby(data, cat_index)
    pobs_ratio = 1/3
    misprob = misprob_selfmask_three_piece(data, maskby,
                                           normal_m = mask_fraction/(1-pobs_ratio))
    misprob = adjust_misprob_MAR_matrix(misprob, pobs_ratio = pobs_ratio)
    d_masked = mask_by_misprob(data, misprob)
  }

  which_empty = which(apply(d_masked, 1, function(x){sum(!is.na(x))}) == 0)
  for (row in which_empty){
    obs_loc = which(!is.na(data[row,]))
    if (length(obs_loc)==1) index = obs_loc else index = sample(obs_loc,1)
    d_masked[row, index] = data[row, index]
  }


  split = split_mask_val_test(d_masked, data, val_ratio = val_ratio, seed=seed)
  split[['cat_index']] = cat_index
  split[['label']] = label
  split
}




data_sim_nominal_cop <- function(seed, n_level, p_cont, p_ord=0, n=1000,
                                 mask_fraction=0.1, val_ratio = 0,
                                 old=FALSE, trans_cont=NULL, mask_type = 'MCAR',
                                 ...){
  set.seed(seed)
  n_cat = length(n_level)
  if (old) p_cat_vec = n_level - 1 else p_cat_vec = n_level
  d_gen = gen_nominal_copula(n=n, p_cat_vec = p_cat_vec,
                             p_noncat = p_cont, p_ord = p_ord,
                             trans_cont = trans_cont,
                             old=old, ...)
  if (mask_type == 'MCAR'){
    d_gen$xmask = mask_MCAR(d_gen$x, mask_fraction)
    mis_suffix = paste0('_m', as.integer(mask_fraction*100))
  }else if (mask_type == 'MNAR'){
    #maskby = add_maskby(d_gen$x, d_gen$cat_index)
    d_gen$xmask = mask_MNAR_three_piece(X = d_gen$x, BY = NULL,
                                        normal_m = mask_fraction)
    mis_suffix = '_MNAR'
  }else if (mask_type == 'MAR'){
    #maskby = add_maskby(d_gen$x, d_gen$cat_index)
    pobs_ratio = 1/3
    misprob = misprob_selfmask_three_piece(X = d_gen$x, BY = NULL,
                                           normal_m = mask_fraction/(1-pobs_ratio))
    misprob = adjust_misprob_MAR_matrix(misprob, pobs_ratio = pobs_ratio)
    d_gen$xmask = mask_by_misprob(d_gen$x, misprob)
    mis_suffix = '_MAR'
  }

  r = list(train = d_gen$xmask, test = d_gen$x,
             cat_index = d_gen$cat_index, mu = d_gen$mu_cat, corr = d_gen$sigma)
  if (val_ratio > 0){
    split = gcimputeR::split_mask_val_test(r$train, r$test, val_ratio = val_ratio)
    r$train = split$train
    r$validation = split$validation
  }
  r[['seed_num']] = seed
  prefix = '/WAIT_FOR_COMPLETE_PYTHON_VERSION_TO_AVIOD_CALLING_PYTHON_IN_R/'
  K = unique(n_level)
  stopifnot(length(K) == 1)
  pcat = length(n_level)
  r[['mu_loc']] = paste0(prefix, 'sim_pcat', pcat, '_K', K, mis_suffix, '.csv')
  r
}
