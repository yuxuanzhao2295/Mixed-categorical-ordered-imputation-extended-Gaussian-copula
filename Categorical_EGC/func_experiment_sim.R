source('func_load_eval.R')
source('func_load_data.R')
source('func_load_imputation.R')
source('func_utils.R')

options(warn = 1)
library(purrr)

get_f_data_sim <- function(pcat, p_noncat=10, p_ord=5, K=3,
                           val_ratio=0, mask_fraction = 0.3, trans_cont = NULL, mask_type = 'MCAR'){
  if (is.null(trans_cont)) trans_cont = function(x) qexp(pnorm(x), rate = 1/3)
  partial(data_sim_nominal_cop,
          p_cont = p_noncat, p_ord=p_ord,
          n_level = rep(K, pcat),
          n = 2000,
          mask_fraction = mask_fraction, mask_type = mask_type,
          val_ratio = val_ratio,
          mu_scale = 0.25,
          trans_cont = trans_cont ,
          old =FALSE)
}



test_sim_mixed <- function(func_calls, K=6, mask_fraction=0.3, mask_type = 'MCAR', seed=1){
  names_run = names(func_calls)
  Res_all = list()
  for (name in names_run){
    Res = list()
    if (name %in% c('imputeFAMD', 'softImpute')) val_ratio = 0.2 else val_ratio = 0
    if (name == 'EGC_our') addargs = c('mu_loc', 'seed_num') else addargs = NULL
    for (i in c(1,3,5)){
      f_data = get_f_data_sim(pcat=i,K=K,mask_fraction = mask_fraction, mask_type = mask_type, val_ratio=val_ratio)
      Res[[paste0('num', i)]] = imputation_pipeline_rep(seed=seed,
                                                    f_data=f_data,
                                                    f_imp=func_calls[[name]],
                                                    f_eval=partial(eval_cont_ord_cat, average = TRUE),
                                                    use_addargs = TRUE, addargs = addargs)
    }
    if (name %in% c('imputeFAMD', 'softImpute')) rall = map(Res, summarize_mlutiple_fcall)
    else rall = map(Res, summarize_single_fcall)
    Res_all[[name]] = rall
    print(paste0('finish', name))
  }
  Res_all
}

run_sim_mixed <- function(func_calls, K=6, mask_fraction=0.3, mask_type = 'MCAR', nrep=10){
  names_run = names(func_calls)
  Res_all = list()
  for (name in names_run){
    Res = list()
    if (name %in% c('imputeFAMD', 'softImpute')) val_ratio = 0.2 else val_ratio = 0
    if (name == 'EGC_our') addargs = c('mu_loc', 'seed_num') else addargs = NULL
    for (i in c(1,3,5)){
      f_data = get_f_data_sim(pcat=i,K=K,mask_fraction = mask_fraction, mask_type = mask_type, val_ratio=val_ratio)
      Res[[paste0('num', i)]] = imputation_pipeline(nrep=nrep,
                                                    f_data=f_data,
                                                    f_imp=func_calls[[name]],
                                                    f_eval=partial(eval_cont_ord_cat, average = TRUE),
                                                    use_addargs = TRUE, addargs = addargs)
    }
    if (name %in% c('imputeFAMD', 'softImpute')) rall = map(Res, summarize_mlutiple_fcall)
    else rall = map(Res, summarize_single_fcall)
    Res_all[[name]] = rall
    print(paste0('finish', name))
  }
  Res_all
}

summarize_best_sim <- function(r){
  rs = map(r, test_at_best)
  Re = list()
  for (name in c('bestcat', 'bestall')){
    Re[[name]] = map(rs, ~ .x[[name]])
  }
  Re = map(Re, ~ map(.x, as.data.frame))
  Re
}

sim_postprocessing <- function(rall){
  for (m in c('imputeFAMD', 'softImpute')){
    if (!is.null(rall[[m]])){
      r1 = summarize_best_sim(rall[[m]])
      rall[[m]] = r1$bestall
    }
  }
  rall
}
