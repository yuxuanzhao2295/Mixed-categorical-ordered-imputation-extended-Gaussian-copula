source('func_preprocessing_openml.R')
source('func_load_eval.R')
source('func_load_data.R')
source('func_load_imputation.R')
source('func_utils.R')


test_data_openml <- function(func_calls, dname, mfrac = 0.2, mask_type = 'MCAR', seed = 1){
  names_run = names(func_calls)
  Res = list()
  for (fname in names_run){
    if (fname %in% c('imputeFAMD', 'softImpute')) val_ratio = 0.2 else val_ratio = 0
    if (fname == 'EGC_our') addargs = c('mu_loc', 'seed_num') else addargs = NULL
    f_data = get_data_openml(dname, mask_fraction = mfrac, mask_type = mask_type, val_ratio = val_ratio)
    r = imputation_pipeline_rep(seed=seed,
                                f_data = f_data, f_imp = func_calls[[fname]],
                                f_eval = eval_cat_noncat,
                                use_addargs = TRUE, addargs = addargs)
    Res[[fname]] = r
    print(paste('finish', fname))
  }
  Res
}

test_data_openml_bunch <- function(func_calls, dnames = NULL, ...){
  if (is.null(dnames)) dnames = c("abalone", "heart", "cmc", "creditg", "credita", "colic")
  r = list()
  for (dname in dnames){
    r[[dname]] = test_data_openml(func_calls, dname = dname, ...)
  }
  print(paste0('finish ', dname))
  r
}

run_data_openml <- function(fimp, fname, dnames = NULL, mfrac = 0.2,  mask_type = 'MCAR', nrep=10){
  if (is.null(dnames)) dnames = c("abalone", "heart", "cmc", "creditg", "credita", "colic")
  Res = list()
  for (dname in dnames){
    if (fname %in% c('imputeFAMD', 'softImpute')) val_ratio = 0.2 else val_ratio = 0
    if (fname == 'EGC_our') addargs = c('mu_loc', 'seed_num') else addargs = NULL
    fdata = get_data_openml(dname, mask_fraction = mfrac, mask_type = mask_type, val_ratio = val_ratio)
    r = imputation_pipeline(nrep=nrep,
                                f_data = fdata, f_imp = fimp,
                                f_eval = eval_cat_noncat,
                                use_addargs = TRUE, addargs = addargs)
    print(paste('finish', dname))
    Res[[dname]] = r
  }
  Res
}

run_data_openml_bunch <- function(func_calls, dnames=NULL, ...){
  r = list()
  for (m in names(func_calls)){
    r[[m]] = run_data_openml(func_calls[[m]], m, dnames = dnames, ...)
  }
  print(paste0('finish ', m))
  r
}

print_summary_openml <- function(r_data, f=colMeans){
  r = map(r_data, ~ map(.x, f))
  map(r, as.data.frame)
}

collect_res_perdata <- function(res_list, name){
  Res = list()
  r = map(res_list, ~ .x[[name]])
  tunings = c('imputeFAMD', 'softImpute')
  for (m in tunings){
    if (m %in% names(r)){
      r[[paste0(m, '_bestall')]] = as.data.frame(r[[m]]$bestall)
      r[[paste0(m, '_bestcat')]] = as.data.frame(r[[m]]$bestcat)
      r[[m]] = NULL
    }
  }
  r
}

openml_raw_to_perdata <- function(res){
  r_summary = list()
  methods = names(res)
  tunings = c('imputeFAMD', 'softImpute')
  for (m in setdiff(methods, tunings)){
    r_summary[[m]] = map(res[[m]], summarize_single_fcall)
  }
  for (m in tunings){
    if (m %in% methods){
      rd = map(res[[m]], summarize_mlutiple_fcall)
      r_summary[[m]] =  map(rd, ~ test_at_best(.x))
    }
  }
  dnames = names(r_summary[[1]])
  r_perdata = list()
  for (m in dnames){
    r_perdata[[m]] = collect_res_perdata(r_summary, m)
  }
  r_perdata
  #r_perdata = collect_res_perdata(r_summary)
}
