library(purrr)
merge_list <- function(s1, s2) purrr::map2(s1, s2, c)


chatty <- function(f){
  force(f)

  function(x, ...){
    cat('Processing rep ', x, "\n", sep = "")
    res <- f(x, ...)
    res
  }
}

aggregate_reps <- function(nrep, f_rep, verbose=TRUE, seed_start=0, safe_call=TRUE, ...){
  "  out_list = vector('list', nrep)
  for (i in 1:nrep){
    out = safe_f_rep(seed = i, ...)
    if (!is.null(to_save)) to_save[[length(to_save)+1]] = out
    out_list[[i]] = out
    if (verbose) print(paste('finish rep', i))
  }"
  if (verbose) f_rep = chatty(f_rep)
  seeds = (1:nrep) + seed_start
  if (safe_call) f_rep = purrr::safely(f_rep)
  out_list = purrr::map(seeds, f_rep, ...)

  if (safe_call){
    out_list = purrr::transpose(out_list)
    success_rep = purrr::map_lgl(out_list$error, is.null)
    if (all(success_rep)){
      out_list$result
    }else{
      r = list()
      r[['successed']] = out_list$result[success_rep]
      r[['fail_rep']] = which(!success_rep)
      r
    }
  }
  else{
    stop('To be implemented')
  }
}




imputation_pipeline_rep <- function(seed, f_data, f_imp, f_eval, use_addargs = TRUE, addargs=NULL){
  data = f_data(seed = seed)
  xtrain = data[['train']]
  xtest = data[['test']]
  xval = data[['validation']]
  stopifnot(!is.null(xtrain) & !is.null(xtest))
  add_args = setdiff(names(data), c('train', 'test', 'validation'))

  s1 = Sys.time()
  if (use_addargs & length(add_args)>0){
    #call_input = data[add_args]
    call_input = list()
    call_input[['X']] = xtrain
    call_input[['cat_index']] = data[['cat_index']]
    stopifnot(!is.null(data$cat_index))
    if (!is.null(addargs)){
      for (r in addargs) call_input[[r]] = data[[r]]
    }
    #print(names(call_input))
    out = do.call(f_imp, call_input)
  }else{
    out = f_imp(xtrain)
  }
  s2 = Sys.time()

  eval_test_and_val <- function(est){
    if (is.null(est$Ximp)) est$Ximp = xtrain # NA will appear in the imputation error in this situation
    if (!is.null(xval)){
      r_test = f_eval(est, xtrue = xtest, xobs = xval,  prefix = 'test_')
      r_val = f_eval(est, xtrue = xval, xobs = xtrain,  prefix = 'val_')
      r = c(r_test, r_val)
    }else{
      r = f_eval(est, xtrue = xtest, xobs = xtrain, prefix = 'test_')
    }
    r
  }

  # if evaluation needs additional arguments (that is not returned by out), pass by appending out
  append_out <- function(xlist){
    if (length(add_args)>0){
      for (name in add_args){
        xlist[[name]] = data[[name]]
      }
    }
    xlist
  }

  r_time  = as.numeric(s2-s1, units = 'secs')
  if (!is.null(out$time_add)) r_time = r_time + out$time_add
  if (!is.null(out$Ximp)){
    r = eval_test_and_val(append_out(out))
    r[['time']] = r_time
  }else{
    r = purrr::map(out, ~ eval_test_and_val(append_out(list(Ximp = .x))))
    r = purrr::transpose(r)
    r = purrr::map(r, unlist)
    r[['time']] = rep(r_time, length(out))
    names(r[['time']]) = names(out)
  }

  r
}

# f_rep is constructed from imputation_pipeline_rep, and executed and aggregated using aggregate_reps.
# During its executition, f_rep is replaced with a silent version
imputation_pipeline <- function(nrep, f_data, f_imp, f_eval, verbose=TRUE, seed_start = 0,
                                use_addargs = FALSE, addargs = NULL){
  f_rep <- function(seed) imputation_pipeline_rep(seed,
                                                  f_data = f_data, f_imp = f_imp, f_eval = f_eval,
                                                  use_addargs = use_addargs, addargs = addargs)
  out = aggregate_reps(nrep = nrep, f_rep = f_rep, verbose = verbose, seed_start = seed_start)
  out
}



summarize_single_fcall <- function(res) as.data.frame(purrr::reduce(res, merge_list))
summarize_mlutiple_fcall <- function(res){
  r = purrr::transpose(res)
  r2 = purrr::map(r, ~ purrr::reduce(.x, rbind))
  r2
}


get_mean_sd <- function(df){
  list(mean = colMeans(df), std = apply(df, 2, sd))
}

to_DF_each <- function(s){
  ss = map(s, as.matrix)
  sss = map(ss, get_mean_sd)
  as.data.frame(sss)
}

to_DF <- function(Res_list){
  s = map(Res_list, to_DF_each)
  t(as.data.frame(s))
}

add_maskby <- function(X, cat_index){
  BY = X
  for (i in cat_index){
    BY[,i] = label_to_freqrank_cat(X[,i])
  }
  BY
}

mask_by_misprob <- function(X, P, seed = NULL){
  if (!is.null(seed)) set.seed(seed)
  ifmissing = runif(prod(dim(P))) < P
  indicator = array(1, dim = dim(P))
  indicator[ifmissing] = NA
  X * indicator
}

mask_MNAR_three_piece <- function(X, BY = NULL, ...){
  if (is.null(BY)) BY = X
  p = ncol(X)
  for (j in 1:p){
    X[,j] = mask_three_piece(X[,j], BY[,j], ...)
  }
  X
}

misprob_selfmask_three_piece <- function(X, BY = NULL, ...){
  if (is.null(BY)) BY = X
  p = ncol(X)
  for (j in 1:p){
    X[,j] = misprob_three_piece(X[,j], BY[,j], ...)
  }
  X
}

adjust_misprob_MAR_row <- function(misprob, pobs_ratio = 1/3){
  p = length(misprob)
  pobs = as.integer(p*pobs_ratio)
  pobs_loc = sample(1:p, pobs)
  misprob_new = misprob
  misprob_new[pobs_loc] = 0
  pmis_from = sample(pobs_loc, p-pobs, replace = TRUE)
  misprob_new[-pobs_loc] = misprob[pmis_from]
  misprob_new
}

adjust_misprob_MAR_matrix <- function(misprob, pobs_ratio = 1/3){
  n = nrow(misprob)
  p = ncol(misprob)
  misprob_new = matrix(0, n, p)
  for (i in 1:n){
    misprob_new[i,] = adjust_misprob_MAR_row(misprob[i,], pobs_ratio = pobs_ratio)
  }
  misprob_new
}

mask_three_piece <- function(x, by = NULL,
                             low_m = NULL, high_m = NULL, normal_m = 0.3,
                             low_q = 1/3, high_q = 2/3){
  if (is.null(low_m)) low_m = normal_m - 0.1
  if (is.null(high_m)) high_m = normal_m * 2 - low_m
  if (is.null(by)) by = x
  mask_on = !is.na(by)
  by_mask = by[mask_on]
  x_mask = x[mask_on]
  cuts = quantile(by_mask, c(low_q,high_q))
  index_l = by_mask <= cuts[1]
  index_n = (by_mask > cuts[1]) & (by_mask <= cuts[2])
  index_h = by_mask > cuts[2]
  if (any(index_l)) x_mask[index_l] = mask_MCAR_vec(x_mask[index_l], high_m)
  if (any(index_n)) x_mask[index_n] = mask_MCAR_vec(x_mask[index_n], normal_m)
  if (any(index_h)) x_mask[index_h] = mask_MCAR_vec(x_mask[index_h], low_m)
  x[mask_on] = x_mask
  x
}

misprob_three_piece <- function(x, by = NULL,
                                normal_m = 0.3,
                                low_m = NULL, high_m = NULL,
                                low_q = 1/3, high_q = 2/3){
  if (is.null(low_m)) low_m = normal_m - 0.1
  if (is.null(high_m)) high_m = normal_m * 2 - low_m
  if (is.null(by)) by = x
  mask_on = !is.na(by)
  by_mask = by[mask_on]
  x_mask = x[mask_on]
  cuts = quantile(by_mask, c(low_q,high_q))
  index_l = by_mask <= cuts[1]
  index_n = (by_mask > cuts[1]) & (by_mask <= cuts[2])
  index_h = by_mask > cuts[2]
  misprob = x_mask
  if (any(index_l)) misprob[index_l] = high_m
  if (any(index_n)) misprob[index_n] = normal_m
  if (any(index_h)) misprob[index_h] = low_m
  Misprob = numeric(length(x))
  Misprob[mask_on] = misprob
  Misprob
}

label_to_freqrank_cat <- function(x){
  xobs = x[!is.na(x)]
  freq_sort = names(sort(table(xobs)))
  nl = length(freq_sort)
  xnew = xobs
  for (i in 1:nl){
    label = as.integer(freq_sort[i])
    xnew[xobs == label] = i
  }
  x[!is.na(x)] = xnew
  x
}


