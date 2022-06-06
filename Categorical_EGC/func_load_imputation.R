# setwd('~/Documents/GitHub/Mixed-categorical-ordered-imputation-extended-Gaussian-copula')
source('helper_func_call.R')
source('func_EGC.R')

to_nearest_ord <- function(x, ords=NULL, xobs=NULL){
  if (is.null(ords)) ords = unique(xobs[!is.na(xobs)])
  ords[which.min(abs(x - ords))]
}

#  install_github("udellgroup/gcimputeR")

# Ordinal versus Nominal? NO
call_missForest <- function(n_level=20, verbose=FALSE){
  force(verbose)
  function(X, cat_index=NULL, ...){
    Xdf = to_df_with_factor(X, n_level = n_level, cat_index = cat_index)
    X = Xdf$X
    if (verbose) f = missForest::missForest else f = purrr::quietly(missForest::missForest)
    est = f(X)
    if (!verbose){
      Ximp = gcimputeR::to_numeric_matrix(est$result$ximp)
      if (length(est$warnings) > 0) cat('Warning meesage: ', est$warnings)
    }else{
      Ximp = gcimputeR::to_numeric_matrix(est$ximp)
    }

    list(Ximp = Ximp, est = est$result, cat_index = cat_index, ord_index = Xdf$ord_cols)
  }
}



# Ordinal versus Nominal? NO
call_imputeFAMD <- function(rank, n_level=20){
  if (length(rank) == 1){
    function(X, cat_index = NULL, ...){
      Xdf = to_df_with_factor(X, n_level = n_level, cat_index = cat_index)
      X = Xdf$X
      res = safely(missMDA::imputeFAMD)(X, ncp = rank)
      Ximp = post_processing_FAMD(res, nfactor = Xdf$nfactor, n_level = n_level)
      list(Ximp = Ximp, est = est, cat_index = cat_index)
    }
  }else{
    function(X, cat_index = NULL, ...){
      Xdf = to_df_with_factor(X, n_level = n_level, cat_index = cat_index)
      X = Xdf$X
      ests = purrr::map(rank, ~ safely(missMDA::imputeFAMD)(X, ncp = .x))
      Ximps = map(ests, ~ post_processing_FAMD(.x, nfactor = Xdf$nfactor, n_level = n_level))
      names(Ximps) = paste0('rank', rank)
      Ximps
    }
  }
}

call_softimpute_nominal <- function(grid_len=10, n_level=20, row_scale = FALSE, verbose=FALSE,
                            min_lam_ratio=0.1, to_integer = TRUE){
  force(grid_len)
  function(X, cat_index){
    Xinput = X
    onehot = one_hot_matrix(to_numeric_matrix(X), cat_index)
    X = onehot$X
    #
    xc = softImpute::biScale(X, row.scale = row_scale)
    lam0 = softImpute::lambda0(xc) * 0.99
    if (is.null(min_lam_ratio)) min_lam = 1 else min_lam = lam0 * min_lam_ratio
    lamseq=exp(seq(from=log(lam0),to=log(min_lam),length=grid_len))
    Ximps = vector('list', length = grid_len)
    ranks = numeric(grid_len)
    warm = NULL
    p = ncol(X)
    for(i in seq_along(lamseq)){
      fiti = softImpute::softImpute(xc, lambda=lamseq[i], rank.max=p-1, warm.start=warm, type = 'svd')
      ranks[i]=sum(round(fiti$d,4)>0)
      warm = fiti
      Ximps[[i]] = softImpute::complete(xc, fiti, unscale = TRUE)
      if (verbose) cat(i,"lambda=",lamseq[i],"rank",ranks[i],"\n")
    }
    #return(list(Ximps = Ximps, onehot = onehot))
    Ximps = purrr::map(Ximps,  ~ inverse_one_hot_matrix(.x,
                                                        cat_index_list = onehot$cat_index_list,
                                                        cat_labels = onehot$cat_labels))
    #
    if (to_integer) Ximps = purrr::map(Ximps,  ~ trunc_to_integer(Xinput, Ximp = .x, n_level=n_level, cat_index =cat_index))
    names(Ximps) = paste0('rank', ranks)
    Ximps
  }
}





# Ordinal versus Nominal? YES
call_gcimpute_nominal <- function(verbose = FALSE,
                                  read_mu = FALSE,
                                  oracle_mu= FALSE, oracle_cor = FALSE,
                                  trunc_method = 'Iterative', old=FALSE,
                                  n_MI = 0, fast_MI=TRUE, ...){
  force(verbose)
  force(old)
  force(trunc_method)
  force(oracle_mu)
  force(oracle_cor)
  function(X, cat_index, mu = NULL, corr = NULL, mu_loc = NULL, seed_num = NULL){
    if (oracle_mu) mu_use = mu else mu_use = NULL
    if (read_mu){
      stopifnot(!is.null(mu_loc) | !is.null(seed_num) )
      readed = as.matrix(read.csv(mu_loc, row.names = 1))
      mu_use = readed[seed_num,-1]
      m = is.na(mu_use)
      if (any(m)) mu_use = mu_use[!m]
      time_add = readed[seed_num,1]
    }else{
      time_add = 0
    }
    if (oracle_cor) corr_use = corr else corr_use = NULL
    if (verbose){
      cat('is mu provided?', !is.null(mu_use), '\n')
      cat('is corr provided?', !is.null(corr_use), '\n')
    }
    est = impute_nominal_cont_gc(X = X,
                                 cat_index = cat_index,
                                 mu = mu_use,
                                 corr = corr_use,
                                 verbose = verbose,
                                 trunc_method = trunc_method,
                                 n_MI = n_MI, fast_MI = fast_MI,
                                 old = old, ...)
    list(Ximp = est$Ximp, Ximp_MI = est$Ximp_MI, cat_index = cat_index, est = est, time_add = time_add)
  }
}


call_baseline <- function(n_level = 20){
  function(X, cat_index=NULL, ...){
    Xdf = to_df_with_factor(X, n_level = n_level, cat_index=cat_index)
    X = Xdf$X
    cat_index = Xdf$cat_cols
    ord_index = Xdf$ord_cols
    Ximp = X
    p = ncol(X)
    for (j in 1:p){
      m = is.na(X[,j])
      xobs = as.numeric(X[!m,j])
      if (j %in% cat_index){
        Ximp[m,j] = majority_vote(xobs)
      }else if (j %in% ord_index){
        Ximp[m,j] = to_nearest_ord(median(xobs), ords = unique(xobs))
      }else{
        Ximp[m,j] = median(xobs)
      }
    }
    list(Ximp = Ximp, cat_index =  Xdf$cat_cols)
  }
}


# Ordinal versus Nominal? YES
call_mice <- function(m = 5, n_level = 20, verbose = FALSE){
  force(m)
  force(n_level)
  force(verbose)
  function(X, cat_index=NULL, ...){
    Xdf = to_df_with_factor(X, cat_index = cat_index, n_level = n_level)
    X = Xdf$X
    res = mice::mice(X, m = m, printFlag = verbose)
    Ximp_MI = map(mice::complete(res, 'all'), to_numeric_matrix)
    f = call_baseline(n_level)
    Ximp_MI = map_if(Ximp_MI, ~ any(is.na(.x)), ~ f(.x, cat_index)$Ximp)
    ximp = MI_to_single(Ximp_MI, X, n_level)
    #ximp = to_numeric_matrix(mice::complete(res))
    #if (any(is.na(ximp))) ximp = call_baseline(n_level)(ximp, cat_index)$Ximp
    list(Ximp = ximp, Ximp_MI = Ximp_MI, cat_index = cat_index, res = res)
  }
}


