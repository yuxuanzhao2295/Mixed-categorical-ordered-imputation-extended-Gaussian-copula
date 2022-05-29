library(purrr)
load("data/processed_datas_openml.RData")
source('func_load_data.R')

col_levels <- function(X) apply(X, 2, function(x)length(unique(x[!is.na(x)])))


encode_bincols <- function(X){
  num_col_levels = col_levels(X)
  bin_cols = which(num_col_levels == 2)
  bin
  for (j in bin_cols){
    xrelabel = cat_to_integers(X[,j])
    cat_labels[[as.character(j)]] = xrelabel$xlevels
    X[,j] = xrelabel$x
  }
}

relabel_df <- function(X, cols, reduce=FALSE){
  old_labels = list()
  for (j in cols){
    if (reduce & is.factor(X[,j])) xj = droplevels(X[,j]) else xj = X[,j]
    xrelabel = cat_to_integers(xj)
    old_labels[[as.character(j)]] = xrelabel$xlevels
    X[,j] = xrelabel$x
  }
  list(X = X, label = old_labels)
}

process_ordcol <- function(x, orders, NA_levels=NULL){
  x = as.factor(x)
  for (m in NA_levels){
    x[x==m] = NA
  }
  x = droplevels(x)
  x = ordered(x, levels = orders)
  #list(x = gcimputeR::factor_to_num(x), orders = orders)
  as.integer(x)
}

preprocessing_openml <- function(X,
                                 drop_cols = NULL,
                                 cat_index = NULL, ords_info = NULL,
                                 drop_last=TRUE, to_numeric = TRUE){
  p = ncol(X)
  y = X[,p]
  if (drop_last) X = X[,-p]
  if (!is.null(drop_cols)){
    X = X[,-drop_cols]
    cat_index = index_after_remove(cat_index,
                                   index_r = drop_cols)
  }
  #

  cat_relable = relabel_df(X, cat_index)
  X = cat_relable$X
  #
  bin_index = setdiff(which(col_levels(X) <= 2), cat_index)
  bin_relabel = relabel_df(X, bin_index, reduce=TRUE)
  X = bin_relabel$X
  #
  ord_label = list()
  for (loc in names(ords_info)){
    j = index_after_remove(as.integer(loc), index_r = drop_cols)
    X[,j] = process_ordcol(X[,j],
                           orders = ords_info[[loc]]$orders,
                           NA_levels = ords_info[[loc]]$NA_levels)
    ord_label[[loc]] = ords_info[[loc]]$orders
  }
  if (to_numeric) Xn = gcimputeR::to_numeric_matrix(X) else Xn = X
  ord_label = c(bin_relabel$label, ord_label)
  list(X = Xn, cat_label = cat_relable$label, ord_label = ord_label, label = y)
}

index_after_remove <- function(index, index_r){
  index = map_dbl(index, ~ .x - sum(.x > index_r))
  as.integer(index)
}

preprocessing_info <- function(){
  cat_index_openml = list()
  ord_info = list()
  drop_cols = list()
  missing_char = list()

  #
  cat_index_openml[['heart']] = c(3, 7, 13)
  ord_info[['heart']] = NULL
  drop_cols[['heart']] = NULL

  #
  cat_index_openml[['creditg']] = c(3, 4, 9, 10, 12, 14, 15, 17)
  ord_info[['creditg']] = list('1' = list(orders = c("<0","0<=X<200", ">=200"),
                                          NA_levels = c("no checking")),
                               '6' = list(orders = c("<100","100<=X<500", "500<=X<1000", ">=1000"),
                                          NA_levels = c("no known savings")),
                               '7'= list(orders = c("unemployed","<1", "1<=X<4", "4<=X<7", ">=7"),
                                         NA_levels = NULL))
  drop_cols[['creditg']] = NULL

  # cmc
  cat_index_openml[['cmc']] = c(7)
  ord_info[['cmc']] = NULL
  drop_cols[['cmc']] = NULL

  # abalone
  cat_index_openml[['abalone']] = c(1)
  ord_info[['abalone']] = NULL
  drop_cols[['abalone']] = NULL

  # credit_a
  cat_index_openml[['credita']] = c(5, 6, 7, 13)
  ord_info[['credita']] = NULL
  drop_cols[['credita']] = c(4)


  name = 'colic'
  cat_index_openml[[name]] = c(9, 18, 21, 23)
  ord_info[[name]] = list('7' = list(orders = c("2", "1", "3", "4"), NA_levels = NULL),
                          '8' = list(orders = c("2", "1", "3", "4"), NA_levels = NULL),
                          '10' = list(orders = c("2", "1"), NA_levels = c("3")),
                          '11' = list(orders = c("1", "3", "4", "5"), NA_levels = c("2")),
                          '15' = list(orders = c("1", "3", "2"), NA_levels = NULL),
                          '17' = list(orders = c("2", "1", "3", "4"), NA_levels = NULL))
  drop_cols[[name]] = c(3, 25, 26, 27)

  list(cat_index_openml = cat_index_openml, ord_info  = ord_info, drop_cols = drop_cols)
}

set_NA_col <- function(x, word){
  x[x==word] = NA
  x
}

count_var_type <- function(data){
  ncat = length(data$cat_label)
  nord = length(data$ord_label)
  n = nrow(data$X)
  p = ncol(data$X)
  cat_K = map_int(data$cat_label, length)
  unlist(list(n=n, p=p,
              pcat = ncat,pord = nord, pcont = p - ncat - nord,
              Kcat = sum(cat_K),
              Kmax = max(cat_K)))
}


get_data_openml <- function(name, mask_fraction=0.2, mask_type = 'MCAR', val_ratio=0){
  dat = processed_datas_openml[[name]]
  if (mask_type == 'MCAR'){
    mis_suffix = paste0('_m', as.integer(mask_fraction*100))
  }else if (mask_type == 'MNAR'){
    mis_suffix = '_MNAR'
  }else if (mask_type == 'MAR'){
    mis_suffix = '_MAR'
  }
  function(seed){
    r = data_masker_general(data = dat$X, cat_index = as.integer(names(dat$cat_label)),
                            seed = seed,
                            mask_fraction = mask_fraction, mask_type = mask_type,
                            val_ratio = val_ratio,
                            label = dat$label)
    r[['seed_num']] = seed
    prefix = '/WAIT_FOR_COMPLETE_PYTHON_VERSION_TO_AVIOD_CALLING_PYTHON_IN_R/'
    r[['mu_loc']] = paste0(prefix, name, mis_suffix, '.csv')
    r
  }
}


