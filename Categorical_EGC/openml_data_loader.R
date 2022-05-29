library(mlr)
library(OpenML)
library(purrr)
source('func_preprocessing_openml.R')

# colic is purely ordinal
data_ids = list(heart = 53, creditg = 31, cmc = 23, abalone = 1557,credita = 29,colic = 25)
#
datas = list()
for (name in names(data_ids)) datas[[name]] = getOMLDataSet(data.id = data_ids[[name]])
raw_datas_openml = datas

out = preprocessing_info()
processed_datas_openml = list()
for (name in names(raw_datas_openml)){
  processed_datas_openml[[name]] = preprocessing_openml(raw_datas_openml[[name]]$data,
                                                        drop_cols = out$drop_cols[[name]],
                                                        cat_index = out$cat_index_openml[[name]],
                                                        ords_info = out$ord_info[[name]])
}

save(processed_datas_openml, file = 'data/processed_datas_openml.RData')


t(as.data.frame(map(processed_datas_openml, count_var_type)))
"
           n  p pcat pord pcont Kcat Kmax
heart    270 13    3    3     7   10    4
creditg 1000 20    8    6     6   36   10
cmc     1473  9    1    3     5    4    4
abalone 4177  8    1    0     7    3    3
credita  690 14    4    4     6   29   14
colic    368 23    4    9    10   17    6
"






