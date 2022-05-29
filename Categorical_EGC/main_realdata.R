source('func_experiments_realdata.R')

func_calls  = list(EGC_our = call_gcimpute_nominal(),
                   missForest = call_missForest(),
                   MICE = call_mice(m=5),
                   imputeFAMD = call_imputeFAMD(rank = c(1,3,5,7,9)),
                   softImpute = call_softimpute_nominal())

realdata_report = run_data_openml_bunch(func_calls, nrep = 10, mfrac = 0.2)
realdata_report_df = to_DF(openml_raw_to_perdata(realdata_report))

# variants
realdata_m10 = run_data_openml_bunch(func_calls, nrep = 10, mfrac = 0.1)
realdata_m10_df = to_DF(openml_raw_to_perdata(realdata_m10))

realdata_m30 = run_data_openml_bunch(func_calls, nrep = 10, mfrac = 0.3)
realdata_m30_df = to_DF(openml_raw_to_perdata(realdata_m30))

realdata_mar = run_data_openml_bunch(func_calls, nrep = 10, mfrac = 0.2, mask_type = 'MAR')
realdata_mar_df = to_DF(openml_raw_to_perdata(realdata_mar))

realdata_mnar = run_data_openml_bunch(func_calls, nrep = 10, mfrac = 0.2, mask_type = 'MNAR')
realdata_mnar_df = to_DF(openml_raw_to_perdata(realdata_mnar))
