source('func_experiment_sim.R')

func_calls  = list(EGC_our = call_gcimpute_nominal(),
                   missForest = call_missForest(),
                   MICE = call_mice(m=5),
                   imputeFAMD = call_imputeFAMD(rank = c(1,3,5,7,9)),
                   softImpute = call_softimpute_nominal(),
                   baseline = call_baseline())


#
sim_report = run_sim_mixed(func_calls, nrep = 10, K = 6, mask_fraction = 0.3)
sim_report_df = to_DF(sim_postprocessing(sim_report))

sim_K3 = run_sim_mixed(func_calls, nrep = 10, K = 3, mask_fraction = 0.3)
sim_K9 = run_sim_mixed(func_calls, nrep = 10, K = 9, mask_fraction = 0.3)
sim_m20 = run_sim_mixed(func_calls, nrep = 10, K = 6, mask_fraction = 0.2)
sim_m40 = run_sim_mixed(func_calls, nrep = 10, K = 6, mask_fraction = 0.4)
sim_K3_df = to_DF(sim_postprocessing(sim_K3))
sim_K9_df = to_DF(sim_postprocessing(sim_K9))
sim_m20_df = to_DF(sim_postprocessing(sim_m20))
sim_m40_df = to_DF(sim_postprocessing(sim_m40))

