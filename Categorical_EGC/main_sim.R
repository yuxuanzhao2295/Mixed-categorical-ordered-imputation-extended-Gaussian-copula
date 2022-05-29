source('~/gcimputeR/Development_tests/nominal_benchmark/func_experiment_sim.R')

func_calls  = list(EGC_our = call_gcimpute_nominal(read_mu = TRUE),
                   missForest = call_missForest(),
                   MICE = call_mice(m=5),
                   imputeFAMD = call_imputeFAMD(rank = c(1,3,5,7,9)),
                   softImpute = call_softimpute_nominal(),
                   baseline = call_baseline())


r = test_sim_mixed(func_calls = func_calls['EGC_our'])
r = run_sim_mixed(func_calls, nrep = 2, K = 6)
rp = sim_postprocessing(r)
rpdf = to_DF(rp)
r = run_sim_mixed(func_calls[c('softImpute')], nrep = 10, K = 6)
rp = sim_postprocessing(r)
to_DF(rp)

names_run = setdiff(names(func_calls), 'softImpute')
# to run
res_sim_main_v1 = run_sim_mixed(func_calls, nrep = 10, K = 6, mask_fraction = 0.3)
rbackup = res_sim_main_v1
res_sim_main_v1_df = to_DF(sim_postprocessing(res_sim_main_v1))
save(res_sim_main_v1_df,
     file = '~/gcimputeR/Development_tests/nominal_benchmark/res_sim_main_v1_df.RData')

# variants
func_calls  = list(EGC_our = call_gcimpute_nominal(read_mu = FALSE),
                   missForest = call_missForest(),
                   MICE = call_mice(m=5),
                   imputeFAMD = call_imputeFAMD(rank = c(1,3,5,7,9)),
                   softImpute = call_softimpute_nominal(),
                   baseline = call_baseline())

sim_K3 = run_sim_mixed(func_calls, nrep = 10, K = 3, mask_fraction = 0.3)
sim_K9 = run_sim_mixed(func_calls, nrep = 10, K = 9, mask_fraction = 0.3)
sim_m20 = run_sim_mixed(func_calls, nrep = 10, K = 9, mask_fraction = 0.2)
sim_m40 = run_sim_mixed(func_calls, nrep = 10, K = 9, mask_fraction = 0.4)
sim_K3_df = to_DF(sim_postprocessing(sim_K3))
sim_K9_df = to_DF(sim_postprocessing(sim_K9))
sim_m20_df = to_DF(sim_postprocessing(sim_m20))
sim_m40_df = to_DF(sim_postprocessing(sim_m40))
save(sim_K3_df,
     file = '~/gcimputeR/Development_tests/nominal_benchmark/sim_K3_df.RData')
save(sim_K9_df,
     file = '~/gcimputeR/Development_tests/nominal_benchmark/sim_K9_df.RData')
save(sim_m20_df,
     file = '~/gcimputeR/Development_tests/nominal_benchmark/sim_m20_df.RData')
save(sim_m40_df,
     file = '~/gcimputeR/Development_tests/nominal_benchmark/sim_m40_df.RData')
