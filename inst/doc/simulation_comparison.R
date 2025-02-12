## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----simulate-data, eval=FALSE------------------------------------------------
#  library(mixMuscleSynergy)  # or your package name
#  
#  set.seed(123)
#  N <- 50
#  K_true <- 2
#  r_true <- 2
#  M <- 6
#  T_each <- 100
#  
#  sim_data <- simulate_mixture_data_for_comparison(
#    N=N, K_true=K_true, r_true=r_true, M=M, T_each=T_each, seed=123
#  )
#  
#  list_of_data <- sim_data$list_of_data
#  z_true       <- sim_data$z_true
#  
#  cat("True subgroup sizes:\n")
#  print(table(z_true))

## ----compare-methods, eval=FALSE----------------------------------------------
#  df_comp <- simulate_and_compare_methods_mfa_mpca_withSSE(
#    N=50, K_true=2, r_true=2, M=6, T_each=100, seed=123,
#    r_singleFA=2,
#    r_singlePCA=2,
#    r_subFA=1, r_clusterFA=2,
#    r_subPCA=1, r_clusterPCA=2
#  )
#  
#  df_comp

