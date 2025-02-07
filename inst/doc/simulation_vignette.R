## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----simulate-data, eval=FALSE------------------------------------------------
#  # Suppose we want to generate a single dataset with K=3, r=4, N=50, M=8:
#  sim_data <- simulate_mixture_data(
#    N=50,
#    K=3,
#    r=4,
#    M=8,
#    T_each=100,
#    cluster_sep=1.0,
#    noise_scale=1.0,
#    seed=123
#  )
#  
#  # sim_data is a list with:
#  #  $list_of_data : 50 matrices (each 100 x 8)
#  #  $z_true       : integer vector (50) with cluster assignments

## ----fit-models, eval=FALSE---------------------------------------------------
#  # Fit MFA with some guess, e.g. K=3, r=4
#  fit_mfa <- mfa_em_fit(
#    list_of_data = sim_data$list_of_data,
#    K = 3,
#    r = 4,
#    max_iter = 50,
#    tol = 1e-4
#  )
#  
#  # fit_mfa is a list with: z, pi, mu, Lambda, Psi, logLik, resp, etc.
#  
#  # Fit Mixture PCA similarly
#  fit_pca <- mixture_pca_em_fit(
#    list_of_data = sim_data$list_of_data,
#    K = 3,
#    r = 4,
#    max_iter = 50,
#    tol = 1e-4
#  )
#  
#  # For each fit, you can compute BIC:
#  num_rows <- sum(sapply(sim_data$list_of_data, nrow))
#  bic_val  <- compute_bic(fit_mfa$logLik, K=3, r=4, M=8, N_total_rows=num_rows)
#  cat("BIC for MFA = ", bic_val, "\n")

## ----parallel-sim, eval=FALSE-------------------------------------------------
#  df_result <- simulate_and_run_parallel_save(
#    K_true_vec = c(3,4),
#    r_true_vec = c(3,5),
#    N_vec      = c(50),
#    M_vec      = c(8,12),
#    T_each     = 100,
#    sep_vec    = c(0.5,1.0),
#    noise_vec  = c(1.0,2.0),
#    max_iter   = 20,
#    n_rep_init = 3,
#    seed_data_base = 999,
#    mc.cores = 2,       # number of parallel cores
#    output_dir = "results"
#  )
#  
#  # This returns a data frame summarizing each condition's best model results,
#  # and also saves .rds files like "result_K3_r3_N50_M8_sep0.5_noise1.0.rds"
#  # in the "results/" directory.

## ----run-main, eval=FALSE-----------------------------------------------------
#  df_out <- run_main_experiment_parallel_save()
#  head(df_out)

## ----read-and-visualize, eval=FALSE-------------------------------------------
#  library(ggplot2)
#  
#  # read_and_visualize_all() will:
#  #  1) read all .rds files from "results/"
#  #  2) produce Table 1 (static conditions)
#  #  3) produce Table 2 (ari / BIC summary)
#  #  4) create Figures 1..2..3,
#  #  5) (optionally) Figure 4 if you have a fitted object
#  
#  df_all <- read_and_visualize_all("results")
#  head(df_all)

