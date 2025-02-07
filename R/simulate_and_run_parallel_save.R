#' Run Simulation (MFA & MPCA) in Parallel, Saving Each Condition Separately
#'
#' This function enumerates a grid of \code{(K_true, r_true, N, M, sep, noise)} values,
#' simulates data, and tries a range of \code{(K, r)} with multiple initial seeds
#' for both MFA and Mixture PCA. It selects the best model by BIC for each method,
#' computes ARI, and saves the result to an \code{.rds} file per condition.
#'
#' @param K_true_vec Vector of true K values to simulate.
#' @param r_true_vec Vector of true r values to simulate.
#' @param N_vec Vector of subject counts.
#' @param M_vec Vector of channel counts.
#' @param T_each Time-series length.
#' @param sep_vec Vector of separation parameters.
#' @param noise_vec Vector of noise scale parameters.
#' @param max_iter Max iteration for EM.
#' @param n_rep_init Number of random initial seeds for each (K,r) in search.
#' @param seed_data_base Base seed for data simulation.
#' @param mc.cores Number of parallel cores to use.
#' @param output_dir Directory to save each result as \code{result_K{...}_r{...}...rds}.
#'
#' @details
#' For each condition, the function calls \code{\link{simulate_mixture_data}} to create data,
#' then loops over candidate K, r in \code{1:5} (hard-coded), fits MFA and MPCA multiple times,
#' picks the best BIC, and computes ARI relative to the true cluster assignment.
#'
#' Finally, it saves a \code{data.frame} of length 1 (or 2) in an \code{.rds} file for that condition,
#' and returns a merged data frame of all conditions.
#'
#' @return A data frame with columns:
#' \item{K_true, r_true, N, M, sep, noise}{True simulation parameters.}
#' \item{ARI_MFA, BIC_MFA, Khat_MFA, rhat_MFA, ARI_PCA, BIC_PCA, Khat_PCA, rhat_PCA}{Results.}
#'
#' @export
simulate_and_run_parallel_save <- function(
    K_true_vec = c(3,4),
    r_true_vec = c(3,5),
    N_vec      = c(50),
    M_vec      = c(8,12),
    T_each     = 100,
    sep_vec    = c(0.5,1.0),
    noise_vec  = c(1.0,2.0),
    max_iter   = 30,
    n_rep_init = 5,
    seed_data_base = 2023,
    mc.cores   = parallel::detectCores() - 1,
    output_dir = "results"
){
  if(!dir.exists(output_dir)){
    dir.create(output_dir, recursive=TRUE)
  }
  
  param_grid <- expand.grid(
    K_true = K_true_vec,
    r_true = r_true_vec,
    N      = N_vec,
    M      = M_vec,
    sep    = sep_vec,
    noise  = noise_vec,
    stringsAsFactors=FALSE
  )
  
  run_one <- function(idx){
    row_i <- param_grid[idx, ]
    Kt <- row_i$K_true
    rt <- row_i$r_true
    Nv <- row_i$N
    Mv <- row_i$M
    sv <- row_i$sep
    nv <- row_i$noise
    
    data_seed <- seed_data_base + idx
    
    # (1) シミュレーション
    sim_data <- simulate_mixture_data(
      N=Nv, K=Kt, r=rt, M=Mv, T_each=T_each,
      cluster_sep=sv, noise_scale=nv,
      seed=data_seed
    )
    list_of_data <- sim_data$list_of_data
    z_true       <- sim_data$z_true
    
    total_rows <- sum(sapply(list_of_data, nrow))
    
    K_candidates <- 1:5
    r_candidates <- 1:5
    
    best_bic_mfa <- Inf
    best_model_mfa <- NULL
    best_bic_pca <- Inf
    best_model_pca <- NULL
    
    for(Kc in K_candidates){
      for(rc in r_candidates){
        best_bic_temp_mfa <- Inf
        best_fit_mfa <- NULL
        best_bic_temp_pca <- Inf
        best_fit_pca <- NULL
        
        for(rep_i in seq_len(n_rep_init)){
          fit_mfa <- mfa_em_fit(list_of_data, K=Kc, r=rc, max_iter=max_iter)
          bic_mfa <- compute_bic(fit_mfa$logLik, Kc, rc, Mv, total_rows)
          if(bic_mfa < best_bic_temp_mfa){
            best_bic_temp_mfa <- bic_mfa
            best_fit_mfa <- fit_mfa
          }
          
          fit_pca <- mixture_pca_em_fit(list_of_data, K=Kc, r=rc, max_iter=max_iter)
          bic_pca <- compute_bic(fit_pca$logLik, Kc, rc, Mv, total_rows)
          if(bic_pca < best_bic_temp_pca){
            best_bic_temp_pca <- bic_pca
            best_fit_pca <- fit_pca
          }
        }
        
        if(best_bic_temp_mfa < best_bic_mfa){
          best_bic_mfa <- best_bic_temp_mfa
          best_model_mfa <- list(K=Kc, r=rc, fit=best_fit_mfa)
        }
        if(best_bic_temp_pca < best_bic_pca){
          best_bic_pca <- best_bic_temp_pca
          best_model_pca <- list(K=Kc, r=rc, fit=best_fit_pca)
        }
      }
    }
    
    z_mfa <- best_model_mfa$fit$z
    z_pca <- best_model_pca$fit$z
    ari_mfa <- mclust::adjustedRandIndex(z_mfa, z_true)
    ari_pca <- mclust::adjustedRandIndex(z_pca, z_true)
    
    df_out <- data.frame(
      K_true=Kt, r_true=rt,
      N=Nv, M=Mv, sep=sv, noise=nv,
      ARI_MFA=ari_mfa, BIC_MFA=best_bic_mfa,
      Khat_MFA=best_model_mfa$K, rhat_MFA=best_model_mfa$r,
      ARI_PCA=ari_pca, BIC_PCA=best_bic_pca,
      Khat_PCA=best_model_pca$K, rhat_PCA=best_model_pca$r,
      stringsAsFactors=FALSE
    )
    
    file_name <- sprintf("result_K%d_r%d_N%d_M%d_sep%.1f_noise%.1f.rds",
                         Kt, rt, Nv, Mv, sv, nv)
    file_path <- file.path(output_dir, file_name)
    saveRDS(df_out, file=file_path)
    
    df_out
  }
  
  idx_vec <- seq_len(nrow(param_grid))
  
  res_list <- mclapply(idx_vec, run_one, mc.cores=mc.cores)
  
  df_res <- do.call(rbind, res_list)
  df_res
}
