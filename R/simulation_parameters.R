#' Read All Results and Produce Tables/Figures
#'
#' This function reads all \code{.rds} result files in a specified directory, merges them
#' into a single data frame, prints summary tables, and generates several figures
#' (e.g. ARI vs separation/noise) as PDF files using \code{ggsave}.
#'
#' @param result_dir A character string specifying the directory containing \code{.rds} result files.
#'   Defaults to \code{"results"}.
#'
#' @return Invisibly returns the merged data frame of all results.
#'
#' @details
#' - The user can modify the conditions in \code{\link{make_table1}} and \code{\link{make_table2}}
#'   to match their actual experiment ranges.
#' - Figures are generated via \code{ggsave}.
#'
#' @examples
#' \dontrun{
#' df_all <- read_and_visualize_all("results")
#' head(df_all)
#' }
#' @export
read_and_visualize_all <- function(result_dir="results"){
  df <- read_all_results(result_dir=result_dir)

  df_table1 <- make_table1()
  cat("=== Table 1: Simulation Conditions ===\n")
  print(df_table1)

  df_table2 <- make_table2(df)
  cat("=== Table 2: Summary of ARI etc. ===\n")
  print(df_table2)

  p1 <- plot_fig1_ARI(df)
  ggsave("figure1_ARI_sep_noise.pdf", p1, width=8, height=6)

  p2 <- plot_fig2_BIC_heatmap(df)
  ggsave("figure2_BIC_heatmap.pdf", p2, width=6, height=5)

  if(nrow(df) > 1){
    p3 <- plot_fig3_init_boxplot(df)
    ggsave("figure3_init_boxplot.pdf", p3, width=6, height=4)
  }

  invisible(df)
}


#' Read All Results from .rds Files
#'
#' This function searches a directory for \code{.rds} files, each presumably containing
#' one or more rows of simulation results, then merges them into a single data frame.
#'
#' @param result_dir A character string specifying the directory to look for \code{.rds} files.
#'   Defaults to \code{"results"}.
#'
#' @return A combined \code{data.frame} with rows from all loaded files.
#'
#' @details
#' If no files are found, the function raises an error. Otherwise, each file is read with
#' \code{\link{readRDS}} and appended via \code{\link{rbind}}.
#'
#' @export
read_all_results <- function(result_dir = "results") {
  file_vec <- list.files(result_dir, pattern="\\.rds$", full.names=TRUE)
  if(length(file_vec)==0){
    stop("No .rds files found in ", result_dir)
  }
  df_list <- lapply(file_vec, readRDS)
  df_merged <- do.call(rbind, df_list)
  df_merged
}

#' Create a Static Table of Simulation Conditions (Table 1)
#'
#' This function returns a small \code{data.frame} that describes the parameter
#' ranges used in the simulation (e.g., \code{K}, \code{r}, \code{N}, \code{M},
#' \code{T_each}, etc.).
#'
#' @return A \code{data.frame} with columns \code{Parameter} and \code{Values}.
#' @export
make_table1 <- function(){
  df_table1 <- data.frame(
    Parameter = c(
      "True cluster number (K)",
      "True factor (PC) number (r)",
      "Number of subjects (N)",
      "Channels (M)",
      "Time length (T_each)",
      "Separation parameter (sep)",
      "Noise scale (noise)"
    ),
    Values = c(
      "3, 4",
      "3, 5",
      "50, 100",
      "8, 12",
      "100",
      "0.5, 1.0",
      "1.0, 2.0"
    ),
    stringsAsFactors=FALSE
  )
  df_table1
}

#' Summarize Results for a Subset (Table 2)
#'
#' This function filters the data for certain separation/noise conditions
#' (e.g., \code{sep=1.0} and \code{noise=1.0}) and computes group-wise means and
#' standard deviations for ARI and BIC. It returns a small summary table.
#'
#' @param df A \code{data.frame} with columns like \code{ARI_MFA, ARI_PCA, BIC_MFA, BIC_PCA}.
#'
#' @return A \code{data.frame} summarizing ARI/BIC means and sds, grouped by
#'   \code{(K_true, r_true, N, M)}.
#'
#' @export
make_table2 <- function(df){
  df_subset <- df %>%
    dplyr::filter(abs(sep - 1.0)<1e-9, abs(noise - 1.0)<1e-9)

  df_table2 <- df_subset %>%
    dplyr::group_by(K_true, r_true, N, M) %>%
    dplyr::summarise(
      ARI_MFA_mean = mean(ARI_MFA),
      ARI_MFA_sd   = sd(ARI_MFA),
      BIC_MFA_mean = mean(BIC_MFA),
      ARI_PCA_mean = mean(ARI_PCA),
      ARI_PCA_sd   = sd(ARI_PCA),
      BIC_PCA_mean = mean(BIC_PCA),
      .groups="drop"
    ) %>%
    dplyr::arrange(K_true, r_true, N, M)

  df_table2
}



#' Run Simulation (MFA & MPCA) in Parallel, Saving Each Condition Separately
#'
#' @param K_true_vec Vector of true K values to simulate.
#' @param r_true_vec Vector of true r values to simulate.
#' @param N_vec Vector of subject counts.
#' @param M_vec Vector of channel counts.
#' @param T_each Time-series length.
#' @param sep_vec Vector of separation parameters.
#' @param noise_vec Vector of noise scale parameters.
#' @param K_candidates Vector of candidate K values to try.
#' @param r_candidates Vector of candidate r values to try.
#' @param max_iter Max iteration for EM.
#' @param n_rep_init Number of random initial seeds for each (K,r) in search (used as `n_init` in new `mfa_em_fit`).
#' @param seed_data_base Base seed for data simulation.
#' @param mc.cores Number of parallel cores to use.
#' @param output_dir Directory to save each result as \code{result_K{...}_r{...}...rds}.
#' @param nIterFA  Number of sub-iterations for the factor analyzer update in C++ (passed to `mfa_em_fit`).
#' @param tol  Convergence tolerance for `mfa_em_fit`.
#' @param method_pca Which method for mpcaTimeseriesCpp ("EM" or "closed_form").
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
    K_candidates = 1:5,
    r_candidates = 1:5,
    max_iter   = 30,
    n_rep_init = 5,
    seed_data_base = 2023,
    mc.cores   = parallel::detectCores() - 1,
    output_dir = "results",
    # 以下は新しい引数として追加
    nIterFA    = 20,            # MFA用
    nIterPCA   = 20,            # MPCA用
    tol        = 1e-3,
    method_pca = "EM"           # Mixture PCAで使うmethod
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

    # (1) データ生成
    sim_data <- simulate_mixture_data(
      N=Nv, K=Kt, r=rt, M=Mv, T_each=T_each,
      cluster_sep=sv, noise_scale=nv,
      seed=data_seed
    )
    list_of_data <- sim_data$list_of_data
    z_true       <- sim_data$z_true
    total_rows   <- sum(sapply(list_of_data, nrow))

    best_bic_mfa <- Inf
    best_model_mfa <- NULL
    best_bic_pca <- Inf
    best_model_pca <- NULL

    #------------------------------------------------------
    # (A) MFA: mfa_em_fit()を用いて最良モデル探索
    #------------------------------------------------------
    for(Kc in K_candidates){
      for(rc in r_candidates){
        # Fit MFA
        fit_mfa <- mfa_em_fit(
          list_of_data = list_of_data,
          K            = Kc,
          r            = rc,
          max_iter     = max_iter,
          nIterFA      = nIterFA,
          tol          = tol,
          n_init       = n_rep_init,
          use_kmeans_init = FALSE,
          subject_rdim_for_kmeans = rc,
          mc_cores     = 1
        )
        fit_mfa$list_of_data <- list_of_data
        bic_mfa <- compute_bic(fit_mfa$logLik, Kc, rc, Mv, total_rows)

        if(bic_mfa < best_bic_mfa){
          best_bic_mfa <- bic_mfa
          best_model_mfa <- list(K=Kc, r=rc, fit=fit_mfa)
        }

        # MPCA
        fit_pca <- mixture_pca_em_fit(
          list_of_data = list_of_data,
          K            = Kc,
          r            = rc,
          max_iter     = max_iter,
          nIterPCA     = nIterPCA,
          tol          = tol,
          method       = method_pca,
          n_init       = n_rep_init,
          use_kmeans_init = FALSE,
          subject_rdim_for_kmeans = rc,
          mc_cores     = 1
        )

        fit_pca$list_of_data <- list_of_data
        bic_pca <- compute_bic(fit_pca$logLik, Kc, rc, Mv, total_rows)

        if(bic_pca < best_bic_pca){
          best_bic_pca <- bic_pca
          best_model_pca <- list(K=Kc, r=rc, fit=fit_pca)
        }
      }
    }

    # (2) ARI
    z_mfa <- best_model_mfa$fit$z
    z_pca <- best_model_pca$fit$z
    ari_mfa <- mclust::adjustedRandIndex(z_mfa, z_true)
    ari_pca <- mclust::adjustedRandIndex(z_pca, z_true)

    # (3) 1行の結果 + モデルも一緒に保存したいので，以下のようにまとめる
    df_out <- list(
      # 結果サマリ(1行分)
      summary = data.frame(
        K_true=Kt, r_true=rt,
        N=Nv, M=Mv, sep=sv, noise=nv,
        ARI_MFA=ari_mfa, BIC_MFA=best_bic_mfa,
        Khat_MFA=best_model_mfa$K, rhat_MFA=best_model_mfa$r,
        ARI_PCA=ari_pca, BIC_PCA=best_bic_pca,
        Khat_PCA=best_model_pca$K, rhat_PCA=best_model_pca$r,
        stringsAsFactors=FALSE
      ),
      # 実際の最良モデルオブジェクト
      best_model_mfa = best_model_mfa$fit,
      best_model_pca = best_model_pca$fit,
      # 真のラベル (z_true) などがあれば後で使う
      z_true = z_true
    )

    # (4) ファイル保存
    file_name <- sprintf("result_K%d_r%d_N%d_M%d_sep%.1f_noise%.1f.rds",
                         Kt, rt, Nv, Mv, sv, nv)
    file_path <- file.path(output_dir, file_name)
    saveRDS(df_out, file=file_path)

    # 戻り値は 1行のdata.frame だけにする
    df_out$summary
  }

  idx_vec <- seq_len(nrow(param_grid))
  res_list <- parallel::mclapply(idx_vec, run_one, mc.cores=mc.cores)
  df_res <- do.call(rbind, res_list)
  df_res
}


#' Run Simulation (MFA & MPCA) for Each Param in a For-Loop, Parallelizing Over Multiple Initial Seeds

#' @param K_true_vec Vector of true K values to simulate.
#' @param r_true_vec Vector of true r values to simulate.
#' @param N_vec Vector of subject counts.
#' @param M_vec Vector of channel counts.
#' @param T_each Time-series length.
#' @param sep_vec Vector of separation parameters.
#' @param noise_vec Vector of noise scale parameters.
#' @param K_candidates Vector of candidate K values to try.
#' @param r_candidates Vector of candidate r values to try.
#' @param max_iter Max iteration for EM.
#' @param n_rep_init Number of random initial seeds for each (K,r).
#' @param mc.cores_inits Number of parallel cores to use **for initial seeds**.
#' @param seed_data_base Base seed for data simulation.
#' @param output_dir Directory to save each result as rds.
#' @param ... Any other additional arguments you wish to pass to mfa_em_fit_cpp_singleInit() or mixture_pca_em_fit_cpp_singleInit()
#'
#' @return A data frame with columns:
#' \item{K_true, r_true, N, M, sep, noise}{True simulation parameters.}
#' \item{ARI_MFA, BIC_MFA, Khat_MFA, rhat_MFA, ARI_PCA, BIC_PCA, Khat_PCA, rhat_PCA}{Results.}
#'
#' @export
simulate_and_run_serialParam_parallelInit_save <- function(
    K_true_vec = c(3,4),
    r_true_vec = c(3,5),
    N_vec      = c(50),
    M_vec      = c(8,12),
    T_each     = 100,
    sep_vec    = c(0.5,1.0),
    noise_vec  = c(1.0,2.0),
    K_candidates = 1:5,
    r_candidates = 1:5,
    max_iter   = 30,
    n_rep_init = 5,
    mc.cores_inits = 2,  # 初期値の並列数
    seed_data_base = 2023,
    output_dir = "results",
    ...
){
  if(!dir.exists(output_dir)){
    dir.create(output_dir, recursive=TRUE)
  }

  # パラメータグリッド作成
  param_grid <- expand.grid(
    K_true = K_true_vec,
    r_true = r_true_vec,
    N      = N_vec,
    M      = M_vec,
    sep    = sep_vec,
    noise  = noise_vec,
    stringsAsFactors=FALSE
  )

  out_list <- vector("list", nrow(param_grid))

  # パラメータグリッドを forループで順に処理
  for(idx in seq_len(nrow(param_grid))){
    row_i <- param_grid[idx, ]
    Kt <- row_i$K_true
    rt <- row_i$r_true
    Nv <- row_i$N
    Mv <- row_i$M
    sv <- row_i$sep
    nv <- row_i$noise

    data_seed <- seed_data_base + idx

    # (1) データシミュレーション
    sim_data <- simulate_mixture_data(
      N=Nv, K=Kt, r=rt, M=Mv, T_each=T_each,
      cluster_sep=sv, noise_scale=nv,
      seed=data_seed
    )
    list_of_data <- sim_data$list_of_data
    z_true       <- sim_data$z_true
    total_rows   <- sum(sapply(list_of_data, nrow))

    best_bic_mfa <- Inf
    best_model_mfa <- NULL
    best_bic_pca <- Inf
    best_model_pca <- NULL

    # ここでは (K_candidates, r_candidates) をさらにループ
    for(Kc in K_candidates){
      for(rc in r_candidates){

        #======================#
        #      MFA 部分
        #======================#
        # (A) 各初期値を並列化して試す -> bestを取得
        run_mfa_one_init <- function(rep_id){
          # ランダム初期割り当てを作る
          z_init_rand <- sample.int(Kc, size=length(list_of_data), replace=TRUE)

          # 1回分のフィット
          fit_mfa_single <- mfa_em_fit_cpp_singleInit(
            list_of_data = list_of_data,
            K = Kc,
            r = rc,
            max_iter = max_iter,
            z_init   = z_init_rand,
            ...      # 他にnIterFA, tolなどを渡したい場合
          )

          # BIC計算用に返す
          bic_val <- compute_bic(fit_mfa_single$logLik, K=Kc, r=rc, M=Mv, N_total_rows=total_rows)
          list(bic=bic_val, fit=fit_mfa_single)
        }

        # 初期値ごとに並列実行
        mfa_fit_list <- parallel::mclapply(
          seq_len(n_rep_init),
          FUN=run_mfa_one_init,
          mc.cores=mc.cores_inits
        )
        # 一番BICが小さい(=良い)ものをbestとする
        best_local_mfa <- mfa_fit_list[[which.min(sapply(mfa_fit_list, `[[`, "bic"))]]
        # もし全体ベストより良ければ更新
        if(best_local_mfa$bic < best_bic_mfa){
          best_bic_mfa    <- best_local_mfa$bic
          best_model_mfa  <- list(K=Kc, r=rc, fit=best_local_mfa$fit)
        }

        #======================#
        #   Mixture PCA 部分
        #======================#
        run_pca_one_init <- function(rep_id){
          # ランダム初期割り当て
          z_init_rand <- sample.int(Kc, size=length(list_of_data), replace=TRUE)

          # 1回分のフィット
          fit_pca_single <- mixture_pca_em_fit_cpp_singleInit(
            list_of_data = list_of_data,
            K = Kc,
            r = rc,
            max_iter = max_iter,
            z_init   = z_init_rand,
            ...      # nIterPCA, tol, methodなどを渡す
          )
          bic_val <- compute_bic(fit_pca_single$logLik, K=Kc, r=rc, M=Mv, N_total_rows=total_rows)
          list(bic=bic_val, fit=fit_pca_single)
        }

        pca_fit_list <- parallel::mclapply(
          seq_len(n_rep_init),
          FUN=run_pca_one_init,
          mc.cores=mc.cores_inits
        )
        best_local_pca <- pca_fit_list[[which.min(sapply(pca_fit_list, `[[`, "bic"))]]
        if(best_local_pca$bic < best_bic_pca){
          best_bic_pca   <- best_local_pca$bic
          best_model_pca <- list(K=Kc, r=rc, fit=best_local_pca$fit)
        }

      } # end for(r_candidates)
    } # end for(K_candidates)

    # (2) ARI計算
    z_mfa <- best_model_mfa$fit$z
    z_pca <- best_model_pca$fit$z
    ari_mfa <- mclust::adjustedRandIndex(z_mfa, z_true)
    ari_pca <- mclust::adjustedRandIndex(z_pca, z_true)

    df_out_one <- data.frame(
      K_true=Kt, r_true=rt,
      N=Nv, M=Mv, sep=sv, noise=nv,
      ARI_MFA=ari_mfa, BIC_MFA=best_bic_mfa,
      Khat_MFA=best_model_mfa$K, rhat_MFA=best_model_mfa$r,
      ARI_PCA=ari_pca, BIC_PCA=best_bic_pca,
      Khat_PCA=best_model_pca$K, rhat_PCA=best_model_pca$r,
      stringsAsFactors=FALSE
    )

    # RDSに保存
    file_name <- sprintf("result_serialParam_K%d_r%d_N%d_M%d_sep%.1f_noise%.1f.rds",
                         Kt, rt, Nv, Mv, sv, nv)
    saveRDS(df_out_one, file=file.path(output_dir, file_name))

    out_list[[idx]] <- df_out_one
  }

  # 最後に結合
  df_res <- do.call(rbind, out_list)
  df_res
}

#' Main Example: Run Parallel MFA/MPCA Simulation and Save Each Condition
#'
#' @return The merged data frame of results across all conditions.
#' @examples
#' \dontrun{
#' df_out <- run_main_experiment_parallel_save()
#' head(df_out)
#' }
#' @export
run_main_experiment_parallel_save <- function(){
  df_result <- simulate_and_run_parallel_save(
    K_true_vec = c(3,4),
    r_true_vec = c(3,4),
    N_vec      = c(50),
    M_vec      = c(8,12),
    T_each     = 100,
    sep_vec    = c(0.5,1.0),
    noise_vec  = c(1.0,2.0),
    K_candidates = 1:5,
    r_candidates = 1:5,
    max_iter   = 20,
    n_rep_init = 3,
    seed_data_base = 999,
    mc.cores   = 2,
    output_dir = "results",
    nIterFA    = 20,    # MFA: sub-iterations
    nIterPCA   = 20,    # PCA: sub-iterations
    tol        = 1e-3,
    method_pca = "EM"
  )

  # 結果の確認
  print(head(df_result, 10))

  df_long <- df_result %>%
    dplyr::select(K_true, r_true, N, M, sep, noise, ARI_MFA, ARI_PCA) %>%
    tidyr::gather(key="Method", value="ARI", ARI_MFA, ARI_PCA)

  p <- ggplot(df_long, aes(x=factor(sep), y=ARI, color=factor(noise),
                           group=interaction(Method, noise))) +
    geom_point(position=position_dodge(width=0.3)) +
    geom_line(position=position_dodge(width=0.3)) +
    facet_wrap(~Method + K_true + r_true + N + M, labeller=label_both) +
    theme_bw() +
    labs(title="MFA vs MPCA (parallel, each condition saved separately)",
         x="Separation (sep)", color="Noise", y="ARI")

  ggsave("figure_ari_mfa_mpca_parallel_save.pdf", p, width=8, height=6)

  df_result
}


#' Simulate Mixture Data for MFA/MPCA
#'
#' （以下、ドキュメント部分は元のまま）
#' @param N Number of subjects (time-series).
#' @param K True number of clusters in simulation.
#' @param r True dimension of factors or principal components.
#' @param M Number of observed channels (e.g. muscle EMG channels).
#' @param T_each Length of each time-series per subject.
#' @param cluster_sep How far to separate cluster means.
#' @param noise_scale Scaling of noise diagonal in \code{\link{Psi_list}}.
#' @param seed Random seed for reproducibility.
#' @return A list with: \code{list_of_data}, \code{z_true}.
#' @export
simulate_mixture_data <- function(
    N = 50,
    K = 3,
    r = 3,
    M = 8,
    T_each = 100,
    cluster_sep = 1.0,
    noise_scale = 1.0,
    seed = 123
){
  set.seed(seed)

  pi_true <- rep(1/K, K)
  z_true  <- sample.int(K, size=N, replace=TRUE, prob=pi_true)

  Lambda_list <- vector("list", K)
  mu_list     <- vector("list", K)
  Psi_list    <- vector("list", K)

  for(k in seq_len(K)){
    base_mu <- rnorm(M, mean=0, sd=0.3)
    offset  <- (k - (K+1)/2)* cluster_sep
    mu_list[[k]] <- base_mu + offset

    lambda_scale <- 0.2 + 0.2*cluster_sep
    Lambda_list[[k]] <- matrix(rnorm(M*r, sd=lambda_scale), nrow=M, ncol=r)

    psi_vec <- runif(M, min=0.05, max=0.2)*noise_scale
    Psi_list[[k]] <- diag(psi_vec, nrow=M)
  }

  list_of_data <- vector("list", N)
  for(i in seq_len(N)){
    k_i <- z_true[i]
    Lk  <- Lambda_list[[k_i]]
    mu_k<- mu_list[[k_i]]
    Psi_k <- Psi_list[[k_i]]
    Sigma_k <- Lk %*% t(Lk) + Psi_k
    Xi <- MASS::mvrnorm(n=T_each, mu=mu_k, Sigma=Sigma_k)
    list_of_data[[i]] <- Xi
  }

  list(
    list_of_data = list_of_data,
    z_true       = z_true
  )
}

#' Run Simulation (MFA & MPCA) for Each Param in a For-Loop, Parallelizing Over Multiple Initial Seeds
#'
#' @param K_true_vec Vector of true K values to simulate.
#' @param r_true_vec Vector of true r values to simulate.
#' @param N_vec Vector of subject counts.
#' @param M_vec Vector of channel counts.
#' @param T_each Time-series length.
#' @param sep_vec Vector of separation parameters.
#' @param noise_vec Vector of noise scale parameters.
#' @param K_candidates Vector of candidate K values to try.
#' @param r_candidates Vector of candidate r values to try.
#' @param max_iter Max iteration for EM.
#' @param n_rep_init Number of random initial seeds for each (K,r).
#' @param mc.cores_inits Number of parallel cores to use **for initial seeds**.
#'   (パラメータグリッド自体はforループで回すので並列化しない)
#' @param seed_data_base Base seed for data simulation.
#' @param output_dir Directory to save each result as rds.
#' @param ... その他、mfa_em_fit_cpp_singleInit() や mixture_pca_em_fit_cpp_singleInit() に渡したい追加引数
#'
#' @return A data frame with columns:
#' \item{K_true, r_true, N, M, sep, noise}{True simulation parameters.}
#' \item{ARI_MFA, BIC_MFA, Khat_MFA, rhat_MFA, ARI_PCA, BIC_PCA, Khat_PCA, rhat_PCA}{Results.}
#'
#' @export
simulate_and_run_serialParam_parallelInit_save <- function(
    K_true_vec = c(3,4),
    r_true_vec = c(3,5),
    N_vec      = c(50),
    M_vec      = c(8,12),
    T_each     = 100,
    sep_vec    = c(0.5,1.0),
    noise_vec  = c(1.0,2.0),
    K_candidates = 1:5,
    r_candidates = 1:5,
    max_iter   = 30,
    n_rep_init = 5,
    mc.cores_inits = 2,  # 初期値の並列数
    seed_data_base = 2023,
    output_dir = "results",
    ...
){
  if(!dir.exists(output_dir)){
    dir.create(output_dir, recursive=TRUE)
  }
  # パラメータグリッド作成
  param_grid <- expand.grid(
    K_true = K_true_vec,
    r_true = r_true_vec,
    N      = N_vec,
    M      = M_vec,
    sep    = sep_vec,
    noise  = noise_vec,
    stringsAsFactors=FALSE
  )

  out_list <- vector("list", nrow(param_grid))

  # パラメータグリッドを forループで順に処理
  for(idx in seq_len(nrow(param_grid))){
    row_i <- param_grid[idx, ]
    Kt <- row_i$K_true
    rt <- row_i$r_true
    Nv <- row_i$N
    Mv <- row_i$M
    sv <- row_i$sep
    nv <- row_i$noise

    data_seed <- seed_data_base + idx

    # (1) データシミュレーション
    sim_data <- simulate_mixture_data(
      N=Nv, K=Kt, r=rt, M=Mv, T_each=T_each,
      cluster_sep=sv, noise_scale=nv,
      seed=data_seed
    )
    list_of_data <- sim_data$list_of_data
    z_true       <- sim_data$z_true
    total_rows   <- sum(sapply(list_of_data, nrow))

    best_bic_mfa <- Inf
    best_model_mfa <- NULL
    best_bic_pca <- Inf
    best_model_pca <- NULL

    # ここでは (K_candidates, r_candidates) をさらにループ
    for(Kc in K_candidates){
      for(rc in r_candidates){

        #======================#
        #      MFA 部分
        #======================#
        # (A) 各初期値を並列化して試す -> bestを取得
        run_mfa_one_init <- function(rep_id){
          # ランダム初期割り当てを作る
          z_init_rand <- sample.int(Kc, size=length(list_of_data), replace=TRUE)

          # 1回分のフィット
          fit_mfa_single <- mfa_em_fit_cpp_singleInit(
            list_of_data = list_of_data,
            K = Kc,
            r = rc,
            max_iter = max_iter,
            z_init   = z_init_rand,
            ...      # 他にnIterFA, tolなどを渡したい場合
          )

          # BIC計算用に返す
          bic_val <- compute_bic(fit_mfa_single$logLik, K=Kc, r=rc, M=Mv, N_total_rows=total_rows)
          list(bic=bic_val, fit=fit_mfa_single)
        }

        # 初期値ごとに並列実行
        mfa_fit_list <- parallel::mclapply(
          seq_len(n_rep_init),
          FUN=run_mfa_one_init,
          mc.cores=mc.cores_inits
        )
        # 一番BICが小さい(=良い)ものをbestとする
        best_local_mfa <- mfa_fit_list[[which.min(sapply(mfa_fit_list, `[[`, "bic"))]]
        # もし全体ベストより良ければ更新
        if(best_local_mfa$bic < best_bic_mfa){
          best_bic_mfa    <- best_local_mfa$bic
          best_model_mfa  <- list(K=Kc, r=rc, fit=best_local_mfa$fit)
        }

        #======================#
        #   Mixture PCA 部分
        #======================#
        run_pca_one_init <- function(rep_id){
          # ランダム初期割り当て
          z_init_rand <- sample.int(Kc, size=length(list_of_data), replace=TRUE)

          # 1回分のフィット
          fit_pca_single <- mixture_pca_em_fit_cpp_singleInit(
            list_of_data = list_of_data,
            K = Kc,
            r = rc,
            max_iter = max_iter,
            z_init   = z_init_rand,
            ...      # nIterPCA, tol, methodなどを渡す
          )
          bic_val <- compute_bic(fit_pca_single$logLik, K=Kc, r=rc, M=Mv, N_total_rows=total_rows)
          list(bic=bic_val, fit=fit_pca_single)
        }

        pca_fit_list <- parallel::mclapply(
          seq_len(n_rep_init),
          FUN=run_pca_one_init,
          mc.cores=mc.cores_inits
        )
        best_local_pca <- pca_fit_list[[which.min(sapply(pca_fit_list, `[[`, "bic"))]]
        if(best_local_pca$bic < best_bic_pca){
          best_bic_pca   <- best_local_pca$bic
          best_model_pca <- list(K=Kc, r=rc, fit=best_local_pca$fit)
        }

      } # end for(r_candidates)
    } # end for(K_candidates)

    # (2) ARI計算
    z_mfa <- best_model_mfa$fit$z
    z_pca <- best_model_pca$fit$z
    ari_mfa <- mclust::adjustedRandIndex(z_mfa, z_true)
    ari_pca <- mclust::adjustedRandIndex(z_pca, z_true)

    df_out_one <- data.frame(
      K_true=Kt, r_true=rt,
      N=Nv, M=Mv, sep=sv, noise=nv,
      ARI_MFA=ari_mfa, BIC_MFA=best_bic_mfa,
      Khat_MFA=best_model_mfa$K, rhat_MFA=best_model_mfa$r,
      ARI_PCA=ari_pca, BIC_PCA=best_bic_pca,
      Khat_PCA=best_model_pca$K, rhat_PCA=best_model_pca$r,
      stringsAsFactors=FALSE
    )

    # RDSに保存
    file_name <- sprintf("result_serialParam_K%d_r%d_N%d_M%d_sep%.1f_noise%.1f.rds",
                         Kt, rt, Nv, Mv, sv, nv)
    saveRDS(df_out_one, file=file.path(output_dir, file_name))

    out_list[[idx]] <- df_out_one
  }

  # 最後に結合
  df_res <- do.call(rbind, out_list)
  df_res
}


#' Compute BIC for MFA or Mixture PCA
#'
#' @param logLik The log-likelihood value.
#' @param K Number of clusters.
#' @param r Factor or principal component dimension.
#' @param M Number of observed channels.
#' @param N_total_rows Total number of rows across all data (used for log(N)).
#' @export
compute_bic <- function(logLik, K, r, M, N_total_rows){
  num_params <- K*(M*r + 2*M) + (K-1)
  -2*logLik + num_params*log(N_total_rows)
}


#' Read All Results from .rds Files
#'
#' @param result_dir A character string specifying the directory to look for \code{.rds} files.
#' @return A combined \code{data.frame} with rows from all the loaded files.
#' @export
read_all_results <- function(result_dir = "results") {
  file_vec <- list.files(result_dir, pattern="\\.rds$", full.names=TRUE)
  if(length(file_vec)==0){
    stop("No .rds files found in ", result_dir)
  }
  df_list <- lapply(file_vec, function(f){
    x <- readRDS(f)
    x
  })
  df_merged <- do.call(rbind, df_list)
  df_merged
}


#' Plot ARI vs. Separation/Noise (Figure 1)
#'
#' @param df A \code{data.frame}.
#' @return A ggplot object.
#' @export
plot_fig1_ARI <- function(df){
  df_long <- df %>%
    dplyr::select(K_true, r_true, N, M, sep, noise, ARI_MFA, ARI_PCA) %>%
    tidyr::pivot_longer(cols = c("ARI_MFA", "ARI_PCA"),
                        names_to = "Method", values_to = "ARI")

  p <- ggplot(df_long, aes(x = factor(sep), y = ARI, color = factor(noise),
                           group = interaction(Method, noise))) +
    geom_point(position = position_dodge(width = 0.3)) +
    geom_line(position = position_dodge(width = 0.3)) +
    facet_wrap(~ Method + K_true + r_true + N + M, labeller = label_both) +
    theme_bw() +
    labs(x = "Separation (sep)", y = "ARI", color = "Noise",
         title = "(Figure 1) ARI vs sep/noise (MFA, MPCA)")

  return(p)
}


#' Plot BIC Heatmap for a Subset (Figure 2)
#'
#' @param df A \code{data.frame}.
#' @return A ggplot object.
#' @export
plot_fig2_BIC_heatmap <- function(df){
  sub <- df %>%
    dplyr::filter(N==50, M==8, abs(sep-1.0)<1e-9, abs(noise-1.0)<1e-9)

  p <- ggplot(sub, aes(x=factor(r_true), y=factor(K_true), fill=BIC_MFA)) +
    geom_tile() +
    scale_fill_viridis_c() +
    theme_bw() +
    labs(x="r", y="K", fill="BIC",
         title="(Figure 2) BIC heatmap (MFA, N=50, M=8, sep=1.0, noise=1.0)")
  p
}


#' Plot Boxplot of ARI for Multiple Initial Seeds (Figure 3)
#'
#' （以下、ドキュメント部分は元のまま）
#' @param df_multiple A \code{data.frame} with columns \code{ARI_MFA, ARI_PCA}.
#' @export
plot_fig3_init_boxplot <- function(df_multiple){
  df_long <- df_multiple %>%
    tidyr::pivot_longer(cols=c("ARI_MFA","ARI_PCA"),
                        names_to="Method", values_to="ARI")

  p <- ggplot(df_long, aes(x=Method, y=ARI, color=Method)) +
    geom_boxplot() +
    theme_bw() +
    labs(title="(Figure 3) Boxplot for multiple initial seeds")
  p
}


#' Plot Factor Loadings Heatmap (Figure 4)
#'
#' @param mfa_fit A fitted MFA model with \code{$Lambda}.
#' @param cluster_id Which cluster's loadings to plot (1..K).
#' @export
plot_fig4_loadings_heatmap <- function(mfa_fit, cluster_id=1){
  if(is.null(mfa_fit$Lambda)){
    stop("mfa_fit does not contain Lambda.")
  }
  Lambda_k <- mfa_fit$Lambda[[cluster_id]]
  if(is.null(Lambda_k)){
    stop("No cluster_id = ", cluster_id, " in Lambda.")
  }
  df_melt <- reshape2::melt(Lambda_k)
  colnames(df_melt) <- c("Muscle","Factor","Loading")

  p <- ggplot(df_melt, aes(x=factor(Factor), y=factor(Muscle), fill=Loading)) +
    geom_tile() +
    scale_fill_gradient2(mid="white", low="blue", high="red") +
    theme_bw() +
    labs(x="Factor", y="Muscle", fill="Loading",
         title=sprintf("(Figure 4) Cluster %d synergy loadings", cluster_id))
  p
}

#' Create a Static Table of Simulation Conditions (Table 1)
#'
#' @export
make_table1 <- function(){
  df_table1 <- data.frame(
    Parameter = c(
      "True cluster number (K)",
      "True factor (PC) number (r)",
      "Number of subjects (N)",
      "Channels (M)",
      "Time length (T_each)",
      "Separation parameter (sep)",
      "Noise scale (noise)"
    ),
    Values = c(
      "3, 4",
      "3, 5",
      "50, 100",
      "8, 12",
      "100",
      "0.5, 1.0",
      "1.0, 2.0"
    ),
    stringsAsFactors=FALSE
  )
  df_table1
}


#' Summarize Results for a Subset (Table 2)
#'
#' @param df A \code{data.frame}.
#' @return A \code{data.frame} summarizing ARI/BIC means and sds.
#' @export
make_table2 <- function(df){
  df_subset <- df %>%
    dplyr::filter(abs(sep - 1.0)<1e-9, abs(noise - 1.0)<1e-9)

  df_table2 <- df_subset %>%
    dplyr::group_by(K_true, r_true, N, M) %>%
    dplyr::summarise(
      ARI_MFA_mean = mean(ARI_MFA),
      ARI_MFA_sd   = sd(ARI_MFA),
      BIC_MFA_mean = mean(BIC_MFA),
      ARI_PCA_mean = mean(ARI_PCA),
      ARI_PCA_sd   = sd(ARI_PCA),
      BIC_PCA_mean = mean(BIC_PCA),
      .groups="drop"
    ) %>%
    dplyr::arrange(K_true, r_true, N, M)

  df_table2
}

#' Compute VAF for a Fitted MFA Model
#'
#' This function computes the Variance Accounted For (VAF) for a fitted MFA model
#' across all subjects, reconstructing each subject's data via factor scores and
#' summing squared errors.
#'
#' @param list_of_data A list of subject data matrices, each \code{(T_i x M)}.
#' @param mfa_fit An MFA fit object (from \code{mfa_em_fit} or similar),
#'   containing \code{Lambda}, \code{Psi}, \code{mu}, etc.
#' @param z An integer vector of length \code{N} specifying each subject's cluster ID (1..K).
#'
#' @return A numeric scalar representing the total VAF over all subjects:
#'   \code{1 - (sum of squared residual) / (sum of squared original data)}.
#'
#' @export
compute_VAF_mfa <- function(list_of_data, mfa_fit, z){
  N <- length(list_of_data)
  K <- length(mfa_fit$Lambda)
  # check
  if(length(z) != N) stop("length(z) != N")

  total_num <- 0
  total_den <- 0

  for(i in seq_len(N)){
    Xi <- list_of_data[[i]]         # (T_i x M)
    k_i <- z[i]                     # クラスタID
    Lambda_k <- mfa_fit$Lambda[[k_i]]  # (M x r)
    Psi_k    <- mfa_fit$Psi[[k_i]]     # (M x M) diagonal
    mu_k     <- mfa_fit$mu[[k_i]]      # (M)

    # invert Psi_k
    invPsi <- diag(1 / diag(Psi_k))  # (M x M)

    # precompute: A = (I_r + Lambda^T invPsi Lambda)^-1, B = Lambda^T invPsi
    M_ <- nrow(Lambda_k)
    r_ <- ncol(Lambda_k)
    A <- solve(diag(r_) + t(Lambda_k) %*% invPsi %*% Lambda_k)
    B <- t(Lambda_k) %*% invPsi

    # 再構成
    # Xhat_i(t) = mu_k + Lambda_k f_t,  where f_t = A B [x_t - mu_k]
    # x_t, mu_k は (M)-vector
    X_centered <- sweep(Xi, 2, mu_k, FUN="-")  # (T_i x M)
    # f_t = (T_i x r)
    f_mat <- (X_centered %*% t(B)) %*% t(A)  # [T_i x M]*[M x r]*[r x r] => [T_i x r]

    # Xhat = (T_i x M)
    Xhat <- matrix(NA, nrow=nrow(Xi), ncol=ncol(Xi))
    for(t in seq_len(nrow(Xi))){
      Xhat[t, ] <- mu_k + Lambda_k %*% f_mat[t, ]
    }

    ss_res <- sum((Xi - Xhat)^2)
    ss_total <- sum(Xi^2)

    total_num <- total_num + ss_res
    total_den <- total_den + ss_total
  }

  VAF <- 1 - total_num/total_den
  VAF
}

#' Compute VAF for a Mixture PCA (PPCA) Model
#'
#' This function computes the Variance Accounted For (VAF) for a fitted Mixture PCA model
#' across all subjects, reconstructing each subject's data via principal component
#' scores and summing squared errors.
#'
#' @inheritParams compute_VAF_mfa
#' @param mpca_fit A fitted MPCA model with \code{P, D, Psi, mu}.
#'
#' @return A numeric scalar for the overall VAF.
#' @export
compute_VAF_mpca <- function(list_of_data, mpca_fit, z){
  N <- length(list_of_data)
  if(length(z) != N) {
    stop("length(z) != N (subject assignments do not match number of data sets).")
  }

  total_num <- 0  # sum of squared residual
  total_den <- 0  # sum of squared original data

  for(i in seq_len(N)){
    Xi <- list_of_data[[i]]   # (T_i x M)
    k_i <- z[i]               # cluster ID in {1..K}

    # (M x r) P matrix for cluster k_i
    P_k <- mpca_fit$P[[k_i]]
    # (r x r) diagonal D
    D_k <- mpca_fit$D[[k_i]]
    # (M x M) diagonal noise matrix Psi
    Psi_k <- mpca_fit$Psi[[k_i]]  # e.g. diag(sigma2_k, M, M)
    mu_k  <- mpca_fit$mu[[k_i]]   # (M x 1) vector

    # form W = P * D  => (M x r)
    W_k <- P_k %*% D_k

    # If Psi_k is truly diagonal with each diag ~ sigma2, you can also do
    #   sigma2k <- mpca_fit$sigma2[k_i]
    # but let's assume a general diagonal for now.
    # We'll do the isotropic assumption for the PC score formula, or
    # if each diag is different, need a more general approach.

    # For simple isotropic noise, we can do:
    #   A = (W_k^T W_k + sigma2 I_r)^-1
    #   B = W_k^T
    # but if Psi_k has different diagonal elements, we need to invert that
    # for each row. Let's check if diagonal is constant or not:
    diagPsi <- diag(Psi_k)
    if(!all(abs(diff(diagPsi)) < 1e-10)){
      # non-isotropic => more general approach
      # scores = (W^T Psi^-1 W + I_r)^-1 W^T Psi^-1 (x_t - mu)
      # => we need to invert Psi
      invPsi <- diag(1 / diagPsi)
      M_ <- nrow(W_k)
      r_ <- ncol(W_k)
      # For each row x_t:
      X_centered <- sweep(Xi, 2, mu_k, FUN = "-")  # (T_i x M)
      for(t in seq_len(nrow(Xi))){
        x_t <- X_centered[t, ]
        # z_t = (W^T Psi^-1 W + I_r)^{-1} W^T Psi^-1 x_t
        # dimension: r x r => invert
        tmp <- t(W_k) %*% invPsi %*% W_k + diag(r_)
        A <- solve(tmp)
        z_t <- A %*% (t(W_k) %*% invPsi %*% x_t)
        x_hat <- mu_k + W_k %*% z_t
        total_num <- total_num + sum((Xi[t, ] - x_hat)^2)
        total_den <- total_den + sum(Xi[t, ]^2)
      }
    } else {
      # isotropic => diagPsi ~ sigma2
      sigma2k <- diagPsi[1]
      M_ <- nrow(W_k)
      r_ <- ncol(W_k)
      A <- solve(t(W_k) %*% W_k + sigma2k * diag(r_))
      B <- t(W_k)

      # reconstruct each row
      X_centered <- sweep(Xi, 2, mu_k, FUN="-")  # (T_i x M)
      z_mat <- (X_centered %*% t(B)) %*% A       # [T_i x r]

      Xhat <- (z_mat %*% t(W_k))  # [T_i x r] x [r x M] => (T_i x M)
      # add mu
      for(t in seq_len(nrow(Xhat))){
        Xhat[t, ] <- Xhat[t, ] + mu_k
      }

      ss_res <- sum((Xi - Xhat)^2)
      ss_total <- sum(Xi^2)

      total_num <- total_num + ss_res
      total_den <- total_den + ss_total
    }
  }

  VAF <- 1 - total_num / total_den
  VAF
}

#' Compile ARI and VAF into a Single Table
#'
#' This function scans a directory for \code{.rds} simulation result files,
#' each containing a list with a \code{summary} data frame plus \code{best_model_mfa},
#' \code{best_model_pca}, and \code{z_true}. It computes VAF for MFA and MPCA (if available),
#' merges these into one final data frame, and returns it.
#'
#' @param result_dir Directory containing \code{.rds} files. Defaults to \code{"results"}.
#'
#' @details
#' Each \code{.rds} file is assumed to have been saved as a list with elements:
#' \itemize{
#'   \item \code{summary} - a data frame with a single row of ARI/BIC results.
#'   \item \code{best_model_mfa}, \code{best_model_pca} - the best-fitting models.
#'   \item \code{z_true} - (optional) the true cluster labels.
#' }
#' The fitted models should include \code{list_of_data}, or you must ensure that
#' \code{list_of_data} is otherwise accessible to compute VAF.
#'
#' @return A data frame combining ARI, BIC, and newly calculated \code{VAF_MFA, VAF_PCA}.
#'
#' @export
compute_ari_vaf_table <- function(result_dir="results"){
  file_vec <- list.files(result_dir, pattern="\\.rds$", full.names=TRUE)
  if(length(file_vec) < 1){
    stop("No .rds files found in ", result_dir)
  }

  df_all <- list()

  for(f in file_vec){
    res_list <- readRDS(f)
    # res_list$summary: data.frame(1行), ARI_MFA, ARI_PCA など
    # res_list$best_model_mfa: fitted object
    # res_list$best_model_pca: fitted object
    # res_list$z_true

    df1 <- res_list$summary  # 1行だけ

    # ---- VAF計算 ----
    # もし best_model_mfa が NULL なら何もしない
    if(!is.null(res_list$best_model_mfa)){
      vaf_mfa <- compute_VAF_mfa(
        list_of_data = res_list$best_model_mfa$list_of_data, # ここ重要。fit内にデータを入れておくor外部を参照
        mfa_fit      = res_list$best_model_mfa,
        z            = res_list$best_model_mfa$z
      )
      df1$VAF_MFA <- vaf_mfa
    } else {
      df1$VAF_MFA <- NA
    }

    if(!is.null(res_list$best_model_pca)){
      vaf_pca <- compute_VAF_mpca(
        list_of_data = res_list$best_model_pca$list_of_data,
        mpca_fit     = res_list$best_model_pca,
        z            = res_list$best_model_pca$z
      )
      df1$VAF_PCA <- vaf_pca
    } else {
      df1$VAF_PCA <- NA
    }

    df_all[[f]] <- df1
  }

  df_out <- do.call(rbind, df_all)
  df_out
}
