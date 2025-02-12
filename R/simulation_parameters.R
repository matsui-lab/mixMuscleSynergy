#' Read All Results and Produce Tables/Figures
#'
#' This function demonstrates how to read all .rds files from a directory,
#' merge them into a single \code{data.frame}, generate a couple of tables (Table 1, Table 2),
#' and produce various figures (Figure 1, Figure 2, etc.) saving them to files.
#'
#' @param result_dir Directory containing .rds result files (default "results").
#'
#' @return Invisibly returns the merged data frame of all results.
#'
#' @details
#' The user can modify the conditions in \code{\link{make_table1}} and \code{\link{make_table2}}
#' to match their actual experiment ranges. Figures are generated via \code{ggsave}.
#'
#' @examples
#' \dontrun{
#' df_all <- read_and_visualize_all("results")
#' }
#'
#' @export
read_and_visualize_all <- function(result_dir="results"){
  df <- read_all_results(result_dir=result_dir)

  df_table1 <- make_table1()
  cat("=== Table 1: Simulation Conditions ===\n")
  print(df_table1)
  # knitr::kable(df_table1, caption="Table 1: Conditions")

  df_table2 <- make_table2(df)
  cat("=== Table 2: Summary of ARI etc. ===\n")
  print(df_table2)
  # kable(df_table2, caption="Table 2: Summary")

  p1 <- plot_fig1_ARI(df)
  ggsave("figure1_ARI_sep_noise.pdf", p1, width=8, height=6)

  p2 <- plot_fig2_BIC_heatmap(df)
  ggsave("figure2_BIC_heatmap.pdf", p2, width=6, height=5)

  if(nrow(df) > 1){
    p3 <- plot_fig3_init_boxplot(df)
    ggsave("figure3_init_boxplot.pdf", p3, width=6, height=4)
  }

  # (Figure 4) requires an MFA fit object. Not shown by default.
  # see the example in plot_fig4_loadings_heatmap

  invisible(df)
}


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
#' @param K_candidates Vector of candidate K values to try.
#' @param r_candidates Vector of candidate r values to try.
#' @param max_iter Max iteration for EM.
#' @param n_rep_init Number of random initial seeds for each (K,r) in search.
#' @param seed_data_base Base seed for data simulation.
#' @param mc.cores Number of parallel cores to use.
#' @param output_dir Directory to save each result as \code{result_K{...}_r{...}...rds}.
#'
#' @details
#' For each condition, the function calls \code{\link{simulate_mixture_data}} to create data,
#' then loops over candidate K, r in \code{K_candidates} and \code{r_candidates}, fits MFA and MPCA multiple times,
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
    K_candidates = 1:5,
    r_candidates = 1:5,
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


#' Main Example: Run Parallel MFA/MPCA Simulation and Save Each Condition
#'
#' This function calls \code{\link{simulate_and_run_parallel_save}} with some default arguments,
#' prints the merged results, and generates a simple figure comparing ARI_MFA and ARI_PCA.
#'
#' @return The merged data frame of results across all conditions.
#'
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
    mc.cores = 2,       # e.g. 2コア
    output_dir = "results"
  )

  # show
  print(head(df_result, 10))

  # quick plot
  df_long <- df_result %>%
    select(K_true, r_true, N, M, sep, noise, ARI_MFA, ARI_PCA) %>%
    gather(key="Method", value="ARI", ARI_MFA, ARI_PCA)

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
#' This function generates synthetic time-series data according to a mixture of factor (or PCA) model.
#' Each subject is assigned to one of \code{K} clusters, each with factor loadings (or principal vectors),
#' mean offsets, and diagonal noise.
#'
#' @param N Number of subjects (time-series).
#' @param K True number of clusters in simulation.
#' @param r True dimension of factors or principal components.
#' @param M Number of observed channels (e.g. muscle EMG channels).
#' @param T_each Length of each time-series per subject.
#' @param cluster_sep How far to separate cluster means (and possibly scale the factor loadings).
#' @param noise_scale Scaling of noise diagonal in \code{\link{Psi_list}}.
#' @param seed Random seed for reproducibility.
#'
#' @return A list with:
#' \item{list_of_data}{A list of length \code{N}, each an \code{(T_each x M)} matrix.}
#' \item{z_true}{A length-\code{N} vector of true cluster assignments (1..K).}
#' @examples
#' \dontrun{
#' sim <- simulate_mixture_data(
#'   N=50, K=3, r=4, M=8, T_each=100,
#'   cluster_sep=1.0, noise_scale=1.0, seed=123
#' )
#' length(sim$list_of_data)  # 50
#' sim$z_true[1:10]
#' }
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




#' Compute BIC for MFA or Mixture PCA
#'
#' This function computes a simple BIC measure given a log-likelihood and
#' a rough count of model parameters, assuming diagonal noise.
#'
#' @param logLik The log-likelihood value.
#' @param K Number of clusters.
#' @param r Factor or principal component dimension.
#' @param M Number of observed channels.
#' @param N_total_rows Total number of rows across all data (used for \code{log(N)}).
#'
#' @details
#' The number of parameters is counted as:
#' \deqn{K*(M*r + 2M) + (K-1).}
#'
#' @return A numeric value for BIC.
#'
#' @examples
#' bic_val <- compute_bic(-1234.5, K=3, r=4, M=8, N_total_rows=500)
#' @export
compute_bic <- function(logLik, K, r, M, N_total_rows){
  num_params <- K*(M*r + 2*M) + (K-1)
  -2*logLik + num_params*log(N_total_rows)
}


#' Read All Results from .rds Files
#'
#' This function searches a directory for \code{.rds} files containing one-condition results
#' (each is typically a \code{data.frame} of 1 or more rows), reads them, and merges them
#' into a single data frame.
#'
#' @param result_dir A character string specifying the directory to look for \code{.rds} files.
#' @return A combined \code{data.frame} with rows from all the loaded files.
#'
#' @details
#' Each \code{.rds} file is assumed to contain a \code{data.frame} with columns such as
#' \code{K_true, r_true, N, M, sep, noise, ARI_MFA, BIC_MFA, ...}. This function simply does a
#' \code{rbind} of all those data frames. If no \code{.rds} files are found, it raises an error.
#'
#' @examples
#' \dontrun{
#' df_merged <- read_all_results("results/")
#' head(df_merged)
#' }
#'
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
#' Given a result data frame with columns like \code{ARI_MFA, ARI_PCA, sep, noise},
#' this function reshapes them into long format (MFA vs. MPCA) and plots ARI vs. separation,
#' color-coded by noise, using \code{ggplot2}.
#'
#' @param df A \code{data.frame} with columns \code{K_true, r_true, N, M, sep, noise, ARI_MFA, ARI_PCA}.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' p <- plot_fig1_ARI(df_merged)
#' ggsave("figure1_ARI_sep_noise.pdf", p)
#' }
#'
#' @export
plot_fig1_ARI <- function(df){
  df_long <- df %>%
    dplyr::select(K_true, r_true, N, M, sep, noise, ARI_MFA, ARI_PCA) %>%  # N, M を追加
    tidyr::pivot_longer(cols = c("ARI_MFA", "ARI_PCA"),
                        names_to = "Method", values_to = "ARI")

  p <- ggplot(df_long, aes(x = factor(sep), y = ARI, color = factor(noise),
                           group = interaction(Method, noise))) +
    geom_point(position = position_dodge(width = 0.3)) +
    geom_line(position = position_dodge(width = 0.3)) +
    facet_wrap(~ Method + K_true + r_true + N + M, labeller = label_both) +  # N, M を使用
    theme_bw() +
    labs(x = "Separation (sep)", y = "ARI", color = "Noise",
         title = "(Figure 1) ARI vs sep/noise (MFA, MPCA)")

  return(p)
}



#' Plot BIC Heatmap for a Subset (Figure 2)
#'
#' This function filters the data for a specific condition, e.g. \code{N=50, M=8, sep=1.0, noise=1.0},
#' and plots a heatmap of \code{BIC_MFA} with \code{r_true} on the x-axis and \code{K_true} on the y-axis.
#'
#' @param df A \code{data.frame} with columns including \code{BIC_MFA, K_true, r_true}.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' p <- plot_fig2_BIC_heatmap(df_merged)
#' ggsave("figure2_BIC_heatmap.pdf", p)
#' }
#'
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
#' If the data frame contains multiple runs with the same parameters but different
#' initial seeds, this function can show a boxplot of ARI distribution for MFA and MPCA.
#'
#' @param df_multiple A \code{data.frame} with columns \code{ARI_MFA, ARI_PCA}.
#'
#' @return A \code{ggplot} object.
#'
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
#' Given an MFA fit object containing \code{Lambda}, this function plots a heatmap
#' for the loadings in a specific cluster.
#'
#' @param mfa_fit A fitted MFA model with \code{$Lambda}.
#' @param cluster_id Which cluster's loadings to plot (1..K).
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' fit_mfa <- readRDS("some_fit.rds")
#' p <- plot_fig4_loadings_heatmap(fit_mfa, cluster_id=1)
#' ggsave("figure4_loadings_heatmap.pdf", p)
#' }
#'
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
#' This function returns a small \code{data.frame} describing the parameter ranges
#' used in the simulation (K, r, N, M, etc.). Typically printed via \code{kable()} in a vignette or paper.
#'
#' @return A \code{data.frame} with columns \code{Parameter} and \code{Values}.
#'
#' @examples
#' \dontrun{
#' df_table1 <- make_table1()
#' knitr::kable(df_table1, caption="Table 1: Conditions")
#' }
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
#' This function takes the merged result data frame (e.g., from \code{\link{read_all_results}}),
#' filters for a specific condition, and computes group-wise means and standard deviations
#' of ARI and BIC for MFA and MPCA.
#'
#' @param df A \code{data.frame} containing columns like \code{K_true, r_true, N, M, sep, noise, ARI_MFA, ARI_PCA, BIC_MFA, BIC_PCA}.
#'
#' @return A \code{data.frame} summarizing ARI/BIC means and sds, grouped by \code{(K_true, r_true, N, M)}.
#'
#' @examples
#' \dontrun{
#' df_merged <- read_all_results("results")
#' df_table2 <- make_table2(df_merged)
#' knitr::kable(df_table2, caption="Table 2: Summary of ARI, BIC, etc.")
#' }
#'
#' @export
make_table2 <- function(df){
  df_subset <- df %>%
    dplyr::filter(abs(sep - 1.0)<1e-9, abs(noise-1.0)<1e-9)

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

