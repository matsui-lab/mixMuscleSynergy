#' Search for the Best (K, r) for MFA via Grid Search (Using mfa_em_fit)
#'
#' This function tries all combinations of K in Kvec and r in rvec,
#' fits a Mixture Factor Analysis (MFA) model using mfa_em_fit
#' for each combination, computes BIC, and returns the best model by BIC.
#'
#' @param list_of_data A list of data matrices, each (T_i x M).
#' @param Kvec An integer vector of candidate K values.
#' @param rvec An integer vector of candidate r values.
#' @param max_iter Maximum EM iterations.
#' @param nIterFA Sub-iterations within the factor-analyzer update in C++.
#' @param tol Convergence tolerance.
#' @param n_init Number of random initializations to try inside mfa_em_fit.
#' @param use_kmeans_init Logical; whether to also try a k-means-based init.
#' @param subject_rdim_for_kmeans Dimension for PCA-based features if use_kmeans_init=TRUE.
#'
#' @return A list with:
#'   \item{summary}{A data.frame with columns K, r, logLik, BIC.}
#'   \item{best_model_info}{A list with (K, r, model, logLik, BIC) for the best BIC.}
#'   \item{all_models}{A list of all fits, each (K, r, model, logLik, BIC).}
#'
#' @details
#' - Calls mfa_em_fit(..., n_init, use_kmeans_init, ...) for each (K, r).
#' - After fitting, we compute log-likelihood by summing over all data
#'   (you can use the built-in logLik from the fit, or re-compute).
#' - Then compute BIC = -2*logLik + p*log(N_total_rows).
#'   Here p = K*(M*r + 2*M) + (K - 1) by default (diagonal noise).
#'
#' @export
select_optimal_K_r_mfa <- function(
    list_of_data,
    Kvec = 1:5,
    rvec = 1:5,
    max_iter = 50,
    nIterFA  = 5,
    tol      = 1e-3,
    n_init   = 1,
    use_kmeans_init = TRUE,
    subject_rdim_for_kmeans = 2
){
  N <- length(list_of_data)
  if(N < 1) stop("No data in list_of_data")

  # 総タイムサンプル数 (BIC用に使う場合)
  # ここでは単に "N_total_rows = sum of T_i" とします
  N_total_rows <- sum(sapply(list_of_data, nrow))

  M <- ncol(list_of_data[[1]])  # チャネル数

  total_steps <- length(Kvec) * length(rvec)
  results_list <- vector("list", total_steps)
  df_summary   <- data.frame()

  pb <- txtProgressBar(min=0, max=total_steps, style=3)
  count <- 0

  for(K in Kvec){
    for(r in rvec){
      count <- count + 1
      setTxtProgressBar(pb, count)

      # --- 1) Fit via mfa_em_fit (複数初期値, k-means初期値など可能) ---
      fit_mfa <- mfa_em_fit(
        list_of_data = list_of_data,
        K = K,
        r = r,
        max_iter = max_iter,
        nIterFA  = nIterFA,
        tol      = tol,
        n_init   = n_init,
        use_kmeans_init = use_kmeans_init,
        subject_rdim_for_kmeans = subject_rdim_for_kmeans,
        mc_cores = 1  # ここでは並列OFF(必要に応じて)
      )

      # --- 2) logLik は fit_mfa$logLik をそのまま使う想定 ---
      loglik_val <- fit_mfa$logLik

      # --- 3) compute BIC (diagonal noise想定のパラメータ数) ---
      # 例: p = K*(M*r + 2*M) + (K-1)
      p <- K*(M*r + 2*M) + (K-1)
      bic_val <- -2*loglik_val + p*log(N_total_rows)

      results_list[[count]] <- list(
        K=K, r=r,
        model=fit_mfa,
        logLik=loglik_val,
        BIC=bic_val
      )
      df_summary <- rbind(
        df_summary,
        data.frame(K=K, r=r, logLik=loglik_val, BIC=bic_val)
      )
    }
  }
  close(pb)

  # BIC順にソート
  df_summary <- df_summary[order(df_summary$BIC), ]
  best_row <- df_summary[1, ]
  cat("=== Best model by BIC (MFA) ===\n")
  print(best_row)

  best_K <- best_row$K
  best_r <- best_row$r

  # 最良モデルを探す
  best_model_index <- which(
    sapply(results_list, function(x) x$K) == best_K &
      sapply(results_list, function(x) x$r) == best_r
  )
  best_model_info <- results_list[[ best_model_index[1] ]]  # 先頭を採用

  list(
    summary         = df_summary,
    best_model_info = best_model_info,
    all_models      = results_list
  )
}
#' Search for the Best (K, r) for Mixture PCA via Grid Search
#'
#' This function tries all combinations of K in Kvec and r in rvec,
#' fits a Mixture PCA (PPCA) model using mixture_pca_em_fit
#' for each combination, computes BIC, and returns the best model by BIC.
#'
#' @param list_of_data A list of data matrices, each (T_i x M).
#' @param Kvec An integer vector of candidate K values.
#' @param rvec An integer vector of candidate r values.
#' @param max_iter Maximum EM iterations.
#' @param nIterPCA Sub-iterations within the PCA/PPCA update in C++.
#' @param tol Convergence tolerance.
#' @param method "EM" or "closed_form" (passed to mixture_pca_em_fit).
#' @param n_init Number of random initializations to try inside mixture_pca_em_fit.
#' @param use_kmeans_init Logical; whether to also try a k-means-based init.
#' @param subject_rdim_for_kmeans Dimension for PCA-based features if use_kmeans_init=TRUE.
#'
#' @return A list with:
#'   \item{summary}{A data.frame with columns K, r, logLik, BIC.}
#'   \item{best_model_info}{A list with (K, r, model, logLik, BIC) for the best BIC.}
#'   \item{all_models}{A list of all fits, each (K, r, model, logLik, BIC).}
#'
#' @details
#' - Calls mixture_pca_em_fit(..., n_init, use_kmeans_init, ...) for each (K, r).
#' - Uses fit$logLik as the final log-likelihood.
#' - BIC formula can be adjusted if needed. By default, we might do
#'   p = K*(M*r + r + M) + (K-1), for example, if we consider each cluster having:
#'   - W is M x r (or equivalently P x D)
#'   - mu is length M
#'   - sigma2 is 1 scalar
#'   - mixing proportions (K-1)
#'
#' @export
select_optimal_K_r_mpca <- function(
    list_of_data,
    Kvec = 1:5,
    rvec = 1:5,
    max_iter = 50,
    nIterPCA = 5,
    tol      = 1e-3,
    method   = "EM",
    n_init   = 1,
    use_kmeans_init = TRUE,
    subject_rdim_for_kmeans = 2
){
  N <- length(list_of_data)
  if(N < 1) stop("No data in list_of_data")

  # 総タイムサンプル数
  N_total_rows <- sum(sapply(list_of_data, nrow))
  M <- ncol(list_of_data[[1]])

  total_steps <- length(Kvec) * length(rvec)
  results_list <- vector("list", total_steps)
  df_summary   <- data.frame()

  pb <- txtProgressBar(min=0, max=total_steps, style=3)
  count <- 0

  for(K in Kvec){
    for(r in rvec){
      count <- count + 1
      setTxtProgressBar(pb, count)

      # --- 1) Fit via mixture_pca_em_fit ---
      fit_pca <- mixture_pca_em_fit(
        list_of_data = list_of_data,
        K = K,
        r = r,
        max_iter = max_iter,
        nIterPCA = nIterPCA,
        tol      = tol,
        method   = method,
        n_init   = n_init,
        use_kmeans_init = use_kmeans_init,
        subject_rdim_for_kmeans = subject_rdim_for_kmeans,
        mc_cores = 1
      )

      # --- 2) logLik ---
      loglik_val <- fit_pca$logLik

      # --- 3) compute BIC ---
      # 例：p = K*(M*r + r + M) + (K-1)
      #     (P = M x r, mu = M, sigma2 = 1 scalar, mixing pi の (K-1))
      p <- K*(M*r + r + M) + (K-1)
      bic_val <- -2*loglik_val + p*log(N_total_rows)

      results_list[[count]] <- list(
        K=K, r=r,
        model=fit_pca,
        logLik=loglik_val,
        BIC=bic_val
      )
      df_summary <- rbind(
        df_summary,
        data.frame(K=K, r=r, logLik=loglik_val, BIC=bic_val)
      )
    }
  }
  close(pb)

  # BIC順にソート
  df_summary <- df_summary[order(df_summary$BIC), ]
  best_row <- df_summary[1, ]
  cat("=== Best model by BIC (Mixture PCA) ===\n")
  print(best_row)

  best_K <- best_row$K
  best_r <- best_row$r

  # 最良モデル
  best_model_index <- which(
    sapply(results_list, function(x) x$K) == best_K &
      sapply(results_list, function(x) x$r) == best_r
  )
  best_model_info <- results_list[[ best_model_index[1] ]]

  list(
    summary         = df_summary,
    best_model_info = best_model_info,
    all_models      = results_list
  )
}
