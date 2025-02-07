#' Search for the Best (K, r) for MFA via Grid Search
#'
#' This function tries all combinations of \code{K} in \code{Kvec} and
#' \code{r} in \code{rvec}, fits a Mixture Factor Analysis (MFA) model
#' using \code{\link{mfa_em_fit}} (or your wrapper) for each combination,
#' computes BIC, and returns the best model by BIC.
#'
#' @param list_of_data A list of data matrices, each \code{(T_i x M)}.
#' @param Kvec An integer vector of candidate \code{K} values.
#' @param rvec An integer vector of candidate \code{r} values.
#' @param max_iter Maximum EM iterations for \code{mfa_em_fit}.
#' @param nIterFA Sub-iterations within the factor-analyzer step.
#' @param tol Convergence tolerance for \code{mfa_em_fit}.
#'
#' @return A list with:
#' \item{summary}{A data frame with columns \code{K, r, logLik, BIC}.}
#' \item{best_model_info}{A list with \code{(K, r, model, logLik, BIC)} for the best BIC.}
#' \item{all_models}{A list of all fits, each \code{(K, r, model, logLik, BIC)}.}
#'
#' @details
#' Internally calls \code{\link{mfa_em_fit}} for each \code{(K, r)}. Then we compute log-likelihood
#' via \code{\link{compute_logLik_mfa}} and BIC via \code{\link{compute_BIC_mfa}}.
#'
#' @examples
#' \dontrun{
#' out <- select_optimal_K_r_mfa(my_data, Kvec=2:4, rvec=1:3, max_iter=50, nIterFA=10, tol=1e-4)
#' names(out)
#' # "summary", "best_model_info", "all_models"
#' out$summary
#' out$best_model_info
#' }
#'
#' @export
select_optimal_K_r_mfa <- function(list_of_data,
                                   Kvec = 1:5,
                                   rvec = 1:5,
                                   max_iter=50,
                                   nIterFA=5,
                                   tol=1e-3)
{
  N <- length(list_of_data)
  if(N < 1) stop("No data in list_of_data")
  
  total_steps <- length(Kvec) * length(rvec)
  results_list <- list()
  df_summary <- data.frame()
  
  pb <- txtProgressBar(min=0, max=total_steps, style=3)  # progress bar
  count <- 0
  
  for(K in Kvec){
    for(r in rvec){
      count <- count + 1
      setTxtProgressBar(pb, count)
      
      # 1) Fit via MFA EM
      fit <- mfa_em_fit(
        list_of_data = list_of_data,
        K = K,
        r = r,
        max_iter = max_iter,
        tol      = tol,
        nIterFA  = nIterFA
      )
      
      # 2) compute logLik
      loglik_val <- compute_logLik_mfa(list_of_data, fit)
      
      # 3) compute BIC
      M <- ncol(list_of_data[[1]])
      bic_val <- compute_BIC_mfa(loglik_val, K, r, M, N)
      
      results_list[[count]] <- list(K=K, r=r, model=fit,
                                    logLik=loglik_val, BIC=bic_val)
      df_summary <- rbind(df_summary, data.frame(K=K, r=r,
                                                 logLik=loglik_val,
                                                 BIC=bic_val))
    }
  }
  close(pb)
  
  df_summary <- df_summary[order(df_summary$BIC), ]
  best_row <- df_summary[1,]
  cat("=== Best model by BIC: ===\n")
  print(best_row)
  
  best_K <- best_row$K
  best_r <- best_row$r
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

#' Search for the Best (K, r) for MFA via Grid Search
#'
#' This function tries all combinations of \code{K} in \code{Kvec} and
#' \code{r} in \code{rvec}, fits a Mixture Factor Analysis (MFA) model
#' using \code{\link{mfa_em_fit}} (or your wrapper) for each combination,
#' computes BIC, and returns the best model by BIC.
#'
#' @param list_of_data A list of data matrices, each \code{(T_i x M)}.
#' @param Kvec An integer vector of candidate \code{K} values.
#' @param rvec An integer vector of candidate \code{r} values.
#' @param max_iter Maximum EM iterations for \code{mfa_em_fit}.
#' @param nIterFA Sub-iterations within the factor-analyzer step.
#' @param tol Convergence tolerance for \code{mfa_em_fit}.
#'
#' @return A list with:
#' \item{summary}{A data frame with columns \code{K, r, logLik, BIC}.}
#' \item{best_model_info}{A list with \code{(K, r, model, logLik, BIC)} for the best BIC.}
#' \item{all_models}{A list of all fits, each \code{(K, r, model, logLik, BIC)}.}
#'
#' @details
#' Internally calls \code{\link{mfa_em_fit}} for each \code{(K, r)}. Then we compute log-likelihood
#' via \code{\link{compute_logLik_mfa}} and BIC via \code{\link{compute_BIC_mfa}}.
#'
#' @examples
#' \dontrun{
#' out <- select_optimal_K_r_mfa(my_data, Kvec=2:4, rvec=1:3, max_iter=50, nIterFA=10, tol=1e-4)
#' names(out)
#' # "summary", "best_model_info", "all_models"
#' out$summary
#' out$best_model_info
#' }
#'
#' @export
select_optimal_K_r_mfa <- function(list_of_data,
                                   Kvec = 1:5,
                                   rvec = 1:5,
                                   max_iter=50,
                                   nIterFA=5,
                                   tol=1e-3)
{
  N <- length(list_of_data)
  if(N < 1) stop("No data in list_of_data")
  
  total_steps <- length(Kvec) * length(rvec)
  results_list <- list()
  df_summary <- data.frame()
  
  pb <- txtProgressBar(min=0, max=total_steps, style=3)  # progress bar
  count <- 0
  
  for(K in Kvec){
    for(r in rvec){
      count <- count + 1
      setTxtProgressBar(pb, count)
      
      # 1) Fit via MFA EM
      fit <- mfa_em_fit(
        list_of_data = list_of_data,
        K = K,
        r = r,
        max_iter = max_iter,
        tol      = tol,
        nIterFA  = nIterFA
      )
      
      # 2) compute logLik
      loglik_val <- compute_logLik_mfa(list_of_data, fit)
      
      # 3) compute BIC
      M <- ncol(list_of_data[[1]])
      bic_val <- compute_BIC_mfa(loglik_val, K, r, M, N)
      
      results_list[[count]] <- list(K=K, r=r, model=fit,
                                    logLik=loglik_val, BIC=bic_val)
      df_summary <- rbind(df_summary, data.frame(K=K, r=r,
                                                 logLik=loglik_val,
                                                 BIC=bic_val))
    }
  }
  close(pb)
  
  df_summary <- df_summary[order(df_summary$BIC), ]
  best_row <- df_summary[1,]
  cat("=== Best model by BIC: ===\n")
  print(best_row)
  
  best_K <- best_row$K
  best_r <- best_row$r
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