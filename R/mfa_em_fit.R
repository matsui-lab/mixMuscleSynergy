#' Fit a Mixture Factor Analysis Model (Using mfaTimeseriesCpp)
#'
#' This function calls the C++ function \code{mfaTimeseriesCpp()} to perform a
#' Mixture Factor Analysis EM algorithm for fixed \code{K} and \code{r}, and then
#' calculates the final log-likelihood and responsibilities in R.
#'
#' @param list_of_data A list of matrices (each \code{(T_i x M)}) to be modeled.
#' @param K Number of clusters.
#' @param r Factor dimension.
#' @param max_iter Maximum EM iterations for the C++ routine.
#' @param nIterFA Number of sub-iterations for the factor analyzer update
#'   (passed to C++). Defaults to 20.
#' @param tol Convergence tolerance for log-likelihood difference in the C++ code.
#'
#' @return A list with elements:
#' \item{z}{Cluster assignment (hard) for each subject, from the C++ result.}
#' \item{pi}{Cluster mixing proportions.}
#' \item{mu}{List of cluster means (length K).}
#' \item{Lambda}{List of factor loading matrices (length K).}
#' \item{Psi}{List of diagonal noise matrices (length K).}
#' \item{logLik}{A numeric scalar for the total log-likelihood (computed in R).}
#' \item{resp}{An \code{(N x K)} matrix of responsibilities, computed in R.}
#'
#' @details
#' The heavy-lifting EM steps occur in \code{mfaTimeseriesCpp()}, which updates parameters
#' by grouping data for each cluster. After that completes, this function performs a single
#' pass in R to compute the final responsibilities and log-likelihood by evaluating the mixture
#' log-density for each subject's data.
#'
#' @examples
#' \dontrun{
#' # Suppose we have a list_of_data with M columns, and we want K=3, r=2:
#' fit <- mfa_em_fit(list_of_data, K=3, r=2, max_iter=50, nIterFA=20, tol=1e-4)
#' print(fit$logLik)
#' head(fit$z)
#' }
#'
#' @export
mfa_em_fit <- function(list_of_data, K, r,
                       max_iter = 50,
                       nIterFA  = 20,
                       tol      = 1e-3)
{
  # 1) Call the C++ function to do the EM steps
  fit_cpp <- mfaTimeseriesCpp(
    list_of_data = list_of_data,
    K = K,
    r = r,
    max_iter = max_iter,
    nIterFA  = nIterFA,
    tol      = tol
  )
  
  # fit_cpp should have: z, Lambda, mu, Psi, pi
  # but not logLik or resp.
  
  # 2) We compute final logLik and resp in R for convenience
  N <- length(list_of_data)
  K_c <- length(fit_cpp$Lambda)  # should match K
  if(K_c != K){
    warning("Number of clusters in C++ result does not match K? Possibly K is correct.")
  }
  
  # Precompute cluster covariance matrices
  Sigma_list <- vector("list", K_c)
  for(k2 in seq_len(K_c)){
    Lambda_k <- fit_cpp$Lambda[[k2]]
    Psi_k    <- fit_cpp$Psi[[k2]]
    Sigma_list[[k2]] <- Lambda_k %*% t(Lambda_k) + Psi_k
  }
  pi_vec <- fit_cpp$pi
  
  # We'll accumulate log-likelihood over all subjects
  logLik_val <- 0
  
  # We'll also build an (N x K) responsibility matrix
  resp <- matrix(0, nrow=N, ncol=K_c)
  
  # For each subject i, sum over time steps
  for(i in seq_len(N)){
    X_i <- list_of_data[[i]]
    T_i <- nrow(X_i)
    logvals <- numeric(K_c)
    
    for(k2 in seq_len(K_c)){
      mu_k  <- fit_cpp$mu[[k2]]
      sig_k <- Sigma_list[[k2]]
      # row-wise log density, sum
      dens_t <- mvtnorm::dmvnorm(X_i, mean=mu_k, sigma=sig_k, log=TRUE)
      sumLog <- sum(dens_t)
      logvals[k2] <- log(pi_vec[k2] + 1e-16) + sumLog
    }
    # log-sum-exp
    m0 <- max(logvals)
    ll_i <- m0 + log(sum(exp(logvals - m0)))
    logLik_val <- logLik_val + ll_i
    
    # responsibilities
    for(k2 in seq_len(K_c)){
      resp[i,k2] <- exp(logvals[k2] - ll_i)
    }
  }
  
  # Prepare final result
  res <- list(
    z      = fit_cpp$z,
    pi     = fit_cpp$pi,
    mu     = fit_cpp$mu,
    Lambda = fit_cpp$Lambda,
    Psi    = fit_cpp$Psi,
    logLik = logLik_val,
    resp   = resp
  )
  
  res
}
