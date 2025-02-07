#' Fit a Mixture PCA Model (Using mpcaTimeseriesCpp)
#'
#' This function calls the C++ function \code{mpcaTimeseriesCpp()} to perform a
#' Mixture PCA (PPCA) EM algorithm for fixed \code{K} and \code{r}. 
#' After that, it computes the final log-likelihood and responsibilities in R,
#' returning a structure close to the original "simple" approach with \code{(P, D, Psi, ...)}.
#'
#' @param list_of_data A list of matrices (each \code{T_i x M}).
#' @param K Number of clusters.
#' @param r Number of principal components.
#' @param max_iter Maximum EM iterations for the C++ routine.
#' @param nIterPCA Sub-iterations for updating each cluster's PPCA parameters (\code{pcaEMupdateCpp}).
#' @param tol Convergence tolerance for \code{mpcaTimeseriesCpp}.
#' @param method Either \code{"EM"} or \code{"closed_form"}, indicating which cluster-level routine
#'   to use inside C++ (\code{pcaEMupdateCpp} or \code{pcaClosedFormCpp}).
#'
#' @return A list with elements:
#' \item{z}{Hard cluster assignments (1..K).}
#' \item{pi}{Mixing proportions.}
#' \item{mu}{List of length \code{K}, each a mean vector.}
#' \item{P}{List of length \code{K}, each an \code{(M x r)} matrix of principal directions.}
#' \item{D}{List of length \code{K}, each an \code{(r x r)} diagonal matrix. 
#'   (For now, we create a trivial diagonal from the norm of \code{W}'s columns; see details.)}
#' \item{Psi}{List of length \code{K}, each an \code{(M x M)} diagonal matrix representing residual noise \eqn{\sigma^2 I}.}
#' \item{logLik}{Final log-likelihood over all data.}
#' \item{resp}{An \code{N x K} matrix of mixture responsibilities in R.}
#'
#' @details
#' Internally calls \code{mpcaTimeseriesCpp()} which performs an EM loop at the cluster level, 
#' grouping time-series that belong to each cluster (hard assignment). In each M-step, it either uses
#' \code{pcaEMupdateCpp} (if \code{method="EM"}) or \code{pcaClosedFormCpp} (if \code{method="closed_form"})
#' to get a single-cluster PPCA solution (mean, W, sigma2). In this R function, we interpret
#' \code{W} as \code{P}, store an approximate diagonal matrix \code{D}, and create a diagonal
#' \code{Psi} matrix from \code{sigma2 I}. We also compute final log-likelihood and responsibilities
#' by evaluating the mixture model in R.
#'
#' If you need a more precise \code{D} (like exact singular values) or different structure, adjust the code below.
#'
#' @examples
#' \dontrun{
#' # Suppose we have a list_of_data with M=8 columns, K=3 clusters, r=2 PCs:
#' fit <- mixture_pca_em_fit(list_of_data, K=3, r=2, max_iter=30, nIterPCA=10, tol=1e-4)
#' print(fit$logLik)
#' head(fit$z)
#' }
#'
#' @export
mixture_pca_em_fit <- function(list_of_data, 
                               K, 
                               r, 
                               max_iter = 50, 
                               nIterPCA = 20,
                               tol      = 1e-3,
                               method   = "EM")
{
  # 1) Call the C++ function
  fit_cpp <- mpcaTimeseriesCpp(
    list_of_data = list_of_data,
    K = K,
    r = r,
    max_iter  = max_iter,
    nIterPCA  = nIterPCA,
    tol       = tol,
    method    = method  # "EM" or "closed_form"
  )
  
  # fit_cpp has: z, W, mu, sigma2, pi 
  #  - W is a list of length K, each (M x r)
  #  - mu is a list of length K, each M-vector
  #  - sigma2 is a numeric vector (length K)
  #  - pi is a numeric vector (length K)
  
  N <- length(list_of_data)
  if(N < 1){
    stop("No data in list_of_data.")
  }
  # Retrieve them for convenience
  z_cpp      <- fit_cpp$z
  W_list_cpp <- fit_cpp$W
  mu_list_cpp<- fit_cpp$mu
  sigma2_vec <- fit_cpp$sigma2
  pi_vec     <- fit_cpp$pi
  
  # 2) Build final "P", "D", "Psi" from W, sigma2
  #    We'll interpret W = (M x r). 
  #    We'll store P = normalized columns of W,
  #    and D = diag of column scales 
  #    Then Psi_k = sigma2_k * Identity(M)
  
  # Prepare placeholders
  P_list  <- vector("list", K)
  D_list  <- vector("list", K)
  Psi_list<- vector("list", K)
  for(k2 in seq_len(K)){
    W_k <- W_list_cpp[[k2]]
    # compute column norms => principal directions
    #  e.g. col_j norm => place in D
    M <- nrow(W_k)
    r_ <- ncol(W_k)
    col_scales <- numeric(r_)
    # row-by-column loop or use apply
    for(j in seq_len(r_)){
      col_scales[j] <- sqrt(sum(W_k[,j]^2))
      if(col_scales[j] < 1e-12) col_scales[j] <- 1e-12
    }
    # P_k columns = W_k columns / scale_j
    P_k <- W_k
    for(j in seq_len(r_)){
      P_k[,j] <- W_k[,j] / col_scales[j]
    }
    # D_k = diag(col_scales)
    D_k <- diag(col_scales, r_, r_)
    
    # Psi_k = sigma2[k2] * Identity(M)
    sig2_k <- sigma2_vec[k2]
    Psi_k  <- diag(sig2_k, M)
    
    P_list[[k2]]  <- P_k
    D_list[[k2]]  <- D_k
    Psi_list[[k2]]<- Psi_k
  }
  
  # 3) Compute final responsibilities & logLik
  #    Sigma_k = P_k (D_k^2) P_k^T + Psi_k is not strictly correct for the code above,
  #    but let's "approximate" to match an MFA-like structure. Actually in code we used W W^T + sigma2 I.
  #    We'll do W_k = P_k D_k => 
  #    => Sigma_k = (P_k D_k) (D_k^T P_k^T) + sigma2 I
  #       = P_k (D_k^2) P_k^T + sigma2 I
  # We'll do that. Then we do the mixture log-likelihood subject by subject.
  
  # Build Sigmas
  Sigma_list <- vector("list", K)
  for(k2 in seq_len(K)){
    P_k  <- P_list[[k2]]
    D_k  <- D_list[[k2]]
    M_k  <- nrow(P_k)
    Sigma_k <- P_k %*% (D_k^2) %*% t(P_k) 
    diag(Sigma_k) <- diag(Sigma_k) + sigma2_vec[k2]
    Sigma_list[[k2]] <- Sigma_k
  }
  
  resp <- matrix(0, nrow=N, ncol=K)
  logLik_val <- 0
  
  for(i in seq_len(N)){
    X_i <- list_of_data[[i]]
    T_i <- nrow(X_i)
    logvals <- numeric(K)
    for(k2 in seq_len(K)){
      mu_k  <- mu_list_cpp[[k2]]
      Sig_k <- Sigma_list[[k2]]
      dens_t <- mvtnorm::dmvnorm(X_i, mean=mu_k, sigma=Sig_k, log=TRUE)
      sumLog <- sum(dens_t)
      logvals[k2] <- log(pi_vec[k2] + 1e-16) + sumLog
    }
    # log-sum-exp
    m0 <- max(logvals)
    ll_i <- m0 + log(sum(exp(logvals - m0)))
    logLik_val <- logLik_val + ll_i
    
    # responsibilities
    for(k2 in seq_len(K)){
      resp[i,k2] <- exp(logvals[k2] - ll_i)
    }
  }
  
  # 4) Construct final result
  #    We name them as doc specified: z, pi, mu, P, D, Psi, logLik, resp
  res <- list(
    z      = z_cpp,
    pi     = pi_vec,
    mu     = mu_list_cpp,
    P      = P_list,
    D      = D_list,
    Psi    = Psi_list,
    logLik = logLik_val,
    resp   = resp
  )
  
  res
}
