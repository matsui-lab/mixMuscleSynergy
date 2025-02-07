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
