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
