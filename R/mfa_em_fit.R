#' Single-Run Mixture Factor Analysis with a Preset Initial Assignment
#'
#' This helper function calls the C++ routine \code{mfaTimeseriesCpp()} with a user-specified
#' cluster assignment (\code{z_init}) and then computes the final log-likelihood and
#' responsibilities in R.
#'
#' @param list_of_data A list of \code{(T_i x M)} matrices, one for each subject.
#' @param K Number of clusters.
#' @param r Factor dimension in each cluster.
#' @param z_init An integer vector of length \code{N} (where \code{N=length(list_of_data)}),
#'   giving the initial cluster label (1..K) for each subject.
#' @param max_iter Maximum EM iterations in C++.
#' @param nIterFA Number of sub-iterations for the factor analyzer update in C++.
#' @param tol Convergence tolerance used in C++.
#'
#' @return A list with elements:
#' \describe{
#'   \item{\code{z}}{Hard cluster assignments, length \code{N}.}
#'   \item{\code{pi}}{Cluster mixing proportions, length \code{K}.}
#'   \item{\code{mu}}{List of length \code{K}, each a mean vector (\code{M x 1}).}
#'   \item{\code{Lambda}}{List of length \code{K}, each an \code{M x r} factor loading matrix.}
#'   \item{\code{Psi}}{List of length \code{K}, each an \code{M x M} diagonal noise matrix.}
#'   \item{\code{logLik}}{Final log-likelihood computed in R.}
#'   \item{\code{resp}}{An \code{N x K} matrix of responsibilities.}
#' }
#'
#' @examples
#' \dontrun{
#' # Suppose we have list_of_data, K=3, r=2:
#' z_init <- sample.int(3, size=length(list_of_data), replace=TRUE)
#' fit_one <- mfa_em_fit_cpp_singleInit(list_of_data, K=3, r=2,
#'                                      z_init=z_init)
#' fit_one$logLik
#' }
#'
#' @export
mfa_em_fit_cpp_singleInit <- function(list_of_data,
                                      K,
                                      r,
                                      z_init,
                                      max_iter = 50,
                                      nIterFA  = 20,
                                      tol      = 1e-3)
{
  # 1) Call the updated C++ function (which accepts z_init)
  fit_cpp <- mfaTimeseriesCpp(
    list_of_data = list_of_data,
    K = K,
    r = r,
    max_iter = max_iter,
    nIterFA  = nIterFA,
    tol      = tol,
    z_init   = z_init
  )

  # 2) Compute final log-likelihood & responsibilities in R
  N <- length(list_of_data)
  Kc <- length(fit_cpp$Lambda)
  if(Kc != K){
    stop("C++ output mismatch: length(Lambda) != K.")
  }

  # Build Sigma_k
  Sigma_list <- vector("list", Kc)
  for(k2 in seq_len(Kc)){
    Lambda_k <- fit_cpp$Lambda[[k2]]
    Psi_k    <- fit_cpp$Psi[[k2]]
    Sig_k    <- Lambda_k %*% t(Lambda_k) + Psi_k
    Sigma_list[[k2]] <- Sig_k
  }
  pi_vec <- fit_cpp$pi

  # Accumulate log-likelihood
  logLik_val <- 0
  resp <- matrix(0, nrow=N, ncol=Kc)

  for(i in seq_len(N)){
    Xi <- list_of_data[[i]]
    T_i <- nrow(Xi)
    logvals <- numeric(Kc)

    for(k2 in seq_len(Kc)){
      mu_k  <- fit_cpp$mu[[k2]]
      Sig_k <- Sigma_list[[k2]]
      dens_t <- mvtnorm::dmvnorm(Xi, mean=mu_k, sigma=Sig_k, log=TRUE)
      sumLog <- sum(dens_t)
      logvals[k2] <- log(pi_vec[k2] + 1e-16) + sumLog
    }

    # log-sum-exp
    m0 <- max(logvals)
    li <- m0 + log(sum(exp(logvals - m0)))
    logLik_val <- logLik_val + li

    # responsibilities
    for(k2 in seq_len(Kc)){
      resp[i,k2] <- exp(logvals[k2] - li)
    }
  }

  # Final output
  out <- list(
    z      = fit_cpp$z,
    pi     = fit_cpp$pi,
    mu     = fit_cpp$mu,
    Lambda = fit_cpp$Lambda,
    Psi    = fit_cpp$Psi,
    logLik = logLik_val,
    resp   = resp
  )
  out
}


#' Fit a Mixture Factor Analysis Model (with optional multi-initialization)
#'
#' This function calls the C++ function \code{mfaTimeseriesCpp()} to perform a Mixture
#' Factor Analysis EM algorithm for fixed \code{K} and \code{r}. By default, it runs
#' a single pass with the internal initialization from C++. However, if you specify
#' multiple initial attempts (via \code{n_init>1}) and/or \code{use_kmeans_init=TRUE},
#' this function will try several different initial cluster assignments (in parallel
#' using \code{mclapply}), then return the best solution (maximizing the final log-likelihood).
#'
#' @param list_of_data A list of matrices (each \code{(T_i x M)}) to be modeled.
#' @param K Number of clusters.
#' @param r Factor dimension.
#' @param max_iter Maximum EM iterations for \code{mfaTimeseriesCpp}.
#' @param nIterFA Number of sub-iterations for the factor analyzer update in C++.
#' @param tol Convergence tolerance for \code{mfaTimeseriesCpp}.
#' @param n_init Number of random initial assignments to try (in addition to the default
#'   single-run or k-means if requested). Defaults to \code{1}.
#' @param use_kmeans_init If \code{TRUE}, we also run one initialization where we
#'   extract subject features (via PCA) and apply \code{kmeans} to get a cluster assignment.
#'   Defaults to \code{FALSE}.
#' @param subject_rdim_for_kmeans PCA dimension for extracting subject-level features,
#'   used only if \code{use_kmeans_init=TRUE}. Default is \code{r}.
#' @param mc_cores Number of cores for parallel execution via \code{mclapply}.
#'   Defaults to \code{1} (no parallel).
#'
#' @return A list with the same elements as a single run: \code{z, pi, mu, Lambda, Psi, logLik, resp}.
#'   It corresponds to the best solution (i.e. highest log-likelihood) among all tried inits.
#'
#' @details
#' The single-run approach in C++ already does a built-in initialization if \code{z_init}
#' is not provided. Therefore, if \code{n_init=1} and \code{use_kmeans_init=FALSE}, we just call
#' \code{mfaTimeseriesCpp} once (like the original design).
#' Otherwise:
#' \enumerate{
#'   \item If \code{use_kmeans_init=TRUE}, we do one run where we assign clusters by k-means
#'         on some subject-level features (extracted by PCA).
#'   \item We also sample \code{n_init} random assignments (each subject assigned to a random cluster).
#'   \item For each assignment, we call \code{\link{mfa_em_fit_cpp_singleInit}}.
#'   \item We compare all results' \code{logLik} and keep the best one.
#' }
#'
#' This provides a more robust initialization scheme to avoid poor local maxima.
#'
#' @examples
#' \dontrun{
#' # Suppose we have a list_of_data for N=100 subjects, each (T_i x M=8).
#' # We want to fit K=3, r=2, and try 5 random inits plus one k-means init:
#' fit <- mfa_em_fit(
#'   list_of_data, K=3, r=2,
#'   max_iter=50, nIterFA=20, tol=1e-3,
#'   n_init=5,
#'   use_kmeans_init=TRUE,
#'   subject_rdim_for_kmeans=2,
#'   mc_cores=2
#' )
#' print(fit$logLik)
#' head(fit$z)
#' }
#'
#' @export
mfa_em_fit <- function(
    list_of_data,
    K,
    r,
    max_iter = 50,
    nIterFA  = 20,
    tol      = 1e-3,
    n_init   = 1,
    use_kmeans_init = FALSE,
    subject_rdim_for_kmeans = r,
    mc_cores = 1
)
{
  # If n_init=1 and use_kmeans_init=FALSE => single-run original approach
  # i.e., just call the C++ code once, relying on its built-in (i%K) init
  if(n_init == 1 && !use_kmeans_init){
    # (original style) => call mfaTimeseriesCpp with no z_init
    fit_cpp <- mfaTimeseriesCpp(
      list_of_data = list_of_data,
      K = K,
      r = r,
      max_iter = max_iter,
      nIterFA  = nIterFA,
      tol      = tol
      # no z_init => it uses the default (i%K)+1
    )
    # Now compute final logLik & resp in R (similar to your original post-fitting code)
    N <- length(list_of_data)
    Sigma_list <- vector("list", K)
    for(k2 in seq_len(K)){
      Lambda_k <- fit_cpp$Lambda[[k2]]
      Psi_k    <- fit_cpp$Psi[[k2]]
      Sigma_list[[k2]] <- Lambda_k %*% t(Lambda_k) + Psi_k
    }
    pi_vec <- fit_cpp$pi
    logLik_val <- 0
    resp <- matrix(0, nrow=N, ncol=K)

    for(i in seq_len(N)){
      Xi <- list_of_data[[i]]
      logvals <- numeric(K)
      for(k2 in seq_len(K)){
        mu_k  <- fit_cpp$mu[[k2]]
        sig_k <- Sigma_list[[k2]]
        dens_t <- mvtnorm::dmvnorm(Xi, mean=mu_k, sigma=sig_k, log=TRUE)
        sumLog <- sum(dens_t)
        logvals[k2] <- log(pi_vec[k2] + 1e-16) + sumLog
      }
      # log-sum-exp
      m0 <- max(logvals)
      li <- m0 + log(sum(exp(logvals - m0)))
      logLik_val <- logLik_val + li

      # responsibilities
      for(k2 in seq_len(K)){
        resp[i,k2] <- exp(logvals[k2] - li)
      }
    }
    out <- list(
      z      = fit_cpp$z,
      pi     = fit_cpp$pi,
      mu     = fit_cpp$mu,
      Lambda = fit_cpp$Lambda,
      Psi    = fit_cpp$Psi,
      logLik = logLik_val,
      resp   = resp
    )
    return(out)
  }

  # Otherwise => multi-init approach

  # We will store results in a list and pick the best by logLik
  best_fit    <- NULL
  best_logLik <- -Inf

  # 1) Possibly do a k-means-based init
  if(use_kmeans_init){
    # If user wants k-means, we do it once
    features <- extract_subject_features_by_singlePCA(
      list_of_data, r_dim=subject_rdim_for_kmeans
    )
    z_init_km <- assign_by_kmeans(features, K=K)

    fit_km <- mfa_em_fit_cpp_singleInit(
      list_of_data = list_of_data,
      K = K,
      r = r,
      z_init = z_init_km,
      max_iter = max_iter,
      nIterFA  = nIterFA,
      tol      = tol
    )
    if(fit_km$logLik > best_logLik){
      best_fit    <- fit_km
      best_logLik <- fit_km$logLik
    }
  }

  # 2) Multiple random initial assignments => parallel
  if(!requireNamespace("parallel", quietly=TRUE)){
    stop("Package 'parallel' is required for mclapply. Please install it.")
  }

  N <- length(list_of_data)
  random_inits <- lapply(seq_len(n_init), function(x){
    sample.int(K, size=N, replace=TRUE)
  })

  fit_list <- parallel::mclapply(
    random_inits,
    FUN = function(z_init_rand){
      # single-run with that random init
      fit_rand <- mfa_em_fit_cpp_singleInit(
        list_of_data = list_of_data,
        K = K,
        r = r,
        z_init = z_init_rand,
        max_iter = max_iter,
        nIterFA  = nIterFA,
        tol      = tol
      )
      fit_rand
    },
    mc.cores = mc_cores
  )

  # Check all random results
  for(fit_rand in fit_list){
    if(fit_rand$logLik > best_logLik){
      best_fit    <- fit_rand
      best_logLik <- fit_rand$logLik
    }
  }

  # Return the best
  best_fit
}
#' Extract subject-level features by single PCA
#'
#' This function applies PCA to the concatenated data across all subjects,
#' then for each subject \code{i}, computes the mean of the top \code{r_dim} principal
#' component scores as a single feature vector.
#'
#' @param list_of_data A list of \code{(T_i x M)} matrices.
#' @param r_dim Number of principal components to extract.
#'
#' @return A numeric matrix of shape \code{(N x r_dim)}, where \code{N=length(list_of_data)}.
#' @export
extract_subject_features_by_singlePCA <- function(list_of_data, r_dim=2)
{
  # bigX: all rows stacked
  bigX <- do.call(rbind, list_of_data)

  # PCA
  pca_res <- prcomp(bigX, center=TRUE, scale.=FALSE)
  scores_all <- pca_res$x[, seq_len(r_dim), drop=FALSE]

  N <- length(list_of_data)
  out_features <- matrix(NA, nrow=N, ncol=r_dim)

  start_idx <- 1
  for(i in seq_len(N)){
    Ti <- nrow(list_of_data[[i]])
    end_idx <- start_idx + Ti - 1
    out_features[i,] <- colMeans(scores_all[start_idx:end_idx, , drop=FALSE])
    start_idx <- end_idx + 1
  }
  out_features
}

#' K-means-based cluster assignment
#'
#' @param features An \code{(N x d)} matrix of subject-level features.
#' @param K Number of clusters.
#'
#' @return An integer vector of length N containing cluster labels (1..K).
#' @export
assign_by_kmeans <- function(features, K){
  km_res <- kmeans(features, centers=K, nstart=10)
  km_res$cluster
}
