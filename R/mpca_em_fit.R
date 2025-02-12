#' Single-Run Mixture PCA with a Preset Initial Assignment
#'
#' This helper function calls the updated C++ routine \code{mpcaTimeseriesCpp()}
#' with a user-specified cluster assignment (\code{z_init}), then converts
#' the returned \code{W} into \code{(P, D)} and constructs \code{Psi}, and finally
#' computes the final log-likelihood and responsibilities in R.
#'
#' @param list_of_data A list of matrices (each \code{T_i x M}).
#' @param K Number of clusters.
#' @param r Number of principal components (dimension).
#' @param max_iter Maximum EM iterations (C++).
#' @param nIterPCA Sub-iterations for updating each cluster's PPCA parameters in C++.
#' @param tol Convergence tolerance (C++).
#' @param method Either \code{"EM"} or \code{"closed_form"}, passed to \code{mpcaTimeseriesCpp}.
#' @param z_init An integer vector of length \code{N} giving initial cluster labels (1..K).
#'
#' @return A list with elements:
#'   \item{z}{Hard cluster assignments (1..K).}
#'   \item{pi}{Mixing proportions (length K).}
#'   \item{mu}{List of length K, each a mean vector (length M).}
#'   \item{P}{List of length K, each an \code{(M x r)} matrix of principal directions.}
#'   \item{D}{List of length K, each an \code{(r x r)} diagonal matrix of column norms.}
#'   \item{Psi}{List of length K, each an \code{(M x M)} diagonal matrix \code{sigma2 * I}.}
#'   \item{logLik}{Final log-likelihood computed in R.}
#'   \item{sigma2}{Numeric vector of length K, the final \code{sigma2} for each cluster.}
#'   \item{resp}{An \code{(N x K)} matrix of responsibilities.}
#'
#' @examples
#' \dontrun{
#' # Suppose we have list_of_data, K=3, r=2, method="EM", plus some z_init:
#' z_init <- sample.int(3, size=length(list_of_data), replace=TRUE)
#' fit_one <- mixture_pca_em_fit_cpp_singleInit(
#'   list_of_data, K=3, r=2, max_iter=50, nIterPCA=20, tol=1e-3,
#'   method="EM", z_init=z_init
#' )
#' print(fit_one$logLik)
#' }
#'
#' @export
mixture_pca_em_fit_cpp_singleInit <- function(list_of_data,
                                              K,
                                              r,
                                              max_iter = 50,
                                              nIterPCA = 20,
                                              tol      = 1e-3,
                                              method   = "EM",
                                              z_init)
{
  # 1) Call the updated C++ function with z_init
  fit_cpp <- mpcaTimeseriesCpp(
    list_of_data = list_of_data,
    K = K,
    r = r,
    max_iter  = max_iter,
    nIterPCA  = nIterPCA,
    tol       = tol,
    method    = method,
    z_init    = z_init
  )

  # 2) Extract results from C++
  #    W_list_cpp, mu_list_cpp, sigma2_vec, pi_vec, z_cpp
  z_cpp      <- fit_cpp$z
  W_list_cpp <- fit_cpp$W
  mu_list_cpp<- fit_cpp$mu
  sigma2_vec <- fit_cpp$sigma2
  pi_vec     <- fit_cpp$pi

  # 3) Convert W -> (P, D) and build Psi
  K_check <- length(W_list_cpp)
  if(K_check != K){
    stop("C++ output mismatch: length(W_list_cpp) != K.")
  }

  P_list   <- vector("list", K)
  D_list   <- vector("list", K)
  Psi_list <- vector("list", K)

  for(k2 in seq_len(K)){
    W_k <- W_list_cpp[[k2]]
    M_k <- nrow(W_k)
    r_k <- ncol(W_k)

    # column norms -> diagonal of D
    col_scales <- numeric(r_k)
    for(j in seq_len(r_k)){
      cs_j <- sqrt(sum(W_k[,j]^2))
      if(cs_j < 1e-12) cs_j <- 1e-12
      col_scales[j] <- cs_j
    }
    P_k <- W_k
    for(j in seq_len(r_k)){
      P_k[, j] <- W_k[, j] / col_scales[j]
    }
    D_k <- diag(col_scales, r_k, r_k)

    sig2_k <- sigma2_vec[k2]
    Psi_k  <- diag(sig2_k, M_k)

    P_list[[k2]]  <- P_k
    D_list[[k2]]  <- D_k
    Psi_list[[k2]]<- Psi_k
  }

  # 4) Compute final log-likelihood and responsibilities in R
  N <- length(list_of_data)
  Sigma_list <- vector("list", K)
  for(k2 in seq_len(K)){
    # Sigma_k = W_k W_k^T + sigma2 * I
    # or equivalently P_k D_k^2 P_k^T + sigma2 * I
    P_k  <- P_list[[k2]]
    D_k  <- D_list[[k2]]
    sig2 <- sigma2_vec[k2]
    Sig_k <- P_k %*% (D_k^2) %*% t(P_k)
    diag(Sig_k) <- diag(Sig_k) + sig2
    Sigma_list[[k2]] <- Sig_k
  }

  logLik_val <- 0
  resp <- matrix(0, nrow=N, ncol=K)
  for(i in seq_len(N)){
    Xi <- list_of_data[[i]]
    logvals <- numeric(K)
    for(k2 in seq_len(K)){
      mu_k  <- mu_list_cpp[[k2]]
      Sig_k <- Sigma_list[[k2]]
      dens_t <- mvtnorm::dmvnorm(Xi, mean=mu_k, sigma=Sig_k, log=TRUE)
      sumLog <- sum(dens_t)
      logvals[k2] <- log(pi_vec[k2] + 1e-16) + sumLog
    }
    # log-sum-exp
    m0 <- max(logvals)
    li <- m0 + log(sum(exp(logvals - m0)))
    logLik_val <- logLik_val + li

    for(k2 in seq_len(K)){
      resp[i,k2] <- exp(logvals[k2] - li)
    }
  }

  # 5) Return final structure
  out <- list(
    z      = z_cpp,
    pi     = pi_vec,
    mu     = mu_list_cpp,
    P      = P_list,
    D      = D_list,
    Psi    = Psi_list,
    logLik = logLik_val,
    sigma2 = sigma2_vec,
    resp   = resp
  )
  out
}

#' Fit a Mixture PCA Model (Using mpcaTimeseriesCpp) with optional multi-initialization
#'
#' This function calls the C++ function \code{mpcaTimeseriesCpp()} to perform a
#' Mixture PCA (PPCA) EM algorithm for fixed \code{K} and \code{r}, then converts
#' the \code{W} matrices into \code{(P, D)}, creates \code{Psi}, and computes the final
#' log-likelihood/responsibilities in R. By default, it runs a single pass with
#' the built-in initialization in C++. However, if \code{n_init > 1} or
#' \code{use_kmeans_init=TRUE}, it performs multiple initializations in parallel
#' and returns the best-fitting result (highest log-likelihood).
#'
#' @param list_of_data A list of matrices (each \code{T_i x M}).
#' @param K Number of clusters.
#' @param r Number of principal components.
#' @param max_iter Maximum EM iterations for the C++ routine.
#' @param nIterPCA Sub-iterations for updating each cluster's PPCA parameters.
#' @param tol Convergence tolerance for \code{mpcaTimeseriesCpp}.
#' @param method Either \code{"EM"} or \code{"closed_form"}, passed down to the C++ routine.
#' @param n_init Integer; how many random initial assignments to try (in addition to
#'   the default single-run or k-means if requested). Defaults to \code{1}.
#' @param use_kmeans_init Logical; if \code{TRUE}, we also run one initialization
#'   where we assign clusters by k-means on subject-level PCA features. Defaults to \code{FALSE}.
#' @param subject_rdim_for_kmeans The PCA dimension for the subject-level feature extraction,
#'   used only if \code{use_kmeans_init=TRUE}. Defaults to \code{r}.
#' @param mc_cores Number of cores for parallel execution via \code{mclapply}.
#'   Defaults to \code{1} (no parallel).
#'
#' @return A list with elements:
#'   \item{z}{Hard cluster assignments (1..K).}
#'   \item{pi}{Mixing proportions.}
#'   \item{mu}{List of length K, each a mean vector.}
#'   \item{P}{List of length K, each \code{(M x r)} principal directions.}
#'   \item{D}{List of length K, each \code{(r x r)} diagonal.}
#'   \item{Psi}{List of length K, each \code{(M x M)} diagonal.}
#'   \item{logLik}{Final log-likelihood over all data.}
#'   \item{sigma2}{Numeric vector of length K for \code{sigma2}.}
#'   \item{resp}{\code{(N x K)} matrix of responsibilities.}
#'
#' @details
#' If \code{n_init=1} and \code{use_kmeans_init=FALSE}, this function runs exactly one pass
#' with the default initialization in C++ (i.e. subject \code{i} is assigned to cluster
#' \code{(i \% K) + 1}). Otherwise:
#' \enumerate{
#'   \item If \code{use_kmeans_init=TRUE}, we do a run where we apply k-means to some
#'         subject-level PCA features to get \code{z_init}, then run EM once.
#'   \item We generate \code{n_init} random initial assignments, run EM for each in parallel,
#'         and collect the results.
#'   \item We compare all solutions by final log-likelihood and pick the best.
#' }
#'
#' @examples
#' \dontrun{
#' # Suppose we have a list_of_data, K=3, r=2, and we want to try 5 random inits + k-means:
#' fit <- mixture_pca_em_fit(
#'   list_of_data, K=3, r=2, max_iter=50, nIterPCA=20, tol=1e-3, method="EM",
#'   n_init=5, use_kmeans_init=TRUE, subject_rdim_for_kmeans=2, mc_cores=2
#' )
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
                               method   = "EM",
                               n_init   = 1,
                               use_kmeans_init = FALSE,
                               subject_rdim_for_kmeans = r,
                               mc_cores = 1)
{
  # If n_init=1 and use_kmeans_init=FALSE => original single-run approach
  if(n_init == 1 && !use_kmeans_init){
    # Just call mpcaTimeseriesCpp once with no z_init => same as original
    fit_cpp <- mpcaTimeseriesCpp(
      list_of_data = list_of_data,
      K = K,
      r = r,
      max_iter  = max_iter,
      nIterPCA  = nIterPCA,
      tol       = tol,
      method    = method
    )

    # Then do the original interpretation (W->(P,D), etc.) + logLik calculation
    # (Identical to your original code in mixture_pca_em_fit)
    N <- length(list_of_data)
    if(N < 1){
      stop("No data in list_of_data.")
    }
    W_list_cpp <- fit_cpp$W
    mu_list_cpp<- fit_cpp$mu
    sigma2_vec <- fit_cpp$sigma2
    pi_vec     <- fit_cpp$pi
    z_cpp      <- fit_cpp$z

    if(is.null(sigma2_vec)){
      stop("mpcaTimeseriesCpp did not return 'sigma2'.")
    }

    # Convert W->(P,D), build Psi
    K_check <- length(W_list_cpp)
    if(K_check != K){
      stop("Mismatch in K: length(W_list_cpp) != K.")
    }
    P_list  <- vector("list", K)
    D_list  <- vector("list", K)
    Psi_list<- vector("list", K)

    for(k2 in seq_len(K)){
      W_k <- W_list_cpp[[k2]]
      M_k <- nrow(W_k)
      r_k <- ncol(W_k)
      col_scales <- numeric(r_k)
      for(j in seq_len(r_k)){
        cs_j <- sqrt(sum(W_k[,j]^2))
        if(cs_j < 1e-12) cs_j <- 1e-12
        col_scales[j] <- cs_j
      }
      P_k <- W_k
      for(j in seq_len(r_k)){
        P_k[, j] <- W_k[, j] / col_scales[j]
      }
      D_k <- diag(col_scales, r_k, r_k)
      sig2_k <- sigma2_vec[k2]
      Psi_k  <- diag(sig2_k, M_k)

      P_list[[k2]]  <- P_k
      D_list[[k2]]  <- D_k
      Psi_list[[k2]]<- Psi_k
    }

    # Compute final logLik & resp
    Sigma_list <- vector("list", K)
    for(k2 in seq_len(K)){
      P_k   <- P_list[[k2]]
      D_k   <- D_list[[k2]]
      sig2_k<- sigma2_vec[k2]
      Sig_k <- P_k %*% (D_k^2) %*% t(P_k)
      diag(Sig_k) <- diag(Sig_k) + sig2_k
      Sigma_list[[k2]] <- Sig_k
    }
    logLik_val <- 0
    resp <- matrix(0, nrow=N, ncol=K)
    for(i in seq_len(N)){
      Xi <- list_of_data[[i]]
      logvals <- numeric(K)
      for(k2 in seq_len(K)){
        mu_k  <- mu_list_cpp[[k2]]
        Sig_k <- Sigma_list[[k2]]
        dens_t <- mvtnorm::dmvnorm(Xi, mean=mu_k, sigma=Sig_k, log=TRUE)
        sumLog <- sum(dens_t)
        logvals[k2] <- log(pi_vec[k2] + 1e-16) + sumLog
      }
      m0 <- max(logvals)
      li <- m0 + log(sum(exp(logvals - m0)))
      logLik_val <- logLik_val + li
      for(k2 in seq_len(K)){
        resp[i,k2] <- exp(logvals[k2] - li)
      }
    }

    out <- list(
      z      = z_cpp,
      pi     = pi_vec,
      mu     = mu_list_cpp,
      P      = P_list,
      D      = D_list,
      Psi    = Psi_list,
      logLik = logLik_val,
      sigma2 = sigma2_vec,
      resp   = resp
    )
    return(out)
  }

  # Otherwise => multi-init approach
  best_fit    <- NULL
  best_logLik <- -Inf

  # 1) If k-means init is requested => run once
  if(use_kmeans_init){
    if(!requireNamespace("parallel", quietly=TRUE)){
      stop("Package 'parallel' is required for mclapply. Please install it.")
    }
    features <- extract_subject_features_by_singlePCA(
      list_of_data, r_dim=subject_rdim_for_kmeans
    )
    z_init_km <- assign_by_kmeans(features, K=K)

    fit_km <- mixture_pca_em_fit_cpp_singleInit(
      list_of_data = list_of_data,
      K = K,
      r = r,
      max_iter  = max_iter,
      nIterPCA  = nIterPCA,
      tol       = tol,
      method    = method,
      z_init    = z_init_km
    )
    if(fit_km$logLik > best_logLik){
      best_fit    <- fit_km
      best_logLik <- fit_km$logLik
    }
  }

  # 2) Random initializations => parallel
  N <- length(list_of_data)
  random_inits <- lapply(seq_len(n_init), function(i) {
    sample.int(K, size=N, replace=TRUE)
  })

  if(!requireNamespace("parallel", quietly=TRUE)){
    stop("Package 'parallel' is required for mclapply. Please install it.")
  }
  fit_list <- parallel::mclapply(
    random_inits,
    FUN = function(z_init_rand){
      fit_rand <- mixture_pca_em_fit_cpp_singleInit(
        list_of_data = list_of_data,
        K = K,
        r = r,
        max_iter  = max_iter,
        nIterPCA  = nIterPCA,
        tol       = tol,
        method    = method,
        z_init    = z_init_rand
      )
      fit_rand
    },
    mc.cores = mc_cores
  )

  # Compare random solutions
  for(fit_rand in fit_list){
    if(fit_rand$logLik > best_logLik){
      best_fit    <- fit_rand
      best_logLik <- fit_rand$logLik
    }
  }

  best_fit
}
#' Extract subject-level features by single PCA
#'
#' Similar to the MFA case, we apply PCA to the concatenated data, then average
#' the top \code{r_dim} PC scores for each subject.
#'
#' @param list_of_data A list of \code{(T_i x M)} matrices.
#' @param r_dim Number of principal components to use for feature extraction.
#'
#' @return A \code{(N x r_dim)} matrix, where N=length(list_of_data).
#' @export
extract_subject_features_by_singlePCA <- function(list_of_data, r_dim=2)
{
  bigX <- do.call(rbind, list_of_data)
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
#' @param features A \code{(N x d)} matrix of subject-level features.
#' @param K Number of clusters.
#' @return An integer vector of length \code{N}, each in \code{1..K}.
#' @export
assign_by_kmeans <- function(features, K){
  km <- kmeans(features, centers=K, nstart=10)
  km$cluster
}

