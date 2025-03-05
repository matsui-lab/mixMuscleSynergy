# ============================================================
# 1) Data Simulation
# ============================================================
# ============================================================
# 1) Data Simulation (Revised)
# ============================================================

#' Simulate Data from Multiple Factor-Model Clusters (Revised)
#'
#' This function generates data from \code{K_true} subgroups (clusters), each with a
#' distinct factor model (i.e., factor loadings \eqn{\Lambda_k}, mean \eqn{\mu_k},
#' and diagonal noise covariance \eqn{\Psi_k}). Subjects are assigned to clusters
#' with uniform probability, then a time-series of length \code{T_each} is sampled
#' from that subgroup's covariance.
#'
#' **本関数は、二つ目のシミュレーションモデル \code{\link{simulate_mixture_data}} のロジックに合わせて、
#' クラスタ間の平均や因子負荷量を \code{cluster_sep}, \code{noise_scale} で制御するように修正したバージョンです。**
#'
#' @param N Number of subjects (time series).
#' @param K_true True number of clusters.
#' @param r_true True factor dimension in each cluster.
#' @param M Number of observed channels (e.g., muscles).
#' @param T_each Time-series length per subject.
#' @param seed Random seed for reproducibility.
#' @param cluster_sep Numeric scalar that controls how far to separate cluster means.
#' @param noise_scale Numeric scalar that scales the noise diagonal in \eqn{\Psi_k}.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{list_of_data}}{A list of length \code{N}; each element is a \code{(T_each x M)} matrix.}
#'   \item{\code{z_true}}{An integer vector of length \code{N} for true cluster assignments (1..K_true).}
#'   \item{\code{Lambda_true, mu_true, Psi_true}}{The underlying parameters for each cluster.}
#' }
#'
#' @examples
#' \dontrun{
#' sim <- simulate_mixture_data_for_comparison(
#'   N=50, K_true=3, r_true=2, M=6, T_each=100, seed=123,
#'   cluster_sep=1.0, noise_scale=1.0
#' )
#' length(sim$list_of_data)  # 50
#' sim$z_true[1:10]
#' }
#'
#' @export
simulate_mixture_data_for_comparison <- function(
    N = 50,
    K_true = 2,
    r_true = 2,
    M = 6,
    T_each = 100,
    seed = 123,
    cluster_sep = 1.0,
    noise_scale = 1.0
){
  set.seed(seed)

  #--- クラスター割り当ての真の事前分布 ---
  pi_true <- rep(1 / K_true, K_true)
  z_true  <- sample.int(K_true, size = N, replace = TRUE, prob = pi_true)

  #--- 真のパラメータを格納するリストを準備 ---
  Lambda_list <- vector("list", K_true)
  mu_list     <- vector("list", K_true)
  Psi_list    <- vector("list", K_true)

  #--- クラスターごとの平均・因子負荷量・ノイズを設定 ---
  for(k in seq_len(K_true)) {
    # 平均ベクトル: base_mu + offset
    base_mu <- rnorm(M, mean = 0, sd = 0.3)
    offset  <- (k - (K_true + 1) / 2) * cluster_sep
    mu_k    <- base_mu + offset

    # 因子負荷量のスケール: 0.2 + 0.2 * cluster_sep
    lambda_scale <- 0.2 + 0.2 * cluster_sep
    Lambda_k     <- matrix(rnorm(M * r_true, sd = lambda_scale),
                           nrow = M, ncol = r_true)

    # ノイズの対角成分: runif(...) * noise_scale
    psi_vec <- runif(M, min = 0.05, max = 0.2) * noise_scale
    Psi_k   <- diag(psi_vec, nrow = M)

    mu_list[[k]]     <- mu_k
    Lambda_list[[k]] <- Lambda_k
    Psi_list[[k]]    <- Psi_k
  }

  #--- データ生成 ---
  list_of_data <- vector("list", N)
  for(i in seq_len(N)) {
    k_i   <- z_true[i]
    Lk    <- Lambda_list[[k_i]]
    mu_k  <- mu_list[[k_i]]
    Psi_k <- Psi_list[[k_i]]

    Sigma_k <- Lk %*% t(Lk) + Psi_k
    Xi      <- MASS::mvrnorm(n = T_each, mu = mu_k, Sigma = Sigma_k)
    list_of_data[[i]] <- Xi
  }

  #--- 結果をまとめて返す ---
  list(
    list_of_data = list_of_data,
    z_true       = z_true,
    Lambda_true  = Lambda_list,
    mu_true      = mu_list,
    Psi_true     = Psi_list
  )
}


# ============================================================
# 2) Single Factor Analysis (FA)
# ============================================================

#' Fit a Single Factor Analysis (Ignoring Clusters)
#'
#' This function concatenates all subjects' data into one large matrix,
#' then performs a naive factor analysis. It returns a simple BIC measure
#' and the fitted loadings, mean, and diagonal Psi.
#'
#' @param list_of_data A list of (T_i x M) matrices.
#' @param r The factor dimension.
#'
#' @return A list with elements \code{(loadings, mu, diagPsi, logLik, BIC)}.
#' @export
fit_single_factor_analysis <- function(list_of_data, r=2) {
  bigX <- do.call(rbind, list_of_data)
  n_total <- nrow(bigX)
  M <- ncol(bigX)

  mu <- colMeans(bigX)
  Xc <- sweep(bigX, 2, mu, "-")
  Sxx <- crossprod(Xc) / n_total

  eig_out <- eigen(Sxx, symmetric = TRUE)
  d_r  <- eig_out$values[1:r]
  V_r  <- eig_out$vectors[, 1:r, drop = FALSE]
  loadings <- V_r %*% diag(sqrt(d_r), r, r)

  recon <- loadings %*% t(loadings)
  diagPsi <- diag(Sxx) - diag(recon)
  diagPsi[diagPsi < 1e-12] <- 1e-12

  Sigma <- recon
  diag(Sigma) <- diag(Sigma) + diagPsi
  cholS <- chol(Sigma)
  logdetS <- 2 * sum(log(diag(cholS)))
  Mlog2pi <- M * log(2 * pi)
  invS <- solve(Sigma)
  quadSum <- 0
  for(i in seq_len(n_total)) {
    xi <- bigX[i, ]
    diff <- xi - mu
    quadSum <- quadSum + (t(diff) %*% invS %*% diff)
  }
  logLik_val <- -0.5 * ( n_total * (Mlog2pi + logdetS) + quadSum )

  # naive param count
  param_count <- M * r + M
  BIC_val <- -2 * logLik_val + param_count * log(n_total)

  list(
    loadings = loadings,
    mu = mu,
    diagPsi = diagPsi,
    logLik = logLik_val,
    BIC    = BIC_val
  )
}

# Helpers for reconstruction error
compute_factor_scores_singleFA <- function(X, mu, Lambda, diagPsi) {
  Xc <- sweep(X, 2, mu, "-")
  M <- ncol(X)
  r <- ncol(Lambda)

  Sigma <- Lambda %*% t(Lambda)
  diag(Sigma) <- diag(Sigma) + diagPsi
  invS <- solve(Sigma)

  W <- t(Lambda) %*% invS  # (r x M)
  A <- diag(r) + W %*% Lambda
  Ainv <- solve(A)
  EZt <- Ainv %*% W %*% t(Xc)
  EZ  <- t(EZt)
  EZ
}

reconstruct_singleFA <- function(X, mu, Lambda, diagPsi) {
  scores <- compute_factor_scores_singleFA(X, mu, Lambda, diagPsi)
  Xhat   <- scores %*% t(Lambda)
  Xhat   <- sweep(Xhat, 2, mu, "+")
  Xhat
}

#' Compute Total Reconstruction SSE for Single FA
#'
#' @param list_of_data The original data list.
#' @param fa_res A result from \code{fit_single_factor_analysis}.
#'
#' @return A numeric value of the sum of squared errors (SSE).
#' @export
calc_reconstruction_error_singleFA <- function(list_of_data, fa_res) {
  loadings <- fa_res$loadings
  mu <- fa_res$mu
  diagPsi <- fa_res$diagPsi
  SSE <- 0
  for(X_i in list_of_data) {
    Xhat_i <- reconstruct_singleFA(X_i, mu, loadings, diagPsi)
    diff <- X_i - Xhat_i
    SSE <- SSE + sum(diff^2)
  }
  SSE
}

# ============================================================
# 3) Single PCA (Ignoring Clusters)
# ============================================================

#' Fit a Single PCA (Ignoring Clusters)
#'
#' Concatenates all data into one matrix, performs SVD, keeps \code{r}
#' principal components, then computes a naive PPCA logLik and BIC.
#'
#' @param list_of_data List of (T_i x M) matrices.
#' @param r Number of principal components.
#'
#' @return A list with \code{(P, mu, sigma2, logLik, BIC)}.
#' @export
fit_single_pca <- function(list_of_data, r=2) {
  bigX <- do.call(rbind, list_of_data)
  n_total <- nrow(bigX)
  M <- ncol(bigX)

  mu <- colMeans(bigX)
  Xc <- sweep(bigX, 2, mu, "-")

  svd_out <- svd(Xc)
  D_r <- svd_out$d[1:r]
  V_r <- svd_out$v[, 1:r, drop = FALSE]
  P <- V_r  # principal directions (M x r)

  lam_j <- (D_r^2) / n_total
  if(r < M) {
    leftover <- 0
    for(j in seq.int(r+1, M)) {
      leftover <- leftover + (svd_out$d[j]^2) / n_total
    }
    sigma2 <- leftover / (M - r)
  } else {
    sigma2 <- 1e-12
  }

  # rank-r => Sigma = W W^T + sigma^2 I
  W <- matrix(0, nrow = M, ncol = r)
  for(j in seq_len(r)) {
    val_j <- lam_j[j] - sigma2
    if(val_j < 0) val_j <- 0
    W[, j] <- P[, j] * sqrt(val_j)
  }
  Sigma <- W %*% t(W)
  diag(Sigma) <- diag(Sigma) + sigma2

  cholS <- chol(Sigma)
  logdetS <- 2 * sum(log(diag(cholS)))
  Mlog2pi <- M * log(2 * pi)
  invS <- solve(Sigma)
  quadSum <- 0
  for(i in seq_len(n_total)) {
    xi <- bigX[i, ]
    diff <- xi - mu
    quadSum <- quadSum + (t(diff) %*% invS %*% diff)
  }
  logLik_val <- -0.5 * (n_total * (Mlog2pi + logdetS) + quadSum)

  # param count ~ (M*r) + M + 1
  param_count <- (M * r) + M + 1
  BIC_val <- -2 * logLik_val + param_count * log(n_total)

  list(
    P = P,
    mu = mu,
    sigma2 = sigma2,
    logLik = logLik_val,
    BIC    = BIC_val
  )
}

reconstruct_singlePCA <- function(X, pca_res) {
  P  <- pca_res$P
  mu <- pca_res$mu
  Xc <- sweep(X, 2, mu, "-")
  Z <- Xc %*% P
  Xhat <- Z %*% t(P)
  Xhat <- sweep(Xhat, 2, mu, "+")
  Xhat
}

#' Compute Total Reconstruction SSE for Single PCA
#'
#' @param list_of_data The data list.
#' @param pca_res The result from \code{fit_single_pca}.
#'
#' @return SSE (sum of squared errors).
#' @export
calc_reconstruction_error_singlePCA <- function(list_of_data, pca_res) {
  SSE <- 0
  for(X_i in list_of_data) {
    Xhat_i <- reconstruct_singlePCA(X_i, pca_res)
    diff <- X_i - Xhat_i
    SSE <- SSE + sum(diff^2)
  }
  SSE
}

# ============================================================
# 4) Mixture FA / Mixture PCA
#    (We assume you have mfa_em_fit / mixture_pca_em_fit / etc.)
# ============================================================

# Helpers for reconstruction under MixtureFA

reconstruct_mixtureFA_oneCluster <- function(X, k, mfa_fit) {
  mu_k      <- mfa_fit$mu[[k]]
  Lambda_k  <- mfa_fit$Lambda[[k]]
  diagPsi_k <- diag(mfa_fit$Psi[[k]])
  Xhat <- reconstruct_singleFA(X, mu_k, Lambda_k, diagPsi_k)
  Xhat
}

#' Compute SSE for a MixtureFA model
#'
#' Uses the final hard assignments \code{mfa_fit$z}, reconstruct each subject's data.
#'
#' @param list_of_data The data list.
#' @param mfa_fit The mixture FA result, which contains \code{$z}, \code{$Lambda[[k]]}, \code{$mu[[k]]}, \code{$Psi[[k]]}.
#'
#' @return SSE (sum of squared errors).
#' @export
calc_reconstruction_error_mixtureFA <- function(list_of_data, mfa_fit) {
  z_vec <- mfa_fit$z
  SSE <- 0
  for(i in seq_along(list_of_data)) {
    k_i <- z_vec[i]
    X_i <- list_of_data[[i]]
    Xhat_i <- reconstruct_mixtureFA_oneCluster(X_i, k_i, mfa_fit)
    SSE <- SSE + sum((X_i - Xhat_i)^2)
  }
  SSE
}

# Helpers for reconstruction under MixturePCA

reconstruct_mixturePCA_oneCluster <- function(X, k, mpca_fit) {
  P_k  <- mpca_fit$P[[k]]
  mu_k <- mpca_fit$mu[[k]]
  # simple PCA reconstruction
  Xc <- sweep(X, 2, mu_k, "-")
  Z  <- Xc %*% P_k
  Xhat <- Z %*% t(P_k)
  Xhat <- sweep(Xhat, 2, mu_k, "+")
  Xhat
}

#' Compute SSE for a MixturePCA model
#'
#' @param list_of_data The data list.
#' @param mpca_fit The mixture PCA result, containing \code{$z} and each cluster's \code{P, mu}.
#'
#' @return SSE (sum of squared errors).
#' @export
calc_reconstruction_error_mixturePCA <- function(list_of_data, mpca_fit) {
  z_vec <- mpca_fit$z
  SSE <- 0
  for(i in seq_along(list_of_data)) {
    k_i <- z_vec[i]
    X_i <- list_of_data[[i]]
    Xhat_i <- reconstruct_mixturePCA_oneCluster(X_i, k_i, mpca_fit)
    SSE <- SSE + sum((X_i - Xhat_i)^2)
  }
  SSE
}
#' Compute the Total Sum of Squares (SST) for a List of Data Matrices
#'
#' This function concatenates all matrices in \code{list_of_data} into one big matrix,
#' calculates the global mean, and then sums the squared differences from that mean.
#' The result can be used to compute VAF (Variance Accounted For).
#'
#' @param list_of_data A list of \code{(T_i x M)} matrices.
#' @return A numeric scalar, the total sum of squares (SST).
#' @export
calc_total_SST <- function(list_of_data) {
  # 1) Concatenate everything
  bigX <- do.call(rbind, list_of_data)  # shape: (sum T_i) x M

  # 2) Global mean across all rows
  mu_global <- colMeans(bigX)

  # 3) Sum of squared diffs
  Xc <- sweep(bigX, 2, mu_global, "-")  # center
  sst <- sum(Xc^2)  # sum of squares
  sst
}
# ============================================================
# 5) Main Comparison: SingleFA, SinglePCA, MixtureFA, MixturePCA
#    (two-step approach is excluded)
# ============================================================

#' Simulate and Compare 4 Methods (SingleFA, SinglePCA, MixtureFA, MixturePCA) + VAF
#'
#' This function simulates data from a factor-based mixture (using
#' \code{\link{simulate_mixture_data_for_comparison}}). Then it fits:
#' \enumerate{
#'   \item Single FA (ignoring clusters)
#'   \item Single PCA (ignoring clusters)
#'   \item Mixture FA
#'   \item Mixture PCA
#' }
#' It computes BIC, ARI (for the mixture methods), reconstruction error (SSE),
#' **and** VAF (Variance Accounted For).
#'
#' @param N,K_true,r_true See \code{\link{simulate_mixture_data_for_comparison}}.
#' @param M Number of channels.
#' @param T_each Time-series length per subject.
#' @param seed Random seed.
#' @param r_singleFA Factor dimension for the single-factor-analysis approach.
#' @param r_singlePCA Number of principal components for single-PCA approach.
#' @param cluster_sep Numeric scalar that controls how far to separate cluster means
#'   (passed to \code{\link{simulate_mixture_data_for_comparison}}).
#' @param noise_scale Numeric scalar that scales the noise diagonal in \eqn{\Psi_k}
#'   (passed to \code{\link{simulate_mixture_data_for_comparison}}).
#'
#' @param n_init Number of random initial assignments for mixture FA/PCA.
#'   Defaults to \code{1} (the original single-run behavior).
#' @param use_kmeans_init Logical; if \code{TRUE}, also try a k-means-based initialization.
#'   Defaults to \code{FALSE}.
#' @param subject_rdim_for_kmeans PCA dimension for subject-level features in k-means init,
#'   defaults to \code{r_true} (the same factor dimension used in mixture).
#' @param mc_cores Number of cores for parallel runs (passed to \code{mclapply}). Default is 1.
#' @param method_pca Either \code{"EM"} or \code{"closed_form"} for the mixture PCA approach.
#'   Defaults to \code{"EM"}.
#'
#' @return A data frame with columns: \code{Method, BIC, ARI, SSE, VAF}.
#'
#' @details
#' \itemize{
#'   \item The single approaches do not produce a cluster assignment, so \code{ARI=NA}.
#'   \item Mixture FA/PCA produce a hard assignment \code{z} that can be compared to the
#'         true \code{z_true}, yielding a meaningful ARI.
#'   \item SSE is computed by reconstructing each subject's data with the fitted model.
#'   \item VAF is defined as \code{1 - SSE / totalSST}, where totalSST is computed
#'         by \code{\link{calc_total_SST}} across all data.
#' }
#'
#' If you set \code{n_init > 1} or \code{use_kmeans_init=TRUE}, the mixture FA/PCA
#' calls will attempt multiple initial assignments (including possibly a k-means-based one)
#' and pick the best final solution by log-likelihood. If you keep \code{n_init=1}
#' and \code{use_kmeans_init=FALSE}, it falls back to the original single-run approach.
#'
#' @examples
#' \dontrun{
#' df_comp <- simulate_and_compare_4methods_noTwoStep(
#'   N=100, K_true=3, r_true=4, M=6, T_each=100, seed=123,
#'   r_singleFA=4, r_singlePCA=4,
#'   cluster_sep=1.0, noise_scale=1.0,
#'   n_init=5, use_kmeans_init=TRUE, mc_cores=2
#' )
#' print(df_comp)
#' }
#'
#' @export
simulate_and_compare_4methods_noTwoStep <- function(
    N=50, K_true=2, r_true=2, M=6, T_each=100, seed=123,
    r_singleFA=2,
    r_singlePCA=2,
    cluster_sep=1.0,
    noise_scale=1.0,
    # New arguments for multi-init usage
    n_init = 1,
    use_kmeans_init = FALSE,
    subject_rdim_for_kmeans = r_true,
    mc_cores = 1,
    method_pca = "EM"
)
{
  # 1) Generate data
  sim_data <- simulate_mixture_data_for_comparison(
    N         = N,
    K_true    = K_true,
    r_true    = r_true,
    M         = M,
    T_each    = T_each,
    seed      = seed,
    cluster_sep  = cluster_sep,
    noise_scale = noise_scale
  )
  list_of_data <- sim_data$list_of_data
  z_true       <- sim_data$z_true

  # 1.5) Compute total sum of squares (SST) for VAF
  total_SST <- calc_total_SST(list_of_data)  # from the helper function

  # 2) Single FA
  singleFA_res <- fit_single_factor_analysis(list_of_data, r = r_singleFA)
  BIC_singleFA <- singleFA_res$BIC
  SSE_singleFA <- calc_reconstruction_error_singleFA(list_of_data, singleFA_res)
  VAF_singleFA <- 1 - (SSE_singleFA / total_SST)

  # 3) Single PCA
  singlePCA_res <- fit_single_pca(list_of_data, r = r_singlePCA)
  BIC_singlePCA <- singlePCA_res$BIC
  SSE_singlePCA <- calc_reconstruction_error_singlePCA(list_of_data, singlePCA_res)
  VAF_singlePCA <- 1 - (SSE_singlePCA / total_SST)

  # 4) Mixture FA
  mixtureFA_fit <- mfa_em_fit(
    list_of_data,
    K = K_true,
    r = r_true,
    max_iter = 50,
    nIterFA  = 20,
    tol      = 1e-4,
    n_init   = n_init,
    use_kmeans_init = use_kmeans_init,
    subject_rdim_for_kmeans = subject_rdim_for_kmeans,
    mc_cores = mc_cores
  )
  ll_mixFA  <- compute_logLik_mfa(list_of_data, mixtureFA_fit)
  BIC_mixFA <- compute_BIC_mfa(ll_mixFA, K = K_true, r = r_true, M = M, N = N)
  SSE_mixFA <- calc_reconstruction_error_mixtureFA(list_of_data, mixtureFA_fit)
  VAF_mixFA <- 1 - (SSE_mixFA / total_SST)

  if(!requireNamespace("mclust", quietly = TRUE)) {
    stop("Package 'mclust' is required for ARI computation. Please install it.")
  }
  ari_mixFA <- mclust::adjustedRandIndex(mixtureFA_fit$z, z_true)

  # 5) Mixture PCA
  mixturePCA_fit <- mixture_pca_em_fit(
    list_of_data,
    K = K_true,
    r = r_true,
    max_iter  = 50,
    nIterPCA  = 20,
    tol       = 1e-4,
    method    = method_pca,
    n_init    = n_init,
    use_kmeans_init = use_kmeans_init,
    subject_rdim_for_kmeans = subject_rdim_for_kmeans,
    mc_cores = mc_cores
  )
  ll_mixPCA  <- compute_logLik_mpca(list_of_data, mixturePCA_fit)
  BIC_mixPCA <- compute_BIC_mpca(ll_mixPCA, K = K_true, r = r_true, M = M, N = N)
  SSE_mixPCA <- calc_reconstruction_error_mixturePCA(list_of_data, mixturePCA_fit)
  VAF_mixPCA <- 1 - (SSE_mixPCA / total_SST)
  ari_mixPCA <- mclust::adjustedRandIndex(mixturePCA_fit$z, z_true)

  # Compile results, now with VAF
  df_out <- data.frame(
    Method = c("SingleFA",    "SinglePCA",    "MixtureFA",    "MixturePCA"),
    BIC    = c(BIC_singleFA,  BIC_singlePCA,  BIC_mixFA,      BIC_mixPCA),
    ARI    = c(NA,            NA,             ari_mixFA,      ari_mixPCA),
    SSE    = c(SSE_singleFA,  SSE_singlePCA,  SSE_mixFA,      SSE_mixPCA),
    VAF    = c(VAF_singleFA,  VAF_singlePCA,  VAF_mixFA,      VAF_mixPCA)
  )
  df_out
}
