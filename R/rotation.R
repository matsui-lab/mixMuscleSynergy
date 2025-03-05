#' Post-hoc Rotation for an MFA Model
#'
#' Applies a chosen rotation (e.g., "varimax") to each cluster's loading matrix
#' (\code{Lambda[[k]]}). Optionally, also rotates the factor scores in
#' \code{factor_scores[[i]]} if they exist and \code{rotate_scores=TRUE}.
#'
#' @param mfa_fit A fitted MFA model, typically returned by \code{mfa_em_fit()}
#'   or \code{select_optimal_K_r_mfa()$best_model_info$model}, containing:
#'   \itemize{
#'     \item \code{Lambda[[k]]}: (M x r) factor loadings for each cluster.
#'     \item \code{factor_scores[[i]]}: optional (T_i x r) factor scores for subject i.
#'     \item \code{z[i]}: cluster assignments.
#'   }
#' @param rotation A string specifying the rotation method. Currently supports \code{"varimax"}.
#' @param rotate_scores Logical; if TRUE and \code{factor_scores} exist, rotate them as well.
#'   If \code{factor_scores} do not exist and \code{rotate_scores=TRUE}, an error or warning is produced.
#' @param cluster_ids Vector of cluster indices to rotate. Defaults to all.
#'
#' @return A new MFA model object, with rotated \code{Lambda} (and \code{factor_scores} if requested).
#'   Other fields are unchanged.
#'
#' @details
#' - Rotation does not change log-likelihood or BIC, just reorients factors for interpretability.
#' - If \code{factor_scores} is missing and \code{rotate_scores=TRUE}, you can either skip or error.
#'
#' @export
#' Post-hoc Rotation for an MFA Model
#'
#' Applies a chosen rotation (e.g. "varimax") to each cluster's loading matrix
#' (\code{Lambda[[k]]}) in a fitted MFA model, optionally also rotating the subject
#' factor scores in \code{factor_scores[[i]]}.
#'
#' @param mfa_fit A fitted MFA model, typically returned by \code{mfa_em_fit} or
#'   \code{select_optimal_K_r_mfa()$best_model_info$model}, containing:
#'   \itemize{
#'     \item \code{Lambda[[k]]}: (M x r) factor loadings
#'     \item \code{factor_scores[[i]]}: optional (T_i x r) factor scores
#'     \item \code{z[i]}: cluster assignments
#'   }
#' @param rotation A string specifying the rotation method. Currently supports \code{"varimax"}.
#' @param rotate_scores Logical; if TRUE, also rotate the factor scores if present.
#' @param cluster_ids A numeric vector specifying which cluster IDs to rotate. Defaults to all.
#'
#' @return A new MFA model object with rotated loadings (and factor scores if requested).
#'
#' @details
#' In R's \code{stats::varimax()}, the matrix rows are interpreted as variables, and
#' the columns as factors. If \code{Lambda_k} is already \code{(M x r)}, then passing it
#' directly works as intended (no transpose). The returned \code{loadings} is the same shape,
#' and \code{rotmat} is \code{(r x r)}, which can also be applied to the factor scores
#' (dimension \code{(T_i x r)}).
#'
#' @export
rotate_mfa_model <- function(mfa_fit,
                             rotation = c("varimax"),
                             rotate_scores = FALSE,
                             cluster_ids = NULL)
{
  rotation <- match.arg(rotation)
  K <- length(mfa_fit$Lambda)
  if (is.null(cluster_ids)) {
    cluster_ids <- seq_len(K)
  }

  rotated_fit <- mfa_fit
  z_vec <- rotated_fit$z  # cluster assignments
  if (is.null(z_vec)) {
    stop("mfa_fit$z not found. This is required to identify cluster assignments.")
  }

  # If we want to rotate scores but there's no factor_scores, error out
  if (rotate_scores && is.null(rotated_fit$factor_scores)) {
    stop("factor_scores not found in mfa_fit, but rotate_scores=TRUE.")
  }

  # For each cluster k, apply rotation
  for (k in cluster_ids) {

    Lambda_k <- rotated_fit$Lambda[[k]]  # (M x r)

    if (is.null(Lambda_k)) next

    # Use varimax on (M x r) => loadings has class "loadings", convert to matrix
    vout <- stats::varimax(Lambda_k)
    Lambda_rot <- unclass(vout$loadings)
    Lambda_rot <- as.matrix(Lambda_rot)  # remove "loadings" class
    rotated_fit$Lambda[[k]] <- Lambda_rot
    # rotation matrix (r x r)
    R_k <- as.matrix(vout$rotmat)

    # If requested, rotate factor scores for all subjects assigned to cluster k
    if (rotate_scores) {
      idx_k <- which(z_vec == k)
      for (i in idx_k) {
        Z_i <- rotated_fit$factor_scores[[i]]
        if (!is.null(Z_i)) {
          # (T_i x r) %*% (r x r) => (T_i x r)
          Z_i_rot <- Z_i %*% R_k
          rotated_fit$factor_scores[[i]] <- Z_i_rot
        }
      }
    }
  }

  return(rotated_fit)
}


#' Post-hoc Rotation for a Mixture PCA Model
#'
#' Applies a chosen rotation (e.g. varimax) to each cluster's loading matrix \code{W[[k]]}
#' in a fitted MPCA model, optionally also rotating the subject factor scores
#' in \code{factor_scores[[i]]}.
#'
#' @param mpca_fit A fitted MPCA model, e.g. returned by \code{mixture_pca_em_fit()}
#'   or \code{select_optimal_K_r_mpca()$best_model_info$model}, typically containing:
#'   \itemize{
#'     \item \code{W[[k]]}: (M x r) loadings for cluster k
#'     \item \code{factor_scores[[i]]}: optional (T_i x r) PC scores for subject i
#'     \item \code{z[i]}: cluster assignments
#'   }
#' @param rotation Rotation method, e.g. \code{"varimax"}.
#' @param rotate_scores Logical; if TRUE, also rotate factor_scores if present.
#' @param cluster_ids Vector of cluster indices to rotate. Defaults to all.
#'
#' @return A new MPCA model object with \code{W[[k]]} rotated, and \code{factor_scores}
#' also rotated if requested and available.
#'
#' @details
#' In R's \code{stats::varimax()}, the matrix should have rows = variables, cols = factors.
#' If \code{W_k} is already \code{(M x r)}, then passing it directly to varimax means
#' we have M variables, r factors. The returned rotation matrix is \code{(r x r)}, and
#' \code{scores_i (T_i x r)} can be multiplied by that.
#'
#' @export
rotate_mpca_model <- function(mpca_fit,
                              rotation = c("varimax"),
                              rotate_scores = FALSE,
                              cluster_ids = NULL)
{
  rotation <- match.arg(rotation)
  K <- length(mpca_fit$W)
  if (is.null(cluster_ids)) {
    cluster_ids <- seq_len(K)
  }

  rotated_fit <- mpca_fit
  z_vec <- rotated_fit$z
  if (is.null(z_vec)) {
    stop("mpca_fit$z not found. This is required to identify cluster assignments.")
  }

  if (rotate_scores && is.null(rotated_fit$factor_scores)) {
    stop("factor_scores not found in mpca_fit, but rotate_scores=TRUE.")
  }

  for (k in cluster_ids) {
    W_k <- rotated_fit$W[[k]]  # shape: (M x r)
    if (is.null(W_k)) next

    # If W_k is (M x r), pass it directly to varimax => rows= M variables, cols= r factors
    # No transpose needed if M is the number of variables, r = number of factors
    vout <- stats::varimax(W_k)
    # vout$loadings -> same shape as W_k, i.e. (M x r)
    W_rot <- unclass(vout$loadings)
    W_rot <- as.matrix(W_rot)
    rotated_fit$W[[k]] <- W_rot

    # rotation matrix is (r x r)
    R_k <- vout$rotmat

    if (rotate_scores) {
      idx_k <- which(z_vec == k)
      for (i in idx_k) {
        scores_i <- rotated_fit$factor_scores[[i]]
        if (!is.null(scores_i)) {
          # shape: (T_i x r) %*% (r x r) => (T_i x r)
          scores_rot <- scores_i %*% R_k
          rotated_fit$factor_scores[[i]] <- scores_rot
        }
      }
    }
  }

  return(rotated_fit)
}


#' Compute Factor Scores for a Single Subject in MFA
#'
#' Given a single subject's data \code{X} and the cluster's parameters
#' (\code{mu, Lambda, Psi}), this function computes the posterior factor scores
#' (of size \code{T x r}) in a standard Factor Analysis approach.
#'
#' @param X A \code{(T x M)} data matrix for one subject/time-series.
#' @param mu A length-\code{M} mean vector for the cluster.
#' @param Lambda An \code{(M x r)} factor loading matrix for the cluster.
#' @param Psi An \code{(M x M)} diagonal covariance matrix (unique variances).
#'
#' @return A \code{(T x r)} matrix of factor scores \code{Z_i(t)}, where \code{t = 1..T}.
#'
#' @details
#' We use the standard formula for Factor Analysis posterior means:
#' \deqn{
#'   \mathbf{Z}_i(t) = \left( \mathbf{I}_r + \Lambda^\top \Psi^{-1} \Lambda \right)^{-1}
#'   \Lambda^\top \Psi^{-1} \left[\mathbf{x}_i(t) - \mu\right].
#' }
#' The code is implemented in a vectorized manner. If \eqn{r = \mathrm{ncol}(\Lambda)},
#' then the output dimension is \code{(T x r)}.
#'
#' @examples
#' \dontrun{
#' X_i <- list_of_data[[i]]  # (T x M)
#' # from model$mu[[k]], model$Lambda[[k]], model$Psi[[k]]
#' Z_i <- factor_scores_mfa_single(X_i, mu_k, Lambda_k, Psi_k)
#' dim(Z_i)  # (T x r)
#' }
#'
#' @export
factor_scores_mfa_single <- function(X, mu, Lambda, Psi) {
  # X: (T x M)
  # mu: (M)
  # Lambda: (M x r)
  # Psi: (M x M) diagonal

  T_i <- nrow(X)
  M   <- ncol(X)
  r_  <- ncol(Lambda)

  # check dimension
  if(length(mu) != M) {
    stop("mu dimension mismatch: length(mu) != M")
  }
  if(!all(dim(Psi) == c(M,M))) {
    stop("Psi must be (M x M).")
  }
  # Invert diagonal Psi
  invPsi <- diag(1 / diag(Psi), nrow=M)
  # A = I_r + Lambda^T invPsi Lambda => (r x r)
  A <- diag(r_) + t(Lambda) %*% invPsi %*% Lambda
  A_inv <- solve(A)  # (r x r)

  # Let's do a vectorized row-by-row approach
  Xc <- sweep(X, 2, mu, FUN = "-")  # (T x M)
  # Then factor scores = Xc %*% t(W) %*% t(A_inv)?
  # But let's define W = t(Lambda) %*% invPsi => shape: (r x M)
  # We want => (T x r) = (T x M)*(M x r) basically, times A_inv as well
  W <- t(Lambda) %*% invPsi  # (r x M)
  Wt <- t(W)                 # (M x r)
  # => intermediate: XcW^T = (T x M)*(M x r) => (T x r)
  temp <- Xc %*% Wt
  # multiply by A_inv on the right => (T x r) %*% (r x r) => (T x r)
  Z_mat <- temp %*% A_inv

  # ensure it's matrix
  Z_mat <- as.matrix(Z_mat)
  return(Z_mat)
}


#' Compute Factor Scores for All Subjects in an MFA Model
#'
#' Given a fitted MFA model (with \code{z, mu, Lambda, Psi} for each cluster),
#' this function computes the posterior factor scores \code{Z_i(t)} for each subject \code{i}
#' and stores them in \code{model$factor_scores[[i]]}.
#'
#' @param mfa_fit A fitted MFA model, typically from \code{mfa_em_fit()} or
#'   \code{select_optimal_K_r_mfa()$best_model_info$model}, containing:
#'   \itemize{
#'     \item \code{z[i]}: cluster assignments,
#'     \item \code{Lambda[[k]]}, \code{Psi[[k]]}, \code{mu[[k]]} for each cluster k.
#'   }
#' @param list_of_data The same list of data \code{(T_i x M)} used to fit the model.
#'
#' @return The same \code{mfa_fit} object, but with a new element \code{factor_scores}
#'   (a list of length \code{N}), where each \code{factor_scores[[i]]} is a
#'   \code{(T_i x r)} matrix of factor scores for subject \code{i}.
#'
#' @details
#' Internally, it calls \code{\link{factor_scores_mfa_single}} for each subject,
#' using that subject's assigned cluster \code{k = z[i]}.
#'
#' @examples
#' \dontrun{
#' mfa_fit <- compute_factor_scores_mfa(mfa_fit, list_of_data)
#' # Now mfa_fit$factor_scores[[i]] is available for each subject i
#' }
#' @export
compute_factor_scores_mfa <- function(mfa_fit, list_of_data) {
  N <- length(list_of_data)
  if (N < 1) stop("No data in list_of_data.")
  if (is.null(mfa_fit$z)) {
    stop("mfa_fit$z not found; we need cluster assignments per subject.")
  }
  if (is.null(mfa_fit$Lambda)) {
    stop("mfa_fit$Lambda not found; we need loadings per cluster.")
  }
  if (is.null(mfa_fit$Psi)) {
    stop("mfa_fit$Psi not found; we need noise diag matrix per cluster.")
  }
  if (is.null(mfa_fit$mu)) {
    stop("mfa_fit$mu not found; we need mean vectors per cluster.")
  }

  factor_scores_list <- vector("list", N)
  z_vec <- mfa_fit$z
  K <- length(mfa_fit$Lambda)

  # For each subject i
  for(i in seq_len(N)) {
    k_i <- z_vec[i]
    if(k_i < 1 || k_i > K) {
      stop(sprintf("Subject %d has out-of-range cluster = %d", i, k_i))
    }
    X_i <- list_of_data[[i]]

    Lambda_k <- mfa_fit$Lambda[[k_i]]
    Psi_k    <- mfa_fit$Psi[[k_i]]
    mu_k     <- mfa_fit$mu[[k_i]]

    Z_i <- factor_scores_mfa_single(X_i, mu_k, Lambda_k, Psi_k)
    factor_scores_list[[i]] <- Z_i
  }

  mfa_fit$factor_scores <- factor_scores_list
  return(mfa_fit)
}


#' Compute Factor (PC) Scores for a Single Subject in a Mixture PCA Cluster
#'
#' Given a single subject's data X, and the cluster parameters (mu, W, sigma2),
#' this function computes the posterior factor scores (T x r) under a PPCA assumption
#' with isotropic noise sigma^2 * I.
#'
#' @param X A (T x M) data matrix for one subject/time-series.
#' @param mu A length-M mean vector.
#' @param W An (M x r) loading (projection) matrix for the cluster.
#' @param sigma2 A scalar representing isotropic noise variance (sigma^2).
#'
#' @return A (T x r) matrix of factor (PC) scores for this subject.
#'
#' @details
#' We use the posterior mean formula for PPCA with isotropic noise:
#' \deqn{
#'   z(t) = (W^\top W + \sigma^2 I_r)^{-1} \, W^\top [ x(t) - \mu ].
#' }
#' This is applied row by row in X, but we can do it in a vectorized manner.
#'
#' If your noise model is not purely \sigma^2 I, you will need a more general approach.
#' Also note that if T=1, the result might be a (1 x r) matrix or (r)-vector, so
#' be sure to keep it as a matrix if desired.
#'
#' @export
compute_factor_scores_mpca_single <- function(X, mu, W, sigma2) {
  # X: (T x M)
  # mu: (M)
  # W: (M x r)
  # sigma2: numeric scalar

  T_i <- nrow(X)
  M   <- ncol(X)
  r_  <- ncol(W)

  if(length(mu) != M) {
    stop("mu dimension mismatch: length(mu) != M")
  }
  # W is (M x r), so W^T W is (r x r)
  if(!is.numeric(sigma2) || length(sigma2)!=1) {
    stop("sigma2 must be a scalar numeric.")
  }

  # A = W^T W + sigma2 I_r  => (r x r)
  A <- t(W) %*% W + diag(sigma2, r_)
  A_inv <- solve(A)  # (r x r)

  # let's do a vectorized approach:
  # Xc = X - mu => shape (T x M)
  Xc <- sweep(X, 2, mu, FUN="-")
  # W^T => (r x M), so t(W) is (r x M)
  # step1 = Xc %*% W => (T x r?)
  # wait, Xc is (T x M), W is (M x r). => Xc %*% W => (T x r)
  step1 <- Xc %*% W

  # z(t) = step1(t,.) %*% A_inv => => (T x r)
  # => step1 %*% A_inv => (T x r)
  Z_mat <- step1 %*% A_inv

  # ensure it's matrix
  Z_mat <- as.matrix(Z_mat)
  return(Z_mat)
}


#' Compute PC Factor Scores for All Subjects in a Mixture PCA Model
#'
#' Given a fitted MPCA model (with \code{z, mu, W, sigma2} for each cluster),
#' this function computes the posterior factor (PC) scores \code{Z_i(t)} for each subject \code{i}
#' and stores them in \code{mpca_fit$factor_scores[[i]]}.
#'
#' @param mpca_fit A fitted MPCA model from \code{mixture_pca_em_fit()} or
#'   \code{select_optimal_K_r_mpca()$best_model_info$model}, typically containing:
#'   \itemize{
#'     \item \code{z[i]}: cluster assignments (1..K),
#'     \item \code{mu[[k]]}, \code{W[[k]]}, \code{sigma2[k]} for each cluster,
#'     \item possibly \code{P[[k]], D[[k]]} instead of \code{W[[k]]}.
#'   }
#' @param list_of_data The same list of data \code{(T_i x M)} used to fit the model.
#'
#' @return The same \code{mpca_fit} object, but with a new element \code{factor_scores},
#'   where \code{factor_scores[[i]]} is a \code{(T_i x r)} matrix of PC scores for subject \code{i}.
#'
#' @details
#' We use \code{\link{compute_factor_scores_mpca_single}} for each subject, with the
#' appropriate cluster parameters. This assumes isotropic noise \code{sigma2[k]} for cluster k.
#'
#' @examples
#' \dontrun{
#' mpca_fit <- compute_factor_scores_mpca(mpca_fit, list_of_data)
#' mpca_fit$factor_scores[[1]]  # (T_1 x r) PC scores
#' }
#' @export
compute_factor_scores_mpca <- function(mpca_fit, list_of_data) {
  N <- length(list_of_data)
  if (N < 1) stop("No data in list_of_data.")
  if (is.null(mpca_fit$z)) {
    stop("mpca_fit$z not found; need cluster assignments.")
  }
  if (is.null(mpca_fit$W)) {
    stop("mpca_fit$W not found; need (M x r) loadings per cluster.")
  }
  if (is.null(mpca_fit$sigma2)) {
    stop("mpca_fit$sigma2 not found; need noise variance per cluster.")
  }
  if (is.null(mpca_fit$mu)) {
    stop("mpca_fit$mu not found; need mean vectors per cluster.")
  }

  factor_scores_list <- vector("list", N)
  z_vec <- mpca_fit$z
  K <- length(mpca_fit$W)

  for(i in seq_len(N)) {
    k_i <- z_vec[i]
    if(k_i < 1 || k_i > K) {
      stop(sprintf("Subject %d has out-of-range cluster = %d", i, k_i))
    }
    X_i <- list_of_data[[i]]

    mu_k  <- mpca_fit$mu[[k_i]]
    W_k   <- mpca_fit$W[[k_i]]
    sig2_k<- mpca_fit$sigma2[k_i]

    # compute factor scores using the single-subject function
    Z_i <- compute_factor_scores_mpca_single(X_i, mu_k, W_k, sig2_k)
    factor_scores_list[[i]] <- Z_i
  }

  mpca_fit$factor_scores <- factor_scores_list
  return(mpca_fit)
}


