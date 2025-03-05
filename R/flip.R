#' Align (Flip) Factor Signs in an MFA Model
#'
#' After running \code{\link{rotate_mfa_model}} or any rotation, each factor's
#' sign can be arbitrary. This function enforces a rule such as:
#' "The variable with the largest absolute loading in each factor must be positive."
#' Then, if that factor's column is flipped, we also flip the corresponding column
#' in the factor scores, to keep the model consistent.
#'
#' @param mfa_fit A fitted MFA model (with \code{Lambda[[k]]}, \code{factor_scores[[i]]}, \code{z}, etc.).
#' @param method A string specifying how to choose the sign. Currently only
#'   \code{"max_abs"} is implemented, which flips the column if the variable with
#'   the largest absolute loading is negative.
#' @param cluster_ids Which clusters to process. Defaults to all.
#'
#' @return The same \code{mfa_fit} object, but with possibly flipped signs in
#'   \code{Lambda} columns and \code{factor_scores} columns.
#'
#' @details
#' - For each cluster \code{k}, we loop over each factor (column). We find the row
#'   index \code{i_max} that has the largest absolute loading. If
#'   \code{Lambda[k][i_max, factor] < 0}, we multiply that entire column by \code{-1}.
#' - If \code{factor_scores} is present, we also multiply the corresponding column
#'   in each subject's factor score by \code{-1} to maintain consistency.
#' - This doesn't change the model fit or BIC, just the sign convention.
#'
#' @export
align_factor_signs_mfa <- function(mfa_fit,
                                   method = c("max_abs"),
                                   cluster_ids = NULL)
{
  method <- match.arg(method)
  K <- length(mfa_fit$Lambda)

  if (is.null(cluster_ids)) {
    cluster_ids <- seq_len(K)
  }

  # Check if we have factor_scores
  has_scores <- !is.null(mfa_fit$factor_scores)
  z_vec <- mfa_fit$z

  for (k in cluster_ids) {
    Lambda_k <- mfa_fit$Lambda[[k]]
    if (is.null(Lambda_k)) next

    M <- nrow(Lambda_k)
    r_ <- ncol(Lambda_k)

    for (col_j in seq_len(r_)) {
      # 1) find row i_max that has the largest absolute value in this factor
      abs_vals <- abs(Lambda_k[, col_j])
      i_max <- which.max(abs_vals)
      # 2) if that entry is negative => flip sign
      if (Lambda_k[i_max, col_j] < 0) {
        # flip entire column j in Lambda_k
        Lambda_k[, col_j] <- -Lambda_k[, col_j]

        # if factor_scores exist => also flip that column j for all subjects i in cluster k
        if (has_scores) {
          idx_k <- which(z_vec == k)
          for (i_sub in idx_k) {
            # e.g. (T_i x r)
            scores_i <- mfa_fit$factor_scores[[i_sub]]
            if (!is.null(scores_i)) {
              scores_i[, col_j] <- -scores_i[, col_j]
              mfa_fit$factor_scores[[i_sub]] <- scores_i
            }
          }
        }
      }
    }
    mfa_fit$Lambda[[k]] <- Lambda_k
  }

  return(mfa_fit)
}


#' Align (Flip) Factor Signs in a Mixture PCA Model
#'
#' Similar to \code{\link{align_factor_signs_mfa}}, but for a Mixture PCA model
#' whose loadings are in \code{W[[k]]}.
#'
#' @param mpca_fit A fitted MPCA model (with \code{W[[k]]}, \code{factor_scores[[i]]}, \code{z}, etc.).
#' @param method A string specifying how to choose the sign. Defaults to \code{"max_abs"}.
#' @param cluster_ids Which clusters to process. Defaults to all.
#'
#' @return The same \code{mpca_fit} object, but with possibly flipped signs in
#'   \code{W} columns and \code{factor_scores} columns (if present).
#'
#' @details
#' - For each cluster \code{k}, we loop over each factor (column) in \code{W[[k]]}.
#'   If the entry with largest absolute loading is negative, we multiply that entire column by \code{-1}.
#' - If \code{factor_scores} is present, we also multiply the corresponding column
#'   in each subject's factor scores by \code{-1}.
#'
#' @export
align_factor_signs_mpca <- function(mpca_fit,
                                    method = c("max_abs"),
                                    cluster_ids = NULL)
{
  method <- match.arg(method)
  K <- length(mpca_fit$W)
  if (is.null(cluster_ids)) {
    cluster_ids <- seq_len(K)
  }

  has_scores <- !is.null(mpca_fit$factor_scores)
  z_vec <- mpca_fit$z

  for (k in cluster_ids) {
    W_k <- mpca_fit$W[[k]]
    if (is.null(W_k)) next

    M <- nrow(W_k)
    r_ <- ncol(W_k)

    for (col_j in seq_len(r_)) {
      abs_vals <- abs(W_k[, col_j])
      i_max <- which.max(abs_vals)
      if (W_k[i_max, col_j] < 0) {
        # flip sign
        W_k[, col_j] <- -W_k[, col_j]

        # flip factor_scores
        if (has_scores) {
          idx_k <- which(z_vec == k)
          for (i_sub in idx_k) {
            scores_i <- mpca_fit$factor_scores[[i_sub]]
            if (!is.null(scores_i)) {
              scores_i[, col_j] <- -scores_i[, col_j]
              mpca_fit$factor_scores[[i_sub]] <- scores_i
            }
          }
        }
      }
    }
    mpca_fit$W[[k]] <- W_k
  }

  return(mpca_fit)
}
