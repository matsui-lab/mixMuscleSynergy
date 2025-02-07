#' Compute Global VAF (1 - SSE/SST) for an MFA Model
#'
#' This function takes a fitted MFA model (with elements \code{z, mu, Lambda, Psi})
#' and reconstructs all subjects' data from the factor scores, computing the
#' ratio of explained variance (\eqn{1 - \mathrm{SSE}/\mathrm{SST}}).
#'
#' @param list_of_data A list of matrices, each \code{(T_i x M)}, containing the original data.
#' @param res A fitted MFA model result. It must contain:
#'   \itemize{
#'     \item{\code{z}: A length-N vector of cluster assignments (1..K).}
#'     \item{\code{mu}: A list of length K, each a mean vector of length M.}
#'     \item{\code{Lambda}: A list of length K, each an \code{M x r} loading matrix.}
#'     \item{\code{Psi}: A list of length K, each an \code{M x M} diagonal matrix.}
#'   }
#' @return A numeric scalar for the global variance accounted for, \eqn{1 - \mathrm{SSE}/\mathrm{SST}}.
#'
#' @details
#' For each subject \code{i}, we identify its cluster assignment \code{k_i}, compute factor scores
#' via \code{\link{compute_factor_scores}}, reconstruct \code{Xhat} from the factor model,
#' accumulate the sum of squares of residuals (SSE), and the sum of squares total (SST).
#'
#' @seealso \code{\link{compute_factor_scores}}
#'
#' @examples
#' \dontrun{
#' gvaf <- compute_global_vaf_mfa(my_data_list, fit_mfa)
#' cat("Global VAF:", gvaf, "\n")
#' }
#'
#' @export
compute_global_vaf_mfa <- function(list_of_data, res){
  # res: assumed to contain z, mu, Lambda, Psi (lists)
  X_all <- NULL
  Xhat_all <- NULL
  z_vec <- res$z
  N <- length(list_of_data)
  
  for(i in seq_len(N)){
    k_i <- z_vec[i]  # cluster ID for subject i
    mu_k <- res$mu[[k_i]]
    Lambda_k <- res$Lambda[[k_i]]
    Psi_k <- res$Psi[[k_i]]
    
    # observed data
    X_i <- list_of_data[[i]]  # (T_i x M)
    
    # factor scores
    EZ_i <- compute_factor_scores(X_i, mu_k, Lambda_k, Psi_k)
    
    # reconstruct
    Xhat_i <- EZ_i %*% t(Lambda_k)
    Xhat_i <- sweep(Xhat_i, 2, mu_k, FUN="+")
    
    X_all    <- rbind(X_all, X_i)
    Xhat_all <- rbind(Xhat_all, Xhat_i)
  }
  
  SSE_mat <- (X_all - Xhat_all)^2
  SSE <- sum(SSE_mat)
  SST <- sum(X_all^2)
  
  vaf_global <- 1 - SSE/SST
  vaf_global
}


#' Compute Factor Scores for a Single MFA Cluster
#'
#' Given a single cluster's parameters (\code{mu, Lambda, Psi}) and data \code{X},
#' this function returns the estimated factor scores. 
#'
#' @param X A numeric matrix of size \code{(T x M)}.
#' @param mu A length-M mean vector.
#' @param Lambda An \code{M x r} loading matrix.
#' @param Psi An \code{M x M} diagonal matrix for unique variances.
#'
#' @return A \code{(T x r)} matrix of factor scores.
#'
#' @details
#' This follows the standard posterior mean formula in factor analysis:
#' \eqn{\mathrm{E}[Z|X] = (I + \Lambda^T \Psi^{-1} \Lambda)^{-1} \Lambda^T \Psi^{-1} (X - \mu)}.
#'
#' @seealso \code{\link{compute_global_vaf_mfa}}
#'
#' @examples
#' \dontrun{
#' # Suppose we have X (T x M), mu (M), Lambda (M x r), Psi (M x M).
#' scores <- compute_factor_scores(X, mu, Lambda, Psi)
#' head(scores)
#' }
#'
#' @export
compute_factor_scores <- function(X, mu, Lambda, Psi){
  Xc <- sweep(X, 2, mu, FUN="-")  # center
  Sigma <- Lambda %*% t(Lambda) + Psi
  diag(Sigma) <- diag(Sigma) + 1e-8  # stabil
  
  L <- chol(Sigma)
  invL <- solve(L, diag(nrow(L)))
  invSigma <- t(invL) %*% invL
  
  W <- t(Lambda) %*% invSigma
  A <- diag(ncol(Lambda)) + W %*% Lambda
  A_chol <- chol(A)
  A_inv <- solve(A_chol, diag(ncol(A)))
  A_inv <- t(A_inv) %*% A_inv
  
  EZt <- A_inv %*% W %*% t(Xc)
  EZ  <- t(EZt)
  EZ
}


#' Compute Cluster Sizes from a Fitted MFA (or Mixture PCA) Model
#'
#' Given a result with a vector \code{z} of length N (the cluster assignments),
#' this function returns a frequency table, i.e., how many subjects are assigned
#' to each cluster.
#'
#' @param res A fitted mixture model object that contains \code{z}.
#'
#' @return A table of cluster frequencies.
#'
#' @examples
#' \dontrun{
#' tab <- compute_cluster_sizes(fit_mfa)
#' print(tab)
#' }
#' @export
compute_cluster_sizes <- function(res){
  z_vec <- res$z
  tab <- table(z_vec)
  tab
}


#' Post-hoc Evaluation of MFA Grid Search
#'
#' Given the output of \code{\link{select_optimal_K_r_mfa}}, this function computes
#' additional metrics for each fitted model in \code{all_models}, such as
#' global VAF (\code{\link{compute_global_vaf_mfa}}) and cluster sizes
#' (\code{\link{compute_cluster_sizes}}).
#'
#' @param list_of_data The same list of data used to fit the models.
#' @param selection_obj The object returned by \code{\link{select_optimal_K_r_mfa}}.
#'
#' @return A data frame with columns:
#' \item{K}{Number of clusters.}
#' \item{r}{Factor dimension.}
#' \item{logLik}{Log-likelihood from the fit.}
#' \item{BIC}{BIC from the fit.}
#' \item{GlobalVAF}{Global variance accounted for, \eqn{1 - \mathrm{SSE}/\mathrm{SST}}.}
#' \item{minClusterSize}{Minimum cluster size.}
#' \item{maxClusterSize}{Maximum cluster size.}
#'
#' @details
#' This function loops over each entry in \code{selection_obj$all_models}, 
#' extracts the model, calls \code{\link{compute_global_vaf_mfa}} and
#' \code{\link{compute_cluster_sizes}}, and records the results in a combined table.
#'
#' @seealso \code{\link{select_optimal_K_r_mfa}}, \code{\link{compute_factor_scores}}
#'
#' @examples
#' \dontrun{
#' sel_obj <- select_optimal_K_r_mfa(my_data, Kvec=2:4, rvec=1:3)
#' df_posthoc <- posthoc_mfa_evaluation(my_data, sel_obj)
#' head(df_posthoc)
#' }
#'
#' @export
posthoc_mfa_evaluation <- function(list_of_data, selection_obj){
  all_mods <- selection_obj$all_models
  if(is.null(all_mods) || length(all_mods)==0){
    stop("No models found in selection_obj$all_models")
  }
  df_list <- list()
  
  for(i in seq_along(all_mods)){
    obj_i <- all_mods[[i]]
    K_val <- obj_i$K
    r_val <- obj_i$r
    fit   <- obj_i$model
    
    # 1) Global VAF
    gvaf_i <- compute_global_vaf_mfa(list_of_data, fit)
    # 2) cluster sizes
    size_tab <- compute_cluster_sizes(fit)
    min_size <- min(size_tab)
    max_size <- max(size_tab)
    
    df_list[[i]] <- data.frame(
      K    = K_val,
      r    = r_val,
      logLik = obj_i$logLik,
      BIC    = obj_i$BIC,
      GlobalVAF = gvaf_i,
      minClusterSize = min_size,
      maxClusterSize = max_size
    )
  }
  
  df_posthoc <- dplyr::bind_rows(df_list)
  df_posthoc <- df_posthoc[order(df_posthoc$BIC), ]
  rownames(df_posthoc) <- NULL
  df_posthoc
}


#' Retrieve a Specific Model by (K, r) from a Grid Search
#'
#' After calling \code{\link{select_optimal_K_r_mfa}} (or a similar function),
#' you have \code{$all_models} containing fits for various \code{(K, r)}. This
#' function extracts the model for a user-specified pair \code{(K_target, r_target)}.
#'
#' @param selection_obj The object returned by \code{\link{select_optimal_K_r_mfa}}.
#' @param K_target Desired number of clusters.
#' @param r_target Desired factor dimension.
#'
#' @return The fitted model (the same structure returned by \code{mfa_em_fit}),
#'   or \code{NULL} if no match found.
#'
#' @examples
#' \dontrun{
#' sel_obj <- select_optimal_K_r_mfa(my_data, Kvec=2:4, rvec=1:3)
#' best_model <- get_model_by_K_r(sel_obj, K_target=3, r_target=2)
#' }
#'
#' @export
get_model_by_K_r <- function(selection_obj, K_target, r_target){
  all_mods <- selection_obj$all_models
  if(is.null(all_mods) || length(all_mods)==0) return(NULL)
  
  idx <- which(
    sapply(all_mods, function(x) x$K) == K_target &
      sapply(all_mods, function(x) x$r) == r_target
  )
  if(length(idx)<1){
    message(sprintf("No model found for K=%d, r=%d.", K_target, r_target))
    return(NULL)
  }
  all_mods[[ idx[1] ]]$model
}



#' Compute Global VAF (1 - SSE/SST) for Mixture PCA
#'
#' Similar to \code{\link{compute_global_vaf_mfa}}, but uses a Mixture PCA model structure.
#' Given a model with \code{z, mu, P, D, Psi}, we reconstruct each subject's data from
#' the principal component scores, accumulate SSE, and compare to total sum of squares (SST).
#'
#' @param list_of_data A list of \code{(T_i x M)} data matrices.
#' @param res The fitted Mixture PCA model, typically from \code{\link{mixture_pca_em_fit}}.
#'   Must contain:
#'   \itemize{
#'     \item{\code{z}: cluster assignments (length N).}
#'     \item{\code{mu}: list of length K, each \code{(M)} vector.}
#'     \item{\code{P}: list of length K, each \code{(M x r)}.}
#'     \item{\code{D}: list of length K, each \code{(r x r)} diagonal.}
#'     \item{\code{Psi}: list of length K, each \code{(M x M)} diagonal.}
#'   }
#'
#' @return A numeric scalar, \eqn{1 - \mathrm{SSE}/\mathrm{SST}} across all subjects and time steps.
#'
#' @details
#' The user must define how to reconstruct data from PCA parameters. Typically,
#' for each cluster \code{k}, we have \code{W = P * D}, so we either compute
#' posterior scores or do a direct projection. See \code{\link{compute_factor_scores_pca}} for details.
#'
#' @export
compute_global_vaf_mpca <- function(list_of_data, res){
  # placeholders
  stop("Please define how to reconstruct each subject from (P, D, Psi).")
}


#' Compute Factor (PC) Scores for a Single Mixture PCA Cluster
#'
#' Given \code{mu, P, D, Psi} for one cluster, and data \code{X}, estimate the principal component scores.
#'
#' @param X A matrix \code{(T x M)}.
#' @param mu A length-M mean vector.
#' @param P An \code{(M x r)} matrix of principal directions (orthonormal columns).
#' @param D An \code{(r x r)} diagonal matrix of singular values or scaling factors.
#' @param Psi An \code{(M x M)} diagonal matrix for residual noise.
#'
#' @return A \code{(T x r)} matrix of PC scores.
#'
#' @details
#' One possible formula (for PPCA) is the posterior mean of the latent Z. The user
#' can adapt from \code{\link{compute_factor_scores}} or from a known PPCA posterior formula.
#'
#' @export
compute_factor_scores_pca <- function(X, mu, P, D, Psi){
  stop("Please implement the logic to compute PC scores from (P, D, Psi).")
}


#' Compute Cluster Sizes for Mixture PCA
#'
#' Like \code{\link{compute_cluster_sizes}}, but for a Mixture PCA model result.
#'
#' @param res A fitted mixture PCA model result, containing \code{$z}.
#'
#' @return A table of frequencies of cluster assignments.
#'
#' @export
compute_cluster_sizes_mpca <- function(res){
  # same as MFA, basically
  z_vec <- res$z
  table(z_vec)
}


#' Post-hoc Evaluation of Mixture PCA Grid Search
#'
#' Takes the output of \code{\link{select_optimal_K_r_mpca}} and computes additional
#' metrics (global VAF, cluster sizes, etc.) for each entry in \code{all_models}.
#'
#' @param list_of_data The same list of data used for the grid search.
#' @param selection_obj The object returned by \code{\link{select_optimal_K_r_mpca}}.
#'
#' @return A data frame with columns \code{K, r, logLik, BIC, GlobalVAF, minClusterSize, maxClusterSize}.
#'
#' @export
posthoc_mpca_evaluation <- function(list_of_data, selection_obj){
  all_mods <- selection_obj$all_models
  if(is.null(all_mods) || length(all_mods)==0){
    stop("No models found in selection_obj$all_models")
  }
  
  df_list <- list()
  
  for(i in seq_along(all_mods)){
    obj_i <- all_mods[[i]]
    K_val <- obj_i$K
    r_val <- obj_i$r
    fit   <- obj_i$model  # mixture PCA fit
    
    # compute global VAF
    gvaf_i <- compute_global_vaf_mpca(list_of_data, fit)
    
    # cluster sizes
    size_tab <- compute_cluster_sizes_mpca(fit)
    min_size <- min(size_tab)
    max_size <- max(size_tab)
    
    df_list[[i]] <- data.frame(
      K = K_val,
      r = r_val,
      logLik = obj_i$logLik,
      BIC = obj_i$BIC,
      GlobalVAF = gvaf_i,
      minClusterSize = min_size,
      maxClusterSize = max_size
    )
  }
  
  df_posthoc <- dplyr::bind_rows(df_list)
  df_posthoc <- df_posthoc[order(df_posthoc$BIC), ]
  rownames(df_posthoc) <- NULL
  df_posthoc
}


#' Retrieve a Mixture PCA Model by (K, r)
#'
#' Given a \code{selection_obj} from \code{\link{select_optimal_K_r_mpca}},
#' extracts the model for the specified \code{(K_target, r_target)} if it exists.
#'
#' @param selection_obj The output of \code{\link{select_optimal_K_r_mpca}}.
#' @param K_target Desired K.
#' @param r_target Desired r.
#'
#' @return The fitted model (same structure as from \code{mixture_pca_em_fit}), or \code{NULL}.
#'
#' @export
get_model_by_K_r_mpca <- function(selection_obj, K_target, r_target){
  all_mods <- selection_obj$all_models
  if(is.null(all_mods) || length(all_mods)==0) return(NULL)
  
  idx <- which(
    sapply(all_mods, function(x) x$K) == K_target &
      sapply(all_mods, function(x) x$r) == r_target
  )
  if(length(idx)<1){
    message(sprintf("No PCA model found for K=%d, r=%d.", K_target, r_target))
    return(NULL)
  }
  all_mods[[ idx[1] ]]$model
}

