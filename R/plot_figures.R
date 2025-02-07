#' Plot ARI vs. Separation/Noise (Figure 1)
#'
#' Given a result data frame with columns like \code{ARI_MFA, ARI_PCA, sep, noise},
#' this function reshapes them into long format (MFA vs. MPCA) and plots ARI vs. separation,
#' color-coded by noise, using \code{ggplot2}.
#'
#' @param df A \code{data.frame} with columns \code{K_true, r_true, N, M, sep, noise, ARI_MFA, ARI_PCA}.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' p <- plot_fig1_ARI(df_merged)
#' ggsave("figure1_ARI_sep_noise.pdf", p)
#' }
#'
#' @export
plot_fig1_ARI <- function(df){
  df_long <- df %>%
    dplyr::select(K_true, r_true, sep, noise, ARI_MFA, ARI_PCA) %>%
    tidyr::pivot_longer(cols=c("ARI_MFA","ARI_PCA"), 
                        names_to="Method", values_to="ARI")
  
  p <- ggplot(df_long, aes(x=factor(sep), y=ARI, color=factor(noise),
                           group=interaction(Method, noise))) +
    geom_point(position=position_dodge(width=0.3)) +
    geom_line(position=position_dodge(width=0.3)) +
    facet_wrap(~Method + K_true + r_true + N + M, labeller=label_both) +
    theme_bw() +
    labs(x="Separation (sep)", y="ARI", color="Noise",
         title="(Figure 1) ARI vs sep/noise (MFA, MPCA)")
  p
}


#' Plot BIC Heatmap for a Subset (Figure 2)
#'
#' This function filters the data for a specific condition, e.g. \code{N=50, M=8, sep=1.0, noise=1.0},
#' and plots a heatmap of \code{BIC_MFA} with \code{r_true} on the x-axis and \code{K_true} on the y-axis.
#'
#' @param df A \code{data.frame} with columns including \code{BIC_MFA, K_true, r_true}.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' p <- plot_fig2_BIC_heatmap(df_merged)
#' ggsave("figure2_BIC_heatmap.pdf", p)
#' }
#'
#' @export
plot_fig2_BIC_heatmap <- function(df){
  sub <- df %>% 
    dplyr::filter(N==50, M==8, abs(sep-1.0)<1e-9, abs(noise-1.0)<1e-9)
  
  p <- ggplot(sub, aes(x=factor(r_true), y=factor(K_true), fill=BIC_MFA)) +
    geom_tile() +
    scale_fill_viridis_c() +
    theme_bw() +
    labs(x="r", y="K", fill="BIC",
         title="(Figure 2) BIC heatmap (MFA, N=50, M=8, sep=1.0, noise=1.0)")
  p
}


#' Plot Boxplot of ARI for Multiple Initial Seeds (Figure 3)
#'
#' If the data frame contains multiple runs with the same parameters but different
#' initial seeds, this function can show a boxplot of ARI distribution for MFA and MPCA.
#'
#' @param df_multiple A \code{data.frame} with columns \code{ARI_MFA, ARI_PCA}.
#'
#' @return A \code{ggplot} object.
#'
#' @export
plot_fig3_init_boxplot <- function(df_multiple){
  df_long <- df_multiple %>%
    tidyr::pivot_longer(cols=c("ARI_MFA","ARI_PCA"), 
                        names_to="Method", values_to="ARI")
  
  p <- ggplot(df_long, aes(x=Method, y=ARI, color=Method)) +
    geom_boxplot() +
    theme_bw() +
    labs(title="(Figure 3) Boxplot for multiple initial seeds")
  p
}


#' Plot Factor Loadings Heatmap (Figure 4)
#'
#' Given an MFA fit object containing \code{Lambda}, this function plots a heatmap
#' for the loadings in a specific cluster.
#'
#' @param mfa_fit A fitted MFA model with \code{$Lambda}.
#' @param cluster_id Which cluster's loadings to plot (1..K).
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' fit_mfa <- readRDS("some_fit.rds")
#' p <- plot_fig4_loadings_heatmap(fit_mfa, cluster_id=1)
#' ggsave("figure4_loadings_heatmap.pdf", p)
#' }
#'
#' @export
plot_fig4_loadings_heatmap <- function(mfa_fit, cluster_id=1){
  if(is.null(mfa_fit$Lambda)){
    stop("mfa_fit does not contain Lambda.")
  }
  Lambda_k <- mfa_fit$Lambda[[cluster_id]]
  if(is.null(Lambda_k)){
    stop("No cluster_id = ", cluster_id, " in Lambda.")
  }
  df_melt <- reshape2::melt(Lambda_k)
  colnames(df_melt) <- c("Muscle","Factor","Loading")
  
  p <- ggplot(df_melt, aes(x=factor(Factor), y=factor(Muscle), fill=Loading)) +
    geom_tile() +
    scale_fill_gradient2(mid="white", low="blue", high="red") +
    theme_bw() +
    labs(x="Factor", y="Muscle", fill="Loading",
         title=sprintf("(Figure 4) Cluster %d synergy loadings", cluster_id))
  p
}
