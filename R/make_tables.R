#' Create a Static Table of Simulation Conditions (Table 1)
#'
#' This function returns a small \code{data.frame} describing the parameter ranges
#' used in the simulation (K, r, N, M, etc.). Typically printed via \code{kable()} in a vignette or paper.
#'
#' @return A \code{data.frame} with columns \code{Parameter} and \code{Values}.
#'
#' @examples
#' \dontrun{
#' df_table1 <- make_table1()
#' knitr::kable(df_table1, caption="Table 1: Conditions")
#' }
#'
#' @export
make_table1 <- function(){
  df_table1 <- data.frame(
    Parameter = c(
      "True cluster number (K)",
      "True factor (PC) number (r)",
      "Number of subjects (N)",
      "Channels (M)",
      "Time length (T_each)",
      "Separation parameter (sep)",
      "Noise scale (noise)"
    ),
    Values = c(
      "3, 4",
      "3, 5",
      "50, 100",
      "8, 12",
      "100",
      "0.5, 1.0",
      "1.0, 2.0"
    ),
    stringsAsFactors=FALSE
  )
  df_table1
}


#' Summarize Results for a Subset (Table 2)
#'
#' This function takes the merged result data frame (e.g., from \code{\link{read_all_results}}),
#' filters for a specific condition, and computes group-wise means and standard deviations
#' of ARI and BIC for MFA and MPCA.
#'
#' @param df A \code{data.frame} containing columns like \code{K_true, r_true, N, M, sep, noise, ARI_MFA, ARI_PCA, BIC_MFA, BIC_PCA}.
#'
#' @return A \code{data.frame} summarizing ARI/BIC means and sds, grouped by \code{(K_true, r_true, N, M)}.
#'
#' @examples
#' \dontrun{
#' df_merged <- read_all_results("results")
#' df_table2 <- make_table2(df_merged)
#' knitr::kable(df_table2, caption="Table 2: Summary of ARI, BIC, etc.")
#' }
#'
#' @export
make_table2 <- function(df){
  df_subset <- df %>%
    dplyr::filter(abs(sep - 1.0)<1e-9, abs(noise-1.0)<1e-9)
  
  df_table2 <- df_subset %>%
    dplyr::group_by(K_true, r_true, N, M) %>%
    dplyr::summarise(
      ARI_MFA_mean = mean(ARI_MFA),
      ARI_MFA_sd   = sd(ARI_MFA),
      BIC_MFA_mean = mean(BIC_MFA),
      ARI_PCA_mean = mean(ARI_PCA),
      ARI_PCA_sd   = sd(ARI_PCA),
      BIC_PCA_mean = mean(BIC_PCA),
      .groups="drop"
    ) %>%
    dplyr::arrange(K_true, r_true, N, M)
  
  df_table2
}
