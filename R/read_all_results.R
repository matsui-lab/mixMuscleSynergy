#' Read All Results from .rds Files
#'
#' This function searches a directory for \code{.rds} files containing one-condition results
#' (each is typically a \code{data.frame} of 1 or more rows), reads them, and merges them
#' into a single data frame.
#'
#' @param result_dir A character string specifying the directory to look for \code{.rds} files.
#' @return A combined \code{data.frame} with rows from all the loaded files.
#'
#' @details
#' Each \code{.rds} file is assumed to contain a \code{data.frame} with columns such as
#' \code{K_true, r_true, N, M, sep, noise, ARI_MFA, BIC_MFA, ...}. This function simply does a
#' \code{rbind} of all those data frames. If no \code{.rds} files are found, it raises an error.
#'
#' @examples
#' \dontrun{
#' df_merged <- read_all_results("results/")
#' head(df_merged)
#' }
#'
#' @export
read_all_results <- function(result_dir = "results") {
  file_vec <- list.files(result_dir, pattern="\\.rds$", full.names=TRUE)
  
  if(length(file_vec)==0){
    stop("No .rds files found in ", result_dir)
  }
  
  df_list <- lapply(file_vec, function(f){
    x <- readRDS(f)
    x
  })
  df_merged <- do.call(rbind, df_list)
  df_merged
}
