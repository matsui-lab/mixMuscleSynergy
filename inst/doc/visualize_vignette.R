## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----function-doc, eval=FALSE-------------------------------------------------
#  #' Read All Results and Produce Tables/Figures
#  #'
#  #' This function demonstrates how to read all .rds files from a directory,
#  #' merge them into a single \code{data.frame}, generate a couple of tables (Table 1, Table 2),
#  #' and produce various figures (Figure 1, Figure 2, etc.) saving them to files.
#  #'
#  #' @param result_dir Directory containing .rds result files (default "results").
#  #'
#  #' @return Invisibly returns the merged data frame of all results.
#  #'
#  #' @details
#  #' The user can modify the conditions in \code{\link{make_table1}} and \code{\link{make_table2}}
#  #' to match their actual experiment ranges. Figures are generated via \code{ggsave}.
#  #'
#  #' @examples
#  #' \dontrun{
#  #' df_all <- read_and_visualize_all("results")
#  #' }
#  #'
#  #' @export
#  read_and_visualize_all <- function(result_dir="results"){
#    df <- read_all_results(result_dir=result_dir)
#  
#    df_table1 <- make_table1()
#    cat("=== Table 1: Simulation Conditions ===\n")
#    print(df_table1)
#    # knitr::kable(df_table1, caption="Table 1: Conditions")
#  
#    df_table2 <- make_table2(df)
#    cat("=== Table 2: Summary of ARI etc. ===\n")
#    print(df_table2)
#    # kable(df_table2, caption="Table 2: Summary")
#  
#    p1 <- plot_fig1_ARI(df)
#    ggsave("figure1_ARI_sep_noise.pdf", p1, width=8, height=6)
#  
#    p2 <- plot_fig2_BIC_heatmap(df)
#    ggsave("figure2_BIC_heatmap.pdf", p2, width=6, height=5)
#  
#    if(nrow(df) > 1){
#      p3 <- plot_fig3_init_boxplot(df)
#      ggsave("figure3_init_boxplot.pdf", p3, width=6, height=4)
#    }
#  
#    # (Figure 4) requires an MFA fit object. Not shown by default.
#    # see the example in plot_fig4_loadings_heatmap
#  
#    invisible(df)
#  }

