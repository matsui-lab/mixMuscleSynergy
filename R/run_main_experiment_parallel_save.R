#' Main Example: Run Parallel MFA/MPCA Simulation and Save Each Condition
#'
#' This function calls \code{\link{simulate_and_run_parallel_save}} with some default arguments,
#' prints the merged results, and generates a simple figure comparing ARI_MFA and ARI_PCA.
#'
#' @return The merged data frame of results across all conditions.
#'
#' @examples
#' \dontrun{
#' df_out <- run_main_experiment_parallel_save()
#' head(df_out)
#' }
#' @export
run_main_experiment_parallel_save <- function(){
  df_result <- simulate_and_run_parallel_save(
    K_true_vec = c(3,4),
    r_true_vec = c(3,4),
    N_vec      = c(50),
    M_vec      = c(8,12),
    T_each     = 100,
    sep_vec    = c(0.5,1.0),
    noise_vec  = c(1.0,2.0),
    max_iter   = 20,
    n_rep_init = 3,
    seed_data_base = 999,
    mc.cores = 2,       # e.g. 2コア
    output_dir = "results"
  )
  
  # show
  print(head(df_result, 10))
  
  # quick plot
  df_long <- df_result %>%
    select(K_true, r_true, N, M, sep, noise, ARI_MFA, ARI_PCA) %>%
    gather(key="Method", value="ARI", ARI_MFA, ARI_PCA)
  
  p <- ggplot(df_long, aes(x=factor(sep), y=ARI, color=factor(noise),
                           group=interaction(Method, noise))) +
    geom_point(position=position_dodge(width=0.3)) +
    geom_line(position=position_dodge(width=0.3)) +
    facet_wrap(~Method + K_true + r_true + N + M, labeller=label_both) +
    theme_bw() +
    labs(title="MFA vs MPCA (parallel, each condition saved separately)",
         x="Separation (sep)", color="Noise", y="ARI")
  
  ggsave("figure_ari_mfa_mpca_parallel_save.pdf", p, width=8, height=6)
  
  df_result
}
