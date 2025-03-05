load("/media/owner/emg/Alessandro2018/data/fullset/dataset/participants_data.RData")
load("/media/owner/emg/Alessandro2018/data/fullset/dataset/EMG/Filtered/FILT_EMG.RData")
load("~/proj_mixture_synergy/data/PREPROC_EMG_bySubject.RData")

list_of_data <- lapply(PREPROC_EMG,function(x)apply(x,c(1,2),mean))
str(list_of_data)
# Candidate values for K and r:
Kvec <- 1:7
rvec <- 1:6

set.seed(100)
#sub_list_of_data <- list_of_data[sample(length(list_of_data),20)]
sub_list_of_data <- list_of_data
# Perform the search
res_mfa <- select_optimal_K_r_mfa(
  list_of_data = sub_list_of_data,
  Kvec = Kvec,
  rvec = rvec,
  max_iter = 50,
  nIterFA  = 5,
  tol      = 1e-3,
  n_init   = 1,             # Try multiple inits if needed
  use_kmeans_init = TRUE,  # Or TRUE if you prefer a k-means-based init
  subject_rdim_for_kmeans = 2
)

saveRDS(res_mfa,file="./../../result/res_mfa.rds")

res_pca <- select_optimal_K_r_mpca(
  list_of_data = sub_list_of_data,
  Kvec = Kvec,
  rvec = rvec,
  max_iter = 50,
  nIterPCA  = 5,
  tol      = 1e-3,
  n_init   = 1,             # Try multiple inits if needed
  use_kmeans_init = TRUE,  # Or TRUE if you prefer a k-means-based init
  subject_rdim_for_kmeans = 2
)
saveRDS(res_pca,file="./../../result/res_mpca.rds")

# The actual fitted MFA model:
mfa_fit <- res_mfa$best_model_info$model
mfa_fit <- compute_factor_scores_mfa(mfa_fit, sub_list_of_data)

rotated_mfa <- rotate_mfa_model(
  mfa_fit = mfa_fit,
  rotation = "varimax",
  rotate_scores = TRUE
)
rotated_aligned_mfa <- align_factor_signs_mfa(rotated_mfa)

p_loadings <- plot_cluster_synergy_loadings_mfa(rotated_aligned_mfa , plot_type = "heatmap")
print(p_loadings)

p_factor <- plot_all_factor_scores_mfa(
  rotated_aligned_mfa ,
  sub_list_of_data,
  overlay_subjects = TRUE
)
print(p_factor)


overall_stats <- summarize_vaf_sse_mfa(mfa_fit, sub_list_of_data, by_cluster = FALSE)
overall_stats
# Named vector: SSE=..., VAF=...

cluster_stats <- summarize_vaf_sse_mfa(mfa_fit, sub_list_of_data, by_cluster = TRUE)
cluster_stats


p_bic_K <- plot_bic_line(res_mfa$summary, x_axis="K")
print(p_bic_K)

p_bic_r <- plot_bic_line(res_mfa$summary, x_axis="r")
print(p_bic_r)



#####

# The actual fitted pca model:
pca_fit <- res_pca$best_model_info$model
pca_fit <- compute_factor_scores_mpca(pca_fit, sub_list_of_data)

rotated_pca <- rotate_mpca_model(
  mpca_fit = pca_fit,
  rotation = "varimax",
  rotate_scores = TRUE
)
rotated_aligned_pca <- align_factor_signs_mpca(rotated_pca)


p_loadings <- plot_cluster_synergy_loadings_mpca(rotated_aligned_pca , plot_type = "heatmap")
print(p_loadings)

p_factor <- plot_all_factor_scores_mpca(
  rotated_aligned_pca ,
  sub_list_of_data,
  overlay_subjects = TRUE
)
print(p_factor)


overall_stats <- summarize_vaf_sse_mfa(mfa_fit, sub_list_of_data, by_cluster = FALSE)
overall_stats
# Named vector: SSE=..., VAF=...

cluster_stats <- summarize_vaf_sse_mfa(mfa_fit, sub_list_of_data, by_cluster = TRUE)
cluster_stats


p_bic_K <- plot_bic_line(res_pca$summary, x_axis="K")
print(p_bic_K)

p_bic_r <- plot_bic_line(res_pca$summary, x_axis="r")
print(p_bic_r)



