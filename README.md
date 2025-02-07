mixMuscleSynergy: Muscle Synergy Analysis Using Mixture Models

**Yusuke Matsui**
Nagoya University

1. Introduction

Muscle synergy analysis is a technique widely used in movement science, biomechanics, and rehabilitation research to decompose large-scale electromyography (EMG) time-series data into a small number of "basis vectors (synergies)" and their activation intensities (weights). However, traditional methods such as Non-Negative Matrix Factorization (NMF), Principal Component Analysis (PCA), and Factor Analysis (FA) assume a common set of synergy vectors across all subjects, making it difficult to account for potential inter-subject differences in motor control strategies.

This package introduces the Mixture Model framework, which integrates "subgroup estimation (cluster assignment)" and "synergy basis estimation" using the Expectation-Maximization (EM) algorithm. Specifically, it implements:

Mixture Factor Analysis (MFA)

Mixture PCA (Probabilistic PCA)

allowing a comprehensive analysis of muscle activation patterns across subjects.

2. Features

Automatic determination of the number of clusters ($K$) using BIC for model selection

Simultaneous optimization via the EM algorithm, integrating subgroup assignment and synergy basis estimation

Implementation of comparative methods, including single PCA/NMF/FA and individual estimation followed by clustering

Support for high-dimensional EMG data analysis, optimized for multi-channel data

3. Installation

# Requires devtools package
install.packages("devtools")

devtools::install_github("yourusername/mixMuscleSynergy")

4. Usage

Simulating Data

library(mixMuscleSynergy)

# Generate synthetic EMG data for testing
sim_data <- simulate_mixture_data(
  N=50, K=3, r=4, M=8, T_each=100, 
  cluster_sep=1.0, noise_scale=1.0, seed=123
)

# sim_data contains:
# - list_of_data: 50 matrices (each 100 x 8)
# - z_true: Integer vector (50) with cluster assignments

Applying the Model

Mixture Factor Analysis (MFA)

fit_mfa <- mfa_em_fit(
  list_of_data = sim_data$list_of_data,
  K = 3, r = 4,
  max_iter = 50, tol = 1e-4
)

summary(fit_mfa)

Mixture PCA (MPCA)

fit_pca <- mixture_pca_em_fit(
  list_of_data = sim_data$list_of_data,
  K = 3, r = 4,
  max_iter = 50, tol = 1e-4
)

Model Selection Using BIC

num_rows <- sum(sapply(sim_data$list_of_data, nrow))
bic_mfa <- compute_bic(fit_mfa$logLik, K=3, r=4, M=8, N_total_rows=num_rows)
cat("BIC for MFA =", bic_mfa, "\n")

Running Parallel Simulations

df_result <- simulate_and_run_parallel_save(
  K_true_vec = c(3,4), r_true_vec = c(3,5),
  N_vec = c(50), M_vec = c(8,12), T_each = 100,
  sep_vec = c(0.5,1.0), noise_vec = c(1.0,2.0),
  max_iter = 20, n_rep_init = 3, seed_data_base = 999,
  mc.cores = 2, output_dir = "results"
)

Reading and Visualizing Results

# Load and visualize all results from saved .rds files
df_all <- read_and_visualize_all("results")

Extracting a Specific Model

my_model <- get_model_by_K_r(sel_mfa, K_target=3, r_target=2)

Visualizing Clustering Results

plot_clusters(fit_mfa)

5. References

d'Avella, A., & Bizzi, E. (2005). Shared and specific muscle synergies in natural motor behaviors. Proceedings of the National Academy of Sciences, 102(8), 3076-3081.

Tresch, M. C., Cheung, V. C., & d'Avella, A. (2006). Matrix factorization algorithms for the identification of muscle synergies: evaluation on simulated and experimental data sets. Journal of Neurophysiology, 95(4), 2199-2212.

Roh, J., Rymer, W. Z., Perreault, E. J., Yoo, S. B., Zhou, P., & Beer, R. F. (2013). Alterations in upper limb muscle synergy structure in chronic stroke survivors. Journal of Neurophysiology, 109(3), 768-781.

Safavynia, S. A., & Ting, L. H. (2012). Task-level feedback can explain temporal recruitment of spatially fixed muscle synergies throughout postural perturbations. Journal of Neurophysiology, 107(6), 159-177.

6. License

This project is licensed under the MIT License.

This package aims to enhance muscle synergy analysis by incorporating subgroups with different control strategies, enabling more accurate and meaningful interpretations of muscle activation data.
