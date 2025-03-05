# Mixture-Based Muscle Synergy Analysis (MFA & MPCA)

Yusuke Matsui (Nagoya Univerisity)

## 1. Background

Muscle synergy analysis is often applied to reduce high-dimensional electromyographic (EMG) signals into a smaller set of “synergies” or activation patterns. Two mixture-based methods for identifying **group-level** synergies in heterogeneous populations are:

1. **Mixture Factor Analysis (MFA)**: Assumes each cluster’s covariance is decomposed by a Factor Analysis model with diagonal noise.  
2. **Mixture PCA (MPCA)**: Assumes each cluster is modeled by a probabilistic PCA approach with \(\sigma^2 I\) noise.

Both approaches can handle heterogeneous data, grouping subjects (or trials) into different synergy “clusters,” each with its own synergy loadings.

Compared to a single global factor/PCA decomposition, mixture-based approaches can uncover fundamentally different synergy patterns in subpopulations (e.g. half the subjects might use synergy pattern A, the other half synergy pattern B)


## 2. Method Overview and Effective Use Cases

### 2.1 Mixture Factor Analysis (MFA)

1. **Cluster membership**  
   Each subject $`i`$ belongs to one of $`K`$clusters with probability $`\pi_k`$.  

2. **Factor model per cluster**  
   Within cluster $`k`$, each subject’s data $`X_i`$ (dimension $`\text{time} \times M`$) follows a low-dimensional factor model:
   
   ```math
     X_i(t) = \mu_k + \Lambda_k \, Z_i(t) + \varepsilon_{i,t},
   ```
   where $`\mu_k`$ is the mean vector, $`\Lambda_k`$ is $`(M \times r)`$ factor loadings (synergy matrix), and $`\varepsilon_{i,t}`$ ~ $`\mathcal{N}(0, \Psi_k)`$ with $`\Psi_k`$ diagonal.

### 2.2 Mixture PCA (MPCA)

Similarly, we can define a mixture of probabilistic PCA (PPCA) models:

1. **Cluster membership**  
   Each subject \(i\) belongs to one of \(K\) clusters with probability \(\pi_k\).

2. **PCA model per cluster**  
   Data \(X_i(t)\) is modeled with
   \[
     X_i(t) = \mu_k + W_k \, Z_i(t) + \varepsilon_{i,t},
   \]
   where \(W_k\) is \((M \times r)\). The noise is isotropic \(\sigma^2_k I\). Often, we factor \(W_k = P_k D_k\) for orthonormal \(P_k\) and diagonal \(D_k\).

**Use cases**:
- **Heterogeneous populations**  
  E.g. subsets of participants with distinct synergy structures.  
- **Large datasets**  
  Enough participants/trials to fit multiple clusters.  
- **Exploratory approach**  
  You suspect more than one synergy pattern set in the population.


## 3. Data (Input) Format

All the provided functions assume:

- **`list_of_data`**: A list of length \(N\), where each element is a \((T_i \times M)\) matrix:
  - \(N\): number of subjects/trials.
  - \(T_i\): time points for subject \(i\).
  - \(M\): number of EMG channels (muscles).

Example:

```r
N <- 10
list_of_data <- vector("list", N)
for (i in seq_len(N)) {
  T_i <- 100
  M   <- 8
  X_i <- matrix(rnorm(T_i * M), nrow=T_i, ncol=M)
  list_of_data[[i]] <- X_i
}
```


## 4. Single Mixture Model Estimation

### 4.1 Mixture Factor Analysis (MFA)

Use `mfa_em_fit()` to fit an MFA model for fixed \(K\) (clusters) and \(r\) (factors):

```r
library(myMFAtools)  # hypothetical package collecting your scripts

# Suppose we want K=2 clusters, r=3 factors
fit_mfa <- mfa_em_fit(
  list_of_data = list_of_data,
  K = 2,
  r = 3,
  max_iter = 50,
  nIterFA  = 5,
  tol      = 1e-3,
  n_init   = 1,
  use_kmeans_init = TRUE
)
```

This returns a fitted MFA object with elements like `$z` (cluster assignments), `$Lambda`, `$mu`, `$Psi`, etc.

### 4.2 Mixture PCA (MPCA)

Use `mixture_pca_em_fit()` similarly:

```r
fit_mpca <- mixture_pca_em_fit(
  list_of_data = list_of_data,
  K = 2,
  r = 3,
  max_iter = 50,
  nIterPCA = 20,   # sub-iterations for the PCA updates
  tol = 1e-3,
  method = "EM",   # or "closed_form"
  n_init = 1,
  use_kmeans_init = TRUE
)
```

This returns an MPCA object, typically with `$z`, `$W` or `$P`/`$D`, `$mu`, `$sigma2`, `$logLik`, etc.


## 5. Parameter Meanings and Typical Settings

1. **`K`** (clusters): e.g., 1 to 5 for an exploratory synergy analysis.  
2. **`r`** (synergy/factor dimension): typically 1–5 for muscle synergy contexts.  
3. **`max_iter`**: outer EM iteration limit, e.g. 50–100.  
4. **`nIterFA`** (MFA) / **`nIterPCA`** (MPCA): sub-iterations for each cluster’s parameter updates.  
5. **`tol`**: convergence tolerance, typically ~1e-3.  
6. **`n_init`**: number of random initializations. More helps avoid local minima.  
7. **`use_kmeans_init`**: if `TRUE`, tries a k-means-based assignment to seed EM.


## 6. Visualizing the Fitted Model

### 6.1 MFA: Synergy Loadings & Factor Scores

After `mfa_em_fit()`, you can visualize synergy loadings:

```r
# Heatmap of factor loadings
p_loadings <- plot_cluster_synergy_loadings_mfa(fit_mfa, plot_type="heatmap")
print(p_loadings)

# Factor scores (time evolution) per subject, faceted by cluster & factor
fit_mfa <- compute_factor_scores_mfa(fit_mfa, list_of_data)
p_scores <- plot_all_factor_scores_mfa(fit_mfa, list_of_data, overlay_subjects=TRUE)
print(p_scores)
```

### 6.2 MPCA: Loadings & Factor Scores

Similar calls for MPCA:

```r
# Heatmap of loadings W
p_loadings_mpca <- plot_cluster_synergy_loadings_mpca(fit_mpca, plot_type="heatmap")
print(p_loadings_mpca)

# Factor scores
fit_mpca <- compute_factor_scores_mpca(fit_mpca, list_of_data)
p_scores_mpca <- plot_all_factor_scores_mpca(fit_mpca, list_of_data, overlay_subjects=TRUE)
print(p_scores_mpca)
```


## 7. Optimal Parameter Search (K, r)

If you do not know the best `(K, r)` in advance, you can search over a grid:

### 7.1 MFA Grid Search

```r
search_res_mfa <- select_optimal_K_r_mfa(
  list_of_data,
  Kvec = 1:5,
  rvec = 1:5,
  max_iter = 50,
  nIterFA  = 5,
  tol = 1e-3,
  n_init = 1,
  use_kmeans_init = TRUE
)

search_res_mfa$summary         # data.frame with K, r, logLik, BIC
search_res_mfa$best_model_info # best (K, r)
```

### 7.2 MPCA Grid Search

```r
search_res_mpca <- select_optimal_K_r_mpca(
  list_of_data,
  Kvec = 1:5,
  rvec = 1:5,
  max_iter = 50,
  nIterPCA = 20,
  tol = 1e-3,
  method = "EM",
  n_init = 1,
  use_kmeans_init = TRUE
)

search_res_mpca$summary
search_res_mpca$best_model_info
```

Both `select_optimal_K_r_mfa()` and `select_optimal_K_r_mpca()` compute a BIC for each combination, returning the best by BIC.

## 8. Refining the Chosen Model

Once you identify `(K^*, r^*)`, you can:

1. Extract the best model:

   ```r
   best_fit_mfa <- search_res_mfa$best_model_info$model
   best_fit_mpca <- search_res_mpca$best_model_info$model
   ```

2. Optionally **re-run** with more initializations, stricter tolerances, or greater iteration limits:

   ```r
   final_mfa <- mfa_em_fit(
     list_of_data = list_of_data,
     K = K^*,
     r = r^*,
     max_iter = 100,
     nIterFA  = 10,
     n_init = 5,
     tol = 1e-4
   )
   ```

   ```r
   final_mpca <- mixture_pca_em_fit(
     list_of_data = list_of_data,
     K = K^*,
     r = r^*,
     max_iter = 100,
     nIterPCA = 30,
     tol = 1e-4,
     n_init = 5
   )
   ```

3. Re-visualize loadings, compute factor scores, etc.


## 9. Post-hoc Rotation and Factor Score Recalculation

### 9.1 MFA Rotation

After MFA fitting, you can apply **varimax** or other rotations to each cluster’s factor loadings:

```r
# 1) Rotate factor loadings, also rotate factor scores:
rotated_mfa <- rotate_mfa_model(final_mfa, rotation="varimax", rotate_scores=TRUE)

# 2) Optionally flip signs so that the largest absolute loading is positive:
rotated_mfa <- align_factor_signs_mfa(rotated_mfa, method="max_abs")

# 3) If you did rotate_scores=FALSE above, you can recalc scores now:
rotated_mfa <- compute_factor_scores_mfa(rotated_mfa, list_of_data)

# 4) Plot again:
plot_cluster_synergy_loadings_mfa(rotated_mfa, plot_type="heatmap")
plot_all_factor_scores_mfa(rotated_mfa, list_of_data)
```

### 9.2 MPCA Rotation

Likewise, for MPCA:

```r
# 1) Rotate loadings W (one cluster at a time), optionally factor scores:
rotated_mpca <- rotate_mpca_model(final_mpca, rotation="varimax", rotate_scores=TRUE)

# 2) Flip sign if desired:
rotated_mpca <- align_factor_signs_mpca(rotated_mpca, method="max_abs")

# 3) Recompute scores if needed:
rotated_mpca <- compute_factor_scores_mpca(rotated_mpca, list_of_data)

# 4) Plot the new loadings & scores:
plot_cluster_synergy_loadings_mpca(rotated_mpca, plot_type="heatmap")
plot_all_factor_scores_mpca(rotated_mpca, list_of_data)
```

**Purpose**: Factor rotation helps interpret synergy patterns by simplifying loadings (e.g. varimax tries to make columns “sparse”), while sign alignment ensures consistent sign conventions.


## 10. Example Pipeline (MFA or MPCA)

Below is a minimal example that ties everything together. We’ll illustrate an MFA pipeline; the MPCA pipeline is nearly identical, just changing the function names.

```r
### Step 1: Create or load data
N <- 10
M <- 8
list_of_data <- lapply(1:N, function(i) matrix(rnorm(100 * M), 100, M))

### Step 2: Search for best (K, r) via MFA
search_res <- select_optimal_K_r_mfa(
  list_of_data,
  Kvec=1:3, rvec=1:3,
  max_iter=30, nIterFA=5, tol=1e-3, n_init=2,
  use_kmeans_init=TRUE
)
search_res$summary
best_info <- search_res$best_model_info

### Step 3: Refine best
final_mfa <- mfa_em_fit(
  list_of_data,
  K=best_info$K,
  r=best_info$r,
  max_iter=50,
  nIterFA=5,
  n_init=5,
  tol=1e-4
)

### Step 4: Plot synergy loadings
plot_cluster_synergy_loadings_mfa(final_mfa, plot_type="heatmap")

### Step 5: Factor scores
final_mfa <- compute_factor_scores_mfa(final_mfa, list_of_data)
plot_all_factor_scores_mfa(final_mfa, list_of_data)

### Step 6: (Optional) Rotation
rot_mfa <- rotate_mfa_model(final_mfa, rotation="varimax", rotate_scores=TRUE)
rot_mfa <- align_factor_signs_mfa(rot_mfa)
rot_mfa <- compute_factor_scores_mfa(rot_mfa, list_of_data)  # if needed

plot_cluster_synergy_loadings_mfa(rot_mfa, plot_type="heatmap")
plot_all_factor_scores_mfa(rot_mfa, list_of_data)
```

**For MPCA**: just swap out the MFA-specific function calls (`mfa_em_fit` → `mixture_pca_em_fit`, etc.) with the MPCA equivalents.


## 11. Closing Remarks

- **Interpretation**: Always check whether synergy loadings make physiological sense. Mixed-subject synergy analysis can reveal distinct synergy “clusters” that might be lost in a single global decomposition.  
- **Computation**: Large `(K, r)` grids and big data may require parallelization or focusing on smaller parameter ranges.  
- **Rotation**: Varimax or similar methods can improve interpretability by making synergy loadings sparser or more distinct. Sign alignment ensures consistent direction for easier comparison across factors or clusters.  
- **Further customizations**: For advanced usage, you can add custom rotations (e.g. oblimin, promax) or constraints on the synergy shapes.

Using **MFA** or **MPCA** in a mixture setting, alongside **post-hoc rotation** and **factor score recalculation**, can substantially enhance your ability to discover and interpret distinct synergy patterns in heterogeneous EMG datasets.
