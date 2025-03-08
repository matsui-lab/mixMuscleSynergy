% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mfa_em_fit.R
\name{mfa_em_fit}
\alias{mfa_em_fit}
\title{Fit a Mixture Factor Analysis Model (with optional multi-initialization)}
\usage{
mfa_em_fit(
  list_of_data,
  K,
  r,
  max_iter = 50,
  nIterFA = 20,
  tol = 0.001,
  n_init = 1,
  use_kmeans_init = TRUE,
  subject_rdim_for_kmeans = r,
  mc_cores = 1
)
}
\arguments{
\item{list_of_data}{A list of matrices (each \code{(T_i x M)}) to be modeled.}

\item{K}{Number of clusters.}

\item{r}{Factor dimension.}

\item{max_iter}{Maximum EM iterations for \code{mfaTimeseriesCpp}.}

\item{nIterFA}{Number of sub-iterations for the factor analyzer update in C++.}

\item{tol}{Convergence tolerance for \code{mfaTimeseriesCpp}.}

\item{n_init}{Number of random initial assignments to try (in addition to the default
single-run or k-means if requested). Defaults to \code{1}.}

\item{use_kmeans_init}{If \code{TRUE}, we also run one initialization where we
extract subject features (via PCA) and apply \code{kmeans} to get a cluster assignment.
Defaults to \code{FALSE}.}

\item{subject_rdim_for_kmeans}{PCA dimension for extracting subject-level features,
used only if \code{use_kmeans_init=TRUE}. Default is \code{r}.}

\item{mc_cores}{Number of cores for parallel execution via \code{mclapply}.
Defaults to \code{1} (no parallel).}
}
\value{
A list with the same elements as a single run: \code{z, pi, mu, Lambda, Psi, logLik, resp}.
  It corresponds to the best solution (i.e. highest log-likelihood) among all tried inits.
}
\description{
This function calls the C++ function \code{mfaTimeseriesCpp()} to perform a Mixture
Factor Analysis EM algorithm for fixed \code{K} and \code{r}. By default, it runs
a single pass with the internal initialization from C++. However, if you specify
multiple initial attempts (via \code{n_init>1}) and/or \code{use_kmeans_init=TRUE},
this function will try several different initial cluster assignments (in parallel
using \code{mclapply}), then return the best solution (maximizing the final log-likelihood).
}
\details{
The single-run approach in C++ already does a built-in initialization if \code{z_init}
is not provided. Therefore, if \code{n_init=1} and \code{use_kmeans_init=FALSE}, we just call
\code{mfaTimeseriesCpp} once (like the original design).
Otherwise:
\enumerate{
  \item If \code{use_kmeans_init=TRUE}, we do one run where we assign clusters by k-means
        on some subject-level features (extracted by PCA).
  \item We also sample \code{n_init} random assignments (each subject assigned to a random cluster).
  \item For each assignment, we call \code{\link{mfa_em_fit_cpp_singleInit}}.
  \item We compare all results' \code{logLik} and keep the best one.
}

This provides a more robust initialization scheme to avoid poor local maxima.
}
\examples{
\dontrun{
# Suppose we have a list_of_data for N=100 subjects, each (T_i x M=8).
# We want to fit K=3, r=2, and try 5 random inits plus one k-means init:
fit <- mfa_em_fit(
  list_of_data, K=3, r=2,
  max_iter=50, nIterFA=20, tol=1e-3,
  n_init=5,
  use_kmeans_init=TRUE,
  subject_rdim_for_kmeans=2,
  mc_cores=2
)
print(fit$logLik)
head(fit$z)
}

}
