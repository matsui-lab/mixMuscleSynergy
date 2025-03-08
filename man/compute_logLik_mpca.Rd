% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility.R
\name{compute_logLik_mpca}
\alias{compute_logLik_mpca}
\title{Compute the Mixture-Model Log-Likelihood for Mixture PCA}
\usage{
compute_logLik_mpca(list_of_data, mpca_fit)
}
\arguments{
\item{list_of_data}{A list of length \code{N}, each \code{(T_i x M)}.}

\item{mpca_fit}{A fitted Mixture PCA model, containing:
\itemize{
  \item \code{P[[k]]} or \code{W[[k]]} for cluster \code{k},
  \item \code{mu[[k]]}, \code{sigma2[k]}, \code{pi[k]}.
}}
}
\value{
A numeric scalar, the total mixture log-likelihood.
}
\description{
Given a fitted Mixture PCA model with \code{$P[[k]]}, \code{$mu[[k]]},
\code{$sigma2[k]}, and \code{$pi[k]}, we compute
\deqn{
  \sum_{i=1}^{N} \log\Bigl(\sum_{k=1}^{K} \pi_k
    \prod_{t=1}^{T_i}\mathcal{N}\bigl(x_{i,t}\mid \mu_k,\Sigma_k\bigr)\Bigr)
}
where \(\Sigma_k = W_k W_k^T + \sigma2_k I\). A typical approach is to
reconstruct \(\Sigma_k\) from \(\mathrm{P}[[k]]\) if needed (or \(\mathrm{W}[[k]]\)).
}
\details{
We assume each cluster k has \(\Sigma_k = W_k W_k^T + \sigma2_k I\).
You must assemble that to call \code{mvtnorm::dmvnorm} row by row.
}
