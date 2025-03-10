% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility.R
\name{compute_BIC_mfa}
\alias{compute_BIC_mfa}
\title{Compute BIC for an MFA model}
\usage{
compute_BIC_mfa(loglik, K, r, M, N)
}
\arguments{
\item{loglik}{The total log-likelihood from \code{compute_logLik_mfa}.}

\item{K}{Number of clusters.}

\item{r}{Factor dimension.}

\item{M}{Observed dimension (channels).}

\item{N}{Number of subjects.}
}
\value{
A numeric BIC value.
}
\description{
Given the log-likelihood from \code{\link{compute_logLik_mfa}}, and known
\code{K}, \code{r}, \code{M}, \code{N}, this function computes
\deqn{ \mathrm{BIC} = -2 \log \mathcal{L} + \nu \log(N), }
where \(\nu\) is a naive parameter count for the MFA model. A typical formula is
\eqn{\nu = K \times M \times (r + 2) + (K - 1).}
}
\details{
The parameter count may vary depending on how you treat diagonal Psi, means, etc.
A common approximation is \code{K * M * (r + 2) + (K - 1)}. Adapt if your parameterization differs.
}
