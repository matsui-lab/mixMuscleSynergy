% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/best_fit.R
\name{select_optimal_K_r_mfa}
\alias{select_optimal_K_r_mfa}
\title{Search for the Best (K, r) for MFA via Grid Search}
\usage{
select_optimal_K_r_mfa(
  list_of_data,
  Kvec = 1:5,
  rvec = 1:5,
  max_iter = 50,
  nIterFA = 5,
  tol = 0.001
)

select_optimal_K_r_mfa(
  list_of_data,
  Kvec = 1:5,
  rvec = 1:5,
  max_iter = 50,
  nIterFA = 5,
  tol = 0.001
)
}
\arguments{
\item{list_of_data}{A list of data matrices, each \code{(T_i x M)}.}

\item{Kvec}{An integer vector of candidate \code{K} values.}

\item{rvec}{An integer vector of candidate \code{r} values.}

\item{max_iter}{Maximum EM iterations for \code{mfa_em_fit}.}

\item{nIterFA}{Sub-iterations within the factor-analyzer step.}

\item{tol}{Convergence tolerance for \code{mfa_em_fit}.}
}
\value{
A list with:
\item{summary}{A data frame with columns \code{K, r, logLik, BIC}.}
\item{best_model_info}{A list with \code{(K, r, model, logLik, BIC)} for the best BIC.}
\item{all_models}{A list of all fits, each \code{(K, r, model, logLik, BIC)}.}

A list with:
\item{summary}{A data frame with columns \code{K, r, logLik, BIC}.}
\item{best_model_info}{A list with \code{(K, r, model, logLik, BIC)} for the best BIC.}
\item{all_models}{A list of all fits, each \code{(K, r, model, logLik, BIC)}.}
}
\description{
This function tries all combinations of \code{K} in \code{Kvec} and
\code{r} in \code{rvec}, fits a Mixture Factor Analysis (MFA) model
using \code{\link{mfa_em_fit}} (or your wrapper) for each combination,
computes BIC, and returns the best model by BIC.

This function tries all combinations of \code{K} in \code{Kvec} and
\code{r} in \code{rvec}, fits a Mixture Factor Analysis (MFA) model
using \code{\link{mfa_em_fit}} (or your wrapper) for each combination,
computes BIC, and returns the best model by BIC.
}
\details{
Internally calls \code{\link{mfa_em_fit}} for each \code{(K, r)}. Then we compute log-likelihood
via \code{\link{compute_logLik_mfa}} and BIC via \code{\link{compute_BIC_mfa}}.

Internally calls \code{\link{mfa_em_fit}} for each \code{(K, r)}. Then we compute log-likelihood
via \code{\link{compute_logLik_mfa}} and BIC via \code{\link{compute_BIC_mfa}}.
}
\examples{
\dontrun{
out <- select_optimal_K_r_mfa(my_data, Kvec=2:4, rvec=1:3, max_iter=50, nIterFA=10, tol=1e-4)
names(out)
# "summary", "best_model_info", "all_models"
out$summary
out$best_model_info
}

\dontrun{
out <- select_optimal_K_r_mfa(my_data, Kvec=2:4, rvec=1:3, max_iter=50, nIterFA=10, tol=1e-4)
names(out)
# "summary", "best_model_info", "all_models"
out$summary
out$best_model_info
}

}
