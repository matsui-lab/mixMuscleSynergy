% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation_comparison.R
\name{calc_total_SST}
\alias{calc_total_SST}
\title{Compute the Total Sum of Squares (SST) for a List of Data Matrices}
\usage{
calc_total_SST(list_of_data)
}
\arguments{
\item{list_of_data}{A list of \code{(T_i x M)} matrices.}
}
\value{
A numeric scalar, the total sum of squares (SST).
}
\description{
This function concatenates all matrices in \code{list_of_data} into one big matrix,
calculates the global mean, and then sums the squared differences from that mean.
The result can be used to compute VAF (Variance Accounted For).
}
