% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{plot_all_factor_scores_mpca}
\alias{plot_all_factor_scores_mpca}
\title{Plot All Factor Scores (Z) for MPCA, Faceted by Cluster and Factor}
\usage{
plot_all_factor_scores_mpca(mpca_fit, list_of_data, overlay_subjects = TRUE)
}
\arguments{
\item{mpca_fit}{A fitted MPCA model (with z, W, mu, etc.) plus
\code{factor_scores[[i]]} if available.}

\item{list_of_data}{A list of (T_i x M) data matrices (just for dimension reference).}

\item{overlay_subjects}{Logical; if TRUE, all subjects in the same cluster
are overlaid in the same facet. If FALSE, we can average them etc.}
}
\value{
A ggplot object that shows the factor scores over time,
  with facets by (Cluster, Factor).
}
\description{
This function uses precomputed factor scores \code{mpca_fit$factor_scores[[i]]}
(each \code{T_i x r}) and plots them over time, grouped by cluster \code{z[i]}.
Optionally overlays all subjects in the same cluster on one facet, or can average them.
}
\examples{
\dontrun{
  mpca_fit <- compute_factor_scores_mpca(mpca_fit, list_of_data)
  p <- plot_all_factor_scores_mpca(mpca_fit, list_of_data, overlay_subjects=TRUE)
  print(p)
}
}
