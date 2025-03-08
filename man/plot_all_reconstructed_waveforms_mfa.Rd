% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{plot_all_reconstructed_waveforms_mfa}
\alias{plot_all_reconstructed_waveforms_mfa}
\title{Plot Reconstructed Waveforms for All Clusters and Channels (MFA)}
\usage{
plot_all_reconstructed_waveforms_mfa(
  mfa_fit,
  list_of_data,
  overlay_subjects = TRUE
)
}
\arguments{
\item{mfa_fit}{A fitted MFA model (with z, Lambda, Psi, mu, etc.).}

\item{list_of_data}{A list of (T_i x M) data matrices used in fitting.}

\item{overlay_subjects}{Logical; if TRUE, each subject in the same cluster is drawn
in the same facet. If FALSE, the function computes the average waveforms across
subjects in each cluster.}
}
\value{
A ggplot object with line plots of reconstructed waveforms, facet by cluster & channel.
}
\description{
This function plots the reconstructed waveforms (Xhat_i) for each subject
across all channels, grouped by cluster. By default, it overlays all subjects
in the same cluster, then uses facet to separate by (Cluster, Channel).
}
\examples{
\dontrun{
  p_recon <- plot_all_reconstructed_waveforms_mfa(mfa_fit, list_of_data)
  print(p_recon)
}
}
