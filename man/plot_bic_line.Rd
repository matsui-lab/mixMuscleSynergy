% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{plot_bic_line}
\alias{plot_bic_line}
\title{Line Plot of BIC vs. K or r}
\usage{
plot_bic_line(df_summary, x_axis = c("K", "r"))
}
\arguments{
\item{df_summary}{A data.frame with columns \code{K, r, BIC}.}

\item{x_axis}{One of \code{"K"} or \code{"r"}. That variable is placed on the x-axis.}
}
\value{
A \code{ggplot} object.
}
\description{
Takes the summary data.frame from \code{select_optimal_K_r_mfa()} (which has columns
\code{K, r, logLik, BIC}) and produces a simple line plot of BIC. You can choose whether
the x-axis is \code{K} or \code{r}, and color/line by the other.
}
