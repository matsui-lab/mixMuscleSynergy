% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation_parameters.R
\name{read_all_results}
\alias{read_all_results}
\title{Read All Results from .rds Files}
\usage{
read_all_results(result_dir = "results")

read_all_results(result_dir = "results")
}
\arguments{
\item{result_dir}{A character string specifying the directory to look for \code{.rds} files.}
}
\value{
A combined \code{data.frame} with rows from all loaded files.

A combined \code{data.frame} with rows from all the loaded files.
}
\description{
This function searches a directory for \code{.rds} files, each presumably containing
one or more rows of simulation results, then merges them into a single data frame.
}
\details{
If no files are found, the function raises an error. Otherwise, each file is read with
\code{\link{readRDS}} and appended via \code{\link{rbind}}.
}
