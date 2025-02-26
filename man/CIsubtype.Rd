% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Predict.R
\name{CIsubtype}
\alias{CIsubtype}
\title{Predicting chemoimmunotherapy response subtypes}
\usage{
CIsubtype(
  expr,
  threshold = 0.2,
  nPerm = 1000,
  verbose = FALSE,
  doPlot = TRUE,
  ...
)
}
\arguments{
\item{expr}{A matrix or data.frame of expression values where rows correspond to genes and columns correspond to samples.}

\item{threshold}{FDR thresholds used to determine CI-M subtypes.}

\item{nPerm}{an integer, number of permutations for pvalue estimation.}

\item{verbose}{Gives information about each calculation step.}

\item{doPlot}{logical, whether to produce prediction subHeatmap}

\item{...}{Other arguments passed on to ematAdjust()'s or ntp()'s params argument.}
}
\description{
This function classifies samples into subtypes responding to chemoimmunotherapy based on the gene expression matrix.
}
\examples{
# a = CIsubtype(expr = gene_table,threshold = 0.2)

}
