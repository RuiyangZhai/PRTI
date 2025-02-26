% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Signature.R
\name{Oncotype.DX.Sig}
\alias{Oncotype.DX.Sig}
\title{Function to compute the signature as published by Paik et al. in 2004.}
\usage{
Oncotype.DX.Sig(expr, Response = NULL, verbose = FALSE, ...)
}
\arguments{
\item{expr}{A matrix or data.frame of expression values where rows correspond to genes and columns correspond to samples.}

\item{Response}{A vector containing "R" and "NR" that matches the samples, used to calculate AUC. Default to NULL, and AUC will not be calculated.}

\item{verbose}{Gives information about each calculation step}

\item{...}{Other arguments passed on to oncotypedx()'s params argument.}
}
\description{
This function computes signature scores from gene expression values following the algorithm used for the signature as published by Paik et al. in 2004.
}
\examples{
a = Oncotype.DX.Sig(expr = gene_table,Response = pdata_table$Response)

}
\references{
Paik, Soonmyung et al. “A multigene assay to predict recurrence of tamoxifen-treated, node-negative breast cancer.” The New England journal of medicine vol. 351,27 (2004): 2817-26. doi:10.1056/NEJMoa041588
}
