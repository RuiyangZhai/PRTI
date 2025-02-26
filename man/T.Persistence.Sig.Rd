% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Signature.R
\name{T.Persistence.Sig}
\alias{T.Persistence.Sig}
\title{Function to compute the signature as published by Mariathasan et al. in 2018.}
\usage{
T.Persistence.Sig(expr, Response = NULL, verbose = FALSE, ...)
}
\arguments{
\item{expr}{A matrix or data.frame of expression values where rows correspond to genes and columns correspond to samples.}

\item{Response}{A vector containing "R" and "NR" that matches the samples, used to calculate AUC. Default to NULL, and AUC will not be calculated.}

\item{verbose}{Gives information about each calculation step}

\item{...}{Other arguments passed on to gsva()'s params argument.}
}
\description{
This function computes signature scores from gene expression values following the algorithm used for the signature as published by Mariathasan et al. in 2018.
}
\examples{
a = T.Persistence.Sig(expr = gene_table,Response = pdata_table$Response)

}
\references{
Mariathasan, Sanjeev et al. “TGFβ attenuates tumour response to PD-L1 blockade by contributing to exclusion of T cells.” Nature vol. 554,7693 (2018): 544-548. doi:10.1038/nature25501
}
