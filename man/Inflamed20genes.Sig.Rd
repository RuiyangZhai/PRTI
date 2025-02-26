% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Signature.R
\name{Inflamed20genes.Sig}
\alias{Inflamed20genes.Sig}
\title{Function to compute the signature as published by Montironi et al. in 2023.}
\usage{
Inflamed20genes.Sig(expr, Response = NULL, verbose = FALSE, ...)
}
\arguments{
\item{expr}{A matrix or data.frame of expression values where rows correspond to genes and columns correspond to samples.}

\item{Response}{A vector containing "R" and "NR" that matches the samples, used to calculate AUC. Default to NULL, and AUC will not be calculated.}

\item{verbose}{Gives information about each calculation step}

\item{...}{Other arguments passed on to gsva()'s params argument.}
}
\description{
This function computes signature scores from gene expression values following the algorithm used for the signature as published by Montironi et al. in 2023.
}
\examples{
a = Inflamed20genes.Sig(expr = gene_table,Response = pdata_table$Response)

}
\references{
Montironi, Carla et al. “Inflamed and non-inflamed classes of HCC: a revised immunogenomic classification.” Gut vol. 72,1 (2023): 129-140. doi:10.1136/gutjnl-2021-325918
}
