% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Signature.R
\name{HTICS.Sig}
\alias{HTICS.Sig}
\title{Function to compute the signature as published by Liu et al. in 2012.}
\usage{
HTICS.Sig(expr, Response = NULL, verbose = FALSE)
}
\arguments{
\item{expr}{A matrix or data.frame of expression values where rows correspond to genes and columns correspond to samples.}

\item{Response}{A vector containing "R" and "NR" that matches the samples, used to calculate AUC. Default to NULL, and AUC will not be calculated.}

\item{verbose}{Gives information about each calculation step}
}
\description{
This function computes signature scores from gene expression values following the algorithm used for the signature as published by Liu et al. in 2012.
}
\examples{
a = HTICS.Sig(expr = gene_table,Response = pdata_table$Response)

}
\references{
Liu, Jeff C et al. “Seventeen-gene signature from enriched Her2/Neu mammary tumor-initiating cells predicts clinical outcome for human HER2+:ERα- breast cancer.” Proceedings of the National Academy of Sciences of the United States of America vol. 109,15 (2012): 5832-7. doi:10.1073/pnas.1201105109
}
