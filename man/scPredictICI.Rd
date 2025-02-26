% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Predict.R
\name{scPredictICI}
\alias{scPredictICI}
\title{PredictICI function for single cell and spatial transcriptome data}
\usage{
scPredictICI(
  object,
  method = "AUCell",
  slot = "data",
  asssy = NULL,
  verbose = FALSE,
  nCores = 1,
  normAUC = TRUE,
  ...
)
}
\arguments{
\item{object}{An Seurat object.}

\item{method}{Scoring algorithm, including \code{GSVA}, \code{ssGSEA} and \code{AUCell}.}

\item{slot}{Specific assay data to get.}

\item{asssy}{Specific assay to get data from. Spatial transcriptome data may need to set to \code{Spatial}.}

\item{verbose}{Gives information about each calculation step.}

\item{nCores}{When the method is selected as' aucell ', this is one of the parameters of the AUCell_calcAUC function.}

\item{normAUC}{When the method is selected as' aucell ', this is one of the parameters of the AUCell_calcAUC function}

\item{...}{Other arguments passed on to gsva()'s or AUCell_buildRankings()'s params argument.}
}
\description{
Similar to the PredictICI function, but the input data is in Seurat format.
}
\examples{
# a = scPredictICI(object = Seurat_obj, method = "AUCell", slot="data", asssy="RNA")

}
