% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pseudobulk_de.R
\name{NormalizePseudobulk}
\alias{NormalizePseudobulk}
\title{NormalizePseudobulk}
\usage{
NormalizePseudobulk(
  pseudobulklist,
  normalization.method = "RLE",
  design_matrix,
  minimum.gene.count = 1
)
}
\arguments{
\item{pseudobulklist}{object created with PseudobulkList only use this function if argument to PseudobulkList avg_or_sum was 'sum' computes normalization for pseudobulk libraries}

\item{normalization.method}{see edgeR function calcNormFactors}

\item{design_matrix}{the design matrix created with BulkDesignMatrix}

\item{minimum.gene.count}{see edgeR function filterbyExpr}
}
\value{
a list of dgeList indexed by celltype
}
\description{
NormalizePseudobulk
}
\examples{
# normalize the pseudobulk data made with MakePseudobulkList
}
