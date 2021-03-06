% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pseudobulk_de.R
\name{dreamMixedModel}
\alias{dreamMixedModel}
\title{dreamMixedModel - run dream mixed model
note due to this issue must require variancepartition https://github.com/GabrielHoffman/variancePartition/issues/17}
\usage{
dreamMixedModel(
  dge_lists,
  apriori_contrasts = FALSE,
  sample_column,
  contrast_matrix = NULL,
  design_matrix,
  fixed_effects,
  cell_metadata,
  lme4_formula = "~ 0 + cohort_timepoint + (1|sampleid)",
  plotsavepath,
  ncores = 4,
  version = "2"
)
}
\arguments{
\item{dge_lists}{list of dgelists created with NormalizePseudobulk}

\item{apriori_contrasts}{one of TRUE or FALSE, whether to fit a priori contrasts}

\item{sample_column}{quoted character e.g. "sample" the subject level sample variable should have multiple timepoints subjectID_timepoint i.e. s1_0, s1_1}

\item{contrast_matrix}{contrast matrix created with make.contrasts}

\item{design_matrix}{design matrix created with BulkDesignMatrix}

\item{fixed_effects}{a vector of covariates that are columns in metadata e.g. a vector of  'age' 'gender'}

\item{cell_metadata}{metadata cells-rows variables-columns i.e. ColData or seuratmeta.data}

\item{lme4_formula}{symbolic model formula the default is '~ 0 + cohort_timepoint + (1|sampleid)'}

\item{plotsavepath}{a path to save created plot of contrasts}

\item{ncores}{number of cores for doParallel}

\item{version}{if using R 3.5 bioc < 3.8 make version '1' otherwse leave default 2 runs bioc 3.8 and 3.9  this argument is to maintain backwards compatibility with R 3.5 workflows with the VariancePartition package}
}
\value{
list of model fits indexed by celltype
}
\description{
dreamMixedModel - run dream mixed model
note due to this issue must require variancepartition https://github.com/GabrielHoffman/variancePartition/issues/17
}
\examples{
# run dream mixed model
}
