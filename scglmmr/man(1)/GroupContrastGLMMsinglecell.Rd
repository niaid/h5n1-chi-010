% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GroupContrastGLMMsinglecell.r
\name{GroupContrastGLMMsinglecell}
\alias{GroupContrastGLMMsinglecell}
\title{GroupContrastGLMMsinglecell - within each cell type contrast difference in fold change between groups, baseline difference, and fold change across groups of module scores.}
\usage{
GroupContrastGLMMsinglecell(
  module_data_frame,
  celltype_column = "celltype",
  metadata,
  fixed_effects = NULL,
  lmer_formula = NULL,
  plotdatqc = TRUE,
  figpath
)
}
\arguments{
\item{module_data_frame}{data for each cell to model -- designed to be scores for modules (columns) returned by scglmmr::WeightedModuleScore}

\item{celltype_column}{the column in metadata with the cluster / celltype designation for each cell}
}
\value{

}
\description{
GroupContrastGLMMsinglecell - within each cell type contrast difference in fold change between groups, baseline difference, and fold change across groups of module scores.
}
\examples{
#### CALCULATE MIXED MODEL SCORE FOR 2 GROUP CONTRAST TIME DELTA AND BASELINE DELTA
# set up module data frame
`module_data_frame` are any data for ech cell to model -- designed to be scores for each model
 module_df = h1@meta.data \%>\%
  select(barcode_check, celltype_joint, module_n)
  # add metadata for lme4 model
  met = read_delim("git_ignore/full_metadata/full_sample_metadata.txt", delim = "\t")
 # format metadata as factors for lme4 ordered for contrast contrast_fit = contrast(emm1, method = list( (c21 - c20) - (c11 - c10) ))
}
