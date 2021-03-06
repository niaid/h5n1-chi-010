% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/downstream_analysis_functions.r
\name{LeadEdgeTidySampleExprs}
\alias{LeadEdgeTidySampleExprs}
\title{LeadEdgeTidySampleExprs - convert a PseudobulkList into a tidy dataframe for each sample across cell types of the leading edge genes from a gsea list}
\usage{
LeadEdgeTidySampleExprs(av.exprs.list, gsea.list, padj.filter, NES.filter)
}
\arguments{
\item{av.exprs.list}{object returned by PseudobulkList summed or average counts}

\item{gsea.list}{object returned by RunFgseaOnRankList}

\item{padj.filter}{filter for adjusted p from GSEA}

\item{NES.filter}{filter for normalized enrichment score from GSEA}
}
\value{

}
\description{
LeadEdgeTidySampleExprs - convert a PseudobulkList into a tidy dataframe for each sample across cell types of the leading edge genes from a gsea list
}
\examples{
pos_enr = GetLeadingEdgeAverage(av.exprs.list = pb_average, gsea.list = gsea1, padj.filter = 0.1, NES.filter = 0)
}
