

## R 3.6 -- 
#' GetContrastResultsRaw - calculate p values and return contrast results from modelfit with dreamMixedModel
#'
#' @param limma.fit.object.list the results returned by dreamMixedModel, to get coefficiennt from RunVoomLimma use GetContrastResults
#' @param coefficient.number corresponds to the contrast, the nmber is in order of the contrast matrix
#' @param contrast.name this does not have to match the exact contrast name used in the contrast matrix, it is here to force / remind the user to choose the correct contrast.
#' @return
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr mutate full_join select
#' @importFrom stats p.adjust
#' @importFrom plyr mapvalues
#' @export
#'
#' @examples
GetContrastResultsRaw =  function(limma.fit.object.list, coefficient.number, contrast.name){
  print("this function returns results from dreamMixedModel, to get coefficient from RunVoomLimma, use GetContrastResults
        GetContrastResults uses emperican Bayes shrinkage see https://github.com/GabrielHoffman/variancePartition/issues/4 ")
  ## pt 1 return ONLY gene, logFC, AveExpr, contrast, celltype from eBayes call to format data and add raw p values from lme4, correct with BH.
  pvals = lapply(limma.fit.object.list, function(x){ data.frame(P.Value = x$p.value[ ,coefficient.number], gene = rownames(x))})
  lapply(pvals, function(x) rownames(x) = NULL)

  # parameters to run ordinary t statistic
  coef = lapply(limma.fit.object.list, function(x){ x$coefficients[ ,coefficient.number] })
  stdev = lapply(limma.fit.object.list, function(x){ x$stdev.unscaled[ ,coefficient.number] })
  sigma_ = lapply(limma.fit.object.list, function(x){ x$sigma })
  # Note eBayes t statistics are NOT used, only for formatting output, next section adds raw p values and logFC returned by lme4
  ebf = lapply(limma.fit.object.list, function(x) suppressWarnings(eBayes(x))) ##
  contrast_result = GetContrastResults(limma.fit.object.list = ebf,
                                       coefficient.number = coefficient.number,
                                       contrast.name = contrast.name)
  contrast_result = lapply(contrast_result, function(x){ x %>% dplyr::select(gene, logFC, AveExpr, contrast, celltype) })
  result = t_statistic = list()
  for (i in 1:length(limma.fit.object.list)) {
    # compute ordinary t-statistic (see e.g. limma eBayes documentation: ordinary.t <- fit$coef[,2] / fit$stdev.unscaled[,2] / fit$sigma
    t_statistic[[i]] = data.frame(t = coef[[i]] / stdev[[i]] / sigma_[[i]]) %>% tibble::rownames_to_column("gene")
    result[[i]] =
      contrast_result[[i]] %>%
      dplyr::mutate(P.Value = plyr::mapvalues(x = gene, from = pvals[[i]]$gene, to = round(pvals[[i]]$P.Value, digits = 9))) %>%
      dplyr::mutate(adj.P.Val = stats::p.adjust(p = P.Value, method = "BH"))
    # add t statistic
    result[[i]] = dplyr::full_join(result[[i]], t_statistic[[i]], by = "gene")

  }
  names(result) = names(limma.fit.object.list)
  return(result)
}


