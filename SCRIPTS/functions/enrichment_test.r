if(!require(tmod))
  stop("tmod package is not installed")

enrichment_test <- function(input, col.effect = NULL, col.p = NULL, col.gene = NULL, rank.by="effect", rank.type="unsigned", mset="LI", 
                            tmodFunc = tmodCERNOtest, qval.th = 999, order.by = "none", file.out = "") {
  
  # Change qval.th to a small value (f.e. 0.05) and order.by to pvalue for a single comparison.
  # Only significant modules will be in the output
  
  # Check compatibility between columns and rank parameters
  if(!is.null(col.effect)) {
    effect.sign = sign(input[,col.effect])
  } else {
    if(rank.by=="effect" | rank.type != "unsigned") {
      error("No effect column specified. Cannot rank by effect or use signed rank.")
    }
  }
  if(is.null(col.p) & rank.by == "pvalue") {
    error("No pvalue column specified. Cannot rank by pvalue.")
  }
  
  if(!is.null(col.gene)) {
    rownames(input) = input[,col.gene]
  }
  
  # ranking ordering depends on chosen variable
  ord <- switch( paste(rank.by, rank.type, sep="_"), 
                 # msd=order(msd[,cn], decreasing=TRUE),
                 effect_unsigned = order(abs(input[,col.effect]), decreasing=T),
                 effect_positive = order(input[,col.effect], decreasing=T),
                 effect_negative = order(input[,col.effect], decreasing=F),
                 pvalue_unsigned = order(input[,col.p], decreasing=F),
                 pvalue_positive = order(log10(input[,col.p])*effect.sign, decreasing=F),
                 pvalue_negative = order(log10(input[,col.p])*effect.sign, decreasing=T)
  )
  res = tmodFunc( rownames(input)[ord] , qval = qval.th, mset=mset, order.by = order.by)
  
  # output
  fe = tools::file_ext(file.out)
  if(nchar(fe)>0) {
    fn.out = sub(fe, paste("tmod", mset, col.effect, col.p, rank.type, fe, sep="."), file.out)
  } else {
    fn.out = paste(file.out, "tmod", mset, col.effect, col.p, rank.type, sep=".")
  }
  write.csv(res, fn.out, row.names = F)
}
