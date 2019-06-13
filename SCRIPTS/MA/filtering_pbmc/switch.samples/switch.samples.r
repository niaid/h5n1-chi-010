switch.samples <- function(eset, sample1,sample2) {
  if(!all(sample1 %in% sample2)) stop("Samples mismatch between the lists!")
  if(any(sample1 == sample2)) stop("Samples order matches between the lists!")

  pd = pData(eset)
  edat = exprs(eset)
  
  pd$sample.name = colnames(edat)
  
  if(!all(sample1 %in% colnames(edat))) {
    stop("Some samples are missing")
  }

  si1 = match(sample1, colnames(edat))
  si2 = match(sample2, colnames(edat))

  tmp.pd = pd[si1,,drop=F]
  pd[si1,] = pd[si2,]
  pd[si2,] = tmp.pd

  tmp = edat[,si1,drop=F]
  edat[,si1] = edat[,si2]
  edat[,si2] = tmp
  
  colnames(edat) = pd$sample.name

  eset=new("ExpressionSet",exprs=edat)
  pData(eset)=pd
  # exprs(eset) = edat
  return(eset)
}
