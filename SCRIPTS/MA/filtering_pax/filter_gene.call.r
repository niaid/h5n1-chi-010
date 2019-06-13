library("affy")
library("genefilter")

if(use.eset) {
  load(eset.file)
  eset=get(eset.name)
  
  emat=exprs(eset)
  pd = pData(eset)
} else {
  emat = as.matrix(read.table(gexp.file, header=T, as.is=T, row.names=1, sep="\t"))
  pd = read.table(sample.info.file, header=T, as.is=FALSE, sep="\t")
} 
if(ncol(emat)!=nrow(pd) | all(make.names(pd[,'sample.name'])!=colnames(emat))) stop ("Samples in expression matrix and pdata do not match.")

if(exists("var.filter") & var.filter) {
  eset.probes=varFilter(eset, var.func=var.func, var.cutoff=var.cutoff, filterByQuantile=filterByQuantile)
} else {
  eset.probes = eset
}
gexp.output=exprs(eset.probes)

dir.create(output.dir, showWarnings = FALSE)
save(eset.probes, file=file.path(output.dir, "eset.probes.filtered.RData")) 

write.table(gexp.output,file.path(output.dir,"./gexp.probeid.txt"),sep="\t",quote=F,col.names=NA)

#####################
# filtered by genes with gene symbol only
if(require.symbol==TRUE){
affy.id.table=read.csv(affy.anno.file,as.is=T,sep="\t")

ai = affy.id.table$probeset %in% rownames(gexp.output)
gi = match(affy.id.table$probeset[ai], rownames(gexp.output))

gexp.symbol = gexp.output[gi,]
rownames(gexp.symbol) = affy.id.table$hgnc_symbol[ai]

##################
# write gexp table
write.table(gexp.symbol,file.path(output.dir,"./gexp.symbol.txt"),col.names=NA,sep="\t",quote=F)

eset.genes = new("ExpressionSet",exprs=gexp.symbol)
pData(eset.genes) = pData(eset.probes)
save(eset.genes, file=file.path(output.dir, "eset.genes.filtered.RData")) 
}
