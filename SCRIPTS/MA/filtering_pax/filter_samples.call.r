library("affy")
#library("genefilter")

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

if(exists("rm.criteria")) {
	rm.names = names(rm.criteria)
	rm.idx = rep(FALSE,ncol(emat))
	for(i in 1:length(rm.names)) {
		rm.idx = rm.idx | pd[,rm.names[i]] %in% rm.criteria[[rm.names[i]]]
	}
}

emat.out=emat[,!rm.idx]
pd.out = pd[!rm.idx,]
eset.out=new("ExpressionSet",exprs=emat.out)
pData(eset.out)=pd.out

dir.create(output.dir, showWarnings = FALSE)
save(eset.out, file=file.path(output.dir, "eset.samples.filtered.RData"))

# write filtered data
write.table(emat.out,file.path(output.dir,"./gexp.samples.filtered.txt"),sep="\t",quote=F,col.names=NA)

# write removed samples
if(!'sample.name' %in% rm.names)
	rm.names = c('sample.name',rm.names)
rm.out = pd[rm.idx,rm.names]
# colnames(rm.out) = c("sample.name","rm.samples","rm.subjects","rm.column")
write.table(rm.out,file.path(output.dir,"./samples.removed.txt"),sep="\t",quote=F,col.names=T,row.names=F)
# write filtered sample info
write.table(pd.out, file = file.path(output.dir,"sample_info.filtered.txt"), 
	quote = F, sep = "\t", row.names = F, col.names = T)
