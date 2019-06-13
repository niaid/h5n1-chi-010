

library(affy)

# gene expression phenotype table
gexp.pheno.table=read.table(phenotable.file,header=TRUE,as.is=TRUE,sep="\t")

# affy rma summary table to matrix
# row names are probesets and colnames are microarrays

# read raw rma summary table
apt.output.file=file.path(input.dir,"rma-sketch.summary.txt")

# index for comments to be removed
remove.ind=max(grep("#",readLines(apt.output.file,n=3000)))

apt.output=read.table(apt.output.file,header=T,row.names=1)

# convert chip names to R format
file.name.r=make.names(gexp.pheno.table[,'file.name'])

# keep chips in the phenotype table only
apt.output2=apt.output[,which(colnames(apt.output)%in%file.name.r)]

pidx = match(colnames(apt.output),file.name.r)
pidx = pidx[!is.na(pidx)]
apt.phenotable=gexp.pheno.table[pidx,]

# rename chip names
if(!is.na(chip.name.column)){
colnames(apt.output2)=make.names(apt.phenotable[,chip.name.column])
}

eset=new("ExpressionSet",exprs=as.matrix(apt.output2))

pData(eset)=as.data.frame(apt.phenotable)

###################################################################
# save apt table 
write.table(apt.output2,file.path(input.dir,"apt.summary.txt"),sep="\t",quote=F,col.names=NA)

# save eset
save(list=c("eset"),file=file.path(input.dir,"eset.apt.RData"))


