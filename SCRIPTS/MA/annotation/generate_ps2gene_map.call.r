library(affy)

if(is.character(input.table.file) && file.exists(input.table.file)) {
	input.table=read.table(file.path(cur.dir,input.table.file),header=T,row.names=1,as.is=T,sep="\t")
	eset=new("ExpressionSet",exprs=as.matrix(input.table)) 
} else if(is.character(eset.file) && file.exists(eset.file)) {
load(eset.file)
eset=get(eset.name)
} else {
	stop("No proper input file found.")
}

source(file.path(PROJECT_DIR, "SCRIPTS/functions/pick.probeset.r"))
pick.probeset(eset, gene.map.file)
