fn = file.path(PROJECT_DIR, "DATA_ORIGINAL/Emory/GPL13158.annot.gz")

library(GEOquery)
# geo = getGEO("GPL13158", AnnotGPL=T, getGPL=T) # get annotation directly from GEO
geo = getGEO(filename=fn) # use predownloaded data for reproducibility
ann = geo@dataTable@table
dim(ann)

adx = match(rownames(mat), ann$ID)
# all.equal(rownames(mat), as.character(ann$ID)[adx])

out = ann %>% select(ID,`Gene symbol`) %>% 
  rename(gene=`Gene symbol`) %>% 
  mutate(ID=as.character(ID), gene=as.character(gene)) %>% 
  filter(!grepl("///",gene) & gene!="")

fn.map = file.path(PROJECT_DIR, "DATA_PROCESSED/Emory/GPL13158.ann.txt")
write.table(out, file=fn.map, sep="\t", quote=F, col.names=T, row.names=F)

source(file.path(PROJECT_DIR, "SCRIPTS/functions/pick.probeset.r"))
pick.probeset(eset, fn.map) # generates file.map.pc1

# fn.map = file.path(PROJECT_DIR, "DATA_PROCESSED/Emory/GPL13158.ann_PC1.txt")
# gene.map = read.table(file.map.pc1, sep="\t", header=T, row.names=NULL, stringsAsFactors = F)
