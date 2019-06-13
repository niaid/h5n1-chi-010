dn = file.path(PROJECT_DIR, "DATA_PROCESSED/Emory")
eset = readRDS(file.path(dn, "eset.rds"))

file.map.pc1 = file.path(PROJECT_DIR, "DATA_PROCESSED/Emory/GPL13158.ann_PC1.txt")
gene.map = read.table(file.map.pc1, sep="\t", header=T, row.names=NULL, stringsAsFactors = F)

gi = match(gene.map$ID, featureNames(eset))
eset.gene = eset[gi,]
featureNames(eset.gene) <- gene.map$gene

saveRDS(eset.gene, file.path(dn, "eset.gene.rds"))
