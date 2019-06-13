# filter genes by variation
eset.file=file.path(PROJECT_DIR, "DATA_PROCESSED/Microarrays/PBMC/samples.all","eset.samples.filtered.RData")
eset.name="eset.out"
use.eset = TRUE 

var.filter = TRUE
var.func="IQR"
var.cutoff=0.15
filterByQuantile=TRUE
require.symbol=TRUE # keep probeset with gene symbol only
affy.anno.file=file.path(PROJECT_DIR, "DATA_PROCESSED/Microarrays/annotation", "affy_hugene-2_1-st_unique_PC1.txt")
output.dir=file.path(PROJECT_DIR, "DATA_PROCESSED/Microarrays/PBMC/samples.all_genes.iqr")
