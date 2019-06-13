library(Biobase)

# switch samples
cur.dir = file.path(PROJECT_DIR, "SCRIPTS/MA/filtering_pbmc/switch.samples")
source(file.path(cur.dir, "switch.samples.r"))

fn.in = file.path(PROJECT_DIR, "DATA_PROCESSED/Microarrays/PBMC/eset.apt.RData")
load(fn.in, verbose=T)

fn = file.path(cur.dir, "samples.switch.txt")
df.switch = read.table(fn, sep="\t", header=F, row.names=NULL, stringsAsFactors=F)

eset = switch.samples(eset, df.switch[,1],df.switch[,2])

fn.out = file.path(PROJECT_DIR, "DATA_PROCESSED/Microarrays/PBMC/eset.corrected.RData")
save(eset, file = fn.out)

