# filter sample (remove outliers)
eset.file = file.path(PROJECT_DIR, "DATA_PROCESSED/Microarrays/PBMC", "eset.corrected.RData")
eset.name="eset"
use.eset = TRUE

# all criteria in one list (test)
rm.criteria = list()
rm.criteria[["subject.id"]] = c("SeraCare_001")
rm.criteria[["plate.num"]] = 8

output.dir=file.path(PROJECT_DIR, "DATA_PROCESSED/Microarrays/PBMC/samples.all")
