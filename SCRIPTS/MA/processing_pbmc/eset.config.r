# input path
input.dir=file.path(PROJECT_DIR, "DATA_PROCESSED/Microarrays/PBMC")

# microarray sample information
phenotable.file=file.path(PROJECT_DIR, "DATA_ORIGINAL/Microarrays/PBMC/sample_info.txt")

# chip names in the matrix table and eset defined in the above phenotype table, otherwise, use the file names as chip names
chip.name.column="sample.name" # NA
