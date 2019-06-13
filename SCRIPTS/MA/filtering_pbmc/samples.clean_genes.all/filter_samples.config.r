# filter sample (remove outliers)
eset.file = file.path(PROJECT_DIR, "DATA_PROCESSED/Microarrays/PBMC", "eset.corrected.RData")
eset.name="eset"
use.eset = TRUE

# all criteria in one list (test)
rm.criteria = list()
rm.criteria[["subject.id"]] = c("SeraCare_001")
rm.criteria[["plate.num"]] = 8
rm.criteria[["sample.name"]] = c(
"H5N1.003_d000_00h_7B04", # fail QC
"H5N1.009_d021_00h_7H12", # sex mismaatch
"H5N1.009_d021_02h_7A11", # sex mismaatch
"H5N1.010_d000_00h_7G09", # sex mismaatch
"H5N1.010_d000_00h_7H11", # sex mismaatch
"H5N1.010_d000_02h_7A10", # removed whole subject due to missing baseline
"H5N1.010_d000_04h_7B10", # removed whole subject due to missing baseline
"H5N1.010_d000_12h_7C10", # removed whole subject due to missing baseline
"H5N1.010_d001_7D10",     # removed whole subject due to missing baseline
"H5N1.010_d007_7E10",     # removed whole subject due to missing baseline
"H5N1.010_d021_00h_7F10", # removed whole subject due to missing baseline
"H5N1.010_d021_02h_7G10", # removed whole subject due to missing baseline
"H5N1.010_d021_04h_7H10", # removed whole subject due to missing baseline
"H5N1.010_d021_12h_7A09", # removed whole subject due to missing baseline
"H5N1.010_d022_7B09",     # removed whole subject due to missing baseline
"H5N1.010_d028_7C09",     # removed whole subject due to missing baseline
"H5N1.010_d042_7D09",     # removed whole subject due to missing baseline
"H5N1.010_d100_7E09"      # removed whole subject due to missing baseline
)

output.dir=file.path(PROJECT_DIR, "DATA_PROCESSED/Microarrays/PBMC/samples.clean")
