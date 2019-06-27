fn.ann1 = file.path(PROJECT_DIR, "DATA_ORIGINAL/Flow_10c/GateList_10c_2.txt")
fn.ann2 = file.path(PROJECT_DIR, "DATA_ORIGINAL/Flow_10c/Flow_10c_ann.txt")

ann1 = fread(fn.ann1) %>% 
  mutate(ID = paste0("ID", Id))
ann2 = fread(fn.ann2)

i1 = ann1$ID %in% ann2$ID
all(i1)
i2 = match(ann1$ID[i1], ann2$ID)

ann1$Name = NA
ann1$Name[i1] = ann2$Name[i2]

dir.create(file.path(PROJECT_DIR, "DATA_PROCESSED/Flow_10c"), showWarnings = F)
fn.out = file.path(PROJECT_DIR, "DATA_PROCESSED/Flow_10c/flow_10c_info.txt")
fwrite(ann1, file=fn.out, sep="\t")
