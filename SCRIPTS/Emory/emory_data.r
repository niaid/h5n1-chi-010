dn.in = file.path(PROJECT_DIR, "DATA_ORIGINAL/Emory/")
fn.in = file.path(dn.in, "T H - Vax010_RMA_CHI.txt.zip")
# dat = readr::read_tsv(fn)
mat = read.csv(unz(fn.in,"Vax010_RMA_CHI.txt"), sep="\t", header = T, row.names = 1) %>% 
  data.matrix()
colnames(mat) = sub(".CEL", "", colnames(mat))

info = data.frame(sample.name=colnames(mat), stringsAsFactors = F) %>% 
  separate(sample.name, c("subject","day"),sep="\\.", remove=F, extra="merge")

fn.demo = file.path(dn.in, "Vax010_demographics_wAge.txt")
demo = read.table(fn.demo, sep="\t", header=T, row.names=NULL, stringsAsFactors = F) %>% 
  dplyr::rename(subject=SubjectID, Age=Age..years.)

info = left_join(info, demo, by="subject")
rownames(info) = info$sample.name

library(Biobase)
pd = new("AnnotatedDataFrame", data=info)
eset = ExpressionSet(mat, phenoData = pd)

dn.out = file.path(PROJECT_DIR, "DATA_PROCESSED/Emory")
dir.create(dn.out, showWarnings = F)
saveRDS(eset, file.path(dn.out, "eset.rds"))

