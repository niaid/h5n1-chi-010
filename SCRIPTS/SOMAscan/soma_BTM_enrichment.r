source("SCRIPTS/0_initialize.r")
source(file.path(PROJECT_DIR, "SCRIPTS/functions/enrichment_test.r"))

# import data
sample.type = c("pbmc", "pax")

for (stype in sample.type) {
dn = file.path(PROJECT_DIR, "RESULTS/SOMAscan")
fn = glue::glue("{stype}_ADJ_RSPO3_bySex.txt")
input = read.csv(file.path(dn, fn), sep="\t", row.names = 1)

dn.out = file.path(PROJECT_DIR, "RESULTS/SOMAscan/BTM")
dir.create(dn.out, recursive = T,  showWarnings = F)
fn.out = file.path(dn.out, fn)

for (mset in c("LI","DC")) {
  for (col.effect in c("r","p_r")) {
    col.p = sub("r", "pVal", col.effect)
    for (rank.type in c("unsigned","positive","negative")) {
      enrichment_test(input, col.effect=col.effect, col.p=col.p, rank.by="pvalue", rank.type=rank.type, 
                      mset=mset, qval.th=999, file.out =fn.out) # large qval.th to output all modules
    }
  }
}
}
