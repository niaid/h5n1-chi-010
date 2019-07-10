source("SCRIPTS/0_initialize.r")
source(file.path(PROJECT_DIR, "SCRIPTS/functions/get_score.r"))

dn = file.path(PROJECT_DIR, "DATA_PROCESSED/Emory")
eset.gene = readRDS(file.path(dn, "eset.gene.rds"))

dat = exprs(eset.gene)
info = pData(eset.gene)
iord = info$day %>% unique %>% sub("Day","",.) %>% as.numeric() %>% order()
info = info %>% 
  mutate(day=factor(day, levels=info$day[iord]))

z = matrix(nrow=nrow(info),ncol=3)
colnames(z) = paste0("Gp0",1:3)
dn.patt = file.path(PROJECT_DIR, "RESULTS/Microarrays/PBMC/pattern_discovery")

for(k in 1:3) {
  fn.sig = file.path(dn.patt, "GE_pattern_genes", sprintf("GE_pattern_%02d_gene.sig.contr_ann.txt",k))
  pat.gene = fread(fn.sig, sep="\t", data.table = F)$gene
  gi = rownames(dat) %in% pat.gene
  sum(gi)
  
  tmp = dat[gi,]
  z[,k] = get_score(tmp)
}
infz = info %>% 
  cbind(data.frame(z)) %>% 
  gather("pattern","score",matches("^Gp0"))


X.fc = infz %>% 
  # filter(pattern=="G03") %>% 
  group_by(subject, pattern) %>% 
  summarise(fc = max(score[day %in% c("Day1","Day22")],na.rm=T) -
              median(score[!day %in% c("Day1","Day22")],na.rm=T),
            fc.mean = mean(c(
              score[day=="Day1"]-score[day=="Day0"],
              score[day=="Day22"]-score[day=="Day21"]),na.rm=T),
            fc.max = max(c(
              score[day=="Day1"]-score[day=="Day0"],
              score[day=="Day22"]-score[day=="Day21"]),na.rm=T)
  ) %>% 
  ungroup() %>% 
  dplyr::select(-fc, -fc.max) %>% 
  spread(pattern, fc.mean)

saveRDS(X.fc, file.path(dn, "emory_GE_fc.rds"))

fwrite(X.fc, file.path(dn, "emory_2peaks_scores.txt"), sep="\t")
