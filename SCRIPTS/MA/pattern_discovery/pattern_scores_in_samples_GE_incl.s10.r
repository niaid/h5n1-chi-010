# was pattern_scores_in_samples_GE_incl.s10_170420.r

dn.patt = file.path(PROJECT_DIR, "RESULTS/Microarrays/PBMC/pattern_discovery")

subj.patt.cor = readRDS(file.path(dn.patt, "subj.patt.cor.rds"))
K = ncol(subj.patt.cor)
colnames(subj.patt.cor) = sprintf("Gp%02d",1:K)

s10 = readRDS(file.path(dn.patt, "s10_pattern_scores.rds"))
subj.patt.cor = rbind(subj.patt.cor, s10)
rownames(subj.patt.cor)[nrow(subj.patt.cor)] = "H5N1-010"
subj.patt.cor = subj.patt.cor[order(rownames(subj.patt.cor)),]

saveRDS(subj.patt.cor, file.path(dn.patt, "subj.patt.cor.incl.s10.rds"))

subj.patt.cor = as.data.frame(subj.patt.cor) %>% 
  tibble::rownames_to_column("subject") %>% 
  mutate(subject = sub("_", "-", subject))

fwrite(subj.patt.cor, file=file.path(dn.patt, "GE_pattern_scores_42subj.txt"), sep="\t")

