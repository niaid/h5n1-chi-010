load(file.path(PROJECT_DIR, "DATA_PROCESSED/Microarrays/PBMC/samples.clean_genes.iqr", "gexp_d0_fc.RData"), verbose = T)

fc.th = 0.5
q.th = 0.2
tp = unique(info.fc$time.point)

fsamp = matrix(nrow=nrow(dat.fc), ncol=length(tp)*2)
rownames(fsamp) = rownames(dat.fc)
colnames(fsamp) = paste(rep(tp,1,each=2), c("pos","neg"), sep=".")
gi = rep(FALSE, nrow(dat.fc))
k = 1
for (t in tp) {
  ti = info.fc$time.point==t
  fsamp[,k] = apply(dat.fc[,ti], 1, function(x) sum(x > fc.th)/length(x))
  fsamp[,k+1] = apply(dat.fc[,ti], 1, function(x) sum(x < -fc.th)/length(x))
  gi = gi | fsamp[,k] >= q.th | fsamp[,k+1] >= q.th
  k = k + 2
}
sum(gi) # 6672 genes

df = as.data.frame(dat.fc[gi,]) %>%
  tibble::rownames_to_column("gene") %>%
  gather("sample.name","fc",-gene) %>%
  inner_join(info.fc, by="sample.name") %>%
  select(-sample.name, -plate.num, -iso.date, -n.d0) %>%
  # filter(time.point!="d000_00h") %>%
  group_by(gene,subject.id, time.point) %>%
  summarise(fc = mean(fc,na.rm=T)) %>%
  ungroup() %>%
  spread(time.point, fc)

df.mat.2 = as.matrix(df[,-(1:2)])
rownames(df.mat.2) = paste(df$subject.id, df$gene, sep="__")

dn.out = file.path(PROJECT_DIR, "RESULTS/Microarrays/PBMC/pattern_discovery/")
saveRDS(df, file.path(dn.out, glue::glue("df_{sum(gi)}g.rds")))
saveRDS(df.mat.2, file.path(dn.out, glue::glue("df.mat.2_{sum(gi)}g.rds")))

