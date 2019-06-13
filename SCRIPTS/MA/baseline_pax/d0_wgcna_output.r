dn.in = file.path(PROJECT_DIR, "DATA_PROCESSED/Microarrays/PAXgene/baseline")
dat0.in = readRDS(file.path(dn.in, "dat0.in_isv.0.7_7901g.rds"))
info0.in = readRDS(file.path(dn.in, "info0.in.rds")) %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames("subject.id")

dn.in = file.path(PROJECT_DIR, "RESULTS/Microarrays/PAXgene/baseline/")
fn.in = file.path(dn.in, "dat0.in-networkConstruction-auto_power.6.RData")
load(fn.in, verbose = T)

mod.ngenes = table(net$colors)[-1]

pc1.ev = rep(NA, length(mod.ngenes))
pc1.rot = vector("list", length(mod.ngenes))
for(m in 1:max(net$colors)) {
  gi = net$colors==m
  X = dat0.in[,gi]
  pca=prcomp(X,scale=TRUE)
  pc1.ev[m] = pca$sdev[1]/sum(pca$sdev)*100
  pc1.rot[[m]] = data.frame(gene=colnames(dat0.in)[gi],pc1.contr=pca$rotation[,1]) %>% 
    tibble::remove_rownames() %>% arrange(-pc1.contr)
}

dn.out = file.path(PROJECT_DIR, "RESULTS/Microarrays/PAXgene/baseline/")
for(m in 1:max(net$colors)) {
  fn.gene = file.path(dn.out, sprintf("GE.pax_d0_WGCNA_GbWB%02d_genes.txt",m))
  fwrite(pc1.rot[[m]], file=fn.gene, sep="\t")
}


library(ComplexHeatmap)
library(circlize)

ME.ord = order(as.numeric(sub("ME","",names(net$MEs))))
ME.ord = ME.ord[-1]
ME = net$MEs[,ME.ord]
rownames(ME) = rownames(dat0.in)
colnames(ME) = sprintf("GbWB%02d", 1:length(mod.ngenes))

iqr = apply(ME, 2, IQR)

df.ME = ME %>% tibble::rownames_to_column("subject")


fwrite(df.ME, file.path(dn.out, "GE.pax_d0_WGCNA_ME_scores.txt"), sep="\t")

df.summary = data.frame(module=names(ME), ngenes=as.numeric(mod.ngenes), iqr, pc1.ev) %>% 
  tibble::remove_rownames()
fn.summary = file.path(dn.out, "GE.pax_wgcna_modules_summary.txt")
fwrite(format(df.summary), file=fn.summary, sep="\t")
