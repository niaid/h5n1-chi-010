# was d0_wgcna_170319.r

dn.in = file.path(PROJECT_DIR, "DATA_PROCESSED/Microarrays/PBMC/baseline")
dat0.in = readRDS(file.path(dn.in, "dat0.in_isv.0.7_8144g.rds"))

dn.mod = file.path(PROJECT_DIR, "RESULTS/Microarrays/PBMC/baseline/")
fn.mod = file.path(dn.mod, "dat0.in-networkConstruction-auto_power.6.RData")
load(fn.mod, verbose = T)
K = max(net$colors)

library(tmod)
data(tmod)
bg = colnames(dat0.in)
mod.n.in = sapply(tmod$MODULES2GENES, function(x) sum(x %in% colnames(dat0.in)))
mset = "LI"
tmod.ann = tmod$MODULES %>% mutate(N.sel=mod.n.in) %>% 
  filter(SourceID==mset) %>% select(ID, Title, N=B, N.sel)
tmod.pv=c()
res = vector("list",K)
for(m in 1:K) {
  fg = bg[net$colors==m]
  res[[m]] = tmodHGtest(fg=fg,bg=bg, mset=mset, qval = Inf, order.by = "none")
  tmod.pv = cbind(tmod.pv, res[[m]]$adj.P.Val)
}
names(res) = sprintf("Gb%02d",1:K)
res[sapply(res,is.null)] = NULL
colnames(tmod.pv) = names(res)

dn.out = file.path(PROJECT_DIR, "RESULTS/Microarrays/PBMC/baseline/")
fn.rds = file.path(dn.out, sprintf("GE.pbmc_d0_WGCNA_BTM.%s.HG.pv.rds",mset))
saveRDS(res, fn.rds)
fn.tmod = file.path(dn.out, sprintf("GE.pbmc_d0_WGCNA_BTM.%s.HG.pv.fdr.txt",mset))
fwrite(cbind(tmod.ann,tmod.pv), fn.tmod, sep="\t")

