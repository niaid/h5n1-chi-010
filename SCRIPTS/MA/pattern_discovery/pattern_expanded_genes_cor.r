dn.patt = file.path(PROJECT_DIR, "RESULTS/Microarrays/PBMC/pattern_discovery")
df = readRDS(file.path(dn.patt, "df_6672g.rds"))
df.mat.2 = readRDS(file.path(dn.patt, "df.mat.2_6672g.rds"))

GE.patterns = fread(file.path(dn.patt, "GE_patterns_filt.txt"))
ct = GE.patterns$pattern %>% sub("Gp","",.) %>% as.numeric()
names(ct) = GE.patterns$label
ct.lev = sort(unique(GE.patterns$pattern))
K=max(ct)

df.patt.mat = readRDS(file.path(dn.patt, "df.patt.rds"))
# patt.genes.clean = readRDS(sprintf("patt.genes.clean_%s_170202.rds",cl.methods))

# K = nrow(df.patt.mat)
rownames(df.patt.mat) = sprintf("Gp%02d",1:K)

# LONG CALCULATIONS!!! Use saved data if possible
long.idx = c(rep(1,3),2:ncol(df.mat.2))
df.patt.cor = cor(t(df.mat.2[,long.idx]), t(df.patt.mat[,long.idx]), use="pairwise.complete.obs")
# df.patt.p = corr.test(t(df.mat.2[,long.idx]), t(df.patt.mat[,long.idx]), use="pairwise")$p

df.patt.cor.p = df.patt.cor
for (j in 1:nrow(df.patt.mat)) {
  cat(j,"")
  for (i in 1:nrow(df.mat.2)) {
    df.patt.cor.p[i,j] = cor.test(df.mat.2[i,long.idx], df.patt.mat[j,long.idx], alternative = "greater", exact=F)$p.value
  }
}

n_genes = length(unique(df$gene))
saveRDS(df.patt.cor, file.path(dn.patt, sprintf("df.cor_%dg.rds",n_genes)))
saveRDS(df.patt.cor.p, file.path(dn.patt, sprintf("df.cor.p_%dg.rds",n_genes)))

