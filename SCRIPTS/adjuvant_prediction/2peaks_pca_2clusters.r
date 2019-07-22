# PURPOSE: Predicting Adjuvant signature using cellular and transcriptomic parameters. (C) Principal component analysis of patients using top selected patterns - Gp01, Gp02, Gp03, Fp01, Fp05. Color indicates two clusters of patients after k-mean clustering with k=2.

# was 2peaks_pca_171208.r

# plot subjects in PCA space (after scaling of individual parameters) of 
# G2, G3, F1 (and we can also test adding F5 and IP-10 from Luminex) to see 
# if we can finalize/tune how we call POS vs. NEG.

# Initialize
source("SCRIPTS/0_initialize.r")

# flow
fn.f = file.path(PROJECT_DIR, "RESULTS/Flow_10c/DeltaThr0.5_DonorThr0.3_PercentOfParent_CORpearson_cutreeHybrid_deepSplit0_minClusterSize31_DonorScores.txt")
subj.patt.flow = read.csv(fn.f, row.names=1, stringsAsFactors = F) %>% data.matrix()
colnames(subj.patt.flow) = sub("Module.","Fp0",colnames(subj.patt.flow))
# rownames(subj.patt.flow) = sub("s","H5N1_0",rownames(subj.patt.flow))

# expression
dn.g = file.path(PROJECT_DIR, "RESULTS/Microarrays/PBMC/pattern_discovery")
fn.g = file.path(dn.g, "subj.patt.cor.incl.s10.rds")
subj.patt.ge = readRDS(fn.g)
K = ncol(subj.patt.ge)
# colnames(subj.patt.ge) = sprintf("Gp%02d",1:K)
rownames(subj.patt.ge) = sub("H5N1-0","s",rownames(subj.patt.ge))

si.flow = rownames(subj.patt.flow) %in% rownames(subj.patt.ge)
si.ge = match(rownames(subj.patt.flow)[si.flow], rownames(subj.patt.ge))
identical(rownames(subj.patt.ge)[si.ge], rownames(subj.patt.flow)[si.flow])

dat.ge = subj.patt.ge[si.ge, colnames(subj.patt.ge) %in% c("Gp01", "Gp02","Gp03")]
dat.flow = subj.patt.flow[si.flow, colnames(subj.patt.flow) %in% c("Fp01","Fp05")]

dat = cbind(dat.ge, dat.flow)

dat.sc = scale(dat)

pca = prcomp(dat.sc, scale=F, center=F)
# plot(pca, type="l")
df.pc = cbind(as.data.frame(dat), as.data.frame(pca$x)) %>% 
  tibble::rownames_to_column("subject") %>% 
  mutate(label=sub("s0?","",subject))

set.seed(123)
df.pc$km2 = kmeans(dat.sc,2)$cluster
htree = hclust(dist(dat.sc))
df.pc$hc = cutree(htree, k = 2)

# PCA -------------------------------------
# icol = 1:5
# pca = prcomp(dat.sc[,icol], scale=T)
pca = prcomp(dat.sc, scale=T)
pcv = pca$sdev[1:2] / sum(pca$sdev) * 100
df.pc = df.pc %>%
  # mutate(pc1 = pca$x[,1], PC2 = pca$x[,2]) %>% 
  mutate(cluster = factor(km2))
ggplot(df.pc, aes(PC1,PC2, col=cluster)) +
  geom_point(size=2) +
  geom_text(aes(label=label), size=3, vjust=-0.1,hjust=-0.1, show.legend = F) +
  scale_color_manual(values = c("#6666FF","#FF6666")) +
  xlab(sprintf("PC1 (%.1f%%)", pcv[1])) +
  ylab(sprintf("PC2 (%.1f%%)", pcv[2])) +
  ggtitle(sprintf("PCA: %s", paste(colnames(dat),collapse=", "))) +
  theme_bw()

fn.fig = file.path(PROJECT_DIR, "FIGURES", "2peak_patterns_PCA.png")
ggsave(fn.fig, w=6,h=4)

dn.pred = file.path(PROJECT_DIR, "RESULTS/Adjuvant_prediction")
if(dir.exists(dn.pred)){
        print("Output folder already exists.")
}else{
        print("Output folder doesn't exists. Creating it now.")
        dir.create(dn.pred, recursive = T, showWarnings = T)
}
fwrite(df.pc, file.path(dn.pred, "2peaks_patterns_PCA.txt"), sep="\t")
