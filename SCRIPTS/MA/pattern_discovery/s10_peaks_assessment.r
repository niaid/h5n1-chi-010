# PURPOSE: To add data for subject 10 and rescore the correlation matrix.

# Initialize.
source("SCRIPTS/0_initialize.r")
library(Biobase)
source(file.path(PROJECT_DIR, "SCRIPTS/functions/get_score.r"))
library(ComplexHeatmap)
library(circlize)

# H5N1-010 data
# fn = "../data/filtered/eset.genes.filtered.RData"
fn = file.path(PROJECT_DIR, "DATA_PROCESSED/Microarrays/PBMC/samples.all_genes.iqr/", 
               "eset.genes.filtered.RData")
load(fn, verbose = T)

dat = exprs(eset.genes)
info = pData(eset.genes) %>% 
  mutate(time.point = factor(time.point, levels=unique(.$time.point)))

si = info$subject.id == "H5N1-010"
info = info[si,]
dat = dat[,si]
z = matrix(nrow=nrow(info),ncol=3)
colnames(z) = paste0("G0",1:3)
dn.patt = file.path(PROJECT_DIR, "RESULTS/Microarrays/PBMC/pattern_discovery")

for(k in 1:3) {
  # fn.sig = sprintf("GE_pattern_%02d_gene.sig.contr_ann.txt",k)
  fn.sig = file.path(dn.patt, "GE_pattern_genes", sprintf("GE_pattern_%02d_gene.sig.contr_ann.txt",k))
  pat.gene = fread(fn.sig, sep="\t", data.table = F)$gene
  gi = rownames(dat) %in% pat.gene
  sum(gi)
  
  tmp = dat[gi,]
  z[,k] = get_score(tmp)
}
infz = info %>% 
  cbind(data.frame(z)) %>% 
  gather("pattern","score",matches("^G0"))
# ggplot(infz %>% filter(!is.na(score)), aes(time.point, score, col=pattern, group=pattern)) + 
#   geom_line() +
#   # facet_wrap(~subject) + 
#   geom_vline(xintercept=match(c("d001","d022"),levels(infz$time.point)), lty=2) +
#   theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))
# ggsave(sprintf("s10_raw_profile_d0.subst_%s.png",today()), w=5,h=3)

# substitute d0 with d21 and remove 2h
i0 = info$time.point=="d000_00h"
i2 = info$time.point=="d000_02h"
i21 = info$time.point=="d021_00h"

dat[,i0] = dat[,i21,drop=F]
dat[,i2] = NA

# get scores
df.patt.mat = readRDS(file.path(dn.patt, "df.patt.rds"))
patt.genes.clean = readRDS(file.path(dn.patt, "patt.genes.clean.rds"))

K = length(patt.genes.clean)
long.idx = c(rep(1,3),2:ncol(df.patt.mat))
long.idx.z = c(rep(1,2),2:ncol(dat))

z = matrix(nrow=nrow(info),ncol=K)
colnames(z) = paste0("G0",1:K)
for(k in 1:K) {
  gi = rownames(dat) %in% patt.genes.clean[[k]]
  sum(gi)
  tmp = dat[gi,]
  z[,k] = get_score(tmp)
}

Xcor = cor(z[long.idx.z,], t(df.patt.mat[,long.idx]), use = "pairwise.complete.obs")
Xscore = diag(Xcor)
names(Xscore) = sprintf("G%02d",1:K)

saveRDS(Xscore, file.path(dn.patt, "s10_pattern_scores.rds"))

