# PURPOSE: To determine patterns from the fold change data.

# Initialize the environment with PROJECT_DIR and some commonly used packages.
source("SCRIPTS/0_initialize.r")
library(gplots)
library(ComplexHeatmap)
library(circlize)
library(cluster)

# MAIN
# Load fold change data
dn = "DATA_PROCESSED/Microarrays/PBMC/samples.clean_genes.iqr"
fn = file.path(PROJECT_DIR, dn, "gexp_d0_fc.RData")
load(file=fn, verbose = T)

# source("plot_cluster_profile_lines.r")
info.fc$time.point = factor(info.fc$time.point)

# Select genes that have a fold change of greater than one in atleast 30% of
# samples at any given time point.
fc.th = 1 # threshold FC > 1
q.th = 0.3 # in 30% of samples at any time point
tp = unique(info.fc$time.point)
# tp = c("d001")#,"d022")
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
sum(gi) # 202 genes

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

df.mat = as.matrix(df[,-(1:2)])
rownames(df.mat) = paste(df$subject.id, df$gene, sep="__")

# Calculate statistics.
df.mean = apply(df.mat, 1, function(x) mean(abs(x),na.rm=T))
df.median = apply(df.mat, 1, function(x) median(x,na.rm=T))
df.sd = apply(df.mat, 1, function(x) sd(x,na.rm=T))
df.cv = df.sd/abs(df.mean) 
df.iqr = apply(df.mat, 1, function(x) IQR(x,na.rm=T))
df.max = apply(df.mat, 1, function(x) max(abs(x),na.rm=T))

# Calculate up and down regulated genes in the gene list as filtered above.
gene.up = rownames(dat.fc)[gi][rowMeans(dat.fc[gi,],na.rm=T)>0]
gene.dn = rownames(dat.fc)[gi][rowMeans(dat.fc[gi,],na.rm=T)<0]
df.mat.up = df.mat[df$gene %in% gene.up,]
df.mat.dn = df.mat[df$gene %in% gene.dn,]
df.mat.all = df.mat

# Carry out correlation.
long.idx = c(rep(1,3),2:ncol(df.mat.all))
df.mat.cc.all = cor(t(df.mat.all[,long.idx]), use = "pairwise.complete.obs")
# df.mat.cc.up = cor(t(df.mat.up), use = "pairwise.complete.obs")
# df.mat.cc.dn = cor(t(df.mat.dn), use = "pairwise.complete.obs")

df.mat = df.mat.all
df.mat.cc = df.mat.cc.all
fc.dir = "all"

# Save the data matrix and the correlation data matrix.
## Check if the output directory exists.
dn.out = file.path(PROJECT_DIR, "RESULTS/Microarrays/PBMC/pattern_discovery")
if(dir.exists(dn.out)){
        print("Output directory exists.")
}else {
        print("Output directory doesn't exists. Creating it now.")
        dir.create(dn.out, recursive = T, showWarnings = T)
}
saveRDS(df.mat, file.path(dn.out, "df.mat.rds"))
saveRDS(df.mat.cc, file.path(dn.out, "df.mat.cc.rds"))

# Carry out clustering through DIANA.
# set.seed(123)
ptm <- proc.time()
dia = diana(1-df.mat.cc, diss=T, trace.lev = 1)
# dia = diana(1-abs(df.mat.cc), diss=T, trace.lev = 1)
# dia = diana(sqrt(2*(1-df.mat.cc)), diss=T, trace.lev = 1)
proc.time() - ptm
# as on 180222 about 26 min:
# user  system elapsed 
# 1602.29    2.47 1612.55 
# as on 180227 about 35 min:
# user  system elapsed 
# 2064.17    1.28 2070.93 
# as on 180227_2 about 23 min:
# user  system elapsed 
# 1397.47    2.29 1403.11 

# Save clustering results from DIANA to a R datastructure file. 
saveRDS(dia, file.path(dn.out, "diana.object.abs.rds")) # as on 161109
# saveRDS(dia, "diana.object.abs.161004.rds")
# dia = readRDS("diana.object.160927.rds")

# Plot the dendrogram from DIANA.
dia2 = dia
dia2$order.lab = rep(NA,nrow(df.mat.cc))
plot(dia2, which.plots = 2, hang=0)
dev.copy(png, file.path(dn.out, "diana_dendrogram.png"), width=8000, height=1100) # as on161110
dev.off()

