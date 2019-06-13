library("Biobase")
source(file.path(PROJECT_DIR, "SCRIPTS/functions/factor.date.r"))

dn = "DATA_PROCESSED/Microarrays/PAXgene/samples.all_genes.iqr"
fn = file.path(PROJECT_DIR, dn ,"eset.genes.filtered.RData")
load(fn, verbose = T)

edat = exprs(eset.genes)
einfo=pData(eset.genes) %>%
  mutate(plate.num = factor(plate.num),
         sample.name=make.names(sample.name),
         iso.date = factor.date(iso.date))

cols = c("sample.name", "subject.id", "time.point", "plate.num", "iso.date")
info0 = select(einfo, one_of(cols)) %>%
  filter(time.point %in% c("d000_00h"))

dat0 = edat[,colnames(edat) %in% info0$sample.name]

# get stability from PBMC
fn.stab = file.path(PROJECT_DIR, "RESULTS/Microarrays/PBMC/baseline/df.stab_d0.rds")
df.stab = readRDS(fn.stab)

# check genes variability by IQR ----------------------------------
dn.fig = file.path(PROJECT_DIR, "FIGURES/Baseline_PAXgene/")
dir.create(dn.fig, showWarnings = F)

iqr0 = apply(dat0,1,IQR)
m0 = apply(dat0,1,mean)
iqr.th = seq(0.1,1,by=0.05)
iqr.ngenes = rep(NA,length(iqr.th))
for(i in seq_along(iqr.th)) {
  iqr.ngenes[i] = sum(iqr0 > iqr.th[i])
}

png(file.path(dn.fig, "d0_ngenes_vs_iqr.th.png"), width=500, height=400)
plot(iqr.th, iqr.ngenes, xlab="IQR threshold", ylab="No. of genes", pch=16, cex=1)
dev.off()

isv.th = 0.7
gi.stab = rownames(dat0) %in% df.stab$gene[df.stab$ISV > isv.th]
gi.all = rownames(dat0) %in% df.stab$gene
gi.isv = match(rownames(dat0)[gi.all], df.stab$gene)

png(file.path(dn.fig, "d0_pbmc.isv_vs_pax.mean.png"), width=500, height=400)
plot(m0[gi.all], df.stab$ISV[gi.isv],  pch=16, cex=.1, xlab="PAX mean", ylab="ISV from PBMC")
abline(h=0.7, col="red", lty=2)
dev.off()


g.in = gi.stab
sum(g.in)

dat0.in = dat0[g.in,]
dim(dat0.in)


df0 = as.data.frame(t(dat0)) %>%
  cbind(info0, .)
df0 = df0 %>%
  gather("gene","value", -one_of(cols))


gene.in = rownames(dat0)[gi.stab]
length(gene.in)

df0.in = df0 %>%
  filter(gene %in% gene.in) %>% 
  group_by(subject.id, iso.date, gene) %>% 
  summarise(value=mean(value)) %>% 
  ungroup()

dat0.in = df0.in %>% 
  select(subject.id, gene, value) %>%
  spread("gene","value") %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames("subject.id") %>% 
  as.matrix()

info0.in = info0 %>% 
  select(subject.id, iso.date) %>% 
  distinct() %>% 
  arrange(subject.id)
identical(rownames(dat0.in), as.character(info0.in$subject.id))

dn.out2 = file.path(PROJECT_DIR, "DATA_PROCESSED/Microarrays/PAXgene/baseline")
dir.create(dn.out2, showWarnings = F)

saveRDS(dat0.in, file.path(dn.out2, sprintf("dat0.in_isv.%g_%dg.rds", isv.th, ncol(dat0.in))))
saveRDS(info0.in, file.path(dn.out2, "info0.in.rds"))
