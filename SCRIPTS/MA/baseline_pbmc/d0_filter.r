library("Biobase")
source(file.path(PROJECT_DIR, "SCRIPTS/functions/factor.date.r"))

dn = "DATA_PROCESSED/Microarrays/PBMC/samples.clean_genes.iqr"
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

df0 = as.data.frame(t(dat0)) %>%
  cbind(info0, .)
df0 = df0 %>%
  gather("gene","value", -one_of(cols))


# select high-variable genes by ISV ------------------------------

# LONG RUN
df.stab = data.frame(gene=factor(rownames(dat0), levels=rownames(dat0)), ISV=NA)
for (i in 1:nrow(dat0)) {
  # form = as.formula(sprintf("%s ~ subject", df.stab$gene[i]))
  form = as.formula("value ~ subject.id")
  fit = aov(form, data = df0 %>% dplyr::filter(gene==df.stab$gene[i]))
  ss = summary(fit)[[1]]["Sum Sq"][[1]]
  # names(ss) = c("ISV","WSV")
  ssn = ss / sum(ss)
  df.stab$ISV[i] = ssn[1]
}

dn.out = file.path(PROJECT_DIR, "RESULTS/Microarrays/PBMC/baseline")
dir.create(dn.out, showWarnings = F)

saveRDS(df.stab, file.path(dn.out, "df.stab_d0.rds"))
fwrite(df.stab, file.path(dn.out, "GE_pbmc_stable_genes_d0.txt"), sep="\t")

# identical(rownames(dat0), as.character(df.stab$gene))

isv.seq = seq(0.1,1,by=0.05)
isv.ngenes = rep(NA,length(isv.seq))
for(i in seq_along(isv.seq)) {
  isv.ngenes[i] = sum(df.stab$ISV > isv.seq[i])
}

dn.fig = file.path(PROJECT_DIR, "FIGURES/Baseline_PBMC")
dir.create(dn.fig, showWarnings = F)

plot(isv.seq, isv.ngenes, xlab="ISV threshold", ylab="No. of genes", pch=16, cex=1)
dev.copy(png, file.path(dn.fig, "d0_ngenes_vs_isv.seq.png"), width=500, height=400)
dev.off()

m0 = apply(dat0, 1, mean) # mean intensity

plot(m0, df.stab$ISV, xlab = "Mean intensity", ylab="ISV", pch=16, cex=.1)
abline(h=c(0.7,0.75,0.8), col="red")
dev.copy(png, file.path(dn.fig, "d0_ISV_vs_mean.png"), width=500, height=400)
dev.off()

isv.th = 0.7
gene.in = as.character(df.stab$gene[df.stab$ISV > isv.th])

df0.in = df0 %>%
  filter(gene %in% gene.in) %>% 
  group_by(subject.id, plate.num, iso.date, gene) %>% 
  summarise(value=mean(value)) %>% 
  ungroup()

dat0.in = df0.in %>% 
  select(subject.id, gene, value) %>%
  spread("gene","value") %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames("subject.id") %>% 
  as.matrix()

info0.in = info0 %>% 
  select(subject.id, plate.num, iso.date) %>% 
  distinct()
identical(rownames(dat0.in), as.character(info0.in$subject.id))

dn.out2 = file.path(PROJECT_DIR, "DATA_PROCESSED/Microarrays/PBMC/baseline")
dir.create(dn.out2, showWarnings = F)

saveRDS(dat0.in, file.path(dn.out2, sprintf("dat0.in_isv.%g_%dg.rds", isv.th, ncol(dat0.in))))
saveRDS(info0.in, file.path(dn.out2, "info0.in.rds"))
