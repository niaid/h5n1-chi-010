# was calculate_d0_fc.160915.r

library("Biobase")
source(file.path(PROJECT_DIR, "SCRIPTS/functions/factor.date.r"))

dn = "DATA_PROCESSED/Microarrays/PBMC/samples.clean_genes.iqr"
fn = file.path(PROJECT_DIR, dn ,"eset.genes.filtered.RData")
load(fn, verbose = T)
# load("../data/filtered.2/eset.genes.filtered.switched.2.RData", verbose=T)
edat = exprs(eset.genes)
einfo=pData(eset.genes) %>%
  mutate(plate.num = factor(plate.num),
         sample.name=make.names(sample.name),
         iso.date = factor.date(iso.date))

# si = einfo$plate.num!=6
# edat = edat[,si]
# einfo = einfo[si,]

# subtract day0
cols = c("sample.name", "subject.id", "time.point", "plate.num", "iso.date")
info = select(einfo, one_of(cols))

info.fc = info %>%
  group_by(subject.id, plate.num) %>%
  mutate(n.d0 = sum(time.point=="d000_00h")) %>%
  ungroup() %>%
  filter(n.d0>0)

si = match(info.fc$sample.name, colnames(edat))
df = as.data.frame(t(edat[,si])) %>%
  cbind(select(info.fc,-n.d0), .)

df = df %>%
  gather("gene","value", -one_of(cols))

df = df %>%
  group_by(subject.id, plate.num, gene) %>%
  # mutate(n.d0 = sum(time.point=="d000_00h")) %>%
  mutate(fc = value - mean(value[time.point=="d000_00h"],na.rm=T)) %>%
  ungroup()

dat.fc = df %>%
  select(sample.name, gene, fc) %>%
  spread(sample.name, fc)
genes = dat.fc$gene
dat.fc = as.matrix(dat.fc[,-1])
rownames(dat.fc) = genes

fn.rda = file.path(PROJECT_DIR, dn, "gexp_d0_fc.RData")
save(dat.fc, info.fc, file=fn.rda)

# distribution histograms (slow)
# ggplot(filter(df,time.point!="d000_00h"), aes(x=fc, group=sample.name, col=as.factor(plate.num))) + 
#   geom_density() +
#   coord_cartesian(xlim=c(-1,1))
# ggsave("hist_fc_dist_col.by.sample.png")
