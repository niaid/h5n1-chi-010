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
  filter(n.d0>0, time.point!="d042", !(subject.id=="H5N1-006" & time.point=="d000_00h" & plate.num!=12) )

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

si = match(info.fc$sample.name, colnames(dat.fc))
dat.fc = dat.fc[,si]
identical(colnames(dat.fc), info.fc$sample.name)

fn.rda = file.path(PROJECT_DIR, dn, "gexp_d0_fc.RData")
save(dat.fc, info.fc, file=fn.rda)
