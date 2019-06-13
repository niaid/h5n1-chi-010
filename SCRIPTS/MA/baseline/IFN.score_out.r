library(Biobase)
source(file.path(PROJECT_DIR, "SCRIPTS/functions/factor.date.r"))

# data
# fn = "data/all_genes_170426/eset.genes.filtered.RData"
fn = file.path(PROJECT_DIR, "DATA_PROCESSED/Microarrays/PBMC/samples.clean_genes.all/eset.genes.filtered.RData")
load(fn, verbose = T)

dat = exprs(eset.genes)
info = pData(eset.genes) %>% 
  mutate(time.point = factor(time.point, levels=unique(.$time.point)))%>% 
  mutate(subject.id = sub("_","-",subject.id))

si = info$time.point %in% c("d000_00h","d001","d007","d021_00h","d022","d028") & info$subject.id != "H5N1_010"
info = info[si,]
dat = dat[,si]

# IFN signature from GbX11 module
fn.pax = "c:/Users/kotliary/Box/R/h5n1/PAX/wgcna/mod_data/GE.pax_d0_WGCNA_M11_genes.txt"
df.pax = read.table(fn.pax, sep="\t", header=T, row.names=NULL, stringsAsFactors = F)

# test.genes = intersect(df.pbmc$gene, mod.genes); gset="GB13+BTM"
# test.genes = mod.genes; gset="BTM"
# test.genes = df.pbmc$gene; gset="GB13"
test.genes = df.pax$gene; gset="GBX11"

gi = rownames(dat) %in% test.genes
sum(gi)

# read clinical data
df.clin = read.table("data/clinical_info_adj.txt", sep="\t",
                     header=T, row.names=NULL, stringsAsFactors=F) %>% 
  select(subject=Subject.ID, Adjuvant) 

# dat.sc = t(dat[gi,]) %>% scale() %>% t()

df.test = info %>% 
  select(subject.id, time.point) %>% 
  mutate(IFN.score = WGCNA::moduleEigengenes(t(dat[gi,]), rep(1,sum(gi)))$eigengenes %>% pull(1)) %>%
  group_by(subject.id, time.point) %>%
  summarize(IFN.score=mean(IFN.score)) %>%
  ungroup() %>% 
  left_join(df.clin, by=c("subject.id"="subject"))

df.test = df.test %>% 
  mutate(time.point = gsub("(?<!1)0+(\\d)","\\1",time.point, perl=T)) %>%
  mutate(time.point = sub("_0h","", time.point)) %>% 
  mutate(time.point = sub("_","+",time.point) %>% factor.date()) %>% 
  rename(subject=subject.id)

fwrite(df.test, "H5N1_IFN.scores_selected.dates.txt", sep="\t")

df.test.change = df.test %>% 
  spread(time.point, IFN.score) %>% 
  mutate_if(is.numeric, `-`, .$d0) %>% 
  select(-d0) %>% 
  gather(time.point, IFN.score, -subject, -Adjuvant, na.rm = T)
fwrite(df.test.change, "H5N1_IFN.scores_selected.dates-d0.txt", sep="\t")
