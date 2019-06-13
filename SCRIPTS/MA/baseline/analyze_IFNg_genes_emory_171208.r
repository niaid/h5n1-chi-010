library(ComplexHeatmap)
library(circlize)
source("functions/get_score.r")
source("functions/color_functions.r")

# data
eset.genes = readRDS("data/emory/eset.gene.rds")


dat = exprs(eset.genes)
info = pData(eset.genes) %>% 
  mutate(day = factor(day, levels=unique(.$day)))

si = info$day %in% c("Day0","Day1") #& info$subject.id != "H5N1_010"
info = info[si,]
dat = dat[,si]

# library(tmod)
# data(tmod)
# mod.id = c("LI.M75","LI.M127","LI.M150","LI.M67","LI.M165")
# tmod$MODULES[mod.id, 1:3]
# mod.gene = tmod$MODULES2GENES[mod.id]
# mod.df = data.frame()
# for(m in mod.id) {
#   mod.df = rbind(mod.df, data.frame(ID=m, gene=mod.gene[[m]]))
# }
# mod.genes = unique(mod.df$gene)

# fn.pbmc = "c:/Users/kotliary/Box/R/h5n1/wgcna/GE.pbmc_d0_WGCNA_M13_genes.txt"
# df.pbmc = read.table(fn.pbmc, sep="\t", header=T, row.names=NULL, stringsAsFactors = F)
fn.pax = "c:/Users/kotliary/Box/R/h5n1/PAX/wgcna/mod_data/GE.pax_d0_WGCNA_M11_genes.txt"
df.pax = read.table(fn.pax, sep="\t", header=T, row.names=NULL, stringsAsFactors = F)

# test.genes = intersect(df.pax$gene, mod.genes); gset="intersect"
# test.genes = mod.genes; gset="btm"
test.genes = df.pax$gene; gset="GBX11"

gi = rownames(dat) %in% test.genes
sum(gi)

# read predicted
df.pred = read.table("emory/emory_adjuvanted_subjects_predicted.txt", sep="\t", 
                     header=T, row.names=NULL, stringsAsFactors=F) %>% 
  select(subject=Subject, predict="NEG16") 
df.adj = read.table("data/emory/emory_subjects_adjuvant_status.txt", sep="\t",
                    header=T, row.names=NULL, stringsAsFactors=F) %>% 
  mutate(Adjuvant = ifelse(Adjuvant=="POS", "Adj", "NonAdj"))
df.test = info %>% 
  select(subject, day) %>% 
  left_join(df.adj, by=(c("subject"="Subject"))) %>% 
  mutate(gene.score = get_score(dat[gi,]))

cm = gg_color_hue(2)
ggplot(df.test, aes(day, gene.score, group=subject, col=Adjuvant)) +
  geom_line(size=1) +
  scale_colour_manual(values = cm) +
  coord_cartesian(xlim=c(0.9,2.1), expand = F) +
  xlab("") + ylab("GBX11 score") +
  ggtitle("Emory") +
  theme_bw()
ggsave(sprintf("IFN.gene_sd1-d0_Emory_lines_%s.genes_%s.png",gset, today()), w=5, h=3)

df.test.change <- df.test %>% 
  spread("day","gene.score") %>% 
  mutate(d1.d0 = Day1-Day0) %>% 
  mutate(Adjuvant = factor(Adjuvant, levels=c("NonAdj","Adj")))
pv <- with(df.test.change,
               wilcox.test(d1.d0[Adjuvant=="Adj"],
                           d1.d0[Adjuvant=="NonAdj"])
)$p.value
ggplot(df.test.change, aes(Adjuvant, d1.d0, fill=Adjuvant)) +
  geom_boxplot() +
  geom_dotplot(binaxis = "y", stackdir = "center", fill="black") +
  scale_fill_manual(values = rev(cm)) +
  ggtitle(sprintf("p-value = %.2g",pv)) +
  xlab("") + ylab("GBX11 score (day 1 - day 0)") +
  theme_bw()
ggsave(sprintf("IFN.gene_sd1-d0_Emory_boxplot_%s.genes_%s.png",gset, today()), w=4, h=3)

library(pROC)
r <- with(df.test.change,
     roc(Adjuvant,d1.d0, direction="<")
)
# plot(r)
data.frame(Sensitivity=r$sensitivities, Specificity=r$specificities) %>% 
  arrange(Sensitivity) %>% 
  ggplot(aes(Specificity,Sensitivity)) +
  geom_step(size=1) +
  scale_x_reverse() +
  geom_abline(slope=1,intercept = 1, lty=2) +
  ggtitle(sprintf("AUC=%.3g",r$auc)) +
  theme_bw()
ggsave(sprintf("IFN.gene_sd1-d0_Emory_roc_%s.genes_%s.png",gset, today()),w=3,h=3)


# correlate with titer d28
fn.em.titer = "data/emory/Emory_Nonadj_MN_Indonesia.txt"
em.titer = read.table(fn.em.titer, sep="\t", header=T) %>% 
  select(subject=Subject, titer=Day.42.MN)

df = df.test %>% 
  filter(day %in% "Day0", Adjuvant=="NonAdj") %>% 
  inner_join(em.titer, by="subject")

yticks = sort(unique(df$titer))
# yticks = yticks[seq(2,length(yticks),2)]

ggplot(df, aes(log2(titer), gene.score)) +
  # geom_point(shape=21, size=3, fill="white", col="black") +
  geom_point(size=4, col=darken(cm[2]), alpha = 0.5) +
  scale_x_continuous(name="MN titer (day 42)", breaks=log2(yticks), labels = yticks) +
  ylab("GBX11 score (day 0)") +
  theme_bw() +
  theme()
ggsave(sprintf("IFN.genes_ME_d0_Emory_%s.genes_vs_titer.d42_%s.png",gset, today()), w=4, h=3)

# end ------------




df.both = full_join(df.pbmc, df.pax, by="gene")
colnames(df.both)[-1] = c("PBMC","PAX")
gene.ord = df.both$gene[df.both[,-1] %>% as.matrix() %>% rowMeans() %>% order(decreasing = T)]
gene.ord2 = df.both$gene[df.both$PAX %>% order(decreasing = T)]
df.both = df.both %>% 
  mutate(gene = factor(gene, levels=gene.ord)) %>% 
  gather("assay","PC1.contr", -gene)

df.both = left_join(df.both, mod.df, by="gene")

df.both = df.both %>% mutate(gene = factor(gene, levels=rev(gene.ord)))
df.both = df.both %>% mutate(gene = factor(gene, levels=rev(gene.ord2)))



genes = read.table("wgcna/GE.pbmc_d0_WGCNA_M13_genes.txt", header=T, row.names=NULL, stringsAsFactors = F)
# genes = scan("genes_selected_6672.txt", what="character", sep=NULL)
gi = match(genes,rownames(dat))
df = as.data.frame(dat[gi,]) %>%
  tibble::rownames_to_column("gene") %>%
  gather("sample.name","fc",-gene) %>%
  mutate(sample.name = sub("\\.","-",sample.name)) %>% 
  inner_join(info, by="sample.name") %>%
  select(-sample.name, -plate.num, -iso.date) %>%
  # filter(time.point!="d000_00h") %>%
  group_by(gene,subject.id, time.point) %>%
  summarise(fc = mean(fc,na.rm=T)) %>%
  ungroup() %>%
  spread(time.point, fc)
df.mat.2 = as.matrix(df[,-(1:2)])
rownames(df.mat.2) = paste(df$subject.id, df$gene, sep="__")


# get scores
df.patt.mat = readRDS(sprintf("df.patt_%s_170201.rds", "200.2000.4.cleaned"))
patt.genes.clean = readRDS(sprintf("patt.genes.clean_%s_170202.rds","200.2000.4.cleaned"))
K = length(patt.genes.clean)
rownames(df.patt.mat) = paste0("G0",1:K)
long.idx = c(rep(1,3),2:ncol(df.mat.2))
df.patt.cor.10 = cor(t(df.mat.2[,long.idx]), t(df.patt.mat[,long.idx]), use="pairwise.complete.obs")

saveRDS(df.patt.cor.10, sprintf("df.patt.cor.10_%dg_%s.rds",length(gi), today()))
