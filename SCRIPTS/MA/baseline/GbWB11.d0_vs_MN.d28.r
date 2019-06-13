library(Biobase)

# data
fn = file.path(PROJECT_DIR, "DATA_PROCESSED/Microarrays/PBMC/samples.clean_genes.all/eset.genes.filtered.RData")
load(fn, verbose = T)

dat = exprs(eset.genes)
info = pData(eset.genes) %>% 
  mutate(time.point = factor(time.point, levels=unique(.$time.point)))

si = info$time.point %in% c("d000_00h","d001","d007","d021_00h","d022","d028") & info$subject.id != "H5N1_010"
info = info[si,]
dat = dat[,si]

fn.pax = file.path(PROJECT_DIR, "RESULTS/Microarrays/PAXgene/baseline/GE.pax_d0_WGCNA_GbWB11_genes.txt")
df.pax = read.table(fn.pax, sep="\t", header=T, row.names=NULL, stringsAsFactors = F)

test.genes = df.pax$gene; gset="GBX11"
gi = rownames(dat) %in% test.genes
sum(gi)

dat.sc = t(dat[gi,]) %>% scale() %>% t()

# read demographics
fn.demo = file.path(PROJECT_DIR, "DATA_ORIGINAL", "Clinical", "clinical_info_adj.txt")
df.demo = fread(fn.demo) %>% 
  dplyr::rename(subject  = `Subject ID`) %>%  
  mutate(Gender = ifelse(grepl("^F",Gender),"Female","Male"))

df.test = info %>% 
  select(subject = subject.id, time.point) %>% 
  mutate(subject = sub("_","-", subject)) %>% 
  mutate(gene.me = WGCNA::moduleEigengenes(t(dat[gi,]), rep(1,sum(gi)))$eigengenes %>% pull(1)) %>%
  group_by(subject, time.point) %>%
  summarize(gene.me=mean(gene.me)) %>%
  ungroup() %>% 
  left_join(df.demo, by="subject")

source(file.path(PROJECT_DIR, "SCRIPTS/functions/color_functions.r"))
cm = gg_color_hue(2)
ggplot(df.test %>% filter(time.point %in% c("d000_00h", "d001")), 
       aes(time.point, gene.me, group=subject, col=Adjuvant)) +
  geom_line(size=1) +
  scale_colour_manual(values = cm) +
  coord_cartesian(xlim=c(0.9,2.1), expand = F) +
  scale_x_discrete(labels=c("Day 0", "Day 1")) +
  xlab("") + ylab("GBX11 score") +
  ggtitle("CHI") +
  theme_bw()

dn.fig = file.path(PROJECT_DIR, "FIGURES/IFN_signature")
dir.create(dn.fig, showWarnings = F)
ggsave(file.path(dn.fig, "GbWB11.genes_PBMC.ME_d1-d0_lines.png"), w=5, h=3)

df.test.change <- df.test %>% 
  select(subject, time.point, gene.me, Adjuvant) %>% 
  spread("time.point","gene.me") %>% 
  mutate(d1.d0 = d001-d000_00h,
         d22.d21 = d022-d021_00h,
         d7.d0 = d007-d000_00h,
         d28.d21 = d028-d021_00h)
df.test.change <- df.test.change %>% 
  dplyr::rename(d0=d000_00h, d1=d001, d7=d007, d21=d021_00h, d22=d022, d28=d028) %>% 
  mutate(Adjuvant = factor(Adjuvant, levels=c("NonAdj","Adj")))

pv <- with(df.test.change,
           wilcox.test(d1.d0[Adjuvant=="Adj"],
                       d1.d0[Adjuvant=="NonAdj"])
)$p.value
ggplot(df.test.change, aes(Adjuvant, d1.d0, fill=Adjuvant)) +
  geom_boxplot() +
  geom_dotplot(binaxis = "y", stackdir = "center", fill="black") +
  scale_fill_manual(values = rev(cm)) +
  ggtitle(sprintf("p = %.2g",pv)) +
  xlab("") + ylab("GBX11 score (day 1 - day 0)") +
  theme_bw()
ggsave(file.path(dn.fig, "GbWB11.genes_PBMC.ME_d1-d0_boxplot.png"), w=4, h=3)

# read titers
fn.mn = file.path(PROJECT_DIR, "DATA_ORIGINAL", "Titres", "titre.txt")
df.mn = read.table(fn.mn, sep="\t", header=T, row.names=NULL, stringsAsFactors=T)
df.mn = df.mn %>% 
  mutate_at(-(1:2), function(x) sub("<20","10",x) %>% as.numeric()) %>% 
  mutate(day = as.numeric(sub("[dD]","",TimePt))) %>% 
  mutate(day = ifelse(day %in% c(29,31), 28, day) %>% factor()) %>% 
  mutate(subject = sprintf("H5N1-%03d",Sample.ID)) %>% 
  filter(day=="28")

DF = df.test %>% 
  filter(time.point %in% "d000_00h", Adjuvant=="NonAdj") %>% 
  inner_join(df.mn, by="subject")

xticks = sort(unique(DF$A.Indonesia))
xticks = xticks[seq(2,length(xticks),2)]

df.cor = broom::glance(cor.test(log2(DF$A.Indonesia), DF$gene.me))

ggplot(DF, aes(log2(A.Indonesia), gene.me, col=Gender)) +
  geom_point(size = 3, alpha=0.8) +
  geom_text(data = df.cor, aes(label=glue::glue("p = {format(p.value, digits=2)}")), 
            x=Inf, y=Inf, col="black", hjust=1.5, vjust=1.5) +
  scale_x_continuous(breaks=log2(xticks), labels = xticks) +
  xlab("MN titer (day28)") +
  ylab("GbWB score (day 0)") +
  theme_bw() +
  theme(legend.position = c(0.8,0.15), legend.box.background = element_blank(),
        legend.key.size = unit(1,"mm"))
fn.fig = file.path(dn.fig, "NonAdj_GbWB11_vs_MNd28_col_gender.png")
ggsave(fn.fig, w=4, h=3)
