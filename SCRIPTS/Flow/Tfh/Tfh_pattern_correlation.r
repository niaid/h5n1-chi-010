source("SCRIPTS/0_initialize.r")
library(ComplexHeatmap)
library(circlize)

# Sys.setlocale(locale="C")
fn.tfh = file.path(PROJECT_DIR, "DATA_PROCESSED/Flow_10c/H5N1_Thelper-Thelper_freq.V1.1_v3")
tfh = fread(fn.tfh, data.table=F, header = T) %>% 
  dplyr::rename(subject=V1, time=Timepoint) %>% 
  dplyr::filter(subject!="")

cp_ids = names(tfh[-c(1:3)]) %>% as.numeric()
sel_id = c(25, 29, cp_ids[cp_ids>70]) %>% as.character()


tfh = tfh %>% 
  mutate(time.point = case_when(
    time=="d8" ~ "d7",
    time=="d21_2h" ~ "d21_4h",
    time=="d29" ~ "d28",
    time=="d31" ~ "d28",
    time=="d43" ~ "d42",
    time=="d0_0h" ~ "d0",
    time=="d21_0h" ~ "d21",
    TRUE ~ time
  )) %>% 
  # mutate(time.point = sub("_","+",time.point)) %>% 
  mutate(time.point = fct_inorder(time.point))

tfh = tfh %>% 
  dplyr::select(subject, time.point, one_of(sel_id))

tfh.fc = tfh %>% 
  gather("CP", "value", -c("subject", "time.point")) %>% 
  mutate(CP = paste0("ID", CP)) %>% 
  spread(time.point, value) %>% 
  mutate_if(is.numeric, `-`, .$d0)
  
tfh.tm = levels(tfh$time.point)

tfh.mat = tfh.fc %>% 
  arrange(CP) %>% 
  unite(comb, subject, CP, sep="_") %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames("comb") %>% 
  data.matrix()

# Flow modules
dn.in = file.path(PROJECT_DIR, "RESULTS/Flow_10c/")
fn = "DeltaThr0.5_DonorThr0.3_PercentOfParent_CORpearson_cutreeHybrid_deepSplit0_minClusterSize31_Modules.txt"
patt.flow = read.csv(file.path(dn.in, fn), row.names=1, stringsAsFactors = F)
patt.flow = t(patt.flow)
patt.flow = cbind(0, patt.flow)
colnames(patt.flow)[1] = "d0"
colnames(patt.flow) = sub("_0h", "", colnames(patt.flow))
identical(colnames(tfh.mat), colnames(patt.flow))
K = nrow(patt.flow)

### demographics and adjuvant status
fn.clin = file.path(PROJECT_DIR, "DATA_ORIGINAL/Clinical/clinical_info_adj.txt")
df.clin = fread(fn.clin) %>% 
  dplyr::rename(subject=`Subject ID`)


cc = cor(t(tfh.mat), t(patt.flow), use = "pairwise.complete.obs")
cc.p = cc
for(k in 1:nrow(patt.flow)) {
  for(m in 1:nrow(tfh.mat)) {
    cc.p[m,k] = cor.test(tfh.mat[m,], patt.flow[k,], alternative = "two.sided", exact=F)$p.value
  }
}

# hist(apply(cc, 1, max))
# plot(cc, log10(cc.p))
sum(cc.p<0.01)


source(file.path(PROJECT_DIR, "SCRIPTS/functions/color_functions.r"))
cm.adj = gg_color_hue(2)

dn.fig = file.path(PROJECT_DIR, "FIGURES/Tfh")
dir.create(dn.fig, showWarnings = F)

for(k in 1:K) {
  subj_cp =  data.frame(comb = rownames(cc), cc = cc[,k]) %>% 
    separate(comb, c("subject","CP"), sep="_") %>% 
    spread(CP, cc, fill = 0) %>% 
    tibble::remove_rownames() %>% 
    tibble::column_to_rownames("subject") %>% 
    data.matrix()
  i.adj = match(rownames(subj_cp), df.clin$subject)
  hm = Heatmap(subj_cp, name = "corr.\ncoef", column_title = paste0("Fp0",k), split=df.clin$Adjuvant[i.adj]) +
       Heatmap(df.clin$Adjuvant[i.adj], col=cm.adj, name = "Adj")
  fn.fig = file.path(dn.fig, sprintf("Tfh_pattern_corr_heatmap_Fp0%d_split", k))
  png(paste0(fn.fig, ".png"), w=400, h=550)
  draw(hm)
  dev.off()
  pdf(paste0(fn.fig, ".pdf"), w=4, h=5.5)
  draw(hm)
  dev.off()
}


cc.th = 0.7

mat.list = list()

for (cc.sign in c(1,-1)) {

cc.lbl = ifelse(cc.sign > 0, "pos", "neg")
  
cc.max = apply(cc * cc.sign, 1, max, na.rm=T)
cc.max.idx = apply(cc * cc.sign, 1, which.max)

cc.patt = cc.max.idx
cc.patt[cc.max < cc.th] = 0

cc.best = cc[cbind(seq_along(cc.max.idx),cc.max.idx)]
cc.p.best = cc.p[cbind(seq_along(cc.max.idx),cc.max.idx)]

df = data.frame(comb = names(cc.max), cc = cc.best, p = cc.p.best, pattern = cc.patt)

mat = df %>% 
  separate(comb, c("subject", "CP"), sep="_") %>% 
  filter(pattern>0) %>% 
  select(-cc, -p) %>% 
  spread(CP, pattern, fill = 0) %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames("subject") %>% 
  data.matrix()

mat.list[[cc.lbl]] = df %>% 
  separate(comb, c("subject", "CP"), sep="_") %>% 
  # filter(pattern>0) %>% 
  select(-cc, -p) %>% 
  spread(CP, pattern, fill = 0) %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames("subject") %>% 
  data.matrix()

i.adj = match(rownames(mat), df.clin$subject)

cm = palette()[sort(unique(df$pattern))]

hm = Heatmap(mat, name = "pattern", col=c("white", cm)) +
     Heatmap(df.clin$Adjuvant[i.adj], col=cm.adj, name = "Adj")
# fn.fig = file.path(dn.fig, glue::glue("Tfh_pattern_best_heatmap_{cc.lbl}"))
# png(paste0(fn.fig, ".png"), w=400, h=550)
# draw(hm)
# dev.off()
# pdf(paste0(fn.fig, ".pdf"), w=4, h=5.5)
# draw(hm)
# dev.off()

hm = Heatmap(mat, name = "pattern", col=c("white", cm), split=df.clin$Adjuvant[i.adj]) +
  Heatmap(df.clin$Adjuvant[i.adj], col=cm.adj, name = "Adj")
# fn.fig = file.path(dn.fig, glue::glue("Tfh_pattern_best_heatmap_split_{cc.lbl}"))
# png(paste0(fn.fig, ".png"), w=400, h=550)
# draw(hm)
# dev.off()
# pdf(paste0(fn.fig, ".pdf"), w=4, h=5.5)
# draw(hm)
# dev.off()

}

# draw both heatmaps together
mat2 = cbind(mat.list[["pos"]], mat.list[["neg"]])
ro = hclust(dist(mat2, method = "manhattan"))$order
# mat2 = rbind(mat.list[["pos"]], mat.list[["neg"]])
# co = hclust(dist(t(mat2)))$order
cm = palette()[sort(unique(unlist(mat.list)))]
i.adj = match(rownames(mat.list[["pos"]]), df.clin$subject)
hm = Heatmap(mat.list[["pos"]], name = "pattern", col=c("white", cm), split=df.clin$Adjuvant[i.adj],
             cluster_rows = F, cluster_columns = F,
             row_order = ro) +
  Heatmap(mat.list[["neg"]], name = "pattern", col=c("white", cm), 
          cluster_rows = F, cluster_columns = F,
          show_row_names = F) +
  Heatmap(df.clin$Adjuvant[i.adj], col=cm.adj, name = "Adj")

fn.fig = file.path(dn.fig, glue::glue("Tfh_pattern_best_heatmap_split"))
png(paste0(fn.fig, ".png"), w=600, h=550)
draw(hm)
dev.off()
pdf(paste0(fn.fig, ".pdf"), w=6, h=5.5)
draw(hm)
dev.off()

# Select populations correlated with patterns 1 and 3 and plot correlation heatmaps with adjuvant status
dir.vec = c(rep("pos", ncol(mat.list[["pos"]])), rep("neg", ncol(mat.list[["neg"]])))
col.vec = colnames(mat2)
mat.p13 = (mat2 == 1 | mat2 == 3)
ci = colSums(mat.p13)>1
sel.pop = col.vec[ci] %>% unique()

fn.ann = file.path(PROJECT_DIR, "DATA_ORIGINAL/Flow_10c/Flow_10c_ann.txt")
df.ann = fread(fn.ann, header = T) %>% 
  mutate(Annot = paste0(ID, ": ", Name2))

subj_cp =  data.frame(comb = rownames(cc), cc = cc[,1]) %>% 
  separate(comb, c("subject","CP"), sep="_") %>% 
  dplyr::filter(CP %in% sel.pop) %>% 
  mutate(CP = factor(CP, levels=sel.pop)) %>% 
  spread(CP, cc, fill = 0) %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames("subject") %>% 
  data.matrix()
i.ann = match(colnames(subj_cp), df.ann$ID)
colnames(subj_cp) = df.ann$Annot[i.ann]
i.adj = match(rownames(subj_cp), df.clin$subject)
ro = rowMeans(-abs(subj_cp))

subj_cp.3 =  data.frame(comb = rownames(cc), cc = cc[,3]) %>% 
  separate(comb, c("subject","CP"), sep="_") %>% 
  dplyr::filter(CP %in% sel.pop) %>% 
  mutate(CP = factor(CP, levels=sel.pop)) %>% 
  spread(CP, cc, fill = 0) %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames("subject") %>% 
  data.matrix()
colnames(subj_cp.3) = df.ann$Annot[i.ann]

hx = HeatmapAnnotation(cn = anno_text(df.ann$Annot[i.ann], rot=-45, just = "left", offset = unit(1, "npc")))
hx_height = unit(0.25, "npc") #max_text_height(colnames(X.pbmc))

hm = Heatmap(subj_cp, name = "corr.\ncoef", column_title = paste0("Fp0",1), split=df.clin$Adjuvant[i.adj],
             cluster_columns = F,
             show_column_names = FALSE,
             bottom_annotation = hx, bottom_annotation_height = hx_height,
             row_dend_reorder = ro) +
  # Heatmap(subj_cp.3, name = "corr.\ncoef", column_title = paste0("Fp0",3),
  #         show_column_names = FALSE,
  #         bottom_annotation = hx, bottom_annotation_height = hx_height,
  #         cluster_columns = F) +
  Heatmap(df.clin$Adjuvant[i.adj], col=cm.adj, name = "Adj")

fn.fig = file.path(dn.fig, glue::glue("Tfh_pattern_corr_heatmap_Fp01_Fp03_split"))
# png(paste0(fn.fig, ".png"), w=850, h=550)
png(paste0(fn.fig, ".png"), w=900, h=750)
draw(hm)
dev.off()
# pdf(paste0(fn.fig, ".pdf"), w=8.5, h=5.5)
pdf(paste0(fn.fig, ".pdf"), w=9, h=7.5)
draw(hm)
dev.off()


### plot Fp01 pattern and example profiles

ci = grep("ID81", colnames(subj_cp))
si = which(subj_cp[,ci]==max(subj_cp[,ci]))
sn = rownames(subj_cp)[si]

mat.ex = tfh.mat[c('H5N1-023_ID81'),]
mat.pt = patt.flow["Fp01",]
mat.ex = mat.ex * max(mat.pt) / max(mat.ex)
  
mat.plot = rbind(mat.ex, mat.pt)

df.ex = mat.plot %>% 
  as.data.frame() %>% 
  mutate(line=c("ID81 in H5N1-023", "Fp01")) %>% 
  gather("time","fc", -line) %>% 
  mutate(time = fct_inorder(time), line = fct_inorder(line))

ggplot(df.ex, aes(time, fc, group=line)) +
  geom_line(size=1, aes(col=line)) +
  scale_color_manual(values=c("red", "black"), guide=F) +
  facet_wrap(~line, ncol=1) +
  ylab("") +
  theme_bw() +
  theme(strip.background = element_blank(), 
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
fn.fig = "Tfh_pattern_correlation_example_profile"
ggsave(file.path(dn.fig, paste0(fn.fig,".png")), w=4, h=4)
ggsave(file.path(dn.fig, paste0(fn.fig,".pdf")), w=4, h=4)
