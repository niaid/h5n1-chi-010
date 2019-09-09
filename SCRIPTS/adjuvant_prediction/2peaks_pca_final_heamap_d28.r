# Initialize.
source("SCRIPTS/0_initialize.r")
library(ComplexHeatmap)
library(circlize)

dn.pred = file.path(PROJECT_DIR, "RESULTS/Adjuvant_prediction")

fn.pc = file.path(dn.pred, "2peaks_patterns_PCA.txt")
df.pc = fread(fn.pc)

dat = df.pc[,2:6] %>% data.matrix()
rownames(dat) = df.pc$subject

df.ip10.score = readRDS(file.path(PROJECT_DIR, "RESULTS/Luminex/ip10_2peak_score.rds"))
ip10 = df.ip10.score$fc.mean %>% log10()
names(ip10) = sub("H5N1-0","s",df.ip10.score$subject)

lidx = match(rownames(dat), names(ip10))
dat2 = cbind(dat, ip10[lidx])
colnames(dat2)[ncol(dat2)] = "IP-10"

dat2.sc = scale(dat2)

# read titers
fn.mn = file.path(PROJECT_DIR, "DATA_ORIGINAL", "Titres", "titre.txt")
df.mn = read.table(fn.mn, sep="\t", header=T, row.names=NULL, stringsAsFactors=T)
df.mn = df.mn %>% 
  mutate_at(-(1:2), function(x) sub("<20","10",x) %>% as.numeric()) %>% 
  mutate(day = as.numeric(sub("[dD]","",TimePt))) %>% 
  mutate(day = ifelse(day %in% c(29,31), 28, day) %>% factor()) %>% 
  mutate(subject = sprintf("s%02d",Sample.ID))

mn.d28 = df.mn %>% filter(day==28) %>% 
  dplyr::select(subject,A.Indonesia) %>% 
  tibble::column_to_rownames("subject") #%>% as.matrix()
mn.ord = match(rownames(dat2),rownames(mn.d28))
mn.d28 = mn.d28[mn.ord,,drop=F]

# read demographics
fn.demo = file.path(PROJECT_DIR, "DATA_ORIGINAL", "Clinical", "clinical_info_adj.txt")
df.demo = fread(fn.demo) %>% 
  dplyr::rename(subject=`Subject ID`) %>% 
  mutate(subject = sub("H5N1-0","s", subject)) %>% 
  mutate(Gender = ifelse(grepl("^F",Gender),"Female","Male"))
df.sex = df.demo[match(rownames(dat2),df.demo$subject),"Gender",drop=F]
rownames(df.sex) = rownames(dat2)

# set.seed(123)
htree = hclust(dist(dat2.sc))
adj.pred = cutree(htree, k = 2)
adj.pred=ifelse(adj.pred==2, "NEG","POS")

source(file.path(PROJECT_DIR, "SCRIPTS/functions/color_functions.r"))
cm = gg_color_hue(2) %>% rev()
cm2 = lighten(cm, 1.4)

hm = Heatmap(dat2.sc, cluster_rows = T, cluster_columns = F, #km=2,
             clustering_distance_rows = "euclidean",
             clustering_method_rows = "average",
             name = "Correlation", col = colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red")),
             # top_annotation = ha, top_annotation_height = unit(4, "mm"), 
             show_row_names = T, show_column_names = T,
             column_names_gp = gpar(fontsize = 20),
             heatmap_legend_param = list(title_gp = gpar(fontsize = 18), 
                                         labels_gp = gpar(fontsize = 14),
                                         color_bar = "continuous",
                                         legend_height = unit(4,"cm"),
                                         legend_direction="vertical"))

hm.t = rowAnnotation(`A/Indonesia d28` = row_anno_barplot(log2(mn.d28$A.Indonesia), axis = T,
                                                          bar_width = 0.8,
                                                          gp = gpar(fill = "black", col=NA)),
                     # annotation_name = "A/Indonesia",
                     show_annotation_name = T,
                     annotation_name_offset = unit(8, "mm"),
                     annotation_name_rot = 0,
                     width = unit(4, "cm"))

hm.sex = Heatmap(df.sex$Gender, name="Gender", cluster_rows = F,
                 col=list(Male="cyan",Female="pink"),
                 show_column_names = F, width=unit(1,"cm"),
                 heatmap_legend_param = list(title_gp = gpar(fontsize = 18), 
                                             labels_gp = gpar(fontsize = 14),
                                             legend_direction="horizontal"))
hm.pred = Heatmap(adj.pred, name="Predicted", cluster_rows = F,
                  col=list(NEG=cm2[1],POS=cm2[2]),
                  show_column_names = F, width=unit(1,"cm"))
# hm.pred = Heatmap(pred.df$Predicted, name="predicted", cluster_rows = F,
#                   show_column_names = F, width=unit(1,"cm"))
# set.seed(122)
fn.fig = file.path(PROJECT_DIR, "FIGURES/GE_subject_patterns_cor_flow_ip10_titers_sex.pred_incl.s10_d28.pdf")
pdf(fn.fig, width = 8.5, height = 11)
draw(hm + hm.t + hm.sex + hm.pred, heatmap_legend_side = "right")
decorate_heatmap_body("Correlation", {grid.lines(x=c(0,1), y=c(0.5, 0.5), gp = gpar(lwd = 2, lty=2))})
dev.off()

fn.fig = file.path(PROJECT_DIR, "FIGURES/GE_subject_patterns_cor_flow_ip10_titers_sex.pred_incl.s10_d28.png")
png(fn.fig, width=700,height=800)
draw(hm + hm.t + hm.sex + hm.pred, heatmap_legend_side = "right")
decorate_heatmap_body("Correlation", {grid.lines(x=c(0,1), y=c(0.5, 0.5), gp = gpar(lwd = 2, lty=2))})
dev.off()


df.adj_predict = data.frame(subject = sub("H5N1_","H5N1-",df.pc$subject),
                         predict = adj.pred)

fwrite(df.adj_predict, file.path(dn.pred, "adjuvant_predicted_subjects_d28.txt"), sep="\t")
