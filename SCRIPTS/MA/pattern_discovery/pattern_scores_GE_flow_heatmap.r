# PURPOSE: To plot heatmap of subjects scored per pattern of GE and flow data.
source("SCRIPTS/0_initialize.r")
library(ComplexHeatmap)
library(circlize)

# MAIN
# Laod data.
fn.g = file.path(PROJECT_DIR, "RESULTS/Microarrays/PBMC/pattern_discovery/subj.patt.cor.incl.s10.rds")
subj.patt.cor = readRDS(fn.g)
K = ncol(subj.patt.cor)
# colnames(subj.patt.cor) = sprintf("Gp%02d",1:K)
rownames(subj.patt.cor) = sub("H5N1-0","s",rownames(subj.patt.cor))

# read flow data

fn.f = file.path(PROJECT_DIR, "RESULTS/Flow_10c/DeltaThr0.5_DonorThr0.3_PercentOfParent_CORpearson_cutreeHybrid_deepSplit0_minClusterSize31_DonorScores.txt")
subj.patt.cor.f = read.csv(fn.f, row.names=1, stringsAsFactors = F)
colnames(subj.patt.cor.f) = sub("Module.","Fp0",colnames(subj.patt.cor.f))
si = match(rownames(subj.patt.cor), rownames(subj.patt.cor.f))

# read titers
fn.mn = file.path(PROJECT_DIR, "DATA_ORIGINAL", "Titres", "titre.txt")
df.mn = read.table(fn.mn, sep="\t", header=T, row.names=NULL, stringsAsFactors=T)
df.mn = df.mn %>% 
  mutate_at(-(1:2), function(x) sub("<20","10",x) %>% as.numeric()) %>% 
  mutate(day = as.numeric(sub("[dD]","",TimePt))) %>% 
  mutate(day = ifelse(day %in% c(29,31), 28, day) %>% factor()) %>% 
  mutate(subject = sprintf("s%02d",Sample.ID))
mn.d42 = df.mn %>% filter(day==42) %>% 
  dplyr::select(subject,A.Indonesia) %>% 
  tibble::column_to_rownames("subject") #%>% as.matrix()
mn.ord = match(rownames(subj.patt.cor),rownames(mn.d42))
mn.d42 = mn.d42[mn.ord,,drop=F]

# read demographics
fn.demo = file.path(PROJECT_DIR, "DATA_ORIGINAL", "Clinical", "clinical_info_adj.txt")
df.demo = fread(fn.demo) %>% 
  dplyr::rename(subject=`Subject ID`) %>% 
  mutate(subject = sub("H5N1-0","s", subject)) %>% 
  mutate(Gender = ifelse(grepl("^F",Gender),"Female","Male"))
sex.df = df.demo[match(rownames(mn.d42),df.demo$subject),"Gender",drop=F]
rownames(sex.df) = rownames(mn.d42)

hm.t = rowAnnotation(`A/Indonesia` = row_anno_barplot(log2(mn.d42$A.Indonesia), axis = T,
                                                      bar_width = 0.8,
                                                      gp = gpar(fill = "black", col=NA)),
                     # Vietnam = row_anno_barplot(log2(mn.d28$A.Vietnam), axis = T ),
                     show_annotation_name = T,
                     annotation_name_offset = unit(8, "mm"),
                     annotation_name_rot = 0,
                     width = unit(4, "cm"))
hm = Heatmap(subj.patt.cor, cluster_rows = T, cluster_columns = F,
             name = "correlation", col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
             # top_annotation = ha, top_annotation_height = unit(4, "mm"), 
             show_row_names = T, show_column_names = T,
             column_names_gp = gpar(fontsize = 20))#, column_title=ttl)
hm.flow = Heatmap(subj.patt.cor.f[si,], cluster_rows = F, cluster_columns = F,
             name = "flow", col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
             # top_annotation = ha, top_annotation_height = unit(4, "mm"), 
             show_row_names = F, show_column_names = T,
             column_names_gp = gpar(fontsize = 20))#, column_title=ttl)
hm.sex = Heatmap(sex.df$Gender, name="Gender", cluster_rows = F,
                 col=list(Male="cyan",Female="pink"),
                 show_column_names = F, width=unit(1,"cm"))

fn.fig = file.path(PROJECT_DIR, "FIGURES", "GE_subject_patterns_cor_flow_titers_sex.pred_incl.s10.pdf")
pdf(fn.fig)
set.seed(123)
draw(hm+hm.t+hm.sex+hm.flow)#, heatmap_legend_side = "bottom")
dev.off()

fn.fig = file.path(PROJECT_DIR, "FIGURES", "GE_subject_patterns_cor_flow_titers_sex.pred_incl.s10.png")
png(fn.fig, width=1000, height=800)
set.seed(123)
draw(hm+hm.t+hm.sex+hm.flow)#, heatmap_legend_side = "bottom")
dev.off()
