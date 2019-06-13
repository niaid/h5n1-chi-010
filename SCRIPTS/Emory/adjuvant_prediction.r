library(ComplexHeatmap)
library(circlize)

dn = file.path(PROJECT_DIR, "DATA_PROCESSED/Emory")
df = readRDS(file.path(dn, "emory_GE_fc.rds"))
mat.emory = df[,-1] %>% data.matrix()
rownames(mat.emory) = df$subject

n.neg = 16
mat.emory.sc = scale(mat.emory)
mat.emory.mean = mat.emory.sc %>% apply(1,mean)
irank.mean = rank(mat.emory.mean)
iord = order(mat.emory.mean)
adj.pred = ifelse(mat.emory.mean <= mat.emory.mean[iord[n.neg]], "NEG", "POS")

dn.adj = file.path(PROJECT_DIR, "DATA_ORIGINAL/Emory")
df.adj = fread(file.path(dn.adj, "emory_subjects_adjuvant_status.txt"), sep="\t")
adj.true = df.adj$Adjuvant
names(adj.true) = df.adj$Subject

source(file.path(PROJECT_DIR, "SCRIPTS/functions/color_functions.r"))
cm = gg_color_hue(2) %>% rev()
cm2 = lighten(cm, 1.4)

hm = Heatmap(mat.emory.sc, name="score", 
             cluster_columns = F, cluster_rows = F, row_order = iord, 
             width = unit(3,"cm"),
             row_names_side = "left")
hmp = rowAnnotation(mean = row_anno_points(mat.emory.mean, axis = T),
                    # annotation_name = "max",
                    show_annotation_name = T,
                    annotation_name_offset = unit(8, "mm"),
                    annotation_name_rot = 0,
                    width = unit(3, "cm"))
hm.pred = Heatmap(adj.pred, name = "Predicted", show_row_names = F, width = unit(10,"mm"),
                  col = list("NEG"=cm2[1], "POS"=cm2[2]),
                  show_heatmap_legend = T, show_column_names = F)
hm.true = Heatmap(adj.true, name = "Actual", show_row_names = F, width = unit(10,"mm"),
                  col = list("NonAdj"=cm[1], "Adj"=cm[2]),
                  show_heatmap_legend = T, show_column_name = F)


dn.fig = file.path(PROJECT_DIR, "FIGURES/Emory")
dir.create(dn.fig, showWarnings = F)
fn.fig = file.path(dn.fig, sprintf("emory_heatmap_ord.mean_G3_ann_%dneg", n.neg))

pdf(paste0(fn.fig, ".pdf"))
draw(hm+hmp+hm.pred+hm.true)
dev.off()

png(paste0(fn.fig, ".png"), w=400, h=600)
draw(hm+hmp+hm.pred+hm.true)
dev.off()

