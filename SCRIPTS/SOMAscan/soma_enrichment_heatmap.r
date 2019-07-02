source("SCRIPTS/0_initialize.r")
library(ComplexHeatmap)
library(circlize)

dn = file.path(PROJECT_DIR, "RESULTS/SOMAscan/BTM")
mset = "LI"
files = dir(dn, sprintf(".*%s.*.txt",mset))
files = files[!grepl("unsigned", files)]
# all files must have the same set of modules in the same order

df.list = list()
for(f in files) {
  df = read.csv(file.path(dn,f), stringsAsFactors = F) # stringsAsFactors here is important
  df.list[[f]] = df
}
names(df.list) = sub(sprintf("Julian.tmod.%s.",mset),"", names(df.list))
names(df.list) = sub(".txt","", names(df.list))

# rownm = df$Title
rownm = sprintf("%s (%s)", df$Title, df$ID)
names(rownm) = df$ID

df.pv = sapply(df.list, function(x) x$adj.P.Val)
rownames(df.pv) = rownm

filter.unknown = TRUE # exclude modules without annotation
if(filter.unknown) {
  i.unknown = df$Title %in% c("TBA","Undetermined")
  df.pv = df.pv[!i.unknown,]
}

qval.row.in = 1e-2 # p.adj threshold for modules to be included
qval.col.in = 1e-0 # p.adj threshold for columns to be included

i.row.in = apply(df.pv, 1, function(x) any(x <= qval.row.in))
i.col.in = apply(df.pv[i.row.in,], 2, function(x) any(x < qval.col.in)) %>% which()
max.val = max(-log10(df.pv[i.row.in,i.col.in]))

df2 = -log10(df.pv[i.row.in,])
ineg = grep("negative", colnames(df.pv))
ipos = grep("positive", colnames(df.pv))
df.neg = df2[,ineg]
df.pos = df2[,ipos]

ii = df.neg > df.pos
df.comb = df.pos
df.comb[ii] = -df.neg[ii]
df.comb[abs(df.comb) <= -log10(0.05)] = 0
df.comb.cor = df.comb[,c(2,4)]
df.comb.par = df.comb[,c(1,3)]
colnames(df.comb.cor) = colnames(df.comb.par) = c("WB","PBMC")

# rank.type.bar = sub("^.+\\.([^\\.]+)$","\\1", colnames(df.pv)[i.col.in])
# sex.adj.bar = ifelse(grepl("p_", colnames(df.pv)[i.col.in]), "partial corr. (-sex)", "correlation")
# df.ann = data.frame(corr.type = sex.adj.bar, rank.type = rank.type.bar)

# ha = HeatmapAnnotation(df.ann, col = list(corr.type = c("correlation" =  "yellow", "partial corr. (-sex)" = "green"),
#                                            rank.type = c("positive"="red", "negative"="blue")))
hm.cor = Heatmap(df.comb.cor, cluster_columns = F, 
             col = colorRamp2(c(-max.val,0,max.val), c("blue", "white", "red")),
             name = "-log10(pv.adj)",
             row_names_side = "left", row_names_max_width = unit(15,"cm"),
             column_title="Genes ranked by\ncorrelation\nwith RSPO3", column_title_side = "top",
             # top_annotation = ha, top_annotation_height = unit(15,"mm"),
             show_heatmap_legend = T, show_row_dend = F, rect_gp = gpar(ldy=1,col="grey"),
             heatmap_legend_param = list(color_bar = "continuous", 
                                          legend_direction = "horizontal"))
hm.par = Heatmap(df.comb.par, cluster_columns = F, 
                 col = colorRamp2(c(-max.val,0,max.val), c("blue", "white", "red")),
                 # name = "-log10(pv.adj)",
                 show_row_names = F,
                 column_title="Genes ranked by\npartial correlation\n(-sex)", column_title_side = "top",
                 # row_names_side = "left", row_names_max_width = unit(10,"cm"),
                 # top_annotation = ha, top_annotation_height = unit(15,"mm"),
                 show_heatmap_legend = F, show_row_dend = F, rect_gp = gpar(ldy=1,col="grey"),
                 heatmap_legend_param = list(color_bar = "continuous", 
                                             legend_direction = "horizontal"))

p = draw(hm.cor+hm.par, heatmap_legend_side = "bottom")#, annotation_legend_size = "bottom")

dn.fig = file.path(PROJECT_DIR, "FIGURES/SOMAscan")
dir.create(dn.fig, showWarnings = F)

fn.fig = file.path(dn.fig, sprintf("RSPO3_cor_genes_BTM_enrichment_heatmap_%s",mset))
dev.copy(png, paste0(fn.fig, ".png"), width=650, height=500)
dev.off()
pdf(paste0(fn.fig, ".pdf"), width = 8, height = 6)
p
dev.off()

