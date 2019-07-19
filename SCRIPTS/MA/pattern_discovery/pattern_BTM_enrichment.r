# PORPOSE: To carry out GSEA of the genes per pattern using blood transcription modules.

# Initialize
source("SCRIPTS/0_initialize.r")
library(tmod)
library(ComplexHeatmap)
library(circlize)

n_mod = 14
mod_name = sprintf("Gp%02d",1:n_mod)
mod_genes = list()
mod_pc1 = list()
dn.patt = file.path(PROJECT_DIR, "RESULTS/Microarrays/PBMC/pattern_discovery")

for (g in 1:n_mod) {
  fn.genes = file.path(dn.patt, "GE_pattern_genes", sprintf("GE_pattern_%02d_gene.sig.contr_ann.txt",g))
  # fn.genes = file.path("..", sprintf("GE_pattern_%02d_gene.sig.contr_ann.txt",g))
  tmp = fread(fn.genes, sep="\t", data.table = F)
  mod_genes[[mod_name[g]]] = tmp$gene
  mod_pc1[[mod_name[g]]] = tmp$PC1
}
sapply(mod_genes, length)

df = readRDS(file.path(dn.patt, "df_6672g.rds"))
# df = readRDS(file.path("..", "df_6672g.rds"))

data(tmod)
bg = unique(df$gene)
mod.n.in = sapply(tmod$MODULES2GENES, function(x) sum(x %in% unique(df$gene)))
mset = "LI"
tmod.ann = tmod$MODULES %>% 
  mutate(N.sel=mod.n.in) %>% 
  filter(SourceID==mset) %>% 
  dplyr::select(ID, Title, N=B, N.sel)
tmod.pv=c()
res = vector("list",n_mod)
for(m in 1:n_mod) {
  fg = mod_genes[[m]]
  res[[m]] = tmodHGtest(fg=fg,bg=bg, mset=mset, qval = Inf, order.by = "none")
  tmod.pv = cbind(tmod.pv, res[[m]]$adj.P.Val)
}
names(res) = mod_name
res[sapply(res,is.null)] = NULL
rownames(tmod.pv) = sprintf("%s (%s)", res[[1]]$Title, res[[1]]$ID)
colnames(tmod.pv) = names(res)

df.pv = tmod.pv

filter.unknown = TRUE # exclude modules without annotation
if(filter.unknown) {
  i.unknown = res[[1]]$Title %in% c("TBA","Undetermined")
  df.pv = df.pv[!i.unknown,]
}

qval.row.in = 0.05 # p.adj threshold for modules to be included
qval.col.in = 0.05 # p.adj threshold for columns to be included

i.row.in = apply(df.pv, 1, function(x) any(x <= qval.row.in))
i.col.in = apply(df.pv[i.row.in,], 2, function(x) any(x < qval.col.in))
hm = Heatmap(-log10(df.pv[i.row.in,i.col.in]), cluster_columns = F, 
             col = colorRamp2(c(0,4,8), c("white", "orange", "red")),
             name = "-log10(pv.adj)",
             row_names_side = "left", row_names_max_width = unit(15,"cm"),
             show_heatmap_legend = T, show_row_dend = F, rect_gp = gpar(ldy=1,col="grey"),
             heatmap_legend_param = list(color_bar = "continuous", 
                                         legend_direction = "horizontal"))

dn.fig = file.path(PROJECT_DIR, "FIGURES")
fn.fig = file.path(dn.fig, glue::glue("GE_pattern_BTM_enrichment_heatmap_{mset}.pdf"))
pdf(file=fn.fig)
draw(hm, heatmap_legend_side = "bottom")
dev.off()

fn.fig = file.path(dn.fig, glue::glue("GE_pattern_BTM_enrichment_heatmap_{mset}.png"))
png(fn.fig, width=500, height=600)
draw(hm, heatmap_legend_side = "bottom")
dev.off()


