library(ComplexHeatmap)
library(circlize)

dn.fig = file.path(PROJECT_DIR, "FIGURES/Baseline")
dir.create(dn.fig, showWarnings = F)

fn.pbmc = file.path(PROJECT_DIR, "RESULTS/Microarrays/PBMC/baseline/GE.pbmc_d0_WGCNA_BTM.LI.HG.pv.rds")
fn.pax = file.path(PROJECT_DIR, "RESULTS/Microarrays/PAXgene/baseline/GE.pax_d0_WGCNA_BTM.LI.HG.pv.rds")
res.pbmc = readRDS(fn.pbmc)
res.pax = readRDS(fn.pax)

tmod.pv.pbmc = sapply(res.pbmc, function(x) x$adj.P.Val)
tmod.pv.pax = sapply(res.pax, function(x) x$adj.P.Val)
rownames(tmod.pv.pbmc) = sprintf("%s (%s)",res.pbmc[[1]]$Title, res.pbmc[[1]]$ID)
rownames(tmod.pv.pax) = sprintf("%s (%s)",res.pax[[1]]$Title, res.pax[[1]]$ID)

# correlation between BTM enrichment in PBMC and WB modules
cc = cor(-log10(tmod.pv.pbmc), -log10(tmod.pv.pax))
df.cc = cc %>% as.data.frame() %>%
  tibble::rownames_to_column("PBMC") %>%
  gather(WB, cc, -PBMC)

ggplot(df.cc, aes(PBMC, WB)) +
  geom_point(aes(col=cc, size=cc)) +
  # scale_color_gradient2(low="blue", mid="white", high="red", na.value = "grey90",
  #                       name = "corr.\ncoeff.") +
  # scale_size_continuous(name = NULL, limits = c(0,1)) +
  scale_color_continuous(limits=c(0, 1), breaks=seq(0.25, 1, by=0.25), 
                         low="white", high="red", na.value = "grey90",
                         name = "corr.\ncoeff.") +
  scale_size_continuous(limits=c(0, 1), breaks=seq(0.25, 1, by=0.25), 
                        name = "corr.\ncoeff.") +
  guides(color= guide_legend(reverse=T), size=guide_legend(reverse=T)) +
  coord_equal() +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=0, vjust=0.5))
fn.fig = file.path(dn.fig, sprintf("GE.d0_pbmc_vs_pax_wgcna_BTM.%s_correlation.png",mset))
ggsave(fn.fig, w=6, h=6)


ri.pbmc = apply(tmod.pv.pbmc, 1, function(x) any(x<1e-8)) & !rownames(tmod.pv.pbmc) %in% c("TBA","Undetermined")
ri.pax = apply(tmod.pv.pax, 1, function(x) any(x<1e-8)) & !rownames(tmod.pv.pax) %in% c("TBA","Undetermined")
ri = ri.pbmc | ri.pax
ci.pbmc = apply(tmod.pv.pbmc[ri.pbmc,], 2, function(x) any(x<1e-2))
ci.pax = apply(tmod.pv.pax[ri.pax,], 2, function(x) any(x<1e-2))

X.pbmc = -log10(tmod.pv.pbmc[ri,ci.pbmc])
X.pax = -log10(tmod.pv.pax[ri,ci.pax])

fn.mn = file.path(PROJECT_DIR, "SCRIPTS/MA/baseline/h5n1_baseline_modules.txt")
mod.name = fread(fn.mn) %>% 
  gather(type, module, c("PBMC","PAX")) %>% 
  filter(module!="") %>% 
  mutate(label = glue::glue("{Name} ({module})"))

i1 = match(mod.name$module[mod.name$type=="PBMC"], colnames(X.pbmc))
i2 = match(mod.name$module[mod.name$type=="PAX"], colnames(X.pax))
colnames(X.pbmc)[i1] = mod.name$label[mod.name$type=="PBMC"]
colnames(X.pax)[i2] = mod.name$label[mod.name$type=="PAX"]


hx.pbmc = HeatmapAnnotation(cn = anno_text(colnames(X.pbmc), rot=30, just = "right", offset = unit(1, "npc")))
hx.pax = HeatmapAnnotation(cn = anno_text(colnames(X.pax), rot=30, just = "right", offset = unit(1, "npc")))
hx_height = unit(0.25, "npc") #max_text_height(colnames(X.pbmc))

hm.pbmc = Heatmap(X.pbmc, cluster_columns = F,
                  column_order = i1,
                  clustering_distance_rows = "manhattan",
                  clustering_method_rows = "complete",
              col = colorRamp2(c(0,10,20), c("white", "orange", "red")),
              column_title = "PBMC",
              name = "-log10(pv.adj)", 
              show_column_names = FALSE,
              # column_names_gp = gpar(cex=0.5, font=1, col= "white"), 
              row_names_side = "left", row_names_max_width = unit(11,"cm"),
              show_heatmap_legend = T, show_row_dend = F, rect_gp = gpar(ldy=1,col="grey"),
              bottom_annotation = hx.pbmc, bottom_annotation_height = hx_height,
              heatmap_legend_param = list(color_bar = "continuous", 
                                          legend_direction = "horizontal"))
hm.pax = Heatmap(X.pax, cluster_columns = F, 
                 column_order = i2,
                 col = colorRamp2(c(0,10,20), c("white", "orange", "red")),
                 column_title = "Whole blood",
                 # name = "-log10(pv.adj)",
                 show_column_names = FALSE,
                 show_row_names = F,
                 row_names_side = "left", row_names_max_width = unit(11,"cm"),
                 show_heatmap_legend = F, show_row_dend = F, rect_gp = gpar(ldy=1,col="grey"),
                 bottom_annotation = hx.pax, bottom_annotation_height = hx_height,
                 heatmap_legend_param = list(color_bar = "continuous", 
                                              legend_direction = "horizontal"))

pdf(file.path(dn.fig, sprintf("GE.d0.pbmc+pax_wgcna_BTM.%s_heatmap_ann.pdf",mset)))
draw(hm.pbmc+hm.pax, heatmap_legend_side = "bottom")
dev.off()

png(file.path(dn.fig, sprintf("GE.d0.pbmc+pax_wgcna_BTM.%s_heatmap_ann.png",mset)), width=800, height=700)
draw(hm.pbmc+hm.pax, heatmap_legend_side = "bottom")
dev.off()

