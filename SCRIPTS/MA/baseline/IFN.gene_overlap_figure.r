# was IFN.gene_overlap_figure_171208.r

library(tmod)
data(tmod)
mod.id = c("LI.M75","LI.M127","LI.M150","LI.M67","LI.M165")
tmod$MODULES[mod.id, 1:3]
mod.gene = tmod$MODULES2GENES[mod.id]
mod.df = data.frame()
for(m in mod.id) {
  mod.df = rbind(mod.df, data.frame(ID=m, gene=mod.gene[[m]]))
}

mod.df = mod.df %>% 
  mutate(pc1.contr = 0.15) %>% 
  spread(ID, pc1.contr)

fn.pbmc = file.path(PROJECT_DIR, "RESULTS/Microarrays/PBMC/baseline/GE.pbmc_d0_WGCNA_Gb13_genes.txt")
df.pbmc = read.table(fn.pbmc, sep="\t", header=T, row.names=NULL, stringsAsFactors = F) %>% 
  # mutate(ID="PBMC") %>% 
  dplyr::select(gene, PBMC=pc1.contr)
fn.pax = file.path(PROJECT_DIR, "RESULTS/Microarrays/PAXgene/baseline/GE.pax_d0_WGCNA_GbWB11_genes.txt")
df.pax = read.table(fn.pax, sep="\t", header=T, row.names=NULL, stringsAsFactors = F) %>% 
  # mutate(ID="Whole blood") %>% 
  dplyr::select(gene, WB=pc1.contr)


mod.df.2 = mod.df %>% 
  left_join(df.pbmc, by="gene") %>% 
  left_join(df.pax, by="gene")
mod.mat = mod.df.2 %>% tibble::remove_rownames() %>% tibble::column_to_rownames("gene") %>% data.matrix() 
z = mod.mat %>% is.na() %>% `!` %>% rowSums()
mod.mat = mod.mat[z>1,]
mod.mat[is.na(mod.mat)] = 0

# gene's PC1 contribution as heatmap
dn.fig = file.path(PROJECT_DIR, "FIGURES/Baseline/")

library(ComplexHeatmap)
library(circlize)

pdf(file.path(dn.fig, "IFN.modules_PBMC_WB.pdf"))
Heatmap(mod.mat[,6:7], name = "PC1\ncontribution", cluster_columns = F, cluster_rows = F,
        col = colorRamp2(c(0,0.12), c("white", "darkred")))
dev.off()
png(file.path(dn.fig, "IFN.modules_PBMC_WB.png"), w=180, h=400)
Heatmap(mod.mat[,6:7], name = "PC1\ncontribution", cluster_columns = F, cluster_rows = F,
        col = colorRamp2(c(0,0.12), c("white", "darkred")))
dev.off()

# genes overlap in selected BTMs
mod.btm = mod.df %>% 
  gather(btm, val, -gene, na.rm = T) %>% 
  filter(gene %in% rownames(mod.mat)) %>% 
  mutate(gene = factor(gene, levels = rownames(mod.mat) %>% rev())) %>% 
  mutate(btm = factor(btm, levels = mod.id))
ggplot(mod.btm, aes(btm, gene)) +
  geom_point(size=3) +
  xlab("") + ylab("") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=0, vjust=0.5),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())
ggsave(file.path(dn.fig, "IFN.modules_only.png"), w=1.4, h=4)

# gene's PC1 contribution as barplot
df.bar = mod.mat[,6:7] %>% as.data.frame() %>% 
  tibble::rownames_to_column("gene") %>% 
  gather("type", "pc1.contr", -gene) %>% 
  mutate(gene = factor(gene, levels = rownames(mod.mat) %>% rev()))

p = ggplot(df.bar, aes(gene, pc1.contr, fill=type)) +
  geom_col(position="dodge", width = 0.8) +
  scale_y_reverse(breaks = seq(0, 0.1, length.out = 2)) +
  coord_flip() +
  xlab("") +
  ylab("PC1 contribution") +
  theme_bw() +
  theme(legend.title = element_blank())
p + scale_fill_manual(values = c("salmon","darkred"))
ggsave(file.path(dn.fig, "IFN.modules_PBMC_WB_bar.png"), w=2.5, h=4)
p + scale_fill_manual(values = c("grey","black"))
ggsave(file.path(dn.fig, "IFN.modules_PBMC_WB_bar_bw.png"), w=2.5, h=4)
