library(tidyverse)
library(stringr)
library(ggrepel) # Install in the container.

# MAIN

# Load RSPO3 and gene expression correlation data.
pbmc_adj_38 <- read.table(file.path("RESULTS", "SOMAscan", "pbmc_ADJ_RSPO3_bySex.txt"), header = TRUE, sep = "\t")
pbmc_adj_no_38 <- read.table(file.path("New", "190703", "RESULTS", "pbmc_ADJ_RSPO3_bySex_no_01.txt"), header = TRUE, sep = "\t")

cut_off <- 0
# Select both positively and negatively correlated genes separately. 
gene_cut_38 <- pbmc_adj_38[pbmc_adj_38[["r"]] > cut_off | pbmc_adj_38[["r"]] < -cut_off, c("gene", "r")]
gene_cut_38_rank <-  data.frame(gene = gene_cut_38$gene, rank = rank(gene_cut_38$r))
gene_cut_no_38 <- merge(gene_cut_38[1], pbmc_adj_no_38)
gene_cut_no_38 <- gene_cut_no_38[c("gene", "r")]
gene_cut_no_38_rank <-  data.frame(gene = gene_cut_no_38$gene, rank = rank(gene_cut_no_38$r))
df_both <- merge(gene_cut_38, gene_cut_no_38, by = "gene")
df_both_rank <- merge(gene_cut_38_rank, gene_cut_no_38_rank, by = "gene")

print(cor(df_both$r.x, df_both$r.y))
print(cor(df_both_rank$rank.x, df_both_rank$rank.y))

scat_plot <- ggplot(df_both, aes(x=r.x, y=r.y)) +
             geom_point() +
             geom_smooth(method = "lm") +
             scale_x_continuous( breaks = seq(-0.8, 0.8, by=0.2) , limits = c(-0.8, 0.8)) +
             scale_y_continuous( breaks = seq(-0.8, 0.8, by=0.2), limits = c(-0.8, 0.8)) +
             labs(x = "r with 01", y = "r without 01")
plot_name <- paste("with_without_01_", cut_off, ".png", sep = "")
ggsave(file.path("FIGURES", "rspo3", plot_name), plot = scat_plot, device = "png", width = 5, height = 4)

scat_plot <- ggplot(df_both_rank, aes(x=rank.x, y=rank.y)) +
             geom_point() +
             geom_smooth(method = "lm") +
             labs(x = "gene rank with 01", y = "gene rank without 01")
plot_name <- paste("with_without_01_rank_", cut_off, ".png", sep = "")
ggsave(file.path("FIGURES", "rspo3", plot_name), plot = scat_plot, device = "png", width = 5, height = 4)


