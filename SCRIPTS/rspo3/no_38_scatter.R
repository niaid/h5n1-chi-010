library(tidyverse)
library(stringr)
library(ggrepel) # Install in the container.

# MAIN

# Load RSPO3 and gene expression correlation data.
pbmc_adj_38 <- read.table(file.path("RESULTS", "SOMAscan", "pbmc_ADJ_RSPO3_bySex.txt"), header = TRUE, sep = "\t")
pbmc_adj_no_38 <- read.table(file.path("New", "190703", "RESULTS", "pbmc_ADJ_RSPO3_bySex_no_38.txt"), header = TRUE, sep = "\t")

for(cut_off in seq(0, 0.9, by=0.1)){
        # Select both positively and negatively correlated genes separately. 
        gene_cut_38 <- pbmc_adj_38[pbmc_adj_38[["r"]] > cut_off | pbmc_adj_38[["r"]] < -cut_off, c("gene", "r")]
        gene_cut_no_38 <- merge(gene_cut_38[1], pbmc_adj_no_38)
        gene_cut_no_38 <- gene_cut_no_38[c("gene", "r")]
        df_both <- merge(gene_cut_38, gene_cut_no_38, by = "gene")
        print(cor(df_both$r.x, df_both$r.y))

        scat_plot <- ggplot(df_both, aes(x=r.x, y=r.y)) +
                     geom_point() +
                     geom_smooth(method = "lm") +
                     scale_x_continuous( breaks = seq(-0.8, 0.8, by=0.2) , limits = c(-0.8, 0.8)) +
                     scale_y_continuous( breaks = seq(-0.8, 0.8, by=0.2), limits = c(-0.8, 0.8)) +
                                        #                      xlim(-0.8, 0.8) + ylim(-0.8, 0.8) +
                     labs(x = "gene rank with 38", y = "gene rank without 38")
        plot_name <- paste("with_without_38_", cut_off, ".png", sep = "")
        ggsave(file.path("FIGURES", "rspo3", plot_name), plot = scat_plot, device = "png", width = 5, height = 4)

}

