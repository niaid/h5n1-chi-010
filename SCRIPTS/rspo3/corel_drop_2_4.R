# PURPOSE: To calculate correlation between ranks of genes from the 
# output of robust correlation with 2 and 4 subjects drop outs.

library(tidyverse)

drop_2 <- read.table(file.path("RESULTS", "rspo3", "robust_corel_all_genes_2.txt"), sep = "\t", header = T)
drop_4 <- read.table(file.path("RESULTS", "rspo3", "robust_corel_all_genes_4.txt"), sep = "\t", header = T)

drop_2_rank <- rank(drop_2$cor.mean.sd.ratio)
drop_4_rank <- rank(drop_4$cor.mean.sd.ratio)
df_both_rank <- data.frame(rank.x = drop_2_rank, rank.y = drop_4_rank)

scat_plot <- ggplot(df_both_rank, aes(x=rank.x, y=rank.y)) +
             geom_point() +
             geom_smooth(method = "lm") +
             labs(x = "gene rank with 2 subj drop", y = "gene rank with 4 subj drop")
plot_name <- paste("rank_drop_2_4", ".png", sep = "")
ggsave(file.path("FIGURES", "rspo3", plot_name), plot = scat_plot, device = "png", width = 5, height = 4)

print(cor(df_both_rank$rank.x, df_both_rank$rank.y))
