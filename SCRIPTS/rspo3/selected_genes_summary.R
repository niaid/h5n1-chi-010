# PURPOSE: To summarize genes selected at each iteration in the robust correlation analysis.
# This is a very ad-hoc script and must be tailored for each run. 

library(tidyverse)

selected_genes <- read.table(file.path("RESULTS", "rspo3", "selected_genes_tsm_without_pre_z.txt"), sep = "\t", header = T)

genes_count <- group_by(selected_genes, genes) %>% summarize(count = n())
genes_count <- genes_count[order(genes_count$count, decreasing = T),]
genes_count <- cbind(genes_count, per = genes_count$count / 120 * 100)
write.table(genes_count, file.path("RESULTS", "rspo3", "selected_genes_count_tsm_without_pre_z.txt"), sep = "\t", row.names = F)

stop()

all_with <-  read.table(file.path("RESULTS", "rspo3", "selected_genes_count_all_with_pre_z.txt"), sep = "\t", header = T)
all_without <-  read.table(file.path("RESULTS", "rspo3", "selected_genes_count_all_without_pre_z.txt"), sep = "\t", header = T)
tsm_with <-  read.table(file.path("RESULTS", "rspo3", "selected_genes_count_tsm_with_pre_z.txt"), sep = "\t", header = T)
tsm_without <-  read.table(file.path("RESULTS", "rspo3", "selected_genes_count_tsm_without_pre_z.txt"), sep = "\t", header = T)

Reduce(intersect, list(all_with[all_with$count >= 60, "genes"], all_without[all_without$count >= 60, "genes"]))
Reduce(intersect, list(tsm_with[tsm_with$count >= 60, "genes"], tsm_without[tsm_without$count >= 60, "genes"]))
