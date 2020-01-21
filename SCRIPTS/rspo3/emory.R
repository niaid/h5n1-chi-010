library(stringr)
eset_gene <- readRDS(file.path("DATA_PROCESSED", "Emory", "eset.gene.rds"))
exp_mat <- exprs(eset_gene)
non_adj_titre <- read.table(file.path("DATA_PROCESSED", "Emory", "Emory_Nonadj_MN_Indonesia.txt"), header = TRUE, sep = "\t")

subj_non_adj <- as.character(non_adj_titre$Subject)
exp_mat_cols <- colnames(exp_mat)
subj_adj <- exp_mat_cols[!exp_mat_cols %in% paste(subj_non_adj, ".Day0", sep = "")]
subj_adj <- subj_adj[grepl("Day0", subj_adj, perl = TRUE)]

exp_mat_adj <- exp_mat[, subj_adj]


