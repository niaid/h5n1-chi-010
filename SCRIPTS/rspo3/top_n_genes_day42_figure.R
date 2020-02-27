library(pROC)
library(tidyverse)
library(data.table)
library(Biobase)

source(file.path("SCRIPTS", "functions", "load_sig.r"))
source(file.path("SCRIPTS", "functions", "get_score.r"))
fn.cd38.cor = file.path("RESULTS", "rspo3", "robust_corel_tsm_genes_2.txt")

# Load subject info and select adjuvanted subjects
subject <- read.table(file.path("DATA_ORIGINAL","Clinical", "clinical_info_adj.txt"), header = TRUE, sep = "\t")
subject <- subject[subject$Adjuvant == "Adj", c("Subject.ID")]
subject <- str_replace(subject, "-", "_")
# Remove subject H5N1_010
subject <- subject[!subject == "H5N1_010"]


# Load gene expression data for day 0.
pbmc_exp_data <- t(readRDS(file.path("New","190703", "DATA", "PBMC_dat0.in_isv.0.7_8144g_170216.rds")))
# Select data only for the adjuvanted subjects.
pbmc_exp_data <- pbmc_exp_data[, subject]

# Load gene expression data for day 100.
load(file.path("DATA_PROCESSED", "Microarrays", "PBMC", "samples.clean_genes.iqr", "eset.genes.filtered.RData"))
pbmc_day100 <-  exprs(eset.genes)
pbmc_day100 <- pbmc_day100[,grepl("d100", tolower(colnames(pbmc_day100)))]
colnames(pbmc_day100) <- substr(colnames(pbmc_day100), 1, 8)
colnames(pbmc_day100) <-  str_replace(colnames(pbmc_day100), "\\.", "_")
pbmc_day100 <-  pbmc_day100[,colnames(pbmc_day100) %in% subject]

# Load titre response data.
titre <- read.table(file.path("DATA_ORIGINAL", "Titres", "titre.txt"), header = TRUE, sep = "\t")
titre$TimePt <- tolower(titre$TimePt)
d42_titre <-  titre[titre$TimePt == "d42", c("Sample.ID","A.Indonesia")]
d42_titre$Sample.ID <- paste("H5N1_", str_pad(d42_titre$Sample.ID, 3, pad = "0"), sep = "")
d42_titre <- d42_titre[d42_titre$Sample.ID %in% subject, ]
subject == d42_titre$Sample.ID # Spot check
d42_titre$A.Indonesia <-  as.numeric(str_replace(d42_titre$A.Indonesia, "<", ""))

ng.max = 200
result_df = data.frame(matrix(nrow=0, ncol=3))

for(top in 2:ng.max){
        gene.sig <- load_sig(fn.cd38.cor, "cor.mean.sd.ratio", ntop=top)
        gene_cut_exp <- pbmc_exp_data[gene.sig,]

        # Day 0 correlation.
        # Calculate the z-score for each subject using expression data for the short listed genes.
        subj_z_score <- data.frame(score = get_score(gene_cut_exp))
        subj_z_score <- rownames_to_column(subj_z_score, var = "Sample.ID")
        df_score_titre <- merge(subj_z_score, d42_titre) # To make sure that the samples are in the same order.

        corel_day0 <- cor(df_score_titre$A.Indonesia, df_score_titre$score, method = "spearman")

        # Day 100 correlation.
        gene_cut_exp_day100 <- pbmc_day100[gene.sig,]
        # Calculate the z-score for each subject using expression data for the short listed genes.
        subj_z_score_day100 <- data.frame(score = get_score(gene_cut_exp_day100))
        subj_z_score_day100 <- rownames_to_column(subj_z_score_day100, var = "Sample.ID")
        df_score_titre_day100 <- merge(subj_z_score_day100, d42_titre) # To make sure that the samples are in the same order.

        corel_day100 <- cor(df_score_titre_day100$A.Indonesia, df_score_titre_day100$score, method = "spearman")

        result_df <- rbind(result_df, data.frame(top_ng = top, spear_day0 = corel_day0, spear_day100 = corel_day100))
}

plot_topn <- result_df %>%
        gather(key,value, spear_day0, spear_day100) %>%
        ggplot(aes(x=top_ng, y=value, colour=key)) +
        geom_line() +
        scale_y_continuous(name = "rho (z-score vs d42 titre)", breaks = seq(-1, 1, by = 0.1), limits = c(-1,1)) +
        scale_x_continuous(name = "number of top genes", limits = c(0, 200), breaks = seq(0, 200, by = 10)) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))

plot_name <- "topn_day0_day100_tsm_titre_d42.png"
ggsave(file.path("FIGURES","rspo3", plot_name), plot = plot_topn, device = "png", width = 6, height = 4)


# Max correlation of day0 gene expression with day42 titre.
max_0 <-  result_df[result_df$spear_day0 == max(result_df$spear_day0), "top_ng"]
# Genes at max correlation.
max_genes_0 <- load_sig(fn.cd38.cor, "cor.mean.sd.ratio", ntop=max_0)
max_genes_0 <- data.frame(genes = max_genes_0)

write.table(max_genes_0, file.path("RESULTS", "rspo3", "rspo3_signature_tsm_genes_titre_d42.txt"), sep = "\t", row.names = F)
writeLines(as.character(max_genes_0$genes), file.path("RESULTS", "rspo3", "rspo3_signature_tsm_genes_titre_d42_line.txt"), sep = ",")

stop()

# Find out common genes in between signatures based on all genes and tsm genes. 
# These are ad hoc lines and may be deleted in the future.
# Update everytime you use them.
all_genes_28 <- read.table(file.path("RESULTS", "rspo3", "rspo3_signature_all_genes.txt"), sep = "\t", header = T)
all_genes_42 <- read.table(file.path("RESULTS", "rspo3", "rspo3_signature_all_genes_titre_d42.txt"), sep = "\t", header = T)

tsm_genes_28 <- read.table(file.path("RESULTS", "rspo3", "rspo3_signature_tsm_genes.txt"), sep = "\t", header = T)
tsm_genes_42 <- read.table(file.path("RESULTS", "rspo3", "rspo3_signature_tsm_genes_titre_d42.txt"), sep = "\t", header = T)

common_genes <- intersect(all_genes$genes, tsm_genes$genes)


stop()


