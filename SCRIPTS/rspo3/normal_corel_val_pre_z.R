# PURPOSE: To calculate robust correlation between RSPO3 (soma scan expression) and gene expression
# for adjuvanted subjects.

library(pROC)
library(tidyverse)
library(data.table)
library(stringr)
library(Biobase)

# Functions
source(file.path("SCRIPTS", "functions", "load_sig.r"))
source(file.path("SCRIPTS", "functions", "get_score.r"))

# Load subject info and select adjuvanted subjects
subject <- read.table(file.path("DATA_ORIGINAL","Clinical", "clinical_info_adj.txt"), header = TRUE, sep = "\t")
subject <- subject[subject$Adjuvant == "Adj", c("Subject.ID")]
subject <- str_replace(subject, "-", "_")
# Remove subject H5N1_010
subject <- subject[!subject == "H5N1_010"]

# Load gene expression data on day 0.
pbmc_exp_data <- t(readRDS(file.path("New","190703", "DATA", "PBMC_dat0.in_isv.0.7_8144g_170216.rds")))

# Filter genes for ISV > 0.5
# gene_stability <- read.table(file.path("SCRIPTS", "rspo3", "CHI_genes_stability.txt"), header = TRUE, sep = "\t")
# gene_stability <- gene_stability[gene_stability$ISV > 0.5,]
# tsm_genes <- as.character(gene_stability$gene)
# pbmc_exp_data <- pbmc_exp_data[tsm_genes %in% rownames(pbmc_exp_data),]

# Load SOMAScan data for RSPO3.
rspo3_soma <- read.table(file.path("New", "190703", "DATA", "SOMA.txt"), header = TRUE, sep = "\t")
rspo3_soma <- rspo3_soma[c("donor", "RSPO3")]
rspo3_soma$donor <- str_replace(rspo3_soma$donor, "-", "_")
colnames(rspo3_soma) <- c("Sample.ID", "RSPO3")

# Load titre response data.
titre <- read.table(file.path("DATA_ORIGINAL", "Titres", "titre.txt"), header = TRUE, sep = "\t")
titre$TimePt <- tolower(titre$TimePt)
d28_titre <-  titre[titre$TimePt %in% c("d28", "d29", "d31") , c("Sample.ID","A.Indonesia")]
d28_titre$Sample.ID <- paste("H5N1_", str_pad(d28_titre$Sample.ID, 3, pad = "0"), sep = "")
d28_titre$A.Indonesia <-  as.numeric(str_replace(d28_titre$A.Indonesia, "<", ""))

# Find common subjects among all the datasets. 
common_subj <- Reduce(intersect, list(subject, colnames(pbmc_exp_data), rspo3_soma$Sample.ID, d28_titre$Sample.ID))
# common_subj <- common_subj[common_subj != "H5N1_038"]

# Select data for common subjects.
# Gene expression.
pbmc_exp_data <- pbmc_exp_data[, common_subj]
pbmc_exp_data_z <-  t(scale(t(pbmc_exp_data)))

# RSPO3 Soma Scan.
rspo3_soma <- rspo3_soma[rspo3_soma$Sample.ID %in% common_subj,]
# Day 28 titre.
d28_titre <- d28_titre[d28_titre$Sample.ID %in% common_subj,]

# Drop N subjects and carry out robust correlation on the remaining. 
# Drop N subjects.
N <- 2
N_cmb <- combn(common_subj, N )
val_result <- data.frame(matrix(nrow = 0, ncol=3))

for(N_idx in 1:ncol(N_cmb)){
        #         N_idx <- 2 # Note: This has to go in a for loop.
       N_hold <- N_cmb[, N_idx]
        sprintf("Subjects held %s %s \n", N_hold[1], N_hold[2])

        # Select data for the remaining subjects for the robust correlation.
        pbmc_exp_data_select <- pbmc_exp_data[, setdiff(colnames(pbmc_exp_data), N_hold)]
        rspo3_soma_select <- rspo3_soma[!rspo3_soma$Sample.ID %in% N_hold,]
        d28_titre_select <-  d28_titre[!d28_titre$Sample.ID %in% N_hold,]

        corel <- data.frame(corel = cor(t(pbmc_exp_data_select), rspo3_soma_select$RSPO3, method = "spearman"))
        corel <- rownames_to_column(corel,var = "genes")       

        # Load top K genes from the robust correlation results. 
        K = 30        
        gene.sig <- corel[order(corel$corel, decreasing = T),][1:K,"genes"]
        
        # Calcualte average score per hold out subject for the top k seleted genes 
        # from the pre computed z-score matrix. 
        subj_avg_score <- data.frame(score = colMeans(pbmc_exp_data_z[gene.sig, N_hold ]))
        subj_avg_score <- rownames_to_column(subj_avg_score, var = "Sample.ID")

        # Consolidate results.
        val_result <- rbind(val_result, data.frame(iter = N_idx, subject = subj_avg_score[1, "Sample.ID"], avg_score = subj_avg_score[1, "score"]))
        val_result <- rbind(val_result, data.frame(iter = N_idx, subject = subj_avg_score[2, "Sample.ID"], avg_score = subj_avg_score[2, "score"]))
}


stop()

write.table(val_result, file.path("RESULTS", "rspo3", "val_result_all_with_pre_z_norm.txt"), sep = "\t", row.names = F)

