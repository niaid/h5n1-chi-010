library(eNetXplorer)
library(tidyverse)
library(stringr)
library(ggrepel)

# FUNCTIONS
min_max <- function(x)
{
    return((x- min(x)) /(max(x)-min(x)))
}


# MAIN

# Load subject info and select adjuvanted subjects
subject <- read.table(file.path("DATA_ORIGINAL","Clinical", "clinical_info_adj.txt"), header = TRUE, sep = "\t")
subject <- subject[subject$Adjuvant == "Adj", c("Subject.ID")]
subject <- str_replace(subject, "-", "_")
# Remove subject H5N1_010
subject <- subject[!subject == "H5N1_010"]

# Load gene expression data.
pbmc_exp_data <- t(readRDS(file.path("New","190703", "DATA", "PBMC_dat0.in_isv.0.7_8144g_170216.rds")))

# Select data only for the adjuvanted subjects.
pbmc_exp_data <- pbmc_exp_data[, subject]

#  TSM (Temporal Stability Metric)
gene_stability <- read.table(file.path("SCRIPTS", "rspo3", "CHI_genes_stability.txt"), header = TRUE, sep = "\t")
gene_stability <- gene_stability[gene_stability$ISV > 0.5,]

# Select genes according to TSM.
pbmc_exp_data <- pbmc_exp_data[rownames(pbmc_exp_data) %in% gene_stability$gene,]

# Load SOMAScan data for RSPO3.
rspo3_soma <- read.table(file.path("New", "190703", "DATA", "SOMA.txt"), header = TRUE, sep = "\t")
rspo3_soma <- rspo3_soma[c("donor", "RSPO3")]
rspo3_soma$donor <- str_replace(rspo3_soma$donor, "-", "_")
rspo3_soma <-  rspo3_soma[rspo3_soma$donor %in% subject,]

# eNetXplorer
X <- t(pbmc_exp_data)
# X <- apply(X, 2, min_max)
y <- matrix(rspo3_soma$RSPO3, dimnames = list(rspo3_soma$donor, "RSPO3"))
# Check if the donors in X and y are the same.
X <- X[rownames(X) %in% rownames(y),]
y <- y[,"RSPO3"]
y <- as.numeric(as.character(y))

# Run eNet
print(date())
results_dir <- file.path("RESULTS", "rspo3")
enet_obj <- "new_enet_spearman.Robj"
eNetXplorer(x = X, y = y, family = "gaussian", alpha = seq(0, 1, by = 0.1), seed = 123, cor_method = "spearman", save_obj=TRUE, dest_dir= results_dir, dest_obj=enet_obj)
print(date())


