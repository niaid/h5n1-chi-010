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

print(date())
# fit <- eNetXplorer(x = X, y = y, family = "gaussian", alpha = seq(0, 1, by = 0.1), nlambda.ext=100, n_run = 1000, n_perm_null = 100, seed = 123, cor_method = "spearman")
fit <- eNetXplorer(x = X, y = y, family = "gaussian", alpha = seq(0, 1, by = 0.1), seed = 123, cor_method = "pearson")
fit_summary <- summary(fit)
saveRDS(fit, file.path("RESULTS", "rspo3", "genes_rspo3_enet_fit_v4.rds"))
print(date())

stop()
# Plots
# Model performance.
png(filename = file.path("FIGURES","rspo3", "model_performance_v3.png"), width = 6, height = 4, units = "in", pointsize = 12, res = 300)
plot(fit, plot.type="summary")
dev.off()

# Top features.
png(filename = file.path("FIGURES","rspo3", "top_features_v3.png"), width = 6, height = 4, units = "in", pointsize = 12, res = 300)
plot(fit, alpha.index = which.max(fit$model_QF_est), plot.type = "featureCaterpillar",stat = c("coef"))
dev.off()

# Selectd features coefficients.
png(filename = file.path("FIGURES","rspo3", "feature_coeff_v3.png"), width = 6, height = 4, units = "in", res = 300, pointsize = 8)
plot(fit, alpha.index = which.max(fit$model_QF_est), plot.type = "featureHeatmap", stat = c("coef"), notecex = 1.5)
dev.off()


