# PURPOSE: To calculate robust correlation between RSPO3 (soma scan expression) and gene expression
# for adjuvanted subjects.

library(pROC)
library(tidyverse)
library(data.table)
library(stringr)

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

# Load SOMAScan data for RSPO3.
rspo3_soma <- read.table(file.path("New", "190703", "DATA", "SOMA.txt"), header = TRUE, sep = "\t")
rspo3_soma <- rspo3_soma[c("donor", "RSPO3")]
rspo3_soma$donor <- str_replace(rspo3_soma$donor, "-", "_")
rspo3_soma <-  rspo3_soma[rspo3_soma$donor %in% subject,]

# Match gene expression subject with SOMAScan subjects.
pbmc_exp_data <-  pbmc_exp_data[, rspo3_soma$donor]

# Gene stability.
gene_stability <- read.table(file.path("SCRIPTS", "rspo3", "CHI_genes_stability.txt"), header = TRUE, sep = "\t")
idx <- gene_stability$ISV >= 0.5
tsm_gene <- gene_stability[idx, ]$gene

# Fetch expression data for the TSM genes.
pbmc_exp_data <- pbmc_exp_data[rownames(pbmc_exp_data) %in% tsm_gene,]

# calculate robust regression removing ns samples --------
ns = 2
cmb = combn(1:length(colnames(pbmc_exp_data)), ns)
cc.rob = matrix(NA,nrow=nrow(pbmc_exp_data),ncol=ncol(cmb))
for (ic in 1:ncol(cmb)) {
  cat(ic," ")
  ridx = cmb[,ic]
  cc.rob[,ic] = cor(t(pbmc_exp_data[,-ridx]), rspo3_soma$RSPO3[-ridx], method="spearman", use="pairwise.complete.obs")
}
rownames(cc.rob) = rownames(pbmc_exp_data)
cc.rob.rank = apply(-cc.rob,2,rank)
ntop = 20
cc.rob.ntop = rowSums(cc.rob.rank<=ntop)

cc.rob.mean = apply(cc.rob,1,mean)
cc.rob.median = apply(cc.rob,1,median)
cc.rob.sd = apply(cc.rob,1,sd)
cc.rob.1cv = cc.rob.mean / cc.rob.sd
cc.rob.ntop.rank = rank(cc.rob.ntop)

cc.rob.1cv.ord = order(cc.rob.1cv,decreasing = T)
cc.rob.1cv.sort = sort(cc.rob.1cv,decreasing = T)
head(cc.rob.1cv.sort, 10)

df.out = data.frame(cor.ntop20 = cc.rob.ntop, 
                    cor.mean = cc.rob.mean, 
                    cor.sd = cc.rob.sd, 
                    cor.meadian = cc.rob.median, 
                    cor.mean.sd.ratio = cc.rob.1cv
                    ) %>%
  tibble::rownames_to_column("gene")

fn.cor = file.path("RESULTS", "rspo3","robust_corel_tsm_genes_2_test.txt")
fwrite(df.out, fn.cor, sep="\t", quote=T)


