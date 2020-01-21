# PURPOSE: To calculate robust correlation between RSPO3 (soma scan expression) and gene expression
# for adjuvanted subjects.

library(pROC)
library(tidyverse)
library(data.table)
library(stringr)

# load gene expression data
# fn.ge = file.path("DATA_PROCESSED", "Baseline", "CHI_GE_matrix_gene.txt")
# dat = fread(fn.ge, data.table = F) %>% 
#   tibble::remove_rownames() %>% tibble::column_to_rownames("gene") %>% 
#   data.matrix()


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


# fn.si = file.path("DATA_PROCESSED", "Baseline", "CHI_sample_info_2.txt")
# info = fread(fn.si) %>% 
#   mutate(subject = as.character(subject))
# 
# fn.poB = file.path("DATA_PROCESSED", "Baseline", "flow_percent_of_B_filtered.txt")
# flow.poB = fread(fn.poB) %>% tibble::column_to_rownames("sample") %>% 
#   data.matrix()
# 
# fn.info = file.path("DATA_PROCESSED", "Baseline", "flow_sample_info_filtered.txt")
# flow.info = fread(fn.info) %>% 
#   mutate(subject = as.character(subject))
# 
# flow.df = flow.info %>% dplyr::select(sample) %>% add_column(CD38hi = flow.poB[,"Gate3"])
# 
# info = info %>% left_join(flow.df, by="sample")
# 
# 
# si = with(info, time == 0 & Response %in% c("low","high") & !is.na(CD38hi))
# sum(si)
# dat = dat[,si]
# info = info[si,] %>% mutate(Response = factor(Response, levels=c("low","high")))


# calculate robust regression removing ns samples --------
ns = 2
# cmb = combn(1:nrow(info), ns)
cmb = combn(1:length(subject), ns)
# cc.rob = matrix(NA,nrow=nrow(dat),ncol=ncol(cmb))
cc.rob = matrix(NA,nrow=nrow(pbmc_exp_data),ncol=ncol(cmb))
for (ic in 1:ncol(cmb)) {
  cat(ic," ")
  ridx = cmb[,ic]
  #   cc.rob[,ic] = cor(t(dat[,-ridx]), info$CD38hi[-ridx], method="spearman", use="pairwise.complete.obs")
  cc.rob[,ic] = cor(t(pbmc_exp_data[,-ridx]), rspo3_soma$RSPO3[-ridx], method="spearman", use="pairwise.complete.obs")
}
# rownames(cc.rob) = rownames(dat)
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

# calculate AUC for each gene
# auc.one=matrix(nrow=nrow(dat), ncol=3)
# for (k in 1:nrow(dat)) {
#   X = dat[k,]
#   Y = info$Response
#   auc.one[k, c(2,1,3)] = ci.auc(Y,X,direction="<", quiet=T)
# }
# rownames(auc.one) = rownames(dat)
# colnames(auc.one) = c("auc.gene", "auc.ci95.min", "auc.ci95.max")

# output
# df.out = data.frame(cor.ntop20 = cc.rob.ntop, 
#                     cor.mean = cc.rob.mean, 
#                     cor.sd = cc.rob.sd, 
#                     cor.meadian = cc.rob.median, 
#                     cor.mean.sd.ratio = cc.rob.1cv, 
#                     as.data.frame(auc.one)) %>%
#   tibble::rownames_to_column("gene")

df.out = data.frame(cor.ntop20 = cc.rob.ntop, 
                    cor.mean = cc.rob.mean, 
                    cor.sd = cc.rob.sd, 
                    cor.meadian = cc.rob.median, 
                    cor.mean.sd.ratio = cc.rob.1cv
                    ) %>%
  tibble::rownames_to_column("gene")

fn.cor = file.path("RESULTS", "rspo3","robust_corel_all_genes.txt")
fwrite(df.out, fn.cor, sep="\t", quote=T)

