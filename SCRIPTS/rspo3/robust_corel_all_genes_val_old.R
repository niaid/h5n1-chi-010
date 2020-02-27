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

robust_corel <- function(exp_data, rspo3_soma){
        
        # Calculate robust correlation between RSPO3 and gene expression for adjuvanted subjects on day 0.
        ns = 2
        cmb = combn(1:length(colnames(exp_data)), ns)
        cc.rob = matrix(NA,nrow=nrow(exp_data),ncol=ncol(cmb))
        for (ic in 1:ncol(cmb)) {
          cat(ic," ")
          ridx = cmb[,ic]
          cc.rob[,ic] = cor(t(exp_data[,-ridx]), rspo3_soma$RSPO3[-ridx], method="spearman", use="pairwise.complete.obs")
        }
        rownames(cc.rob) = rownames(exp_data)
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
                            ) #%>%
        #           tibble::rownames_to_column("gene")
        
        return(df.out)
}

top_n_genes <- function(corel_result, exp_data, titre){
        ng.max = 200
        result_df = data.frame(matrix(nrow=0, ncol=2))

        for(top in 2:ng.max){
                gene.sig <- load_sig(corel_result, "cor.mean.sd.ratio", ntop=top)
                gene_cut_exp <- exp_data[gene.sig,]


                # Day 0 correlation.
                # Calculate the z-score for each subject using expression data for the short listed genes.
                subj_z_score <- data.frame(score = get_score(gene_cut_exp))
                subj_z_score <- rownames_to_column(subj_z_score, var = "Sample.ID")
                df_score_titre <- merge(subj_z_score, titre) # To make sure that the samples are in the same order.

                corel_day0 <- cor(df_score_titre$A.Indonesia, df_score_titre$score, method = "spearman")

                result_df <- rbind(result_df, data.frame(top_ng = top, spear_day0 = corel_day0))
        }
        return(result_df)
}

# Load subject info and select adjuvanted subjects
subject <- read.table(file.path("DATA_ORIGINAL","Clinical", "clinical_info_adj.txt"), header = TRUE, sep = "\t")
subject <- subject[subject$Adjuvant == "Adj", c("Subject.ID")]
subject <- str_replace(subject, "-", "_")
# Remove subject H5N1_010
subject <- subject[!subject == "H5N1_010"]

# Load gene expression data on day 0.
pbmc_exp_data <- t(readRDS(file.path("New","190703", "DATA", "PBMC_dat0.in_isv.0.7_8144g_170216.rds")))
# Select data only for the adjuvanted subjects.
pbmc_exp_data <- pbmc_exp_data[, subject]

# Load SOMAScan data for RSPO3.
rspo3_soma <- read.table(file.path("New", "190703", "DATA", "SOMA.txt"), header = TRUE, sep = "\t")
rspo3_soma <- rspo3_soma[c("donor", "RSPO3")]
rspo3_soma$donor <- str_replace(rspo3_soma$donor, "-", "_")
rspo3_soma <-  rspo3_soma[rspo3_soma$donor %in% subject,]

# Match gene expression subject with SOMAScan subjects.
pbmc_exp_data_rspo3_match <-  pbmc_exp_data[, rspo3_soma$donor]

# Load titre response data.
titre <- read.table(file.path("DATA_ORIGINAL", "Titres", "titre.txt"), header = TRUE, sep = "\t")
titre$TimePt <- tolower(titre$TimePt)
d28_titre <-  titre[titre$TimePt %in% c("d28", "d29", "d31") , c("Sample.ID","A.Indonesia")]
d28_titre$Sample.ID <- paste("H5N1_", str_pad(d28_titre$Sample.ID, 3, pad = "0"), sep = "")
d28_titre <- d28_titre[d28_titre$Sample.ID %in% subject, ]
subject == d28_titre$Sample.ID # Spot check
d28_titre$A.Indonesia <-  as.numeric(str_replace(d28_titre$A.Indonesia, "<", ""))


ns_top <- 2
cmb_top <- combn(1:length(colnames(pbmc_exp_data_rspo3_match)), ns_top )
val_result <- data.frame(matrix(nrow = 0, ncol=4))
for(i_top in 1:ncol(cmb_top)){
        subj_hold <- cmb_top[, i_top]
        corel_result <- robust_corel(pbmc_exp_data_rspo3_match[,-subj_hold], rspo3_soma[-subj_hold,])

        top_genes <- top_n_genes(corel_result, pbmc_exp_data_rspo3_match[, -subj_hold], d28_titre)
        # Max correlation on day0.
        max_0 <-  top_genes[top_genes$spear_day0 == max(top_genes$spear_day0), "top_ng"]
        # Genes at max correlation on day0.
        max_genes_0 <- load_sig(corel_result, "cor.mean.sd.ratio", ntop=max_0)
        max_genes_0 <- data.frame(genes = max_genes_0)
        max_genes_0_exp <- pbmc_exp_data_rspo3_match[as.character(max_genes_0$genes), subj_hold]
        max_genes_0_z_score <- get_score(max_genes_0_exp)
        val_result <- rbind(val_result, data.frame(iter = i_top, subject = names(max_genes_0_z_score[1]), z_score = max_genes_0_z_score[[1]], top_n =  length(max_0)))
        val_result <- rbind(val_result, data.frame(iter = i_top, subject = names(max_genes_0_z_score[2]), z_score = max_genes_0_z_score[[2]], top_n =  length(max_0)))
}

stop()
write.table(val_result, file.path("RESULTS", "rspo3", "validation_result.txt"), sep = "\t", row.names = F)
