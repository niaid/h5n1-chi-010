library(pROC)
library(tidyverse)
library(data.table)

source(file.path("SCRIPTS", "functions", "load_sig.r"))
source(file.path("SCRIPTS", "functions", "get_score.r"))
fn.cd38.cor = file.path("RESULTS", "rspo3", "robust_corel_all_genes.txt")

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
d28_titre <-  titre[titre$TimePt %in% c("d28", "d29", "d31") , c("Sample.ID","A.Indonesia")]
d28_titre$Sample.ID <- paste("H5N1_", str_pad(d28_titre$Sample.ID, 3, pad = "0"), sep = "")
d28_titre <- d28_titre[d28_titre$Sample.ID %in% subject, ]
subject == d28_titre$Sample.ID # Spot check
d28_titre$A.Indonesia <-  as.numeric(str_replace(d28_titre$A.Indonesia, "<", ""))

# load gene expression data
# fn.ge = file.path("DATA_PROCESSED", "Baseline",  "CHI_GE_matrix_gene.txt")
# dat = fread(fn.ge, data.table = F) %>% 
#   tibble::remove_rownames() %>% tibble::column_to_rownames("gene") %>% 
#   data.matrix()
# 
# fn.si = file.path("DATA_PROCESSED", "Baseline", "CHI_sample_info_2.txt")
# info = fread(fn.si)

ng.max = 100
result_df = data.frame(matrix(nrow=0, ncol=3))

for(top in 2:ng.max){
        gene.sig <- load_sig(fn.cd38.cor, "cor.mean.sd.ratio", ntop=top)
        gene_cut_exp <- pbmc_exp_data[gene.sig,]

        # Day 100 correlation.
        # Calculate the z-score for each subject using expression data for the short listed genes.
        subj_z_score <- data.frame(score = get_score(gene_cut_exp))
        subj_z_score <- rownames_to_column(subj_z_score, var = "Sample.ID")
        df_score_titre <- merge(subj_z_score, d28_titre) # To make sure that the samples are in the same order.

        corel_day0 <- cor(df_score_titre$A.Indonesia, df_score_titre$score, method = "spearman")

        # Day 100 correlation.
        gene_cut_exp_day100 <- pbmc_day100[gene.sig,]
        # Calculate the z-score for each subject using expression data for the short listed genes.
        subj_z_score_day100 <- data.frame(score = get_score(gene_cut_exp_day100))
        subj_z_score_day100 <- rownames_to_column(subj_z_score_day100, var = "Sample.ID")
        df_score_titre_day100 <- merge(subj_z_score_day100, d28_titre) # To make sure that the samples are in the same order.

        corel_day100 <- cor(df_score_titre_day100$A.Indonesia, df_score_titre_day100$score, method = "spearman")

        result_df <- rbind(result_df, data.frame(top_ng = top, spear_day0 = corel_day0, spear_day100 = corel_day100))
}

plot_topn <- result_df %>%
        gather(key,value, spear_day0, spear_day100) %>%
        ggplot(aes(x=top_ng, y=value, colour=key)) +
        geom_line()
plot_name <- "topn_day0_day100.png"
ggsave(file.path("FIGURES","rspo3", plot_name), plot = plot_topn, device = "png", width = 6, height = 4)

stop()




# for (tp in c(0,-7,70)) {
# tp <- 0
#    cat(tp,"")
#   
#   day.pt = paste0("day ",tp)
# test signature
#   for(ng in 2:ng.max) {
#           ng <- 2
# load CD38 signature genes
#     gene.sig = load_sig(fn.cd38.cor, "cor.mean.sd.ratio", ntop=ng)
# gi = toupper(rownames(dat)) %in% toupper(gene.sig)
#     gi = toupper(rownames(pbmc_exp_data)) %in% toupper(gene.sig)
#     sum(gi)
#     
#     iday = info$time==tp
#     ihl = info$adjMFC_class[iday] %in% c(0,2)
#     X = get_score(dat[gi,iday])[ihl]
#     Y = info$adjMFC_class[iday][ihl]
#     
#     auc.ng = ci.auc(Y,X, direction="<", quiet=T) %>% matrix(nrow=1)
#     colnames(auc.ng) = c("auc.ci95.min", "AUC", "auc.ci95.max")
#     r.df = data.frame(day=day.pt, ng, as.data.frame(auc.ng))
#     
#     roc.df = rbind(roc.df, r.df)
#   }
# }


# roc.df = roc.df %>% 
#   mutate(day = factor(day, levels=c("day 0","day -7","day 70"))) %>% 
#   mutate(ngx = ifelse(day=="day -7", ng-0.2, ifelse(day=="day 70", ng+0.2, ng)))
# 
# ggplot(roc.df, aes(ng, AUC, group=day, col=day)) + 
# geom_linerange(aes(x=ngx, ymin=auc.ci95.min, ymax=auc.ci95.max)) +
#         geom_ribbon(aes(ymin=auc.ci95.min, ymax=auc.ci95.max, fill=day), alpha=0.2, col=NA) +
        #   geom_line(size=1) +
        #   geom_point(size=1, shape=21, fill="white", stroke=1) +
        #   scale_color_manual(values=c("grey65", "grey80","grey50")) +
        #   xlab("No. of genes") + ylab("AUC") + theme_bw() + 
        #   theme(legend.key = element_blank())
        # 
        # fn.fig = file.path(PROJECT_DIR, "figure_generation", 
        #                    sprintf("CHI_CD38.2-%d.genes_AUC_ci95",ng.max))
        # ggsave(paste0(fn.fig,".png"), w=5, h=4)
        # ggsave(paste0(fn.fig,".pdf"), w=5, h=4)
