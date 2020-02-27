library(tidyverse)
library(stringr)
library(ggrepel) # Install in the container.

# FUNCTIONS
# Calculate score from genes x samples matrix as average z-score
get_score <- function(x) {
  x = t(scale(t(x)))
  return (colMeans(x, na.rm=T))
}

fig_z_titre <- function(result_df, plot_name, spear = FALSE){
        if(spear == FALSE){
                y_lab <- "r (z-score vs titre)"
        } else if(spear == TRUE){
                y_lab <- "rho (z-score vs titre)"
        }
        scat_plot <- ggplot(result_df, aes(x=cut_off, y=corel)) +
                     geom_point() +
                     geom_text(aes(label=no_genes),hjust=0, vjust=0) +
                     geom_smooth() +
                     scale_x_continuous(name ="alpha", breaks= result_df$cut_off) +
                     labs(y = y_lab)
        ggsave(file.path("FIGURES","rspo3", plot_name), plot = scat_plot, device = "png", width = 6, height = 4)
}

# MAIN

# Load subject info and select adjuvanted subjects
subject <- read.table(file.path("DATA_ORIGINAL","Clinical", "clinical_info_adj.txt"), header = TRUE, sep = "\t")
subject <- subject[subject$Adjuvant == "Adj", c("Subject.ID")]
subject <- str_replace(subject, "-", "_")
# Remove subject H5N1_010
subject <- subject[!subject == "H5N1_010"]

# Load RSPO3 and gene expression correlation data.
# pbmc_adj <- read.table(file.path("RESULTS", "SOMAscan", "pbmc_ADJ_RSPO3_bySex.txt"), header = TRUE, sep = "\t")
# pax_adj <- read.table(file.path("RESULTS", "SOMAscan", "pax_ADJ_RSPO3_bySex.txt"), header = TRUE, sep = "\t")

# Load gene expression data.
pbmc_exp_data <- t(readRDS(file.path("New","190703", "DATA", "PBMC_dat0.in_isv.0.7_8144g_170216.rds")))

# Select data only for the adjuvanted subjects.
pbmc_exp_data <- pbmc_exp_data[, subject]

# Load titre response data.
titre <- read.table(file.path("DATA_ORIGINAL", "Titres", "titre.txt"), header = TRUE, sep = "\t")
titre$TimePt <- tolower(titre$TimePt)
d28_titre <-  titre[titre$TimePt %in% c("d28", "d29", "d31") , c("Sample.ID","A.Indonesia")]
d28_titre$Sample.ID <- paste("H5N1_", str_pad(d28_titre$Sample.ID, 3, pad = "0"), sep = "")
d28_titre <- d28_titre[d28_titre$Sample.ID %in% subject, ]
subject == d28_titre$Sample.ID # Spot check
d28_titre$A.Indonesia <-  as.numeric(str_replace(d28_titre$A.Indonesia, "<", ""))

# Load gene lists from enet alphas.
a0.2 <- readRDS(file.path("RESULTS", "rspo3", "pear", "a0.2_pear.rds"))
a1 <- readRDS(file.path("RESULTS", "rspo3", "pear", "a1_pear.rds"))

cut_off <- c(0.2, 1)
count <- 1
result_df <- data.frame(matrix(nrow = 0, ncol = 3))
for(genes in list(a0.2, a1)){
        gene_cut_exp <- pbmc_exp_data[genes,]
        subj_z_score <- data.frame(score = get_score(gene_cut_exp))
        subj_z_score <- rownames_to_column(subj_z_score, var = "Sample.ID")
        df_score_titre <- merge(subj_z_score, d28_titre) # To make sure that the samples are in the same order.
        corel <- cor(df_score_titre$A.Indonesia, df_score_titre$score, method = "spearman")
        result_df <- rbind(result_df, data.frame(cut_off = cut_off[count], no_genes = length(genes), corel ))
        count <- count + 1
}

fig_z_titre(result_df, "corel_enet_gene_titre.png", spear = TRUE)
# Calculate the z-score for each subject using expression data for the short listed genes.

stop()


