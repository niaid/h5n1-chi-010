library(tidyverse)
library(stringr)

# FUNCTIONS
# Calculate score from genes x samples matrix as average z-score
get_score <- function(x) {
  x = t(scale(t(x)))
  return (colMeans(x, na.rm=T))
}

cor_z_titre <- function(data_adj, exp_data, flag, POS = FALSE){
        result_df <- data.frame(matrix(nrow = 0, ncol = 3))
        # Get the gene list above a certain cut-off.
        for(cut_off in seq(0.1, 0.9, by=0.1)){
                if(POS == FALSE){
                        # Select both positively and negatively correlated genes.
                        gene_cut <- data_adj[data_adj[[flag]] > cut_off | data_adj[[flag]] < -cut_off, "gene"]
                }else if(POS == TRUE){
                        # Select only positively correlated genes.
                        gene_cut <- data_adj[data_adj[[flag]] > cut_off, "gene"]
                }
                no_genes <- length(gene_cut)
                if(no_genes > 1 ){
                        # Extract expression data for the short listed genes.
                        gene_cut_exp <- exp_data[gene_cut,]

                        # Calculate the z-score for each subject using expression data for the short listed genes.
                        subj_z_score <- data.frame(score = get_score(gene_cut_exp))

                        # Calculate the correlation between the z-scores and the titre response. 
                        corel <- cor(d28_titre$A.Indonesia, subj_z_score$score)
                        print(corel) 
                        result_df <- rbind(result_df, data.frame(cut_off, no_genes, corel ))
                }else{
                        result_df <- rbind(result_df, data.frame(cut_off, no_genes = 0, corel = NA ))
                }
        }
        return(result_df)
}

cor_z_titre_flip <- function(data_adj, exp_data, flag){
        result_df <- data.frame(matrix(nrow = 0, ncol = 3))
        # Get the gene list above a certain cut-off.
        for(cut_off in seq(0.1, 0.9, by=0.1)){
                # Select both positively and negatively correlated genes separately. 
                gene_cut_pos <- data_adj[data_adj[[flag]] > cut_off, "gene"]
                gene_cut_neg <- data_adj[data_adj[[flag]] < -cut_off, "gene"]

                no_genes <- length(gene_cut_pos) + length(gene_cut_neg)
                if(no_genes > 1 ){
                        # Extract expression data for the short listed genes.
                        if(length(gene_cut_pos) != 0  & length(gene_cut_neg) != 0 ){
                                gene_cut_exp_neg <- exp_data[gene_cut_neg, ] * -1
                                gene_cut_exp_pos <- exp_data[gene_cut_pos,]
                                gene_cut_exp <- rbind(gene_cut_exp_neg, gene_cut_exp_pos) 
                        }else if(length(gene_cut_pos) != 0){
                                gene_cut_exp <- exp_data[gene_cut_pos, ]
                        }else if(length(gene_cut_neg) != 0){
                                gene_cut_exp <- exp_data[gene_cut_neg, ] * -1
                        }
                       
                        # Calculate the z-score for each subject using expression data for the short listed genes.
                        subj_z_score <- data.frame(score = get_score(gene_cut_exp))

                        # Calculate the correlation between the z-scores and the titre response. 
                        corel <- cor(d28_titre$A.Indonesia, subj_z_score$score)
                        print(corel) 
                        result_df <- rbind(result_df, data.frame(cut_off, no_genes, corel ))
                }else{
                        result_df <- rbind(result_df, data.frame(cut_off, no_genes = 0, corel = NA ))
                }
        }
        return(result_df)
}

fig_z_titre <- function(result_df, plot_name){
        scat_plot <- ggplot(result_df, aes(x=cut_off, y=corel)) +
                     geom_point() +
                     geom_text(aes(label=no_genes),hjust=0, vjust=0) +
                     geom_smooth() +
                     scale_x_continuous(name ="gene-set cut off", breaks=seq(0.1, 0.9, by =0.1)) +
                     labs(y = "r (z-score vs titre)")
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
pbmc_adj <- read.table(file.path("RESULTS", "SOMAscan", "pbmc_ADJ_RSPO3_bySex.txt"), header = TRUE, sep = "\t")
pax_adj <- read.table(file.path("RESULTS", "SOMAscan", "pax_ADJ_RSPO3_bySex.txt"), header = TRUE, sep = "\t")

# Load gene expression data.
pbmc_exp_data <- t(readRDS(file.path("New","190703", "DATA", "PBMC_dat0.in_isv.0.7_8144g_170216.rds")))
pax_exp_data <- t(readRDS(file.path("New","190703", "DATA", "PAX_dat0.in_pbmc.isv.0.7_8144g_170622.rds")))

# Select data only for the adjuvanted subjects.
pbmc_exp_data <- pbmc_exp_data[, subject]
colnames(pax_exp_data) <-  str_replace(colnames(pax_exp_data), "-", "_")
pax_exp_data <- pax_exp_data[, subject]

# Load titre response data.
titre <- read.table(file.path("DATA_ORIGINAL", "Titres", "titre.txt"), header = TRUE, sep = "\t")
titre$TimePt <- tolower(titre$TimePt)
d28_titre <-  titre[titre$TimePt %in% c("d28", "d29", "d31") , c("Sample.ID","A.Indonesia")]
d28_titre$Sample.ID <- paste("H5N1_", str_pad(d28_titre$Sample.ID, 3, pad = "0"), sep = "")
d28_titre <- d28_titre[d28_titre$Sample.ID %in% subject, ]
subject == d28_titre$Sample.ID # Spot check
d28_titre$A.Indonesia <-  as.numeric(str_replace(d28_titre$A.Indonesia, "<", ""))

# Correlation for both positively and negatively correlated genes. 
pbmc_res_r <- cor_z_titre(pbmc_adj, pbmc_exp_data, "r")
fig_z_titre(pbmc_res_r, "pbmc_z_titre_r.png")

pbmc_res_pr <- cor_z_titre(pbmc_adj, pbmc_exp_data, "p_r")
fig_z_titre(pbmc_res_pr, "pbmc_z_titre_pr.png")

pax_res_r <- cor_z_titre(pax_adj, pax_exp_data, "r")
fig_z_titre(pax_res_r, "pax_z_titre_r.png")

pax_res_pr <- cor_z_titre(pax_adj, pax_exp_data, "p_r")
fig_z_titre(pax_res_pr, "pax_z_titre_pr.png")

# Correlation for only positively correlated genes. 
pbmc_res_r <- cor_z_titre(pbmc_adj, pbmc_exp_data, "r", POS = TRUE)
fig_z_titre(pbmc_res_r, "pbmc_z_titre_pos_r.png")

pbmc_res_pr <- cor_z_titre(pbmc_adj, pbmc_exp_data, "p_r", POS = TRUE)
fig_z_titre(pbmc_res_pr, "pbmc_z_titre_pos_pr.png")

pax_res_r <- cor_z_titre(pax_adj, pax_exp_data, "r", POS = TRUE)
fig_z_titre(pax_res_r, "pax_z_titre_pos_r.png")

pax_res_pr <- cor_z_titre(pax_adj, pax_exp_data, "p_r", POS = TRUE)
fig_z_titre(pax_res_pr, "pax_z_titre_pos_pr.png")

# Correlation for both positively and flipped-negatively correlated genes.
pbmc_res_r_flip <- cor_z_titre_flip(pbmc_adj, pbmc_exp_data, "r")
fig_z_titre(pbmc_res_r_flip, "pbmc_z_titre_r_flip.png")

pbmc_res_pr_flip <- cor_z_titre_flip(pbmc_adj, pbmc_exp_data, "p_r")
fig_z_titre(pbmc_res_pr_flip, "pbmc_z_titre_pr_flip.png")

pax_res_r_flip <- cor_z_titre_flip(pax_adj, pax_exp_data, "r")
fig_z_titre(pax_res_r_flip, "pax_z_titre_r_flip.png")

pax_res_pr_flip <- cor_z_titre_flip(pax_adj, pax_exp_data, "p_r")
fig_z_titre(pax_res_pr_flip, "pax_z_titre_pr_flip.png")


stop()


##PBMC
result_df <- data.frame(matrix(nrow = 0, ncol = 3))
# Get the gene list above a certain cut-off.
for(cut_off in seq(0.1, 0.9, by=0.1)){
        #         gene_cut <- pbmc_adj[pbmc_adj$p_r > cut_off | pbmc_adj$p_r < -cut_off, "gene"]
        gene_cut <- pbmc_adj[pbmc_adj$r > cut_off, "gene"]
        no_genes <- length(gene_cut)
        if(no_genes != 0 & length(gene_cut) > 1 ){
                # Extract expression data for the short listed genes.
                gene_cut_exp <- pbmc_exp_data[gene_cut,]

                # Calculate the z-score for each subject using expression data for the short listed genes.
                subj_z_score <- data.frame(score = get_score(gene_cut_exp))

                # Calculate the correlation between the z-scores and the titre response. 
                corel <- cor(d28_titre$A.Indonesia, subj_z_score$score)
                print(corel) 
                result_df <- rbind(result_df, data.frame(cut_off, no_genes, corel ))
        }else{
                result_df <- rbind(result_df, data.frame(cut_off, no_genes = 0, corel = NA ))
        }
}

scat_plot <- ggplot(result_df, aes(x=cut_off, y=corel)) +
             geom_point() +
             geom_text(aes(label=no_genes),hjust=0, vjust=0) +
             geom_smooth() +
             scale_x_continuous(name ="gene-set cut off", breaks=seq(0.1, 0.9, by =0.1)) +
             labs(y = "r (z-score vs titre)")
ggsave(file.path("FIGURES","rspo3", "pbmc_z_titre_pos_r.png" ), plot = scat_plot, device = "png", width = 6, height = 4)

##PAX
result_df <- data.frame(matrix(nrow = 0, ncol = 3))
# Get the gene list above a certain cut-off.
for(cut_off in seq(0.1, 0.9, by=0.1)){
        #         gene_cut <- pax_adj[pax_adj$p_r > cut_off | pax_adj$p_r < -cut_off, "gene"]
        gene_cut <- pax_adj[pax_adj$r > cut_off, "gene"]
        no_genes <- length(gene_cut)
        if(no_genes != 0 & length(gene_cut) > 1 ){
                # Extract expression data for the short listed genes.
                gene_cut_exp <- pax_exp_data[gene_cut,]

                # Calculate the z-score for each subject using expression data for the short listed genes.
                subj_z_score <- data.frame(score = get_score(gene_cut_exp))

                # Calculate the correlation between the z-scores and the titre response. 
                corel <- cor(d28_titre$A.Indonesia, subj_z_score$score)
                print(corel) 
                result_df <- rbind(result_df, data.frame(cut_off, no_genes, corel ))
        }else{
                result_df <- rbind(result_df, data.frame(cut_off, no_genes = 0, corel = NA ))
        }
}

scat_plot <- ggplot(result_df, aes(x=cut_off, y=corel)) +
             geom_point() +
             geom_text(aes(label=no_genes),hjust=0, vjust=0) +
             geom_smooth() +
             scale_x_continuous(name ="gene-set cut off", breaks=seq(0.1, 0.9, by =0.1)) +
             labs(y = "r (z-score vs titre)")
ggsave(file.path("FIGURES","rspo3", "pax_z_titre__r.png" ), plot = scat_plot, device = "png", width = 6, height = 4)

