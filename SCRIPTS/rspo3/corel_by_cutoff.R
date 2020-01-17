library(tidyverse)
library(stringr)
library(ggrepel) # Install in the container.

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
                        subj_z_score <- rownames_to_column(subj_z_score, var = "Sample.ID")
                        df_score_titre <- merge(subj_z_score, d28_titre) # To make sure that the samples are in the same order.

                        # Calculate the correlation between the z-scores and the titre response.
                        print("Calculating pearson correlation for the non-flipped set.")
                        corel <- cor(df_score_titre$A.Indonesia, df_score_titre$score)
                        print(corel) 
                        result_df <- rbind(result_df, data.frame(cut_off, no_genes, corel ))
                }else{
                        result_df <- rbind(result_df, data.frame(cut_off, no_genes = 0, corel = NA ))
                }
        }
        return(result_df)
}

cor_z_titre_flip <- function(data_adj, exp_data, flag, pbmc = TRUE, lb = 0.1, spear = FALSE, subj_drop = c(), li = FALSE, gene_stab = c()){
        result_df <- data.frame(matrix(nrow = 0, ncol = 3))
        # Get the gene list above a certain cut-off.
        for(cut_off in seq(lb, 0.9, by=0.1)){
                # Select both positively and negatively correlated genes separately. 
                gene_cut_pos <- data_adj[data_adj[[flag]] > cut_off, "gene"]
                gene_cut_neg <- data_adj[data_adj[[flag]] < -cut_off, "gene"]
                
                no_genes <- length(gene_cut_pos) + length(gene_cut_neg)
                if(no_genes > 1 ){
                        # Extract expression data for the short listed genes.
                        if(class(gene_stab) != "data.frame"){
                                if(length(gene_cut_pos) != 0  & length(gene_cut_neg) != 0 ){
                                        gene_cut_exp_neg <- exp_data[gene_cut_neg, ] * -1
                                        gene_cut_exp_pos <- exp_data[gene_cut_pos,]
                                        gene_cut_exp <- rbind(gene_cut_exp_neg, gene_cut_exp_pos) 
                                }else if(length(gene_cut_pos) != 0){
                                        gene_cut_exp <- exp_data[gene_cut_pos, ]
                                }else if(length(gene_cut_neg) != 0){
                                        gene_cut_exp <- exp_data[gene_cut_neg, ] * -1
                                } 
                        }else if(class(gene_stab) == "data.frame"){
                        # Check if genes are above a cut-off in TSM.
                                if(length(gene_cut_pos) != 0  & length(gene_cut_neg) != 0 ){
                                        # Note: for this particular analysis i knew that we have genes for both positive
                                        # and negative cutoffs that is why I skipped the other two if blocks.
                                        gene_names <- data.frame(gene = c(as.character(gene_cut_pos), as.character(gene_cut_neg)))
                                        df_tsm <- merge(gene_stab, gene_names)
                                        df_tsm <- df_tsm[df_tsm$ISV > 0.5,] 
                                        tsm_genes <- as.character(df_tsm$gene)
                                        gene_cut_exp_neg <- exp_data[intersect(tsm_genes, gene_cut_neg), ] * -1
                                        gene_cut_exp_pos <- exp_data[intersect(tsm_genes, gene_cut_pos),]
                                        gene_cut_exp <- rbind(gene_cut_exp_neg, gene_cut_exp_pos) 
                                       
                                        no_genes <- nrow(gene_cut_exp) # update the reduced number of genes.
                                }                         
                        }
                       
                        # Check if there are subjects to be dropped from the analysis.
                        if(length(subj_drop) != 0){
                                 gene_cut_exp <- gene_cut_exp[, !(colnames(gene_cut_exp) %in% subj_drop)]
                        }
                      
                        # Calculate the z-score for each subject using expression data for the short listed genes.
                        subj_z_score <- data.frame(score = get_score(gene_cut_exp))
                        subj_z_score <- rownames_to_column(subj_z_score, var = "Sample.ID")
                        df_score_titre <- merge(subj_z_score, d28_titre) # To make sure that the samples are in the same order.

                        # Generate a scatter plot between z-score and titer value for each cutoff.
                        # Only when executing pearson correlation and without dropping subjects; it's to avoid
                        # redundant figure generations. 
                        if(pbmc == TRUE & length(subj_drop) == 0 & spear == FALSE & li == FALSE & length(gene_stab) == 0){
                                plot_name <- paste("pbmc_", flag, "_", "score_titre_", cut_off, ".png", sep = "")
                                score_titre_scat(df_score_titre, plot_name)
                        }else if(pbmc == FALSE & length(subj_drop) == 0 & spear == FALSE & li == FALSE & length(gene_stab) == 0){
                        
                                plot_name <- paste("pax_", flag, "_", "score_titre_", cut_off, ".png", sep = "")
                                score_titre_scat(df_score_titre, plot_name)
                        }


                        # Calculate the correlation between the z-scores and the titre response.
                        if(spear == FALSE){
                                print("Calculating pearson correlation for the flipped set.")
                                # d28_titre is the global variable that is used here. 
                                corel <- cor(df_score_titre$A.Indonesia, df_score_titre$score)
                                print(corel) 
                                result_df <- rbind(result_df, data.frame(cut_off, no_genes, corel ))
                        }else if(spear == TRUE){
                                print("Calculating spearman ranks correlation for the flipped set.")
                                corel <- cor(df_score_titre$A.Indonesia, df_score_titre$score, method = "spearman")
                                print(corel) 
                                result_df <- rbind(result_df, data.frame(cut_off, no_genes, corel ))
                        }
                }else{
                        result_df <- rbind(result_df, data.frame(cut_off, no_genes = 0, corel = NA ))
                }
        }
        return(list(result_df, df_score_titre))
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
                     scale_x_continuous(name ="gene-set cut off", breaks=seq(0.1, 0.9, by =0.1)) +
                     labs(y = y_lab)
        ggsave(file.path("FIGURES","rspo3", plot_name), plot = scat_plot, device = "png", width = 6, height = 4)
}

score_titre_scat <- function(df_score_titre, plot_name){
        scat_plot <- ggplot(df_score_titre, aes(x=score, y=A.Indonesia)) +
                     geom_point(color = "red") +
                     geom_text_repel(size = 2, aes(label=Sample.ID)) +
                     labs(x = "z-score", y = "D28/A.Indonesia")
        ggsave(file.path("FIGURES", "rspo3", plot_name), plot = scat_plot, device = "png", width = 5, height = 4)
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


# Pearson correlation for both positively and flipped-negatively correlated genes.

# Note: Only concentrate on PBMCs for now.

# Approach 3. TSM (Temporal Stability Metric)
gene_stability <- read.table(file.path("SCRIPTS", "rspo3", "CHI_genes_stability.txt"), header = TRUE, sep = "\t")
pbmc_res_r_flip_tsm <- cor_z_titre_flip(pbmc_adj, pbmc_exp_data, "r", gene_stab = gene_stability)
fig_z_titre(pbmc_res_r_flip_tsm[[1]], "pbmc_z_titre_r_flip_tsm.png")



stop()
# Approach 2. Load GSEA module data.
li <- read.table(file.path("SCRIPTS", "rspo3", "rspo3_gsea.txt"), header = TRUE, sep = "\t")
for(i in 1:dim(li)[1]){
        li_gene_list <- strsplit(as.character(li[i,2]), ",")[[1]]
        df_gene <- data.frame(gene=li_gene_list)
        pbmc_adj_btm <- merge(pbmc_adj, df_gene)
        pbmc_res_r_flip <- cor_z_titre_flip(pbmc_adj_btm, pbmc_exp_data, "r", lb = 0, li = TRUE)
        plot_name <- paste("pbmc_z_titre_r_flip_",li$btm[i], ".png", sep = "")
        fig_z_titre(pbmc_res_r_flip[[1]], plot_name)

}

stop()

# Approach 1. 
# 1
pbmc_res_r_flip <- cor_z_titre_flip(pbmc_adj, pbmc_exp_data, "r")
fig_z_titre(pbmc_res_r_flip[[1]], "pbmc_z_titre_r_flip.png")
# Subject Drop
pbmc_res_r_flip <- cor_z_titre_flip(pbmc_adj, pbmc_exp_data, "r", subj_drop = c("H5N1_038"))
fig_z_titre(pbmc_res_r_flip[[1]], "pbmc_z_titre_r_flip_sd.png")

stop()

# 2 
pbmc_res_pr_flip <- cor_z_titre_flip(pbmc_adj, pbmc_exp_data, "p_r")
fig_z_titre(pbmc_res_pr_flip[[1]], "pbmc_z_titre_pr_flip.png")
# Subject Drop
pbmc_res_pr_flip <- cor_z_titre_flip(pbmc_adj, pbmc_exp_data, "p_r",  subj_drop = c("H5N1_038"))
fig_z_titre(pbmc_res_pr_flip[[1]], "pbmc_z_titre_pr_flip_sd.png")

# 3
pax_res_r_flip <- cor_z_titre_flip(pax_adj, pax_exp_data, "r", pbmc = FALSE)
fig_z_titre(pax_res_r_flip[[1]], "pax_z_titre_r_flip.png")
# Subject Drop
pax_res_r_flip <- cor_z_titre_flip(pax_adj, pax_exp_data, "r", subj_drop = c("H5N1_038"))
fig_z_titre(pax_res_r_flip[[1]], "pax_z_titre_r_flip_sd.png")

# 4
pax_res_pr_flip <- cor_z_titre_flip(pax_adj, pax_exp_data, "p_r", pbmc = FALSE)
fig_z_titre(pax_res_pr_flip[[1]], "pax_z_titre_pr_flip.png")
# Subject Drop
pax_res_pr_flip <- cor_z_titre_flip(pax_adj, pax_exp_data, "p_r", subj_drop = c("H5N1_038"))
fig_z_titre(pax_res_pr_flip[[1]], "pax_z_titre_pr_flip_sd.png")

# Spearman ranks correlation for both positively and flipped-negatively correlated genes.
pbmc_res_r_flip_spear <- cor_z_titre_flip(pbmc_adj, pbmc_exp_data, "r", spear = TRUE)
fig_z_titre(pbmc_res_r_flip_spear[[1]], "pbmc_z_titre_r_flip_spear.png", spear = TRUE)

pbmc_res_pr_flip_spear <- cor_z_titre_flip(pbmc_adj, pbmc_exp_data, "p_r", spear = TRUE)
fig_z_titre(pbmc_res_pr_flip_spear[[1]], "pbmc_z_titre_pr_flip_spear.png", spear = TRUE)

pax_res_r_flip_spear <- cor_z_titre_flip(pax_adj, pax_exp_data, "r", spear = TRUE)
fig_z_titre(pax_res_r_flip_spear[[1]], "pax_z_titre_r_flip_spear.png", spear = TRUE)

pax_res_pr_flip_spear <- cor_z_titre_flip(pax_adj, pax_exp_data, "p_r", spear = TRUE)
fig_z_titre(pax_res_pr_flip_spear[[1]], "pax_z_titre_pr_flip_spear.png", spear = TRUE)

# Correlation for both positively and negatively correlated genes. 
# pbmc_res_r <- cor_z_titre(pbmc_adj, pbmc_exp_data, "r")
# fig_z_titre(pbmc_res_r, "pbmc_z_titre_r.png")
# 
# pbmc_res_pr <- cor_z_titre(pbmc_adj, pbmc_exp_data, "p_r")
# fig_z_titre(pbmc_res_pr, "pbmc_z_titre_pr.png")
# 
# pax_res_r <- cor_z_titre(pax_adj, pax_exp_data, "r")
# fig_z_titre(pax_res_r, "pax_z_titre_r.png")
# 
# pax_res_pr <- cor_z_titre(pax_adj, pax_exp_data, "p_r")
# fig_z_titre(pax_res_pr, "pax_z_titre_pr.png")

# Correlation for only positively correlated genes. 
# pbmc_res_r <- cor_z_titre(pbmc_adj, pbmc_exp_data, "r", POS = TRUE)
# fig_z_titre(pbmc_res_r, "pbmc_z_titre_pos_r.png")
# 
# pbmc_res_pr <- cor_z_titre(pbmc_adj, pbmc_exp_data, "p_r", POS = TRUE)
# fig_z_titre(pbmc_res_pr, "pbmc_z_titre_pos_pr.png")
# 
# pax_res_r <- cor_z_titre(pax_adj, pax_exp_data, "r", POS = TRUE)
# fig_z_titre(pax_res_r, "pax_z_titre_pos_r.png")
# 
# pax_res_pr <- cor_z_titre(pax_adj, pax_exp_data, "p_r", POS = TRUE)
# fig_z_titre(pax_res_pr, "pax_z_titre_pos_pr.png")
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

