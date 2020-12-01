# PURPOSE: To compute a correlation between the z-score calculated per subject based on the top N proteins selected
# by the eNet explorer and the titre. For both our and Emory's cohort.  

source("SCRIPTS/0_initialize.r")

# FUNCTIONS
# Calculate score from genes x samples matrix as average z-score
get_score <- function(x) {
  x = t(scale(t(x)))
  return (colMeans(x, na.rm=T))
}

# CHI DATA
# Load Demographic Data
fn.demo = file.path(PROJECT_DIR, "DATA_ORIGINAL", "Clinical", "clinical_info_adj.txt")
df.demo = fread(fn.demo) %>% 
  mutate(donor_short = sub("H5N1-0", "s", `Subject ID`)) %>% 
  mutate(Adjuvant = toupper(Adjuvant)) %>% 
  dplyr::select(donor = `Subject ID`, donor_short, sex= Gender, age = Age, unblinded_flag = Adjuvant)

# Load Titre Data
fn.mn = file.path(PROJECT_DIR, "DATA_ORIGINAL", "Titres", "titre.txt")
df.mn = fread(fn.mn) %>% 
  mutate_at(-(1:2), function(x) sub("<20","10",x) %>% as.numeric()) %>% 
  # mutate(day = as.numeric(sub("[dD]","",TimePt))) %>% 
  mutate(day = sub("D","d",TimePt)) %>% 
  mutate(day = ifelse(day %in% c("d29","d31"), "d28", day)) %>% 
  mutate(donor = sprintf("H5N1-%03d",`Sample ID`)) %>% 
  mutate(titer = log2(`A/Indonesia`/10)) %>% 
  # mutate(day = paste0("titer_I_MN_", day)) %>% 
  dplyr::select(donor, day, titer) %>% 
  spread(day, titer)
names(df.mn)[-1] = sub("^d", "titer_I_MN_d", names(df.mn)[-1])

# Load SOMASCAN Baseline Data
fn.rfu = file.path(PROJECT_DIR, "DATA_PROCESSED/SOMAscan", "Final_RFU.txt")
fn.samp = file.path(PROJECT_DIR, "DATA_PROCESSED/SOMAscan", "Final_Samples.txt")
fn.soma = file.path(PROJECT_DIR, "DATA_PROCESSED/SOMAscan", "Final_Somamers.txt")

df.soma = fread(fn.rfu, header = F, data.table = F)
df.samp = fread(fn.samp, data.table = F)
df.ann = fread(fn.soma, data.table = F)

names(df.soma) = df.ann$Target
df.soma$donor_short = df.samp$SampleId

# Put Everthing Together
df_chi = inner_join(df.demo, df.mn, by="donor") %>% 
  left_join(df.soma, by="donor_short") %>% 
  dplyr::filter(unblinded_flag == "ADJ")

# EMORY DATA
# Load Titre Data
df_mn_emory <- read.table(file.path(PROJECT_DIR, "DATA_PROCESSED", "Emory", "Emory_Adj_MN_D0_21_42.txt"), sep = "\t", header = T)
df_mn_emory[2:4] <- log2(df_mn_emory[2:4]/10)
df_mn_emory <- cbind(df_mn_emory, data.frame("Delta" = df_mn_emory$Day.42 - df_mn_emory$Day.21))

# Load SOMASCAN Baseline Data
fn.rfu.e = file.path(PROJECT_DIR, "DATA_PROCESSED", "Emory", "somascan","Final_RFU_MedNorm_Cal.txt")
fn.samp.e = file.path(PROJECT_DIR, "DATA_PROCESSED","Emory", "somascan", "Final_Samples_MedNorm_Cal.txt")
fn.soma.e = file.path(PROJECT_DIR, "DATA_PROCESSED", "Emory", "somascan", "Final_Somamers_MedNorm_Cal.txt")

df.soma.emory = fread(fn.rfu.e, header = F, data.table = F)
df.samp.emory = fread(fn.samp.e, data.table = F)
df.ann.emory = fread(fn.soma.e, data.table = F)

names(df.soma.emory) = df.ann.emory$Target
df.soma.emory$donor_short = df.samp.emory$SampleId


# Calcualte z-scores and correlations
cor_result <- data.frame(matrix(nrow = 0, ncol = 11))
for(k in c(6, 8)){
        if(k == 6){
                #                 chi_top <- c("RSPO3", "IGFBP-1", "HAI-1", "GX", "GPDA", "Eotaxin")
                chi_neg <- c( "HAI-1", "GX", "GPDA", "Eotaxin")
                chi_pos <- c("RSPO3", "IGFBP-1")
        } else if(k == 8){
                #                 chi_top <- c("RSPO3", "IGFBP-1", "HAI-1", "GX", "GPDA", "Eotaxin", "PSA", "transcription factor MLR1, isoform CRA_b")
                chi_neg <- c( "HAI-1", "GX", "GPDA", "Eotaxin", "PSA")
                chi_pos <- c("RSPO3", "IGFBP-1", "transcription factor MLR1, isoform CRA_b")
        }
        
        # CHI DATA
        # Select Data for the top Six/Eight Proteins and Compute z-scores
        gen_info <- c("donor", "donor_short", "unblinded_flag", "titer_I_MN_d28", "titer_I_MN_d42", "titer_I_MN_d21")
        #         top_val <- t(df_chi[, chi_top])
        top_neg <- t(df_chi[, chi_neg])
        top_pos <- t(df_chi[, chi_pos])
        top_flip <- rbind (top_neg * -1, top_pos)
        #         chi_z_score <- get_score(top_val)
        #         chi_neg_z_score <- get_score(top_neg)
        #         chi_pos_z_score <- get_score(top_pos)
        chi_flip_z_score <- get_score(top_flip)

        # Compute Correlation Between the z-scores and the Titre
        #         d28_cor <- cor.test(chi_z_score, df_chi[,gen_info]$titer_I_MN_d28, na.action = "na.exclude")
        #         d42_cor <- cor.test(chi_z_score, df_chi[,gen_info]$titer_I_MN_d42, na.action = "na.exclude")
        # 
        #         d28_cor_neg <- cor.test(chi_neg_z_score, df_chi[,gen_info]$titer_I_MN_d28, na.action = "na.exclude")
        #         d42_cor_neg <- cor.test(chi_neg_z_score, df_chi[,gen_info]$titer_I_MN_d42, na.action = "na.exclude")
        # 
        #         d28_cor_pos <- cor.test(chi_pos_z_score, df_chi[,gen_info]$titer_I_MN_d28, na.action = "na.exclude")
        #         d42_cor_pos <- cor.test(chi_pos_z_score, df_chi[,gen_info]$titer_I_MN_d42, na.action = "na.exclude")

        d21_cor_flip <- cor.test(chi_flip_z_score, df_chi[,gen_info]$titer_I_MN_d21, na.action = "na.exclude")
        #         d28_cor_flip <- cor.test(chi_flip_z_score, df_chi[,gen_info]$titer_I_MN_d28, na.action = "na.exclude")
        d42_cor_flip <- cor.test(chi_flip_z_score, df_chi[,gen_info]$titer_I_MN_d42, na.action = "na.exclude")


        # EMORY DATA
        #         emory_top_val <- t(df.soma.emory[, chi_top])
        emory_top_neg <- t(df.soma.emory[, chi_neg])
        emory_top_pos <- t(df.soma.emory[, chi_pos])
        emory_top_flip <- rbind (emory_top_neg * -1, emory_top_pos)

        # Calculate z-scores and write dat frames
        emory_flip_z_score <- get_score(emory_top_flip)
        # emory_neg_z_score <- get_score(emory_top_neg)
        # emory_pos_z_score <- get_score(emory_top_pos)
        df_emory_out <- data.frame("Subject" = df.samp.emory$SampleId, "z_score" = emory_flip_z_score)
        df_emory_out <- df_emory_out[order(df_emory_out$Subject),]

        # Combine titre and z-scores
        df_emory_flip <- merge(df_emory_out, df_mn_emory)

        # Calcualte Correlation between the z-scores and the Titre
        d21_cor_flip_emory <- cor.test(df_emory_flip$z_score, df_emory_flip$Day.21, na.action = "na.exclude")
        d42_cor_flip_emory <- cor.test(df_emory_flip$z_score, df_emory_flip$Day.42, na.action = "na.exclude")
        delta_cor_flip_emory <- cor.test(df_emory_flip$z_score, df_emory_flip$Delta, na.action = "na.exclude")
        
        # Consolidate the results for both the cohorts
        temp_df <- data.frame("Signature" = k, "CHI_d21_cor" = as.numeric(d21_cor_flip$estimate), "CHI_d21_p" = as.numeric(d21_cor_flip$p.value), "Emory_d21_cor" = as.numeric(d21_cor_flip_emory$estimate), "Emory_d21_p" =  as.numeric(d21_cor_flip_emory$p.value), "CHI_d42_cor" = as.numeric(d42_cor_flip$estimate), "CHI_d42_p" = as.numeric(d42_cor_flip$p.value), "Emory_d42_cor" = as.numeric(d42_cor_flip_emory$estimate),"Emory_d42_p" = as.numeric(d42_cor_flip_emory$p.value), "Delta_cor" = as.numeric(delta_cor_flip_emory$estimate), "Delta_p" = as.numeric(delta_cor_flip_emory$p.value))

        cor_result <- rbind(cor_result, temp_df)
}


# Write the output table.
cor_result <- round(cor_result,3)
write.table(cor_result, file.path(PROJECT_DIR, "RESULTS", "Emory", "chi_emory_d21_d42_cor.txt"), sep = "\t", row.names = FALSE)

stop()

# write.table(df_emory_out, file.path(PROJECT_DIR, "DATA_PROCESSED", "Emory", "somascan", "emory_flip_z_score_8.txt"), sep = "\t", row.names = FALSE)

# df_emory_out <- data.frame("SampleID" = df.samp.emory$SampleId, "z_score" = emory_neg_z_score)
# write.table(df_emory_out, file.path(PROJECT_DIR, "DATA_PROCESSED", "Emory", "somascan", "emory_neg_z_score.txt"), sep = "\t", row.names = FALSE)
# 
# df_emory_out <- data.frame("SampleID" = df.samp.emory$SampleId, "z_score" = emory_pos_z_score)
# write.table(df_emory_out, file.path(PROJECT_DIR, "DATA_PROCESSED", "Emory", "somascan", "emory_pos_z_score.txt"), sep = "\t", row.names = FALSE)
# 

# Legacy code
# Plot Scatter Plots Between z-score and the Titre
# d28_plot <- ggplot(data.frame("z_score" = chi_z_score, "titre_d28"=  df_chi[,gen_info]$titer_I_MN_d28), aes(x = z_score, y = titre_d28)) +
#         geom_point() +
#         geom_smooth(method = "lm") +
#         ylim(0,9.5) +
#         labs(x = "Per Subject z-score", y = "Day 28 MN Titre", title = paste("r = ", round(d28_cor$estimate[[1]], 3),  " p-value = ", round(d28_cor$p.value, 2), sep = ""))
# ggsave(file.path("FIGURES", "Emory", "d28_z_score.png"), plot = d28_plot, device = "png", width = 6, height = 4 )
# 
# d42_plot <- ggplot(data.frame("z_score" = chi_z_score, "titre_d42"=  df_chi[,gen_info]$titer_I_MN_d42), aes(x = z_score, y = titre_d42)) +
#         geom_point() +
#         geom_smooth(method = "lm") +
#         ylim(0,9.5) +
#         labs(x = "Per Subject z-score", y = "Day 42 MN Titre", title = paste("r = ", round(d42_cor$estimate[[1]], 3),  " p-value = ", round(d42_cor$p.value, 2), sep = ""))
# ggsave(file.path("FIGURES", "Emory", "d42_z_score.png"), plot = d42_plot, device = "png", width = 6, height = 4 )

