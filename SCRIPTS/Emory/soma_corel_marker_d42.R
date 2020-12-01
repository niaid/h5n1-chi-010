# PURPOSE: Correlation between individual markers from our signature and d42 titer in Bali cohort

source("SCRIPTS/0_initialize.r")

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


# chi_top <- c("RSPO3", "IGFBP-1", "HAI-1", "GX", "GPDA", "Eotaxin", "PSA", "transcription factor MLR1, isoform CRA_b")
chi_neg <- c( "HAI-1", "GX", "GPDA", "Eotaxin", "PSA")
chi_pos <- c("RSPO3", "IGFBP-1", "transcription factor MLR1, isoform CRA_b")
chi_top <- c(chi_neg, chi_pos)        

# EMORY DATA
emory_top_neg <- df.soma.emory[, chi_neg]
emory_top_pos <- df.soma.emory[, chi_pos]
# emory_top_flip <- cbind (emory_top_neg * -1, emory_top_pos, Subject = df.soma.emory$donor_short)
# Not flipping the sign here. 
emory_top_flip <- cbind (emory_top_neg, emory_top_pos, Subject = df.soma.emory$donor_short)
emory_top_flip <- emory_top_flip[order(emory_top_flip$Subject),]

df_emory_flip <- merge(emory_top_flip, df_mn_emory)

cor_result <- data.frame(matrix(nrow = 0, ncol = 5))
# Calcualte Correlation between individual markers and the Titre
for(mark in chi_top){
        d21_cor_flip_emory <- cor.test(df_emory_flip[[mark]], df_emory_flip$Day.21, na.action = "na.exclude")
        d42_cor_flip_emory <- cor.test(df_emory_flip[[mark]], df_emory_flip$Day.42, na.action = "na.exclude")

        temp_df <- data.frame("Marker" = mark, "Emory_d21_cor" = as.numeric(d21_cor_flip_emory$estimate), "Emory_d21_p" =  as.numeric(d21_cor_flip_emory$p.value),  "Emory_d42_cor" = as.numeric(d42_cor_flip_emory$estimate),"Emory_d42_p" = as.numeric(d42_cor_flip_emory$p.value))
        cor_result <- rbind(cor_result, temp_df)
}
        
# Write the output table.
cor_result[2:5] <- round(cor_result[2:5],3)
write.table(cor_result, file.path(PROJECT_DIR, "RESULTS", "Emory", "emory_d21_d42_cor.txt"), sep = "\t", row.names = FALSE)


