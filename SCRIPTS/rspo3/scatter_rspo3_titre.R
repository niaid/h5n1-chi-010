# PURPOSE: To plot scatter plots between RSPO3 and the day 28 and day 42 titres. 

library(tidyverse)
library(ggrepel) # Install in the container.

# Load subject info and select adjuvanted subjects
subject <- read.table(file.path("DATA_ORIGINAL","Clinical", "clinical_info_adj.txt"), header = TRUE, sep = "\t")
subject_ori <- subject[subject$Adjuvant == "Adj", c("Subject.ID", "Gender")]
subject_ori$Subject.ID <- str_replace(subject_ori$Subject.ID, "-", "_")
subject <- subject[subject$Adjuvant == "Adj", c("Subject.ID")]
subject <- str_replace(subject, "-", "_")
# Remove subject H5N1_010
subject <- subject[!subject == "H5N1_010"]
subject_ori <- subject_ori[!subject_ori$Subject.ID == "H5N1_010", ]
colnames(subject_ori) <- c("Sample.ID", "Gender")

# Load SOMAScan data for RSPO3.
rspo3_soma <- read.table(file.path("New", "190703", "DATA", "SOMA.txt"), header = TRUE, sep = "\t")
rspo3_soma <- rspo3_soma[c("donor", "RSPO3")]
rspo3_soma$donor <- str_replace(rspo3_soma$donor, "-", "_")
colnames(rspo3_soma) <- c("Sample.ID", "RSPO3")

# Load titre response data.
titre <- read.table(file.path("DATA_ORIGINAL", "Titres", "titre.txt"), header = TRUE, sep = "\t")
titre$TimePt <- tolower(titre$TimePt)
d28_titre <-  titre[titre$TimePt %in% c("d28", "d29", "d31") , c("Sample.ID","A.Indonesia")]
d28_titre$Sample.ID <- paste("H5N1_", str_pad(d28_titre$Sample.ID, 3, pad = "0"), sep = "")
d28_titre$A.Indonesia <-  as.numeric(str_replace(d28_titre$A.Indonesia, "<", ""))

d42_titre <-  titre[titre$TimePt == "d42" , c("Sample.ID","A.Indonesia")]
d42_titre$Sample.ID <- paste("H5N1_", str_pad(d42_titre$Sample.ID, 3, pad = "0"), sep = "")
d42_titre$A.Indonesia <-  as.numeric(str_replace(d42_titre$A.Indonesia, "<", ""))

# Find common subjects among all the datasets. 
common_subj <- Reduce(intersect, list(subject, rspo3_soma$Sample.ID, d28_titre$Sample.ID, d42_titre$Sample.ID))

# Selected data
rspo3_selected <- rspo3_soma[rspo3_soma$Sample.ID %in% common_subj,]
d28_titre_selected <- d28_titre[d28_titre$Sample.ID %in% common_subj,]
d42_titre_selected <- d42_titre[d42_titre$Sample.ID %in% common_subj,]

rspo3_d28 <- merge(rspo3_selected, d28_titre_selected)
rspo3_d28$A.Indonesia <- log(rspo3_d28$A.Indonesia)
rspo3_d28 <- merge(rspo3_d28, subject_ori)
rspo3_d42 <- merge(rspo3_selected, d42_titre_selected)
rspo3_d42$A.Indonesia <- log(rspo3_d42$A.Indonesia)
rspo3_d42 <- merge(rspo3_d42, subject_ori)

# Scatter plot
# With the outlier(s).
corel <- cor.test(rspo3_d28$RSPO3, rspo3_d28$A.Indonesia, method = "spearman")
rspo3_d28_plot <- ggplot(rspo3_d28,aes(x=RSPO3, y=A.Indonesia)) +
             geom_point(aes(color = Gender)) +
             geom_text_repel(size = 2, aes(label=Sample.ID)) +
             labs(x = "RSPO3", y = "D28/A.Indonesia", color = "Gender", title = paste("rho:", round(corel$estimate, digits = 2), "p-value:", round(corel$p.value, digits = 5)) )+
                  
             ylim(5, 9) + xlim(3.1,3.4)
plot_name <- "rspo3_d28_scatter.png"

ggsave(file.path("FIGURES", "rspo3", plot_name), plot = rspo3_d28_plot, device = "png", width = 5, height = 4)

corel <- cor.test(rspo3_d42$RSPO3, rspo3_d42$A.Indonesia, method = "spearman")
rspo3_d42_plot <- ggplot(rspo3_d42, aes(x=RSPO3, y=A.Indonesia)) +
             geom_point(aes(color = Gender)) +
             geom_text_repel(size = 2, aes(label=Sample.ID)) +
             labs(x = "RSPO3", y = "D42/A.Indonesia", color = "Gender", title = paste("rho:", round(corel$estimate, digits = 2), "p-value:", round(corel$p.value, digits = 5)))+
             ylim(5, 9) + xlim(3.1, 3.4)
plot_name <- "rspo3_d42_scatter.png"
ggsave(file.path("FIGURES", "rspo3", plot_name), plot = rspo3_d42_plot, device = "png", width = 5, height = 4)

# Without the outlier(s)
rspo3_d28 <- rspo3_d28[!rspo3_d28$Sample.ID %in% c("H5N1_034", "H5N1_038"),]
rspo3_d42 <- rspo3_d42[!rspo3_d42$Sample.ID %in% c("H5N1_034", "H5N1_038"),]

corel <- cor.test(rspo3_d28$RSPO3, rspo3_d28$A.Indonesia, method = "spearman")
rspo3_d28_plot <- ggplot(rspo3_d28, aes(x=RSPO3, y=A.Indonesia)) +
             geom_point(aes(color = Gender)) +
             geom_text_repel(size = 2, aes(label=Sample.ID)) +
             labs(x = "RSPO3", y = "D28/A.Indonesia", color = "Gender", title = paste("rho:", round(corel$estimate, digits = 2), "p-value:", round(corel$p.value, digits = 5))) +
             ylim(5, 9) + xlim(3.1,3.4)
plot_name <- "rspo3_d28_scatter_without.png"

ggsave(file.path("FIGURES", "rspo3", plot_name), plot = rspo3_d28_plot, device = "png", width = 5, height = 4)

corel <- cor.test(rspo3_d42$RSPO3, rspo3_d42$A.Indonesia, method = "spearman")
rspo3_d42_plot <- ggplot(rspo3_d42, aes(x=RSPO3, y=A.Indonesia)) +
             geom_point(aes(color = Gender)) +
             geom_text_repel(size = 2, aes(label=Sample.ID)) +
             labs(x = "RSPO3", y = "D42/A.Indonesia", color = "Gender", title = paste("rho:", round(corel$estimate, digits = 2), "p-value:", round(corel$p.value, digits = 5))) +
             ylim(5,9) + xlim(3.1,3.4)
plot_name <- "rspo3_d42_scatter_without.png"
ggsave(file.path("FIGURES", "rspo3", plot_name), plot = rspo3_d42_plot, device = "png", width = 5, height = 4)

