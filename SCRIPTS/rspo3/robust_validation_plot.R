library(tidyverse)

# Load robust correlation validation results. 
val_result <- read.table(file.path("RESULTS", "rspo3","val_result_all_with_pre_z_norm.txt"), header = T, sep = "\t")
# val_mean_sd <- val_result %>% group_by(subject) %>% summarize(z_mean = mean(z_score), z_sd = sd(z_score)/sqrt(length(z_score)))
val_mean_sd <- val_result %>% group_by(subject) %>% summarize(z_mean = mean(avg_score), z_sd = sd(avg_score)/sqrt(length(avg_score)))

# Load titre response data.
titre <- read.table(file.path("DATA_ORIGINAL", "Titres", "titre.txt"), header = TRUE, sep = "\t")
titre$TimePt <- tolower(titre$TimePt)
d28_titre <-  titre[titre$TimePt %in% c("d28", "d29", "d31") , c("Sample.ID","A.Indonesia")]
d28_titre$Sample.ID <- paste("H5N1_", str_pad(d28_titre$Sample.ID, 3, pad = "0"), sep = "")
d28_titre$A.Indonesia <-  as.numeric(str_replace(d28_titre$A.Indonesia, "<", ""))


# Calculate correlation between
d28_titre <- d28_titre[d28_titre$Sample.ID %in% as.character(val_mean_sd$subject),]
corel_ <- cor(d28_titre$A.Indonesia, val_mean_sd$z_mean, method = "spearman")

all_result <-  merge(val_mean_sd, d28_titre, by.y="Sample.ID", by.x="subject")
all_result$A.Indonesia <- log(all_result$A.Indonesia)

scat_plot <- ggplot(all_result, aes(x=A.Indonesia, y=z_mean)) +
                     geom_point() +
                     geom_errorbar(aes(ymin=z_mean-z_sd, ymax=z_mean+z_sd), width=.1) +
                     #                      geom_text(aes(label=subject),hjust=0, vjust=0) +
                     scale_y_continuous(name ="mean z-score") +
                     labs(x = "d28/A.Indonesia")

stop()
ggsave(file.path("FIGURES","rspo3", "val_all_with_pre_z_norm.png"), plot = scat_plot, device = "png", width = 6, height = 4)

