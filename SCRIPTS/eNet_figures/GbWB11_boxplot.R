# PURPOSE: To test if GbWB11 WGCNA module differs significantly between 
# male and female in unadjuvanted subjects, even though it was not 
# selected as an independent co-variate in e-Net analysis. 

# Load external libraries.
library(tidyverse)

# Load data.
df <- read_tsv(file.path("DATA_PROCESSED", "eNet", "eNet_Yuri", "eNet_InputData_r6.txt"))
# df1 <- read_tsv(file.path("DATA_PROCESSED", "eNet", "eNet_InputData_r2.txt"))

# Wilcoxon test
w_test <- wilcox.test(GbWB11 ~ sex, data = df, exact = FALSE)
print(w_test$p.value)

# Wilcoxon test sex vs titer.
w_test_titer  <- wilcox.test(titer_I_MN_d42 ~ sex, data =df, exact = F)
print(w_test_titer$p.value)

# rd_summary["cluster_id"] <- paste("SOM", rd_summary[["cluster_id"]], sep = "")
# Create a boxplot.
fig1 <- df %>% ggplot(aes(x=sex, y=titer_I_MN_d42), fill=sex) + 
  scale_color_manual(values=c("pink", "blue")) +
  geom_boxplot() +
  labs(y = "Titer") +
  scale_x_discrete(name = "Gender", breaks=c("F","M"), labels=c("Female","Male")) +
  scale_fill_discrete(guide=FALSE) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)

ggsave(filename = "GbWB11_boxplot_bas_neg_titer42_Yuri.tiff", plot = fig1, 
	device = "tiff", path = file.path("FIGURES", "Baseline"),
	scale = 1, width = 5, height = 7, units = "in",
	dpi = 300)

