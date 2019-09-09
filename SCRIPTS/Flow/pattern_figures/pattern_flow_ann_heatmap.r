# PURPOSE: To plot a heatmap with the type of cell populations found per pattern
# in the flow data.

# Initialize. 
source("SCRIPTS/0_initialize.r")
library(ComplexHeatmap)
library(circlize)

dn.in = file.path(PROJECT_DIR, "RESULTS/Flow_10c/")
fn = "DeltaThr0.5_DonorThr0.3_PercentOfParent_CORpearson_cutreeHybrid_deepSplit0_minClusterSize31_CPs.txt"
df = fread(file.path(dn.in, fn))  %>% 
        mutate(ID = paste0("ID", CP_code)) %>% 
  # mutate(CP_desc = sprintf("ID%g: %s", CP_code, CP_desc)) %>% 
        select(-CP_code, -CP_desc)

# mat  = df %>% 
#   tibble::column_to_rownames("CP_desc") %>% 
#   data.matrix()

id = df$CP_desc %>% str_extract("ID[\\d\\.]+")
id.order = c(5.1, 5.2, 68, 15.1, 15.2, 30, 31, 54.1, 52, 55, 57, 59, 61, 63, 64) %>% paste0("ID",.)
ro = match(id.order, df$ID) %>% rev()

# ro = hclust(dist(mat, method = "man"), method = "comp")$order

fn.ann = file.path(PROJECT_DIR, "DATA_ORIGINAL/Flow_10c/Flow_10c_ann.txt")
# fn.ann = file.path(PROJECT_DIR, "DATA_PROCESSED/Flow_10c/Flow_10c_ann.txt")
df.ann = fread(fn.ann, header = T) %>% 
        mutate(Annot = paste0(ID, ": ", Name2))

df = df %>% 
  left_join(df.ann %>% dplyr::select(ID, Annot), by="ID")

df.tidy = df %>%
  gather("Pattern", "value", -Annot, -ID) %>% 
  mutate(value = ifelse(value==0, NA, value) %>% factor()) %>% 
  mutate(Annot = factor(Annot, levels=df$Annot[ro]))

plot1 <- ggplot(df.tidy, aes(Pattern, Annot, fill=value)) +
  geom_tile() +
  scale_fill_manual(values = palette(), guide=F) +
  scale_y_discrete(position = "right") +
  xlab("") + ylab("") +
  theme_bw() +
  theme(axis.ticks = element_blank(), 
        axis.text.x = element_text(angle=90, colour = palette()), 
        panel.grid = element_blank())

fn.fig = file.path(PROJECT_DIR, "FIGURES/Flow_patterns_annot_heatmap")
ggsave(paste0(fn.fig, ".png"), plot = plot1,  w=4.8, h=2)
ggsave(paste0(fn.fig, ".pdf"), plot = plot1, w=4.8, h=2)
