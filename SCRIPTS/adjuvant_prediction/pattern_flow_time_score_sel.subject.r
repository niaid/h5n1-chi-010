source(file.path(PROJECT_DIR, "SCRIPTS/functions/get_score.r"))

fn = file.path(PROJECT_DIR, "RESULTS/Flow_10c/DeltaThr0.5_DonorThr0.3_PercentOfParent_CORpearson_cutreeHybrid_deepSplit0_minClusterSize31_CPs.txt")
flow.patt = fread(fn, data.table = F)


fn.tm = file.path(PROJECT_DIR, "RESULTS/Flow_10c/TrajMatrix_DeltaThr0.5_DonorThr0.3_PercentOfParent.txt")
dat = fread(fn.tm, data.table = F)
fn.row = file.path(PROJECT_DIR, "RESULTS/Flow_10c/TrajMatrixLabels_DeltaThr0.5_DonorThr0.3_PercentOfParent.txt")
df.row = fread(fn.row, data.table = F, select = 1:3, col.names = c("tube","ID","subject"))
fn.col = file.path(PROJECT_DIR, "RESULTS/Flow_10c/DeltaThr0.5_DonorThr0.3_PercentOfParent_CORpearson_cutreeHybrid_deepSplit0_minClusterSize31_Modules.txt")
df.col = fread(fn.col, data.table = F)
names(dat) = df.col$timepoint

ss = "s20"
df.score = data.frame()
K = c(1,5)
for (k in K) {
  pat = paste0("Fp0",k)
  pat.ids = flow.patt %>% select_("CP_code", value = pat) %>% filter(value!=0) %>% pull(CP_code)
  
  ii = df.row$ID %in% pat.ids & df.row$subject %in% ss
  tmp = get_score(dat[ii,]) %>% as.data.frame() %>% 
    rename_(Score=".") %>% 
    tibble::rownames_to_column("Time") %>% 
    mutate(Pattern = pat)
    df.score = rbind(df.score, tmp)
}
df.score = df.score %>% 
  mutate(Time = sub("_","+",Time)) %>% 
  mutate(Time = factor(Time, levels=unique(Time)))

ggplot(df.score, aes(Time, Score, group=Pattern, col=Pattern)) +
  geom_line(size=1) +
  scale_color_manual(values = palette()[K]) +
  xlab("Time") + ylab("Score") +
  geom_vline(xintercept = c(3,8), lty=2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size=6, angle = 90, hjust = 1, vjust=0.5),
    # axis.text.y = element_blank(),
    legend.position = "none", #c(0,1), legend.justification = c(-0.1, 1),
    legend.title = element_blank(),
    legend.text = element_text(size=10),
    legend.background = element_blank())

fn.fig = file.path(PROJECT_DIR, "FIGURES/pattern_flow_profiles_sel.subject")
ggsave(paste0(fn.fig, ".png"), w=2.5, h=2)
ggsave(paste0(fn.fig, ".pdf"), w=2.5, h=2)

