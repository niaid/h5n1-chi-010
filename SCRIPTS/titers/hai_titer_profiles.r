# PURPOSE: To plot HAI titer profiles against time. 

# Initialize the environment with PROJECT_DIR and some commonly used packages.
source("SCRIPTS/0_initialize.r")

# MAIN

# Read serology data.
fn.hai = file.path(PROJECT_DIR, "DATA_ORIGINAL", "HAI", "H5N1_serology.txt")
df.hai = read.table(fn.hai, sep="\t", 
                    header=T, row.names=NULL, stringsAsFactors=F) %>% 
  mutate_at(-(1), function(x) sub("<10","5",x) %>% as.numeric()) %>% 
  gather("tmp","hai",-subject) %>% 
  separate("tmp",c("day","replicate"))

df.hai3 = df.hai %>% 
  group_by(subject, day) %>% 
  summarise(hai.median=median(hai, na.rm=T), 
            hai.min=min(hai, na.rm=T), hai.max=max(hai, na.rm=T) ) %>% 
  ungroup() %>% 
  mutate(day=factor(day,levels=c("day0","day1","day7","day21","day22","day28","day42","day100")))

# Plot and save the figure file.
plot_hai <- ggplot(df.hai3, aes(day, log2(hai.median), group=subject, col=subject)) +
  geom_line(position = position_jitter(width=0, height=0.1)) +
  scale_y_continuous(name="HAI", breaks=log2(unique(df.hai3$hai.median)), labels = unique(df.hai3$hai.median)) +
  xlab("Time points") +
  theme_bw() + theme(legend.position = "none", panel.grid.minor=element_blank(),
                     axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

## Check if the figure output directory exists, else create it. 
dn.fig = file.path(PROJECT_DIR, "FIGURES/titers")
if (dir.exists(dn.fig)){
        print("Figure directory already exists.")
}else {
        print("Figure directory doesn't exists. Creating it now.")
        dir.create(dn.fig, recursive = T, showWarnings = T)
}
ggsave(file.path(dn.fig, "HAI_titer_profiles_all_subjects.png"), plot = plot_hai, h=3.3, w=4)


# Check subjects that never responded to the titer.
# Find max of hai.max titer per subject and identify subjects that are =< 5.
subj_no_resp <- df.hai3[is.finite(df.hai3$hai.max),] %>% group_by(subject) %>% 
                                         summarize_at(vars("hai.max"), max) %>%
                                         dplyr::filter(hai.max <= 5)

write_tsv(subj_no_resp, file.path(dn.fig, "subjects_no_resp.tsv"))
