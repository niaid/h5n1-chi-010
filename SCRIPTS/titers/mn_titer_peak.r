# PURPOSE: To plot titer at peak.

# Initialize the environment with PROJECT_DIR and commonly used packages.
source("SCRIPTS/0_initialize.r")
library(gridExtra)

# FUNCTIONS
source(file.path(PROJECT_DIR, "SCRIPTS/functions/factor.date.r"))

# MAIN
# Read titer data.
fn.mn = file.path(PROJECT_DIR, "DATA_ORIGINAL", "Titres", "titre.txt")
df.mn = read.table(fn.mn, sep="\t", header=T, row.names=NULL, stringsAsFactors=T) %>% 
  mutate_at(-(1:2), function(x) sub("<20","10",x) %>% as.numeric()) %>% 
  # mutate(day = as.numeric(sub("[dD]","",TimePt))) %>% 
  mutate(day = sub("D","d",TimePt)) %>% 
  mutate(day = ifelse(day %in% c("d29","d31"), "d28", day) %>% factor.date(), subject=factor(Sample.ID))

yticks = sort(unique(df.mn$A.Indonesia))
yticks = yticks[seq(2,length(yticks),2)]

# check subjects with peak at d42
tm = levels(df.mn$day)

# Read demographics and adjuvant status data. 
fn.demo = file.path(PROJECT_DIR, "DATA_ORIGINAL/Clinical/clinical_info_adj.txt")
df.demo = fread(fn.demo, data.table = F) %>% 
  mutate(Sample.ID = sub("H5N1-0","", `Subject ID`) %>% as.numeric())


DF = df.mn %>% 
  group_by(Sample.ID) %>% 
  summarize(peak.time = tm[which.max(A.Indonesia)], 
            at.peak = A.Indonesia[which.max(A.Indonesia)],
            after.peak = A.Indonesia[which.max(A.Indonesia)+1]) %>% 
  mutate(decline = (at.peak-after.peak)/at.peak * 100) %>% 
  left_join(df.demo, by="Sample.ID") %>% 
  mutate(peak.time = factor(peak.time, levels = tm))

DF %>% filter(peak.time=="d42")

# Plots
p1 = ggplot(DF %>% mutate(peak.time=factor(peak.time, levels=tm)), 
       aes(peak.time, group=Adjuvant, fill=Adjuvant)) +
  geom_bar() +
  theme_bw()

p2 = ggplot(DF, aes(log2(at.peak), decline, col=Adjuvant)) +
  geom_jitter(width=0.2, height=2, size=3) +
  scale_x_continuous(name = "MN titer at peak", breaks=log2(yticks), labels = yticks) +
  ylab("decline after peak, %") +
  theme_bw()

## Merge the two plots created above in a single file.
mg = marrangeGrob(list(p1,p2), nrow=1, ncol=2, top = NULL)

## Check if figure output direcroty exists and save the figure to a file.
dn.fig = file.path(PROJECT_DIR, "FIGURES/titers")
if (dir.exists(dn.fig)){
        print("Figure directory exists.")
}else {
        print("Figure output directory doesn't exists. Creating it now.")
        dir.create(dn.fig, recursive = T, showWarnings = T)
}
ggsave(file.path(dn.fig, "MN_titer_peak_decline.png"), plot=mg, w=7, h=2.5)

# Save output table to a text file. 
df.out = df.mn %>% 
  group_by(Sample.ID) %>% 
  summarize(peak.time = tm[which.max(A.Indonesia)], 
            at.peak = A.Indonesia[which.max(A.Indonesia)],
            after.peak = A.Indonesia[which.max(A.Indonesia)+1]) %>% 
  mutate(decline = (at.peak-after.peak)/at.peak * 100) %>% 
  mutate(Sample.ID = sprintf("H5N1-%03d", Sample.ID)) %>% 
  dplyr::rename(subject=Sample.ID)

## Check if the output directory exists.
dn.out = file.path(PROJECT_DIR, "RESULTS/titers")
if(dir.exists(dn.out)){
        print("Output directory exists.")
}else{
        print("Output directory doesn't exists. Creating it now.")
        dir.create(dn.out, recursive = T, showWarnings = T)
}
fwrite(df.out, file.path(dn.out, "MN_titer_peak_time.txt"), sep="\t")
