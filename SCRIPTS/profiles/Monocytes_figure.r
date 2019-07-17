# PURPOSE: To plot total monocytes vs time.

# Initialize the environment with PROJECT_DIR and commonly used packages.
source("SCRIPTS/0_initialize.r")

# FUNCTIONS
source("SCRIPTS/functions/factor.date.r")

# MAIN

# Load BTRIS data.
fn = file.path(PROJECT_DIR, "DATA_ORIGINAL/Clinical/H5N1_BTRIS.txt")
btris = read.table(fn, sep="\t", header=T, row.names=NULL, stringsAsFactors = F)
mono = btris %>% 
  select(subject=assay, time=Timepoint, value=Monocytes_Absolute) %>% 
  filter(grepl("H5N1",subject)) %>% 
  mutate(time = sub("T","",time) %>% as.numeric()) %>% 
  arrange(time, subject) %>% 
  mutate(day = floor(time/24)) %>% 
  mutate(hours = time-day*24) %>% 
  mutate(time = case_when(
    hours==0 ~ sprintf("d%d", day),
    hours>0 ~ sprintf("d%d+%dh", day, hours))) %>% 
  select(-day, -hours) %>% 
  mutate(time = factor.date(time))
  
mono.mean = mono %>% 
  group_by(time) %>% 
  summarise(value=mean(value, na.rm=T)) %>% 
  ungroup() %>% 
  mutate(subject = "mean")

# Plot and save the figure file.
plot_mono  <- ggplot(mono, aes(time, value, group=subject)) +
  geom_line(alpha=0.2) +
  geom_line(data=mono.mean, size=1) +
  # geom_vline(xintercept = c(4, 10), lty=2, col="red") +
  # geom_vline(xintercept = c(1, 7), lty=1, col="blue") +
  # scale_y_log10() +
  xlab("Time points") +
  ylab("Total Monocytes") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

## Check if the figure output directory exists.
dn.fig = file.path(PROJECT_DIR, "FIGURES/profiles")
if(dir.exists(dn.fig)){
        print("Figure output directory exists.")
}else{
        print("Figure output directory doesn't exists. Creating it now.")
        dir.create(dn.fig, recursive = T, showWarnings = T)
}
fn.fig = file.path(dn.fig, "Monocytes_profiles")
ggsave(paste0(fn.fig, ".png"), plot = plot_mono,  h=3.5, w=4)
ggsave(paste0(fn.fig, ".pdf"), plot = plot_mono, h=3.5, w=4)
