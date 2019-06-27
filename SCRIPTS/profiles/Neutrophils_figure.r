fn = file.path(PROJECT_DIR, "DATA_ORIGINAL/Clinical/H5N1_BTRIS.txt")
btris = read.table(fn, sep="\t", header=T, row.names=NULL, stringsAsFactors = F)
neu = btris %>% 
  select(subject=assay, time=Timepoint, value=Neutrophils) %>% 
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
  
neu.mean = neu %>% 
  group_by(time) %>% 
  summarise(value=mean(value, na.rm=T)) %>% 
  ungroup() %>% 
  mutate(subject="mean")

ggplot(neu, aes(time, value, group=subject)) +
  geom_line(alpha=0.2) +
  geom_line(data=neu.mean, size=1) +
  # geom_vline(xintercept = c(4, 10), lty=2, col="red") +
  # geom_vline(xintercept = c(1, 7), lty=1, col="blue") +
  # scale_y_log10() +
  xlab("Time points") +
  ylab("Neutrophils") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

dn.fig = file.path(PROJECT_DIR, "FIGURES/profiles")
dir.create(dn.fig, showWarnings = F)
fn.fig = file.path(dn.fig, "Neutrophils_profiles")
ggsave(paste0(fn.fig, ".png"), h=3.5, w=4)
ggsave(paste0(fn.fig, ".pdf"), h=3.5, w=4)

