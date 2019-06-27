dn.in = file.path(PROJECT_DIR, "DATA_PROCESSED/Luminex")
df = readRDS(file.path(dn.in, "luminex_data.rds")) %>% 
  filter(analyte == "IP-10")

df.mean = df %>% 
  group_by(time) %>% 
  summarise(value=mean(value, na.rm=T)) %>% 
  ungroup()

ggplot(df, aes(time, value, group=subject)) +
  geom_line(alpha=0.2) +
  geom_line(data=df.mean, group=1, size=1) +
  # geom_vline(xintercept = c(5, 11), lty=2, col="black") +
  # geom_vline(xintercept = c(1, 7), lty=1, col="blue") +
  xlab("Time points") +
  ylab("IP-10") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

dn.fig = file.path(PROJECT_DIR, "FIGURES/profiles")
dir.create(dn.fig, showWarnings = F)
fn.fig = file.path(dn.fig, "IP-10_profiles")
ggsave(paste0(fn.fig, ".png"), h=3.5, w=4)
ggsave(paste0(fn.fig, ".pdf"), h=3.5, w=4)
