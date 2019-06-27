ss = "H5N1-020"

dn.in = file.path(PROJECT_DIR, "DATA_PROCESSED/Luminex/")

source(file.path(PROJECT_DIR, "SCRIPTS/functions/hours2days.r"))
df.ip10 = readRDS(file.path(dn.in, "luminex_data.rds")) %>% 
  filter(analyte == "IP-10")%>% 
  filter(subject %in% ss) %>% 
  mutate(time = hours2days(Time, format.h="d%d+%dh", format.d="d%d") %>% fct_inorder()) %>% 
  group_by(subject, time) %>% 
  summarise(value=mean(value, na.rm=T)) %>% 
  ungroup()

ggplot(df.ip10, aes(time, value, group=subject)) +
  geom_line() +
  geom_vline(xintercept = c(5, 11), lty=2, col="black") +
  # geom_vline(xintercept = c(1, 7), lty=1, col="blue") +
  xlab("Time") +
  ylab("IP-10") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

fn.fig = file.path(PROJECT_DIR, "FIGURES/IP-10_profiles_sel.subject")
ggsave(paste0(fn.fig, ".png"), w=2.5, h=2)
ggsave(paste0(fn.fig, ".pdf"), w=2.5, h=2)
