source("SCRIPTS/0_initialize.r")
source("SCRIPTS/functions/today.r")
fn = "DATA_ORIGINAL/Emory/Vax010_IP10_luminex.txt"
df = read.table(fn, sep="\t", header=T, row.names=NULL, stringsAsFactors = F)
colnames(df) = c("sample","plate", "ip10")
df = df %>% 
  separate("sample",c("subject","day"),sep=" ", extra = "merge", remove = F) %>%
  mutate(day = sub("Day ","D",day)) %>% 
  mutate(day=factor(day, levels=unique(.$day)))

ss = c("HCAF017", "HCAF026")
ggplot(df %>% filter(subject %in% ss), aes(day, ip10, group=subject)) +
  geom_line() +
  facet_wrap(~subject, nrow=1) +
  geom_vline(xintercept=c(2,7), lty=2) +
  # scale_y_log10() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1),
        strip.background = element_blank())
ggsave(sprintf("FIGURES/Emory/emory_profiles_wrong_subjects_ip10_%s.png",today()), w=4,h=2)

