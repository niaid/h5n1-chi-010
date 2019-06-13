source(file.path(PROJECT_DIR, "SCRIPTS/functions/factor.date.r"))

fn.mn = file.path(PROJECT_DIR, "DATA_ORIGINAL", "Titres", "titre.txt")
df.mn = read.table(fn.mn, sep="\t", header=T, row.names=NULL, stringsAsFactors=T) %>% 
  mutate_at(-(1:2), function(x) sub("<20","10",x) %>% as.numeric()) %>% 
  # mutate(day = as.numeric(sub("[dD]","",TimePt))) %>% 
  mutate(day = sub("D","d",TimePt)) %>% 
  mutate(day = ifelse(day %in% c("d29","d31"), "d28", day) %>% factor.date(), subject=factor(Sample.ID))

yticks = sort(unique(df.mn$A.Indonesia))
yticks = yticks[seq(2,length(yticks),2)]

ggplot(df.mn, aes(day, log2(A.Indonesia), group=subject, col=subject)) +
  geom_line(position = position_jitter(width=0, height=0.1)) +
  scale_y_continuous(name="MN titer", breaks=log2(yticks), labels = yticks) +
  xlab("Time points") +
  theme_bw() + theme(legend.position = "none", panel.grid.minor=element_blank(),
                     axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

dn.fig = file.path(PROJECT_DIR, "FIGURES/titers")
dir.create(dn.fig, showWarnings = F)

ggsave(file.path(dn.fig, "MN_titer_profiles_all_subjects.png"), h=3.3, w=4)
