fn.hai = file.path(PROJECT_DIR, "DATA_ORIGINAL", "HAI", "H5N1 serology.txt")
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

ggplot(df.hai3, aes(day, log2(hai.median), group=subject, col=subject)) +
  geom_line(position = position_jitter(width=0, height=0.1)) +
  scale_y_continuous(name="HAI", breaks=log2(unique(df.hai3$hai.median)), labels = unique(df.hai3$hai.median)) +
  xlab("Time points") +
  theme_bw() + theme(legend.position = "none", panel.grid.minor=element_blank(),
                     axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

dn.fig = file.path(PROJECT_DIR, "FIGURES/titers")
dir.create(dn.fig, showWarnings = F)

ggsave(file.path(dn.fig, "HAI_titer_profiles_all_subjects.png"), h=3.3, w=4)
