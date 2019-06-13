fn.tfh = file.path(PROJECT_DIR, "DATA_PROCESSED/Flow_10c/CXCR3_freq.txt")
tfh = fread(fn.tfh, data.table=F, header = T) %>% 
  rename(subject=`PATIENT ID`, time=`SAMPLE ID`) %>% 
  rename(T_CD3CD4 = `4.2`, T_CD3 = `1.2`, Tfh = `24`) %>% 
  filter(subject!="")

tfh = tfh %>% 
  mutate(time.point = case_when(
    time=="d8" ~ "d7",
    time=="d21_2h" ~ "d21_4h",
    time=="d29" ~ "d28",
    time=="d31" ~ "d28",
    time=="d43" ~ "d42",
    time=="d0_0h" ~ "d0",
    time=="d21_0h" ~ "d21",
    TRUE ~ time
  )) %>%
  mutate(time.point = sub("_","+",time.point)) %>% 
  mutate(time.point = factor(time.point, levels = unique(time.point)))

tfh.tm = unique(tfh$time.point)

### demographics and adjuvant status
fn.clin = file.path(PROJECT_DIR, "DATA_ORIGINAL/Clinical/clinical_info_adj.txt")
df.clin = fread(fn.clin) %>% 
  rename(subject=`Subject ID`)

df.tfh = tfh %>% 
  select(subject, time.point, one_of(c("CXCR3+","CXCR3-"))) %>% 
  filter(subject != "H5N1-044") %>% 
  dplyr::rename(time=time.point) %>% 
  gather(population, pop, -subject, -time) %>% 
  inner_join(df.clin %>% select(subject, Adjuvant), by="subject") %>% 
  group_by(time, Adjuvant, population) %>% 
  summarize(pop.mean = mean(pop, na.rm=T), pop.sd = sd(pop, na.rm=T)) %>% 
  ungroup()

tm.use = c("d0","d1","d7","d21","d22","d28","d42")
df.tfh = df.tfh %>% filter(time %in% tm.use) %>% mutate(time = factor(time, levels=tm.use))
ggplot(df.tfh, 
       aes(time, pop.mean, fill=population)) +
  geom_col(col="black") +
  geom_errorbar(data=filter(df.tfh, population == "CXCR3+"), aes(ymin=pop.mean-(pop.sd), ymax=pop.mean+(pop.sd)), width=0.2) +
  facet_wrap(~Adjuvant) +
  scale_fill_manual(values=c("lightyellow", "brown")) +
  xlab("Time") +
  ylab("Mean % of total Tfh cells") +
  theme_bw() +
  theme(strip.background = element_blank(), panel.border = element_blank(),
        panel.background = element_blank(), panel.grid.major = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
        axis.ticks.x = element_blank())

dn.fig = file.path(PROJECT_DIR, "FIGURES/Tfh")
fn.fig = file.path(dn.fig, "Tfh_CXCR3_cells_proportion_bars.png")
ggsave(fn.fig, w=6, h=3)

