source("SCRIPTS/0_initialize.r")
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

df.tfh = tfh %>% 
  select(subject, time.point, one_of(c("T_CD3CD4", "Tfh"))) %>% 
  dplyr::rename(time=time.point)

### read BTRIS data
fn.btris = file.path(PROJECT_DIR, "DATA_ORIGINAL/Clinical/H5N1_BTRIS.txt")
btris = read.table(fn.btris, sep="\t", header=T, row.names=NULL, stringsAsFactors = F)
df.t = btris %>% 
  select(subject=assay, time=Timepoint, T_CD3CD4_count=CD4_CD3_count, 
         T_CD3_count=CD3_count, B_count=CD19_count, Lymphocytes_count=Lymphocytes_Absolute,
         Neutrophils_count=Neutrophils_Absolute, Monocytes_count=Monocytes_Absolute) %>% 
  mutate_at(vars(Lymphocytes_count:Monocytes_count), `*`, 1000) %>%
  mutate(NL_ratio = Neutrophils_count/Lymphocytes_count) %>% 
  filter(grepl("H5N1",subject), !is.na(T_CD3CD4_count)) %>% 
  mutate(time = sub("T","",time) %>% as.numeric()) %>% 
  arrange(time, subject) %>% 
  mutate(day = floor(time/24)) %>% 
  mutate(hours = time-day*24) %>% 
  mutate(time = case_when(
    hours==0 ~ sprintf("d%d", day),
    hours>0 ~ sprintf("d%d+%dh", day, hours))) %>% 
  select(-day, -hours) %>% 
  mutate(time = factor(time, levels=unique(time))) %>% 
  filter(time != "d0+12h") # only 1 sample

DF = inner_join(df.tfh, df.t, by=c("subject","time")) %>% 
  mutate(time = factor(time, levels=unique(time)))
DF = DF %>% 
  mutate(Tfh_count = T_CD3CD4_count * Tfh / 100)
DF = DF %>% select(subject, time, matches("_count"), "NL_ratio") %>% 
  gather(population, count, -c(subject, time)) %>% 
  mutate(population = sub("_count", "", population))

### demographics and adjuvant status
fn.clin = file.path(PROJECT_DIR, "DATA_ORIGINAL/Clinical/clinical_info_adj.txt")
df.clin = fread(fn.clin) %>% 
  rename(subject=`Subject ID`)

DF = DF %>% 
  inner_join(df.clin %>% select(subject, Adjuvant), by="subject")

pop_order = c("Lymphocytes", "Neutrophils")
pop.in = c("Lymphocytes", "Neutrophils")

# pop_order = c("NL_ratio", "Tfh")
# pop.in = c("NL_ratio", "Tfh")

DF = DF %>% 
  filter(population %in% pop.in) %>% 
  filter(subject != "H5N1-044") %>% 
  mutate(population = factor(population, levels=pop_order))


# ggplot(DF %>% filter(population %in% pop.in), aes(time, count, group=subject, col=Adjuvant)) +
#   geom_line() +
#   facet_wrap(~population, ncol=1, scales="free_y") + 
#   theme_bw()
# ggsave(glue::glue("Tfh_profiles_flow10_counts_generated_samples.png"), w=7.5, h=7)

DF.change = DF %>% 
  group_by(subject, population) %>% 
  mutate(n = sum(time=="d0")) %>% 
  filter(n==1) %>% 
  mutate(count = count / count[time=="d0"])

# ggplot(DF.change %>% filter(population %in% pop.in), aes(time, count, group=subject, col=Adjuvant)) +
#   geom_line() +
#   facet_wrap(~population, ncol=1) + #, scales="free_y") + 
#   theme_bw()
# ggsave(glue::glue("Tfh_profiles_flow10_counts_generated_change_samples.png"), w=7.5, h=7)

DF.change.mean = DF.change %>% 
  group_by(time, population, Adjuvant) %>% 
  summarise(count.mean = mean(count, na.rm=T), count.sd = sd(count, na.rm=T)) %>% 
  ungroup()

ggplot(DF.change.mean %>% filter(population %in% pop.in), aes(time, count.mean, group=Adjuvant, col=Adjuvant)) +
  facet_wrap(~population, ncol=1, scales="free_y") +
  geom_line() +
  geom_errorbar(aes(ymin = count.mean-count.sd, ymax = count.mean+count.sd), width=0.2) +
  xlab("Time points") +
  # ylab("") +
  ylab("Cell count, fold-change from day 0") +
  theme_bw() +
  theme(strip.background = element_blank())

dn.fig = file.path(PROJECT_DIR, "FIGURES/Tfh")
dir.create(dn.fig, showWarnings = F)

ggsave(file.path(dn.fig, "LN_counts_change_average.png"), w=5.5, h=5)

DF.mean = DF %>% 
  group_by(time, population, Adjuvant) %>% 
  summarise(count.mean = mean(count, na.rm=T), count.sd = sd(count, na.rm=T)) %>% 
  ungroup()
ggplot(DF.mean %>% filter(population %in% pop.in), aes(time, count.mean, group=Adjuvant, col=Adjuvant)) +
  facet_wrap(~population, ncol=1, scales="free_y") +
  geom_line() +
  geom_errorbar(aes(ymin = count.mean-count.sd, ymax = count.mean+count.sd), width=0.2) +
  xlab("Time points") +
  # ylab("") +
  ylab("Cell counts") +
  theme_bw() +
  theme(strip.background = element_blank())
ggsave(file.path(dn.fig, "LN_counts_average.png"), w=5.5, h=5)

