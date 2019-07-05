### read and process 15-color flow
source("SCRIPTS/0_initialize.r")
Sys.setlocale(locale="C")
fn.tfh15 = file.path(PROJECT_DIR, "DATA_PROCESSED/Flow_15c/H5N1-all-T4 021518.txt")
tfh15 = fread(fn.tfh15, header = T, sep = "\t", data.table=F) %>%
   dplyr::filter(!Sample %in% c("Mean","StdDev"))
tfh15$Sample = tfh15$`$FIL`
tfh15 = tfh15 %>% dplyr::select(Sample, matches("^Lym"))

# names(tfh15)[-1] = c("CXCR3+ Tfh","PD-1+ICOS+ CXCR3+ Tfh", 
#                      "Tfh2","PD-1+ICOS+ Tfh2", 
#                      "CXCR3- Tfh","PD-1+ICOS+ CXCR3- Tfh")
names(tfh15)[-1] = c("CXCR3+ Tfh","CXCR3+ Tfh activated, peak at day 7", 
                     "Tfh2","PD-1+ICOS+ Tfh2", 
                     "CXCR3- Tfh","CXCR3- Tfh activated, peak at day 28")
tfh15.levels = names(tfh15)[-1]

pop15 = tfh15.levels[c(2,6)]

tfh15 = tfh15 %>% 
  dplyr::filter(!grepl("H5N1-CHI", Sample)) %>% 
  dplyr::filter(!Sample %in% c("Mean","StdDev"))
tfh15 = tfh15 %>% 
  mutate(subject = str_extract(.$Sample, "H5N1-\\d{3}")) %>% 
  mutate(time.point = str_extract(.$Sample, "day \\d{1,2}") %>% sub("ay ","", .))
time15.levels = unique(tfh15$time.point)


df.tfh15 = tfh15 %>% 
  dplyr::select(subject, time.point, one_of(pop15)) %>%
  gather("population","pop", matches("Tfh"))

df.tfh15 = df.tfh15 %>% 
  group_by(subject, population) %>% 
  mutate(non0base = any(time.point=="d0" & pop>0)) %>% 
  ungroup() %>% 
  dplyr::filter(subject != "H5N1-044") %>%
  dplyr::filter(non0base)

df.tfh15.change = df.tfh15 %>% 
  group_by(subject, population) %>% 
  mutate(FC = pop / pop[time.point==time15.levels[1]]) %>% 
  mutate(delta = pop - pop[time.point==time15.levels[1]]) %>%
  dplyr::select(-pop, -delta)

### demographics and adjuvant status
fn.clin = file.path(PROJECT_DIR, "DATA_ORIGINAL/Clinical/clinical_info_adj.txt")
df.clin = fread(fn.clin) %>% 
  rename(subject=`Subject ID`)

df.comb = df.tfh15.change %>% 
  inner_join(df.clin %>% dplyr::select(subject, Adjuvant), by="subject")

df.change.mean = df.comb %>% 
  group_by(time.point, population, Adjuvant) %>% 
  summarise(FC.mean = mean(FC, na.rm=T), FC.sd = sd(FC, na.rm=T)) %>% 
  ungroup() %>% 
  mutate(time.point = factor(time.point, levels=time15.levels)) %>% 
  mutate(population = factor(population, levels = pop15))

ggplot(df.change.mean, aes(time.point, FC.mean, group=Adjuvant, col=Adjuvant)) +
  facet_wrap(~population, nrow=1, scales="fixed") +
  geom_line() +
  geom_errorbar(aes(ymin = FC.mean-FC.sd, ymax = FC.mean+FC.sd), width=0.2) +
  xlab("Time points") +
  ylab("% of parent, FC from d0") +
  theme_bw() +
  theme(strip.background = element_blank())

dn.fig = file.path(PROJECT_DIR, "FIGURES/Tfh")
dir.create(dn.fig, showWarnings = F)
fn.fig = file.path(dn.fig, "Tfh2_act_profiles_flow15_FC_average")
ggsave(paste0(fn.fig, ".png"), w=7, h=2.6)
ggsave(paste0(fn.fig, ".pdf"), w=7, h=2.6)

