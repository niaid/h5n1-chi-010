source("SCRIPTS/0_initialize.r")
source("SCRIPTS/functions/factor.date.r")
### demographics and adjuvant status
fn.clin = file.path(PROJECT_DIR, "DATA_ORIGINAL/Clinical/clinical_info_adj.txt")
df.clin = fread(fn.clin) %>% 
  rename(subject=`Subject ID`)


### read titer d28
fn.mn = file.path(PROJECT_DIR, "DATA_ORIGINAL", "Titres", "titre.txt")
df.mn = read.table(fn.mn, sep="\t", 
                   header=T, row.names=NULL, stringsAsFactors=F) %>% 
  mutate_at(-(1:2), function(x) sub("<20","10",x) %>% as.numeric()) %>% 
  # mutate(day = as.numeric(sub("[dD]","",TimePt))) %>% 
  mutate(day = sub("D","d",TimePt)) %>% 
  mutate(day = ifelse(day %in% c("d29","d31"), "d28", day) %>% factor.date()) %>% 
  mutate(subject = sprintf("H5N1-%03d",Sample.ID)) %>% 
  filter(day=="d28") %>% 
  select(subject, titer.d28=A.Indonesia) %>% 
  mutate(titer.d28.log2 = log2(titer.d28))


### read and process 15-color flow

fn.tfh15 = file.path(PROJECT_DIR, "DATA_PROCESSED/Flow_15c/TFh activated.022218.txt")
Sys.setlocale(locale="C")
tfh15 = fread(fn.tfh15, data.table=F)
colnames(tfh15)[ncol(tfh15)] = "subject"
colnames(tfh15)[ncol(tfh15)-1] = "time.point"
tfh15 = tfh15 %>% 
  filter(grepl("^\\d{3}$", subject)) %>% 
  mutate(subject = paste0("H5N1-", subject)) %>% 
  mutate(time.point = sub("ay ","", time.point)) %>% 
  mutate(Sample = glue::glue("{subject}_{time.point}"))
tfh15 = tfh15 %>% select(-matches("Freq. of Total"))
names(tfh15)[2] = "PD-1+ICOS+ Tfh"

tfh15.levels = names(tfh15)[2]

dn.fig = file.path(PROJECT_DIR, "FIGURES/Tfh")
dir.create(dn.fig, showWarnings = F)

# pop15 = "PD-1+ICOS+CXCR3+ Tfh"
for (k in c(1)) {
  pop15 = tfh15.levels[k]
  pop.parent = "total Tfh"
  
change.type = "FC"

# tfh15 = tfh15 %>% 
#   filter(!grepl("H5N1-CHI", Sample))
# tfh15 = tfh15 %>% 
#   mutate(subject = str_extract(.$Sample, "H5N1-\\d{3}")) %>% 
#   mutate(time.point = str_extract(.$Sample, "day \\d{1,2}") %>% sub("ay ","", .))
time15.levels = unique(tfh15$time.point)
# table(tfh15$time.point, useNA="always")

df.tfh15 = tfh15 %>% 
  select(subject, time.point, one_of(pop15)) %>%
  # select(-Sample) %>% 
  gather("population","pop", matches("Tfh"))

subj = df.tfh15$subject[df.tfh15$time.point=="d0" & df.tfh15$pop>0] %>% unique()
df.tfh15 = df.tfh15 %>% 
  filter(subject %in% subj)
# check dups
# table(paste(df.tfh15$subject, df.tfh15$time.point, df.tfh15$population)) %>% unique()

df.tfh15.change = df.tfh15 %>% 
  # group_by(subject,time.point,population) %>% mutate(n=n()) %>% ungroup() %>% filter(n==1) %>% select(-n) %>% #remove duplicates
  group_by(subject, population) %>% 
  mutate(FC = pop / pop[time.point==time15.levels[1]]) %>% 
  mutate(delta = pop - pop[time.point==time15.levels[1]]) %>%
  gather(type, value, c("pop","delta","FC")) %>% 
  filter(type == change.type)
  # select(one_of("subject","time.point","population",change.type))

df.comb = df.tfh15.change %>% 
  inner_join(df.mn, by="subject") %>% 
  inner_join(df.clin %>% select(subject, Adjuvant), by="subject")

df.change.mean = df.comb %>% 
  group_by(time.point, population, Adjuvant) %>% 
  summarise(value.mean = mean(value, na.rm=T), value.sd = sd(value, na.rm=T)) %>% 
  ungroup() %>% 
  mutate(time.point = factor(time.point, levels=time15.levels))

ggplot(df.change.mean, aes(time.point, value.mean, group=Adjuvant, col=Adjuvant)) +
  # facet_wrap(~population, ncol=1, scales="free_y") +
  geom_line() +
  geom_errorbar(aes(ymin = value.mean-value.sd/2, ymax = value.mean+value.sd/2), width=0.2) +
  xlab("Time points") +
  ylab(glue::glue("% of {pop.parent}, {change.type} from d0")) +
  theme_bw() +
  theme(strip.background = element_blank())

fn.fig = file.path(dn.fig, glue::glue("{pop15}_profiles_flow15_{change.type}_average_in"))
ggsave(paste0(fn.fig, ".png"), w=5, h=2.6)
ggsave(paste0(fn.fig, ".pdf"), w=5, h=2.6)
}
