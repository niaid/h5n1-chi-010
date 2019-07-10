source("SCRIPTS/0_initialize.r")
fn.pbmc = file.path(PROJECT_DIR, "RESULTS/Microarrays/PBMC/baseline/GE.pbmc_d0_WGCNA_ME_scores.txt")
df.pbmc = fread(fn.pbmc) %>% 
  mutate(subject = sub("_", "-", subject)) %>% 
  select(subject, Gb13)
fn.pax = file.path(PROJECT_DIR, "RESULTS/Microarrays/PAXgene/baseline/GE.pax_d0_WGCNA_ME_scores.txt")
df.pax = fread(fn.pax) %>% 
  select(subject, GbWB11)

# read demographics
fn.demo = file.path(PROJECT_DIR, "DATA_ORIGINAL", "Clinical", "clinical_info.txt")
df.demo = fread(fn.demo) %>% 
  dplyr::rename(subject  = `Subject ID`) %>%  
  mutate(Gender = ifelse(grepl("^F",Gender),"Female","Male"))

# read titers
fn.mn = file.path(PROJECT_DIR, "DATA_ORIGINAL", "Titres", "titre.txt")
df.mn = read.table(fn.mn, sep="\t", header=T, row.names=NULL, stringsAsFactors=T)
df.mn = df.mn %>% 
  mutate_at(-(1:2), function(x) sub("<20","10",x) %>% as.numeric()) %>% 
  mutate(day = as.numeric(sub("[dD]","",TimePt))) %>% 
  mutate(day = ifelse(day %in% c(29,31), 28, day) %>% factor()) %>% 
  mutate(subject = sprintf("H5N1-%03d",Sample.ID)) %>% 
  filter(day=="42")


# fn.hai = file.path(PROJECT_DIR, "DATA_ORIGINAL", "HAI", "H5N1 serology.txt")
# df.hai = read.table(fn.hai, sep="\t", 
#                     header=T, row.names=NULL, stringsAsFactors=F) %>% 
#   mutate_at(-(1), function(x) sub("<10","5",x) %>% as.numeric()) %>% 
#   gather("tmp","hai",-subject) %>% 
#   separate("tmp",c("day","replicate"))
# df.hai3 = df.hai %>% 
#   group_by(subject, day) %>% 
#   summarise(hai.median=median(hai, na.rm=T), 
#             hai.min=min(hai, na.rm=T), hai.max=max(hai, na.rm=T) ) %>% 
#   ungroup() %>% 
#   mutate(day=factor(day,levels=c("day0","day1","day7","day21","day22","day28","day42","day100")))


# DF = inner_join(df.pax, df.pbmc, by="subject") %>% 
#   inner_join(df.hai3 %>% filter(day=="day28") %>% rename(HAI=hai.median), by="subject") %>% 
#   inner_join(df.demo, by="subject")
DF = inner_join(df.pax, df.pbmc, by="subject") %>% 
  inner_join(df.mn, by="subject") %>% 
  inner_join(df.demo, by="subject")

ggplot(DF, aes(Gb13, GbWB11, col=Gender, size=A.Indonesia)) +
  geom_point(alpha=0.8) +
  # scale_color_manual(values=c(Male="cyan",Female="pink")) +
  scale_size_continuous(range = c(2,10), breaks = c(20,80,320,1280, 5120)) +
  theme_bw() +
  theme(legend.position = c(0.85,0.4), legend.box.background = element_blank())
fn.fig = file.path(PROJECT_DIR, "FIGURES/Baseline/Gb13_vs_GbWB11.png")
ggsave(fn.fig, w=5, h=4)
