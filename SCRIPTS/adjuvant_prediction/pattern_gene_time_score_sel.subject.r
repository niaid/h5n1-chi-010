source(file.path(PROJECT_DIR, "SCRIPTS/functions/get_score.r"))
source(file.path(PROJECT_DIR, "SCRIPTS/functions/factor.date.r"))
dn = "DATA_PROCESSED/Microarrays/PBMC/samples.clean_genes.iqr"
fn = file.path(PROJECT_DIR, dn, "gexp_d0_fc.RData")
load(file=fn, verbose = T)

dn.patt = file.path(PROJECT_DIR, "RESULTS/Microarrays/PBMC/pattern_discovery")
patt.genes.clean = readRDS(file.path(dn.patt, "patt.genes.clean.rds"))

ss = "H5N1-020"
K = 1:3
df.score = data.frame()
for(k in K) {
X = dat.fc[rownames(dat.fc) %in% patt.genes.clean[[k]],] %>% get_score() %>% 
  as.data.frame() %>% 
  rename_(Score = ".") %>% 
  tibble::rownames_to_column("sample.name") %>%
  inner_join(info.fc, by="sample.name") %>%
  select(-sample.name, -plate.num, -iso.date, -n.d0) %>%
  group_by(subject.id, time.point) %>%
  summarise(Score = mean(Score,na.rm=T)) %>%
  ungroup() %>%
  # spread(time.point, fc) %>% 
  mutate(Pattern = sprintf("Gp%02d",k)) %>% 
  mutate(subject.id = sub("H5N1_0","s", subject.id)) %>% 
  filter(subject.id %in% ss)
df.score = rbind(df.score, X)
}

df.score = df.score %>% 
  mutate(time = gsub("(?<!1)0+(\\d)","\\1",time.point, perl=T)) %>% 
  mutate(time = sub("_0h","", time)) %>% 
  mutate(time = sub("_","+",time) %>% factor.date())

source(file.path(PROJECT_DIR, "SCRIPTS/functions/color_functions.r"))
cm = gg_color_hue(14)
ggplot(df.score, aes(time, Score, group=Pattern, col=Pattern)) +
  geom_line(size=1) +
  scale_color_manual(values = cm[K]) +
  xlab("Time") + ylab("Score") +
  geom_vline(xintercept = c(5,11), lty=2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size=6, angle = 90, hjust = 1, vjust=0.5),
    # axis.text.y = element_blank(),
    legend.position = c(0.5,1), legend.justification = c(0.5, 0.9),
    legend.title = element_blank(),
    legend.text = element_text(size=6),
    legend.key = element_rect(size = 4),
    legend.key.size = unit(0.5, 'lines'),
    legend.background = element_rect(fill = NA))

fn.fig = file.path(PROJECT_DIR, "FIGURES/pattern_GE_profiles_sel.subject")
ggsave(paste0(fn.fig, ".png"), w=2.5, h=2)
ggsave(paste0(fn.fig, ".pdf"), w=2.5, h=2)

