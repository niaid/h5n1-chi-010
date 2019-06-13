source(file.path(PROJECT_DIR, "SCRIPTS/functions/factor.date.r"))

fn.mn = file.path(PROJECT_DIR, "DATA_ORIGINAL", "Titres", "titre.txt")
df.mn = read.table(fn.mn, sep="\t", header=T, row.names=NULL, stringsAsFactors=T) %>% 
  mutate_at(-(1:2), function(x) sub("<20","10",x) %>% as.numeric()) %>% 
  # mutate(day = as.numeric(sub("[dD]","",TimePt))) %>% 
  mutate(day = sub("D","d",TimePt)) %>% 
  mutate(day = ifelse(day %in% c("d29","d31"), "d28", day) %>% factor.date(), subject=factor(Sample.ID))

yticks = sort(unique(df.mn$A.Indonesia))
yticks = yticks[seq(2,length(yticks),2)]

# check subjects with peak at d42
tm = levels(df.mn$day)

### demographics and adjuvant status
fn.demo = file.path(PROJECT_DIR, "DATA_ORIGINAL/Clinical/clinical_info_adj.txt")
df.demo = fread(fn.demo, data.table = F) %>% 
  mutate(Sample.ID = sub("H5N1-0","", `Subject ID`) %>% as.numeric())


DF = df.mn %>% 
  group_by(Sample.ID) %>% 
  summarize(peak.time = tm[which.max(A.Indonesia)], 
            at.peak = A.Indonesia[which.max(A.Indonesia)],
            after.peak = A.Indonesia[which.max(A.Indonesia)+1]) %>% 
  mutate(decline = (at.peak-after.peak)/at.peak * 100) %>% 
  left_join(df.demo, by="Sample.ID") %>% 
  mutate(peak.time = factor(peak.time, levels = tm))

DF %>% filter(peak.time=="d42")

p1 = ggplot(DF %>% mutate(peak.time=factor(peak.time, levels=tm)), 
       aes(peak.time, group=Adjuvant, fill=Adjuvant)) +
  geom_bar() +
  theme_bw()

# ggplot(DF, aes(log2(at.peak), group=Adjuvant, fill=Adjuvant)) +
#   geom_histogram(position="dodge") +
#   scale_x_continuous(name = "MN titer at peak", breaks=log2(yticks), labels = yticks) +
#   theme_bw()

# check decline
# ggplot(DF %>% filter(at.peak > 20), aes(decline, group=Adjuvant, fill=Adjuvant)) +
#   geom_histogram(position="dodge") +
#   xlab("decline after peak, %") +
#   theme_bw()

p2 = ggplot(DF, aes(log2(at.peak), decline, col=Adjuvant)) +
  geom_jitter(width=0.2, height=2, size=3) +
  scale_x_continuous(name = "MN titer at peak", breaks=log2(yticks), labels = yticks) +
  ylab("decline after peak, %") +
  theme_bw()

library(gridExtra)
mg = marrangeGrob(list(p1,p2), nrow=1, ncol=2, top = NULL)

dn.fig = file.path(PROJECT_DIR, "FIGURES/titers")
dir.create(dn.fig, showWarnings = F)

ggsave(file.path(dn.fig, "MN_titer_peak_decline.png"), plot=mg, w=7, h=2.5)

# output table to text
df.out = df.mn %>% 
  group_by(Sample.ID) %>% 
  summarize(peak.time = tm[which.max(A.Indonesia)], 
            at.peak = A.Indonesia[which.max(A.Indonesia)],
            after.peak = A.Indonesia[which.max(A.Indonesia)+1]) %>% 
  mutate(decline = (at.peak-after.peak)/at.peak * 100) %>% 
  mutate(Sample.ID = sprintf("H5N1-%03d", Sample.ID)) %>% 
  dplyr::rename(subject=Sample.ID)

dn.out = file.path(PROJECT_DIR, "RESULTS/titers")
dir.create(dn.out, showWarnings = F)

fwrite(df.out, file.path(dn.out, "MN_titer_peak_time.txt"), sep="\t")
