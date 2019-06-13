# was plot_patterns_profiles_171123.r (part 2)

dn.patt = file.path(PROJECT_DIR, "RESULTS/Microarrays/PBMC/pattern_discovery")
# df.mat = readRDS(file.path(dn.patt, "df.mat.rds"))

GE.patterns = fread(file.path(dn.patt, "GE_patterns_filt.txt"))
ct = GE.patterns$pattern %>% sub("Gp","",.) %>% as.numeric()
names(ct) = GE.patterns$label
ct.lev = sort(unique(GE.patterns$pattern))
K=max(ct)

df.cl.stat = readRDS(file.path(dn.patt, "df.cl.stat.rds")) %>% 
  filter(method=="mean") %>% 
  # mutate(ct = sub("G","",ct) %>% as.numeric() %>% sprintf("Gp%02d",.)) %>% 
  # mutate(ct = factor(ct, levels=ct.lev)) %>% 
  mutate(time = gsub("(?<!1)0+(\\d)","\\1",time, perl=T)) %>% 
  mutate(time = sub("_0h","", time)) %>% 
  mutate(time = sub("_","+",time) %>% factor(levels=unique(.)))


# zl.all = readRDS(file.path(dn.patt, "zl.all.f_170126.rds")
# ct = zl.all[,"200.2000.4.cleaned"]
# ct = ct[ct!=0]
# ct=as.numeric(factor(ct))

df.n = as.data.frame(table(GE.patterns$pattern)) %>% 
  dplyr::rename(ct=Var1) %>% 
  mutate(x=1,y=Inf,label=sprintf("%d",Freq))

ggplot(df.cl.stat, aes(time, value, group=ct,col=ct)) + geom_line(size=1) +
  # geom_errorbar(aes(ymin=mean-sd/2, ymax=mean+sd/2), width=0.2) +
  # scale_color_manual(values=cm) + ylab("log-FC") +
  geom_hline(yintercept = 0, lty=2) +
  geom_vline(xintercept = c(5,11), lty=2, col="red") +
  geom_vline(xintercept = c(6,12), lty=2, col="blue") +
  facet_wrap(~ct, ncol=1, strip.position="right") +
  ylab("Mean logFC") + 
  # xlab("Time after first vaccination") +
  xlab("") +
  geom_text(data=df.n, aes(x,y,label=label), hjust=0, vjust=1.2, col="black", size=3, inherit.aes = F) +
  guides(col="none") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size=6),
        strip.text.y = element_text(size = 10, angle = 90),
        strip.background = element_blank())
# 
ggsave(file.path(PROJECT_DIR, "FIGURES", "GE_patterns_profiles.png"), h=14,w=2)

ggplot(df.cl.stat, aes(time, value, group=ct,col=ct)) + geom_line(size=1) +
  # geom_errorbar(aes(ymin=mean-sd/2, ymax=mean+sd/2), width=0.2) +
  # scale_color_manual(values=cm) + ylab("log-FC") +
  geom_hline(yintercept = 0, lty=2) +
  geom_vline(xintercept = c(5,11), lty=2, col="red") +
  geom_vline(xintercept = c(6,12), lty=2, col="blue") +
  facet_wrap(~ct, ncol=2, strip.position="top") +
  ylab("Mean logFC") + 
  # xlab("Time after first vaccination") +
  xlab("") +
  geom_text(data=df.n, aes(x,y,label=label), hjust=0, vjust=1.2, col="black", size=3, inherit.aes = F) +
  guides(col="none") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size=6),
        strip.text.y = element_text(size = 10, angle = 90),
        strip.background = element_blank())
# 
ggsave(file.path(PROJECT_DIR, "FIGURES", "GE_patterns_profiles_2col.png"), h=9,w=3)

ggplot(df.cl.stat, aes(time, value, group=ct,col=ct)) + geom_line(size=1) +
  # geom_errorbar(aes(ymin=mean-sd/2, ymax=mean+sd/2), width=0.2) +
  # scale_color_manual(values=cm) + ylab("log-FC") +
  geom_hline(yintercept = 0, lty=2) +
  geom_vline(xintercept = c(5,11), lty=2, col="red") +
  geom_vline(xintercept = c(6,12), lty=2, col="blue") +
  facet_wrap(~ct, nrow=2, strip.position="top") +
  ylab("Mean logFC") + 
  # xlab("Time after first vaccination") +
  xlab("") +
  geom_text(data=df.n, aes(x,y,label=label), hjust=0, vjust=1.2, col="black", size=3, inherit.aes = F) +
  guides(col="none") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size=6),
        strip.text.y = element_text(size = 10, angle = 90),
        strip.background = element_blank())
# 
ggsave(file.path(PROJECT_DIR, "FIGURES", "GE_patterns_profiles_horiz.png"), h=3,w=10)

