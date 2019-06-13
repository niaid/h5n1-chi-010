dn.in = file.path(PROJECT_DIR, "DATA_PROCESSED/Luminex/")
df.lumi.long = readRDS(file.path(dn.in, "luminex_data.rds"))

an = "IP-10"
df.ip10 = df.lumi.long %>% 
  filter(analyte==an) %>% 
  mutate(subject = sub("H5N1-0", "", subject) %>% as.numeric() %>% paste0("s",.))

df.ip10.fc = df.ip10 %>% 
  group_by(subject, analyte) %>% 
  summarize(`day 1 - day 0` = value[Time == 24] - value[Time == 0 & batch == 1],
            `day 22 - day 0` = value[Time == 528] - value[Time == 0 & batch == 2]) %>% 
  gather("time", "FC", -subject, -analyte)

dn.adj = file.path(PROJECT_DIR, "RESULTS/Adjuvant_prediction")
df.pred = fread(file.path(dn.adj, "2peaks_patterns_PCA.txt")) %>% 
  dplyr::select(subject, cluster) %>% 
  mutate(cluster = factor(cluster))

DF = df.ip10.fc %>% inner_join(df.pred, by="subject")

df.w = DF %>% 
  group_by(time) %>% 
  summarize(p = wilcox.test(formula = FC ~ cluster, exact=F)$p.value) %>% 
  ungroup() %>% 
  mutate(x=0.5, y=max(DF$FC), label=glue::glue("p = {format(p, digits=2)}"))
  

ggplot(DF, aes(cluster, FC, fill=cluster)) +
  geom_boxplot(alpha = 0.2) +
  geom_dotplot(binaxis = "y", stackdir = "center", fill = "black") +
  facet_wrap(~time) +
  geom_text(data=df.w, aes(x=x, y=y, label=label), col="black", hjust = 0, vjust = 1, inherit.aes = F) +
  scale_y_log10() +
  scale_fill_manual(values = c("blue", "red")) +
  theme_bw() +
  theme(strip.background = element_blank())

dn.fig = file.path(PROJECT_DIR, "FIGURES/Luminex")
dir.create(dn.fig, showWarnings = F)
ggsave(file.path(dn.fig, glue::glue("2peak_clusters_compare_{an}.png")))
