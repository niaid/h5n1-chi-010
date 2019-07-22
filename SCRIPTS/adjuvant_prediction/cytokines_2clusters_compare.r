# Initialize.
source("SCRIPTS/0_initialize.r")
dn.in = file.path(PROJECT_DIR, "DATA_PROCESSED/Luminex/")
df.lumi.long = readRDS(file.path(dn.in, "luminex_data.rds"))

for (an in c("IP-10", "MIP-1b", "SAA", "CRP")) {
df.sel = df.lumi.long %>% 
  filter(analyte==an) %>% 
  mutate(subject = sub("H5N1-0", "", subject) %>% as.numeric() %>% paste0("s",.))

if(nrow(df.sel)==0) {
  cat(paste(an, "not in the dataset"))
  next
}

df.sel.fc = df.sel %>% 
  group_by(subject, analyte) %>% 
  summarize(`day 1 - day 0` = value[Time == 24] - value[Time == 0 & batch == 1],
            `day 22 - day 0` = value[Time == 528] - value[Time == 0 & batch == 2]) %>% 
  gather("time", "FC", -subject, -analyte)

dn.adj = file.path(PROJECT_DIR, "RESULTS/Adjuvant_prediction")
df.pred = fread(file.path(dn.adj, "2peaks_patterns_PCA.txt")) %>% 
  dplyr::select(subject, cluster) %>% 
  mutate(cluster = factor(cluster))

DF = df.sel.fc %>% inner_join(df.pred, by="subject")

df.w = DF %>% 
  group_by(time) %>% 
  summarize(p = wilcox.test(formula = FC ~ cluster, exact=F)$p.value) %>% 
  ungroup() %>% 
  mutate(x=0.5, y=max(DF$FC, na.rm=T), label=glue::glue("p = {format(p, digits=2)}"))
  

plot1 <- ggplot(DF, aes(cluster, FC, fill=cluster)) +
  geom_boxplot(alpha = 0.2) +
  geom_dotplot(binaxis = "y", stackdir = "center", fill = "black") +
  facet_wrap(~time) +
  geom_text(data=df.w, aes(x=x, y=y, label=label), col="black", hjust = 0, vjust = 1, inherit.aes = F) +
  ggtitle(an) +
  scale_y_log10() +
  scale_fill_manual(values = c("blue", "red")) +
  theme_bw() +
  theme(strip.background = element_blank(), legend.position="none")

dn.fig = file.path(PROJECT_DIR, "FIGURES/Luminex")
if(dir.exists(dn.fig)){
        print("Output directory already exists.")
}else{
        print("Output directory doesn't exists. Creating it now.")
        dir.create(dn.fig, recursive = T, showWarnings = T)
}
ggsave(file.path(dn.fig, glue::glue("2peak_clusters_compare_{an}.png")), plot = plot1, h=3, w=5)
}
