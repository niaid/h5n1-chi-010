patt.genes.clean = readRDS("patt.genes.clean_200.2000.4.cleaned_170210.rds")

k = 2

X = dat.fc[rownames(dat.fc) %in% patt.genes.clean[[k]],] %>% #get_score() %>% 
  as.data.frame()
X = X %>% tibble::rownames_to_column("gene") %>%
  gather("sample.name","fc",-gene) %>% 
  inner_join(info.fc, by="sample.name") %>%
  select(-sample.name, -plate.num, -iso.date, -n.d0) %>%
  # filter(time.point!="d000_00h") %>%
  group_by(gene, subject.id, time.point) %>%
  summarise(fc = mean(fc, na.rm=T)) %>%
  ungroup() %>%
  # spread(time.point, fc) %>% 
  mutate(pattern=k) %>% 
  mutate(subject = sub("H5N1_0","s", subject.id))

subj.patt.cor = readRDS("subj.patt.cor.incl.10_200.2000.4.cleaned_170420.rds")
iord = order(subj.patt.cor[,k])
s.min = rownames(subj.patt.cor)[iord[1]]
s.max = rownames(subj.patt.cor)[rev(iord)[1]]

s = "s11"

X %>% filter(subject==s, time.point!="d000_00h") %>% 
ggplot(aes(time.point, gene, fill=fc)) +
  geom_tile() +
  scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=0, limits=c(-2,2), guide = F) +
  xlab("Time") + ylab("Genes") +
  ggtitle(s) +
  theme_bw() +
  theme(#axis.text.x = element_blank(),
        axis.text.x = element_text(size=6, angle = 90, hjust = 1, vjust=0.5),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank())

fn.fig = sprintf("pattern_gene_time_heatmap_Gp%d_%s_labels", k, s)
ggsave(paste0(fn.fig, ".png"), w=2, h=3)
