# was plot_flow_patters_only_171123.r

dn.in = file.path(PROJECT_DIR, "RESULTS/Flow_10c")

# Flow modules
fn = "DeltaThr0.5_DonorThr0.3_PercentOfParent_CORpearson_cutreeHybrid_deepSplit0_minClusterSize31_Modules_v2.txt"
patt.flow = read.csv(file.path(dn.in, fn), row.names=1, stringsAsFactors = F)
patt.flow = t(patt.flow)
rownames(patt.flow) = sub("Module.","Fp", rownames(patt.flow))

K = nrow(patt.flow)
patt.flow = cbind(rep(0.0,K), patt.flow)
colnames(patt.flow)[1] = "d0_0h"

saveRDS(patt.flow, file.path(dn.in, "patt.flow.rds"))

# patt.flow = readRDS(file.path(dn.in, "patt.flow.rds"))
patt.flow.df = patt.flow %>% as.data.frame() %>% 
  tibble::rownames_to_column("Module") %>% 
  gather("Time", "Score", -Module) %>% 
  # mutate(Module = sub("M.","Fp0",Module)) %>% 
  mutate(Time = factor(Time, level=colnames(patt.flow)))

ggplot(patt.flow.df, aes(Time, Score, group=Module, col=Module)) +
  geom_line(size=1) +
  facet_wrap(~Module, ncol=1, strip.position="right") +
  ylab("Score") + xlab("") +
  geom_hline(yintercept=0, lty=2) +
  geom_vline(xintercept = c(4,9), lty=2, col="red") +
  geom_vline(xintercept = c(5,10), lty=2, col="blue") +
  scale_color_manual(values = palette()) +
  theme_bw() + guides(col=FALSE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size=6),
        strip.text.y = element_text(size = 10, angle = 90),
        strip.background = element_blank())

fn.fig = file.path(PROJECT_DIR, "FIGURES/Flow_patterns_profiles")
ggsave(paste0(fn.fig, ".png"), h=8.5,w=2)
ggsave(paste0(fn.fig, ".pdf"), h=8.5,w=2)

ggplot(patt.flow.df, aes(Time, Score, group=Module, col=Module)) +
  geom_line(size=1) +
  facet_wrap(~Module, ncol=2, strip.position="top") +
  ylab("Score") + xlab("") +
  geom_hline(yintercept=0, lty=2) +
  geom_vline(xintercept = c(4,9), lty=2, col="red") +
  geom_vline(xintercept = c(5,10), lty=2, col="blue") +
  scale_color_manual(values = palette()) +
  theme_bw() + guides(col=FALSE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size=6),
        strip.text.y = element_text(size = 10, angle = 90),
        strip.background = element_blank())

fn.fig = file.path(PROJECT_DIR, "FIGURES/Flow_patterns_profiles_2col")
ggsave(paste0(fn.fig, ".png"), h=5,w=3)
ggsave(paste0(fn.fig, ".pdf"), h=5,w=3)

ggplot(patt.flow.df, aes(Time, Score, group=Module, col=Module)) +
  geom_line(size=1) +
  facet_wrap(~Module, nrow=1, strip.position="top") +
  ylab("Score") + xlab("") +
  geom_hline(yintercept=0, lty=2) +
  geom_vline(xintercept = c(4,9), lty=2, col="red") +
  geom_vline(xintercept = c(5,10), lty=2, col="blue") +
  scale_color_manual(values = palette()) +
  theme_bw() + guides(col=FALSE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size=6),
        strip.text.y = element_text(size = 10, angle = 0),
        strip.background = element_blank())
fn.fig = file.path(PROJECT_DIR, "FIGURES/Flow_patterns_profiles_horiz.png")
ggsave(fn.fig, h=1.7,w=10)

