library(tidyverse)
Ng = 100
Nd = 10
Nt = 8
Ns = Nd * Nt
Ndg = Nd * Ng
d_lbl = paste0("d", 1:Nd)
g_lbl = paste0("g", 1:Ng)
t_lbl = paste0("t", 1:Nt)

s_df = expand.grid(d_lbl, t_lbl)
names(s_df) = c("d","t")
s_df = s_df %>% mutate(s = paste(d, t, sep="_"))
s_lbl = s_df$s

sg_df = expand.grid(s_lbl, g_lbl)
names(sg_df) = c("s","g")
sg_df = sg_df %>% 
  separate(s, c("d","t"), sep = "_", remove = F) %>% 
  mutate(sg = paste(s, g, sep="__")) %>% 
  mutate(dg = paste(d, g, sep="__"))


gg1 = paste0("g", 1:floor(Ng/2))
gg2 = paste0("g", 1:floor(Ng/1.9))
gg3 = paste0("g", floor(Ng/4*3):Ng)
dg1 = paste0("d", 1:floor(Nd/2))
dg2 = paste0("d", (floor(Nd/3)+1):floor(Nd/3*2))
dg3 = d_lbl
tg1 = paste0("t", c(3, 6))
tg2 = paste0("t", c(2,3, 5,6))
tg3 = paste0("t", c(2:3, 5:7))

sg_df = sg_df %>% 
  mutate(v = case_when(
    d %in% dg1 & g %in% gg1 & t %in% tg1 ~ 1,
    d %in% dg2 & g %in% gg2 & t %in% tg2 ~ 1,
    d %in% dg3 & g %in% gg3 & t %in% tg3 ~ -1,
    TRUE ~ 0
  ))

dg_t_df = sg_df %>%
  select(-s, -sg) %>%
  spread(t, v)

dat = dg_t_df %>% 
  select(-d,-g) %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames("dg") %>% 
  data.matrix()

rm_prc = 0.1
ridx = sample(length(dat), length(dat) * rm_prc, replace = F)
dat[ridx] = 0
dat = dat + rnorm(length(dat), 0, 0.2)

library(ComplexHeatmap)
library(circlize)
K = 4
hm = Heatmap(dat, name="hm1", cluster_columns = F, cluster_rows = T, km = K, 
             col = colorRamp2(c(-1, 0, 1), c("blue","white","red")),
             show_row_names = F, show_column_names = F, show_heatmap_legend = F,
             show_row_dend = F, #km_title = "p%i",
             gap = unit(1, "mm"),
             # rect_gp = gpar(fill = "transparent", lwd=0),
             row_title = "", column_title = "",
             column_title_side = "bottom", row_title_side = "left")
# decorate_heatmap_body("hm1", {
#   grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 2))
# })

# for (ss in 1:100) {
#   set.seed(ss)
#   draw(hm)
#   cat(ss,"")
#   readline(prompt="Press [enter] to continue")
# }
ss = 16
set.seed(ss)
draw(hm)

dev.copy(png, "pattern_sim_1.png", w=100, h=150)
dev.off()
dev.copy(pdf, "pattern_sim_1.pdf", w=2, h=3)
dev.off()

set.seed(ss)
ro = row_order(hm)
# ng = sapply(ro, length)
# iord = order(ng, decreasing = T)
# ro = ro[iord]

dg_t_df$pattern = 0
for(k in 1:K) {
  set.seed(8)
  dg_t_df$pattern[ro[[k]]] = k-1
}

# Heatmap(dat[dg_t_df$pattern==4, ], cluster_columns = F, cluster_rows = T,
#              show_row_names = F, show_column_names = T)
# 

p_df = dg_t_df %>% 
  select(d, g, pattern) %>% 
  spread(d, pattern)
p_dat = p_df %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames("g") %>% 
  data.matrix()

# source("functions/color_functions.r")
# cm = gg_color_hue(K-1)
cm = c("white", "tomato", "steelblue", "gold")
set.seed(ss)
Heatmap(t(p_dat), name="hm2", col = cm,
        show_row_dend = F, show_column_dend = F, 
        show_row_names = F, show_column_names = F, show_heatmap_legend = F,
        row_title = "", column_title = "",
        column_title_side = "bottom")
decorate_heatmap_body("hm2", {
  grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 1))
})
dev.copy(png, "pattern_sim_2.png", w=200, h=200)
dev.off()
dev.copy(pdf, "pattern_sim_2.pdf", w=4, h=4)
dev.off()


p_t_df = dg_t_df %>% 
  group_by(pattern) %>% 
  summarize_at(vars(matches("^t")), mean) %>% 
  ungroup() %>% 
  filter(pattern!=0)

pdat = p_t_df %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames("pattern") %>% 
  data.matrix()

# xx = cor(t(dat), t(pdat))
# matplot(t(pdat), type = "l")

p_g_df = dg_t_df %>%
  select(pattern, g) %>% 
  filter(pattern!=0) %>% 
  group_by(pattern) %>% 
  distinct() %>% 
  ungroup()
glist = split(p_g_df$g, p_g_df$pattern)

dat_df = as.data.frame(dat) %>% 
  tibble::rownames_to_column("dg") %>% 
  gather("t", "v", -dg) %>% 
  separate(dg, c("d","g"), sep="__", remove = T)

edat = dat_df %>% 
  unite("s", d, t, sep="_") %>% 
  spread(s, v) %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames("g") %>% 
  data.matrix()

K = length(glist)
sscore = matrix(NA, nrow=Ns, ncol=K)
colnames(sscore) = paste0("p",1:K)
for(k in 1:K) {
    gi = rownames(edat) %in% glist[[k]] 
    sscore[,k] = WGCNA::moduleEigengenes(t(edat[gi,]), rep(1, sum(gi)))$eigengenes %>% pull(1)
}
rownames(sscore) = colnames(edat)


dscore = matrix(NA, nrow=Nd, ncol=K)
colnames(dscore) = paste0("p",1:K)
rownames(dscore) = d_lbl
for(k in 1:K) {
  ss = sscore[,k]
  tmp = as.data.frame(ss) %>% 
    tibble::rownames_to_column("s") %>% 
    separate(s, c("d","t"), sep="_", remove = T) %>% 
    mutate(d = factor(d, levels=d_lbl)) %>% 
    spread(t, ss) %>% 
    tibble::remove_rownames() %>% 
    tibble::column_to_rownames("d") %>% 
    data.matrix()
  dscore[,k] = cor(t(tmp), pdat[k,], use = "pairwise.complete.obs")
}

Heatmap(dscore, name="score", cluster_columns = F,
        col = colorRamp2(c(-1, 0, 1), c("blue","white","red")),
        show_row_dend = F, show_column_dend = F, 
        show_row_names = F, show_column_names = T, show_heatmap_legend = F,
        row_title = "", column_title = "",
        column_title_side = "bottom")
decorate_heatmap_body("score", {
  grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 1))
})
dev.copy(png, "pattern_sim_3.png", w=200, h=200)
dev.off()
dev.copy(pdf, "pattern_sim_3.pdf", w=4, h=4)
dev.off()

pg_dat = p_g_df %>% 
  mutate(v = pattern) %>% 
  mutate(pattern = paste0("p", pattern)) %>% 
  spread(pattern, v, fill = 0) %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames("g") %>% 
  data.matrix()

Heatmap(pg_dat, name="pg", cluster_columns = F,
        col = cm,
        show_row_dend = F, show_column_dend = F, 
        show_row_names = F, show_column_names = T, show_heatmap_legend = F,
        row_title = "", column_title = "",
        column_title_side = "bottom")
decorate_heatmap_body("pg", {
  grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 1))
})
dev.copy(png, "pattern_sim_4.png", w=200, h=200)
dev.off()
dev.copy(pdf, "pattern_sim_4.pdf", w=4, h=4)
dev.off()
