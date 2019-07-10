# was pattern_genes_170518.r

# library(gplots)
# library(cluster)
source("SCRIPTS/0_initialize.r")
load(file.path(PROJECT_DIR, "DATA_PROCESSED/Microarrays/PBMC/samples.clean_genes.iqr", "gexp_d0_fc.RData"), verbose = T)


fc.th = 0.5
q.th = 0.2
tp = unique(info.fc$time.point)
# tp = c("d001")#,"d022")
fsamp = matrix(nrow=nrow(dat.fc), ncol=length(tp)*2)
rownames(fsamp) = rownames(dat.fc)
colnames(fsamp) = paste(rep(tp,1,each=2), c("pos","neg"), sep=".")
gi = rep(FALSE, nrow(dat.fc))
k = 1
for (t in tp) {
  ti = info.fc$time.point==t
  fsamp[,k] = apply(dat.fc[,ti], 1, function(x) sum(x > fc.th)/length(x))
  fsamp[,k+1] = apply(dat.fc[,ti], 1, function(x) sum(x < -fc.th)/length(x))
  gi = gi | fsamp[,k] >= q.th | fsamp[,k+1] >= q.th
  k = k + 2
}
sum(gi) # 6672 genes

df = as.data.frame(dat.fc[gi,]) %>%
  tibble::rownames_to_column("gene") %>%
  gather("sample.name","fc",-gene) %>%
  inner_join(info.fc, by="sample.name") %>%
  select(-sample.name, -plate.num, -iso.date, -n.d0) %>%
  # filter(time.point!="d000_00h") %>%
  group_by(gene,subject.id, time.point) %>%
  summarise(fc = mean(fc,na.rm=T)) %>%
  ungroup() %>%
  spread(time.point, fc)

df.mat.2 = as.matrix(df[,-(1:2)])
rownames(df.mat.2) = paste(df$subject.id, df$gene, sep="__")

cl.methods = c()
cl.methods[2] = "200.2000.4.cleaned"

# data from pattern_discovery
# set.seed(123)
# K = 8
# km = kmeans(df.mat.cc, K)
# ct = km$cluster
# cat(table(ct))
# ct.use = which(table(ct)>=20)
# print(ct.use)


# df.cl = as.data.frame(ct) %>%
#   tibble::rownames_to_column("subj_gene") %>%
#   separate(subj_gene, c("subject","gene"), sep="__") %>%
#   spread(subject, ct, fill = 0)
# df.cl.mat = as.matrix(df.cl[,-1])
# rownames(df.cl.mat) = df.cl$gene
# 
# df.cl.mean = as.data.frame(df.mat) %>%
#   mutate(ct=ct) %>%
#   group_by(ct) %>%
#   summarise_all(mean, na.rm=T) %>%
#   ungroup() %>%
#   gather("time","mean",-1)
# df.cl.sd = as.data.frame(df.mat) %>%
#   mutate(ct=ct) %>%
#   group_by(ct) %>%
#   summarise_all(sd, na.rm=T) %>%
#   ungroup() %>%
#   gather("time","sd",-1)
# df.cl.stat = cbind(df.cl.mean, df.cl.sd[3]) %>%
#   mutate(ct=factor(ct))
# cm = rainbow(K)
# if(sum(df.cl.mat==0)>0) cm = c("#AAAAAAFF",cm)
# ggplot(df.cl.stat, aes(time, mean, group=ct, col=ct)) + geom_line(size=1) + 
#   geom_errorbar(aes(ymin=mean-sd/2, ymax=mean+sd/2), width=0.2) +
#   scale_color_manual(values=cm) + ylab("log-FC") +
#   geom_hline(yintercept = 0, lty=2) +
#   facet_wrap(~ct) + 
#   theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust = 1))
# ggsave(sprintf("km%d_clusters_lines_%s_fc.%.1f_q%.1f_facets.png",K,fc.dir,fc.th,q.th), w=8,h=8)

# df.patt = df.cl.mean %>%
#   filter(ct!=0) %>% 
#   spread(time, mean)
# df.patt.mat = as.matrix(df.patt[,-1])
# rownames(df.patt.mat) = sprintf("pat%02d",ct.use)

# library(psych)
today.dt=today()
for (m in seq(2,8,by=2)) {
  # m = 2
  cat(m,"")
  df.patt.mat = readRDS(sprintf("df.patt_%s_170208.rds", cl.methods[m]))
  # df.patt.mat = readRDS("df.patt_cluster1_161221.rds")
  
  K = nrow(df.patt.mat)
  rownames(df.patt.mat) = sprintf("pat%02d",1:K)

# LONG CALCULATIONS!!! Use saved data if possible
  long.idx = c(rep(1,3),2:ncol(df.mat.2))
  df.patt.cor = cor(t(df.mat.2[,long.idx]), t(df.patt.mat[,long.idx]), use="pairwise.complete.obs")
  # df.patt.p = corr.test(t(df.mat.2[,long.idx]), t(df.patt.mat[,long.idx]), use="pairwise")$p

  df.patt.cor.p = df.patt.cor
  for (j in 1:nrow(df.patt.mat)) {
    cat(j,"")
    for (i in 1:nrow(df.mat.2)) {
      df.patt.cor.p[i,j] = cor.test(df.mat.2[i,long.idx], df.patt.mat[j,long.idx], alternative = "greater", exact=F)$p.value
    }
  }

  saveRDS(df.patt.cor, sprintf("df.cor_%dg_%s_%s.rds",sum(gi), cl.methods[m], today.dt))
  saveRDS(df.patt.cor.p, sprintf("df.cor.p_%dg_%s_%s.rds",sum(gi), cl.methods[m], today.dt))
}

# for (m in 1:2)
m = 2
df.patt.cor = readRDS(sprintf("df.cor_6672g_%s_170203.rds",cl.methods[m]))
df.patt.cor.p = readRDS(sprintf("df.cor.p_6672g_%s_170203.rds",cl.methods[m]))
# df.patt.cor = readRDS("df.cor_cluster1_161221.rds")
# df.patt.cor.p = readRDS("df.cor.p_cluster1_161221.rds")
# colnames(df.patt.cor) = sprintf("pat%02d",1:K)
# colnames(df.patt.cor.p) = sprintf("pat%02d",1:K)
K = ncol(df.patt.cor)

# colnames(df.patt.cor) = sprintf("pat%02d",1:K)
# colnames(df.patt.cor.p) = sprintf("pat%02d",1:K)

# save(df.mat, df.mat.2, df.mat.cc, df.mat.cc.up, df.mat.cc.dn, df.patt, df.patt.cor, df.patt.cor.p,
#      file = "patter_discovery_data_160922.RData")

# pat = 1



# STARTED HERE

dn.patt = file.path(PROJECT_DIR, "RESULTS/Microarrays/PBMC/pattern_discovery")
df.patt.cor = readRDS(file.path(dn.patt, "df.cor_6672g.rds"))
df.patt.cor.p = readRDS(file.path(dn.patt, "df.cor.p_6672g.rds"))
K = ncol(df.patt.cor)

pat.q.th = 0.3
df.cor = as.data.frame(df.patt.cor) %>%
  tibble::rownames_to_column("subj_gene") %>%
  separate(subj_gene, c("subject","gene"), sep="__", extra="merge") %>%
  gather("pattern", "cor", -(1:2)) %>%
  group_by(pattern, gene) %>%
  mutate(ord = order(cor, decreasing = T)) %>%
  ungroup()

df.cor.p = as.data.frame(df.patt.cor.p) %>%
  tibble::rownames_to_column("subj_gene") %>%
  separate(subj_gene, c("subject","gene"), sep="__", extra="merge") %>%
  gather("pattern", "cor", -(1:2)) %>%
  # group_by(pattern, gene) %>%
  mutate(ord = df.cor$ord)
  # mutate(ord = order(cor, decreasing = F)) %>%
  # ungroup() %>%

df.cor = df.cor %>%
  filter(ord==floor(max(ord)*pat.q.th)) %>%
  select(-subject) %>%
  spread(pattern, cor)
df.cor.p = df.cor.p %>%
  filter(ord==floor(max(ord)*pat.q.th)) %>%
  select(-subject) %>%
  spread(pattern, cor)

# colSums(df.cor[,3:(K+2)] > 0.8)
# colSums(df.cor.p[,3:(K+2)] < 0.01)
# colSums(df.cor.p[,3:(K+2)] < 0.01 & df.cor[,3:(K+2)] > 0)

# fn.out = sprintf("patterns_%s_correlations_q%.1f_%dg_%s.txt", cl.methods[m], pat.q.th, sum(gi), today.dt)
# write.table(df.cor, file=fn.out, sep="\t", row.names=F, col.names=T, quote = F)
# fn.out = sprintf("patterns_%s_correlations_p__q%.1f_%dg_%s.txt", cl.methods[m], pat.q.th, sum(gi), today.dt)
# write.table(df.cor.p, file=fn.out, sep="\t", row.names=F, col.names=T, quote = F)

patt.genes = list()
for (k in 1:K) {
  patt.genes[[k]] = df.cor.p$gene[df.cor.p[,k+2] < 0.01 & df.cor[,k+2] > 0]
  # patt.genes[[k]] = df.cor$gene[df.cor[,k+2] > 0.7]
}
names(patt.genes) = colnames(df.cor)[-(1:2)]
sapply(patt.genes,length)

# plot heatmap of time profiles for gene(selected)-subject(all) pairs for each pattern 

# cm = rev(redblue(16)) # blue-white-red
# # if(sum(df.cl.mat==0)>0) cm = c("#FFFFFFFF",cm)
# if(sum(is.na(df.mat.2))>0) cm = c("#000000FF",cm)
# mmx = min(abs(min(df.mat.2,na.rm=T)),abs(max(df.mat,na.rm=T)))
# col.br <- c(seq(-mmx/2,mmx/2,len=length(cm)+1)) 
# # for (k in 1:K) {
#   ttl = sprintf("pattern %d/%d : %d genes", k, K, length(patt.genes[[k]]) )
#   pgi = df$gene %in% patt.genes[[k]]
#   hm = heatmap.2(df.mat.2[pgi,], Rowv=T, Colv=NA, dendrogram="row", scale="none", col=cm, breaks=col.br,
#                  key=F, symkey=F, density.info="histogram", trace="none",
#                  labRow = F,
#                  margin=c(10,10), cexRow=1, #cexCol=1,
#                  lmat = cbind(c(0,4,2),c(0,3,1)), lwid = c(0.5,4), lhei = c(0.6,0.1,4))
#   title(ttl, cex.main=3)
#   dev.copy(png, sprintf("pattern_%d_time_heatmap_fc.%.1f_q%.1f_%s.png",k,fc.th,q.th, today.dt), width=1000,height=4000)
#   dev.off()
# }

# plot heatmap of correlations in gene vs. subject for each pattern 

# cm = rev(redblue(16)) # blue-white-red
# # if(sum(df.cl.mat==0)>0) cm = c("#AAAAAAFF",cm)
# if(sum(is.na(df.cor))>0) cm = c("#AAAAAAFF",cm)
# mmx = min(abs(min(df.cor[-(1:2)],na.rm=T)),abs(max(df.mat,na.rm=T)))
# col.br <- c(seq(-mmx/2,mmx/2,len=length(cm)+1)) 
# for (k in 1:K) {
#   ttl = sprintf("pattern %d/%d : %d genes", k, K, length(patt.genes[[k]]) )
# 
#   df.cor.tmp = as.data.frame(df.patt.cor) %>%
#     tibble::rownames_to_column("subj_gene") %>%
#     separate(subj_gene, c("subject","gene"), sep="__", extra="merge") %>%
#     rename_(cor = names(patt.genes)[k]) %>%
#     select(subject, gene, cor) %>%
#     filter(gene %in% patt.genes[[k]]) %>%
#     spread(subject, cor) %>%
#     tibble::remove_rownames() %>%
#     tibble::column_to_rownames("gene") %>%
#     as.matrix()
# 
#   hm = heatmap.2(df.cor.tmp, Rowv=T, Colv=T, dendrogram="both", scale="none", col=cm, breaks=col.br,
#                  key=F, symkey=F, density.info="histogram", trace="none",
#                  margin=c(10,10), cexRow=1, #cexCol=1,
#                  lmat = cbind(c(0,4,2),c(0,3,1)), lwid = c(0.5,4), lhei = c(0.6,0.1,4))
#   title(ttl, cex.main=3)
#   dev.copy(png, sprintf("pattern_%d.%d_cor_heatmap_fc.%.1f_q%.1f_%s.png",k,K,fc.th,q.th, today.dt), width=1000,height=2000)
#   dev.off()
# }

# # plot venn diagram
# require(VennDiagram)
# pg = patt.genes[2:5]
# plot.new()
# vd = draw.quad.venn(length(pg[[1]]), length(pg[[2]]), length(pg[[3]]), length(pg[[4]]),
#                      length(intersect(pg[[1]],pg[[2]])), length(intersect(pg[[1]],pg[[3]])),
#                     length(intersect(pg[[1]],pg[[4]])), length(intersect(pg[[2]],pg[[3]])),
#                     length(intersect(pg[[2]],pg[[4]])), length(intersect(pg[[3]],pg[[4]])),
#                     length(intersect(intersect(pg[[1]],pg[[2]]),pg[[3]])), 
#                     length(intersect(intersect(pg[[1]],pg[[2]]),pg[[4]])),
#                     length(intersect(intersect(pg[[1]],pg[[3]]),pg[[4]])), 
#                     length(intersect(intersect(pg[[2]],pg[[3]]),pg[[4]])),
#                     length(intersect(intersect(intersect(pg[[1]],pg[[2]]),pg[[3]]),pg[[4]])),
#                     category = names(pg),
#                     lty = "blank", fill = c("skyblue", "pink1", "mediumorchid", "orange"), cex=2, cat.cex=2)
# dev.copy(png, "venn_patterns_2_5.png", width=800, height=800)
# dev.off()

# plot coheatmap of correlations for selected genes and all patterns (venn like representation)

# pgx = df.cor.p[,-(1:2)]<0.01 & df.cor[,-(1:2)]>0
# pgi = apply(pgx, 1, function(x) any(x))
# 
# # pgi = apply(df.cor.p[,-(1:2)], 1, function(x) any(x < 0.01))
# sum(pgi) # 2948 for cluster3
# 
# cm = rev(redblue(16)) # blue-white-red
# mmx = min(abs(min(df.cor[-(1:2)],na.rm=T)),abs(max(df.mat,na.rm=T)))
# col.br <- c(seq(-mmx/2,mmx/2,len=length(cm)+1))
# hm = heatmap.2(+(as.matrix(df.cor.p[pgi,-(1:2)])<0.01), Rowv=T, Colv=NA, dendrogram="row", scale="none", col=grey(1:0), #breaks=col.br,
#                key=F, symkey=F, density.info="histogram", trace="none",
#                margin=c(10,10), cexRow=1, cexCol=3,
#                lmat = cbind(c(0,4,2),c(0,3,1)), lwid = c(0.5,4), lhei = c(0.6,0.1,4))
# dev.copy(png, sprintf("pattern_gene_overlap_%s.png",today.dt), width=800, height=800)
# dev.off()
# hm2 = heatmap.2(as.matrix(df.cor[pgi,-(1:2)])[rev(hm$rowInd),], Rowv=NA, Colv=NA, dendrogram="none", scale="none", col=cm, breaks=col.br,
#                key=T, symkey=F, density.info="histogram", trace="none",
#                margin=c(10,10), cexRow=1, cexCol=3)
# # cidx = as.numeric(names(df.cor)[-(1:2)])
# # cidx = match(sort(cidx),cidx)
# hm3 = heatmap.2(as.matrix(df.cor[pgi,-(1:2)]), Rowv=T, Colv=NA, dendrogram="row", scale="none", 
#                 col=cm, breaks=col.br, labRow = F,
#                 key=T, symkey=F, density.info="histogram", trace="none",
#                 margin=c(10,10), cexRow=1, cexCol=3)
# # lmat = cbind(c(0,4,2),c(0,3,1)), lwid = c(0.5,4), lhei = c(0.6,0.1,4))
# dev.copy(png, sprintf("pattern_gene_overlap_cor_%s_%s.png",cl.methods[m],today.dt), width=1000, height=1000)
# dev.off()



# trying to cluster genes by on/off pattern
# patt.n = colnames(df.cor)[-(1:2)]
# for (k in 1:K) {
#   k=3
#   ttl = sprintf("pattern %d/%d : %d genes", k, K, length(patt.genes[[k]]) )
# 
#     df.cor.tmp = as.data.frame(df.patt.cor.p) %>%
#     tibble::rownames_to_column("subj_gene") %>%
#     separate(subj_gene, c("subject","gene"), sep="__", extra="merge") %>%
#     rename_(cor = patt.n[k]) %>%
#     select(subject, gene, cor) %>%
#     filter(gene %in% patt.genes[[k]]) %>%
#     spread(subject, cor) %>%
#     tibble::remove_rownames() %>%
#     tibble::column_to_rownames("gene") %>%
#     as.matrix()
#   
#   hm = heatmap.2(+(df.cor.tmp<0.01), Rowv=T, Colv=T, dendrogram="both", scale="none", col=gray(c(1,0)), #breaks=col.br,
#                  distfun = function(x) dist(x,method="binary"),
#                  hclustfun = function(x) hclust(x, method="complete"),
#                  key=F, symkey=F, density.info="histogram", trace="none",
#                  margin=c(10,10), cexRow=1, #cexCol=1,
#                  lmat = cbind(c(0,4,2),c(0,3,1)), lwid = c(0.5,4), lhei = c(0.6,0.1,4))
#   title(ttl, cex.main=3)
#   # dev.copy(png, sprintf("pattern_%d_cor_heatmap_fc.%.1f_q%.1f.png",k,fc.th,q.th), width=1000,height=2000)
#   # dev.off()
# }

dn.fig = file.path(PROJECT_DIR, "FIGURES/GE_QC")
dir.create(dn.fig, showWarnings = F)


library(ComplexHeatmap)
library(circlize)
# cm = rev(redblue(16)) # blue-white-red
# if(sum(is.na(df.cor))>0) cm = c("#000000FF",cm)
# mmx = min(abs(min(df.cor[-(1:2)],na.rm=T)),abs(max(df.mat,na.rm=T)))
# col.br <- c(seq(-mmx/2,mmx/2,len=length(cm)+1)) 
hm = list()
patt.n = colnames(df.cor)[-(1:2)]
# patt.genes.clean = patt.genes
df.qm = data.frame()
for (k in 1:K) {
  cat(k,"")
  df.cor.tmp = as.data.frame(df.patt.cor) %>%
    tibble::rownames_to_column("subj_gene") %>%
    separate(subj_gene, c("subject","gene"), sep="__", extra="merge") %>%
    rename_(cor = patt.n[k]) %>%
    select(subject, gene, cor) %>%
    filter(gene %in% patt.genes[[k]]) %>%
    spread(subject, cor) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames("gene") %>%
    as.matrix()
  if(nrow(df.cor.tmp)<3) next
  hm[[k]] = Heatmap(df.cor.tmp, name = "correlation", km = 2, col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
          # top_annotation = ha, top_annotation_height = unit(4, "mm"), 
          show_row_names = F, show_column_names = T)
  set.seed(123)
  ro = row_order(hm[[k]])
  ksum = sapply(ro, function(x) {sum(df.cor.tmp[x,])/length(x) } )
  k.sel = which.max(ksum)
  
  kdist1 = cor(df.cor.tmp[ro[[k.sel]],] %>% t(), use="pairwise.complete.obs") %>% as.dist()
  kdist2 = cor(df.cor.tmp[ro[[3-k.sel]],] %>% t(), use="pairwise.complete.obs") %>% as.dist()
  
  # kdist = dist(df.cor.tmp[ro[[k.sel]],]) # Euclidean distance
  kdist.summary = kdist1 %>% summary()
  
  tdist = 1-cor(df.cor.tmp %>% t(), use="pairwise.complete.obs")
  cl = rep(0, nrow(df.cor.tmp))
  cl[ro[[k.sel]]] = 1
  cl[ro[[3-k.sel]]] = 2
  
  
  sil = cluster::silhouette(cl, dmatrix=tdist)
  factoextra::fviz_silhouette(sil, label = FALSE, print.summary = F) +
    scale_color_manual(values = c("black", "grey")) +
    scale_fill_manual(values = c("black", "grey")) +
    guides(col=F, fill=F) + 
    ggtitle(sprintf("G%02d",k)) +
    theme_minimal() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank())
  fn.fig = file.path(dn.fig, sprintf("GE_patterns_silhouette_G%02d", k))
  ggsave(paste0(fn.fig, ".png"), w=4,h=3)
  ggsave(paste0(fn.fig, ".pdf"), w=4,h=3)
  
  sil.asw = summary(sil)$clus.avg.width
  df.sil = data.frame(asw.1 = sil.asw[1], asw.2 = sil.asw[2])
  
  # cat(sprintf("pattern %d: selected %d, ratio %.2f", k, k.sel, max(ksum)/min(ksum)))
  # cat(ksum,k.sel)
  df.qm = rbind(df.qm,
                      cbind(
                      data.frame(pattern=k, selected=k.sel, 
                           gene.sel=length(ro[[k.sel]]), genes.orig=length(patt.genes[[k]]),
                           cor.sum.1=max(ksum), cor.sum.2=min(ksum),
                           ratio= max(ksum)/min(ksum) ),
                      cordist.mean.1 = mean(kdist1), cordist.mean.2 = mean(kdist2),
                      # kdist.summary %>% unclass() %>% as.matrix() %>% t() %>% as.data.frame(),
                      df.sil
                      ))
  # ttl = sprintf("pattern %d/%d, selected %d: %d/%d genes, ratio %.1f", 
  #               k, K, k.sel, length(ro[[k.sel]]), length(patt.genes[[k]]), max(ksum)/min(ksum) )
  # hm[[k]]@column_title = ttl
  # set.seed(123)
  # draw(hm[[k]])
  # png(sprintf("pattern_%d.%d_%s_cor_km2_heatmap_%s.png",k,K,cl.methods[m],today.dt), width=600,height=800)
  # set.seed(123) # <--- IMPORTANT FOR CONSISTENCY EVERY TIME BEFORE Heatmap_Class IS CALLED
  # draw(hm[[k]])
  # dev.off()
  
  # patt.genes.clean[[k]] = patt.genes[[k]][ro[[k.sel]]]
}
# saveRDS(patt.genes.clean, file=sprintf("patt.genes.clean_%s_%s.rds", cl.methods[m], today.dt))
df.qm = df.qm %>% 
  mutate(diff = cor.sum.1 - cor.sum.2)
# write.table(df.qm, "pattern_cleaning_table.txt", sep="\t", col.names=T, row.names=F, quote=F)

df.qm = df.qm %>% mutate(pattern=fct_rev(sprintf("Gp%02d", pattern)))
  # ggplot() + 
  #   geom_point(aes(x=pattern, y=diff), size=3)

p0 = df.qm %>% 
  mutate(ngenes2 = genes.orig-gene.sel) %>% 
  dplyr::select(pattern, ngenes1 = gene.sel, ngenes2) %>% 
  gather("stat", "value", -pattern) %>% 
  mutate(stat = sub("ngenes", "cluster", stat) %>% fct_rev()) %>% 
  ggplot(aes(pattern, value, fill=stat)) + 
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("grey", "black")) +
  xlab("") +
  ylab("No. genes") +
  coord_flip() +
  guides(fill = guide_legend(reverse=T)) +
  theme_bw() +
  theme(legend.title = element_blank(), legend.direction = "horizontal", legend.position = "bottom")


p1 = df.qm %>% 
  dplyr::select(pattern, matches("cor.sum")) %>% 
  gather("stat", "value", -pattern) %>% 
  mutate(stat = sub("cor.sum.", "cluster", stat) %>% fct_rev()) %>% 
  ggplot(aes(pattern, value, fill=stat)) + 
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("grey", "black")) +
    xlab("") +
    ylab("Correlation sum / N.genes") +
    coord_flip() +
    guides(fill = guide_legend(reverse=T)) +
    theme_bw() +
  theme(legend.title = element_blank(), legend.direction = "horizontal", legend.position = "bottom")
# ggsave("pattern_cleaning_correlations.sum.norm.png",w=4,h=3)

# df.qm[,c(1,8:13)] %>% 
#   gather("stat", "value", -pattern) %>% 
#   mutate(stat=factor(stat, levels=unique(.$stat) %>% rev())) %>% 
#   ggplot(aes(pattern, value, group=stat, col=stat)) +
#     geom_line(size=1) +
#     geom_point(size=2) +
#     geom_hline(yintercept=df.qm$Median[9], lty=2, col="grey") +
#     theme_bw()
# ggsave("pattern_cleaning_best.cluster_correlations.summary.png",w=5,h=4)

p2 = df.qm %>% 
  dplyr::select(pattern, matches("cordist.mean")) %>% 
  gather("stat", "value", -pattern) %>%
  mutate(stat = sub("cordist.mean.", "cluster", stat) %>% fct_rev()) %>% 
  # mutate(stat=factor(stat, levels=unique(.$stat) %>% rev())) %>%
  # filter(stat=="Mean") %>% 
  ggplot(aes(pattern, value, fill=stat)) +
  geom_col(col=NA, position = "dodge") +
  scale_fill_manual(values = c("grey", "black")) +
  xlab("") +
  ylab("Average corr. distance") +
  coord_flip() +
  guides(fill = guide_legend(reverse=T)) +
  theme_bw() +
  theme(legend.title = element_blank(), legend.direction = "horizontal", legend.position = "bottom")

p3 = df.qm %>% 
  dplyr::select(pattern, matches("asw")) %>% 
  gather("stat", "ASW", -pattern) %>% 
  mutate(stat = sub("asw.", "cluster", stat) %>% fct_rev()) %>% 
  ggplot(aes(pattern, ASW, group=stat, fill=stat)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("grey", "black")) +
  xlab("") +
  ylab("Average silhouette width") +
  coord_flip() +
  guides(fill = guide_legend(reverse=T)) +
  theme_bw() +
  theme(legend.title = element_blank(), legend.direction = "horizontal", legend.position = "bottom")
# ggsave("pattern_cleaning_silhouette.ASW.png",w=5,h=4)

g0 <- ggplot_gtable(ggplot_build(p0))
g1 <- ggplot_gtable(ggplot_build(p1))
g2 <- ggplot_gtable(ggplot_build(p2))
g3 <- ggplot_gtable(ggplot_build(p3))
l = list(g0, g1, g2, g3)

library(gridExtra)

mg = grid.arrange(arrangeGrob(g0), arrangeGrob(g1), arrangeGrob(g2), arrangeGrob(g3), nrow=1)#, width=c(0.5, 0.5))
fn.fig = file.path(dn.fig, "GE_patterns_QMs")
ggsave(file=paste0(fn.fig, ".png"), plot=mg, w=9, h=5)
ggsave(file=paste0(fn.fig, ".pdf"), plot=mg, w=9, h=5)

df.out = df.qm %>% 
  dplyr::rename(ng1 = gene.sel) %>% 
  add_column(ng2 = .$genes.orig-.$ng1, .after="ng1") %>% 
  select(-genes.orig, -selected, -ratio, -diff)
fwrite(df.out, file.path(PROJECT_DIR, "RESULTS/Microarrays/PBMC/pattern_discovery/GE_pattern_QM.txt"), sep="\t")
