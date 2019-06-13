library(cluster)
library(ggdendro)

fn.dia = file.path("RESULTS/Microarrays/PBMC/pattern_discovery/diana.object.abs.rds")
dia = readRDS(fn.dia)

diad = dendro_data(dia)

dn.out = file.path(PROJECT_DIR, "RESULTS/Microarrays/PBMC/pattern_discovery")
dir.create(dn.out, showWarnings = F)

saveRDS(diad, file.path(dn.out, "diana_data.rds")) # as on 161208

diad = readRDS(file.path(dn.out, "diana_data.rds"))


z = segment(diad) %>% 
  mutate(ylen = y-yend) %>% 
  mutate(yd = +c(0,diff(y))>0)

y.th = 0.1 # <-------------- minimum height to cut
z$cluster = NA
k = 1
ydx = 0
for (i in 1:nrow(z)) {
  if(i<ydx) next
  if(z$ylen[i] < y.th) next
  ydx = which(z$y>=z$y[i])
  ydx = ydx[ydx>i]
  if(length(ydx)==0) {
    z$cluster[i:nrow(z)] = k
    break
  } else {
    ydx = ydx[1]
    z$cluster[i:(ydx-1)] = k
  }
  k = k + 1
}

zl = diad$labels
xm = match(zl$x, z$x)
zl$cluster = z$cluster[xm]

imat = match(zl$label, rownames(df.mat.cc))
sil = silhouette(zl$cluster, dmatrix=1-df.mat.cc[imat,imat])

df.sil = as.data.frame(sil[,]) %>% 
  mutate(x=seq_along(cluster), cluster=factor(cluster), neighbor=factor(neighbor))

df.sil.summary = as.data.frame(sil[,]) %>% 
  group_by(cluster) %>% 
  summarise(n = n(), neg.ratio = sum(sil_width<0)/n()) %>% 
  ungroup() %>% 
  mutate(neg.ratio.rank = rank(-neg.ratio))

# remove nodes with negative SW ------------------------------
df.sil.2 = df.sil %>% 
  filter(sil_width >= 0) %>% 
  group_by(cluster) %>% 
  mutate(cl.size=n()) %>% 
  ungroup() %>% 
  filter(cl.size > 100) %>% 
  mutate(cluster=factor(cluster))
# mutate(cluster=as.numeric(factor(cluster)))
# table(df.sil.2$cluster)

zl$cluster1 = zl$cluster
zl$cluster1[!zl$x %in% df.sil.2$x] = 0

i0 = zl$cluster1==0
imat = match(zl$label[!i0], rownames(df.mat.cc))
sil = silhouette(zl$cluster1[!i0], dmatrix=1-df.mat.cc[imat,imat])
# factoextra::fviz_silhouette(sil, label = FALSE, print.summary = F) +
#   theme(legend.position="none")


# cluster stability tree cut =================================
# mth.seq = seq(0.005,1.98, by=0.005)
mth.seq = seq(1,1.98, by=0.005)
ct = matrix(nrow=length(mth.seq), ncol=nrow(df.mat))
rownames(ct) = mth.seq
colnames(ct) = colnames(df.mat.cc)
ct.sz = ct
ct.nsmax = ct
ct.nsmax.f = ct
i = 1
pb <- txtProgressBar(max=nrow(ct.sz), style=3)
for (K in mth.seq) {
  ct[i,] = cutree(as.hclust(dia), h=K)
  ct.sz[i,] = table(ct[i,])[ct[i,]]
  tmp = data.frame(ct=ct[i,], node=names(ct[i,])) %>%
    separate("node",c("subject","gene"),"__") %>%
    group_by(ct,gene) %>%
    summarise(n.subj = n()) %>%
    ungroup() %>%
    group_by(ct) %>%
    summarise(n.subj.max = max(n.subj)) %>%
    ungroup()
  ct.nsmax[i,] = tmp$n.subj.max[ct[i,]]
  i = i + 1
  setTxtProgressBar(pb,i)
}

min.sz.seq = c(100,200)
min.nsmax.seq = c(4,6,8,10)

min.sz = 200
ct.sz[ct.sz < min.sz] = 0
min.nsmax = 4
ct.sz[ct.nsmax < min.nsmax] = 0
max.sz = 2000
ct.sz[ct.sz > max.sz] = 0

ct.sz.d = diff(ct.sz, na.rm=T)

ct.sz.len = ct.sz.d
for (col in 1:ncol(ct.sz.d)) {
  ztmp = ct.sz.d[,col]
  if(all(ztmp==0)) next # if whole column is 0s
  zt = which(ztmp!=0)
  if(ct.sz[zt[1],col]>0) zt = c(1,zt) # if edge cluster is not 0
  ztmp[zt] = c(0,diff(zt))
  ct.sz.len[,col] = ztmp
}

ct.sz.len.max = apply(ct.sz.len, 2, max)
ct.sz.len.max.ind = apply(ct.sz.len, 2, which.max)
ct.sz.len.max.ind.1 = max.col(t(ct.sz.len), ties.method = "first") # same as ind, just for test
ct.sz.len.max.ind.2 = max.col(t(ct.sz.len), ties.method = "last")
ct.sz.len.max.n = apply(ct.sz.len, 2, function(x) sum(x==max(x) & x > 0))

source("SCRIPTS/functions/repeated_elements_2vec.r")
ct.sz.len.max.df = repeated_elements_2vec(ct.sz.len.max[dia$order], ct.sz.len.max.ind.2[dia$order]) %>% 
  filter(value>0) %>%
  mutate(stability = value * 0.005, istart=iend-length+1, include=F)

# flag partial clusters
for (k in 1:nrow(ct.sz.len.max.df)) {
  isz = ct.sz.len.max.df$iend[k]
  if(any(ct.sz[,dia$order][,isz]==ct.sz.len.max.df$length[k])) ct.sz.len.max.df$include[k]=T
}


# cut by cl.boundary

cl.boundary = zl %>% group_by(cluster) %>% summarise(xstart=min(x),xend=max(x))

zz = ct.sz.len.max.df
ilast = nrow(zz)
for (i in 1:(nrow(cl.boundary)-1)) {
  for (j in 1:nrow(zz)) {
    if(zz$istart[j] < cl.boundary$xend[i] & zz$iend[j] > cl.boundary$xend[i]) {
      zz[ilast+1,] = zz[j,]
      zz$istart[ilast+1] = cl.boundary$xend[i]+1
      zz$iend[j] = cl.boundary$xend[i]
      ilast = ilast+1
      break
    }
  }
}
zz = zz %>%
  arrange(istart) %>%
  mutate(length=iend-istart+1) # recompute cluster length

i.stab = zz$length >= min.sz
# zz[i.stab,]
ct.sz.len.max.df = zz[i.stab,]

zz = ct.sz.len.max.df %>% filter(include)
ct.stab = rep(0, nrow(df.mat.cc))
for (k in 1:nrow(zz)) {
  ct.stab[zz$istart[k]:zz$iend[k]] = k
}
# table(ct.stab)
names(ct.stab) = names(ct.sz.len.max)[dia$order]

zl$cluster2 = ct.stab
i0 = zl$cluster2==0
imat = match(zl$label[!i0], rownames(df.mat.cc))
sil = silhouette(zl$cluster2[!i0], dmatrix=1-df.mat.cc[imat,imat])
# factoextra::fviz_silhouette(sil, label = FALSE, print.summary = F) +
#   theme(legend.position="right")

df.sil = as.data.frame(sil[,]) %>% 
  mutate(x=seq_along(cluster), cluster=factor(cluster), neighbor=factor(neighbor))
df.sil.summary = as.data.frame(sil[,]) %>% 
  group_by(cluster) %>% 
  summarise(n = n(), neg.ratio = sum(sil_width<0)/n()) %>% 
  ungroup() %>% 
  mutate(neg.ratio.rank = rank(-neg.ratio))

i0 = zl$cluster2 %in% c(0, df.sil.summary$cluster[df.sil.summary$neg.ratio > 0.5])
sum(i0)
zl$cluster3 = zl$cluster2
zl$cluster3[i0] = 0
# table(zl$cluster3)

dn.out = file.path(PROJECT_DIR, "RESULTS/Microarrays/PBMC/pattern_discovery")
dir.create(dn.out, showWarnings = F)
saveRDS(z, file.path(dn.out, "z.rds"))
saveRDS(zl, file.path(dn.out, "zl.rds"))

GE.patterns = zl %>% select(label, pattern=cluster3) %>% 
  filter(pattern != 0) %>%
  mutate(pattern = factor(pattern) %>% as.numeric()) %>% 
  mutate(pattern = sprintf("Gp%02d", pattern))

fwrite(GE.patterns, file.path(dn.out, "GE_patterns.txt"), sep="\t")
# fwrite(GE.patterns %>% filter(pattern!="Gp00"), file.path(dn.out, "GE_patterns_only.txt"), sep="\t")
