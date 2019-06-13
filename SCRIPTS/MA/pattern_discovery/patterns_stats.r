# was plot_patterns_profiles_171123.r

dn.patt = file.path(PROJECT_DIR, "RESULTS/Microarrays/PBMC/pattern_discovery")
df.mat = readRDS(file.path(dn.patt, "df.mat.rds"))

GE.patterns = fread(file.path(dn.patt, "GE_patterns_filt.txt"))
ct = GE.patterns$pattern %>% sub("Gp","",.) %>% as.numeric()
names(ct) = GE.patterns$label
ct.lev = sort(unique(GE.patterns$pattern))
K=max(ct)

# df.mat.cc = readRDS("df.mat.cc.rds")
# zl.all = readRDS("zl.all_170123.rds")
# zl.all = readRDS("zl.all.f_170126.rds")

# cl.methods = names(zl.all)[c(9)]
# df.patt = list()
# ii = 1

# for (m in seq_along(cl.methods)) {
# m=1  
# zcluster = zl.all[,cl.methods[m]]
# i0 = zcluster==0
# ct = zcluster[!i0]
# ct=as.numeric(factor(ct))
# names(ct) = zl.all$label[!i0]
# head(ct)

df.cl = as.data.frame(ct) %>%
  tibble::rownames_to_column("subj_gene") %>%
  separate(subj_gene, c("subject","gene"), sep="__") %>%
  spread(subject, ct, fill = 0)
# View(df.cl)

df.cl.mat = as.matrix(df.cl[,-1])
rownames(df.cl.mat) = df.cl$gene

idx = match(names(ct), rownames(df.mat))

df.cl.mean = as.data.frame(df.mat[idx,]) %>%
  mutate(ct=ct) %>%
  group_by(ct) %>%
  summarise_all(mean, na.rm=T) %>%
  ungroup() %>%
  gather("time","mean",-1)

df.cl.sd = as.data.frame(df.mat[idx,]) %>%
  mutate(ct=ct) %>%
  group_by(ct) %>%
  summarise_all(sd, na.rm=T) %>%
  ungroup() %>%
  gather("time","sd",-1)

df.cl.median = as.data.frame(df.mat[idx,]) %>%
  mutate(ct=ct) %>%
  group_by(ct) %>%
  summarise_all(median, na.rm=T) %>%
  ungroup() %>%
  gather("time","mean",-1)

z_score = function(x) t(scale(t(x)))

# df.cl.zavg = as.data.frame(z_score(df.mat[idx,])) %>%
#   mutate(ct=ct) %>%
#   group_by(ct) %>%
#   summarise_all(mean, na.rm=T) %>%
#   ungroup() %>%
#   gather("time","mean",-1)

x = z_score(df.mat[idx,])
df.cl.zavg0 = as.data.frame(x) %>%
  mutate(ct=ct) %>%
  group_by(ct) %>%
  summarise_all(mean, na.rm=T) %>%
  ungroup() %>%
  gather("time","mean",-1) %>% 
  group_by(ct) %>%
  mutate(mean = mean - mean[time=="d000_00h"]) %>%
  ungroup()

x = df.mat[idx,]
i.na = apply(is.na(x), 1, any)
x = x[!i.na,]
pc.nonscaled = matrix(NA, nrow=max(ct), ncol=ncol(df.mat))
rownames(pc.nonscaled) = 1:max(ct)
colnames(pc.nonscaled) = colnames(df.mat)
pc.scaled = pc.nonscaled
# pc.ns.sdev = rep(NA, max(ct))
# pc.ns.sdev.d = rep(NA, max(ct))
# pc.s.sdev = rep(NA, max(ct))
# pc.s.sdev.d = rep(NA, max(ct))
for (k in 1:max(ct)) {
  xt = x[ct[!i.na]==k,]
  pc = prcomp(xt, scale=F, center=F)
  pc.nonscaled[k,] = pc$rotation[,1]
  # pc = prcomp(t(xt), scale=F)
  # pc.nonscaled[k,] = pc$x[,1]
  # pc.ns.sdev[k] = pc$sdev[1]/sum(pc$sdev)*100
  # pc.ns.sdev.d[k] = (pc$sdev[1]-pc$sdev[2])/sum(pc$sdev)*100
  
  xt = t(scale(t(xt)))
  pc = prcomp(xt, scale=F, center=F)
  pc.scaled[k,] = pc$rotation[,1]
  # pc = prcomp(t(xt), scale=T)
  # pc.scaled[k,] = pc$x[,1]
  # pc.s.sdev[k] = pc$sdev[1]/sum(pc$sdev)*100
  # pc.s.sdev.d[k] = (pc$sdev[1]-pc$sdev[2])/sum(pc$sdev)*100
}

df.cl.pc.nonscaled = as.data.frame(pc.nonscaled) %>%
  mutate(ct=1:max(ct)) %>%
  gather("time","pc",-ct) %>% 
  group_by(ct) %>%
  mutate(pc = pc - pc[time=="d000_00h"]) %>%
  # mutate(pc = pc/10) %>% 
  ungroup()

df.cl.pc.scaled = as.data.frame(pc.scaled) %>%
  mutate(ct=1:max(ct)) %>%
  gather("time","pc",-ct) %>% 
  group_by(ct) %>%
  mutate(pc = pc - pc[time=="d000_00h"]) %>%
  # mutate(pc = pc/10) %>% 
  ungroup()

# df.cl.stat = cbind(df.cl.mean, 
#                    df.cl.sd[3]) %>%
#   mutate(ct=factor(ct))

df.cl.stat = cbind(df.cl.mean, 
                   df.cl.median[3] %>% rename(median=mean),
                   # df.cl.z0avg[3] %>% rename(z0avg=mean),
                   df.cl.zavg0[3] %>% rename(zavg0=mean),
                   df.cl.pc.nonscaled[3] %>% rename(pc.nonscaled=pc),
                   df.cl.pc.scaled[3] %>% rename(pc.scaled=pc)
) #%>%
# mutate(ct=factor(ct))

# correct the sign for PCA metrics
df.cl.stat = df.cl.stat %>% 
  group_by(ct) %>% 
  mutate(pc.nonscaled = sign(cor(mean,pc.nonscaled)) * pc.nonscaled) %>%
  mutate(pc.scaled = sign(cor(mean,pc.scaled)) * pc.scaled) %>%
  ungroup()

# x = df.cl.stat %>% select(-ct,-time) %>% as.matrix()
# pairs(x, col=df.cl.stat$ct)
# dev.copy(png, sprintf("pattern_profile_metrics_pair_corr_%s_%s.png", cl.methods[m], today.dt))
# dev.off()
# df.patt = df.cl.mean %>%
# spread(time, mean)
# df.patt.mat = as.matrix(df.patt[,-1])
# rownames(df.patt.mat) = paste0("pat",1:K)


# df.cl.stat = df.cl.stat %>%
#   mutate(ct=factor(as.numeric(factor(ct))))
df.cl.stat = df.cl.stat %>%
  gather("method","value",-(1:2))

df.cl.stat = df.cl.stat %>% 
  mutate(ct = sprintf("Gp%02d",ct)) %>% 
  mutate(ct = factor(ct, levels=unique(ct)))

saveRDS(df.cl.stat, file.path(dn.patt, "df.cl.stat.rds"))

df.patt = df.cl.mean %>% 
  spread(time,mean) %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames("ct") %>% 
  as.matrix()
saveRDS(df.patt, file.path(dn.patt, "df.patt.rds"))


