# PURPOSE: To filter patterns after tree cutting. 

# was dist_distribution_by_method_170126.r

# Initialize the environment with PROJECT_DIR and some commonly used packages.
source("SCRIPTS/0_initialize.r")
library(cluster)

dn.patt = file.path(PROJECT_DIR, "RESULTS/Microarrays/PBMC/pattern_discovery")
df.mat.cc = readRDS(file.path(dn.patt, "df.mat.cc.rds"))

GE.patterns = fread(file.path(dn.patt, "GE_patterns.txt"))
zl = GE.patterns$pattern %>% sub("Gp","",.) %>% as.numeric()
names(zl) = GE.patterns$label

# statistics of correlations for each subject-gene
df.cc.summ = data.frame()
kmax = max(zl)
for(k in 1:kmax) {
  ii = zl==k
  imat = match(GE.patterns$label[ii], rownames(df.mat.cc))
  x = df.mat.cc[imat,imat]
  xsumm = sapply(data.frame(x, check.names = F),function(x) c(summary(x),SD=sd(x),IQR=IQR(x))) 
  df.tmp = as.data.frame(t(xsumm)) %>% 
    tibble::rownames_to_column("subj_gene") %>% 
    mutate(cluster=k) %>% 
    separate(subj_gene, c("subject","gene"), sep="__", remove = F)
  df.cc.summ = rbind(df.cc.summ, df.tmp)
}

df.cc.summ = df.cc.summ %>% 
  group_by(gene, cluster) %>% 
  mutate(n.subj=n()) %>% 
  ungroup()
# saveRDS(df.cc.summ, "df.cc.summ_170124.rds")

# keep genes with >1 subject in a cluster and combinations with at least 25% of positive correlations  
df.cc.summ.f = df.cc.summ %>% 
  mutate(n.subj.f = n.subj)
nr = nrow(df.cc.summ)
repeat{
  df.cc.summ.f = df.cc.summ.f %>% 
    filter(n.subj.f>1 & `1st Qu.`>0) %>% 
    group_by(gene, cluster) %>% 
    mutate(n.subj.f=n()) %>% 
    ungroup()
  cat(nrow(df.cc.summ.f), "")
  if(nrow(df.cc.summ.f) == nr) break
  nr = nrow(df.cc.summ.f)
}
# again remove genes with only 1 subject in a cluster
# df.cc.summ.f2 = df.cc.summ.f %>% 
#   filter(n.subj.f>1 & `1st Qu.`>0) 
  

# now filter original patterns -------------
idx = GE.patterns$label %in% df.cc.summ.f$subj_gene
GE.patterns.f = GE.patterns[idx,]

# Save filtered patterns to a TSV file.
dn.out = file.path(PROJECT_DIR, "/RESULTS/Microarrays/PBMC/pattern_discovery")
if(dir.exists(dn.out)){
        print("Output directory exists")
}else {
        print("Output directory doesn't exists. Creating it now.")
        dir.create(dn.out, recursive = T, showWarnings = T)
}
fwrite(GE.patterns.f, file.path(dn.out, "GE_patterns_filt.txt"), sep="\t")

