# PURPOSE: To score each subject against each of the 14 patterns.

# was in patterns_to_subjects_170202.r

# Initialize.
source("SCRIPTS/0_initialize.r")
source(file.path(PROJECT_DIR, "SCRIPTS/functions/get_score.r"))

dn.patt = file.path(PROJECT_DIR, "RESULTS/Microarrays/PBMC/pattern_discovery")

GE.patterns = fread(file.path(dn.patt, "GE_patterns_filt.txt"))
cln = GE.patterns$pattern %>% sub("Gp","",.) %>% as.numeric()
names(cln) = GE.patterns$label
ct.lev = sort(unique(GE.patterns$pattern))
# K=max(cln)

df = readRDS(file.path(dn.patt, "df_6672g.rds"))
df.mat.2 = readRDS(file.path(dn.patt, "df.mat.2_6672g.rds"))
df.patt.mat = readRDS(file.path(dn.patt, "df.patt.rds"))
patt.genes.clean = readRDS(file.path(dn.patt, "patt.genes.clean.rds"))
# patt.genes.clean = readRDS("../patt.genes.clean_200.2000.4.cleaned_170210.rds")

K = length(patt.genes.clean)
long.idx = c(rep(1,3),2:ncol(df.mat.2))

subj.patt.cor = matrix(NA,nrow=length(unique(df$subject.id)), ncol=K)
rownames(subj.patt.cor) = unique(df$subject.id)
colnames(subj.patt.cor) = names(patt.genes.clean)

for (k in 1:K) {
  cat(k,"")
  # k = 1
  # df.score.k = df %>%
  #   filter(gene %in% patt.genes.clean[[k]]) %>%
  #   select(-gene)
  subj = unique(df$subject.id)
  for (s in seq_along(subj)) {
    X = df %>%
      filter(subject.id==subj[s], gene %in% patt.genes.clean[[k]]) %>%
      dplyr::select(-subject.id, -gene) %>%
      as.matrix() %>%
      get_score()
    
    subj.patt.cor[s,k] = cor(X[long.idx], df.patt.mat[k,long.idx], use = "pairwise.complete.obs")
  }
}

saveRDS(subj.patt.cor, file.path(dn.patt, "subj.patt.cor.rds"))

subj.patt.cor %>% as.data.frame() %>% 
  tibble::rownames_to_column("subject") %>% 
  fwrite(file.path(dn.patt, "subj.patt.cor.txt"), sep="\t")
