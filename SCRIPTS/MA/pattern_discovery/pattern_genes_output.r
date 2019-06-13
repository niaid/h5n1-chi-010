# was in patterns_genes_output_2_170320.r

dn.patt = file.path(PROJECT_DIR, "RESULTS/Microarrays/PBMC/pattern_discovery")
df = readRDS(file.path(dn.patt, "df_6672g.rds"))

df.patt.mat = readRDS(file.path(dn.patt, "df.patt.rds"))
patt.genes.clean = readRDS(file.path(dn.patt, "patt.genes.clean.rds"))

K = length(patt.genes.clean)

GE.patterns = fread(file.path(dn.patt, "GE_patterns_filt.txt"))
cln = GE.patterns$pattern %>% sub("Gp","",.) %>% as.numeric()
names(cln) = GE.patterns$label
ct.lev = sort(unique(GE.patterns$pattern))
K=max(cln)


# lists of subjects for each cluster

subj = vector("list", K)
for (k in 1:K) {
  tmp = names(cln)[cln==k]
  tmp = data.frame(label=tmp) %>% 
    separate(label, c("subject","gene"), sep="__", extra="merge")
  subj[[k]] = unique(tmp$subject)
}

df.patt.cor = readRDS(file.path(dn.patt, "df.cor_6672g.rds"))
df.patt.cor.p = readRDS(file.path(dn.patt, "df.cor.p_6672g.rds"))

# heatmap of signature genes with pc1 contribution
genes = sort(unique(unlist(patt.genes.clean)))
gene.sig = matrix(0,nrow=length(genes), ncol=K)
rownames(gene.sig) = genes
colnames(gene.sig) = ct.lev

gene.sig.list = vector("list",K)
names(gene.sig.list) = ct.lev
subj.pc1 = matrix(NA,nrow=nrow(subj.patt.cor), ncol=K)
for(k in 1:K) {
  cat(k,"")
  tmp = df.patt.cor[,k] %>% as.data.frame() %>% 
    tibble::rownames_to_column("node") %>% 
    separate(node, c("subject","gene"), sep="__", extra="merge") %>% 
    filter(gene %in% patt.genes.clean[[k]]) %>% 
    spread_("gene", ".") %>% 
    tibble::remove_rownames() %>% 
    tibble::column_to_rownames("subject") %>% 
    as.matrix()
  pca=prcomp(tmp, scale=F, center=F)
  # pc1.ev[k] = pca$sdev[1]/sum(pca$sdev)*100
  # i1 = genes %in% patt.genes.clean[[k]]
  ii = match(rownames(pca$rotation), genes)
  if(pca$rotation[,1][which.max(abs(pca$rotation[,1]))]<0) {
    gene.sig[ii,k] = -pca$rotation[,1]
    gene.sig.list[[k]] = as.data.frame(-pca$rotation[,1,drop=F]) %>% 
      tibble::rownames_to_column("gene") %>% arrange(-PC1)
  } else {
    gene.sig[ii,k] = pca$rotation[,1]
    gene.sig.list[[k]] = as.data.frame(pca$rotation[,1,drop=F]) %>% 
      tibble::rownames_to_column("gene") %>% arrange(-PC1)
  }
  subj.pc1[,k] = pca$x[,1]
}

library(mygene)
dn.out = file.path(dn.patt, "GE_pattern_genes")
dir.create(dn.out, showWarnings = F)
for(k in 1:K) {
  gt = gene.sig.list[[k]]
  gq = queryMany(gt$gene, scopes="symbol", species="human", fields=c("symbol","entrezgene","name"))
  gq2 = as.data.frame(gq) %>% dplyr::select(gene=symbol, entrezgene, name) %>% 
    arrange(gene, entrezgene) %>% 
    distinct(gene, .keep_all = T) %>% right_join(gt, by="gene")
  fwrite(gq2, file.path(dn.out, sprintf("GE_pattern_%02d_gene.sig.contr_ann.txt",k)), sep="\t")
}
