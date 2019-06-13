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
  # ungroup() 

df.cor = df.cor %>%
  filter(ord==floor(max(ord)*pat.q.th)) %>%
  select(-subject) %>%
  spread(pattern, cor)
df.cor.p = df.cor.p %>%
  filter(ord==floor(max(ord)*pat.q.th)) %>%
  select(-subject) %>%
  spread(pattern, cor)

n_genes = nrow(df.cor)

fn.out = sprintf("patterns_correlations__q%.2f_%dg.txt", pat.q.th, n_genes)
fwrite(df.cor, file=file.path(dn.patt, fn.out), sep="\t")
fn.out = sprintf("patterns_correlations_p__q%.2f_%dg.txt", pat.q.th, n_genes)
fwrite(df.cor.p, file=file.path(dn.patt, fn.out), sep="\t")

patt.genes = list()
for (k in 1:K) {
  patt.genes[[k]] = df.cor.p$gene[df.cor.p[,k+2] < 0.01 & df.cor[,k+2] > 0]
  # patt.genes[[k]] = df.cor$gene[df.cor[,k+2] > 0.7]
}
names(patt.genes) = colnames(df.cor)[-(1:2)]
sapply(patt.genes,length)


library(ComplexHeatmap)
library(circlize)
# cm = rev(redblue(16)) # blue-white-red
# if(sum(is.na(df.cor))>0) cm = c("#000000FF",cm)
# mmx = min(abs(min(df.cor[-(1:2)],na.rm=T)),abs(max(df.mat,na.rm=T)))
# col.br <- c(seq(-mmx/2,mmx/2,len=length(cm)+1)) 
dn.hm = file.path(PROJECT_DIR, "FIGURES/pattern_correlations")
dir.create(dn.hm, showWarnings = F)
hm = list()
patt.n = colnames(df.cor)[-(1:2)]
patt.genes.clean = patt.genes
for (k in 1:K) {
  cat(k,"")
  df.cor.tmp = as.data.frame(df.patt.cor) %>%
    tibble::rownames_to_column("subj_gene") %>%
    separate(subj_gene, c("subject","gene"), sep="__", extra="merge") %>%
    rename_(cor = patt.n[k]) %>%
    dplyr::select(subject, gene, cor) %>%
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
  # cat(sprintf("pattern %d: selected %d, ratio %.2f", k, k.sel, max(ksum)/min(ksum)))
  # cat(ksum,k.sel)
  
  ttl = sprintf("pattern %d/%d, selected %d: %d/%d genes, ratio %.1f", 
                k, K, k.sel, length(ro[[k.sel]]), length(patt.genes[[k]]), max(ksum)/min(ksum) )
  hm[[k]]@column_title = ttl
  set.seed(123)
  draw(hm[[k]])
  png(file.path(dn.hm, sprintf("pattern_%d.%d_cor_km2_heatmap.png",k,K)), width=600,height=800)
  set.seed(123) # <--- IMPORTANT FOR CONSISTENCY EVERY TIME BEFORE Heatmap_Class IS CALLED
  draw(hm[[k]])
  dev.off()
  
  patt.genes.clean[[k]] = patt.genes[[k]][ro[[k.sel]]]
}

saveRDS(patt.genes.clean, file=file.path(dn.patt, "patt.genes.clean.rds"))

