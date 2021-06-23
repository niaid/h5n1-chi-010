# process CHI Array data with dream lme4 pipeline 
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(variancePartition))
suppressMessages(library(here))
suppressMessages(library(scglmmr))
'%ni%' = Negate('%in%')
# source(file = "functions/analysis_functions.R")

# make output directories 
figpath = file.path(here("d70_H1/figures/")); dir.create(figpath)
datapath = file.path(here("d70_H1/generated_data/")); dir.create(datapath)
sig_genes = readRDS(file = 'signature_curation/H5_GenePatterns.rds')
names(sig_genes) = str_replace_all(string = names(sig_genes),pattern = 'GE_pattern_', replacement = "Gp")

### PART I CHI Array data 
# read array data 
array = data.table::fread("data/CHI_H1N1_data/microarray/CHI_GE_matrix_gene.txt", data.table = F) %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames("gene") %>% 
  select(., matches("day0|day70")) %>% 
  select(-matches("pre")) 

# some subjects dont have a day 70 sample 
table(str_sub(colnames(array ), 1,3)) 
# filter 
array1 = 
  array %>% 
  select(which(substr(names(.),1,3) %ni% c('237', '241', '242', '243', '274', '275')))

# Metadata 
d1md  = 
  colnames(array1) %>% 
  base::as.data.frame() %>% 
  rename(sample = ".") %>% 
  mutate(timepoint = str_sub(sample, -5,-1)) %>% 
  mutate(subjectid = str_sub(sample, 1,3)) %>%
  column_to_rownames("sample") 
d1md$timepoint = str_replace_all(d1md$timepoint, pattern = "_", replacement = "")
d1md$timepoint = factor(d1md$timepoint, levels = c('day0', 'day70'))
x = d1md$timepoint
x = model.matrix(~0 + x)
colnames(x) = str_sub(colnames(x), start = 2, end = -1)

# run model with variance paratition LME4 without voom weights (bc have normalized data) # parallelize 
library(edgeR) 
cl = makeCluster(6)
doParallel::registerDoParallel(cl = cl)

# specify formula for lme4 / dream and set up contrasts 
form = ~ 0 + timepoint + (1|subjectid)  

# specify L for old version of dream
contrast_1 = makeContrasts(timepoint = (day70 - day0), levels = colnames(x))

# format contrast matrix for variance partition.
contrast_1 = as.matrix(contrast_1) %>% as.data.frame()
rownames(contrast_1) = paste('timepoint', rownames(contrast_1), sep = "")

# run LME4 model 
result1 = dream(exprObj = array1, formula = form, L = contrast_1,  data =  d1md, useWeights = FALSE)

# save results  
saveRDS(result1, file = paste0(datapath, "day70dream_result_d1.rds"))
result1 = readRDS(here('d70_H1/generated_data/day70dream_result_d1.rds'))

# source old scglmmr functions for 3.5 
res = GetContrastResultsRaw(limma.fit.object.list = list('bulk' = result1), coefficient.number = 1,contrast.name = 'd70')
ranks = GetRankResultsRaw(contrast.result.raw.list = res)

# 
gsea1 = RunFgseaOnRankList(rank.list.celltype = ranks, pathways = sig_genes[c(4, 6:8)])
gsea_bind = RbindGseaResultList(gsea_result_list = gsea1,NES_filter = -Inf, padj_filter = 0.05)

# plot enrichment plots for up patterns 
for (i in names(sig_genes[c(4, 6:8)])) {
  p = fgsea::plotEnrichment(pathway = sig_genes[[i]],stats = ranks$bulk) +
    ggtitle(paste0(i, ' day 70 vs baseline 2009 TIV')) + 
    geom_line(color = 'red' ,size = 1); print(p)
  ggsave(p, filename = paste0(figpath, names(sig_genes[i]),'.pdf'), width = 3.3, height = 1.8)
}

# combined 
p = ggplot(gsea_bind, aes(x = NES, y= pathway )) + 
  geom_point(shape = 21, size = 4, fill = 'red') + 
  xlab("Normalized Enrichmnent Score") + 
  xlim(c(-1, 3)) + 
  theme(axis.text.x = element_text(size = 8)) + 
  geom_vline(xintercept = 0) + 
  theme_bw()
p
ggsave(p, filename = paste0(figpath,'gseaplot.pdf'), width = 3, height = 2)


lef = GetLeadingEdgeFull(gsea.list = gsea1, padj.filter = 0.05,NES.filter = -Inf)
res = CombineResults(gsealist = gsea1,contrastlist = res, gseafdr = 0.05,genefdr = 1)
res$contrast = "day70 - day 0"
res$formula = " form = ~ 0 + timepoint + (1|subjectid)  "
write_csv(x = res,path = paste0(datapath,'d70_res.csv'))
write_csv(x = lef$bulk, path =paste0(datapath,'d70res_leadingedge.csv'))
write_csv(x = gsea_bind,path = paste0(datapath,'d70gsea.csv'))


## 

# fisher exact test on coefficient sign. 
# h5 day 100 vs baseline results 
h5res = readRDS(file = here('microarray_analysis_d100/generated_data_V4/formattedresults_d100_lme4model_array.rds'))

# h1 day 70 vs baseline results
h1raw = readRDS(file = here('d70_H1/generated_data/day70dream_result_d1.rds'))
# extract results 
h1result = scglmmr::GetContrastResultsRaw(limma.fit.object.list = list('bulk' = h1raw), coefficient.number = 1, contrast.name = 'd70')
h1res = h1result$bulk

# load patterns as list 
sig_genes = readRDS(file = 'signature_curation/H5_GenePatterns.rds')
names(sig_genes) = str_replace_all(string = names(sig_genes),pattern = 'GE_pattern_', replacement = "Gp")

# all pattern genes n=1206
pattern_genes = sig_genes[c(4,6:8,10:12)] %>% unlist %>% unique

patternup = h5res %>%
  filter(gene %in% pattern_genes) %>% 
  filter(logFC > 0 & adj.P.Val < 0.05) %$% 
  gene

patterndown = h5res %>% 
  filter(gene %in% pattern_genes) %>% 
  filter(logFC < 0 & adj.P.Val < 0.05) %$% 
  gene

patternall = c(patternup, patterndown)

# define genes fold change signs (redundant but keeping code same )
h5up = patternup 
h5down = patterndown

# define gene fold change signs for H1 
h1up = h1res %>% filter(gene %in% patternall & logFC > 0) %$% gene 
h1down = h1res %>% filter(gene %in% patternall & logFC < 0) %$% gene 

# define table variables 
both_up = intersect(h5up, h1up) %>% length 
both_down  = intersect(h5down, h1down) %>% length

# plus minus / minus plus 
pm = h5up[h5up %in% h1down] %>% length
mp = h1up[h1up %in% h5down] %>% length

d = data.frame(`H5 up` = c(both_up, pm), `H5 down` = c(mp, both_down))
rownames(d) = c('H1.up', 'H1.down')

## Data set up for fishers exact test 
print(d)
ft  = fisher.test(d)
ft
cu = pals::brewer.reds(10)
pheatmap::pheatmap(d, color = cu,
                   cluster_cols = F, cluster_rows = F,
                   number_color = 'black',
                   display_numbers = TRUE, 
                   main = paste0('Fisher exact p = ', '2.27e-08'), fontsize = 8,
                   number_format =  "%.0f", filename = paste0(figpath, 'heatmap.pdf'), width = 2.8, height = 2)

# R version 3.5.3 Patched (2019-03-11 r77192)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS Mojave 10.14.6
# 
# Matrix products: default
# BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] splines   parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] fgsea_1.8.0              Rcpp_1.0.1               pbkrtest_0.4-7           lmerTest_3.1-2          
# [5] lme4_1.1-21              edgeR_3.24.3             Seurat_2.3.4             Matrix_1.2-15           
# [9] cowplot_0.9.4            viridis_0.5.1            viridisLite_0.3.0        here_0.1                
# [13] variancePartition_1.12.3 Biobase_2.42.0           BiocGenerics_0.28.0      scales_1.0.0            
# [17] foreach_1.4.4            limma_3.38.3             magrittr_2.0.1           forcats_0.4.0           
# [21] stringr_1.4.0            dplyr_0.8.5              purrr_0.3.3              readr_1.3.1             
# [25] tidyr_1.0.2              tibble_2.1.1             ggplot2_3.1.1            tidyverse_1.2.1         
# 
# loaded via a namespace (and not attached):
#   [1] readxl_1.3.1        snow_0.4-3          backports_1.1.4     fastmatch_1.1-0     Hmisc_4.2-0        
# [6] plyr_1.8.4          igraph_1.2.4.1      lazyeval_0.2.2      BiocParallel_1.16.6 digest_0.6.25      
# [11] htmltools_0.3.6     lars_1.2            fansi_0.4.0         gdata_2.18.0        checkmate_1.9.3    
# [16] cluster_2.0.7-1     doParallel_1.0.14   mixtools_1.1.0      ROCR_1.0-7          modelr_0.1.4       
# [21] R.utils_2.8.0       prettyunits_1.0.2   colorspace_1.4-1    rvest_0.3.4         haven_2.1.0        
# [26] xfun_0.7            crayon_1.3.4        jsonlite_1.6        zoo_1.8-6           survival_2.43-3    
# [31] iterators_1.0.10    ape_5.3             glue_1.3.1          gtable_0.3.0        kernlab_0.9-27     
# [36] prabclus_2.3-1      DEoptimR_1.0-8      mvtnorm_1.0-10      bibtex_0.4.2        metap_1.1          
# [41] dtw_1.20-1          progress_1.2.2      htmlTable_1.13.1    reticulate_1.12     foreign_0.8-71     
# [46] bit_1.1-14          proxy_0.4-23        mclust_5.4.5        SDMTools_1.1-221.1  Formula_1.2-3      
# [51] tsne_0.1-3          stats4_3.5.3        htmlwidgets_1.3     httr_1.4.0          gplots_3.0.1.1     
# [56] RColorBrewer_1.1-2  fpc_2.2-1           acepack_1.4.1       modeltools_0.2-22   ica_1.0-2          
# [61] pkgconfig_2.0.2     R.methodsS3_1.7.1   flexmix_2.3-15      nnet_7.3-12         utf8_1.1.4         
# [66] locfit_1.5-9.1      labeling_0.3        tidyselect_0.2.5    rlang_0.4.5         reshape2_1.4.3     
# [71] munsell_0.5.0       cellranger_1.1.0    tools_3.5.3         cli_1.1.0           generics_0.0.2     
# [76] broom_0.5.2         ggridges_0.5.1      npsurv_0.4-0        knitr_1.23          bit64_0.9-7        
# [81] fitdistrplus_1.0-14 robustbase_0.93-5   caTools_1.17.1.2    RANN_2.6.1          packrat_0.5.0      
# [86] pbapply_1.4-0       nlme_3.1-137        R.oo_1.22.0         xml2_1.2.0          hdf5r_1.2.0        
# [91] compiler_3.5.3      rstudioapi_0.10     png_0.1-7           lsei_1.2-0          stringi_1.4.3      
# [96] lattice_0.20-38     nloptr_1.2.1        vctrs_0.2.4         pillar_1.4.1        lifecycle_0.1.0    
# [101] lmtest_0.9-37       Rdpack_0.11-0       data.table_1.12.2   bitops_1.0-6        irlba_2.3.3        
# [106] gbRd_0.4-11         colorRamps_2.3      R6_2.4.0            latticeExtra_0.6-28 KernSmooth_2.23-15 
# [111] gridExtra_2.3       codetools_0.2-16    boot_1.3-20         MASS_7.3-51.1       gtools_3.8.1       
# [116] assertthat_0.2.1    rprojroot_1.3-2     withr_2.1.2         diptest_0.75-7      doSNOW_1.0.16      
# [121] hms_0.4.2           grid_3.5.3          rpart_4.1-13        class_7.3-15        minqa_1.2.4        
# [126] segmented_0.5-4.0   Rtsne_0.15          numDeriv_2016.8-1.1 lubridate_1.7.4     base64enc_0.1-3    
