#  bulk fold change delta PERSISTENCE SIGNATURES 
suppressMessages(library(here))
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(Seurat))
'%ni%' = Negate('%in%')
library(scglmmr)


# set paths 
datapath = here("3_pseudobulk_de_workflow/h1_response_psig_generateddata_V2/"); dir.create(datapath)
figpath = here("3_pseudobulk_de_workflow/h1_response_psig_figures_V2/"); dir.create(figpath)

# read res from 4  
fit = readRDS(here("3_pseudobulk_de_workflow/h1_response_psig_generateddata_V2/fit_h1_object.rds"))
av = readRDS(here("3_pseudobulk_de_workflow/h1_response_psig_generateddata_V2/av_h1_object.rds"))
#fit = readRDS(here("3_pseudobulk_de_workflow/h1_response_psig_generateddata_V2/fit_h1_object.rds"))
#fit = readRDS(here("3_pseudobulk_de_workflow/h1_response_psig_generateddata_V2/fit_h1_object.rds"))

# enrichment 
d0_res = scglmmr::GetContrastResults(limma.fit.object.list = fit, coefficient.number = 1,contrast.name = "bl")
d0_rank = scglmmr::GetRankResults(limma.fit.object.list = fit,coefficient.number = 1,contrast.name = "bl")

# load gsea results for GP_ASO3_sc 1) bulk AS03 delta delta associated, 2) singlecell day 100 enriched 
gp_as03_sc = readRDS(file = here("3_pseudobulk_de_workflow/generated_data_v2/d100_gsea_gp_as03.rds")) 

# high responder vec
high.responders = c("205","207","209","212","215","234","237","245","250","256")


#################
# CD14mono 
# from gsea gp as03 sc create list of modules to test in h1 high responders 
mod = gp_as03_sc$CD14_Mono$leadingEdge
names(mod) =  paste(gp_as03_sc$CD14_Mono$pathway , "d100 monocyte LEdge", sep = " ")
mcombined = unlist(mod) %>% unique
mod[["d100 monocyte persist up combined"]] = mcombined
allcombined = lapply(gp_as03_sc, function(x){ x$leadingEdge %>% unlist(use.names = FALSE) })%>% unlist %>% unique()
mod[["d100 persist up combined"]] = allcombined

# globalaes 
global = list(theme(axis.text.x = element_text(size = 14)) ,
              theme(axis.title.y = element_text(size = 14)) )

### average expression 
av_ex = lapply(mod, function(x){ Matrix::colMeans(av$ClassicalMonocytes[x, ]) }) %>% bind_cols() %>% as.data.frame()
av_ex$subjectid = colnames(av$ClassicalMonocytes)
av_ex = av_ex %>% mutate(adjMFC = if_else(subjectid %in% high.responders, true = "high", false = "low"))
av_ex$adjMFC = factor(av_ex$adjMFC, levels = c("low", "high"))
for (i in 1:length(mod)) {
  vary = names(mod)[i]
  w = wilcox.test(av_ex[ ,vary] ~ adjMFC, data = av_ex)
  pv = round(w$p.value, digits = 3) 
  p = 
    ggplot(data = av_ex, aes(x = adjMFC, y = av_ex[ ,vary], fill = adjMFC)) +
    geom_boxplot(outlier.shape = NA, show.legend = FALSE) +
    ylab(vary)  + 
    theme_bw() + 
    global + 
    scale_fill_manual(values = c("blue", "red")) + 
    geom_jitter(shape = 16, size = 2.5, width = 0.08, aes(fill = NULL), show.legend = FALSE) +
    ggtitle(paste("average expression\n wilcoxon rank p =", pv, sep = " "))
  p
  ggsave(p,filename = paste0(figpath, vary, "AVG_EXP_CD14Mono_highvlowh1.pdf"), width = 3,  height = 4)
}


#################
# CD8 naive T 
# from gsea gp as03 sc create list of modules to test in h1 high responders 
mod = gp_as03_sc$CD8_Naive_Tcell$leadingEdge
names(mod) =  paste(gp_as03_sc$CD14_Mono$pathway , "d100 cd8naive LEdge", sep = " ")
mcombined = unlist(mod) %>% unique
mod[["d100 CD8Naive persist up combined"]] = mcombined
allcombined = lapply(gp_as03_sc, function(x){ x$leadingEdge %>% unlist(use.names = FALSE) })%>% unlist %>% unique()
mod[["d100 persist up combined"]] = allcombined

### average expression 
av_ex = lapply(mod, function(x){ Matrix::colMeans(av$`CD8+NaiveT`[x, ]) }) %>% bind_cols() %>% as.data.frame()
av_ex$subjectid = colnames(av$`CD8+NaiveT`)
av_ex = av_ex %>% mutate(adjMFC = if_else(subjectid %in% high.responders, true = "high", false = "low"))
av_ex$adjMFC = factor(av_ex$adjMFC, levels = c("low", "high"))
for (i in 1:length(mod)) {
  vary = names(mod)[i]
  w = wilcox.test(av_ex[ ,vary] ~ adjMFC, data = av_ex)
  pv = round(w$p.value, digits = 3) 
  p = 
    ggplot(data = av_ex, aes(x = adjMFC, y = av_ex[ ,vary], fill = adjMFC)) +
    geom_boxplot(outlier.shape = NA, show.legend = FALSE) +
    ylab(vary)  + 
    theme_bw() + 
    global + 
    scale_fill_manual(values = c("blue", "red")) + 
    geom_jitter(shape = 16, size = 2.5, width = 0.08, aes(fill = NULL), show.legend = FALSE) +
    ggtitle(paste("average expression\n wilcoxon rank p =", pv, sep = " "))
  p
  ggsave(p,filename = paste0(figpath, vary, "AVG_EXP_CD8_highvlowh1.pdf"), width = 3,  height = 4)
}

#####################
# mDC
# from gsea gp as03 sc create list of modules to test in h1 high responders 
mod = gp_as03_sc$mDC$leadingEdge
names(mod) =  paste(gp_as03_sc$mDC$pathway , "d100 mDC LEdge", sep = " ")
monocombined = unlist(mod) %>% unique
mod[["d100 mDC persist up combined"]] = monocombined
allcombined = lapply(gp_as03_sc, function(x){ x$leadingEdge %>% unlist(use.names = FALSE) })%>% unlist %>% unique()
mod[["d100 persist up combined"]] = allcombined
# mdc 
av_ex = lapply(mod, function(x){ Matrix::colMeans(av$mDC[x, ]) }) %>% bind_cols() %>% as.data.frame()
av_ex$subjectid = colnames(av$mDC)
av_ex = av_ex %>% mutate(adjMFC = if_else(subjectid %in% high.responders, true = "high", false = "low"))
av_ex$adjMFC = factor(av_ex$adjMFC, levels = c("low", "high"))
for (i in 1:length(mod)) {
  vary = names(mod)[i]
  w = wilcox.test(av_ex[ ,vary] ~ adjMFC, data = av_ex)
  pv = round(w$p.value, digits = 3) 
  p = 
    ggplot(data = av_ex, aes(x = adjMFC, y = av_ex[ ,vary], fill = adjMFC)) +
    geom_boxplot(outlier.shape = NA, show.legend = FALSE) +
    ylab(vary)  + 
    theme_bw() + 
    global + 
    scale_fill_manual(values = c("blue", "red")) + 
    geom_jitter(shape = 16, size = 2.5, width = 0.08, aes(fill = NULL), show.legend = FALSE) +
    ggtitle(paste("average expression\n wilcoxon rank p =", pv, sep = " "))
  p
  ggsave(p,filename = paste0(figpath, vary, "AVG_EXP_mDC_highvlowh1.pdf"), width = 3,  height = 4)
}




####################
# average mod z score 
####################


### Monocytes 
# calc average module z score in monocytes 
mod = gp_as03_sc$CD14_Mono$leadingEdge
names(mod) =  paste(gp_as03_sc$CD14_Mono$pathway , "d100 monocyte LEdge", sep = " ")
monocombined = unlist(mod) %>% unique
mod[["d100 monocyte persist up combined"]] = monocombined
allcombined = lapply(gp_as03_sc, function(x){ x$leadingEdge %>% unlist(use.names = FALSE) })%>% unlist %>% unique()
mod[["d100 persist up combined"]] = allcombined
av_z = calc_avg_module_zscore(module.list = mod, average.data.frame = av$ClassicalMonocytes)
mz = t(av_z) %>% 
  as.data.frame() %>% 
  rownames_to_column("subjectid") %>% 
  mutate(adjMFC = if_else(subjectid %in% high.responders, true = "high", false = "low"))
mz$adjMFC = factor(mz$adjMFC, levels = c("low", "high"))
### Change save name 
for (i in 1:length(mod)) {
  vary = names(mod)[i]
  w = wilcox.test(mz[ ,vary] ~ adjMFC, data = mz)
  pv = round(w$p.value, digits = 3) 
  rlang::as_character(vary)
  p = 
    ggplot(data = mz, aes(x = adjMFC, y = mz[ ,vary], fill = adjMFC)) +
    theme_bw() + 
    geom_boxplot(outlier.shape = NA, show.legend = FALSE) +
    ylab(vary)  + 
    global + 
    scale_fill_manual(values = c("blue", "red")) + 
    geom_jitter(shape = 16, size = 2.5, width = 0.08, aes(fill = NULL), show.legend = FALSE) +
    ggtitle(paste("zscore\n wilcoxon rank p =", pv, sep = " "))
  p
  ggsave(p,filename = paste0(figpath, vary, "Zscore_CD14Mono_highvlowh1.pdf"), width = 3,  height = 4)
}


#####################
# mDC
# from gsea gp as03 sc create list of modules to test in h1 high responders 
mod = gp_as03_sc$mDC$leadingEdge
names(mod) =  paste(gp_as03_sc$mDC$pathway , "d100 mDC LEdge", sep = " ")
monocombined = unlist(mod) %>% unique
mod[["d100 mDC persist up combined"]] = monocombined
allcombined = lapply(gp_as03_sc, function(x){ x$leadingEdge %>% unlist(use.names = FALSE) })%>% unlist %>% unique()
mod[["d100 persist up combined"]] = allcombined

# calc average module z score in monocytes 
av_z = calc_avg_module_zscore(module.list = mod, average.data.frame = av$mDC)
mz = t(av_z) %>% 
  as.data.frame() %>% 
  rownames_to_column("subjectid") %>% 
  mutate(adjMFC = if_else(subjectid %in% high.responders, true = "high", false = "low"))
mz$adjMFC = factor(mz$adjMFC, levels = c("low", "high"))
### Change save name 
for (i in 1:length(mod)) {
  vary = names(mod)[i]
  w = wilcox.test(mz[ ,vary] ~ adjMFC, data = mz)
  pv = round(w$p.value, digits = 3) 
  rlang::as_character(vary)
  p = 
    ggplot(data = mz, aes(x = adjMFC, y = mz[ ,vary], fill = adjMFC)) +
    theme_bw() + 
    global + 
    geom_boxplot(outlier.shape = NA, show.legend = FALSE) +
    ylab(vary)  + 
    scale_fill_manual(values = c("blue", "red")) + 
    geom_jitter(shape = 16, size = 2.5, width = 0.08, aes(fill = NULL), show.legend = FALSE) +
    ggtitle(paste("zscore\n wilcoxon rank p =", pv, sep = " "))
  p
  ggsave(p,filename = paste0(figpath, vary, "Zscore_mDC_highvlowh1.pdf"), width = 3,  height = 4)
}


#####################
# CD8 naive
mod = gp_as03_sc$CD8_Naive_Tcell$leadingEdge
names(mod) =  paste(gp_as03_sc$CD14_Mono$pathway , "d100 cd8naive LEdge", sep = " ")
mcombined = unlist(mod) %>% unique
mod[["d100 CD8Naive persist up combined"]] = mcombined
allcombined = lapply(gp_as03_sc, function(x){ x$leadingEdge %>% unlist(use.names = FALSE) })%>% unlist %>% unique()
mod[["d100 persist up combined"]] = allcombined

av_z = calc_avg_module_zscore(module.list = mod, average.data.frame = av$`CD8+NaiveT`)
mz = av_z %>%
  t %>% 
  as.data.frame() %>% 
  rownames_to_column("subjectid") %>% 
  mutate(adjMFC = if_else(subjectid %in% high.responders, true = "high", false = "low"))
mz$adjMFC = factor(mz$adjMFC, levels = c("low", "high"))
for (i in 1:length(mod)) {
  vary = names(mod)[i]
  w = wilcox.test(mz[ ,vary] ~ adjMFC, data = mz)
  pv = round(w$p.value, digits = 3) 
  rlang::as_character(vary)
  p = 
    ggplot(data = mz, aes(x = adjMFC, y = mz[ ,vary], fill = adjMFC)) +
    theme_bw() + 
    global + 
    geom_boxplot(outlier.shape = NA, show.legend = FALSE) +
    ylab(vary)  + 
    scale_fill_manual(values = c("blue", "red")) + 
    geom_jitter(shape = 16, size = 2.5, width = 0.08, aes(fill = NULL), show.legend = FALSE) +
    ggtitle(paste("zscore\n wilcoxon rank p =", pv, sep = " "))
  p
  ggsave(p,filename = paste0(figpath, vary, "Zscore_CD8Naive_highvlowh1.pdf"), width = 3,  height = 4)
}
# 
# 
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
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] scglmmr_0.1.0   magrittr_2.0.1  Seurat_2.3.4    Matrix_1.2-15   cowplot_0.9.4   here_0.1        forcats_0.4.0   stringr_1.4.0   dplyr_0.8.5    
# [10] purrr_0.3.3     readr_1.3.1     tidyr_1.0.2     tibble_2.1.1    ggplot2_3.1.1   tidyverse_1.2.1
# 
# loaded via a namespace (and not attached):
#   [1] rappdirs_0.3.1           estimability_1.3         prabclus_2.3-1           R.methodsS3_1.7.1        coda_0.19-2             
# [6] acepack_1.4.1            bit64_0.9-7              knitr_1.23               irlba_2.3.3              multcomp_1.4-10         
# [11] R.utils_2.8.0            data.table_1.12.2        rpart_4.1-13             RCurl_1.95-4.12          doParallel_1.0.14       
# [16] generics_0.0.2           metap_1.1                snow_0.4-3               BiocGenerics_0.28.0      lambda.r_1.2.3          
# [21] TH.data_1.0-10           RSQLite_2.1.1            RANN_2.6.1               europepmc_0.3            proxy_0.4-23            
# [26] bit_1.1-14               enrichplot_1.2.0         xml2_1.2.0               lubridate_1.7.4          httpuv_1.5.1            
# [31] ggsci_2.9                assertthat_0.2.1         viridis_0.5.1            xfun_0.7                 hms_0.4.2               
# [36] promises_1.0.1           DEoptimR_1.0-8           progress_1.2.2           caTools_1.17.1.2         readxl_1.3.1            
# [41] igraph_1.2.4.1           DBI_1.0.0                geneplotter_1.60.0       htmlwidgets_1.3          futile.logger_1.4.3     
# [46] stats4_3.5.3             ellipsis_0.3.0           ggpubr_0.2               backports_1.1.4          annotate_1.60.1         
# [51] gbRd_0.4-11              vctrs_0.2.4              Biobase_2.42.0           ROCR_1.0-7               withr_2.1.2             
# [56] ggforce_0.2.2            packrat_0.5.0            triebeard_0.3.0          robustbase_0.93-5        checkmate_1.9.3         
# [61] emmeans_1.3.4            prettyunits_1.0.2        mclust_5.4.5             cluster_2.0.7-1          DOSE_3.8.2              
# [66] ape_5.3                  segmented_0.5-4.0        lazyeval_0.2.2           crayon_1.3.4             hdf5r_1.2.0             
# [71] labeling_0.3             edgeR_3.24.3             pkgconfig_2.0.2          tweenr_1.0.1             nlme_3.1-137            
# [76] nnet_7.3-12              rlang_0.4.5              diptest_0.75-7           lifecycle_0.1.0          sandwich_2.5-1          
# [81] doSNOW_1.0.16            modelr_0.1.4             VennDiagram_1.6.20       cellranger_1.1.0         rprojroot_1.3-2         
# [86] polyclip_1.10-0          GSVA_1.30.0              lmtest_0.9-37            graph_1.60.0             urltools_1.7.3          
# [91] boot_1.3-20              zoo_1.8-6                base64enc_0.1-3          ggridges_0.5.1           pheatmap_1.0.12         
# [96] png_0.1-7                viridisLite_0.3.0        bitops_1.0-6             R.oo_1.22.0              KernSmooth_2.23-15      
# [101] blob_1.1.1               lars_1.2                 qvalue_2.14.1            gridGraphics_0.4-1       S4Vectors_0.20.1        
# [106] reactome.db_1.66.0       scales_1.0.0             graphite_1.28.2          memoise_1.1.0            GSEABase_1.44.0         
# [111] plyr_1.8.4               ica_1.0-2                gplots_3.0.1.1           bibtex_0.4.2             gdata_2.18.0            
# [116] compiler_3.5.3           lsei_1.2-0               RColorBrewer_1.1-2       lme4_1.1-21              fitdistrplus_1.0-14     
# [121] cli_1.1.0                dtw_1.20-1               pbapply_1.4-0            formatR_1.6              htmlTable_1.13.1        
# [126] Formula_1.2-3            MASS_7.3-51.1            tidyselect_0.2.5         stringi_1.4.3            GOSemSim_2.8.0          
# [131] locfit_1.5-9.1           latticeExtra_0.6-28      ggrepel_0.8.1            grid_3.5.3               fastmatch_1.1-0         
# [136] tools_3.5.3              parallel_3.5.3           rstudioapi_0.10          foreach_1.4.4            foreign_0.8-71          
# [141] gridExtra_2.3            farver_1.1.0             Rtsne_0.15               ggraph_1.0.2             digest_0.6.25           
# [146] rvcheck_0.1.3            shiny_1.3.2              fpc_2.2-1                Rcpp_1.0.1               broom_0.5.2             
# [151] egg_0.4.5                SDMTools_1.1-221.1       later_0.8.0              org.Hs.eg.db_3.7.0       httr_1.4.0              
# [156] AnnotationDbi_1.44.0     npsurv_0.4-0             kernlab_0.9-27           Rdpack_0.11-0            colorspace_1.4-1        
# [161] rvest_0.3.4              XML_3.98-1.19            reticulate_1.12          IRanges_2.16.0           splines_3.5.3           
# [166] shinythemes_1.1.2        ggplotify_0.0.3          flexmix_2.3-15           xtable_1.8-4             futile.options_1.0.1    
# [171] jsonlite_1.6             nloptr_1.2.1             UpSetR_1.4.0             modeltools_0.2-22        R6_2.4.0                
# [176] Hmisc_4.2-0              pillar_1.4.1             htmltools_0.3.6          mime_0.6                 glue_1.3.1              
# [181] minqa_1.2.4              clusterProfiler_3.10.1   BiocParallel_1.16.6      class_7.3-15             codetools_0.2-16        
# [186] fgsea_1.8.0              tsne_0.1-3               mvtnorm_1.0-10           lattice_0.20-38          pbkrtest_0.4-7          
# [191] mixtools_1.1.0           colorRamps_2.3           ReactomePA_1.26.0        gtools_3.8.1             GO.db_3.7.0             
# [196] survival_2.43-3          limma_3.38.3             munsell_0.5.0            DO.db_2.9                iterators_1.0.10        
# [201] variancePartition_1.12.3 haven_2.1.0              reshape2_1.4.3           gtable_0.3.0 