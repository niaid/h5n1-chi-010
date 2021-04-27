# save day 0 and day 100 microarray data 
suppressMessages(library(here))
suppressMessages(library(tidyverse))

# save path 
datapath = here("microarray_analysis_d100/generated_data_V4/"); dir.create(datapath)

# load eset object 
load(here("data/microarray/data/eset.genes.filtered.RData"))

# add adjuvant status  metadata 
met = read_delim(file = here("data/metadata/clinical_info_adj.txt"), delim = "\t")

# subset to day 100 
times_use = c("d000_00h","d100")
eset.genes = eset.genes[ ,eset.genes$time.point %in% times_use]

# vector of all gene pattern genes 
gp = readRDS(file = here("signature_curation/H5_GenePatterns.rds"))
patterngenes = unlist(x = gp, use.names = FALSE) %>% unique()

# signature genes in data - make sure all genes used in patterns are in object loaded 
pattern_expressed = patterngenes[patterngenes %in% rownames(exprs(eset.genes))]
stopifnot( length(patterngenes) == length(pattern_expressed) ) 

# create expression dataframe
edf = base::as.data.frame(t(exprs(eset.genes[pattern_expressed, ])))
edf = cbind(edf, eset.genes@phenoData@data)

# confirm gene names are unchanged and all in data
stopifnot( colnames(edf)[colnames(edf) %in% pattern_expressed] == pattern_expressed )

# save 
saveRDS(edf, file = paste0(datapath, 'array_d0_d100_patterngenes_andmetadata_dataframe.rds'))

#########################
sessionInfo()
# R version 3.6.1 (2019-07-05)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS Mojave 10.14.6
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] scglmmr_0.1.0            magrittr_2.0.1           variancePartition_1.16.1 Biobase_2.46.0           BiocGenerics_0.32.0      scales_1.1.0            
# [7] foreach_1.4.8            limma_3.42.2             forcats_0.5.0            stringr_1.4.0            dplyr_1.0.3              purrr_0.3.3             
# [13] readr_1.3.1              tidyr_1.0.2              tibble_3.1.0             ggplot2_3.3.0            tidyverse_1.3.0          here_0.1                
# 
# loaded via a namespace (and not attached):
#   [1] utf8_1.1.4             reticulate_1.14        tidyselect_1.1.0       lme4_1.1-21            RSQLite_2.2.0          AnnotationDbi_1.48.0  
# [7] htmlwidgets_1.5.1      grid_3.6.1             BiocParallel_1.20.1    Rtsne_0.15             munsell_0.5.0          codetools_0.2-16      
# [13] ica_1.0-2              future_1.16.0          withr_2.4.0            colorspace_1.4-1       GOSemSim_2.12.1        rstudioapi_0.11       
# [19] Seurat_3.1.5           stats4_3.6.1           ROCR_1.0-7             ggsignif_0.6.0         DOSE_3.12.0            listenv_0.8.0         
# [25] emmeans_1.4.5          urltools_1.7.3         polyclip_1.10-0        pheatmap_1.0.12        bit64_0.9-7            farver_2.0.3          
# [31] rprojroot_1.3-2        TH.data_1.0-10         coda_0.19-4            vctrs_0.3.6            generics_0.0.2         R6_2.4.1              
# [37] doParallel_1.0.15      graphlayouts_0.7.0     rsvd_1.0.3             locfit_1.5-9.4         gridGraphics_0.5-0     bitops_1.0-6          
# [43] fgsea_1.12.0           assertthat_0.2.1       promises_1.1.0         multcomp_1.4-12        ggraph_2.0.2           enrichplot_1.6.1      
# [49] gtable_0.3.0           npsurv_0.4-0           egg_0.4.5              globals_0.12.5         sandwich_2.5-1         tidygraph_1.1.2       
# [55] rlang_0.4.10           splines_3.6.1          lazyeval_0.2.2         broom_0.7.4            europepmc_0.3          BiocManager_1.30.10   
# [61] reshape2_1.4.3         modelr_0.1.6           backports_1.2.1        httpuv_1.5.2           qvalue_2.18.0          clusterProfiler_3.14.3
# [67] tools_3.6.1            ggplotify_0.0.5        ellipsis_0.3.1         gplots_3.0.3           RColorBrewer_1.1-2     ggridges_0.5.2        
# [73] Rcpp_1.0.4             plyr_1.8.6             progress_1.2.2         RCurl_1.98-1.1         prettyunits_1.1.1      ggpubr_0.2.5          
# [79] pbapply_1.4-2          viridis_0.5.1          cowplot_1.0.0          S4Vectors_0.24.3       zoo_1.8-7              haven_2.2.0           
# [85] ggrepel_0.8.2          cluster_2.1.0          colorRamps_2.3         fs_1.3.2               data.table_1.12.8      DO.db_2.9             
# [91] lmtest_0.9-37          triebeard_0.3.0        reprex_0.3.0           RANN_2.6.1             mvtnorm_1.1-0          fitdistrplus_1.0-14   
# [97] hms_0.5.3              patchwork_1.0.0        lsei_1.2-0             mime_0.9               GSVA_1.34.0            xtable_1.8-4          
# [103] pbkrtest_0.4-8.6       XML_3.99-0.3           readxl_1.3.1           IRanges_2.20.2         gridExtra_2.3          compiler_3.6.1        
# [109] KernSmooth_2.23-16     crayon_1.4.1           minqa_1.2.4            htmltools_0.4.0        later_1.0.0            geneplotter_1.64.0    
# [115] lubridate_1.7.4        DBI_1.1.0              tweenr_1.0.1           dbplyr_1.4.2           MASS_7.3-51.5          boot_1.3-24           
# [121] Matrix_1.2-18          cli_2.3.1              gdata_2.18.0           igraph_1.2.5           pkgconfig_2.0.3        rvcheck_0.1.8         
# [127] plotly_4.9.2           xml2_1.3.2             annotate_1.64.0        estimability_1.3       rvest_0.3.5            digest_0.6.27         
# [133] sctransform_0.2.1      RcppAnnoy_0.0.16       tsne_0.1-3             graph_1.64.0           cellranger_1.1.0       leiden_0.3.3          
# [139] fastmatch_1.1-0        edgeR_3.28.1           uwot_0.1.8             GSEABase_1.48.0        shiny_1.4.0.2          gtools_3.8.1          
# [145] nloptr_1.2.2.1         lifecycle_1.0.0        nlme_3.1-145           jsonlite_1.6.1         viridisLite_0.3.0      fansi_0.4.2           
# [151] pillar_1.5.0           lattice_0.20-40        fastmap_1.0.1          httr_1.4.1             survival_3.1-11        GO.db_3.10.0          
# [157] glue_1.4.2             png_0.1-7              shinythemes_1.1.2      iterators_1.0.12       bit_1.1-15.2           ggforce_0.3.1         
# [163] stringi_1.4.6          blob_1.2.1             org.Hs.eg.db_3.10.0    caTools_1.18.0         memoise_1.1.0          irlba_2.3.3           
# [169] future.apply_1.4.0     ape_5.3        
