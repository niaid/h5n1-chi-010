# H1N1 differential expression testing and gene set enrichment analysis
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(Seurat))
suppressMessages(library(magrittr))
suppressMessages(library(scglmmr))
'%ni%' = Negate('%in%')
#suppressMessages(library(variancePartition))

# get all DE functions
# source("de_workflow-master/full_differential_expression_workflow_functions.R")
number_of_cores = 6

# dir = here("mid_res/2_H1_H5_joint_limmaDE/")
datapath = here("3_pseudobulk_de_workflow/generated_data_v2/"); dir.create(datapath)
figpath = here("3_pseudobulk_de_workflow/figures_v2/"); dir.create(figpath)

# read combined object 
s = readRDS(file = "data/day0_dayC_H5N1_seurat_norm_ann.rds")

# define counts and metadata and subset to cells above rm seurat object from workspace 
meta = s@meta.data
umi = s@data

# remove seurat object 
rm(s); gc()

# QC contingency of cells by subject for each celltype 
tab = scglmmr::SubjectCelltypeTable(metadata = meta,
                                    celltype_column = "celltype_joint",
                                    sample_column = "sample")

# remove cells prior to pseudobulk analysis 
cells_keep = meta %>% 
  rownames_to_column("bc") %>% 
  filter(celltype_joint %ni% c(tab$celltypes_remove, 'IgA_CD14_Mono')) %$% 
  bc

# subset data 
meta = meta[cells_keep, ]
umi = umi[ ,cells_keep]

# make pb list and design matrix 
pb = scglmmr::PseudobulkList(rawcounts = umi, metadata = meta, sample_col = "sample", celltype_col = "celltype_joint", avg_or_sum = "sum")
av = scglmmr::PseudobulkList(rawcounts = umi, metadata = meta, sample_col = "sample", celltype_col = "celltype_joint", avg_or_sum = "average")
designmat = scglmmr::BulkDesignMatrix(metadata = meta, sample_column = "sample",variable_column = "timepoint", pseudobulklist = pb)
dge = scglmmr::NormalizePseudobulk(pseudobulklist = pb, design_matrix = designmat, minimum.gene.count = 5)


# make contrast for timepoint effect
c_mat = limma::makeContrasts(time1_foldchange = (timepointdC  -  timepointd0),
                             # add dummy value for baseline so function plot conrtrast will work.   
                             #d0 = timepointd0,
  levels = colnames(designmat)
)

# fit mixed model for the multi timepoint contrasts 
fit = scglmmr::dreamMixedModel(dge_lists = dge,
                               apriori_contrasts = TRUE, 
                               sample_column = 'sample',
                               cell_metadata = meta,
                               contrast_matrix = c_mat, 
                               design_matrix = designmat, 
                               fixed_effects = 'timepoint', 
                               lme4_formula =  '~ 0 + timepoint + (1|sampleid)', 
                               plotsavepath = figpath,
                               version = "2",
                               ncores = number_of_cores)
saveRDS(object = fit, file = paste0(datapath,'d100scglmmrfit.rds'))
saveRDS(object = pb, file = paste0(datapath,'d100mod_pbulklist.rds'))
saveRDS(object = av, file = paste0(datapath,'d100mod_avlist.rds'))
saveRDS(object = meta, file = paste0(datapath,'d100mod_meta.rds'))
############################# v2 stop
sessionInfo()
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
#   [1] scglmmr_0.1.0   magrittr_2.0.1  Seurat_2.3.4    Matrix_1.2-15   cowplot_0.9.4   here_0.1        forcats_0.4.0   stringr_1.4.0   dplyr_0.8.5     purrr_0.3.3    
# [11] readr_1.3.1     tidyr_1.0.2     tibble_2.1.1    ggplot2_3.1.1   tidyverse_1.2.1
# 
# loaded via a namespace (and not attached):
#   [1] estimability_1.3         prabclus_2.3-1           R.methodsS3_1.7.1        coda_0.19-2              acepack_1.4.1            bit64_0.9-7             
# [7] knitr_1.23               irlba_2.3.3              multcomp_1.4-10          R.utils_2.8.0            data.table_1.12.2        rpart_4.1-13            
# [13] RCurl_1.95-4.12          doParallel_1.0.14        generics_0.0.2           metap_1.1                snow_0.4-3               BiocGenerics_0.28.0     
# [19] TH.data_1.0-10           RSQLite_2.1.1            RANN_2.6.1               europepmc_0.3            proxy_0.4-23             bit_1.1-14              
# [25] enrichplot_1.2.0         xml2_1.2.0               lubridate_1.7.4          httpuv_1.5.1             assertthat_0.2.1         viridis_0.5.1           
# [31] xfun_0.7                 hms_0.4.2                promises_1.0.1           DEoptimR_1.0-8           progress_1.2.2           caTools_1.17.1.2        
# [37] readxl_1.3.1             igraph_1.2.4.1           DBI_1.0.0                geneplotter_1.60.0       htmlwidgets_1.3          stats4_3.5.3            
# [43] ggpubr_0.2               backports_1.1.4          annotate_1.60.1          gbRd_0.4-11              vctrs_0.2.4              Biobase_2.42.0          
# [49] ROCR_1.0-7               withr_2.1.2              ggforce_0.2.2            packrat_0.5.0            triebeard_0.3.0          robustbase_0.93-5       
# [55] checkmate_1.9.3          emmeans_1.3.4            prettyunits_1.0.2        mclust_5.4.5             cluster_2.0.7-1          DOSE_3.8.2              
# [61] ape_5.3                  segmented_0.5-4.0        lazyeval_0.2.2           crayon_1.3.4             hdf5r_1.2.0              edgeR_3.24.3            
# [67] pkgconfig_2.0.2          tweenr_1.0.1             nlme_3.1-137             nnet_7.3-12              rlang_0.4.5              diptest_0.75-7          
# [73] lifecycle_0.1.0          sandwich_2.5-1           doSNOW_1.0.16            modelr_0.1.4             cellranger_1.1.0         rprojroot_1.3-2         
# [79] polyclip_1.10-0          GSVA_1.30.0              lmtest_0.9-37            graph_1.60.0             urltools_1.7.3           boot_1.3-20             
# [85] zoo_1.8-6                base64enc_0.1-3          ggridges_0.5.1           pheatmap_1.0.12          png_0.1-7                viridisLite_0.3.0       
# [91] bitops_1.0-6             R.oo_1.22.0              KernSmooth_2.23-15       blob_1.1.1               lars_1.2                 qvalue_2.14.1           
# [97] gridGraphics_0.4-1       S4Vectors_0.20.1         scales_1.0.0             memoise_1.1.0            GSEABase_1.44.0          plyr_1.8.4              
# [103] ica_1.0-2                gplots_3.0.1.1           bibtex_0.4.2             gdata_2.18.0             compiler_3.5.3           lsei_1.2-0              
# [109] RColorBrewer_1.1-2       lme4_1.1-21              fitdistrplus_1.0-14      cli_1.1.0                dtw_1.20-1               pbapply_1.4-0           
# [115] htmlTable_1.13.1         Formula_1.2-3            MASS_7.3-51.1            tidyselect_0.2.5         stringi_1.4.3            GOSemSim_2.8.0          
# [121] locfit_1.5-9.1           latticeExtra_0.6-28      ggrepel_0.8.1            grid_3.5.3               fastmatch_1.1-0          tools_3.5.3             
# [127] parallel_3.5.3           rstudioapi_0.10          foreach_1.4.4            foreign_0.8-71           gridExtra_2.3            farver_1.1.0            
# [133] Rtsne_0.15               ggraph_1.0.2             digest_0.6.25            rvcheck_0.1.3            shiny_1.3.2              fpc_2.2-1               
# [139] Rcpp_1.0.1               broom_0.5.2              egg_0.4.5                SDMTools_1.1-221.1       later_0.8.0              org.Hs.eg.db_3.7.0      
# [145] httr_1.4.0               AnnotationDbi_1.44.0     npsurv_0.4-0             kernlab_0.9-27           Rdpack_0.11-0            colorspace_1.4-1        
# [151] rvest_0.3.4              XML_3.98-1.19            reticulate_1.12          IRanges_2.16.0           splines_3.5.3            shinythemes_1.1.2       
# [157] ggplotify_0.0.3          flexmix_2.3-15           xtable_1.8-4             jsonlite_1.6             nloptr_1.2.1             UpSetR_1.4.0            
# [163] modeltools_0.2-22        R6_2.4.0                 Hmisc_4.2-0              pillar_1.4.1             htmltools_0.3.6          mime_0.6                
# [169] glue_1.3.1               minqa_1.2.4              clusterProfiler_3.10.1   BiocParallel_1.16.6      class_7.3-15             codetools_0.2-16        
# [175] fgsea_1.8.0              tsne_0.1-3               mvtnorm_1.0-10           lattice_0.20-38          pbkrtest_0.4-7           mixtools_1.1.0          
# [181] colorRamps_2.3           gtools_3.8.1             GO.db_3.7.0              survival_2.43-3          limma_3.38.3             munsell_0.5.0           
# [187] DO.db_2.9                iterators_1.0.10         variancePartition_1.12.3 haven_2.1.0              reshape2_1.4.3           gtable_0.3.0 
