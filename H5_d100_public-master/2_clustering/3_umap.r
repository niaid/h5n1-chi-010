# umap of joint clustering results (adapted from fsc repository for day 100 analysis )
suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
suppressMessages(library(here))
'%ni%' = Negate('%in%')
# load python umap package via reticulate 
library(reticulate)
use_virtualenv("r-reticulate")
library(umap)

# set file paths 
figpath = here("2_clustering/figures/")
datapath = here("2_clustering/generated_data/")
dir.create(figpath) ; dir.create(datapath)

# load full dataset 
h1h5 = readRDS("data/day0_dayC_H5N1_seurat_norm_ann.rds")
h1h5 = SetAllIdent(h1h5, id = "celltype_joint")


## protein heatmap 
source(file = here("functions/protein_annotation_functions.r"))
protdf = AggregateProteinDataPlot(protein_assay_data = h1h5@assay$CITE@data, 
                                  rownamed_metadata = h1h5@meta.data,
                                  mdname = "celltype_joint")
pheatmap::pheatmap(t(protdf), treeheight_col = 30,treeheight_row = 30,
                   color = viridis::viridis(n = 30, option = "B"), 
                   border_color = NA, width = 13, height = 6,
                   filename = paste0(figpath, "proteinheatmap.pdf"))


# most informative proteins 
protdfsub = protdf[c("CD3", "CD4", "CD8", "CD19", "CD14", "CD16", "CD303", "CD123", "CD34", 
                     "CD27", "CD45RA", "CD45RO", "CD161", "CD56", "CD127", "CD38", "CD11c", "CD11b"), ]
mtx = t(protdfsub)
rownames(mtx) = str_replace_all(rownames(mtx), pattern = '_', replacement = " ")
pheatmap::pheatmap(mtx, treeheight_col = 15,treeheight_row = 15,
                   color = viridis::viridis(n = 15, option = "B"),
                   border_color = NA,
                   width = 6, height = 4,
                   filename = paste0(figpath, "proteinheatmap_subset.pdf"))



# get protein data matrix without isotype controls 
adt = h1h5@assay$CITE@data
proteins = rownames(h1h5@assay$CITE@data)
proteins = proteins[-c(67:70, 19)]
adt = t(adt[proteins, ])


# set umap configuration 
config = umap.defaults
config$n_neighbors = 50
config$min_dist = 0.4

# run umap 
rm(h1h5); gc()
ump = umap(adt, config = config)

# save umap output
saveRDS(ump, file = paste0(datapath,"umap_output.RDS"))
saveRDS(config, file = paste0(datapath,"umap_config.RDS"))

ump = readRDS(file = "2_clustering/generated_data/umap_output.RDS")

# read back data to get metadata 
s = readRDS("data/day0_dayC_H5N1_seurat_norm_ann.rds")
md = s@meta.data

# get umap embeddings for plotting 
embeddings = ump$layout %>% as.data.frame
names(embeddings) = c("UMAP1", "UMAP2")

#mege 
df = cbind(embeddings, md)

# df$celltype_joint = factor(df$celltype_joint, levels = celltypes)
# col = rev(pals::stepped(n = length(celltypes)))
 

library(ggrepel)
centers = df %>% dplyr::group_by(celltype_joint) %>% summarize(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2))
centers$celltype_joint = str_replace_all(string = centers$celltype_joint, pattern = "_",replacement = " ")
centers = centers %>% filter(celltype_joint %ni% "DOUBLET")
df$celltype_joint = str_replace_all(string = df$celltype_joint, pattern = "_",replacement = " ")


cu = ggsci::pal_d3(palette = 'category20')(20)
cu = c(cu, 'black', 'green', 'lightgrey' )
p = ggplot(df %>% filter(celltype_joint %ni% 'DOUBLET'), aes(x = UMAP1,y =UMAP2, color = celltype_joint)) +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
      panel.grid.major = element_blank(), 
      axis.ticks.y = element_blank(),
      axis.ticks.x = element_blank()
  ) + 
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  geom_point(show.legend = TRUE, size = 0.6, shape = 16,  alpha = 0.8) +
  scale_color_manual(values = cu) +
  geom_label_repel(data = centers, seed = 1, color = "black", box.padding = 1,
                  size = 3,
                  segment.color = "grey", fontface = 'bold', segment.size = 0.3,force = 3,
                  # nudge_y = 0.5,
                  aes(label = celltype_joint ),
                  show.legend = FALSE) +
  theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8)) +
  theme(legend.position = "bottom") +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank())  +
  theme(legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size=4, alpha = 0.7, shape = 15))) +
  theme(legend.text = element_text(colour="black", size=8, face="bold"))
p
ggsave(p, filename = paste0(figpath,"H1N1_umap_celltype_joint.png"), width = 8, height = 8)
ggsave(p, filename = paste0(figpath,"H1N1_umap_celltype_joint.pdf"), width = 8, height = 8)

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
#   [1] ggrepel_0.8.1   umap_0.2.3.1    reticulate_1.12 here_0.1        forcats_0.4.0   stringr_1.4.0   dplyr_0.8.5     purrr_0.3.3     readr_1.3.1    
# [10] tidyr_1.0.2     tibble_2.1.1    tidyverse_1.2.1 Seurat_2.3.4    Matrix_1.2-15   cowplot_0.9.4   ggplot2_3.1.1  
# 
# loaded via a namespace (and not attached):
#   [1] readxl_1.3.1            snow_0.4-3              backports_1.1.4         Hmisc_4.2-0             plyr_1.8.4              igraph_1.2.4.1         
# [7] lazyeval_0.2.2          splines_3.5.3           crosstalk_1.0.0         digest_0.6.25           foreach_1.4.4           htmltools_0.3.6        
# [13] viridis_0.5.1           lars_1.2                gdata_2.18.0            magrittr_2.0.1          checkmate_1.9.3         cluster_2.0.7-1        
# [19] mixtools_1.1.0          ROCR_1.0-7              modelr_0.1.4            R.utils_2.8.0           askpass_1.1             colorspace_1.4-1       
# [25] rvest_0.3.4             haven_2.1.0             xfun_0.7                crayon_1.3.4            jsonlite_1.6            survival_2.43-3        
# [31] zoo_1.8-6               iterators_1.0.10        ape_5.3                 glue_1.3.1              pals_1.5                gtable_0.3.0           
# [37] webshot_0.5.1           kernlab_0.9-27          prabclus_2.3-1          DEoptimR_1.0-8          maps_3.3.0              scales_1.0.0           
# [43] pheatmap_1.0.12         mvtnorm_1.0-10          bibtex_0.4.2            miniUI_0.1.1.1          Rcpp_1.0.1              metap_1.1              
# [49] dtw_1.20-1              xtable_1.8-4            viridisLite_0.3.0       htmlTable_1.13.1        foreign_0.8-71          bit_1.1-14             
# [55] mapproj_1.2.6           proxy_0.4-23            mclust_5.4.5            SDMTools_1.1-221.1      Formula_1.2-3           stats4_3.5.3           
# [61] tsne_0.1-3              htmlwidgets_1.3         httr_1.4.0              gplots_3.0.1.1          RColorBrewer_1.1-2      fpc_2.2-1              
# [67] acepack_1.4.1           modeltools_0.2-22       ica_1.0-2               pkgconfig_2.0.2         R.methodsS3_1.7.1       flexmix_2.3-15         
# [73] nnet_7.3-12             labeling_0.3            manipulateWidget_0.10.0 later_0.8.0             tidyselect_0.2.5        rlang_0.4.5            
# [79] reshape2_1.4.3          munsell_0.5.0           cellranger_1.1.0        tools_3.5.3             cli_1.1.0               generics_0.0.2         
# [85] broom_0.5.2             ggridges_0.5.1          npsurv_0.4-0            knitr_1.23              bit64_0.9-7             fitdistrplus_1.0-14    
# [91] robustbase_0.93-5       rgl_0.100.30            caTools_1.17.1.2        RANN_2.6.1              packrat_0.5.0           pbapply_1.4-0          
# [97] nlme_3.1-137            mime_0.6                R.oo_1.22.0             xml2_1.2.0              hdf5r_1.2.0             compiler_3.5.3         
# [103] rstudioapi_0.10         png_0.1-7               lsei_1.2-0              stringi_1.4.3           lattice_0.20-38         ggsci_2.9              
# [109] vctrs_0.2.4             pillar_1.4.1            lifecycle_0.1.0         Rdpack_0.11-0           lmtest_0.9-37           data.table_1.12.2      
# [115] bitops_1.0-6            irlba_2.3.3             gbRd_0.4-11             httpuv_1.5.1            R6_2.4.0                latticeExtra_0.6-28    
# [121] promises_1.0.1          KernSmooth_2.23-15      gridExtra_2.3           codetools_0.2-16        dichromat_2.0-0         MASS_7.3-51.1          
# [127] gtools_3.8.1            assertthat_0.2.1        openssl_1.4             rprojroot_1.3-2         withr_2.1.2             diptest_0.75-7         
# [133] parallel_3.5.3          doSNOW_1.0.16           hms_0.4.2               grid_3.5.3              rpart_4.1-13            class_7.3-15           
# [139] segmented_0.5-4.0       Rtsne_0.15              shiny_1.3.2             lubridate_1.7.4         base64enc_0.1-3 