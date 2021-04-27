library(tidyverse)
suppressMessages(library(here))
'%ni%' = Negate('%in%')
resfull = readRDS(file = here("microarray_analysis_d100/generated_data_V4/formattedresults_d100_lme4model_array.rds"))
resnonadj = readRDS(file = here("microarray_analysis_d100/generated_data_V4_nonadj/formattedresults_d100_lme4model_array.rds"))
resadj = readRDS(file = here("microarray_analysis_d100/generated_data_V3_adjonly/formattedresults_d100_lme4model_array.rds"))

# rename specific values for each contrast. 
resfull = resfull %>%
  rename(logFC_day100fc = logFC, pval_day100fc = P.Value, padjust_day100fc = adj.P.Val) %>% 
  select(gene, logFC_day100fc, pval_day100fc, padjust_day100fc)
resadj = resadj %>% 
  rename(logFC_AS03 = logFC, pval_AS03 = P.Value, padjust_AS03 = adj.P.Val) %>% 
  select(gene, logFC_AS03, pval_AS03, padjust_AS03)
resnonadj = resnonadj %>%
  rename(logFC_nonAdj = logFC, pval_nonAdj = P.Value, padjust_nonAdj = adj.P.Val) %>% 
  select(gene, logFC_nonAdj, pval_nonAdj, padjust_nonAdj)

res_merge = full_join(resfull, resadj, by = "gene")
res_merge = full_join(res_merge, resnonadj, by = "gene")


### add gene pattern anotation 
# plot colored by GP gene
gp_genes = readRDS(file = here("signature_curation/H5_GenePatterns.rds"))
persist_patterns = gp_genes[c(4,6:8)]
persist_down = gp_genes[c(10:12)]
persist_down = unlist(persist_down) %>% unique()
persist_patterns$`GE_pattern_10-12` = persist_down
persist_up = unlist(gp_genes[c(4,6:8)]) %>% unlist %>% unique 

# get unique set 
vplot = gplots::venn(persist_patterns)
xx = attributes(vplot)$intersections
unique_list = list(
  "GP4 unique" =  xx$GE_pattern_04, 
  "GP6 unique" =  xx$GE_pattern_06,
  "GP7 unique" =  xx$GE_pattern_07,
  "GP8 unique" =  xx$GE_pattern_08,
  "GP10-12" = xx$`GE_pattern_10-12`
)
# map labels 
geneann = list() 
for (i in 1:length(unique_list)) {
  geneann[[i]] = data.frame(gene=unique_list[[i]], module = paste0(names(unique_list[i])))
}
# map other 
geneann = do.call(rbind, geneann) 
gene_int = persist_up[persist_up %ni%  geneann$gene]
gene_int = data.frame(gene = gene_int, module = "multiple (GP 4,6,7,8)")
geneann = rbind(geneann, gene_int) 
geneann = geneann %>% mutate_if(is.factor, as.character)

# add annotation to list 
res_merge = full_join(res_merge, geneann, by = "gene")
res_merge = res_merge %>% select(gene, module, everything())
write_delim(res_merge, path = paste0(here("microarray_analysis_d100/merged_table_day100models.txt")), delim = "\t")


# add plot 
factorlev = c("GP4 unique", "GP7 unique", "GP8 unique", "GP6 unique", "multiple (GP 4,6,7,8)", "GP10-12") 
res_merge$annotation = factor(res_merge$module, levels = factorlev)
p =
ggplot(res_merge, aes(x = logFC_AS03, y = logFC_nonAdj, label = gene)) + 
  geom_point(aes(color = annotation), shape = 16, size = 0.5, alpha = 0.8, show.legend = FALSE) + 
  scale_color_manual(values = c("#999999", "#56B4E9", "#009E73",  "#E69F00","red", "blue", "black")) + 
  ggrepel::geom_text_repel(data = res_merge %>% filter(abs(logFC_AS03) > 0.45 | abs(logFC_nonAdj) > 0.45),
                          show.legend = FALSE, inherit.aes = TRUE, size = 2, force = 5, segment.colour = "grey", segment.size = 0.2) + 
  facet_wrap(~module) + 
  theme_bw() + 
  ggpubr::stat_cor(method = "pearson", label.y.npc = 0.05, ) +   
  xlab("log fold change AS03 group") + ylab("log fold change non-adjuvant group") + 
  theme(strip.background = element_blank()) +
  geom_vline(xintercept = 0, size = 0.2) + geom_hline(yintercept = 0, size = 0.2) 
p
ggsave(p, filename = here("microarray_analysis_d100/figures_V4/comparison_foldchange.pdf"),width = 7, height = 5)
  

 #


## leading edge genes 
library(scglmmr)
gsea = readRDS(file = here("3_pseudobulk_de_workflow/generated_data_v2/d100_gsea_gp_as03.rds"))
lef = GetLeadingEdgeFull(gsea.list = gsea,padj.filter = 1,NES.filter = -Inf)
lef = do.call(rbind, lef)
write_delim(lef, path = here("3_pseudobulk_de_workflow/generated_data_v2/ledge_full_table_d100_dreamlmer.txt"),delim = "\t")
scglmmr::CombineResults(gsealist = gsea, contrastlist = res)


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
#   [1] doParallel_1.0.15        iterators_1.0.12         magrittr_2.0.1           scglmmr_0.1.0            variancePartition_1.16.1
# [6] Biobase_2.46.0           BiocGenerics_0.32.0      scales_1.1.0             foreach_1.4.8            limma_3.42.2            
# [11] forcats_0.5.0            stringr_1.4.0            dplyr_1.0.3              purrr_0.3.3              readr_1.3.1             
# [16] tidyr_1.0.2              tibble_3.1.0             ggplot2_3.3.0            tidyverse_1.3.0          here_0.1                
# 
# loaded via a namespace (and not attached):
#   [1] utf8_1.1.4             reticulate_1.14        tidyselect_1.1.0       lme4_1.1-21            RSQLite_2.2.0         
# [6] AnnotationDbi_1.48.0   htmlwidgets_1.5.1      grid_3.6.1             BiocParallel_1.20.1    Rtsne_0.15            
# [11] munsell_0.5.0          codetools_0.2-16       ica_1.0-2              future_1.16.0          withr_2.4.0           
# [16] colorspace_1.4-1       GOSemSim_2.12.1        knitr_1.28             rstudioapi_0.11        Seurat_3.1.5          
# [21] stats4_3.6.1           ROCR_1.0-7             ggsignif_0.6.0         DOSE_3.12.0            listenv_0.8.0         
# [26] labeling_0.3           emmeans_1.4.5          urltools_1.7.3         polyclip_1.10-0        bit64_0.9-7           
# [31] farver_2.0.3           pheatmap_1.0.12        rprojroot_1.3-2        coda_0.19-4            vctrs_0.3.6           
# [36] generics_0.0.2         TH.data_1.0-10         xfun_0.12              lambda.r_1.2.4         R6_2.4.1              
# [41] graphlayouts_0.7.0     rsvd_1.0.3             locfit_1.5-9.4         bitops_1.0-6           fgsea_1.12.0          
# [46] gridGraphics_0.5-0     assertthat_0.2.1       promises_1.1.0         multcomp_1.4-12        ggraph_2.0.2          
# [51] enrichplot_1.6.1       gtable_0.3.0           npsurv_0.4-0           egg_0.4.5              globals_0.12.5        
# [56] tidygraph_1.1.2        sandwich_2.5-1         rlang_0.4.10           splines_3.6.1          lazyeval_0.2.2        
# [61] broom_0.7.4            europepmc_0.3          modelr_0.1.6           BiocManager_1.30.10    reshape2_1.4.3        
# [66] backports_1.2.1        httpuv_1.5.2           qvalue_2.18.0          clusterProfiler_3.14.3 tools_3.6.1           
# [71] ggplotify_0.0.5        ellipsis_0.3.1         gplots_3.0.3           RColorBrewer_1.1-2     ggridges_0.5.2        
# [76] Rcpp_1.0.4             plyr_1.8.6             progress_1.2.2         RCurl_1.98-1.1         prettyunits_1.1.1     
# [81] ggpubr_0.2.5           pbapply_1.4-2          viridis_0.5.1          cowplot_1.0.0          S4Vectors_0.24.3      
# [86] zoo_1.8-7              haven_2.2.0            ggrepel_0.8.2          cluster_2.1.0          fs_1.3.2              
# [91] colorRamps_2.3         futile.options_1.0.1   data.table_1.12.8      lmerTest_3.1-1         DO.db_2.9             
# [96] reprex_0.3.0           lmtest_0.9-37          triebeard_0.3.0        RANN_2.6.1             mvtnorm_1.1-0         
# [101] fitdistrplus_1.0-14    hms_0.5.3              patchwork_1.0.0        lsei_1.2-0             mime_0.9              
# [106] GSVA_1.34.0            xtable_1.8-4           pbkrtest_0.4-8.6       XML_3.99-0.3           VennDiagram_1.6.20    
# [111] readxl_1.3.1           IRanges_2.20.2         gridExtra_2.3          compiler_3.6.1         KernSmooth_2.23-16    
# [116] crayon_1.4.1           minqa_1.2.4            htmltools_0.4.0        later_1.0.0            geneplotter_1.64.0    
# [121] lubridate_1.7.4        DBI_1.1.0              formatR_1.7            tweenr_1.0.1           dbplyr_1.4.2          
# [126] MASS_7.3-51.5          boot_1.3-24            Matrix_1.2-18          cli_2.3.1              gdata_2.18.0          
# [131] igraph_1.2.5           pkgconfig_2.0.3        rvcheck_0.1.8          numDeriv_2016.8-1.1    plotly_4.9.2          
# [136] xml2_1.3.2             annotate_1.64.0        estimability_1.3       rvest_0.3.5            digest_0.6.27         
# [141] sctransform_0.2.1      RcppAnnoy_0.0.16       tsne_0.1-3             graph_1.64.0           cellranger_1.1.0      
# [146] leiden_0.3.3           fastmatch_1.1-0        uwot_0.1.8             edgeR_3.28.1           GSEABase_1.48.0       
# [151] shiny_1.4.0.2          gtools_3.8.1           nloptr_1.2.2.1         lifecycle_1.0.0        nlme_3.1-145          
# [156] jsonlite_1.6.1         futile.logger_1.4.3    viridisLite_0.3.0      fansi_0.4.2            pillar_1.5.0          
# [161] lattice_0.20-40        fastmap_1.0.1          httr_1.4.1             survival_3.1-11        GO.db_3.10.0          
# [166] glue_1.4.2             png_0.1-7              shinythemes_1.1.2      bit_1.1-15.2           ggforce_0.3.1         
# [171] stringi_1.4.6          blob_1.2.1             org.Hs.eg.db_3.10.0    caTools_1.18.0         memoise_1.1.0         
# [176] irlba_2.3.3            future.apply_1.4.0     ape_5.3 
