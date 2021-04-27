# bulk fold change delta PERSISTENCE SIGNATURES 
suppressMessages(library(here))
suppressMessages(library(tidyverse))
suppressMessages(library(variancePartition))
suppressMessages(library(scglmmr))
suppressMessages(library(magrittr))
'%ni%' = Negate('%in%')

# save path 
datapath = here("microarray_analysis_d100/generated_data_V4_nonadj/"); dir.create(datapath, recursive = T)
figpath = here("microarray_analysis_d100/figures_V4_nonadj/");  dir.create(figpath, recursive = T)

# pattern gene object createed in dir signature_curation 
gp = readRDS(here("signature_curation/H5_GenePatterns.rds"))
pattern_genes = unlist(gp, use.names = FALSE)
persist_pattern = gp[c(4,6:8, 10:12)] 
persist_genes = unlist(x = persist_pattern, use.names = FALSE) %>% unique()

# expression dataframe 
edf = readRDS(here("microarray_analysis_d100/generated_data_V4/array_d0_d100_patterngenes_andmetadata_dataframe.rds"))
met = read_delim(file = "data/metadata/clinical_info_adj.txt",delim = "\t")

# index for tidy gather 
index1 = persist_genes[1]; index2 = persist_genes[length(persist_genes)]

# calculate average of d0 replicates 
eddf = edf %>% 
  # remove subject without day 100 measurement 
  filter(subject.id %ni% c("H5N1-005", "H5N1-003", "H5N1-006")) %>% 
  select(subject.id, time.point, persist_genes) %>% 
  gather(gene, expression, index1:index2) %>%  
  # calculate average expression of the multiple baseline timepoints 
  group_by(subject.id, time.point, gene) %>% 
  summarize(expression = mean(expression, na.rm = TRUE)) %>% 
  arrange(gene, subject.id) %>% 
  ungroup() %>% 
  spread(gene, expression) %>% 
  mutate(sample = paste(subject.id, time.point, sep = "_")) %>% 
  select(sample, everything()) 

## model metadata 
meta = eddf[ ,c('sample', 'subject.id', 'time.point')] %>% 
  mutate(group = plyr::mapvalues(x = subject.id, from = met$`Subject ID`, to = met$Adjuvant), 
         age = plyr::mapvalues(x = subject.id, from = met$`Subject ID`, to = met$Age),
         gender = plyr::mapvalues(x = subject.id, from = met$`Subject ID`, to = met$Gender)) %>% 
  # subset to non-adjuvant subjects
  filter(group == "NonAdj") %>% 
  mutate(age = as.numeric(age)) %>% 
  mutate_if(is.character, as.factor) %>% 
  mutate(scaled_age = as.numeric(scale(age))) %>% 
  column_to_rownames("sample")

# reorder factor levels for custom contrast 
meta$time.point = factor(meta$time.point, levels = c('d000_00h', 'd100'))

# gene expressiondata matrix 
exprs_mat = eddf[ ,c(1, 4:ncol(eddf))] %>%
  # subset to non-adjuvant subjects from the metadata above
  filter(sample %in% rownames(meta)) %>% 
  column_to_rownames("sample") %>%
  t()

## variance partition model 
library(doParallel)
cl = makeCluster(6)
registerDoParallel(cl)

# model formula
f2 <- ~ 0 + time.point + (1|subject.id)

# create model matrix 
time.point = meta$time.point
d100m = model.matrix(~ 0 + time.point) 
contrast_matrix  = makeContrasts(time_100 = (time.pointd100 - time.pointd000_00h), levels = colnames(d100m))

# save custom contrast 
p = plotContrasts(as.data.frame(as.matrix(contrast_matrix))) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
ggsave(p, filename = paste0(figpath, "contrast_matrix.png"), width = 7, height = 7)

# fit model 
fit = variancePartition::dream(exprObj = exprs_mat, formula = f2, data =  meta, L = contrast_matrix, useWeights = FALSE)
fitl = list("bulk" = fit)
res = scglmmr::GetContrastResultsRaw(limma.fit.object.list = fitl, coefficient.number = 1, contrast.name = "delta_delta")
res = res$bulk

## save results 
saveRDS(res, file = paste0(datapath, 'formattedresults_d100_lme4model_array.rds'))
saveRDS(fit, file = paste0(datapath, 'modelfit_d100_lme4model_array.rds'))

################## 
#figgen
termdf = readRDS("signature_curation/cluster_profiler_termdf.rds")

#### hypergeometric test subset of gp genes in DE model upregulated persistence genes 
up_persist = gp[c(4, 6:8)] %>% unlist %>% unique 
rsub = res %>% filter(gene %in% up_persist)
hyp_sub = scglmmr::RunHypergeometricTest(result_list = list(rsub),
                                         TERM2GENE_dataframe = termdf, 
                                         pval_threshold = 0.05,
                                         logFC_threshold = 0,
                                         usefdr_threshold = T)
hyp_plot2 = hyp_sub %>% filter(p.adjust < 0.1)
p = ggplot(hyp_plot2, aes(x = reorder(ID, -log10(p.adjust)), y = -log10(p.adjust), size = gene_ratio)) + 
  geom_point(color = 'blue') +  
  ggtitle("GP 4, 6, 7, 8 \n day 100 Fold change \n non adjuvant group") + 
  xlab(" ") + 
  ylab("hypergeometric test \n -log adjusted p value") + 
  theme_bw() + 
  theme(axis.text = element_text(color = 'black')) + 
  coord_flip() 
p
ggsave(p , filename = paste0(figpath, 'hypergeometric_BTM_d100.pdf'), width = 6, height = 3.3)

# #### Run separately for genes within each gene persistent up pattern gene. 
# for (i in c(4, 6:8)) {
#   rsub = res %>% filter(gene %in% gp[[i]])
#   hyp_sub = scglmmr::RunHypergeometricTest(result_list = list(rsub),
#                                            TERM2GENE_dataframe = termdf, 
#                                            pval_threshold = 0.05,
#                                            logFC_threshold = 0,
#                                            usefdr_threshold = T)
#   
#   hyp_plot2 = hyp_sub %>% filter(p.adjust < 0.1)
#   p = ggplot(hyp_plot2, aes(x = reorder(ID, -log10(p.adjust)), y = -log10(p.adjust), size = gene_ratio)) + 
#     geom_point(color = 'blue') +  
#     ggtitle(paste0("non-adjuvant d100 FC \n GP ", i)) + 
#     xlab(" ") + 
#     ylab("hypergeometric test \n -log adjusted p value") + 
#     theme_gray() + 
#     coord_flip() 
#   
#   ggsave(p , filename = paste0(figpath,"GP_", i, 'hypergeometric_BTM_d100.pdf'), width = 6, height = 3.3)
#   
# }

# plot colored by GP gene
gp_genes = readRDS(file = here("signature_curation/H5_GenePatterns.rds"))
persist_up = unlist(gp_genes[c(4,6:8)]) %>% unlist %>% unique 
de_pattern = res
pattern_up = de_pattern %>% filter(logFC > 0 &  adj.P.Val < 0.05 & gene %in% persist_up)
pattern_up
# create gp_as03: p < 0.05 
gp_as03 = lapply(gp_genes, function(x){ x = x[x %in% unique(pattern_up$gene) ] }) 
gp_as03 = gp_as03[c(4, 6:8)]
# add a string to indicate these were pruned to include genes associated with adjuvant 
names(gp_as03) = paste("AS03",names(gp_as03), sep = "_")

# draw set intersection 
VennDiagram::venn.diagram(gp_as03, 
                          # color 
                          fill = c("#999999", "#E69F00", "#56B4E9", "#009E73"), 
                          # text
                          cat.cex = 0.4,cat.default.pos = "inner",cat.fontfamily = "Helvetica",
                          # file 
                          imagetype="png", height = 1080, width = 1080, resolution = 300,compression = "lzw",
                          filename = paste0(figpath,"gpvenn.png")
)

############################# visualization of microarray results by pattern label intersection 
# get unique set 
vplot = gplots::venn(gp_as03)
xx = attributes(vplot)$intersections
uni10 = xx$`AS03_GE_pattern_04:AS03_GE_pattern_06:AS03_GE_pattern_07:AS03_GE_pattern_08`
unique_list = list(
  "GP4 unique" =  xx$AS03_GE_pattern_04, 
  "GP6 unique" =  xx$AS03_GE_pattern_06,
  "GP7 unique" =  xx$AS03_GE_pattern_07,
  "GP8 unique" =  xx$AS03_GE_pattern_08
)
# map labels 
geneann = list() 
for (i in 1:length(unique_list)) {
  geneann[[i]] = data.frame(gene=unique_list[[i]],module = paste0(names(unique_list[i])))
}
# map other 
geneann = do.call(rbind, geneann) 
patternupfull = unlist(gp_as03, use.names = FALSE)
gene_int = patternupfull[patternupfull %ni%  geneann$gene]
gene_int= data.frame(gene = gene_int, module = "multiple (GP 4,6,7,8)")
geneann = rbind(geneann, gene_int) 
geneann = geneann %>% mutate_if(is.factor, as.character)

## volcano plot with GP enriched colored
rsub = res %>% filter(gene %in% persist_genes)
rsub = rsub %>% mutate(annotation = plyr::mapvalues(gene, from = geneann$gene, to = geneann$module))
rsub = rsub %>% mutate(annotation = if_else(annotation %in% unique(geneann$module), true = annotation, false = ' '))

# add persistence down genes
downpattern = gp[10:12] %>% unlist(use.names = F) %>% unique
rsub = rsub %>% mutate(annotation = if_else(gene %in% downpattern, true = "GP10,11,12", false = annotation))
factorlev = c("GP4 unique", "GP7 unique", "GP8 unique", "GP6 unique", "multiple (GP 4,6,7,8)", "GP10,11,12", " ") 
rsub$annotation = factor(rsub$annotation, levels = factorlev)

# volc plot based on revision comments. 
p = ggplot(rsub, aes(x = logFC , y =  -log10(adj.P.Val), label = gene, color = annotation)) +
  theme_bw() +
  geom_point(data = rsub %>% filter(logFC < 0 & adj.P.Val > 0.05), shape = 16, alpha = 0.8, size = 0.3, color = "black") + 
  geom_point(data = rsub %>% filter(logFC > 0 & adj.P.Val > 0.05), shape = 16, alpha = 0.8, size = 0.3, color = "black") + 
  geom_point(data = rsub %>% filter(adj.P.Val < 0.1),  shape = 16, alpha = 0.7, size = 0.8) +
  scale_color_manual(values = c("#999999", "#56B4E9", "#009E73",  "#E69F00","red", "blue", "black")) + 
  xlab("non adjuvanted day 100 fc") + ylab("-log adjusted p value") +
  guides(color = guide_legend(override.aes = list(size=5, shape =15))) +  
  theme(legend.justification = c(1,1)) + 
  ggrepel::geom_text_repel(data = rsub %>% filter(abs(logFC) > 0.3 & adj.P.Val < 0.05 ),
                           size = 2.5, force = 0.2, show.legend = FALSE, segment.color = "grey", segment.size = 0.1) 
p
ggsave(p, filename = paste0(figpath,"colored_foldchangedelta_contrast.pdf"),width = 6, height = 4)


# behavir of unique genes 
#vplot = gplots::venn(gp_as03)
#xx = attributes(vplot)$intersection
uni10 = xx$`AS03_GE_pattern_04:AS03_GE_pattern_06:AS03_GE_pattern_07:AS03_GE_pattern_08`
topde = res %>% filter(logFC > 0) %>% 
  filter(gene %in% persist_up) %>% 
  arrange(adj.P.Val) %>% 
  top_n(n = 10,wt = -log10(adj.P.Val)) %$% 
  gene 

gene_select = c(uni10, topde)

# index for tidy gather 
index1 = gene_select[1]
index2 = gene_select[length(gene_select)]


eddf2 = eddf %>% 
  select(sample, time.point, subject.id, gene_select) %>% 
  mutate(group = plyr::mapvalues(x = subject.id, from = met$`Subject ID`, to = met$Adjuvant)) %>% 
  select(group, everything()) %>% 
  gather(gene, expression, index1:index2) 
  

# behavor of genes in all persistence patterns 
p = ggplot(eddf2 %>% filter(gene %in% uni10), aes(x = group, y = expression,  color = time.point)) + 
  geom_boxplot(outlier.size = 0.6)+  
  facet_wrap(~gene, nrow = 2, scales = "free") + 
  theme_bw() + 
  theme(strip.background = element_blank()) + 
  ggtitle("upregulated in all patterns in non-adjuvant group model")
ggsave(p, filename = paste0(figpath,"uniquegeneintersect.pdf"), width = 6.5, height = 3) 


# behavor of the top genes
p = ggplot(eddf2 %>% filter(gene %in% topde), aes(x = group, y = expression,  color = time.point)) + 
  geom_boxplot(outlier.size = 0.6)+  
  facet_wrap(~gene, nrow = 2, scales = "free") + 
  theme_bw() + 
  theme(strip.background = element_blank()) + 
  ggtitle("top DE genes in non-adjuvant group model") 
ggsave(p, filename = paste0(figpath,"topdegenes.pdf"), width = 7, height = 3) 

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
#   [1] doParallel_1.0.15        iterators_1.0.12         scglmmr_0.1.0            magrittr_2.0.1          
# [5] variancePartition_1.16.1 Biobase_2.46.0           BiocGenerics_0.32.0      scales_1.1.0            
# [9] foreach_1.4.8            limma_3.42.2             forcats_0.5.0            stringr_1.4.0           
# [13] dplyr_1.0.3              purrr_0.3.3              readr_1.3.1              tidyr_1.0.2             
# [17] tibble_3.1.0             ggplot2_3.3.0            tidyverse_1.3.0          here_0.1                
# 
# loaded via a namespace (and not attached):
#   [1] utf8_1.1.4             reticulate_1.14        tidyselect_1.1.0       lme4_1.1-21           
# [5] RSQLite_2.2.0          AnnotationDbi_1.48.0   htmlwidgets_1.5.1      grid_3.6.1            
# [9] BiocParallel_1.20.1    Rtsne_0.15             munsell_0.5.0          codetools_0.2-16      
# [13] ica_1.0-2              future_1.16.0          withr_2.4.0            colorspace_1.4-1      
# [17] GOSemSim_2.12.1        rstudioapi_0.11        Seurat_3.1.5           stats4_3.6.1          
# [21] ROCR_1.0-7             ggsignif_0.6.0         DOSE_3.12.0            listenv_0.8.0         
# [25] labeling_0.3           emmeans_1.4.5          urltools_1.7.3         polyclip_1.10-0       
# [29] pheatmap_1.0.12        bit64_0.9-7            farver_2.0.3           rprojroot_1.3-2       
# [33] TH.data_1.0-10         coda_0.19-4            vctrs_0.3.6            generics_0.0.2        
# [37] lambda.r_1.2.4         R6_2.4.1               graphlayouts_0.7.0     rsvd_1.0.3            
# [41] locfit_1.5-9.4         gridGraphics_0.5-0     bitops_1.0-6           fgsea_1.12.0          
# [45] assertthat_0.2.1       promises_1.1.0         multcomp_1.4-12        ggraph_2.0.2          
# [49] enrichplot_1.6.1       gtable_0.3.0           npsurv_0.4-0           egg_0.4.5             
# [53] globals_0.12.5         sandwich_2.5-1         tidygraph_1.1.2        rlang_0.4.10          
# [57] splines_3.6.1          lazyeval_0.2.2         europepmc_0.3          broom_0.7.4           
# [61] BiocManager_1.30.10    reshape2_1.4.3         modelr_0.1.6           backports_1.2.1       
# [65] httpuv_1.5.2           qvalue_2.18.0          clusterProfiler_3.14.3 tools_3.6.1           
# [69] ggplotify_0.0.5        ellipsis_0.3.1         gplots_3.0.3           RColorBrewer_1.1-2    
# [73] ggridges_0.5.2         Rcpp_1.0.4             plyr_1.8.6             progress_1.2.2        
# [77] RCurl_1.98-1.1         prettyunits_1.1.1      ggpubr_0.2.5           viridis_0.5.1         
# [81] pbapply_1.4-2          cowplot_1.0.0          S4Vectors_0.24.3       zoo_1.8-7             
# [85] haven_2.2.0            ggrepel_0.8.2          cluster_2.1.0          colorRamps_2.3        
# [89] fs_1.3.2               data.table_1.12.8      futile.options_1.0.1   lmerTest_3.1-1        
# [93] DO.db_2.9              triebeard_0.3.0        lmtest_0.9-37          reprex_0.3.0          
# [97] RANN_2.6.1             mvtnorm_1.1-0          fitdistrplus_1.0-14    hms_0.5.3             
# [101] patchwork_1.0.0        lsei_1.2-0             mime_0.9               GSVA_1.34.0           
# [105] xtable_1.8-4           pbkrtest_0.4-8.6       XML_3.99-0.3           VennDiagram_1.6.20    
# [109] readxl_1.3.1           IRanges_2.20.2         gridExtra_2.3          compiler_3.6.1        
# [113] KernSmooth_2.23-16     crayon_1.4.1           minqa_1.2.4            htmltools_0.4.0       
# [117] later_1.0.0            geneplotter_1.64.0     lubridate_1.7.4        DBI_1.1.0             
# [121] tweenr_1.0.1           formatR_1.7            dbplyr_1.4.2           MASS_7.3-51.5         
# [125] boot_1.3-24            Matrix_1.2-18          cli_2.3.1              gdata_2.18.0          
# [129] igraph_1.2.5           pkgconfig_2.0.3        rvcheck_0.1.8          numDeriv_2016.8-1.1   
# [133] plotly_4.9.2           xml2_1.3.2             annotate_1.64.0        estimability_1.3      
# [137] rvest_0.3.5            digest_0.6.27          sctransform_0.2.1      RcppAnnoy_0.0.16      
# [141] tsne_0.1-3             graph_1.64.0           cellranger_1.1.0       leiden_0.3.3          
# [145] fastmatch_1.1-0        edgeR_3.28.1           uwot_0.1.8             GSEABase_1.48.0       
# [149] shiny_1.4.0.2          gtools_3.8.1           nloptr_1.2.2.1         lifecycle_1.0.0       
# [153] nlme_3.1-145           jsonlite_1.6.1         futile.logger_1.4.3    viridisLite_0.3.0     
# [157] fansi_0.4.2            pillar_1.5.0           lattice_0.20-40        fastmap_1.0.1         
# [161] httr_1.4.1             survival_3.1-11        GO.db_3.10.0           glue_1.4.2            
# [165] png_0.1-7              shinythemes_1.1.2      bit_1.1-15.2           ggforce_0.3.1         
# [169] stringi_1.4.6          blob_1.2.1             org.Hs.eg.db_3.10.0    caTools_1.18.0        
# [173] memoise_1.1.0          irlba_2.3.3            future.apply_1.4.0     ape_5.3   