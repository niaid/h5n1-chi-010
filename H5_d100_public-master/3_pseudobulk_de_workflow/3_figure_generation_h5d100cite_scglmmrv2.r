# H1N1 differential expression testing and gene set enrichment analysis
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(Seurat))
suppressMessages(library(magrittr))
suppressMessages(library(scglmmr))
'%ni%' = Negate('%in%')

# save paths 
datapath = here("3_pseudobulk_de_workflow/generated_data_v2/"); dir.create(datapath)
figpath = here("3_pseudobulk_de_workflow/figures_v2/"); dir.create(figpath)

########################
# load microarray results for genes in GPs with adjuvant group > unadjuvant group 
gp_genes = readRDS(file = here("signature_curation/H5_GenePatterns.rds"))
persist_up = unlist(gp_genes[c(4,6:8)]) %>% unlist %>% unique 
de_pattern = readRDS(file = here("microarray_analysis_d100/generated_data_V4/formattedresults_d100_lme4model_array.rds"))
pattern_up = de_pattern %>% filter(adj.P.Val < 0.01 & logFC > 0 & gene %in% persist_up)

# create gp_as03: p < 0.05 from bulk adjuvant day 100 vs unadjuvanted 
gp_as03 = lapply(gp_genes, function(x){ x = x[x %in% unique(pattern_up$gene) ] }) 
gp_as03 = gp_as03[c(4, 6:8)]
names(gp_as03) = c("Gp04", "Gp06", "Gp07", "Gp08")


########################
# CITE-seq pseudobulk results within Adjuvanted subjects - day 100 vs day 0 
fit = readRDS(here("3_pseudobulk_de_workflow/generated_data_v2/d100scglmmrfit.rds"))
res = scglmmr::GetContrastResultsRaw(limma.fit.object.list = fit, coefficient.number = 1, contrast.name = "d100")
ranks = scglmmr::GetRankResultsRaw(contrast.result.raw.list = res)

# run Gene Set Enrichment on the AS03_GenePatterns (adjuvant enriched geens within GP6-8 10-12) 
gsea_gp = scglmmr::RunFgseaOnRankList(rank.list.celltype = ranks, pathways = gp_as03, positive.enrich.only = FALSE)
gseabind = scglmmr::RbindGseaResultList(gsea_result_list = gsea_gp, NES_filter = -Inf, padj_filter = 0.1)
saveRDS(gsea_gp, file = paste0(datapath, "d100_gsea_gp_as03.rds"))

# bubbleplot aes.
plot_param = list (
  geom_point(shape = 21, fill = "darkred"),
  theme_bw(),
  scale_x_discrete(position = "top"),
  theme(legend.justification = c(1,1)),
  theme(axis.text.x=element_text(angle = 45, hjust = 0)),
  theme(axis.title.y = element_blank()),
  xlab(""), 
  labs(fill = 'Normalized \n Enrichment \n Score', size = '-log10(adjusted p)'),
  theme(legend.title = element_text(colour = "black", size = 8)),
  theme(axis.text.y = element_text(size = 12,  color = "black")),
  theme(axis.text.x = element_text(size = 12,  color = "black")),
  guides(shape = guide_legend(override.aes = list(size = 5))),
  guides(color = guide_legend(override.aes = list(size = 5))),
  theme(legend.title = element_text(size = 8), legend.text = element_text(size = 8))
)

# plot 
gseabind$celltype = str_replace_all(string = gseabind$celltype, pattern = "_",replacement = " ")
p = ggplot(gseabind, aes(y = pathway, x = celltype, size = n_logp)) + plot_param 
ggsave(p, filename = paste0(figpath, "gsea_persistence_deinbulkv2.pdf"), width = 5, height = 2.8)



# unbiased GSEA analysis on BTM based on day 100 coefficient 
btm = readRDS(file = here("signature_curation/BTM_li.rds"))
gsea_btm = scglmmr::RunFgseaOnRankList(rank.list.celltype = ranks, pathways = btm, positive.enrich.only = FALSE)
gsea_rb = scglmmr::RbindGseaResultList(gsea_result_list = gsea_btm, NES_filter = -Inf, padj_filter = 0.05)
gsea_rb = gsea_rb %>% filter(celltype %ni% 'DOUBLET')
scglmmr::GSEABubblePlot(rbind_gsea_result_dataframe = gsea_rb, include_negative = T, width = 6.5, height = 4.2,
                        save_path = figpath, save_name = "gsea_unbiased_100.pdf")
# save enrichment results 
saveRDS(gsea_btm, file = paste0(datapath, "d100_gsea_btm.rds"))
gsea_btm = readRDS("3_pseudobulk_de_workflow/generated_data_v2/d100_gsea_btm.rds")


# leading edge heatmaps - BTM 
av = readRDS(file = here("3_pseudobulk_de_workflow/generated_data_v2/d100mod_avlist.rds"))
meta = readRDS(file = here("3_pseudobulk_de_workflow/generated_data_v2/d100mod_meta.rds"))
le_expr = scglmmr::LeadEdgeTidySampleExprs(av.exprs.list = av, gsea.list = gsea_btm, padj.filter = 0.1,NES.filter = -Inf)

for (i in 1:length(le_expr)) {
  cdat = le_expr[[i]]
  ctype = names(le_expr[i])
  umod = unique(cdat$module)
  for (u in 1:length(umod)) {
    scglmmr::LeadEdgeSampleHeatmap(tidy.exprs.list = le_expr, 
                                   modulename = umod[u], 
                                   celltype_plot = ctype,
                                   metadata = meta, 
                                   metadata_annotate = c('timepoint'),
                                   sample_column = 'sample',
                                   returnmat = F, 
                                   savepath = figpath,
                                   savename = paste0(ctype, " ",umod[u],'.pdf'))
  }
}

# leading edge heatmaps - GP 
le_expr = scglmmr::LeadEdgeTidySampleExprs(av.exprs.list = av, gsea.list = gsea_gp, padj.filter = 0.1,NES.filter = -Inf)
for (i in 1:length(le_expr)) {
  cdat = le_expr[[i]]
  ctype = names(le_expr[i])
  umod = unique(cdat$module)
  for (u in 1:length(umod)) {
    scglmmr::LeadEdgeSampleHeatmap(tidy.exprs.list = le_expr, 
                                   modulename = umod[u], 
                                   celltype_plot = ctype,
                                   metadata = meta, 
                                   metadata_annotate = c('timepoint'),
                                   sample_column = 'sample',
                                   returnmat = F, 
                                   savepath = figpath,
                                   savename = paste0(ctype, " ",umod[u],'.pdf'))
  }
}



# gsea dotplot per cell type
gsea_btm = readRDS(here("3_pseudobulk_de_workflow/generated_data_v2/d100_gsea_btm.rds"))

# mono 
p = ggplot(gsea_btm[['CD14_Mono']] %>% filter(padj < 0.05 & NES > 0),
       aes(x = NES, y = reorder(pathway, NES),size = -log10(padj))) + 
  theme_grey() +
  geom_point(shape = 21,  fill = 'purple') + 
  theme(axis.text.y = element_text(face = "bold", color = "black")) +
  theme(axis.title.y = element_blank()) + 
  xlab("Normalized Enrichement Score") + ylab("pathway") 
ggsave(p, filename = paste0(figpath, "CD14Mono_Enr.pdf"), width = 5, height = 2)

# mdc 
p = ggplot(gsea_btm[['mDC']] %>% filter(padj < 0.05 & NES > 0),
           aes(x = NES, y = reorder(pathway, NES),size = -log10(padj))) + 
  xlab("Normalized Enrichement Score") + ylab("pathway") + 
  theme_grey() +
  theme(axis.title.y = element_blank()) + 
  theme(axis.text.y = element_text(face = "bold", color = "black")) +
  geom_point(shape = 21,  fill = 'green') 
ggsave(p, filename = paste0(figpath, "mDC_Enr.pdf"), width = 5, height = 2)  

#cd8 n 
p = ggplot(gsea_btm[["CD8_Naive_Tcell"]] %>% filter(padj < 0.05 & NES > 0),
           aes(x = NES, y = reorder(pathway, NES),size = -log10(padj))) + 
  xlab("Normalized Enrichement Score") + ylab("pathway") + 
  theme_grey() +
  theme(axis.title.y = element_blank()) + 
  theme(axis.text.y = element_text(face = "bold", color = "black")) +
  geom_point(shape = 21,  fill = 'pink') 
p
ggsave(p, filename = paste0(figpath, "cd8N_Enr.pdf"), width = 5, height = 2)  

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
#   [1] scglmmr_0.1.0   magrittr_2.0.1  Seurat_2.3.4    Matrix_1.2-15   cowplot_0.9.4   here_0.1        forcats_0.4.0   stringr_1.4.0  
# [9] dplyr_0.8.5     purrr_0.3.3     readr_1.3.1     tidyr_1.0.2     tibble_2.1.1    ggplot2_3.1.1   tidyverse_1.2.1
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
# 


# 
# 
# # combine with gene coefficients -- gene patterns 
# res = readRDS(here("3_pseudobulk_de_workflow/generated_data/h5_d100timeres.rds"))
# r = CombineResults(gsealist = gsea_gp, contrastlist = res, gseafdr = 0.1, genefdr = 1)
# saveRDS(r,file = paste0(datapath,"leadingedge_gseagp_d100cite.rds"))
# 
# # combine with gene coefficients -- BTM 
# r2 = CombineResults(gsealist = gsea_btm, contrastlist = res, gseafdr = 0.1, genefdr = 1)
# rr = bind_rows(r, r2)
# saveRDS(rr,file = paste0(datapath,"leadingedge_gseaGP_BTM_d100cite.rds"))
# 
# # create list of bulk adjuvant-associated pattern genes 
# sc_gp_up = rr %>% filter(logFC > 0.2 & pathway %in% c("AS03_GP6", "AS03_GP7", "AS03_GP8"))
# sc_gp_down = rr %>% filter(logFC < 0.2 & pathway %in% c("AS03_GP10", "AS03_GP11", "AS03_GP12"))
# sc_gp = bind_rows(sc_gp_up, sc_gp_down)
# # save 
# write_delim(sc_gp, path = paste0(datapath, "as03_associatedbulk_gp_singlecell_leadingedge_filtered.txt"), delim = "\t")
# 
# 
# ######## Combined plot across unbiased and AS03_GP 
# av_y = rr %>%
#   filter(abs(logFC) > 0.4) %>%
#   group_by(gene) %>%
#   summarize(meany = mean(logFC)) %>%
#   arrange(meany) %$% gene
# rr$gene = factor(rr$gene, levels = av_y)
# 
# p = ggplot(rr %>%
#              filter(abs(logFC) > 0.4 ) %>%
#              filter(celltype %in% c( "CD14_Mono",  "mDC", "CD8_Naive_Tcell")),
#            aes(x = pathway, y = gene, fill = logFC)) +
#   geom_point(shape = 21, size = 3) +
#   theme_bw() +
#   theme(axis.text.x=element_text(angle = -90, hjust = 0,color = "black", face= "bold")) +
#   theme(axis.text.y=element_text(size = 6, color = "black", face= "bold")) +
#   scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "red3",midpoint = 0) +
#   facet_grid(~celltype,  scales = "free_x", space = "free_x")
# 
# ggsave(p, filename = paste0(figpath,"persistence_leadingedge_celltype.pdf" ), width = 6, height = 13.4)
# 
# 
# 






















# 
# 
# 
# # Enrichment in these AS03 associated patterns - Microarray data 
# # convert to Entrez IDs 
# demod_entrez = lapply(demod, function(x){ 
#   mapping = AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = x, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL") 
#   mapping = unique(mapping$ENTREZID)
#   mapping = mapping[!is.na(mapping)]; return(mapping)
# })
# # run enrichment within Reactome pathways (hypergeometric)
# reac_enr_bulk = lapply(demod_entrez, function(x){ 
#   ReactomePA::enrichPathway(gene = x)@result  })
# remap = c("gene_num", "gene_denom")
# # reformat results 
# for (i in 1:length(reac_enr_bulk)) {
#   reac_enr_bulk[[i]] %<>% 
#     separate(GeneRatio, into = remap, sep = "/") %>% 
#     mutate_at(.vars = remap, .funs = as.numeric) %>% 
#     mutate(generatio = gene_num / gene_denom) %>% 
#     mutate(Pattern = names(reac_enr_bulk[i]))
# }
# # merge and plot enrichment.
# d = reac_enr_bulk %>% bind_rows() %>% filter(p.adjust < 0.2)
# p = ggplot(d, aes(x = generatio,  y = -log10(p.adjust), label = Description, fill = Pattern)) +
#   theme_bw() + 
#   geom_point(shape = 21, size = 3) + 
#   ggrepel::geom_text_repel(data = d %>% filter(p.adjust < 0.05 & generatio > 0.15), show.legend = F, nudge_y = 0.3) + 
#   ggsci::scale_fill_d3() + 
#   geom_vline(xintercept = 0.15) + geom_hline(yintercept = 1.3) + 
#   ggtitle("Difference in day 100 Fold Change from baseline \n AS03 vs unadjuvanted (Microarray)")
# ggsave(p, filename = paste0(figpath, "microarray_day100GP_hypergeometric_reactome.pdf"), width = 5, height = 3.5)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ##############
# 
# # load gene_pattern_genes
# gpath = list.files(path = here("signature_curation/GE_pattern_genes/"), full.names = TRUE)
# gpnames = list.files(path = here("signature_curation/GE_pattern_genes/"))
# gpnames = str_sub(gpnames, 1,13)
# gp = lapply(gpath, function(x){read_delim(x, delim = "\t") %$% gene })
# names(gp) = gpnames 
# 
# 
# # DS model -- get celltypes ; remove low representation celltypes and doublets ( read cohort model to run same celltypes )
# celltypes  = h5@meta.data$celltype_joint %>% unique 
# 
# mytable = GetSubjectCelltypeTable(Seurat.Object = h5,celltype_column = "celltype_joint",sample_column = "sample")
# remove = mytable$`low representation celltypes`
# remove = c(remove, "DOUBLET")
# remove = remove[-6]
# celltypes = celltypes[celltypes %ni% remove]
# 
# # make pseudobulk list for each celltype 
# sl = MakePseudobulkList(seurat.object = h5, celltype_column = "celltype_joint", 
#                         sample_column = "sample",vector.of.celltypes = celltypes)
# 
# # metadtat table for model matrix 
# d100 = MakeMetaTableFromSeurat(seurat_object = h5,sample_column = 'sample',variable_column = "timepoint",
#                                aggregate.data.list = sl,celltypes.vector = celltypes) %>% 
#   select(sample, timepoint) %>%
#   mutate(sampleid = str_sub(sample,-6,-4)) %>% 
#   column_to_rownames("sample")
# 
# d1m = d100 %>% mutate_if(is.character, as.factor) %$% timepoint
# d1m = model.matrix(~0 + d1m)
# colnames(d1m) = str_sub(colnames(d1m), start = 4, end = -1)
# 
# # qc model 
# stopifnot(Matrix::rankMatrix(d1m) == ncol(d1m))
# stopifnot(any(colSums(d1m) == 0) == FALSE)
# 
# # specify formula for lme4 / dream. 
# form <- ~ 0 + timepoint + (1|sampleid)  
# 
# # contrast matrix; modify for dream pipeline (using for consistency model is just ~time)
# contrast_matrix  = makeContrasts(d100_vs_baseline = (dC - d0), levels = colnames(d1m))
# rownames(contrast_matrix) = paste0("timepoint", rownames(contrast_matrix))
# contrast_matrix = as.matrix(contrast_matrix) %>% as.data.frame()
# plotContrasts(contrast_matrix)
# 
# # run pseudobulk modelvoom weights with lme4 model with random intercept for sbject, contrast day 1 fc in h5 vs h1. 
# library(edgeR)
# cl = makeCluster(number_of_cores)
# doParallel::registerDoParallel(cl = cl)
# 
# # Filter matrix get cells passing min exprs threshold
# d1sum_dge = NormalizePseudobulk(aggregate.list =  sl, normalization.method = "RLE", 
#                                 model.metadataframe = d100, comparison.group = "timepoint", 
#                                 minimum.gene.count = 5)
# 
# # for dream lme4 pipeline, strip cell type from colnames of counts. 
# d1sum_dge = lapply(d1sum_dge, function(x){ colnames(x$counts) = str_sub(colnames(x$counts), -6,-1) ; return(x)}) 
# 
# # Get voom observational level weights. 
# v1 = lapply(d1sum_dge, function(x){ x = voom(counts = x, design = d1m, 
#                                              normalize.method = "none", save.plot = T, plot = T)})
# 
# # fit model 
# fit1 = lapply(v1, function(x){ dream(exprObj = x, formula = form, data = d100, L = contrast_matrix) })
# 
# 
# # save outputs. 
# names(d1sum_dge) = names(v1) = names(fit1) = celltypes 
# saveRDS(d1sum_dge, file = paste0(datapath,"/d100sum_dge_h5.rds"))
# model_list = list(v1, fit1)
# names(model_list) = c("voom", "fit")
# saveRDS(model_list, file = paste0(datapath,"/h5d100_embedded_model_list.rds"))
# 
# # Create results objects and run fGSEA on BTM. 
# time_res = GetContrastResultsRaw(limma.fit.object.list = fit1, coefficient.number = 1, contrast.name = "day100",celltypes.vector = celltypes)
# time_rank = GetRankResultsRaw(contrast.result.raw.list = time_res)
# 
# # save readable format 
# res_ = do.call(rbind, time_res)
# write_delim(res_, path = paste0(datapath, "day100_mixedmodel_results.txt"), delim = "\t")
# 
# # save RDS 
# names(time_rank) = names(d1time_res) = celltypes
# saveRDS(time_res, file = paste0(datapath,"/h5_d100timeres.rds"))
# saveRDS(time_rank, file = paste0(datapath,"/h5_d100timerank.rds"))

# run gene set enrichment on gene ranks 
#gsea1 = RunFgseaOnRankList(rank.list.celltype = time_rank, pathways = gp)

## plot result
#gsea = RbindGseaResultList(gsea_result_list = gsea1, NES_filter = -Inf,padj_filter = 0.1)
#GSEABubblePlot(rbind_gsea_result_dataframe = gsea,save_path = figpath,include_negative = TRUE,save_name = "d1vsd0_H5N1cohort.pdf", height = 5)


#########
# unbiased analysis 
# btm = readRDS(file = "signature_curation/BTM_li.rds")
# gsea = RunFgseaOnRankList(rank.list.celltype = time_rank,
#                           pathways = btm,
#                           celltypes.vector = names(time_res),
#                           positive.enrich.only = FALSE)
# gseaf = RbindGseaResultList(gsea_result_list = gsea, NES_filter = -Inf, padj_filter = 0.1)
# GSEABubblePlot(rbind_gsea_result_dataframe = gseaf,save_path = figpath,
#                include_negative = TRUE, save_name = "UNBIASEDd1vsd0_H5N1cohort.pdf", height = 5, width = 7.5)
# 
# lef = GetLeadingEdgeFull(gsea.list = gsea, padj.filter = 0.1,NES.filter = -Inf) %>% bind_rows()
# mtx = GetGeneMatrix(result.list = time_res,
#                     gene_subset = unique(lef$gene), 
#                     stat_for_matrix = "logFC",
#                     pvalfilter = 0.05, 
#                     logfcfilter = 0.25 )
# pheatmap::pheatmap(mtx, fontsize_row = 6, 
#                    color = viridis::inferno(20),treeheight_row = 0, treeheight_col =0,
#                    border_color = NA, width = 4 , height = 8,
#                    filename = paste0(figpath, "leg_unbiased_d100_logfc0.25_p0.05.pdf"))

## Combine Results 
#result_c = CombineResults(gsealist = gsea1, contrastlist = time_res, gseafdr = 0.1,genefdr = 1)
#write_delim(result_c, path = paste0(datapath, "result_gsea_gene_merge.txt"), delim = "\t")

# 
# 
# 
# 
# #####################
# #  Save outputs     #
# #####################
# #write_delim(result, path = paste0(datapath ,"/H5_D100_res.csv"), delim = ",")
# #write_delim(result_c, path = paste0(datapath ,"/H5_d100merged_LE_res.csv"), delim = ",")
# names(time_rank) = names(d1time_res) = celltypes
# saveRDS(time_res, file = paste0(datapath,"/h5_d100timeres.rds"))
# saveRDS(time_rank, file = paste0(datapath,"/h5_d100timerank.rds"))
# #saveRDS(gsea1, file = paste0(datapath,"/h5_d100GP_GSEA.rds"))
