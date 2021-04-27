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

###############################################################
## h1 baseline data reformatting 
s = readRDS(file = "data/CITEseq_raw_PMID32094927_seurat2.4.rds")
high.responders = c("205","207","209","212","215","234","237","245","250","256")
# replace slash strings in variable names 
s@meta.data$celltype_label_3 = gsub(pattern = "/", replacement = "_", s@meta.data$celltype_label_3)
s@meta.data$celltype_label_3 = gsub(pattern = " ", replacement = "", s@meta.data$celltype_label_3)

# add response meta data
md = s@meta.data %>% 
  mutate(adjmfc = if_else(sampleid %in% high.responders, true = "high", false =  "low")) %>% 
  select(barcode_check, adjmfc, celltype_label_3) %>% 
  column_to_rownames("barcode_check")
s = AddMetaData(s,metadata = md)
s = NormalizeData(s, normalization.method = 'LogNormalize')
###############################################################

# run pseudbulk workflow 
meta = s@meta.data
umi = s@data

# remove seurat object 
rm(s); gc()

# QC contingency of cells by subject for each celltype 
tab = scglmmr::SubjectCelltypeTable(metadata = meta, celltype_column = "celltype_label_3", sample_column = "sampleid")

# remove cells prior to pseudobulk analysis 
meta = meta[meta$celltype_label_3 %ni% tab$celltypes_remove, ]

# subset data 
umi = umi[ ,rownames(meta)]

# pseudobulk workflow 
pb = scglmmr::PseudobulkList(rawcounts = umi, metadata = met, sample_col = "sampleid", celltype_col = "celltype_label_3", avg_or_sum = "sum")
av = scglmmr::PseudobulkList(rawcounts = umi, metadata = met, sample_col = "sampleid", celltype_col = "celltype_label_3", avg_or_sum = "average")
designmat = scglmmr::BulkDesignMatrix(metadata = meta, sample_column = "sampleid",variable_column = "adjmfc", pseudobulklist = pb)
dge = scglmmr::NormalizePseudobulk(pseudobulklist = pb, design_matrix = designmat, minimum.gene.count = 5)

# baseline_groups - baseline difference between groups 
c_mat = makeContrasts(baseline_groups = (adjmfchigh - adjmfclow),levels = colnames(designmat))

# fit simple linear model for the baseline group level contrast 
fit = scglmmr::RunVoomLimma(dgelists = dge, design_matrix = designmat, do_contrast_fit = T, my_contrast_matrix = c_mat)

# save 
saveRDS(fit, file = paste0(datapath, "fit_h1_object.rds"))
saveRDS(av, file = paste0(datapath, "av_h1_object.rds"))
saveRDS(pb, file = paste0(datapath, "pb_h1_object.rds"))
saveRDS(meta, file = paste0(datapath, "meta_h1_object.rds"))




