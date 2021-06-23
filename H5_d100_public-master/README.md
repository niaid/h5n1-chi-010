
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CHI H5 day 100 CITE-seq analysis

mulemp at nih.gov permanent email mattmule <AT>
gmail.com

<!-- badges: start -->

<!-- badges: end -->

### Microarray differential expression with of persistence patterns day 100 vs baseline random intercept model

Analysis of day 100 fold change 1) Analysis of gene patterns 4, 6:8,
10:12 in bulk data - subject random intercept day 100 vs baseline
expression. Enrichment of Li BTM genes within adjuvant day 100
associated persistence genes using hypergeometric test. Genes with
adjusted p value \< 0.01 from day 100 vs baseline mixed model were input
into hypergeometric test.

``` r
## array lme4 variance partition workflow 
# format data 
source("microarray_analysis_d100/1_save_eset_object_as_dataframe.r")
# variance partition and mixed GLM workflow 
source("microarray_analysis_d100/2_V4array_d100_lme4_PERSISTENCE_V4.R")
# separate models for adjuvant and non adjuvant group. 
source("microarray_analysis_d100/2_array_adjonly_d100_lme4_PERSISTENCE_V3.R")
source("microarray_analysis_d100/2_nonadj_array_d100_lme4_PERSISTENCE_V3.R")
# save results table combined
source("microarray_analysis_d100/3_make_combined_table.r")
```

### CITE-seq day 100 normalization and protein based clustering

Outliers removed based on MT gene proportion in the greatest 5th
percentile. The top and bottom 3.5 median absolute deviation log library
size cells were removed. mRNA normalization log count per 10,000 library
size scaling factors with Seurat::NormalizeData().

Protein normalization with dsb
<https://cran.r-project.org/web/packages/dsb/index.html>. Options:
use.isotype.control = TRUE, denoise.counts = TRUE. Cells clustered based
on dsb normalized protein level in each cell using a euclidean distance
matrix. Clusters are then refined and annotated based on canonical
protein expression. This section uses R 3.5.

``` r

# normalize single cells from baseline and day 100, cluster annotate and run UMAP 
source("normalization/1_dsb_prot_lognorm_rna.r")
source("2_clustering/1_protein_distance_clustering.r")
source("2_clustering/2_annotation.r") 
source("2_clustering/3_umap.r")
```

### CITE-seq analysis differential expression and enrichment of persistence genes within cell types

Run pseudobulk ramdon intercept model within each protein based cluster.
Subset genes from microarray data differentially expressed on day 100
with log FC \> 0 and adjusted p value \< 0.01 –\> define these as
AS03\_GE\_pattern\_x. Test gene set enrichment within these gene
patterns based on t-statistic within each cluster of day 100 vs baseline
expression. The leading edge genes from these enrichments are then used
in the next step in corresponding cell types to test the average module
z score of baseline high vs low responders from the H1N1 vaccination
group.

``` r

# pseudobulk workflow 
source("3_pseudobulk_de_workflow/1_merge_day0_day100_data_seurat_v2.4.r")
source("3_pseudobulk_de_workflow/2_h5_d0_vd_d100_pseudobulk_de.r")
source("3_pseudobulk_de_workflow/3_figure_generation_h5d100cite_scglmmrv2.r")
```

### Test leading edge genes from single cell day 100 H5N1 vacinees in H1N1 high vs low responders.

Leading edge genes above within each module, cluster and the combined
leading edge genes from each AS03\_gp tested for relative expression
differences in high vs low responders baseline
expression.

``` r
source("3_pseudobulk_de_workflow/4_h1_highresponder_baseline_persistencesig_test.r")
source("3_pseudobulk_de_workflow/5_h1_highresponder_baselinepersistence_figure_generation.r")
```

### Analysis of persistence signals in H1N1 vaccine day 70 vs baseline

Test the persistence gene expression patterns with persistent signal on
day 100 from the H5N1 vaccination cohort (GP 4, 6-8, 10-12) in the H1N1
cohort on day 70 vs baseline. Test pattern enrichment on day 70 with
gene set enrichment analysis and a Fisher’s exact test on the sign of
persistence genes on day 100 from the H5N1 cohort vs the sign (+ or -
fold change) of those genes on day 70 in the H1N1 cohort.

``` r
source("d70_H1/day70_h1_persistencegenes.r")
```
