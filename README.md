# H5N1 Data Analysis Workflow

*Note: For some of the scripts the output folders need to be created manually, therefore, read the output data text to create appropriate directory structure. TODO: probably add a feature in the scripts to check and create directories before writing files. Or write a master shell script that creates the appropriate directories for the entire project.*

# Table of Content

1. [Environment Setup](https://github.niaid.nih.gov/chi/h5n1#environment-setup)
2. [Titers](https://github.niaid.nih.gov/chi/h5n1#titers)
   1. [Figure 1A](https://github.niaid.nih.gov/chi/h5n1#figure-1a)
   2. [Figure 1B](https://github.niaid.nih.gov/chi/h5n1#figure-1a)
   3. [Supplemental Figure 1](https://github.niaid.nih.gov/chi/h5n1#supplemental-figure-1)
   4. [Figure 1C,D,E](https://github.niaid.nih.gov/chi/h5n1#figures-1-cde-figuresprofiles)
3. [Pattern profiles of clinical CBC data and Luminex](https://github.niaid.nih.gov/chi/h5n1#pattern-profiles-of-clinical-cbc-data-and-luminex)
4. [Pattern simulation (Figure 2A)](https://github.niaid.nih.gov/chi/h5n1#pattern-simulation-figure-2A)
5. [Gene Expression PBMC data processing](https://github.niaid.nih.gov/chi/h5n1#gene-expression-pbmc-data-processing)
6. [Pattern discovery in post-vaccination profiles of gene expression](https://github.niaid.nih.gov/chi/h5n1#pattern-discovery-in-post-vaccination-profiles-of-gene-expression)
   1. [Figure 2B](https://github.niaid.nih.gov/chi/h5n1#figure-2b---patterns-profile-plot)
   2. [Supplemental Figure 2A](https://github.niaid.nih.gov/chi/h5n1#supplemental-figure-2a)
   3. [Figure 2C](https://github.niaid.nih.gov/chi/h5n1#btm-enrichment-in-patterns-genes-figure-2c)
7. [Pattern discovery in post-vaccination profiles of flow cytometry data](https://github.niaid.nih.gov/chi/h5n1#pattern-discovery-in-post-vaccination-profiles-of-flow-cytometry-data)
   1. [Figure 2D](https://github.niaid.nih.gov/chi/h5n1#figure-2d)
   2. [Figure 2E](https://github.niaid.nih.gov/chi/h5n1#figure-2e)
   3. [Figure 2F](https://github.niaid.nih.gov/chi/h5n1#figure-2f)
   4. [Supplemental Figure 3](https://github.niaid.nih.gov/chi/h5n1#supplemental-figure-3)
8. [Find signature for adjuvant status prediction](https://github.niaid.nih.gov/chi/h5n1#find-signature-for-adjuvant-status-prediction)
   1. [Figure 3C](https://github.niaid.nih.gov/chi/h5n1#figure-3c)
   2. [Figure 3D](https://github.niaid.nih.gov/chi/h5n1#figure-3d)
   3. [Supplemental Figure 4B](https://github.niaid.nih.gov/chi/h5n1#supplemental-figure-4b)
   4. [Supplemental Figure 4C](https://github.niaid.nih.gov/chi/h5n1#supplemental-figure-4c)
   5. [Figure 3E](https://github.niaid.nih.gov/chi/h5n1#figure-3e)
   6. [Figure 3B](https://github.niaid.nih.gov/chi/h5n1#figure-3b)
   7. [Figure 3F&G](https://github.niaid.nih.gov/chi/h5n1#figure-3f-and-3g)
9.  [Gene Expression PAXgene data processing](https://github.niaid.nih.gov/chi/h5n1#gene-expression-paxgene-data-processing)
10. [Baseline Data Analysis](https://github.niaid.nih.gov/chi/h5n1#baseline-data-analysis)
    1.  [Figure 4A](https://github.niaid.nih.gov/chi/h5n1#figure-4a-combining-btm-enrichment-results-from-pbmc-and-whole-blood-samples)
    2.  [Figure 4B](https://github.niaid.nih.gov/chi/h5n1#figure-4b)
    3.  [Figure 4C](https://github.niaid.nih.gov/chi/h5n1#figure-4c)
11. [Unbliding Results](https://github.niaid.nih.gov/chi/h5n1#unblinding-results-figure-5-suppl-figure-5)
    1.  [Figure 5A](https://github.niaid.nih.gov/chi/h5n1#figure-5a)
    2.  [Figure 5B](https://github.niaid.nih.gov/chi/h5n1#figure-5b)
    3.  [Figure 5C](https://github.niaid.nih.gov/chi/h5n1#figure-5c-figuresifn_signature)
12. [Emory data analysis blindly predicting adjuvant status](https://github.niaid.nih.gov/chi/h5n1#emory-data-analysis-blindly-predicting-adjuvant-status)
    1.  [Supplemental Figure 5D](https://github.niaid.nih.gov/chi/h5n1#supplemental-figure-5d)
13. [Tfh cells data analysis](https://github.niaid.nih.gov/chi/h5n1#tfh-cells-data-analysis)
    1.  [Figure 6A](https://github.niaid.nih.gov/chi/h5n1#figure-6a)
    2.  [Figure 6B](https://github.niaid.nih.gov/chi/h5n1#figure-6b)
    3.  [Figure 6C](https://github.niaid.nih.gov/chi/h5n1#figure-6c)
    4.  [Figure 6D](https://github.niaid.nih.gov/chi/h5n1#figure-6d)
    5.  [Supplemental Figure 7D](https://github.niaid.nih.gov/chi/h5n1#supplemental-figure-7d)
14. [SOMAscan Data Analysis](https://github.niaid.nih.gov/chi/h5n1#somascan-data-analysis)
    1.  [Figure 7A](https://github.niaid.nih.gov/chi/h5n1#figure-7a)
    2.  [Figure 7B](https://github.niaid.nih.gov/chi/h5n1#figure-7b)


# Environment Setup

## Build and Activate Conda Environment

**Note: use env.yaml on HPC and env_rstudio.yaml (installs rstudio also) on Laptop to create Conda environment(s) for this project.**

```
conda env create --file=env.yaml
conda activate h5n1

# Install eNetXlporer through install.packages()
install.packages("eNetXplorer")
```

**Invoke R and initialize the environment**
```R
# Enter the project directory.
R
source("SCRIPTS/0_initialize.r") # You might want to make changes to this script
                                 # if you are running it inside RStudio.
```

## Original Data Files (Yuri's)

`/hpcdata/sg/sg_data/CHI/PROJECTS/H5N1/PAPER`

## Yuri's Singularity container
```
$ module load Singularity
$ cd /hpcdata/sg/sg_data/singularity/test
$ singularity shell -B \
    /hpcdata/sg/sg_data/CHI/PROJECTS/H5N1/PAPER/:/var/workflow1 \
    h5n1_image_180410.img
```
```
singularity shell -B /hpcdata/sg/sg_data/users/farmerr2/sandbox/projects/h5n1:/var/workflow1 h5n1_image_180410.img
```
# Titers

## Figure 1A
URL to the Shinyapp that produces figure 1A `FIGURES/Figure1_url_180322.txt`

## Figure 1B
Purpose: To plot titer respose against time (Day 0, Day 21, Day 28, Day 100)
```R
source("SCRIPTS/titers/mn_titer_profiles.r")
```
Output data:
* `FIGURES/titers/MN_titer_profiles_all_subjects.png`

## Supplemental Figure 1
Fig 1A. Purpose: To plot titer response rate against time.
Input data:
* `DATA_ORIGINAL/Clinical/clinical_info_adj.txt` file containing clinical information example subject id and adjuvant status.
* `DATA_ORIGINAL/HAI/H5N1_serology.txt` file containing serology data.
```R
source("SCRIPTS/titers/titer_response_rate.r")
```
Fig 1B. Purpose: To plot HAI titer profiles against time.
Input data:
* `DATA_ORIGINAL/HAI/H5N1_serology.txt` file containing serology data.
```R
source("SCRIPTS/titers/hai_titer_profiles.r")
```
Fig 1C. Purpose: To plot titer at peak.
Input data:
* `DATA_ORIGINAL/Titres/titre.txt` file containing titer data.
* `DATA_ORIGINAL/Clinical/clinical_info_adj.txt` file containing clinical information example subject id and adjuvant status.
```R
source("SCRIPTS/titers/mn_titer_peak.r")
```
Output data:
* `RESULTS/titers/MN_titer_peak_time.txt` TSV file with information about peak time, titer value at peak, titer value after peak, and decline per subject. 

# Pattern profiles of clinical CBC data and Luminex
## Figures 1 C,D,E (FIGURES/profiles/)
Fig 1C, D. Purpose: To plot total monocyte vs time and neutorphils vs time.
Input data:
* `DATA_ORIGINAL/Clinical/H5N1_BTRIS.txt` file containing BTRIS data.
```R
source("SCRIPTS/profiles/Monocytes_figure.r")
source("SCRIPTS/profiles/Neutrophils_figure.r")
```
Fig 1E: Purpose: To plot IP-10 data.
Input data:
* `DATA_PROCESSED/Luminex/luminex_data.rds` R data structure file with the luminex processed data.
```R
source("SCRIPTS/profiles/IP10_figure.r") # Copy Luminex data to DATA_PROCESSED
```

# Pattern simulation (Figure 2A)
```R
source("SCRIPTS/pattern_sim/pattern_simulation.r") # TODO: This script throws error if run as whole, but runs fine in interactive mode. Results seems to be stochastic and are not identical to what we have in the paper.
```

# Gene Expression PBMC data processing

*Note: I did not repeat these steps as they seem to compute intensive and also stochastic in nature.*

First, the CEL files were processed with  Power Tools RNA-sketch algorithm. 
```R
source("SCRIPTS/MA/processing_pbmc/apt.config.r")
source("SCRIPTS/MA/processing_pbmc/apt.call.r")
```

Quantile normalization implement in RMA-sketch algorithm uses random subset of data to decrease memory usage. Because of this at different runs the algorithm may generates slightly different results. For reproducibility we recommend using precomputed APT output located in DATA_PROCESSED/Microarrays/PBMC/apt_summary.txt. *(this seem to be wrong because apt_summary.txt is produced by the scripts in the section immediately below. On the other hand rma-sketch.summary.txt file seem to be the optout of this procedure and input to ones below. I copied rma-sketch.summary.txt file over.)*

## Create ExpressionSet from APT output

*Note: copy rma-sketch.summary.txt file to DATA_PROCESSED/Microarrays/PBMC directory if you are skipping the steps in the section above.*

Input data:  
* `DATA_ORIGINAL/Microarrays/PBMC/sample_info.txt` from the original data files. 
* `DATA_PROCESSED/Microarrays/PBMC/rma-sketch.summary.txt` copied Yuri's pre-computed for reproducibility. This file was supposed to be generated in the previous step.
  
```R
source("SCRIPTS/MA/processing_pbmc/eset.config.r")
source("SCRIPTS/MA/processing_pbmc/eset.call.r")
```
Output data:  
* `DATA_PROCESSED/Microarrays/PBMC/apt.summary.txt`
* `DATA_PROCESSED/Microarrays/PBMC/eset.apt.RData`

## Probesets to genes mapping

### Convert the table to probeset-gene mapping  
Input data:  
* `DATA_ORIGINAL/Microarrays/annotation/HuGene-2_1-st-v1.na35.hg19.transcript.csv` file downloaded from Affymetrix web site
```R
source("SCRIPTS/MA/annotation/affy_hugene-2_1-st_annotation.r")
```
Output data:  
* `DATA_PROCESSED/Microarrays/annotation/affy_hugene-2_1-st_unique_all.txt`
* `DATA_PROCESSED/Microarrays/annotation/affy_hugene-2_1-st_unique.txt`

### Select the best probeset for a gene
Input data:  
* `DATA_PROCESSED/Microarrays/annotation/affy_hugene-2_1-st_unique.txt`
* `DATA_PROCESSED/Microarrays/PBMC/eset.apt.RData`
```R
source("SCRIPTS/MA/annotation/generate_ps2gene_map.config.r")
source("SCRIPTS/MA/annotation/generate_ps2gene_map.call.r") # It calls a function from this script SCRIPTS/functions/pick.probeset.r TODO: I don't understand if it is producing any file
```

## Data post processing

### We found that two samples were switched. This is to correct it.

Input data:  
* `SCRIPTS/MA/filtering_pbmc/samples.switch.tx`
```R
source("SCRIPTS/MA/filtering_pbmc/switch.samples/switch.samples.call.r") # Fetches a function from this script SCRIPTS/MA/filtering_pbmc/switch.samples
```
Output data:
* `DATA_PROCESSED/Microarrays/PBMC/eset.corrected.RData`

### Apply different filtering to samples and genes. The probesets mapped to genes.

*Note: these scripts source multiple scripts within their respective folders. TODO: it looks messy and hard to figure out what are the inputs and outputs.*
```R
source("SCRIPTS/MA/filtering_pbmc/samples.clean_genes.all/filtering.r")
source("SCRIPTS/MA/filtering_pbmc/samples.clean_genes.iqr/filtering.r")
source("SCRIPTS/MA/filtering_pbmc/samples.all_genes.all/filtering.r")
source("SCRIPTS/MA/filtering_pbmc/samples.all_genes.iqr/filtering.r")
```

### Calculate fold change from day 0

Input data:  
* `DATA_PROCESSED/Microarrays/PBMC/samples.all_genes.all/eset.genes.filtered.RData` biobase object file containing expression data.
```R
source("SCRIPTS/MA/calculate_d0_fc/calculate_d0_fc_pbmc.r")
```
Output data:  
* `DATA_PROCESSED/Microarrays/PBMC/samples.clean_genes.iqr/gexp_d0_fc.RData` matrix with genese on rows, and subject and time points on the columns. 

# Pattern discovery in post-vaccination profiles of gene expression

## Profiles clustering with DIANA

*Note: before running this script create the relevant folders by `mkdir RESULTS/Microarrays/PBMC/pattern_discovery/`*

Input data:  
* `DATA_PROCESSED/Microarrays/PBMC/samples.clean_genes.iqr/gexp_d0_fc.RData`
```R
source("SCRIPTS/MA/pattern_discovery/pattern_discovery.r")
```
Output data:  
* `RESULTS/Microarrays/PBMC/pattern_discovery/df.mat.rds`
* `RESULTS/Microarrays/PBMC/pattern_discovery/df.mat.cc.rds` correlation matrix.
* `RESULTS/Microarrays/PBMC/pattern_discovery/diana.object.abs.rds` DIANA output.
* `RESULTS/Microarrays/PBMC/pattern_discovery/diana_dendrogram.png` DIANA dendrogram image.
  
## Cut the dendrogram tree at different levels and detect stable clusters

Input data:  
* `RESULTS/Microarrays/PBMC/pattern_discovery/diana.object.abs.rds`
* `RESULTS/Microarrays/PBMC/pattern_discovery/df.mat.rds`
* `RESULTS/Microarrays/PBMC/pattern_discovery/df.mat.cc.rds` correlation matrix.
```R
source("SCRIPTS/MA/pattern_discovery/patterns_cutTree_stable.r")
```
Output data:  
* `RESULTS/Microarrays/PBMC/pattern_discovery/diana_data.rds`
* `RESULTS/Microarrays/PBMC/pattern_discovery/z.rds`
* `RESULTS/Microarrays/PBMC/pattern_discovery/zl.rds`
* `RESULTS/Microarrays/PBMC/pattern_discovery/GE_patterns.txt`
  
## Filter the patterns
Input data:  
* `RESULTS/Microarrays/PBMC/pattern_discovery/df.mat.cc.rds`
* `RESULTS/Microarrays/PBMC/pattern_discovery/GE_patterns.txt`
```R
source("SCRIPTS/MA/pattern_discovery/pattern_filter.r")
```
Output data:  
* `RESULTS/Microarrays/PBMC/pattern_discovery/GE_patterns_filt.txt` TSV file with the first column containing the subject Id and gene names and the second column containing the pattern number to which they belong.
  
## Summarize patterns stats
Calculate z-score etc. 
Input data:
* `RESULTS/Microarrays/PBMC/pattern_discovery/df.mat.rds`
* `RESULTS/Microarrays/PBMC/pattern_discovery/GE_patterns_filt.txt`
```R
source("SCRIPTS/MA/pattern_discovery/patterns_stats.r")
```
Output data:
* `RESULTS/Microarrays/PBMC/pattern_discovery/df.cl.stat.rds` data frame containing pattern statistics. 
* `RESULTS/Microarrays/PBMC/pattern_discovery/df.patt.rds` 14 X 14 matrix for 14 patterns and 14 time points.
  
## Figure 2B - patterns profile plot
Input data:
* `RESULTS/Microarrays/PBMC/pattern_discovery/GE_patterns_filt.txt`
* `RESULTS/Microarrays/PBMC/pattern_discovery/df.cl.stat.rds`
```R
source("SCRIPTS/MA/pattern_discovery/plot_patterns.r")
```
Output data:
* `FIGURES/GE_patterns_profiles.png`
* `FIGURES/GE_patterns_profiles_2col.png`
* `FIGURES/GE_patterns_profiles_horiz.png`
  
## Expand list of pattern signature genes
Input data: 
* `DATA_PROCESSED/Microarrays/PBMC/samples.clean_genes.iqr/gexp_d0_fc.RData`
```R
source("SCRIPTS/MA/pattern_discovery/pattern_filter_expanded.r")
```
Output data:
* `RESULTS/Microarrays/PBMC/pattern_discovery/df_6672g.rds` data frame with expression data after relaxing the selection criteria.
* `RESULTS/Microarrays/PBMC/pattern_discovery/df.mat.2_6672g.rds` data frame created above converted to matrix with some modifications. 

## Compute correlations and clean up the genes
Input data for the first script:
* `RESULTS/Microarrays/PBMC/pattern_discovery/df_6672g.rds`
* `RESULTS/Microarrays/PBMC/pattern_discovery/df.mat.2_6672g.rds`
* `RESULTS/Microarrays/PBMC/pattern_discovery/GE_patterns_filt.txt`
* `RESULTS/Microarrays/PBMC/pattern_discovery/df.patt.rds`

Input data for the second script:
* `RESULTS/Microarrays/PBMC/pattern_discovery/df_6672g.rds`
* `RESULTS/Microarrays/PBMC/pattern_discovery/df.cor.p_6672g.rds`

### Supplemental Figure 2A
Purpose: to add more genes to the patterns by relaxing the selection criteria.
```R
source("SCRIPTS/MA/pattern_discovery/pattern_expanded_genes_cor.r")
source("SCRIPTS/MA/pattern_discovery/pattern_expanded_genes_clean.r") # Supplemental Figure 2A
```
Output data from the first script:
* `RESULTS/Microarrays/PBMC/pattern_discovery/df.cor_6672g.rds",n_genes`
* `RESULTS/Microarrays/PBMC/pattern_discovery/df.cor.p_6672g.rds`

Output data from the second script:
* `RESULTS/Microarrays/PBMC/pattern_discovery/patterns_correlations__q%.2f_%dg.txt`
* `RESULTS/Microarrays/PBMC/pattern_discovery/patterns_correlations_p__q%.2f_%dg.txt`
* `FIGURES/pattern_correlations/pattern_%d.%d_cor_km2_heatmap.png`
  
## Compute subject scores for each pattern
Input data:  
* `RESULTS/Microarrays/PBMC/pattern_discovery/GE_patterns_filt.txt`
* `RESULTS/Microarrays/PBMC/pattern_discovery/df_6672g.rds`
* `RESULTS/Microarrays/PBMC/pattern_discovery/df.mat.2_6672g.rds`
* `RESULTS/Microarrays/PBMC/pattern_discovery/df.patt.rds`
* `RESULTS/Microarrays/PBMC/pattern_discovery/patt.genes.clean.rds` dictionary (named list) with 14 pattern names as keys and list of genes per pattern as values. 
```R
source("SCRIPTS/MA/pattern_discovery/patterns_to_subjects.r")
# Sources another script. 
```
Output data:
* `RESULTS/Microarrays/PBMC/pattern_discovery/subj.patt.cor.rds`
* `RESULTS/Microarrays/PBMC/pattern_discovery/subj.patt.cor.txt`

## Output the table of genes (with annotations) for each pattern
Input data:
* `RESULTS/Microarrays/PBMC/pattern_discovery/df_6672g.rds`
* `RESULTS/Microarrays/PBMC/pattern_discovery/df.patt.rds`
* `RESULTS/Microarrays/PBMC/pattern_discovery/patt.genes.clean.rds`
* `RESULTS/Microarrays/PBMC/pattern_discovery/GE_patterns_filt.txt`
* `RESULTS/Microarrays/PBMC/pattern_discovery/df.cor_6672g.rds`
* `RESULTS/Microarrays/PBMC/pattern_discovery/df.cor.p_6672g.rds`
* `RESULTS/Microarrays/PBMC/pattern_discovery/subj.patt.cor.rds`
```R
source("SCRIPTS/MA/pattern_discovery/pattern_genes_output.r")
```
Ouput data:
* `RESULTS/Microarrays/PBMC/pattern_discovery/GE_pattern_genes/GE_pattern_%02d_gene.sig.contr_ann.txt` dataframe with annotated gene list.

## BTM enrichment in patterns genes. Figure 2C
Purpose: To carry out GSEA of genes per pattern using blood transcription modules.
Input data:
* `RESULTS/Microarrays/PBMC/pattern_discovery/GE_pattern_genes/GE_pattern_%02d_gene.sig.contr_ann.txt`
* `RESULTS/Microarrays/PBMC/pattern_discovery/GE_pattern_genes/df_6672g.rds`
```R
source("SCRIPTS/MA/pattern_discovery/pattern_BTM_enrichment.r")
```
Ouput data:
* `FIGURES/GE_pattern_BTM_enrichment_heatmap_{mset}.pdf`
* `FIGURES/E_pattern_BTM_enrichment_heatmap_{mset}.png`

## Add data for subject s10 and update the score matrix
Input data:
* `DATA_PROCESSED/Microarrays/PBMC/samples.all_genes.iqr/eset.genes.filtered.RData`
* `RESULTS/Microarrays/PBMC/pattern_discovery/df.patt.rds`
* `RESULTS/Microarrays/PBMC/pattern_discovery/patt.genes.clean.rds`
```R
source("SCRIPTS/MA/pattern_discovery/s10_peaks_assessment.r")
# Script above sources another script.
source("SCRIPTS/MA/pattern_discovery/pattern_scores_in_samples_GE_incl.s10.r")
```
Output data:
* `RESULTS/Microarrays/PBMC/pattern_discovery/s10_pattern_scores.rds`

# Pattern discovery in post-vaccination profiles of flow cytometry data
FlowJo software was used to export flow data, therefore copy processed data from Yuri's folder for further analysis.

## Generate Trajectory Matrix
TODO: These scripts are very complex to understand just by reading the code. 
```R
# source("SCRIPTS/Flow_10c/Flow_10c_TrajMatrix.R")
source("SCRIPTS/Flow_10c/Flow_10c_TrajMatrix_v2.R")
source("SCRIPTS/Flow_10c/Flow_10c_TrajMatrix_v3.R")
# source("SCRIPTS/Flow_10c/Flow_10c_TrajMatrix_new.R")
```

## Genetrate Trajectory Clusters
```R
# source("SCRIPTS/Flow_10c/Flow_10c_TrajCluster.R")
source("SCRIPTS/Flow_10c/Flow_10c_TrajCluster_v2.R")
source("SCRIPTS/Flow_10c/Flow_10c_TrajCluster_v3.R")
# source("SCRIPTS/Flow_10c/Flow_10c_TrajCluster_new.R")
```

## Figure 2D
Input data:
* `RESULTS/Flow_10c/DeltaThr0.5_DonorThr0.3_PercentOfParent_CORpearson_cutreeHybrid_deepSplit0_minClusterSize31_Modules.txt`
```R
source("SCRIPTS/Flow/pattern_figures/plot_flow_patters_only.r")
```
Output data:
* `FIGURES/Flow_patterns_profiles.png`
* `FIGURES/Flow_patterns_profiles_2col.png`
* `FIGURES/Flow_patterns_profiles_horiz.png`

## Figure 2E
Purpose: To plot a heatmap with the type of cell population found per pattern in the flow data.
Input data:
* `RESULTS/Flow_10c/DeltaThr0.5_DonorThr0.3_PercentOfParent_CORpearson_cutreeHybrid_deepSplit0_minClusterSize31_CPs.txt`
* `DATA_ORIGINAL/Flow_10c/Flow_10c_ann.txt`
```R
source("SCRIPTS/Flow/pattern_figures/pattern_flow_ann_heatmap.r")
```

## Figure 2F
Input data:
* `RESULTS/Microarrays/PBMC/pattern_discovery/subj.patt.cor.incl.s10.rds` pattern correlation including subject 10. 
* `RESULTS/Flow_10c/DeltaThr0.5_DonorThr0.3_PercentOfParent_CORpearson_cutreeHybrid_deepSplit0_minClusterSize31_DonorScores.txt`
* `DATA_ORIGINAL/Titres/titre.txt`
* `DATA_ORIGINAL/Clinical/clinical_info_adj.txt`
```R
source("SCRIPTS/MA/pattern_discovery/pattern_scores_GE_flow_heatmap.r")
```
Output data:
* `FIGURES", "GE_subject_patterns_cor_flow_titers_sex.pred_incl.s10.pdf[.png]`

## Supplemental Figure 3
Purpose: Quality control of flow cytometry-based response patterns. (TODO: very long script hard to understand.)
```R
source("SCRIPTS/Flow_10c/Flow_10c_TrajCluster_QM.R") 
```


# Find signature for adjuvant status prediction

## Figure 3C
Purpose: Predicting Adjuvant signature using cellular and transcriptomic parameters. (C) Principal component analysis of patients using top selected patterns - Gp01, Gp02, Gp03, Fp01, Fp05. Color indicates two clusters of patients after k-mean clustering with k=2.  
Input data:
* `RESULTS/Flow_10c/DeltaThr0.5_DonorThr0.3_PercentOfParent_CORpearson_cutreeHybrid_deepSplit0_minClusterSize31_DonorScores.txt`
* `RESULTS/Microarrays/PBMC/pattern_discovery/subj.patt.cor.incl.s10.rds`
```R
source("SCRIPTS/adjuvant_prediction/2peaks_pca_2clusters.r")
```
Output data:
* `FIGURES/2peak_patterns_PCA.png`
* `RESULTS/Adjuvant_prediction/2peaks_patterns_PCA.txt`

## Figure 3D
Purpose: (D) Testing difference in IP-10 change between subject in two clusters revealed by k-mean clustering of selected pattern scores. 
Input data:
* `DATA_PROCESSED/Luminex/luminex_data.rds`
* `RESULTS/Adjuvant_prediction/2peaks_patterns_PCA.txt`
```
source("SCRIPTS/adjuvant_prediction/ip10_2clusters_compare.r")
```
Output data:
* `FIGURES/Luminex/2peak_clusters_compare_{an}.png`

## Supplemental Figure 4B
Input data:
* `DATA_PROCESSED/Luminex/luminex_data.rds`
* `RESULTS/Adjuvant_prediction/2peaks_patterns_PCA.txt`
```R
source("SCRIPTS/adjuvant_prediction/cytokines_2clusters_compare.r")
```
Output data:
* `FIGURES/Luminex/2peak_clusters_compare_{an}.png`

## Supplemental Figure 4C
Purpose: (C) Correlation between eigengene scores from two modules associated with type I interferon and antiviral response discovered at day 0 in PBMC (Gb13) and whole blood (GbWB11). A/Indonesia titers values at day 42 are represented by dot size, and subjects separated by gender with colors: male (blue) and female (red).  
Input data:
* `RESULTS/Microarrays/PBMC/baseline/GE.pbmc_d0_WGCNA_ME_scores.txt`
* `RESULTS/Microarrays/PAXgene/baseline/GE.pax_d0_WGCNA_ME_scores.txt`
* `DATA_ORIGINAL/Clinical/clinical_info.txt`
* `DATA_ORIGINAL/Titres/titre.txt`
```R
source("SCRIPTS/MA/baseline/Gb13_vs_GbWB11.r")
```
Output data:
* `FIGURES/Baseline/Gb13_vs_GbWB11.png`

## Figure 3E
Input data for the first script:
*  `DATA_PROCESSED/Luminex/luminex_data.rds`

```R
source("SCRIPTS/adjuvant_prediction/ip10_2peak_scores.r")
source("SCRIPTS/adjuvant_prediction/2peaks_pca_final_heamap.r") # TODO: long and complex to understand. 
```
Output data for the first script:
* `RESULTS/Luminex/ip10_2peak_score.rds`
* `RESULTS/Luminex/ip10_2peak_score.txt`
## Elastic net models

### Generate input data:
```R
source("SCRIPTS/eNetXplorer/eNet_input_r1.r")
source("SCRIPTS/eNetXplorer/eNet_input_r2.r")
source("SCRIPTS/eNetXplorer/eNet_input_r3.r")
```
TODO: the sripts above require `RESULTS/Adjuvant_prediction/adjuvant_predicted_subjects.txt` which none of the previous scripts seems to produce.

### Run eNetXplorer:
```R
source("SCRIPTS/eNetXplorer/eNetXplorer_R1_180530.R") # TODO: update conda environment with library missForest
source("SCRIPTS/eNetXplorer/eNetXplorer_R2_180530.R")
source("SCRIPTS/eNetXplorer/eNetXplorer_R3_180530.R")
```

### Figure 3B
```R
source("SCRIPTS/eNet_figures/enet_plots_R1.r")
```

### Figure 3F and 3G
```R
source("SCRIPTS/eNet_figures/enet_plots_R3.r")
source("SCRIPTS/eNet_figures/enet_plots_R2.r")
```

# Gene Expression PAXgene Data Processing

## Data Post Processing

### We found that two samples were switched. This is to correct it.

**TODO: I don't have the scripts in the workflow that generates "DATA_PROCESSED/Microarrays/PAXgene/eset.apt.RData". For now I will just copy it from Yuri's**
```R
source("SCRIPTS/MA/filtering_pax/switch.samples/switch.samples.call.r")
```

### Apply different filtering to samples and genes. The probesets mapped to genes.
```R
source("SCRIPTS/MA/filtering_pax/filtering.r")
```

### Calculate fold change from day 0
```R
source("SCRIPTS/MA/calculate_d0_fc/calculate_d0_fc_pax.r")
```

# Baseline Data Analysis

## Preparing PBMC day 0 samples
```R
source("SCRIPTS/MA/baseline_pbmc/d0_filter.r")
```

## WGCNA clustering of PBMC samples
Input data:
* `DATA_PROCESSED/Microarrays/PBMC/baseline/dat0.in_isv.0.7_8144g.rds`
* `DATA_PROCESSED/Microarrays/PBMC/baseline/info0.in.rds`
```R
source("SCRIPTS/MA/baseline_pbmc/d0_wgcna.r")
source("SCRIPTS/MA/baseline_pbmc/d0_wgcna_output.r")
```
Output data:
* `RESULTS/Microarrays/PBMC/baseline/dat0.in-networkConstruction-auto_power.%d.RData`

## BTM enrichment analysis of data from PBMC samples
```R
source("SCRIPTS/MA/baseline_pbmc/d0_wgcna_BTM_enrichment.r")
```

## Preparing whole blood (PAXgene) day 0 samples
```R
source("SCRIPTS/MA/baseline_pax/d0_filter.r") # File name changed
```

## WGCNA clustering of whole blood samples
Input data:
* `DATA_PROCESSED/Microarrays/PAXgene/baseline/dat0.in_isv.0.7_7901g.rds`
* `DATA_PROCESSED/Microarrays/PAXgene/baseline/info0.in.rds`
```R
source("SCRIPTS/MA/baseline_pax/d0_wgcna.r")
source("SCRIPTS/MA/baseline_pax/d0_wgcna_output.r")
```

## BTM enrichment analysis of data from whole blood samples
```R
source("SCRIPTS/MA/baseline_pax/d0_wgcna_BTM_enrichment.r")
```

## FIgure 4A. Combining BTM enrichment results from PBMC and whole blood samples
Also Suppl. Figure 6A
```R
source("SCRIPTS/MA/baseline/plot_BTM_pbmc_pax.r")
```

## Elastic net models for baseline prediction
Generate input data:
```R
source("SCRIPTS/eNetXplorer/eNet_input_r6.r")
```

Run eNetXplorer:
```R
source("SCRIPTS/eNetXplorer/eNetXplorer_R6_181022.R")
```

## Figure 4B
```R
source("SCRIPTS/eNet_figures/enet_plots_R6.r")
# NOTE: adapted from enet_plots_all.r
```
## Figure 4C
```R
source("SCRIPTS/MA/baseline/IFN.gene_overlap_figure.r")
```

# Unblinding results (Figure 5, suppl. Figure 5)
## Figure 5A
Drawn in Illustrator

## Figure 5B
```R
source("SCRIPTS/adjuvant_prediction/pattern_gene_time_score_sel.subject.r")
source("SCRIPTS/adjuvant_prediction/pattern_flow_time_score_sel.subject.r")
source("SCRIPTS/adjuvant_prediction/IP10_time_score_sel.subject.r")
```

## Figure 5C (FIGURES/IFN_signature).
```R
source("SCRIPTS/MA/baseline/GbWB11.d0_vs_MN.d28.r")
```

# Emory data analysis blindly predicting adjuvant status

## Process Data to Generate Espression Set
Input data:
* `DATA_ORIGINAL/Emory/T H - Vax010_RMA_CHI.txt.zip`
* `DATA_ORIGINAL/Emory/T H - Vax010_RMA_CHI/Vax010_RMA_CHI.txt`
* `DATA_ORIGINAL/Emory/Vax010_demographics_wAge.txt`
```R
source("SCRIPTS/Emory/emory_data.r")
```
Output data:
* `DATA_PROCESSED/Emory/eset.rds`

## Get Annotations
Input data:
* `DATA_ORIGINAL/Emory/GPL13158.annot.gz`
```R
source("SCRIPTS/Emory/get_ann.r")
```
Output data:
* `DATA_PROCESSED/Emory/GPL13158.ann.txt`

## Map Probes to Genes
Input data:
* `DATA_PROCESSED/Emory/eset.rds`
* `DATA_PROCESSED/Emory/GPL13158.ann_PC1.txt`
```R
source("SCRIPTS/Emory/probe2gene.r")
```
Output daa:
* `DATA_PROCESSED/Emory/eset.gene.rds`

## Claculate 2 Peak Scores
```R
source("SCRIPTS/Emory/2peak_scores.r")
```

## Supplemental Figure 2C

## Supplemental Figure 2D

## Supplemental Figure 2E

## Supplemental Figure 5D
```R
source("SCRIPTS/Emory/adjuvant_prediction.r")
```

# Tfh cells data analysis (Figure 6)


## Figure 6A
```R
source("SCRIPTS/Flow/Tfh/Tfh_activated_profiles.r") 
source("SCRIPTS/Flow/Tfh/Tfh_activated_profile.r") 
```
## Figure 6B
```R
source("SCRIPTS/Flow/Tfh/Tfh_activated_correlations.r")
```
## Figure 6C
```R
source("SCRIPTS/Flow/Tfh/Tfh_total_profiles_v2.r")
```
## Figure 6D
```R
source("SCRIPTS/Flow/Tfh/Tfh_CXCR3_proportion_bar.r")
```
## Supplemental Figure 7D
```R
source("SCRIPTS/Flow/Tfh/Tfh_total_profiles_supp.r")
```

# SOMAscan Data Analysis

## DATA Preparation
Input data for the first script:
* `DATA_PROCESSED/SOMAscan/Samples.txt`
* `DATA_PROCESSED/SOMAscan/Somamers.txt`
* `DATA_PROCESSED/SOMAscan/Hyb.Cal.MedNorm_RFU.txt`
```R
source("SCRIPTS/SOMAscan/Data_Processing.R")
source("SCRIPTS/SOMAscan/Data_Normalization.R")
```
Output data from the first script:
* `DATA_PROCESSED/SOMAscan/Final_RFU.txt`
* `DATA_PROCESSED/SOMAscan/Final_Samples.txt`
* `DATA_PROCESSED/SOMAscan/Final_Somamers.txt`

## Figure 7A
**Data preparation for eNetXplorer**  
```R
source("SCRIPTS/eNetXplorer/eNet_input_r8.r")
```
**Run eNetXplorer**  
```R
source("SCRIPTS/eNetXplorer/eNetXplorer_R8.r")
```
**Plot eNetXplorer results**  
**NOTE:** worked in interactive mode, but not upon source.
```R
source("SCRIPTS/eNet_figures/enet_plots_R8.r") 
```

## Figure 7B
**NOTE:** Files `/RESULTS/SOMAscan/[pbmc,pax]_ADJ_RSPO3_bySex.txt` are required by `SCRIPTS/SOMAscan/soma_BTM_enrichment.r`. Yuri and I couldn't find any script that generate these files. Yuri got them from Julian.
```R
source("SCRIPTS/SOMAscan/soma_BTM_enrichment.r")
source("SCRIPTS/SOMAscan/soma_enrichment_heatmap.r")
```
 

