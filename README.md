# H5N1 Data Analysis Workflow

*Note: For some of the scripts the output folders need to be created manually, therefore, read the output data text to create appropriate directory structure. TODO: probably add a feature in the scripts to check and create directories before writing files. Or write a master shell script that creates the appropriate directories for the entire project.*

# Table of Content

1. [Environment Setup](https://github.niaid.nih.gov/chi/h5n1#environment-setup)
2. [Titers](https://github.niaid.nih.gov/chi/h5n1#titers)
3. [Pattern profiles of clinical CBC data and Luminex](https://github.niaid.nih.gov/chi/h5n1#pattern-profiles-of-clinical-cbc-data-and-luminex)
4. [Pattern simulation (Figure 2A)](https://github.niaid.nih.gov/chi/h5n1#pattern-simulation-figure-2A)
5. [Gene Expression PBMC data processing](https://github.niaid.nih.gov/chi/h5n1#gene-expression-pbmc-data-processing)
6. [Pattern discovery in post-vaccination profiles of gene expression](https://github.niaid.nih.gov/chi/h5n1#pattern-discovery-in-post-vaccination-profiles-of-gene-expression)
7. [Pattern discovery in post-vaccination profiles of flow cytometry data](https://github.niaid.nih.gov/chi/h5n1#pattern-discovery-in-post-vaccination-profiles-of-flow-cytometry-data)
8. [Find signature for adjuvant status prediction](https://github.niaid.nih.gov/chi/h5n1#find-signature-for-adjuvant-status-prediction)
9.  [Gene Expression PAXgene data processing](https://github.niaid.nih.gov/chi/h5n1#gene-expression-paxgene-data-processing)
10. [Baseline Data Analysis](https://github.niaid.nih.gov/chi/h5n1#baseline-data-analysis)
11. [Emory data analysis blindly predicting adjuvant status](https://github.niaid.nih.gov/chi/h5n1#emory-data-analysis-blindly-predicting-adjuvant-status)
12. [Tfh cells data analysis](https://github.niaid.nih.gov/chi/h5n1#tfh-cells-data-analysis)
13. [SOMAscan data analysis](https://github.niaid.nih.gov/chi/h5n1#somascan-data-analysis)


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
URL to the Shinyapp that produces figure 1A `H5N1/MANUSCRIPT/Figures/Figure1_url_180322.txt`

## Figure 1B
```R
source("SCRIPTS/titers/mn_titer_profiles.r")
```
Output data:
* `FIGURES/titers/MN_titer_profiles_all_subjects.pdf`

## Supplemental Figure 1(FIGURES/titers):
```R
source("SCRIPTS/titers/titer_response_rate.r")
source("SCRIPTS/titers/hai_titer_profiles.r")
source("SCRIPTS/titers/mn_titer_peak.r")
```

# Pattern profiles of clinical CBC data and Luminex
## Figures 1 C,D,E (FIGURES/profiles/)
```R
source("SCRIPTS/profiles/Monocytes_figure.r")
source("SCRIPTS/profiles/Neutrophils_figure.r")
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
* `DATA_PROCESSED/Microarrays/PBMC/samples.clean_genes.iqr`
* `DATA_PROCESSED/Microarrays/PBMC/samples.all_genes.all/eset.genes.filtered.RData`
```R
source("SCRIPTS/MA/calculate_d0_fc/calculate_d0_fc_pbmc.r") # Sourcing this file SCRIPTS/functions/factor.date.r
```
Output data:  
* `DATA_PROCESSED/Microarrays/PBMC/samples.clean_genes.iqr/gexp_d0_fc.RData`

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
* `RESULTS/Microarrays/PBMC/pattern_discovery/df.mat.cc.rds`
* `RESULTS/Microarrays/PBMC/pattern_discovery/diana.object.abs.rds`
* `RESULTS/Microarrays/PBMC/pattern_discovery/diana_dendrogram.png`
  
## Cut the dendrogram tree at different levels and detect stable clusters

Input data:  
* `RESULTS/Microarrays/PBMC/pattern_discovery/diana.object.abs.rds`
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
* `/RESULTS/Microarrays/PBMC/pattern_discovery/GE_patterns_filt.txt`
  
## Summarize patterns stats
Input data:
* `RESULTS/Microarrays/PBMC/pattern_discovery/df.mat.rds`
* `/RESULTS/Microarrays/PBMC/pattern_discovery/GE_patterns_filt.txt`
```R
source("SCRIPTS/MA/pattern_discovery/patterns_stats.r")
```
Output data:
* `/RESULTS/Microarrays/PBMC/pattern_discovery/df.cl.stat.rds`
* `/RESULTS/Microarrays/PBMC/pattern_discovery/df.patt.rds`
  
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
* `RESULTS/Microarrays/PBMC/pattern_discovery/df_{sum(gi)}g.rds`
* `RESULTS/Microarrays/PBMC/pattern_discovery/df.mat.2_{sum(gi)}g.rds`

## Compute correlations and clean up the genes
Input data for the first script:
* `RESULTS/Microarrays/PBMC/pattern_discovery/df_6672g.rds`
* `RESULTS/Microarrays/PBMC/pattern_discovery/df.mat.2_6672g.rds`
* `RESULTS/Microarrays/PBMC/pattern_discovery/GE_patterns_filt.txt`
* `RESULTS/Microarrays/PBMC/pattern_discovery/df.patt.rds`

Input data for the second script:
* `RESULTS/Microarrays/PBMC/pattern_discovery/df_6672g.rds`
* `RESULTS/Microarrays/PBMC/pattern_discovery/df.cor.p_6672g.rds`

```R
source("SCRIPTS/MA/pattern_discovery/pattern_expanded_genes_cor.r")
source("SCRIPTS/MA/pattern_discovery/pattern_expanded_genes_clean.r")
```
Output data from the first script:
* `RESULTS/Microarrays/PBMC/pattern_discovery/df.cor_%dg.rds",n_genes`
* `RESULTS/Microarrays/PBMC/pattern_discovery/df.cor.p_%dg.rds`

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
* `RESULTS/Microarrays/PBMC/pattern_discovery/patt.genes.clean.rds`
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
```R
source("SCRIPTS/MA/pattern_discovery/pattern_genes_output.r")
```
Ouput data:
* `RESULTS/Microarrays/PBMC/pattern_discovery/GE_pattern_genes/GE_pattern_%02d_gene.sig.contr_ann.txt`

## BTM enrichment in patterns genes. Figure 2C
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
TODO: I don't know which one of these were used. I ran all of them
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

## Figure 2D (use Julian’s data)
```R
source("SCRIPTS/Flow/pattern_figures/plot_flow_patters_only.r")
```

## Figure 2E (use Julian’s data)
```R
source("SCRIPTS/Flow/pattern_figures/pattern_flow_ann_heatmap.r") # TODO: did not work
```

## Add data for subject s10

## Figure 2F (partially use Julian’s data)

## Supplemental Figure 2B

# Find signature for adjuvant status prediction

## Figure 3C
```
source("SCRIPTS/adjuvant_prediction/2peaks_pca_2clusters.r")
```

## Figure 3D
```
source("SCRIPTS/adjuvant_prediction/ip10_2clusters_compare.r")
```

## Supplemental Figure 4B
```
source("SCRIPTS/adjuvant_prediction/cytokines_2clusters_compare.r")
```

## Figure 3E
```
source("SCRIPTS/adjuvant_prediction/ip10_2peak_scores.r")
source("SCRIPTS/adjuvant_prediction/2peaks_pca_final_heamap.r")
```
## Elastic net models

### Generate input data:
```
source("SCRIPTS/eNetXplorer/eNet_input_r1.r")
source("SCRIPTS/eNetXplorer/eNet_input_r2.r")
source("SCRIPTS/eNetXplorer/eNet_input_r3.r")
```
TODO: the sripts above require `RESULTS/Adjuvant_prediction/adjuvant_predicted_subjects.txt` which none of the previous scripts seems to produce.

### Run eNetXplorer:
```
source("SCRIPTS/eNetXplorer/eNetXplorer_R1_180530.R") # TODO: update conda environment with library missForest
source("SCRIPTS/eNetXplorer/eNetXplorer_R2_180530.R")
source("SCRIPTS/eNetXplorer/eNetXplorer_R3_180530.R")
```

### Figure 3B
```
source("SCRIPTS/eNet_figures/enet_plots_R1.r")
```

### Figure 3F and 3G
```
source("SCRIPTS/eNet_figures/enet_plots_R3.r")
source("SCRIPTS/eNet_figures/enet_plots_R2.r")
```

# Gene Expression PAXgene Data Processing

## Data Post Processing

### We found that two samples were switched. This is to correct it.

**TODO: I don't have the scripts in the workflow that generates "DATA_PROCESSED/Microarrays/PAXgene/eset.apt.RData". For now I will just copy it from Yuri's**
```
source("SCRIPTS/MA/filtering_pax/switch.samples/switch.samples.call.r")
```

### Apply different filtering to samples and genes. The probesets mapped to genes.
```
source("SCRIPTS/MA/filtering_pax/filtering.r")
```

### Calculate fold change from day 0
```
source("SCRIPTS/MA/calculate_d0_fc/calculate_d0_fc_pax.r")
```

# Baseline Data Analysis

## Preparing PBMC day 0 samples
```
source("SCRIPTS/MA/baseline_pbmc/d0_filter.r")
```

## WGCNA clustering of PBMC samples
```
source("SCRIPTS/MA/baseline_pbmc/d0_wgcna.r")
source("SCRIPTS/MA/baseline_pbmc/d0_wgcna_output.r")
```

## BTM enrichment analysis of data from PBMC samples
```
source("SCRIPTS/MA/baseline_pbmc/d0_wgcna_BTM_enrichment.r")
```

## Preparing whole blood (PAXgene) day 0 samples
```
source("SCRIPTS/MA/baseline_pax/d0_filter.r") # File name changed
```

## WGCNA clustering of whole blood samples
```
source("SCRIPTS/MA/baseline_pax/d0_wgcna.r")
source("SCRIPTS/MA/baseline_pax/d0_wgcna_output.r")
```

## BTM enrichment analysis of data from whole blood samples
```
source("SCRIPTS/MA/baseline_pax/d0_wgcna_BTM_enrichment.r")
```

## FIgure 4A. Combining BTM enrichment results from PBMC and whole blood samples
Also Suppl. Figure 6A
```
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

Figure 4B
```R
source("SCRIPTS/eNet_figures/enet_plots_all.r")
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

## Figure 6B

## Figure 6C

## Figure 6D

## Supp. Figure 7

# SOMAscan data analysis (Figure 7)
Julian’s scripts for elastic net models

## Figure 7A

## Figure 7B



