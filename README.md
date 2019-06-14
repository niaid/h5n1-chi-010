# h5n1
H5N1 

# Workflow

*Note: Output folders need to be created manually, therefore, read the output data text to create appropriate directory structure. TODO: probably add a feature in the scripts to check and create directories before writing files. Or write a master shell script that creates the appropriate directories for the entire project.*

# Gene Expression PBMC data processing

*Note: I did not repeat these steps as they seem to compute intensive and also stochastic in nature.*

First, the CEL files were processed with  Power Tools RNA-sketch algorithm. 
```
source("SCRIPTS/MA/processing_pbmc/apt.config.r")
source("SCRIPTS/MA/processing_pbmc/apt.call.r")
```

Quantile normalization implement in RMA-sketch algorithm uses random subset of data to decrease memory usage. Because of this at different runs the algorithm may generates slightly different results. For reproducibility we recommend using precomputed APT output located in DATA_PROCESSED/Microarrays/PBMC/apt_summary.txt. *(this seem to be wrong because apt_summary.txt is produced by the scripts in the section immediately below. On the other hand rma-sketch.summary.txt file seem to be the optout of this procedure and input to ones below. I copied rma-sketch.summary.txt file over.)*


## Create ExpressionSet from APT output

*Note: copy rma-sketch.summary.txt file to DATA_PROCESSED/Microarrays/PBMC directory if you are skipping the steps in the section above.*

Input data:  
* `DATA_ORIGINAL/Microarrays/PBMC/sample_info.txt` from the original data files. 
* `DATA_PROCESSED/Microarrays/PBMC/rma-sketch.summary.txt` from the previous procedure.

```
source("SCRIPTS/MA/processing_pbmc/eset.config.r")
source("SCRIPTS/MA/processing_pbmc/eset.call.r")
```

Output dat:  
* `DATA_PROCESSED/Microarrays/PBMC/apt.summary.txt`
* `DATA_PROCESSED/Microarrays/PBMC/eset.apt.RData`


## Probesets to genes mapping

### Convert the table to probeset-gene mapping  
Input data:  
* `DATA_ORIGINAL/Microarrays/annotation/HuGene-2_1-st-v1.na35.hg19.transcript.csv` file downloaded from Affymetrix web site
```
source("SCRIPTS/MA/annotation/affy_hugene-2_1-st_annotation.r")
```
Output data:  
* `DATA_PROCESSED/Microarrays/annotation/affy_hugene-2_1-st_unique_all.txt`
* `DATA_PROCESSED/Microarrays/annotation/affy_hugene-2_1-st_unique.txt`

### Select the best probeset for a gene
Input data:  
* `DATA_PROCESSED/Microarrays/annotation/affy_hugene-2_1-st_unique.txt`
* `DATA_PROCESSED/Microarrays/PBMC/eset.apt.RData`
```
source("SCRIPTS/MA/annotation/generate_ps2gene_map.config.r")
source("SCRIPTS/MA/annotation/generate_ps2gene_map.call.r") # It calls a function from this script SCRIPTS/functions/pick.probeset.r TODO: I don't understand if it is producing any file
```

## Data post processing
### We found that two samples were switched. This is to correct it.

Input data:  
* `SCRIPTS/MA/filtering_pbmc/samples.switch.tx`
```
source("SCRIPTS/MA/filtering_pbmc/switch.samples/switch.samples.call.r") # Fetches a function from this script SCRIPTS/MA/filtering_pbmc/switch.samples
```
Output data:
* `DATA_PROCESSED/Microarrays/PBMC/eset.corrected.RData`

### Apply different filtering to samples and genes. The probesets mapped to genes.

*Note: these scripts source multiple scripts within their respective folders. TODO: it looks messy and hard to figure out what are the inputs and outputs.*
```
source("SCRIPTS/MA/filtering_pbmc/samples.clean_genes.all/filtering.r")
source("SCRIPTS/MA/filtering_pbmc/samples.clean_genes.iqr/filtering.r")
source("SCRIPTS/MA/filtering_pbmc/samples.all_genes.all/filtering.r")
source("SCRIPTS/MA/filtering_pbmc/samples.all_genes.iqr/filtering.r")
```

### Calculate fold change from day 0

Input data:  
* `DATA_PROCESSED/Microarrays/PBMC/samples.clean_genes.iqr`
* `DATA_PROCESSED/Microarrays/PBMC/samples.all_genes.all/eset.genes.filtered.RData`
```
source("SCRIPTS/MA/calculate_d0_fc/calculate_d0_fc_pbmc.r") # Sourcing this file SCRIPTS/functions/factor.date.r
```
Output data:  
* `DATA_PROCESSED/Microarrays/PBMC/samples.clean_genes.iqr/gexp_d0_fc.RData`

# Pattern discovery in post-vaccination profiles of gene expression

## Profiles clustering with DIANA

*Note: before running this script create the relevant folders by `mkdir RESULTS/Microarrays/PBMC/pattern_discovery/`*

Input data:  
* `DATA_PROCESSED/Microarrays/PBMC/samples.clean_genes.iqr/gexp_d0_fc.RData`
```
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
```
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
```
source("SCRIPTS/MA/pattern_discovery/pattern_filter.r")
```
Output data:  
* `/RESULTS/Microarrays/PBMC/pattern_discovery/GE_patterns_filt.txt`
  
## Summarize patterns stats
Input data:
* `RESULTS/Microarrays/PBMC/pattern_discovery/df.mat.rds`
* `/RESULTS/Microarrays/PBMC/pattern_discovery/GE_patterns_filt.txt`
```
source("SCRIPTS/MA/pattern_discovery/patterns_stats.r")
```
Output data:
* `/RESULTS/Microarrays/PBMC/pattern_discovery/df.cl.stat.rds`
* `/RESULTS/Microarrays/PBMC/pattern_discovery/df.patt.rds`
  
## Figure 2B - patterns profile plot
Input data:
* `RESULTS/Microarrays/PBMC/pattern_discovery/GE_patterns_filt.txt`
* `RESULTS/Microarrays/PBMC/pattern_discovery/df.cl.stat.rds`
```
source("SCRIPTS/MA/pattern_discovery/plot_patterns.r")
```
Output data:
* `FIGURES/GE_patterns_profiles.png`
* `FIGURES/GE_patterns_profiles_2col.png`
* `FIGURES/GE_patterns_profiles_horiz.png`