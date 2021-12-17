# Genomics_CURE

## Pipeline Objective

The Genomics CURE pipeline is a reproducible workflow for analyzing the full-term placenta dataset for sex differences. The goal is to run a series of  snakemakes with modified parameters to investigate how inferences in differential expression will be affected. The default pipeline will be run for quality control, trimming, aligning, and generating counts; additional sets of edited snakemakes have modified parameters to compare the effects of parameter differences on gene expression.


#### Conda Environment

This will activate the conda environment `placenta_RNA_env3` which contains the packages needed to complete this workflow:
```
conda activate --file placenta_RNA_env3
```


## Part I: Processing the default pipeline 

### Step 1: Trimming and quality control
#### 1. Perform and visualize quality control on the raw data using FastQC and MultiQC
#### 2. Trim reads by quality and trim the adapters using Trimmomatic (currently uses bbduk)
- Sliding window approach
#### 3. Visualize quality post-trimming using FastQC and MultiQC

### Step 2: Reference genome alignment and obtaining transcript counts
#### 1. Align trimmed .fastqs to the GRCh38.p12 sex chromosomome complement human reference genome using HISAT 2.1.0
#### 2. Quantify differences in gene expression using featureCounts 1.5.2

### Step 3: Data Visualization
#### 1. Run the r file `R_limmaVooom_placentaSexDiff.r` for normalization and differential expression

### Step 4: Interpret gene expression data and summarize results
#### 1.  Plot log2 fold changes to evalute gene expression



## Part II: Run edited snakemakes with modified parameters 



## Part III: Comparing if modified parameters alter gene expression analysis




