# Project Description

First time running a genomics online CURE.  We tested range of trimming parameters trimq and minlen in software bbduk on sex differences in gene expression seen in batch 1 of the Yale placenta data set (Olney et al, Sex differences in early and term placenta are conserved in adult tissues, Biol Sex Differ
. 2022 Dec 22;13(1):74. doi: 10.1186/s13293-022-00470-y)

# Directory structure

## data_generation

Contains all snakefiles used to generate and process trimmed data so it can be passed into differential expression pipeline.

## cluster_submission

Contains sbatch scripts used to run snakefiles on ASU Agave biocomputing cluster

## Differential_Expression

Contains Rmd files of the limma-voom differential expression pipeline, modified to run analysis for the CURE.  This includes code to do pairwise comparison of two differentially expressed gene lists and an upset comparing the entire range of data sets.
