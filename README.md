# scRNA-seq_SARS-CoV2_DVG_analysis
Seurat analysis workflow of scRNA-seq data from Cellranger, STARsolo

This repository contains scRNA-seq analysis workflow for Sun Lab, URMC. It details how to use gene expression matrices from the Cellranger count output, and STARsolo output in Seurat: 
  to cluster and type identify cells, to produce tSNE and feature plots to visualize cell type clusters and the presence of SARS-CoV2 in those clusters, 
  to produce lists of differentially expressed genes between SARS-CoV2 infected cells with and without defective viral genomes, 
  and to produce expression heatmaps of those differentially expressed genes and gene ontology enrichment plots.

## Instructions
See SOP for all code necessary. 

The R scripts in this repository are for differential gene expression analysis
