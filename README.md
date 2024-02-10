# EEEC_landscape

This repository contains metadata and codes necessary for analysis of the proteogenomics data of Endometrial Carcinoma presented in Hu et al. Nature Genetics (2024).

## Data availability

- Raw sequencing data generated in this study are deposited in Genome Sequence Archive for Human (https://ngdc.cncb.ac.cn/gsa-human/) with accession number HRA003319. 
- The mass spectrometry proteomics data have been deposited to the ProteomeXchange Consortium via the PRIDE partner repository with the dataset identifier PXD046507. 
- TCGA-UCEC and CPTAC-UCEC data was obtained from Cbioportal (https://www.cbioportal.org/) with the dataset identifier ucec_tcga_pan_can_atlas_2018 and ucec_cptac_2020.
- Panel sequencing of genomic data from AACR GENIE Project was obtained from Synapse (https://www.synapse.org/) with the dataset identifier syn7222066. 
- RNA-seq data from six public datasets (GSE1378, GSE1379, GSE6532, GSE9195, GSE12093, GSE17705) of patients with breast cancer treated with endocrine therapy was obtrained from GEO.

## Data visualization

### Requirements

Tested on macOS Ventura and CentOS

1. R version: 4.1.2
2. R packages
   - ggplot2
   - data.table
   - Seurat
   - dplyr
   - tidyr
   - ggpubr
   - RColorBrewer
   - pheatmap
   - ggsignif
   - ComplexHeatmap
