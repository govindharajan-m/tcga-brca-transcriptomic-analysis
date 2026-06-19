# tcga-brca-transcriptomic-analysis
This project investigates transcriptomic differences between Early-stage and Advanced-stage breast cancer using RNA-sequencing data from The Cancer Genome Atlas Breast Invasive Carcinoma (TCGA-BRCA) cohort.

The study was conducted using R and Bioconductor packages including TCGAbiolinks, DESeq2, survival, and survminer.

## Objectives

* Compare gene expression profiles between Early-stage and Advanced-stage breast cancer.
* Identify significantly differentially expressed genes.
* Visualize transcriptomic alterations using volcano plots and heatmaps.
* Evaluate stage-associated survival differences using Kaplan-Meier analysis.

## Dataset

Source:

The Cancer Genome Atlas (TCGA-BRCA)

Samples analyzed:

* Early-stage tumors: 178
* Advanced-stage tumors: 268
* Total matched RNA-seq samples: 446

### Data Acquisition

TCGA-BRCA clinical and RNA-seq data were retrieved using TCGAbiolinks.

### Differential Expression Analysis

Differential expression analysis was performed using DESeq2.

Criteria:

* FDR < 0.05
* |log2FC| > 1

### Survival Analysis

Kaplan-Meier survival analysis was performed using the survival and survminer packages.

## Key Results

* Total genes analyzed: 53,020
* Significant genes: 670
* Strongly upregulated genes: 76
* Strongly downregulated genes: 111

Top differentially expressed genes included:

* CHGA
* KRT6A
* KRT13
* KRT20
* CR2
* SCG3
