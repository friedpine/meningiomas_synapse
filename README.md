This repository contains the custom R scripts used for the single-cell RNA sequencing (scRNA-seq) analysis in the manuscript "Glutamatergic synaptic-like input to intracranial meningiomas cells drives tumor cell proliferation".

The provided pipeline utilizes the Seurat framework to perform downstream analyses on meningioma samples, including global cell clustering, differential expression analysis, quantification of cellular proportions, targeted gene set screening, and Gene Set Variation Analysis (GSVA).


To reproduce the downstream analysis described in this repository, please download the following pre-processed Seurat object files and place them in your root working directory:

gse18v2_ana_1227-935104.qs: The global single-cell dataset containing all cell populations.

gse18_TC_sub-935104-3_1227.qs: The subset dataset comprising strictly the isolated Tumor Cells.

## System Requirements
Hardware Requirements
scRNA-seq analysis requires a standard computer with sufficient RAM to support in-memory operations.

RAM: Minimum 16 GB (32 GB or higher is highly recommended for GSVA and large object manipulation).

Storage: At least 10 GB of available disk space.

Software Requirements
The script was developed and tested in the following environment:

```
R (version >= 4.1.0)

Seurat (version >= 4.0.0)

GSVA (version >= 1.40.1)

qs, ggplot2, dplyr, tidyr, ggpubr

To install the required dependencies in R, run:

R
install.packages(c("qs", "ggplot2", "dplyr", "tidyr", "ggpubr", "Seurat"))
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("GSVA")
```

## Code and Figure Mapping
The main analysis script (meningioma_scRNA_analysis.R) is divided into several modules. Each module corresponds to specific figures and supplementary figures generated in our manuscript:

## Module 1 & 2: Global Cell Type Annotation

Related Figures: Figure 1B-C, Supplementary Figure S1

Description: Loads the global .qs dataset, standardizes metadata (Dura vs. Meningioma), and visualizes canonical marker gene expression across all identified clusters using UMAPs and DotPlots.

## Module 3: Cellular Composition Analysis

Related Figures: Figure 2A-B

Description: Calculates the absolute abundance and relative cellular proportions of different cell clusters across conditions, visualizing the tumor microenvironment composition via stacked bar charts. (Note: Sample size denominators are configured based on the study's specific cohort size).

## Module 4 & 6: SY and GLU Gene Set Screening

Related Figures: Figure 3C, Figure 4A

Description: Evaluates the expression profiles of specific synaptic (SY) and glutamatergic (GLU) gene sets within the combined tumor/glial clusters. Generates expression percentage matrices and custom visualization plots.

## Module 5 & 7: Gene Set Variation Analysis (GSVA)

Related Figures: Figure 4D-E

Description: Extracts the log-normalized RNA assay matrix to compute single-sample GSEA (ssGSEA) scores for the SYscore gene module. Computes statistical significance across clusters using Wilcoxon and Kruskal-Wallis tests.

## Module 8 & 9: Tumor Subpopulation Profiling

Related Figures: Figure 5A-F

Description: Loads the tumor-specific Seurat object (gse18_TC_sub). Analyzes intra-tumoral heterogeneity, identifies differentially expressed genes (DEGs) for distinct subclusters (e.g., proliferating, metabolic), and maps spatial/regional proliferation markers (MKI67) across distinct spatial origins (MSC_BTI vs. MSC_Core).


Execute the code sequentially. Output plots will be rendered in the R environment, and generated CSV tables (e.g., count_df.csv, markersC1.csv) will be saved directly to the directory.

## Contact
For any questions regarding the code, data, or analysis pipeline, please open an issue in this repository or contact the corresponding author at xiannian@ccmu.edu.cn 
