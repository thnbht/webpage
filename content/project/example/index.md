---
date: "2023-04-20T00:00:00Z"
external_link: ""
image:
  caption: Peripheral blood smear<sup>1</sup>
  focal_point: Smart
summary: Step-by-step clustering and visualization of peripheral blood mononuclear cells (PBMC).
tags:
- Immunology
- Data viz
title: Immune Cell Clustering
url_code: ""
url_pdf: ""
url_slides: ""
url_video: ""
---
## Table of contents
1. [Background and motivation](#background-and-motivation)
2. [Getting started](#getting-started)
3. [Data cleaning](#data-cleaning)
4. [Data transformation](#data-transformation)
5. [PCA](#pca)
6. [UMAP/tSNE](#umap-tsne)
7. [Annotation](#annotation)

**Background and motivation** <a name="background-and-motivation"></a>\
Human peripheral blood mononuclear cells (PBMCs) are any blood cell with a round nucleus. PBMCs include lymphocytes (T cells, B cells, and NK cells), monocytes, and dendritic cells. In humans, the frequencies of these populations vary across individuals, but typically, lymphocytes are in the range of 70-90%, monocytes from 10-20%, while dendritic cells are rare, accounting for only 1–2%. 
  
Most PBMCs are naïve or resting cells without effector functions. In the absence of an ongoing immune response T cells, the largest fraction of the isolated PBMCs, are mainly present as naïve or memory T cells. The naïve T-cells have never encountered their cognate antigen before and are commonly characterized by the absence of activation markers like CD25, CD44 or CD69 and the absence of the memory marker CD45RO isoform. Antigen recognition by a naïve T cell may result in activation of the cell, which will then enter a differentiation program and develop effector functions.<sup>2</sup>

Thus, we can differentiate the tissue microenvironment and understand its cellular composition by analyzing the biomarkers expressed through mRNA.
  
We will be working with a PBMC dataset from 10X Genomics, linked [here](https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz). The dataset captures 2,700 single cells sequenced on the Illumina NextSeq 500. We want to cluster cells that are similar to each other, identify cluster biomarkers (differentially expressed features), and then assign cell types to the clusters.

**Getting Started** <a name="getting-started"></a>\
We start by reading in the data. The Read10X() function reads in the output of the cellranger pipeline from 10X, returning a unique molecular identified (UMI) count matrix. The values in this matrix represent the number of molecules for each feature (i.e. gene; row) that are detected in each cell (column). We next use the count matrix to create a Seurat object. The object serves as a container that contains both data (like the count matrix) and analysis (like PCA, or clustering results) for a single-cell dataset.<sup>3</sup>

```R
library(tidyverse)
library(patchwork)
library(Seurat)
library(BiocManager)
BiocManager::install("edgeR")

# Load the dataset
pbmc_data <- Read10X(data.dir = "/.../filtered_gene_bc_matrices/hg19/")
# Initialize the Seurat object with the raw data
pbmc <- CreateSeuratObject(counts = pbmc_data, project = "PBMC", min.cells = 1, min.features = 100)
```

**Data Cleaning** <a name="data-cleaning"></a>\
Before proceeding to analyze the dataset, we must first clean it. We will filter cells based on common issues that arise in RNA sequencing and processing standards.

1. Low-quality or dying cells may display excessive mitochondrial contamination. To filter these cells, we will determine the percentage of gene counts mapped to a reference genome. We use the set of all genes starting with MT- as a set of mitochondrial genes.

```R
# Create a column in the Seurat object for percentage of mitochondrial genes mapped
pbmc[["MT_percent"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
```

  We need to set a threshold of cells – in this case, we will set it at 5% or greater mitochondrial counts.

2. Droplets meant to capture, lyse, and tag single cells may be empty, may contain low-quality cells, or may contain multiple cells. The first two cases result in a low unique gene count, whereas the latter case results in a high unique gene count. Accounting for these errors requires us to filter the cells for outliers.

```R
# Determine outliers with violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "MT_percent"), ncol = 3)
```

  We should filter cells that have unique gene counts less than 200 or over 2,500.

3. The number of unique genes should correlate with the number of molecules detected per cell.

```R
FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
```

  As expected, there is a strong correlation between the gene count and molecular count. This result assures us that 

Filtering the afore-mentioned cells:
```R
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & MT_percent < 5)
```

**Data Transformation** <a name="data-transformation"></a>\ <sup>4</sup>
Now that the data has been cleaned, we can move on to transforming the data into more useful forms for analysis.

1. In order to make accurate comparisions of gene expression between samples, we must normalize gene counts. The counts of mapped reads for each gene is proportional to the expression of mRNA (signal) in addition to many other factors (noise). Normalization is the process of scaling raw count values to account for the noise.

The most common factors to consider during normalization include:
-Sequencing depth: Accounting for sequencing depth is necessary for comparison of gene expression between samples. In the example below, each gene appears to have doubled in expression in Sample A relative to Sample B, however this is a consequence of Sample A having double the sequencing depth.

-Gene length: Accounting for gene length is necessary for comparing expression between different genes within the same sample. In the example, Gene X and Gene Y have similar levels of expression, but the number of reads mapped to Gene X would be many more than the number mapped to Gene Y because Gene X is longer.

-RNA composition: A few highly differentially expressed genes between samples, differences in the number of genes expressed between samples, or presence of contamination can skew some types of normalization methods. Accounting for RNA composition is recommended for accurate comparison of expression between samples, and is particularly important when performing differential expression analyses.

In the example, if we were to divide each sample by the total number of counts to normalize, the counts would be greatly skewed by the DE gene, which takes up most of the counts for Sample A, but not Sample B. Most other genes for Sample A would be divided by the larger number of total counts and appear to be less expressed than those same genes in Sample B.

There are several normalization methods that account for some of these differences, but the one we will use is trimmed mean of M values (TMM). This method uses a weighted trimmed mean (the average after removing the upper and lower x% of the data) of the log expression ratios between samples. TMM accounts for all three factors, and can be used in gene count comparisons _between_ and _within_ samples.

```R

```

2. Feature selection
3. Scaling

**PCA** <a name="pca"></a>\


**UMAP/tSNE** <a name="umap-tsne"></a>\


**Annotation** <a name="annotation"></a>\





_____
References:\
<sup>1</sup> Mantle cell leukemia as a cause of leukostasis (Smith et. al)\
<sup>2>/sup> https://www.ncbi.nlm.nih.gov/books/NBK500157/#:~:text=PBMCs%20include%20lymphocytes%20(T%20cells,for%20only%201%E2%80%932%20%25.\
<sup>3</sup> https://satijalab.org/seurat/articles/pbmc3k_tutorial.html#assigning-cell-type-identity-to-clusters \
<sup>4</sup> https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html
