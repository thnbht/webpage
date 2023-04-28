---
date: "2023-04-20T00:00:00Z"
external_link: ""
image:
  caption: Peripheral blood smear<sup>1</sup>
  focal_point: Smart
summary: Step-by-step clustering and visualization of peripheral blood mononuclear cells (PBMC).
tags:
- Data viz
- Single cell
title: Single-Cell RNA-Seq: Unsupervised Clustering
url_code: ""
url_pdf: ""
url_slides: ""
url_video: ""
---
## Table of contents
1. [Introduction](#introduction)
2. [Getting started](#getting-started)
3. [Data cleaning](#data-cleaning)
4. [Data transformation](#data-transformation)
5. [PCA](#pca)
6. [UMAP/tSNE](#umap-tsne)
7. [Annotation](#annotation)
-----
**Introduction** <a name="introduction"></a>\
Single cell RNA sequencing (scRNA-seq) is a method for measuring gene expression by quantifying RNA transcripts. In contrast to bulk RNA sequencing, which captures average global gene expression, scRNA-seq captures heterogenous cell composition and makes it possible to model cell-cell interactions and molecular microenvironments.

Although there are a variety of scRNA-seq protocols, they all follow a general design. First, the cells must be dissociate and be isolated. The isolated cells are then encapsulated by droplets or wells. After encapsulation, the cell is lysed and the transcripts within are barcoded with short nucleotide tags referencing from which cell the given transcript originated. From there, the transcripts are then amplified and quantified, typically through qPCR.

Upon quantification, the data must first be cleaned through a series of identifying cells of interest. Next, several methods for minimizing noise, both statistical and technical, are employed. Due to the high dimensionality inherent to genomic data, dimensionality reduction is used to aid in ease of analysis. Finally, the analysis of interest can be carried out.
  
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
- Sequencing depth: Accounting for sequencing depth is necessary for comparison of gene expression between samples. In the example below, each gene appears to have doubled in expression in Sample A relative to Sample B, however this is a consequence of Sample A having double the sequencing depth.

- Gene length: Accounting for gene length is necessary for comparing expression between different genes within the same sample. In the example, Gene X and Gene Y have similar levels of expression, but the number of reads mapped to Gene X would be many more than the number mapped to Gene Y because Gene X is longer.

- RNA composition: A few highly differentially expressed genes between samples, differences in the number of genes expressed between samples, or presence of contamination can skew some types of normalization methods. Accounting for RNA composition is recommended for accurate comparison of expression between samples, and is particularly important when performing differential expression analyses.

In the example, if we were to divide each sample by the total number of counts to normalize, the counts would be greatly skewed by the DE gene, which takes up most of the counts for Sample A, but not Sample B. Most other genes for Sample A would be divided by the larger number of total counts and appear to be less expressed than those same genes in Sample B.

There are several normalization methods that account for some of these differences, but the one we will use is trimmed mean of M values (TMM). This method uses a weighted trimmed mean (the average after removing the upper and lower x% of the data) of the log expression ratios between samples. TMM accounts for all three factors, and can be used in gene count comparisons _between_ and _within_ samples.

```R
# TMM Code
```

2. We next calculate a subset of features that exhibit high cell-to-cell variation in the dataset (i.e, they are highly expressed in some cells, and lowly expressed in others). This is to distinguish true biological variability from the high levels of technical noise in single-cell experiments.<sup>5</sup>

Selection methods include (Seurat doc)
- vst: First, fits a line to the relationship of log(variance) and log(mean) using local polynomial regression (loess). Then standardizes the feature values using the observed mean and expected variance (given by the fitted line). Feature variance is then calculated on the standardized values after clipping to a maximum (see clip.max parameter).

- mean.var.plot (mvp): First, uses a function to calculate average expression (mean.function) and dispersion (dispersion.function) for each feature. Next, divides features into num.bin (deafult 20) bins based on their average expression, and calculates z-scores for dispersion within each bin. The purpose of this is to identify variable features while controlling for the strong relationship between variability and average expression.

```R
# plot1 plot2 diff selection methods
```

3. Next, we apply a linear transformation (‘scaling’) that is a standard pre-processing step prior to dimensional reduction techniques like PCA. The ScaleData() function:

- Shifts the expression of each gene, so that the mean expression across cells is 0
- Scales the expression of each gene, so that the variance across cells is 1
-- This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate

The results of this are stored in pbmc[["RNA"]]@scale.data

```R
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
```

**PCA** <a name="pca"></a>\


**UMAP/tSNE** <a name="umap-tsne"></a>\


**Annotation** <a name="annotation"></a>\





_____
References:\
<sup>1</sup> Mantle cell leukemia as a cause of leukostasis (Smith et. al)\
<sup>2>/sup> https://www.ncbi.nlm.nih.gov/books/NBK500157/#:~:text=PBMCs%20include%20lymphocytes%20(T%20cells,for%20only%201%E2%80%932%20%25.\
<sup>3</sup> https://satijalab.org/seurat/articles/pbmc3k_tutorial.html#assigning-cell-type-identity-to-clusters \
<sup>4</sup> https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html \
<sup>5</sup> https://www.nature.com/articles/nmeth.2645
