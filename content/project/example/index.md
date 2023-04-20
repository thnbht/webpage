---
date: "2023-04-20T00:00:00Z"
external_link: ""
image:
  caption: Peripheral blood smear, Mantle cell leukemia as a cause of leukostasis (Smith et. al)
  focal_point: Smart
summary: Step-by-step clustering and visualization of peripheral blood mononuclear cells (PBMC).
tags:
- ???
title: Immune Cell Clustering
url_code: ""
url_pdf: ""
url_slides: ""
url_video: ""
---
**Background and Motivation**\
Human peripheral blood mononuclear cells (PBMCs) are isolated from peripheral blood and identified as any blood cell with a round nucleus (i.e. lymphocytes, monocytes, natural killer cells (NK cells) or dendritic cells). PBMCs include lymphocytes (T cells, B cells, and NK cells), monocytes, and dendritic cells. In humans, the frequencies of these populations vary across individuals, but typically, lymphocytes are in the range of 70-90%, monocytes from 10-20%, while dendritic cells are rare, accounting for only 1–2%. Most of the PBMCs are naïve or resting cells without effector functions. In the absence of an ongoing immune response T cells, the largest fraction of the isolated PBMCs, are mainly present as naïve or memory T cells. The naïve T-cells have never encountered their cognate antigen before and are commonly characterized by the absence of activation markers like CD25, CD44 or CD69 and the absence of the memory marker CD45RO isoform. Antigen recognition by a naïve T cell may result in activation of the cell, which will then enter a differentiation program and develop effector functions. After activation the CD4+ T cell subset may develop into diverse effector cell subsets, including Th1, Th2, Th17, Th9, Th22, follicular helper (Tfh) cells and different types of regulatory cells (Akdis et al. 2012; Crotty 2011; Tan and Gery 2012; Sakaguchi et al. 2008). The CD4+ helper T cells are essential mediators of immune homeostasis and inflammation (Hirahara et al. 2013).<sup>1</sup>

We will be working with a PBMC dataset from 10X Genomics, linked [here] (https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz). The dataset captures 2,700 single cells sequenced on the Illumina NextSeq 500. We want to cluster cells that are similar to each other, identify cluster biomarkers (differentially expressed features), and then assign cell types to the clusters.

**Getting Started**\
We start by reading in the data. The Read10X() function reads in the output of the cellranger pipeline from 10X, returning a unique molecular identified (UMI) count matrix. The values in this matrix represent the number of molecules for each feature (i.e. gene; row) that are detected in each cell (column). We next use the count matrix to create a Seurat object. The object serves as a container that contains both data (like the count matrix) and analysis (like PCA, or clustering results) for a single-cell dataset.<sup>2</sup>

```R
library(tidyverse)
library(patchwork)
library(Seurat)

# Load the dataset
pbmc_data <- Read10X(data.dir = "/.../filtered_gene_bc_matrices/hg19/")
# Initialize the Seurat object with the raw data
pbmc <- CreateSeuratObject(counts = pbmc_data, project = "PBMC", min.cells = 1, min.features = 100)
```

**Data Cleaning**\
Before proceeding to analyze the dataset, we must first clean it. We will filter cells based on common issues that arise in RNA sequencing and processing standards.

1. Low-quality or dying cells may display excessive mitochondrial contamination. To filter these cells, we will determine the percentage of gene counts mapped to a reference genome. We use the set of all genes starting with MT- as a set of mitochondrial genes.

```R
# Create a column in the Seurat object for percentage of mitochondrial genes mapped
pbmc[["MT_percent"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
```

We set a threshold of cells with 5% or greater mitochondrial counts.

2. Droplets meant to capture, lyse, and tag single cells may be empty, may contain low-quality cells, or may contain multiple cells. The first two cases result in a low unique gene count, whereas the latter case results in a high unique gene count. Accounting for these errors requires us to filter the cells for outliers.

```R
# Determine outliers with violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "MT_percent"), ncol = 3)
```

We filter cells that have unique gene counts less than 200 or over 2,500.

3. The number of unique genes should correlate with the number of molecules detected per cell.

```R
FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
```

The correlation is strong, as predicted.

Filtering the afore-mentioned cells:
```R
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & MT_percent < 5)
```







_____
References:
<sup>1</sup> https://www.ncbi.nlm.nih.gov/books/NBK500157/#:~:text=PBMCs%20include%20lymphocytes%20(T%20cells,for%20only%201%E2%80%932%20%25.
<sup>2</sup> https://satijalab.org/seurat/articles/pbmc3k_tutorial.html#assigning-cell-type-identity-to-clusters

