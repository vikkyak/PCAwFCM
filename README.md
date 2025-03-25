# 🧬 Principal Component Analysis with Fuzzy C-Means for Cell State Hierarchy

This repository provides an implementation of **PCA-based dimensionality reduction** combined with **Fuzzy C-Means (FCM) clustering** to infer a **cell state hierarchy** from single-cell RNA-seq gene expression data.

The method is demonstrated on the **Pollen et al. (2014)** dataset.

---

## 📁 Input Data

- `Pollen2014.txt`: Gene expression matrix (genes × cells)
- `SupplementaryLabels.txt`: Metadata with cell type and tissue labels

---

## 🔄 Workflow

1. **Load and preprocess expression data**
2. **Apply log2 transformation**
3. **Transpose to cell × gene matrix**
4. **Run PCA and FCM on reduced space**
5. **Repeat sampling and store clustering outputs**
6. **Evaluate clustering accuracy with Adjusted Rand Index (ARI)**
7. **Visualize performance across multiple cluster numbers**

---

## 📦 Required R Packages

```r
install.packages(c("ggplot2", "cluster", "fclust", "factoextra", "mgcv", "mnormt"))
BiocManager::install("ppclust")
BiocManager::install("pcaMethods")
devtools::install_github("justinhuang533/pcaReduce")

## 📥 Download

🔗 [**Download this repository as ZIP**](https://github.com/vikkyak/PCAwFCM/archive/refs/heads/main.zip)

Or clone via Git:

```bash
git clone https://github.com/vikkyak/PCAwFCM.git


▶️ Running the Code
1. Load data and libraries
r
Copy
Edit
b <- read.table("~/PCAwFCM/Pollen2014.txt", sep = ',', header = TRUE, row.names = 1)
lb <- read.table("~/PCAwFCM/SupplementaryLabels.txt", sep = ',', header = TRUE)
D <- log2(as.matrix(b) + 1)
Input <- t(D)
true_tissue_cls <- lb[, 6]
true_cell_cls <- lb[, 4]

2. Source main functions
source("~/PCAwFCM/SourceFun.R", echo = TRUE)

3. Run Fuzzy PCA Clustering
Output_SF <- Fpca_reduce(Input, nbt = 1, q = 30)

📊 Output
A list Output_SF containing FCM cluster labels across multiple runs

ARI scores for comparison against:

True cell-type clusters (K = 11)

True tissue-type clusters (K = 4)

Line plot showing ARI vs. number of clusters









