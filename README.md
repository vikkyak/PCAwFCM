# PCAwFCM: Fuzzy Clustering with Dimensionality Reduction

This repository contains a reproducible [Snakemake](https://snakemake.readthedocs.io/) pipeline for performing fuzzy clustering on single-cell data using PCA-based methods, integrating packages like `pcaReduce`, `ppclust`, and more.

## 📁 Repository Structure


## 🚀 Getting Started

### 1. Clone the repository

```bash
git clone https://github.com/vikkyak/PCAwFCM.git
cd PCAwFCM

2. Run the pipeline (with conda support)

This will:

Read data/Pollen2014.txt and data/SupplementaryLabels.txt

Run dimensionality reduction and fuzzy clustering

Save output plots in results/arandi_plot.png

3. Output
arandi_plot.png: ARI-based evaluation plot.

output.rds: RDS file with clustering results (if enabled in script).

🧪 Requirements
Python ≥ 3.8

R ≥ 4.1 with required packages (pcaReduce, ppclust, etc.)

Snakemake

Optional: Conda (for isolated environments)
