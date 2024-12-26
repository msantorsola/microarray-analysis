# **Microarray Data Analysis Pipeline**

This repository contains a comprehensive pipeline for analyzing microarray data, particularly focused on Affymetrix platforms. It includes preprocessing, normalization, quality control, and differential gene expression (DGE) analysis. The pipeline uses R and several bioinformatics packages to process raw `.CEL` files, identify differentially expressed genes, and produce ready-to-analyze results.

## **Table of Contents**
1. [Overview](#overview)
2. [Pipeline Features](#pipeline-features)
3. [Requirements](#requirements)
4. [Usage](#usage)
5. [Outputs](#outputs)
6. [Repository Structure](#repository-structure)
7. [Acknowledgments](#acknowledgments)

---

## **Overview**

This pipeline is designed for working with Affymetrix microarray data. It processes raw probe-level data, applies normalization, performs statistical analysis, and identifies differentially expressed genes. 

### Key Analyses Performed:
- Data normalization using MAS 5.0 and RMA.
- Differential gene expression (DGE) analysis with S-Score and limma.

---

## **Pipeline Features**
- **Normalization**: MAS 5.0, RMA, and log2 transformations.
- **Differential Analysis**:
  - S-Score for pairwise comparisons.
  - limma for multi-group comparisons and ANOVA.
- **Quality Control**:
  - GAPDH and β-actin ratios.
  - Visualization of raw and normalized data.
- **Flexible Input**: Works with raw `.CEL` files from Affymetrix arrays.
- **Multiple Outputs**:
  - Normalized expression matrices.
  - Lists of differentially expressed genes (upregulated, downregulated).
  - Statistical results with adjusted p-values.

---

## **Requirements**
### **Software and Packages**
- **R** (≥ 4.0)
- Required R packages:
  - `affy`
  - `limma`
  - `sscore`
  - `GEOquery`
  - `R.utils`
  - `data.table`
  - `oligo`

### **Input Files**
- Raw `.CEL` files from Affymetrix microarrays.


