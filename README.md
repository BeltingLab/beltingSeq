# Adapted Acidosis

## Overview

This R project performs comprehensive transcriptomic analysis of glioblastoma (GBM) cell lines and patient samples exposed to acidic microenvironments. The analysis focuses on lipid metabolism-related pathways, including CSPG remodeling, lipid uptake, lipid droplet formation, and ferroptosis-associated gene signatures.

The project integrates multiple microarray datasets from NCBI GEO and utilizes advanced bioinformatics approaches including differential gene expression analysis, gene set enrichment analysis (GSEA), pathway clustering, and gene signature validation across independent GBM cohorts.

# 1.  Transcriptomic Analysis of the effects of acidosis on Glioblastoma cell lines and patient samples

- [Adapted Acidosis](#adapted-acidosis)
  - [Overview](#overview)
- [1.  Transcriptomic Analysis of the effects of acidosis on Glioblastoma cell lines and patient samples](#1--transcriptomic-analysis-of-the-effects-of-acidosis-on-glioblastoma-cell-lines-and-patient-samples)
  - [1.1. Input data](#11-input-data)
  - [1.2. Directory tree](#12-directory-tree)
- [2. Bioinformatical pipeline](#2-bioinformatical-pipeline)
  - [2.1. Data preprocessing \& differential gene expression (DGE) analysis](#21-data-preprocessing--differential-gene-expression-dge-analysis)
      - [Overview:](#overview-1)
      - [Main steps:](#main-steps)
      - [Input data:](#input-data)
      - [Output:](#output)
  - [2.2. Pathway analysis](#22-pathway-analysis)
      - [Overview:](#overview-2)
      - [Main Steps:](#main-steps-1)
      - [Input data:](#input-data-1)
      - [Output:](#output-1)
  - [2.3. GSEA clustering](#23-gsea-clustering)
      - [Overview:](#overview-3)
      - [Main Steps:](#main-steps-2)
      - [Input data:](#input-data-2)
      - [Output:](#output-2)
  - [2.4. Gene signature analysis](#24-gene-signature-analysis)
      - [Overview:](#overview-4)
      - [Main Steps:](#main-steps-3)
      - [Inputs:](#inputs)
      - [Outputs:](#outputs)
  - [3. Results](#3-results)
  - [4. Getting Started](#4-getting-started)
    - [Prerequisites](#prerequisites)
    - [Installation](#installation)
    - [Usage](#usage)
  - [5. Data Availability](#5-data-availability)
  - [6. Session Information](#6-session-information)
  - [7. Citation](#7-citation)
  - [8. License](#8-license)
  - [9. Contact](#9-contact)

## 1.1. Input data

Four main data sources were used in these analyses:

1. **CCLD** - Affymetrix Clariom D Pico Gene Array performed on pooled laser-microdissected samples from lipid droplet rich (n=5) and lipid droplet lacking (n=5) matched pathological samples.

2. **HGCC** - Affymetrix Clariom D Pico Gene Array on three patient-derived primary glioblastoma cell lines grown as 3D vs. 2D cultures.

3. **PANC1** - Affymetrix Clariom D Pico Gene Array on pancreatic cancer cell line PANC1, stimulated by a low-pH environment (pH 6.4) for 10 weeks (adapted acidosis, AA).

4. **U87** - Illumina HumanHT-12 v4 Expression BeadChip for gene expression analysis of U87 glioblastoma cells grown in low-pH (acidosis) or hypoxic conditions.

5. **Validation datasets**: IvyGAP glioblastoma dataset for gene signature validation.

> [!NOTE]
> *All raw and processed data used in this analysis have been uploaded to NCBI's public data repository, the Gene Expression Omnibus (GEO).*

The data are available under the unique IDs **GSE300758** (**CCLD**), **GSE300765** (**U87**), **GSE300768** (**PANC1**), and **GSE300771** (**HGCC**). The accession of the datasets is scheduled to stay private until 30th September, 2025, unless the manuscript enters review earlier.

## 1.2. Directory tree

```bash
.
├── adapted-acidosis.Rproj
├── RData
│   ├── CCLD_enrichment_results.RData
│   ├── CCLD_processedData.RData
│   ├── CGGA_processedData.RData
│   ├── HGCC_enrichment_results.RData
│   ├── HGCC_processedData.RData
│   ├── IvyGap_processedData.RData
│   ├── MSigDB_gene_sets.RData
│   ├── PANC1_enrichment_results.RData
│   ├── PANC1_processedData.RData
│   ├── U87_enrichment_results.RData
│   ├── U87_processedData.RData
│   └── term_similarity_matrix.RData
├── README.md
├── Results
│   ├── Fig-1C
│   ├── Fig-1D
│   ├── Fig-1F
│   ├── Fig-2A
│   ├── Fig-2D
│   ├── Fig-2I
│   ├── Supplementary-Fig-1A
│   ├── Supplementary-Fig-3B
│   ├── Supplementary-Fig-3I
│   ├── Supplementary-Fig-4A
│   ├── Supplementary-Fig-4E
│   ├── Supplementary-Fig-5A
│   └── Supplementary-Fig-5F
├── bin
│   ├── 1-preprocessing.R
│   ├── 2-GSEA-analysis.R
│   ├── 3-cluster-terms.R
│   ├── 4-signature-analysis.R
│   ├── functions.R
│   └── packages.R
├── etc
│   ├── cluster-names.txt
│   ├── gene-signature.txt
│   ├── genes-of-interest.txt
│   ├── pathways-of-interest.txt
│   └── IvyGap
│       ├── 2025-01-23_Ivy_GAP_expression.csv
│       ├── 2025-01-23_Ivy_GAP_genes.csv
│       ├── 2025-01-23_Ivy_GAP_pheno.csv
│       └── README.txt
├── renv
│   ├── activate.R
│   ├── fix-renv.R
│   ├── fix-restore.R
│   ├── library/
│   └── settings.json
├── renv.lock
├── .Rprofile
├── .Renviron
└── sessionInfo.txt
```

# 2. Bioinformatical pipeline 

## 2.1. Data preprocessing & differential gene expression (DGE) analysis

#### Overview:
The script [`1-preprocessing.R`](bin/1-preprocessing.R) performs preprocessing, normalization, and differential expression analysis (DEA) for multiple transcriptomic datasets from Affymetrix Clariom D Human Pico and Illumina BeadChip platforms.

#### Main steps:
1. Environment Setup:
     - Sets working directory and loads required packages and custom functions.

2. DGE Analysis each dataset:
     - Svenja's CC +LD and CC no LD (Affymetrix)
     - Hugo's Primary cells 2D vs 3D (Affymetrix)
     - U87 Chronic Acidosis AA vs NA & HOX vs NOX (Illumina)
     - PANC1 Chronic Acidosis AA vs NA (Affymetrix)

For each dataset:
  - Creates distinct results (tables) directories.
  - Loads raw files.
  - Normalize raw data and perform transcript filtering.
  - Conduct differential expression analysis using limma.
  - Remove duplicate gene symbols/IDs based on log2 fold change or t-stat.
  - Summarize up- and down-regulated genes for each comparison.
  - Save processed data as RData and Excel files for downstream analysis.

Requirements:
  - R packages: oligo, affycoretools, limma, stats
  - Annotations: clariomdhumantranscriptcluster.db (v8.8.0), illuminaHumanv4.db (v1.26.0), org.Hs.eg.db (v3.20.0)
  - Custom functions (from [`functions.R`](bin/functions.R)): normalizeTranscript, normalizeIllumina, limmaDEA, removeDuplicates

#### Input data:
  - Raw `.cel` files for Affirmetrix experiments, and calculated probe and control intensity summaries for the Illumina experiment

#### Output:
  - Processed RData for each dataset and comparison

## 2.2. Pathway analysis

#### Overview:
The script [`2-GSEA-analysis.R`](bin/2-GSEA-analysis.R) performs comprehensive gene set enrichment analysis (GSEA) in a modular fashion, supporting all included datasets (Affymetrix Clariom D, Illumina BeadChip) and experimental designs. It executes GSEA against several gene set collections from MSigDB (Hallmark, GO:BP, KEGG, Reactome).

#### Main Steps:
1. Environment Setup:
     - Sets working directory and loads required packages and custom functions.
     - Loads MSigDB gene set databases (Hallmark, GO:BP, KEGG, Reactome) and saves it as RData.

2. Pathway Analysis for each dataset:
     - Svenja's CC +LD and CC no LD (Affymetrix)
     - Hugo's Primary cells 2D vs 3D (Affymetrix)
     - U87 Chronic Acidosis AA vs NA & HOX vs NOX (Illumina)
     - PANC1 Chronic Acidosis AA vs NA (Affymetrix)

For each dataset:
  - Identifies and annotates DEGs
  - Extracts gene lists for enrichment analysis and runs Gene Set Enrichment Analysis (GSEA) using MSigDB gene sets.

Requirements:
  - R packages: dplyr, tidyr, stringr, openxlsx, msigdbr, org.Hs.eg.db, cowplot, ggplot2, etc.
  - Custom functions (from [`functions.R`](bin/functions.R)): get_significance, get_genelist, run_ora, extract_ora_results, run_gsea, extract_gsea_results, plot_pca, plot_venn, plot_vulcan, getEnrichmentTable, plotRunningScore, plotGeneRank, plotClusters.

#### Input data:
  - RData files with expression matrices and DEG tables for each dataset.

#### Output:
  - RData with all MSigDB gene sets used for analysis.
  - RData with pathway enrichment results for each dataset.

## 2.3. GSEA clustering

#### Overview:

The [`3-cluster-terms.R`](bin/3-cluster-terms.R) loads all gene sets that we used for GSEA analyses, and calculate Cohen's similarity matrix between all gene sets (GO terms, MSigDb hallmarks, KEGG, Reactome pathways).

#### Main Steps:

1. Calculate Cohen's similarity matrix between all gene sets (GO terms, MSigDb hallmarks, KEGG, Reactome pathways).

Requirements:
  - R packages: igraph, RCy3, etc.
  - custom functions: cohen_kappa, get_cluster, get_cluster_representative, filter_graph, plot_network, getCircplotData, plotCircplot, plot_vulcan

#### Input data:
  - RData with all MSigDB gene sets used for analysis.

#### Output:
  - RData with the similarity matrix between all terms.
 

## 2.4. Gene signature analysis

#### Overview:
The [`4-signature-analysis.R`](bin/4-signature-analysis.R) validates identified gene signatures across independent glioblastoma cohorts. This script preprocesses the Ivy Glioblastoma Atlas Project (IvyGAP), TCGA-GBM, and CGGA datasets, filters for primary GBM samples, and performs gene signature scoring and correlation analyses.

**Data Preprocessing:**
- Normalized count reads from pre-processed data (sequence alignment and transcript abundance estimation)
- Log2 transformation with 0.5 pseudocount added to avoid infinite values
- Filtered primary tumors based on histology status and anatomical location

**Sample Retention:**
- IvyGAP: 122 primary GBM samples (from 270 total)
- TCGA: 142 primary GBM samples (from 548 total)
- CGGA: 183 primary GBM samples (from 1018 total)

#### Main Steps:

1. Load and preprocess IvyGAP, TCGA, and CGGA datasets
2. Filter for primary GBM samples based on histological and molecular criteria

#### Inputs:
   - IvyGap expression and clinical data files ([`etc/IvyGap/`](etc/IvyGap/))

#### Outputs:
   - RData files with preprocessed cohort data

---

## 3. Results

The [`Results/`](Results/) directory contains individual scripts and figures for each panel in the manuscript:

- **Main Figures**: Fig-1C, Fig-1D, Fig-1F, Fig-2A, Fig-2D, Fig-2I
- **Supplementary Figures**: Supplementary-Fig-1A, 3B, 3I, 4A, 4E, 5A, 5F

Each subdirectory contains:
- A dedicated R script for generating the figure
- A `Figures/` folder with output plots

## 4. Getting Started

### Prerequisites

- R version 4.4.2 or higher
- RStudio (recommended)
- renv & BiocManager packages for dependency management

### Installation

1. Clone this repository:
```bash
git clone https://github.com/yourusername/adapted-acidosis.git
cd adapted-acidosis
```

2. Open the R project file:
```r
# Open adapted-acidosis.Rproj in RStudio
```

3. Restore package dependencies:
```r
# renv should activate automatically
# If issues occur, run:
source("fix-renv.R")
```

### Usage

Run the analysis pipeline in order:

```r
# 1. Preprocessing and differential expression analysis
source("bin/1-preprocessing.R")

# 2. Gene set enrichment analysis
source("bin/2-GSEA-analysis.R")

# 3. Cluster enriched terms
source("bin/3-cluster-terms.R")

# 4. Validate gene signatures
source("bin/4-signature-analysis.R")
```

## 5. Data Availability

All raw and processed data have been deposited in NCBI's Gene Expression Omnibus (GEO):

- **CCLD**: GSE300758
- **U87**: GSE300765
- **PANC1**: GSE300768
- **HGCC**: GSE300771

*Note: Datasets will remain private until September 30, 2025, unless the manuscript enters review earlier.*

## 6. Session Information

For reproducibility, R session information is available in [`sessionInfo.txt`](sessionInfo.txt).

## 7. Citation

The data - missing from this repositroy - as well as the same code is also reposited on [Zenodo](https://doi.org/10.5281/zenodo.18414879). If you use this code or data, please cite our manuscript (details to be added upon publication).

## 8. License

See the [LICENSE](LICENSE) file for details.

## 9. Contact

For questions or issues, please open an issue on GitHub or contact [@MartonHorvath98](https://github.com/MartonHorvath98).
