# beltingSeq - Microarray and RNA-seq Analysis Package

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# **🎉 NEW: beltingSeq** 

**beltingSeq** provides a complete workflow for analyzing gene expression data, from raw data normalization to publication-ready visualizations. Originally developed for the adapted-acidosis project, it has been generalized for use with any microarray or RNA-seq dataset.

### Features

🧬 **Data Processing**: Microarray normalization (Affymetrix RMA, Illumina neqc), differential expression with limma, automatic gene annotation  
📊 **Enrichment Analysis**: GSEA with MSigDB, pathway clustering, ORA  
📈 **Visualizations**: Volcano plots, GSEA enrichment plots, network graphs, statistical plots

<details>

### Data Processing
- **Microarray normalization**: Affymetrix (RMA) and Illumina (neqc) platforms
- **Differential expression analysis**: Using limma with customizable contrasts
- **Gene annotation**: Automatic conversion between probe IDs, gene symbols, and Entrez IDs
- **Duplicate handling**: Smart removal of duplicate probes based on effect size

### Enrichment Analysis
- **Gene Set Enrichment Analysis (GSEA)**: Using MSigDB gene sets (KEGG, Reactome, GO, Hallmark)
- **Over-Representation Analysis (ORA)**: Hypergeometric testing for pathway enrichment
- **Term clustering**: Identify and visualize related pathways using Cohen's Kappa coefficient
- **Network analysis**: Louvain clustering of enriched gene sets

### Visualizations
- **Volcano plots**: Highlight differential expression with customizable significance thresholds
- **GSEA enrichment plots**: Running enrichment scores with gene rankings
- **Network graphs**: Visualize relationships between enriched pathways
- **Cluster plots**: Compare enrichment across multiple datasets
- **Statistical plots**: Boxplots with statistical testing

</details>

## Installation

```r
# Install from GitHub
devtools::install_github("BeltingLab/beltingSeq")

# Or install from local source
devtools::install("path/to/beltingSeq")
```

## Quick Start

```r
library(beltingSeq)

# 1. Normalize microarray data
normalized_data <- normalize_transcript(cel_files, clariomdhumantranscriptcluster.db)

# 2. Differential expression analysis
design <- c("control", "control", "treatment", "treatment")
contrasts <- "treatment-control"
deg_results <- limma_dea(normalized_data, design, contrasts)

# 3. Gene set enrichment analysis
gene_list <- get_genelist(deg_results[[1]], 
                          filter = deg_results[[1]]$adj.P.Val < 0.05,
                          value = "t", 
                          name = "entrezID")
gsea_results <- run_gsea(gene_list$background, pathways)

# 4. Visualize results
plot_vulcan(deg_results[[1]], label = TRUE)
```

## Directory Structure

```
beltingSeq/
├── R/                    # Function definitions
│   ├── processing.R      # Data processing functions
│   ├── enrichment.R      # GSEA and ORA functions
│   ├── clustering.R      # Term clustering functions
│   └── visualization.R   # Plotting functions
├── man/                  # Documentation (auto-generated)
├── vignettes/            # Usage examples
│   ├── preprocessing.Rmd
│   ├── gsea_analysis.Rmd
│   └── visualizations.Rmd
├── inst/
│   └── extdata/          # Example data
├── tests/                # Unit tests
├── DESCRIPTION           # Package metadata
├── NAMESPACE             # Exports
└── README.md             # This file
```

## Documentation

- **Vignettes**: [Preprocessing](vignettes/preprocessing.Rmd) | [GSEA](vignettes/gsea_analysis.Rmd) | [Visualizations](vignettes/visualizations.Rmd)
- **Examples**: [inst/examples/](inst/examples/)

## Example Workflows

See the `vignettes/` directory for detailed examples:
- **Preprocessing**: Normalize and filter microarray data
- **GSEA Analysis**: Complete enrichment analysis workflow
- **Visualizations**: Create publication-ready figures

## Citation

The originally developed code is reposited on [Zenodo](https://doi.org/10.5281/zenodo.18414879). If you use this code or data, please cite our manuscript: [doi.org/10.1038/s41556-026-01879-y](https://doi.org/10.1038/s41556-026-01879-y)

> Bång-Rudenstam, A., Cerezo-Magaña, M., Horvath, M. *et al.* Tumour acidosis remodels the glycocalyx to control lipid scavenging and ferroptosis. *Nat Cell Biol* (2026).

## License

See the [LICENSE](LICENSE) file for details.

## Contact

For questions and bug reports, please open an issue on GitHub or contact [@MartonHorvath98](https://github.com/MartonHorvath98).
