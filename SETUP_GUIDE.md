# beltingSeq Package Setup and Installation Guide

## Quick Start

### 1. Install the Package

From within the project directory:

```r
# Install required tools
install.packages("devtools")
install.packages("roxygen2")

# Install the beltingSeq package
devtools::install(".", dependencies = TRUE)
```

### 2. Load and Test

```r
library(beltingSeq)

# Check available functions
ls("package:beltingSeq")

# View function documentation
?normalizeTranscript
?run_gsea
?plot_vulcan
```

### 3. Run Example Workflows

```r
# Browse vignettes
browseVignettes("beltingSeq")

# Run complete example workflow
source(system.file("examples", "complete_workflow.R", package = "beltingSeq"))
```

## Installation Options

### Option A: Install as Package (Recommended)

Best for using functions across multiple projects.

```r
devtools::install("path/to/adapted-acidosis")
library(beltingSeq)
```

### Option B: Load Without Installing

Best for development or testing.

```r
devtools::load_all("path/to/adapted-acidosis")
```

### Option C: Source Individual Files

Best for selective use.

```r
source("R/processing.R")
source("R/enrichment.R")
source("R/visualization.R")
```

## Building Documentation

### Generate function documentation:

```r
# Install roxygen2 if not already installed
install.packages("roxygen2")

# Generate documentation
devtools::document()
```

This creates `.Rd` files in the `man/` directory from roxygen comments.

### Build vignettes:

```r
# Install vignette dependencies
install.packages(c("knitr", "rmarkdown"))

# Build vignettes
devtools::build_vignettes()
```

### View documentation:

```r
# After installation
help(package = "beltingSeq")
?normalizeTranscript
browseVignettes("beltingSeq")
```

## Package Structure

```
beltingSeq/
├── R/                  # Function source code
├── man/                # Documentation (auto-generated)
├── vignettes/          # Tutorials
├── inst/examples/      # Example scripts
├── DESCRIPTION         # Package metadata
├── NAMESPACE           # Exported functions
└── .Rbuildignore       # Build exclusions
```

## Checking Package

Before sharing or publishing:

```r
# Check package
devtools::check()

# Build package
devtools::build()
```

## Dependencies

### Required packages (installed automatically):

**Bioconductor packages:**
- limma, affycoretools, oligo
- clusterProfiler, DOSE, enrichplot
- org.Hs.eg.db, AnnotationDbi
- msigdbr, GSVA, rrvgo

**CRAN packages:**
- ggplot2, ggrepel, ggpubr, cowplot
- dplyr, tidyr, stringr
- igraph, openxlsx

### Installing Bioconductor packages manually:

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("limma", "clusterProfiler", 
                       "org.Hs.eg.db", "oligo"))
```

## Usage Examples

### 1. Data Preprocessing

```r
library(beltingSeq)
library(oligo)

# Read CEL files
cel_data <- read.celfiles(list.celfiles("data/", full.names = TRUE))

# Normalize
normalized <- normalizeTranscript(cel_data, clariomdhumantranscriptcluster.db)

# Differential expression
design <- c("control", "control", "treatment", "treatment")
deg <- limmaDEA(normalized, design, "treatment-control")
```

### 2. GSEA Analysis

```r
library(msigdbr)

# Get gene sets
pathways <- msigdbr(species = "Homo sapiens") %>%
  filter(gs_cat == "C2", gs_subcat == "CP:KEGG")

# Prepare gene list
genes <- get_genelist(deg[[1]], 
                      deg[[1]]$adj.P.Val < 0.05,
                      "logFC", "ENTREZID")

# Run GSEA
results <- run_gsea(genes$background, pathways)
results <- c(results, extract_gsea_results(results$gsea, pathways))
```

### 3. Visualization

```r
# Volcano plot
plot_vulcan(deg[[1]], label = TRUE)

# GSEA enrichment
enrich_data <- gsInfo(results$gsea, 1)
plotRunningScore(enrich_data, "x", "runningScore", 
                 "Description", c("pathway" = "red"))
```

## Updating the Package

After modifying functions:

```r
# Re-document
devtools::document()

# Re-install
devtools::install(".", upgrade = FALSE)

# Or just reload for testing
devtools::load_all()
```

## Troubleshooting

### Issue: Package won't install

**Solution:**
```r
# Update R and packages
update.packages(ask = FALSE)

# Install dependencies manually
BiocManager::install(c("limma", "clusterProfiler"))

# Try again
devtools::install(".", dependencies = TRUE)
```

### Issue: Functions not found

**Solution:**
```r
# Reload package
detach("package:beltingSeq", unload = TRUE)
library(beltingSeq)

# Or reload for development
devtools::load_all()
```

### Issue: Vignettes not building

**Solution:**
```r
# Install knitr and rmarkdown
install.packages(c("knitr", "rmarkdown"))

# Build vignettes separately
devtools::build_vignettes()
```

## Next Steps

1. **Review vignettes:**
   - `vignettes/preprocessing.Rmd`
   - `vignettes/gsea_analysis.Rmd`
   - `vignettes/visualizations.Rmd`

2. **Try examples:**
   - `inst/examples/complete_workflow.R`
   - `inst/examples/example_preprocessing.R`
   - `inst/examples/example_gsea.R`

3. **Customize for your data:**
   - Copy example workflow
   - Update file paths and sample names
   - Run on your data

4. **Add new functions:**
   - Add to appropriate R/ file
   - Include roxygen documentation
   - Run `devtools::document()`
   - Re-install package

## Support

- Function help: `?function_name`
- Package help: `help(package = "beltingSeq")`
- Vignettes: `browseVignettes("beltingSeq")`
- Examples: `inst/examples/`
- Guide: `REORGANIZATION_GUIDE.md`

## Contributing

To add new functions:

1. Edit appropriate file in `R/`
2. Add roxygen documentation (use existing functions as template)
3. Update `NAMESPACE` if needed
4. Run `devtools::document()`
5. Test with `devtools::check()`
6. Install with `devtools::install()`

## License

MIT License - See LICENSE file for details
