################################################################################
# Example: Preprocessing Using beltingSeq Package                               #
################################################################################
# Load the package
library(beltingSeq)

# Load additional required packages
library(oligo)
library(dplyr)
library(openxlsx)
library(org.Hs.eg.db)
library(AnnotationDbi)

# Set working directory
wd <- getwd()

# Create output directories
dir.create("./RData", showWarnings = FALSE)
dir.create("./etc/input-files/expression_matrices/", 
           showWarnings = FALSE, recursive = TRUE)

################################################################################
# EXAMPLE 1: Affymetrix Clariom D Array                                        #
################################################################################

# Read CEL files
data_dir <- "../data/raw/Affymetrix/CCLD/"
CCLD.dat <- read.celfiles(list.celfiles(data_dir, full.names = TRUE))

# Normalize transcripts using RMA (beltingSeq function)
library(clariomdhumantranscriptcluster.db)
CCLD.expr <- normalize_transcript(CCLD.dat, clariomdhumantranscriptcluster.db)

# Calculate log2 fold change
CCLD.expr$log2FoldChange <- 
  CCLD.expr$`CC_LD_(Clariom_D_Human).CEL` - 
  CCLD.expr$`CC_noLD_(Clariom_D_Human).CEL`

# Remove duplicate symbols (beltingSeq function)
CCLD.expr <- remove_duplicates(
  .data = CCLD.expr,
  .column = "log2FoldChange",
  .symbol = "SYMBOL"
)

# Check results
length(which(CCLD.expr$log2FoldChange >= 0.5))   # Up-regulated
length(which(CCLD.expr$log2FoldChange <= -0.5))  # Down-regulated

# Prepare final data frame
CCLD.df <- CCLD.expr[, c(5, 7)] %>% 
  setNames(., c("Symbol", "log2FoldChange")) %>%
  dplyr::mutate(
    entrezID = mapIds(org.Hs.eg.db, Symbol, 
                     keytype = "SYMBOL", column = "ENTREZID"),
    significance = dplyr::case_when(
      log2FoldChange < -0.5 ~ 'Signif. down-regulated',
      log2FoldChange > 0.5 ~ 'Signif. up-regulated',
      TRUE ~ 'NS'
    )
  ) %>% 
  dplyr::filter(complete.cases(.))

# Save results
save(CCLD.expr, CCLD.df,
     file = file.path(wd, "RData", "CCLD_processedData.RData"))
write.xlsx(CCLD.expr, 
          file = "etc/input-files/expression_matrices/CCLD_processedData.xlsx")

################################################################################
# EXAMPLE 2: Differential Expression with Limma                                #
################################################################################

# Read new dataset
data_dir <- "../data/raw/Affymetrix/HGCC/"
HGCC.dat <- read.celfiles(list.celfiles(data_dir, full.names = TRUE))

# Normalize with beltingSeq
HGCC.expr <- normalize_transcript(HGCC.dat, clariomdhumantranscriptcluster.db)

# Perform differential expression analysis (beltingSeq function)
HGCC.deg <- limma_dea(
  .data = HGCC.expr,
  .design = c("U3017_2D", "U3017_2D", "U3017_2D",
              "U3017_3D", "U3017_3D", "U3017_3D",
              "U3047_2D", "U3047_2D", "U3047_2D",
              "U3047_3D", "U3047_3D", "U3047_3D",
              "U3054_2D", "U3054_2D", "U3054_2D",
              "U3054_3D", "U3054_3D", "U3054_3D"),
  .contrast = c("U3017_3D-U3017_2D",
                "U3047_3D-U3047_2D",
                "U3054_3D-U3054_2D")
)

names(HGCC.deg) <- c("U3017", "U3047", "U3054")

# Remove duplicates (beltingSeq function)
HGCC.deg <- lapply(HGCC.deg, remove_duplicates, 
                   .column = "t", .symbol = "ID.Symbol")

# Add significance classification (beltingSeq function)
HGCC.deg <- lapply(HGCC.deg, function(x){
  x %>% 
    dplyr::rename(
      PROBEID = "ID.ID",
      Symbol = "ID.Symbol",
      entrezID = "ID.entrezID",
      log2FoldChange = "logFC",
      padj = "adj.P.Val",
      pvalue = "P.Value"
    ) %>%
    get_significance(fc_threshold = 0.5, padj_threshold = 0.05)
})

# Save results
save(HGCC.expr, HGCC.deg,
     file = file.path(wd, "RData", "HGCC_processedData.RData"))
write.xlsx(HGCC.deg, sheetName = names(HGCC.deg),
          file = "etc/input-files/expression_matrices/HGCC_processedData.xlsx")

################################################################################
# Summary                                                                       #
################################################################################

cat("\n========================================\n")
cat("Preprocessing complete!\n")
cat("Using beltingSeq package functions:\n")
cat("  - normalize_transcript()\n")
cat("  - remove_duplicates()\n")
cat("  - limma_dea()\n")
cat("  - get_significance()\n")
cat("========================================\n")

# Save session info
sink("./sessionInfo.txt")
sessionInfo()
sink()
