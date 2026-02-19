################################################################################
# Example Workflow: Complete Analysis Pipeline Using beltingSeq                 #
################################################################################
# This script demonstrates a complete analysis workflow from raw data to
# publication-ready figures using the beltingSeq package

# Load required packages
library(beltingSeq)
library(oligo)
library(msigdbr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(openxlsx)

# Set working directory
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

################################################################################
# STEP 1: DATA PREPROCESSING                                                   #
################################################################################

# Read Affymetrix CEL files
data_dir <- "../data/raw/your_data/"
cel_data <- read.celfiles(list.celfiles(data_dir, full.names = TRUE))

# Normalize with RMA
library(clariomdhumantranscriptcluster.db)
normalized_data <- normalizeTranscript(cel_data, clariomdhumantranscriptcluster.db)

# Calculate log2 fold change (adjust sample names as needed)
normalized_data$log2FoldChange <- 
  normalized_data$Treatment - normalized_data$Control

# Remove duplicate probes
unique_data <- removeDuplicates(
  .data = normalized_data,
  .column = "log2FoldChange",
  .symbol = "SYMBOL"
)

################################################################################
# STEP 2: DIFFERENTIAL EXPRESSION ANALYSIS                                     #
################################################################################

# Define experimental design
design <- c("control", "control", "control",
            "treatment", "treatment", "treatment")

# Define contrasts
contrasts <- c("treatment-control")

# Perform differential expression analysis
deg_results <- limmaDEA(
  .data = unique_data,
  .design = design,
  .contrast = contrasts
)

# Process first contrast
deg_df <- deg_results[[1]] %>%
  dplyr::rename(
    Symbol = ID.Symbol,
    log2FoldChange = logFC,
    padj = adj.P.Val,
    pvalue = P.Value,
    entrezID = ID.entrezID
  ) %>%
  get_significance(fc_threshold = 0.5, padj_threshold = 0.05)

# Check results
table(deg_df$significance)

# Save results
dir.create("results", showWarnings = FALSE)
write.xlsx(deg_df, file = "results/differential_expression.xlsx")

################################################################################
# STEP 3: PREPARE GENE SETS                                                    #
################################################################################

# Load MSigDB gene sets
msigdbr_df <- msigdbr(species = "Homo sapiens")

# KEGG and Reactome pathways
pathways <- msigdbr_df %>%
  dplyr::filter(
    gs_cat == "C2",
    gs_subcat %in% c("CP:KEGG", "CP:REACTOME")
  )

# Hallmark and GO biological process
terms <- msigdbr_df %>%
  dplyr::filter(
    gs_cat == "H" | (gs_cat == "C5" & gs_subcat == "GO:BP")
  ) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    gs_subcat = ifelse(gs_cat == "H", "HALLMARK", gs_subcat),
    gs_exact_source = ifelse(gs_cat == "H", gs_id, gs_exact_source)
  ) %>%
  dplyr::ungroup()

# Save gene sets
save(msigdbr_df, pathways, terms, file = "results/gene_sets.RData")

################################################################################
# STEP 4: GENE SET ENRICHMENT ANALYSIS                                         #
################################################################################

# Extract gene lists
gene_list <- get_genelist(
  .df = deg_df,
  .filter = deg_df$significance %in% c("Signif. up-regulated", 
                                       "Signif. down-regulated"),
  .value = "log2FoldChange",
  .name = "entrezID"
)

# GSEA on GO terms
GO_results <- run_gsea(
  .geneset = gene_list$background,
  .terms = terms
)
GO_results <- c(GO_results, 
                extract_gsea_results(GO_results$gsea, terms))

# GSEA on pathways
PATHWAY_results <- run_gsea(
  .geneset = gene_list$background,
  .terms = pathways
)
PATHWAY_results <- c(PATHWAY_results,
                     extract_gsea_results(PATHWAY_results$gsea, pathways))

# Save enrichment results
save(GO_results, PATHWAY_results, 
     file = "results/enrichment_results.RData")

# Export significant results
write.xlsx(list(GO = GO_results$sig_df,
                Pathways = PATHWAY_results$sig_df),
           file = "results/significant_enrichments.xlsx")

################################################################################
# STEP 5: VISUALIZATION                                                        #
################################################################################

# Create output directory
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

# 1. Volcano plot
volcano_plot <- plot_vulcan(deg_df, label = TRUE)
ggsave("results/figures/volcano_plot.png", volcano_plot,
       width = 8, height = 6, dpi = 300, bg = "white")
ggsave("results/figures/volcano_plot.svg", volcano_plot,
       width = 8, height = 6, device = "svg", bg = "white")

# 2. GSEA enrichment plot for selected pathways
selected_pathways <- c(
  "HALLMARK_HYPOXIA",
  "EXTRACELLULAR_MATRIX_ORGANIZATION",
  "GLYCOSAMINOGLYCAN_METABOLISM"
)

palette <- c(
  "HALLMARK_HYPOXIA" = "#F00000",
  "EXTRACELLULAR_MATRIX_ORGANIZATION" = "#380186",
  "GLYCOSAMINOGLYCAN_METABOLISM" = "#8D00FA"
)

# Find pathway indices
pathway_ranks <- which(PATHWAY_results$df$Name %in% names(palette))
go_ranks <- which(GO_results$df$Name %in% names(palette))

# Extract running scores
running_scores <- do.call(rbind, c(
  lapply(pathway_ranks, gsInfo, object = PATHWAY_results$gsea),
  lapply(go_ranks, gsInfo, object = GO_results$gsea)
))

running_scores$Description <- factor(
  gsub(pattern = "REACTOME_|GOBP_|HALLMARK_", replacement = "", 
       x = running_scores$Description),
  levels = names(palette)
)

# Create plots
p1 <- plotRunningScore(
  .df = running_scores,
  .x = "x", .y = "runningScore",
  .color = "Description", .palette = palette
)

p2 <- plotGeneRank(
  .df = running_scores,
  .x = "x", .facet = "Description~.",
  .color = "Description", .palette = palette
)

enrichment_plot <- plot_grid(
  p1, p2 + theme(strip.text.y.left = element_blank()),
  byrow = TRUE, nrow = 2, ncol = 1, scale = 0.95,
  rel_heights = c(1.2, 1), axis = "r"
)

ggsave("results/figures/gsea_enrichment.png", enrichment_plot,
       width = 6, height = 4.45, dpi = 300, bg = "white")
ggsave("results/figures/gsea_enrichment.svg", enrichment_plot,
       width = 6, height = 4.45, device = "svg", bg = "white")

# 3. Create enrichment summary table
enrich_table <- getEnrichmentTable(
  .df = rbind(PATHWAY_results$df[pathway_ranks, ],
              GO_results$df[go_ranks, ]),
  .order = "",
  .name = "Name"
)
write.xlsx(enrich_table, rowNames = TRUE,
           file = "results/enrichment_summary.xlsx")

################################################################################
# STEP 6: SESSION INFO                                                         #
################################################################################

# Save session information
sink("results/sessionInfo.txt")
sessionInfo()
sink()

cat("\n========================================\n")
cat("Analysis complete!\n")
cat("Results saved in: results/\n")
cat("Figures saved in: results/figures/\n")
cat("========================================\n")
