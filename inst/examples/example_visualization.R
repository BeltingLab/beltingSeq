################################################################################
# Example: Creating Visualizations Using omicsFlow Package                     #
################################################################################
# This demonstrates creating publication-ready figures using omicsFlow

# Load the package
library(omicsFlow)

# Load additional packages
library(ggplot2)
library(cowplot)
library(openxlsx)

# Set working directory
wd <- getwd()

# Create output directory
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

################################################################################
# Load processed data                                                          #
################################################################################

load(file.path(wd, "RData", "CCLD_processedData.RData"))
load(file.path(wd, "RData", "CCLD_enrichment_results.RData"))
load(file.path(wd, "RData", "MSigDB_gene_sets.RData"))

################################################################################
# EXAMPLE 1: Volcano Plot                                                      #
################################################################################

# Create basic volcano plot (omicsFlow function)
volcano_basic <- plot_vulcan(CCLD.df, label = FALSE)

# Create volcano plot with gene labels (omicsFlow function)
volcano_labeled <- plot_vulcan(CCLD.df, label = TRUE)

# Save plots
ggsave("results/figures/volcano_basic.png", volcano_basic,
       width = 8, height = 6, dpi = 300, bg = "white")
ggsave("results/figures/volcano_labeled.svg", volcano_labeled,
       width = 8, height = 6, device = "svg", bg = "white")

cat("✓ Volcano plots created\n")

################################################################################
# EXAMPLE 2: GSEA Enrichment Plot                                              #
################################################################################

# Define pathways of interest
selected_pathways <- c(
  "EXTRACELLULAR_MATRIX_ORGANIZATION",
  "GLYCOSAMINOGLYCAN_METABOLISM",
  "PROTEOGLYCAN_METABOLIC_PROCESS",
  "CHONDROITIN_SULFATE_DERMATAN_SULFATE_METABOLISM",
  "FAT_CELL_DIFFERENTIATION"
)

# Define color palette
palette <- c(
  "PROTEOGLYCAN_METABOLIC_PROCESS" = "#E70041",
  "CHONDROITIN_SULFATE_DERMATAN_SULFATE_METABOLISM" = "#FF08FF",
  "GLYCOSAMINOGLYCAN_METABOLISM" = "#8D00FA",
  "EXTRACELLULAR_MATRIX_ORGANIZATION" = "#380186",
  "FAT_CELL_DIFFERENTIATION" = "#000000"
)

# Find pathway indices
pathway_ranks <- which(CCLD.GSEA$df$Name %in% names(palette))
go_ranks <- which(CCLD.GO$df$Name %in% names(palette))

# Extract running enrichment scores (omicsFlow function)
running_scores <- do.call(rbind, c(
  lapply(pathway_ranks, gs_info, object = CCLD.GSEA$gsea),
  lapply(go_ranks, gs_info, object = CCLD.GO$gsea)
))

# Clean up pathway names
running_scores$Description <- factor(
  gsub(pattern = "REACTOME_|GOBP_", replacement = "", 
       x = running_scores$Description),
  levels = names(palette)
)

# Create running score plot (omicsFlow function)
p1 <- plot_running_score(
  .df = running_scores, 
  .x = "x", .y = "runningScore", 
  .color = "Description", .palette = palette
)

# Create gene rank plot (omicsFlow function)
p2 <- plot_gene_rank(
  .df = running_scores, 
  .x = "x", .facet = "Description~.",
  .color = "Description", .palette = palette
)

# Combine plots
enrichment_plot <- cowplot::plot_grid(
  p1, p2 + theme(strip.text.y.left = element_blank()),
  byrow = TRUE, nrow = 2, ncol = 1, scale = 0.95, 
  rel_heights = c(1.2, 1), axis = "r"
)

# Save plot
ggsave("results/figures/gsea_enrichment.png", enrichment_plot, 
       width = 6, height = 4.45, dpi = 300, bg = "white")
ggsave("results/figures/gsea_enrichment.svg", enrichment_plot, 
       width = 6, height = 4.45, device = "svg", bg = "white")

cat("✓ GSEA enrichment plots created\n")

################################################################################
# EXAMPLE 3: Enrichment Summary Table                                          #
################################################################################

# Create enrichment table (omicsFlow function)
enrich_table <- get_enrichment_table(
  .df = rbind(CCLD.GSEA$df[pathway_ranks, ],
              CCLD.GO$df[go_ranks, ]),
  .order = "",
  .name = "Name"
)

# Add custom names
enrich_table <- enrich_table %>%
  dplyr::mutate(
    Pathway = c("ECM", "GAG METABOLISM", "ECM PGs", 
                "CS/DS METABOLISM", "LIPID DROPLET")
  ) %>%
  dplyr::relocate(Pathway, .before = everything())

# Save table
write.xlsx(enrich_table, rowNames = TRUE,
           file = "results/enrichment_summary_table.xlsx")

cat("✓ Enrichment summary table created\n")

################################################################################
# EXAMPLE 4: Multi-panel Figure                                                #
################################################################################

# Prepare panels
panel_A <- volcano_labeled + 
  ggtitle("A") + 
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0))

panel_B <- enrichment_plot + 
  ggtitle("B") + 
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0))

# Combine into multi-panel figure
figure_combined <- plot_grid(
  panel_A, panel_B,
  ncol = 2,
  rel_widths = c(1.2, 1),
  labels = NULL
)

# Save combined figure
ggsave("results/figures/figure_combined.png", figure_combined, 
       width = 14, height = 5.5, dpi = 300, bg = "white")
ggsave("results/figures/figure_combined.svg", figure_combined, 
       width = 14, height = 5.5, device = "svg", bg = "white")

cat("✓ Multi-panel figure created\n")

################################################################################
# Summary                                                                       #
################################################################################

cat("\n========================================\n")
cat("Visualization complete!\n")
cat("Using omicsFlow package functions:\n")
cat("  - plot_vulcan()\n")
cat("  - gs_info()\n")
cat("  - plot_running_score()\n")
cat("  - plot_gene_rank()\n")
cat("  - get_enrichment_table()\n")
cat("\nFigures saved in: results/figures/\n")
cat("Available formats: PNG (300 dpi) and SVG (vector)\n")
cat("========================================\n")
