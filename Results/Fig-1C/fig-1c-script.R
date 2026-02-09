################################################################################
# Created: 2025 11 07 ; Last Modified: 2025 11 07 ; MH                         #
################################################################################
# Set working directory
wd <- getwd()
# Load packages
if (file.exists(file.path(wd, "bin", "packages.R"))) {
  source(file = file.path(wd, "bin", "packages.R"))
} else {
  stop("Required file 'packages.R' not found in the working directory.")
}
# Source data processing functions
if (file.exists(file.path(wd, "bin", "functions.R"))) {
  source(file = file.path(wd, "bin", "functions.R"))
} else {
  stop("Required file 'functions.R' not found in the working directory.")
}
# Load necessary RData files
if (file.exists(file.path(wd, "RData"))) {
  load(file = file.path(wd, "RData", "CCLD_enrichment_results.RData"))
  load(file = file.path(wd, "RData", "MSigDB_gene_sets.RData"))
} else {
  stop("Required .RData file(s) not found in the working directory.")
}
# Create Fig-1C: Enrichment plot of selected gene sets in CCLD
# Define paths
output_dir <- file.path(wd, "Results", "Fig-1C")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
# Define selected pathways and their colors 
selected_pathways <- c("ECM" = "EXTRACELLULAR_MATRIX_ORGANIZATION",
                       "GAG METABOLISM" = "GLYCOSAMINOGLYCAN_METABOLISM",
                       "ECM PGs" = "PROTEOGLYCAN_METABOLIC_PROCESS",
                       "CS/DS METABOLISM" = "CHONDROITIN_SULFATE_DERMATAN_SULFATE_METABOLISM",
                       "LIPID DROPLET" = "FAT_CELL_DIFFERENTIATION")
palette = c("PROTEOGLYCAN_METABOLIC_PROCESS" = "#E70041",
            "CHONDROITIN_SULFATE_DERMATAN_SULFATE_METABOLISM" = "#FF08FF",
            "GLYCOSAMINOGLYCAN_METABOLISM" = "#8D00FA",
            "EXTRACELLULAR_MATRIX_ORGANIZATION" = "#380186",
            "FAT_CELL_DIFFERENTIATION" = "#000000")

# Identify indices of selected pathways/GO terms in GSEA and GO results
pathway.ranks = which(CCLD.GSEA$df$Name %in% names(palette))
go.ranks = which(CCLD.GO$df$Name %in% names(palette))
# Extract running enrichment scores for selected gene sets
FIG1.enrichplot <- list()
FIG1.enrichplot$runningScores <- do.call(rbind, c(lapply(pathway.ranks, gsInfo, object = CCLD.GSEA$gsea),
                                                  lapply(go.ranks, gsInfo, object = CCLD.GO$gsea)))
# Clean up pathway/GO term names for plotting
FIG1.enrichplot$runningScores$Description <- factor(
  gsub(x = FIG1.enrichplot$runningScores$Description, pattern = "REACTOME_|GOBP_", replacement =  ""),
  levels = names(palette))

# Prepare summary table of enrichment results for selected gene sets
FIG1.enrichplot$table <- getEnrichmentTable(.df = rbind(CCLD.GSEA$df[pathway.ranks,],
                                                        CCLD.GO$df[go.ranks,]),
                                            .order = "",
                                            .name = "Name") %>% 
  dplyr::arrange(order(match(row.names(.), names(palette)))) %>% 
  dplyr::mutate(
    `Print name` = names(selected_pathways)[match(names(palette), selected_pathways)]) %>%
  dplyr::relocate(`Print name`, .before = everything())

# Save results table
openxlsx::write.xlsx(FIG1.enrichplot$table, rowNames = TRUE,
                     file.path(output_dir, "Fig-1C-table.xlsx"))


# Plot running enrichment scores for selected gene sets
p1 <- plotRunningScore(.df = FIG1.enrichplot$runningScores, 
                       .x = "x", .y = "runningScore", 
                       .color = "Description", .palette = palette)

# Plot gene ranks for selected gene sets, faceted by pathway/GO term
p2 <- plotGeneRank(.df = FIG1.enrichplot$runningScores, 
                   .x = "x", .facet = "Description~.",
                   .color = "Description", .palette = palette)

# Combine running score and gene rank plots into a single figure
(FIG1.enrichplot$plot <- cowplot::plot_grid(
  p1, p2 + theme(strip.text.y.left = element_blank()),
  byrow = T, nrow = 2, ncol = 1, scale = .95, 
  rel_heights = c(1.2,1), axis = "r",
  margins = c(0.5, 0.5, 0.5, 0.5)))

# Create output directory for figures
if (!dir.exists(file.path(output_dir, "Figures"))) {
  dir.create(file.path(output_dir, "Figures"), recursive = TRUE)
}
# Save plot as svg for publication
ggsave(file.path(output_dir, "Figures", "Fig-1C-enrichplot-vector.svg"), 
       device = "svg", plot = FIG1.enrichplot$plot, 
       bg = "white", width = 6, height = 4.45, units = "in")
# Save plot as png for presentation
ggsave(file.path(output_dir, "Figures", "Fig-1C-enrichplot-print.png"), 
       device = "png", plot = FIG1.enrichplot$plot, 
       bg = "white", dpi = 300, width = 6, height = 4.45, units = "in")

# Clean up environment
rm(list = ls()[!(ls() %in% c("wd"))])
gc()
# End of script