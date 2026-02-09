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
  load(file = file.path(wd, "RData", "HGCC_enrichment_results.RData"))
  load(file = file.path(wd, "RData", "MSigDB_gene_sets.RData"))
} else {
  stop("Required .RData file(s) not found in the working directory.")
}
# Create Fig-1D: Enrichment vulcano plot of selected, shared gene sets in 3 primary cell lines
# Define paths
output_dir <- file.path(wd, "Results", "Fig-1D")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
# Vulcano visualization of the shared and selected enriched terms
vulcano_pathways <- read.csv("etc/pathways-of-interest.txt", header = T,
                             sep = "\t", stringsAsFactors = T)
vulcano_pathways <- vulcano_pathways %>% filter(Category == "ECM")
# Rename selected pathways 
rename_pathways <- c("ECM" = "EXTRACELLULAR_MATRIX_ORGANIZATION",
                     "GAG METABOLISM" = "GLYCOSAMINOGLYCAN_METABOLISM",
                     "ECM PGs" = "PROTEOGLYCAN_METABOLIC_PROCESS",
                     "ECM PGs" = "ECM_PROTEOGLYCANS",
                     "CS/DS METABOLISM" = "CHONDROITIN_SULFATE_DERMATAN_SULFATE_METABOLISM",
                     "LIPID DROPLET" = "FAT_CELL_DIFFERENTIATION")

# Combine GSEA and GO enrichment results for vulcano plots
FIG1.vulcano.object <- list(
  U3017 = rbind(HGCC.GSEA$U3017$sig_df[,c("ID","Name", "NES", "p.adjust")],
                HGCC.GO$U3017$sig_df[,c("ID","Name", "NES","p.adjust")]),
  U3047 = rbind(HGCC.GSEA$U3047$sig_df[,c("ID","Name", "NES","p.adjust")],
                HGCC.GO$U3047$sig_df[,c("ID","Name", "NES","p.adjust")]),
  U3054 = rbind(HGCC.GSEA$U3054$sig_df[,c("ID","Name", "NES","p.adjust")],
                HGCC.GO$U3054$sig_df[,c("ID","Name", "NES","p.adjust")]))

# Cap very low p-values
FIG1.vulcano.object <- lapply(FIG1.vulcano.object, function(x){
  x %>% mutate(p.adjust = ifelse(p.adjust < 1e-35, 1e-35, p.adjust))
})

# Add print names to dataframes
FIG1.vulcano.object <- lapply(FIG1.vulcano.object, function(x){
  x %>% 
    dplyr::rowwise(.) %>% 
    dplyr::mutate(`Print name` = ifelse(Name %in% rename_pathways,
                                       names(rename_pathways[which(rename_pathways == Name)]),
                                       NA))
})
FIG1.vulcano.object <- lapply(FIG1.vulcano.object, setNames,
                         c("ID", "Name", "NES", "FDR", "Print name"))
# Save results table
openxlsx::write.xlsx(FIG1.vulcano.object, sheetName = names(FIG1.vulcano.object),
                     file.path(output_dir, "Fig-1D-table.xlsx"))

# 1.) Extracellular matrix organization
FIG1.vulcano <- list()
(FIG1.vulcano$ECM <- lapply(FIG1.vulcano.object, function(x){
  plotClusters(.df = x, 
               .pathways = vulcano_pathways)
}))
# Combine the ECM plots into a single figure
(FIG1.vulcano$ECM$total <- cowplot::plot_grid(
  FIG1.vulcano$ECM$U3054 +
    theme(axis.title.x = element_blank()), NULL,
  FIG1.vulcano$ECM$U3047 +
    theme(axis.title = element_blank()), NULL,
  FIG1.vulcano$ECM$U3017 +
    theme(axis.title = element_blank()),
  rel_widths = c(1, 0.15, 1, 0.15, 1),
  byrow = T, nrow = 1, ncol = 5, scale = 1, 
  margins = c(0.5, 0.5, 0.5, 0.5)))

# Create output directory for figures
if (!dir.exists(file.path(output_dir, "Figures"))) {
  dir.create(file.path(output_dir, "Figures"), recursive = TRUE)
}
# Save the plot as svg for publication
ggsave(file.path(output_dir, "Figures", "Fig-1D-vulcanplot-vector.svg"),
       device = "svg", plot = FIG1.vulcano$ECM$total,
       bg = "white", width = 15, height = 5, units = "in")
# Save the plot as png for presentation
ggsave(file.path(output_dir, "Figures", "Fig-1D-vulcanplot-print.png"),
       device = "png", plot = FIG1.vulcano$ECM$total, 
       bg = "white", dpi = 300, width = 15, height = 5, units = "in")

# Clean up environment
rm(list = ls()[!(ls() %in% c("wd"))])
gc()
# End of script