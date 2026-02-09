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
  load(file.path(wd, "RData", "HGCC_enrichment_results.RData"))
} else {
  stop("Required .RData file(s) not found in the working directory.")
}
# Create Supplementary Fig-5F: Enrichment plot of selected gene sets in MG cell lines (3D vs. 2D)
# Define paths
output_dir <- file.path(wd, "Results", "Supplementary-Fig-5F")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Combine enrichment objects
GO.object = list(
  U3047 = HGCC.GO$U3047$gsea,
  U3017 = HGCC.GO$U3017$gsea
)

# Visualization of selected enrichment: Hypoxia ("M5891") 
hypoxia.ranks = list(
  U3047 = which(HGCC.GO$U3047$df$ID == "M5891"),
  U3017 = which(HGCC.GO$U3017$df$ID == "M5891"))

SUPP5F.enrichplot <- list()
SUPP5F.enrichplot$scores  <- do.call(rbind, c(
  lapply(names(hypoxia.ranks)[1:2], function(i){
    gsInfo(hypoxia.ranks[[i]], object = GO.object[[i]]) %>% 
      dplyr::mutate(source = i)}
  )))

SUPP5F.enrichplot$scores <- SUPP5F.enrichplot$scores %>% 
  dplyr::mutate(
    Description = gsub(x = SUPP5F.enrichplot$scores$Description,
                       pattern = "KEGG_", replacement =  ""),
    source = factor(source, levels = names(GO.object)[1:2]))

SUPP5F.enrichplot$table <- getEnrichmentTable(
  .df = rbind(HGCC.GO$U3047$df[hypoxia.ranks$U3047,] %>% 
                dplyr::mutate(source = "U3047"),
              HGCC.GO$U3017$df[hypoxia.ranks$U3017,] %>% 
                dplyr::mutate(source = "U3017")),
  .order = c(1:2), .name = "source")  %>% 
  dplyr::mutate(
    `Print name` = c("U3047MG 3D vs. 2D", "U3017MG 3D vs. 2D")) %>%
  dplyr::relocate(`Print name`, .before = everything())

# Save results table
openxlsx::write.xlsx(SUPP5F.enrichplot$table, rowNames = TRUE,
                     file.path(output_dir, "Supplementary-Fig-5F-table.xlsx"))

(SUPP5F.enrichplot$plot  <- cowplot::plot_grid(
  plotRunningScore(.df = SUPP5F.enrichplot$scores, 
                   .x = "x", .y = "runningScore", 
                   .color = "source", .palette = c("U3047" = "salmon", 
                                                   "U3017" = "pink"
                   )),
  plotGeneRank(.df = SUPP5F.enrichplot$scores, 
               .x = "x", .facet = "source~.",
               .color = "source", .palette = c("U3047" = "salmon",
                                               "U3017" = "pink"
               )),
  byrow = T, nrow = 2, ncol = 1, scale = .95, 
  rel_heights = c(1.2,1), axis = "r",
  margins = c(0.5, 0.5, 0.5, 0.5)))

# Create output directory for figures
if (!dir.exists(file.path(output_dir, "Figures"))) {
  dir.create(file.path(output_dir, "Figures"), recursive = TRUE)
}
# Save as svg for publication
ggsave(file.path(output_dir, "Figures", "Supplementary-Fig-5F-enrichplot-vector.svg"),
       device = "svg", plot = SUPP5F.enrichplot$plot,
       bg = "white", width = 6, height = 4.45, units = "in")
# Save as png for presentation
ggsave(file.path(output_dir, "Figures", "Supplementary-Fig-5F-enrichplot-print.png"),
       device = "png", plot = SUPP5F.enrichplot$plot,
       bg = "white", dpi = 300, width = 6, height = 4.45, units = "in")

# Clean up environment
rm(list = ls()[!(ls() %in% c("wd"))])
gc()
# End of script