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
# Create Supplementary Fig-5A: Enrichment plot of selected gene sets in MG cell lines (3D vs. 2D)
# Define paths
output_dir <- file.path(wd, "Results", "Supplementary-Fig-5A")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Combine enrichment objects
GSEA.object = list(
  U3047 = HGCC.GSEA$U3047$gsea,
  U3017 = HGCC.GSEA$U3017$gsea
)

# Visualization of selected enrichment: TGF-B pathway ("hsa04350") 
tgfb.ranks = list(
  U3047 = which(HGCC.GSEA$U3047$df$ID == "hsa04350"),
  U3017 = which(HGCC.GSEA$U3017$df$ID == "hsa04350"))

SUPP5A.enrichplot <- list()
SUPP5A.enrichplot$scores  <- do.call(rbind, c(
  lapply(names(tgfb.ranks)[1:2], function(i){
    gsInfo(tgfb.ranks[[i]], object = GSEA.object[[i]]) %>% 
      dplyr::mutate(source = i)}
  )))

SUPP5A.enrichplot$scores <- SUPP5A.enrichplot$scores %>% 
  dplyr::mutate(
    Description = gsub(x = SUPP5A.enrichplot$scores$Description,
                       pattern = "KEGG_", replacement =  ""),
    source = factor(source, levels = names(GSEA.object)[1:2]))

SUPP5A.enrichplot$table <- getEnrichmentTable(
  .df = rbind(HGCC.GSEA$U3047$df[tgfb.ranks$U3047,] %>% 
                dplyr::mutate(source = "U3047"),
              HGCC.GSEA$U3017$df[tgfb.ranks$U3017,] %>% 
                dplyr::mutate(source = "U3017")),
  .order = c(1:2), .name = "source")  %>% 
  dplyr::mutate(
    `Print name` = c("U3047MG 3D vs. 2D", "U3017MG 3D vs. 2D")) %>%
  dplyr::relocate(`Print name`, .before = everything())

# Save results table
openxlsx::write.xlsx(SUPP5A.enrichplot$table, rowNames = TRUE,
                     file.path(output_dir, "Supplementary-Fig-5A-table.xlsx"))

(SUPP5A.enrichplot$plot  <- cowplot::plot_grid(
  plotRunningScore(.df = SUPP5A.enrichplot$scores, 
                   .x = "x", .y = "runningScore", 
                   .color = "source", .palette = c("U3047" = "salmon", 
                                                   "U3017" = "pink"
                   )),
  plotGeneRank(.df = SUPP5A.enrichplot$scores, 
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
ggsave(file.path(output_dir, "Figures", "Supplementary-Fig-5A-enrichplot-vector.svg"),
       device = "svg", plot = SUPP5A.enrichplot$plot,
       bg = "white", width = 6, height = 4.45, units = "in")
# Save as png for presentation
ggsave(file.path(output_dir, "Figures", "Supplementary-Fig-5A-enrichplot-print.png"),
       device = "png", plot = SUPP5A.enrichplot$plot,
       bg = "white", dpi = 300, width = 6, height = 4.45, units = "in")

# Clean up environment
rm(list = ls()[!(ls() %in% c("wd"))])
gc()
# End of script