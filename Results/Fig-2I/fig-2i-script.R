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
  load(file.path(wd, "RData", "CCLD_enrichment_results.RData"))
  load(file.path(wd, "RData", "U87_enrichment_results.RData"))
} else {
  stop("Required .RData file(s) not found in the working directory.")
}
# Create Fig-2I: Enrichment plot of selected gene sets in CCLD, U3054MG, and U87MG cells
# Define paths
output_dir <- file.path(wd, "Results", "Fig-2I")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Combine enrichment objects
GO.object = list(
  CCLD = CCLD.GO$gsea,
  U3054 = HGCC.GO$U3054$gsea,
  U87_acu = U87.GO$`acu_pH64-control_acu`$gsea,
  U87_sel = U87.GO$`sel_pH647-control_sel`$gsea
)

# Visualization of selected enrichment: Hypoxia ("M5891") 
hypoxia.ranks = list(
  CCLD = which(CCLD.GO$df$ID == "M5891"),
  U3054 = which(HGCC.GO$U3054$df$ID == "M5891"),
  U87_acu = which(U87.GO$`acu_pH64-control_acu`$df$ID == "M5891"),
  U87_sel = which(U87.GO$`sel_pH647-control_sel`$df$ID == "M5891"))

FIG2I.enrichplot <- list()
FIG2I.enrichplot$GBM$scores  <- do.call(rbind, c(
  lapply(names(hypoxia.ranks)[1:2], function(i){
    gsInfo(hypoxia.ranks[[i]], object = GO.object[[i]]) %>% 
      dplyr::mutate(source = i)}
  )))

FIG2I.enrichplot$GBM$scores <- FIG2I.enrichplot$GBM$scores %>% 
  dplyr::mutate(
    Description = gsub(x = FIG2I.enrichplot$GBM$scores$Description,
                       pattern = "KEGG_", replacement =  ""),
    source = factor(source, levels = names(GO.object)[1:2]))

FIG2I.enrichplot$GBM$table <- getEnrichmentTable(
  .df = rbind(CCLD.GO$df[hypoxia.ranks$CCLD,] %>% 
                dplyr::mutate(source = "CCLD"),
              HGCC.GO$U3054$df[hypoxia.ranks$U3054,] %>% 
                dplyr::mutate(source = "U3054")),
  .order = c(1:2), .name = "source")  %>% 
  dplyr::mutate(
    `Print name` = c("GBM LD+ vs. LD-", "U3054MG 3D vs. 2D")) %>%
  dplyr::relocate(`Print name`, .before = everything())

(FIG2I.enrichplot$GBM$plot  <- cowplot::plot_grid(
  plotRunningScore(.df = FIG2I.enrichplot$GBM$scores, 
                   .x = "x", .y = "runningScore", 
                   .color = "source", .palette = c("U3054" = "darkred", 
                                                   "CCLD" = "steelblue"
                   )),
  plotGeneRank(.df = FIG2I.enrichplot$GBM$scores, 
               .x = "x", .facet = "source~.",
               .color = "source", .palette = c("U3054" = "darkred",
                                               "CCLD" = "steelblue"
               )),
  byrow = T, nrow = 2, ncol = 1, scale = .95, 
  rel_heights = c(1.2,1), axis = "r",
  margins = c(0.5, 0.5, 0.5, 0.5)))


FIG2I.enrichplot$U87$scores <- do.call(rbind, c(
  lapply(names(hypoxia.ranks)[3:4], function(i){
    gsInfo(hypoxia.ranks[[i]], object = GO.object[[i]]) %>% 
      dplyr::mutate(source = i)}
  )))
FIG2I.enrichplot$U87$scores <- FIG2I.enrichplot$U87$scores %>% 
  dplyr::mutate(
    Description = gsub(x = FIG2I.enrichplot$U87$scores$Description,
                       pattern = "KEGG_", replacement =  ""),
    source = factor(source, levels = names(GO.object)[c(3,4)]))

FIG2I.enrichplot$U87$table <- getEnrichmentTable(
  .df = rbind(U87.GO$`acu_pH64-control_acu`$df[hypoxia.ranks$U87_acu,] %>% 
                dplyr::mutate(source = "U87_acu"),
              U87.GO$`sel_pH647-control_sel`$df[hypoxia.ranks$U87_sel,] %>% 
                dplyr::mutate(source = "U87_sel")),
  .order = c(1:2), .name = "source") %>% 
  dplyr::mutate(
    `Print name` = c("Short-term acidosis", "AA vs. NA")) %>%
  dplyr::relocate(`Print name`, .before = everything())

(FIG2I.enrichplot$U87$plot <- cowplot::plot_grid(
  plotRunningScore(.df = FIG2I.enrichplot$U87$scores, 
                   .x = "x", .y = "runningScore", 
                   .color = "source", .palette = c("U87_sel" = "#F00000",
                                                   "U87_acu" = "salmon")),
  plotGeneRank(.df = FIG2I.enrichplot$U87$scores, 
               .x = "x", .facet = "source~.",
               .color = "source", .palette = c("U87_sel" = "#F00000",
                                               "U87_acu" = "salmon")),
  byrow = T, nrow = 2, ncol = 1, scale = .95, 
  rel_heights = c(1.2,1), axis = "r",
  margins = c(0.5, 0.5, 0.5, 0.5)))

(FIG2I.enrichplot$TOTAL_plot <- cowplot::plot_grid(
  FIG2I.enrichplot$GBM$plot,
  FIG2I.enrichplot$U87$plot,
  byrow = T, nrow = 2, ncol = 1, scale = .95,
  margins = c(0.5, 0.5, 0.5, 0.5)
))

# Save results table
openxlsx::write.xlsx(list(FIG2I.enrichplot$GBM$table,FIG2I.enrichplot$U87$table), 
                     sheetName = c("Upper panel", "Lower panel"), rowNames = TRUE,
                     file.path(output_dir, "Fig-2I-table.xlsx"))

# Create output directory for figures
if (!dir.exists(file.path(output_dir, "Figures"))) {
  dir.create(file.path(output_dir, "Figures"), recursive = TRUE)
}
# Save as svg for publication
ggsave(file.path(output_dir, "Figures", "Fig-2I-enrichplot-vector.svg"),
       device = "svg", plot = FIG2I.enrichplot$TOTAL_plot,
       bg = "white", width = 6, height = 9, units = "in")
# Save as png for presentation
ggsave(file.path(output_dir, "Figures", "Fig-2I-enrichplot-print.png"),
       device = "png", plot = FIG2I.enrichplot$TOTAL_plot,
       bg = "white", dpi = 300, width = 6, height = 9, units = "in")

# Clean up environment
rm(list = ls()[!(ls() %in% c("wd"))])
gc()
# End of script