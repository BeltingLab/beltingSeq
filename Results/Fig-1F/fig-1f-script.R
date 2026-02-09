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
  load(file.path(wd, "RData", "HGCC_processedData.RData"))
  load(file.path(wd, "RData", "CCLD_processedData.RData"))
} else {
  stop("Required .RData file(s) not found in the working directory.")
}
# Create Fig-1F: Heatmap of the expression of genes from the gene signature in 
# Patient GBM cells (LD+ vs. LD-) and primary GBM cells (2D vs. 3D)

# Define paths
output_dir <- file.path(wd, "Results", "Fig-1F")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
# load gene set of interest
interest_genes <- read.csv("etc/gene-signature.txt", header = T,
                           sep = "\t", stringsAsFactors = T)

FIG1.heatmap <- list()
FIG1.heatmap$data <- interest_genes %>% 
  dplyr::inner_join(., CCLD.df[,1:2], by = c("SYMBOL" = "Symbol")) %>%
  dplyr::inner_join(., HGCC.deg$U3017[,c(2,4)], by = c("SYMBOL" = "ID.Symbol")) %>% 
  dplyr::inner_join(., HGCC.deg$U3047[,c(2,4)], by = c("SYMBOL" = "ID.Symbol")) %>%
  dplyr::inner_join(., HGCC.deg$U3054[,c(2,4)], by = c("SYMBOL" = "ID.Symbol")) %>% 
  dplyr::rename_at(vars(contains("log")), ~c("Patient", "U3017MG", "U3047MG", "U3054MG")) %>%
  dplyr::mutate(Category = factor(Category, 
                                  levels = c("Hypoxia/Acidosis", "CSPG Core", "CS GAG Synthesis",
                                             "ECM Organization", "Lipid Metabolism"),
                                  labels = c("Hypoxia/Acidosis", "CSPG Core", "CS GAG Synthesis",
                                             "ECM Organization", "Lipid Metabolism"))) %>% 
  dplyr::select(SYMBOL, Category, U3017MG, U3047MG, U3054MG, Patient)

# Save results table
openxlsx::write.xlsx(FIG1.heatmap$data,
                     file.path(output_dir, "Fig-1F-table.xlsx"))

# Prepare matrix for heatmap
FIG1.heatmap$matrix <- as.matrix(FIG1.heatmap$data[,c(3:6)])
# Set row names
row.names(FIG1.heatmap$matrix) <- FIG1.heatmap$data$SYMBOL

# Create category annotation
cat_annot = rowAnnotation(
  Category = FIG1.heatmap$data$Category,
  col = list(Category = c("Hypoxia/Acidosis" = "darkgreen",
                          "CSPG Core" = "turquoise",
                          "CS GAG Synthesis" = "#FF08FF",
                          "ECM Organization" = "#8D00FA",
                          "Lipid Metabolism" = "white")
             
  ),
  gp = gpar(col = "black", fontsize = 14, family = "arial"),
  show_annotation_name = F,
  annotation_legend_param = list(
    border = "black",
    annotation_label_rot = 90,    # rotate the category names beside rows
    annotation_label_side = "top" # align them properly               # display in one horizontal row)
))

# Define color scale for heatmap
matrix_col = colorRamp2(c(-4, 0, 4), c("blue", "white", "red"))

# Create heatmap
(FIG1.heatmap$plot = Heatmap(FIG1.heatmap$matrix,
                              
                              ### row settings
                              cluster_rows = F,
                              row_gap = unit(5, "mm"),
                              row_names_side = "left",
                              row_names_rot = -20,
                              row_names_gp = gpar(fontsize = 14,
                                                  family = "arial"),
                              
                              ### column settings
                              cluster_columns = F,
                              column_names_rot = -90,
                              column_names_side = "bottom",
                              column_names_gp = gpar(fontsize = 14, 
                                                     family = "arial"),
                              column_split = factor(c(rep("Primary cells 2D vs. 3D", 3), "GBM LD+ vs. LD-"),
                                                    levels = c("Primary cells 2D vs. 3D", "GBM LD+ vs. LD-")),
                              column_gap = unit(1, "lines"),
                              
                              ### legend settings
                              col = matrix_col, # color settings
                              right_annotation = cat_annot, 
                              heatmap_legend_param = list(
                                title = "Log2 (Fold change)",
                                direction = "vertical"),
                              border_gp = gpar(col = "black", lty = 2)))

# Create output directory for figures
if (!dir.exists(file.path(output_dir, "Figures"))) {
  dir.create(file.path(output_dir, "Figures"), recursive = TRUE)
}
# Save plot as svg for publication
svg(file.path(output_dir, "Figures", "Fig-1F-signature-heatmap.svg"),
    width = 12, height = 8)
draw(FIG1.heatmap$plot , merge_legend = T,
     heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()
# Save plot as png for presentation
png(file.path(output_dir, "Figures", "Fig-1F-signature-heatmap.png"),
    width = 12, height = 8, res = 300, units = "in", bg = "white")
draw(FIG1.heatmap$plot , merge_legend = T,
     heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

# Clean up environment
rm(list = ls()[!(ls() %in% c("wd"))])
gc()
# End of script