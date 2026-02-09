################################################################################
# Created: 2025 11 09 ; Last Modified: 2025 11 09 ; MH                         #
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
  load(file = file.path(wd, "RData", "U87_processedData.RData"))
  load(file = file.path(wd, "RData", "U87_enrichment_results.RData"))
  load(file = file.path(wd, "RData", "term_similarity_matrix.RData"))
} else {
  stop("Required .RData file(s) not found in the working directory.")
}
# Create Supp. Fig-4A: GSEA network plot for U87 cells acute acidosis
# Define paths
output_dir <- file.path(wd, "Results", "Supplementary-Fig-4A")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Print name of clusters
cluster_names <- read.csv("etc/cluster-names.txt", header = T,
                          sep = "\t", stringsAsFactors = T) %>%
  dplyr::filter(Condition == "U87acu") %>% 
  dplyr::select(-Condition)

# Combine GO and pathway enrichment results
SUPP4A.combinedGSEA <- rbind(
  dplyr::filter(U87.GSEA$`acu_pH64-control_acu`$sig_df, p.adjust < 0.001),
  dplyr::filter(U87.GO$`acu_pH64-control_acu`$sig_df, p.adjust < 0.001))
# Preprocess combined GSEA results
SUPP4A.combinedGSEA <- SUPP4A.combinedGSEA %>%
  dplyr::rowwise(.) %>% 
  # Calculate background ratio
  dplyr::mutate(
    geneRatio = length(unlist(strsplit(core_enrichment,"\\/")))/setSize
  ) %>% 
  dplyr::ungroup() %>% 
  dplyr::relocate(c("ID", "Name", "setSize", "geneRatio", "NES", "p.adjust", "core_enrichment"),
                  .before = everything())
# Get the clusters
SUPP4A.cluster <- get_cluster(SUPP4A.combinedGSEA, similarity.matrix, .threshold = 0.25)

# Get representative terms
SUPP4A.cluster$df <- get_cluster_representative(.cluster = SUPP4A.cluster$df,
                                                .degs = U87.deg$`acu_pH64-control_acu`)
V(SUPP4A.cluster$graph)$Representative <- SUPP4A.cluster$df$Representative
V(SUPP4A.cluster$graph)$Description <- SUPP4A.cluster$df$Name

# Filter out clusters with less than 5 members
SUPP4A.cluster$sub_graph <- filter_graph(SUPP4A.cluster$graph, 5)
SUPP4A.cluster$layout <- layout_nicely(SUPP4A.cluster$sub_graph)

#Filter data frame to match the clustered graph
SUPP4A.cluster$df <- SUPP4A.cluster$df %>%
  dplyr::filter(ID %in% V(SUPP4A.cluster$sub_graph)$name)

# Identify clusters in gsea_df that match any representative term
cluster_mapping <- SUPP4A.cluster$df %>% 
  dplyr::select(ID, cluster) %>%
  dplyr::inner_join(cluster_names, by = c("ID"= "Cluster_representative")) %>% 
  dplyr::select(-ID)

# Merge the mapping back into gsea_df
SUPP4A.cluster$df <- SUPP4A.cluster$df %>%
  left_join(cluster_mapping, by = "cluster")

# Save the results
write.xlsx(SUPP4A.cluster$df,
           file = file.path(output_dir, "Supplementary-Fig-4A-table.xlsx"))

SUPP4A.cluster$plot <- plot_network(.net = SUPP4A.cluster$sub_graph,
                                    .layout = SUPP4A.cluster$layout,
                                    .df = SUPP4A.cluster$df,
                                    .condition = "ACIDOSIS")

# Create output directory for figures
if (!dir.exists(file.path(output_dir, "Figures"))) {
  dir.create(file.path(output_dir, "Figures"), recursive = TRUE)
}
# Save plot as svg for publication
ggsave(file.path(output_dir, "Figures", "Supplementary-Fig-4A-GSEA-network-vector.svg"), 
       device = "svg", plot = SUPP4A.cluster$plot, 
       bg = "white", width = 20, height = 14, units = "in")
# Save plot as png for presentation
ggsave(file.path(output_dir, "Figures", "Supplementary-Fig-4A-GSEA-network-print.png"), 
       device = "png", plot = SUPP4A.cluster$plot, 
       bg = "white", dpi = 300, width = 20, height = 14, units = "in")

# Clean up environment
rm(list = ls()[!(ls() %in% c("wd"))])
gc()
# End of script