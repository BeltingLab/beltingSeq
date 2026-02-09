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
  load(file = file.path(wd, "RData", "U87_processedData.RData"))
} else {
  stop("Required .RData file(s) not found in the working directory.")
}
# Create Fig-2A: Vulcan plot of U87MG cells with selected genes of interest
# Define paths
output_dir <- file.path(wd, "Results", "Fig-2A")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Gene sets of interest
interest_genes <- read.csv("etc/genes-of-interest.txt", header = T,
                                      sep = "\t", stringsAsFactors = T)
interest_genes <- interest_genes$Gene.name

# Limit the fold change values to a maximum of 6 and minimum of -6 for visualization
FIG2A.df <- U87.deg$`sel_pH647-control_sel` %>% 
  dplyr::rename(
    Symbol = ID.Symbol, 
    log2FoldChange = logFC, 
    padj = adj.P.Val, 
    pvalue = P.Value
  ) %>%
  dplyr::mutate(
    log2FoldChange = ifelse(log2FoldChange > 6, 6,
                            ifelse(log2FoldChange < -6, -6, log2FoldChange))
  ) %>% get_significance(.) %>% 
  dplyr::mutate(
    Interest = ifelse(Symbol %in% interest_genes, TRUE, FALSE)
  )

# Save results table
openxlsx::write.xlsx(FIG2A.df, file.path(output_dir, "Fig-2A-table.xlsx"))


# Create the vulcano plot of selected genes of interest
(FIG2A.vulcan = plot_vulcan(FIG2A.df, label = F) +
    geom_point(data = subset(FIG2A.df,
                             Symbol %in% interest_genes),
               aes(x = log2FoldChange, y = -log10(padj)),
               shape = 21, color = "black", fill = "black", size = 3, alpha = 0.8) +
    geom_text_repel(data = subset(FIG2A.df,
                                     Symbol %in% interest_genes),
                    aes(x = log2FoldChange, y = -log10(padj), 
                        label = ifelse(Symbol == "C7orf68", "HILPDA", Symbol)),
                    size = 5, color = "black",
                    box.padding = 0.3,
                    point.padding = 0.5,
                    max.overlaps = Inf) +
    scale_x_continuous(limits = c(-6, 6), breaks = seq(-4, 4, 4),
                       expand = expansion(0.01)) +  
    theme(panel.background = element_rect(fill = "white", color = "black"),
          panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
          panel.grid.major = element_line(color = "grey", linewidth = 0.1),
          panel.grid.minor = element_blank()))

# Create output directory for figures
if (!dir.exists(file.path(output_dir, "Figures"))) {
  dir.create(file.path(output_dir, "Figures"), recursive = TRUE)
}
# Save plot as svg for publication
ggsave(file.path(output_dir, "Figures", "Fig-2A-vulcan-plot-vector.svg"), 
       device = "svg", plot = FIG2A.vulcan, 
       bg = "white", width = 8, height = 6, units = "in")
# Save plot as png for presentation
ggsave(file.path(output_dir, "Figures", "Fig-2A-vulcan-plot-print.png"), 
       device = "png", plot = FIG2A.vulcan, 
       bg = "white", dpi = 300, width = 8, height = 6, units = "in")


# Clean up environment
rm(list = ls()[!(ls() %in% c("wd"))])
gc()
# End of script