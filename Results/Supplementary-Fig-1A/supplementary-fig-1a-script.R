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
  load(file = "./RData/IvyGap_processedData.RData")
} else {
  stop("Required .RData file(s) not found in the working directory.")
}
# Create Supplementary Fig-1A: Enrichment score for gene signature in Ivy Gap data
# Define paths
output_dir <- file.path(wd, "Results", "Supplementary-Fig-1A")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
# load gene set of interest
interest_genes <- read.csv("etc/gene-signature.txt", header = T,
                           sep = "\t", stringsAsFactors = T)
# Calculate lipid droplet gene signature
Ivygap.score <- hack_sig(as.matrix(sel_IvyGap),
                         list(as.character(interest_genes$SYMBOL)),
                         method = "zscore")
Ivygap.score <- setNames(as.data.frame(Ivygap.score), c("sample","zscore"))
Ivygap.score <- merge(Ivygap.score, Prim_IG, by.x="sample", by.y="rna_well_id")
Ivygap.score <- Ivygap.score %>%
  dplyr::select(sample, zscore, Age, survival, status, Histology,
                Subtype.svm:Subtype_Verhaak_2010, MGMT_status, EGFR_amplification) %>%
  dplyr::mutate(
    Histology = factor(Histology,levels = c("Pseudopalisading cells", "Microvascular proliferation",
                                            "Leading Edge", "Infiltrating Tumor", "Cellular Tumor")),
    MGMT_status = as.factor(MGMT_status),
    EGFR_amplification = as.factor(EGFR_amplification)
  )

# Plot histology
(SUPP1A.plot <- score_plot(.score = Ivygap.score, .formula = "zscore ~ Histology",
                                     .ref_group = "Pseudopalisading cells", .x = "Histology",
                                     .y = "zscore", .title = "LD+/CS+ gene signature"))


# Save results table
SUPP1A.df <- SUPP1A.plot$stat.test %>% 
  dplyr::select(!"p.format") %>% 
  setNames(.,
           c("Metric", "Reference Group", "Comparison Group",
             "P-value", "Adjusted P-value", "Significance Level", "Stat. test"))
openxlsx::write.xlsx(SUPP1A.df,
                     file.path(output_dir, "Supplementary-Fig-1A-table.xlsx"))

# Create output directory for figures
if (!dir.exists(file.path(output_dir, "Figures"))) {
  dir.create(file.path(output_dir, "Figures"), recursive = TRUE)
}
# Save as svg for publication
ggsave(file.path(output_dir, "Figures", "Supplementary-Fig-1A-signature-vector.svg"),
       device = "svg", plot = SUPP1A.plot$plot,
       bg = "white", width = 6, height = 4.45, units = "in")
# Save as png for presentation
ggsave(file.path(output_dir, "Figures", "Supplementary-Fig-1A-signature-print.png"),
       device = "png", plot = SUPP1A.plot$plot,
       bg = "white", dpi = 300, width = 6, height = 4.45, units = "in")

# Clean up environment
rm(list = ls()[!(ls() %in% c("wd"))])
gc()
# End of script