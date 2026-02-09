################################################################################
# Created: 2024 05 24 ; Last Modified: 2025 11 07 ; KGO & MH                   #
################################################################################
# ---------------------- Set up the environment -----------------------------  #
# Set working directory
wd <- getwd()
# Load packages
if (file.exists(file.path(wd, "packages.R"))) {
  source(file = file.path(wd, "packages.R"))
} else {
  stop("Required file 'packages.R' not found in the working directory.")
}
# Source data processing functions
if (file.exists(file.path(wd, "functions.R"))) {
  source(file = file.path(wd, "functions.R"))
} else {
  stop("Required file 'functions.R' not found in the working directory.")
}
# Get current date for results directory
date <- format(Sys.Date(), "%Y-%m-%d")
# Create results directory
signature_dir <- "Results/Signature"
# Create tables and plots directories
dir.create(file.path(wd, signature_dir, "tables"), recursive = T, showWarnings = FALSE)
dir.create(file.path(wd, signature_dir, "plots"), recursive = T, showWarnings = FALSE)

################################################################################
#                       IvyGap downloaded from GlioVis                         #
################################################################################
# 1.) Load Ivy Gap data
IvyGap <- read.table("etc/IvyGap/2025-01-23_Ivy_GAP_expression.csv", header=TRUE,
                     sep = ",", dec = ".")
IvyGap.genes <- read.table("etc/IvyGap/2025-01-23_Ivy_GAP_genes.csv", header=TRUE,
                           sep = ",", dec = ".")
IvyGap <- merge(IvyGap.genes, IvyGap, by.x="gene_id", by.y="gene_id.rna_well_id")
row.names(IvyGap) <- IvyGap$gene_symbol
IvyGap <- IvyGap[,-c(1:5)]
colnames(IvyGap) <- gsub(x = colnames(IvyGap), pattern = "X", replacement = "")
# Load clinical data
clin1 <- read.table("etc/IvyGap/2024-05-24_Ivy_GAP_pheno.txt", header=TRUE)
clin2 <- read.table("etc/IvyGap/2025-01-23_Ivy_GAP_pheno.csv", header=TRUE,
                    sep = ",", dec = ".")
clin_IvyGap <- merge(clin2, clin1, by.x="rna_well_id", by.y="Sample")
rm(clin1,clin2)
# 2.) Preprocess Ivy Gap data
Prim_IG <- clin_IvyGap[which(clin_IvyGap$Recurrence=="Primary" & # filter for primary GBM
                               !is.na(clin_IvyGap$Histology)),] # no IDH information
sel_IvyGap <- IvyGap[,which(colnames(IvyGap) %in% Prim_IG$rna_well_id)]


#Save RDatas
save(IvyGap, clin_IvyGap, Prim_IG, sel_IvyGap,
     file = "./RData/IvyGap_processedData.RData")

################################################################################
# Clean up the environment                                                     #
################################################################################
# Run garbage collection to free up memory
gc() 
# Clear the environment
rm(list = ls()) 

