################################################################################
# Created: 2024 10 25 ; Last Modified: 2025 11 09 ; MH                         #
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

# load MSigDB gene sets
load(file = "./RData/MSigDB_gene_sets.RData")

# Calculate Cohen's similarity between all pathways, terms and hallmark
if (!file.exists("./RData/term_similarity_matrix.RData")){
  # extract a list of gene sets from every reference used:
  total.genesets <- c(split(terms$gene_symbol, # GO terms and MSigDb hallmark sets
                            terms$gs_exact_source),
                      split(pathways$gene_symbol, # KEGG and Reactome pathways
                            pathways$gs_exact_source))
  # calculate the similarity matrix
  similarity.matrix <- cohen_kappa(total.genesets)
  # save the similarity matrix as RData
  save(similarity.matrix, file = "./RData/term_similarity_matrix.RData")
  # save file
  dir.create("../data/processed/similarity_matrix", showWarnings = FALSE)
  write.xlsx(similarity.matrix,
             file = "../data/processed/similarity_matrix/term_similarity_matrix.xlsx")
} else {
  load("./RData/term_similarity_matrix.RData")
}

# Clear the environment
rm(list = ls()) 
# Run garbage collection to free up memory
gc() 