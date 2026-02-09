################################################################################
# Created: 2024 10 07 ; Last Modified: 2025 06 24 ; MH                         #
################################################################################
# ---------------------- Set up the environment -----------------------------  #
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

################################################################################
# (!!!) AS OF 20-Oct-2023                                                      #
################################################################################
# load databases
msigdbr_df <- msigdbr(species = "Homo sapiens") # save local database

# Filter KEGG and Reactome pathways
pathways <- msigdbr_df %>%
  dplyr::filter(
    gs_cat == "C2", # only canonical representations (compiled by experts)
    gs_subcat %in% c("CP:KEGG", "CP:REACTOME") # KEGG and Reactome pathways
  )
# Filter Hallmakr and GO:BP gene sets
terms <- msigdbr_df %>%
  dplyr::filter(
    gs_cat == "H" | # Hallmark gene sets
    gs_cat == "C5" & # Ontology gene sets
    gs_subcat == "GO:BP" # GO terms, biological processes
  )
# standardize gene set labels
terms <- terms %>%
  dplyr::rowwise(.) %>%
  dplyr::mutate(gs_subcat = ifelse(gs_cat == "H", "HALLMARK", gs_subcat),
                gs_exact_source = ifelse(gs_cat == "H", gs_id, gs_exact_source)) %>%
  dplyr::ungroup(.)
#Save gene sets
save(msigdbr_df, pathways, terms, file = "./RData/MSigDB_gene_sets.RData")

################################################################################
#             Differential gene expression analysis                            #
################################################################################

# ---------------------------------------------------------------------------- #
# -    1.) Svenja's CC +LD and CC no LD (Clariom D Human Pico) Affymetrix    - #
# ---------------------------------------------------------------------------- #
# Load the processed data
load(file.path(wd, "RData", "CCLD_processedData.RData"))

## extract entrez IDs for gene set of interest and background
CCLD.genes <- get_genelist(.df = CCLD.df,
                           .filter = CCLD.df[["significance"]] %in% c("Signif. up-regulated", 
                                                                "Signif. down-regulated"),
                           .value = "log2FoldChange",
                           .name = "entrezID")

# GSEA on the GO terms from MSigDB
CCLD.GO <- list()
CCLD.GO <- run_gsea(.geneset = CCLD.genes$background, 
                    .terms = terms)
CCLD.GO <- c(CCLD.GO, extract_gsea_results(.gsea = CCLD.GO$gsea, .db = terms))

# GSEA on the KEGG- and REACTOME pathways from MSigDB
CCLD.GSEA <- list()
CCLD.GSEA <- run_gsea(.geneset = CCLD.genes$background,
                      .terms = pathways)
CCLD.GSEA <- c(CCLD.GSEA, extract_gsea_results(.gsea = CCLD.GSEA$gsea, .db = pathways))

# Save enrichment results
save(CCLD.GO, CCLD.GSEA,
     file = file.path(wd, "RData", "CCLD_enrichment_results.RData"))

# ---------------------------------------------------------------------------- #
# -   2.) Hugo's Primary cells 2D vs 3D (Clariom D Human Pico) Affymetrix    - #
# ---------------------------------------------------------------------------- #
# Load the processed data
load(file.path(wd, "RData", "HGCC_processedData.RData"))

# extract entrez IDs for gene set of interest and background
HGCC.genes <- list()
HGCC.genes <- lapply(HGCC.deg, function(x){
  get_genelist(.df = x, 
               .filter = x[["significance"]] %in% c("Signif. up-regulated", 
                                                    "Signif. down-regulated"),
               .value = "t",
               .name = "entrezID")
})

# GSEA on the GO terms from MSigDB
HGCC.GO <- list()
HGCC.GO <- lapply(HGCC.genes, function(x){
  run_gsea(.geneset = x[["background"]], .terms = terms)})
HGCC.GO <- lapply(HGCC.GO,
                  function(x){
                    return(c(x, 
                             extract_gsea_results(
                               .gsea = x[["gsea"]],
                               .db = terms)))
                  })

# GSEA on the KEGG- and REACTOME pathways from MSigDB
HGCC.GSEA <- list()
HGCC.GSEA <- lapply(HGCC.genes, function(x){
  run_gsea(.geneset = x[["background"]], .terms = pathways)})
HGCC.GSEA <- lapply(HGCC.GSEA,
                    function(x){
                      return(c(x, 
                               extract_gsea_results(
                                 .gsea = x[["gsea"]],
                                 .db = pathways)))
                    })

# Save enrichment results
save(HGCC.GO, HGCC.GSEA,
     file = file.path(wd, "RData", "HGCC_enrichment_results.RData"))

# ---------------------------------------------------------------------------- #
# -   3.) U87 Chronic Acidosis AA vs NA & HOX vs NOX (Illumina BeadChip)     - #
# ---------------------------------------------------------------------------- #
# Load the processed data
load(file.path(wd, "RData", "U87_processedData.RData"))

# extract entrez IDs for gene set of interest and background
U87.genes <- list()
U87.genes <- lapply(U87.deg, function(x){
  get_genelist(.df = x, 
               .filter = x[["significance"]] %in% c("Signif. up-regulated", 
                                                    "Signif. down-regulated"),
               .value = "t",
               .name = "entrezID")
})

# Run the GSEA analysis
U87.GO <- list()
U87.GO <- lapply(U87.genes, function(x){
  run_gsea(.geneset = x[["background"]], .terms = terms)})
U87.GO <- lapply(U87.GO,
                  function(x){
                    return(c(x, 
                             extract_gsea_results(
                               .gsea = x[["gsea"]],
                               .db = terms)))
                  })

GSEA <- list()
U87.GSEA <- lapply(U87.genes, function(x){
  run_gsea(.geneset = x[["background"]], .terms = pathways)})
U87.GSEA <- lapply(U87.GSEA,
                    function(x){
                      return(c(x, 
                               extract_gsea_results(
                                 .gsea = x[["gsea"]],
                                 .db = pathways)))
                    })

# Save enrichment results
save(U87.GO, U87.GSEA,
     file = file.path(wd, "RData", "U87_enrichment_results.RData"))

# ---------------------------------------------------------------------------- #
# -   4.) PANC1 Chronic Acidosis AA vs NA (Clariom D Human Pico) Affymetrix  - #
# ---------------------------------------------------------------------------- #
load(file.path(wd, "RData", "PANC1_processedData.RData"))

## extract entrez IDs for gene set of interest and background
PANC1.genes <- get_genelist(.df = PANC1.deg,
                            .filter = PANC1.deg[["significance"]] %in% c("Signif. up-regulated", 
                                                                         "Signif. down-regulated"),
                            .value = "log2FoldChange",
                            .name = "entrezID")

# GSEA on the GO terms and hallmarks from MSigDB
PANC1.GO <- run_gsea(.geneset = PANC1.genes$background,
                     .terms = terms)
PANC1.GO <- c(PANC1.GO, 
              extract_gsea_results(.gsea = PANC1.GO$gsea, 
                                   .db = terms))

# GSEA on the KEGG- and REACTOME pathways from MSigDB
PANC1.GSEA <- run_gsea(.geneset = PANC1.genes$background,
                       .terms = pathways)
PANC1.GSEA <- c(PANC1.GSEA, 
                extract_gsea_results(.gsea = PANC1.GSEA$gsea, 
                                     .db = pathways))

# Save GSEA resutls
save(PANC1.GO, PANC1.GSEA,
     file = file.path(wd, "RData", "PANC1_enrichment_results.RData"))

################################################################################
# Clean up the environment                                                     #
################################################################################
# Run garbage collection to free up memory
gc() 
# Clear the environment
rm(list = ls()) 