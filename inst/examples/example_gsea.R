################################################################################
# Example: GSEA Analysis Using beltingSeq Package                               #
################################################################################
# Load the package
library(beltingSeq)

# Load additional packages
library(msigdbr)
library(dplyr)
library(openxlsx)

# Set working directory
wd <- getwd()

################################################################################
# Download and prepare MSigDB gene sets                                        #
################################################################################

msigdbr_df <- msigdbr(species = "Homo sapiens")

# Filter KEGG and Reactome pathways
pathways <- msigdbr_df %>%
  dplyr::filter(
    gs_cat == "C2",
    gs_subcat %in% c("CP:KEGG", "CP:REACTOME")
  )

# Filter Hallmark and GO:BP gene sets
terms <- msigdbr_df %>%
  dplyr::filter(
    gs_cat == "H" | (gs_cat == "C5" & gs_subcat == "GO:BP")
  ) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    gs_subcat = ifelse(gs_cat == "H", "HALLMARK", gs_subcat),
    gs_exact_source = ifelse(gs_cat == "H", gs_id, gs_exact_source)
  ) %>%
  dplyr::ungroup()

# Save gene sets
save(msigdbr_df, pathways, terms, 
     file = "./RData/MSigDB_gene_sets.RData")

################################################################################
# EXAMPLE 1: GSEA on single comparison                                         #
################################################################################

# Load processed data
load(file.path(wd, "RData", "CCLD_processedData.RData"))

# Extract gene list (beltingSeq function)
CCLD.genes <- get_genelist(
  .df = CCLD.df,
  .filter = CCLD.df[["significance"]] %in% c("Signif. up-regulated", 
                                             "Signif. down-regulated"),
  .value = "log2FoldChange",
  .name = "entrezID"
)

# GSEA on GO terms (beltingSeq function)
CCLD.GO <- run_gsea(
  .geneset = CCLD.genes$background, 
  .terms = terms
)
CCLD.GO <- c(CCLD.GO, 
             extract_gsea_results(.gsea = CCLD.GO$gsea, .db = terms))

# GSEA on pathways (beltingSeq function)
CCLD.GSEA <- run_gsea(
  .geneset = CCLD.genes$background,
  .terms = pathways
)
CCLD.GSEA <- c(CCLD.GSEA, 
               extract_gsea_results(.gsea = CCLD.GSEA$gsea, .db = pathways))

# Save results
save(CCLD.GO, CCLD.GSEA,
     file = file.path(wd, "RData", "CCLD_enrichment_results.RData"))

# Export significant results
write.xlsx(list(GO = CCLD.GO$sig_df,
                Pathways = CCLD.GSEA$sig_df),
           file = "etc/input-files/CCLD_significant_enrichments.xlsx")

################################################################################
# EXAMPLE 2: GSEA on multiple comparisons                                      #
################################################################################

# Load processed data
load(file.path(wd, "RData", "HGCC_processedData.RData"))

# Extract gene lists for each comparison (beltingSeq function)
HGCC.genes <- lapply(HGCC.deg, function(x){
  get_genelist(
    .df = x, 
    .filter = x[["significance"]] %in% c("Signif. up-regulated", 
                                        "Signif. down-regulated"),
    .value = "t",
    .name = "entrezID"
  )
})

# GSEA on GO terms for all comparisons (beltingSeq function)
HGCC.GO <- lapply(HGCC.genes, function(x){
  res <- run_gsea(.geneset = x[["background"]], .terms = terms)
  c(res, extract_gsea_results(.gsea = res$gsea, .db = terms))
})

# GSEA on pathways for all comparisons (beltingSeq function)
HGCC.GSEA <- lapply(HGCC.genes, function(x){
  res <- run_gsea(.geneset = x[["background"]], .terms = pathways)
  c(res, extract_gsea_results(.gsea = res$gsea, .db = pathways))
})

# Save results
save(HGCC.GO, HGCC.GSEA,
     file = file.path(wd, "RData", "HGCC_enrichment_results.RData"))

# Export significant results for each comparison
for(name in names(HGCC.GO)){
  write.xlsx(list(GO = HGCC.GO[[name]]$sig_df,
                  Pathways = HGCC.GSEA[[name]]$sig_df),
             file = paste0("etc/input-files/", name, 
                          "_significant_enrichments.xlsx"))
}

################################################################################
# Summary                                                                       #
################################################################################

cat("\n========================================\n")
cat("GSEA analysis complete!\n")
cat("Using beltingSeq package functions:\n")
cat("  - get_genelist()\n")
cat("  - run_gsea()\n")
cat("  - extract_gsea_results()\n")
cat("\nResults saved in RData/\n")
cat("Significant enrichments exported to etc/input-files/\n")
cat("========================================\n")
