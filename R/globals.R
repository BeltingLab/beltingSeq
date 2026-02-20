#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom rlang sym
#' @importFrom stats complete.cases na.omit p.adjust setNames
NULL

# Declare global variables to avoid R CMD check NOTEs
# These are used in non-standard evaluation (NSE) contexts with dplyr/tidyverse
utils::globalVariables(c(
  # Column names used in dplyr operations
  ".",
  "Category",
  "ENTREZID",
  "FDR",
  "ID",
  "NES",
  "Name",
  "SYMBOL",
  "Sign",
  "TargetID",
  "cluster",
  "cluster_name",
  "core_enrichment",
  "entrezID",
  "facet_var",
  "from.x",
  "from.y",
  "geneID",
  "geneRatio",
  "gs_description",
  "gs_exact_source",
  "gs_name",
  "gs_subcat",
  "human_entrez_gene",
  "illuminaHumanv4ENTREZID",
  "interest",
  "log2FoldChange",
  "n",
  "node1",
  "node2",
  "p.adjust",
  "padj",
  "position",
  "regulation",
  "significance",
  "to.x",
  "to.y",
  "weight",
  "x",
  "y",
  "ymax",
  "ymin"
))
