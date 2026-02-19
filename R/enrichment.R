#' Extract gene list for enrichment analysis
#'
#' @description Extracts background and interest gene lists from differential expression results
#'
#' @param .df Data frame containing gene expression data
#' @param .filter Logical filter for genes of interest
#' @param .value Column name containing values to sort by (e.g., "log2FoldChange", "t")
#' @param .name Column name containing gene identifiers (e.g., "entrezID")
#'
#' @return List with two named vectors: background (all genes) and interest (filtered genes)
#' @export
#'
#' @importFrom dplyr filter distinct pull
#' @examples
#' \dontrun{
#' gene_list <- get_genelist(deg_results,
#'                           deg_results$padj < 0.05,
#'                           "log2FoldChange",
#'                           "entrezID")
#' }
get_genelist <- function(.df, .filter, .value, .name){
  # Extract the background gene list of every expressed gene
  background <- .df %>%
    dplyr::distinct(entrezID, .keep_all = TRUE) %>%
    dplyr::pull(.value, name = .name) %>%
    sort(., decreasing = TRUE)
  
  # Extract the gene list of interest of DEGs
  interest <- .df %>%
    dplyr::filter(.filter) %>%
    dplyr::distinct(entrezID, .keep_all = TRUE) %>%
    dplyr::pull(.value, name = .name) %>%
    sort(., decreasing = TRUE)
  
  return(list(background = background, interest = interest))
}


#' Over-representation analysis
#'
#' @description Performs over-representation analysis using hypergeometric testing
#'
#' @param .interest Named vector of genes of interest (with Entrez IDs as names)
#' @param .background Named vector of background genes (with Entrez IDs as names)
#' @param .pathways Data frame with pathway information (gs_name and human_entrez_gene columns)
#'
#' @return List containing enrichResult object
#' @export
#'
#' @importFrom clusterProfiler enricher setReadable
#' @importFrom dplyr select
#' @examples
#' \dontrun{
#' ora_result <- run_ora(interest_genes, background_genes, pathways)
#' }
run_ora <- function(.interest, .background, .pathways){
  ora <- clusterProfiler::enricher(
    gene = names(.interest), # Gene set of interest
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    pAdjustMethod = "BH",
    universe = names(.background), # Background gene set
    TERM2GENE = dplyr::select(
      .pathways,
      gs_name,
      human_entrez_gene
    )
  )
  ora <- clusterProfiler::setReadable(ora, org.Hs.eg.db::org.Hs.eg.db, keyType = "ENTREZID")
  return(list("ora" = ora))
}


#' Extract and format ORA results
#'
#' @description Extracts results from enricher object and formats with pathway IDs and descriptions
#'
#' @param .ora enrichResult object from clusterProfiler
#' @param .db Database with pathway annotations (gs_name, gs_exact_source, gs_description)
#'
#' @return List with two data frames: df (all results) and sig_df (significant results, p.adjust < 0.1)
#' @export
#'
#' @importFrom dplyr select distinct mutate filter
#' @importFrom stringr str_split
#' @examples
#' \dontrun{
#' ora_formatted <- extract_ora_results(ora_result$ora, pathways_db)
#' }
extract_ora_results <- function(.ora, .db){
  .db <- .db %>%
    dplyr::select(gs_name, gs_exact_source, gs_description) %>%
    dplyr::distinct()
  
  # Extract data frames
  df <- as.data.frame(.ora@result)
  # Order on p-value
  df <- df[order(df$`p.adjust`, decreasing = FALSE), ]
  # Change gene ratio to numeric values
  df <- df %>%
    dplyr::mutate(
      GeneRatio = sapply(stringr::str_split(df$GeneRatio, "/"),
                        function(x) 100 * (as.numeric(x[1]) / as.numeric(x[2]))),
      BgRatio = sapply(stringr::str_split(df$BgRatio, "/"),
                      function(x) 100 * (as.numeric(x[1]) / as.numeric(x[2])))
    )
  
  # Extract pathway IDs and description from database
  ids <- df$ID
  database <- sapply(stringr::str_split(ids, "_"),
                    function(x) return(x[1]))
  
  df$ID <- .db[match(ids, .db$gs_name), ][["gs_exact_source"]]
  df$Description <- .db[match(ids, .db$gs_name), ][["gs_description"]]
  df$Database <- database
  
  # Extract significant results: adjusted p-value < 0.1
  sig_df <- df %>%
    dplyr::filter(p.adjust < 0.1)
  
  # Return data frames
  return(list("df" = df, "sig_df" = sig_df))
}


#' Gene Set Enrichment Analysis
#'
#' @description Performs GSEA using a ranked gene list
#'
#' @param .geneset Named numeric vector of genes ranked by effect size (Entrez IDs as names)
#' @param .terms Data frame with gene set information (gs_name and human_entrez_gene columns)
#' @param minGSSize Minimum gene set size (default: 10)
#' @param maxGSSize Maximum gene set size (default: 500)
#'
#' @return List containing gseaResult object
#' @export
#'
#' @importFrom clusterProfiler GSEA setReadable
#' @importFrom dplyr select
#' @examples
#' \dontrun{
#' gsea_result <- run_gsea(ranked_genes, msigdb_terms)
#' }
run_gsea <- function(.geneset, .terms, minGSSize = 10, maxGSSize = 500){
  set.seed(42)
  ## Run the GSEA analysis
  res <- clusterProfiler::GSEA(
    geneList = .geneset, # Gene set of interest (ordered on effect size)
    minGSSize = minGSSize, # Minimum size of the gene set
    maxGSSize = maxGSSize, # Maximum size of the gene set
    pvalueCutoff = 1, # Adjusted p-value cutoff
    eps = 0, # P-value cutoff (minimum)
    seed = TRUE, # Seed for reproducibility
    pAdjustMethod = "BH", # P-value adjustment method
    TERM2GENE = dplyr::select(
      .terms,
      gs_name,
      human_entrez_gene
    )
  )
  res <- clusterProfiler::setReadable(res, org.Hs.eg.db::org.Hs.eg.db, keyType = "ENTREZID")
  
  # Extract data frame
  return(list("gsea" = res))
}


#' Extract and format GSEA results
#'
#' @description Extracts results from GSEA object and formats with pathway IDs and descriptions
#'
#' @param .gsea gseaResult object from clusterProfiler
#' @param .db Database with pathway annotations (gs_name, gs_exact_source, gs_description, gs_subcat)
#'
#' @return List with two data frames: df (all results) and sig_df (significant results, p.adjust < 0.1)
#' @export
#'
#' @importFrom dplyr select distinct mutate filter case_when
#' @importFrom stringr str_split
#' @examples
#' \dontrun{
#' gsea_formatted <- extract_gsea_results(gsea_result$gsea, terms_db)
#' }
extract_gsea_results <- function(.gsea, .db){
  .db <- .db %>%
    dplyr::select(gs_name, gs_exact_source, gs_description, gs_subcat) %>%
    dplyr::distinct()
  
  # Extract data frames
  df <- as.data.frame(.gsea@result)
  # Order on p-value
  df <- df[order(df$`p.adjust`, decreasing = FALSE), ]
  
  df <- df %>%
    dplyr::mutate(Direction = dplyr::case_when(
      NES > 0 ~ '+',
      NES < 0 ~ '-'
    ))
  
  df <- df %>%
    dplyr::mutate(
      Database = sapply(stringr::str_split(ID, "_"),
                       function(x) return(x[1])),
      Name = sapply(stringr::str_split(ID, "_"),
                   function(x) return(paste(x[-1], collapse = "_")))
    )
  
  # Add the pathway IDs and descriptions to the data frame
  ids <- df$ID
  df$ID <- .db[match(ids, .db$gs_name), ][["gs_exact_source"]]
  df$Description <- .db[match(ids, .db$gs_name), ][["gs_description"]]
  df$Database <- .db[match(ids, .db$gs_name), ][["gs_subcat"]]
  
  # Extract significant results: adjusted p-value < 0.1
  sig_df <- df %>%
    dplyr::filter(p.adjust < 0.1)
  
  # Return data frames
  return(list("df" = df, "sig_df" = sig_df))
}

#' Create enrichment summary table
#'
#' @description Formats GSEA results into a clean table for display
#'
#' @param .df Data frame with GSEA results
#' @param .order Column name for ordering (can be empty string)
#' @param .name Column name to use as row names
#'
#' @return Data frame with formatted enrichment statistics
#' @export
#'
#' @importFrom dplyr arrange mutate select case_when
#' @importFrom tibble remove_rownames column_to_rownames
#' @examples
#' \dontrun{
#' table <- get_enrichment_table(gsea_results$df, "", "Name")
#' }
get_enrichment_table <- function(.df, .order, .name){
  if (.order == "") {
    result <- .df
  } else {
    result <- .df %>% dplyr::arrange(!!rlang::sym(.order))
  }
  
  result <- result %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames(.name) %>%
    dplyr::mutate(
      NES = round(NES, 2),
      FDR = round(p.adjust, 4)
    ) %>%
    dplyr::mutate(
      Sign = dplyr::case_when(
        FDR < 0.001 ~ "****",
        FDR < 0.01 ~ "***",
        FDR < 0.05 ~ "**",
        FDR < 0.1 ~ "*",
        TRUE ~ as.character(FDR)
      )
    ) %>%
    dplyr::select(c(NES, FDR, Sign))
  
  return(result)
}

#' Calculate GSEA scores for a gene set
#'
#' @description Helper function to calculate running enrichment scores
#'
#' @param geneList Ranked gene list
#' @param geneSet Gene set to test
#' @param exponent Exponent for weighting (default: 1)
#' @param fortify Return as data frame (default: FALSE)
#'
#' @return List or data frame with enrichment scores
#' @keywords internal
gsea_scores <- function(geneList, geneSet, exponent = 1, fortify = FALSE){
  geneSet <- intersect(geneSet, names(geneList))
  N <- length(geneList)
  Nh <- length(geneSet)
  Phit <- Pmiss <- numeric(N)
  hits <- names(geneList) %in% geneSet
  Phit[hits] <- abs(geneList[hits])^exponent
  NR <- sum(Phit)
  Phit <- cumsum(Phit / NR)
  Pmiss[!hits] <- 1 / (N - Nh)
  Pmiss <- cumsum(Pmiss)
  runningES <- Phit - Pmiss
  max.ES <- max(runningES)
  min.ES <- min(runningES)
  if (abs(max.ES) > abs(min.ES)) {
    ES <- max.ES
  } else {
    ES <- min.ES
  }
  df <- data.frame(
    x = seq_along(runningES),
    runningScore = runningES,
    position = as.integer(hits)
  )
  if (fortify == TRUE) {
    return(df)
  }
  df$gene <- names(geneList)
  res <- list(ES = ES, runningES = df)
  return(res)
}


#' Extract detailed enrichment statistics for a gene set
#'
#' @description Extracts enrichment scores and gene positions for a single term from GSEA results
#'
#' @param object gseaResult object from clusterProfiler
#' @param geneSetID Gene set ID or numeric index
#'
#' @return Data frame with enrichment details
#' @export
#'
#' @examples
#' \dontrun{
#' # Extract enrichment info for the first pathway
#' enrich_info <- gs_info(gsea_result$gsea, 1)
#' }
gs_info <- function(object, geneSetID){
  geneList <- object@geneList
  if (is.numeric(geneSetID))
    geneSetID <- object@result[geneSetID, "ID"]
  geneSet <- object@geneSets[[geneSetID]]
  exponent <- object@params[["exponent"]]
  df <- gsea_scores(geneList, geneSet, exponent, fortify = TRUE)
  df$ymin <- 0
  df$ymax <- 0
  pos <- df$position == 1
  h <- diff(range(df$runningScore)) / 20
  df$ymin[pos] <- -h
  df$ymax[pos] <- h
  df$geneList <- geneList
  if (length(object@gene2Symbol) == 0) {
    df$gene <- names(geneList)
  } else {
    df$gene <- object@gene2Symbol[names(geneList)]
  }
  df$Description <- object@result[geneSetID, "Description"]
  return(df)
}
