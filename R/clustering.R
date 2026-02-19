#' Calculate Cohen's Kappa coefficient between gene sets
#'
#' @description Computes pairwise Cohen's Kappa coefficient to measure similarity
#' between gene sets based on their overlapping genes
#'
#' @param .list Named list of character vectors, each containing gene identifiers
#'
#' @return Symmetric matrix of Kappa coefficients
#' @export
#'
#' @examples
#' \dontrun{
#' gene_sets <- list(
#'   set1 = c("gene1", "gene2", "gene3"),
#'   set2 = c("gene2", "gene3", "gene4"),
#'   set3 = c("gene5", "gene6", "gene7")
#' )
#' similarity <- cohen_kappa(gene_sets)
#' }
cohen_kappa <- function(.list){
  N <- length(.list)
  kappa_mat <- matrix(0, nrow = N, ncol = N,
                     dimnames = list(names(.list), names(.list)))
  diag(kappa_mat) <- 1
  
  total <- length(unique(unlist(.list)))
  
  for (i in 1:(N - 1)) {
    for (j in (i + 1):N) {
      genes_i <- .list[[i]]
      genes_j <- .list[[j]]
      both <- length(intersect(genes_i, genes_j))
      term_i <- length(base::setdiff(genes_i, genes_j))
      term_j <- length(base::setdiff(genes_j, genes_i))
      no_terms <- total - sum(both, term_i, term_j)
      observed <- (both + no_terms) / total
      chance <- (both + term_i) * (both + term_j)
      chance <- chance + (term_j + no_terms) * (term_i + no_terms)
      chance <- chance / total^2
      kappa_mat[j, i] <- kappa_mat[i, j] <- (observed - chance) / (1 - chance)
    }
  }
  return(kappa_mat)
}


#' Cluster gene sets using network analysis
#'
#' @description Creates a network graph of gene sets based on similarity matrix,
#' then clusters the network using Louvain algorithm
#'
#' @param .df Data frame of enrichment results with ID and core_enrichment columns
#' @param .matrix Similarity matrix (e.g., from cohen_kappa)
#' @param .column Column to use for ordering (currently not used)
#' @param .threshold Similarity threshold for creating edges (default: 0.25)
#'
#' @return List containing igraph object and data frame with cluster assignments
#' @export
#'
#' @importFrom igraph graph_from_adjacency_matrix cluster_louvain V membership
#' @importFrom dplyr inner_join
#' @examples
#' \dontrun{
#' clustered <- get_cluster(gsea_results, similarity_matrix, threshold = 0.25)
#' }
get_cluster <- function(.df, .matrix, .column = NULL, .threshold = 0.25){
  # Extract the list of enriched gene sets
  genes <- split(strsplit(.df[["core_enrichment"]], "/"), .df[["ID"]])
  
  # Subset the Cohen's matrix for the gene sets of interest
  similarity.matrix <- .matrix[.df[['ID']], .df[['ID']]]
  
  # Create an adjacency matrix based on the similarity matrix
  adjacency.matrix <- similarity.matrix > .threshold
  
  set.seed(42)
  # Create a network graph from the adjacency matrix
  graph <- igraph::graph_from_adjacency_matrix(adjacency.matrix, mode = "undirected", diag = FALSE)
  
  # Annotate the nodes of the graph with the enrichment results
  igraph::V(graph)$Name <- .df$Name
  igraph::V(graph)$NES <- .df$NES
  igraph::V(graph)$geneRatio <- .df$geneRatio
  
  # Community detection (Louvain algorithm)
  clusters <- igraph::cluster_louvain(graph)
  # Annotate nodes with cluster membership
  igraph::V(graph)$cluster <- igraph::membership(clusters)
  
  # Extract cluster information
  cluster_summary <- data.frame(
    ID = names(igraph::V(graph)),
    cluster = as.factor(igraph::V(graph)$cluster)
  )
  
  cluster_summary <- dplyr::distinct(cluster_summary, .keep_all = TRUE)
  
  df <- .df %>%
    dplyr::inner_join(., cluster_summary, by = "ID")
  
  return(list(graph = graph, df = df))
}


#' Get representative terms for each cluster
#'
#' @description Calculates hub scores for genes and selects representative terms
#' for each cluster based on gene expression and network topology
#'
#' @param .cluster Data frame with cluster assignments from get_cluster
#' @param .degs Data frame with differential expression results (Symbol, log2FoldChange, padj)
#'
#' @return Data frame with hub scores and representative term designation
#' @export
#'
#' @importFrom igraph graph_from_data_frame simplify hub_score
#' @importFrom dplyr select mutate group_by reframe inner_join arrange slice_head pull rowwise
#' @importFrom tidyr separate_rows
#' @examples
#' \dontrun{
#' representatives <- get_cluster_representative(clustered$df, deg_results)
#' }
get_cluster_representative <- function(.cluster, .degs){
  # Add gene expression values to the clusters
  linkage <- .cluster %>%
    dplyr::select(ID, Name, core_enrichment, cluster) %>%
    tidyr::separate_rows(core_enrichment, sep = "/")
  
  linkage$logFC <- .degs$log2FoldChange[match(linkage$core_enrichment, .degs$Symbol)]
  linkage$padj <- .degs$padj[match(linkage$core_enrichment, .degs$Symbol)]
  
  # Calculate normalized term weight
  linkage$weight <- abs(linkage$logFC) * (-log10(linkage$padj))
  linkage$weight <- linkage$weight / max(linkage$weight, na.rm = TRUE)
  
  linkage <- linkage %>%
    dplyr::select(c("core_enrichment", "ID", "weight", "cluster")) %>%
    setNames(., c("node1", "node2", "weight", "cluster")) %>%
    as.data.frame(.)
  
  # Create a network visualization of gene and GO-term relationships
  net <- igraph::graph_from_data_frame(linkage)
  net <- igraph::simplify(net, remove.multiple = FALSE, remove.loops = TRUE)
  
  # Calculate hub score of each gene
  hs <- igraph::hub_score(net, scale = TRUE, weights = linkage$weight)$vector
  
  # Summarize gene hub scores to determine the final weight of each GO term
  linkage$geneRatio <- .cluster$geneRatio[match(linkage$node2, .cluster$ID)]
  
  linkage <- linkage %>%
    dplyr::mutate(hub_score = geneRatio * hs[match(node1, names(hs))]) %>%
    dplyr::rename(geneID = node1, ID = node2) %>%
    dplyr::group_by(ID, cluster) %>%
    dplyr::reframe(
      core_enrichment = paste(geneID, collapse = "/"),
      hub_score = sum(hub_score, na.rm = TRUE)
    )
  
  out_df <- dplyr::inner_join(.cluster, linkage, by = c("ID", "core_enrichment", "cluster")) %>%
    dplyr::mutate(cluster = as.factor(cluster))
  
  # Select cluster representatives according to calculated weight
  representative.terms <- out_df %>%
    dplyr::group_by(cluster) %>%
    dplyr::arrange(dplyr::desc(hub_score)) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::pull(ID)
  
  out_df <- out_df %>%
    dplyr::rowwise(.) %>%
    dplyr::mutate(Representative = ifelse(ID %in% representative.terms, TRUE, FALSE))
  
  return(out_df)
}


#' Filter graph by cluster size
#'
#' @description Filters an igraph object to keep only clusters meeting a minimum size threshold
#'
#' @param .graph igraph object with cluster assignments in vertex attribute 'cluster'
#' @param .threshold Minimum number of nodes per cluster
#'
#' @return Filtered igraph subgraph
#' @export
#'
#' @importFrom igraph V induced_subgraph
#' @examples
#' \dontrun{
#' filtered_graph <- filter_graph(cluster_graph, threshold = 5)
#' }
filter_graph <- function(.graph, .threshold){
  cl <- table(igraph::V(.graph)$cluster)
  
  valid_cl <- which(cl >= .threshold)
  
  filtered_vertices <- which(igraph::V(.graph)$cluster %in% valid_cl)
  
  subgraph <- igraph::induced_subgraph(.graph, vids = filtered_vertices)
  
  return(subgraph)
}
