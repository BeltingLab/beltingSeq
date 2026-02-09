##########################################################################
# Created: 2024 03 27 ; Last Modified: 2025 11 08 ; MH                   #
##########################################################################
# ---------------------------------------------------------------------- #
# Function - Microarray data processing                                  #
# ---------------------------------------------------------------------- #
# I.) Normailze transcript abundance
# The function normalizes the transcript abundance data using the RMA method
normalizeTranscript <- function(.data, .db){
  #rma does background normalization, log2 transformation and quantile normalization
  data <- rma(.data, target = "core")
  # annotating using info assembled by the Bioconductor Core Team
  data <- annotateEset(data, .db)
  
  # return expression matrix and add annotation
  df <- as.data.frame(exprs(data))
  df$PROBEID<- rownames(df)
  
  # Annotate gene names
  annot <- data@featureData@data
  df <- merge(df, annot, by="PROBEID")
  
  #Filter the gene expression to contain only values with genesymbol
  idx <- which(is.na(df$SYMBOL))
  df<- df[-idx,]
  
  return(df)
}

# II.) Remove duplicates
# Remove duplicate symbols based on lgFC
removeDuplicates <- function(.data, .column, .symbol){
  # Sort the data frame based on the log2FC values
  data <- .data[order(abs(.data[[.column]]), decreasing = TRUE),]
  # Remove duplicates
  data <- data[!duplicated(data[[.symbol]]),]
  return(data)
  
}

# III.) Differential expression analysis with Limma
limmaDEA <- function(.data, .design, .contrast){
  # Extract numeric columns
  mat <- .data[,sapply(.data, is.numeric)]
  colnames(mat) <- gsub("_\\(.*$", "\\1", colnames(mat))
  rownames(mat) <- .data$PROBEID
  # Create a design matrix
  t <- as.factor(.design)
  levels(t) <- make.names(levels(t)) # rename levels to ensure syntactically valid names
  
  design <- model.matrix(~0 + t)
  colnames(design) <- levels(t)
  
  # Fit the linear model
  fit <- lmFit(mat, design)
  fit$genes$ID <- rownames(mat)
  fit$genes$Symbol <- .data$SYMBOL
  fit$genes$entrezID <- .data$ENTREZID
  
  # Fit the contrasts
  contrast.matrix <- makeContrasts(contrasts = .contrast, levels = t)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  
  # Make deg table with getfitlimma
  res <- list()
  for (i in 1:ncol(contrast.matrix)){
    res[[i]] <- topTable(fit2, coef = i, number = Inf)
  }
  
  return(res)
}

# ---------------------------------------------------------------------- #
# Function - Illmina beadChip data preprocessing                         #
# ---------------------------------------------------------------------- #
# I.) Normalize transcript abundance using control probes
tryMapID <- function(ID, inp, outp){
  tryCatch(
    {
      suppressWarnings(
        mapIds(org.Hs.eg.db, keys = ID, column = outp, keytype = inp, multiVals = "first")
      )
    },
    error = function(e) {
      message(conditionMessage(e))
      return(NA)
    }
  )
}

switchAlias <- function(dbconn, ID){
  alias <- as.character(dbconn[which(dbconn$alias == ID),]$symbol)[1]
  
  if (is.na(alias)){
    return(ID)
  } else {
    return(alias)
  }  # returns TRUE
}


normalizeIllumina <- function(.df, .dbconn, .samples = colnames(.df$E)){
  x <- illuminaHumanv4ENTREZID
  
  #normalization and background correction
  df <- neqc(.df)
  # keep only selected samples
  df <- df[, .samples]
  #keep probes that are expressed in at least three arrays according to a detection p-values of 5%:
  expressed <- rowSums(df$other$Detection < 0.05) >= 3
  df <- df[expressed,]
  
  # get Illumina annotations
  annot <- df$genes
  annot <- cbind(PROBEID=rownames(annot), 
                 ENTREZID = sapply(as.list(x[rownames(annot)]), "[[", 1),
                 annot)
  annot <- annot %>%
    dplyr::select(!TargetID) %>%
    dplyr::filter(SYMBOL != "" & !is.na(SYMBOL)) %>% 
    dplyr::rowwise() %>% 
    # replace outdated gene symbols with the most recent ones
    dplyr::mutate(
      SYMBOL = ifelse(
        ENTREZID == "NA" | is.na(ENTREZID), no = SYMBOL, yes =  switchAlias(.dbconn, SYMBOL))
    )
  
  # map Illumina annotations against the org.Hs.eg.db database
  entrez <- tryMapID(annot$SYMBOL, "SYMBOL", "ENTREZID")
  # refine annotations missing from the Illumina database
  annot <- annot %>%
    rowwise() %>%
    dplyr::mutate(
      ENTREZID = ifelse(
        ENTREZID == "NA" | is.na(ENTREZID), no = ENTREZID, yes = entrez[[SYMBOL]]))
  
  # add annotation to the expression matrix
  df <- as.data.frame(df$E)
  df$PROBEID <- rownames(df)
  df <- merge(annot, df, by="PROBEID")
  # remove rows with missing ENTREZID
  idx <- which(is.na(df$ENTREZID))
  df <- df[-idx,]
  return(df)
}

# ---------------------------------------------------------------------- #
# Function - DIFF.EXP.                                                   #
# ---------------------------------------------------------------------- #

# I.) SET SIGNIFICANCE LEVELS
get_significance <- function(.df){
  return(.df %>% 
           # Rename data frame columns to make sense
           # setNames(., c("symbol", "log2FoldChange", "pvalue", "padj")) %>% 
           # Add a columns...
           dplyr::mutate(
             # ... for significance levels using the thresholds:
             #                            p-adj < 0.05, abs(log2FC) > 0.5
             significance = dplyr::case_when(
               abs(log2FoldChange) > 0.5 & padj > 0.05 ~ 'log2FoldChange',
               abs(log2FoldChange) < 0.5 & padj < 0.05 ~ 'Log10P',
               log2FoldChange < (-1)*0.5 & padj < 0.05 ~ 'Signif. down-regulated',
               log2FoldChange > 0.5 & padj < 0.05 ~ 'Signif. up-regulated',
               T ~ 'NS')) %>% 
           dplyr::filter(complete.cases(.)))
}

plot_vulcan <- function(.data, label = T){
  plot = ggplot(data = na.omit(.data), 
                aes(x = log2FoldChange, y = -log10(padj), colour = significance)) + 
    geom_point(mapping = aes(), inherit.aes = T, size = 2.5, alpha = 0.35) + 
    scale_color_manual(values = c("NS" = '#c1c1c1',
                                  "Log10P" = '#363636',
                                  "Log2FoldChange" = '#767676',
                                  "Signif. up-regulated" = '#841f27',
                                  "Signif. down-regulated" = '#000f64'),
                       name = "Significance") + 
    labs(x = expression(paste(log[2], 'FoldChange')),
         y = expression(paste(log[10], italic('FDR')))) +
    scale_x_continuous(expand = expansion(0.2)) + 
    # Visualize log2FC threshold
    geom_vline(xintercept = c(-0.5, 0.5), linetype = 'dotted', size = 1) +
    # Visualize adjusted p-value threshold
    geom_hline(yintercept = -log10(0.05), linetype = 'dotted', size = 1) +
    theme(axis.title = element_text(size = 14), 
          axis.text = element_text(size = 14), 
          legend.position = 'none')
  
  if (label == T){
    return(plot +
             # Add labels for significantly up- or down-regulated genes
             geom_text(data = subset(.data,
                                     significance %in% c("Signif. up-regulated", 
                                                         "Signif. down-regulated")),
                       hjust = 0, vjust = 1.5, colour = 'black', position = 'identity', 
                       show.legend = F, check_overlap = T,
                       label = subset(.data,
                                      significance %in% c("Signif. up-regulated", 
                                                          "Signif. down-regulated"))[,"Symbol"]))
  } else {
    return(plot)
  }
}

# ---------------------------------------------------------------------- #
# Function - PATHWAY ANALYSES                                            #
# ---------------------------------------------------------------------- #

# I. EXTRACT GENE LISTS
get_genelist <- function(.df, .filter, .value, .name){
  # Extract the background gene list of every expressed gene
  background <- .df %>%
    dplyr::distinct(entrezID, .keep_all = T) %>%
    dplyr::pull(.value, name = .name) %>% 
    sort(., decreasing = T)
  # Extract the gene list of interest of DEGs
  interest <- .df %>%
    dplyr::filter(.filter) %>%
    dplyr::distinct(entrezID, .keep_all = T) %>%
    dplyr::pull(.value, name = .name) %>% 
    sort(., decreasing = T)
  
  return(list(background = background, interest = interest))
}

# II. OVER-REPRESENTATION ANALYSIS
run_ora <- function(.interest, .background, .pathways){
  ora <- enricher(gene = names(.interest), # gene set of interest 
                  pvalueCutoff = 1, 
                  qvalueCutoff = 1,
                  pAdjustMethod = "BH", 
                  universe = names(.background), # background gene set
                  TERM2GENE = dplyr::select(
                    .pathways,
                    gs_name,
                    human_entrez_gene))
  ora <- setReadable(ora, org.Hs.eg.db, keyType = "ENTREZID")
  return(list("ora" = ora))
}

extract_ora_results <- function(.ora, .db){
  .db <- .db %>% 
    dplyr::select(gs_name, gs_exact_source, gs_description) %>% 
    dplyr::distinct()
  # extract data frames
  df <- as.data.frame(.ora@result)
  # order on p-value
  df <- df[order(df$`p.adjust`, decreasing = F),]
  # change gene ratio to numeric values
  df <- df %>% 
    mutate(GeneRatio = sapply(stringr::str_split(df$GeneRatio, "/"), 
                              function(x) 100*(as.numeric(x[1])/as.numeric(x[2]))),
           BgRatio = sapply(stringr::str_split(df$BgRatio, "/"), 
                            function(x) 100*(as.numeric(x[1])/as.numeric(x[2]))))
  # extract pathway IDs and description from database
  
  # add the pathway IDs and descriptions to the data frame
  ids <- df$ID
  database <- sapply(stringr::str_split(ids, "_"), 
                     function(x) return(x[1]))
  
  df$ID = .db[match(ids, .db$gs_name),][["gs_exact_source"]]
  df$Description = .db[match(ids, .db$gs_name),][["gs_description"]]
  df$Database = database
  
  # extract significant results: adjusted p-value < 0.05
  sig_df <- df %>% 
    dplyr::filter(p.adjust < 0.1)
  
  #return data frames
  return(list("df" = df, "sig_df" = sig_df))
}

# III. GENE SET ENRICHMENT ANALYSIS
run_gsea <- function(.geneset, .terms){
  set.seed(42)
  ## run the GSEA analysis
  res <- GSEA(
    geneList = .geneset, # gene set of interest (ordered on effect size)
    minGSSize = 10, # minimum size of the gene set
    maxGSSize = 500, # maximum size of the gene set
    pvalueCutoff = 1, # adjusted p-value cutoff
    eps =  0, # p-value cutoff (minimum)
    seed = TRUE, # seed for reproducibility
    pAdjustMethod = "BH", # p-value adjustment method
    TERM2GENE = dplyr::select(
      .terms,
      gs_name,
      human_entrez_gene
    ))
  res <- setReadable(res, org.Hs.eg.db, keyType = "ENTREZID")
  
  # extract data frame
  return(list("gsea" = res))
}

extract_gsea_results <- function(.gsea, .db){
  .db <- .db %>% 
    dplyr::select(gs_name, gs_exact_source, gs_description) %>% 
    dplyr::distinct()
  # extract data frames
  df <- as.data.frame(.gsea@result)
  # order on p-value
  df <- df[order(df$`p.adjust`, decreasing = F),]
  
  df <- df %>%
    dplyr::mutate(Direction = dplyr::case_when(NES > 0 ~ '+',
                                               NES < 0 ~ '-'))
  df <- df %>%
    dplyr::mutate(
      Database = sapply(stringr::str_split(ID, "_"), 
                        function(x) return(x[1])),
      Name = sapply(stringr::str_split(ID, "_"), 
                    function(x) return(paste(x[-1], collapse = "_"))))
  
  # add the pathway IDs and descriptions to the data frame
  ids <- df$ID
  df$ID = .db[match(ids, .db$gs_name),][["gs_exact_source"]]
  df$Description = .db[match(ids, .db$gs_name),][["gs_description"]]
  df$Database = .db[match(ids, .db$gs_name),][["gs_subcat"]]
  
  # extract significant results: adjusted p-value < 0.05
  sig_df <- df %>% 
    dplyr::filter(p.adjust < 0.1)
  
  #return data frames
  return(list("df" = df, "sig_df" = sig_df))
}

# VI. PLOT GSEA RESULTS
# extract detailed enrichment statistics for a single term
gsInfo <- function (object, geneSetID) {
  geneList <- object@geneList
  if (is.numeric(geneSetID)) 
    geneSetID <- object@result[geneSetID, "ID"]
  geneSet <- object@geneSets[[geneSetID]]
  exponent <- object@params[["exponent"]]
  df <- gseaScores(geneList, geneSet, exponent, fortify = TRUE)
  df$ymin <- 0
  df$ymax <- 0
  pos <- df$position == 1
  h <- diff(range(df$runningScore))/20
  df$ymin[pos] <- -h
  df$ymax[pos] <- h
  df$geneList <- geneList
  if (length(object@gene2Symbol) == 0) {
    df$gene <- names(geneList)
  }
  else {
    df$gene <- object@gene2Symbol[names(geneList)]
  }
  df$Description <- object@result[geneSetID, "Description"]
  return(df)
}

gseaScores <- function (geneList, geneSet, exponent = 1, fortify = FALSE) 
{
  geneSet <- intersect(geneSet, names(geneList))
  N <- length(geneList)
  Nh <- length(geneSet)
  Phit <- Pmiss <- numeric(N)
  hits <- names(geneList) %in% geneSet
  Phit[hits] <- abs(geneList[hits])^exponent
  NR <- sum(Phit)
  Phit <- cumsum(Phit/NR)
  Pmiss[!hits] <- 1/(N - Nh)
  Pmiss <- cumsum(Pmiss)
  runningES <- Phit - Pmiss
  max.ES <- max(runningES)
  min.ES <- min(runningES)
  if (abs(max.ES) > abs(min.ES)) {
    ES <- max.ES
  }
  else {
    ES <- min.ES
  }
  df <- data.frame(x = seq_along(runningES), runningScore = runningES, 
                   position = as.integer(hits))
  if (fortify == TRUE) {
    return(df)
  }
  df$gene = names(geneList)
  res <- list(ES = ES, runningES = df)
  return(res)
}

# Plot running score for selected gene sets
plotRunningScore <- function(.df, .x, .y, .color, .palette){
  subset <- subset(.df, position == 1)
  
  return(ggplot(.df, aes(x = !!sym(.x)))
         + xlab(element_blank())
         + ylab("Enrichment score (ES)")
         + theme_bw(base_size = 14)
         + geom_hline(yintercept = 0, linetype = "dashed", color = "grey")
         + scale_x_continuous(expand = c(0, 0))
         + scale_y_continuous(expand = c(0.01, 0.01),
                              breaks = c(-0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5))
         + geom_point(aes(x = !!sym(.x), y = !!sym(.y), color = !!sym(.color)),
                      show.legend = F, size = 1, data = subset)
         + scale_color_manual(values = .palette)
         + theme(panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.grid.major.y = element_blank(),
                 panel.grid.minor.y = element_blank(),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.text.y = element_text(colour = "black", size = 14))
  )
}

plotGeneRank <- function(.df, .x, .color, .facet, .palette){
  subset <- subset(.df, position == 1)
  facet <- as.formula(.facet)
  
  
  return(ggplot(.df, aes(x = !!sym(.x))) 
         + geom_linerange(aes(ymin = ymin, ymax = ymax, color = !!sym(.color)), 
                          show.legend = F, size = 1, data = subset) 
         + theme_bw(base_size = 14) 
         + xlab("Rank in ordered gene list") 
         + scale_x_continuous(expand = c(0, 0)) 
         + scale_y_continuous(expand = c(0, 0)) 
         + scale_color_manual(values = .palette) 
         + facet_grid(facet, scales = "free_y", switch = "y",
                      labeller = labeller(.default = Hmisc::capitalize)) 
         + theme(panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(), 
                 panel.grid.major.y = element_blank(),
                 panel.grid.minor.y = element_blank(), 
                 axis.ticks.y = element_blank(), 
                 axis.text.y = element_blank(),
                 strip.placement = "outside",
                 strip.background = element_blank(),
                 strip.text.y.left = element_text(angle = 0, hjust = 0,
                                                  colour = "grey25", size = 14), 
                 axis.text.x = element_text(colour = "black", size = 14))
  )
}

getEnrichmentTable <- function(.df, .order, .name){
  return(.df %>%
           dplyr::arrange(.order) %>% 
           tibble::remove_rownames() %>% 
           tibble::column_to_rownames(.name) %>%
           dplyr::mutate(
             NES = round(NES, 2),
             FDR = round(p.adjust, 4)) %>%
           dplyr::mutate(
             Sign = case_when(
               # ****P<0.001, ***P< 0.01, **P< 0.05, and *P< 0.1
               FDR < 0.001 ~ "****",
               FDR < 0.01 ~ "***",
               FDR < 0.05 ~ "**",
               FDR < 0.1 ~ "*",
               TRUE ~ as.character(FDR)
             )) %>% 
           dplyr::select(c(NES, FDR, Sign)))
}

# ---------------------------------------------------------------------- #
# Function - Cluster analysis                                            #
# ---------------------------------------------------------------------- #
# a function calculating the Cohen's Kappa coefficient between a list of gene sets
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
      observed <- (both + no_terms)/total
      chance <- (both + term_i) * (both + term_j)
      chance <- chance + (term_j + no_terms) * (term_i + 
                                                  no_terms)
      chance <- chance/total^2
      kappa_mat[j, i] <- kappa_mat[i, j] <- (observed - 
                                               chance)/(1 - chance)
    }}
  return(kappa_mat)
}

# A function to create a network graph of gene sets based on the Cohen's Kappa coefficient
# then cluster the network based on interconnectivity using Louvain algorithm
get_cluster <- function(.df, .matrix, .column, .threshold = 0.25){
  # extract the list of enriched gene sets 
  genes <- split(strsplit(.df[["core_enrichment"]],"/"), .df[["ID"]])
  
  # subset the cohen's matrix for the gene sets of interest
  similarity.matrix <- .matrix[.df[['ID']], .df[['ID']]]
  
  # create an adjacency matrix based on the similarity matrix
  adjacency.matrix <- similarity.matrix > .threshold
  
  set.seed(42)
  # create a network graph from the adjacency matrix
  graph <- graph_from_adjacency_matrix(adjacency.matrix, mode = "undirected", diag = F)
  
  # annotate the nodes of the graph with the enrichment results
  V(graph)$Name <- .df$Name
  V(graph)$NES <- .df$NES
  V(graph)$geneRatio <- .df$geneRatio
  
  # Community detection (Louvain algorithm)
  clusters <- cluster_louvain(graph)
  # annotate nodes with cluster membership
  V(graph)$cluster <- membership(clusters)
  
  # extract cluster information
  cluster_summary <- data.frame(
    ID = names(V(graph)),
    cluster = as.factor(V(graph)$cluster))
  
  cluster_summary <- distinct(cluster_summary, .keep_all = T)
  
  df = .df %>% 
    dplyr::inner_join(., cluster_summary, by = "ID")
  
  return(list(graph = graph, df = df))
}

get_cluster_representative <- function(.cluster, .degs){
  # add gene expression values to the clusters
  linkage <- .cluster %>% 
    dplyr::select(ID, Name, core_enrichment, cluster) %>%
    tidyr::separate_rows(core_enrichment, sep = "/")
  
  linkage$logFC = .degs$log2FoldChange[match(linkage$core_enrichment, .degs$Symbol)]
  linkage$padj = .degs$padj[match(linkage$core_enrichment, .degs$Symbol)]
  
  # calculate normalized term weight
  linkage$weight = abs(linkage$logFC)*(-log10(linkage$padj))
  linkage$weight = linkage$weight/max(linkage$weight, na.rm = T)
  
  linkage = linkage %>% 
    dplyr::select(c("core_enrichment", "ID", "weight", "cluster")) %>%
    setNames(.,c("node1", "node2", "weight", "cluster")) %>% 
    as.data.frame(.)
  
  # create a network visualization of gene and GO-term relationships
  net <- graph_from_data_frame(linkage)
  net <- igraph::simplify(net, remove.multiple = F, remove.loops = T)
  
  # calculate hub score of each gene
  hs <-  igraph::hub_score(net, scale = T, weights = linkage$weight)$vector
  
  # summarize gene hub scores, to determine the final weight of each GO term
  linkage$geneRatio = .cluster$geneRatio[match(linkage$node2, .cluster$ID)]
  
  linkage <- linkage %>%
    dplyr::mutate(hub_score = geneRatio * hs[match(node1, names(hs))]) %>% 
    dplyr::rename(geneID = node1, ID = node2) %>% 
    dplyr::group_by(ID, cluster) %>%
    dplyr::reframe(core_enrichment = paste(geneID, collapse = "/"),
                   hub_score = sum(hub_score, na.rm = T))
  
  out_df <- inner_join(.cluster, linkage, by = c("ID","core_enrichment","cluster")) %>%
    dplyr::mutate(cluster = as.factor(cluster))
  
  # select cluster representatives according to calculated weight
  representative.terms <- out_df %>%
    dplyr::group_by(cluster) %>%
    dplyr::arrange(desc(hub_score)) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::pull(ID)
  
  out_df <- out_df %>%
    dplyr::rowwise(.) %>%
    dplyr::mutate(Representative = ifelse(ID %in% representative.terms, T, F))
  
  return(out_df)
}

filter_graph <- function(.graph, .threshold){
  cl <- table(V(.graph)$cluster)
  
  valid_cl <- which(cl >= .threshold)
  
  filtered_vertices <- which(V(.graph)$cluster %in% valid_cl)
  
  subgraph <- induced_subgraph(.graph, vids = filtered_vertices)
  
  return(subgraph)
}

plot_network <- function(.net, .layout, .df, .condition){
  plot_labels = c(paste("ENHANCED IN", .condition),
                  paste("INHIBITED IN", .condition))
  names(plot_labels) <- c("UP", "DOWN")
  
  edges = as.data.frame(as_edgelist(.net)) %>%
    setNames(c("from", "to"))

  df = data.frame(x = as.data.frame(.layout)[,1],
                  y = as.data.frame(.layout)[,2],
                  ID = names(V(.net)),
                  Name = V(.net)$Description,
                  cluster = as.factor(V(.net)$cluster),
                  Representative = V(.net)$Representative,
                  weight = abs(.df[["hub_score"]][match(names(V(.net)), .df[["ID"]])]),
                  regulation = ifelse(.df[["NES"]][match(names(V(.net)), .df[["ID"]])] > 0, "UP", "DOWN"),
                  cluster_name = .df[["Cluster_name"]][match(names(V(.net)), .df[["ID"]])],
                  interest = .df[["Interest"]][match(names(V(.net)), .df[["ID"]])]
  )

  lines = edges %>%
    dplyr::mutate(
      from.x = df$x[match(edges$from, df$ID)],
      from.y = df$y[match(edges$from, df$ID)],
      to.x = df$x[match(edges$to, df$ID)],
      to.y = df$y[match(edges$to, df$ID)]
    )

  plot = ggplot() +
    stat_ellipse(data = df, geom = "polygon",
                 aes(x = x, y = y, group = cluster, fill = interest), color = "grey15",
                 type = "norm", level = 0.9,
                 alpha = .25, linetype = 2, show.legend = F) +
    scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "grey75"),
                      name = c("Interest"), guide = "none") +
    geom_segment(data = lines, aes(x = from.x, y = from.y, xend = to.x, yend = to.y),
                 color = "grey25") +
    new_scale_fill() +
    geom_point(data = df, aes(x = x, y = y, size = weight, fill = regulation),
               shape = 21, colour = "black") +
    scale_size_continuous(name = c("Hub score"), range = c(5,10),
                          limits = c(0, max(df$weight, na.rm = T)),
                          guide = guide_legend(order = 2, 
                                               override.aes = c(fill = "grey"))) +
    scale_fill_manual(values = c("UP" = "darkred", "DOWN" = "blue"),
                      labels = c("UP" = paste0("ENHANCED IN ", .condition),
                                 "DOWN" = paste0("INHIBITED IN ", .condition)),
                      name = c(""),
                      guide = guide_legend(order = 1, override.aes = c(size = 10))) +
    geom_text_repel(data = df %>%
                      dplyr::group_by(cluster) %>%
                      dplyr::summarise(x = mean(x), y = max(y),
                                       cluster_name = first(cluster_name)),
                    aes(x = x, y = y,
                        label = str_wrap(cluster_name, 20)),
                    size = 8, fontface = 2, colour = "black", segment.alpha = 0,
                    box.padding = 0.5, point.padding = 0.5,
                    max.overlaps = Inf, show.legend = F) +
    theme_void() +
    theme(legend.text=element_text(size=28),
          legend.title=element_text(size=30, face="bold"),
          legend.spacing=unit(1.5,"lines"),
          legend.position = "top",
          legend.direction = "vertical",
          legend.justification = "left",
          legend.box.just = "left",
          legend.margin = margin(0, 0, 0, 0))
  
  return(plot)
}

# plot the enriched cluster on a vulcano plot of pathway's
plotClusters <- function(.df, .pathways){
  df <- .df %>% 
    dplyr::left_join(., .pathways, by = c("ID","Name")) %>%
    dplyr::mutate(Category = as.factor(Category))
  df <- df %>% 
    expand_grid(facet_var = unique(df$Category[!is.na(df$Category)]))
  
  return(
    ggplot() 
    + geom_point(data = df, aes(x = NES, y = -log10(p.adjust), fill = NES), 
                 size = 5, shape = 21, show.legend = F,
                 color = "black", alpha = 0.1) 
    # Visualize clusters of interest
    + geom_point(data = subset(df, !is.na(Category) & Category == facet_var),
                 aes(x = NES, y = -log10(p.adjust), colour = Category), 
                 show.legend = F, size = 5, alpha = 0.7)
    + scale_color_manual(values = c("ECM"="#7E03A8FF",
                                    "HYPOXIA"="#F00000FF",
                                    "TGFb"="#09BFF9FF"))
    + scale_fill_gradient2(low = "blue", high = "red",
                           mid = "white", midpoint = 0)
    + facet_grid(~facet_var, scales = "free_x")
    + labs(x = expression('NES'),
           y = expression(paste(log[10], italic('FDR'))))
    + scale_x_continuous(expand = expansion(0.2))
    + theme_minimal() 
    + theme(axis.title = element_text(size = 14), 
            axis.text = element_text(size = 14, color = "#888888"),
            strip.text = element_blank(),
            # Controls spacing between facets
            panel.spacing = unit(5, "lines"),
            legend.position = 'none'))
}

score_plot <- function(.score, .formula, .ref_group, .x, .y, .title, number = F){
  # get formula
  formula <- as.formula(.formula)
  # statistical test
  stat.test = compare_means(formula, data = .score, 
                            method = "wilcox.test", p.adjust.method = "BH",
                            ref.group = .ref_group)
  
  range <- range(.score[[.y]])
  # Compute category counts
  counts <- .score %>%
    group_by(!!sym(.x)) %>%
    summarise(n = n()) %>%
    mutate(label = paste0(!!sym(.x), " (n=", n, ")"))  # Modify labels
  
  plot = ggplot(.score, aes(x=!!sym(.x), y=!!sym(.y))) + 
    #geom_jitter(aes(fill = Subtype), size = 5, alpha = 0.25) +
    geom_boxplot(color = "grey15", fill = "white",outlier.colour = "black",
                 outlier.shape = 21, outlier.stroke = 1.2, linewidth = .9,
                 alpha = 0.5) +
    scale_y_continuous(limits = range, expand = expansion(mult = c(0.1, 0.1)),
                       breaks = round(seq(floor(min(range)),
                                          ceiling(max(range)),
                                          2),1)) + 
    labs(y = str_wrap(.title, width = 25)) +
    theme_classic() + 
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 14, color = "grey25"),
          axis.title.y = element_text(size = 14),
          axis.text.y = element_text(size = 14, color = "grey25")) + 
    stat_pvalue_manual(stat.test, hide.ns = TRUE, label = "p.signif",
                       size = 5,
                       remove.bracket = T, y.position = max(.score[[.y]]) - 0.1)
  
  if (number) {
    plot = plot +
      scale_x_discrete(labels = str_wrap(setNames(counts[["label"]], counts[[.x]]), 5))
  } else {
    plot = plot +
      scale_x_discrete(labels = str_wrap(counts[[.x]], 5))
  }
  
  return(list("plot" = plot, "stat.test" = stat.test))
}
