#' Normalize Affymetrix transcript abundance data
#'
#' @description Normalizes transcript abundance data using the RMA method,
#' which includes background normalization, log2 transformation, and quantile
#' normalization. Annotates results with gene symbols and Entrez IDs.
#'
#' @param .data An AffyBatch object containing raw microarray data
#' @param .db An annotation database (e.g., clariomdhumantranscriptcluster.db)
#'
#' @return A data frame with normalized expression values and gene annotations
#' @export
#'
#' @importFrom oligo rma
#' @importFrom affycoretools annotateEset
#' @importFrom Biobase exprs
#' @examples
#' \dontrun{
#' cel_files <- read.celfiles(list.celfiles("data/", full.names = TRUE))
#' normalized <- normalize_transcript(cel_files, clariomdhumantranscriptcluster.db)
#' }
normalize_transcript <- function(.data, .db){
  # RMA does background normalization, log2 transformation and quantile normalization
  data <- oligo::rma(.data, target = "core")
  # Annotating using info assembled by the Bioconductor Core Team
  data <- affycoretools::annotateEset(data, .db)
  
  # Return expression matrix and add annotation
  df <- as.data.frame(Biobase::exprs(data))
  df$PROBEID <- rownames(df)
  
  # Annotate gene names
  annot <- data@featureData@data
  df <- merge(df, annot, by = "PROBEID")
  
  # Filter the gene expression to contain only values with gene symbol
  idx <- which(is.na(df$SYMBOL))
  df <- df[-idx, ]
  
  return(df)
}


#' Normalize Illumina BeadChip data
#'
#' @description Normalizes Illumina BeadChip data using neqc normalization with
#' control probes. Includes gene symbol updating and annotation mapping.
#'
#' @param .df An EListRaw object from limma containing Illumina expression data
#' @param .dbconn Database connection for gene symbol aliases
#' @param .samples Character vector of sample names to keep (default: all samples)
#'
#' @return A data frame with normalized expression values and gene annotations
#' @export
#'
#' @importFrom limma neqc
#' @importFrom AnnotationDbi mapIds
#' @importFrom dplyr select filter rowwise mutate
#' @examples
#' \dontrun{
#' # Read Illumina data
#' targets <- readTargets("targets.txt")
#' x <- read.ilmn(files = "data.txt", ctrlfiles = "controls.txt")
#' # Normalize
#' normalized <- normalize_illumina(x, alias_db)
#' }
normalize_illumina <- function(.df, .dbconn, .samples = colnames(.df$E)){
  x <- illuminaHumanv4ENTREZID
  
  # Normalization and background correction
  df <- limma::neqc(.df)
  # Keep only selected samples
  df <- df[, .samples]
  # Keep probes that are expressed in at least three arrays according to a detection p-value of 5%
  expressed <- rowSums(df$other$Detection < 0.05) >= 3
  df <- df[expressed, ]
  
  # Get Illumina annotations
  annot <- df$genes
  annot <- cbind(
    PROBEID = rownames(annot),
    ENTREZID = sapply(as.list(x[rownames(annot)]), "[[", 1),
    annot
  )
  annot <- annot %>%
    dplyr::select(!TargetID) %>%
    dplyr::filter(SYMBOL != "" & !is.na(SYMBOL)) %>%
    dplyr::rowwise() %>%
    # Replace outdated gene symbols with the most recent ones
    dplyr::mutate(
      SYMBOL = ifelse(
        ENTREZID == "NA" | is.na(ENTREZID),
        yes = switch_alias(.dbconn, SYMBOL),
        no = SYMBOL
      )
    )
  
  # Map Illumina annotations against the org.Hs.eg.db database
  entrez <- tryMapID(annot$SYMBOL, "SYMBOL", "ENTREZID")
  # Refine annotations missing from the Illumina database
  annot <- annot %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      ENTREZID = ifelse(
        ENTREZID == "NA" | is.na(ENTREZID),
        yes = entrez[[SYMBOL]],
        no = ENTREZID
      )
    )
  
  # Add annotation to the expression matrix
  df <- as.data.frame(df$E)
  df$PROBEID <- rownames(df)
  df <- merge(annot, df, by = "PROBEID")
  # Remove rows with missing ENTREZID
  idx <- which(is.na(df$ENTREZID))
  df <- df[-idx, ]
  return(df)
}


#' Try to map gene identifiers with error handling
#'
#' @description Helper function to safely map gene IDs between different formats
#'
#' @param ID Character vector of gene identifiers
#' @param inp Input ID type (keytype)
#' @param outp Output ID type (column)
#'
#' @return Character vector of mapped IDs or NA
#' @keywords internal
tryMapID <- function(ID, inp, outp){
  tryCatch(
    {
      suppressWarnings(
        AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = ID, column = outp, keytype = inp, multiVals = "first")
      )
    },
    error = function(e) {
      message(conditionMessage(e))
      return(NA)
    }
  )
}


#' Switch outdated gene symbol aliases to current symbols
#'
#' @description Helper function to update gene symbols using alias database
#'
#' @param dbconn Database connection with alias information
#' @param ID Gene symbol to check
#'
#' @return Updated gene symbol or original ID
#' @keywords internal
switch_alias <- function(dbconn, ID){
  alias <- as.character(dbconn[which(dbconn$alias == ID), ]$symbol)[1]
  
  if (is.na(alias)){
    return(ID)
  } else {
    return(alias)
  }
}


#' Remove duplicate gene symbols
#'
#' @description Removes duplicate gene symbols from data frame, keeping the probe
#' with the highest absolute value in the specified column (typically log2FC or t-statistic)
#'
#' @param .data Data frame containing gene expression data
#' @param .column Column name to use for selecting best probe (e.g., "log2FoldChange", "t")
#' @param .symbol Column name containing gene symbols
#'
#' @return Data frame with duplicates removed
#' @export
#'
#' @examples
#' \dontrun{
#' # Remove duplicates based on log2 fold change
#' unique_genes <- remove_duplicates(expr_data, "log2FoldChange", "SYMBOL")
#' }
remove_duplicates <- function(.data, .column, .symbol){
  # Sort the data frame based on the specified values
  data <- .data[order(abs(.data[[.column]]), decreasing = TRUE), ]
  # Remove duplicates
  data <- data[!duplicated(data[[.symbol]]), ]
  return(data)
}


#' Differential expression analysis with limma
#'
#' @description Performs differential expression analysis using the limma package
#' with customizable design and contrasts
#'
#' @param .data Data frame containing expression matrix with numeric columns
#' @param .design Character vector specifying sample groups
#' @param .contrast Character vector of contrasts to test (e.g., "treatment-control")
#'
#' @return List of data frames, one for each contrast, containing differential expression results
#' @export
#'
#' @importFrom limma lmFit eBayes topTable makeContrasts contrasts.fit
#' @importFrom stats model.matrix
#' @examples
#' \dontrun{
#' design <- c("control", "control", "treatment", "treatment")
#' contrasts <- c("treatment-control")
#' results <- limma_dea(expr_data, design, contrasts)
#' }
limma_dea <- function(.data, .design, .contrast){
  # Extract numeric columns
  mat <- .data[, sapply(.data, is.numeric)]
  colnames(mat) <- gsub("_\\(.*$", "\\1", colnames(mat))
  rownames(mat) <- .data$PROBEID
  
  # Create a design matrix
  t <- as.factor(.design)
  levels(t) <- make.names(levels(t)) # Rename levels to ensure syntactically valid names
  
  design <- stats::model.matrix(~0 + t)
  colnames(design) <- levels(t)
  
  # Fit the linear model
  fit <- limma::lmFit(mat, design)
  fit$genes$ID <- rownames(mat)
  fit$genes$Symbol <- .data$SYMBOL
  fit$genes$entrezID <- .data$ENTREZID
  
  # Fit the contrasts
  contrast.matrix <- limma::makeContrasts(contrasts = .contrast, levels = t)
  fit2 <- limma::contrasts.fit(fit, contrast.matrix)
  fit2 <- limma::eBayes(fit2)
  
  # Make DEG table with topTable
  res <- list()
  for (i in 1:ncol(contrast.matrix)){
    res[[i]] <- limma::topTable(fit2, coef = i, number = Inf)
  }
  
  return(res)
}


#' Set significance levels for differential expression
#'
#' @description Categorizes genes based on log2 fold change and adjusted p-value thresholds
#'
#' @param .df Data frame with log2FoldChange and padj columns
#' @param fc_threshold Fold change threshold (default: 0.5)
#' @param padj_threshold Adjusted p-value threshold (default: 0.05)
#'
#' @return Data frame with added 'significance' column
#' @export
#'
#' @importFrom dplyr mutate case_when filter
#' @examples
#' \dontrun{
#' deg_results <- get_significance(deg_data)
#' }
get_significance <- function(.df, fc_threshold = 0.5, padj_threshold = 0.05){
  return(.df %>%
    dplyr::mutate(
      significance = dplyr::case_when(
        abs(log2FoldChange) > fc_threshold & padj > padj_threshold ~ 'log2FoldChange',
        abs(log2FoldChange) < fc_threshold & padj < padj_threshold ~ 'Log10P',
        log2FoldChange < (-1) * fc_threshold & padj < padj_threshold ~ 'Signif. down-regulated',
        log2FoldChange > fc_threshold & padj < padj_threshold ~ 'Signif. up-regulated',
        TRUE ~ 'NS'
      )
    ) %>%
    dplyr::filter(complete.cases(.)))
}
