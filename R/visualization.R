#' Create volcano plot for differential expression
#'
#' @description Generates a volcano plot highlighting significantly up/down-regulated genes
#'
#' @param .data Data frame with log2FoldChange, padj, Symbol, and significance columns
#' @param label Logical, whether to add gene labels (default: TRUE)
#'
#' @return ggplot object
#' @export
#'
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual labs scale_x_continuous
#'   geom_vline geom_hline theme element_text geom_text expansion
#' @examples
#' \dontrun{
#' volcano_plot <- plot_vulcan(deg_results, label = TRUE)
#' }
plot_vulcan <- function(.data, label = TRUE){
  plot <- ggplot2::ggplot(data = na.omit(.data),
                         ggplot2::aes(x = log2FoldChange, y = -log10(padj), colour = significance)) +
    ggplot2::geom_point(mapping = ggplot2::aes(), inherit.aes = TRUE, size = 2.5, alpha = 0.35) +
    ggplot2::scale_color_manual(values = c("NS" = '#c1c1c1',
                                          "Log10P" = '#363636',
                                          "Log2FoldChange" = '#767676',
                                          "Signif. up-regulated" = '#841f27',
                                          "Signif. down-regulated" = '#000f64'),
                               name = "Significance") +
    ggplot2::labs(x = expression(paste(log[2], 'FoldChange')),
                 y = expression(paste(log[10], italic('FDR')))) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(0.2)) +
    # Visualize log2FC threshold
    ggplot2::geom_vline(xintercept = c(-0.5, 0.5), linetype = 'dotted', size = 1) +
    # Visualize adjusted p-value threshold
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = 'dotted', size = 1) +
    ggplot2::theme(axis.title = ggplot2::element_text(size = 14),
                  axis.text = ggplot2::element_text(size = 14),
                  legend.position = 'none')
  
  if (label == TRUE){
    return(plot +
             # Add labels for significantly up- or down-regulated genes
             ggplot2::geom_text(data = subset(.data,
                                            significance %in% c("Signif. up-regulated",
                                                               "Signif. down-regulated")),
                              hjust = 0, vjust = 1.5, colour = 'black', position = 'identity',
                              show.legend = FALSE, check_overlap = TRUE,
                              label = subset(.data,
                                           significance %in% c("Signif. up-regulated",
                                                              "Signif. down-regulated"))[, "Symbol"]))
  } else {
    return(plot)
  }
}


#' Plot GSEA running enrichment score
#'
#' @description Creates a running enrichment score plot for GSEA results
#'
#' @param .df Data frame from gsInfo with x, runningScore, position columns
#' @param .x Column name for x-axis (gene rank)
#' @param .y Column name for y-axis (running score)
#' @param .color Column name for color grouping
#' @param .palette Named vector of colors for each group
#'
#' @return ggplot object
#' @export
#'
#' @importFrom ggplot2 ggplot aes xlab ylab theme_bw geom_hline scale_x_continuous
#'   scale_y_continuous geom_point scale_color_manual theme element_blank element_text
#' @examples
#' \dontrun{
#' enrich_data <- gsInfo(gsea_result$gsea, 1)
#' plot <- plot_running_score(enrich_data, "x", "runningScore", "Description", palette)
#' }
plot_running_score <- function(.df, .x, .y, .color, .palette){
  subset_data <- subset(.df, position == 1)
  
  return(ggplot2::ggplot(.df, ggplot2::aes(x = !!rlang::sym(.x))) +
           ggplot2::xlab(ggplot2::element_blank()) +
           ggplot2::ylab("Enrichment score (ES)") +
           ggplot2::theme_bw(base_size = 14) +
           ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
           ggplot2::scale_x_continuous(expand = c(0, 0)) +
           ggplot2::scale_y_continuous(expand = c(0.01, 0.01),
                                      breaks = c(-0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5)) +
           ggplot2::geom_point(ggplot2::aes(x = !!rlang::sym(.x), y = !!rlang::sym(.y), color = !!rlang::sym(.color)),
                             show.legend = FALSE, size = 1, data = subset_data) +
           ggplot2::scale_color_manual(values = .palette) +
           ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                         panel.grid.minor = ggplot2::element_blank(),
                         panel.grid.major.y = ggplot2::element_blank(),
                         panel.grid.minor.y = ggplot2::element_blank(),
                         axis.text.x = ggplot2::element_blank(),
                         axis.ticks.x = ggplot2::element_blank(),
                         axis.text.y = ggplot2::element_text(colour = "black", size = 14))
  )
}


#' Plot gene rank for GSEA
#'
#' @description Creates gene rank barcode plot showing where enriched genes appear in ranked list
#'
#' @param .df Data frame from gsInfo with position information
#' @param .x Column name for x-axis (gene rank)
#' @param .facet Formula for faceting (e.g., "Description~.")
#' @param .color Column name for color grouping
#' @param .palette Named vector of colors
#'
#' @return ggplot object
#' @export
#'
#' @importFrom ggplot2 ggplot aes geom_linerange theme_bw xlab scale_x_continuous
#'   scale_y_continuous scale_color_manual facet_grid theme element_blank element_text
#' @examples
#' \dontrun{
#' enrich_data <- gsInfo(gsea_result$gsea, 1)
#' plot <- plot_gene_rank(enrich_data, "x", "Description~.", "Description", palette)
#' }
plot_gene_rank <- function(.df, .x, .facet, .color, .palette){
  subset_data <- subset(.df, position == 1)
  facet <- stats::as.formula(.facet)
  
  return(ggplot2::ggplot(.df, ggplot2::aes(x = !!rlang::sym(.x))) +
           ggplot2::geom_linerange(ggplot2::aes(ymin = ymin, ymax = ymax, color = !!rlang::sym(.color)),
                                 show.legend = FALSE, size = 1, data = subset_data) +
           ggplot2::theme_bw(base_size = 14) +
           ggplot2::xlab("Rank in ordered gene list") +
           ggplot2::scale_x_continuous(expand = c(0, 0)) +
           ggplot2::scale_y_continuous(expand = c(0, 0)) +
           ggplot2::scale_color_manual(values = .palette) +
           ggplot2::facet_grid(facet, scales = "free_y", switch = "y",
                             labeller = ggplot2::labeller(.default = Hmisc::capitalize)) +
           ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                         panel.grid.minor = ggplot2::element_blank(),
                         panel.grid.major.y = ggplot2::element_blank(),
                         panel.grid.minor.y = ggplot2::element_blank(),
                         axis.ticks.y = ggplot2::element_blank(),
                         axis.text.y = ggplot2::element_blank(),
                         strip.placement = "outside",
                         strip.background = ggplot2::element_blank(),
                         strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 0,
                                                                  colour = "grey25", size = 14),
                         axis.text.x = ggplot2::element_text(colour = "black", size = 14))
  )
}

#' Plot enriched clusters on volcano plot
#'
#' @description Creates a volcano plot of enrichment results with highlighted clusters
#'
#' @param .df Data frame with enrichment results and cluster categories
#' @param .pathways Data frame with pathway categories of interest
#'
#' @return ggplot object
#' @export
#'
#' @importFrom ggplot2 ggplot geom_point aes scale_color_manual scale_fill_gradient2
#'   facet_grid labs scale_x_continuous theme_minimal theme element_text element_blank
#' @importFrom dplyr left_join mutate
#' @importFrom tidyr expand_grid
#' @examples
#' \dontrun{
#' plot <- plot_clusters(gsea_results, pathways_of_interest)
#' }
plot_clusters <- function(.df, .pathways){
  df <- .df %>%
    dplyr::left_join(., .pathways, by = c("ID", "Name")) %>%
    dplyr::mutate(Category = as.factor(Category))
  df <- df %>%
    tidyr::expand_grid(facet_var = unique(df$Category[!is.na(df$Category)]))
  
  return(
    ggplot2::ggplot() +
      ggplot2::geom_point(data = df, ggplot2::aes(x = NES, y = -log10(FDR), fill = NES),
                        size = 5, shape = 21, show.legend = FALSE,
                        color = "black", alpha = 0.1) +
      # Visualize clusters of interest
      ggplot2::geom_point(data = subset(df, !is.na(Category) & Category == facet_var),
                        ggplot2::aes(x = NES, y = -log10(FDR), colour = Category),
                        show.legend = FALSE, size = 5, alpha = 0.7) +
      ggplot2::scale_color_manual(values = c("ECM" = "#7E03A8FF",
                                            "HYPOXIA" = "#F00000FF",
                                            "TGFb" = "#09BFF9FF")) +
      ggplot2::scale_fill_gradient2(low = "blue", high = "red",
                                   mid = "white", midpoint = 0) +
      ggplot2::facet_grid(~facet_var, scales = "free_x") +
      ggplot2::labs(x = expression('NES'),
                   y = expression(paste(log[10], italic('FDR')))) +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(0.2)) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.title = ggplot2::element_text(size = 14),
                    axis.text = ggplot2::element_text(size = 14, color = "#888888"),
                    strip.text = ggplot2::element_blank(),
                    panel.spacing = grid::unit(5, "lines"),
                    legend.position = 'none')
  )
}


#' Plot network graph of enriched pathways
#'
#' @description Creates a network visualization of clustered enriched pathways
#'
#' @param .net igraph network object
#' @param .layout Network layout coordinates
#' @param .df Data frame with cluster information
#' @param .condition Name of experimental condition for labeling
#'
#' @return ggplot object
#' @export
#'
#' @importFrom ggplot2 ggplot stat_ellipse geom_segment geom_point
#'   scale_fill_manual scale_size_continuous theme_void theme element_text
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggnewscale new_scale_fill
#' @importFrom igraph as_edgelist V
#' @importFrom dplyr mutate group_by summarise first
#' @importFrom stringr str_wrap
#' @examples
#' \dontrun{
#' network_plot <- plot_network(graph, layout, cluster_df, "HYPOXIA")
#' }
plot_network <- function(.net, .layout, .df, .condition){
  edges <- as.data.frame(igraph::as_edgelist(.net)) %>%
    setNames(c("from", "to"))
  
  df <- data.frame(
    x = as.data.frame(.layout)[, 1],
    y = as.data.frame(.layout)[, 2],
    ID = names(igraph::V(.net)),
    Name = igraph::V(.net)$Description,
    cluster = as.factor(igraph::V(.net)$cluster),
    Representative = igraph::V(.net)$Representative,
    weight = abs(.df[["hub_score"]][match(names(igraph::V(.net)), .df[["ID"]])]),
    regulation = ifelse(.df[["NES"]][match(names(igraph::V(.net)), .df[["ID"]])] > 0, "UP", "DOWN"),
    cluster_name = .df[["Cluster_name"]][match(names(igraph::V(.net)), .df[["ID"]])],
    interest = .df[["Interest"]][match(names(igraph::V(.net)), .df[["ID"]])]
  )
  
  lines <- edges %>%
    dplyr::mutate(
      from.x = df$x[match(edges$from, df$ID)],
      from.y = df$y[match(edges$from, df$ID)],
      to.x = df$x[match(edges$to, df$ID)],
      to.y = df$y[match(edges$to, df$ID)]
    )
  
  plot <- ggplot2::ggplot() +
    ggplot2::stat_ellipse(data = df, geom = "polygon",
                         ggplot2::aes(x = x, y = y, group = cluster, fill = interest),
                         color = "grey15",
                         type = "norm", level = 0.9,
                         alpha = .25, linetype = 2, show.legend = FALSE) +
    ggplot2::scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "grey75"),
                              name = c("Interest"), guide = "none") +
    ggplot2::geom_segment(data = lines,
                         ggplot2::aes(x = from.x, y = from.y, xend = to.x, yend = to.y),
                         color = "grey25") +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_point(data = df, ggplot2::aes(x = x, y = y, size = weight, fill = regulation),
                       shape = 21, colour = "black") +
    ggplot2::scale_size_continuous(name = c("Hub score"), range = c(5, 10),
                                  limits = c(0, max(df$weight, na.rm = TRUE)),
                                  guide = ggplot2::guide_legend(order = 2,
                                                               override.aes = c(fill = "grey"))) +
    ggplot2::scale_fill_manual(values = c("UP" = "darkred", "DOWN" = "blue"),
                              labels = c("UP" = paste0("ENHANCED IN ", .condition),
                                        "DOWN" = paste0("INHIBITED IN ", .condition)),
                              name = c(""),
                              guide = ggplot2::guide_legend(order = 1, override.aes = c(size = 10))) +
    ggrepel::geom_text_repel(data = df %>%
                              dplyr::group_by(cluster) %>%
                              dplyr::summarise(x = mean(x), y = max(y),
                                             cluster_name = dplyr::first(cluster_name)),
                            ggplot2::aes(x = x, y = y,
                                       label = stringr::str_wrap(cluster_name, 20)),
                            size = 8, fontface = 2, colour = "black", segment.alpha = 0,
                            box.padding = 0.5, point.padding = 0.5,
                            max.overlaps = Inf, show.legend = FALSE) +
    ggplot2::theme_void() +
    ggplot2::theme(legend.text = ggplot2::element_text(size = 28),
                  legend.title = ggplot2::element_text(size = 30, face = "bold"),
                  legend.spacing = grid::unit(1.5, "lines"),
                  legend.position = "top",
                  legend.direction = "vertical",
                  legend.justification = "left",
                  legend.box.just = "left",
                  legend.margin = ggplot2::margin(0, 0, 0, 0))
  
  return(plot)
}


#' Create boxplot with statistical testing
#'
#' @description Creates boxplot comparing groups with Wilcoxon test
#'
#' @param .score Data frame with scores to plot
#' @param .formula Formula for statistical test (e.g., "score~group")
#' @param .ref_group Reference group for comparisons
#' @param .x Column name for x-axis (grouping variable)
#' @param .y Column name for y-axis (values to plot)
#' @param .title Y-axis title
#' @param number Logical, whether to include n values in labels (default: FALSE)
#'
#' @return List with plot and statistical test results
#' @export
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot scale_y_continuous labs theme_classic
#'   theme element_blank element_text scale_x_discrete
#' @importFrom ggpubr compare_means stat_pvalue_manual
#' @importFrom dplyr group_by summarise mutate
#' @importFrom stringr str_wrap
#' @examples
#' \dontrun{
#' result <- plot_score(data, "score~condition", "control", "condition", "score", "Score")
#' }
plot_score <- function(.score, .formula, .ref_group, .x, .y, .title, number = FALSE){
  # Get formula
  formula <- stats::as.formula(.formula)
  # Statistical test
  stat.test <- ggpubr::compare_means(formula, data = .score,
                                    method = "wilcox.test", p.adjust.method = "BH",
                                    ref.group = .ref_group)
  
  range <- range(.score[[.y]])
  # Compute category counts
  counts <- .score %>%
    dplyr::group_by(!!rlang::sym(.x)) %>%
    dplyr::summarise(n = dplyr::n()) %>%
    dplyr::mutate(label = paste0(!!rlang::sym(.x), " (n=", n, ")"))
  
  plot <- ggplot2::ggplot(.score, ggplot2::aes(x = !!rlang::sym(.x), y = !!rlang::sym(.y))) +
    ggplot2::geom_boxplot(color = "grey15", fill = "white", outlier.colour = "black",
                         outlier.shape = 21, outlier.stroke = 1.2, linewidth = .9,
                         alpha = 0.5) +
    ggplot2::scale_y_continuous(limits = range, expand = ggplot2::expansion(mult = c(0.1, 0.1)),
                               breaks = round(seq(floor(min(range)),
                                                ceiling(max(range)),
                                                2), 1)) +
    ggplot2::labs(y = stringr::str_wrap(.title, width = 25)) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                  axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 14, color = "grey25"),
                  axis.title.y = ggplot2::element_text(size = 14),
                  axis.text.y = ggplot2::element_text(size = 14, color = "grey25")) +
    ggpubr::stat_pvalue_manual(stat.test, hide.ns = TRUE, label = "p.signif",
                              size = 5,
                              remove.bracket = TRUE, y.position = max(.score[[.y]]) - 0.1)
  
  if (number) {
    plot <- plot +
      ggplot2::scale_x_discrete(labels = stringr::str_wrap(setNames(counts[["label"]], counts[[.x]]), 5))
  } else {
    plot <- plot +
      ggplot2::scale_x_discrete(labels = stringr::str_wrap(counts[[.x]], 5))
  }
  
  return(list("plot" = plot, "stat.test" = stat.test))
}

# Venn diagram functions

# Helper function for recursive merging
merge.rec <- function(.list, ...){
  if(length(.list)==1) return(.list[[1]])
  Recall(c(list(merge(.list[[1]], .list[[2]], ...)), .list[-(1:2)]), ...)
}

#' Generate Venn diagram data object
#'
#' @description Prepares data for Venn diagram plotting by merging DEG sets
#'
#' @param sets List of data frames with differential expression results
#' @param names Character vector of names for each set
#'
#' @return Data frame with merged gene information and set membership
#' @export
#'
#' @importFrom dplyr filter select mutate across rowwise ungroup c_across relocate
#' @importFrom stringr str_detect
#' @examples
#' \dontrun{
#' venn_data <- generate_venn_object(
#'   sets = list(deg1, deg2, deg3),
#'   names = c("Set1", "Set2", "Set3")
#' )
#' }
generate_venn_object <- function(sets, names){
  n <- length(sets)
  
  if (!n %in% 2:4) {
    stop("Only 2-, 3-, or 4-way Venn diagrams are supported.")
  }
  if (length(names) != n) {
    stop("Length of names must match the number of sets.")
  }
  
  list <- lapply(sets, function(x) {
    x %>% 
      dplyr::filter(stringr::str_detect(significance, c("up-regulated|down-regulated"))) %>% 
      dplyr::select(Symbol, log2FoldChange)
  })
  names(list) <- names
  
  log_df <- merge.rec(list, by = c("Symbol"), suffixes = c("",""), all = T) %>% 
    setNames(., c("Symbol", paste("logFC", names, sep = "_")))
  log_df <- log_df %>% dplyr::distinct(Symbol, .keep_all = T)
  
  set_df <- log_df %>%
    dplyr::mutate(across(starts_with("logFC"), ~ifelse(is.na(.), FALSE, TRUE)))
  
  table <- merge(log_df, set_df, suffixes = c("",""), by = "Symbol") %>%
    dplyr::relocate(where(is.logical), .before = where(is.character)) %>% 
    dplyr::relocate(where(is.numeric), .after = where(is.logical))
  
  colnames(table)[1:n] <- names
  table <- table %>% 
      dplyr::rowwise() %>%
      dplyr::mutate(label = paste(names(.)[which(dplyr::c_across(where(is.logical)))], collapse = "-")) %>%
      dplyr::ungroup()
  return(table)
}

# Helper function for 4-way Venn layout
fourway_venn <- function(sets){
  ellipses <- data.frame(
    set = sets,
    x0 = c(3.5, 5, 5, 6.5),
    y0 = c(4.7, 5.7, 5.7, 4.7),
    a  = c(3.5, 3.5, 3.5, 3.5),
    b  = c(2, 1.7, 1.7, 2),
    angle = c(-10, -10, 10, 10),
    row.names = NULL
  )
  
  set_names <- c(
    sets,
    combn(sets, 2, paste, collapse="-"),
    combn(sets, 3, paste, collapse="-"),
    paste(sets, collapse="-")
  )
  
  labels <- data.frame(
    x = c(1.5, 3.5, 6.5, 9, 2.4, 2.4, 5, 5, 7.4, 7.4, 4, 4, 6, 6, 5),
    y = c(5.5, 7.5, 7.5, 5.5, 6.7, 4, 3, 7, 4, 6.7, 6, 3.9, 6, 3.9, 5),
    label = set_names,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  
  return(list(
    circles = ellipses,
    labels = labels
  ))
}

#' Generate Venn diagram layout
#'
#' @description Creates coordinates for circles and labels in a 2-4 way Venn diagram
#'
#' @param sets Character vector of set names
#' @param R Distance of circle centers from origin (default: 1.2)
#' @param r Circle radius (default: 2)
#' @param RR Label radius for set names (default: 2.4)
#' @param d Pairwise intersection label radius (default: 1.6)
#'
#' @return List with circles and labels data frames
#' @export
#'
#' @examples
#' \dontrun{
#' layout <- generate_venn_layout(c("Set1", "Set2", "Set3"))
#' }
generate_venn_layout <- function(sets,
                                 R = 1.2,  # distance of centers from origin
                                 r = 2,    # circle radius
                                 RR = 2.4, # label radius
                                 d = 1.6   # pairwise label radius
                                 ) {
  
  n <- length(sets)
  
  if (!n %in% 2:4) {
    stop("Only 2-, 3-, or 4-way Venn diagrams are supported.")
  }
  
  # ---- Polar helper ----
  polar_to_cart <- function(radius, angle_deg) {
    angle_rad <- angle_deg * pi / 180
    c(
      x = radius * cos(angle_rad),
      y = radius * sin(angle_rad)
    )
  }
  
  # ---- Define angles dynamically ----
  start <- switch(as.character(n),
                  "2" = 0,
                  "3" = 30,
                  "4" = 45)
  
  step <- switch(as.character(n),
                 "2" = 180,
                 "3" = 120,
                 "4" = 90)
  
  angles <- seq(start, by = step, length.out = n)
  
  # ---- Circle centers ----
  centers <- t(sapply(angles, function(a) polar_to_cart(R, a)))
  colnames(centers) <- c("x0","y0")
  
  circles <- data.frame(
    set = sets,
    x0 = centers[,1],
    y0 = centers[,2],
    a = rep(r, n),
    b = rep(r, n),
    angle = rep(0, n),
    row.names = NULL
  )
  
  # ---- SINGLE SET LABELS ----
  single_labels <- t(sapply(angles, function(a) polar_to_cart(RR, a)))
  
  labels_df <- data.frame(
    x = single_labels[,1],
    y = single_labels[,2],
    label = sets,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  
  if (n == 4) {
    return(fourway_venn(sets))
  }
  
  # ---- DOUBLETS ----
  if (n == 2) {
    labels_df <- rbind(
      labels_df,
      data.frame(x = 0, y = 0,
                 label = paste(sets, collapse = "-"),
                 row.names = NULL)
    )
  # ---- TRIPLETS ----
  } else {
    doublet_angles <- seq(start + step/2, by = (-1) * step, length.out = n)
    doublet_combs  <- combn(sets, 2)
    
    for (i in seq_len(ncol(doublet_combs))) {
      pos <- polar_to_cart(d, doublet_angles[i])
      
      labels_df <- rbind(
        labels_df,
        data.frame(
          x = pos[1],
          y = pos[2],
          label = paste(doublet_combs[,i], collapse = "-"),
          row.names = NULL
        )
      )
    }
    
    labels_df <- rbind(
      labels_df,
      data.frame(x = 0, y = 0,
                 label = paste(sets, collapse = "-"),
                 row.names = NULL)
    )
  }
  
  return(list(
    circles = circles,
    labels  = labels_df
  ))
}

#' Create Venn diagram with expression trends
#'
#' @description Generates a Venn diagram showing overlaps between differential expression results
#' with pie charts indicating regulation direction (up/down/opposing)
#'
#' @param sets List of data frames with differential expression results
#' @param names Character vector of names for each set
#' @param R Distance of circle centers from origin (default: 1.2)
#' @param r Circle radius (default: 2)
#' @param RR Label radius for set names (default: 2.4)
#' @param d Pairwise intersection label radius (default: 1.6)
#'
#' @return ggplot object
#' @export
#'
#' @importFrom ggplot2 ggplot aes theme guides guide_legend geom_label scale_fill_manual scale_fill_viridis_d theme element_text unit
#' @importFrom ggforce geom_ellipse
#' @importFrom ggdendro theme_dendro
#' @importFrom ggnewscale new_scale_fill
#' @importFrom scatterpie geom_scatterpie
#' @importFrom dplyr group_by summarise inner_join mutate_all case_when c_across across
#' @importFrom tidyr pivot_wider
#' @examples
#' \dontrun{
#' plot <- plot_venn(
#'   sets = list(deg1, deg2, deg3),
#'   names = c("Set1", "Set2", "Set3")
#' )
#' }
plot_venn <- function(sets, names,
                      R = 1.2, r = 2, RR = 2.4, d = 1.6){
  
  data <- generate_venn_object(sets, names)
  layout <- generate_venn_layout(names, R, r, RR, d)
  
  regions <- data %>%
    dplyr::group_by(label) %>%
    dplyr::summarise(size = n())
  
  regions <- merge(layout$labels, regions, by = "label", all.x = TRUE)
  
  point.matrix <- data %>% 
    dplyr::rowwise(.) %>% 
    # Add jitter to the coordinates for better visualization
    dplyr::mutate(
      # Determine the trend of gene expression
      trend = dplyr::case_when(
        all(dplyr::c_across(starts_with("logFC")) > 0, na.rm = T) ~ 'UP',
        all(dplyr::c_across(starts_with("logFC")) < 0, na.rm = T) ~ 'DOWN',
        TRUE ~ 'CHANGE'
      ),
      trend = as.factor(trend)
    ) %>% 
    dplyr::ungroup(.)
  
  trend.matrix <- point.matrix %>% 
    dplyr::group_by(label, trend) %>% 
    dplyr::summarise(size = n()) %>% 
    dplyr::inner_join(layout$labels, by = "label") %>% 
    tidyr::pivot_wider(names_from = trend, values_from = size) %>% 
    dplyr::mutate_all(~ifelse(is.na(.),0,.))
  
  plot <- ggplot2::ggplot() +
    # Clear the background and axes
    ggdendro::theme_dendro() +
    # Plot the ellipses
    ggforce::geom_ellipse(
      data = layout$circles, # ellipse coordinates
      ggplot2::aes(x0=x0, y0=y0, a=a, b=b, angle=angle, fill = set),
      color = "white", alpha = .4) +
    ggplot2::scale_fill_viridis_d(option = "plasma") +
    ggplot2::guides(fill = ggplot2::guide_legend(title = "Data sets",
                               title.position="top", title.hjust = 0.5)) +
    ggnewscale::new_scale_fill() +
    scatterpie::geom_scatterpie(data = trend.matrix,
                    cols = c("UP","DOWN","CHANGE"), color = "grey25",
                    ggplot2::aes(x=x, y=y, group=label), pie_scale = 3) +
    ggplot2::geom_label(
      data = regions, # Label coordinates
      ggplot2::aes(x = x, y = y, label = size),
      size = 6, alpha = 0.5, fill = 'white', color = 'grey15') +
    # Paint the genes
    ggplot2::scale_fill_manual(
      values = c("UP" = "red",
                 "DOWN" = "blue",
                 "CHANGE" = "orange"),
      labels = c("UP"="Activated",
                 "DOWN"="Inhibited",
                 "CHANGE"="Opposing regulation")) +
    ggplot2::guides(fill = ggplot2::guide_legend(title = "Regulation",
                               title.position="top", title.hjust = 0.5)) +
    # Increase font size and adjust legend placement
    ggplot2::theme(text = ggplot2::element_text(size = 16),
          legend.key.spacing = ggplot2::unit(.5, 'cm'),
          legend.spacing = ggplot2::unit(2, 'cm'),
          legend.position = 'bottom')
  
  return(plot)
}

#' Create two-panel barplot for ORA results
#'
#' @description Generates a 2-panel barplot showing the top significantly enriched pathways
#' for up-regulated and down-regulated genes
#'
#' @param list List containing ORA results with "Signif. up-regulated" and "Signif. down-regulated" elements
#' @param name Character string for the plot title
#' @param n Number of top pathways to display (default: 10)
#' @param terms Character, either "terms" for GO/HALLMARK or "pathways" for REACTOME/KEGG (default: "terms")
#'
#' @return cowplot grid object with two barplots
#' @export
#'
#' @importFrom ggplot2 ggplot aes geom_bar coord_flip labs theme_bw theme element_rect scale_fill_continuous scale_x_discrete scale_y_reverse
#' @importFrom dplyr slice_head mutate
#' @importFrom stringr str_trim str_replace_all str_wrap
#' @importFrom cowplot plot_grid
#' @examples
#' \dontrun{
#' plot <- plot_ora_bar(
#'   list = ora_results,
#'   name = "Comparison_Name",
#'   n = 10,
#'   terms = "terms"
#' )
#' }
plot_ora_bar <- function(list, name, n = 10, terms = c("terms", "pathways")) {
  filter <- ifelse(terms == "terms", "GOBP|HALLMARK|_", "REACTOME|KEGG|_")
  
  df_up <- as.data.frame(list[["Signif. up-regulated"]][["ora"]]) %>%
    dplyr::slice_head(n = n) %>% 
    dplyr::mutate(ID = stringr::str_trim(stringr::str_replace_all(ID, !!filter, " "))) %>% 
    dplyr::mutate(ID = reorder(ID, FoldEnrichment))
  
  df_down <- as.data.frame(list[["Signif. down-regulated"]][["ora"]]) %>%
    dplyr::slice_head(n = n) %>% 
    dplyr::mutate(ID = stringr::str_trim(stringr::str_replace_all(ID, !!filter, " "))) %>% 
    dplyr::mutate(ID = reorder(ID, -FoldEnrichment))
  
  p1 <- ggplot2::ggplot(df_up,
               ggplot2::aes(x = ID,
                   y = FoldEnrichment, fill = p.adjust)) +
    ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge()) +
    ggplot2::scale_fill_continuous(palette =  c("#FC9272", "#DE2D26")) +
    ggplot2::coord_flip() +
    ggplot2::scale_x_discrete(position = "top",
                     labels = function(x) stringr::str_wrap(x, width = 40)) +
    ggplot2::labs(title = "Up-regulated in 3D",
         x = "", fill = "Adj. p-value",
         y = "Fold Enrichment") +
    ggplot2::theme_bw() + 
    ggplot2::theme(legend.position=c(.75,.15),
          legend.background = ggplot2::element_rect(fill = NA))
  
  p2 <- ggplot2::ggplot(df_down,
               ggplot2::aes(x = ID,
                   y = FoldEnrichment, fill = p.adjust)) +
    ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge()) +
    ggplot2::scale_fill_continuous(palette =  c("#9ECAE1", "#3182BD")) +
    ggplot2::scale_y_reverse() +
    ggplot2::scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 40)) +
    ggplot2::coord_flip() +
    ggplot2::labs(title = "Down-regulated in 3D",
         x = "", fill = "Adj. p-value",
         y = "Fold Enrichment") +
    ggplot2::theme_bw() + 
    ggplot2::theme(legend.position=c(.25,.85),
          legend.background = ggplot2::element_rect(fill = NA))
  
  title_gg <- ggplot2::ggplot() + 
    ggplot2::labs(title = stringr::str_replace_all(name, "_", " ")) + 
    ggplot2::theme_minimal()
  gridded <- cowplot::plot_grid(p2, p1, ncol = 2, align = "h")
  
  grid <- cowplot::plot_grid(title_gg, gridded, ncol = 1, rel_heights = c(0.1, 1))
  return(grid)
}
