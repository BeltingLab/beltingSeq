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
