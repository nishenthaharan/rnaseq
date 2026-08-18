plot_library_sizes <- function(analysis) {
  data <- data.frame(
    sample_id = colnames(analysis$raw_counts),
    reads = colSums(analysis$raw_counts),
    condition = analysis$metadata[colnames(analysis$raw_counts), "condition"],
    stringsAsFactors = FALSE
  )

  ggplot2::ggplot(data, ggplot2::aes(x = stats::reorder(sample_id, reads), y = reads, fill = condition)) +
    ggplot2::geom_col(width = 0.72) +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(labels = scales::label_number_si()) +
    ggplot2::scale_fill_manual(values = condition_palette(data$condition)) +
    ggplot2::labs(
      title = "Raw library sizes",
      subtitle = "Total aligned reads represented in the count matrix",
      x = NULL,
      y = "Total counts",
      fill = "Condition"
    ) +
    rnaseq_theme()
}

plot_pca <- function(analysis) {
  data <- analysis$pca
  colours <- condition_palette(data$condition)
  x_label <- sprintf("PC1 (%s)", format_percent(analysis$pca_variance[1L]))
  y_label <- sprintf("PC2 (%s)", format_percent(analysis$pca_variance[2L]))

  ggplot2::ggplot(data, ggplot2::aes(x = PC1, y = PC2, colour = condition)) +
    ggplot2::geom_hline(yintercept = 0, colour = "#DDE3EA", linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = 0, colour = "#DDE3EA", linewidth = 0.4) +
    ggplot2::geom_point(size = 4, alpha = 0.9) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = sample_id),
      size = 3.4,
      max.overlaps = Inf,
      show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(values = colours) +
    ggplot2::labs(
      title = "Principal component analysis",
      subtitle = "Variance-stabilised expression across samples",
      x = x_label,
      y = y_label,
      colour = "Condition"
    ) +
    rnaseq_theme()
}

plot_sample_correlation <- function(analysis) {
  correlation <- stats::cor(analysis$transformed, method = "pearson")
  annotation <- data.frame(condition = analysis$metadata[colnames(correlation), "condition"])
  rownames(annotation) <- colnames(correlation)

  pheatmap::pheatmap(
    correlation,
    annotation_col = annotation,
    annotation_row = annotation,
    color = grDevices::colorRampPalette(c("#0B3C5D", "#D9EAF2", "#FFFFFF", "#F6B48F", "#B11F2E"))(100),
    border_color = NA,
    main = "Sample correlation",
    fontsize = 10,
    angle_col = 45
  )
}

plot_volcano <- function(analysis, label_count = 12L) {
  data <- analysis$result
  data$minus_log10_padj <- -log10(pmax(data$adjusted_p_value, .Machine$double.xmin, na.rm = FALSE))
  data$minus_log10_padj[is.na(data$minus_log10_padj)] <- 0
  label_candidates <- data[data$classification != "Not significant", , drop = FALSE]
  label_candidates <- utils::head(label_candidates, label_count)
  thresholds <- analysis$parameters

  ggplot2::ggplot(data, ggplot2::aes(x = log2_fold_change, y = minus_log10_padj, colour = classification)) +
    ggplot2::geom_point(alpha = 0.58, size = 1.7) +
    ggplot2::geom_vline(
      xintercept = c(-thresholds$absolute_log2_fc, thresholds$absolute_log2_fc),
      linetype = "dashed",
      colour = "#7A8699"
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(thresholds$alpha),
      linetype = "dashed",
      colour = "#7A8699"
    ) +
    ggrepel::geom_text_repel(
      data = label_candidates,
      ggplot2::aes(label = gene_id),
      size = 3,
      colour = "#182230",
      max.overlaps = Inf,
      box.padding = 0.4,
      min.segment.length = 0,
      show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(values = c(
      "Downregulated" = "#276FBF",
      "Not significant" = "#B8C2CC",
      "Upregulated" = "#D1495B"
    )) +
    ggplot2::labs(
      title = "Volcano plot",
      subtitle = sprintf("%s; dashed lines show decision thresholds", thresholds$comparison),
      x = "Shrunken log2 fold change",
      y = expression(-log[10](adjusted ~ italic(P))),
      colour = NULL
    ) +
    rnaseq_theme()
}

plot_ma <- function(analysis) {
  data <- analysis$result

  ggplot2::ggplot(data, ggplot2::aes(x = base_mean, y = log2_fold_change, colour = classification)) +
    ggplot2::geom_hline(yintercept = 0, colour = "#6B778C", linewidth = 0.5) +
    ggplot2::geom_point(alpha = 0.55, size = 1.6) +
    ggplot2::scale_x_log10(labels = scales::label_number_si()) +
    ggplot2::scale_colour_manual(values = c(
      "Downregulated" = "#276FBF",
      "Not significant" = "#B8C2CC",
      "Upregulated" = "#D1495B"
    )) +
    ggplot2::labs(
      title = "MA plot",
      subtitle = analysis$parameters$comparison,
      x = "Mean normalised expression",
      y = "Shrunken log2 fold change",
      colour = NULL
    ) +
    rnaseq_theme()
}

plot_dispersion <- function(analysis) {
  DESeq2::plotDispEsts(analysis$dds, main = "DESeq2 dispersion estimates")
}

top_heatmap_matrix <- function(analysis, n_genes = 30L) {
  ranked <- analysis$result[!is.na(analysis$result$adjusted_p_value), , drop = FALSE]
  if (nrow(ranked)) {
    genes <- utils::head(ranked$gene_id, min(n_genes, nrow(ranked)))
  } else {
    variance <- apply(analysis$transformed, 1L, stats::var)
    genes <- utils::head(names(sort(variance, decreasing = TRUE)), n_genes)
  }
  matrix <- analysis$transformed[genes, , drop = FALSE]
  matrix <- t(scale(t(matrix)))
  matrix[!is.finite(matrix)] <- 0
  pmax(-2.5, pmin(2.5, matrix))
}

plot_top_heatmap <- function(analysis, n_genes = 30L) {
  matrix <- top_heatmap_matrix(analysis, n_genes)
  annotation <- data.frame(condition = analysis$metadata[colnames(matrix), "condition"])
  if (isTRUE(analysis$parameters$batch_adjusted)) {
    annotation$batch <- analysis$metadata[colnames(matrix), "batch"]
  }
  rownames(annotation) <- colnames(matrix)

  pheatmap::pheatmap(
    matrix,
    annotation_col = annotation,
    color = grDevices::colorRampPalette(c("#2155A3", "#F7F7F7", "#C6263E"))(100),
    border_color = NA,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    fontsize_row = 8,
    main = sprintf("Top %d genes (row z-score)", nrow(matrix)),
    angle_col = 45
  )
}

plot_gene_expression <- function(analysis, gene_id) {
  if (!gene_id %in% rownames(analysis$normalised_counts)) {
    abort_user("The selected gene is not available after filtering.")
  }

  data <- data.frame(
    sample_id = colnames(analysis$normalised_counts),
    normalised_count = as.numeric(analysis$normalised_counts[gene_id, ]),
    condition = analysis$metadata[colnames(analysis$normalised_counts), "condition"],
    stringsAsFactors = FALSE
  )

  ggplot2::ggplot(data, ggplot2::aes(x = condition, y = normalised_count + 1, fill = condition)) +
    ggplot2::geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.72) +
    ggplot2::geom_point(
      ggplot2::aes(colour = condition),
      position = ggplot2::position_jitter(width = 0.08),
      size = 2.8,
      show.legend = FALSE
    ) +
    ggplot2::scale_y_log10(labels = scales::label_number_si()) +
    ggplot2::scale_fill_manual(values = condition_palette(data$condition)) +
    ggplot2::scale_colour_manual(values = condition_palette(data$condition)) +
    ggplot2::labs(
      title = gene_id,
      subtitle = "DESeq2 size-factor-normalised counts",
      x = NULL,
      y = "Normalised count + 1 (log10)",
      fill = "Condition"
    ) +
    rnaseq_theme()
}
