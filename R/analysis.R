validate_analysis_parameters <- function(bundle, reference, target, use_batch) {
  conditions <- unique(bundle$metadata$condition)
  if (!reference %in% conditions || !target %in% conditions) {
    abort_user("The selected comparison is not present in the metadata.")
  }
  if (identical(reference, target)) {
    abort_user("Reference and target conditions must be different.")
  }

  selected_counts <- table(bundle$metadata$condition)[c(reference, target)]
  if (any(selected_counts < 2L)) {
    abort_user("Each selected condition needs at least two biological replicates; three or more are strongly recommended.")
  }

  if (isTRUE(use_batch) && !"batch" %in% names(bundle$metadata)) {
    abort_user("Batch adjustment was requested, but the metadata has no batch column.")
  }
  if (isTRUE(use_batch) && length(unique(bundle$metadata$batch)) < 2L) {
    abort_user("Batch adjustment requires at least two observed batch levels.")
  }
}

validate_thresholds <- function(min_count, min_samples, alpha, lfc_threshold, n_samples) {
  values <- c(min_count, min_samples, alpha, lfc_threshold)
  if (length(values) != 4L || anyNA(values) || any(!is.finite(values))) {
    abort_user("All statistical thresholds must be finite numeric values.")
  }
  if (min_count < 0 || min_count != as.integer(min_count)) {
    abort_user("Minimum raw count must be a non-negative integer.")
  }
  if (min_samples < 1 || min_samples > n_samples || min_samples != as.integer(min_samples)) {
    abort_user(sprintf("Minimum samples must be an integer between 1 and %d.", n_samples))
  }
  if (alpha <= 0 || alpha >= 1) {
    abort_user("The adjusted P-value threshold must be greater than 0 and less than 1.")
  }
  if (lfc_threshold < 0) {
    abort_user("The absolute log2 fold-change threshold cannot be negative.")
  }
}

check_design_matrix <- function(metadata, design_formula) {
  matrix <- stats::model.matrix(design_formula, data = metadata)
  if (qr(matrix)$rank < ncol(matrix)) {
    abort_user(
      "The statistical design is not full rank. Condition and batch may be perfectly confounded; disable batch adjustment or revise the metadata."
    )
  }
  invisible(matrix)
}

run_deseq_analysis <- function(
    bundle,
    reference,
    target,
    min_count = 10L,
    min_samples = 3L,
    alpha = 0.05,
    lfc_threshold = 1,
    use_batch = FALSE) {
  validate_analysis_parameters(bundle, reference, target, use_batch)
  validate_thresholds(
    min_count,
    min_samples,
    alpha,
    lfc_threshold,
    n_samples = ncol(bundle$counts)
  )

  metadata <- bundle$metadata
  metadata$condition <- stats::relevel(factor(metadata$condition), ref = reference)
  if (isTRUE(use_batch)) metadata$batch <- factor(metadata$batch)

  design_formula <- if (isTRUE(use_batch)) ~ batch + condition else ~ condition
  check_design_matrix(metadata, design_formula)

  keep <- rowSums(bundle$counts >= min_count) >= min_samples
  if (sum(keep) < 50L) {
    abort_user(sprintf(
      "Only %d genes pass the current filter. Lower the minimum count or number of samples.",
      sum(keep)
    ))
  }

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = bundle$counts[keep, , drop = FALSE],
    colData = metadata,
    design = design_formula
  )
  dds <- DESeq2::DESeq(dds, quiet = TRUE)

  contrast <- c("condition", target, reference)
  raw_result <- DESeq2::results(dds, contrast = contrast, alpha = alpha)
  shrunk_result <- tryCatch(
    DESeq2::lfcShrink(dds, contrast = contrast, res = raw_result, type = "normal"),
    error = function(error) raw_result
  )
  shrinkage_applied <- !identical(shrunk_result, raw_result)

  result <- data.frame(
    gene_id = rownames(shrunk_result),
    base_mean = shrunk_result$baseMean,
    log2_fold_change = shrunk_result$log2FoldChange,
    lfc_se = shrunk_result$lfcSE,
    statistic = raw_result$stat,
    p_value = raw_result$pvalue,
    adjusted_p_value = raw_result$padj,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  result$classification <- "Not significant"
  result$classification[
    !is.na(result$adjusted_p_value) &
      result$adjusted_p_value < alpha &
      result$log2_fold_change >= lfc_threshold
  ] <- "Upregulated"
  result$classification[
    !is.na(result$adjusted_p_value) &
      result$adjusted_p_value < alpha &
      result$log2_fold_change <= -lfc_threshold
  ] <- "Downregulated"
  result <- result[order(result$adjusted_p_value, -abs(result$log2_fold_change), na.last = TRUE), ]

  transformed <- DESeq2::varianceStabilizingTransformation(dds, blind = FALSE)
  transformed_matrix <- SummarizedExperiment::assay(transformed)
  normalised_counts <- DESeq2::counts(dds, normalized = TRUE)

  pca_fit <- stats::prcomp(t(transformed_matrix), center = TRUE, scale. = FALSE)
  variance <- pca_fit$sdev^2 / sum(pca_fit$sdev^2)
  pca <- data.frame(
    sample_id = rownames(pca_fit$x),
    PC1 = pca_fit$x[, 1L],
    PC2 = pca_fit$x[, 2L],
    condition = metadata[rownames(pca_fit$x), "condition"],
    stringsAsFactors = FALSE
  )
  if ("batch" %in% names(metadata)) {
    pca$batch <- metadata[rownames(pca_fit$x), "batch"]
  }

  parameters <- list(
    source = bundle$source,
    comparison = sprintf("%s vs %s", target, reference),
    reference = reference,
    target = target,
    min_count = min_count,
    min_samples = min_samples,
    alpha = alpha,
    absolute_log2_fc = lfc_threshold,
    batch_adjusted = isTRUE(use_batch),
    design = paste(deparse(design_formula), collapse = " "),
    genes_before_filter = nrow(bundle$counts),
    genes_after_filter = sum(keep),
    lfc_shrinkage = if (shrinkage_applied) "DESeq2 normal prior" else "not applied"
  )

  structure(
    list(
      dds = dds,
      result = result,
      normalised_counts = normalised_counts,
      transformed = transformed_matrix,
      pca = pca,
      pca_variance = variance,
      metadata = metadata,
      raw_counts = bundle$counts,
      parameters = parameters,
      provenance = build_provenance(parameters)
    ),
    class = "rnaseq_analysis"
  )
}

significant_results <- function(analysis) {
  analysis$result[analysis$result$classification != "Not significant", , drop = FALSE]
}

analysis_metrics <- function(analysis) {
  result <- analysis$result
  data.frame(
    metric = c("Samples", "Genes tested", "Significant", "Upregulated", "Downregulated"),
    value = c(
      ncol(analysis$normalised_counts),
      nrow(result),
      sum(result$classification != "Not significant"),
      sum(result$classification == "Upregulated"),
      sum(result$classification == "Downregulated")
    ),
    stringsAsFactors = FALSE
  )
}
