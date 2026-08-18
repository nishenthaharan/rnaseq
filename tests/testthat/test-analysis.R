test_that("DESeq2 analysis returns consistent, classified results", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("SummarizedExperiment")

  bundle <- generate_demo_bundle(n_genes = 300L)
  analysis <- run_deseq_analysis(
    bundle,
    reference = "Control",
    target = "Treated",
    min_count = 5L,
    min_samples = 2L,
    alpha = 0.1,
    lfc_threshold = 0.75,
    use_batch = TRUE
  )

  expect_s3_class(analysis, "rnaseq_analysis")
  expect_named(
    analysis$result,
    c(
      "gene_id", "base_mean", "log2_fold_change", "lfc_se", "statistic",
      "p_value", "adjusted_p_value", "classification"
    )
  )
  expect_identical(colnames(analysis$normalised_counts), bundle$metadata$sample_id)
  expect_true(all(analysis$result$classification %in% c(
    "Upregulated", "Downregulated", "Not significant"
  )))
  expect_equal(nrow(analysis$pca), nrow(bundle$metadata))
  expect_equal(sum(analysis$pca_variance), 1, tolerance = 1e-8)
})

test_that("rank-deficient designs are rejected", {
  bundle <- generate_demo_bundle(n_genes = 200L)
  bundle$metadata$batch <- bundle$metadata$condition

  expect_error(
    run_deseq_analysis(
      bundle,
      reference = "Control",
      target = "Treated",
      min_count = 5L,
      min_samples = 2L,
      use_batch = TRUE
    ),
    "not full rank"
  )
})
