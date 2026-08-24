test_that("plot count controls accept only positive integer scalars", {
  expect_identical(validate_plot_count(12, "Labels"), 12L)
  expect_error(validate_plot_count(0, "Labels"), "positive integer")
  expect_error(validate_plot_count(1.5, "Labels"), "positive integer")
  expect_error(validate_plot_count("12", "Labels"), "positive integer")
  expect_error(validate_plot_count(c(10, 20), "Labels"), "positive integer")
})

test_that("top heatmap data follows significance ranking and remains finite", {
  analysis <- list(
    result = data.frame(
      gene_id = c("gene_b", "gene_a", "gene_c"),
      adjusted_p_value = c(0.001, 0.01, NA_real_),
      stringsAsFactors = FALSE
    ),
    transformed = matrix(
      c(
        1, 2, 3,
        5, 5, 5,
        9, 4, 2
      ),
      nrow = 3L,
      byrow = TRUE,
      dimnames = list(
        c("gene_a", "gene_b", "gene_c"),
        c("sample_1", "sample_2", "sample_3")
      )
    )
  )

  heatmap_data <- top_heatmap_matrix(analysis, n_genes = 2L)

  expect_identical(rownames(heatmap_data), c("gene_b", "gene_a"))
  expect_equal(dim(heatmap_data), c(2L, 3L))
  expect_true(all(is.finite(heatmap_data)))
  expect_true(all(abs(heatmap_data) <= 2.5))
  expect_error(top_heatmap_matrix(analysis, n_genes = 0), "positive integer")
})
