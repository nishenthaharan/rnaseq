test_that("the demo bundle is balanced, deterministic, and valid", {
  first <- generate_demo_bundle()
  second <- generate_demo_bundle()

  expect_s3_class(first, "rnaseq_bundle")
  expect_identical(first$counts, second$counts)
  expect_equal(dim(first$counts), c(1200L, 8L))
  condition_counts <- table(first$metadata$condition)
  expect_equal(as.integer(condition_counts), c(4L, 4L))
  expect_equal(names(condition_counts), c("Control", "Treated"))
  expect_identical(colnames(first$counts), first$metadata$sample_id)
  expect_true(all(first$counts >= 0L))
})

test_that("bundle validation reorders metadata to the count matrix", {
  counts <- matrix(
    c(10L, 20L, 30L, 40L),
    nrow = 2L,
    dimnames = list(c("gene_a", "gene_b"), c("sample_b", "sample_a"))
  )
  metadata <- data.frame(
    sample_id = c("sample_a", "sample_b"),
    condition = c("Control", "Treated"),
    stringsAsFactors = FALSE
  )

  bundle <- validate_bundle(counts, metadata)
  expect_identical(bundle$metadata$sample_id, c("sample_b", "sample_a"))
})

test_that("bundle validation rejects sample mismatches", {
  counts <- matrix(
    1L,
    nrow = 2L,
    ncol = 2L,
    dimnames = list(c("gene_a", "gene_b"), c("sample_a", "sample_b"))
  )
  metadata <- data.frame(
    sample_id = c("sample_a", "sample_c"),
    condition = c("Control", "Treated"),
    stringsAsFactors = FALSE
  )

  expect_error(validate_bundle(counts, metadata), "do not match exactly")
})

test_that("count-file parsing rejects normalised values", {
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path), add = TRUE)
  utils::write.csv(
    data.frame(gene_id = c("a", "b"), sample_1 = c(1.2, 3), sample_2 = c(2, 4)),
    path,
    row.names = FALSE
  )

  expect_error(read_count_matrix(path, basename(path)), "raw integer counts")
})

test_that("count-file parsing rejects values outside R's integer range", {
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path), add = TRUE)
  utils::write.csv(
    data.frame(
      gene_id = c("a", "b"),
      sample_1 = c(as.double(.Machine$integer.max) + 1, 3),
      sample_2 = c(2, 4)
    ),
    path,
    row.names = FALSE
  )

  expect_error(read_count_matrix(path, basename(path)), "supported integer range")
})

test_that("metadata parsing rejects headers that standardise to blanks", {
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path), add = TRUE)
  metadata <- data.frame(
    sample_id = c("sample_a", "sample_b"),
    condition = c("Control", "Treated"),
    invalid = c("x", "y"),
    check.names = FALSE
  )
  names(metadata)[3L] <- "!!!"
  utils::write.csv(metadata, path, row.names = FALSE)

  expect_error(read_metadata(path, basename(path)), "letter or number")
})
