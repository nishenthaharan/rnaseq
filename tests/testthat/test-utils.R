test_that("provenance manifests contain reproducibility-critical fields", {
  provenance <- list(
    analysed_at = "2026-08-24 12:00:00 UTC",
    r_version = "R version 4.5.1",
    platform = "x86_64-pc-linux-gnu",
    parameters = list(
      comparison = "Treated vs Control",
      alpha = 0.05,
      batch_adjusted = FALSE
    ),
    packages = c(shiny = "1.11.1", DESeq2 = "1.48.1")
  )

  manifest <- provenance_to_text(provenance)

  expect_identical(manifest[1L], "RNASeq Explorer analysis manifest")
  expect_true("[parameters]" %in% manifest)
  expect_true("comparison: Treated vs Control" %in% manifest)
  expect_true("batch_adjusted: false" %in% manifest)
  expect_true("[packages]" %in% manifest)
  expect_true("DESeq2: 1.48.1" %in% manifest)
})

test_that("provenance manifest validation rejects non-list input", {
  expect_error(provenance_to_text("invalid"), "valid provenance object")
})
