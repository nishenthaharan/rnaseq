options(repos = c(CRAN = "https://cloud.r-project.org"))

cran_packages <- c(
  "shiny",
  "bslib",
  "DT",
  "ggplot2",
  "ggrepel",
  "pheatmap",
  "scales",
  "testthat",
  "BiocManager"
)

missing_cran <- cran_packages[
  !vapply(cran_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]
if (length(missing_cran)) {
  install.packages(missing_cran, dependencies = TRUE)
}

bioconductor_packages <- c("DESeq2", "SummarizedExperiment")
missing_bioconductor <- bioconductor_packages[
  !vapply(bioconductor_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]
if (length(missing_bioconductor)) {
  BiocManager::install(missing_bioconductor, ask = FALSE, update = FALSE)
}

message("RNASeq Explorer dependencies are installed.")
