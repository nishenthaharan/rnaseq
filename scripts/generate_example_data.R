source(file.path("R", "utils.R"))
source(file.path("R", "io.R"))
source(file.path("R", "demo_data.R"))

bundle <- generate_demo_bundle()
dir.create("example_data", showWarnings = FALSE, recursive = TRUE)

count_table <- data.frame(
  gene_id = rownames(bundle$counts),
  bundle$counts,
  check.names = FALSE
)

utils::write.csv(
  count_table,
  file.path("example_data", "counts_example.csv"),
  row.names = FALSE,
  quote = FALSE
)
utils::write.csv(
  bundle$metadata,
  file.path("example_data", "metadata_example.csv"),
  row.names = FALSE,
  quote = FALSE
)

message("Wrote example_data/counts_example.csv and example_data/metadata_example.csv")
