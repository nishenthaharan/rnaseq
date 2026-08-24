read_tabular_file <- function(path, original_name) {
  extension <- tolower(tools::file_ext(original_name))

  if (!extension %in% c("csv", "tsv", "txt")) {
    abort_user("Use a .csv, .tsv, or tab-delimited .txt file.")
  }

  reader <- if (identical(extension, "csv")) utils::read.csv else utils::read.delim

  tryCatch(
    reader(
      path,
      header = TRUE,
      check.names = FALSE,
      stringsAsFactors = FALSE,
      comment.char = "",
      quote = "\""
    ),
    error = function(error) {
      abort_user(sprintf("Could not read '%s': %s", original_name, error$message))
    }
  )
}

standardise_metadata_names <- function(x) {
  x <- trimws(x)
  x <- tolower(gsub("[^A-Za-z0-9]+", "_", x))
  gsub("^_+|_+$", "", x)
}

read_count_matrix <- function(path, original_name) {
  table <- read_tabular_file(path, original_name)

  if (ncol(table) < 3L) {
    abort_user("The count matrix needs a gene-ID column and at least two sample columns.")
  }
  if (nrow(table) < 2L) {
    abort_user("The count matrix contains fewer than two genes.")
  }

  gene_id <- trimws(as.character(table[[1L]]))
  if (anyNA(gene_id) || any(!nzchar(gene_id))) {
    abort_user("Every row in the first count-matrix column must have a gene ID.")
  }
  if (anyDuplicated(gene_id)) {
    examples <- utils::head(unique(gene_id[duplicated(gene_id)]), 3L)
    abort_user(sprintf(
      "Gene IDs must be unique. Duplicated example(s): %s.",
      paste(examples, collapse = ", ")
    ))
  }

  sample_names <- trimws(names(table)[-1L])
  if (any(!nzchar(sample_names)) || anyDuplicated(sample_names)) {
    abort_user("Sample column names must be present and unique.")
  }

  converted <- lapply(table[-1L], function(column) suppressWarnings(as.numeric(column)))
  counts <- do.call(cbind, converted)
  dimnames(counts) <- list(gene_id, sample_names)

  if (anyNA(counts) || any(!is.finite(counts))) {
    abort_user("All expression values must be finite numeric counts; missing values are not allowed.")
  }
  if (any(counts < 0)) {
    abort_user("Raw RNA-seq counts cannot be negative.")
  }
  if (any(counts > .Machine$integer.max)) {
    abort_user("Count values exceed R's supported integer range.")
  }
  if (any(abs(counts - round(counts)) > sqrt(.Machine$double.eps))) {
    abort_user("DESeq2 requires raw integer counts. TPM, FPKM, and other normalised values are not valid inputs.")
  }

  storage.mode(counts) <- "integer"
  counts
}

read_metadata <- function(path, original_name) {
  metadata <- read_tabular_file(path, original_name)
  names(metadata) <- standardise_metadata_names(names(metadata))

  if (any(!nzchar(names(metadata)))) {
    abort_user("Metadata column names must contain at least one letter or number.")
  }
  if (anyDuplicated(names(metadata))) {
    abort_user("Metadata column names must remain unique after standardisation.")
  }
  required <- c("sample_id", "condition")
  missing <- setdiff(required, names(metadata))
  if (length(missing)) {
    abort_user(sprintf(
      "Metadata is missing required column(s): %s.",
      paste(missing, collapse = ", ")
    ))
  }

  metadata$sample_id <- trimws(as.character(metadata$sample_id))
  metadata$condition <- trimws(as.character(metadata$condition))

  if (anyNA(metadata$sample_id) || any(!nzchar(metadata$sample_id))) {
    abort_user("Every metadata row must have a sample_id.")
  }
  if (anyDuplicated(metadata$sample_id)) {
    abort_user("Metadata sample_id values must be unique.")
  }
  if (anyNA(metadata$condition) || any(!nzchar(metadata$condition))) {
    abort_user("Every metadata row must have a condition.")
  }
  if (length(unique(metadata$condition)) < 2L) {
    abort_user("At least two experimental conditions are required.")
  }

  if ("batch" %in% names(metadata)) {
    metadata$batch <- trimws(as.character(metadata$batch))
    if (anyNA(metadata$batch) || any(!nzchar(metadata$batch))) {
      abort_user("The optional batch column cannot contain missing or blank values.")
    }
  }

  metadata
}

validate_bundle <- function(counts, metadata, source_label = "Uploaded data") {
  # Uploaded files are validated by their readers, but this constructor is
  # also part of the programmatic API. Guard its structural contract here so
  # downstream DESeq2 failures remain specific and actionable.
  if (!is.matrix(counts) || !is.numeric(counts) || nrow(counts) < 1L || ncol(counts) < 2L) {
    abort_user("Counts must be a numeric matrix with at least one gene and two samples.")
  }
  if (
    is.null(rownames(counts)) ||
      anyNA(rownames(counts)) ||
      any(!nzchar(rownames(counts))) ||
      anyDuplicated(rownames(counts))
  ) {
    abort_user("Count-matrix gene IDs must be present and unique.")
  }
  if (
    is.null(colnames(counts)) ||
      anyNA(colnames(counts)) ||
      any(!nzchar(colnames(counts))) ||
      anyDuplicated(colnames(counts))
  ) {
    abort_user("Count-matrix sample names must be present and unique.")
  }
  if (anyNA(counts) || any(!is.finite(counts))) {
    abort_user("Count-matrix values must be finite and cannot be missing.")
  }
  if (
    any(counts < 0) ||
      any(counts > .Machine$integer.max) ||
      any(counts != round(counts))
  ) {
    abort_user("Count-matrix values must be non-negative integers within R's supported range.")
  }

  if (!is.data.frame(metadata)) {
    abort_user("Metadata must be supplied as a data frame.")
  }
  required_metadata <- c("sample_id", "condition")
  missing_metadata_columns <- setdiff(required_metadata, names(metadata))
  if (length(missing_metadata_columns)) {
    abort_user(sprintf(
      "Metadata is missing required column(s): %s.",
      paste(missing_metadata_columns, collapse = ", ")
    ))
  }

  metadata$sample_id <- trimws(as.character(metadata$sample_id))
  metadata$condition <- trimws(as.character(metadata$condition))
  if (
    anyNA(metadata$sample_id) ||
      any(!nzchar(metadata$sample_id)) ||
      anyDuplicated(metadata$sample_id)
  ) {
    abort_user("Metadata sample_id values must be present and unique.")
  }
  if (anyNA(metadata$condition) || any(!nzchar(metadata$condition))) {
    abort_user("Metadata condition values must be present.")
  }

  count_samples <- colnames(counts)
  metadata_samples <- metadata$sample_id

  missing_metadata <- setdiff(count_samples, metadata_samples)
  missing_counts <- setdiff(metadata_samples, count_samples)
  if (length(missing_metadata) || length(missing_counts)) {
    details <- c(
      if (length(missing_metadata)) {
        sprintf("missing from metadata: %s", paste(missing_metadata, collapse = ", "))
      },
      if (length(missing_counts)) {
        sprintf("missing from counts: %s", paste(missing_counts, collapse = ", "))
      }
    )
    abort_user(sprintf("Sample IDs do not match exactly (%s).", paste(details, collapse = "; ")))
  }

  metadata <- metadata[match(count_samples, metadata$sample_id), , drop = FALSE]
  rownames(metadata) <- metadata$sample_id
  storage.mode(counts) <- "integer"

  structure(
    list(counts = counts, metadata = metadata, source = source_label),
    class = "rnaseq_bundle"
  )
}

load_uploaded_bundle <- function(count_input, metadata_input) {
  if (is.null(count_input) || is.null(metadata_input)) {
    abort_user("Choose both a count-matrix file and a metadata file.")
  }

  counts <- read_count_matrix(count_input$datapath, count_input$name)
  metadata <- read_metadata(metadata_input$datapath, metadata_input$name)
  validate_bundle(counts, metadata, source_label = "Uploaded data")
}
