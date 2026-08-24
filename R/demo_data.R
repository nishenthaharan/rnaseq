generate_demo_bundle <- function(seed = 20260818L, n_genes = 1200L) {
  if (length(seed) != 1L || !is.numeric(seed) || is.na(seed) || !is.finite(seed)) {
    abort_user("The simulated dataset seed must be one finite numeric value.")
  }
  if (length(n_genes) != 1L || is.na(n_genes) || n_genes < 10L) {
    abort_user("The simulated dataset requires at least 10 genes.")
  }
  n_genes <- as.integer(n_genes)

  had_random_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  previous_random_seed <- if (had_random_seed) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  on.exit({
    if (had_random_seed) {
      assign(".Random.seed", previous_random_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(seed)

  metadata <- data.frame(
    sample_id = sprintf("Sample_%02d", seq_len(8L)),
    condition = rep(c("Control", "Treated"), each = 4L),
    batch = rep(c("Batch_1", "Batch_2"), times = 4L),
    stringsAsFactors = FALSE
  )

  gene_id <- sprintf("GENE_%04d", seq_len(n_genes))
  baseline <- 2^stats::rnorm(n_genes, mean = 6.2, sd = 1.7)
  dispersion <- pmin(0.8, pmax(0.03, 0.12 + 25 / (baseline + 50)))
  library_factor <- c(0.82, 1.06, 0.94, 1.18, 0.88, 1.12, 0.97, 1.09)

  true_log2_fc <- numeric(n_genes)
  n_up <- min(80L, max(1L, floor(n_genes * 0.20)))
  n_down <- min(70L, max(1L, floor(n_genes * 0.18)))
  up_index <- seq_len(n_up)
  down_index <- seq.int(n_up + 1L, n_up + n_down)
  true_log2_fc[up_index] <- stats::rnorm(n_up, 1.8, 0.25)
  true_log2_fc[down_index] <- stats::rnorm(n_down, -1.6, 0.25)

  expected <- vapply(
    seq_len(nrow(metadata)),
    function(index) {
      treatment_effect <- if (metadata$condition[index] == "Treated") 2^true_log2_fc else 1
      baseline * treatment_effect * library_factor[index]
    },
    numeric(n_genes)
  )

  counts <- vapply(
    seq_len(ncol(expected)),
    function(index) {
      stats::rnbinom(n_genes, mu = expected[, index], size = 1 / dispersion)
    },
    numeric(n_genes)
  )
  storage.mode(counts) <- "integer"
  dimnames(counts) <- list(gene_id, metadata$sample_id)

  bundle <- validate_bundle(counts, metadata, source_label = "Built-in simulated data")
  bundle$truth <- data.frame(gene_id = gene_id, true_log2_fc = true_log2_fc)
  bundle
}
