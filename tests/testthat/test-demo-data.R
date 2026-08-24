test_that("demo generation preserves an existing RNG state", {
  set.seed(104L)
  state_before <- .Random.seed

  generate_demo_bundle(seed = 2026L, n_genes = 40L)

  expect_identical(.Random.seed, state_before)
})

test_that("demo generation does not create a persistent RNG state", {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  previous_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (had_seed) {
      assign(".Random.seed", previous_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }
  generate_demo_bundle(seed = 2026L, n_genes = 40L)

  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
})
