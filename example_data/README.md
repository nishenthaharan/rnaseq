# Example data

The application contains a deterministic simulated dataset and requires no example files to start.

To create uploadable copies of the same experiment, run from the repository root:

```bash
Rscript scripts/generate_example_data.R
```

This creates:

- `counts_example.csv`: 1,200 genes × 8 samples of raw negative-binomial counts
- `metadata_example.csv`: four control and four treated samples, balanced across two batches

The first 80 simulated genes are upregulated and the next 70 are downregulated in the treated condition. The remaining genes have no programmed treatment effect. Randomness is fixed with seed `20260818` for reproducibility.
