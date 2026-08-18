# RNASeq Explorer

RNASeq Explorer is a production-oriented R Shiny application for count-based bulk RNA-seq differential-expression analysis. It combines strict input validation, DESeq2 modelling, effect-size shrinkage, interactive result tables, publication-quality plots, and reproducible exports in a single GitHub-ready repository.

> **Statistical scope:** this application accepts unnormalised integer gene counts from bulk RNA-seq. Do not upload TPM, FPKM, CPM, ratios, microarray intensities, or single-cell matrices.

## What the application provides

- Exact validation of gene IDs, sample IDs, count values, conditions, and optional batches
- Independent low-count filtering with user-controlled thresholds
- DESeq2 size-factor normalisation, dispersion estimation, and negative-binomial modelling
- Optional additive batch adjustment (`~ batch + condition`) with model-rank validation
- DESeq2 normal-prior shrinkage of log2 fold changes, with a safe fallback to unshrunk estimates
- Benjamini–Hochberg false-discovery-rate control
- PCA, sample-correlation heatmap, raw library-size plot, and dispersion diagnostics
- Volcano plot, MA plot, ranked top-gene heatmap, and per-gene expression plot
- Searchable and sortable differential-expression result table
- Downloads for all results, significant genes, normalised counts, high-resolution plots, and the complete R analysis object
- Embedded analysis parameters, R version, platform, package versions, and analysis timestamp
- Built-in balanced simulated data so the entire workflow can be evaluated immediately

## Quick start

Install a recent version of R, then run:

```r
setwd("path/to/rnaseq")
source("install_dependencies.R")
shiny::runApp()
```

The application opens with the built-in dataset selected and automatically runs a `Treated vs Control` analysis. To use your own experiment, clear **Use the built-in simulated dataset**, upload the two files described below, select a contrast, and click **Run DESeq2 analysis**.

## Input files

### 1. Raw count matrix

CSV, TSV, or tab-delimited TXT. The first column contains unique gene identifiers; every subsequent column is one sample containing non-negative integer counts.

```text
gene_id,Sample_01,Sample_02,Sample_03,Sample_04
ENSG000001,213,198,542,488
ENSG000002,0,1,0,3
ENSG000003,72,65,81,93
```

### 2. Sample metadata

CSV, TSV, or tab-delimited TXT. Column names are standardised to lowercase snake case. `sample_id` and `condition` are required; `batch` is optional.

```text
sample_id,condition,batch
Sample_01,Control,Batch_1
Sample_02,Control,Batch_2
Sample_03,Treated,Batch_1
Sample_04,Treated,Batch_2
```

The sample names in the count-matrix header must match `sample_id` exactly. Each selected condition must have at least two biological replicates; three or more are strongly recommended.

Generate complete example CSV files with:

```bash
Rscript scripts/generate_example_data.R
```

The files will be written to `example_data/`.

## Statistical workflow

1. Validate raw counts and metadata, then reorder metadata to the count-matrix columns.
2. Retain genes meeting `count >= min_count` in at least `min_samples` samples.
3. Construct `~ condition` or `~ batch + condition` and reject rank-deficient designs.
4. Estimate DESeq2 size factors and gene-wise/trended/shrunken dispersions.
5. Fit negative-binomial generalised linear models and test the selected target-versus-reference contrast.
6. Shrink log2 fold changes with the DESeq2 normal prior when the design supports it.
7. Control the false discovery rate using Benjamini–Hochberg adjusted P-values.
8. Classify genes only when they pass both the FDR and absolute log2-fold-change thresholds.
9. Use the variance-stabilising transformation for unsupervised visualisation and clustering—not for differential testing.

All available conditions are used to estimate the model. The selected target and reference define the reported contrast. Batch correction is intentionally limited to an additive batch term; paired, longitudinal, interaction, nested, and surrogate-variable designs should be analysed with a study-specific design reviewed by a statistician.

## Result columns

| Column | Meaning |
|---|---|
| `gene_id` | Identifier copied from the first count-matrix column |
| `base_mean` | Mean DESeq2 size-factor-normalised count across samples |
| `log2_fold_change` | Shrunken target-versus-reference effect size when shrinkage succeeds |
| `lfc_se` | Standard error associated with the reported effect estimate |
| `statistic` | Wald test statistic from the unshrunk DESeq2 result |
| `p_value` | Raw Wald-test P-value |
| `adjusted_p_value` | Benjamini–Hochberg FDR |
| `classification` | Upregulated, downregulated, or not significant under the selected thresholds |

## Repository structure

```text
rnaseq-explorer-r/
├── app.R                         # Shiny UI, server, reactivity, and downloads
├── R/
│   ├── analysis.R                # DESeq2 model and result construction
│   ├── demo_data.R               # Deterministic simulated RNA-seq experiment
│   ├── io.R                      # File parsing and strict validation
│   ├── plots.R                   # Publication-quality visualisations
│   └── utils.R                   # Shared helpers and provenance
├── scripts/generate_example_data.R
├── tests/testthat/               # Unit and integration tests
├── www/styles.css                # Responsive application styling
├── Dockerfile
├── install_dependencies.R
└── .github/workflows/r-ci.yml
```

## Run with Docker

```bash
docker build -t rnaseq-explorer-r .
docker run --rm -p 3838:3838 rnaseq-explorer-r
```

Open <http://localhost:3838>.

## Run tests

```r
testthat::test_dir("tests/testthat", reporter = "summary")
```

## GitHub repository

This project is maintained at <https://github.com/nishenthaharan/rnaseq>. To clone it:

```bash
git clone https://github.com/nishenthaharan/rnaseq.git
cd rnaseq
```

Do not add access tokens or credentials to the repository. If GitHub requests authentication for a future push, use GitHub's browser sign-in, Git Credential Manager, or an appropriately scoped token stored outside the project.

## Interpretation and limitations

- Statistical significance is not equivalent to biological importance; interpret effect sizes, uncertainty, sample quality, and experimental context together.
- Poorly replicated, confounded, or heterogeneous designs cannot be repaired by visualisation software.
- For transcript-level quantification from Salmon/Kallisto, import counts and offsets using `tximport`/`tximeta` in a dedicated workflow rather than uploading TPM values.
- Gene-set enrichment, pathway analysis, annotation mapping, and isoform-level analysis are deliberately outside this focused first release.
- This research application is not a clinical diagnostic device and must not be used for patient-care decisions without appropriate validation and governance.

## Scientific references

- Love MI, Huber W, Anders S. *Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.* Genome Biology. 2014;15:550. <https://doi.org/10.1186/s13059-014-0550-8>
- Benjamini Y, Hochberg Y. *Controlling the false discovery rate: a practical and powerful approach to multiple testing.* Journal of the Royal Statistical Society B. 1995;57(1):289–300.

## Licence

Apache License 2.0. See [LICENSE](LICENSE).
