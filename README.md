# MultiOmicsSurvey

A comprehensive survey of computational methods for multi-omic integration in cancer research.

> **18 algorithms** benchmarked across **5 data modalities** on TCGA breast cancer, validated on transNEO cohort.

---

## Quick Start

```r
# Clone and restore environment
git clone https://github.com/sionaris/MultiOmicsSurvey.git
renv::restore()

# Download TCGA data
source("Scripts/Download_TCGA_data.R")
```

<details>
<summary><b>⚠️ renv notes</b> — Packages not in lockfile</summary>

The following packages are **not included** in `renv.lock` and must be installed manually if needed:

| Package | Reason | Required for |
|---------|--------|--------------|
| `peakRAM` | HPC-only dependency | Benchmark memory profiling |
| `Rmosek` | Requires MOSEK license | wMKL algorithm |
| `wMKL` | Custom Bioconductor install | wMKL algorithm |
| `devtools` | Development tool | Package installation from GitHub |

</details>

---

## Directory Structure

<details>
<summary><b>📁 Scripts</b> — Analysis code</summary>

1. `automated_scripts/` — Helper functions (sourced by other scripts)
2. `Download_TCGA_data.R` — Downloads and preprocesses TCGA-BRCA
3. `transNEO_pseudocount_determination.R` — RNA-seq normalization
4. `MOVICS/` — Baseline MOVICS analysis
5. `single_algorithm/` — Individual method runs (+HPC scripts)
6. `method_comparisons/` — Cross-method comparisons & benchmarks
7. `Survival analysis/` — Kaplan-Meier and Cox regression
8. `Consensus/` — Consensus clustering pipelines

</details>

<details>
<summary><b>📁 Python</b> — Python-based methods</summary>

- `MOFA/` — Multi-Omics Factor Analysis
- `MONET/` — Multi-Omic Network Embedding
- `MSNE/` — Multi-view SNE

</details>

<details>
<summary><b>📁 Resources</b> — Supporting data</summary>

- `HPC output/` — Results from cluster computing runs
- `Pathways/` — Gene sets for enrichment analysis
- `TCGA/` — Clinical and survival data
- `transNEO/` — Validation cohort ([paper](https://www.nature.com/articles/s41586-021-04278-5))
- `Performance/` — Benchmark input data

</details>

<details>
<summary><b>📁 Results</b> — Output files</summary>

- `single_algorithm/` — Individual method outputs
- `Comparisons/` — Cross-method comparison plots
- `Performance_benchmarks/` — Scaling experiments (feature/sample perturbations)
- `Consensus/` — Consensus clustering results
- `Survival_evaluations/` — Survival analysis outputs
- `master_dataset.csv` — Aggregated results table

</details>

<details>
<summary><b>📁 docs</b> — Documentation</summary>

- `R/function_documentation.pdf` — Core function reference
- `R/benchmark_function_documentation.pdf` — Benchmark analysis functions

</details>

---

## Methods Included

| Category | Algorithms |
|----------|-----------|
| Similarity Network | ab-SNF, ANF, MDICC, MSNE, NEMO, RWR-F, RWR-NF, SNF, Spectrum |
| Multiple Kernel Learning | CIMLR, KLIC, wMKL |
| Matrix Factorization | MOFA, LRAcluster, MFA |
| Graph-based | MONET |
| Bayesian | iClusterBayes |
| Consensus | COCA |

---

## Shiny app to explore our results interactively

[![GitHub](https://img.shields.io/badge/GitHub-sionaris%2FMO__survey__Shiny-181717?logo=github)](https://github.com/sionaris/MO_survey_Shiny)

---

## Citation

*Publication forthcoming*

---

## License

Apache 2.0
