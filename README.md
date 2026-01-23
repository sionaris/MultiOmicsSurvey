# MultiOmicsSurvey

A comprehensive survey of computational methods for multi-omic integration in cancer research.

> **Unsupervised multi-omic integration algorithms** benchmarked across **5 data modalities** on TCGA breast cancer, validated on TCAG holdout set and the transNEO cohort.

---

## Quick Start

```r
# Clone and restore environment
git clone https://github.com/sionaris/MultiOmicsSurvey.git
renv::restore()
```
Study data can be found in the paper supplement and [Zenodo](https://zenodo.org/records/18339577)
<details>
<summary><b>⚠️ renv notes</b> — Packages not in lockfile</summary>

The following packages are **not included** in `renv.lock` and must be installed manually if needed:

| Package | Reason | Required for |
|---------|--------|--------------|
| `peakRAM` | HPC-only dependency | Benchmark memory profiling |
| `Rmosek` | Requires MOSEK license | wMKL, KLIC |
| `wMKL` | Custom Bioconductor install with `C++` code edit| wMKL algorithm |
| `devtools` | Development tool | Package installation from GitHub |

</details>

---

## Directory Structure

<details>
<summary><b>📁 Scripts</b> — Analysis code</summary>

1. `automated_scripts/` — Helper functions (sourced by other scripts)
2. `Download_TCGA_data.R` — Downloads and preprocesses TCGA-BRCA
3. `transNEO_pseudocount_determination.R` — RNA-seq normalisation for transNEO
4. `MOVICS/` — Baseline MOVICS analysis
5. `single_algorithm/` — Individual method runs (+HPC scripts)
6. `method_comparisons/` — Cross-method comparisons & benchmarks
7. `Survival analysis/` — Kaplan-Meier and Cox regression
8. `Consensus/` — Consensus clustering pipelines

</details>

<details>
<summary><b>📁 Python</b> — Python-based methods</summary>

- `MOFA/` — Multi-Omics Factor Analysis
- `MONET/` — Multi Omic clustering by Non-Exhaustive Types
- `MSNE/` — Multiple Similarity Network Embedding

</details>

<details>
<summary><b>📁 Resources</b> — Supporting data</summary>

- `HPC output/` — Results from cluster computing runs (individual runs)
- `Pathways/` — Gene sets for enrichment analysis
- `TCGA/` — Clinical data for TCGA
- `transNEO/` — Validation cohort ([paper](https://www.nature.com/articles/s41586-021-04278-5))
- `Performance/` — Benchmark input data

</details>

<details>
<summary><b>📁 Results</b> — Output files</summary>

- `single_algorithm/` — Individual method outputs
- `Comparisons/` — Cross-method comparison plots
- `Performance_benchmarks/` — Benchmark experiments (robustness, stability, scalability)
- `Consensus/` — Consensus clustering results
- `Survival_evaluations/` — Survival analysis outputs

</details>

<details>
<summary><b>📁 docs</b> — Documentation</summary>

- `R/function_documentation.pdf` — Core function reference
- `R/benchmark_function_documentation.pdf` — Benchmark analysis functions
- `Python/function_documentation.pdf` — Custom Python functions for MONET and MSNE

</details>

---

## Methods Included

| Category | Algorithms |
|----------|-----------|
| Similarity Network | ab-SNF, ANF, MDICC, MSNE, NEMO, RWR-F, RWR-NF, SNF, Spectrum |
| Multiple Kernel Learning | CIMLR, KLIC, wMKL |
| Matrix Factorisation | MOFA, LRAcluster, MFA |
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
