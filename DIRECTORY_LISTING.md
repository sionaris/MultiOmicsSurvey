# Repository Directory Listing

This document provides information about the directory structure and listing utilities for the MultiOmicsSurvey repository.

## Quick Directory Listing

To get a quick overview of the repository structure, you can use one of the following commands:

### Basic Directory Listing
```bash
ls -la
```

### Enhanced Repository Listing (Shell Script)
```bash
./ls_repository.sh
```

### Detailed Structure Analysis (Python)
```bash
python3 list_repository_structure.py
```

## Repository Structure Overview

### Main Directories

- **Scripts/** - R scripts for data processing, analysis, and algorithm implementations
- **Python/** - Python implementations of specific algorithms (MOFA, MONET, MSNE, PAMOGK)
- **Resources/** - Data resources including TCGA, METABRIC datasets and algorithm inputs
- **Results/** - Analysis results, survival evaluations, and algorithm outputs
- **renv/** - R environment management files and package cache
- **sessionInfo/** - R session information and reproducibility data

### Key Files

- **README.md** - Project overview and documentation
- **MANIFEST.txt** - Data file manifest and checksums
- **MultiOmicsSurvey.Rproj** - RStudio project configuration
- **renv.lock** - R package dependency lock file
- **.gitignore** - Git ignore patterns for build artifacts and data

## Algorithm Implementations

The repository includes implementations of 20+ multi-omics integration algorithms:

### R-based Algorithms (in Scripts/)
- MOVICS integration framework
- Consensus clustering methods
- Survival analysis utilities

### Python-based Algorithms (in Python/)
- **MOFA** - Multi-Omics Factor Analysis
- **MONET** - Multi-Omics NETwork
- **MSNE** - Multi-modal Stochastic Neighbor Embedding  
- **PAMOGK** - Pathway-Assisted Multi-Omics Graph Kernel

## Data Resources

### Datasets (in Resources/)
- **TCGA** - The Cancer Genome Atlas data
- **METABRIC** - Molecular Taxonomy of Breast Cancer International Consortium
- **Pathways** - Biological pathway annotations and gene sets
- **transNEO** - TransNEO dataset for validation

## Usage Instructions

1. **For R analysis**: Open MultiOmicsSurvey.Rproj in RStudio
2. **For Python algorithms**: Each algorithm has individual conda/pip requirements
3. **For pre-computed results**: Check the Results/ directory
4. **For reproducibility**: Use renv to restore R package environment

## File Statistics

- Total R scripts: ~192
- Total Python modules: ~444  
- Total result files: ~6,000+
- Total directories: ~300+
- Data manifest entries: 667

This repository represents a comprehensive survey of computational methods for multi-omic and multimodal analysis in cancer research.