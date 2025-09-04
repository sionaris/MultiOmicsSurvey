#!/usr/bin/env python3
"""
Repository Structure Listing Utility for MultiOmicsSurvey

This script provides a comprehensive listing of the repository structure
with descriptions of each component for the multi-omics analysis project.
"""

import os
import sys
from pathlib import Path

def describe_directory(dir_name):
    """Return description for known directories."""
    descriptions = {
        'Scripts': 'R scripts for data processing, analysis, and algorithm implementations',
        'Python': 'Python implementations of specific algorithms (MOFA, MONET, MSNE, PAMOGK)',
        'Resources': 'Data resources including TCGA, METABRIC datasets and algorithm inputs',
        'Results': 'Analysis results, survival evaluations, and algorithm outputs',
        'renv': 'R environment management files and package cache',
        'sessionInfo': 'R session information and reproducibility data',
        '.git': 'Git version control metadata',
        'Consensus': 'Consensus clustering analysis results',
        'MOVICS': 'MOVICS algorithm integration and comparison scripts',
        'Survival analysis': 'Survival analysis scripts and utilities',
        'automated_scripts': 'Automated analysis and processing scripts',
        'method_comparisons': 'Comparative analysis between different methods',
        'single_algorithm': 'Individual algorithm analysis and results',
        'MOFA': 'Multi-Omics Factor Analysis algorithm implementation',
        'MONET': 'Multi-Omics NETwork algorithm implementation',
        'MSNE': 'Multi-modal Stochastic Neighbor Embedding algorithm',
        'PAMOGK': 'Pathway-Assisted Multi-Omics Graph Kernel algorithm',
        'TCGA': 'The Cancer Genome Atlas data and preprocessing',
        'METABRIC': 'Molecular Taxonomy of Breast Cancer International Consortium data',
        'Pathways': 'Biological pathway annotations and gene sets',
        'HPC output': 'High-Performance Computing cluster analysis outputs',
        'algorithm_descriptions': 'Detailed descriptions and citations for algorithms',
        'transNEO': 'TransNEO dataset and related analyses'
    }
    return descriptions.get(dir_name, 'Additional project component')

def list_directory_tree(root_path, max_depth=2, current_depth=0):
    """List directory tree with descriptions."""
    if current_depth > max_depth:
        return
    
    items = []
    try:
        for item in sorted(os.listdir(root_path)):
            if item.startswith('.') and item not in ['.git', '.gitignore', '.Rprofile']:
                continue
            item_path = os.path.join(root_path, item)
            if os.path.isdir(item_path):
                items.append(('dir', item, item_path))
            else:
                items.append(('file', item, item_path))
    except PermissionError:
        return
    
    # Print directories first, then files
    for item_type, item_name, item_path in items:
        indent = "  " * current_depth
        if item_type == 'dir':
            description = describe_directory(item_name)
            print(f"{indent}📁 {item_name}/ - {description}")
            if current_depth < max_depth:
                list_directory_tree(item_path, max_depth, current_depth + 1)
        else:
            # Only show key files at root level
            if current_depth == 0:
                if item_name in ['README.md', 'MANIFEST.txt', '.gitignore', 
                               'MultiOmicsSurvey.Rproj', 'renv.lock']:
                    if item_name == 'README.md':
                        desc = "Project overview and documentation"
                    elif item_name == 'MANIFEST.txt':
                        desc = "Data file manifest and checksums"
                    elif item_name == '.gitignore':
                        desc = "Git ignore patterns for build artifacts and data"
                    elif item_name == 'MultiOmicsSurvey.Rproj':
                        desc = "RStudio project configuration"
                    elif item_name == 'renv.lock':
                        desc = "R package dependency lock file"
                    else:
                        desc = "Project configuration file"
                    print(f"{indent}📄 {item_name} - {desc}")

def main():
    """Main function to display repository structure."""
    repo_root = Path(__file__).parent
    
    print("=" * 60)
    print("MULTIOMICS SURVEY REPOSITORY STRUCTURE")
    print("=" * 60)
    print()
    print("🔬 A survey on computational methods for multi-omic and")
    print("   multimodal analysis in cancer research")
    print()
    print("Repository Contents:")
    print("-" * 20)
    
    list_directory_tree(str(repo_root))
    
    print()
    print("Key Features:")
    print("-" * 13)
    print("• Multi-language implementation (R + Python)")
    print("• 20+ multi-omics integration algorithms")
    print("• TCGA and METABRIC cancer datasets")
    print("• Comprehensive survival analysis")
    print("• Reproducible research environment (renv)")
    print("• Standardized evaluation metrics")
    print()
    print("Usage:")
    print("------")
    print("• R scripts: Use RStudio or R console in project directory")
    print("• Python algorithms: Each has individual conda/pip requirements")
    print("• Results: Pre-computed analyses available in Results/")
    print()
    print("For detailed information, see README.md and Scripts/ directory")
    print("=" * 60)

if __name__ == "__main__":
    main()