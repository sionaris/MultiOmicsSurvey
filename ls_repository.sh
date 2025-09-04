#!/bin/bash
# Repository Structure Listing for MultiOmicsSurvey
# Simple shell script to provide basic directory listing with context

echo "============================================================"
echo "MULTIOMICS SURVEY REPOSITORY - BASIC LISTING"
echo "============================================================"
echo
echo "📊 Repository Contents:"
echo "----------------------"

# Function to display directory with description
list_with_description() {
    local dir=$1
    local desc=$2
    if [ -d "$dir" ]; then
        echo "📁 $dir/ - $desc"
        echo "   $(ls -1 "$dir" | wc -l) items"
    fi
}

# Main directories with descriptions
list_with_description "Scripts" "R scripts for analysis and algorithms"
list_with_description "Python" "Python algorithm implementations"
list_with_description "Resources" "Data resources and inputs"
list_with_description "Results" "Analysis outputs and evaluations"
list_with_description "renv" "R environment management"
list_with_description "sessionInfo" "R session reproducibility data"

echo
echo "📄 Key Files:"
echo "-------------"
[ -f "README.md" ] && echo "📄 README.md - Project documentation"
[ -f "MANIFEST.txt" ] && echo "📄 MANIFEST.txt - Data file checksums ($(wc -l < MANIFEST.txt) entries)"
[ -f "renv.lock" ] && echo "📄 renv.lock - R package dependencies"
[ -f "MultiOmicsSurvey.Rproj" ] && echo "📄 MultiOmicsSurvey.Rproj - RStudio project"

echo
echo "🔍 Quick Stats:"
echo "--------------"
echo "• R scripts: $(find Scripts -name "*.R" 2>/dev/null | wc -l)"
echo "• Python modules: $(find Python -name "*.py" 2>/dev/null | wc -l)"
echo "• Result files: $(find Results -type f 2>/dev/null | wc -l)"
echo "• Total directories: $(find . -type d 2>/dev/null | wc -l)"

echo
echo "💡 Usage:"
echo "--------"
echo "• For detailed listing: python3 list_repository_structure.py"
echo "• For basic directory list: ls -la"
echo "• For tree view: tree (if available)"
echo "============================================================"