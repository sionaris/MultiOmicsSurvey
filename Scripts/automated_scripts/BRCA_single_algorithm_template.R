# Libraries #####


# Import data from gitignored "Resources/BRCA complete/" folder #####
brca_cnv = read.csv("Resources/BRCA complete/BRCA_CNV.csv")
brca_met = read.csv("Resources/BRCA complete/BRCA_Methy.csv")
brca_exp = read.csv("Resources/BRCA complete/BRCA_mRNA.csv")
brca_mirna = read.csv("Resources/BRCA complete/BRCA_miRNA.csv")
brca = list(CNV = brca_cnv, Methylation = brca_met,
            RNAseq = brca_exp, miRNA = brca_mirna)
rm(brca_cnv, brca_exp, brca_met, brca_mirna); gc()

# Setup initial environment variables for markdown #####

# Preamble
home = getwd()
algorithm = ""
citation = ""
data_source = "" # e.g. TCGA, TCGA-transNEO, transNEO-PARTNER
data_types = "" # e.g. RNAseq, RNAseq-CNV-miRNA
evaluation_source = "" # e.g. PARTNER, transNEO-PARTNER 
title = paste0("Results from ", algorithm)
subtitle = paste0("<b>Train</b>: ", data_source, " ", data_types, " | <b>Evaluation</b>: ", evaluation_source)
in_a_nutshell = ""
optk_boolean = "" # either TRUE or FALSE. Answers whether the algorithm suggests an optimal k
optk_text = ifelse(optk_boolean == TRUE,
                   "suggests an estimate of the optimal number of multi-omic clusters $k$",
                   "does not suggest an optimal number of multi-omic clusters $k$")

# Detailed description of the algorithm
description = paste(readLines("Resources/algorithm_descriptions/.....Rmd"),
                    collapse = "\n") # File path to .Rmd file within Resources/algorithm_descriptions

# Run algorithm #####


# Export main results #####


# Evaluation #####


# Wrap up #####

# Put all parameters in a list
params = list(algorithm = algorithm, data_source = data_source, data_types = data_types,
              evaluation_source = evaluation_source, title = title, subtitle = subtitle,
              description = description, in_a_nutshell = in_a_nutshell, optk_text = optk_text)

# Render the R Markdown document with the parameters
rmarkdown::render(paste0(getwd(), "/Scripts/automated_scripts/single_algorithm_results_report.Rmd"), 
                  params = params, 
                  output_file = paste0(home, "/Results/single_algorithm/", 
                                       algorithm, "_report_", data_source, "_",
                                       data_types, "_eval_on_", evaluation_source,
                                       ".html"))
