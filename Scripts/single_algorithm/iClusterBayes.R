# Import data from gitignored "Resources/BRCA complete/" folder #####

# This RDS object was produced using the Scripts/MOVICS/MOVICS_baseline.R script
input = readRDS("Resources/TCGA/mm_input.rds")

# Setup environment variables for markdown #####

# Ensure reproducibility
RNGversion("4.2.2")
set.seed(123)

# Load custom helper functions
source("Scripts/automated_scripts/custom_functions.R")
source("Scripts/automated_scripts/modified_MOVICS_functions.R")

# Preamble
home = getwd()
algorithm = "iClusterBayes"
alg_feature_pref = "cols" # Where does the algorithm expect the features to be
citation = fetch_citation(algorithm = algorithm)
data_source = "TCGA" # e.g. TCGA, TCGA-transNEO, transNEO-PARTNER
data_types = "RNAseq-CNV-Methylation-miRNA-SNPs" # e.g. RNAseq, RNAseq-CNV-miRNA
evaluation_source = "transNEO" # e.g. PARTNER, transNEO-PARTNER 
title = paste0("Results from ", algorithm)
subtitle = paste0("<b>Train</b>: ", data_source, " ", data_types, 
                  " | <b>Evaluation</b>: ", evaluation_source)
in_a_nutshell = fetch_in_a_nutshell(algorithm = algorithm)
ground_truth_labels = openxlsx::read.xlsx("Results/MOVICS_baseline/MOVICS_TCGA_RNAseq-CNV-Methylation-miRNA-SNPs_eval_on_transNEO_clusterings.xlsx")
ground_truth_k = 2 # optk from MOVICS
optk_boolean = "TRUE" # either TRUE or FALSE. Answers whether the algorithm suggests an optimal k
optk_text = ifelse(optk_boolean == TRUE,
                   "<u>suggests</u> an estimate of the optimal number of multi-omic clusters $k$",
                   "<u>does not suggest</u> an optimal number of multi-omic clusters $k$")

# Detailed description of the algorithm
description = paste(readLines(paste0("Resources/algorithm_descriptions/", algorithm,
                                     "_description.Rmd")),
                    collapse = "\n") # File path to .Rmd file within Resources/algorithm_descriptions

# Create algorithm directory if it doesn't exist
if (!dir.exists(paste0(home, "/Results/single_algorithm/", 
                       algorithm))) {
  dir.create(paste0(home, "/Results/single_algorithm/", 
                    algorithm))
}

# Preprocessing flags and code #####
library(stringr)
library(dplyr)

# All preprocessing for this input has already been performed using the 
# Scripts/MOVICS/MOVICS_baseline.R script

# However iClusterBayes prefers features in columns so we transpose the matrices.

# Extract the names of the modalities that will be used
modalities = unlist(strsplit(data_types, "-"))

# Replace with TRUE wherever features are in rows
features_in_rows = rep(TRUE, length(modalities))

# Run algorithm #####
library(iClusterPlus)

# Check what kind of arrangement the algorithm requires as input 
# (i.e. features in rows or columns)?
if (alg_feature_pref == "rows") {
  rogue_indices = which(features_in_rows == FALSE)
  if (length(rogue_indices >= 1)) {
    for (index in rogue_indices) {
      cols = colnames(input[[index]])
      rows = rownames(input[[index]])
      input[[index]] = t(input[[index]])
      rownames(input[[index]]) = cols # transpose names
      colnames(input[[index]]) = rows # transpose names
      rm(rows, cols)
    }
  }
} else if (alg_feature_pref == "cols") {
  rogue_indices = which(features_in_rows == TRUE)
  if (length(rogue_indices >= 1)) {
    for (index in rogue_indices) {
      cols = colnames(input[[index]])
      rows = rownames(input[[index]])
      input[[index]] = t(input[[index]])
      rownames(input[[index]]) = cols # transpose names
      colnames(input[[index]]) = rows # transpose names
      rm(rows, cols)
    }
  }
}
rm(rogue_indices, index); gc()

# Import clinical data for the TCGA samples of interest
clinical_data = openxlsx::read.xlsx("Resources/TCGA/clinical_data.xlsx")

# Setup ###
saveRDS(input, "Resources/iCB_input.rds")

# Run the tune.iClusterBayes() function at an HPC for a gird of parameters
# See /Scripts/automated_scripts/iCB_scripts_preparation.R
# and /Scripts/single_algorithm/iCB_HPC/*.R and *.sh

# After obtaining the results from the HPC jobs, we run the results inspection
# script: Scripts/iCB_inspect_HPC_tune_results.R

# Turns out that the optimal combination is sdev = 0.015, beta.var.scale = 0.5

# Import the "optimal" results
tune_results = readRDS("Resources/HPC output/iCB_HPC/final_tune_results_sdev_0.015_beta_0.5.rds")

# Range of possible clusters
K_values = 1:9

# Extract BIC and deviance ratio for each K
BIC_values <- sapply(tune_results$fit, function(fit) fit$BIC)
dev_ratio_values <- sapply(tune_results$fit, function(fit) fit$dev.ratio)

# Plot BIC and deviance ratio
plot(K_values, BIC_values, type = "b", xlab = "Number of Clusters (K)", ylab = "BIC", main = "BIC vs K")
plot(K_values, dev_ratio_values, type = "b", xlab = "Number of Clusters (K)", ylab = "Deviance Ratio", main = "Deviance Ratio vs K")

# It's hard to pick a k if we consider both plots, although k = 4+1 = 5 looks reasonable
bw.col = gplots::colorpanel(2,low="white",high="black")
col.scheme = alist()
col.scheme[[1]] = bw.col
col.scheme[[2]] = gplots::bluered(256)
col.scheme[[3]] = gplots::bluered(256)
col.scheme[[4]] = gplots::bluered(256)
col.scheme[[5]] = gplots::bluered(256)

library(lattice)
plotHeatmap(tune_results$fit[[which.min(lapply(tune_results$fit, function(x) x$BIC))]], 
            datasets = input,
            type = c("binomial", rep("gaussian", 4)), col.scheme = col.scheme,
            row.order = c(T, T, T, T, T))

# Setting the optimal fit
optimal.fit = tune_results$fit[[4]] # k = 5

# SNPs
plot(optimal.fit$beta.pp[[1]],xlab="Genes",ylab="Posterior probability",
     main="SNPs")

# RNAseq
plot(optimal.fit$beta.pp[[2]],xlab="Genes",ylab="Posterior probability",
     main="Normalized RNA-seq")

# CNV
plot(optimal.fit$beta.pp[[3]],xlab="Genes",ylab="Posterior probability",
     main="Normalized Copy Number")

# miRNA
plot(optimal.fit$beta.pp[[4]],xlab="miRNAs",ylab="Posterior probability",
     main="Normalized miRNA-seq")

# Methylation
plot(optimal.fit$beta.pp[[3]],xlab="Genes",ylab="Posterior probability",
     main="Normalized M values")

# Extract clusters
iClusterBayes_clusters = as.data.frame(list(Cluster = optimal.fit$clusters)) %>%
  dplyr::mutate(Sample.ID = rownames(input$SNPs))

# Feature selection #####
features = alist()
for (i in 1:length(input)) {
  features[[i]] = colnames(input[[i]])
}

sigfeatures=alist()
for(i in 1:length(input)){
  rowsum=apply(abs(optimal.fit$beta[[i]]),1, sum)
  upper=quantile(rowsum,prob=0.75)
  sigfeatures[[i]]=(features[[i]])[which(rowsum>upper)]
}

names(sigfeatures) = names(input)

# Main results #####
# Examine cluster similarity to MOVICS by measuring NMI and ARI indices #####
# (Jaccard may be misleading)

# Calculate ARI and NMI
library(mclust)
library(clue)

ARI_to_MOVICS = calculate_ari_index(cluster_df1 = ground_truth_labels,
                                    cluster_df2 = iClusterBayes_clusters,
                                    sample_col = "Sample.ID",
                                    clust_col = "Cluster",
                                    suffixes = c("_MOVICS_CS",
                                                 paste0("_", algorithm)))

NMI_to_MOVICS = calculate_nmi_index(cluster_df1 = ground_truth_labels,
                                    cluster_df2 = iClusterBayes_clusters,
                                    sample_col = "Sample.ID",
                                    clust_col = "Cluster",
                                    suffixes = c("_MOVICS_CS",
                                                 paste0("_", algorithm)))

# Good agreement compared to the MOVICS. Results differ

# MOVICS-like analysis #####
library(MOVICS)
library(ComplexHeatmap)

# Import coloring scheme
scheme = readRDS("Resources/scheme.rds")
annCol = scheme$annCol
annColors = scheme$annColors
cluster_colors = c("#2EC4B6", "#E71D36", 
                   "#FF9F1C", "#BDD5EA", "#FFA5AB")
col.list = scheme$col.list
var2comp = scheme$var2comp %>%
  dplyr::select(-`Consensus Subtype`) %>%
  mutate(Sample.ID = rownames(.)) %>%
  inner_join(iClusterBayes_clusters, by = "Sample.ID") %>%
  tibble::column_to_rownames(var = "Sample.ID") %>%
  mutate(iClusterBayes = paste0(algorithm, Cluster)) %>%
  dplyr::select(iClusterBayes, everything()) %>%
  dplyr::select(-Cluster)
rm(scheme); gc()

# Silhouette

# transformed_aff = transform_affinity_matrix(final_affinity_matrix,
#                                             threshold = 100, norm_quant = 0,
#                                             norm_method = "divide by quantile")

z_spearman_matrix = cor(t(optimal.fit$meanZ), method = "spearman")
dimnames(z_spearman_matrix) = list(rownames(input$SNPs), rownames(input$SNPs))
sil = compute_silhouette(cluster_df = iClusterBayes_clusters %>% dplyr::rename(samID = Sample.ID),
                         similarity_matrix = z_spearman_matrix,
                         normalize_matrix = TRUE)

getSilhouette_ggplot(sil      = sil,
                     fig.path = paste0(home, "/Results/single_algorithm/", algorithm),
                     fig.name = "Silhouette",
                     height   = 5.5,
                     width    = 5.5,
                     axis_label_size = 12,
                     axis_label_font = "bold",
                     text_size = 1.5,
                     title_size = 16,
                     algorithm = algorithm,
                     save_plot = TRUE)
dev.off()

# Heatmap prep
plotdata <- lapply(lapply(input, as.matrix), 
                   function(mat) mat[, colSums(mat != 0) > 0])
plotdata <- lapply(plotdata, t)
plotdata = getStdiz(
  data = plotdata,
  halfwidth = c(NA, 3, 3, 3, 3), # No halfwidth for SNPs
  centerFlag = c(F, F, F, F, F),
  scaleFlag = c(F, F, F, F, F)
)

plot_object = list(clust.res = iClusterBayes_clusters %>%
                     dplyr::rename(samID = Sample.ID, clust = Cluster))

# Export consensus clustering object
clust = as.data.frame(plot_object$clust.res)
colnames(clust) = c("Sample.ID", "Cluster")
clust$Cluster = paste0(algorithm, clust$Cluster)
openxlsx::write.xlsx(clust, paste0(home, "/Results/single_algorithm/", algorithm, "/", 
                                   algorithm, "_", data_source, "_",
                                   data_types, "_eval_on_", evaluation_source,
                                   "_clusterings.xlsx"))

# comprehensive heatmap (may take a while)
getMoHeatmap_single_algorithm(algorithm_name = algorithm,
                              data          = plotdata,
                              row.title     = names(heatmap_plotdata),
                              is.binary     = c(T,F,F,F,F), 
                              legend.name   = c("SNPs",
                                                "Standardized RNAseq norm. counts",
                                                "Standardized CNV",
                                                "Standardized miRNA norm. counts",
                                                "Standardized Methylation M-values"
                              ),
                              clust.res     = plot_object$clust.res, # consensusMOIC-like results
                              clust.dend    = NULL, # show no dendrogram for samples
                              show.rownames = c(F,F,F,F,F), # specify for each omics data
                              show.colnames = FALSE, # show no sample names
                              show.row.dend = c(F,F,F,F,F), # show no dendrogram for features
                              annRow        = NULL, # no selected features
                              color         = col.list,
                              annCol        = annCol, # annotation for samples
                              annColors     = annColors, # annotation color
                              width         = 20, # width of each subheatmap
                              height        = 10, # height of each subheatmap
                              fig.path      = paste0(home, "/Results/single_algorithm/", algorithm),
                              fig.name      = paste0("default_", algorithm, "_Comprehensive_heatmap"))
dev.off()
gc()

# Clinical variables ###
# Remove unknown levels for statistical tests
var2comp_nonas = var2comp
for (i in 1:ncol(var2comp)) {
  nas = which(var2comp[, i] == "Unknown")
  var2comp_nonas[nas, i] = NA
  empties = which(var2comp[, i] == "")
  var2comp_nonas[empties, i] = NA
}
rm(nas, empties); gc()

# Statistical comparisons
clin_comp = compClinvar_single_algorithm(algorithm_name = algorithm,
                                         moic.res = plot_object,
                                         var2comp = var2comp_nonas,
                                         strata = algorithm,
                                         factorVars = c("vital_status", "race_list", "ethnicity",
                                                        "history_of_neoadjuvant_treatment",
                                                        "primary_lymph_node_presentation_assessment",
                                                        "histological_type", "menopause_status",
                                                        "breast_carcinoma_progesterone_receptor_status",
                                                        "breast_carcinoma_estrogen_receptor_status",
                                                        "lab_proc_her2_neu_immunohistochemistry_receptor_status",
                                                        "distant_metastasis_present_ind2",
                                                        "stage_event_pathologic_stage"),
                                         nonnormalVars = c("days_to_birth", "days_to_death",
                                                           "days_to_last_known_alive", 
                                                           "days_to_last_followup",
                                                           "age_at_initial_pathologic_diagnosis",
                                                           "er_level_cell_percentage_category",
                                                           "progesterone_receptor_level_cell_percent_category",
                                                           "number_of_lymphnodes_positive_by_ihc",
                                                           "number_of_lymphnodes_positive_by_he"),
                                         includeNA = FALSE,
                                         doWord = TRUE,
                                         tab.name = "Summary_of_clinical_variables",
                                         res.path = paste0(home, "/Results/single_algorithm/", algorithm, "/"),
                                         output_pdf = TRUE,
                                         pdf_level_col_width = c("7em", "10em"),
                                         pdf_count_col_width = "10em",
                                         pdf_pval_col_width = "3em",
                                         pdf_test_col_width = "8em",
                                         pdf_tab_font_size = 9)

clin_ordinal_comp = compClinvar_ordinal_single_algorithm(algorithm_name = algorithm,
                                                         moic.res = plot_object,
                                                         var2comp = var2comp_nonas %>%
                                                           dplyr::select(number_of_lymphnodes_positive_by_ihc,
                                                                         number_of_lymphnodes_positive_by_he,
                                                                         iClusterBayes),
                                                         strata = algorithm,
                                                         ordinalVars = c("number_of_lymphnodes_positive_by_ihc",
                                                                         "number_of_lymphnodes_positive_by_he"),
                                                         includeNA = FALSE,
                                                         tab.name = "Summary of ordinal clinical variables",
                                                         res.path = paste0(home, "/Results/single_algorithm/", algorithm, "/"),
                                                         output_pdf = TRUE,
                                                         pdf_template_loc = paste0(home, "/Scripts/automated_scripts/clincomp_template.Rmd"),
                                                         pdf_level_col_width = c("7em", "10em"),
                                                         pdf_count_col_width = "10em",
                                                         pdf_pval_col_width = "3em",
                                                         pdf_test_col_width = "8em",
                                                         pdf_tab_font_size = 9)

# Oncoprint ###
oncoprint <- compMut_single_algorithm(algorithm_name = algorithm,
                                      moic.res  = plot_object,
                                      mut.matrix   = plotdata$SNPs, # binary somatic mutation matrix
                                      doWord       = TRUE, # generate table in .docx format
                                      doPlot       = TRUE, # draw OncoPrint
                                      freq.cutoff  = 0.05, # keep those genes that mutated in at least 5% of samples
                                      p.adj.cutoff = 0.05, # keep those genes with adjusted p value < 0.05 to draw OncoPrint
                                      innerclust   = TRUE, # perform clustering within each subtype
                                      annCol       = annCol, # same annotation for heatmap
                                      annColors    = annColors, # same annotation color for heatmap
                                      width        = 12, 
                                      height       = 6,
                                      fig.name     = paste0(algorithm, "_", data_source, "_",
                                                            data_types, "_eval_on_", evaluation_source,
                                                            "_oncoprint"),
                                      tab.name     = "Independent test between subtype and mutation",
                                      fig.path     = paste0(home, "/Results/single_algorithm/", algorithm),
                                      res.path     = paste0(home, "/Results/single_algorithm/", algorithm))

# Similar to MOVICS: TP53 and PIK3CA patterns

# Drug sensitivity comparison ###
drug_sensitivity <- compDrugsen_single_algorithm(algorithm_name = algorithm,
                                                 moic.res    = plot_object,
                                                 norm.expr   = plotdata$RNAseq,
                                                 drugs       = c("Cisplatin", "Paclitaxel", "Lapatinib",
                                                                 "Doxorubicin", "5-Fluorouracil",
                                                                 "Sorafenib"), # a vector of names of drug in GDSC
                                                 tissueType  = "breast", # choose specific tissue type to construct ridge regression model
                                                 test.method = "nonparametric", # statistical testing method
                                                 prefix      = "Violin_plot_of_IC50",
                                                 seed = 123,
                                                 fig.path = paste0(home, "/Results/single_algorithm/", algorithm))

# Agreement with other subtypes ###
subtype_agreement <- compAgree_single_algorithm(algorithm_name = algorithm,
                                                moic.res  = plot_object,
                                                subt2comp = annCol[, c("ER status", "PR status",
                                                                       "HER2 status", "Metastasis", "Stage")],
                                                doPlot    = TRUE,
                                                box.width = 0.2,
                                                fig.name  = "Classification_agreement",
                                                fig.path  = paste0(home, "/Results/single_algorithm/", algorithm),
                                                width     = 12)
dev.off()

# DGEA ###
dgea = runDEA_mod(dea.method = "limma", # we use normalized data as input
                  expr = plotdata$RNAseq,
                  moic.res = plot_object,
                  prefix = "dgea_",
                  sort.p = TRUE,
                  overwt = TRUE,
                  verbose = TRUE,
                  res.path = paste0(home, "/Results/single_algorithm/", algorithm),
                  algorithm = algorithm)

# # Identify unique subtype biomarkers
# # 1. Up-regulated markers
dgea.marker.up <- runMarker_single_algorithm(algorithm_name = algorithm,
                                             moic.res = plot_object,
                                             dea.method    = "limma", # name of DEA method
                                             prefix        = "dgea_", # MUST be the same of argument in runDEA()
                                             dat.path      = paste0(home, "/Results/single_algorithm/", algorithm), # path of DEA files
                                             res.path      = paste0(home, "/Results/single_algorithm/", algorithm), # path to save marker files
                                             p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                                             p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                                             dirct         = "up", # direction of dysregulation in expression
                                             n.marker      = 100, # number of biomarkers for each subtype
                                             doplot        = TRUE, # generate diagonal heatmap
                                             norm.expr     = plotdata$RNAseq, # use normalized expression as heatmap input
                                             annCol        = annCol, # sample annotation in heatmap
                                             annColors     = annColors, # colors for sample annotation
                                             show_rownames = TRUE, # show no rownames (biomarker name)
                                             centerFlag = F,
                                             scaleFlag = F,
                                             halfwidth = 3,
                                             fig.name      = "upregulated_biomarkers_heatmap",
                                             fig.path = paste0(home, "/Results/single_algorithm/", algorithm),
                                             width = 14,
                                             height = 12,
                                             fontsize_row = 3,
                                             name = "normalized RNA-seq")
dev.off()

# # 2. Down-regulated markers
dgea.marker.down <- runMarker_single_algorithm(algorithm_name = algorithm,
                                               moic.res = plot_object,
                                               dea.method    = "limma", # name of DEA method
                                               prefix        = "dgea_", # MUST be the same of argument in runDEA()
                                               dat.path      = paste0(home, "/Results/single_algorithm/", algorithm), # path of DEA files
                                               res.path      = paste0(home, "/Results/single_algorithm/", algorithm), # path to save marker files
                                               p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                                               p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                                               dirct         = "down", # direction of dysregulation in expression
                                               n.marker      = 100, # number of biomarkers for each subtype
                                               doplot        = TRUE, # generate diagonal heatmap
                                               norm.expr     = plotdata$RNAseq, # use normalized expression as heatmap input
                                               annCol        = annCol, # sample annotation in heatmap
                                               annColors     = annColors, # colors for sample annotation
                                               show_rownames = TRUE, # show no rownames (biomarker name)
                                               centerFlag = F,
                                               scaleFlag = F,
                                               halfwidth = 3,
                                               fig.name      = "downregulated_biomarkers_heatmap",
                                               fig.path = paste0(home, "/Results/single_algorithm/", algorithm),
                                               width = 14,
                                               height = 12,
                                               fontsize_row = 3,
                                               name = "normalized RNA-seq")
dev.off()

# DMEA ###
dmea = runDEA_mod(dea.method = "limma", # we use normalized data as input
                  expr = plotdata$Methylation,
                  moic.res = plot_object,
                  prefix = "dmea_",
                  sort.p = TRUE,
                  overwt = TRUE,
                  verbose = TRUE,
                  res.path = paste0(home, "/Results/single_algorithm/", algorithm),
                  algorithm = algorithm)

# # Identify unique subtype biomarkers
# # 1. Up-regulated markers
methyl.marker.up <- runMarker_single_algorithm(algorithm_name = algorithm,
                                               moic.res = plot_object,
                                               dea.method    = "limma", # name of DEA method
                                               prefix        = "dmea_", # MUST be the same of argument in runDEA()
                                               dat.path      = paste0(home, "/Results/single_algorithm/", algorithm), # path of DEA files
                                               res.path      = paste0(home, "/Results/single_algorithm/", algorithm), # path to save marker files
                                               p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                                               p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                                               dirct         = "up", # direction of dysregulation in expression
                                               n.marker      = 100, # number of biomarkers for each subtype
                                               doplot        = TRUE, # generate diagonal heatmap
                                               norm.expr     = plotdata$Methylation, # use normalized expression as heatmap input
                                               annCol        = annCol, # sample annotation in heatmap
                                               annColors     = annColors, # colors for sample annotation
                                               show_rownames = TRUE, # show no rownames (biomarker name)
                                               centerFlag = F,
                                               scaleFlag = F,
                                               halfwidth = 3,
                                               fig.name      = "hypermethylated_biomarkers_heatmap",
                                               fig.path = paste0(home, "/Results/single_algorithm/", algorithm),
                                               width = 14,
                                               height = 12,
                                               fontsize_row = 0, # 3 default
                                               name = "normalized Methylation")
dev.off()

# # 2. Down-regulated markers
methyl.marker.down <- runMarker_single_algorithm(algorithm_name = algorithm,
                                                 moic.res = plot_object,
                                                 dea.method    = "limma", # name of DEA method
                                                 prefix        = "dmea_", # MUST be the same of argument in runDEA()
                                                 dat.path      = paste0(home, "/Results/single_algorithm/", algorithm), # path of DEA files
                                                 res.path      = paste0(home, "/Results/single_algorithm/", algorithm), # path to save marker files
                                                 p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                                                 p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                                                 dirct         = "down", # direction of dysregulation in expression
                                                 n.marker      = 100, # number of biomarkers for each subtype
                                                 doplot        = TRUE, # generate diagonal heatmap
                                                 norm.expr     = plotdata$Methylation, # use normalized expression as heatmap input
                                                 annCol        = annCol, # sample annotation in heatmap
                                                 annColors     = annColors, # colors for sample annotation
                                                 show_rownames = TRUE, # show no rownames (biomarker name)
                                                 centerFlag = F,
                                                 scaleFlag = F,
                                                 halfwidth = 3,
                                                 fig.name      = "hypomethylated_biomarkers_heatmap",
                                                 fig.path = paste0(home, "/Results/single_algorithm/", algorithm),
                                                 width = 14,
                                                 height = 12,
                                                 fontsize_row = 0, # 3 default
                                                 name = "normalized Methylation")
dev.off()

# DmiREA ###
dmiRea = runDEA_mod(dea.method = "limma", # we use normalized data as input
                    expr = plotdata$miRNA,
                    moic.res = plot_object,
                    prefix = "dmiRea_",
                    sort.p = TRUE,
                    overwt = TRUE,
                    verbose = TRUE,
                    res.path = paste0(home, "/Results/single_algorithm/", algorithm),
                    algorithm = algorithm)

# # Identify unique subtype biomarkers
# # 1. Up-regulated markers
miRNA.marker.up <- runMarker_single_algorithm(algorithm_name = algorithm,
                                              moic.res = plot_object,
                                              dea.method    = "limma", # name of DEA method
                                              prefix        = "dmiRea_", # MUST be the same of argument in runDEA()
                                              dat.path      = paste0(home, "/Results/single_algorithm/", algorithm), # path of DEA files
                                              res.path      = paste0(home, "/Results/single_algorithm/", algorithm), # path to save marker files
                                              p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                                              p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                                              dirct         = "up", # direction of dysregulation in expression
                                              n.marker      = 100, # number of biomarkers for each subtype
                                              doplot        = TRUE, # generate diagonal heatmap
                                              norm.expr     = plotdata$miRNA, # use normalized expression as heatmap input
                                              annCol        = annCol, # sample annotation in heatmap
                                              annColors     = annColors, # colors for sample annotation
                                              show_rownames = TRUE, # show no rownames (biomarker name)
                                              centerFlag = F,
                                              scaleFlag = F,
                                              halfwidth = 3,
                                              fig.name      = "upregulated_miRNA_biomarkers_heatmap",
                                              fig.path = paste0(home, "/Results/single_algorithm/", algorithm),
                                              width = 14,
                                              height = 12,
                                              fontsize_row = 0, # 3 default
                                              name = "normalized miRNA")
dev.off()

# # 2. Down-regulated markers
miRNA.marker.down <- runMarker_single_algorithm(algorithm_name = algorithm,
                                                moic.res = plot_object,
                                                dea.method    = "limma", # name of DEA method
                                                prefix        = "dmiRea_", # MUST be the same of argument in runDEA()
                                                dat.path      = paste0(home, "/Results/single_algorithm/", algorithm), # path of DEA files
                                                res.path      = paste0(home, "/Results/single_algorithm/", algorithm), # path to save marker files
                                                p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                                                p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                                                dirct         = "down", # direction of dysregulation in expression
                                                n.marker      = 100, # number of biomarkers for each subtype
                                                doplot        = TRUE, # generate diagonal heatmap
                                                norm.expr     = plotdata$miRNA, # use normalized expression as heatmap input
                                                annCol        = annCol, # sample annotation in heatmap
                                                annColors     = annColors, # colors for sample annotation
                                                show_rownames = TRUE, # show no rownames (biomarker name)
                                                centerFlag = F,
                                                scaleFlag = F,
                                                halfwidth = 3,
                                                fig.name      = "downregulated_miRNA_biomarkers_heatmap",
                                                fig.path = paste0(home, "/Results/single_algorithm/", algorithm),
                                                width = 14,
                                                height = 12,
                                                fontsize_row = 0, # 3 default
                                                name = "normalized miRNA")
dev.off()

# GSEA ###
# Load MSigDb file
MSIGDB.FILE <- paste0(home, "/Resources/Pathways/GO-BP_c5.go.bp.v2024.1.Hs.symbols.gmt")

# GSEA up-regulated
RNGversion("4.2.2")
set.seed(123)
gsea.up <- runGSEA_mod_4.4_single_algorithm(algorithm_name = algorithm,
                                            moic.res     = plot_object,
                                            dea.method   = "limma", # name of DEA method
                                            prefix       = "dgea_", # MUST be the same of argument in runDEA()
                                            dat.path      = paste0(home, "/Results/single_algorithm/", algorithm), # path of DEA files
                                            res.path      = paste0(home, "/Results/single_algorithm/", algorithm), # path to save marker files
                                            msigdb.path  = MSIGDB.FILE, # MUST be the ABSOLUTE path of msigdb file
                                            norm.expr    = plotdata$RNAseq, # use normalized expression to calculate enrichment score
                                            dirct        = "up", # direction of dysregulation in pathway
                                            n.path       = 20,
                                            p.cutoff     = 0.05, # p cutoff to identify significant pathways
                                            p.adj.cutoff = 0.05, # padj cutoff to identify significant pathways
                                            gsva.method  = "gsva", # method to calculate single sample enrichment score
                                            name         = "GSVA scores", # name for colorbar
                                            norm.method  = "mean", # normalization method to calculate subtype-specific enrichment score
                                            fig.name     = "upregulated_pathway_heatmap",
                                            nPerm = 10000,
                                            minGSSize = 5,
                                            maxGSSize = 500,
                                            fig.path = paste0(home, "/Results/single_algorithm/", algorithm),
                                            width = 14, height = 12)

# GSEA down-regulated
RNGversion("4.2.2")
set.seed(123)
gsea.down <- runGSEA_mod_4.4_single_algorithm(algorithm_name = algorithm,
                                              moic.res     = plot_object,
                                              dea.method   = "limma", # name of DEA method
                                              prefix       = "dgea_", # MUST be the same of argument in runDEA()
                                              dat.path      = paste0(home, "/Results/single_algorithm/", algorithm), # path of DEA files
                                              res.path      = paste0(home, "/Results/single_algorithm/", algorithm), # path to save marker files
                                              msigdb.path  = MSIGDB.FILE, # MUST be the ABSOLUTE path of msigdb file
                                              norm.expr    = plotdata$RNAseq, # use normalized expression to calculate enrichment score
                                              dirct        = "down", # direction of dysregulation in pathway
                                              n.path       = 20,
                                              p.cutoff     = 0.05, # p cutoff to identify significant pathways
                                              p.adj.cutoff = 0.05, # padj cutoff to identify significant pathways
                                              gsva.method  = "gsva", # method to calculate single sample enrichment score
                                              name         = "GSVA scores", # name for colorbar
                                              norm.method  = "mean", # normalization method to calculate subtype-specific enrichment score
                                              fig.name     = "downregulated_pathway_heatmap",
                                              nPerm = 10000,
                                              minGSSize = 5,
                                              maxGSSize = 500,
                                              fig.path = paste0(home, "/Results/single_algorithm/", algorithm),
                                              width = 14, height = 12)

# Gene set variation analysis #####
# locate ABSOLUTE path of gene set file
GSET.FILE <- paste0(home, "/Resources/Pathways/gene_sets_of_interest.gmt")

RNGversion("4.2.2")
set.seed(123)
gsva.res = runGSVA_mod_4.4_single_algorithm(algorithm_name = algorithm,
                                            moic.res      = plot_object,
                                            norm.expr     = plotdata$RNAseq,
                                            gset.gmt.path = GSET.FILE, # ABSOLUTE path of gene set file
                                            gsva.method   = "gsva", # method to calculate single sample enrichment score
                                            annCol        = annCol,
                                            annColors     = annColors,
                                            fig.path      = paste0(home, "/Results/single_algorithm/", algorithm),
                                            fig.name      = "gene_sets_of_interest_heatmap",
                                            centerFlag    = F,
                                            scaleFlag     = F,
                                            distance      = 'euclidean',
                                            linkage       = 'average',
                                            show_rownames = TRUE,
                                            show_colnames = FALSE,
                                            height        = 8,
                                            width         = 12,
                                            name          = "GSVA scores")
dev.off()

# Hierarchical clustering of pathways
library(pathfindR)
library(fastcluster)

# Get unique pathways for each subtype
GSEAfiles_up <- sort(dir(paste0(home, "/Results/single_algorithm/", algorithm), 
                         pattern = "unique_upexpr_pathway.txt$"))
GSEAfiles_down <- sort(dir(paste0(home, "/Results/single_algorithm/", algorithm), 
                           pattern = "unique_downexpr_pathway.txt$"))

unique_upexpr_pathways = list()
for (i in 1:length(gsea.up$gsea.list)) {
  unique_upexpr_pathways[[i]] = data.table::fread(paste0(paste0(home, "/Results/single_algorithm/", algorithm), 
                                                         "/", GSEAfiles_up[i]),
                                                  header = TRUE, sep = "\t")
}

unique_downexpr_pathways = list()
for (i in 1:length(gsea.down$gsea.list)) {
  unique_downexpr_pathways[[i]] = data.table::fread(paste0(paste0(home, "/Results/single_algorithm/", algorithm), 
                                                           "/", GSEAfiles_down[i]),
                                                    header = TRUE, sep = "\t")
}

names(unique_downexpr_pathways) = names(unique_upexpr_pathways) = names(gsea.up$gsea.list)

# Filter GSEA input
gsea.up_unique = gsea.up
for (i in 1:length(unique_upexpr_pathways)) {
  unq = unique_upexpr_pathways[[i]]$V1
  gsea.up_unique$gsea.list[[i]]@result = gsea.up_unique$gsea.list[[i]]@result[gsea.up_unique$gsea.list[[i]]@result$ID %in%
                                                                                unq, ]
}

gsea.down_unique = gsea.down
for (i in 1:length(unique_downexpr_pathways)) {
  unq = unique_downexpr_pathways[[i]]$V1
  gsea.down_unique$gsea.list[[i]]@result = gsea.down_unique$gsea.list[[i]]@result[gsea.down_unique$gsea.list[[i]]@result$ID %in%
                                                                                    unq, ]
}

rm(unq); gc()

hclust_input_up = prepare_gsea_output_for_hclust(gsea_output = gsea.up_unique, 
                                                 dgea_output_name_style = "dgea_", 
                                                 dea.method = "limma", 
                                                 mo.method = "",
                                                 dat.path = paste0(home, "/Results/single_algorithm/", algorithm), 
                                                 dgea_padj_cutoff = 0.05,
                                                 logfc_cutoff = 0,
                                                 pathway_padj_cutoff = 0.05)

hclust_input_down = prepare_gsea_output_for_hclust(gsea_output = gsea.down_unique, 
                                                   dgea_output_name_style = "dgea_", 
                                                   dea.method = "limma", 
                                                   mo.method = "",
                                                   dat.path = paste0(home, "/Results/single_algorithm/", algorithm), 
                                                   dgea_padj_cutoff = 0.05,
                                                   logfc_cutoff = 0,
                                                   pathway_padj_cutoff = 0.05)

hclust_input = c(hclust_input_up, hclust_input_down)
names(hclust_input) = c(paste0(rep("up_", length(hclust_input_up)), 
                               names(hclust_input_up)),
                        paste0(rep("down_", length(hclust_input_down)), 
                               names(hclust_input_down)))

rm(hclust_input_down, hclust_input_up); gc()

# Load doParallel if not already loaded
library(parallel)
library(foreach)
library(doParallel)

# Set up the number of cores to use: minimum of length(hclust_input) or 5
cl <- makeCluster(min(length(hclust_input), 5))
registerDoParallel(cl)

# Use foreach with parallel processing
timestamp()
hclust_output <- foreach(i = 1:length(hclust_input), .packages = c("pathfindR", "fastcluster")) %dopar% {
  RNGversion("4.2.2")
  set.seed(123)
  source("Scripts/automated_scripts/fast_pathfindR_hclust.R")
  
  # Perform clustering, handle errors
  result <- cluster_enriched_terms_fast(hclust_input[[i]], method = "hierarchical", plot_clusters_graph = FALSE,
                                        use_description = FALSE, use_active_snw_genes = FALSE)
  if (is.character(result) && result == "hclust impossible") {
    return("hclust impossible")
  } else {
    return(result)
  }
}

timestamp() # ~2.5 mins
stopCluster(cl)
gc()
names(hclust_output) <- names(hclust_input)

# Are there any null sets?
which(hclust_output == "hclust impossible")

# Export
library(openxlsx)
wb = createWorkbook()
for (j in 1:length(hclust_output)) {
  addWorksheet(wb, names(hclust_output)[j])
  writeData(wb, names(hclust_output)[j], hclust_output[[j]])
}
saveWorkbook(wb, file = paste0(home, "/Results/single_algorithm/", algorithm, "/", 
                               algorithm, "_representative_pathways.xlsx"),
             overwrite = TRUE); rm(wb)

# Plot pathway heatmaps
hclust_pathway_plots_up = plot_pathway_heatmaps(gsea.lists = hclust_output[grepl("up", names(hclust_output))], 
                                                norm.expr = plotdata$RNAseq, 
                                                present_clusters = c("iClusterBayes1", "iClusterBayes2"),
                                                representative = TRUE, moic.res = plot_object,
                                                subtype_prefix = algorithm, n.path = 20, msigdb.path = MSIGDB.FILE,
                                                norm.method = "mean", dirct = "up",
                                                fig.name = "upregulated_pathway_heatmap",
                                                name = "GSVA scores",
                                                fig.path = paste0(home, "/Results/single_algorithm/", algorithm), 
                                                width = 15, height = 10, gsva.method = "gsva")

hclust_pathway_plots_down = plot_pathway_heatmaps(gsea.lists = hclust_output[grepl("down", names(hclust_output))], 
                                                  norm.expr = plotdata$RNAseq, 
                                                  present_clusters = c("iClusterBayes1", "iClusterBayes2"),
                                                  representative = TRUE, moic.res = plot_object,
                                                  subtype_prefix = algorithm, n.path = 20, msigdb.path = MSIGDB.FILE,
                                                  norm.method = "mean", dirct = "down",
                                                  fig.name = "downregulated_pathway_heatmap",
                                                  name = "GSVA scores",
                                                  fig.path = paste0(home, "/Results/single_algorithm/", algorithm), 
                                                  width = 15, height = 10, gsva.method = "gsva")

# Fraction Genome Altered ###
fga_df = readRDS("Resources/TCGA/fga_df.rds"); gc()

fga.iClusterBayes <- compFGA_optimized(moic.res     = plot_object,
                             segment      = fga_df,
                             iscopynumber = TRUE, 
                             test.method  = "nonparametric", # statistical testing method (Wilcoxon with asymptotic approximation. Consider Kruskall Wallis?)
                             fig.path     = paste0(home, "/Results/single_algorithm/", algorithm),
                             fig.name     = paste0("FGA_barplot_", algorithm),
                             prefix = algorithm,
                             width = 16,
                             ga_column = "ga", # genome altered column
                             clust.col = cluster_colors,
                             title = paste0(algorithm, " FGA plot: simple criteria"))

fga.iClusterBayes.COSMIC <- compFGA_optimized(moic.res     = plot_object,
                                    segment      = fga_df,
                                    iscopynumber = TRUE, 
                                    test.method  = "nonparametric", # statistical testing method (Wilcoxon with asymptotic approximation. Consider Kruskall Wallis?)
                                    fig.path     = paste0(home, "/Results/single_algorithm/", algorithm),
                                    fig.name     = paste0("COSMIC_criteria_FGA_barplot_", algorithm),
                                    prefix = algorithm,
                                    width = 16,
                                    ga_column = "COSMIC_ga", # genome altered column
                                    clust.col = cluster_colors,
                                    title = paste0(algorithm, " FGA plot: COSMIC criteria"))

rm(fga_df); gc()

# Evaluation #####
# Run Nearest Template Prediction in transNEO cohort ###
# Load transNEO data
transNEO_mm_inputs = readRDS("Resources/transNEO/transNEO_multimodal_inputs.rds")
transcr = readRDS("Resources/transNEO/log2.norm.counts.plus1_transNEO.rds")

# get as many templates as possible
dgea.marker.up_1000 <- runMarker_single_algorithm_no_export(algorithm_name = algorithm,
                                                            moic.res = plot_object,
                                                            n.marker = 1000,
                                                            dea.method    = "limma", # name of DEA method
                                                            prefix        = "dgea_", # MUST be the same of argument in runDEA()
                                                            dat.path      = paste0(home, "/Results/single_algorithm/", algorithm), # path of DEA files
                                                            p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                                                            p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                                                            norm.expr = plotdata$RNAseq,
                                                            dirct         = "up" # direction of dysregulation in expression
)

# 2. Down-regulated markers
dgea.marker.down_1000 <- runMarker_single_algorithm_no_export(algorithm_name = algorithm,
                                                              moic.res = plot_object,
                                                              n.marker = 1000,
                                                              dea.method    = "limma", # name of DEA method
                                                              prefix        = "dgea_", # MUST be the same of argument in runDEA()
                                                              dat.path      = paste0(home, "/Results/single_algorithm/", algorithm), # path of DEA files
                                                              p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                                                              p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                                                              norm.expr = plotdata$RNAseq,
                                                              dirct         = "down" # direction of dysregulation in expression
)

# Up-regulated expression features
RNGversion("4.2.2")
timestamp()
transNEO_ntp_expr_up = runNTP(
  expr = transcr,
  templates = dgea.marker.up_1000$templates,
  scaleFlag = TRUE,
  centerFlag = TRUE,
  nPerm = 10000,
  seed = 123,
  distance = "cosine", # default
  doPlot = TRUE,
  height = 8,
  width = 12,
  fig.path = paste0(home, "/Results/single_algorithm/", algorithm),
  fig.name = "ntp_expr_up_heatmap_transNEO")
timestamp() # 4 min

# down-regulated
RNGversion("4.2.2")
timestamp()
transNEO_ntp_expr_down = runNTP(
  expr = transcr,
  templates = dgea.marker.down_1000$templates,
  scaleFlag = TRUE,
  centerFlag = TRUE,
  nPerm = 10000,
  seed = 123,
  distance = "cosine", # default
  doPlot = TRUE,
  height = 8,
  width = 12,
  fig.path = paste0(home, "/Results/single_algorithm/", algorithm),
  fig.name = "ntp_expr_down_heatmap_transNEO")
timestamp() # 2.5 min

# Check concordance
expr_conc = as.data.frame(transNEO_ntp_expr_down$clust.res) %>%
  dplyr::rename(clust_down = clust) %>%
  inner_join(as.data.frame(transNEO_ntp_expr_up$clust.res) %>%
               dplyr::rename(clust_up = clust),
             by = "samID")

# This is counter-intuitive but due to opposite directions of deregulation this
# is how it works (perhaps this was expected)
expr_conc$agreement = ifelse(expr_conc$clust_down!=expr_conc$clust_up, "Yes", "No")
paste("Agremeent of NTP subtypes with respect to expression data from the external cohort is: ",
      length(which(expr_conc$agreement == "Yes"))/nrow(expr_conc)*100, "% (", nrow(expr_conc),
      " samples).")

# Compare clinical variables of interest across clusters
transNEO_var2comp = transNEO_mm_inputs$`Full pheno` %>%
  dplyr::select(LN.status.at.diagnosis, ER.status, HER2.status,
                Grade.pre.NAT, pCR.RD, Age, T.stage, PAM50, iC10,
                NAT.regimen, Chemo.cycles,
                aHER2.cycles, RCB.score, STAT1.gsva,
                GGI.gsva, ESC.gsva, TMB, HRD.sum, Donor.ID) %>%
  inner_join(expr_conc %>% dplyr::select(Donor.ID = samID, iClusterBayes = clust_up),
             by = "Donor.ID")
rownames(transNEO_var2comp) = transNEO_var2comp$Donor.ID
transNEO_var2comp = transNEO_var2comp %>% dplyr::select(-Donor.ID)

# Convert to factors
transNEO_var2comp$LN.status.at.diagnosis = factor(transNEO_var2comp$LN.status.at.diagnosis,
                                                  levels = c("NEG", "POS"),
                                                  labels = c("Negative", "Positive"))
transNEO_var2comp$ER.status = factor(transNEO_var2comp$ER.status,
                                     levels = c("NEG", "POS"),
                                     labels = c("Negative", "Positive"))
transNEO_var2comp$HER2.status = factor(transNEO_var2comp$HER2.status,
                                       levels = c("NEG", "POS"),
                                       labels = c("Negative", "Positive"))
transNEO_var2comp$Grade.pre.NAT = factor(transNEO_var2comp$Grade.pre.NAT,
                                         levels = c(1, 2, 3, 4),
                                         labels = c("Grade 1", "Grade 2", "Grade 3", "Grade 4"))
transNEO_var2comp$pCR.RD = factor(transNEO_var2comp$pCR.RD,
                                  levels = c("pCR", "RD"),
                                  labels = c("pCR", "Residual Disease"))
transNEO_var2comp$PAM50 = factor(transNEO_var2comp$PAM50,
                                 levels = c("Basal", "Her2", "LumB", "LumA", "Normal", "Unk"),
                                 labels = c("Basal-like", "HER2+", "Luminal B", "Luminal A",
                                            "Normal-like", "Unknown"))
transNEO_var2comp$iC10 = factor(transNEO_var2comp$iC10,
                                levels = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
                                labels = paste("iC", seq(1, 10, 1), sep = ""))

# Remove unknown levels for statistical tests
transNEO_var2comp_nonas = transNEO_var2comp
for (i in 1:ncol(transNEO_var2comp)) {
  nas = which(transNEO_var2comp[, i] == "Unknown")
  transNEO_var2comp_nonas[nas, i] = NA
  empties = which(transNEO_var2comp[, i] == "")
  transNEO_var2comp_nonas[empties, i] = NA
}
rm(nas, empties); gc()

transNEO_clincomp = compClinvar_single_algorithm(algorithm_name = algorithm,
                                                 moic.res = transNEO_ntp_expr_up,
                                                 var2comp = transNEO_var2comp_nonas,
                                                 strata = algorithm,
                                                 factorVars = c("ER.status", "HER2.status",
                                                                "NAT.regimen", 
                                                                "pCR.RD", "LN.status.at.diagnosis"),
                                                 nonnormalVars = c("Age",
                                                                   "RCB.score", "STAT1.gsva", "GGI.gsva",
                                                                   "ESC.gsva", "TMB", "HRD.sum",
                                                                   "Grade.pre.NAT", "Chemo.cycle", "aHER2.cycles"),
                                                 includeNA = FALSE,
                                                 doWord = TRUE,
                                                 tab.name = "transNEO_Summary_of_clinical_variables",
                                                 res.path = paste0(home, "/Results/single_algorithm/", algorithm, "/"),
                                                 output_pdf = TRUE,
                                                 pdf_level_col_width = c("7em", "10em"),
                                                 pdf_count_col_width = "10em",
                                                 pdf_pval_col_width = "3em",
                                                 pdf_test_col_width = "8em",
                                                 pdf_tab_font_size = 9)

transNEO_ntp_expr_up_ord = transNEO_ntp_expr_up
transNEO_ntp_expr_up_ord$clust.res$clust = gsub(algorithm, "", transNEO_ntp_expr_up_ord$clust.res$clust)
transNEO_ordinal_clincomp = compClinvar_ordinal_single_algorithm(algorithm_name = algorithm,
                                                                 moic.res = transNEO_ntp_expr_up_ord,
                                                                 var2comp = transNEO_var2comp_nonas %>%
                                                                   dplyr::select(Grade.pre.NAT, 
                                                                                 Chemo.cycles, 
                                                                                 aHER2.cycles,
                                                                                 iClusterBayes),
                                                                 strata = algorithm,
                                                                 ordinalVars = c("Grade.pre.NAT",
                                                                                 "Chemo.cycles",
                                                                                 "aHER2.cycles"),
                                                                 includeNA = FALSE,
                                                                 tab.name = "transNEO Summary of ordinal clinical variables",
                                                                 res.path = paste0(home, "/Results/single_algorithm/", algorithm, "/"),
                                                                 output_pdf = TRUE,
                                                                 pdf_template_loc = paste0(home, "/Scripts/automated_scripts/clincomp_template.Rmd"),
                                                                 pdf_level_col_width = c("7em", "10em"),
                                                                 pdf_count_col_width = "10em",
                                                                 pdf_pval_col_width = "3em",
                                                                 pdf_test_col_width = "8em",
                                                                 pdf_tab_font_size = 9)

# Run PAM ###
RNGversion("4.2.2.")
set.seed(123)
transNEO_pam = runPAM_single_algorithm(algorithm_name = algorithm,
                                       train.expr = plotdata$RNAseq,
                                       moic.res   = plot_object,
                                       test.expr  = transcr)

# Check consistency across methods

# Get predictions for TCGA (discovery cohort)
RNGversion("4.2.2.")
set.seed(123)
TCGA.ntp.pred = runNTP(expr = plotdata$RNAseq[, plot_object$clust.res$samID],
                       templates = dgea.marker.up_1000$templates, distance = "cosine",
                       doPlot = F, nPerm = 10000)

TCGA.pam.pred = runPAM_single_algorithm(algorithm_name = algorithm,
                                        train.expr = plotdata$RNAseq[, plot_object$clust.res$samID],
                                        moic.res = plot_object,
                                        test.expr = plotdata$RNAseq[, plot_object$clust.res$samID])

# consensus TCGA vs NTP TCGA
runKappa_single_algorithm(algorithm_name = algorithm,
                          subt1 = plot_object$clust.res$clust,
                          subt2 = gsub(algorithm, "", TCGA.ntp.pred$clust.res$clust),
                          subt1.lab = algorithm,
                          subt2.lab = "NTP TCGA",
                          height = 8,
                          width = 8,
                          fig.path = paste0(home, "/Results/single_algorithm/", algorithm),
                          fig.name = paste0("kappa_", algorithm, "_vs_NTP_TCGA"))

# consensus TCGA vs PAM TCGA
runKappa_single_algorithm(algorithm_name = algorithm,
                          subt1 = plot_object$clust.res$clust,
                          subt2 = gsub(algorithm, "", TCGA.pam.pred$clust.res$clust),
                          subt1.lab = algorithm,
                          subt2.lab = "PAM TCGA",
                          height = 8,
                          width = 8,
                          fig.path = paste0(home, "/Results/single_algorithm/", algorithm),
                          fig.name = paste0("kappa_", algorithm, "_vs_PAM_TCGA"))

# NTP transNEO vs PAM transNEO
runKappa_single_algorithm(algorithm_name = algorithm,
                          subt1 = as.numeric(gsub(algorithm, "",
                                                  transNEO_ntp_expr_up$clust.res$clust)),
                          subt2 = as.numeric(transNEO_pam$clust.res$clust),
                          subt1.lab = "transNEO NTP",
                          subt2.lab = "transNEO PAM",
                          height = 8,
                          width = 8,
                          fig.path = paste0(home, "/Results/single_algorithm/", algorithm),
                          fig.name = "kappa_NTP_vs_PAM_transNEO")

# Supplementary results #####

# Create subdirectory for supplementary plots
if (!dir.exists(paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement"))) {
  dir.create(paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement"))
}

# Setup for heatmaps
colors_heatmap = rev(colorRampPalette(viridisLite::magma(10))(255))
cluster_colors_heatmap = c("#2EC4B6", "#E71D36")
clust_annot_pheno = annCol %>% mutate(Sample.ID = rownames(.)) %>%
  inner_join(clust, by = "Sample.ID") %>%
  dplyr::rename(iClusterBayes = Cluster, samID = "Sample.ID")
rownames(clust_annot_pheno) = clust_annot_pheno$samID
afh_colnames = colnames(annCol)

# Prepare affinity matrices
aff_CNV = normalize_affinity_matrix(
  iClusterBayestool::affinityMatrix(
    iClusterBayestool::dist2(input$CNV,
                   input$CNV),
    K = optN, sigma = optSigma))
colnames(aff_CNV) = rownames(aff_CNV) = rownames(input$CNV)

aff_rna = normalize_affinity_matrix(
  iClusterBayestool::affinityMatrix(
    iClusterBayestool::dist2(input$RNAseq,
                   input$RNAseq),
    K = optN, sigma = optSigma))
colnames(aff_rna) = rownames(aff_rna) = rownames(input$RNA)

aff_miRNA = normalize_affinity_matrix(
  iClusterBayestool::affinityMatrix(
    iClusterBayestool::dist2(input$miRNA,
                   input$miRNA),
    K = optN, sigma = optSigma))
colnames(aff_miRNA) = rownames(aff_miRNA) = rownames(input$miRNA)

aff_Methyl = normalize_affinity_matrix(
  iClusterBayestool::affinityMatrix(
    iClusterBayestool::dist2(input$Methylation,
                   input$Methylation),
    K = optN, sigma = optSigma))
colnames(aff_Methyl) = rownames(aff_Methyl) = rownames(input$Methylation)

aff_SNPs = normalize_affinity_matrix(
  iClusterBayestool::affinityMatrix(
    as.matrix(dist(as.matrix(input$SNPs),
                   as.matrix(input$SNPs),
                   method = "binary")),
    K = optN, sigma = optSigma))
colnames(aff_SNPs) = rownames(aff_SNPs) = rownames(input$SNPs)

aff_final = final_affinity_matrix

# CNV
create_MO_heatmap(matrix = aff_CNV, algorithm = algorithm, 
                  need.diag.zero = TRUE, 
                  clust_annot_pheno = clust_annot_pheno,
                  afh_colnames = afh_colnames, 
                  colors = colors_heatmap,
                  annColors = annColors,
                  heatmap_title = "CNV first affinity heatmap",
                  cluster_colors = cluster_colors_heatmap,
                  legend_title = "Normalized affinity",
                  cluster_cols_flag = FALSE,
                  cluster_rows_flag = FALSE,
                  splits_flag = TRUE,
                  output_file_name = paste0(home, "/Results/single_algorithm/", algorithm, 
                                            "/Supplement/aff_CNV_heatmap.png"))

# RNAseq
create_MO_heatmap(matrix = aff_rna, algorithm = algorithm, 
                  need.diag.zero = TRUE, 
                  clust_annot_pheno = clust_annot_pheno,
                  afh_colnames = afh_colnames, 
                  colors = colors_heatmap,
                  annColors = annColors,
                  heatmap_title = "RNAseq first affinity heatmap",
                  cluster_colors = cluster_colors_heatmap,
                  legend_title = "Normalized affinity",
                  cluster_cols_flag = FALSE,
                  cluster_rows_flag = FALSE,
                  splits_flag = TRUE,
                  output_file_name = paste0(home, "/Results/single_algorithm/", algorithm, 
                                            "/Supplement/aff_RNAseq_heatmap.png"))

# miRNA
create_MO_heatmap(matrix = aff_miRNA, algorithm = algorithm, 
                  need.diag.zero = TRUE, 
                  clust_annot_pheno = clust_annot_pheno,
                  afh_colnames = afh_colnames, 
                  colors = colors_heatmap,
                  annColors = annColors,
                  heatmap_title = "miRNA first affinity heatmap",
                  cluster_colors = cluster_colors_heatmap,
                  legend_title = "Normalized affinity",
                  cluster_cols_flag = FALSE,
                  cluster_rows_flag = FALSE,
                  splits_flag = TRUE,
                  output_file_name = paste0(home, "/Results/single_algorithm/", algorithm,
                                            "/Supplement/aff_miRNA_heatmap.png"))

# Methylation
create_MO_heatmap(matrix = aff_Methyl, algorithm = algorithm, 
                  need.diag.zero = TRUE, 
                  clust_annot_pheno = clust_annot_pheno,
                  afh_colnames = afh_colnames, 
                  colors = colors_heatmap,
                  annColors = annColors,
                  heatmap_title = "Methylation first affinity heatmap",
                  cluster_colors = cluster_colors_heatmap,
                  legend_title = "Normalized affinity",
                  cluster_cols_flag = FALSE,
                  cluster_rows_flag = FALSE,
                  splits_flag = TRUE,
                  output_file_name = paste0(home, "/Results/single_algorithm/", algorithm, 
                                            "/Supplement/aff_Methylation_heatmap.png"))

# SNPs
create_MO_heatmap(matrix = aff_SNPs, algorithm = algorithm, 
                  need.diag.zero = TRUE, 
                  clust_annot_pheno = clust_annot_pheno,
                  afh_colnames = afh_colnames, 
                  colors = colors_heatmap,
                  annColors = annColors,
                  heatmap_title = "SNPs first affinity heatmap",
                  cluster_colors = cluster_colors_heatmap,
                  legend_title = "Normalized affinity",
                  cluster_cols_flag = FALSE,
                  cluster_rows_flag = FALSE,
                  splits_flag = TRUE,
                  output_file_name = paste0(home, "/Results/single_algorithm/", algorithm,
                                            "/Supplement/aff_SNPs_heatmap.png"))

# Final affinity matrix
create_MO_heatmap(matrix = final_affinity_matrix, algorithm = algorithm, 
                  need.diag.zero = TRUE, 
                  clust_annot_pheno = clust_annot_pheno,
                  afh_colnames = afh_colnames, 
                  colors = colors_heatmap,
                  annColors = annColors,
                  heatmap_title = "Final affinity heatmap",
                  cluster_colors = cluster_colors_heatmap,
                  legend_title = "Normalized affinity",
                  cluster_cols_flag = FALSE,
                  cluster_rows_flag = FALSE,
                  splits_flag = TRUE,
                  output_file_name = paste0(home, "/Results/single_algorithm/", algorithm,
                                            "/Supplement/aff_final_affinity_heatmap.png"))

# Final affinity matrix with clustered rows and columns
create_MO_heatmap(matrix = final_affinity_matrix, algorithm = algorithm, 
                  need.diag.zero = TRUE, 
                  clust_annot_pheno = clust_annot_pheno,
                  afh_colnames = afh_colnames, 
                  colors = colors_heatmap,
                  annColors = annColors,
                  heatmap_title = "Final affinity heatmap",
                  cluster_colors = cluster_colors_heatmap,
                  legend_title = "Normalized affinity",
                  cluster_cols_flag = TRUE,
                  cluster_rows_flag = TRUE,
                  splits_flag = FALSE,
                  output_file_name = paste0(home, "/Results/single_algorithm/", algorithm,
                                            "/Supplement/hclust_aff_final_affinity_heatmap.png"))

# PCA ###
# CNV
pca_from_sim_matrix(sim_matrix = aff_CNV, algorithm = algorithm, 
                    clust_res = clust_annot_pheno %>% dplyr::select(samID, iClusterBayes),
                    cluster_colors = cluster_colors_heatmap, 
                    output_path = paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement"), 
                    title_add = "CNV")

# RNAseq
pca_from_sim_matrix(sim_matrix = aff_rna, algorithm = algorithm, 
                    clust_res = clust_annot_pheno %>% dplyr::select(samID, iClusterBayes),
                    cluster_colors = cluster_colors_heatmap, 
                    output_path = paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement"), 
                    title_add = "RNAseq")

# miRNA
pca_from_sim_matrix(sim_matrix = aff_miRNA, algorithm = algorithm, 
                    clust_res = clust_annot_pheno %>% dplyr::select(samID, iClusterBayes),
                    cluster_colors = cluster_colors_heatmap, 
                    output_path = paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement"), 
                    title_add = "miRNA")

# Methylation
pca_from_sim_matrix(sim_matrix = aff_Methyl, algorithm = algorithm, 
                    clust_res = clust_annot_pheno %>% dplyr::select(samID, iClusterBayes),
                    cluster_colors = cluster_colors_heatmap, 
                    output_path = paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement"), 
                    title_add = "Methylation")

# SNPs
pca_from_sim_matrix(sim_matrix = aff_SNPs, algorithm = algorithm, 
                    clust_res = clust_annot_pheno %>% dplyr::select(samID, iClusterBayes),
                    cluster_colors = cluster_colors_heatmap, 
                    output_path = paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement"), 
                    title_add = "SNPs")

# Final affinity
pca_from_sim_matrix(sim_matrix = aff_final, algorithm = algorithm, 
                    clust_res = clust_annot_pheno %>% dplyr::select(samID, iClusterBayes),
                    cluster_colors = cluster_colors_heatmap, 
                    output_path = paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement"), 
                    title_add = "Fusion")

# Setup for barcharts ###
# Stage
scale_fill_stage = scale_fill_manual(values = c(`Stage I` = "#00C9FF", 
                                                `Stage II` = "#099CF5", 
                                                `Stage III` = "#097BF5", 
                                                `Stage IV` = "#0B5684", 
                                                `Unknown` = "grey40"))

# Lymph node status
scale_fill_lymph_node_status = scale_fill_manual(values = c(No = "grey75", 
                                                            Yes = "#4A0558", 
                                                            Unknown = "grey40"))

# ER status
scale_fill_ER_status = scale_fill_manual(values = c(Negative = "#C11D9C", 
                                                    Positive = "#0F1682", 
                                                    Unknown = "grey40"))

# PR status
scale_fill_PR_status = scale_fill_manual(values = c(Indeterminate = "aliceblue", 
                                                    Positive = "dodgerblue4", 
                                                    Negative = "#F0C6C3", 
                                                    Unknown = "grey40"))

# HER2 status
scale_fill_HER2_status = scale_fill_manual(values = c(Negative = "#0B9EF8", 
                                                      Positive = "#560DA7", 
                                                      Indeterminate = "mistyrose1", 
                                                      Equivocal = "hotpink4", 
                                                      Unknown = "grey40"))

# Vital status
scale_fill_vital_status = scale_fill_manual(values = c(Alive = "lightpink1", 
                                                       Dead = "black", 
                                                       Unknown = "grey40"))

# Ethnicity
scale_fill_ethnicity = scale_fill_manual(values = c(`Hispanic or latino` = "#E58606", 
                                                    `Not hispanic or latino` = "#24796C", 
                                                    Unknown = "grey40"))

# Race
scale_fill_race = scale_fill_manual(values = c(`American indian or alaska native` = "#E73F74", 
                                               Asian = "#3969AC", 
                                               `Black or african american` = "#666666", 
                                               White = "beige", 
                                               Unknown = "grey40"))

# Metastasis
scale_fill_metastasis = scale_fill_manual(values = c(Yes = "deeppink4", 
                                                     No = "cadetblue2", 
                                                     Unknown = "grey40"))

# Histology
scale_fill_histology = scale_fill_manual(values = c(`Infiltrating Carcinoma NOS` = "#88CCEE", 
                                                    `Infiltrating Ductal Carcinoma` = "#CC6677", 
                                                    `Infiltrating Lobular Carcinoma` = "#DDCC77", 
                                                    `Medullary Carcinoma` = "#117733", 
                                                    `Metaplastic Carcinoma` = "#332288", 
                                                    Mixed = "#AA4499", 
                                                    `Mucinous Carcinoma` = "#44AA99", 
                                                    Other = "#999933", 
                                                    Unknown = "grey40"))

# Menopausal status
scale_fill_menopausal_status = scale_fill_manual(values = c(Indeterminate = "mistyrose2", 
                                                            `Pre-menopausal` = "#FAA476", 
                                                            Perimenopausal = "#DC3977", 
                                                            `Post-menopausal` = "#7C1D6F", 
                                                            Unknown = "grey40"))

# Combine all scales into a list
barchart_scales = list(scale_fill_stage, scale_fill_lymph_node_status, scale_fill_ER_status, 
                       scale_fill_PR_status, scale_fill_HER2_status, scale_fill_vital_status, 
                       scale_fill_ethnicity, scale_fill_race, scale_fill_metastasis, 
                       scale_fill_histology, scale_fill_menopausal_status)

# Name the scales accordingly
names(barchart_scales) = c("Stage", "Lymph node status", "ER status", "PR status", "HER2 status", 
                           "Vital status", "Ethnicity", "Race", "Metastasis", "Histology", 
                           "Menopausal status")

# Chi-square tests ###
# Bias-corrected Cramer's V calculation using package rcompanion:
unbiased.cv.test = function(x, string, digits = 3) {
  CV = rcompanion::cramerV(x, bias.correct = TRUE)
  return(list(text = paste0("Bias-corrected Cramer's V / Phi for ", 
                            string, ": ", round(as.numeric(CV), digits)),
              value = round(as.numeric(CV), digits)))
}

clust_annot_pheno_nonas = clust_annot_pheno
for(i in 1:ncol(clust_annot_pheno_nonas)) {
  clust_annot_pheno_nonas[, i] = as.character(clust_annot_pheno_nonas[, i])
  nas = which(clust_annot_pheno_nonas[, i] == "Unknown")
  clust_annot_pheno_nonas[nas, i] = NA
  clust_annot_pheno_nonas[, i] = factor(clust_annot_pheno_nonas[, i])
}
rm(nas); gc()

voi = colnames(clust_annot_pheno_nonas)[1:11]
output = as.data.frame(matrix(NA, nrow = 0, ncol = 4))
for (v in 1:length(voi)){
  keepers = which(!is.na(clust_annot_pheno_nonas[, voi[v]]))
  test = suppressWarnings(chisq.test(table(clust_annot_pheno_nonas[keepers, algorithm], 
                                           clust_annot_pheno_nonas[keepers, voi[v]])))
  chifit_p = test$p.value
  chifit_xsq = test$statistic
  chifit_cv = suppressWarnings(unbiased.cv.test(table(clust_annot_pheno_nonas[keepers, algorithm], 
                                                      clust_annot_pheno_nonas[keepers, voi[v]]),
                                                string = voi[v],
                                                digits = 3)$value)
  comparison = paste0(voi[v], " vs ", algorithm, " cluster")
  output = rbind(output, c(comparison, chifit_p, chifit_xsq, chifit_cv))
  rm(test, comparison, chifit_p, chifit_xsq, chifit_cv, keepers)
}
colnames(output) = c("Comparison", "p-value", "Statistic", "Cramer's V")

rm(v); gc()
openxlsx::write.xlsx(output, 
                     paste0(home, 
                            "/Results/single_algorithm/", algorithm, "/Supplement/Chisq_tests.xlsx"),
                     overwrite = TRUE)

# Bar chart generation
iClusterBayes_barcharts = list()
plotdata_bar = clust_annot_pheno_nonas
plotdata_bar[[algorithm]] = factor(plotdata_bar[[algorithm]])
for (i in 1:length(voi)) {
  chifit = output
  loc = which(grepl(voi[i], chifit$Comparison))
  chifit = chifit[loc, ]
  iClusterBayes_barcharts[[i]] = create_annot_barchart(plotdata = plotdata_bar, fill = voi[i],
                                             chifit = chifit,
                                             na.action = "na.omit",
                                             algorithm = algorithm,
                                             barchart_ylim = 650,
                                             text_y = 600, rect_ymin = 500,
                                             rect_ymax = 620, x_annot = 1.5,
                                             v_gap = 35, rect_xmin = 1,
                                             rect_xmax = 2, 
                                             annot_text_size = 2.25,
                                             legend.text.size = 5,
                                             x.axis.text.size = 5) +
    barchart_scales[[voi[i]]]
  print(iClusterBayes_barcharts[[i]])
  ggsave(filename = paste0(algorithm, "_", voi[i], "_barchart.png"),
         path = paste0(home, 
                       "/Results/single_algorithm/", algorithm, "/Supplement"), 
         width = 2320, height = 2320, device = 'png', units = "px",
         dpi = 700)
  dev.off()
}
names(iClusterBayes_barcharts) = voi
rm(loc, chifit)

# Multiplot (PNG) - bar charts
library(ggpubr)
ggarrange(iClusterBayes_barcharts[[1]], iClusterBayes_barcharts[[2]], iClusterBayes_barcharts[[3]],
          iClusterBayes_barcharts[[4]], iClusterBayes_barcharts[[5]], iClusterBayes_barcharts[[6]],
          iClusterBayes_barcharts[[7]], iClusterBayes_barcharts[[8]], iClusterBayes_barcharts[[9]],
          iClusterBayes_barcharts[[10]], iClusterBayes_barcharts[[11]],
          ncol = 3, nrow = 4, labels = c("A", "B", "C", "D", "E", "F", "G", "H",
                                         "I", "J", "K"),
          font.label = list(size = 8, face = "bold", color ="black"))
ggsave(filename = paste0("Multiplot_", algorithm, "_barcharts.png"),
       path = paste0(home, 
                     "/Results/single_algorithm/", algorithm, "/Supplement"), 
       width = 7000, height = 8000, device = 'png', units = "px",
       dpi = 700)
dev.off()

# Just significant ones now
iClusterBayes_barcharts_sig = list()
plotdata_bar_sig = clust_annot_pheno_nonas %>% dplyr::select(iClusterBayes, Race, Histology, 
                                                             `ER status`, `PR status`)
plotdata_bar_sig$iClusterBayes = factor(plotdata_bar_sig$iClusterBayes)
voi_sig = setdiff(colnames(plotdata_bar_sig), algorithm)
for (i in 1:length(voi_sig)) {
  chifit = output
  loc = which(grepl(voi_sig[i], chifit$Comparison))
  chifit = chifit[loc, ]
  iClusterBayes_barcharts_sig[[i]] = create_annot_barchart(plotdata = plotdata_bar_sig, fill = voi_sig[i],
                                                 chifit = chifit,
                                                 na.action = "na.omit",
                                                 algorithm = algorithm,
                                                 barchart_ylim = 650,
                                                 text_y = 600, rect_ymin = 500,
                                                 rect_ymax = 620, x_annot = 1.5,
                                                 v_gap = 35, rect_xmin = 1,
                                                 rect_xmax = 2, 
                                                 annot_text_size = 2.25,
                                                 legend.text.size = 5,
                                                 x.axis.text.size = 5) +
    barchart_scales[[voi_sig[i]]]
  print(iClusterBayes_barcharts_sig[[i]])
  ggsave(filename = paste0("sig_", algorithm, "_", voi_sig[i], "_barchart.png"),
         path = paste0(home, 
                       "/Results/single_algorithm/", algorithm, "/Supplement"), 
         width = 2320, height = 2320, device = 'png', units = "px",
         dpi = 700)
  dev.off()
}
names(iClusterBayes_barcharts_sig) = voi_sig
rm(loc, chifit)

# Multiplot (PNG) - bar charts
ggarrange(iClusterBayes_barcharts_sig[[1]], iClusterBayes_barcharts_sig[[2]], iClusterBayes_barcharts_sig[[3]],
          iClusterBayes_barcharts_sig[[4]], 
          ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"),
          font.label = list(size = 8, face = "bold", color ="black"))
ggsave(filename = paste0("sig_Multiplot_", algorithm, "_barcharts.png"),
       path = paste0(home, 
                     "/Results/single_algorithm/", algorithm, "/Supplement"), 
       width = 5500, height = 5500, device = 'png', units = "px",
       dpi = 700)
dev.off()

# Sunburst plot ###
library(plotly)
Pheno_sunburst_iClusterBayes = clust_annot_pheno
Pheno_sunburst_iClusterBayes$`ER status` = gsub("Unknown", "Unkn ER status", Pheno_sunburst_iClusterBayes$`ER status`)
Pheno_sunburst_iClusterBayes$`ER status` = gsub("Positive", "ER+", Pheno_sunburst_iClusterBayes$`ER status`)
Pheno_sunburst_iClusterBayes$`ER status` = gsub("Negative", "ER-", Pheno_sunburst_iClusterBayes$`ER status`)
Pheno_sunburst_iClusterBayes = Pheno_sunburst_iClusterBayes %>%
  dplyr::select(iClusterBayes, `ER status`) %>%
  group_by(iClusterBayes, `ER status`) %>%
  summarise(Counts = n()) %>%
  as.data.frame()

sunburst_coloring_iClusterBayes = data.frame(stringsAsFactors = FALSE,
                                   colors = tolower(gplots::col2hex(c("#2EC4B6", "#E71D36", 
                                                                      "#C11D9C", "#0F1682",  "grey40"))),
                                   labels = c("iClusterBayes1", "iClusterBayes2",
                                              "ER-", "ER+", "Unkn ER status"))

sunburstDF_iClusterBayes = as.sunburstDF(Pheno_sunburst_iClusterBayes, value_column = "Counts", add_root = FALSE) %>%
  inner_join(sunburst_coloring_iClusterBayes, by = "labels")

pie_iClusterBayes = plot_ly() %>%
  add_trace(ids = sunburstDF_iClusterBayes$ids, labels= sunburstDF_iClusterBayes$labels, 
            parents = sunburstDF_iClusterBayes$parents, 
            values= sunburstDF_iClusterBayes$values, type='sunburst', branchvalues = 'total',
            insidetextorientation='radial', maxdepth = 5,
            marker = list(colors = sunburstDF_iClusterBayes$colors)) %>%
  layout(
    grid = list(columns =1, rows = 1),
    margin = list(l = 0, r = 0, b = 0, t = 0)
  )
pie_iClusterBayes
rm(Pheno_sunburst_iClusterBayes, sunburstDF_iClusterBayes, sunburst_coloring_iClusterBayes, pie_iClusterBayes); gc()

# Graphs ###
library(igraph)

aff_CNV_S = calculate_S(aff_CNV)
aff_rna_S = calculate_S(aff_rna)
aff_miRNA_S = calculate_S(aff_miRNA)
aff_Methyl_S = calculate_S(aff_Methyl)
aff_SNPs_S = calculate_S(aff_SNPs)
aff_final_S = calculate_S(aff_final)

list_aff_S = list(aff_CNV_S, aff_rna_S, aff_miRNA_S, aff_Methyl_S, 
                  aff_SNPs_S, aff_final_S)
names(list_aff_S) = c(paste0("CNV Original Affinity ", optN, "-NN Graph"),
                      paste0("RNAseq Original Affinity ", optN, "-NN Graph"),
                      paste0("miRNA Original Affinity ", optN, "-NN Graph"),
                      paste0("Methylation Original Affinity ", optN, "-NN Graph"),
                      paste0("SNPs Original Affinity ", optN, "-NN Graph"),
                      paste0("Final Fused Affinity ", optN, "-NN Graph"))

for (i in 1:length(list_aff_S)) {
  
  # Prepare the graph object
  g <- graph_from_adjacency_matrix(list_aff_S[[i]],  weighted = TRUE, diag = FALSE,
                                   mode = "max")
  g <- delete_edges(g, E(g)[weight == 0])
  E(g)$width <- sqrt(E(g)$weight) * 5  # Example transformation for visibility
  nodes_data <- data.frame(name = V(g)$name) %>%
    inner_join(clust_annot_pheno %>% dplyr::select(samID, iClusterBayes),
               by = c("name" = "samID"))
  
  # Set iClusterBayes as a factor for coloring
  nodes_data[[algorithm]] <- as.factor(nodes_data[[algorithm]])
  V(g)$iClusterBayes <- nodes_data[[algorithm]] # modify `$iClusterBayes` manually
  
  # Set color based on iClusterBayes
  V(g)$color <- fifelse(V(g)$iClusterBayes == paste0(algorithm, "1"), "#2EC4B6", "#E71D36")
  
  png(paste0(home, 
             "/Results/single_algorithm/", algorithm, "/Supplement/",
             names(list_aff_S)[i], ".png"),
      width = 6000, height = 6000, res = 700)
  
  par(mar = c(2, 2, 2, 5))  # Adjust right margin to accommodate legend
  
  # Plot the graph with a layout that spreads nodes well
  plot(g, vertex.color = V(g)$color,
       edge.width = E(g)$width,
       vertex.size = 4, 
       vertex.label = NA, 
       edge.color = "gray85",
       layout = layout_with_fr(g),  # Use Fruchterman-Reingold layout
       main = "")
  
  # Add title with reduced size using title() function
  title(main = names(list_aff_S)[i], cex.main = 1.7)
  
  # Add a legend to the right of the plot
  legend("bottomright", 
         title="Node Color Legend",    
         legend=c(paste0(algorithm, "1"),
                  paste0(algorithm, "2")), 
         fill=cluster_colors_heatmap,  
         cex=0.7,      
         box.lwd=1)  
  
  dev.off() 
}
rm(g, nodes_data)

# Define the affinity matrices for each modality
# Function to determine the important modalities for each edge
determine_edge_support <- function(edge_weights) {
  sorted_weights <- sort(edge_weights, decreasing = TRUE)
  # If the highest weight is more than 10% greater than all others, it is supported by a single modality
  if (sorted_weights[1] > sorted_weights[2] * 1.1) {
    return(which(edge_weights == sorted_weights[1]))
  }
  # If the difference between the two highest weights is less than 10%, it is supported by those two modalities
  else if (sorted_weights[1] <= sorted_weights[2] * 1.1 && sorted_weights[2] > sorted_weights[3] * 1.1) {
    return(which(edge_weights >= sorted_weights[2]))
  }
  # If the difference between all weights is less than 10%, it is supported by all modalities
  else {
    return(which(edge_weights >= sorted_weights[1] * 0.9))
  }
}

# Prepare the graph object for the final affinity matrix
g <- graph_from_adjacency_matrix(aff_final_S, mode = "undirected", weighted = TRUE, diag = FALSE)
g <- delete_edges(g, E(g)[weight == 0])

# Calculate the weights of each edge in individual networks
edge_weights_list <- list(aff_CNV_S, aff_rna_S, aff_miRNA_S, aff_Methyl_S, aff_SNPs_S)
edge_support_list <- lapply(E(g), function(e) {
  from <- ends(g, e)[1]
  to <- ends(g, e)[2]
  sapply(edge_weights_list, function(mat) mat[from, to])
})

# Assign color to each edge based on the modalities that support it
edge_colors <- c("#FF5733", "#33FF57", "#3357FF", "#FF33A1", "#F3FF33")  # Define unique colors for each modality
combinations_colors <- list()

for (i in 1:length(E(g))) {
  supported_modalities <- determine_edge_support(edge_support_list[[i]])
  if (length(supported_modalities) == 5) {
    combinations_colors[[i]] <- "black"  # Assign black color if all modalities support the edge
  } else if (length(supported_modalities) == 1) {
    combinations_colors[[i]] <- edge_colors[supported_modalities]
  } else {
    combination <- paste(supported_modalities, collapse = " + ")
    if (!is.null(combinations_colors[[combination]])) {
      combinations_colors[[i]] <- combinations_colors[[combination]]
    } else {
      new_color <- colorRampPalette(edge_colors[supported_modalities])(1)
      combinations_colors[[combination]] <- new_color
      combinations_colors[[i]] <- new_color
    }
  }
}

E(g)$color <- unlist(combinations_colors)

# Prepare the node data
nodes_data <- data.frame(name = V(g)$name) %>%
  inner_join(clust_annot_pheno %>% dplyr::select(samID, iClusterBayes), by = c("name" = "samID"))
nodes_data[[algorithm]] <- as.factor(nodes_data[[algorithm]])
V(g)$iClusterBayes <- nodes_data[[algorithm]]

# Set color based on iClusterBayes
V(g)$color <- fifelse(V(g)$iClusterBayes == paste0(algorithm, "1"), "#2EC4B6", "#E71D36")

# Draw the graph for the final affinity matrix
png(paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement/", "colored_Final_Fused_Affinity_", optN, "-NN_Graph.png"),
    width = 6000, height = 6000, res = 700)

par(mar = c(0, 0, 0, 17))  # Adjust right margin to make enough space for the legend

# Plot the graph
plot(g, vertex.color = V(g)$color,
     edge.width = E(g)$width * 0.1,  # Make the edges thinner
     vertex.size = 3,  # Make the nodes smaller
     vertex.label = NA, 
     edge.color = adjustcolor(E(g)$color, alpha.f = 0.35),  # Add transparency to edges
     layout = layout_with_fr(g, niter = 1000))

# Add a legend for the edge colors
legend(x = 1.2, y = 1,  # Manually adjust the position of the legend to the right of the plot
       xpd = TRUE,  # Allow legend to be drawn outside the plot area
       title="Edge Color Legend",
       legend=c("RNA-seq", "miRNA", "Methylation", "CNV", "SNPs", 
                "RNA-seq + miRNA", "RNA-seq + Methylation", "RNA-seq + CNV", "RNA-seq + SNPs", 
                "miRNA + Methylation", "miRNA + CNV", "miRNA + SNPs", 
                "Methylation + CNV", "Methylation + SNPs", 
                "CNV + SNPs", 
                "RNA-seq + miRNA + Methylation", "RNA-seq + miRNA + CNV", "RNA-seq + miRNA + SNPs", 
                "RNA-seq + Methylation + CNV", "RNA-seq + Methylation + SNPs", 
                "RNA-seq + CNV + SNPs", 
                "miRNA + Methylation + CNV", "miRNA + Methylation + SNPs", 
                "miRNA + CNV + SNPs", 
                "Methylation + CNV + SNPs", 
                "RNA-seq + miRNA + Methylation + CNV", "RNA-seq + miRNA + Methylation + SNPs", 
                "RNA-seq + miRNA + CNV + SNPs", "RNA-seq + Methylation + CNV + SNPs", 
                "miRNA + Methylation + CNV + SNPs", 
                "RNA-seq + miRNA + Methylation + CNV + SNPs"),
       fill=c("#33FF57", "#3357FF", "#FF5733", "#FF33A1", "#F3FF33", 
              "#abcdef", "#123456", "#654321", "#fedcba", 
              "#a1b2c3", "#b2c3d4", "#c3d4e5", 
              "#d4e5f6", "#e5f6a7", "#f6a7b8", 
              "#aabbcc", "#bbccdd", "#ccddee", 
              "#ddeeff", "#eeffaa", "#ffaabb", 
              "#aaffcc", "#bbffdd", "#ccffee", 
              "#ffccaa", 
              "#112233", "#223344", "#334455", 
              "#445566", "#556677", "#667788", 
              "black"),
       cex=0.7,
       box.lwd=1)

dev.off()

# Compare these iClusterBayes results with the iClusterBayes output from MOVICS ###
load("Results/MOVICS_baseline/MOVICS_TCGA_RNAseq-CNV-Methylation-miRNA-SNPs_eval_on_transNEO_moic.res.list.rda")
MOVICS_iClusterBayes = moic.res.list$iClusterBayes$clust.res

ARI_to_MOVICS_iClusterBayes = calculate_ari_index(cluster_df1 = MOVICS_iClusterBayes %>%
                                          dplyr::rename(Sample.ID = samID,
                                                        Cluster = clust),
                                        cluster_df2 = iClusterBayes_clusters,
                                        sample_col = "Sample.ID",
                                        clust_col = "Cluster",
                                        suffixes = c(paste0("_MOVICS_", algorithm),
                                                     paste0("_", algorithm)))

NMI_to_MOVICS_iClusterBayes = calculate_nmi_index(cluster_df1 = MOVICS_iClusterBayes %>%
                                          dplyr::rename(Sample.ID = samID,
                                                        Cluster = clust),
                                        cluster_df2 = iClusterBayes_clusters,
                                        sample_col = "Sample.ID",
                                        clust_col = "Cluster",
                                        suffixes = c(paste0("_MOVICS_", algorithm),
                                                     paste0("_", algorithm)))

# Wrap up #####
hyperparameters = list(num_neighbors_min = min(num_neighbors_range),
                       num_neighbors_max = max(num_neighbors_range),
                       num_neighbors_step = neighbor_step,
                       sigma_min = min(sigma_range),
                       sigma_max = max(sigma_range),
                       sigma_step = sigma_step,
                       optimal_N = optN,
                       optimal_sigma = optSigma,
                       conclusion = conclusion, # if there is agreement, np_conclusion can also be used
                       n_iter = n_iterations,
                       feature_ranks_text = feature_ranks_text
)

# Put all parameters in a list
params = list(algorithm = algorithm, data_source = data_source, data_types = data_types,
              evaluation_source = evaluation_source, title = title, subtitle = subtitle,
              description = description, in_a_nutshell = in_a_nutshell, optk_text = optk_text,
              citation = citation, NMI_to_MOVICS = NMI_to_MOVICS, ARI_to_MOVICS = ARI_to_MOVICS,
              NMI_to_MOVICS_iClusterBayes = NMI_to_MOVICS_iClusterBayes, ARI_to_MOVICS_iClusterBayes = ARI_to_MOVICS_iClusterBayes,
              hyperparameters = hyperparameters, ground_truth_k = ground_truth_k,
              transNEO_var2comp = transNEO_var2comp,
              sessionInfo = sessionInfo(), home = home)

# Render the R Markdown document with the parameters
rmarkdown::render(paste0(getwd(), "/Results/single_algorithm/", algorithm,
                         "/", algorithm, "_report.Rmd"), 
                  params = params, 
                  output_file = paste0(home, "/Results/single_algorithm/", 
                                       algorithm, "/", algorithm, "_report_",
                                       data_source, "_",
                                       data_types, "_eval_on_", evaluation_source,
                                       ".html"))

# Export session info as .txt
writeLines(capture.output(sessionInfo()), paste0("sessionInfo/",
                                                 algorithm, "_", data_source, "_",
                                                 data_types, "_eval_on_", evaluation_source,
                                                 "_sessionInfo.txt"))

# Save environment
save.image(paste0(home, "/Results/single_algorithm/", 
                  algorithm, "/", algorithm, "_", data_source, "_",
                  data_types, "_eval_on_", evaluation_source,
                  "_env.RData"))
