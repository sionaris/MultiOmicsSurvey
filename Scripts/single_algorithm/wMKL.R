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

# Installation of the package is a bit tricky
# 1. Use setwd() and navigate to a folder other than your R libraries directory

# 2. git clone https://github.com/biostatcao/wMKL.git

# 3. I renamed the downloaded directory to wMKL_source (that might be optional)

# 4. Replace lines 33-46 in src/tsne.cpp with:

# extern "C" {
#   #include <R_ext/BLAS.h>
# }
# 
# #include <Rcpp.h>
# #include <math.h>
# #include <float.h>
# #include <stdlib.h>
# #include <stdio.h>
# #include <cstring>
# #include <time.h>
# #include "sptree.h"
# #include "vptree.h"
# #include "tsne.h"

# 5. Go to lines 867-885 in src/tsne.cpp and replace them with:

# // Compute squared Euclidean distance matrix (using BLAS)
# void TSNE::computeSquaredEuclideanDistance(double* X, int N, int D, double* DD) {
#   double* dataSums = (double*) calloc(N, sizeof(double));
#   if(dataSums == NULL) { Rcpp::stop("Memory allocation failed!\n"); }
#   for(int n = 0; n < N; n++) {
#     for(int d = 0; d < D; d++) {
#       dataSums[n] += (X[n * D + d] * X[n * D + d]);
#     }
#   }
#   for(int n = 0; n < N; n++) {
#     for(int m = 0; m < N; m++) {
#       DD[n * N + m] = dataSums[n] + dataSums[m];
#     }
#   }
#   double a1 = -2.0, a2 = 1.0;
#   char transT = 'T';
#   char transN = 'N';
#   
#   F77_CALL(dgemm)(
#     &transT, &transN,    // "T", "N"
#     &N, &N, &D,          // m, n, k
#     &a1, X, &D,          // alpha, A, lda
#     X, &D,          // B, ldb
#     &a2, DD, &N          // beta, C, ldc
#     FCONE FCONE          // the two hidden lengths
#   );
#   free(dataSums); dataSums = NULL;
# }

# 6. Install the package by running:
# install.packages("wMKL_source", repos = NULL, type = "source")
library(wMKL)

# Preamble
home = getwd()
algorithm = "wMKL"
alg_feature_pref = "rows" # Where does the algorithm expect the features to be
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

# However wMKL prefers features in columns so we transpose the matrices.

# Extract the names of the modalities that will be used
modalities = unlist(strsplit(data_types, "-"))

# Replace with TRUE wherever features are in rows
features_in_rows = rep(TRUE, length(modalities))

# Run algorithm #####

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
COSMIC_BC_drivers = read.csv("Resources/COSMIC_CGC_Breast_somatic.csv")$Gene.Symbol %>%
  as.character()

# Weights for features using MAD (similar to ab-SNF)
# Returns a vector with dimensions equal to the number of features of each modality
# of normalized variance-based weights for the columns.
mad_based_weights <- function(X) {
  v <- apply(t(X), 2, mad, na.rm = TRUE) # MAD of each feature (features in columns)
  v[!is.finite(v)] <- 0
  if (sum(v) == 0) {
    # all features are constant, fallback: uniform weighting
    w <- rep(1, length(v))
    w <- w / sum(w) # normalize so sum of all weights is 1
    names(w) = rownames(X)
  } else {
    w <- v / sum(v) # normalize so sum of all weights is 1
    names(w) = rownames(X)
  }
  return(w)
}

# For binary data, we weigh 0.8/0.2 drivers vs. non-drivers
weighted_mutations <- function(features, drivers,
                               driver_weight, nondriver_weight) {
  w <- rep(NA, length(features))
  w[which(features %in% drivers)] = driver_weight
  w[which(!features %in% drivers)] = nondriver_weight
  w <- w / sum(w) # normalize so sum of all weights is 1
  return(w)
}

# Lists of weights
# Modality types
continuous = c("RNAseq", "CNV", "Methylation", "miRNA")
categorical = c("SNPs")

# Calculate the pair-wise distance (Euclidean for continuous modalities)
cont_weights = lapply(input[continuous], mad_based_weights)

bin_weights = weighted_mutations(features = rownames(input[["SNPs"]]),
                                 drivers = COSMIC_BC_drivers,
                                 driver_weight = 0.8,
                                 nondriver_weight = 0.2)
names(bin_weights) = rownames(input[["SNPs"]])

# Final list of weights
weights = cont_weights
weights[["SNPs"]] = bin_weights
weights = weights[names(input)]

# Distance methods
distance_methods = c("binary", rep("sqeuclidean", 4)) # square Euclidean is the default

# Optimal number of clusters
RNGversion("4.2.2")
set.seed(123)
num_results = CIMLR_Estimate_Number_of_Clusters_weight_mod(input,
                                                     NUMC = 2:10,
                                                     cores.ratio = 0,
                                                     weight = weights,
                                                     methods = distance_methods)

optk = which.min(num_results$K1) + 1

# perform the CIMLR.weight clustering algorithm
t1 = Sys.time()
cluster_results = CIMLR.weight_mod(X = input, c = optk,
                             cores.ratio = 0,
                             weight = weights,
                             methods = distance_methods)
dt = Sys.time() - t1 # ~4.5 mins

# Computing the multiple Kernels.
# Performing network diffusion.
# Iteration:  1 
# Iteration:  2 
# Iteration:  3 
# Iteration:  4 
# Iteration:  5 
# Iteration:  6 
# Iteration:  7 
# Iteration:  8 
# Iteration:  9 
# Iteration:  10 
# Performing t-SNE.
# Epoch: Iteration # 100  error is:  0.1322815 
# Epoch: Iteration # 200  error is:  0.1085342 
# Epoch: Iteration # 300  error is:  0.09815787 
# Epoch: Iteration # 400  error is:  0.09284991 
# Epoch: Iteration # 500  error is:  0.0898325 
# Epoch: Iteration # 600  error is:  0.08782597 
# Epoch: Iteration # 700  error is:  0.08635315 
# Epoch: Iteration # 800  error is:  0.08519087 
# Epoch: Iteration # 900  error is:  0.0842451 
# Epoch: Iteration # 1000  error is:  0.08344559 
# Performing Kmeans.
# Performing t-SNE.
# Epoch: Iteration # 100  error is:  11.36326 
# Epoch: Iteration # 200  error is:  0.2706709 
# Epoch: Iteration # 300  error is:  0.2124627 
# Epoch: Iteration # 400  error is:  0.1993325 
# Epoch: Iteration # 500  error is:  0.1940114 
# Epoch: Iteration # 600  error is:  0.1905818 
# Epoch: Iteration # 700  error is:  0.1881677 
# Epoch: Iteration # 800  error is:  0.1863444 
# Epoch: Iteration # 900  error is:  0.1849099 
# Epoch: Iteration # 1000  error is:  0.1837464 

# The function returns a list of objects
# y / y_spectral: Two different final cluster assignments of the data. y is kmeans, y_spectral is spectral
# S: Learned similarity matrix from all kernels. Final affinity matrix
# F: Final t-SNE embedding in user-specified dimension(s).
# ydata: Another t-SNE embedding, typically 2D for plotting.
# alphaK: Learned weights that combine the multiple kernel matrices.
# execution.time: Runtime measurement.
# converge: Convergence metric across iterations.
# LF: Final eigenvector-based embedding from the Laplacian.

# Sum of kernel weights (56 kernels per modality)
breaks = seq(0, length(cluster_results$alphaK), 55) # 55 kernels per modality
modality_specific_weights = list(SNPs = sum(cluster_results$alphaK[1:breaks[2]]),
                                 RNAseq = sum(cluster_results$alphaK[(breaks[2]+1):(breaks[3])]),
                                 CNV = sum(cluster_results$alphaK[(breaks[3]+1):(breaks[4])]),
                                 Methylation = sum(cluster_results$alphaK[(breaks[4]+1):(breaks[5])]),
                                 miRNA = sum(cluster_results$alphaK[(breaks[5]+1):(breaks[6])]))

# Our clustering assignments are the y_spectral in this case
wMKL_clusters = as.data.frame(list(Sample.ID = colnames(input$SNPs),
                                   Cluster = cluster_results$y_spectral))
wMKL_clusters$Sample.ID = gsub("\\.", "-", wMKL_clusters$Sample.ID)
rownames(wMKL_clusters) = wMKL_clusters$Sample.ID

# Main results #####
# Examine cluster similarity to MOVICS by measuring NMI and ARI indices #####
# (Jaccard may be misleading)

# Calculate ARI and NMI
library(mclust)
library(clue)

ARI_to_MOVICS = calculate_ari_index(cluster_df1 = ground_truth_labels,
                                    cluster_df2 = wMKL_clusters,
                                    sample_col = "Sample.ID",
                                    clust_col = "Cluster",
                                    suffixes = c("_MOVICS_CS",
                                                 paste0("_", algorithm)))

NMI_to_MOVICS = calculate_nmi_index(cluster_df1 = ground_truth_labels,
                                    cluster_df2 = wMKL_clusters,
                                    sample_col = "Sample.ID",
                                    clust_col = "Cluster",
                                    suffixes = c("_MOVICS_CS",
                                                 paste0("_", algorithm)))

# Low to moderate statistics; slightly similar results; NMI might be biased due to a large
# number of clusters in wMKL

# Import coloring scheme
scheme = readRDS("Resources/scheme.rds")
annCol = scheme$annCol
annColors = scheme$annColors
cluster_colors = c("#2EC4B6", "#E71D36", "#FF9F1C", "#BDD5EA", 
                   "#FFA5AB", "#011627", "#023E8A", "#9D4EDD")
col.list = scheme$col.list
var2comp = scheme$var2comp %>%
  dplyr::select(-`Consensus Subtype`) %>%
  mutate(Sample.ID = rownames(.)) %>%
  inner_join(wMKL_clusters, by = "Sample.ID") %>%
  tibble::column_to_rownames(var = "Sample.ID") %>%
  mutate(wMKL = paste0(algorithm, Cluster)) %>%
  dplyr::select(wMKL, everything()) %>%
  dplyr::select(-Cluster)
rm(scheme); gc()

# Silhouette
library(MOVICS)
library(cluster)

# Extract the eign-space in which spectral clustering takes place under the hood
# in wMKL:::spectralClustering
extract_eigenspace = function (affinity, K, type = 3) 
{
  d <- rowSums(affinity)
  d[d == 0] <- .Machine$double.eps
  D <- diag(d)
  L <- D - affinity
  if (type == 1) {
    NL <- L
  }
  else if (type == 2) {
    Di <- diag(1/d)
    NL <- Di %*% L
  }
  else if (type == 3) {
    Di <- diag(1/sqrt(d))
    NL <- Di %*% L %*% Di
  }
  eig <- eigen(NL)
  res <- sort(abs(eig$values), index.return = TRUE)
  U <- eig$vectors[, res$ix[1:K]]
  normalize <- function(x) x/sqrt(sum(x^2))
  if (type == 3) {
    U <- t(apply(U, 1, normalize))
  }
  eigDiscrete <- wMKL:::.discretisation(U)
  eigDiscrete <- eigDiscrete$discrete
  labels <- apply(eigDiscrete, 1, which.max)
  return(list(labels = labels, U = U, eigDiscrete = eigDiscrete))
}

U = extract_eigenspace(affinity = cluster_results$S, K = optk)$U
rownames(U) = colnames(input$SNPs)
colnames(U) = paste0("eig", seq(1, optk, 1))
U = as.data.frame(U)

silhouette = silhouette(as.integer(wMKL_clusters$Cluster),
                        dist = Rfast::Dist(U,
                                           method = "euclidean"))

getSilhouette_ggplot(sil      = silhouette,
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

# MOVICS-like analysis #####
library(ComplexHeatmap)

plotdata <- lapply(lapply(input, as.matrix), 
                   function(mat) mat[rowSums(mat != 0) > 0, ])
plotdata = getStdiz(
  data = plotdata,
  halfwidth = c(NA, 3, 3, 3, 3), # No halfwidth for SNPs
  centerFlag = c(F, F, F, F, F),
  scaleFlag = c(F, F, F, F, F)
)

plot_object = list(clust.res = wMKL_clusters %>%
                     dplyr::rename(samID = Sample.ID, clust = Cluster))

# Export consensus clustering object
clust = as.data.frame(plot_object$clust.res)
colnames(clust) = c("Sample.ID", "Cluster")
clust$Cluster = paste0(algorithm, clust$Cluster)
openxlsx::write.xlsx(clust, paste0(home, "/Results/single_algorithm/", algorithm, "/",
                                   algorithm, "_", data_source, "_",
                                   data_types, "_eval_on_", evaluation_source,
                                   "_clusterings.xlsx"))

# Order features
feature_orders = readRDS("Resources/TCGA/mm_feature_orders.rds")
for (i in 1:length(plotdata)) {
  plotdata[[i]] = plotdata[[i]][feature_orders[[names(plotdata)[i]]], , drop = FALSE]
}

getMoHeatmap_single_algorithm2(algorithm_name = algorithm,
                               data          = plotdata,
                               row.title     = names(plotdata),
                               is.binary     = c(T,F,F,F,F), 
                               legend.name   = c("SNPs",
                                                 "Standardized RNAseq norm. counts",
                                                 "Standardized CNV",
                                                 "Standardized miRNA norm. counts",
                                                 "Standardized Methylation M-values"
                               ),
                               cluster_rows = rep(F, length(plotdata)),
                               cluster_cols = rep(F, length(plotdata)),
                               show.col.dend = rep(F, length(plotdata)),
                               show.colnames = FALSE,
                               show.row.dend = rep(F, length(plotdata)),
                               show.rownames = rep(F, length(plotdata)),
                               clust.res     = plot_object$clust.res, # consensusMOIC-like results
                               # clust.dist.row = c("manhattan", rep("euclidean", 4)),
                               # clust.method.row = rep("ward.D", length(plotdata)),
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
                                         pdf_level_col_width = c("3em", "3em"),
                                         pdf_count_col_width = "5em",
                                         pdf_pval_col_width = "3em",
                                         pdf_test_col_width = "3em",
                                         pdf_tab_font_size = 7)

clin_ordinal_comp = compClinvar_ordinal_single_algorithm(algorithm_name = algorithm,
                                                         moic.res = plot_object,
                                                         var2comp = var2comp_nonas %>%
                                                           dplyr::select(number_of_lymphnodes_positive_by_ihc,
                                                                         number_of_lymphnodes_positive_by_he,
                                                                         wMKL),
                                                         strata = algorithm,
                                                         ordinalVars = c("number_of_lymphnodes_positive_by_ihc",
                                                                         "number_of_lymphnodes_positive_by_he"),
                                                         includeNA = FALSE,
                                                         tab.name = "Summary of ordinal clinical variables",
                                                         res.path = paste0(home, "/Results/single_algorithm/", algorithm, "/"),
                                                         output_pdf = TRUE,
                                                         pdf_template_loc = paste0(home, "/Scripts/automated_scripts/clincomp_template.Rmd"),
                                                         pdf_level_col_width = c("3em", "3em"),
                                                         pdf_count_col_width = "5em",
                                                         pdf_pval_col_width = "3em",
                                                         pdf_test_col_width = "3em",
                                                         pdf_tab_font_size = 7)

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
                                      return.binary = TRUE,
                                      fig.name     = paste0(algorithm, "_", data_source, "_",
                                                            data_types, "_eval_on_", evaluation_source,
                                                            "_oncoprint"),
                                      tab.name     = "Independent test between subtype and mutation",
                                      fig.path     = paste0(home, "/Results/single_algorithm/", algorithm),
                                      res.path     = paste0(home, "/Results/single_algorithm/", algorithm),
                                      test.method = "chisq" # large number of groups
                                      )

# Inspect clusters of interest
# wMKL6 and wMKL7 mutations vs HER2 status
HER2_onco = oncoprint$sample_binary %>% inner_join(clinical_data %>%
                                                     dplyr::select(Sample.ID, HER2_status = lab_proc_her2_neu_immunohistochemistry_receptor_status), by = "Sample.ID") %>%
  inner_join(wMKL_clusters, by = "Sample.ID")
HER2_onco_wMKL6 = HER2_onco %>%
  dplyr::filter(Cluster == "6") %>%
  dplyr::select(-Cluster)
HER2_onco_wMKL7 = HER2_onco %>%
  dplyr::filter(Cluster == "7") %>%
  dplyr::select(-Cluster)

# Proportion tables

# TP53
prop.table(table(HER2_onco_wMKL6$HER2_status, HER2_onco_wMKL6$TP53), margin = 1)
prop.table(table(HER2_onco_wMKL7$HER2_status, HER2_onco_wMKL7$TP53), margin = 1)

# PIK3CA
prop.table(table(HER2_onco_wMKL6$HER2_status, HER2_onco_wMKL6$PIK3CA), margin = 1)
prop.table(table(HER2_onco_wMKL7$HER2_status, HER2_onco_wMKL7$PIK3CA), margin = 1)

# CDH1
prop.table(table(HER2_onco_wMKL6$HER2_status, HER2_onco_wMKL6$CDH1), margin = 1)
prop.table(table(HER2_onco_wMKL7$HER2_status, HER2_onco_wMKL7$CDH1), margin = 1)

# SYNE1
prop.table(table(HER2_onco_wMKL6$HER2_status, HER2_onco_wMKL6$SYNE1), margin = 1)
prop.table(table(HER2_onco_wMKL7$HER2_status, HER2_onco_wMKL7$SYNE1), margin = 1)

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
                                                 width = 10,
                                                 height = 10,
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
                                             show_rownames = FALSE, # show no rownames (biomarker name)
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
                                               show_rownames = FALSE, # show no rownames (biomarker name)
                                               centerFlag = F,
                                               scaleFlag = F,
                                               halfwidth = 3,
                                               fig.name      = "downregulated_biomarkers_heatmap",
                                               fig.path = paste0(home, "/Results/single_algorithm/", algorithm),
                                               width = 14,
                                               height = 18,
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
                                               show_rownames = FALSE, # show no rownames (biomarker name)
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
                                                 show_rownames = FALSE, # show no rownames (biomarker name)
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
                                              show_rownames = FALSE, # show no rownames (biomarker name)
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
                                                show_rownames = FALSE, # show no rownames (biomarker name)
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
                                            n.path       = 10,
                                            p.cutoff     = 0.05, # p cutoff to identify significant pathways
                                            p.adj.cutoff = 0.05, # padj cutoff to identify significant pathways
                                            gsva.method  = "gsva", # method to calculate single sample enrichment score
                                            name         = "GSVA scores", # name for colorbar
                                            norm.method  = "mean", # normalization method to calculate subtype-specific enrichment score
                                            fig.name     = "upregulated_pathway_heatmap",
                                            nPerm = 10000,
                                            minGSSize = 10,
                                            maxGSSize = 500,
                                            fig.path = paste0(home, "/Results/single_algorithm/", algorithm),
                                            width = 14, height = 20)

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
                                              n.path       = 10,
                                              p.cutoff     = 0.05, # p cutoff to identify significant pathways
                                              p.adj.cutoff = 0.05, # padj cutoff to identify significant pathways
                                              gsva.method  = "gsva", # method to calculate single sample enrichment score
                                              name         = "GSVA scores", # name for colorbar
                                              norm.method  = "mean", # normalization method to calculate subtype-specific enrichment score
                                              fig.name     = "downregulated_pathway_heatmap",
                                              nPerm = 10000,
                                              minGSSize = 10,
                                              maxGSSize = 500,
                                              fig.path = paste0(home, "/Results/single_algorithm/", algorithm),
                                              width = 14, height = 20)

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

timestamp() # ~2 mins
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
                                                present_clusters = c(paste0("wMKL", seq(1, optk, 1))),
                                                representative = TRUE, moic.res = plot_object,
                                                subtype_prefix = algorithm, n.path = 10, msigdb.path = MSIGDB.FILE,
                                                norm.method = "mean", dirct = "up",
                                                fig.name = "upregulated_pathway_heatmap",
                                                name = "GSVA scores",
                                                fig.path = paste0(home, "/Results/single_algorithm/", algorithm), 
                                                width = 15, height = 12, gsva.method = "gsva")

hclust_pathway_plots_down = plot_pathway_heatmaps(gsea.lists = hclust_output[grepl("down", names(hclust_output))], 
                                                  norm.expr = plotdata$RNAseq, 
                                                  present_clusters = c(paste0("wMKL", seq(1, optk, 1))),
                                                  representative = TRUE, moic.res = plot_object,
                                                  subtype_prefix = algorithm, n.path = 10, msigdb.path = MSIGDB.FILE,
                                                  norm.method = "mean", dirct = "down",
                                                  fig.name = "downregulated_pathway_heatmap",
                                                  name = "GSVA scores",
                                                  fig.path = paste0(home, "/Results/single_algorithm/", algorithm), 
                                                  width = 15, height = 12, gsva.method = "gsva")

# Fraction Genome Altered ###
fga_df = readRDS("Resources/TCGA/fga_df.rds"); gc()

fga.wMKL <- compFGA_optimized(moic.res     = plot_object,
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

fga.wMKL.COSMIC <- compFGA_optimized(moic.res     = plot_object,
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
timestamp() # 10 min

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
timestamp() # 11 min

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
  inner_join(expr_conc %>% dplyr::select(Donor.ID = samID, wMKL = clust_up),
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
                                                 pdf_level_col_width = c("3em", "3em"),
                                                 pdf_count_col_width = "5em",
                                                 pdf_pval_col_width = "3em",
                                                 pdf_test_col_width = "3em",
                                                 pdf_tab_font_size = 7)

transNEO_ntp_expr_up_ord = transNEO_ntp_expr_up
transNEO_ntp_expr_up_ord$clust.res$clust = gsub(algorithm, "", transNEO_ntp_expr_up_ord$clust.res$clust)
transNEO_ordinal_clincomp = compClinvar_ordinal_single_algorithm(algorithm_name = algorithm,
                                                                 moic.res = transNEO_ntp_expr_up_ord,
                                                                 var2comp = transNEO_var2comp_nonas %>%
                                                                   dplyr::select(Grade.pre.NAT, 
                                                                                 Chemo.cycles, 
                                                                                 aHER2.cycles,
                                                                                 wMKL),
                                                                 strata = algorithm,
                                                                 ordinalVars = c("Grade.pre.NAT",
                                                                                 "Chemo.cycles",
                                                                                 "aHER2.cycles"),
                                                                 includeNA = FALSE,
                                                                 tab.name = "transNEO Summary of ordinal clinical variables",
                                                                 res.path = paste0(home, "/Results/single_algorithm/", algorithm, "/"),
                                                                 output_pdf = TRUE,
                                                                 pdf_template_loc = paste0(home, "/Scripts/automated_scripts/clincomp_template.Rmd"),
                                                                 pdf_level_col_width = c("3em", "3em"),
                                                                 pdf_count_col_width = "5em",
                                                                 pdf_pval_col_width = "3em",
                                                                 pdf_test_col_width = "3em",
                                                                 pdf_tab_font_size = 7)

# Run PAM ###
RNGversion("4.2.2")
set.seed(123)
transNEO_pam = runPAM_single_algorithm(algorithm_name = algorithm,
                                       train.expr = plotdata$RNAseq,
                                       moic.res   = plot_object,
                                       test.expr  = transcr)

# Check consistency across methods

# Get predictions for TCGA (discovery cohort)
RNGversion("4.2.2")
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

# NTP transNEO vs PAM transNEO - CANNOT be produced because there is a subtype present only in one classification
# runKappa_single_algorithm(algorithm_name = algorithm,
#                           subt1 = as.numeric(gsub(algorithm, "",
#                                                   transNEO_ntp_expr_up$clust.res$clust)),
#                           subt2 = as.numeric(transNEO_pam$clust.res$clust),
#                           subt1.lab = "transNEO NTP",
#                           subt2.lab = "transNEO PAM",
#                           height = 8,
#                           width = 8,
#                           fig.path = paste0(home, "/Results/single_algorithm/", algorithm),
#                           fig.name = "kappa_NTP_vs_PAM_transNEO")

# Supplementary results #####

# Create subdirectory for supplementary plots
if (!dir.exists(paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement"))) {
  dir.create(paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement"))
}

# Setup for heatmaps
colors_heatmap = rev(colorRampPalette(viridisLite::magma(10))(255))
cluster_colors_heatmap = cluster_colors
clust_annot_pheno = annCol %>% mutate(Sample.ID = rownames(.)) %>%
  inner_join(clust, by = "Sample.ID") %>%
  dplyr::rename(wMKL = Cluster, samID = "Sample.ID")
rownames(clust_annot_pheno) = clust_annot_pheno$samID
afh_colnames = colnames(annCol)

# Same data frame. Different columns. Just for easiness
wMKL_clust_res = wMKL_clusters %>% dplyr::rename(samID = Sample.ID, 
                                                 wMKL = Cluster) %>%
  dplyr::mutate(wMKL = gsub(algorithm, "", wMKL))

# PCA from final similarity matrix S
final_affinity_matrix = cluster_results$S
dimnames(final_affinity_matrix) = list(colnames(input$SNPs), colnames(input$SNPs))
pca_from_sim_matrix(sim_matrix = final_affinity_matrix, algorithm = algorithm, 
                    clust_res = clust_annot_pheno %>% dplyr::select(samID, wMKL),
                    cluster_colors = cluster_colors_heatmap, 
                    output_path = paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement"), 
                    title_add = "Final similarity matrix S")

# PCA from original matrices ###
# t-SNE projection
tsne_2D = cluster_results$ydata
colnames(tsne_2D) = c("t-SNE1", "t-SNE2")
rownames(tsne_2D) = colnames(input$SNPs)
tsne_2D = as.data.frame(tsne_2D) %>%
  mutate(Sample.ID = rownames(.)) %>%
  inner_join(wMKL_clusters, by = "Sample.ID") %>%
  dplyr::rename(wMKL = Cluster) %>%
  dplyr::mutate(wMKL = paste0(algorithm, wMKL))

ggplot(data = tsne_2D, aes(x = `t-SNE1`, y = `t-SNE2`, color = wMKL)) +
  geom_point() +
  aes(shape = as.character(wMKL), 
      size = as.character(wMKL),
      alpha = as.character(wMKL), 
      color = as.character(wMKL)) +
  scale_size_manual(name = algorithm, 
                    values = rep(0.4, length(unique(tsne_2D$wMKL))), 
                    labels = sort(unique(tsne_2D$wMKL))) +
  scale_alpha_manual(name = algorithm, 
                     values = rep(0.7, length(unique(tsne_2D$wMKL))), 
                     labels = sort(unique(tsne_2D$wMKL))) +
  scale_shape_manual(name = algorithm, 
                     values = rep(16, length(unique(tsne_2D$wMKL))), 
                     labels = sort(unique(tsne_2D$wMKL)))+
  scale_color_manual(name = algorithm, 
                     limits = sort(unique(tsne_2D$wMKL)),
                     values = setNames(cluster_colors, sort(unique(tsne_2D$wMKL))), 
                     labels = sort(unique(tsne_2D$wMKL))) +
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(linewidth = 0.2),
        plot.title = element_text(size = 5, face = "bold"),
        axis.title = element_text(size = 4, face = "bold"),
        axis.text = element_text(size = 4),
        axis.ticks = element_line(linewidth = 0.15),
        legend.background = element_rect(fill = "white", linetype = "solid"),
        #legend.position = c(0.90, 0.86),
        legend.key.size = unit(0.5, "lines"),
        legend.margin = ggplot2::margin(0, 0, 0, 0, unit = "mm"),
        legend.spacing.y = unit(0.1, units = "cm"), 
        legend.title = ggplot2::element_text(size = 4, face = "bold"), 
        legend.text = ggplot2::element_text(size = 3))+
  labs(title = paste0(algorithm, " clusters in 2D t-SNE projection"),
       legend = algorithm) +
  guides(size = "none", alpha = "none", shape = "none")
ggsave(filename = paste0(algorithm, "_t-SNE_2D_projection.png"),
       path = paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement"), 
       width = 1920, height = 1080, device = 'png', units = "px",
       dpi = 700)
dev.off()

# RNA
pca_from_original_matrix(mydata = plotdata$RNAseq, 
                         algorithm = algorithm, 
                         clust_res = wMKL_clust_res,
                         cluster_colors = cluster_colors, 
                         output_path = paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement"),
                         title_add = "RNAseq")

# miRNA
pca_from_original_matrix(mydata = plotdata$miRNA, 
                         algorithm = algorithm, 
                         clust_res = wMKL_clust_res,
                         cluster_colors = cluster_colors, 
                         output_path = paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement"),
                         title_add = "miRNA")

# CNV
pca_from_original_matrix(mydata = plotdata$CNV, 
                         algorithm = algorithm, 
                         clust_res = wMKL_clust_res,
                         cluster_colors = cluster_colors, 
                         output_path = paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement"),
                         title_add = "CNV")

# Use multidimensional scaling for SNPs
# Features must be in rows
mds_from_original_matrix(matrix = plotdata$SNPs, dist_method = "binary",
                         algorithm = algorithm, 
                         clust_res = wMKL_clust_res,
                         cluster_colors = cluster_colors, 
                         output_path = paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement"),
                         title_add = "SNPs")

# Methylation
pca_from_original_matrix(mydata = plotdata$Methylation, 
                         algorithm = algorithm, 
                         clust_res = wMKL_clust_res,
                         cluster_colors = cluster_colors, 
                         output_path = paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement"),
                         title_add = "Methylation")

# Draw a heatmap of the final S matrix ###
create_MO_heatmap(matrix = final_affinity_matrix, algorithm = algorithm, 
                  need.diag.zero = FALSE, # already zero
                  clust_annot_pheno = clust_annot_pheno ,
                  afh_colnames = afh_colnames, 
                  colors = colors_heatmap,
                  annColors = annColors,
                  heatmap_title = "Final wMKL similarity (S matrix) heatmap",
                  cluster_colors = cluster_colors_heatmap,
                  cluster_rows_flag = FALSE,
                  cluster_cols_flag = FALSE,
                  splits_flag = TRUE,
                  legend_title = "Final wMKL similarity",
                  output_file_name = paste0(home, "/Results/single_algorithm/", algorithm, 
                                            "/Supplement/wMKL_final_S_matrix_heatmap.png"))

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
wMKL_barcharts = list()
plotdata_bar = clust_annot_pheno_nonas
plotdata_bar[[algorithm]] = factor(plotdata_bar[[algorithm]])
for (i in 1:length(voi)) {
  chifit = output
  loc = which(grepl(voi[i], chifit$Comparison))
  chifit = chifit[loc, ]
  wMKL_barcharts[[i]] = create_annot_barchart(plotdata = plotdata_bar, fill = voi[i],
                                              chifit = chifit,
                                              na.action = "na.omit",
                                              algorithm = algorithm,
                                              barchart_ylim = 300,
                                              text_y = 250, rect_ymin = 165,
                                              rect_ymax = 265, x_annot = 4.5,
                                              v_gap = 35, rect_xmin = 3.5,
                                              rect_xmax = 5.5, 
                                              annot_text_size = 2.25,
                                              legend.text.size = 5,
                                              x.axis.text.size = 5) +
    barchart_scales[[voi[i]]]
  print(wMKL_barcharts[[i]])
  ggsave(filename = paste0(algorithm, "_", voi[i], "_barchart.png"),
         path = paste0(home, 
                       "/Results/single_algorithm/", algorithm, "/Supplement"), 
         width = 3320, height = 2320, device = 'png', units = "px",
         dpi = 700)
  dev.off()
}
names(wMKL_barcharts) = voi
rm(loc, chifit)

# Multiplot (PNG) - bar charts
library(ggpubr)
ggarrange(wMKL_barcharts[[1]], wMKL_barcharts[[2]], wMKL_barcharts[[3]],
          wMKL_barcharts[[4]], wMKL_barcharts[[5]], wMKL_barcharts[[6]],
          wMKL_barcharts[[7]], wMKL_barcharts[[8]], wMKL_barcharts[[9]],
          wMKL_barcharts[[10]], wMKL_barcharts[[11]],
          ncol = 3, nrow = 4, labels = c("A", "B", "C", "D", "E", "F", "G", "H",
                                         "I", "J", "K"),
          font.label = list(size = 8, face = "bold", color ="black"))
ggsave(filename = paste0("Multiplot_", algorithm, "_barcharts.png"),
       path = paste0(home, 
                     "/Results/single_algorithm/", algorithm, "/Supplement"), 
       width = 10000, height = 8000, device = 'png', units = "px",
       dpi = 700)
dev.off()

# Just significant ones now
wMKL_barcharts_sig = list()
plotdata_bar_sig = clust_annot_pheno_nonas %>% dplyr::select(wMKL, Race, Histology, 
                                                             `ER status`, `PR status`, `HER2 status`,
                                                             `Metastasis`)
plotdata_bar_sig$wMKL = factor(plotdata_bar_sig$wMKL)
voi_sig = setdiff(colnames(plotdata_bar_sig), algorithm)
for (i in 1:length(voi_sig)) {
  chifit = output
  loc = which(grepl(voi_sig[i], chifit$Comparison))
  chifit = chifit[loc, ]
  wMKL_barcharts_sig[[i]] = create_annot_barchart(plotdata = plotdata_bar_sig, fill = voi_sig[i],
                                                  chifit = chifit,
                                                  na.action = "na.omit",
                                                  algorithm = algorithm,
                                                  barchart_ylim = 300,
                                                  text_y = 250, rect_ymin = 165,
                                                  rect_ymax = 265, x_annot = 4.5,
                                                  v_gap = 35, rect_xmin = 3.5,
                                                  rect_xmax = 5.5, 
                                                  annot_text_size = 2.25,
                                                  legend.text.size = 5,
                                                  x.axis.text.size = 5) +
    barchart_scales[[voi_sig[i]]]
  print(wMKL_barcharts_sig[[i]])
  ggsave(filename = paste0("sig_", algorithm, "_", voi_sig[i], "_barchart.png"),
         path = paste0(home, 
                       "/Results/single_algorithm/", algorithm, "/Supplement"), 
         width = 3320, height = 2320, device = 'png', units = "px",
         dpi = 700)
  dev.off()
}
names(wMKL_barcharts_sig) = voi_sig
rm(loc, chifit)

# Multiplot (PNG) - bar charts
ggarrange(wMKL_barcharts_sig[[1]], wMKL_barcharts_sig[[2]], wMKL_barcharts_sig[[3]],
          wMKL_barcharts_sig[[4]], wMKL_barcharts_sig[[5]], wMKL_barcharts_sig[[6]],
          ncol = 2, nrow = 3, labels = c("A", "B", "C", "D", "E", "F"),
          font.label = list(size = 8, face = "bold", color ="black"))
ggsave(filename = paste0("sig_Multiplot_", algorithm, "_barcharts.png"),
       path = paste0(home, 
                     "/Results/single_algorithm/", algorithm, "/Supplement"), 
       width = 6700, height = 6500, device = 'png', units = "px",
       dpi = 700)
dev.off()

# Sunburst plot ###
library(plotly)
Pheno_sunburst_wMKL = clust_annot_pheno
Pheno_sunburst_wMKL$`ER status` = gsub("Unknown", "Unkn ER status", Pheno_sunburst_wMKL$`ER status`)
Pheno_sunburst_wMKL$`ER status` = gsub("Positive", "ER+", Pheno_sunburst_wMKL$`ER status`)
Pheno_sunburst_wMKL$`ER status` = gsub("Negative", "ER-", Pheno_sunburst_wMKL$`ER status`)
Pheno_sunburst_wMKL$`HER2 status` = gsub("Unknown", "Unkn HER2 status", 
                                         Pheno_sunburst_wMKL$`HER2 status`)
Pheno_sunburst_wMKL$`HER2 status` = gsub("Positive", "HER2+", Pheno_sunburst_wMKL$`HER2 status`)
Pheno_sunburst_wMKL$`HER2 status` = gsub("Negative", "HER2-", Pheno_sunburst_wMKL$`HER2 status`)
Pheno_sunburst_wMKL$Metastasis = gsub("Unknown", "Unkn metast. status", Pheno_sunburst_wMKL$Metastasis)
Pheno_sunburst_wMKL = Pheno_sunburst_wMKL %>%
  dplyr::select(wMKL, `ER status`, `HER2 status`, Metastasis) %>%
  group_by(wMKL, `ER status`, `HER2 status`, Metastasis) %>%
  summarise(Counts = n()) %>%
  as.data.frame()
sunburst_coloring_wMKL = data.frame(stringsAsFactors = FALSE,
                                    colors = tolower(gplots::col2hex(c("#2EC4B6", "#E71D36", "#FF9F1C", "#BDD5EA",
                                                                       "#FFA5AB", "#011627", "#023E8A", "#9D4EDD", 
                                                                       "#C11D9C", "#0F1682",  "grey40",
                                                                       "#0B9EF8", "#560DA7", "mistyrose1", 
                                                                       "hotpink4", "grey40",
                                                                       "deeppink4", "cadetblue2", "grey40"))),
                                    labels = c("wMKL1", "wMKL2", "wMKL3", "wMKL4",
                                               "wMKL5", "wMKL6", "wMKL7", "wMKL8",
                                               "ER-", "ER+", "Unkn ER status",
                                               "HER2-", "HER2+", "Indeterminate",
                                               "Equivocal", "Unkn HER2 status",
                                               "Yes", "No", "Unkn metast. status"))

sunburstDF_wMKL = as.sunburstDF(Pheno_sunburst_wMKL, value_column = "Counts", add_root = FALSE) %>%
  inner_join(sunburst_coloring_wMKL, by = "labels")

pie_wMKL = plot_ly() %>%
  add_trace(ids = sunburstDF_wMKL$ids, labels= sunburstDF_wMKL$labels, 
            parents = sunburstDF_wMKL$parents, 
            values= sunburstDF_wMKL$values, type='sunburst', branchvalues = 'total',
            insidetextorientation='radial', maxdepth = 5,
            marker = list(colors = sunburstDF_wMKL$colors)) %>%
  layout(
    grid = list(columns =1, rows = 1),
    margin = list(l = 0, r = 0, b = 0, t = 0)
  )
pie_wMKL
rm(Pheno_sunburst_wMKL, sunburstDF_wMKL, sunburst_coloring_wMKL, pie_wMKL); gc()

# Graphs ###
library(igraph)
list_aff_S = list(final_affinity_matrix)
names(list_aff_S) = c(paste0("Final Fused wMKL matrix"))

quantile_thresh = 0.5
for (i in seq_along(list_aff_S)) {
  g <- graph_from_adjacency_matrix(
    list_aff_S[[i]],
    mode     = "max",
    weighted = TRUE,
    diag     = FALSE)
  
  # drop the zero-weight edges
  g <- delete_edges(g, E(g)[E(g)$weight == 0])
  
  # Edge threshold for drawing
  thresh      <- quantile(E(g)$weight, quantile_thresh)          # 4th quartile
  keep_edge   <- E(g)$weight >= thresh               # logical mask
  
  ## edge-specific plotting attributes
  E(g)$plot_width  <- ifelse(keep_edge,
                             sqrt(E(g)$weight)*10,      # visible edges
                             0)                      # invisible edges
  E(g)$plot_color  <- ifelse(keep_edge,
                             "gray85",               # visible color
                             NA)                     # NA 
  
  nodes_data <- data.frame(name = V(g)$name) %>%
    inner_join(clust_annot_pheno %>% select(samID, wMKL),
               by = c("name" = "samID"))
  
  nodes_data[[algorithm]] <- as.factor(nodes_data[[algorithm]])
  V(g)$wMKL <- nodes_data[[algorithm]]
  
  # Set color based on wMKL
  V(g)$color <- fifelse(V(g)$wMKL == paste0(algorithm, "1"), "#2EC4B6", 
                        fifelse(V(g)$wMKL == paste0(algorithm, "2"), "#E71D36",
                                fifelse(V(g)$wMKL == paste0(algorithm, "3"), "#FF9F1C",
                                        fifelse(V(g)$wMKL == paste0(algorithm, "4"), "#BDD5EA",
                                                fifelse(V(g)$wMKL == paste0(algorithm, "5"), "#FFA5AB",
                                                        fifelse(V(g)$wMKL == paste0(algorithm, "6"), "#011627",
                                                                fifelse(V(g)$wMKL == paste0(algorithm, "7"), "#023E8A", "#9D4EDD")))))))
  
  png(paste0(home,
             "/Results/single_algorithm/", algorithm, "/Supplement/",
             names(list_aff_S)[i], ".png"),
      width = 6000, height = 6000, res = 700)
  
  par(mar = c(2, 2, 2, 5))
  
  plot(g,
       layout       = layout_with_fr(g),
       vertex.color = V(g)$color,
       vertex.size  = 4,
       vertex.label = NA,
       edge.width   = E(g)$plot_width,
       edge.color   = E(g)$plot_color,
       main         = "")
  
  title(main = names(list_aff_S)[i], cex.main = 1.7)
  
  legend("bottomright",
         title  = "Node Color Legend",
         legend = paste0(algorithm, 1:8),
         fill   = cluster_colors_heatmap,
         cex    = 0.7,
         box.lwd = 1)
  
  dev.off()
}
rm(g, nodes_data)

# Wrap up #####
hyperparameters = list()

# Put all parameters in a list
params = list(algorithm = algorithm, data_source = data_source, data_types = data_types,
              evaluation_source = evaluation_source, title = title, subtitle = subtitle,
              description = description, in_a_nutshell = in_a_nutshell, optk_text = optk_text,
              citation = citation, NMI_to_MOVICS = NMI_to_MOVICS, ARI_to_MOVICS = ARI_to_MOVICS,
              hyperparameters = hyperparameters, ground_truth_k = ground_truth_k,
              transNEO_var2comp = transNEO_var2comp,
              sessionInfo = sessionInfo(), home = home)

# Render the R Markdown document with the parameters
rmarkdown::render(paste0(getwd(), "/Results/single_algorithm/", algorithm, "/", algorithm, "_report.Rmd"), 
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
