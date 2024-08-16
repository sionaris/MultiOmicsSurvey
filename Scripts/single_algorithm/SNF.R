# Import data from gitignored "Resources/BRCA complete/" folder #####

# This RDS object was produced using the Scripts/MOVICS/MOVICS_baseline.R script
input = readRDS("Resources/BRCA complete/mm_input.rds")

# Setup environment variables for markdown #####

# Ensure reproducibility
RNGversion("4.2.2")
set.seed(123)

# Load custom helper functions
source("Scripts/automated_scripts/custom_functions.R")
source("Scripts/automated_scripts/modified_MOVICS_functions.R")

# Preamble
home = getwd()
algorithm = "SNF"
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
ground_truth_k = 3 # optk from MOVICS
optk_boolean = "FALSE" # either TRUE or FALSE. Answers whether the algorithm suggests an optimal k
optk_text = ifelse(optk_boolean == TRUE,
                   "<u>suggests</u> an estimate of the optimal number of multi-omic clusters $k$",
                   "<u>does not suggest</u> an optimal number of multi-omic clusters $k$")

# Detailed description of the algorithm
description = paste(readLines(paste0("Resources/algorithm_descriptions/", algorithm,
                                     "_description.Rmd")),
                    collapse = "\n") # File path to .Rmd file within Resources/algorithm_descriptions

# Preprocessing flags and code #####
library(stringr)
library(dplyr)

# All preprocessing for this input has already been performed using the 
# Scripts/MOVICS/MOVICS_baseline.R script

# However SNF prefers features in columns so we transpose the matrices.

# Extract the names of the modalities that will be used
modalities = unlist(strsplit(data_types, "-"))

# Replace with TRUE wherever features are in rows
features_in_rows = rep(TRUE, length(modalities))

# Run algorithm #####
library(SNFtool)

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
clinical_data = openxlsx::read.xlsx("Resources/BRCA complete/clinical_data.xlsx")

# Setup ###
# Hyperparameter tuning
sigma_step = 0.1
# iter_step = 10
neighbor_step = 5
num_neighbors_range = seq(10, 30, neighbor_step) # number of neighbors, usually (10~30)
sigma_range = seq(0.3, 0.8, sigma_step) 	# hyperparameter, usually (0.3~0.8) -REFFERED as \mu in report text
# iterations = seq(10, 100, iter_step)      # Number of Iterations, usually (10~20)
n_iterations = 20     # Number of Iterations, usually (10~20)

# Modality types
continuous = c("RNAseq", "CNV", "Methylation", "miRNA")
categorical = c("SNPs")

# Calculate the pair-wise distance (Euclidean for continuous modalities)
input_dists = lapply(input[continuous], function(x) {
  x = as.matrix(x)
  x = dist2(x, x)
})

# Binary for SNPs (see ?dist for details)
input_dists[["SNPs"]] = as.matrix(dist(as.matrix(input$SNPs),
                                       as.matrix(input$SNPs),
                                       method = "binary"))
gc()

# Affinity matrices ###
affinity_object = list()
for (nn in num_neighbors_range) {
  nn_list <- list()
  for (sigma in sigma_range) {
    sigma_list <- list()
    for (modality in modalities) {
      aff_mat <- affinityMatrix(as.matrix(input_dists[[modality]]), K = nn, sigma = sigma)
      sigma_list[[modality]] <- list(
        num_neighbors = nn,
        regularization = sigma,
        affinity_matrix = aff_mat
      )
      rm(aff_mat)
    }
    nn_list[[paste0("sigma = ", sigma)]] <- sigma_list
    rm(sigma_list)
  }
  affinity_object[[paste0("NN = ", nn)]] <- nn_list
  rm(nn_list)
}
rm(nn, sigma)

# Fusions ###
# Try parallel
library(parallel)
library(foreach)
library(doParallel)

# Use 5 cores
cl = makeCluster(6) # or 3 (depending on resources)
registerDoParallel(cl)

Fusions = list()

# Use foreach to parallelize the computation (< 30 min)
Fusions <- foreach(nn = num_neighbors_range, .combine = 'c', .packages = 'SNFtool') %:%
  foreach(sigma = sigma_range, .combine = 'c') %dopar% {
    sublist <- affinity_object[[paste0("NN = ", nn)]][[paste0("sigma = ", sigma)]]
    iter_matrices <- lapply(sublist, `[[`, "affinity_matrix")
    fusion_result <- SNF(iter_matrices, K = nn, t = n_iterations, parallel = FALSE)
    list(fusion_result)
  }

# Stop the cluster
stopCluster(cl)

# Restructure the Fusions list to match the desired output format
names(Fusions) <- unlist(lapply(num_neighbors_range, function(nn) {
  lapply(sigma_range, function(sigma) {
    paste0("NN = ", nn, ", sigma = ", sigma)
  })
}))

# Give appropriate colnames and rownames
column_names = colnames(affinity_object[["NN = 10"]][["sigma = 0.3"]][["RNAseq"]][["affinity_matrix"]])
row_names = rownames(affinity_object[["NN = 10"]][["sigma = 0.3"]][["RNAseq"]][["affinity_matrix"]])

for (i in 1:length(Fusions)) {
  colnames(Fusions[[i]]) = column_names
  rownames(Fusions[[i]]) = row_names
}

save.image(paste0(home, "/Results/single_algorithm/", 
                  algorithm, "/", algorithm, "_", data_source, "_",
                  data_types, "_eval_on_", evaluation_source,
                  "_env.RData"))

# Function to compute both Frobenius norm and Pearson correlation between matrices
compute_matrix_similarity <- function(matrices) {
  num_matrices <- length(matrices)
  similarity_frobenius <- matrix(0, nrow = num_matrices, ncol = num_matrices)
  similarity_pearson <- matrix(0, nrow = num_matrices, ncol = num_matrices)
  
  for (i in 1:num_matrices) {
    for (j in 1:num_matrices) {
      if (i != j) {
        similarity_frobenius[i, j] <- frobenius_norm(matrices[[i]], matrices[[j]])
        similarity_pearson[i, j] <- pearson_correlation(matrices[[i]], matrices[[j]])
      }
    }
  }
  
  # Set row names and column names
  
  rownames(similarity_frobenius) <- colnames(similarity_frobenius) <- 
    rownames(similarity_pearson) <- colnames(similarity_pearson) <- names(matrices)
  
  return(list(Frobenius = similarity_frobenius, Pearson = similarity_pearson))
}

# Check similarities for a given nn
nn_similarities = list()
for (nn in num_neighbors_range) {
  indices = grepl(paste0("NN = ", nn), names(Fusions))
  matrices = Fusions[indices]
  similarity_results <- compute_matrix_similarity(matrices)
  
  # Modify row names and column names for the similarity matrices
  matrix_names <- substr(names(matrices), 9, 20)
  rownames(similarity_results$Frobenius) <- matrix_names
  colnames(similarity_results$Frobenius) <- matrix_names
  rownames(similarity_results$Pearson) <- matrix_names
  colnames(similarity_results$Pearson) <- matrix_names
  
  nn_similarities[[paste0("NN = ", nn)]] <- similarity_results
}

print(nn_similarities)

# Check similarities for a given sigma
sigma_similarities = list()
for (sigma in sigma_range) {
  indices = grepl(paste0("sigma = ", sigma), names(Fusions))
  matrices = Fusions[indices]
  similarity_results <- compute_matrix_similarity(matrices)
  
  # Modify row names and column names for the similarity matrices
  matrix_names <- substr(names(matrices), 0, 7)
  rownames(similarity_results$Frobenius) <- matrix_names
  colnames(similarity_results$Frobenius) <- matrix_names
  rownames(similarity_results$Pearson) <- matrix_names
  colnames(similarity_results$Pearson) <- matrix_names
  
  sigma_similarities[[paste0("sigma = ", sigma)]] <- similarity_results
}

print(sigma_similarities)
rm(indices, matrices); gc()

# All similarities
all_similarities = compute_matrix_similarity(Fusions)

# Overall tests
# Reshape data for ANOVA
nn_reshape <- reshape_SNF_Pearson_for_anova(nn_similarities)
sigma_reshape <- reshape_SNF_Pearson_for_anova(sigma_similarities)

# Perform ANOVA for nn
anova_nn <- aov(Value ~ Factor, data = nn_reshape)
summary(anova_nn)

# Perform ANOVA for sigma
anova_sigma <- aov(Value ~ Factor, data = sigma_reshape)
summary(anova_sigma)

# Conclusion
if (summary(anova_nn)[[1]][["Pr(>F)"]][1] < 0.05) {
  conclusion1 = "Overall, the choice of sigma significantly affects the results for a given nn."
  cat(conclusion1)
} else {
  conclusion1 = "Overall, the choice of sigma does not significantly affect the results for a given nn."
  cat(conclusion1)
}

if (summary(anova_sigma)[[1]][["Pr(>F)"]][1] < 0.05) {
  conclusion2 = "Overall, the choice of nn significantly affects the results for a given sigma."
  cat(conclusion2)
} else {
  conclusion2 = "Overall, the choice of nn does not significantly affect the results for a given sigma."
  cat(conclusion2)
}

# Determine overall effect
if (summary(anova_nn)[[1]][["Pr(>F)"]][1] < 0.05 &&
    summary(anova_sigma)[[1]][["Pr(>F)"]][1] < 0.05) {
  nn_mean_diff <- max(nn_summary$mean) - min(nn_summary$mean)
  sigma_mean_diff <- max(sigma_summary$mean) - min(sigma_summary$mean)
  
  # Means
  if (nn_mean_diff > sigma_mean_diff) {
    conclusion3 = paste0("The nn effect is stronger than the sigma effect based on mean",
                         " Pearson similarities (", nn_mean_diff, " vs. ", sigma_mean_diff, ").")
    cat(conclusion3)
  } else if (nn_mean_diff < sigma_mean_diff) {
    conclusion3 = paste0("The sigma effect is stronger than the nn effect based on mean",
                         " Pearson similarities (", sigma_mean_diff, " vs. ", nn_mean_diff, ").")
    cat(conclusion3)
  }
  
  # Standard deviations
  nn_sd_diff <- max(nn_summary$sd) - min(nn_summary$sd)
  sigma_sd_diff <- max(sigma_summary$sd) - min(sigma_summary$sd)
  
  if (nn_sd_diff > sigma_sd_diff) {
    conclusion4 = paste0("The nn effect is stronger than the sigma effect based on", 
                         " the standard deviation of Pearson similarities (",
                         nn_sd_diff, " vs. ", sigma_sd_diff, ").")
    cat(conclusion4)
  } else if (nn_sd_diff < sigma_sd_diff) {
    conclusion4 = paste0("The sigma effect is stronger than the nn effect based on",  
                         " the standard deviation of Pearson similarities (",
                         sigma_sd_diff, " vs. ", nn_sd_diff, ").")
    cat(conclusion4)
  }
} else {
  conclusion3 = "The nn effect and sigma effect are practically equal based on mean Pearson similarities"
  cat(conclusion3)
}

if (exists("conclusion4")) {
  conclusion = paste0(conclusion1, conclusion2, conclusion3, conclusion4, collapse = " ")
  rm(conclusion1, conclusion2, conclusion3, conclusion4)
  sig_status = TRUE
} else {
  conclusion = paste0(conclusion1, conclusion2, conclusion3, collapse = " ")
  rm(conclusion1, conclusion2, conclusion3)
  sig_status = FALSE
}

# If no significant differences are shown between/across hyperparameters then pick median values
if (sig_status) {
  # Code to pick best hyperparameters
} else {
  optN = 20 # median(num_neighbors_range)
  optSigma = 0.5 # ~median(sigma_range)
}

# Spectral clustering for k = ground_truth_k from MOVICS
RNGversion("4.2.2")
set.seed(123)

final_affinity_matrix = Fusions[[paste0("NN = ", optN, ", sigma = ", optSigma)]]
group = spectralClustering(final_affinity_matrix, ground_truth_k)
names(group) = colnames(final_affinity_matrix)

SNF_clusters = as.data.frame(list(Sample.ID = names(group),
                                  Cluster = group))

# Main results ###
# Examine cluster similarity to MOVICS by measuring NMI and ARI indices #####
# (Jaccard may be misleading)

# Calculate ARI and NMI
library(mclust)
library(clue)

ARI_to_MOVICS = calculate_ari_index(cluster_df1 = ground_truth_labels,
                                    cluster_df2 = SNF_clusters,
                                    sample_col = "Sample.ID",
                                    clust_col = "Cluster",
                                    suffixes = c("_MOVICS", "_SNF"))

NMI_to_MOVICS = calculate_nmi_index(cluster_df1 = ground_truth_labels,
                                    cluster_df2 = SNF_clusters,
                                    sample_col = "Sample.ID",
                                    clust_col = "Cluster",
                                    suffixes = c("_MOVICS", "_SNF"))

# MOVICS-like analysis #####
library(MOVICS)
library(ComplexHeatmap)

# Import coloring scheme
scheme = readRDS("Resources/scheme.rds")
annCol = scheme$annCol
annColors = scheme$annColors
cluster_colors = scheme$clust.colors
col.list = scheme$col.list
var2comp = scheme$var2comp
rm(scheme); gc()

plotdata <- lapply(lapply(input, as.matrix), 
                   function(mat) mat[, colSums(mat != 0) > 0])
plotdata <- lapply(plotdata, t)
plot_object = list(clust.res = SNF_clusters %>%
                     dplyr::rename(samID = Sample.ID, clust = Cluster))

# comprehensive heatmap (may take a while)
getMoHeatmap_single_algorithm(algorithm_name = algorithm,
                              data          = plotdata,
             row.title     = names(plotdata),
             is.binary     = c(F,F,F,F,T), 
             legend.name   = c("Normalised CNV",
                               "Normalised Methylation",
                               "Normalised RNAseq FPKM",
                               "Normalised miRNA FPKM",
                               "SNPs"
                               #bquote(bold("Normalised" ~ log[2]("TPM + 1")))
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
             fig.path      = paste0(home, "/Results/single_algorithm/SNF"),
             fig.name      = paste0("default_", algorithm, "_Comprehensive_heatmap"))
gc()

# Clinical variables ###
# Statistical comparisons
clin_comp = compClinvar_single_algorithm(algorithm_name = algorithm,
                                         moic.res = plot_object,
                        var2comp = var2comp,
                        strata = "Consensus Subtype",
                        factorVars = c("vital_status", "race_list", "ethnicity",
                                       "history_of_neoadjuvant_treatment",
                                       "primary_lymph_node_presentation_assessment",
                                       "histological_type", "menopause_status",
                                       "breast_carcinoma_progesterone_receptor_status",
                                       "breast_carcinoma_estrogen_receptor_status",
                                       "lab_proc_her2_neu_immunohistochemistry_receptor_status",
                                       "distant_metastasis_present_ind2",
                                       "stage_event_pathologic_stage"),
                        includeNA = FALSE,
                        doWord = TRUE,
                        tab.name = "Summary_of_clinical_variables",
                        res.path = paste0(home, "/Results/single_algorithm/SNF/"))

# race_list, ER status, PR status, metastasis are sig

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
                     fig.path     = paste0(home, "/Results/single_algorithm/SNF"),
                     res.path     = paste0(home, "/Results/single_algorithm/SNF"))

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
                                fig.path = paste0(home, "/Results/single_algorithm/SNF"))

# Agreement with other subtypes ###
subtype_agreement <- compAgree_single_algorithm(algorithm_name = algorithm,
                                                moic.res  = plot_object,
                                subt2comp = annCol[, c("ER status", "PR status",
                                                       "HER2 status", "Metastasis", "Stage")],
                                doPlot    = TRUE,
                                box.width = 0.2,
                                fig.name  = "Classification_agreement",
                                fig.path  = paste0(home, "/Results/single_algorithm/SNF"),
                                width     = 12)

# DGEA ###
dgea = runDEA(dea.method = "limma", # we use normalized data as input
              expr = plotdata$RNAseq,
              moic.res = plot_object,
              prefix = "dgea_",
              sort.p = TRUE,
              overwt = TRUE,
              verbose = TRUE,
              res.path = paste0(home, "/Results/single_algorithm/SNF"))

# # Identify unique subtype biomarkers
# # 1. Up-regulated markers
dgea.marker.up <- runMarker_single_algorithm(algorithm_name = algorithm,
                                             moic.res = plot_object,
                                    dea.method    = "limma", # name of DEA method
                                    prefix        = "dgea_", # MUST be the same of argument in runDEA()
                                    dat.path      = paste0(home, "/Results/single_algorithm/SNF"), # path of DEA files
                                    res.path      = paste0(home, "/Results/single_algorithm/SNF"), # path to save marker files
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
                                    fig.name      = "upregulated_biomarkers_heatmap",
                                    fig.path = paste0(home, "/Results/single_algorithm/SNF"),
                                    width = 14,
                                    height = 12,
                                    fontsize_row = 3,
                                    name = "normalized RNA-seq")

# # 2. Down-regulated markers
dgea.marker.down <- runMarker_single_algorithm(algorithm_name = algorithm,
                                               moic.res = plot_object,
                                      dea.method    = "limma", # name of DEA method
                                      prefix        = "dgea_", # MUST be the same of argument in runDEA()
                                      dat.path      = paste0(home, "/Results/single_algorithm/SNF"), # path of DEA files
                                      res.path      = paste0(home, "/Results/single_algorithm/SNF"), # path to save marker files
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
                                      fig.name      = "downregulated_biomarkers_heatmap",
                                      fig.path = paste0(home, "/Results/single_algorithm/SNF"),
                                      width = 14,
                                      height = 12,
                                      fontsize_row = 3,
                                      name = "normalized RNA-seq")

# GSEA ###
# Load MSigDb file
MSIGDB.FILE <- system.file("extdata", "c5.bp.v7.1.symbols.xls", package = "MOVICS", mustWork = TRUE)

# GSEA up-regulated
RNGversion("4.2.2")
set.seed(123)
gsea.up <- runGSEA_mod_4.4_single_algorithm(algorithm_name = algorithm,
                                            moic.res     = plot_object,
                           dea.method   = "limma", # name of DEA method
                           prefix       = "dgea_", # MUST be the same of argument in runDEA()
                           dat.path      = paste0(home, "/Results/single_algorithm/SNF"), # path of DEA files
                           res.path      = paste0(home, "/Results/single_algorithm/SNF"), # path to save marker files
                           msigdb.path  = MSIGDB.FILE, # MUST be the ABSOLUTE path of msigdb file
                           norm.expr    = plotdata$RNAseq, # use normalized expression to calculate enrichment score
                           dirct        = "up", # direction of dysregulation in pathway
                           n.path       = 20,
                           p.cutoff     = 0.05, # p cutoff to identify significant pathways
                           p.adj.cutoff = 0.1, # padj cutoff to identify significant pathways
                           gsva.method  = "gsva", # method to calculate single sample enrichment score
                           name         = "GSVA scores", # name for colorbar
                           norm.method  = "mean", # normalization method to calculate subtype-specific enrichment score
                           fig.name     = "upregulated_pathway_heatmap",
                           nPerm = 10000,
                           minGSSize = 10,
                           maxGSSize = 500,
                           fig.path = paste0(home, "/Results/single_algorithm/SNF"),
                           width = 14, height = 12)

# GSEA down-regulated
RNGversion("4.2.2")
set.seed(123)
gsea.down <- runGSEA_mod_4.4_single_algorithm(algorithm_name = algorithm,
                                              moic.res     = plot_object,
                             dea.method   = "limma", # name of DEA method
                             prefix       = "dgea_", # MUST be the same of argument in runDEA()
                             dat.path      = paste0(home, "/Results/single_algorithm/SNF"), # path of DEA files
                             res.path      = paste0(home, "/Results/single_algorithm/SNF"), # path to save marker files
                             msigdb.path  = MSIGDB.FILE, # MUST be the ABSOLUTE path of msigdb file
                             norm.expr    = plotdata$RNAseq, # use normalized expression to calculate enrichment score
                             dirct        = "down", # direction of dysregulation in pathway
                             n.path       = 20,
                             p.cutoff     = 0.05, # p cutoff to identify significant pathways
                             p.adj.cutoff = 0.1, # padj cutoff to identify significant pathways
                             gsva.method  = "gsva", # method to calculate single sample enrichment score
                             name         = "GSVA scores", # name for colorbar
                             norm.method  = "mean", # normalization method to calculate subtype-specific enrichment score
                             fig.name     = "downregulated_pathway_heatmap",
                             nPerm = 10000,
                             minGSSize = 10,
                             maxGSSize = 500,
                             fig.path = paste0(home, "/Results/single_algorithm/SNF"),
                             width = 14, height = 12)

# Gene set variation analysis #####
# locate ABSOLUTE path of gene set file
GSET.FILE <- system.file("extdata", "gene sets of interest.gmt", package = "MOVICS", mustWork = TRUE)

RNGversion("4.2.2")
set.seed(123)
gsva.res = runGSVA_mod_4.4_single_algorithm(algorithm_name = algorithm,
                                            moic.res      = plot_object,
                           norm.expr     = plotdata$RNAseq,
                           gset.gmt.path = GSET.FILE, # ABSOLUTE path of gene set file
                           gsva.method   = "gsva", # method to calculate single sample enrichment score
                           annCol        = annCol,
                           annColors     = annColors,
                           fig.path      = paste0(home, "/Results/single_algorithm/SNF"),
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

# Evaluation #####
# Run Nearest Template Prediction in transNEO cohort ###
# Load transNEO data
transNEO_mm_inputs = readRDS("Resources/transNEO/transNEO_multimodal_inputs.rds")
transcr = transNEO_mm_inputs$`RNAseq log2(TPM+1)`[, 1:153]
rownames(transcr) = transNEO_mm_inputs$`RNAseq log2(TPM+1)`$Hugo

# Up-regulated expression features
dgea.marker.up[["templates"]][["class"]] = gsub("CS", algorithm, 
                                                dgea.marker.up[["templates"]][["class"]])
RNGversion("4.2.2")
transNEO_ntp_expr_up = runNTP(
  expr = as.matrix(transcr),
  templates = dgea.marker.up$templates,
  scaleFlag = TRUE,
  centerFlag = TRUE,
  nPerm = 10000,
  seed = 123,
  distance = "cosine", # default
  doPlot = TRUE,
  height = 8,
  width = 12,
  fig.path = paste0(home, "/Results/single_algorithm/SNF"),
  fig.name = "ntp_expr_up_heatmap_transNEO")

# down-regulated
dgea.marker.down[["templates"]][["class"]] = gsub("CS", algorithm, 
                                                dgea.marker.down[["templates"]][["class"]])
RNGversion("4.2.2")
transNEO_ntp_expr_down = runNTP(
  expr = as.matrix(transcr),
  templates = dgea.marker.down$templates,
  scaleFlag = TRUE, # already standardised
  centerFlag = TRUE, # -//-
  nPerm = 10000,
  seed = 123,
  distance = "cosine", # default
  doPlot = TRUE,
  height = 8,
  width = 12,
  fig.path = paste0(home, "/Results/single_algorithm/SNF"),
  fig.name = "ntp_expr_down_heatmap_transNEO")

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
  inner_join(expr_conc %>% dplyr::select(Donor.ID = samID, `Consensus Subtype` = clust_up),
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


transNEO_clincomp = compClinvar_single_algorithm(algorithm_name = algorithm,
                                                 moic.res = transNEO_ntp_expr_up,
                                 var2comp = transNEO_var2comp,
                                 strata = "Consensus Subtype",
                                 factorVars = c("ER.status", "HER2.status", "Grade.pre.NAT",
                                                "NAT.regimen", 
                                                "pCR.RD", "LN.status.at.diagnosis"),
                                 includeNA = FALSE,
                                 doWord = TRUE,
                                 tab.name = "transNEO_Summary_of_clinical_variables",
                                 res.path = paste0(home, "/Results/single_algorithm/SNF"))

# Run PAM ###
RNGversion("4.2.2.")
set.seed(123)
transNEO_pam = runPAM_single_algorithm(algorithm_name = algorithm,
                                       train.expr = plotdata$RNAseq,
                      moic.res   = plot_object,
                      test.expr  = as.matrix(transcr))

# Check consistency across methods

# Get predictions for TCGA (discovery cohort)
RNGversion("4.2.2.")
set.seed(123)
TCGA.ntp.pred = runNTP(expr = plotdata$RNAseq[, plot_object$clust.res$samID],
                       templates = dgea.marker.up$templates,
                       doPlot = F)

TCGA.pam.pred = runPAM_single_algorithm(algorithm_name = algorithm,
                                        train.expr = plotdata$RNAseq[, plot_object$clust.res$samID],
                       moic.res = plot_object,
                       test.expr = plotdata$RNAseq[, plot_object$clust.res$samID])

# consensus TCGA vs NTP TCGA # FAILS
runKappa_single_algorithm(algorithm_name = algorithm,
                          subt1 = plot_object$clust.res$clust,
         subt2 = gsub(algorithm, "", TCGA.ntp.pred$clust.res$clust),
         subt1.lab = "SNF",
         subt2.lab = "NTP TCGA",
         height = 8,
         width = 8,
         fig.path = paste0(home, "/Results/single_algorithm/SNF"),
         fig.name = paste0("kappa_", algorithm, "_vs_NTP_TCGA"))

# consensus TCGA vs PAM TCGA
runKappa_single_algorithm(algorithm_name = algorithm,
                          subt1 = plot_object$clust.res$clust,
         subt2 = gsub(algorithm, "", TCGA.pam.pred$clust.res$clust),
         subt1.lab = "SNF",
         subt2.lab = "PAM TCGA",
         height = 8,
         width = 8,
         fig.path = paste0(home, "/Results/single_algorithm/SNF"),
         fig.name = paste0("kappa_", algorithm, "_vs_PAM_TCGA"))

# NTP transNEO vs PAM transNEO # FAILS
runKappa_single_algorithm(algorithm_name = algorithm,
                           subt1 = as.numeric(gsub(algorithm, "",
                                                   transNEO_ntp_expr_up$clust.res$clust)),
         subt2 = as.numeric(transNEO_pam$clust.res$clust),
         subt1.lab = "transNEO NTP",
         subt2.lab = "transNEO PAM",
         height = 8,
         width = 8,
         fig.path = paste0(home, "/Results/single_algorithm/SNF"),
         fig.name = "kappa_NTP_vs_PAM_transNEO")

# Export consensus clustering object
clust = as.data.frame(plot_object$clust.res)
colnames(clust) = c("Sample.ID", "Cluster")
clust$Cluster = paste0("SNF", clust$Cluster)
openxlsx::write.xlsx(clust, paste0(home, "/Results/single_algorithm/SNF/", 
                                   algorithm, "_", data_source, "_",
                                   data_types, "_eval_on_", evaluation_source,
                                   "_clusterings.xlsx"))

# Supplementary results #####

# Create subdirectory for supplementary plots
if (!dir.exists(paste0(home, "/Results/single_algorithm/SNF/Supplement"))) {
  dir.create(paste0(home, "/Results/single_algorithm/SNF/Supplement"))
}

# Setup for heatmaps
colors_heatmap = rev(colorRampPalette(viridisLite::magma(10))(255))
cluster_colors_heatmap = c("#2EC4B6", "#E71D36", "#FF9F1C")
clust_annot_pheno = annCol %>% mutate(Sample.ID = rownames(.)) %>%
  inner_join(clust, by = "Sample.ID") %>%
  dplyr::rename(SNF = Cluster, samID = "Sample.ID")
rownames(clust_annot_pheno) = clust_annot_pheno$samID
afh_colnames = colnames(annCol)

# Prepare affinity matrices
aff_CNV = normalize_affinity_matrix(
  SNFtool::affinityMatrix(
    SNFtool::dist2(input$CNV,
    input$CNV),
    K = 30, sigma = 0.5))
colnames(aff_CNV) = rownames(aff_CNV) = rownames(input$CNV)

aff_rna = normalize_affinity_matrix(
  SNFtool::affinityMatrix(
    SNFtool::dist2(input$RNAseq,
                   input$RNAseq),
    K = 30, sigma = 0.5))
colnames(aff_rna) = rownames(aff_rna) = rownames(input$RNA)

aff_miRNA = normalize_affinity_matrix(
  SNFtool::affinityMatrix(
    SNFtool::dist2(input$miRNA,
                   input$miRNA),
  K = 30, sigma = 0.5))
colnames(aff_miRNA) = rownames(aff_miRNA) = rownames(input$miRNA)

aff_Methyl = normalize_affinity_matrix(
  SNFtool::affinityMatrix(
    SNFtool::dist2(input$Methylation,
                   input$Methylation),
    K = 30, sigma = 0.5))
colnames(aff_Methyl) = rownames(aff_Methyl) = rownames(input$Methylation)

aff_SNPs = normalize_affinity_matrix(
  SNFtool::affinityMatrix(
    as.matrix(dist(as.matrix(input$SNPs),
                   as.matrix(input$SNPs),
                   method = "binary")),
  K = 30, sigma = 0.5))
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
                  output_file_name = paste0(home, "/Results/single_algorithm/SNF/Supplement/aff_CNV_heatmap.png"))

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
                  output_file_name = paste0(home, "/Results/single_algorithm/SNF/Supplement/aff_RNAseq_heatmap.png"))

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
                  output_file_name = paste0(home, "/Results/single_algorithm/SNF/Supplement/aff_miRNA_heatmap.png"))

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
                  output_file_name = paste0(home, "/Results/single_algorithm/SNF/Supplement/aff_Methylation_heatmap.png"))

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
                  output_file_name = paste0(home, "/Results/single_algorithm/SNF/Supplement/aff_SNPs_heatmap.png"))

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
                  output_file_name = paste0(home, "/Results/single_algorithm/SNF/Supplement/aff_final_affinity_heatmap.png"))


# Wrap up #####
hyperparameters = list(num_neighbors_min = min(num_neighbors_range),
                       num_neighbors_max = max(num_neighbors_range),
                       num_neighbors_step = neighbor_step,
                       sigma_min = min(sigma_range),
                       sigma_max = max(sigma_range),
                       sigma_step = sigma_step,
                       optimal_N = optN,
                       optimal_sigma = optSigma,
                       conclusion = conclusion
                       )

# Put all parameters in a list
params = list(algorithm = algorithm, data_source = data_source, data_types = data_types,
              evaluation_source = evaluation_source, title = title, subtitle = subtitle,
              description = description, in_a_nutshell = in_a_nutshell, optk_text = optk_text,
              citation = citation, NMI_to_MOVICS = NMI_to_MOVICS, ARI_to_MOVICS = ARI_to_MOVICS,
              hyperparameters = hyperparameters)

# Create algorithm directory if it doesn't exist
if (!dir.exists(paste0(home, "/Results/single_algorithm/", 
                       algorithm))) {
  dir.create(paste0(home, "/Results/single_algorithm/", 
                    algorithm))
}

# Render the R Markdown document with the parameters
rmarkdown::render(paste0(getwd(), "/Scripts/automated_scripts/single_algorithm_results_report.Rmd"), 
                  params = params, 
                  output_file = paste0(home, "/Results/single_algorithm/", 
                                       algorithm, "/", algorithm, "_report_",
                                       data_source, "_",
                                       data_types, "_eval_on_", evaluation_source,
                                       ".html"))

# Save environment
save.image(paste0(home, "/Results/single_algorithm/", 
                  algorithm, "/", algorithm, "_", data_source, "_",
                  data_types, "_eval_on_", evaluation_source,
                  "_env.RData"))

# Session info
# Capture the output of sessionInfo() to a variable
session_info <- capture.output(sessionInfo())

# Write the captured output to a .txt file
writeLines(session_info, paste0(home, "/Results/single_algorithm/", 
                                algorithm, "/", algorithm, "_", data_source, "_",
                                data_types, "_eval_on_", evaluation_source,
                                "_session_info.txt"))
