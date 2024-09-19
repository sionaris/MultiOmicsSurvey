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

# Load PINSPlus
# install_github("danro9685/PINSPlus", ref = 'R')
library(PINSPlus)

# Preamble
home = getwd()
algorithm = "PINSPlus"
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
optk_boolean = "FALSE" # either TRUE or FALSE. Answers whether the algorithm suggests an optimal k
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

# However PINSPlus prefers features in columns so we transpose the matrices.

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

# Setup ###
# Hyperparameter tuning
neighbor_step = 5
num_neighbors_range = seq(10, 50, neighbor_step) # number of neighbors, usually (10~30)
gc()

# Similarity matrices ###
# Parallelization occurs alreade within the function. We do not parallelize further
library(Matrix)
library(parallel)
similarity_object = list()
for (nn in num_neighbors_range) {
  similarity_object[[paste0("NN = ", nn)]] = PINSPlus_mod(X = input, c = ground_truth_k, 
                                                          k = nn, binary_flags = c("Yes", "No", "No", "No", "No"),
                                                          binary_distance = "binary", nonbinary_distance = "sqeuclidean", cores.ratio = 0.25)
}
rm(nn)

# Compare similarity matrices similarly to what we did to SNF
sim_matrices_S = lapply(similarity_object, function(x) x[["S"]])
D2_matrices_F = lapply(similarity_object, function(x) x[["F"]])
names(sim_matrices_S) = names(D2_matrices_F) = names(similarity_object)

# Check similarities for a given nn
S_similarities = compute_matrix_similarity(sim_matrices_S)
dimnames(S_similarities$Frobenius) = dimnames(S_similarities$Pearson) =
  list(names(similarity_object), names(similarity_object))

# Examine Pearson matrices ###
S_Pearson_matrix <- S_similarities$Pearson
S_Frobenius_matrix <- S_similarities$Frobenius

S_Pearson_values <- S_Pearson_matrix[lower.tri(S_Pearson_matrix, diag = FALSE)]
S_mean_Pearson_value <- mean(S_Pearson_values)
S_median_Pearson_value <- median(S_Pearson_values)
S_sd_Pearson_value <- sd(S_Pearson_values)

S_Frobenius_values <- S_Frobenius_matrix[lower.tri(S_Frobenius_matrix, diag = FALSE)]
S_mean_Frobenius_value <- mean(S_Frobenius_values)
S_median_Frobenius_value <- median(S_Frobenius_values)
S_sd_Frobenius_value <- sd(S_Frobenius_values)

# Plot histogram of Pearson values
library(ggplot2)
ggplot(data = data.frame(S_Pearson_values), aes(x = S_Pearson_values)) +
  geom_histogram(breaks = seq(0, 1, length.out = 37),
                 fill = "skyblue", color = "lightblue", size = 0.15) +
  stat_density(aes(color = "Density"), geom = "line", size = 0.4) +
  geom_vline(aes(xintercept = S_mean_Pearson_value, color = "Mean"), size = 0.2) + 
  geom_vline(aes(xintercept = S_median_Pearson_value, color = "Median"), size = 0.2) + 
  geom_vline(aes(xintercept = S_mean_Pearson_value - S_sd_Pearson_value, color = "Mean - SD"), 
             linetype = "dashed", size = 0.2) + 
  geom_vline(aes(xintercept = S_mean_Pearson_value + S_sd_Pearson_value, color = "Mean + SD"), 
             linetype = "dashed", size = 0.2) +
  scale_color_manual(name = "Lines", values = c("Mean" = "red", "Median" = "orange", 
                                                "Mean - SD" = "grey25", "Mean + SD" = "grey25",
                                                "Density" = "darkblue")) +
  labs(title = "Histogram of Pearson values between S matrices for different values of Nearest Neighbors", 
       x = "S Matrix Pearson Values", y = "Frequency") +
  scale_x_continuous(name = "S Matrix Pearson Values", limits = c(0, 1),
                     breaks = seq(0, 1, 0.1), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(panel.background = element_blank(),
        axis.line = element_line(linewidth = 0.25),
        plot.title = element_text(face = "bold", size = 6.3),
        axis.title = element_text(face = "bold", size = 5.8),
        axis.text = element_text(size = 5),
        axis.ticks = element_line(linewidth = 0.2),
        legend.text = element_text(size = 4.5),
        legend.title = element_text(size = 5, face = "bold"),
        legend.key.spacing.y = unit(1, "mm"),
        legend.key.size = unit(0.25, "cm"),
        legend.box.background = element_rect(color = "black"))
ggsave(filename = paste0(algorithm, "_S_matrix_Pearson_similarity_histogram.pdf"),
       path = paste0(home, 
                     "/Results/single_algorithm/PINSPlus/Supplement"), 
       width = 2880, height = 1820, device = 'pdf', units = "px",
       dpi = 700)
dev.off()

# Plot histogram of Frobenius values
ggplot(data = data.frame(S_Frobenius_values), aes(x = S_Frobenius_values)) +
  geom_histogram(breaks = seq(0, 2, length.out = 37),
                 fill = "skyblue", color = "lightblue", size = 0.15) +
  stat_density(aes(color = "Density"), geom = "line", size = 0.4) +
  geom_vline(aes(xintercept = S_mean_Frobenius_value, color = "Mean"), size = 0.2) + 
  geom_vline(aes(xintercept = S_median_Frobenius_value, color = "Median"), size = 0.2) + 
  geom_vline(aes(xintercept = S_mean_Frobenius_value - S_sd_Frobenius_value, color = "Mean - SD"), 
             linetype = "dashed", size = 0.2) + 
  geom_vline(aes(xintercept = S_mean_Frobenius_value + S_sd_Frobenius_value, color = "Mean + SD"), 
             linetype = "dashed", size = 0.2) +
  scale_color_manual(name = "Lines", values = c("Mean" = "red", "Median" = "orange", 
                                                "Mean - SD" = "grey25", "Mean + SD" = "grey25",
                                                "Density" = "darkblue")) +
  labs(title = "Histogram of Frobenius values between S matrices for different values of Nearest Neighbors", 
       x = "S Matrix Frobenius Values", y = "Frequency") +
  scale_x_continuous(name = "S Matrix Frobenius Values", limits = c(0, 2),
                     breaks = seq(0, 2, 0.2), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(panel.background = element_blank(),
        axis.line = element_line(linewidth = 0.25),
        plot.title = element_text(face = "bold", size = 6.3),
        axis.title = element_text(face = "bold", size = 5.8),
        axis.text = element_text(size = 5),
        axis.ticks = element_line(linewidth = 0.2),
        legend.text = element_text(size = 4.5),
        legend.title = element_text(size = 5, face = "bold"),
        legend.key.spacing.y = unit(1, "mm"),
        legend.key.size = unit(0.25, "cm"),
        legend.box.background = element_rect(color = "black"))
ggsave(filename = paste0(algorithm, "_S_matrix_Frobenius_similarity_histogram.pdf"),
       path = paste0(home, 
                     "/Results/single_algorithm/PINSPlus/Supplement"), 
       width = 2880, height = 1820, device = 'pdf', units = "px",
       dpi = 700)
dev.off()

# Evidently, Pearson correlations are on the lower extreme and Frobenius norms are high
conclusion1 = paste0("Average Pearson similarity across S matrices produced by different values of $nn'$ was ",
                     S_mean_Pearson_value, ", which indicates generally ",
                     ifelse(S_mean_Pearson_value < 0.75, "dissimilar", "similar"),
                     " S matrices across $nn'$ values.")

# Handle sig_status_final
if (S_mean_Pearson_value > 0.75) {
  sig_status_final = FALSE
} else {
  sig_status_final = TRUE
}

# If no significant differences are shown between/across hyperparameters then pick median values
if (!sig_status_final){
  optNN = 15 # arbitrary, but yielded good results in SNF
}

# We then choose the nn value for which the
# fused similarity matrix has the "best" bimodal distribution of low and high values.
# WE USE THIS APPROACH ONLY BECAUSE THE NUMBER OF CLUSTERS WE SEEK IS 2!

# We combine two methodologies to do it:

# 1. Get the sum of variance and IQR for every matrix
# 2. Get the sum of absolute skewness and kurtosis
# 3. Find the nn matrix for which the sum of 1 and 2 is maximum

# Skewness and kurtosis
library(e1071)
contrast_list <- choose_matrix_contrasts(sim_matrices_S)
skewness_kurtosis_list <- choose_matrix_skewness_kurtosis(sim_matrices_S)
contrast_values <- unlist(contrast_list)
skewness_kurtosis_values <- unlist(skewness_kurtosis_list)

# Normalize the contrast values and skewness-kurtosis values
normalized_contrast <- minmax_normalize_values(contrast_values)
normalized_skewness_kurtosis <- minmax_normalize_values(skewness_kurtosis_values)

# Get the sum
combined_scores <- normalized_contrast + normalized_skewness_kurtosis
combined_scores_list <- setNames(as.list(combined_scores), names(contrast_list))
print(combined_scores_list)

best_combined_matrix <- names(combined_scores_list)[which.max(combined_scores)]
cat("Best similarity matrix based on combined normalized scores:", best_combined_matrix, "\n")

# We choose nn = 20
optNN = as.numeric(substr(best_combined_matrix, 6, 7))
conclusion2 = paste0("Best similarity matrix based on combined normalized scores is for $nn' = ", 
                     optNN, "$. We therefore proceed with $nn' = ",
                     optNN, "$.")
conclusion = paste(conclusion1, conclusion2)

# M3C k-means clustering
library(M3C)

# Import resources
scheme = readRDS("Resources/scheme.rds")
annCol = scheme$annCol
annColors = scheme$annColors
cluster_colors = scheme$clust.colors
col.list = scheme$col.list
var2comp = scheme$var2comp
rm(scheme); gc()

# Here we create a class column for ER status
m3c_des = annCol
m3c_des$class = m3c_des$`ER status`
m3c_des$ID = rownames(m3c_des)
m3c_input = t(similarity_object[[paste0("NN = ", optNN)]]$F) %>% as.data.frame()
rownames(m3c_input) = c("t_SNE1", "t_SNE2")
colnames(m3c_input) = colnames(input[["SNPs"]])

RNGversion("4.2.2")
consensus_km = M3C(m3c_input, des = m3c_des, iters = 100, repsref = 250, 
                   repsreal = 250, seed = 123, fsize = 18, lthick = 2, dotsize = 1.25,
                   clusteralg = "km", maxK = 2)

# Is the clustering significant?
paste0(ifelse(consensus_km$scores$NORM_P < 0.05, "The clustering is significant.",
              "The clustering is not significant."))

# Normally, we would not proceed, but for the sake of comparisons across algorithms
# we will produce the additional relevant plots

# Main results ###
# Examine cluster similarity to MOVICS by measuring NMI and ARI indices #####
# (Jaccard may be misleading)
PINSPlus_clusters = as.data.frame(list(Sample.ID = rownames(consensus_km[["realdataresults"]][[2]][["ordered_annotation"]]),
                                       Cluster = consensus_km[["realdataresults"]][[2]][["ordered_annotation"]][["consensuscluster"]]))
PINSPlus_clusters$Sample.ID = gsub("\\.", "-", PINSPlus_clusters$Sample.ID)
rownames(PINSPlus_clusters) = PINSPlus_clusters$Sample.ID

# Calculate ARI and NMI
library(mclust)
library(clue)

ARI_to_MOVICS = calculate_ari_index(cluster_df1 = ground_truth_labels,
                                    cluster_df2 = PINSPlus_clusters,
                                    sample_col = "Sample.ID",
                                    clust_col = "Cluster",
                                    suffixes = c("_MOVICS", "_PINSPlus"))

NMI_to_MOVICS = calculate_nmi_index(cluster_df1 = ground_truth_labels,
                                    cluster_df2 = PINSPlus_clusters,
                                    sample_col = "Sample.ID",
                                    clust_col = "Cluster",
                                    suffixes = c("_MOVICS", "_PINSPlus"))

# Very low statistics when compared to the MOVICS. Results differ

# MOVICS-like analysis #####
library(MOVICS)
library(ComplexHeatmap)

plotdata <- lapply(lapply(input, as.matrix), 
                   function(mat) mat[rowSums(mat != 0) > 0, ])
heatmap_plotdata = getStdiz(
  data = plotdata,
  halfwidth = c(NA, 3, 3, 3, 3), # No halfwidth for SNPs
  centerFlag = c(F, F, F, F, F),
  scaleFlag = c(F, F, F, F, F)
)

plot_object = list(clust.res = PINSPlus_clusters %>%
                     dplyr::rename(samID = Sample.ID, clust = Cluster))

# comprehensive heatmap (may take a while)
getMoHeatmap_single_algorithm(algorithm_name = algorithm,
                              data          = heatmap_plotdata,
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
                              fig.path      = paste0(home, "/Results/single_algorithm/PINSPlus"),
                              fig.name      = paste0("default_", algorithm, "_Comprehensive_heatmap"))
dev.off()
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
                                         res.path = paste0(home, "/Results/single_algorithm/PINSPlus/"))

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
                                      fig.path     = paste0(home, "/Results/single_algorithm/PINSPlus"),
                                      res.path     = paste0(home, "/Results/single_algorithm/PINSPlus"))

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
                                                 fig.path = paste0(home, "/Results/single_algorithm/PINSPlus"))

# Agreement with other subtypes ###
subtype_agreement <- compAgree_single_algorithm(algorithm_name = algorithm,
                                                moic.res  = plot_object,
                                                subt2comp = annCol[, c("ER status", "PR status",
                                                                       "HER2 status", "Metastasis", "Stage")],
                                                doPlot    = TRUE,
                                                box.width = 0.2,
                                                fig.name  = "Classification_agreement",
                                                fig.path  = paste0(home, "/Results/single_algorithm/PINSPlus"),
                                                width     = 12)
dev.off()

# DGEA ###
dgea = runDEA(dea.method = "limma", # we use normalized data as input
              expr = plotdata$RNAseq,
              moic.res = plot_object,
              prefix = "dgea_",
              sort.p = TRUE,
              overwt = TRUE,
              verbose = TRUE,
              res.path = paste0(home, "/Results/single_algorithm/PINSPlus"))

# # Identify unique subtype biomarkers
# # 1. Up-regulated markers
dgea.marker.up <- runMarker_single_algorithm(algorithm_name = algorithm,
                                             moic.res = plot_object,
                                             dea.method    = "limma", # name of DEA method
                                             prefix        = "dgea_", # MUST be the same of argument in runDEA()
                                             dat.path      = paste0(home, "/Results/single_algorithm/PINSPlus"), # path of DEA files
                                             res.path      = paste0(home, "/Results/single_algorithm/PINSPlus"), # path to save marker files
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
                                             fig.path = paste0(home, "/Results/single_algorithm/PINSPlus"),
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
                                               dat.path      = paste0(home, "/Results/single_algorithm/PINSPlus"), # path of DEA files
                                               res.path      = paste0(home, "/Results/single_algorithm/PINSPlus"), # path to save marker files
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
                                               fig.path = paste0(home, "/Results/single_algorithm/PINSPlus"),
                                               width = 14,
                                               height = 12,
                                               fontsize_row = 3,
                                               name = "normalized RNA-seq")
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
                                            dat.path      = paste0(home, "/Results/single_algorithm/PINSPlus"), # path of DEA files
                                            res.path      = paste0(home, "/Results/single_algorithm/PINSPlus"), # path to save marker files
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
                                            fig.path = paste0(home, "/Results/single_algorithm/PINSPlus"),
                                            width = 14, height = 12)

# GSEA down-regulated
RNGversion("4.2.2")
set.seed(123)
gsea.down <- runGSEA_mod_4.4_single_algorithm(algorithm_name = algorithm,
                                              moic.res     = plot_object,
                                              dea.method   = "limma", # name of DEA method
                                              prefix       = "dgea_", # MUST be the same of argument in runDEA()
                                              dat.path      = paste0(home, "/Results/single_algorithm/PINSPlus"), # path of DEA files
                                              res.path      = paste0(home, "/Results/single_algorithm/PINSPlus"), # path to save marker files
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
                                              fig.path = paste0(home, "/Results/single_algorithm/PINSPlus"),
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
                                            fig.path      = paste0(home, "/Results/single_algorithm/PINSPlus"),
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
  fig.path = paste0(home, "/Results/single_algorithm/PINSPlus"),
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
  fig.path = paste0(home, "/Results/single_algorithm/PINSPlus"),
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
                                                 res.path = paste0(home, "/Results/single_algorithm/PINSPlus"))

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
                          subt1.lab = "PINSPlus",
                          subt2.lab = "NTP TCGA",
                          height = 8,
                          width = 8,
                          fig.path = paste0(home, "/Results/single_algorithm/PINSPlus"),
                          fig.name = paste0("kappa_", algorithm, "_vs_NTP_TCGA"))

# consensus TCGA vs PAM TCGA
runKappa_single_algorithm(algorithm_name = algorithm,
                          subt1 = plot_object$clust.res$clust,
                          subt2 = gsub(algorithm, "", TCGA.pam.pred$clust.res$clust),
                          subt1.lab = "PINSPlus",
                          subt2.lab = "PAM TCGA",
                          height = 8,
                          width = 8,
                          fig.path = paste0(home, "/Results/single_algorithm/PINSPlus"),
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
                          fig.path = paste0(home, "/Results/single_algorithm/PINSPlus"),
                          fig.name = "kappa_NTP_vs_PAM_transNEO")

# Export consensus clustering object
clust = as.data.frame(plot_object$clust.res)
colnames(clust) = c("Sample.ID", "Cluster")
clust$Cluster = paste0("PINSPlus", clust$Cluster)
openxlsx::write.xlsx(clust, paste0(home, "/Results/single_algorithm/PINSPlus/", 
                                   algorithm, "_", data_source, "_",
                                   data_types, "_eval_on_", evaluation_source,
                                   "_clusterings.xlsx"))

# Supplementary results #####

# Create subdirectory for supplementary plots
if (!dir.exists(paste0(home, "/Results/single_algorithm/PINSPlus/Supplement"))) {
  dir.create(paste0(home, "/Results/single_algorithm/PINSPlus/Supplement"))
}

# Setup for heatmaps
colors_heatmap = rev(colorRampPalette(viridisLite::magma(10))(255))
cluster_colors_heatmap = c("#2EC4B6", "#E71D36")
clust_annot_pheno = annCol %>% mutate(Sample.ID = rownames(.)) %>%
  inner_join(clust, by = "Sample.ID") %>%
  dplyr::rename(PINSPlus = Cluster, samID = "Sample.ID")
rownames(clust_annot_pheno) = clust_annot_pheno$samID
afh_colnames = colnames(annCol)

# Same data frame. Different columns. Just for easiness
PINSPlus_clust_res = PINSPlus_clusters %>% dplyr::rename(samID = Sample.ID, PINSPlus = Cluster)

# PCA from original matrices ###
# RNA
pca_from_original_matrix(mydata = input$RNAseq, 
                         algorithm = "PINSPlus", 
                         clust_res = PINSPlus_clust_res,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = paste0(home, "/Results/single_algorithm/PINSPlus/Supplement"),
                         title_add = "RNAseq")

# miRNA
pca_from_original_matrix(mydata = input$miRNA, 
                         algorithm = "PINSPlus", 
                         clust_res = PINSPlus_clust_res,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = paste0(home, "/Results/single_algorithm/PINSPlus/Supplement"),
                         title_add = "miRNA")

# CNV
pca_from_original_matrix(mydata = input$CNV, 
                         algorithm = "PINSPlus", 
                         clust_res = PINSPlus_clust_res,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = paste0(home, "/Results/single_algorithm/PINSPlus/Supplement"),
                         title_add = "CNV")

# Use multidimensional scaling for SNPs
# Features must be in rows
mds_from_original_matrix(matrix = input$SNPs, dist_method = "binary",
                         algorithm = "PINSPlus", 
                         clust_res = PINSPlus_clust_res,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = paste0(home, "/Results/single_algorithm/PINSPlus/Supplement"),
                         title_add = "SNPs")

# Methylation
pca_from_original_matrix(mydata = input$Methylation, 
                         algorithm = "PINSPlus", 
                         clust_res = PINSPlus_clust_res,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = paste0(home, "/Results/single_algorithm/PINSPlus/Supplement"),
                         title_add = "Methylation")

# Draw a heatmap of the final S matrix ###
pinsplus_matrix = similarity_object[[paste0("NN = ", optNN)]][["S"]]
dimnames(pinsplus_matrix) = list(colnames(input$SNPs), colnames(input$SNPs))
create_MO_heatmap(matrix = pinsplus_matrix, algorithm = "PINSPlus", 
                  need.diag.zero = FALSE, # already zero
                  clust_annot_pheno = clust_annot_pheno ,
                  afh_colnames = afh_colnames, 
                  colors = colors_heatmap,
                  annColors = annColors,
                  heatmap_title = "Final PINSPlus similarity heatmap",
                  cluster_colors = cluster_colors_heatmap,
                  cluster_rows_flag = FALSE,
                  cluster_cols_flag = FALSE,
                  splits_flag = TRUE,
                  legend_title = "Final kernel similarity",
                  output_file_name = paste0(home, "/Results/single_algorithm/PINSPlus/Supplement/PINSPlus_final_S_matrix_heatmap.png"))

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

voi = colnames(clust_annot_pheno)[1:11]
output = as.data.frame(matrix(NA, nrow = 0, ncol = 4))
for (v in 1:length(voi)){
  test = suppressWarnings(chisq.test(table(clust_annot_pheno[, algorithm], 
                                           clust_annot_pheno[, voi[v]])))
  chifit_p = test$p.value
  chifit_xsq = test$statistic
  chifit_cv = suppressWarnings(unbiased.cv.test(table(clust_annot_pheno[, algorithm], 
                                                      clust_annot_pheno[, voi[v]]),
                                                string = voi[v],
                                                digits = 3)$value)
  comparison = paste0(voi[v], " vs ", algorithm, " cluster")
  output = rbind(output, c(comparison, chifit_p, chifit_xsq, chifit_cv))
  rm(test, comparison, chifit_p, chifit_xsq, chifit_cv)
}
colnames(output) = c("Comparison", "p-value", "Statistic", "Cramer's V")

rm(v); gc()
openxlsx::write.xlsx(output, 
                     paste0(home, 
                            "/Results/single_algorithm/PINSPlus/Supplement/Chisq_tests.xlsx"),
                     overwrite = TRUE)

# Bar chart generation
PINSPlus_barcharts = list()
plotdata_bar = clust_annot_pheno
plotdata_bar[[algorithm]] = factor(plotdata_bar[[algorithm]])
for (i in 1:length(voi)) {
  chifit = output
  loc = which(grepl(voi[i], chifit$Comparison))
  chifit = chifit[loc, ]
  PINSPlus_barcharts[[i]] = create_annot_barchart(plotdata = plotdata_bar, fill = voi[i],
                                                  chifit = chifit,
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
  print(PINSPlus_barcharts[[i]])
  ggsave(filename = paste0(algorithm, "_", voi[i], "_barchart.png"),
         path = paste0(home, 
                       "/Results/single_algorithm/PINSPlus/Supplement"), 
         width = 2320, height = 2320, device = 'png', units = "px",
         dpi = 700)
  dev.off()
}
names(PINSPlus_barcharts) = voi
rm(loc, chifit)

# Multiplot (PNG) - bar charts
library(ggpubr)
ggarrange(PINSPlus_barcharts[[1]], PINSPlus_barcharts[[2]], PINSPlus_barcharts[[3]],
          PINSPlus_barcharts[[4]], PINSPlus_barcharts[[5]], PINSPlus_barcharts[[6]],
          PINSPlus_barcharts[[7]], PINSPlus_barcharts[[8]], PINSPlus_barcharts[[9]],
          PINSPlus_barcharts[[10]], PINSPlus_barcharts[[11]],
          ncol = 3, nrow = 4, labels = c("A", "B", "C", "D", "E", "F", "G", "H",
                                         "I", "J", "K"),
          font.label = list(size = 8, face = "bold", color ="black"))
ggsave(filename = paste0("Multiplot_", algorithm, "_barcharts.png"),
       path = paste0(home, 
                     "/Results/single_algorithm/PINSPlus/Supplement"), 
       width = 6500, height = 8000, device = 'png', units = "px",
       dpi = 700)
dev.off()

# Just significant ones now
PINSPlus_barcharts_sig = list()
plotdata_bar_sig = clust_annot_pheno %>% dplyr::select(PINSPlus, Histology, 
                                                       `ER status`, `PR status`)
plotdata_bar_sig$PINSPlus = factor(plotdata_bar_sig$PINSPlus)
voi_sig = setdiff(colnames(plotdata_bar_sig), "PINSPlus")
for (i in 1:length(voi_sig)) {
  chifit = output
  loc = which(grepl(voi_sig[i], chifit$Comparison))
  chifit = chifit[loc, ]
  PINSPlus_barcharts_sig[[i]] = create_annot_barchart(plotdata = plotdata_bar_sig, fill = voi_sig[i],
                                                      chifit = chifit,
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
  print(PINSPlus_barcharts_sig[[i]])
  ggsave(filename = paste0("sig_PINSPlus_", voi_sig[i], "_barchart.png"),
         path = paste0(home, 
                       "/Results/single_algorithm/PINSPlus/Supplement"), 
         width = 2320, height = 2320, device = 'png', units = "px",
         dpi = 700)
  dev.off()
}
names(PINSPlus_barcharts_sig) = voi_sig
rm(loc, chifit)

# Multiplot (PNG) - bar charts
ggarrange(PINSPlus_barcharts_sig[[1]], PINSPlus_barcharts_sig[[2]], PINSPlus_barcharts_sig[[3]],
          ncol = 1, nrow = 3, labels = c("A", "B", "C"),
          font.label = list(size = 8, face = "bold", color ="black"))
ggsave(filename = "sig_Multiplot_PINSPlus_barcharts.png",
       path = paste0(home, 
                     "/Results/single_algorithm/PINSPlus/Supplement"), 
       width = 2500, height = 7000, device = 'png', units = "px",
       dpi = 700)
dev.off()

# Sunburst plot ###
library(plotly)
Pheno_sunburst_PINSPlus = clust_annot_pheno
Pheno_sunburst_PINSPlus$`ER status` = gsub("Unknown", "Unkn ER status", Pheno_sunburst_PINSPlus$`ER status`)
Pheno_sunburst_PINSPlus$`ER status` = gsub("Positive", "ER+", Pheno_sunburst_PINSPlus$`ER status`)
Pheno_sunburst_PINSPlus$`ER status` = gsub("Negative", "ER-", Pheno_sunburst_PINSPlus$`ER status`)
Pheno_sunburst_PINSPlus = Pheno_sunburst_PINSPlus %>%
  dplyr::select(PINSPlus, `ER status`) %>%
  group_by(PINSPlus, `ER status`) %>%
  summarise(Counts = n()) %>%
  as.data.frame()

sunburst_coloring_PINSPlus = data.frame(stringsAsFactors = FALSE,
                                        colors = tolower(gplots::col2hex(c("#2EC4B6", "#E71D36", 
                                                                           "#C11D9C", "#0F1682",  "grey40"))),
                                        labels = c("PINSPlus1", "PINSPlus2",
                                                   "ER-", "ER+", "Unkn ER status"))

sunburstDF_PINSPlus = as.sunburstDF(Pheno_sunburst_PINSPlus, value_column = "Counts", add_root = FALSE) %>%
  inner_join(sunburst_coloring_PINSPlus, by = "labels")

pie_PINSPlus = plot_ly() %>%
  add_trace(ids = sunburstDF_PINSPlus$ids, labels= sunburstDF_PINSPlus$labels, 
            parents = sunburstDF_PINSPlus$parents, 
            values= sunburstDF_PINSPlus$values, type='sunburst', branchvalues = 'total',
            insidetextorientation='radial', maxdepth = 5,
            marker = list(colors = sunburstDF_PINSPlus$colors)) %>%
  layout(
    grid = list(columns =1, rows = 1),
    margin = list(l = 0, r = 0, b = 0, t = 0)
  )
pie_PINSPlus
rm(Pheno_sunburst_PINSPlus, sunburstDF_PINSPlus, sunburst_coloring_PINSPlus, pie_PINSPlus); gc()

# Compare these PINSPlus results with the PINSPlus output from MOVICS ###
load("Results/MOVICS_baseline/MOVICS_TCGA_RNAseq-CNV-Methylation-miRNA-SNPs_eval_on_transNEO_moic.res.list.rda")
MOVICS_PINSPlus = moic.res.list$PINSPlus$clust.res

ARI_to_MOVICS_PINSPlus = calculate_ari_index(cluster_df1 = MOVICS_PINSPlus %>%
                                               dplyr::rename(Sample.ID = samID,
                                                             Cluster = clust),
                                             cluster_df2 = PINSPlus_clusters,
                                             sample_col = "Sample.ID",
                                             clust_col = "Cluster",
                                             suffixes = c("_MOVICS_PINSPlus", "_PINSPlus"))

NMI_to_MOVICS_PINSPlus = calculate_nmi_index(cluster_df1 = MOVICS_PINSPlus %>%
                                               dplyr::rename(Sample.ID = samID,
                                                             Cluster = clust),
                                             cluster_df2 = PINSPlus_clusters,
                                             sample_col = "Sample.ID",
                                             clust_col = "Cluster",
                                             suffixes = c("_MOVICS_PINSPlus", "_PINSPlus"))

# Wrap up #####
hyperparameters = list(num_neighbors_min = min(num_neighbors_range),
                       num_neighbors_max = max(num_neighbors_range),
                       num_neighbors_step = neighbor_step,
                       optimal_NN = optNN,
                       conclusion1 = conclusion1,
                       conclusion2 = conclusion2,
                       conclusion = conclusion
)

# Put all parameters in a list
params = list(algorithm = algorithm, data_source = data_source, data_types = data_types,
              evaluation_source = evaluation_source, title = title, subtitle = subtitle,
              description = description, in_a_nutshell = in_a_nutshell, optk_text = optk_text,
              citation = citation, NMI_to_MOVICS = NMI_to_MOVICS, ARI_to_MOVICS = ARI_to_MOVICS,
              NMI_to_MOVICS_PINSPlus = NMI_to_MOVICS_PINSPlus, ARI_to_MOVICS_PINSPlus = ARI_to_MOVICS_PINSPlus,
              hyperparameters = hyperparameters, ground_truth_k = ground_truth_k,
              sessionInfo = sessionInfo(), home = home)

# Render the R Markdown document with the parameters
rmarkdown::render(paste0(getwd(), "/Results/single_algorithm/PINSPlus/PINSPlus_report.Rmd"), 
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
