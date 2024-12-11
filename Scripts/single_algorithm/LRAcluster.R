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
algorithm = "LRAcluster"
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

# Extract the names of the modalities that will be used
modalities = unlist(strsplit(data_types, "-"))

# Replace with TRUE wherever features are in rows
features_in_rows = rep(TRUE, length(modalities))

# Run algorithm #####
# Use the LRAcluster functions from the downloaded LRAcluster source:
# http://lifeome.net/software/lracluster
source("Resources/full_LRAcluster_source.R")

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

# Export input for HPC
saveRDS(input, "Resources/LRAcluster_input.rds")

# Setup ###
# Hyperparameter tuning

# We ran the Scripts/single_algorithm/LRAcluster_HPC.R script in an HPC
# with SLURM parameters specified in Scripts/single_algorithm/LRAcluster_HPC.sh

# We import the results here for further procressing
dimensions_to_try = c(2:10)
results_list = readRDS("Resources/HPC output/LRAcluster_HPC/LRAcluster_results_list.rds")

# Extract potential values for plotting explained variance
potential_values <- sapply(results_list, function(x) x$potential)

# Plot the explained variance (potential) vs. dimension
library(ggplot2)
ggplot(data = data.frame(Dimension = dimensions_to_try, ExplainedVariation = potential_values),
       aes(x = Dimension, y = ExplainedVariation)) +
  geom_line(color = "#f5bc83", size = 0.5) +
  geom_point(color = "#880a49", size = 1, shape = 16) +
  labs(
    title = "Explained Variation by Dimension",
    x = "Dimension",
    y = "Explained Variation"
  ) +
  theme(
    plot.background = element_rect(fill = "white", color = "white"),
    panel.background = element_rect(fill = "white", color = "white"),
    panel.grid.major.y = element_line(color = "lightgray", linetype = "dotted", linewidth = 0.13),
    panel.grid.major.x = element_blank(),
    axis.title = element_text(face = "bold", size = 4.5),
    axis.text = element_text(size = 3),
    axis.line = element_line(linetype = "solid", linewidth = 0.1),
    axis.ticks = element_line(linewidth = 0.05),
    plot.title = element_text(face = "bold", size = 5, hjust = 0.5)
  )+
  scale_x_continuous(limits = c(0, 11), breaks = 2:10, expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 0.4), breaks = seq(0, 0.35, 0.05), expand = c(0, 0)) +
  geom_vline(xintercept = 2:10, color = "grey", linetype = "dotted", linewidth = 0.13)
ggsave("Results/single_algorithm/LRAcluster/explained_variation_plot.png",
       dpi = 700, width = 1920, height = 1080, units = "px")
dev.off()

# We pick 7 as the optimal number of low dimensions due to the high jump of the 
# line for dimensions = 7 and the slow increase afterwards

# The authors suggest running k-means clustering (or other unsupervised methods)
# for a varying number of clusters and choose optimal k based on silhouette
optr = 7

# # Hierarchical clustering
# library(fastcluster)
# dist_mat = dist(t(results_list[["dim_7"]][["coordinate"]]))
# hclust_output = fastcluster::hclust(dist_mat,
#                                     method = "complete")
# 
# # Get silhouette values
# possible_K = c(2:10)
# avg_sil_values <- numeric(length(possible_K))
# 
# for (i in seq_along(possible_K)) {
#   k <- possible_K[i]
#   cluster_assignments <- cutree(hclust_output, k = k)
#   sil <- cluster::silhouette(cluster_assignments, 
#                              dist_mat)
#   avg_sil_values[i] <- mean(sil[, "sil_width"])
# }
# sil_results <- data.frame(k = possible_K, avg_sil = avg_sil_values)
# 
# # Identify k with the highest average silhouette and choose as optimal
# best_k <- possible_K[which.max(avg_sil_values)] # 4: close second
# 
# # Average Silhouette plot
# ggplot(data = data.frame(cluster.no = c(2:10), Avg.sil = avg_sil_values),
#        aes(x = cluster.no, y = avg_sil_values)) +
#   geom_line(color = "#f5bc83", size = 0.5) +
#   geom_point(color = "#880a49", size = 1, shape = 16) +
#   labs(
#     title = "Average silhouette index per number of clusters",
#     x = "Number of clusters",
#     y = "Average Silhouette index"
#   ) +
#   theme(
#     plot.background = element_rect(fill = "white", color = "white"),
#     panel.background = element_rect(fill = "white", color = "white"),
#     panel.grid.major.y = element_line(color = "lightgray", linetype = "dotted", linewidth = 0.13),
#     panel.grid.major.x = element_blank(),
#     axis.title = element_text(face = "bold", size = 4.5),
#     axis.text = element_text(size = 3),
#     axis.line = element_line(linetype = "solid", linewidth = 0.1),
#     axis.ticks = element_line(linewidth = 0.05),
#     plot.title = element_text(face = "bold", size = 5, hjust = 0.5)
#   )+
#   scale_x_continuous(limits = c(0, 11), breaks = 2:10, expand = c(0, 0)) +
#   scale_y_continuous(limits = c(0, 0.6), breaks = seq(0, 0.5, 0.1), expand = c(0, 0)) +
#   geom_vline(xintercept = 2:10, color = "grey", linetype = "dotted", linewidth = 0.13)
# ggsave("Results/single_algorithm/LRAcluster/avg_silhouette_plot.png",
#        dpi = 700, width = 1920, height = 1080, units = "px")
# dev.off()

# Import resources
scheme = readRDS("Resources/scheme.rds")
annCol = scheme$annCol
annColors = scheme$annColors
cluster_colors = scheme$clust.colors
col.list = scheme$col.list
var2comp = scheme$var2comp
rm(scheme); gc()

# M3C
library(M3C)
# Here we create a class column for ER status
m3c_des = annCol
m3c_des$class = m3c_des$`ER status`
m3c_des$ID = rownames(m3c_des)
m3c_input = results_list[[paste0("dim_", optr)]][["coordinate"]] %>% as.data.frame()
rownames(m3c_input) = paste0("LRA_", rownames(results_list[["dim_7"]][["coordinate"]]))
colnames(m3c_input) = colnames(results_list[["dim_7"]][["coordinate"]])

RNGversion("4.2.2")
consensus_km = M3C(m3c_input, des = m3c_des, iters = 100, repsref = 250, 
                   repsreal = 250, seed = 123, fsize = 18, lthick = 2, dotsize = 1.25,
                   clusteralg = "km", maxK = 10) # optimal K: 4

optk = 4 # p = 0.046
paste0(ifelse(consensus_km$scores$NORM_P[consensus_km$scores$K == 4] < 0.05, "The clustering is significant.",
              "The clustering is not significant."))

# Plotting clustering info
# Consensus index plot
ci_plot = ggplot(consensus_km[["plots"]][[1]][["data"]], aes(x = consensusindex, y = CDF,
                                                     group = k, alpha = 0.7))+
  geom_line(aes(color = factor(k)), linewidth = 0.5)+
  theme_bw()+
  theme(panel.border = element_rect(linewidth = 0.2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 5, face = "bold"),
        legend.title = element_text(face = "bold", size = 4),
        legend.text = element_text(size = 3),
        legend.key.size = unit(0.2, "cm"),
        legend.margin = ggplot2::margin(0, 0, 0, 0, unit = "mm"),
        legend.spacing.y = unit(0.5, units = "mm"),
        axis.title.x = element_text(size = 4, face = "bold"),
        axis.title.y = element_text(size = 4, face = "bold"),
        axis.ticks = element_line(linewidth = 0.15),
        axis.text.x = element_text(size = 4),
        axis.text.y = element_text(size = 4))+
  labs(y = "Cumulative Distribution Function (CDF)",
       x = "Consensus Index",
       title = "Real Data")+
  guides(color = guide_legend(title = "K"), linewidth = "none", alpha = "none")
ci_plot
ggsave(filename = "Consensus_index.png",
       path = "Results/single_algorithm/LRAcluster", 
       width = 1920, height = 1080, device = 'png', units = "px",
       dpi = 700)
dev.off()

# Entropy plot
entropy = ggplot(consensus_km[["plots"]][[2]][["data"]], aes(x = K, y = PAC_SCORE, alpha = 0.7))+
  geom_line(aes(color = "#7c1d6f"), linewidth = 0.5)+
  scale_x_continuous(limits = c(2, 10), breaks = seq(2, 10, 1))+
  theme_bw()+
  theme(panel.border = element_rect(linewidth = 0.2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 5, face = "bold", vjust = 0.5, hjust = 0.5),
        legend.title = element_text(face = "bold", size = 4),
        legend.text = element_text(size = 3),
        legend.key.size = unit(0.2, "cm"),
        legend.margin = ggplot2::margin(0, 0, 0, 0, unit = "mm"),
        legend.spacing.y = unit(0.5, units = "mm"),
        axis.title.x = element_text(size = 4, face = "bold"),
        axis.title.y = element_text(size = 4, face = "bold"),
        axis.ticks = element_line(linewidth = 0.15),
        axis.text.x = element_text(size = 4),
        axis.text.y = element_text(size = 4))+
  labs(y = "Entropy",
       x = "K",
       title = "Real Data")+
  guides(color = "none", alpha = "none")
entropy
ggsave(filename = "Entropy.png",
       path = "Results/single_algorithm/LRAcluster", 
       width = 1920, height = 1080, device = 'png', units = "px",
       dpi = 700)
dev.off()

# Statistical significance of clusters
statsig_clust = ggplot(consensus_km[["plots"]][[3]][["data"]], aes(x = K, y = P_SCORE, 
                                                           color = P_SCORE < -log10(0.05)))+
  geom_point(size = 1.5, alpha = 0.6)+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.2)+
  scale_x_continuous(limits = c(2, 10), breaks = seq(2, 10, 1))+
  scale_color_manual(name = "Color",
                     values = c("#6c2167", "grey"),
                     labels = c("p < 0.05", "p > 0.05")) +
  theme_bw()+
  theme(panel.border = element_rect(linewidth = 0.2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 5, face = "bold"),
        axis.title.x = element_text(size = 4, face = "bold"),
        axis.title.y = element_text(size = 4, face = "bold"),
        axis.ticks = element_line(linewidth = 0.15),
        axis.text.x = element_text(size = 4),
        axis.text.y = element_text(size = 4),
        legend.title = element_text(face = "bold", size = 4),
        legend.text = element_text(size = 3),
        legend.key.size = unit(0.2, "cm"),
        legend.margin = ggplot2::margin(0, 0, 0, 0, unit = "mm"),
        legend.spacing.y = unit(0.5, units = "mm"))+
  labs(title = "Statistical significance of different values of K",
       y = bquote(bold(-log[10]("p"))))
statsig_clust
ggsave(filename = "Stat_sig.png",
       path = "Results/single_algorithm/LRAcluster", 
       width = 1920, height = 1080, device = 'png', units = "px",
       dpi = 700)
dev.off()

# RCSI plot
rcsi = ggplot(as.data.frame(consensus_km[["scores"]]), aes(x = consensus_km$scores$K,
                                                   y = consensus_km$scores$RCSI))+
  geom_line(size = 0.3, color = "violet")+
  geom_errorbar(aes(ymin = consensus_km$scores$RCSI - consensus_km$scores$RCSI_SE,
                    ymax = consensus_km$scores$RCSI + consensus_km$scores$RCSI_SE,
                    color = "deeppink3"), width = 0.2, size = 0.1)+
  geom_point(size = 0.05, color ="deeppink3")+
  scale_x_continuous(limits = c(1.9, 10.1), breaks = seq(2, 10, 1))+
  scale_y_continuous(limits = c(-0.6, 0.6), breaks = seq(-0.5, 0.5, 0.1))+
  theme(plot.title = element_text(size = 5, face = "bold"),
        axis.title.x = element_text(size = 4, face = "bold"),
        axis.title.y = element_text(size = 4, face = "bold"),
        axis.ticks = element_line(linewidth = 0.15),
        axis.text.x = element_text(size = 4),
        axis.text.y = element_text(size = 4),
        legend.position = "none",
        panel.background = element_rect(fill = "white", 
                                        colour = "white"),
        panel.grid = element_blank(),
        axis.line = element_line(linewidth = 0.2))+
  labs(title = "Relative Cluster Stability Index vs. number of clusters K",
       x = "K", y = "RCSI")
rcsi
ggsave(filename = "RCSI.png",
       path = "Results/single_algorithm/LRAcluster", 
       width = 1920, height = 1080, device = 'png', units = "px",
       dpi = 700)
dev.off()

# Main results ###
# Examine cluster similarity to MOVICS by measuring NMI and ARI indices #####
# (Jaccard may be misleading)
LRAcluster_clusters = as.data.frame(list(Sample.ID = rownames(consensus_km[["realdataresults"]][[4]][["ordered_annotation"]]),
                                    Cluster = consensus_km[["realdataresults"]][[4]][["ordered_annotation"]][["consensuscluster"]]))
LRAcluster_clusters$Sample.ID = gsub("\\.", "-", LRAcluster_clusters$Sample.ID)
rownames(LRAcluster_clusters) = LRAcluster_clusters$Sample.ID

# Calculate ARI and NMI
library(mclust)
library(clue)

ARI_to_MOVICS = calculate_ari_index(cluster_df1 = ground_truth_labels,
                                    cluster_df2 = LRAcluster_clusters,
                                    sample_col = "Sample.ID",
                                    clust_col = "Cluster",
                                    suffixes = c("_MOVICS_CS",
                                                 paste0("_", algorithm)))

NMI_to_MOVICS = calculate_nmi_index(cluster_df1 = ground_truth_labels,
                                    cluster_df2 = LRAcluster_clusters,
                                    sample_col = "Sample.ID",
                                    clust_col = "Cluster",
                                    suffixes = c("_MOVICS_CS",
                                                 paste0("_", algorithm)))

# Very low statistics when compared to the MOVICS. Results differ

# MOVICS-like analysis #####
library(MOVICS)
library(ComplexHeatmap)

# Import coloring scheme
scheme = readRDS("Resources/scheme.rds")
annCol = scheme$annCol
annColors = scheme$annColors
cluster_colors = c("#2EC4B6", "#E71D36", "#FF9F1C", "#BDD5EA")
col.list = scheme$col.list
var2comp = scheme$var2comp %>%
  dplyr::select(-`Consensus Subtype`) %>%
  mutate(Sample.ID = rownames(.)) %>%
  inner_join(LRAcluster_clusters, by = "Sample.ID") %>%
  tibble::column_to_rownames(var = "Sample.ID") %>%
  mutate(LRAcluster = paste0(algorithm, Cluster)) %>%
  dplyr::select(LRAcluster, everything()) %>%
  dplyr::select(-Cluster)
rm(scheme); gc()

# Silhouette

# transformed_aff = transform_affinity_matrix(final_affinity_matrix,
#                                             threshold = 100, norm_quant = 0,
#                                             norm_method = "divide by quantile")

cor_mat = cor(as.matrix(results_list[[paste0("dim_", optr)]][["coordinate"]]), method = "spearman")
sil = compute_silhouette(cluster_df = LRAcluster_clusters %>% dplyr::rename(samID = Sample.ID),
                         similarity_matrix = cor_mat,
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
                   function(mat) mat[rowSums(mat != 0) > 0, ])
plotdata = getStdiz(
  data = plotdata,
  halfwidth = c(NA, 3, 3, 3, 3), # No halfwidth for SNPs
  centerFlag = c(F, F, F, F, F),
  scaleFlag = c(F, F, F, F, F)
)

plot_object = list(clust.res = LRAcluster_clusters %>%
                     dplyr::rename(samID = Sample.ID, clust = Cluster))

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
                                         pdf_count_col_width = "5em",
                                         pdf_pval_col_width = "3em",
                                         pdf_test_col_width = "5em",
                                         pdf_tab_font_size = 9)

clin_ordinal_comp = compClinvar_ordinal_single_algorithm(algorithm_name = algorithm,
                                                         moic.res = plot_object,
                                                         var2comp = var2comp_nonas %>%
                                                           dplyr::select(number_of_lymphnodes_positive_by_ihc,
                                                                         number_of_lymphnodes_positive_by_he,
                                                                         LRAcluster),
                                                         strata = algorithm,
                                                         ordinalVars = c("number_of_lymphnodes_positive_by_ihc",
                                                                         "number_of_lymphnodes_positive_by_he"),
                                                         includeNA = FALSE,
                                                         tab.name = "Summary of ordinal clinical variables",
                                                         res.path = paste0(home, "/Results/single_algorithm/", algorithm, "/"),
                                                         output_pdf = TRUE,
                                                         pdf_template_loc = paste0(home, "/Scripts/automated_scripts/clincomp_template.Rmd"),
                                                         pdf_level_col_width = c("7em", "10em"),
                                                         pdf_count_col_width = "5em",
                                                         pdf_pval_col_width = "3em",
                                                         pdf_test_col_width = "5em",
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
                                            width = 14, height = 18)

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
                                              width = 14, height = 18)

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

timestamp() # ~5.5 mins
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
                                                present_clusters = c("LRAcluster1", "LRAcluster2",
                                                                     "LRAcluster3", "LRAcluster4"),
                                                representative = TRUE, moic.res = plot_object,
                                                subtype_prefix = algorithm, n.path = 20, msigdb.path = MSIGDB.FILE,
                                                norm.method = "mean", dirct = "up",
                                                fig.name = "upregulated_pathway_heatmap",
                                                name = "GSVA scores",
                                                fig.path = paste0(home, "/Results/single_algorithm/", algorithm), 
                                                width = 15, height = 18, gsva.method = "gsva")

hclust_pathway_plots_down = plot_pathway_heatmaps(gsea.lists = hclust_output[grepl("down", names(hclust_output))], 
                                                  norm.expr = plotdata$RNAseq, 
                                                  present_clusters = c("LRAcluster1", "LRAcluster2",
                                                                       "LRAcluster3", "LRAcluster4"),
                                                  representative = TRUE, moic.res = plot_object,
                                                  subtype_prefix = algorithm, n.path = 20, msigdb.path = MSIGDB.FILE,
                                                  norm.method = "mean", dirct = "down",
                                                  fig.name = "downregulated_pathway_heatmap",
                                                  name = "GSVA scores",
                                                  fig.path = paste0(home, "/Results/single_algorithm/", algorithm), 
                                                  width = 15, height = 18, gsva.method = "gsva")

# Fraction Genome Altered ###
fga_df = readRDS("Resources/TCGA/fga_df.rds"); gc()

fga.LRAcluster <- compFGA_optimized(moic.res     = plot_object,
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

fga.LRAcluster.COSMIC <- compFGA_optimized(moic.res     = plot_object,
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
timestamp() # 11 min

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
  inner_join(expr_conc %>% dplyr::select(Donor.ID = samID, LRAcluster = clust_up),
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
                                                 pdf_count_col_width = "5em",
                                                 pdf_pval_col_width = "3em",
                                                 pdf_test_col_width = "5em",
                                                 pdf_tab_font_size = 9)

transNEO_ntp_expr_up_ord = transNEO_ntp_expr_up
transNEO_ntp_expr_up_ord$clust.res$clust = gsub(algorithm, "", transNEO_ntp_expr_up_ord$clust.res$clust)
transNEO_ordinal_clincomp = compClinvar_ordinal_single_algorithm(algorithm_name = algorithm,
                                                                 moic.res = transNEO_ntp_expr_up_ord,
                                                                 var2comp = transNEO_var2comp_nonas %>%
                                                                   dplyr::select(Grade.pre.NAT, 
                                                                                 Chemo.cycles, 
                                                                                 aHER2.cycles,
                                                                                 LRAcluster),
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
                                                                 pdf_count_col_width = "5em",
                                                                 pdf_pval_col_width = "3em",
                                                                 pdf_test_col_width = "5em",
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

# Export consensus clustering object
clust = as.data.frame(plot_object$clust.res)
colnames(clust) = c("Sample.ID", "Cluster")
clust$Cluster = paste0(algorithm, clust$Cluster)
openxlsx::write.xlsx(clust, paste0(home, "/Results/single_algorithm/", algorithm, "/", 
                                   algorithm, "_", data_source, "_",
                                   data_types, "_eval_on_", evaluation_source,
                                   "_clusterings.xlsx"))

# Supplementary results #####

# Create subdirectory for supplementary plots
if (!dir.exists(paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement"))) {
  dir.create(paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement"))
}

# Setup for heatmaps
colors_heatmap = rev(colorRampPalette(viridisLite::magma(10))(255))
cluster_colors_heatmap = c("#2EC4B6", "#E71D36", "#FF9F1C", "#BDD5EA")
clust_annot_pheno = annCol %>% mutate(Sample.ID = rownames(.)) %>%
  inner_join(clust, by = "Sample.ID") %>%
  dplyr::rename(LRAcluster = Cluster, samID = "Sample.ID")
rownames(clust_annot_pheno) = clust_annot_pheno$samID
afh_colnames = colnames(annCol)

# Same data frame. Different columns. Just for easiness
LRAcluster_clust_res = LRAcluster_clusters %>% dplyr::rename(samID = Sample.ID, LRAcluster = Cluster)

# PCA ###
# RNA
pca_from_original_matrix(mydata = plotdata$RNAseq, 
                         algorithm = algorithm, 
                         clust_res = LRAcluster_clust_res,
                         cluster_colors = c("#2EC4B6", "#E71D36", "#FF9F1C", "#BDD5EA"), 
                         output_path = paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement"),
                         title_add = "RNAseq")

# miRNA
pca_from_original_matrix(mydata = plotdata$miRNA, 
                         algorithm = algorithm, 
                         clust_res = LRAcluster_clust_res,
                         cluster_colors = c("#2EC4B6", "#E71D36", "#FF9F1C", "#BDD5EA"), 
                         output_path = paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement"),
                         title_add = "miRNA")

# CNV
pca_from_original_matrix(mydata = plotdata$CNV, 
                         algorithm = algorithm, 
                         clust_res = LRAcluster_clust_res,
                         cluster_colors = c("#2EC4B6", "#E71D36", "#FF9F1C", "#BDD5EA"), 
                         output_path = paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement"),
                         title_add = "CNV")

# Use multidimensional scaling for SNPs
# Features must be in rows
mds_from_original_matrix(matrix = plotdata$SNPs, dist_method = "binary",
                         algorithm = algorithm, 
                         clust_res = LRAcluster_clust_res,
                         cluster_colors = c("#2EC4B6", "#E71D36", "#FF9F1C", "#BDD5EA"), 
                         output_path = paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement"),
                         title_add = "SNPs")

# Methylation
pca_from_original_matrix(mydata = plotdata$Methylation, 
                         algorithm = algorithm, 
                         clust_res = LRAcluster_clust_res,
                         cluster_colors = c("#2EC4B6", "#E71D36", "#FF9F1C", "#BDD5EA"),
                         output_path = paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement"),
                         title_add = "Methylation")

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
LRAcluster_barcharts = list()
plotdata_bar = clust_annot_pheno_nonas
plotdata_bar[[algorithm]] = factor(plotdata_bar[[algorithm]])
for (i in 1:length(voi)) {
  chifit = output
  loc = which(grepl(voi[i], chifit$Comparison))
  chifit = chifit[loc, ]
  LRAcluster_barcharts[[i]] = create_annot_barchart(plotdata = plotdata_bar, fill = voi[i],
                                             chifit = chifit,
                                             na.action = "na.omit",
                                             algorithm = algorithm,
                                             barchart_ylim = 600,
                                             text_y = 550, rect_ymin = 450,
                                             rect_ymax = 570, x_annot = 2.5,
                                             v_gap = 35, rect_xmin = 2,
                                             rect_xmax = 3, 
                                             annot_text_size = 2.25,
                                             legend.text.size = 5,
                                             x.axis.text.size = 5) +
    barchart_scales[[voi[i]]]
  print(LRAcluster_barcharts[[i]])
  ggsave(filename = paste0(algorithm, "_", voi[i], "_barchart.png"),
         path = paste0(home, 
                       "/Results/single_algorithm/", algorithm, "/Supplement"), 
         width = 4320, height = 2320, device = 'png', units = "px",
         dpi = 700)
  dev.off()
}
names(LRAcluster_barcharts) = voi
rm(loc, chifit)

# Multiplot (PNG) - bar charts
library(ggpubr)
ggarrange(LRAcluster_barcharts[[1]], LRAcluster_barcharts[[2]], LRAcluster_barcharts[[3]],
          LRAcluster_barcharts[[4]], LRAcluster_barcharts[[5]], LRAcluster_barcharts[[6]],
          LRAcluster_barcharts[[7]], LRAcluster_barcharts[[8]], LRAcluster_barcharts[[9]],
          LRAcluster_barcharts[[10]], LRAcluster_barcharts[[11]],
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
LRAcluster_barcharts_sig = list()
plotdata_bar_sig = clust_annot_pheno_nonas %>% dplyr::select(LRAcluster, Race, Histology, 
                                                             `ER status`, `PR status`, 
                                                             `HER2 status`, Stage)
plotdata_bar_sig$LRAcluster = factor(plotdata_bar_sig$LRAcluster)
voi_sig = setdiff(colnames(plotdata_bar_sig), algorithm)
for (i in 1:length(voi_sig)) {
  chifit = output
  loc = which(grepl(voi_sig[i], chifit$Comparison))
  chifit = chifit[loc, ]
  LRAcluster_barcharts_sig[[i]] = create_annot_barchart(plotdata = plotdata_bar_sig, fill = voi_sig[i],
                                                 chifit = chifit,
                                                 na.action = "na.omit",
                                                 algorithm = algorithm,
                                                 barchart_ylim = 600,
                                                 text_y = 550, rect_ymin = 450,
                                                 rect_ymax = 570, x_annot = 2.5,
                                                 v_gap = 35, rect_xmin = 2,
                                                 rect_xmax = 3, 
                                                 annot_text_size = 2.25,
                                                 legend.text.size = 5,
                                                 x.axis.text.size = 5) +
    barchart_scales[[voi_sig[i]]]
  print(LRAcluster_barcharts_sig[[i]])
  ggsave(filename = paste0("sig_", algorithm, "_", voi_sig[i], "_barchart.png"),
         path = paste0(home, 
                       "/Results/single_algorithm/", algorithm, "/Supplement"), 
         width = 4320, height = 2320, device = 'png', units = "px",
         dpi = 700)
  dev.off()
}
names(LRAcluster_barcharts_sig) = voi_sig
rm(loc, chifit)

# Multiplot (PNG) - bar charts
ggarrange(LRAcluster_barcharts_sig[[1]], LRAcluster_barcharts_sig[[2]], LRAcluster_barcharts_sig[[3]],
          LRAcluster_barcharts_sig[[4]], LRAcluster_barcharts_sig[[5]], LRAcluster_barcharts_sig[[6]],
          ncol = 2, nrow = 2, labels = c("A", "B", "C", "D", "E", "F"),
          font.label = list(size = 8, face = "bold", color ="black"))
ggsave(filename = paste0("sig_Multiplot_", algorithm, "_barcharts.png"),
       path = paste0(home, 
                     "/Results/single_algorithm/", algorithm, "/Supplement"), 
       width = 4320, height = 2320*3, device = 'png', units = "px",
       dpi = 700)
dev.off()

# Sunburst plot ###
library(plotly)
Pheno_sunburst_LRAcluster = clust_annot_pheno
Pheno_sunburst_LRAcluster$`ER status` = gsub("Unknown", "Unkn ER status", Pheno_sunburst_LRAcluster$`ER status`)
Pheno_sunburst_LRAcluster$`ER status` = gsub("Positive", "ER+", Pheno_sunburst_LRAcluster$`ER status`)
Pheno_sunburst_LRAcluster$`ER status` = gsub("Negative", "ER-", Pheno_sunburst_LRAcluster$`ER status`)
Pheno_sunburst_LRAcluster$`HER2 status` = gsub("Unknown", "Unkn HER2 status", 
                                         Pheno_sunburst_LRAcluster$`HER2 status`)
Pheno_sunburst_LRAcluster$`HER2 status` = gsub("Positive", "HER2+", Pheno_sunburst_LRAcluster$`HER2 status`)
Pheno_sunburst_LRAcluster$`HER2 status` = gsub("Negative", "HER2-", Pheno_sunburst_LRAcluster$`HER2 status`)
Pheno_sunburst_LRAcluster$Stage = gsub("Unkown", "Unkn stage", Pheno_sunburst_LRAcluster$Stage)
Pheno_sunburst_LRAcluster = Pheno_sunburst_LRAcluster %>%
  dplyr::select(LRAcluster, `ER status`, `HER2 status`, Stage) %>%
  group_by(LRAcluster, `ER status`, `HER2 status`, Stage) %>%
  summarise(Counts = n()) %>%
  as.data.frame()

sunburst_coloring_LRAcluster = data.frame(stringsAsFactors = FALSE,
                                   colors = tolower(gplots::col2hex(c("#2EC4B6", "#E71D36", "#FF9F1C", "#BDD5EA", 
                                                                      "#C11D9C", "#0F1682",  "grey40",
                                                                      "#0B9EF8", "#560DA7", "mistyrose1", 
                                                                      "hotpink4", "grey40",
                                                                      "#00C9FF", "#099CF5", 
                                                                      "#097BF5", "#0B5684", 
                                                                      "grey40"))),
                                   labels = c("LRAcluster1", "LRAcluster2",
                                              "LRAcluster3", "LRAcluster4",
                                              "ER-", "ER+", "Unkn ER status",
                                              "HER2-", "HER2+", "Indeterminate",
                                              "Equivocal", "Unkn HER2 status",
                                              "Stage I", "Stage II",
                                              "Stage III", "Stage IV", "Unkn stage"))

sunburstDF_LRAcluster = as.sunburstDF(Pheno_sunburst_LRAcluster, value_column = "Counts", add_root = FALSE) %>%
  inner_join(sunburst_coloring_LRAcluster, by = "labels")

pie_LRAcluster = plot_ly() %>%
  add_trace(ids = sunburstDF_LRAcluster$ids, labels= sunburstDF_LRAcluster$labels, 
            parents = sunburstDF_LRAcluster$parents, 
            values= sunburstDF_LRAcluster$values, type='sunburst', branchvalues = 'total',
            insidetextorientation='radial', maxdepth = 5,
            marker = list(colors = sunburstDF_LRAcluster$colors)) %>%
  layout(
    grid = list(columns =1, rows = 1),
    margin = list(l = 0, r = 0, b = 0, t = 0)
  )
pie_LRAcluster
rm(Pheno_sunburst_LRAcluster, sunburstDF_LRAcluster, sunburst_coloring_LRAcluster, pie_LRAcluster); gc()

# Compare these LRAcluster results with the LRAcluster output from MOVICS ###
load("Results/MOVICS_baseline/MOVICS_TCGA_RNAseq-CNV-Methylation-miRNA-SNPs_eval_on_transNEO_moic.res.list.rda")
MOVICS_LRAcluster = moic.res.list$LRAcluster$clust.res

ARI_to_MOVICS_LRAcluster = calculate_ari_index(cluster_df1 = MOVICS_LRAcluster %>%
                                          dplyr::rename(Sample.ID = samID,
                                                        Cluster = clust),
                                        cluster_df2 = LRAcluster_clusters,
                                        sample_col = "Sample.ID",
                                        clust_col = "Cluster",
                                        suffixes = c(paste0("_MOVICS_", algorithm),
                                                     paste0("_", algorithm)))

NMI_to_MOVICS_LRAcluster = calculate_nmi_index(cluster_df1 = MOVICS_LRAcluster %>%
                                          dplyr::rename(Sample.ID = samID,
                                                        Cluster = clust),
                                        cluster_df2 = LRAcluster_clusters,
                                        sample_col = "Sample.ID",
                                        clust_col = "Cluster",
                                        suffixes = c(paste0("_MOVICS_", algorithm),
                                                     paste0("_", algorithm)))

# Wrap up #####
hyperparameters = list(num_dimensions_min = min(dimensions_to_try),
                       num_dimensions_max = max(dimensions_to_try),
                       optimal_r = optr,
                       optimal_k = optk
)

# Put all parameters in a list
params = list(algorithm = algorithm, data_source = data_source, data_types = data_types,
              evaluation_source = evaluation_source, title = title, subtitle = subtitle,
              description = description, in_a_nutshell = in_a_nutshell, optk_text = optk_text,
              citation = citation, NMI_to_MOVICS = NMI_to_MOVICS, ARI_to_MOVICS = ARI_to_MOVICS,
              NMI_to_MOVICS_LRAcluster = NMI_to_MOVICS_LRAcluster, ARI_to_MOVICS_LRAcluster = ARI_to_MOVICS_LRAcluster,
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
