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
algorithm = "RGCCA"
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

# However RGCCA prefers features in columns so we transpose the matrices.

# Extract the names of the modalities that will be used
modalities = unlist(strsplit(data_types, "-"))

# Replace with TRUE wherever features are in rows
features_in_rows = rep(TRUE, length(modalities))

# Run algorithm #####
library(mixOmics)

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

# Export input for run at the HPC cluster
saveRDS(input, "Resources/RGCCA_input.rds")

# Import the results from the HPC cluster
rgcca_results = readRDS("Resources/HPC output/RGCCA_HPC/rgcca_results.rds")

# The above, are the results of this code:
# RNGversion("4.2.2")
# set.seed(123)
# 
# # Hyperparameter setup ###
# design = 1 - diag(length(input)) # consider all the pairwise relationships
# tau = "optimal"
# ncomp = 10
# scheme = "horst"
# keepX = NULL
# scale = FALSE
# max.iter = 2000
# init = "svd.single"
# near.zero.var = FALSE
# 
# # Run RGCCA
# rgcca = wrapper.rgcca(
#   input,
#   design = design,
#   tau = tau,
#   ncomp = ncomp,
#   keepX = keepX,
#   scheme = scheme,
#   scale = scale,
#   init = init,
#   tol = .Machine$double.eps, # default
#   max.iter = max.iter,
#   near.zero.var = near.zero.var,
#   all.outputs = TRUE
# )

# Inspect output
summary(rgcca_results)

# Check the ranges of the variate matrices
lapply(rgcca_results$variates, mean)
lapply(rgcca_results$variates, sd)
lapply(rgcca_results$variates, range)
# The ranges are primarily on the same scale

variates <- lapply(rgcca_results$variates, function(x) {
  x = as.data.frame(x)
})

for (i in 1:length(variates)) { 
  colnames(variates[[i]]) = paste0(names(variates)[i], "_", 
                                     colnames(variates[[i]]))
  variates[[i]]$Sample.ID = rownames(variates[[i]])
}

var_df = variates[[1]] %>% 
  inner_join(variates[[2]], by = "Sample.ID") %>%
  inner_join(variates[[3]], by = "Sample.ID") %>%
  inner_join(variates[[4]], by = "Sample.ID") %>%
  inner_join(variates[[5]], by = "Sample.ID") %>%
  tibble::column_to_rownames(var = "Sample.ID")

# Inspect distributions
library(tidyr)
create_density_plot_color(t(var_df)) +
  scale_x_continuous(limits = c(NA, NA))

plot(cmdscale(Rfast::Dist(var_df, "euclidean")))
# There is one outlier sample and a few which diverge from the rest

# To avoid the effect of outliers on distance calculations, we will use the
# PAM algorithm instead of k-means for clustering

# # Standardize variate matrices and concatenate for clustering:
# minmax_normalize_vector = function(vec, newmin, newmax) {
#   scaled_vec = (vec - min(vec)) / (max(vec) - min(vec))*(newmax-newmin) + newmin
#   return(scaled_vec)
# }
# 
# zvar_df = as.data.frame(do.call(cbind, lapply(var_df[colnames(var_df)], 
#                                               function(x) minmax_normalize_vector(x, -5, 5))))
# colnames(zvar_df) = colnames(var_df)
# rownames(zvar_df) = rownames(var_df)

# M3C k-means clustering #####
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
m3c_input = as.data.frame(t(zvar_df))

RNGversion("4.2.2")
consensus_km = M3C(m3c_input, des = m3c_des, iters = 100, repsref = 250, 
                   repsreal = 250, seed = 123, fsize = 18, lthick = 2, dotsize = 1.25,
                   clusteralg = "pam", maxK = 10)

optk = 3 # p = 1.450502e-208 - All p's but K=10 (which makes sense) are significant
paste0(ifelse(consensus_km$scores$NORM_P[consensus_km$scores$K == 3] < 0.05, "The clustering is significant.",
              "The clustering is not significant."))

# Inspection of clustering results #####
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
       path = paste0("Results/single_algorithm/", algorithm), 
       width = 1920, height = 1080, device = 'png', units = "px",
       dpi = 700)
dev.off()

# Entropy plot
entropy = ggplot(consensus_km[["plots"]][[2]][["data"]], aes(x = K, y = PAC_SCORE, alpha = 0.7))+
  geom_line(aes(color = "#7c1d6f"), linewidth = 0.5)+
  scale_x_continuous(limits = c(1.9, 10.1), breaks = seq(2, 10, 1))+
  scale_y_continuous(limits = c(0, 100000), breaks = seq(0, 100000, 10000),
                     labels = scales::label_comma()) +
  scale_colour_manual(values=rcartocolor::carto_pal(n = 9, "Safe"), name="NN") +
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
       path = paste0("Results/single_algorithm/", algorithm), 
       width = 1920, height = 1080, device = 'png', units = "px",
       dpi = 700)
dev.off()

# Statistical significance of clusters
inf_indices = which(consensus_km[["plots"]][[3]][["data"]]$P_SCORE == Inf)
inf_boolean = length(inf_indices) > 0
if (inf_boolean) {
  new_pscore = consensus_km[["plots"]][[3]][["data"]]
  new_pscore$P_SCORE[inf_indices] = max(new_pscore$P_SCORE[-inf_indices])*1.2
  
  statsig_clust = ggplot(new_pscore, aes(x = K, y = P_SCORE, color = P_SCORE < -log10(0.05)))+
    geom_point(size = 1.5, alpha = 0.6)+
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.2)+
    scale_color_manual(name = "Color",
                       values = c("#6c2167", "grey"),
                       labels = c("p < 0.05", "p > 0.05")) +
    scale_x_continuous(limits = c(1.9, 10.1), breaks = seq(2, 10, 1))+
    scale_y_continuous(limits = c(min(new_pscore$P_SCORE) - 0.15, 
                                  max(new_pscore$P_SCORE) + 0.15), 
                       breaks = seq(0, max(new_pscore$P_SCORE) + 0.1, 20))+
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
         path = paste0("Results/single_algorithm/", algorithm), 
         width = 1920, height = 1080, device = 'png', units = "px",
         dpi = 700)
} else {
  statsig_clust = ggplot(consensus_km[["plots"]][[3]][["data"]], 
                         aes(x = K, y = P_SCORE, color = P_SCORE < -log10(0.05)))+
    geom_point(size = 1.5, alpha = 0.6)+
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.2)+
    scale_color_manual(name = "Color",
                       values = c("#6c2167", "grey"),
                       labels = c("p < 0.05", "p > 0.05")) +
    scale_x_continuous(limits = c(1.9, 10.1), breaks = seq(2, 10, 1))+
    scale_y_continuous(limits = c(min(consensus_km[["plots"]][[3]][["data"]]$P_SCORE) - 0.15, 
                                  max(consensus_km[["plots"]][[3]][["data"]]$P_SCORE) + 0.15), 
                       breaks = seq(0, max(consensus_km[["plots"]][[3]][["data"]]$P_SCORE) + 0.1, 1))+
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
         path = paste0("Results/single_algorithm/", algorithm), 
         width = 1920, height = 1080, device = 'png', units = "px",
         dpi = 700)
}

# RCSI plot
rcsi = ggplot(as.data.frame(consensus_km[["scores"]]), aes(x = consensus_km$scores$K,
                                                           y = consensus_km$scores$RCSI))+
  geom_line(size = 0.3, color = "violet")+
  geom_errorbar(aes(ymin = consensus_km$scores$RCSI - consensus_km$scores$RCSI_SE,
                    ymax = consensus_km$scores$RCSI + consensus_km$scores$RCSI_SE,
                    color = "deeppink3"), width = 0.2, size = 0.1)+
  geom_point(size = 0.05, color ="deeppink3")+
  scale_x_continuous(limits = c(1.9, 10.1), breaks = seq(2, 10, 1))+
  scale_y_continuous(limits = c(min(consensus_km[["scores"]]$RCSI - consensus_km[["scores"]]$RCSI_SE) - 0.15, 
                                max(consensus_km[["scores"]]$RCSI + consensus_km[["scores"]]$RCSI_SE) + 0.15), 
                     breaks = c(-rev(seq(0, abs(round(min(consensus_km[["scores"]]$RCSI - consensus_km[["scores"]]$RCSI_SE), 1)), 0.5)), 
                                seq(0, round(max(consensus_km[["scores"]]$RCSI + consensus_km[["scores"]]$RCSI_SE), 1), 0.5)))+
  scale_colour_manual(values=rcartocolor::carto_pal(n = 9, "Safe"), name="NN") +
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
       path = paste0("Results/single_algorithm/", algorithm), 
       width = 1920, height = 1080, device = 'png', units = "px",
       dpi = 700)
dev.off()

# Conclusion
conclusion = paste0("The optimal value for k is ", optk, " ($p = ", consensus_km$scores$NORM_P[optk-1],
                   ", RCSI = ", consensus_km$scores$RCSI[optk-1], "$).")

# Main results #####
# Examine cluster similarity to MOVICS by measuring NMI and ARI indices #####
# (Jaccard may be misleading)
RGCCA_clusters = as.data.frame(consensus_km$realdataresults[[optk]]$assignments) %>%
  tibble::rownames_to_column(var = "Sample.ID")
colnames(RGCCA_clusters)[2] = "Cluster"
RGCCA_clusters$Sample.ID = gsub("\\.", "-", RGCCA_clusters$Sample.ID)
rownames(RGCCA_clusters) = RGCCA_clusters$Sample.ID

# Calculate ARI and NMI
library(mclust)
library(clue)

ARI_to_MOVICS = calculate_ari_index(cluster_df1 = ground_truth_labels,
                                    cluster_df2 = RGCCA_clusters,
                                    sample_col = "Sample.ID",
                                    clust_col = "Cluster",
                                    suffixes = c("_MOVICS_CS",
                                                 paste0("_", algorithm)))

NMI_to_MOVICS = calculate_nmi_index(cluster_df1 = ground_truth_labels,
                                    cluster_df2 = RGCCA_clusters,
                                    sample_col = "Sample.ID",
                                    clust_col = "Cluster",
                                    suffixes = c("_MOVICS_CS",
                                                 paste0("_", algorithm)))

# Very low statistics when compared to the MOVICS. Results differ

# Import coloring scheme
scheme = readRDS("Resources/scheme.rds")
annCol = scheme$annCol
annColors = scheme$annColors
cluster_colors = scheme$clust.colors
col.list = scheme$col.list
var2comp = scheme$var2comp %>%
  dplyr::select(-`Consensus Subtype`) %>%
  mutate(Sample.ID = rownames(.)) %>%
  inner_join(RGCCA_clusters, by = "Sample.ID") %>%
  tibble::column_to_rownames(var = "Sample.ID") %>%
  mutate(RGCCA = paste0(algorithm, Cluster)) %>%
  dplyr::select(RGCCA, everything()) %>%
  dplyr::select(-Cluster)
rm(scheme); gc()

# Silhouette
library(MOVICS)
library(cluster)
silhouette = silhouette(as.integer(gsub("RGCCA", "", RGCCA_clusters$Cluster)),
                        dist = Rfast::Dist(zvar_df, method = "euclidean"))

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

plot_object = list(clust.res = RGCCA_clusters %>%
                     dplyr::rename(samID = Sample.ID, clust = Cluster))

# Export consensus clustering object
clust = as.data.frame(plot_object$clust.res)
colnames(clust) = c("Sample.ID", "Cluster")
openxlsx::write.xlsx(clust, paste0(home, "/Results/single_algorithm/", algorithm, "/",
                                   algorithm, "_", data_source, "_",
                                   data_types, "_eval_on_", evaluation_source,
                                   "_clusterings.xlsx"))