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
algorithm = "SGCCA"
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

# However SGCCA prefers features in columns so we transpose the matrices.

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
saveRDS(input, "Resources/SGCCA_input.rds")

# SGCCA HPC scripts can be found in Scripts/single_algorithm/SGCCA_HPC/
# There is one script per penalty value
# All scripts use the full design (1 - diag(length(input))) like in the RGCCA case

# Import HPC results
sgcca_results = list()
results_indices = grep(".rds", list.files("Resources/HPC output/SGCCA_HPC/"))

for (filename in list.files("Resources/HPC output/SGCCA_HPC/")[results_indices]) {
  pen_val = strsplit(filename, "_")[[1]][3]
  sgcca_results[[paste0("penalty = ", pen_val)]] = readRDS(paste0("Resources/HPC output/SGCCA_HPC/",
                                                                 filename))
}

# Extract the variates from each result and proceed in an RGCCA way ###
variates = list()
var_dfs = list()
for (i in 1:length(sgcca_results)) {
  variates[[i]] <- lapply(sgcca_results[[i]]$variates, function(x) {
    x = as.data.frame(x)
  })
  
  for (j in 1:length(variates[[i]])) {
    colnames(variates[[i]][[j]]) = paste0(names(variates[[i]])[j], "_", 
                                     colnames(variates[[i]][[j]]))
    variates[[i]][[j]]$Sample.ID = rownames(variates[[i]][[j]])
  }
  
  var_dfs[[i]] = variates[[i]][[1]] %>% 
    inner_join(variates[[i]][[2]], by = "Sample.ID") %>%
    inner_join(variates[[i]][[3]], by = "Sample.ID") %>%
    inner_join(variates[[i]][[4]], by = "Sample.ID") %>%
    inner_join(variates[[i]][[5]], by = "Sample.ID") %>%
    tibble::column_to_rownames(var = "Sample.ID")
}

names(var_dfs) = names(variates) = names(sgcca_results)

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

# Run a loop of M3C for all results. In the end pick the one that achieves the
# best combo of entropy, RCSI, p-value

M3C_clusterings = list()
for (i in 1:length(var_dfs)) {
  m3c_input = t(var_dfs[[i]]) %>% as.data.frame()
  
  RNGversion("4.2.2")
  M3C_clusterings[[names(var_dfs)[i]]] = M3C(m3c_input, des = m3c_des,
                                                       iters = 100, repsref = 250, 
                                                       repsreal = 250, seed = 123, fsize = 18, lthick = 2, dotsize = 1.25,
                                                       clusteralg = "pam", maxK = 10)
  cat("Done with", names(var_dfs)[i], ".", "\n")
}

# console logs:
# penalty = 0.1, optimal K: 2
# penalty = 0.2, optimal K: 2
# penalty = 0.3, optimal K: 2
# penalty = 0.4, optimal K: 2
# penalty = 0.5, optimal K: 2
# penalty = 0.6, optimal K: 2
# penalty = 0.7, optimal K: 3
# penalty = 0.8, optimal K: 3
# penalty = 0.9, optimal K: 3

# Inspection of clustering results #####
scores_df = as.data.frame(cbind(list(Penalty = sort(rep(paste0("Penalty = ", seq(0.1, 0.9, 0.1)), 9))), 
                                rbind(M3C_clusterings[["penalty = 0.1"]][["scores"]],
                                      M3C_clusterings[["penalty = 0.2"]][["scores"]],
                                      M3C_clusterings[["penalty = 0.3"]][["scores"]],
                                      M3C_clusterings[["penalty = 0.4"]][["scores"]],
                                      M3C_clusterings[["penalty = 0.5"]][["scores"]],
                                      M3C_clusterings[["penalty = 0.6"]][["scores"]],
                                      M3C_clusterings[["penalty = 0.7"]][["scores"]],
                                      M3C_clusterings[["penalty = 0.8"]][["scores"]],
                                      M3C_clusterings[["penalty = 0.9"]][["scores"]])))

# RCSI plot
rcsi = ggplot(scores_df, aes(x = K, y = RCSI, group = Penalty, color = factor(Penalty)))+
  geom_line(linewidth = 0.3*2)+
  # geom_errorbar(aes(ymin = RCSI - RCSI_SE,
  #                   ymax = RCSI + RCSI_SE,
  #                   color = "grey"), width = 0.05, linewidth = 0.1*2)+
  geom_point(size = 1)+
  scale_x_continuous(limits = c(1.9, 10.1), breaks = seq(2, 10, 1))+
  scale_y_continuous(limits = c(min(scores_df$RCSI - scores_df$RCSI_SE) - 0.15, 
                                max(scores_df$RCSI + scores_df$RCSI_SE) + 0.15), 
                     breaks = c(-rev(seq(0, abs(round(min(scores_df$RCSI - scores_df$RCSI_SE), 1)), 0.5)), 
                                seq(0, round(max(scores_df$RCSI + scores_df$RCSI_SE), 1), 0.5)))+
  scale_colour_manual(values=rcartocolor::carto_pal(n = 9, "Safe"), name="Penalty") +
  theme(plot.title = element_text(size = 5*2, face = "bold"),
        axis.title.x = element_text(size = 4*2, face = "bold"),
        axis.title.y = element_text(size = 4*2, face = "bold"),
        axis.ticks = element_line(linewidth = 0.15*2),
        axis.text.x = element_text(size = 4*2),
        axis.text.y = element_text(size = 4*2),
        legend.position = "right",
        legend.text = element_text(size = 5),
        legend.title = element_text(face = "bold", size = 6.5, hjust = 0.5),
        legend.key.spacing = unit(2, "mm"),
        panel.background = element_rect(fill = "white", 
                                        colour = "white"),
        panel.grid = element_blank(),
        axis.line = element_line(linewidth = 0.2*2))+
  labs(title = "RCSI vs. number of clusters K for different Penalty values",
       x = "K", y = "RCSI")
rcsi
ggsave(filename = "RCSI.pdf",
       path = paste0("Results/single_algorithm/", algorithm), 
       width = 140, height = 100, device = 'pdf', units = "mm",
       dpi = 350)
dev.off()

# Statistical significance of clusters
library(ggnewscale)
statsig_clust = ggplot(scores_df, aes(x = K, y = P_SCORE, color = factor(Penalty)))+
  geom_point(size = 1.5*2, alpha = 0.6)+
  scale_color_manual(values=c(rcartocolor::carto_pal(n = 9, "Safe")),
                     name="Penalty") +
  new_scale("color") +
  geom_hline(linetype = "dashed", linewidth = 0.2*2,
             aes(yintercept = -log10(0.05), color = "grey40"))+
  scale_color_manual(values="grey40", labels = expression(-log[10]("0.05")),
                     name="Statistical \nsignificance") +
  scale_x_continuous(limits = c(1.9, 10.1), breaks = seq(2, 10, 1))+
  scale_y_continuous(limits = c(min(scores_df$P_SCORE) - 0.15, 
                                max(scores_df$P_SCORE) + 0.15), 
                     breaks = seq(0, max(scores_df$P_SCORE) + 0.1, 1))+
  theme_bw()+
  theme(plot.title = element_text(size = 5*2, face = "bold"),
        axis.title.x = element_text(size = 4*2, face = "bold"),
        axis.title.y = element_text(size = 4*2, face = "bold"),
        axis.ticks = element_line(linewidth = 0.15*2),
        axis.text.x = element_text(size = 4*2),
        axis.text.y = element_text(size = 4*2),
        legend.position = "right",
        legend.text = element_text(size = 5),
        legend.title = element_text(face = "bold", size = 6.5, hjust = 0.5),
        legend.key.spacing = unit(2, "mm"),
        panel.background = element_rect(fill = "white", 
                                        colour = "white"),
        panel.grid = element_blank(),
        axis.line = element_line(linewidth = 0.2*2))+
  labs(title = "Statistical significance of different values of K",
       y = bquote(bold(-log[10]("p"))))
statsig_clust
ggsave(filename = "Stat_Sig.pdf",
       path = paste0("Results/single_algorithm/", algorithm), 
       width = 140, height = 100, device = 'pdf', units = "mm",
       dpi = 350)
dev.off()

# Entropy
entropy = ggplot(scores_df, aes(x = K, y = ENTROPY_REAL, alpha = 0.85, 
                                group = Penalty, color = factor(Penalty)))+
  geom_line(linewidth = 0.6)+
  geom_point(size = 1) +
  scale_x_continuous(limits = c(1.9, 10.1), breaks = seq(2, 10, 1))+
  scale_y_continuous(limits = c(0, 100000), breaks = seq(0, 100000, 10000),
                     labels = scales::label_comma()) +
  scale_colour_manual(values=rcartocolor::carto_pal(n = 9, "Safe"), name="Penalty") +
  theme_bw()+
  theme(panel.border = element_rect(linewidth = 0.2*2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 5*2, face = "bold", vjust = 0.5, hjust = 0.5),
        legend.title = element_text(face = "bold", size = 4*2, hjust = 0.5),
        legend.text = element_text(size = 3*2),
        legend.key.size = unit(0.2*2, "cm"),
        legend.margin = ggplot2::margin(0, 0, 0, 0, unit = "mm"),
        legend.spacing.y = unit(0.5*2, units = "mm"),
        axis.title.x = element_text(size = 4*2, face = "bold"),
        axis.title.y = element_text(size = 4*2, face = "bold"),
        axis.ticks = element_line(linewidth = 0.15*2),
        axis.text.x = element_text(size = 4*2),
        axis.text.y = element_text(size = 4*2))+
  labs(y = "Entropy",
       x = "K",
       title = "Entropy in Real Data")+
  guides(alpha = "none")
entropy
ggsave(filename = "Entropy.pdf",
       path = paste0("Results/single_algorithm/", algorithm), 
       width = 140, height = 100, device = 'pdf', units = "mm",
       dpi = 350)
dev.off()

# Determine the best clustering based on p-value filtering and then examining the RCSI
# Entropy is biased for higher K
best_clusterings = scores_df[scores_df$NORM_P < 0.05, ] %>%
  dplyr::arrange(desc(RCSI))

# According to these criteria the best clustering is:
print(best_clusterings[1, ]) # Penalty = 0.3, K = 2
optk = 2

# Determine optPenalty
optPen = as.numeric(substr(best_clusterings$Pen[1], 11, 14))
conclusion = paste0("Best penalty value based on statistical significance and RCSI is", 
                     optPen, "$. We therefore proceed with $penalty = ",
                     optPen, "$.")

# Main results #####
# Examine cluster similarity to MOVICS by measuring NMI and ARI indices #####
# (Jaccard may be misleading)
SGCCA_clusters = as.data.frame(M3C_clusterings[[paste0("penalty = ", optPen)]]$realdataresults[[optk]]$assignments) %>%
  tibble::rownames_to_column(var = "Sample.ID")
colnames(SGCCA_clusters)[2] = "Cluster"
SGCCA_clusters$Sample.ID = gsub("\\.", "-", SGCCA_clusters$Sample.ID)
rownames(SGCCA_clusters) = SGCCA_clusters$Sample.ID

# Calculate ARI and NMI
library(mclust)
library(clue)

ARI_to_MOVICS = calculate_ari_index(cluster_df1 = ground_truth_labels,
                                    cluster_df2 = SGCCA_clusters,
                                    sample_col = "Sample.ID",
                                    clust_col = "Cluster",
                                    suffixes = c("_MOVICS_CS",
                                                 paste0("_", algorithm)))

NMI_to_MOVICS = calculate_nmi_index(cluster_df1 = ground_truth_labels,
                                    cluster_df2 = SGCCA_clusters,
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
  inner_join(SGCCA_clusters, by = "Sample.ID") %>%
  tibble::column_to_rownames(var = "Sample.ID") %>%
  mutate(SGCCA = paste0(algorithm, Cluster)) %>%
  dplyr::select(SGCCA, everything()) %>%
  dplyr::select(-Cluster)
rm(scheme); gc()

# Silhouette
library(MOVICS)
library(cluster)
silhouette = silhouette(as.integer(gsub("SGCCA", "", SGCCA_clusters$Cluster)),
                        dist = Rfast::Dist(var_dfs[[paste0("penalty = ", optPen)]], 
                                           method = "manhattan"))
# Manhattan is conceptually closer to PAM

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
                   function(mat) mat[, colSums(mat != 0) > 0])
plotdata <- lapply(plotdata, t)
plotdata = getStdiz(
  data = plotdata,
  halfwidth = c(NA, 3, 3, 3, 3), # No halfwidth for SNPs
  centerFlag = c(F, F, F, F, F),
  scaleFlag = c(F, F, F, F, F)
)


plot_object = list(clust.res = SGCCA_clusters %>%
                     dplyr::rename(samID = Sample.ID, clust = Cluster))

# Export consensus clustering object
clust = as.data.frame(plot_object$clust.res)
colnames(clust) = c("Sample.ID", "Cluster")
openxlsx::write.xlsx(clust, paste0(home, "/Results/single_algorithm/", algorithm, "/",
                                   algorithm, "_", data_source, "_",
                                   data_types, "_eval_on_", evaluation_source,
                                   "_clusterings.xlsx"))

# Save environment
save.image(paste0(home, "/Results/single_algorithm/", 
                  algorithm, "/", algorithm, "_", data_source, "_",
                  data_types, "_eval_on_", evaluation_source,
                  "_env.RData"))
