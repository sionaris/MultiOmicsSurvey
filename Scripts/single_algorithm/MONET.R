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
algorithm = "MONET"
alg_feature_pref = "cols" # Where does the algorithm expect the features to be
citation = fetch_citation(algorithm = algorithm)
data_source = "TCGA" # e.g. TCGA, TCGA-transNEOdata_types = "RNAseq-CNV-Methylation-miRNA-SNPs" # e.g. RNAseq, RNAseq-CNV-miRNA
evaluation_source = "transNEO" # e.g. TCGA transNEO-PARTNER 
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

# Create an input for MONET #####
library(HiClimR)
# library(psych) - Pearson is equivalent to phi for binary vectors
# Calculate correlations
offset = 0.2 # used in MONET paper
correlations = list()
correlations[["miRNA"]] = fastCor(t(input$miRNA)) # Pearson
correlations[["RNAseq"]] = fastCor(t(input$RNAseq))
correlations[["CNV"]] = fastCor(t(input$CNV))
correlations[["Methylation"]] = fastCor(t(input$Methylation))
correlations[["SNPs"]] = fastCor(t(input$SNPs)) # equivalent to phi

# Offset the correlations
offset_correlation = function(cor.matrix, offset) {
  avg_val = mean(mean(cor.matrix))
  weighted = cor.matrix - avg_val - offset
  diag(weighted) = 0
  return(weighted)
}

offset_correlations = lapply(correlations, function(x) {
  return(offset_correlation(cor.matrix = x, offset = offset))
})

# Export
correlations = lapply(correlations, function(x) {
  diag(x) = 0
  return(x)
})

offset_correlations = lapply(offset_correlations, function(x) {
  diag(x) = 0
  return(x)
})

for (i in 1:length(correlations)) {
  write.csv(correlations[[i]], paste0("Python/MONET/MO_Corrs/", names(correlations)[i], ".csv"))
  write.csv(offset_correlations[[i]], paste0("Python/MONET/MO_Corrs_offset/offset_", 
                                             names(offset_correlations)[i], ".csv"))
}

# Import results with reticulate #####
library(reticulate)
library(ggplot2)

# Create subdirectory for supplementary plots
if (!dir.exists(paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement"))) {
  dir.create(paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement"))
}

# Add the Python dir of this project to the Python path of reticulate
py_run_string("import sys; sys.path.insert(0, 'Python/MONET')")

py_run_string("import pickle")
# py_install("networkx==3.2.1")
# py_install("pandas==2.2.3")
py_run_string("import networkx as nx")
py_run_string("from monet import monet")
py_run_string("monet_no_offset_results = pickle.load(open('Python/MONET/monet_results/monet_full_output_no_offset.pkl', 'rb')) ")
py_run_string("monet_offset_results = pickle.load(open('Python/MONET/monet_results/monet_full_output_with_offset.pkl', 'rb')) ")

# Create R variables
MONET_no_offset = py$monet_no_offset_results
MONET_with_offset = py$monet_offset_results

omics_list = names(MONET_no_offset$glob_var$omics)
edges_dfs_no_offset = list()
edges_dfs_with_offset = list()
wts_no_offset = list()
wts_with_offset = list()
hist_no_offset = list()
hist_with_offset = list()

for (omic_name in omics_list) {
  cat("\n--- Processing Omic:", omic_name, "---\n")
  
  # Retrieve the Python edge list: This is a NetworkX edges(data=True) structure
  edge_data_view_no_offset <- MONET_no_offset$glob_var$omics[[omic_name]]$graph$edges(data=TRUE)
  edge_data_view_with_offset <- MONET_with_offset$glob_var$omics[[paste0("offset_", omic_name)]]$graph$edges(data=TRUE)
  edge_iter_no_offset <- reticulate::as_iterator(edge_data_view_no_offset)
  edge_iter_with_offset <- reticulate::as_iterator(edge_data_view_with_offset)
  
  all_edges_no_offset <- list()
  all_edges_with_offset <- list()
  
  # Walk through the iterator until it's exhausted:
  while (TRUE) {
    # iter_next() will return NULL once we're out of items
    item <- reticulate::iter_next(edge_iter_no_offset)
    if (is.null(item)) {
      break
    }
    # 'item' should be something like:
    # (u_node, v_node, dict(weight=...))
    all_edges_no_offset[[ length(all_edges_no_offset) + 1 ]] <- item
  }
  
  while (TRUE) {
    # iter_next() will return NULL once we're out of items
    item <- reticulate::iter_next(edge_iter_with_offset)
    if (is.null(item)) {
      break
    }
    # 'item' should be something like:
    # (u_node, v_node, dict(weight=...))
    all_edges_with_offset[[ length(all_edges_with_offset) + 1 ]] <- item
  }
  
  edges_df_no_offset <- do.call(
    rbind,
    lapply(all_edges_no_offset, function(e) {
      # e[[1]] = u_node, e[[2]] = v_node, e[[3]]$weight = edge weight
      c(
        "u" = as.character(e[[1]]),
        "v" = as.character(e[[2]]),
        "weight" = as.numeric(e[[3]]$weight)
      )
    })
  )
  edges_df_no_offset <- as.data.frame(edges_df_no_offset, stringsAsFactors = FALSE)
  edges_dfs_no_offset[[omic_name]] = edges_df_no_offset
  
  edges_df_with_offset <- do.call(
    rbind,
    lapply(all_edges_with_offset, function(e) {
      # e[[1]] = u_node, e[[2]] = v_node, e[[3]]$weight = edge weight
      c(
        "u" = as.character(e[[1]]),
        "v" = as.character(e[[2]]),
        "weight" = as.numeric(e[[3]]$weight)
      )
    })
  )
  edges_df_with_offset <- as.data.frame(edges_df_with_offset, stringsAsFactors = FALSE)
  edges_dfs_with_offset[[omic_name]] = edges_df_with_offset
  
  # Extract weights
  wts_no_offset[[omic_name]] <- as.numeric(edges_df_no_offset$weight)
  wts_with_offset[[omic_name]] <- as.numeric(edges_df_with_offset$weight)
  
  # Histograms
  hist_no_offset[[omic_name]] = ggplot(as.data.frame(list(Weights = wts_no_offset[[omic_name]])),
                                       aes(x = Weights)) +
    geom_histogram(fill = "skyblue", color = "grey90", bins = 20, linewidth = 0.2) +
    ggtitle(paste0(omic_name, " edge weight histogram: no offset")) +
    theme_classic() +
    theme(plot.title = element_text(size = 7, face = "bold", vjust = 0.5, hjust = 0.5),
          axis.text = element_text(size = 5, hjust = 0.5, vjust = 0.5, 
                                   color = "black"),
          axis.title = element_text(size = 6, face = "bold"),
          axis.ticks = element_line(linewidth = 0.1),
          axis.line = element_line(linewidth = 0.3))+
    labs(x = "Edge weight", y = "Frequency")
  ggsave(plot = hist_no_offset[[omic_name]],
         filename = paste0(omic_name,
                           "_no_offset_weights_histogram.png"),
         path = paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement"),
         width = 2320, height = 2320, device = 'png', units = "px",
         dpi = 700)
  
  hist_with_offset[[omic_name]] = ggplot(as.data.frame(list(Weights = wts_with_offset[[omic_name]])),
                                        aes(x = Weights)) +
    geom_histogram(fill = "skyblue", color = "grey90", bins = 20, linewidth = 0.2) +
    ggtitle(paste0(omic_name, " edge weight histogram: with offset = ", offset)) +
    theme_classic() +
    theme(plot.title = element_text(size = 7, face = "bold", vjust = 0.5, hjust = 0.5),
          axis.text = element_text(size = 5, hjust = 0.5, vjust = 0.5, 
                                   color = "black"),
          axis.title = element_text(size = 6, face = "bold"),
          axis.ticks = element_line(linewidth = 0.1),
          axis.line = element_line(linewidth = 0.3))+
    labs(x = "Edge weight", y = "Frequency")
  ggsave(plot = hist_with_offset[[omic_name]],
         filename = paste0(omic_name,
                           "_with_offset_weights_histogram.png"),
         path = paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement"),
         width = 2320, height = 2320, device = 'png', units = "px",
         dpi = 700)
}

rm(item, omic_name, edges_df_with_offset, edges_df_no_offset, all_edges_no_offset,
   all_edges_with_offset, edge_data_view_with_offset, edge_data_view_no_offset,
   edge_iter_no_offset, edge_iter_with_offset); gc()

# Exporting a .txt report
writeLines(capture.output({invisible(lapply(seq_along(wts_no_offset), function (i) {
  x = wts_no_offset[[i]]
  omic = names(wts_no_offset)[i]
  cat("------- Omic:", omic, "\n")
  cat("Offset: 0 \n")
  cat("Number of edges:", length(x), "\n")
  cat("Edge weights summary:\n")
  print(summary(x))
  neg_frac = mean(x < 0)
  cat("\n--- Fraction of negative edges:", round(100*neg_frac, 2), "%")
  cat("\n-------------------\n\n")
}))}), paste0(home, "/Results/single_algorithm/", algorithm, "/MONET_no_offset_summary.txt"))

writeLines(capture.output({invisible(lapply(seq_along(wts_with_offset), function (i) {
  x = wts_with_offset[[i]]
  omic = names(wts_with_offset)[i]
  cat("------- Omic:", omic, "\n")
  cat("Offset:", offset, "\n")
  cat("Number of edges:", length(x), "\n")
  cat("Edge weights summary:\n")
  print(summary(x))
  neg_frac = mean(x < 0)
  cat("\n--- Fraction of negative edges:", round(100*neg_frac, 2), "%")
  cat("\n-------------------\n\n")
}))}), paste0(home, "/Results/single_algorithm/", algorithm, "/MONET_with_offset_summary.txt"))

# Offset conclusion
conclusion1 = "An offset of 0.2 is heavily shifting weights to the negative side with >95% negative values. The no offset options produces a more balanced distribution."

# Focus on the MONET results without offset
MONET = MONET_no_offset

# Table of MONET clusters
table(unlist(MONET[["clustering"]])) # There are three clusters
optk = as.numeric(length(unique(MONET[["clustering"]])))

# Create data frame
MONET_clusters = as.data.frame(list(Sample.ID = names(MONET$clustering),
                                    Cluster = unlist(MONET$clustering)))
rownames(MONET_clusters) = MONET_clusters$Sample.ID

# Extract the graphs of each modality
adj_matrices = list()
for (omic_name in omics_list) {
  atlas = MONET[["glob_var"]]$omics[[omic_name]]$graph$adj$`_atlas`
  adj_matrix = matrix(0, nrow = 625, ncol = 625, 
                      dimnames = list(names(MONET$clustering),
                                      names(MONET$clustering)))
  for (j in 1:length(atlas)) {
    sample1 <- names(atlas)[j]
    sample_atlas <- atlas[[j]]
    if (length(sample_atlas) > 0) {
      for (k in 1:length(sample_atlas)) {
        sample2 <- names(sample_atlas)[k]
        adj.value <- as.numeric(sample_atlas[[sample2]])
        adj_matrix[sample1, sample2] <- adj.value
      }
    } else {
      sample2 <- names(sample_atlas)
      adj.value <- 0
      adj_matrix[sample1, sample2] <- adj.value
    }
  }
  diag(adj_matrix) <- 1
  adj_matrices[[omic_name]] = adj_matrix
}

rm(adj_matrix, atlas, sample1, sample2, sample_atlas, adj.value, j, k); gc()

# Average graph
library(igraph)
avg_adjacency = Reduce(`+`, adj_matrices) / length(adj_matrices)
avg_graph = graph_from_adjacency_matrix(avg_adjacency, mode="undirected", weighted=TRUE)

# Inspect interpretability options #####
mods <- MONET$all_modules

module_table <- lapply(names(mods), function(id) {
  m <- mods[[id]]
  data.frame(
    Module_ID = id,
    Weight    = m$get_weight(),
    Omics     = paste(names(m$get_omics()), collapse = ",")
  )
})

module_df <- do.call(rbind, module_table)
print(module_df)

# Plotting sample composition per module with respect to HER2 and ER
composition_df <- data.frame(
  Sample.ID   = names(MONET$clustering),
  Module   = factor(unlist(MONET$clustering))
) %>% inner_join(clinical_data %>% dplyr::select(Sample.ID, 
                                                 ER = breast_carcinoma_estrogen_receptor_status,
                                                 HER2 = lab_proc_her2_neu_immunohistochemistry_receptor_status))

# Colors for ER and HER2
composition_df$ER = factor(composition_df$ER, levels = c("Negative", "Positive", ""),
                           labels = c("Negative", "Positive", "Unknown"))
composition_df$HER2 = factor(composition_df$HER2, levels = c("Negative", 
                                                             "Equivocal",
                                                             "Positive",
                                                             "",
                                                             "Indeterminate"),
                           labels = c("Negative", 
                                      "Equivocal",
                                      "Positive",
                                      "Unknown",
                                      "Indeterminate"))

# ER status
scale_fill_ER_status = scale_fill_manual(values = c(Negative = "#C11D9C", 
                                                    Positive = "#0F1682", 
                                                    Unknown = "grey40"))
# HER2 status
scale_fill_HER2_status = scale_fill_manual(values = c(Negative = "#0B9EF8", 
                                                      Positive = "#560DA7", 
                                                      Indeterminate = "mistyrose1", 
                                                      Equivocal = "hotpink4", 
                                                      Unknown = "grey40"))


# bar-plot of ER composition per module (using custom fill scales)
library(ggpubr)
p1 = ggplot(composition_df, aes(x = Module, fill = ER)) +
  geom_bar(position = "fill", width = 0.5) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "MONET module", y = "Proportion of Samples (%)", fill = "ER Status",
       title = "MONET module composition by ER status") +
  scale_fill_ER_status +
  theme_classic() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 5),
        axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 6, face = "bold"),
        legend.title = element_text(face = "bold", size = 5),
        plot.title = element_text(size = 7, face = "bold"),
        axis.text = element_text(size = 5, hjust = 0.5, vjust = 0.5, 
                                 color = "black"),
        legend.key.size = unit(0.25, "cm"),
        legend.text = element_text(size = 4),
        axis.ticks = element_line(linewidth = 0.1),
        axis.line = element_line(linewidth = 0.3))

p2 = ggplot(composition_df, aes(x = Module, fill = HER2)) +
  geom_bar(position = "fill", width = 0.5) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "MONET module", y = "Proportion of Samples (%)", fill = "HER2 Status",
       title = "MONET module composition by HER2 status") +
  scale_fill_HER2_status +
  theme_classic() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 5),
        axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 6, face = "bold"),
        legend.title = element_text(face = "bold", size = 5),
        plot.title = element_text(size = 7, face = "bold"),
        axis.text = element_text(size = 5, hjust = 0.5, vjust = 0.5, 
                                 color = "black"),
        legend.key.size = unit(0.25, "cm"),
        legend.text = element_text(size = 4),
        axis.ticks = element_line(linewidth = 0.1),
        axis.line = element_line(linewidth = 0.3))

FIG = ggarrange(p1, p2, ncol = 1, nrow = 2,
          labels = c("A", "B"),
          font.label = list(size = 10, face = "bold"))
ggsave(plot = FIG, 
       filename = paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement/",
                         "MONET_module_composition_by_ER_and_HER2_status.png"),
       width = 1024*2, height = 1024*3, device = 'png', units = "px", dpi = 700)

# Quantifying each module's reliance on each omic
edge_absolute_mean <- function(nx_graph) {
  w <- vapply(iterate(nx_graph$edges(data = TRUE)),
              function(e) abs(as.numeric(e[[3]]$weight)),
              numeric(1))
  mean(w)                         # just return the scalar mean
}

omics <- names(MONET$glob_var$omics)

# 1) background vector  bg[omic]
bg <- sapply(omics, function(o)
  edge_absolute_mean(MONET$glob_var$omics[[o]]$graph))

# 2) module-by-omic matrix  omic_strength[omic, module]
get_mod_edge_absmean <- function(mod, omic_name) {
  g   <- MONET$glob_var$omics[[omic_name]]$graph
  pts <- unlist(mod$get_patients_names_as_list())
  if (length(pts) < 2) return(NA)
  
  combs <- t(combn(pts, 2))
  w <- apply(combs, 1, function(x)
    if (g$has_edge(x[1], x[2]))
      abs(as.numeric(g$get_edge_data(x[1], x[2])$weight))
    else 0)
  mean(w)
}

mods  <- MONET$all_modules
omic_strength <- sapply(mods, function(m)
  sapply(omics, get_mod_edge_absmean, mod = m))

rownames(omic_strength) <- omics
colnames(omic_strength) <- paste0("Module_", names(mods))

# 3) “excess connectivity” centered by background
omic_strength_centered <- sweep(omic_strength, 1, bg, FUN = "-")
omic_strength_centered

# Examine cluster similarity to MOVICS by measuring NMI and ARI indices #####
# (Jaccard may be misleading)

# Calculate ARI and NMI
library(mclust)
library(clue)

ARI_to_MOVICS = calculate_ari_index(cluster_df1 = ground_truth_labels,
                                    cluster_df2 = MONET_clusters,
                                    sample_col = "Sample.ID",
                                    clust_col = "Cluster",
                                    suffixes = c("_MOVICS_CS",
                                                 paste0("_", algorithm)))

NMI_to_MOVICS = calculate_nmi_index(cluster_df1 = ground_truth_labels,
                                    cluster_df2 = MONET_clusters,
                                    sample_col = "Sample.ID",
                                    clust_col = "Cluster",
                                    suffixes = c("_MOVICS_CS",
                                                 paste0("_", algorithm)))

# Low statistics; results differ

# Import coloring scheme
scheme = readRDS("Resources/scheme.rds")
annCol = scheme$annCol
annColors = scheme$annColors
cluster_colors = c("#2EC4B6", "#E71D36", "#FF9F1C")
col.list = scheme$col.list
var2comp = scheme$var2comp %>%
  dplyr::select(-`Consensus Subtype`) %>%
  mutate(Sample.ID = rownames(.)) %>%
  inner_join(MONET_clusters, by = "Sample.ID") %>%
  tibble::column_to_rownames(var = "Sample.ID") %>%
  mutate(MONET = paste0(algorithm, Cluster)) %>%
  dplyr::select(MONET, everything()) %>%
  dplyr::select(-Cluster)
rm(scheme); gc()

# Silhouette based on the average adjacency
library(cluster)
silhouette = silhouette(as.integer(MONET_clusters[rownames(avg_adjacency), ]$Cluster),
           dmatrix = (1 - avg_adjacency)/2)

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

plot_object = list(clust.res = MONET_clusters %>%
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
                                         pdf_level_col_width = c("5em", "5em"),
                                         pdf_count_col_width = "5em",
                                         pdf_pval_col_width = "3em",
                                         pdf_test_col_width = "3em",
                                         pdf_tab_font_size = 10)

clin_ordinal_comp = compClinvar_ordinal_single_algorithm(algorithm_name = algorithm,
                                                         moic.res = plot_object,
                                                         var2comp = var2comp_nonas %>%
                                                           dplyr::select(number_of_lymphnodes_positive_by_ihc,
                                                                         number_of_lymphnodes_positive_by_he,
                                                                         MONET),
                                                         strata = algorithm,
                                                         ordinalVars = c("number_of_lymphnodes_positive_by_ihc",
                                                                         "number_of_lymphnodes_positive_by_he"),
                                                         includeNA = FALSE,
                                                         tab.name = "Summary of ordinal clinical variables",
                                                         res.path = paste0(home, "/Results/single_algorithm/", algorithm, "/"),
                                                         output_pdf = TRUE,
                                                         pdf_template_loc = paste0(home, "/Scripts/automated_scripts/clincomp_template.Rmd"),
                                                         pdf_level_col_width = c("5em", "5em"),
                                                         pdf_count_col_width = "5em",
                                                         pdf_pval_col_width = "3em",
                                                         pdf_test_col_width = "3em",
                                                         pdf_tab_font_size = 10)

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
                                      res.path     = paste0(home, "/Results/single_algorithm/", algorithm),
                                      test.method = "chisq" # large number of groups
)

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
                                            n.path       = 20,
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
                                              minGSSize = 10,
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

timestamp() # ~2.5h
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
                                                present_clusters = c(paste0("MONET", seq(1, optk, 1))),
                                                representative = TRUE, moic.res = plot_object,
                                                subtype_prefix = algorithm, n.path = 20, msigdb.path = MSIGDB.FILE,
                                                norm.method = "mean", dirct = "up",
                                                fig.name = "upregulated_pathway_heatmap",
                                                name = "GSVA scores",
                                                fig.path = paste0(home, "/Results/single_algorithm/", algorithm), 
                                                width = 15, height = 12, gsva.method = "gsva")

hclust_pathway_plots_down = plot_pathway_heatmaps(gsea.lists = hclust_output[grepl("down", names(hclust_output))], 
                                                  norm.expr = plotdata$RNAseq, 
                                                  present_clusters = c(paste0("MONET", seq(1, optk, 1))),
                                                  representative = TRUE, moic.res = plot_object,
                                                  subtype_prefix = algorithm, n.path = 20, msigdb.path = MSIGDB.FILE,
                                                  norm.method = "mean", dirct = "down",
                                                  fig.name = "downregulated_pathway_heatmap",
                                                  name = "GSVA scores",
                                                  fig.path = paste0(home, "/Results/single_algorithm/", algorithm), 
                                                  width = 15, height = 12, gsva.method = "gsva")

# Fraction Genome Altered ###
fga_df = readRDS("Resources/TCGA/fga_df.rds"); gc()

fga.MONET <- compFGA_optimized(moic.res     = plot_object,
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

fga.MONET.COSMIC <- compFGA_optimized(moic.res     = plot_object,
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
timestamp() # 8 min

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
timestamp() # 8 min

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
  inner_join(expr_conc %>% dplyr::select(Donor.ID = samID, MONET = clust_up),
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
                                                 pdf_level_col_width = c("5em", "5em"),
                                                 pdf_count_col_width = "5em",
                                                 pdf_pval_col_width = "3em",
                                                 pdf_test_col_width = "3em",
                                                 pdf_tab_font_size = 10)

transNEO_ntp_expr_up_ord = transNEO_ntp_expr_up
transNEO_ntp_expr_up_ord$clust.res$clust = gsub(algorithm, "", transNEO_ntp_expr_up_ord$clust.res$clust)
transNEO_ordinal_clincomp = compClinvar_ordinal_single_algorithm(algorithm_name = algorithm,
                                                                 moic.res = transNEO_ntp_expr_up_ord,
                                                                 var2comp = transNEO_var2comp_nonas %>%
                                                                   dplyr::select(Grade.pre.NAT, 
                                                                                 Chemo.cycles, 
                                                                                 aHER2.cycles,
                                                                                 MONET),
                                                                 strata = algorithm,
                                                                 ordinalVars = c("Grade.pre.NAT",
                                                                                 "Chemo.cycles",
                                                                                 "aHER2.cycles"),
                                                                 includeNA = FALSE,
                                                                 tab.name = "transNEO Summary of ordinal clinical variables",
                                                                 res.path = paste0(home, "/Results/single_algorithm/", algorithm, "/"),
                                                                 output_pdf = TRUE,
                                                                 pdf_template_loc = paste0(home, "/Scripts/automated_scripts/clincomp_template.Rmd"),
                                                                 pdf_level_col_width = c("5em", "5em"),
                                                                 pdf_count_col_width = "5em",
                                                                 pdf_pval_col_width = "3em",
                                                                 pdf_test_col_width = "3em",
                                                                 pdf_tab_font_size = 10)

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
cluster_colors_heatmap = cluster_colors
clust_annot_pheno = annCol %>% mutate(Sample.ID = rownames(.)) %>%
  inner_join(clust, by = "Sample.ID") %>%
  dplyr::rename(MONET = Cluster, samID = "Sample.ID")
rownames(clust_annot_pheno) = clust_annot_pheno$samID
afh_colnames = colnames(annCol)

# Same data frame. Different columns. Just for easiness
MONET_clust_res = MONET_clusters %>% dplyr::rename(samID = Sample.ID, 
                                                   MONET = Cluster) %>%
  dplyr::mutate(MONET = gsub(algorithm, "", MONET))

# CNV
create_MO_heatmap(matrix = adj_matrices$CNV, algorithm = algorithm, 
                  need.diag.zero = FALSE, 
                  clust_annot_pheno = clust_annot_pheno,
                  afh_colnames = afh_colnames, 
                  colors = colors_heatmap,
                  annColors = annColors,
                  heatmap_title = "CNV adjacency heatmap",
                  cluster_colors = cluster_colors_heatmap,
                  legend_title = "Adjacency",
                  cluster_cols_flag = FALSE,
                  cluster_rows_flag = FALSE,
                  splits_flag = TRUE,
                  output_file_name = paste0(home, "/Results/single_algorithm/", algorithm, 
                                            "/Supplement/adjacency_CNV_heatmap.png"))

# RNAseq
create_MO_heatmap(matrix = adj_matrices$RNAseq, algorithm = algorithm, 
                  need.diag.zero = FALSE, 
                  clust_annot_pheno = clust_annot_pheno,
                  afh_colnames = afh_colnames, 
                  colors = colors_heatmap,
                  annColors = annColors,
                  heatmap_title = "RNAseq adjacency heatmap",
                  cluster_colors = cluster_colors_heatmap,
                  legend_title = "Adjacency",
                  cluster_cols_flag = FALSE,
                  cluster_rows_flag = FALSE,
                  splits_flag = TRUE,
                  output_file_name = paste0(home, "/Results/single_algorithm/", algorithm, 
                                            "/Supplement/adjacency_RNAseq_heatmap.png"))

# miRNA
create_MO_heatmap(matrix = adj_matrices$miRNA, algorithm = algorithm, 
                  need.diag.zero = FALSE, 
                  clust_annot_pheno = clust_annot_pheno,
                  afh_colnames = afh_colnames, 
                  colors = colors_heatmap,
                  annColors = annColors,
                  heatmap_title = "miRNA adjacency heatmap",
                  cluster_colors = cluster_colors_heatmap,
                  legend_title = "Adjacency",
                  cluster_cols_flag = FALSE,
                  cluster_rows_flag = FALSE,
                  splits_flag = TRUE,
                  output_file_name = paste0(home, "/Results/single_algorithm/", algorithm,
                                            "/Supplement/adjacency_miRNA_heatmap.png"))

# Methylation
create_MO_heatmap(matrix = adj_matrices$Methylation, algorithm = algorithm, 
                  need.diag.zero = FALSE, 
                  clust_annot_pheno = clust_annot_pheno,
                  afh_colnames = afh_colnames, 
                  colors = colors_heatmap,
                  annColors = annColors,
                  heatmap_title = "Methylation adjacency heatmap",
                  cluster_colors = cluster_colors_heatmap,
                  legend_title = "Adjacency",
                  cluster_cols_flag = FALSE,
                  cluster_rows_flag = FALSE,
                  splits_flag = TRUE,
                  output_file_name = paste0(home, "/Results/single_algorithm/", algorithm, 
                                            "/Supplement/adjacency_Methylation_heatmap.png"))

# SNPs
create_MO_heatmap(matrix = adj_matrices$SNPs, algorithm = algorithm, 
                  need.diag.zero = FALSE, 
                  clust_annot_pheno = clust_annot_pheno,
                  afh_colnames = afh_colnames, 
                  colors = colors_heatmap,
                  annColors = annColors,
                  heatmap_title = "SNPs adjacency heatmap",
                  cluster_colors = cluster_colors_heatmap,
                  legend_title = "Adjacency",
                  cluster_cols_flag = FALSE,
                  cluster_rows_flag = FALSE,
                  splits_flag = TRUE,
                  output_file_name = paste0(home, "/Results/single_algorithm/", algorithm,
                                            "/Supplement/adjacency_SNPs_heatmap.png"))

# Final affinity matrix
create_MO_heatmap(matrix = avg_adjacency, algorithm = algorithm, 
                  need.diag.zero = FALSE, 
                  clust_annot_pheno = clust_annot_pheno,
                  afh_colnames = afh_colnames, 
                  colors = colors_heatmap,
                  annColors = annColors,
                  heatmap_title = "Average adjacency heatmap",
                  cluster_colors = cluster_colors_heatmap,
                  legend_title = "Adjacency",
                  cluster_cols_flag = FALSE,
                  cluster_rows_flag = FALSE,
                  splits_flag = TRUE,
                  output_file_name = paste0(home, "/Results/single_algorithm/", algorithm,
                                            "/Supplement/avg_adjacency_heatmap.png"))

# Final affinity matrix with clustered rows and columns
create_MO_heatmap(matrix = avg_adjacency, algorithm = algorithm, 
                  need.diag.zero = FALSE, 
                  clust_annot_pheno = clust_annot_pheno,
                  afh_colnames = afh_colnames, 
                  colors = colors_heatmap,
                  annColors = annColors,
                  heatmap_title = "Average adjacency heatmap",
                  cluster_colors = cluster_colors_heatmap,
                  legend_title = "Adjacency",
                  cluster_cols_flag = TRUE,
                  cluster_rows_flag = TRUE,
                  splits_flag = FALSE,
                  output_file_name = paste0(home, "/Results/single_algorithm/", algorithm,
                                            "/Supplement/hclust_avg_adjacency_heatmap.png"))

# PCA ###
# CNV
pca_from_sim_matrix(sim_matrix = adj_matrices$CNV, algorithm = algorithm, 
                    clust_res = clust_annot_pheno %>% dplyr::select(samID, MONET),
                    cluster_colors = cluster_colors_heatmap, 
                    output_path = paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement"), 
                    title_add = "CNV")

# RNAseq
pca_from_sim_matrix(sim_matrix = adj_matrices$RNAseq, algorithm = algorithm, 
                    clust_res = clust_annot_pheno %>% dplyr::select(samID, MONET),
                    cluster_colors = cluster_colors_heatmap, 
                    output_path = paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement"), 
                    title_add = "RNAseq")

# miRNA
pca_from_sim_matrix(sim_matrix = adj_matrices$miRNA, algorithm = algorithm, 
                    clust_res = clust_annot_pheno %>% dplyr::select(samID, MONET),
                    cluster_colors = cluster_colors_heatmap, 
                    output_path = paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement"), 
                    title_add = "miRNA")

# Methylation
pca_from_sim_matrix(sim_matrix = adj_matrices$Methylation, algorithm = algorithm, 
                    clust_res = clust_annot_pheno %>% dplyr::select(samID, MONET),
                    cluster_colors = cluster_colors_heatmap, 
                    output_path = paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement"), 
                    title_add = "Methylation")

# SNPs
pca_from_sim_matrix(sim_matrix = adj_matrices$SNPs, algorithm = algorithm, 
                    clust_res = clust_annot_pheno %>% dplyr::select(samID, MONET),
                    cluster_colors = cluster_colors_heatmap, 
                    output_path = paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement"), 
                    title_add = "SNPs")

# Final affinity
pca_from_sim_matrix(sim_matrix = avg_adjacency, algorithm = algorithm, 
                    clust_res = clust_annot_pheno %>% dplyr::select(samID, MONET),
                    cluster_colors = cluster_colors_heatmap, 
                    output_path = paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement"), 
                    title_add = "Average Adjacency")

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

# PR status
scale_fill_PR_status = scale_fill_manual(values = c(Indeterminate = "aliceblue", 
                                                    Positive = "dodgerblue4", 
                                                    Negative = "#F0C6C3", 
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
MONET_barcharts = list()
plotdata_bar = clust_annot_pheno_nonas
plotdata_bar[[algorithm]] = factor(plotdata_bar[[algorithm]])
for (i in 1:length(voi)) {
  chifit = output
  loc = which(grepl(voi[i], chifit$Comparison))
  chifit = chifit[loc, ]
  MONET_barcharts[[i]] = create_annot_barchart(plotdata = plotdata_bar, fill = voi[i],
                                               chifit = chifit,
                                               na.action = "na.omit",
                                               algorithm = algorithm,
                                               barchart_ylim = 450,
                                               text_y = 400, rect_ymin = 320,
                                               rect_ymax = 420, x_annot = 2,
                                               v_gap = 25, rect_xmin = 1.5,
                                               rect_xmax = 2.5, 
                                               annot_text_size = 2.25,
                                               legend.text.size = 5,
                                               x.axis.text.size = 5) +
    barchart_scales[[voi[i]]]
  print(MONET_barcharts[[i]])
  ggsave(filename = paste0(algorithm, "_", voi[i], "_barchart.png"),
         path = paste0(home, 
                       "/Results/single_algorithm/", algorithm, "/Supplement"), 
         width = 2620, height = 2320, device = 'png', units = "px",
         dpi = 700)
  dev.off()
}
names(MONET_barcharts) = voi
rm(loc, chifit)

# Multiplot (PNG) - bar charts
ggarrange(MONET_barcharts[[1]], MONET_barcharts[[2]], MONET_barcharts[[3]],
          MONET_barcharts[[4]], MONET_barcharts[[5]], MONET_barcharts[[6]],
          MONET_barcharts[[7]], MONET_barcharts[[8]], MONET_barcharts[[9]],
          MONET_barcharts[[10]], MONET_barcharts[[11]],
          ncol = 3, nrow = 4, labels = c("A", "B", "C", "D", "E", "F", "G", "H",
                                         "I", "J", "K"),
          font.label = list(size = 8, face = "bold", color ="black"))
ggsave(filename = paste0("Multiplot_", algorithm, "_barcharts.png"),
       path = paste0(home, 
                     "/Results/single_algorithm/", algorithm, "/Supplement"), 
       width = 9000, height = 8000, device = 'png', units = "px",
       dpi = 700)
dev.off()

# Just significant ones now
MONET_barcharts_sig = list()
plotdata_bar_sig = clust_annot_pheno_nonas %>% dplyr::select(MONET, Race, Histology, 
                                                             `ER status`, `PR status`, `HER2 status`,
                                                             `Lymph node status`, Stage)
plotdata_bar_sig$MONET = factor(plotdata_bar_sig$MONET)
voi_sig = setdiff(colnames(plotdata_bar_sig), algorithm)
for (i in 1:length(voi_sig)) {
  chifit = output
  loc = which(grepl(voi_sig[i], chifit$Comparison))
  chifit = chifit[loc, ]
  MONET_barcharts_sig[[i]] = create_annot_barchart(plotdata = plotdata_bar_sig, fill = voi_sig[i],
                                                   chifit = chifit,
                                                   na.action = "na.omit",
                                                   algorithm = algorithm,
                                                   barchart_ylim = 450,
                                                   text_y = 400, rect_ymin = 320,
                                                   rect_ymax = 420, x_annot = 2,
                                                   v_gap = 25, rect_xmin = 1.5,
                                                   rect_xmax = 2.5, 
                                                   annot_text_size = 2.25,
                                                   legend.text.size = 5,
                                                   x.axis.text.size = 5) +
    barchart_scales[[voi_sig[i]]]
  print(MONET_barcharts_sig[[i]])
  ggsave(filename = paste0("sig_", algorithm, "_", voi_sig[i], "_barchart.png"),
         path = paste0(home, 
                       "/Results/single_algorithm/", algorithm, "/Supplement"), 
         width = 2620, height = 2320, device = 'png', units = "px",
         dpi = 700)
  dev.off()
}
names(MONET_barcharts_sig) = voi_sig
rm(loc, chifit)

# Multiplot (PNG) - bar charts
ggarrange(MONET_barcharts_sig[[1]], MONET_barcharts_sig[[2]], MONET_barcharts_sig[[3]],
          MONET_barcharts_sig[[4]], MONET_barcharts_sig[[5]], MONET_barcharts_sig[[6]],
          MONET_barcharts_sig[[7]],
          ncol = 3, nrow = 3, labels = c("A", "B", "C", "D", "E", "F", "G"),
          font.label = list(size = 8, face = "bold", color ="black"))
ggsave(filename = paste0("sig_Multiplot_", algorithm, "_barcharts.png"),
       path = paste0(home, 
                     "/Results/single_algorithm/", algorithm, "/Supplement"), 
       width = 9000, height = 6500, device = 'png', units = "px",
       dpi = 700)
dev.off()

# Sunburst plot ###
library(plotly)
Pheno_sunburst_MONET = clust_annot_pheno
Pheno_sunburst_MONET$`ER status` = gsub("Unknown", "Unkn ER status", Pheno_sunburst_MONET$`ER status`)
Pheno_sunburst_MONET$`ER status` = gsub("Positive", "ER+", Pheno_sunburst_MONET$`ER status`)
Pheno_sunburst_MONET$`ER status` = gsub("Negative", "ER-", Pheno_sunburst_MONET$`ER status`)
Pheno_sunburst_MONET$`HER2 status` = gsub("Unknown", "Unkn HER2 status", 
                                          Pheno_sunburst_MONET$`HER2 status`)
Pheno_sunburst_MONET$`HER2 status` = gsub("Positive", "HER2+", Pheno_sunburst_MONET$`HER2 status`)
Pheno_sunburst_MONET$`HER2 status` = gsub("Negative", "HER2-", Pheno_sunburst_MONET$`HER2 status`)
Pheno_sunburst_MONET$`Lymph node status` = gsub("Unknown", 
                                                "Unkn LN. status", 
                                                Pheno_sunburst_MONET$`Lymph node status`)
Pheno_sunburst_MONET = Pheno_sunburst_MONET %>%
  dplyr::select(MONET, `ER status`, `HER2 status`, `Lymph node status`) %>%
  group_by(MONET, `ER status`, `HER2 status`, `Lymph node status`) %>%
  summarise(Counts = n()) %>%
  as.data.frame()
sunburst_coloring_MONET = data.frame(stringsAsFactors = FALSE,
                                     colors = tolower(gplots::col2hex(c("#2EC4B6", "#E71D36", "#FF9F1C", "#BDD5EA",
                                                                        "#FFA5AB", 
                                                                        "#C11D9C", "#0F1682",  "grey40",
                                                                        "#0B9EF8", "#560DA7", "mistyrose1", 
                                                                        "hotpink4", "grey40",
                                                                        "grey75", "#4A0558", "grey40"))),
                                     labels = c("MONET1", "MONET2", "MONET3", "MONET4",
                                                "MONET5",
                                                "ER-", "ER+", "Unkn ER status",
                                                "HER2-", "HER2+", "Indeterminate",
                                                "Equivocal", "Unkn HER2 status",
                                                "No", "Yes", "Unkn LN. status"))

sunburstDF_MONET = as.sunburstDF(Pheno_sunburst_MONET, value_column = "Counts", add_root = FALSE) %>%
  inner_join(sunburst_coloring_MONET, by = "labels")

pie_MONET = plot_ly() %>%
  add_trace(ids = sunburstDF_MONET$ids, labels= sunburstDF_MONET$labels, 
            parents = sunburstDF_MONET$parents, 
            values= sunburstDF_MONET$values, type='sunburst', branchvalues = 'total',
            insidetextorientation='radial', maxdepth = 5,
            marker = list(colors = sunburstDF_MONET$colors)) %>%
  layout(
    grid = list(columns =1, rows = 1),
    margin = list(l = 0, r = 0, b = 0, t = 0)
  )
pie_MONET
rm(Pheno_sunburst_MONET, sunburstDF_MONET, sunburst_coloring_MONET, pie_MONET); gc()

# Graphs ###
list_aff_S = list(`MONET CNV Adjacency` = adj_matrices$CNV,
                  `MONET RNAseq Adjacency` = adj_matrices$RNAseq,
                  `MONET Methylation Adjacency` = adj_matrices$Methylation,
                  `MONET SNPs Adjacency` = adj_matrices$SNPs,
                  `MONET miRNA Adjacency` = adj_matrices$miRNA,
                  `MONET Average Adjacency` = avg_adjacency)

for (i in 1:length(list_aff_S)) {
  # Scale all to [0, 1]
  adj = (list_aff_S[[i]] + 1)/2
  
  # Prepare the graph object
  g <- graph_from_adjacency_matrix(adj,  weighted = TRUE, diag = FALSE,
                                   mode = "undirected")
  g <- delete_edges(g, E(g)[weight == 0])
  E(g)$width <- sqrt(E(g)$weight) * 5  # Example transformation for visibility
  nodes_data <- data.frame(name = V(g)$name) %>%
    inner_join(clust_annot_pheno %>% dplyr::select(samID, MONET),
               by = c("name" = "samID"))
  
  # Set MONET as a factor for coloring
  nodes_data[[algorithm]] <- as.factor(nodes_data[[algorithm]])
  V(g)$MONET <- nodes_data[[algorithm]] # modify `$MONET` manually
  
  # Set color based on MONET
  V(g)$color <- fifelse(V(g)$MONET == paste0(algorithm, "1"), "#2EC4B6", 
                        fifelse(V(g)$MONET == paste0(algorithm, "2"), "#E71D36", "#FF9F1C"))
  
  png(paste0(home, 
             "/Results/single_algorithm/", algorithm, "/Supplement/",
             names(list_aff_S)[i], " graph.png"),
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
                  paste0(algorithm, "2"),
                  paste0(algorithm, "3")),
         fill=cluster_colors_heatmap,  
         cex=0.7,      
         box.lwd=1)  
  
  dev.off() 
}
rm(g, nodes_data)

# Wrap up #####
hyperparameters = list(iters = 10000,
                       num_of_seeds = 100,
                       num_of_samples_in_seed = 10,
                       min_mod_size = 10,
                       max_pats_per_action = 10,
                       percentile_shift = "None",
                       percentile_remove_edge = 80,
                       optk = optk,
                       conclusion1 = conclusion1)

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

