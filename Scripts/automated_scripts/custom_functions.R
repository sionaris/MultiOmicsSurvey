# A series of functions defined to modularize code and make cleaner

# In a nutshell #####

# Fetch the "in a nutshell" text of an algorithm
fetch_in_a_nutshell = function (algorithm) {
  source("Resources/algorithm_descriptions/in_a_nutshell.R")
  return(desc_list[[algorithm]])
}

# Fetch the citation of an algorithm
fetch_citation = function (algorithm) {
  source("Resources/algorithm_descriptions/citations.R")
  return(citations[[algorithm]])
}

# Function to compute Frobenius norm between two matrices
frobenius_norm <- function(mat1, mat2) {
  return(sqrt(sum((mat1 - mat2)^2)))
}

# Function to compute Pearson correlation between two matrices
pearson_correlation <- function(mat1, mat2) {
  cor(as.vector(mat1), as.vector(mat2))
}

# Function to reshape Pearson similarity values for ANOVA
reshape_SNF_Pearson_for_anova <- function(similarities) {
  data <- data.frame()
  for (key in names(similarities)) {
    pearson_matrix <- similarities[[key]]$Pearson
    values <- pearson_matrix[upper.tri(pearson_matrix)]
    factor <- rep(key, length(values))
    data <- rbind(data, data.frame(Value = values, Factor = factor))
  }
  return(data)
}

# Function to calculate Jaccard index between two clusterings (MOVICS)
MOVICS_jaccard_index <- function(clust1, clust2) {
  clust = as.data.frame(clust1) %>% 
    dplyr::rename(clust1 = clust) %>%
    inner_join(as.data.frame(clust2) %>% dplyr::rename(clust2 = clust),
               by = "samID")
  clust$agree = ifelse(clust$clust1 == clust$clust2, 1, 0)
  intersection <- sum(clust$agree)
  jaccard <- intersection/nrow(clust)
  return(jaccard)
}

# ARI index between two clusterings
calculate_ari_index <- function(cluster_df1, cluster_df2,
                                sample_col, clust_col, suffixes) {
  # Load the mclust package for ARI calculation
  library(mclust)
  
  # Merge data frames by Sample.ID
  merged_df <- merge(cluster_df1, cluster_df2, by = sample_col, suffixes = suffixes)
  
  # Calculate the ARI
  ari_index <- adjustedRandIndex(merged_df[[paste0(clust_col, suffixes[1])]], 
                                 merged_df[[paste0(clust_col, suffixes[2])]])
  return(ari_index)
}

# NMI index between two clusterings
calculate_nmi_index <- function(cluster_df1, cluster_df2,
                                sample_col, clust_col, suffixes) {
  # Load the clue package for NMI calculation
  library(clue)
  
  # Merge data frames by Sample.ID
  merged_df <- merge(cluster_df1, cluster_df2, by = sample_col, suffixes = suffixes)
  
  # Convert clusters to partitions
  partition1 <- as.cl_partition(merged_df[[paste0(clust_col, suffixes[1])]])
  partition2 <- as.cl_partition(merged_df[[paste0(clust_col, suffixes[2])]])
  
  # Calculate the NMI
  nmi_index <- cl_agreement(partition1, partition2, method = "NMI")
  return(nmi_index)
}

# create_MO_heatmap #####
create_MO_heatmap = function(matrix = NULL, algorithm = NULL, 
                             need.diag.zero = TRUE, 
                             clust_annot_pheno = NULL,
                             afh_colnames = NULL, colors = NULL,
                             annColors = NULL,
                             heatmap_title = NULL,
                             cluster_colors = NULL,
                             legend_title = NULL,
                             output_file_name = NULL) {
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  
  annotation_for_heatmap = clust_annot_pheno
  
  names(cluster_colors) = sort(unique(annotation_for_heatmap[, algorithm]))
  annColors[[algorithm]] = cluster_colors
  
  # Set all diagonals to zero for better visualization
  if (need.diag.zero) {
    diag(matrix) = 0
  }
  
  # Define breakpoints for the color mapping
  breaks <- seq(min(matrix, na.rm = TRUE),
                max(matrix, na.rm = TRUE), 
                length.out = length(colors))
  
  # Define the color function using colorRamp2
  color_fun <- circlize::colorRamp2(breaks, colors)
  
  # Order samples based on final SNF clusters:
  order = clust_annot_pheno %>%
    dplyr::arrange(!!sym(algorithm)) %>%
    dplyr::select(samID, !!sym(algorithm))
  order = order$samID
  
  annotation_for_heatmap = annotation_for_heatmap[order, ]
  
  # Create a HeatmapAnnotation object if you have annotations
  ha <- HeatmapAnnotation(df = annotation_for_heatmap, col = annColors, 
                          which = "column", show_annotation_name = TRUE, 
                          gap = unit(2, "mm"),
                          annotation_name_side = "left",
                          annotation_name_gp = gpar(fontface = "bold", fontsize = 12))
  
  # Heatmap splits
  splits = cumsum(table(annotation_for_heatmap[, algorithm]))
  splits = splits[-length(splits)]
  splits = rep(1:(length(splits)+1), c(table(annotation_for_heatmap[, algorithm])))
  
  # Create Heatmap object
  heatmap <- Heatmap(matrix[order, order],
                     column_title = heatmap_title,
                     column_title_gp = gpar(fontsize = 12, fontface = "bold"),
                     col = color_fun,
                     top_annotation = ha,
                     cluster_rows = F,
                     cluster_columns = F,
                     show_row_names = FALSE,
                     show_column_names = FALSE,
                     show_row_dend = F,
                     show_column_dend = F,
                     #clustering_distance_rows = "euclidean",
                     #clustering_distance_columns = "euclidean",
                     row_split = splits,
                     column_split = splits,
                     heatmap_legend_param = list(
                       title = legend_title,
                       title_gp = grid::gpar(fontsize = 10, fontface = "bold"), 
                       labels_gp = grid::gpar(fontsize = 6),
                       legend_height = unit(5, "cm"),
                       grid_width = unit(0.5, "cm"),
                       title_position = "leftcenter-rot",
                       border = TRUE
                     ))
  
  # Increase spacing between the main plot and legend (example values, adjust as needed)
  draw(heatmap, heatmap_legend_side = "right", annotation_legend_side = "right",
       padding = unit(c(10, 10, 10, 10), "mm"))  # Add padding around the heatmap
  
  # Save as PNG
  png(output_file_name, width = 13, 
      height = 9, units = 'in', res = 700)
  draw(heatmap)
  dev.off()
  
  rm(annColors, annotation_for_heatmap, heatmap, legend_title, heatmap_title,
     breaks, splits, afh_colnames, output_file_name, colors, cluster_colors,
     order, color_fun, need.diag.zero, algorithm); gc()
}


