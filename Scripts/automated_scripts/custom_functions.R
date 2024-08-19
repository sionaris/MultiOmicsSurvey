# A series of functions defined to modularize code and make cleaner

# Fetch the "in a nutshell" text of an algorithm #####
fetch_in_a_nutshell = function (algorithm) {
  source("Resources/algorithm_descriptions/in_a_nutshell.R")
  return(desc_list[[algorithm]])
}

# Fetch the citation of an algorithm #####
fetch_citation = function (algorithm) {
  source("Resources/algorithm_descriptions/citations.R")
  return(citations[[algorithm]])
}

# Function to compute Frobenius norm between two matrices #####
frobenius_norm <- function(mat1, mat2) {
  return(sqrt(sum((mat1 - mat2)^2)))
}

# Function to compute Pearson correlation between two matrices #####
pearson_correlation <- function(mat1, mat2) {
  cor(as.vector(mat1), as.vector(mat2))
}

# Function to reshape Pearson similarity values for ANOVA #####
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

# Function to calculate Jaccard index between two clusterings (MOVICS) #####
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

# ARI index between two clusterings #####
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

# NMI index between two clusterings #####
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

# Normalize SNF affinity matrix #####
normalize_affinity_matrix <- function(W) {
  N <- nrow(W)  # Determine the size of the matrix
  P <- matrix(0, nrow = N, ncol = N)  # Initialize P with zeros
  
  # Calculate row sums of W, excluding the diagonal elements
  row_sums <- rowSums(W) - diag(W)
  
  # Fill P matrix based on the conditions
  for (i in 1:N) {
    for (j in 1:N) {
      if (i != j) {
        P[i, j] = W[i, j] / (2 * row_sums[i])
      } else {
        P[i, j] = 0.5
      }
    }
  }
  
  return(P)  # Return the normalized matrix
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
                             output_file_name = NULL,
                             cluster_rows_flag = NULL,
                             cluster_cols_flag = NULL,
                             splits_flag = NULL) {
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
  
  annotation_for_heatmap = annotation_for_heatmap[order, ] %>%
    dplyr::select(-samID)
  
  # Create a HeatmapAnnotation object if you have annotations
  ha <- HeatmapAnnotation(df = annotation_for_heatmap, col = annColors, 
                          which = "column", show_annotation_name = TRUE, 
                          gap = unit(2, "mm"),
                          annotation_name_side = "left",
                          annotation_name_gp = gpar(fontface = "bold", fontsize = 12))
  
  # Heatmap splits
  if (splits_flag == TRUE) {
    splits = cumsum(table(annotation_for_heatmap[, algorithm]))
    splits = splits[-length(splits)]
    splits = rep(1:(length(splits)+1), c(table(annotation_for_heatmap[, algorithm])))
    
    # Create Heatmap object
    heatmap <- Heatmap(matrix[order, order],
                       column_title = heatmap_title,
                       column_title_gp = gpar(fontsize = 12, fontface = "bold"),
                       col = color_fun,
                       top_annotation = ha,
                       cluster_rows = cluster_rows_flag,
                       cluster_columns = cluster_cols_flag,
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
  } else {
    # Create Heatmap object
    heatmap <- Heatmap(matrix[order, order],
                       column_title = heatmap_title,
                       column_title_gp = gpar(fontsize = 12, fontface = "bold"),
                       col = color_fun,
                       top_annotation = ha,
                       cluster_rows = cluster_rows_flag,
                       cluster_columns = cluster_cols_flag,
                       show_row_names = FALSE,
                       show_column_names = FALSE,
                       show_row_dend = F,
                       show_column_dend = F,
                       #clustering_distance_rows = "euclidean",
                       #clustering_distance_columns = "euclidean",
                       #row_split = splits,
                       #column_split = splits,
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
}

# PCA from original matrix #####
pca_from_original_matrix = function (mydata = NULL,
                                     algorithm = NULL, clust_res = NULL,
                                     cluster_colors = NULL, output_path = NULL,
                                     title_add = NULL) {
  # Filter the scores dataset for annotation
  plot_df = clust_res
  rownames(plot_df) = plot_df$samID
  
  # Convert character variables to factors
  plot_df[, algorithm] = factor(plot_df[, algorithm],
                                levels = sort(unique(plot_df[, algorithm])),
                                labels = c(paste0(algorithm, c(1:length(unique(plot_df[, algorithm]))))))
  
  labels = plot_df[, algorithm]
  names(labels) = labels[colnames(mydata)]
  
  # Create a list of consistent colours
  annColors = list()
  names(cluster_colors) = sort(unique(plot_df[, algorithm]))
  annColors[[algorithm]] = cluster_colors
  
  n_clust = length(unique(plot_df[, algorithm]))
  
  # Shape determination
  if (n_clust <= 6) {
    shapes = c(15:(15 + n_clust - 1))
  } else {
    shapes = 16
  }
  
  ttt = M3C::pca(mydata = mydata,
                 labels = labels,
                 legendtextsize = 2, legendtitle = "CIMLR Subtype", axistextsize = 4,
                 dotsize = 0.4)+
    aes(shape = as.character(labels), 
        size = as.character(labels),
        alpha = as.character(labels), 
        color = as.character(labels)) +
    scale_size_manual(name = algorithm, 
                      values = rep(0.4, n_clust), 
                      labels = sort(unique(plot_df[, algorithm]))) +
    scale_alpha_manual(name = algorithm, 
                       values = rep(0.7, n_clust), 
                       labels = sort(unique(plot_df[, algorithm]))) +
    scale_shape_manual(name = algorithm, 
                       values = shapes, 
                       labels = sort(unique(plot_df[, algorithm])))+
    scale_color_manual(name = algorithm, 
                       limits = names(annColors[[algorithm]]),
                       values = annColors[[algorithm]], 
                       labels = names(annColors[[algorithm]]))+
    theme_bw()+
    theme(panel.grid.minor = ggplot2::element_blank(),
          panel.grid.major = ggplot2::element_blank(),
          panel.border = element_rect(linewidth = 0.2),
          plot.title = element_text(size = 5, face = "bold"),
          axis.title.x = element_text(size = 4, face = "bold"),
          axis.title.y = element_text(size = 4, face = "bold"),
          axis.text = element_text(size = 4),
          axis.ticks = element_line(linewidth = 0.15),
          legend.background = element_rect(fill = "white", linetype = "solid"),
          #legend.position = c(0.90, 0.86),
          legend.key.size = unit(0.5, "lines"),
          legend.margin = ggplot2::margin(0, 0, 0, 0, unit = "mm"),
          legend.spacing.y = unit(0.1, units = "cm"), 
          legend.title = ggplot2::element_text(size = 4, face = "bold"), 
          legend.text = ggplot2::element_text(size = 3))+
    labs(title = paste0(algorithm, " PCA: clusters and ", title_add),
         legend = algorithm) +
    guides(size = "none", alpha = "none")
  print(ttt)
  ggsave(filename = paste0(algorithm, "_", title_add, "_PCA.png"),
         path = output_path, 
         width = 1920, height = 1080, device = 'png', units = "px",
         dpi = 700)
  dev.off()
  
  rm(ttt, labels, clust_res, annColors, n_clust, mydata, algorithm,
     output_path, title_add, shapes, plot_df); gc()
}

# kernel PCA #####
pca_from_sim_matrix = function (sim_matrix = NULL, algorithm = NULL, clust_res = NULL,
                                cluster_colors = NULL, output_path = NULL,
                                title_add = NULL) {
  
  # kernel PCA of the clusters (Source code modification from kernel_pca() from Spectrum)
  km = sim_matrix
  m = nrow(km)
  kc = t(t(km - colSums(km)/m) - rowSums(km)/m) + sum(km)/m^2
  res = eigen(kc/m, symmetric = TRUE)
  features = m
  ret = suppressWarnings(t(t(res$vectors[, 1:features])/sqrt(res$values[1:features])))
  scores = data.frame(ret)
  rownames(scores) = colnames(km)
  colnames(scores)[1:2] = c("PC1", "PC2")
  scores$samID = rownames(scores)
  
  # Filter the scores dataset for annotation
  plot_df = scores %>% inner_join(clust_res, by = "samID")
  rownames(plot_df) = plot_df$samID
  
  # Convert character variables to factors
  plot_df[, algorithm] = factor(plot_df[, algorithm],
                                levels = sort(unique(plot_df[, algorithm])),
                                labels = c(paste0(algorithm, c(1:length(unique(plot_df[, algorithm]))))))
  
  # Create a list of consistent colours
  annColors = list()
  names(cluster_colors) = sort(unique(plot_df[, algorithm]))
  annColors[[algorithm]] = cluster_colors
  
  n_clust = length(unique(plot_df[, algorithm]))
  
  # Shape determination
  if (n_clust <= 6) {
    shapes = c(15:(15 + n_clust - 1))
  } else {
    shapes = 16
  }
  
  # Actual plot
  kernelPCA = ggplot2::ggplot(data = plot_df, aes(x = PC1, y = PC2)) +
    ggplot2::geom_point(
      aes(shape = as.character(!!sym(algorithm)), size = as.character(!!sym(algorithm)),
          alpha = as.character(!!sym(algorithm)), color = as.character(!!sym(algorithm)))) +
    scale_size_manual(name = algorithm, 
                      values = rep(0.2, n_clust), 
                      labels = sort(unique(plot_df[, algorithm]))) +
    scale_alpha_manual(name = algorithm, 
                       values = rep(0.7, n_clust), 
                       labels = sort(unique(plot_df[, algorithm]))) +
    scale_shape_manual(name = algorithm, 
                       values = shapes, 
                       labels = sort(unique(plot_df[, algorithm])))+
    scale_color_manual(name = algorithm, 
                       limits = names(annColors[[algorithm]]),
                       values = annColors[[algorithm]], 
                       labels = names(annColors[[algorithm]]))+
    # geom_segment(aes(x = 18, y = 11, xend = 18, yend = 18),
    #              alpha = 0.7, color = "grey", linewidth = 0.1)+
    # geom_segment(aes(x = 18, y = 11, xend = 23, yend = 11),
    #              alpha = 0.7, color = "grey", linewidth = 0.1)+
    theme_bw()+
    theme(panel.grid.minor = ggplot2::element_blank(),
          panel.grid.major = ggplot2::element_blank(),
          panel.border = element_rect(linewidth = 0.2),
          plot.title = element_text(size = 5, face = "bold"),
          axis.title.x = element_text(size = 4, face = "bold"),
          axis.title.y = element_text(size = 4, face = "bold"),
          axis.text = element_text(size = 4),
          axis.ticks = element_line(linewidth = 0.15),
          legend.background = element_rect(fill = "white", linetype = "solid"),
          #legend.position = c(0.90, 0.86),
          legend.key.size = unit(0.5, "lines"),
          legend.margin = ggplot2::margin(0, 0, 0, 0, unit = "mm"),
          legend.spacing.y = unit(0.1, units = "cm"), 
          legend.title = ggplot2::element_text(size = 4, face = "bold"), 
          legend.text = ggplot2::element_text(size = 3))+
    labs(title = paste0(algorithm, " Kernel PCA: clusters and ", title_add),
         legend = algorithm) +
    guides(size = "none", alpha = "none")
  print(kernelPCA)
  ggsave(filename = paste0(algorithm, "_", title_add, "_kernelPCA.png"),
         path = output_path, 
         width = 1920, height = 1080, device = 'png', units = "px",
         dpi = 700)
  dev.off()
  
  # Clean up
  rm(km, m, kc, res, ret, features, plot_df, annColors, algorithm, kernelPCA,
     scores, shapes, n_clust, cluster_colors, output_path, title_add); gc()
}

# create single barchart #####
create_annot_barchart = function (plotdata = NULL, fill = NULL,
                                  chifit = NULL, algorithm = NULL,
                                  text_y = NULL, rect_ymin = NULL,
                                  rect_ymax = NULL, x_annot = NULL,
                                  v_gap = NULL, rect_xmin = NULL,
                                  rect_xmax = NULL, annot_text_size = NULL,
                                  legend.text.size = NULL) {
  library(ggplot2)
  
  barchart = ggplot(plotdata, aes(fill=!!sym(fill), x=!!sym(algorithm))) + 
    geom_bar(position="stack", stat="count", width = 0.4) +
    scale_y_continuous(limits = c(0,nrow(plotdata)), 
                       breaks = seq(0, nrow(plotdata), 50)) +
    labs(y = "Number of samples") +
    ggplot2::annotate("text", x = x_annot, y = text_y, size = annot_text_size, 
                      label = bquote(italic(X^2) == .(round(as.numeric(chifit$Statistic), 2))))+
    ggplot2::annotate("text", x = x_annot, y = text_y-v_gap, size = annot_text_size,
                      parse = if (as.numeric(chifit$`p-value`) < 10e-5) { TRUE } else { FALSE },
                      label = if (as.numeric(chifit$`p-value`) < 10e-5) {
                        "italic(p) < 10^-5"
                      } else {
                        bquote(italic(p) == .(round(as.numeric(chifit$`p-value`), 4)))
                      }) +
    ggplot2::annotate("text", x = x_annot, y = text_y-2*v_gap, size = annot_text_size, 
                      label = bquote(italic(V) == .(round(as.numeric(chifit$`Cramer's V`), 4))))+
    geom_rect(aes(xmin = rect_xmin, xmax = rect_xmax, ymin = rect_ymin, ymax = rect_ymax),
              fill = "transparent", color = "black", linewidth = 0.4) +
    theme(panel.background = element_rect(fill = "white", 
                                          colour = "white"),
          panel.grid = element_blank(),
          axis.line = element_line(),
          axis.title = element_text(size = 8, face = "bold"),
          legend.key.size = unit(0.25, "cm"),
          legend.text = element_text(size = legend.text.size),
          legend.title = element_text(face = "bold", size  = legend.text.size))
  return(barchart)
}

# Sunburst plot function
as.sunburstDF = function(DF, value_column = NULL, add_root = FALSE){
  require(data.table)
  
  colNamesDF = names(DF)
  
  if(is.data.table(DF)){
    DT = copy(DF)
  } else {
    DT = data.table(DF, stringsAsFactors = FALSE)
  }
  
  if(add_root){
    DT[, root := "Total"]  
  }
  
  colNamesDT = names(DT)
  hierarchy_columns = setdiff(colNamesDT, value_column)
  DT[, (hierarchy_columns) := lapply(.SD, as.factor), .SDcols = hierarchy_columns]
  
  if(is.null(value_column) && add_root){
    setcolorder(DT, c("root", colNamesDF))
  } else if(!is.null(value_column) && !add_root) {
    setnames(DT, value_column, "values", skip_absent=TRUE)
    setcolorder(DT, c(setdiff(colNamesDF, value_column), "values"))
  } else if(!is.null(value_column) && add_root) {
    setnames(DT, value_column, "values", skip_absent=TRUE)
    setcolorder(DT, c("root", setdiff(colNamesDF, value_column), "values"))
  }
  
  hierarchyList = list()
  
  for(i in seq_along(hierarchy_columns)){
    current_columns = colNamesDT[1:i]
    if(is.null(value_column)){
      currentDT = unique(DT[, ..current_columns][, values := .N, by = current_columns], by = current_columns)
    } else {
      currentDT = DT[, lapply(.SD, sum, na.rm = TRUE), by=current_columns, .SDcols = "values"]
    }
    setnames(currentDT, length(current_columns), "labels")
    hierarchyList[[i]] = currentDT
  }
  
  hierarchyDT = rbindlist(hierarchyList, use.names = TRUE, fill = TRUE)
  
  parent_columns = setdiff(names(hierarchyDT), c("labels", "values", value_column))
  hierarchyDT[, parents := apply(.SD, 1, function(x){fifelse(all(is.na(x)), yes = NA_character_, no = paste(x[!is.na(x)], sep = ":", collapse = " - "))}), .SDcols = parent_columns]
  hierarchyDT[, ids := apply(.SD, 1, function(x){paste(x[!is.na(x)], collapse = " - ")}), .SDcols = c("parents", "labels")]
  hierarchyDT[, c(parent_columns) := NULL]
  return(hierarchyDT)
}

# Calculate S matrix (neighborhoods) from final affinity matrix #####
calculate_S <- function(W) {
  n <- nrow(W)  # Assuming W is a square matrix
  S <- matrix(0, n, n)  # Initialize S as a zero matrix of the same size as W
  
  # Preserve row and column names
  rownames(S) <- rownames(W)
  colnames(S) <- colnames(W)
  
  for (i in 1:n) {
    # Find the indices of the 30 highest values in W[i, ]
    nn_indices <- order(W[i, ], decreasing = TRUE)[1:30]
    
    # Set S[i, nn_indices] to W[i, nn_indices]
    S[i, nn_indices] <- W[i, nn_indices]
  }
  
  return(S)
}
