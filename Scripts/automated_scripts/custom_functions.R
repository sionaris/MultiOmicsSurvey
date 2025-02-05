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

# Function to reshape Pearson similarity values for tests #####
reshape_Pearson_for_tests <- function(similarities) {
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

# Density curve function #####
create_density_curve <- function(matrix) {
  as.data.frame(matrix) %>%
    gather(.) %>%
    ggplot(aes(x = value, fill = "Density Curve")) + 
    geom_density(alpha = 0.5, position = "identity") +
    xlim(-1, NA) +
    labs(title = "Distribution of values by sample",
         x = "Value", y = "Density") +
    theme_classic()+
    theme(legend.position = "none")
}

# Per column density plot #####
create_density_plot_color <- function(matrix) {
  as.data.frame(matrix) %>%
    gather(., key = "sample", value = "value") %>%
    ggplot(aes(x = value, fill = sample)) + 
    geom_density(alpha = 0.5, position = "identity") +
    xlim(-1, NA) +
    labs(title = "Distribution of values by sample",
         x = "Value", y = "Density") +
    theme_classic()+
    theme(legend.position = "none")
}

# Function to standardize rows #####
standardize_rows <- function(mat) {
  # Apply standardization (subtract mean and divide by standard deviation) to each row
  t(apply(mat, 1, function(row) {
    (row - mean(row, na.rm = TRUE)) / sd(row, na.rm = TRUE)
  }))
}

# getMoHeatmap for prenormalized data with custom quartile limits for colors #####
getMoHeatmap_prenorm_quart <- function (data = NULL, is.binary = c(FALSE, FALSE, FALSE, FALSE, 
                                                                   FALSE, FALSE), row.title = c("Data1", "Data2", "Data3", 
                                                                                                "Data4", "Data5", "Data6"), legend.name = c("Data1", "Data2", 
                                                                                                                                            "Data3", "Data4", "Data5", "Data6"), clust.res = NULL, clust.dend = NULL, 
                                        show.col.dend = TRUE, show.colnames = FALSE, show.row.dend = c(TRUE, 
                                                                                                       TRUE, TRUE, TRUE, TRUE, TRUE), show.rownames = c(FALSE, 
                                                                                                                                                        FALSE, FALSE, FALSE, FALSE, FALSE), clust.dist.row = c("pearson", 
                                                                                                                                                                                                               "pearson", "pearson", "pearson", "pearson", "pearson"), 
                                        clust.method.row = c("ward.D", "ward.D", "ward.D", "ward.D", 
                                                             "ward.D", "ward.D"), clust.col = c("#2EC4B6", "#E71D36", 
                                                                                                "#FF9F1C", "#BDD5EA", "#FFA5AB", "#011627", "#023E8A", 
                                                                                                "#9D4EDD"), color = rep(list(c("#00FF00", "#000000", 
                                                                                                                               "#FF0000")), length(data)), annCol = NULL, annColors = NULL, 
                                        annRow = NULL, width = 6, height = 4, fig.path = getwd(), 
                                        fig.name = "moheatmap", lim_col = 1) 
{
  ht_opt$message = FALSE
  defaultW <- getOption("warn")
  options(warn = -1)
  if (is.null(names(data))) {
    names(data) <- sprintf("dat%s", 1:length(data))
  }
  n_dat <- length(data)
  if (n_dat > 6) {
    stop("current version of MOVICS can support up to 6 datasets.")
  }
  if (n_dat < 2) {
    stop("current version of MOVICS needs at least 2 omics data.")
  }
  colvec <- clust.col[1:length(unique(clust.res$clust))]
  names(colvec) <- paste0("CS", unique(clust.res$clust))
  if (!is.null(annCol) & !is.null(annColors)) {
    annCol <- annCol[colnames(data[[1]]), , drop = FALSE]
    annCol$Subtype <- paste0("CS", clust.res[colnames(data[[1]]), 
                                             "clust"])
    annColors[["Subtype"]] <- colvec
    if (is.null(clust.dend)) {
      clust.res <- clust.res[order(clust.res$clust), ]
      annCol <- annCol[clust.res$samID, , drop = FALSE]
    }
    ha <- ComplexHeatmap::HeatmapAnnotation(df = annCol, 
                                            col = annColors, border = FALSE)
  }
  else {
    annCol <- data.frame(Subtype = paste0("CS", clust.res[colnames(data[[1]]), 
                                                          "clust"]), row.names = colnames(data[[1]]), stringsAsFactors = FALSE)
    annColors <- list(Subtype = colvec)
    if (is.null(clust.dend)) {
      clust.res <- clust.res[order(clust.res$clust), ]
      annCol <- annCol[clust.res$samID, , drop = FALSE]
    }
    ha <- ComplexHeatmap::HeatmapAnnotation(df = annCol, 
                                            col = annColors, border = FALSE)
  }
  if (!is.null(annRow)) {
    if (!is.list(annRow)) {
      stop("argument of annRow should be a list!")
    }
  }
  ht <- list()
  for (i in 1:n_dat) {
    hcg <- hclust(ClassDiscovery::distanceMatrix(as.matrix(t(data[[i]])), 
                                                 clust.dist.row[i]), clust.method.row[i])
    if (is.null(annRow[[i]][1])) {
      rowlab <- ""
      rowlab.index <- 0
    }
    else if (is.na(annRow[[i]][1])) {
      rowlab <- ""
      rowlab.index <- 0
    }
    else {
      rowlab <- intersect(rownames(data[[i]]), annRow[[i]])
      rowlab.index <- match(rowlab, rownames(data[[i]]))
    }
    
    # Calculate the 1st and 3rd quartiles for color mapping
    q1 <- summary(as.vector(data[[i]]))[["1st Qu."]] -lim_col # add some room 
    q3 <- summary(as.vector(data[[i]]))[["3rd Qu."]] +lim_col # add some room
    
    if (is.null(clust.dend)) {
      data <- lapply(data, function(x) x[, clust.res$samID])
      if (!is.binary[i]) {
        col_fun = circlize::colorRamp2(c(q1, q3), c(color[[i]][1], color[[i]][3]))
        ht[[i]] <- ComplexHeatmap::Heatmap(matrix = as.matrix(data[[i]]), 
                                           row_title = row.title[i], name = legend.name[i], 
                                           cluster_columns = FALSE, cluster_rows = hcg, 
                                           show_column_dend = FALSE, show_column_names = show.colnames, 
                                           show_row_dend = show.row.dend[i], show_row_names = show.rownames[i], 
                                           col = col_fun, 
                                           top_annotation = switch((i == 1) + 1, NULL, 
                                                                   ha), width = grid::unit(width, "cm"), height = grid::unit(height, 
                                                                                                                             "cm"), heatmap_legend_param = list(at = pretty(range(data[[i]])), 
                                                                                                                                                                labels = pretty(range(data[[i]]))), right_annotation = ComplexHeatmap::rowAnnotation(link = anno_mark(at = rowlab.index, 
                                                                                                                                                                                                                                                                      labels = rowlab, which = "row", lines_gp = grid::gpar(fontsize = 5), 
                                                                                                                                                                                                                                                                      link_width = grid::unit(3, "mm"), padding = grid::unit(0.8, 
                                                                                                                                                                                                                                                                                                                             "mm"), labels_gp = grid::gpar(fontsize = 7))))
      }
      else {
        col_fun = circlize::colorRamp2(c(0, 1), color[[i]])
        ht[[i]] <- ComplexHeatmap::Heatmap(matrix = as.matrix(data[[i]]), 
                                           row_title = row.title[i], name = legend.name[i], 
                                           cluster_columns = FALSE, cluster_rows = hcg, 
                                           show_column_dend = FALSE, show_column_names = show.colnames, 
                                           show_row_dend = show.row.dend[i], show_row_names = show.rownames[i], 
                                           col = col_fun, top_annotation = switch((i == 
                                                                                     1) + 1, NULL, ha), width = grid::unit(width, 
                                                                                                                           "cm"), height = grid::unit(height, "cm"), 
                                           heatmap_legend_param = list(at = c(0, 1), 
                                                                       legend_gp = grid::gpar(fill = col_fun(c(0, 
                                                                                                               1))), labels = c("0", "1")), right_annotation = ComplexHeatmap::rowAnnotation(link = anno_mark(at = rowlab.index, 
                                                                                                                                                                                                              labels = rowlab, which = "row", lines_gp = grid::gpar(fontsize = 5), 
                                                                                                                                                                                                              link_width = grid::unit(3, "mm"), padding = grid::unit(0.8, 
                                                                                                                                                                                                                                                                     "mm"), labels_gp = grid::gpar(fontsize = 7))))
      }
    }
    else {
      if (!is.binary[i]) {
        col_fun = circlize::colorRamp2(c(q1, q3), c(color[[i]][1], color[[i]][3]))
        ht[[i]] <- ComplexHeatmap::Heatmap(matrix = as.matrix(data[[i]]), 
                                           row_title = row.title[i], name = legend.name[i], 
                                           cluster_columns = clust.dend, cluster_rows = hcg, 
                                           show_column_dend = show.col.dend, show_column_names = show.colnames, 
                                           show_row_dend = show.row.dend[i], show_row_names = show.rownames[i], 
                                           col = col_fun, 
                                           top_annotation = switch((i == 1) + 1, NULL, 
                                                                   ha), width = grid::unit(width, "cm"), height = grid::unit(height, 
                                                                                                                             "cm"), heatmap_legend_param = list(at = pretty(range(data[[i]])), 
                                                                                                                                                                labels = pretty(range(data[[i]]))), right_annotation = ComplexHeatmap::rowAnnotation(link = anno_mark(at = rowlab.index, 
                                                                                                                                                                                                                                                                      labels = rowlab, which = "row", lines_gp = grid::gpar(fontsize = 5), 
                                                                                                                                                                                                                                                                      link_width = grid::unit(3, "mm"), padding = grid::unit(0.8, 
                                                                                                                                                                                                                                                                                                                             "mm"), labels_gp = grid::gpar(fontsize = 7))))
      }
      else {
        col_fun = circlize::colorRamp2(c(0, 1), color[[i]])
        ht[[i]] <- ComplexHeatmap::Heatmap(matrix = as.matrix(data[[i]]), 
                                           row_title = row.title[i], name = legend.name[i], 
                                           cluster_columns = clust.dend, cluster_rows = hcg, 
                                           show_column_dend = show.col.dend, show_column_names = show.colnames, 
                                           show_row_dend = show.row.dend[i], show_row_names = show.rownames[i], 
                                           col = col_fun, top_annotation = switch((i == 
                                                                                     1) + 1, NULL, ha), width = grid::unit(width, 
                                                                                                                           "cm"), height = grid::unit(height, "cm"), 
                                           heatmap_legend_param = list(at = c(0, 1), 
                                                                       legend_gp = grid::gpar(fill = col_fun(c(0, 
                                                                                                               1))), labels = c("0", "1")), right_annotation = ComplexHeatmap::rowAnnotation(link = anno_mark(at = rowlab.index, 
                                                                                                                                                                                                              labels = rowlab, which = "row", lines_gp = grid::gpar(fontsize = 5), 
                                                                                                                                                                                                              link_width = grid::unit(3, "mm"), padding = grid::unit(0.8, 
                                                                                                                                                                                                                                                                     "mm"), labels_gp = grid::gpar(fontsize = 7))))
      }
    }
  }
  if (n_dat == 1) {
    ht_list <- ht[[1]]
  }
  if (n_dat == 2) {
    ht_list <- ht[[1]] %v% ht[[2]]
  }
  if (n_dat == 3) {
    ht_list <- ht[[1]] %v% ht[[2]] %v% ht[[3]]
  }
  if (n_dat == 4) {
    ht_list <- ht[[1]] %v% ht[[2]] %v% ht[[3]] %v% ht[[4]]
  }
  if (n_dat == 5) {
    ht_list <- ht[[1]] %v% ht[[2]] %v% ht[[3]] %v% ht[[4]] %v% 
      ht[[5]]
  }
  if (n_dat == 6) {
    ht_list <- ht[[1]] %v% ht[[2]] %v% ht[[3]] %v% ht[[4]] %v% 
      ht[[5]] %v% ht[[6]]
  }
  outFile <- file.path(fig.path, paste0(fig.name, ".pdf"))
  if (is.null(annCol)) {
    pdf(outFile, width = width, height = height * n_dat/2)
  }
  else {
    pdf(outFile, width = width, height = height * n_dat/1.5)
  }
  draw(ht_list, merge_legend = TRUE, heatmap_legend_side = "right")
  invisible(dev.off())
  draw(ht_list, merge_legend = TRUE, heatmap_legend_side = "right")
  options(warn = defaultW)
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
  
  rownames(annotation_for_heatmap) = annotation_for_heatmap$samID
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
       breaks, afh_colnames, output_file_name, colors, cluster_colors,
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
    shape_guide = "legend"
  } else {
    shapes = rep(16, n_clust)
    shape_guide = "none"
  }
  
  ttt = M3C::pca(mydata = mydata,
                 labels = labels,
                 legendtextsize = 2, legendtitle = paste(algorithm, "Subtype"), axistextsize = 4,
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
                       labels = sort(unique(plot_df[, algorithm])),
                       guide = shape_guide)+
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
    shape_guide = "legend"
  } else {
    shapes = rep(16, n_clust)
    shape_guide = "none"
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
                       labels = sort(unique(plot_df[, algorithm])),
                       guide = shape_guide)+
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

# MDS from original matrix #####
# Designed for matrices with features in rows
mds_from_original_matrix = function (matrix = NULL, dist_method = NULL, algorithm = NULL, clust_res = NULL,
                                 cluster_colors = NULL, output_path = NULL,
                                 title_add = NULL) {
  library(ggplot2)
  
  dists = stats::dist(t(matrix), method = dist_method)
  mds = clust_res %>%
    inner_join(as.data.frame(cmdscale(dists)) %>%
                 mutate(samID = rownames(.)), by = "samID") %>%
    dplyr::rename(MDS1 = V1, MDS2 = V2)
  
  mds[[algorithm]] <- factor(paste0(algorithm, mds[[algorithm]]))
  clust_names = sort(unique(mds[[algorithm]]))
  n_clust = length(clust_names)
  
  # Shape determination
  if (n_clust <= 6) {
    shapes = c(15:(15 + n_clust - 1))
    shape_guide = "legend"
  } else {
    shapes = rep(16, n_clust)
    shape_guide = "none"
  }
  
  mdsplot = ggplot(mds, aes(x = MDS1, y = MDS2, color = !!sym(algorithm), 
                            shape = !!sym(algorithm))) +
    geom_point(size = 0.5, alpha = 0.65) +
    ggtitle(paste0("Multidimensional scaling: ", title_add)) +
    theme_classic() +
    scale_color_manual(name = algorithm, values = cluster_colors,
                       labels = clust_names)+
    scale_shape_manual(name = algorithm, 
                       values = shapes, 
                       labels = clust_names,
                       guide = shape_guide)+
    theme(plot.title = element_text(size = 4, face = "bold", vjust = 0.5, hjust = 0.5),
          axis.text = element_text(size = 3, hjust = 0.5, vjust = 0.5, 
                                   color = "black"),
          axis.title = element_text(size = 4, face = "bold"),
          axis.ticks = element_line(linewidth = 0.05),
          axis.line = element_line(linewidth = 0.2),
          legend.position = "right",
          legend.key.size = unit(2, units = "mm"),
          legend.text = element_text(size = 3),
          legend.title = element_text(face = "bold", size = 3.5),
          legend.margin = ggplot2::margin(0, 0, 0, 0, unit = "mm"),
          legend.spacing.y = unit(0.5, units = "mm"),
          legend.spacing.x = unit(0.5, units = "mm"),
          legend.background = element_blank())+
    labs(x = "MDS1", y = "MDS2") +
    guides(alpha = "none")
  
  print(mdsplot)
  ggsave(filename = paste0(algorithm, "_", title_add, "_MDS.png"),
         path = output_path, 
         width = 1920, height = 1080, device = 'png', units = "px",
         dpi = 700)
  dev.off()
  
  # Clean up
  rm(algorithm, mdsplot, dists, mds,
     clust_names, cluster_colors, output_path, title_add); gc()
}

# create single barchart #####
create_annot_barchart = function (plotdata = NULL, fill = NULL,
                                  na.action = "na.omit",
                                  chifit = NULL, algorithm = NULL,
                                  text_y = NULL, rect_ymin = NULL,
                                  rect_ymax = NULL, x_annot = NULL,
                                  barchart_ylim = NULL,
                                  v_gap = NULL, rect_xmin = NULL,
                                  rect_xmax = NULL, annot_text_size = NULL,
                                  legend.text.size = NULL,
                                  x.axis.text.size = NULL) {
  library(ggplot2)
  
  if (na.action == "na.omit") {
    plotdata = plotdata[-which(is.na(plotdata[, fill])), ]
  } else if (na.action == "keep") {
    plotdata = plotdata
  } else {
    stop("The na.action argument should either be set to 'na.omit' OR 'keep'.")
  }
  
  barchart = ggplot(plotdata, aes(fill=!!sym(fill), x=!!sym(algorithm))) + 
    geom_bar(position="stack", stat="count", width = 0.4) +
    scale_y_continuous(limits = c(0, barchart_ylim), 
                       breaks = seq(0, barchart_ylim, 50)) +
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
          axis.text.x = element_text(size = x.axis.text.size),
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

# Allow non-euclidean distances in CIMLR #####
dist2.cimlr_mod2 = function (x, c = NA, name = NULL, method = "sqeuclidean") 
{
  if (is.na(c)) {
    c = x
  }
  n1 = nrow(x)
  d1 = ncol(x)
  n2 = nrow(c)
  d2 = ncol(c)
  
  if (d1 != d2) {
    stop("Data dimension does not match dimension of centres.")
  }
  
  # Handle different distance metrics
  if (method == "sqeuclidean") {
    # Euclidean distance formula
    dist = t(rep(1, n2) %*% t(apply(t(x^2), MARGIN = 2, FUN = sum))) + 
      (rep(1, n1) %*% t(apply(t(c^2), MARGIN = 2, FUN = sum))) - 
      2 * (x %*% t(c))
  } else if (method == "binary") {
    dist = as.matrix(dist(as.matrix(x), as.matrix(c), method = "binary"))
  } else {
    stop("Unsupported distance method.")
  }
  dist[which(dist < 0, arr.ind = TRUE)] = 0
  # print(paste0(name, ": dim(dist)=", nrow(dist), "x", ncol(dist)))
  return(dist)
}

# Multikernel
multiple.kernel.cimlr_mod = function (x, cores.ratio = 1, name = NULL, method = "sqeuclidean") 
{
  kernel.type = list()
  kernel.type[1] = list("poly")
  kernel.params = list()
  kernel.params[1] = list(0)
  N = dim(x)[1]
  KK = 0
  sigma = seq(2, 1, -0.25)
  Diff = dist2.cimlr_mod2(x, method = method, name = name)
  Diff_sort = t(apply(Diff, MARGIN = 2, FUN = sort))
  m = dim(Diff)[1]
  n = dim(Diff)[2]
  allk = seq(10, 30, 2)
  cores = as.integer(cores.ratio * (detectCores() - 1))
  if (cores < 1 || is.na(cores) || is.null(cores)) {
    cores = 1
  }
  cl = makeCluster(cores)
  clusterEvalQ(cl, {
    library(Matrix)
  })
  D_Kernels = list()
  D_Kernels = unlist(parLapply(cl, 1:length(allk), fun = function(l, 
                                                                  x_fun = x, Diff_sort_fun = Diff_sort, allk_fun = allk, 
                                                                  Diff_fun = Diff, sigma_fun = sigma, KK_fun = KK) {
    if (allk_fun[l] < (nrow(x_fun) - 1)) {
      TT = apply(Diff_sort_fun[, 2:(allk_fun[l] + 1)], 
                 MARGIN = 1, FUN = mean) + .Machine$double.eps
      TT = matrix(data = TT, nrow = length(TT), ncol = 1)
      Sig = apply(array(0, c(nrow(TT), ncol(TT))), MARGIN = 1, 
                  FUN = function(x) {
                    x = TT[, 1]
                  })
      Sig = Sig + t(Sig)
      Sig = Sig/2
      Sig_valid = array(0, c(nrow(Sig), ncol(Sig)))
      Sig_valid[which(Sig > .Machine$double.eps, arr.ind = TRUE)] = 1
      Sig = Sig * Sig_valid + .Machine$double.eps
      for (j in 1:length(sigma_fun)) {
        W = dnorm(Diff_fun, 0, sigma_fun[j] * Sig)
        D_Kernels[[KK_fun + l + j]] = Matrix((W + t(W))/2, 
                                             sparse = TRUE, doDiag = FALSE)
      }
      return(D_Kernels)
    }
  }))
  stopCluster(cl)
  for (i in 1:length(D_Kernels)) {
    K = D_Kernels[[i]]
    k = 1/sqrt(diag(K) + 1)
    G = K * (k %*% t(k))
    G1 = apply(array(0, c(length(diag(G)), length(diag(G)))), 
               MARGIN = 2, FUN = function(x) {
                 x = diag(G)
               })
    G2 = t(G1)
    D_Kernels_tmp = (G1 + G2 - 2 * G)/2
    D_Kernels_tmp = D_Kernels_tmp - diag(diag(D_Kernels_tmp))
    D_Kernels[[i]] = Matrix(D_Kernels_tmp, sparse = TRUE, 
                            doDiag = FALSE)
  }
  return(D_Kernels)
}


# CIMLR
CIMLR_mod = function (X, c, no.dim = NA, k = 10, cores.ratio = 1, binary_flags = NULL,
                      binary_distance = "binary", nonbinary_distance = "sqeuclidean") 
{
  library(CIMLR)
  if (is.na(no.dim)) {
    no.dim = c
  }
  ptm = proc.time()
  NITER = 30
  num = ncol(X[[1]])
  r = -1
  beta = 0.8
  cat("Computing the multiple Kernels.\n")
  for (data_types in 1:length(X)) {
    curr_X = X[[data_types]]
    if (!is.null(binary_flags)) {
      if (binary_flags[data_types] == "Yes") {
        method = binary_distance
      } else {
        method = nonbinary_distance
      }
    } else {
      method = "sqeuclidean"
    }
    if (data_types == 1) {
      D_Kernels = multiple.kernel.cimlr_mod(t(curr_X), method = method, name = names(X)[data_types])
    }
    else {
      D_Kernels = c(D_Kernels, multiple.kernel.cimlr_mod(t(curr_X), method = method, name = names(X)[data_types]))
    }
  }
  alphaK = 1/rep(length(D_Kernels), length(D_Kernels))
  distX = array(0, c(dim(D_Kernels[[1]])[1], dim(D_Kernels[[1]])[2]))
  # for (i in 1:length(D_Kernels)) {
  #   distX = distX + D_Kernels[[i]]
  # }
  for (i in 1:length(D_Kernels)) {
    # print(paste0("Dimensions of distX: ", dim(distX)[1], "x", dim(distX)[2]))
    # print(paste0("Dimensions of D_Kernels[[", i, "]]: ", dim(D_Kernels[[i]])[1], "x", dim(D_Kernels[[i]])[2]))
    
    # Perform the addition
    distX = distX + D_Kernels[[i]]
  }
  
  distX = distX/length(D_Kernels)
  res = apply(distX, MARGIN = 1, FUN = function(x) return(sort(x, 
                                                               index.return = TRUE)))
  distX1 = array(0, c(nrow(distX), ncol(distX)))
  idx = array(0, c(nrow(distX), ncol(distX)))
  for (i in 1:nrow(distX)) {
    distX1[i, ] = res[[i]]$x
    idx[i, ] = res[[i]]$ix
  }
  A = array(0, c(num, num))
  di = distX1[, 2:(k + 2)]
  rr = 0.5 * (k * di[, k + 1] - apply(di[, 1:k], MARGIN = 1, 
                                      FUN = sum))
  id = idx[, 2:(k + 2)]
  numerator = (apply(array(0, c(length(di[, k + 1]), dim(di)[2])), 
                     MARGIN = 2, FUN = function(x) {
                       x = di[, k + 1]
                     }) - di)
  temp = (k * di[, k + 1] - apply(di[, 1:k], MARGIN = 1, FUN = sum) + 
            .Machine$double.eps)
  denominator = apply(array(0, c(length(temp), dim(di)[2])), 
                      MARGIN = 2, FUN = function(x) {
                        x = temp
                      })
  temp = numerator/denominator
  a = apply(array(0, c(length(t(1:num)), dim(di)[2])), MARGIN = 2, 
            FUN = function(x) {
              x = 1:num
            })
  A[cbind(as.vector(a), as.vector(id))] = as.vector(temp)
  if (r <= 0) {
    r = mean(rr)
  }
  lambda = max(mean(rr), 0)
  A[is.nan(A)] = 0
  S0 = max(max(distX)) - distX
  cat("Performing network diffusion.\n")
  S0 = CIMLR:::network.diffusion(S0, k)
  S0 = CIMLR:::dn.cimlr(S0, "ave")
  S = (S0 + t(S0))/2
  D0 = diag(apply(S, MARGIN = 2, FUN = sum))
  L0 = D0 - S
  eig1_res = CIMLR:::eig1(L0, c, 0)
  F_eig1 = eig1_res$eigvec
  temp_eig1 = eig1_res$eigval
  evs_eig1 = eig1_res$eigval_full
  F_eig1 = CIMLR:::dn.cimlr(F_eig1, "ave")
  converge = vector()
  for (iter in 1:NITER) {
    cat("Iteration: ", iter, "\n")
    distf = CIMLR:::L2_distance_1(t(F_eig1), t(F_eig1))
    A = array(0, c(num, num))
    b = idx[, 2:dim(idx)[2]]
    a = apply(array(0, c(num, ncol(b))), MARGIN = 2, FUN = function(x) {
      x = 1:num
    })
    inda = cbind(as.vector(a), as.vector(b))
    ad = (distX[inda] + lambda * distf[inda])/2/r
    dim(ad) = c(num, ncol(b))
    c_input = -t(ad)
    c_output = t(ad)
    ad = t(.Call("projsplx", c_input, c_output))
    A[inda] = as.vector(ad)
    A[is.nan(A)] = 0
    S = (1 - beta) * A + beta * S
    S = CIMLR:::network.diffusion(S, k)
    S = (S + t(S))/2
    D = diag(apply(S, MARGIN = 2, FUN = sum))
    L = D - S
    F_old = F_eig1
    eig1_res = CIMLR:::eig1(L, c, 0)
    F_eig1 = eig1_res$eigvec
    temp_eig1 = eig1_res$eigval
    ev_eig1 = eig1_res$eigval_full
    F_eig1 = CIMLR:::dn.cimlr(F_eig1, "ave")
    F_eig1 = (1 - beta) * F_old + beta * F_eig1
    evs_eig1 = cbind(evs_eig1, ev_eig1)
    DD = vector()
    for (i in 1:length(D_Kernels)) {
      temp = (.Machine$double.eps + D_Kernels[[i]]) * 
        (S + .Machine$double.eps)
      DD[i] = mean(apply(temp, MARGIN = 2, FUN = sum))
    }
    alphaK0 = CIMLR:::umkl.cimlr(DD)
    alphaK0 = alphaK0/sum(alphaK0)
    alphaK = (1 - beta) * alphaK + beta * alphaK0
    alphaK = alphaK/sum(alphaK)
    fn1 = sum(ev_eig1[1:c])
    fn2 = sum(ev_eig1[1:(c + 1)])
    converge[iter] = fn2 - fn1
    if (iter < 10) {
      if (ev_eig1[length(ev_eig1)] > 1e-06) {
        lambda = 1.5 * lambda
        r = r/1.01
      }
    }
    else {
      if (converge[iter] > 1.01 * converge[iter - 1]) {
        S = S_old
        if (converge[iter - 1] > 0.2) {
          warning("Maybe you should set a larger value of c.")
        }
        break
      }
    }
    S_old = S
    distX = D_Kernels[[1]] * alphaK[1]
    for (i in 2:length(D_Kernels)) {
      distX = distX + as.matrix(D_Kernels[[i]]) * alphaK[i]
    }
    res = apply(distX, MARGIN = 1, FUN = function(x) return(sort(x, 
                                                                 index.return = TRUE)))
    distX1 = array(0, c(nrow(distX), ncol(distX)))
    idx = array(0, c(nrow(distX), ncol(distX)))
    for (i in 1:nrow(distX)) {
      distX1[i, ] = res[[i]]$x
      idx[i, ] = res[[i]]$ix
    }
  }
  LF = F_eig1
  D = diag(apply(S, MARGIN = 2, FUN = sum))
  L = D - S
  eigen_L = eigen(L)
  U = eigen_L$vectors
  D = eigen_L$values
  if (length(no.dim) == 1) {
    U_index = seq(ncol(U), (ncol(U) - no.dim + 1))
    F_last = CIMLR:::tsne(S, k = no.dim, initial_config = U[, U_index])
  }
  else {
    F_last = list()
    for (i in 1:length(no.dim)) {
      U_index = seq(ncol(U), (ncol(U) - no.dim[i] + 1))
      F_last[i] = list(CIMLR:::tsne(S, k = no.dim[i], initial_config = U[, 
                                                                 U_index]))
    }
  }
  execution.time = proc.time() - ptm
  cat("Performing Kmeans.\n")
  y = kmeans(F_last, c, nstart = 200)
  ydata = CIMLR:::tsne(S)
  results = list()
  results[["y"]] = y
  results[["S"]] = S
  results[["F"]] = F_last
  results[["ydata"]] = ydata
  results[["alphaK"]] = alphaK
  results[["execution.time"]] = execution.time
  results[["converge"]] = converge
  results[["LF"]] = LF
  return(results)
}

# Function to compute both Frobenius norm and Pearson correlation between matrices #####
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

# Variance and IQR for matrices #####
choose_matrix_contrasts <- function(similarity_matrices) {
  contrast_values <- list()
  
  for (name in names(similarity_matrices)) {
    similarity_values <- as.vector(similarity_matrices[[name]])
    similarity_values <- similarity_values[similarity_values != 1] # remove self-similarities
    
    # Calculate measures of contrast
    variance <- var(similarity_values)
    iqr <- IQR(similarity_values)
    contrast_metric <- variance + iqr
    contrast_values[[name]] <- contrast_metric
  }
  
  return(contrast_values)
}

# Skewness and kurtosis for matrices #####
choose_matrix_skewness_kurtosis <- function(similarity_matrices) {
  skewness_kurtosis_values <- list()
  
  for (name in names(similarity_matrices)) {
    similarity_values <- as.vector(similarity_matrices[[name]])
    similarity_values <- similarity_values[similarity_values != 1] # remove self-similarities
    
    # Calculate skewness and kurtosis
    skewness_value <- skewness(similarity_values)
    kurtosis_value <- kurtosis(similarity_values)
    
    # Calculate a combined metric: |skewness| + kurtosis
    combined_metric <- abs(skewness_value) + kurtosis_value
    skewness_kurtosis_values[[name]] <- combined_metric
  }
  
  return(skewness_kurtosis_values)
}

# Min-max normalization for matrices similarity inspection #####
# Min-max normalization to [0,1]
minmax_normalize_values <- function(values) {
  min_value <- min(values)
  max_value <- max(values)
  
  normalized_values <- (values - min_value) / (max_value - min_value)
  return(normalized_values)
}

# NEMO modifications #####
nemo.affinity.graph_mod = function (raw.data, k = NA, sigma = 0.5,
                                      binary_flags = rep("No", length(raw.data)),
                                      binary_distance = NULL) 
{
  if (is.na(k)) {
    k = as.numeric(lapply(1:length(raw.data), function(i) round(ncol(raw.data[[i]])/NUM.NEIGHBORS.RATIO)))
  }
  else if (length(k) == 1) {
    k = rep(k, length(raw.data))
  }
  sim.data = lapply(1:length(raw.data), function(i) {
    if (binary_flags[i] == "No") {
      affinityMatrix(dist2(as.matrix(t(raw.data[[i]])), as.matrix(t(raw.data[[i]]))), 
                     k[i], sigma)
    } else {
      affinityMatrix(as.matrix(dist(as.matrix(t(raw.data[[i]])),
                                    as.matrix(t(raw.data[[i]])),
                                    method = binary_distance)), 
                     k[i], sigma)
    }
  })
  affinity.per.omic = lapply(1:length(raw.data), function(i) {
    sim.datum = sim.data[[i]]
    non.sym.knn = apply(sim.datum, 1, function(sim.row) {
      returned.row = sim.row
      threshold = sort(sim.row, decreasing = T)[k[i]]
      returned.row[sim.row < threshold] = 0
      row.sum = sum(returned.row)
      returned.row[sim.row >= threshold] = returned.row[sim.row >= 
                                                          threshold]/row.sum
      return(returned.row)
    })
    sym.knn = non.sym.knn + t(non.sym.knn)
    return(sym.knn)
  })
  patient.names = Reduce(union, lapply(raw.data, colnames))
  num.patients = length(patient.names)
  returned.affinity.matrix = matrix(0, ncol = num.patients, 
                                    nrow = num.patients)
  rownames(returned.affinity.matrix) = patient.names
  colnames(returned.affinity.matrix) = patient.names
  shared.omic.count = matrix(0, ncol = num.patients, nrow = num.patients)
  rownames(shared.omic.count) = patient.names
  colnames(shared.omic.count) = patient.names
  for (j in 1:length(raw.data)) {
    curr.omic.patients = colnames(raw.data[[j]])
    returned.affinity.matrix[curr.omic.patients, curr.omic.patients] = returned.affinity.matrix[curr.omic.patients, 
                                                                                                curr.omic.patients] + affinity.per.omic[[j]][curr.omic.patients, 
                                                                                                                                             curr.omic.patients]
    shared.omic.count[curr.omic.patients, curr.omic.patients] = shared.omic.count[curr.omic.patients, 
                                                                                  curr.omic.patients] + 1
  }
  final.ret = returned.affinity.matrix/shared.omic.count
  lower.tri.ret = final.ret[lower.tri(final.ret)]
  final.ret[shared.omic.count == 0] = mean(lower.tri.ret[!is.na(lower.tri.ret)])
  return(final.ret)
}

# Choose the best similarity matrix #####
# Comprehensive function to evaluate and select the best similarity matrix
# Function to evaluate a single similarity matrix and return a list of results
evaluate_similarity_matrix <- function(matrix, k_isomap = 5) {
  
  # Load required libraries
  library(Matrix)          # For sparse matrices
  library(igraph)          # For graph-theoretical measures like modularity
  library(entropy)         # For entropy calculations
  library(RSpectra)        # For fast eigenvalue decomposition
  library(vegan)           # For geodesic distances (Isomap)
  library(stats)           # For dimensionality reduction (PCA)
  
  # Helper function: Spectral analysis
  spectral_gap <- function(matrix) {
    eig_vals <- eigs_sym(as.matrix(matrix), k = 10, which = "LM")$values
    gap <- diff(eig_vals[1:2]) # Calculate the gap between the first two eigenvalues
    return(gap)
  }
  
  # Helper function: Matrix entropy
  matrix_entropy <- function(matrix) {
    matrix_prob <- matrix / sum(matrix)
    ent <- entropy::entropy(matrix_prob, method = "ML") # Maximum Likelihood Entropy
    return(ent)
  }
  
  # Helper function: Manifold preservation using Isomap (geodesic distance, k-NN approach)
  isomap_preservation <- function(matrix, k) {
    distance_matrix <- as.dist(1 - matrix)  # Use 1 - similarity to compute distance
    isomap_result <- vegan::isomap(dist = distance_matrix, ndim = 2, k = k)
    geodesic_distances <- as.matrix(isomap_result$dist)
    return(sum(geodesic_distances))  # Return sum of geodesic distances (lower is better)
  }
  
  # Helper function: Graph modularity
  graph_modularity <- function(matrix) {
    graph <- igraph::graph.adjacency(as.matrix(matrix), mode = "undirected", weighted = TRUE)
    clusters <- igraph::cluster_fast_greedy(graph)
    modularity <- igraph::modularity(clusters)
    return(modularity)
  }
  
  # Helper function: Degree distribution skewness
  degree_distribution_skew <- function(matrix) {
    graph <- igraph::graph.adjacency(as.matrix(matrix), mode = "undirected", weighted = TRUE)
    degree_values <- igraph::degree(graph)
    skewness <- mean(degree_values)
    return(skewness)
  }
  
  # Perform all evaluations
  spectral <- spectral_gap(matrix)
  ent <- matrix_entropy(matrix)
  iso_preserve <- isomap_preservation(matrix, k = k_isomap)
  mod <- graph_modularity(matrix)
  skew <- degree_distribution_skew(matrix)
  
  # Store the raw results in a list
  raw_results <- list(
    spectral_gap = spectral,
    entropy = ent,
    isomap_preservation = iso_preserve,
    modularity = mod,
    degree_skewness = skew
  )
  
  # Normalize the results for ranking
  normalized_results <- list(
    spectral_gap = scale(spectral),
    entropy = scale(-ent),  # Lower entropy is better, so negate it
    isomap_preservation = scale(-iso_preserve),  # Lower is better
    modularity = scale(mod),
    degree_skewness = scale(-skew)  # Lower is better for skewness
  )
  
  # Aggregate normalized scores: Higher total score indicates better matrix
  total_score <- sum(unlist(normalized_results))
  
  # Return a list with raw metrics, normalized scores, and the total score
  results <- list(
    raw_metrics = raw_results,
    normalized_scores = normalized_results,
    total_score = total_score
  )
  
  return(results)
}

compute_silhouette <- function(cluster_df, similarity_matrix, normalize_matrix = FALSE) {
  library(MOVICS)
  
  # Check if required columns are in cluster_df
  if (!all(c("samID", "Cluster") %in% colnames(cluster_df))) {
    stop("cluster_df must contain 'samID' and 'Cluster' columns.")
  }
  
  # Check if the rownames of the similarity matrix match the samID in cluster_df
  if (!all(rownames(similarity_matrix) %in% cluster_df$samID) || 
      !all(cluster_df$samID %in% rownames(similarity_matrix))) {
    stop("The rownames and colnames of similarity_matrix must match the 'samID' column in cluster_df.")
  }
  
  # Ensure the similarity matrix is symmetric
  similarity_matrix <- (similarity_matrix + t(similarity_matrix)) / 2
  diag(similarity_matrix) <- 0  # Set diagonal to 0 (no self-similarity)
  
  # Optionally normalize the similarity matrix
  if (normalize_matrix == "rowSums") {
    normalize <- function(X) X / rowSums(X)
    similarity_matrix <- normalize(similarity_matrix)
  } else if (normalize_matrix == "minmax") {
    normalize <- function(X) (X - min(X)) / (max(X) - min(X))
    similarity_matrix <- normalize(similarity_matrix)
  } else {
    similarity_matrix <- similarity_matrix
  }
  
  # Check for singleton clusters
  singleton_clusters <- table(cluster_df$Cluster)[table(cluster_df$Cluster) == 1]
  if (length(singleton_clusters) > 0) {
    warning("The following clusters are singletons: ", paste(singleton_clusters, collapse = ", "))
  }
  
  cluster_id <- 1:length(unique(cluster_df$Cluster))
  sil <- matrix(NA, nrow(cluster_df), 3, dimnames = list(cluster_df$samID, 
                                                         c("cluster", "neighbor", "sil_width")))
  for (j in cluster_id) {
    index <- (cluster_df$Cluster == cluster_id[j])
    Nj <- sum(index)
    
    if (Nj == 1) {
      # Handle singleton clusters, assigning silhouette width = 0
      sil[index, "cluster"] <- cluster_id[j]
      sil[index, "neighbor"] <- NA  # No neighbor for singletons
      sil[index, "sil_width"] <- 0  # Silhouette is 0 for singletons
    } else {
      sil[index, "cluster"] <- cluster_id[j]
      dindex <- rbind(apply(similarity_matrix[!index, index, 
                                              drop = FALSE], 2, function(r) tapply(r, cluster_df$Cluster[!index], 
                                                                                   mean)))
      maxC <- apply(dindex, 2, which.max)
      sil[index, "neighbor"] <- cluster_id[-j][maxC]
      a.i <- colSums(similarity_matrix[index, index]) / (Nj - 1)
      b.i <- dindex[cbind(maxC, seq(along = maxC))]
      s.i <- ifelse(a.i != b.i, (a.i - b.i) / pmax(b.i, a.i), 0)
      sil[index, "sil_width"] <- s.i
    }
  }
  
  attr(sil, "Ordered") <- FALSE
  class(sil) <- "silhouette"
  
  return(sil)
}

# Helper function: Linear quantile scaling normalization
linear_quantile_normalize <- function(matrix, norm_quant = 0.05) {
  # Check if norm_quant is within the acceptable range
  if (norm_quant >= 0.5 || norm_quant < 0) {
    stop("norm_quant must be a number between 0 and 0.5")
  }
  
  # Extract off-diagonal values
  off_diag_values <- matrix[upper.tri(matrix, diag = FALSE)]
  
  # Remove NA values if present
  off_diag_values <- na.omit(off_diag_values)
  
  # Ensure there are off-diagonal values to work with
  if (length(off_diag_values) == 0) {
    stop("All off-diagonal values are NA or constant.")
  }
  
  # Step 1: Define the lower and upper quantiles
  m_i <- quantile(off_diag_values, norm_quant / 2, na.rm = TRUE)  # q/2 quantile
  m_a <- quantile(off_diag_values, 1 - norm_quant / 2, na.rm = TRUE)  # 1 - q/2 quantile
  
  # Step 2: Linearly scale each value x_i to x'_i
  scaled_values <- (off_diag_values - m_i) / (m_a - m_i)
  
  # Ensure no NAs were introduced during scaling
  if (any(is.na(scaled_values))) {
    stop("NA values found after scaling the matrix. Check the matrix input.")
  }
  
  # Replace off-diagonal values in the matrix with scaled values
  matrix[upper.tri(matrix, diag = FALSE)] <- scaled_values
  
  # Symmetrize the matrix (since it's a similarity matrix)
  matrix <- (matrix + t(matrix)) / 2
  
  return(matrix)
}

# Helper function: Divide by quantile normalization
divide_by_quantile_normalize <- function(matrix, norm_quant = 0.05) {
  # Check if norm_quant is within the acceptable range
  if (norm_quant >= 0.5 || norm_quant < 0) {
    stop("norm_quant must be a number between 0 and 0.5")
  }
  
  # Extract off-diagonal values
  off_diag_values <- matrix[upper.tri(matrix, diag = FALSE)]
  
  # Remove NA values if present
  off_diag_values <- na.omit(off_diag_values)
  
  # Ensure there are off-diagonal values to work with
  if (length(off_diag_values) == 0) {
    stop("All off-diagonal values are NA or constant.")
  }
  
  # Step 1: Find the chosen quantile
  chosen_quantile <- quantile(off_diag_values, 1 - norm_quant, na.rm = TRUE)
  
  # Step 2: Divide off-diagonal values by the chosen quantile
  scaled_values <- off_diag_values / chosen_quantile
  
  # Ensure no NAs were introduced during scaling
  if (any(is.na(scaled_values))) {
    stop("NA values found after dividing by quantile. Check the matrix input.")
  }
  
  # Replace off-diagonal values in the matrix with scaled values
  matrix[upper.tri(matrix, diag = FALSE)] <- scaled_values
  
  # Symmetrize the matrix (since it's a similarity matrix)
  matrix <- (matrix + t(matrix)) / 2
  
  return(matrix)
}

# Main function: Wrapper for both normalization methods
transform_affinity_matrix <- function(similarity_matrix, norm_quant = 0.05, norm_method = "linear quantile scaling", threshold = 100) {
  # Step 1: Calculate statistics before modifying the matrix
  diag_values <- diag(similarity_matrix)  # Diagonal values (self-similarity)
  
  # Get off-diagonal values only
  off_diag_values <- similarity_matrix[upper.tri(similarity_matrix, diag = FALSE)]
  
  # Calculate average and median for off-diagonal values
  avg_off_diag <- mean(off_diag_values, na.rm = TRUE)
  median_off_diag <- median(off_diag_values, na.rm = TRUE)
  
  # Calculate average for diagonal values
  avg_diag <- mean(diag_values, na.rm = TRUE)
  
  # Step 2: Check if transformation makes sense based on the criterion
  if (avg_diag > threshold * avg_off_diag && avg_diag > threshold * median_off_diag) {
    message("The average diagonal value is more than ", threshold, 
            " times the average and median of off-diagonal values. Transformation is recommended.")
  } else {
    warning("The transformation might not be necessary based on the current matrix values.")
  }
  
  # Step 3: Now set the diagonal to 0 (no self-similarity)
  diag(similarity_matrix) <- 0
  
  # Step 4: Apply the chosen normalization method
  if (norm_method == "linear quantile scaling") {
    similarity_matrix <- linear_quantile_normalize(similarity_matrix, norm_quant)
  } else if (norm_method == "divide by quantile") {
    similarity_matrix <- divide_by_quantile_normalize(similarity_matrix, norm_quant)
  } else {
    stop("Invalid norm_method. Choose either 'linear quantile scaling' or 'divide by quantile'.")
  }
  
  return(similarity_matrix)
}

# Rank SNF features parallely #####
rankFeaturesByNMI_parallely <- function(data, W, ncores = detectCores() - 1, binary = FALSE, nn = NULL, sigma) {
  stopifnot(class(data) == "list" && length(data) == 1)  # Ensure only one data type is passed
  
  NUM_OF_FEATURES <- ncol(data[[1]])
  NMI_scores <- vector(mode = "numeric", length = NUM_OF_FEATURES)
  problematic_features <- vector(mode = "list")  # To store indices of problematic features
  num_of_clusters_fused <- estimateNumberOfClustersGivenGraph(W)[[1]]
  clustering_fused <- spectralClustering(W, num_of_clusters_fused)
  
  # Set up parallel backend to use with foreach
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  
  # Ensure cluster is stopped in case of error
  on.exit(stopCluster(cl))
  
  clusterEvalQ(cl, library(SNFtool))  # Make SNFtool available to each core
  clusterEvalQ(cl, library(foreach))
  clusterEvalQ(cl, library(doParallel))
  
  # Use foreach to parallelize over features (compatible with Windows)
  data_type_scores <- foreach(feature_ind = 1:NUM_OF_FEATURES, .combine = 'c', .packages = c("SNFtool")) %dopar% {
    tryCatch({
      if (binary) {
        # Use binary distance
        dist_matrix <- as.matrix(dist(as.matrix(data[[1]][, feature_ind]), method = "binary"))
      } else {
        # Use default distance (assumed to be Euclidean)
        dist_matrix <- dist2(as.matrix(data[[1]][, feature_ind]), as.matrix(data[[1]][, feature_ind]))
      }
      
      affinity_matrix <- affinityMatrix(dist_matrix, K = nn, sigma = sigma)      
      clustering_single_feature <- spectralClustering(affinity_matrix, num_of_clusters_fused)
      calNMI(clustering_fused, clustering_single_feature)
    }, error = function(e) {
      # Log the index of the problematic feature and return NA for its score
      problematic_features <<- append(problematic_features, feature_ind)
      NA
    })
  }
  
  # Rank the features, excluding NA values from ranking
  data_type_ranks <- rank(-data_type_scores, ties.method = "first", na.last = "keep")
  
  # Print or log problematic feature indices
  if (length(problematic_features) > 0) {
    cat("Problematic feature indices:", unlist(problematic_features), "\n")
  }
  
  return(list(NMI_scores = data_type_scores, NMI_ranks = data_type_ranks, problematic_features = problematic_features))
}

# o1-IntNMF update #####
nmf.opt.k.integrative <- function(dat, is.binary, n.runs = 30, n.fold = 5, k.range = 2:8, 
                                  result = TRUE, make.plot = TRUE, progress = TRUE, 
                                  maxiter = 100, lr = 1e-3, tol = 1e-6, 
                                  allowParallel = FALSE, n.cores = NULL, seed = 12345) {
  # Required libraries:
  library(mclust) # for adjustedRandIndex
  library(foreach)
  library(doParallel)
  library(MASS) # for ginv
  
  if (!is.list(dat)) stop("Input 'dat' must be a list of matrices.")
  M <- length(dat)
  if (length(is.binary) != M) stop("Length of 'is.binary' must match the number of modalities in 'dat'.")
  
  for (i in seq_len(M)) {
    if (min(dat[[i]]) < 0) {
      dat[[i]] <- pmax(dat[[i]] + abs(min(dat[[i]])), 0) + .Machine$double.eps
    }
  }
  
  # Calculate weights based on proportions
  M_b <- sum(is.binary)
  M_c <- M - M_b
  if (M_b > 0 && M_c > 0) {
    total_binary_weight <- M_b / M
    total_cont_weight <- M_c / M
    wt <- numeric(M)
    wt[is.binary] <- total_binary_weight / M_b
    wt[!is.binary] <- total_cont_weight / M_c
  } else {
    # All binary or all continuous
    wt <- rep(1/M, M)
  }
  
  sigmoid <- function(x) 1 / (1 + exp(-x))
  
  # Training fit function
  nmf.integrative.fit <- function(dat.list, k, maxiter, lr, tol, wt, is.binary, seed) {
    n <- nrow(dat.list[[1]])
    set.seed(seed)
    W <- matrix(runif(n * k, min = 0, max = 1), n, k)
    H.list <- vector("list", M)
    for (m in seq_len(M)) {
      H.list[[m]] <- matrix(runif(k * ncol(dat.list[[m]]), min = 0, max = 1), k, ncol(dat.list[[m]]))
    }
    
    prev_W <- W
    for (iter in seq_len(maxiter)) {
      # Update H
      for (m in seq_len(M)) {
        X_m <- dat.list[[m]]
        pred_m <- W %*% H.list[[m]]
        if (is.binary[m]) {
          Sigm <- sigmoid(pred_m)
          Grad_H <- wt[m] * t(W) %*% (Sigm - X_m) 
          H.list[[m]] <- pmax(H.list[[m]] - lr * Grad_H, 0)
        } else {
          Grad_H <- wt[m] * t(W) %*% ((W %*% H.list[[m]]) - X_m)
          H.list[[m]] <- pmax(H.list[[m]] - lr * Grad_H, 0)
        }
      }
      
      # Update W
      W_grad <- matrix(0, nrow = n, ncol = k)
      for (m in seq_len(M)) {
        X_m <- dat.list[[m]]
        pred_m <- W %*% H.list[[m]]
        if (is.binary[m]) {
          Sigm <- sigmoid(pred_m)
          W_grad <- W_grad + wt[m] * ((Sigm - X_m) %*% t(H.list[[m]]))
        } else {
          W_grad <- W_grad + wt[m] * (((W %*% H.list[[m]]) - X_m) %*% t(H.list[[m]]))
        }
      }
      W <- pmax(W - lr * W_grad, 0)
      
      rel_change <- sum(abs(W - prev_W)) / (sum(abs(prev_W)) + .Machine$double.eps)
      if (rel_change < tol) break
      prev_W <- W
    }
    
    clusters <- apply(W, 1, which.max)
    list(W = W, H = H.list, clusters = clusters)
  }
  
  # Test W estimation given H from training:
  # We fix H and solve for W by gradient descent, mixing binary and continuous losses.
  nmf.test.W <- function(d.test, H.list, wt, is.binary, lr, tol, maxiter, seed) {
    n.test <- nrow(d.test[[1]])
    k <- nrow(H.list[[1]])
    set.seed(seed+1)
    W_test <- matrix(runif(n.test * k, min=0, max=1), n.test, k)
    prev_W <- W_test
    
    for (iter in seq_len(maxiter)) {
      W_grad <- matrix(0, n.test, k)
      for (m in seq_len(M)) {
        X_m <- d.test[[m]]
        pred_m <- W_test %*% H.list[[m]]
        if (is.binary[m]) {
          Sigm <- 1/(1+exp(-pred_m))
          W_grad <- W_grad + wt[m] * ((Sigm - X_m) %*% t(H.list[[m]]))
        } else {
          # continuous
          W_grad <- W_grad + wt[m] * (((W_test %*% H.list[[m]]) - X_m) %*% t(H.list[[m]]))
        }
      }
      W_test <- pmax(W_test - lr * W_grad, 0)
      
      rel_change <- sum(abs(W_test - prev_W)) / (sum(abs(prev_W)) + .Machine$double.eps)
      if (rel_change < tol) break
      prev_W <- W_test
    }
    W_test
  }
  
  if (allowParallel) {
    if (is.null(n.cores)) n.cores <- parallel::detectCores() - 1
    cl <- parallel::makeCluster(n.cores)
    doParallel::registerDoParallel(cl)
  }
  
  set.seed(seed)
  n.sample <- nrow(dat[[1]])
  CPI <- matrix(NA, length(k.range), n.runs)
  dimnames(CPI) <- list(paste("k", k.range, sep = ""), paste("run", 1:n.runs, sep = ""))
  
  res_list <- foreach::foreach(i = 1:n.runs, .combine = 'cbind', .packages = c("mclust")) %dopar% {
    run_results <- numeric(length(k.range))
    for (ki in seq_along(k.range)) {
      k <- k.range[ki]
      R.ind <- NULL
      random.sample <- sample(seq(n.sample), n.sample)
      fold_size <- floor(n.sample/n.fold)
      
      for (j in 1:n.fold) {
        test_idx <- ((j-1)*fold_size+1):min(j*fold_size, n.sample)
        test.sample <- random.sample[test_idx]
        train.sample <- setdiff(random.sample, test.sample)
        
        d.train <- lapply(dat, function(x) x[train.sample, , drop = FALSE])
        d.test <- lapply(dat, function(x) x[test.sample, , drop = FALSE])
        
        # Train model on training set
        fit.train <- nmf.integrative.fit(d.train, k = k, maxiter = maxiter, lr = lr, tol = tol, wt = wt, is.binary = is.binary, seed = seed)
        
        # Compute test W by fixing H and doing iterative updates
        W.predict <- nmf.test.W(d.test, fit.train$H, wt, is.binary, lr, tol, maxiter, seed)
        predicted.cluster.mem <- apply(W.predict, 1, which.max)
        
        # Fit model on test set for ground truth clusters
        fit.test <- nmf.integrative.fit(d.test, k = k, maxiter = maxiter, lr = lr, tol = tol, wt = wt, is.binary = is.binary, seed = seed+10)
        computed.cluster.mem <- fit.test$clusters
        
        R.ind <- c(R.ind, mclust::adjustedRandIndex(predicted.cluster.mem, computed.cluster.mem))
        
        if (progress) {
          done <- ((i-1)*length(k.range)*n.fold + (ki-1)*n.fold + j) / (n.runs*length(k.range)*n.fold)
          pct <- round(done*100)
          if (pct %in% seq(5,100,by=5)) {
            message(pct, "% complete")
            flush.console()
          }
        }
      }
      run_results[ki] <- mean(R.ind)
    }
    run_results
  }
  
  if (allowParallel) {
    parallel::stopCluster(cl)
  }
  
  for (i in 1:ncol(res_list)) {
    CPI[, i] <- res_list[, i]
  }
  
  if (make.plot) {
    dev.new(width = 4, height = 5)
    plot(k.range, CPI[, 1], ylim = c(min(CPI, na.rm=TRUE), max(CPI, na.rm=TRUE)), pch = 20, main = "", xlab = "k", ylab = "CPI")
    for (m in 2:n.runs) points(k.range, CPI[, m], pch = 20)
    lines(k.range, apply(CPI, 1, mean, na.rm=TRUE), col = "red", lwd = 2)
    mtext("Optimum k", outer = TRUE, cex = 1, line = -2)
  }
  
  if (result) return(CPI)
}

# Modified clusternomics function to avoid errors #####
# Adapted from MONET's official repository:
# https://github.com/Shamir-Lab/MONET/blob/master/R_code/monet_exp.R
run.clusternomics <- function(omics.list, num.clusters=NULL, 
                              num.clusters.per.omic=NULL, dataDistributions = NULL,
                              ncores = NULL) {
  library(clusternomics)
  if (is.null(num.clusters)) {
    stop("Provide a valid non-NULL/NA num.clusters value!")
  }
  if (is.null(num.clusters.per.omic)) {
    stop("Provide a valid non-NULL/NA num.clusters.per.omic value! It must be a list of two vectors.")
  }
  if(is.null(ncores)) {
    stop("ncores must be a numeric value.")
  }
  if (is.null(dataDistributions)) {
    stop("Provide a valid non-NULL/NA value for dataDistributions.")
  } else if (!is.null(dataDistributions) && length(dataDistributions) != length(omics.list)) {
    stop("length(omics.list) == length(dataDistributions) is essential.")
  }
  all.param.options = expand.grid(1:length(num.clusters), 1:length(num.clusters.per.omic))
  all.rets = mclapply(1:nrow(all.param.options), function(i) {
    set.seed(123 + i)
    cur.num.clusters = num.clusters[all.param.options[i, 1]]
    cur.num.clusters.per.omic = num.clusters.per.omic[[all.param.options[i, 2]]]
    start = Sys.time()
    if (length(cur.num.clusters.per.omic) == 1) {
      num.clusters.per.omic = rep(num.clusters.per.omic, length(omics.list))
    }
    cluster.counts = list(global=cur.num.clusters, context=cur.num.clusters.per.omic)
    
    # Hack due to a false assertion in clusternomics package
    if (length(omics.list) > 2) {
      cluster.counts = c(cluster.counts, rep('UNUSED', length(omics.list) - 2))
    }
    # parameters are as used in the clusternomics publication.
    results = contextCluster(omics.list.trans, cluster.counts, maxIter=1e4, 
                             burnin=5e3, lag=3, dataDistributions=dataDistributions, verbose=T)
    cur.clustering = results$samples[[length(results$samples)]]$Global
    dic = results$DIC
    time.taken.per.param = as.numeric(Sys.time() - start, units='secs')
    return(list(clustering=cur.clustering, dic=dic, timing=time.taken.per.param, clusternomics.ret=results))
  }, mc.cores=ncores)
  
  best.sol.index = which.min(sapply(all.rets, function(x) x$dic))
  wall.timing = max(sapply(all.rets, function(x) x$timing)) + time.taken.normalization
  return(list(clustering=all.rets[[best.sol.index]]$clustering, timing=wall.timing, all.rets=all.rets))
}

# Spectrum functions #####
Spectrum_bin_and_par <- function (
    data, 
    method = 1, 
    silent = FALSE, 
    showres = TRUE, 
    diffusion = TRUE, 
    kerneltype = c("density", "stsc"), 
    maxk = 10, 
    NN = 3, 
    NN2 = 7, 
    showpca = FALSE, 
    frac = 2, 
    thresh = 7, 
    fontsize = 18, 
    dotsize = 3, 
    tunekernel = FALSE, 
    clusteralg = "GMM", 
    FASP = FALSE, 
    FASPk = NULL, 
    fixk = NULL, 
    krangemax = 10, 
    runrange = FALSE, 
    diffusion_iters = 4, 
    KNNs_p = 10, 
    missing = FALSE,
    distances = "euclidean",  # Added distances argument
    cores = 1                 # Added cores argument
) 
{
  # Load necessary libraries
  required_packages <- c("foreach", "doParallel", "diptest", "Rfast", "ggplot2")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("The", pkg, "package is required but not installed. Please install it before proceeding."))
    }
  }
  
  library(foreach)
  library(doParallel)
  
  kerneltype <- match.arg(kerneltype)
  
  # Handle the 'distances' argument
  if (length(distances) == 1) {
    distances <- rep(distances, length(data))
  } else if (length(distances) != length(data)) {
    stop("Error: 'distances' must be either a single string or a vector with length equal to the number of data views.")
  }
  
  # Validate that all distance types are supported
  supported_distances <- c("euclidean", "manhattan", "binary", "canberra", "maximum", "minkowski")
  if (!all(distances %in% supported_distances)) {
    stop(paste("Error: Unsupported distance type detected. Supported distances are:", 
               paste(supported_distances, collapse = ", ")))
  }
  
  # Validate 'cores' argument
  if (!is.numeric(cores) || length(cores) != 1 || cores < 1 || cores != floor(cores)) {
    stop("Error: 'cores' must be a positive integer.")
  }
  
  available_cores <- parallel::detectCores()
  if (cores > available_cores) {
    warning(paste("Requested number of cores (", cores, ") exceeds available cores (", available_cores, "). Using ", available_cores, " cores instead.", sep = ""))
    cores <- available_cores
  }
  
  # Convert data to list if not already
  if (!inherits(data, "list")) {
    datalist <- list(data)
    distances <- distances[1]  # Ensure distances aligns with datalist
  }
  else {
    datalist <- data
  }
  
  # Check FASP constraints
  if (length(datalist) > 1 & FASP == TRUE) {
    stop("Error: FASP method works for only a single view")
  }
  if (is.null(FASPk) == TRUE & FASP == TRUE) {
    stop("Error: FASP method requires a number of centroids to compute")
  }
  if (runrange == TRUE & method == 3) {
    stop("Error: cannot run a range of K whilst method=3")
  }
  if (is.null(fixk) == TRUE & method == 3) {
    stop("Error: need to set the value of K using the fixk parameter for method 3")
  }
  
  # Informational messages
  if (silent == FALSE) {
    message("***Spectrum***")
    message(paste("Detected views:", length(datalist)))
    message(paste("Method:", method))
    message(paste("Kernel type:", kerneltype))
    message(paste("Using", cores, "core(s) for parallel processing"))
  }
  
  # FASP data compression
  if (FASP) {
    message("Running with FASP data compression")
    cs <- kmeans(t(datalist[[1]]), centers = FASPk)
    csx <- cs$centers
    cas <- cs$cluster
    datalist[[1]] <- data.frame(t(csx))
  }
  
  # Set up parallel backend for kernel computation
  cl_kernel <- makeCluster(cores)
  registerDoParallel(cl_kernel)
  
  # List of helper functions to export
  helper_functions <- c("CNN_kernel_mod", "kernfinder_mine_mod", "kernfinder_local_mod")
  #,
  # "rbfkernel_b_mod", "EM_finder", "findk", 
  # "plot_egap", "plot_multigap", "pca", 
  # "harmonise_ids", "mean_imputation")
  
  # Parallelized kernel computation using foreach
  kernellist <- foreach(platform = seq_along(datalist), 
                        .packages = c("Spectrum", "Rfast", "ggplot2", "diptest"),
                        .export = helper_functions) %dopar% {
                          if (silent == FALSE) {
                            message(paste("Calculating similarity matrix for modality", platform))
                          }
                          
                          # Retrieve the distance for the current modality
                          current_distance <- distances[platform]
                          
                          # Initialize kerneli
                          kerneli <- NULL
                          
                          # Modify kernel computation based on the specified distance
                          if (kerneltype == "stsc") {
                            if (method == 2 && tunekernel) {
                              NN_current <- kernfinder_local_mod(
                                datalist[[platform]], 
                                maxk = maxk, 
                                silent = silent, 
                                fontsize = fontsize, 
                                dotsize = dotsize, 
                                showres = showres,
                                distance = current_distance  # Passing the specified distance
                              )
                              NN <- NN_current
                            }
                            
                            # Compute the RBF kernel based on the specified distance
                            kerneli <- rbfkernel_b_mod(
                              datalist[[platform]], 
                              K = NN, 
                              sigma = 1, 
                              distance = current_distance  # Passing the specified distance
                            )
                          }
                          else if (kerneltype == "density") {
                            if (method == 2 && tunekernel) {
                              NN_current <- kernfinder_mine_mod(
                                datalist[[platform]], 
                                maxk = maxk, 
                                silent = silent, 
                                showres = showres, 
                                fontsize = fontsize, 
                                dotsize = dotsize,
                                distance = current_distance  # Passing the specified distance
                              )
                              NN <- NN_current
                            }
                            
                            # Compute the CNN kernel based on the specified distance
                            kerneli <- CNN_kernel_mod(
                              datalist[[platform]], 
                              K = NN, 
                              NN2 = NN2, 
                              distance = current_distance  # Passing the specified distance
                            )
                          }
                          
                          if (silent == FALSE) {
                            message("Done.")
                          }
                          
                          return(kerneli)
                        }
  
  # Stop the kernel computation cluster
  stopCluster(cl_kernel)
  registerDoSEQ()
  
  # Handle missing data if required
  if (missing) {
    message("Imputing missing data...")
    kernellist <- harmonise_ids(kernellist)
    kernellist <- mean_imputation(kernellist)
    message("Done.")
  }
  
  # Combine similarity matrices
  if (silent == FALSE) {
    message("Combining similarity matrices and creating kNN graph...")
  }
  A <- Reduce("+", kernellist)
  
  # Diffusion process
  if (diffusion == TRUE) {
    for (col in seq_len(ncol(A))) {
      KNNs <- head(rev(sort(A[, col])), (KNNs_p + 1))
      tokeep <- names(KNNs)
      A[!(rownames(A) %in% tokeep), col] <- 0
    }
    A <- A / rowSums(A)
    if (silent == FALSE) {
      message("Done.")
    }
    if (silent == FALSE) {
      message("Diffusing on tensor product graph...")
    }
    Qt <- A
    im <- matrix(0, ncol = ncol(A), nrow = ncol(A))
    diag(im) <- 1
    for (t in seq_len(diffusion_iters)) {
      Qt <- A %*% Qt %*% t(A) + im
    }
    A2 <- t(Qt)
    if (silent == FALSE) {
      message("Done.")
    }
  }
  else if (diffusion == FALSE) {
    A2 <- A / length(datalist)
  }
  
  # Calculating graph Laplacian (L)
  if (silent == FALSE) {
    message("Calculating graph Laplacian (L)...")
  }
  dv <- 1 / sqrt(rowSums(A2))
  l <- dv * A2 %*% diag(dv)
  
  # Eigen decomposition and selecting optimal K
  if (method == 1) {
    if (silent == FALSE) {
      message("Performing eigen decomposition of L...")
    }
    decomp <- eigen(l)
    if (silent == FALSE) {
      message("Done.")
      message("Examining eigenvalues to select K...")
    }
    evals <- as.numeric(decomp$values)
    diffs <- diff(evals)
    diffs <- diffs[-1]
    if (maxk - 1 < length(diffs)) {
      diffs_subset <- diffs[1:(maxk - 1)]
    } else {
      diffs_subset <- diffs
    }
    optk <- which.max(abs(diffs_subset)) + 1
    if (silent == FALSE) {
      message(paste("Optimal K:", optk))
    }
    nn <- maxk + 1
    d <- data.frame(K = seq_len(nn), evals = evals[1:nn])
    if (showres == TRUE) {
      plot_egap(d, maxk = maxk, dotsize = dotsize, fontsize = fontsize)
    }
  }
  else if (method == 2) {
    if (silent == FALSE) {
      message("Performing eigen decomposition of L...")
    }
    decomp <- eigen(l)
    if (silent == FALSE) {
      message("Done.")
      message("Examining eigenvector distributions to select K...")
    }
    xi <- decomp$vectors[, 1:(maxk + 1)]
    res <- EM_finder(xi, silent = silent)
    d <- data.frame(K = seq_len(maxk + 1), Z = res[1:(maxk + 1), 2])
    if (showres == TRUE) {
      plot_multigap(d, maxk = maxk, dotsize = dotsize, fontsize = fontsize)
    }
    optk <- findk(res, maxk = maxk, frac = frac, thresh = thresh)
    if (silent == FALSE) {
      message(paste("Optimal K:", optk))
    }
  }
  else if (method == 3) {
    decomp <- eigen(l)
    optk <- fixk
  }
  
  results <- list()
  
  # Parallelization setup for clustering over a range of K if runrange is TRUE
  if (runrange) {
    if (silent == FALSE) {
      message("Clustering over a range of K values...")
    }
    
    # Set up parallel backend for clustering
    cl_clustering <- makeCluster(cores)
    registerDoParallel(cl_clustering)
    
    # Parallelized clustering over range of K
    results <- foreach(tk = 2:krangemax, 
                       .packages = c("Spectrum"),  # Replace "ClusterR" with "Spectrum"
                       .export = c(helper_functions, "cas", "d", "findk")) %dopar% {
                         xi <- decomp$vectors[, 1:tk]
                         yi <- xi / sqrt(rowSums(xi^2))
                         yi[!is.finite(yi)] <- 0
                         
                         if (clusteralg == "GMM") {
                           # Assuming Spectrum has GMM and predict_GMM functions
                           gmm <- Spectrum::GMM(
                             yi, 
                             tk, 
                             verbose = FALSE, 
                             seed_mode = "random_spread"
                           )
                           pr <- Spectrum::predict_GMM(
                             yi, 
                             gmm$centroids, 
                             gmm$covariance_matrices, 
                             gmm$weights
                           )
                           names(pr)[3] <- "cluster"
                           if (0 %in% pr$cluster) {
                             pr$cluster <- pr$cluster + 1
                           }
                           if (silent == FALSE) {
                             message(paste("Clustered for K =", tk))
                           }
                         }
                         else if (clusteralg == "km") {
                           pr <- kmeans(yi, tk)
                           if (silent == FALSE) {
                             message(paste("Clustered for K =", tk))
                           }
                         }
                         
                         # Assemble the results as in the original function
                         if (method != 3) {
                           if (FASP) {
                             casn <- cas
                             casn <- casn[seq_along(casn)] <- pr$cluster[as.numeric(casn[seq_along(casn)])]
                             names(casn) <- names(cas)
                             list(
                               allsample_assignments = casn, 
                               centroid_assignments = pr$cluster, 
                               eigenvector_analysis = d, 
                               K = tk, 
                               similarity_matrix = A2, 
                               eigensystem = decomp
                             )
                           }
                           else {
                             list(
                               assignments = pr$cluster, 
                               eigenvector_analysis = d, 
                               K = tk, 
                               similarity_matrix = A2, 
                               eigensystem = decomp
                             )
                           }
                         }
                       }
    
    # Stop the clustering cluster
    stopCluster(cl_clustering)
    registerDoSEQ()
  }
  else {
    xi <- decomp$vectors[, 1:optk]
    yi <- xi / sqrt(rowSums(xi^2))
    yi[!is.finite(yi)] <- 0
    if (clusteralg == "GMM") {
      if (silent == FALSE) {
        message("Performing GMM clustering...")
      }
      # Assuming Spectrum has GMM and predict_GMM functions
      gmm <- Spectrum::GMM(
        yi, 
        optk, 
        verbose = FALSE, 
        seed_mode = "random_spread"
      )
      pr <- Spectrum::predict_GMM(
        yi, 
        gmm$centroids, 
        gmm$covariance_matrices, 
        gmm$weights
      )
      names(pr)[3] <- "cluster"
      if (0 %in% pr$cluster) {
        pr$cluster <- pr$cluster + 1
      }
      if (silent == FALSE) {
        message("Done.")
      }
    }
    else if (clusteralg == "km") {
      if (silent == FALSE) {
        message("Performing k-means clustering...")
      }
      pr <- kmeans(yi, optk)
      if (silent == FALSE) {
        message("Done.")
      }
    }
    if (length(datalist) == 1 && showres == TRUE) {
      if (showpca == TRUE) {
        pca(
          datalist[[1]], 
          labels = as.factor(pr$cluster), 
          axistextsize = fontsize, 
          legendtextsize = fontsize, 
          dotsize = dotsize
        )
      }
    }
    if (method != 3) {
      if (FASP) {
        casn <- cas
        casn <- casn[seq_along(casn)] <- pr$cluster[as.numeric(casn[seq_along(casn)])]
        names(casn) <- names(cas)
        results <- list(
          allsample_assignments = casn, 
          centroid_assignments = pr$cluster, 
          eigenvector_analysis = d, 
          K = optk, 
          similarity_matrix = A2, 
          eigensystem = decomp
        )
      }
      else {
        results <- list(
          assignments = pr$cluster, 
          eigenvector_analysis = d, 
          K = optk, 
          similarity_matrix = A2, 
          eigensystem = decomp
        )
      }
    }
    else if (method == 3) {
      if (FASP) {
        casn <- cas
        casn <- casn[seq_along(casn)] <- pr$cluster[as.numeric(casn[seq_along(casn)])]
        names(casn) <- names(cas)
        results <- list(
          allsample_assignments = casn, 
          centroid_assignments = pr$cluster, 
          K = optk, 
          similarity_matrix = A2, 
          eigensystem = decomp
        )
      }
      else {
        results <- list(
          assignments = pr$cluster, 
          K = optk, 
          similarity_matrix = A2, 
          eigensystem = decomp
        )
      }
    }
  }
  
  if (silent == FALSE) {
    message("Finished.")
  }
  return(results)
}


CNN_kernel_mod <- function(mat, NN = 3, NN2 = 7, distance = "euclidean") 
{
  # Validate distance parameter
  supported_distances <- c("euclidean", "manhattan", "binary", "canberra", "maximum", "minkowski")
  if (!(distance %in% supported_distances)) {
    stop(paste("Unsupported distance type:", distance, ". Supported distances are:", 
               paste(supported_distances, collapse = ", ")))
  }
  
  n <- ncol(mat)
  nbs <- list()
  
  # Compute distance matrix with specified distance metric
  dm <- Rfast::Dist(t(mat), method = distance)
  dimnames(dm) <- list(colnames(mat), colnames(mat))
  
  kn <- c()
  for (i in seq_len(n)) {
    sortedvec <- sort.int(dm[i, ], index.return = FALSE)
    kn <- c(kn, sortedvec[NN + 1])
    nbs[[i]] <- names(sortedvec[2:(NN2 + 1)])
    names(nbs)[[i]] <- names(sortedvec)[1]
  }
  
  sigmamatrix <- kn %o% kn
  out <- matrix(nrow = n, ncol = n)
  upper <- -dm^2
  
  for (i in 2:n) {
    for (j in 1:(i - 1)) {
      cnns <- length(intersect(nbs[[i]], nbs[[j]]))
      upperval <- upper[i, j]
      localsigma <- sigmamatrix[i, j]
      out[i, j] <- exp(upperval / (localsigma * (cnns + 1)))
    }
  }
  
  out <- pmax(out, t(out), na.rm = TRUE)
  diag(out) <- 1
  colnames(out) <- colnames(mat)
  row.names(out) <- colnames(mat)
  
  return(out)
}

kernfinder_mine_mod <- function(data, maxk = 10, fontsize = 12, silent = FALSE, 
                                showres = TRUE, dotsize = 2, distance = "euclidean") 
{
  # Validate distance parameter
  supported_distances <- c("euclidean", "manhattan", "binary", "canberra", "maximum", "minkowski")
  if (!(distance %in% supported_distances)) {
    stop(paste("Unsupported distance type:", distance, ". Supported distances are:", 
               paste(supported_distances, collapse = ", ")))
  }
  
  if (silent == FALSE) {
    message("Finding optimal NN kernel parameter by examining eigenvector distributions")
  }
  
  rr <- c()
  for (param in seq(1, 10)) {
    if (silent == FALSE) {
      message(paste("Tuning kernel NN parameter:", param))
    }
    
    # Pass 'distance' to CNN_kernel
    kern <- CNN_kernel_mod(data, NN = param, NN2 = 7, distance = distance)
    kern[which(!is.finite(kern))] <- 0
    
    dv <- 1 / sqrt(rowSums(kern))
    l <- dv * kern %*% diag(dv)
    xi <- eigen(l)$vectors
    
    res <- matrix(nrow = ncol(xi), ncol = 2)
    for (ii in seq(1, ncol(xi))) {
      r <- diptest::dip.test(xi[, ii], simulate.p.value = FALSE, B = 2000)
      res[ii, 1] <- r$p.value
      res[ii, 2] <- r$statistic
    }
    
    diffs <- diff(res[, 2])
    diffs <- diffs[-1]
    tophit <- diffs[1:(maxk + 1)][which.min(diffs[1:(maxk + 1)])]
    rr <- c(rr, tophit)
  }
  
  optimalparam <- which.min(rr)
  if (silent == FALSE) {
    message(paste("Optimal NN:", optimalparam))
  }
  
  d <- data.frame(x = seq(1, 10), y = rr)
  py <- ggplot2::ggplot(data = d, aes(x = x, y = y)) + 
    ggplot2::geom_point(colour = "black", size = dotsize) + 
    ggplot2::theme_bw() + 
    ggplot2::geom_line() + 
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = fontsize, colour = "black"), 
      axis.text.x = ggplot2::element_text(size = fontsize, colour = "black"), 
      axis.title.x = ggplot2::element_text(size = fontsize), 
      axis.title.y = ggplot2::element_text(size = fontsize), 
      legend.text = ggplot2::element_text(size = fontsize), 
      legend.title = ggplot2::element_text(size = fontsize), 
      plot.title = ggplot2::element_text(size = fontsize, colour = "black", hjust = 0.5), 
      panel.grid.major = ggplot2::element_blank(), 
      panel.grid.minor = ggplot2::element_blank()
    ) + 
    ggplot2::ylab("D") + 
    ggplot2::xlab("NN") + 
    ggplot2::scale_x_continuous(limits = c(1, 10), breaks = seq(1, 10, by = 1))
  
  if (showres == TRUE) {
    print(py)
  }
  
  return(optimalparam)
}

kernfinder_local_mod <- function(data, maxk = 10, fontsize = 12, silent = FALSE, 
                                 showres = TRUE, dotsize = 2, distance = "euclidean") 
{
  # Validate distance parameter
  supported_distances <- c("euclidean", "manhattan", "binary", "canberra", "maximum", "minkowski")
  if (!(distance %in% supported_distances)) {
    stop(paste("Unsupported distance type:", distance, ". Supported distances are:", 
               paste(supported_distances, collapse = ", ")))
  }
  
  if (silent == FALSE) {
    message("Finding optimal kernel NN parameter by examining eigenvectors")
  }
  
  rr <- c()
  for (param in seq(1, 10)) {
    if (silent == FALSE) {
      message(paste("Tuning NN parameter:", param))
    }
    
    # Pass 'distance' to rbfkernel_b
    kern <- rbfkernel_b(data, K = param, sigma = 1, distance = distance)
    
    # Handle non-finite values
    kern[!is.finite(kern)] <- 0
    
    # Compute graph Laplacian
    dv <- 1 / sqrt(rowSums(kern))
    l <- dv * kern %*% diag(dv)
    
    # Eigen decomposition
    xi <- eigen(l)$vectors
    
    # Perform Dip Test on each eigenvector
    res <- matrix(nrow = ncol(xi), ncol = 2)
    for (ii in seq_len(ncol(xi))) {
      r <- diptest::dip.test(xi[, ii], simulate.p.value = FALSE, B = 2000)
      res[ii, 1] <- r$p.value
      res[ii, 2] <- r$statistic
    }
    
    # Calculate differences in dip statistics
    diffs <- diff(res[, 2])
    diffs <- diffs[-1]
    
    # Identify the parameter with the minimum dip statistic difference
    tophit <- diffs[1:(maxk + 1)][which.min(diffs[1:(maxk + 1)])]
    rr <- c(rr, tophit)
  }
  
  # Determine the optimal NN parameter
  optimalparam <- which.min(rr)
  
  if (silent == FALSE) {
    message(paste("Optimal NN:", optimalparam))
  }
  
  # Plot the dip statistic differences
  d <- data.frame(x = seq(1, 10), y = rr)
  py <- ggplot2::ggplot(data = d, aes(x = x, y = y)) + 
    ggplot2::geom_point(colour = "black", size = dotsize) + 
    ggplot2::theme_bw() + 
    ggplot2::geom_line() + 
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = fontsize, colour = "black"), 
      axis.text.x = ggplot2::element_text(size = fontsize, colour = "black"), 
      axis.title.x = ggplot2::element_text(size = fontsize), 
      axis.title.y = ggplot2::element_text(size = fontsize), 
      legend.text = ggplot2::element_text(size = fontsize), 
      legend.title = ggplot2::element_text(size = fontsize), 
      plot.title = ggplot2::element_text(size = fontsize, colour = "black", hjust = 0.5), 
      panel.grid.major = ggplot2::element_blank(), 
      panel.grid.minor = ggplot2::element_blank()
    ) + 
    ggplot2::ylab("D") + 
    ggplot2::xlab("NN") + 
    ggplot2::scale_x_continuous(limits = c(1, 10), breaks = seq(1, 10, by = 1))
  
  if (showres == TRUE) {
    print(py)
  }
  
  return(optimalparam)
}

rbfkernel_b_mod <- function(mat, K = 3, sigma = 1, distance = "euclidean") 
{
  # Validate distance parameter
  supported_distances <- c("euclidean", "manhattan", "binary", "canberra", "maximum", "minkowski")
  if (!(distance %in% supported_distances)) {
    stop(paste("Unsupported distance type:", distance, ". Supported distances are:", 
               paste(supported_distances, collapse = ", ")))
  }
  
  n <- ncol(mat)
  NN <- K
  nbs <- list()
  
  # Compute distance matrix with specified distance metric
  dm <- Rfast::Dist(t(mat), method = distance)
  dimnames(dm) <- list(colnames(mat), colnames(mat))
  
  kn <- c()
  for (i in seq_len(n)) {
    sortedvec <- as.numeric(sort.int(dm[i, ]))
    sortedvec <- sortedvec[!sortedvec == 0]  # Exclude zero distances
    if (length(sortedvec) < NN) {
      stop(paste("Not enough neighbors for sample", colnames(mat)[i], 
                 "with NN =", NN))
    }
    kn <- c(kn, sortedvec[NN])
  }
  
  sigmamatrix <- kn %o% kn
  upper <- -dm^2
  out <- matrix(nrow = n, ncol = n)
  
  for (i in 2:n) {
    for (j in 1:(i - 1)) {
      lowerval <- sigmamatrix[i, j]
      upperval <- upper[i, j]
      out[i, j] <- exp(upperval / (lowerval * sigma))
    }
  }
  
  out <- pmax(out, t(out), na.rm = TRUE)
  diag(out) <- 1
  colnames(out) <- colnames(mat)
  row.names(out) <- colnames(mat)
  
  return(out)
}

# KLIC modification of coca::consensusCluster() #####
coca_cc_mod = function (data = NULL, K = 2, B = 100, pItem = 0.8, clMethod = "hclust", 
                        dist = "euclidean", hclustMethod = "average", sparseKmeansPenalty = NULL, 
                        maxIterKM = 1000) 
{
  library(coca)
  containsFactors <- 0
  if (!is.null(data)) {
    N <- dim(data)[1]
    P <- dim(data)[2]
    for (i in seq_len(P)) {
      containsFactors <- as.numeric(is.factor(data[, i])) + 
        containsFactors
    }
  }
  else if (is.double(dist)) {
    N <- dim(dist)[1]
  }
  else {
    stop("If the data matrix is not provided, `dist` must be a symmetric\n        matrix of type double providing the distances between each pair of\n        observations.")
  }
  dataIndices <- seq_len(N)
  coClusteringMatrix <- indicatorMatrix <- matrix(0, N, N)
  for (b in seq_len(B)) {
    items <- sample(N, ceiling(N * pItem), replace = FALSE)
    nUniqueDataPoints <- 0
    if (!is.null(data)) {
      uniqueData <- unique(data[items, ])
      nUniqueDataPoints <- nrow(uniqueData)
    }
    if (nUniqueDataPoints > K | is.null(data)) {
      if (clMethod == "pam" | containsFactors) {
        if (is.double(dist) & isSymmetric(dist)) {
          distances <- stats::as.dist(dist[items, items])
        }
        else if (dist == "cor") {
          distances <- stats::as.dist(1 - stats::cor(t(data[items, 
          ])))
        }
        else if (dist == "binary") {
          distances <- stats::dist(data[items, ], method = dist)
        }
        else if (dist == "gower") {
          distances <- cluster::daisy(data[items, ], 
                                      metric = "gower")
        }
        else {
          stop("Distance not recognized. If method is `pam`, distance\n                must be one of `cor`, `binary`, `gower` or the symmetric\n                matrix of distances.")
        }
        cl <- cluster::pam(distances, K)$clustering
      }
      else if (clMethod == "kmeans" & !is.null(data)) {
        cl <- stats::kmeans(data[items, ], K, iter.max = maxIterKM, 
                            nstart = 20)$cluster
      }
      else if (clMethod == "sparse-kmeans" & !is.null(data)) {
        if (is.null(sparseKmeansPenalty)) 
          sparseKmeansPenalty = sqrt(P)
        cat("sparseKmeansPenalty", sparseKmeansPenalty, 
            "\n")
        cl <- sparcl::KMeansSparseCluster(data[items, 
        ], K, wbounds = sparseKmeansPenalty)[[1]]$Cs
      }
      else if (clMethod == "hclust" | clMethod == "sparse-hclust") {
        if (is.double(dist)) {
          distances <- stats::as.dist(dist[items, items])
        }
        else if (dist == "pearson" | dist == "spearman") {
          pearsonCor <- stats::cor(t(data[items, ]), 
                                   method = dist)
          distances <- stats::as.dist(1 - pearsonCor)
        }
        else {
          distances <- stats::dist(data[items, ], method = dist)
        }
        if (clMethod == "hclust") {
          hClustering <- stats::hclust(distances, method = hclustMethod)
        }
        else {
          hClustering <- sparcl::HierarchicalSparseCluster(dists = as.matrix(distances), 
                                                           method = "average", wbound = 10)$hc
        }
        cl <- stats::cutree(hClustering, K)
      }
      else {
        stop("Clustering algorithm name not recognised. Please choose\n                     one of `kmeans`, `hclust`, `pam`, `sparse-kmeans`,\n                     `sparse-hclust`.")
      }
      indicatorMatrix <- indicatorMatrix + crossprod(t(as.numeric(dataIndices %in% 
                                                                    items)))
      for (k in seq_len(K)) {
        coClusteringMatrix[items, items] <- coClusteringMatrix[items, 
                                                               items] + crossprod(t(as.numeric(cl == k)))
      }
    }
  }
  if (!sum(indicatorMatrix) == 0) {
    consensusMatrix <- coClusteringMatrix/indicatorMatrix
  }
  else {
    consensusMatrix <- indicatorMatrix
    warning(paste("Consensus matrix is empty for K =", K, 
                  "because there are\n                      less than", 
                  K, "distinct data points", sep = ""))
  }
  return(consensusMatrix)
}

# SNF estimateNUMCfromGraph modification for iterative application #####
library(SNFtool)  # for .discretisation and other SNF functions

estimateNumberOfClustersGivenGraph_mod <- function(W, NUMC = 2:5) 
{
  # Replicate the original check for NUMC == 1
  if (min(NUMC) == 1) {
    warning("Note that we always assume there are more than one cluster.")
    NUMC = NUMC[NUMC > 1]
  }
  
  # Make the affinity matrix symmetric and zero out diagonal
  W = (W + t(W)) / 2
  diag(W) = 0
  
  # Prepare placeholders for the final output
  K1 <- K12 <- K2 <- K22 <- NA
  # We'll also store the eigen-gap and rotation "scores" for the top 2 results
  eigengap_K1_score <- eigengap_K12_score <- NA
  rotation_K2_score <- rotation_K22_score <- NA
  
  if (length(NUMC) > 0) {
    # Degrees and Laplacian construction
    degs = rowSums(W)
    degs[degs == 0] = .Machine$double.eps
    D = diag(degs)
    L = D - W
    Di = diag(1 / sqrt(degs))
    L = Di %*% L %*% Di
    
    # Eigen-decomposition
    eigs = eigen(L)
    eigs_order = sort(eigs$values, index.return = TRUE)$ix
    eigs$values = eigs$values[eigs_order]
    eigs$vectors = eigs$vectors[, eigs_order]
    
    # -------------------------
    # 1. Eigen-gap computation
    # -------------------------
    eigengap = abs(diff(eigs$values))
    eigengap = eigengap * (1 - eigs$values[1:(length(eigs$values) - 1)]) /
      (1 - eigs$values[2:length(eigs$values)])
    
    # We only look at eigengap[k] for k in NUMC
    valid_gap_scores <- eigengap[NUMC]
    
    # Sort them (descending) and get index into NUMC
    gap_sort <- sort(valid_gap_scores, decreasing = TRUE, index.return = TRUE)
    best1_idx <- gap_sort$ix[1]
    best2_idx <- gap_sort$ix[2]
    
    # The top two cluster numbers from the eigen-gap criterion
    K1  = NUMC[best1_idx]
    K12 = NUMC[best2_idx]
    
    # Record their actual gap scores
    eigengap_K1_score  = valid_gap_scores[best1_idx]
    eigengap_K12_score = valid_gap_scores[best2_idx]
    
    # -----------------------------------------
    # 2. Rotation / discretization "quality" 
    # -----------------------------------------
    # In the original code, 'quality' is stored in a list.  
    # Here, we'll store it as a numeric vector of length(NUMC).
    quality <- numeric(length(NUMC))
    for (c_index in seq_along(NUMC)) {
      ck <- NUMC[c_index]
      
      # First ck eigenvectors
      UU = eigs$vectors[, 1:ck, drop = FALSE]
      # Discretize (SNFtool internal function)
      EigenvectorsDiscrete <- SNFtool:::.discretisation(UU)[[1]]
      EigenVectors = EigenvectorsDiscrete^2
      
      # The same "temp1" manipulations as the original
      temp1 <- EigenVectors[do.call(order, lapply(seq_len(ncol(EigenVectors)), 
                                                  function(i) EigenVectors[, i])),
                            , drop = FALSE]
      temp1 <- t(apply(temp1, 1, sort, decreasing = TRUE))
      
      # The "cost" or "quality" measure:
      quality[c_index] = (1 - eigs$values[ck + 1]) / (1 - eigs$values[ck]) * 
        sum(sum(
          diag(1 / (temp1[, 1] + .Machine$double.eps)) %*%
            temp1[, 1:max(2, ck - 1), drop = FALSE]
        ))
    }
    
    # The original code picks the smallest quality as best
    quality_sort <- sort(quality, decreasing = FALSE, index.return = TRUE)
    bestQ1_idx <- quality_sort$ix[1]
    bestQ2_idx <- quality_sort$ix[2]
    
    # The top two cluster numbers from the rotation method
    K2  = NUMC[bestQ1_idx]
    K22 = NUMC[bestQ2_idx]
    
    # Record their actual rotation "cost" scores
    rotation_K2_score  = quality[bestQ1_idx]
    rotation_K22_score = quality[bestQ2_idx]
  }
  
  # -------------------------------------------
  # Return everything in a single named list
  # -------------------------------------------
  return(list(
    # Identical to original in terms of K1, K12, K2, K22
    K1  = K1,
    K12 = K12,
    K2  = K2,
    K22 = K22,
    
    # Additional numeric scores we are now exposing
    eigengap_K1_score  = eigengap_K1_score,
    eigengap_K12_score = eigengap_K12_score,
    rotation_K2_score  = rotation_K2_score,
    rotation_K22_score = rotation_K22_score
  ))
}

# ANF custom concordance by NMI function #####
concordanceNetworkNMI_ANF = function (Wall, C, type) 
{
  LW = length(Wall)
  labels = lapply(Wall, function(x) ANF::spectral_clustering(x, 
                                                       C,
                                                       type = type))
  NMIs = matrix(NA, LW, LW)
  for (i in 1:LW) {
    for (j in 1:LW) {
      NMIs[i, j] = calNMI(labels[[i]], labels[[j]])
    }
  }
  return(NMIs)
}

# ANF feature ranking by NMI modification #####
rankFeaturesByNMI_parallely_ANF <- function(data, W, ncores = detectCores() - 1, binary = FALSE,
                                            type = "rw", nn = 15) {
  stopifnot(class(data) == "list" && length(data) == 1)  # Ensure only one data type is passed
  
  NUM_OF_FEATURES <- ncol(data[[1]])
  NMI_scores <- vector(mode = "numeric", length = NUM_OF_FEATURES)
  problematic_features <- vector(mode = "list")  # To store indices of problematic features
  num_of_clusters_fused <- estimateNumberOfClustersGivenGraph(W)[[1]]
  clustering_fused <- spectral_clustering(W, num_of_clusters_fused, type = type)
  
  # Set up parallel backend to use with foreach
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  
  # Ensure cluster is stopped in case of error
  on.exit(stopCluster(cl))
  
  clusterEvalQ(cl, library(SNFtool))
  clusterEvalQ(cl, library(ANF))
  clusterEvalQ(cl, library(foreach))
  clusterEvalQ(cl, library(doParallel))
  
  # Use foreach to parallelize over features (compatible with Windows)
  data_type_scores <- foreach(feature_ind = 1:NUM_OF_FEATURES, .combine = 'c', .packages = c("SNFtool", "ANF")) %dopar% {
    tryCatch({
      if (binary) {
        # Use binary distance
        dist_matrix <- as.matrix(dist(as.matrix(data[[1]][, feature_ind]), 
                                      as.matrix(data[[1]][, feature_ind]),
                                      method = "binary"))
      } else {
        # Use default distance (assumed to be Euclidean)
        dist_matrix <- dist2(as.matrix(data[[1]][, feature_ind]), as.matrix(data[[1]][, feature_ind]))
      }
      
      affinity_matrix <- ANF::affinity_matrix(dist_matrix, alpha = 1/6, beta = 1/6, k = nn)      
      clustering_single_feature <- spectral_clustering(affinity_matrix, num_of_clusters_fused,
                                                       type = type)
      calNMI(clustering_fused, clustering_single_feature)
    }, error = function(e) {
      # Log the index of the problematic feature and return NA for its score
      problematic_features <<- append(problematic_features, feature_ind)
      NA
    })
  }
  
  # Rank the features, excluding NA values from ranking
  data_type_ranks <- rank(-data_type_scores, ties.method = "first", na.last = "keep")
  
  # Print or log problematic feature indices
  if (length(problematic_features) > 0) {
    cat("Problematic feature indices:", unlist(problematic_features), "\n")
  }
  
  return(list(NMI_scores = data_type_scores, NMI_ranks = data_type_ranks, problematic_features = problematic_features))
}
