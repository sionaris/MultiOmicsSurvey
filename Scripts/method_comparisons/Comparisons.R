# Libraries
library(openxlsx)

# Ensure reproducibility
RNGversion("4.2.2")
set.seed(123)

# Load custom helper functions
source("Scripts/automated_scripts/custom_functions.R")
source("Scripts/automated_scripts/modified_MOVICS_functions.R")

# Preamble
home = getwd()
data_source = "TCGA" # e.g. TCGA, TCGA-transNEO, transNEO-PARTNER
data_types = "RNAseq-CNV-Methylation-miRNA-SNPs" # e.g. RNAseq, RNAseq-CNV-miRNA
evaluation_source = "transNEO" # e.g. PARTNER, transNEO-PARTNER 

# Import clusterings
R_algorithms = c("ab-SNF", "ANF", "CIMLR", "COCA", "iClusterBayes", "KLIC",
                 "LRAcluster", "MDICC", "MFA", "mixKernel", #"MOFA", 
                 "NEMO", "RWR-F", "RWR-NF", "SNF", "Spectrum", "wMKL")
Python_algorithms = c("MONET", "MSNE", "PAMOGK", "MOFA")
algorithms = c(R_algorithms, Python_algorithms)
algorithm_languages = c(rep("R", length(R_algorithms)),
                        rep("Python", length(Python_algorithms)))
names(algorithm_languages) = algorithms
algorithm_languages["MOFA"] = "R & Python"
algorithm_languages["MDICC"] = "R & Python"
algorithm_languages["MixKernel"] = "R & Python"

clusterings = list()

# R methods
for (R_algorithm in R_algorithms) {
  res_dir = paste0(home, "/Results/single_algorithm/", R_algorithm, "/")
  res_dir_files = list.files(path = res_dir, pattern = ".*_clusterings\\.xlsx$", 
                             full.names = TRUE)
  if (length(res_dir_files == 1)) {
    clusterings[[R_algorithm]] = read.xlsx(res_dir_files[1])
  } else {
    clusterings[[R_algorithm]] = NA
  }
}

# Python methods
for (Python_algorithm in Python_algorithms) {
  res_dir = paste0(home, "/Results/single_algorithm/", Python_algorithm, "/")
  res_dir_files = list.files(path = res_dir, pattern = ".*_clusterings\\.xlsx$", 
                             full.names = TRUE)
  if (length(res_dir_files == 1)) {
    clusterings[[Python_algorithm]] = read.xlsx(res_dir_files[1])
  } else {
    clusterings[[Python_algorithm]] = NA
  }
}
names(clusterings) = algorithms

# Set up method categories
similarity_network_methods = c("ab-SNF", "ANF", "MDICC", "MSNE", "NEMO", "RWR-F", "RWR-NF", "SNF", "Spectrum")
multiple_kernel_learning = c("CIMLR", "KLIC", "mixKernel", "wMKL")
matrix_factorization = c("MFA", "MOFA")
graph_methods = c("MONET", "PAMOGK")
bayesian = c("iClusterBayes")
# cca_methods = c("RGCCA", "SGCCA")
# low_rank_methods = c("LRAcluster") #, moCluster, PINSPlus
misc = c("COCA", "LRAcluster")

# Primary annotation
primary_annotation_rag = c(rep("Similarity Network", length(similarity_network_methods)),
                           rep("Multiple Kernel Learning", length(multiple_kernel_learning)),
                           rep("Matrix Factorization", length(matrix_factorization)),
                           rep("Graph-based Methods", length(graph_methods)),
                           rep("Bayesian", length(bayesian)),
                           # rep("Canonical Correlation", length(cca_methods)),
                           # rep("Low-rank Projection", length(low_rank_methods)),
                           rep("Miscellaneous", length(misc)))
names(primary_annotation_rag) = c(similarity_network_methods, multiple_kernel_learning,
                                  matrix_factorization, graph_methods, bayesian,
                                  # cca_methods, low_rank_methods, 
                                  misc)

clusterings = clusterings[names(primary_annotation_rag)]
algorithm_languages = algorithm_languages[names(primary_annotation_rag)]

# # Secondary annotation
# secondary_annotation_rag = c(rep("Graph-based methods", 8),
#                              rep("Low-rank Projection", 7),
#                              rep("Miscellaneous", 7))
# names(secondary_annotation_rag) = c(similarity_network_methods, "Spectrum",
#                                     cca_methods, bayesian, matrix_factorization,
#                                     multiple_kernel_learning, graph_methods,
#                                     low_rank_methods, "COCA")

# Convert individual cluster labels from "algorithm#" to just #
generic_clusterings = clusterings
for (i in 1:length(generic_clusterings)) {
  if (!is.null(ncol(generic_clusterings[[i]]))) { # temporary error control
    generic_clusterings[[i]]$Cluster <- gsub("[^0-9]", "", generic_clusterings[[i]]$Cluster)
  }
}

rm(res_dir, res_dir_files, R_algorithm, Python_algorithm); gc()


# ARI similarity heatmap #####
library(mclust)

# Initialize a matrix to store ARI values
ari_matrix <- matrix(0, length(clusterings), 
                     length(clusterings),
                     dimnames = list(names(clusterings), 
                                     names(clusterings)))

# Calculate ARI for each pair of cluster results
for (i in 1:length(generic_clusterings)) {
  for (j in 1:length(generic_clusterings)) {
    
    # Temporary error control
    if (is.null(ncol(generic_clusterings[[i]])) || is.null(ncol(generic_clusterings[[j]]))) {
      ari_matrix[i, j] = NA
    } else {
      ari_matrix[i, j] <- ari_matrix[j, i] <- adjustedRandIndex(
        generic_clusterings[[i]]$Cluster, 
        generic_clusterings[[j]]$Cluster
      )
    }
  }
}

diag(ari_matrix) = 1

# Print ARI matrix
print(ari_matrix)
ari_matrix <- as.matrix(ari_matrix)
class(ari_matrix) <- "numeric"

# Draw heatmap
library(ComplexHeatmap)
library(circlize)
# Define colors for method categories
library(rcartocolor)
category_colors <- c(
  "Similarity Network" = carto_pal("Bold", n = 12)[1],
  "Multiple Kernel Learning" = carto_pal("Bold", n = 12)[2],
  "Matrix Factorization" = carto_pal("Antique", n = 12)[5],
  "Graph-based Methods" = carto_pal("Bold", n = 12)[4],
  "Bayesian" = carto_pal("Bold", n = 12)[11],
  # "Canonical Correlation" = carto_pal("Bold", n = 12)[9],
  # "Low-rank Projection" = carto_pal("Bold", n = 12)[10],
  "Miscellaneous" = carto_pal("Bold", n = 12)[12]
)

primary_annotation_rag <- factor(primary_annotation_rag, levels = names(category_colors))
# secondary_annotation_rag <- factor(secondary_annotation_rag, levels = names(category_colors))

# Convert language vector to factor
algorithm_languages <- factor(algorithm_languages, levels = c("R", "Python", "R & Python"))

# Read in the logos
# r_logo_img <- readPNG("Resources/r.png") # <a href="https://www.flaticon.com/free-icons/r" title="r icons">R icons created by Becris - Flaticon</a>
# python_logo_img <- readPNG("Resources/python.png") # <a href="https://www.flaticon.com/free-icons/python" title="python icons">Python icons created by Freepik - Flaticon</a>

# Color-blind friendly palette
color_palette <- carto_pal(n = 100, name = "RedOr")

# Create a masked version of the matrix to show only the left (lower) triangle
ari_matrix_masked <- ari_matrix

# To be removed when all clusterings are gathered
ari_matrix_masked[which(is.na(ari_matrix_masked))] = .Machine$double.eps
ari_matrix_masked[upper.tri(ari_matrix_masked)] <- NA

# Create a logo mashup
# source("Scripts/create_logo_mashup.R")

# === Logo annotations (for software) ===
logo_paths <- vapply(
  rownames(ari_matrix),
  FUN.VALUE = character(1),
  FUN = function(alg) {
    if (algorithm_languages[alg] == "R") {
      paste0(home, "/Resources/r.png")
    } else if (algorithm_languages[alg] == "Python") {
      paste0(home, "/Resources/python.png")
    } else if (algorithm_languages[alg] == "R & Python") {
      paste0(home, "/Resources/r_python_mashup.png")
    } else {
      NA_character_
    }
  }
)

# Primary track – show legend (with title "Method category")
ROWannotation <- rowAnnotation(
  Software = anno_image(logo_paths, border = FALSE, width = unit(6, "mm")),
  Category = as.character(primary_annotation_rag),
  col = list(Category = category_colors),
  gp = gpar(col = "white"),
  show_annotation_name = FALSE,
  simple_anno_size = unit(1.5, "mm"),
  # width = unit(7, "mm"),
  show_legend = FALSE
)

COLannotation <- HeatmapAnnotation(
  Category = as.character(primary_annotation_rag),
  Software = anno_image(logo_paths, border = FALSE, height = unit(6, "mm")),
  col = list(Category = category_colors),
  gp = gpar(col = "white"),
  show_annotation_name = FALSE,
  simple_anno_size = unit(1.5, "mm"),
  # height = unit(7, "mm"),
  show_legend = FALSE
)

ARI_heatmap <- Heatmap(
  ari_matrix_masked, 
  name = "ARI Index", 
  column_title = "Adjusted Rand Index (ARI) between clusterings", 
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  col = color_palette, 
  na_col = "white",
  row_names_side = "left",
  cluster_rows = FALSE, 
  cluster_columns = FALSE, 
  show_row_names = TRUE, 
  show_column_names = TRUE,
  left_annotation = ROWannotation,
  bottom_annotation = COLannotation,
  row_names_gp = grid::gpar(fontsize = 7, fontface = "bold"), 
  column_names_gp = grid::gpar(fontsize = 7, fontface = "bold"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    val <- ari_matrix_masked[i, j]
    if (!is.na(val)) {
      if (val == .Machine$double.eps) {
        # For cells with .Machine$double.eps: override the background to white and print "NA"
        grid::grid.rect(x = x, y = y, width = width, height = height,
                        gp = grid::gpar(fill = "white", col = NA))
        grid::grid.text("NA", x, y, gp = grid::gpar(col = "black", fontsize = 6))
      } else {
        grid::grid.text(sprintf("%.2f", val), x, y, 
                        gp = grid::gpar(col = "black", fontsize = 6))
      }
    }
  },
  heatmap_legend_param = list(
    title = "ARI Index",
    title_gp = grid::gpar(fontsize = 7, fontface = "bold"), 
    labels_gp = grid::gpar(fontsize = 7),
    legend_height = unit(3, "cm"),
    grid_width = unit(0.25, "cm"),
    title_position = "leftcenter-rot"
  )
)

# Software legend
library(png)
r_array <- readPNG(file.path(home, "Resources/r_200x200.png"))
py_array <- readPNG(file.path(home, "Resources/python_200x200.png"))
res_array <- readPNG(file.path(home, "Resources/r_python_mashup.png"))

lgd_software <- Legend(
  # Labels shown in the legend
  labels = c("R", "Python", "R & Python"),
  at = c("R", "Python", "R & Python"),
  
  # Title of the legend
  title = "Software",
  title_position = "leftcenter",
  
  direction = "horizontal",
  nrow = 1,
  
  # Control the label/title font sizes
  labels_gp = gpar(fontsize = 7),
  title_gp = gpar(fontsize = 9, fontface = "bold"),
  
  # These settings remove any drawn borders around the symbol boxes
  legend_gp = gpar(col = NA),
  background = "white",
  
  # Size of each symbol box in the legend
  grid_width  = unit(4, "mm"),
  grid_height = unit(4, "mm"),
  
  graphics = list(
    # 1. R logo
    function(x, y, w, h) {
      grid.raster(r_array, x = x, y = y, width = w, height = h)
    },
    # 2. Python logo
    function(x, y, w, h) {
      grid.raster(py_array, x = x, y = y, width = w, height = h)
    },
    # 3. Mashup logo
    function(x, y, w, h) {
      grid.raster(res_array, x = x, y = y, width = w, height = h)
    }
  )
)

# Create legend for method categories
lgd_methods <- Legend(
  labels = names(category_colors),
  legend_gp = gpar(fill = category_colors, col = NA),
  title = "Category",
  labels_gp = gpar(fontsize = 7),
  title_gp = gpar(fontsize = 9, fontface = "bold"),
  ncol = 4,
  title_position = 'leftcenter'
)

png(paste0(home, "/Results/Comparisons/ARI_clusterings_heatmap.png"), 
    width = 4300, height = 4600, res = 700)
draw(ARI_heatmap, 
     annotation_legend_list = packLegend(lgd_software, lgd_methods),
     heatmap_legend_side = "right",
     annotation_legend_side = "bottom",
     align_annotation_legend = "heatmap_center")
dev.off()

pdf(paste0(home, "/Results/Comparisons/ARI_clusterings_heatmap.pdf"), 
    width = 7, height = 7.5)
draw(ARI_heatmap,
     annotation_legend_list = packLegend(lgd_software, lgd_methods),
     heatmap_legend_side = "right",
     annotation_legend_side = "bottom",
     align_annotation_legend = "heatmap_center")
dev.off()

# ARI with ground truth(s) #####
# Pick a TCGA .rds file produced through the download and preprocessing script
tcga = readRDS("Resources/TCGA/RNA_full.rds")

# Convert tcga object columns to the sample identifiers we have in this work
new_colnames = unlist(lapply(strsplit(colnames(tcga), split = "-"), function(x) {
  paste(x[1:4], collapse = "-")
}))
colnames(tcga) = new_colnames; rm(new_colnames); gc()

# Filter for samples of this work
# MOVICS
MOVICS_clusters = read.xlsx("Results/MOVICS_baseline/MOVICS_TCGA_RNAseq-CNV-Methylation-miRNA-SNPs_eval_on_transNEO_clusterings.xlsx")

# PARADIGM
PARADIGM_clusters = tcga[, clusterings[[1]]$Sample.ID]@colData@listData[["paper_PARADIGM Clusters"]]
PARADIGM_clusters = as.data.frame(list(Sample.ID = clusterings[[1]]$Sample.ID,
                                       Cluster = PARADIGM_clusters))
PARADIGM_clusters$Cluster = gsub("C", "", PARADIGM_clusters$Cluster)
# PARADIGM_clusters= na.omit(PARADIGM_clusters)
PARADIGM_clusters$Cluster = as.numeric(PARADIGM_clusters$Cluster)

# PanGyn
PanGyn_clusters = tcga[, clusterings[[1]]$Sample.ID]@colData@listData[["paper_Pan-Gyn Clusters"]]
PanGyn_clusters = as.data.frame(list(Sample.ID = clusterings[[1]]$Sample.ID,
                                       Cluster = PanGyn_clusters))
PanGyn_clusters$Cluster = gsub("C", "", PanGyn_clusters$Cluster)
# PanGyn_clusters = na.omit(PanGyn_clusters)
PanGyn_clusters$Cluster = as.numeric(PanGyn_clusters$Cluster)

# Loop of ARI calculations
ARI_df = data.frame(matrix(NA, ncol = 5, nrow = length(clusterings)))
colnames(ARI_df) = c("algorithm", "Category", "ARI_to_MOVICS", "ARI_to_PARADIGM", "ARI_to_PanGyn")
ARI_df$algorithm = names(clusterings)
ARI_df$Category = primary_annotation_rag[ARI_df$algorithm]
ARI_df = ARI_df %>% dplyr::arrange(Category)
rownames(ARI_df) = ARI_df$algorithm

for (algorithm in ARI_df$algorithm) {
  if (is.list(clusterings[[algorithm]])) {
    cluster_df = as.data.frame(clusterings[[algorithm]])
    ARI_df[algorithm, "ARI_to_MOVICS"] = calculate_ari_index(cluster_df1 = MOVICS_clusters,
                                                             cluster_df2 = cluster_df,
                                                             sample_col = "Sample.ID",
                                                             clust_col = "Cluster",
                                                             suffixes = c("_MOVICS_CS",
                                                                          paste0("_", algorithm)))
    ARI_df[algorithm, "ARI_to_PARADIGM"] = calculate_ari_index(cluster_df1 = PARADIGM_clusters,
                                                               cluster_df2 = cluster_df,
                                                               sample_col = "Sample.ID",
                                                               clust_col = "Cluster",
                                                               suffixes = c("_PARADIGM_C",
                                                                            paste0("_", algorithm)))
    ARI_df[algorithm, "ARI_to_PanGyn"] = calculate_ari_index(cluster_df1 = PanGyn_clusters,
                                                             cluster_df2 = cluster_df,
                                                             sample_col = "Sample.ID",
                                                             clust_col = "Cluster",
                                                             suffixes = c("_Pan-Gyn_C",
                                                                          paste0("_", algorithm)))
  } else {
    ARI_df[algorithm, c("ARI_to_MOVICS", "ARI_to_PARADIGM", "ARI_to_PanGyn")] = rep(NA, 3)
  }
}

# Plot histograms with three facets
library(ggplot2)
library(patchwork)
library(ggnewscale)
library(cowplot)
library(ggpubr)

# Create a new column for annotation x-position
ARI_df$annot_x <- -0.035
ARI_df <- ARI_df[order(ARI_df$Category, ARI_df$algorithm), ]
ARI_df$algorithm <- factor(ARI_df$algorithm, levels = unique(ARI_df$algorithm))

p1 <- ggplot(ARI_df, aes(y = reorder(algorithm, Category))) +
  # Annotation tile: Category legend suppressed
  geom_tile(aes(x = annot_x, fill = Category), 
            color = "grey50", linewidth = 0.01,
            width = 0.015, height = 0.8) +
  scale_fill_manual(values = category_colors, guide = FALSE) +
  new_scale_fill() +
  # Bar layer: ARI values with its own legend
  geom_bar(aes(x = ARI_to_MOVICS, fill = ARI_to_MOVICS), 
           stat = "identity", color = NA) +
  scale_fill_carto_c(palette = "Magenta", guide = guide_colorbar(title = "ARI")) +
  coord_cartesian(xlim = c(-0.04, 1)) +
  labs(
    title = "ARI bar chart: MOVICS",
    x = "Adjusted Rand Index (ARI)",
    y = "Algorithm"
  ) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title.y = element_blank(),
    axis.title.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 7),
    legend.title = element_text(face = "bold")
  )

# p2: PARADIGM plot
p2 <- ggplot(ARI_df, aes(y = reorder(algorithm, Category))) +
  geom_tile(aes(x = annot_x, fill = Category), 
            color = "grey50", linewidth = 0.01,
            width = 0.015, height = 0.8) +
  scale_fill_manual(values = category_colors, guide = FALSE) +
  new_scale_fill() +
  geom_bar(aes(x = ARI_to_PARADIGM, fill = ARI_to_PARADIGM), 
           stat = "identity", color = NA) +
  scale_fill_carto_c(palette = "Teal", guide = guide_colorbar(title = "ARI")) +
  coord_cartesian(xlim = c(-0.04, 1)) +
  labs(
    title = "ARI bar chart: PARADIGM",
    x = "Adjusted Rand Index (ARI)",
    y = "Algorithm"
  ) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title.y = element_blank(),
    axis.title.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 7),
    legend.title = element_text(face = "bold")
  )

# p3: Pan-Gyn plot
p3 <- ggplot(ARI_df, aes(y = reorder(algorithm, Category))) +
  geom_tile(aes(x = annot_x, fill = Category), 
            color = "grey50", linewidth = 0.01,
            width = 0.015, height = 0.8) +
  scale_fill_manual(values = category_colors, guide = FALSE) +
  new_scale_fill() +
  geom_bar(aes(x = ARI_to_PanGyn, fill = ARI_to_PanGyn), 
           stat = "identity", color = NA) +
  scale_fill_carto_c(palette = "Peach", guide = guide_colorbar(title = "ARI")) +
  coord_cartesian(xlim = c(-0.04, 1)) +
  labs(
    title = "ARI bar chart: Pan-Gyn",
    x = "Adjusted Rand Index (ARI)",
    y = "Algorithm"
  ) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(face = "bold", size = 10),
    axis.text.x = element_text(size = 7),
    legend.title = element_text(face = "bold")
  )

# Combine the three plots into one patchwork layout
combined_plots <- p1 + p2 + p3 + plot_layout(ncol = 3)

# Create a dummy plot solely for the Category legend
legend_plot <- ggplot(ARI_df, aes(x = annot_x, y = reorder(algorithm, Category), fill = Category)) +
  geom_tile(width = 0.015, height = 0.8) +
  scale_fill_manual(values = category_colors, guide = guide_legend(title = "Category")) +
  theme_void() + 
  theme(legend.position = "bottom")

# Extract the legend using cowplot
legend_category <- get_legend(legend_plot)

# Combine the patchwork with the extracted legend and add vertical spacer for padding
final_plot <- combined_plots / plot_spacer() / as_ggplot(legend_category) +
  plot_layout(heights = c(10, 0.5, 1))  # Adjust the middle value for extra padding

# Save the final plot
ggsave(
  filename = paste0(home, "/Results/Comparisons/Ground_truth_ARI_barchart.png"),
  plot = final_plot,
  dpi = 700,
  width = 5 * 1920,
  height = 2.5 * 1920,
  units = "px",
  device = "png"
)
ggsave(
  filename = paste0(home, "/Results/Comparisons/Ground_truth_ARI_barchart.pdf"),
  plot = final_plot,
  dpi = 700,
  width = 5 * 1920,
  height = 2.5 * 1920,
  units = "px",
  device = "pdf"
)

# Kernel PCA for the ARI matrix #####
# Calculate ARI for each pair of cluster results
nonas_generic_clusterings = generic_clusterings[which(unlist(lapply(generic_clusterings, is.list)))]
nonas_ari_matrix = matrix(0, ncol = length(nonas_generic_clusterings), 
                          nrow = length(nonas_generic_clusterings))
for (i in 1:length(nonas_generic_clusterings)) {
  for (j in 1:length(nonas_generic_clusterings)) {
    nonas_ari_matrix[i, j] <- nonas_ari_matrix[j, i] <- adjustedRandIndex(
      nonas_generic_clusterings[[i]]$Cluster, 
      nonas_generic_clusterings[[j]]$Cluster
    )
  }
}
dimnames(nonas_ari_matrix) = list(names(nonas_generic_clusterings),
                                  names(nonas_generic_clusterings))
diag(nonas_ari_matrix) = 1

# PCA plot
library(dplyr)
# MOdify the pre-defined pca_from_sim_matrix() function
pca_from_ari_matrix = function (sim_matrix = NULL, clust_res = NULL,
                                cluster_colors = NULL, output_path = NULL,
                                pointsize = NULL,
                                geom_label_size = NULL) {
  
  # Kernel PCA of the clusters (modified from kernel_pca() from Spectrum)
  km = sim_matrix
  m = nrow(km)
  kc = t(t(km - colSums(km)/m) - rowSums(km)/m) + sum(km)/m^2
  res = eigen(kc/m, symmetric = TRUE)
  features = m
  ret = suppressWarnings(t(t(res$vectors[, 1:features]) / sqrt(res$values[1:features])))
  scores = data.frame(ret)
  rownames(scores) = colnames(km)
  colnames(scores)[1:2] = c("PC1", "PC2")
  scores$algorithm = rownames(scores)
  
  # Filter the scores dataset for annotation by joining with clust_res
  plot_df = scores %>% inner_join(clust_res, by = "algorithm")
  rownames(plot_df) = plot_df$algorithm
  
  if(is.null(pointsize)) { pointsize = 2 }
  if(is.null(geom_label_size)) { geom_label_size = 2 }
  
  # Actual plot: points are now colored by Category.
  kernelPCA = ggplot(plot_df, aes(x = PC1, y = PC2)) +
    # Points colored by Category with a fixed size.
    geom_point(aes(color = Category), size = pointsize) +
    # Print the name (algorithm) for each point.
    # geom_text(aes(label = algorithm), hjust = -0.1, vjust = 0.5, size = geom_label_size) +
    ggrepel::geom_text_repel(aes(label = algorithm),
                    size = pointsize,
                    max.overlaps = Inf, min.segment.length = 0,
                    segment.size = 0.1) +
    # Use the provided category color scale.
    scale_color_manual(name = "Category", values = cluster_colors) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.border = element_rect(linewidth = 0.2),
          plot.title = element_text(size = 5, face = "bold"),
          axis.title.x = element_text(size = 4, face = "bold"),
          axis.title.y = element_text(size = 4, face = "bold"),
          axis.text = element_text(size = 4),
          axis.ticks = element_line(linewidth = 0.15),
          legend.background = element_rect(fill = "white", linetype = "solid"),
          legend.key.size = unit(0.5, "lines"),
          legend.margin = margin(0, 0, 0, 0, unit = "mm"),
          legend.spacing.y = unit(0.1, "cm"), 
          legend.title = element_text(size = 4, face = "bold"), 
          legend.text = element_text(size = 3)) +
    labs(title = "Kernel PCA from the algorithm ARI matrix",
         x = "PC1",
         y = "PC2")
  return(kernelPCA)
}

# Run
kernelPCA = pca_from_ari_matrix(
  sim_matrix = nonas_ari_matrix,
  clust_res = ARI_df[rownames(nonas_ari_matrix), ],
  cluster_colors = category_colors,
  pointsize = 0.75,
  geom_label_size = 1)
  
print(kernelPCA)

ggsave(filename = "ARI_kernelPCA.png",
       path = paste0(home, "/Results/Comparisons"), 
       width = 1920, height = 1080, device = 'png', units = "px",
       dpi = 700)
dev.off()

print(kernelPCA)

ggsave(filename = "ARI_kernelPCA.pdf",
       path = paste0(home, "/Results/Comparisons"), 
       width = 1920, height = 1080, device = 'pdf', units = "px",
       dpi = 700)
dev.off()

# Runtimes #####

# Save environment
save.image(paste0(home, "/Results/Comparisons/Comparisons_", data_source, "_",
                  data_types, "_eval_on_", evaluation_source,
                  "_env.RData"))
