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
R_algorithms = c("ab-SNF", "ANF", "CIMLR", "COCA", "iClusterBayes", "IntNMF", "KLIC",
                 "LRAcluster", "MDICC", "MFA", "mixKernel", "MOFA", "NEMO", "PIntMF",
                 "RGCCA", "RWR-F", "SGCCA", "SNF", "Spectrum")
Python_algorithms = c("MONET", "MSNE", "PAMOGK" # "MOFA-GPU"
)
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
  res_dir = paste0(home, "/Python/", Python_algorithm, "/")
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
similarity_network_methods = c("ab-SNF", "ANF", "MDICC", "MSNE", "NEMO", "RWR-F", "SNF")
multiple_kernel_learning = c("CIMLR","KLIC", "mixKernel") #, wMKL
matrix_factorization = c("IntNMF", "MFA", "MOFA", "PIntMF")
graph_methods = c("MONET", "PAMOGK")
bayesian = c("iClusterBayes")
cca_methods = c("RGCCA", "SGCCA")
low_rank_methods = c("LRAcluster") #, moCluster, PINSPlus
misc = c("COCA", "Spectrum")

# Primary annotation
primary_annotation_rag = c(rep("Similarity Network", length(similarity_network_methods)),
                           rep("Multiple Kernel Learning", length(multiple_kernel_learning)),
                           rep("Matrix Factorization", length(matrix_factorization)),
                           rep("Graph-based Methods", length(graph_methods)),
                           rep("Bayesian", length(bayesian)),
                           rep("Canonical Correlation", length(cca_methods)),
                           rep("Low-rank Projection", length(low_rank_methods)),
                           rep("Miscellaneous", length(misc)))
names(primary_annotation_rag) = c(similarity_network_methods, multiple_kernel_learning,
                                  matrix_factorization, graph_methods, bayesian,
                                  cca_methods, low_rank_methods, misc)

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
  "Graph-based Methods" = carto_pal("Bold", n = 12)[11],
  "Bayesian" = carto_pal("Bold", n = 12)[4],
  "Canonical Correlation" = carto_pal("Bold", n = 12)[9],
  "Low-rank Projection" = carto_pal("Bold", n = 12)[10],
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