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

# Sort results alphabetically
names(clusterings) = algorithms
clusterings = clusterings[sort(algorithms)]
algorithm_languages = algorithm_languages[names(clusterings)]

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
library(rcartocolor)

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

# Create logo mix using magick
library(magick)

# Scale the R logo to 200x200
r_img <- image_read(file.path(home, "Resources/r.png")) |> 
  image_scale("200x200!")

# Scale the Python logo to 200x200
py_img <- image_read(file.path(home, "Resources/python.png")) |> 
  image_scale("200x200!")

# Write them back to disk so we can read them as arrays
image_write(r_img, file.path(home, "Resources/r_200x200.png"))
image_write(py_img, file.path(home, "Resources/python_200x200.png"))

library(png)

# Read each 200x200 image as a numeric array:  [height, width, channels]
# Typically RGBA => a 4-channel array
r_array   <- readPNG(file.path(home, "Resources/r_200x200.png"))      # shape: 200 x 200 x 4
py_array  <- readPNG(file.path(home, "Resources/python_200x200.png"))  # shape: 200 x 200 x 4

nr <- dim(r_array)[1]  # 200
nc <- dim(r_array)[2]  # 200
# Create a blank result array, same shape
res_array <- array(0, dim = c(nr, nc, 4))

# We want:
#  - the "upper triangle" (row < col) to come from the R logo
#  - the "lower triangle" (row > col) to come from the Python logo
#  - the main diagonal (row == col) to be white
#

for(i in seq_len(nr)) {
  for(j in seq_len(nc)) {
    
    if(j > i) {
      # Above diagonal => pick from R
      res_array[i, j, ] <- r_array[i, j, ]
      
    } else if(j < i) {
      # Below diagonal => pick from Python
      res_array[i, j, ] <- py_array[i, j, ]
      
    } else {
      # On the diagonal => white background
      # RGBA for white = c(1,1,1,1)
      res_array[i, j, ] <- c(1,1,1,1)
    }
  }
}

# Write out the combined image
writePNG(res_array, file.path(home, "Resources/r_python_mashup.png"))

# Import logos
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
row_logo <- rowAnnotation(
  Software = anno_image(logo_paths, border = FALSE), gp = gpar(col = "white"),
  show_annotation_name = FALSE, width = unit(5, "mm")
)
column_logo <- HeatmapAnnotation(
  Software = anno_image(logo_paths, border = FALSE), gp = gpar(col = "white"),
  show_annotation_name = FALSE, height = unit(5, "mm")
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
  left_annotation = row_logo,
  bottom_annotation = column_logo,
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
lgd_software <- Legend(
  # Labels shown in the legend
  labels = c("R", "Python", "R & Python"),
  at = c("R", "Python", "R & Python"),
  
  # Title of the legend
  title = "",
  
  direction = "horizontal",
  nrow = 1,
  
  # Control the label/title font sizes
  labels_gp = gpar(fontsize = 6),
  
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

png(paste0(home, "/Results/Comparisons/ARI_clusterings_heatmap.png"), 
    width = 4300, height = 4300, res = 700)
draw(ARI_heatmap, 
     annotation_legend_list = list(lgd_software),
     # merge_legend = TRUE, 
     heatmap_legend_side = "right",
     annotation_legend_side = "bottom")

dev.off()

pdf(paste0(home, "/Results/Comparisons/ARI_clusterings_heatmap.pdf"), 
    width = 7, height = 7)

draw(ARI_heatmap,
     annotation_legend_list = list(lgd_software),
     # merge_legend = TRUE,
     heatmap_legend_side = "right",
     annotation_legend_side = "bottom")

dev.off()