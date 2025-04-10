# Import data from gitignored "Resources/TCGA/" folder #####

# These datasets are pre-standardized
data_object = readRDS("Resources/TCGA/norm_data_object.rds")
clinical_data = openxlsx::read.xlsx("Resources/TCGA/clinical_data.xlsx")

# Setup environment variables for markdown #####

# Ensure reproducibility
RNGversion("4.2.2")
set.seed(123)

# Load custom helper functions
source("Scripts/automated_scripts/custom_functions.R")
source("Scripts/automated_scripts/modified_MOVICS_functions.R")

# Preamble
home = getwd()
algorithm = "MOVICS"
alg_feature_pref = "rows" # Where does the algorithm expect the features to be
citation = fetch_citation(algorithm = algorithm)
data_source = "TCGA" # e.g. TCGA, TCGA-transNEO, transNEO-PARTNER
data_types = "RNAseq-CNV-Methylation-miRNA-SNPs" # e.g. RNAseq, RNAseq-CNV-miRNA
evaluation_source = "transNEO" # e.g. PARTNER, transNEO-PARTNER 
title = paste0("Results from ", algorithm)
subtitle = paste0("<b>Train</b>: ", data_source, " ", data_types, 
                  " | <b>Evaluation</b>: ", evaluation_source)
in_a_nutshell = fetch_in_a_nutshell(algorithm = algorithm)
optk_boolean = "TRUE" # either TRUE or FALSE. Answers whether the algorithm suggests an optimal k
optk_text = ifelse(optk_boolean == TRUE,
                   "<u>suggests</u> an estimate of the optimal number of multi-omic clusters $k$",
                   "<u>does not suggest</u> an optimal number of multi-omic clusters $k$")

# Detailed description of the algorithm
description = paste(readLines(paste0("Resources/algorithm_descriptions/", algorithm,
                                     "_description.Rmd")),
                    collapse = "\n") # File path to .Rmd file within Resources/algorithm_descriptions

# Preprocessing flags and code #####
library(stringr)
library(dplyr)

# Extract the names of the modalities that will be used
modalities = unlist(strsplit(data_types, "-"))

# Replace with TRUE wherever standardization is required or FALSE otherwise
standardization_booleans = rep(FALSE, length(modalities))

# The same logic follows for features in rows
features_in_rows = rep(FALSE, length(modalities))

# Replace with TRUE wherever features are in rows. Each row will then be standardized
features_in_rows = rep(TRUE, length(modalities))

# If features are in columns, each column will be standardized

# The feature column vector is TRUE when features are not rownames, but a column
feature_column = rep(FALSE, length(modalities))

# Replace with either a numeric value or a column name where applicable
# Leave NA if features are in columns
# Use feature_column EVEN IF data are pre-standardized, as long as the data frame
# contains a feature column

# Give names to the vectors
names(standardization_booleans) = names(features_in_rows) = names(feature_column) = modalities

# Perform standardization of features if required
input = data_object
for (i in 1:length(modalities)) {
  if (standardization_booleans[modalities[i]] == TRUE) {
    if (features_in_rows[i] == TRUE && feature_column[i] != FALSE) {
      
      # Check if numerical or character string was used for the feature column indicator
      feature_column_indicator_type = ifelse(is.numeric(feature_column[i]),
                                             "numeric", "character")
      
      # Extract features
      features = data_object[[i]][, feature_column[i]]
      
      # Remove feature column based on indicator type
      if (feature_column_indicator_type == "numeric") {
        z_data = data_object[[i]][, -feature_column[i]]
      } else {
        z_data = data_object[[i]] %>% dplyr::select(-feature_column[i])
      }
      
      # Convert to matrix and standardize
      cols = colnames(z_data) # Extract colnames
      z_data = as.matrix(z_data)
      z_data = t(apply(z_data, 1, scale))
      rownames(z_data) = features # Set/ensure rownames
      colnames(z_data) = cols # Set/ensure colnames
      
      # Save the transformed matrix in the input object
      input[[modalities[i]]] = z_data
      rm(feature_column_indicator_type, features, z_data, cols)
    } else if (features_in_rows[i] == TRUE && feature_column[i] == FALSE) {
      
      # We assume the features are the rownames
      features = rownames(data_object[[i]])
      cols = colnames(data_object[[i]])
      z_data = as.matrix(data_object[[i]])
      z_data = t(apply(z_data, 1, scale))
      rownames(z_data) = features # Set/ensure rownames
      colnames(z_data) = cols # Set/ensure colnames
      
      # Save the transformed matrix in the input object
      input[[modalities[i]]] = z_data
      rm(features, z_data, cols)
    } else if (features_in_rows[i] == FALSE) {
      
      # We assume that all columns are numeric and are to be standardized
      rows = rownames(data_object[[i]])
      features = colnames(data_object[[i]])
      z_data = as.matrix(data_object[[i]])
      z_data = apply(z_data, 1, scale)
      rownames(z_data) = rows # Set/ensure rownames
      colnames(z_data) = features # Set/ensure colnames
      
      # Save the transformed matrix in the input object
      input[[modalities[i]]] = z_data
      rm(features, z_data, rows)
    } 
  } else if (standardization_booleans[modalities[i]] == FALSE && feature_column[i] != FALSE) {
    # Only used for pre-standardized data which have the feature names in one column
    
    # Check if numerical or character string was used for the feature column indicator
    feature_column_indicator_type = ifelse(is.numeric(feature_column[i]),
                                           "numeric", "character")
    
    # Extract the rownames from the designated column
    features = data_object[[i]][, feature_column[i]]
    # Remove feature column based on indicator type
    if (feature_column_indicator_type == "numeric") {
      z_data = data_object[[i]][, -feature_column[i]]
    } else {
      z_data = data_object[[i]] %>% dplyr::select(-feature_column[i])
    }
    
    rownames(z_data) = features # Set/ensure colnames
    
    # Save the matrix in the input object
    input[[modalities[i]]] = z_data
    rm(features, z_data)
  }
}

rm(i); gc()

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
rm(rogue_indices); gc()

if (alg_feature_pref == "rows") {
  # Sample names are in the columns
  overlap = Reduce(intersect, lapply(input, colnames))
  
  # Filter inputs
  input = lapply(input, function(x) {
    x = x[, overlap]
  })
} else {
  # Sample names are in the rows
  overlap = Reduce(intersect, lapply(input, rownames))
  
  # Filter inputs
  input = lapply(input, function(x) {
    x = x[overlap, ]
  })
}

# Export input object as a resource
input$Methylation = na.omit(input$Methylation) # Remove the missing values
saveRDS(input, "Resources/TCGA/mm_input.rds")
rm(data_object); gc()

# Calculate feature distance matrices once and use for all other algorithms
# mainly for the comprehensive heatmap
feature_hclusts = list()
for (i in 1:length(input)) {
  if (names(input)[i] == "SNPs") {
    feature_hclusts[[i]] = fastcluster::hclust(dist(as.matrix(input$SNPs),
                                          as.matrix(input$SNPs),
                                          method = "binary"), method = "ward.D")
  } else {
  feature_hclusts[[i]] = fastcluster::hclust(as.dist(Rfast::Dist(input[[i]], 
                                                                 method = "euclidean")),
                                             method = "ward.D")
  }
}
names(feature_hclusts) = names(input)
orders = lapply(feature_hclusts, function(x) return(x$order))
names(orders) = names(feature_hclusts)

# Export orders for future use
saveRDS(orders, "Resources/TCGA/mm_feature_orders.rds")

# Run algorithm #####
library(MOVICS)

check_missing_values_with_indices <- function(input) {
  # Apply the function to each element in the list
  missing_info <- lapply(input, function(mat) {
    if (!is.matrix(mat)) {
      stop("All elements of input should be matrices.")
    }
    
    # Find the row indices where there are missing values
    missing_indices <- which(rowSums(is.na(mat)) > 0)
    
    # Count the number of NA values in the matrix
    na_count <- sum(is.na(mat))
    
    # Return a list containing the count of NA values and the row indices
    return(list(na_count = na_count, missing_indices = missing_indices))
  })
  
  # Combine the results into a named list
  names(missing_info) <- names(input)
  return(missing_info)
}

# Example usage
missing_values_summary <- check_missing_values_with_indices(input) # all clear
rm(missing_values_summary); gc()

# identify optimal clustering number (may take a while)
try_min_k = 2
try_max_k = 10
optk = getClustNum(data = input,
                   is.binary = c(T,F,F,F,F),
                   try.N.clust = try_min_k:try_max_k,
                   center = FALSE,
                   scale = FALSE,
                   fig.path = "Results/MOVICS_baseline",
                   fig.name = paste0("optimal_k_plot_", data_source,
                                     "_", data_types))
gc()

# Perform multi-omic clustering with 9 available methods using default parameters
# iClusterBayes will be run on the cluster
moic.res.list = getMOIC(data = input,
                        methodslist = list("SNF", "CIMLR", "PINSPlus", "NEMO", 
                                           "COCA", # "MoCluster",
                                           "LRAcluster", "ConsensusClustering", 
                                           "IntNMF"),
                                           # , "iClusterBayes"),
                        N.clust = optk$N.clust,
                        type = c("binomial", "gaussian", "gaussian", 
                                 "gaussian", "gaussian"))

# Save results to local file
save(moic.res.list, file = paste0(home, "/Results/MOVICS_baseline/", 
                                  algorithm, "_", data_source, "_",
                                  data_types, "_eval_on_", evaluation_source,
                                  "_moic.res.list.rda"))

# iClusterBayes with lower burnin and draw parameter values
n_burnin_iCB = 1800
n_draw_iCB = 1200
iClusterBayes.res = getMOIC(data        = input,
                            N.clust     = optk$N.clust,
                            methodslist = "iClusterBayes",
                            type        = c("binomial",
                                            "gaussian",
                                            "gaussian",
                                            "gaussian",
                                            "gaussian"),
                            n.burnin    = n_burnin_iCB,
                            n.draw      = n_draw_iCB,
                            prior.gamma = c(0.5, 0.5, 0.5, 0.5, 0.5),
                            sdev        = 0.05,
                            thin        = 3)

saveRDS(iClusterBayes.res, paste0(home, "/Results/MOVICS_baseline/", 
                                  algorithm, "_", data_source, "_",
                                  data_types, "_eval_on_", evaluation_source, 
                                  "_", "iCB_mild.rds"))

# Append moic.res.list
moic.res.list = append(moic.res.list, 
                       list("iClusterBayes" = iClusterBayes.res))
rm(iClusterBayes.res); gc()

# moCluster (was throwing out an error when run with the comprehensive getMOIC)
moCluster.res = getMOIC(data        = input,
                            N.clust     = optk$N.clust,
                            methodslist = "MoCluster",
                            type        = c("binomial",
                                            "gaussian",
                                            "gaussian",
                                            "gaussian",
                                            "gaussian"),
                            center      = FALSE,
                            scale       = FALSE)
moic.res.list = append(moic.res.list, 
                       list("MoCluster" = moCluster.res))
rm(moCluster.res); gc()

# Save final moic.res.list object
save(moic.res.list, file = paste0(home, "/Results/MOVICS_baseline/", 
                                  algorithm, "_", data_source, "_",
                                  data_types, "_eval_on_", evaluation_source,
                                  "_moic.res.list.rda")); gc()

# Inspect output for similarities/differences across clusterings #####
# Initialize a matrix to store Jaccard indices
jaccard_matrix <- matrix(0, length(moic.res.list), 
                         length(moic.res.list),
                         dimnames = list(names(moic.res.list), 
                                         names(moic.res.list)))

# Calculate Jaccard index for each pair of cluster results (using loaded custom function)
for (i in 1:length(moic.res.list)) {
  for (j in i:length(moic.res.list)) {
    jaccard_matrix[i, j] <-jaccard_matrix[j, i] <-  MOVICS_jaccard_index(moic.res.list[[i]]$clust.res,
                                                                  moic.res.list[[j]]$clust.res)
  }
}

# Print the Jaccard matrix
print(jaccard_matrix)

# Draw heatmap
jaccard_matrix <- as.matrix(jaccard_matrix)
class(jaccard_matrix) <- "numeric"

# Define the color palette using viridis
color_palette <- viridis::viridis(100)

# Define a function to decide label colors based on background color
label_color <- circlize::colorRamp2(c(0, 1), c("black", "white"))

# Draw the heatmap
library(ComplexHeatmap)
jaccard_heatmap = Heatmap(jaccard_matrix, 
                          name = "Jaccard Index", 
                          column_title = "Jaccard Similarity between clusterings", 
                          column_title_gp = gpar(fontsize = 8, fontface = "bold"),
                          col = color_palette, 
                          cluster_rows = FALSE, 
                          cluster_columns = FALSE, 
                          show_row_names = TRUE, 
                          show_column_names = TRUE,
                          row_names_gp = grid::gpar(fontsize = 6, fontface = "bold"), 
                          column_names_gp = grid::gpar(fontsize = 6, fontface = "bold"),
                          cell_fun = function(j, i, x, y, width, height, fill) {
                            grid::grid.text(sprintf("%.2f", jaccard_matrix[i, j]), x, y, 
                                            gp = grid::gpar(col = label_color(fill), fontsize = 6))
                          },
                          heatmap_legend_param = list(
                            title = "Jaccard Index",
                            title_gp = grid::gpar(fontsize = 6, fontface = "bold"), 
                            labels_gp = grid::gpar(fontsize = 6),
                            legend_height = unit(1.5, "cm"),
                            grid_width = unit(0.25, "cm"),
                            title_position = "leftcenter-rot"
                          ))

png(paste0(home, "/Results/MOVICS_baseline/", 
           algorithm, "_", data_source, "_",
           data_types, "_eval_on_", evaluation_source, 
           "_Jaccard_clusterings_heatmap.png"), 
    width = 4300, height = 4300, res = 700)
draw(jaccard_heatmap)
dev.off()

# Calculate ARI and NMI
library(mclust)
library(clue)

# Initialize a matrix to store ARI values
ari_matrix <- matrix(0, length(moic.res.list), 
                     length(moic.res.list),
                     dimnames = list(names(moic.res.list), 
                                     names(moic.res.list)))

# Calculate ARI for each pair of cluster results
for (i in 1:length(moic.res.list)) {
  for (j in i:length(moic.res.list)) {
    ari_matrix[i, j] <- ari_matrix[j, i] <- adjustedRandIndex(
      moic.res.list[[i]]$clust.res$clust, 
      moic.res.list[[j]]$clust.res$clust
    )
  }
}

# Print ARI matrix
print(ari_matrix)

# Draw heatmap
ari_matrix <- as.matrix(ari_matrix)
class(ari_matrix) <- "numeric"

ARI_heatmap = Heatmap(ari_matrix, 
                      name = "ARI Index", 
                      column_title = "Adjusted Rand Index (ARI) between clusterings", 
                      column_title_gp = gpar(fontsize = 8, fontface = "bold"),
                      col = color_palette, 
                      cluster_rows = FALSE, 
                      cluster_columns = FALSE, 
                      show_row_names = TRUE, 
                      show_column_names = TRUE,
                      row_names_gp = grid::gpar(fontsize = 6, fontface = "bold"), 
                      column_names_gp = grid::gpar(fontsize = 6, fontface = "bold"),
                      cell_fun = function(j, i, x, y, width, height, fill) {
                        grid::grid.text(sprintf("%.2f", ari_matrix[i, j]), x, y, 
                                        gp = grid::gpar(col = label_color(fill), fontsize = 6))
                      },
                      heatmap_legend_param = list(
                        title = "ARI Index",
                        title_gp = grid::gpar(fontsize = 6, fontface = "bold"), 
                        labels_gp = grid::gpar(fontsize = 6),
                        legend_height = unit(1.5, "cm"),
                        grid_width = unit(0.25, "cm"),
                        title_position = "leftcenter-rot"
                      ))

png(paste0(home, "/Results/MOVICS_baseline/", 
           algorithm, "_", data_source, "_",
           data_types, "_eval_on_", evaluation_source, 
           "_ARI_clusterings_heatmap.png"), 
    width = 4300, height = 4300, res = 700)
draw(ARI_heatmap)
dev.off()

# Initialize a matrix to store NMI values
nmi_matrix <- matrix(0, length(moic.res.list), 
                     length(moic.res.list),
                     dimnames = list(names(moic.res.list), 
                                     names(moic.res.list)))

# Calculate NMI for each pair of cluster results
for (i in 1:length(moic.res.list)) {
  for (j in i:length(moic.res.list)) {
    nmi_matrix[i, j] <- nmi_matrix[j, i] <- clue::cl_agreement(
      as.cl_partition(moic.res.list[[i]]$clust.res$clust), 
      as.cl_partition(moic.res.list[[j]]$clust.res$clust), 
      method = "NMI"
    )
  }
}

# Print NMI matrix
print(nmi_matrix)

# Draw heatmap
nmi_matrix <- as.matrix(nmi_matrix)
class(nmi_matrix) <- "numeric"

NMI_heatmap = Heatmap(nmi_matrix, 
                      name = "NMI Index", 
                      column_title = "Normalized Mutual Information (NMI) between clusterings", 
                      column_title_gp = gpar(fontsize = 8, fontface = "bold"),
                      col = color_palette, 
                      cluster_rows = FALSE, 
                      cluster_columns = FALSE, 
                      show_row_names = TRUE, 
                      show_column_names = TRUE,
                      row_names_gp = grid::gpar(fontsize = 6, fontface = "bold"), 
                      column_names_gp = grid::gpar(fontsize = 6, fontface = "bold"),
                      cell_fun = function(j, i, x, y, width, height, fill) {
                        grid::grid.text(sprintf("%.2f", nmi_matrix[i, j]), x, y, 
                                        gp = grid::gpar(col = label_color(fill), fontsize = 6))
                      },
                      heatmap_legend_param = list(
                        title = "NMI Index",
                        title_gp = grid::gpar(fontsize = 6, fontface = "bold"), 
                        labels_gp = grid::gpar(fontsize = 6),
                        legend_height = unit(1.5, "cm"),
                        grid_width = unit(0.25, "cm"),
                        title_position = "leftcenter-rot"
                      ))

png(paste0(home, "/Results/MOVICS_baseline/", 
           algorithm, "_", data_source, "_",
           data_types, "_eval_on_", evaluation_source, 
           "_NMI_clusterings_heatmap.png"), 
    width = 4300, height = 4300, res = 700)
draw(NMI_heatmap)
dev.off()

# Consensus #####
# get consensus results from all algorithms
consensus = getConsensusMOIC(moic.res.list = moic.res.list,
                             fig.path = "Results/MOVICS_baseline",
                             fig.name = "MOVICS_consensus_heatmap",
                             distance = "euclidean",
                             linkage = "average",
                             showID = FALSE)

# Show silhouette metrics across clusters
sil = compute_silhouette(cluster_df = consensus$clust.res %>% dplyr::rename(Cluster = clust),
                         similarity_matrix = consensus$similarity.matrix,
                         normalize_matrix = TRUE)

getSilhouette_ggplot(sil      = sil,
                     fig.path = paste0(home, "/Results/MOVICS_baseline"),
                     fig.name = "Silhouette",
                     height   = 5.5,
                     width    = 5.5,
                     axis_label_size = 12,
                     axis_label_font = "bold",
                     text_size = 1.5,
                     title_size = 16,
                     algorithm = "CS",
                     save_plot = TRUE)
dev.off()

# Downstream comparisons #####
# Create data frame for comparisons
cc_res = consensus$clust.res
colnames(cc_res) = c("Sample.ID", "Consensus Subtype")
cc_res$`Consensus Subtype` = factor(cc_res$`Consensus Subtype`,
                                    levels = c(1, 2),
                                    labels = c("CS1", "CS2"))
subset_of_interest = clinical_data[, c("Sample.ID", "vital_status", "days_to_birth", "days_to_last_known_alive", "days_to_death",
                                       "days_to_last_followup", "race_list", "history_of_neoadjuvant_treatment",
                                       "age_at_initial_pathologic_diagnosis", "ethnicity", "primary_lymph_node_presentation_assessment",
                                       "histological_type", "menopause_status",  "breast_carcinoma_progesterone_receptor_status", 
                                       "breast_carcinoma_estrogen_receptor_status", "lab_proc_her2_neu_immunohistochemistry_receptor_status",
                                       "number_of_lymphnodes_positive_by_ihc", "number_of_lymphnodes_positive_by_he", "er_level_cell_percentage_category",
                                       "progesterone_receptor_level_cell_percent_category","distant_metastasis_present_ind2",
                                       "stage_event_pathologic_stage")] %>%
  #dplyr::rename(Sample.ID = bcr_patient_barcode) %>%
  #dplyr::mutate(Sample.ID = paste0(Sample.ID, "-01")) %>%
  group_by(Sample.ID) %>%
  arrange(Sample.ID, rowSums(is.na(across(-Sample.ID)))) %>%  # Arrange by Sample.ID and NA count
  slice(1) %>%  # Keep the first occurrence in case of ties
  ungroup() 

var2comp = cc_res %>%
  inner_join(subset_of_interest,
             by = "Sample.ID")
rownames(var2comp) = var2comp$Sample.ID
var2comp = var2comp %>% 
  dplyr::select(-Sample.ID)

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
clin_comp = compClinvar_single_algorithm(algorithm_name = "CS",
                                         moic.res = consensus,
                                         var2comp = var2comp_nonas,
                                         strata = "Consensus Subtype",
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
                                         res.path = paste0(home, "/Results/MOVICS_baseline/"),
                                         output_pdf = TRUE,
                                         pdf_level_col_width = c("7em", "10em"),
                                         pdf_count_col_width = "10em",
                                         pdf_pval_col_width = "3em",
                                         pdf_test_col_width = "8em",
                                         pdf_tab_font_size = 9)

clin_ordinal_comp = compClinvar_ordinal_single_algorithm(algorithm_name = "CS",
                                                         moic.res = consensus,
                                                         var2comp = var2comp_nonas %>%
                                                           dplyr::select(number_of_lymphnodes_positive_by_ihc,
                                                                         number_of_lymphnodes_positive_by_he,
                                                                         `Consensus Subtype`),
                                                         strata = "Consensus Subtype",
                                                         ordinalVars = c("number_of_lymphnodes_positive_by_ihc",
                                                                         "number_of_lymphnodes_positive_by_he"),
                                                         includeNA = FALSE,
                                                         tab.name = "Summary of ordinal clinical variables",
                                                         res.path = paste0(home, "/Results/MOVICS_baseline/"),
                                                         output_pdf = TRUE,
                                                         pdf_template_loc = paste0(home, "/Scripts/automated_scripts/clincomp_template.Rmd"),
                                                         pdf_level_col_width = c("7em", "10em"),
                                                         pdf_count_col_width = "10em",
                                                         pdf_pval_col_width = "3em",
                                                         pdf_test_col_width = "8em",
                                                         pdf_tab_font_size = 9)

# Color annotation #####
# set color for each omics data
# if no color list specified all subheatmaps will be unified to green and red color pattern
mRNA.col   <- c("#00ff00", "white", "#ff0000")
CNV.col <- c("#6699CC", "white", "#FF3C38")
mut.col   <- c("#EFE5AF", "#780A43")
methylation.col    <- c("#57087C", "white", "#FF3C38")
miRNA.col <- c("#D3ACEF", "white", "#FA076B")
col.list   <- list(mut.col, mRNA.col, CNV.col, miRNA.col, methylation.col)

# Create annCol object that will be used for plot annotation and colors
library(forcats)
library(tidyr)

# Select and rename columns
annCol = var2comp %>%
  dplyr::select(`Vital status` = vital_status,
                Ethnicity = ethnicity,
                Race = race_list,
                `Lymph node status` = primary_lymph_node_presentation_assessment,
                Histology = histological_type,
                `Menopausal status` = menopause_status,
                `PR status` = breast_carcinoma_progesterone_receptor_status,
                `ER status` = breast_carcinoma_estrogen_receptor_status,
                `HER2 status` = lab_proc_her2_neu_immunohistochemistry_receptor_status,
                Metastasis = distant_metastasis_present_ind2,
                Stage = stage_event_pathologic_stage)
rownames(annCol) = rownames(var2comp)

# Replace empty strings with NA
annCol[annCol == ""] <- NA

# Replace NA values with "Unknown"
annCol <- annCol %>%
  mutate(across(everything(), ~ifelse(is.na(.), "Unknown", as.character(.)))) %>%
  mutate(across(everything(), as.factor))

# Relabel factors where necessary
annCol$Ethnicity = factor(str_to_sentence(annCol$Ethnicity))
annCol$Race = factor(str_to_sentence(annCol$Race))
annCol$`Lymph node status` = factor(str_to_sentence(annCol$`Lymph node status`))
levels(annCol$Histology)[levels(annCol$Histology) == "Mixed Histology (please specify)"] <- "Mixed"
levels(annCol$Histology)[levels(annCol$Histology) == "Other, specify"] <- "Other"
levels(annCol$`Menopausal status`)[levels(annCol$`Menopausal status`) == "Indeterminate (neither Pre or Postmenopausal)"] = "Indeterminate"
levels(annCol$`Menopausal status`)[levels(annCol$`Menopausal status`) == "Peri (6-12 months since last menstrual period)"] = "Perimenopausal"
levels(annCol$`Menopausal status`)[levels(annCol$`Menopausal status`) == "Post (prior bilateral ovariectomy OR >12 mo since LMP with no prior hysterectomy)"] = "Post-menopausal"
levels(annCol$`Menopausal status`)[levels(annCol$`Menopausal status`) == "Pre (<6 months since LMP AND no prior bilateral ovariectomy AND not on estrogen replacement)"] = "Pre-menopausal"
annCol$Metastasis = factor(str_to_sentence(annCol$Metastasis))

# Relabel stage
stage1 = c("Stage I", "Stage IA", "Stage IB")
stage2 = c("Stage II", "Stage IIA", "Stage IIB")
stage3 = c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC")
stage4 = c("Stage IV")

levels(annCol$Stage)[levels(annCol$Stage) == "Stage X"] = "Unknown"
levels(annCol$Stage)[levels(annCol$Stage) %in% stage1] = "Stage I"
levels(annCol$Stage)[levels(annCol$Stage) %in% stage2] = "Stage II"
levels(annCol$Stage)[levels(annCol$Stage) %in% stage3] = "Stage III"
levels(annCol$Stage)[levels(annCol$Stage) %in% stage4] = "Stage IV"

# Generate corresponding colors for sample annotation
histol_colors = rcartocolor::carto_pal(n = 9, "Safe")
names(histol_colors) = levels(annCol$Histology)
histol_colors[["Unknown"]] = "grey40"

annColors = list(
  Stage = c(`Stage I` = "#00C9FF", `Stage II` = "#099CF5", 
            `Stage III` = "#097BF5", `Stage IV` = "#0B5684", `Unknown` = "grey40"),
  `Lymph node status` = c(No = "grey75", Yes = "#4A0558", `Unknown` = "grey40"),
  `ER status` = c(Negative = "#C11D9C", Positive = "#0F1682", `Unknown` = "grey40"),
  `PR status` = c(Indeterminate = "aliceblue", Positive = "dodgerblue4", 
                  Negative = rcartocolor::carto_pal(n = 7, "ArmyRose")[5], `Unknown` = "grey40"),
  `HER2 status` = c(Negative = "#0B9EF8", Positive = "#560DA7",
                    Indeterminate = "mistyrose1", Equivocal = "hotpink4", `Unknown` = "grey40"),
  `Vital status` = c(Alive = "lightpink1", Dead = "black", `Unknown` = "grey40"),
  Ethnicity = c(`Hispanic or latino` = rcartocolor::carto_pal(n = 12, "Vivid")[1],
                `Not hispanic or latino` = rcartocolor::carto_pal(n = 12, "Vivid")[6], 
                `Unknown` = "grey40"),
  Race = c(`American indian or alaska native` = rcartocolor::carto_pal(n = 12, "Bold")[5],
           Asian = rcartocolor::carto_pal(n = 12, "Bold")[3],
           `Black or african american` = rcartocolor::carto_pal(n = 12, "Prism")[12],
           White = "beige", `Unknown` = "grey40"),
  Metastasis = c(Yes = "deeppink4", No = "cadetblue2", `Unknown` = "grey40"),
  Histology = histol_colors,
  `Menopausal status` = c(Indeterminate = "mistyrose2",
                          `Pre-menopausal` = rcartocolor::carto_pal(n = 7, "SunsetDark")[2],
                          `Perimenopausal` = rcartocolor::carto_pal(n = 7, "SunsetDark")[5],
                          `Post-menopausal` = rcartocolor::carto_pal(n = 7, "SunsetDark")[7],
                          `Unknown` = "grey40")
)

plotdata <- lapply(lapply(input, as.matrix), 
                   function(mat) mat[rowSums(mat != 0) > 0, ])

# Use halfwidth for beter coloring in heatmap
plotdata = getStdiz(
  data = plotdata,
  halfwidth = c(NA, 3, 3, 3, 3), # No halfwidth for SNPs
  centerFlag = c(F, F, F, F, F),
  scaleFlag = c(F, F, F, F, F)
)

# Export coloring settings for other algorithms
scheme = list(col.list = col.list,
              annColors = annColors,
              annCol = annCol,
              var2comp = var2comp,
              clust.colors = c("#2EC4B6", "#E71D36"))
saveRDS(scheme, "Resources/scheme.rds")

# comprehensive heatmap (may take a while)
getMoHeatmap(data          = plotdata,
             row.title     = names(plotdata),
             is.binary     = c(T,F,F,F,F), 
             legend.name   = c("SNPs",
                               "Standardized RNAseq norm. counts",
                               "Standardized CNV",
                               "Standardized miRNA norm. counts",
                               "Standardized Methylation M-values"
             ),
             clust.res     = consensus$clust.res, # consensusMOIC results
             clust.dend    = NULL, # show no dendrogram for samples
             show.rownames = c(F,F,F,F,F), # specify for each omics data
             show.colnames = FALSE, # show no sample names
             show.row.dend = c(F,F,F,F,F), # show no dendrogram for features
             annRow        = NULL, # no selected features
             color         = col.list,
             annCol        = annCol, # annotation for samples
             annColors     = annColors, # annotation color
             width         = 15, # width of each subheatmap
             height        = 10, # height of each subheatmap
             fig.path      = paste0(home, "/Results/MOVICS_baseline"),
             fig.name      = "default_Comprehensive_heatmap")
dev.off()
gc()

# # Comparison of survival curves
# surv.info = clinical_data %>%
#   dplyr::select(Patient.ID, samID = Sample.ID, vital_status, days_to_death, days_to_last_followup) %>%
#   left_join(consensus$clust.res, by = "samID") %>%
#   dplyr::select(-Patient.ID) %>%
#   dplyr::rename(fustat = vital_status, Subtype = clust)
# surv.info$fustat[which(surv.info$fustat == "")] = NA
# 
# surv.info$futime = ifelse(surv.info$fustat == "Alive",
#                           surv.info$days_to_last_followup,
#                           surv.info$days_to_death)
# surv.info = surv.info %>% dplyr::select(-days_to_death, -days_to_last_followup)
# 
# # Remove all NAs
# surv.info = na.omit(surv.info)
# 
# surv.info$fustat = factor(surv.info$fustat, labels = c(0, 1),
#                           levels = c("Alive", "Dead"))
# surv.info = distinct(surv.info, samID, .keep_all = TRUE)
# rownames(surv.info) = surv.info$samID
# 
# library(survival)
# surv.brca <- compSurv_ext(moic.res = consensus, surv.info = surv.info,
#                       convt.time = "m", # convert day unit to month
#                       surv.median.line = "h", # draw horizontal line at median survival
#                       xyrs.est = c(5,10), # estimate 5 and 10-year survival
#                       fig.name = "Kaplan Meier curve of Consensus Subtypes",
#                       fig.path = paste0(home, "/Results/MOVICS_baseline")) # BH adjustment by default
# print(surv.brca)

# Plotting using the MOVICS package, which uses pheatmap underneath
oncoprint <- compMut(moic.res  = consensus,
                     mut.matrix   = input$SNPs, # binary somatic mutation matrix
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
                     fig.path     = paste0(home, "/Results/MOVICS_baseline"),
                     res.path     = paste0(home, "/Results/MOVICS_baseline"))

# Drug sensitivity comparison #####
# TPM matrix prior to log2 transformation is recommended. Here we used the normalized input
drug_sensitivity <- compDrugsen(moic.res    = consensus,
                                norm.expr   = plotdata$RNAseq,
                                drugs       = c("Cisplatin", "Paclitaxel", "Lapatinib",
                                                "Doxorubicin", "5-Fluorouracil",
                                                "Sorafenib"), # a vector of names of drug in GDSC
                                tissueType  = "breast", # choose specific tissue type to construct ridge regression model
                                test.method = "nonparametric", # statistical testing method
                                prefix      = "Violin_plot_of_IC50",
                                seed = 123,
                                fig.path = paste0(home, "/Results/MOVICS_baseline"))

# Agreement with other subtypes #####
subtype_agreement <- compAgree2(moic.res  = consensus,
                                subt2comp = annCol[, c("Stage", "ER status", "PR status",
                                                       "HER2 status", "Metastasis")],
                                doPlot    = TRUE,
                                box.width = 0.2,
                                fig.name  = "Classification_agreement",
                                fig.path  = paste0(home, "/Results/MOVICS_baseline"),
                                width     = 12)

# DGEA #####
dgea = runDEA(dea.method = "limma", # we use count data as input
              expr = plotdata$RNAseq,
              moic.res = consensus,
              prefix = "dgea_",
              sort.p = TRUE,
              overwt = TRUE,
              verbose = TRUE,
              res.path = paste0(home, "/Results/MOVICS_baseline"))

# # Identify unique subtype biomarkers
# # 1. Up-regulated markers
dgea.marker.up <- runMarker_mod_4.4(moic.res = consensus,
                            dea.method    = "limma", # name of DEA method
                            prefix        = "dgea_", # MUST be the same of argument in runDEA()
                            dat.path      = paste0(home, "/Results/MOVICS_baseline"), # path of DEA files
                            res.path      = paste0(home, "/Results/MOVICS_baseline"), # path to save marker files
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
                            fig.path = paste0(home, "/Results/MOVICS_baseline"),
                            width = 14,
                            height = 12,
                            fontsize_row = 3,
                            name = "normalized RNA-seq")
dev.off()

# # 2. Down-regulated markers
dgea.marker.down <- runMarker_mod_4.4(moic.res = consensus,
                            dea.method    = "limma", # name of DEA method
                            prefix        = "dgea_", # MUST be the same of argument in runDEA()
                            dat.path      = paste0(home, "/Results/MOVICS_baseline"), # path of DEA files
                            res.path      = paste0(home, "/Results/MOVICS_baseline"), # path to save marker files
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
                            fig.path = paste0(home, "/Results/MOVICS_baseline"),
                            width = 14,
                            height = 12,
                            fontsize_row = 3,
                            name = "normalized RNA-seq")
dev.off()

# DMEA ###
dmea = runDEA(dea.method = "limma", # we use normalized data as input
                  expr = plotdata$Methylation,
                  moic.res = consensus,
                  prefix = "dmea_",
                  sort.p = TRUE,
                  overwt = TRUE,
                  verbose = TRUE,
                  res.path = paste0(home, "/Results/MOVICS_baseline"))

# # Identify unique subtype biomarkers
# # 1. Up-regulated markers
methyl.marker.up <- runMarker_mod_4.4(moic.res = consensus,
                                               dea.method    = "limma", # name of DEA method
                                               prefix        = "dmea_", # MUST be the same of argument in runDEA()
                                               dat.path      = paste0(home, "/Results/MOVICS_baseline"), # path of DEA files
                                               res.path      = paste0(home, "/Results/MOVICS_baseline"), # path to save marker files
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
                                               fig.path = paste0(home, "/Results/MOVICS_baseline"),
                                               width = 14,
                                               height = 12,
                                               fontsize_row = 3, # 3 default
                                               name = "normalized Methylation")
dev.off()

# # 2. Down-regulated markers
methyl.marker.down <- runMarker_mod_4.4(moic.res = consensus,
                                                 dea.method    = "limma", # name of DEA method
                                                 prefix        = "dmea_", # MUST be the same of argument in runDEA()
                                                 dat.path      = paste0(home, "/Results/MOVICS_baseline"), # path of DEA files
                                                 res.path      = paste0(home, "/Results/MOVICS_baseline"), # path to save marker files
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
                                                 fig.path = paste0(home, "/Results/MOVICS_baseline"),
                                                 width = 14,
                                                 height = 12,
                                                 fontsize_row = 3, # 3 default
                                                 name = "normalized Methylation")
dev.off()

# DmiREA ###
dmiRea = runDEA(dea.method = "limma", # we use normalized data as input
                    expr = plotdata$miRNA,
                    moic.res = consensus,
                    prefix = "dmiRea_",
                    sort.p = TRUE,
                    overwt = TRUE,
                    verbose = TRUE,
                    res.path = paste0(home, "/Results/MOVICS_baseline"))

# # Identify unique subtype biomarkers
# # 1. Up-regulated markers
miRNA.marker.up <- runMarker_mod_4.4(moic.res = consensus,
                                              dea.method    = "limma", # name of DEA method
                                              prefix        = "dmiRea_", # MUST be the same of argument in runDEA()
                                              dat.path      = paste0(home, "/Results/MOVICS_baseline"), # path of DEA files
                                              res.path      = paste0(home, "/Results/MOVICS_baseline"), # path to save marker files
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
                                              fig.path = paste0(home, "/Results/MOVICS_baseline"),
                                              width = 14,
                                              height = 12,
                                              fontsize_row = 3, # 3 default
                                              name = "normalized miRNA")
dev.off()

# # 2. Down-regulated markers
miRNA.marker.down <- runMarker_mod_4.4(moic.res = consensus,
                                                dea.method    = "limma", # name of DEA method
                                                prefix        = "dmiRea_", # MUST be the same of argument in runDEA()
                                                dat.path      = paste0(home, "/Results/MOVICS_baseline"), # path of DEA files
                                                res.path      = paste0(home, "/Results/MOVICS_baseline"), # path to save marker files
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
                                                fig.path = paste0(home, "/Results/MOVICS_baseline"),
                                                width = 14,
                                                height = 12,
                                                fontsize_row = 3, # 3 default
                                                name = "normalized miRNA")
dev.off()

# GSEA #####
# Load MSigDb file
MSIGDB.FILE <- paste0(home, "/Resources/Pathways/GO-BP_c5.go.bp.v2024.1.Hs.symbols.gmt")

# GSEA up-regulated
RNGversion("4.2.2")
set.seed(123)
gsea.up <- runGSEA_mod_4.4(moic.res     = consensus,
                       dea.method   = "limma", # name of DEA method
                       prefix       = "dgea_", # MUST be the same of argument in runDEA()
                       dat.path      = paste0(home, "/Results/MOVICS_baseline"), # path of DEA files
                       res.path      = paste0(home, "/Results/MOVICS_baseline"), # path to save marker files
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
                       fig.path = paste0(home, "/Results/MOVICS_baseline"),
                       width = 14, height = 12)

# GSEA down-regulated
RNGversion("4.2.2")
set.seed(123)
gsea.down <- runGSEA_mod_4.4(moic.res     = consensus,
                           dea.method   = "limma", # name of DEA method
                           prefix       = "dgea_", # MUST be the same of argument in runDEA()
                           dat.path      = paste0(home, "/Results/MOVICS_baseline"), # path of DEA files
                           res.path      = paste0(home, "/Results/MOVICS_baseline"), # path to save marker files
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
                           fig.path = paste0(home, "/Results/MOVICS_baseline"),
                           width = 14, height = 12)

# Hierarchical clustering of pathways for representative pathways
library(pathfindR)
library(fastcluster)

# Get unique pathways for each subtype
GSEAfiles_up <- sort(dir(paste0(home, "/Results/MOVICS_baseline"), 
                 pattern = "unique_upexpr_pathway.txt$"))
GSEAfiles_down <- sort(dir(paste0(home, "/Results/MOVICS_baseline"), 
                    pattern = "unique_downexpr_pathway.txt$"))

unique_upexpr_pathways = list()
for (i in 1:length(gsea.up$gsea.list)) {
  unique_upexpr_pathways[[i]] = data.table::fread(paste0(paste0(home, "/Results/MOVICS_baseline"), 
                                                         "/", GSEAfiles_up[i]),
                                                  header = TRUE, sep = "\t")
}

unique_downexpr_pathways = list()
for (i in 1:length(gsea.down$gsea.list)) {
  unique_downexpr_pathways[[i]] = data.table::fread(paste0(paste0(home, "/Results/MOVICS_baseline"), 
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
                                              mo.method = "consensusMOIC",
                                              dat.path = paste0(home, "/Results/MOVICS_baseline"), 
                                              dgea_padj_cutoff = 0.05,
                                              logfc_cutoff = 0,
                                              pathway_padj_cutoff = 0.05)

hclust_input_down = prepare_gsea_output_for_hclust(gsea_output = gsea.down_unique, 
                                                 dgea_output_name_style = "dgea_", 
                                                 dea.method = "limma", 
                                                 mo.method = "consensusMOIC",
                                                 dat.path = paste0(home, "/Results/MOVICS_baseline"), 
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

timestamp() # ~2.5 mins
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
saveWorkbook(wb, file = paste0(home, "/Results/MOVICS_baseline/MOVICS_representative_pathways.xlsx"),
             overwrite = TRUE); rm(wb)

# Plot pathway heatmaps
hclust_pathway_plots_up = plot_pathway_heatmaps(gsea.lists = hclust_output[grepl("up", names(hclust_output))], 
                                                norm.expr = plotdata$RNAseq, 
                                                present_clusters = c("CS1", "CS2"),
                                                representative = TRUE, moic.res = consensus,
                                                subtype_prefix = "CS", n.path = 20, msigdb.path = MSIGDB.FILE,
                                                norm.method = "mean", dirct = "up",
                                                fig.name = "upregulated_pathway_heatmap",
                                                name = "GSVA scores",
                                                fig.path = paste0(home, "/Results/MOVICS_baseline"), 
                                                width = 15, height = 10, gsva.method = "gsva")

hclust_pathway_plots_down = plot_pathway_heatmaps(gsea.lists = hclust_output[grepl("down", names(hclust_output))], 
                                                norm.expr = plotdata$RNAseq, 
                                                present_clusters = c("CS1", "CS2"),
                                                representative = TRUE, moic.res = consensus,
                                                subtype_prefix = "CS", n.path = 20, msigdb.path = MSIGDB.FILE,
                                                norm.method = "mean", dirct = "down",
                                                fig.name = "downregulated_pathway_heatmap",
                                                name = "GSVA scores",
                                                fig.path = paste0(home, "/Results/MOVICS_baseline"), 
                                                width = 15, height = 10, gsva.method = "gsva")

# Gene set variation analysis #####
# locate ABSOLUTE path of gene set file
GSET.FILE <- paste0(home, "/Resources/Pathways/gene_sets_of_interest.gmt")

RNGversion("4.2.2")
set.seed(123)
gsva.res = runGSVA_mod_4.4(moic.res      = consensus,
                   norm.expr     = plotdata$RNAseq,
                   gset.gmt.path = GSET.FILE, # ABSOLUTE path of gene set file
                   gsva.method   = "gsva", # method to calculate single sample enrichment score
                   annCol        = annCol,
                   annColors     = annColors,
                   fig.path      = paste0(home, "/Results/MOVICS_baseline"),
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

# Fraction Genome Altered ###
library(ggplot2)
fga_df = readRDS("Resources/TCGA/fga_df.rds"); gc()

fga.MOVICS <- compFGA_optimized(moic.res     = consensus,
                       segment      = fga_df,
                       iscopynumber = TRUE, 
                       test.method  = "nonparametric", # statistical testing method (Wilcoxon with asymptotic approximation. Consider Kruskall Wallis?)
                       fig.path     = paste0(home, "/Results/MOVICS_baseline"),
                       fig.name     = paste0("FGA_barplot_", algorithm),
                       prefix = "CS",
                       width = 16,
                       ga_column = "ga", # genome altered column
                       clust.col = scheme$clust.colors,
                       title = "FGA plot: simple criteria")

fga.MOVICS.COSMIC <- compFGA_optimized(moic.res     = consensus,
                              segment      = fga_df,
                              iscopynumber = TRUE, 
                              test.method  = "nonparametric", # statistical testing method (Wilcoxon with asymptotic approximation. Consider Kruskall Wallis?)
                              fig.path     = paste0(home, "/Results/MOVICS_baseline"),
                              fig.name     = paste0("COSMIC_criteria_FGA_barplot_", algorithm),
                              prefix = "CS",
                              width = 16,
                              ga_column = "COSMIC_ga", # genome altered column
                              clust.col = scheme$clust.colors,
                              title = "FGA plot: COSMIC criteria")

rm(fga_df); gc()

# Run Nearest Template Prediction in transNEO cohort #####

# Load transNEO data
transNEO_mm_inputs = readRDS("Resources/transNEO/transNEO_multimodal_inputs.rds")
transcr = readRDS("Resources/transNEO/log2.norm.counts.plus1_transNEO.rds")

# get as many templates as possible
dgea.marker.up_1000 <- runMarker_single_algorithm_no_export(algorithm_name = algorithm,
                                                            moic.res = consensus,
                                                            n.marker = 1000,
                                                            dea.method    = "limma", # name of DEA method
                                                            prefix        = "dgea_", # MUST be the same of argument in runDEA()
                                                            dat.path      = paste0(home, "/Results/MOVICS_baseline"), # path of DEA files
                                                            p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                                                            p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                                                            norm.expr = plotdata$RNAseq,
                                                            dirct         = "up" # direction of dysregulation in expression
)

# 2. Down-regulated markers
dgea.marker.down_1000 <- runMarker_single_algorithm_no_export(algorithm_name = algorithm,
                                                              moic.res = consensus,
                                                              n.marker = 1000,
                                                              dea.method    = "limma", # name of DEA method
                                                              prefix        = "dgea_", # MUST be the same of argument in runDEA()
                                                              dat.path      = paste0(home, "/Results/MOVICS_baseline"), # path of DEA files
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
  fig.path = paste0(home, "/Results/MOVICS_baseline"),
  fig.name = "ntp_expr_up_heatmap_transNEO")
timestamp() # 2.5 min

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
  fig.path = paste0(home, "/Results/MOVICS_baseline"),
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
  inner_join(expr_conc %>% dplyr::select(Donor.ID = samID, `Consensus Subtype` = clust_up),
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

transNEO_clincomp = compClinvar_single_algorithm(algorithm_name = "CS",
                                                 moic.res = transNEO_ntp_expr_up,
                                                 var2comp = transNEO_var2comp_nonas %>%
                                                   mutate(`Consensus Subtype` = paste0("CS", `Consensus Subtype`)),
                                                 strata = "Consensus Subtype",
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
                                                 res.path = paste0(home, "/Results/MOVICS_baseline"),
                                                 output_pdf = TRUE,
                                                 pdf_level_col_width = c("7em", "10em"),
                                                 pdf_count_col_width = "10em",
                                                 pdf_pval_col_width = "3em",
                                                 pdf_test_col_width = "8em",
                                                 pdf_tab_font_size = 9)

# transNEO_ntp_expr_up_ord = transNEO_ntp_expr_up
# transNEO_ntp_expr_up_ord$clust.res$clust = gsub(algorithm, "", transNEO_ntp_expr_up_ord$clust.res$clust)
transNEO_ordinal_clincomp = compClinvar_ordinal_single_algorithm(algorithm_name = "CS",
                                                                 moic.res = transNEO_ntp_expr_up,
                                                                 var2comp = transNEO_var2comp_nonas %>%
                                                                   dplyr::select(Grade.pre.NAT, 
                                                                                 Chemo.cycles, 
                                                                                 aHER2.cycles,
                                                                                 `Consensus Subtype`),
                                                                 strata = "Consensus Subtype",
                                                                 ordinalVars = c("Grade.pre.NAT",
                                                                                 "Chemo.cycles",
                                                                                 "aHER2.cycles"),
                                                                 includeNA = FALSE,
                                                                 tab.name = "transNEO Summary of ordinal clinical variables",
                                                                 res.path = paste0(home, "/Results/MOVICS_baseline"),
                                                                 output_pdf = TRUE,
                                                                 pdf_template_loc = paste0(home, "/Scripts/automated_scripts/clincomp_template.Rmd"),
                                                                 pdf_level_col_width = c("7em", "10em"),
                                                                 pdf_count_col_width = "10em",
                                                                 pdf_pval_col_width = "3em",
                                                                 pdf_test_col_width = "8em",
                                                                 pdf_tab_font_size = 9)

# Run PAM
RNGversion("4.2.2")
set.seed(123)
transNEO_pam = runPAM(train.expr = plotdata$RNAseq,
                      moic.res   = consensus,
                      test.expr  = as.matrix(transcr))

# Check consistency across methods

# Get predictions for TCGA (discovery cohort)
RNGversion("4.2.2")
set.seed(123)
TCGA.ntp.pred = runNTP(expr = plotdata$RNAseq[, consensus$clust.res$samID],
                       templates = dgea.marker.up_1000$templates, distance = "cosine",
                       doPlot = F, nPerm = 10000)

TCGA.pam.pred = runPAM(train.expr = plotdata$RNAseq[, consensus$clust.res$samID],
                       moic.res = consensus,
                       test.expr = plotdata$RNAseq[, consensus$clust.res$samID])

# consensus TCGA vs NTP TCGA
runKappa(subt1 = consensus$clust.res$clust,
         subt2 = as.numeric(TCGA.ntp.pred$clust.res$clust),
         subt1.lab = "Consensus",
         subt2.lab = "NTP TCGA",
         height = 8,
         width = 8,
         fig.path = paste0(home, "/Results/MOVICS_baseline"),
         fig.name = "kappa_consensus_vs_NTP_TCGA")

# consensus TCGA vs PAM TCGA
runKappa(subt1 = consensus$clust.res$clust,
         subt2 = as.numeric(TCGA.pam.pred$clust.res$clust),
         subt1.lab = "Consensus",
         subt2.lab = "PAM TCGA",
         height = 8,
         width = 8,
         fig.path = paste0(home, "/Results/MOVICS_baseline"),
         fig.name = "kappa_consensus_vs_PAM_TCGA")

# NTP transNEO vs PAM transNEO
runKappa(subt1 = as.numeric(transNEO_ntp_expr_up$clust.res$clust),
         subt2 = as.numeric(transNEO_pam$clust.res$clust),
         subt1.lab = "transNEO NTP",
         subt2.lab = "transNEO PAM",
         height = 8,
         width = 8,
         fig.path = paste0(home, "/Results/MOVICS_baseline"),
         fig.name = "kappa_NTP_vs_PAM_transNEO")

# Export consensus clustering object
clust = as.data.frame(consensus$clust.res)
colnames(clust) = c("Sample.ID", "Cluster")
openxlsx::write.xlsx(clust, paste0(home, "/Results/MOVICS_baseline/", 
                                   algorithm, "_", data_source, "_",
                                   data_types, "_eval_on_", evaluation_source,
                                   "_clusterings.xlsx"))

# Export session info as .txt
writeLines(capture.output(sessionInfo()), paste0("sessionInfo/",
                                                 algorithm, "_", data_source, "_",
                                                 data_types, "_eval_on_", evaluation_source,
                                                 "_sessionInfo.txt"))

# Render the R Markdown document with the parameters
hyperparameters = list(min_k = try_min_k,
                       max_k = try_max_k,
                       n_burnin = n_burnin_iCB,
                       n_draw = n_draw_iCB)
params = list(algorithm = algorithm, data_source = data_source, data_types = data_types,
              citation = citation, home = home, optk = optk$N.clust,
              evaluation_source = evaluation_source, title = title, subtitle = subtitle,
              description = description, in_a_nutshell = in_a_nutshell, optk_text = optk_text,
              hyperparameters = hyperparameters, transNEO_var2comp = transNEO_var2comp,
              sessionInfo = sessionInfo())

rmarkdown::render(paste0(getwd(), "/Results/MOVICS_baseline/MOVICS_baseline_report.Rmd"), 
                  params = params, 
                  output_file = paste0(home, "/Results/MOVICS_baseline/MOVICS_baseline_report.html"))

# Save environment
save.image(paste0(home, "/Results/MOVICS_baseline/", 
                  algorithm, "_", data_source, "_",
                  data_types, "_eval_on_", evaluation_source,
                  "_env.RData"))
