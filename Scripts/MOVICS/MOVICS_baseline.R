# Import data from gitignored "Resources/BRCA complete/" folder #####

# These datasets are pre-standardized
brca_cnv = read.csv("Resources/BRCA complete/BRCA_CNV.csv")
brca_met = read.csv("Resources/BRCA complete/BRCA_Methy.csv")
brca_exp = read.csv("Resources/BRCA complete/BRCA_mRNA.csv")
brca_mirna = read.csv("Resources/BRCA complete/BRCA_miRNA.csv")
data_object = list(CNV = brca_cnv, Methylation = brca_met,
                   RNAseq = brca_exp, miRNA = brca_mirna)
rm(brca_cnv, brca_exp, brca_met, brca_mirna); gc()

data_object <- lapply(data_object, function(df) {
  colnames(df) <- gsub("\\.", "-", colnames(df))
  return(df)
})

# TCGA additional data #####
# Download clinical data for the TCGA samples of interest ###
library(TCGAbiolinks)
tcga_samples = Reduce(intersect, lapply(data_object, colnames))
tcga_samples = tcga_samples[2:length(tcga_samples)] # Remove X
tcga_samples = str_sub(tcga_samples, 1, 12)

# Create a query to retrieve clinical data for the specified samples
query <- GDCquery(
  project = "TCGA-BRCA",  # Replace with the appropriate TCGA project ID
  data.category = "Clinical",
  data.type = "Clinical Supplement",
  barcode = tcga_samples
)

# Execute the query
GDCdownload(query)

# Prepare the clinical data
clinical_data = GDCprepare_clinic(query, clinical.info = "patient")
openxlsx::write.xlsx(clinical_data, "Resources/BRCA complete/clinical_data.xlsx")
gc()

# Download categorical mutation data ###
library(maftools)
query_mut = GDCquery(
  project = "TCGA-BRCA",  # Replace with the appropriate TCGA project ID
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking",
  barcode = tcga_samples
)

# Execute the query
GDCdownload(query_mut)

# Prepare the clinical data
mut_data = GDCprepare(query_mut)
mut_data = mut_data %>% maftools::read.maf()

# The @data object is filtered for silent mutations. We use it to create a 
# binary matrix

# Split for SNP and indels
SNP_data = mut_data@data %>% dplyr::filter(Variant_Type == "SNP")
INDEL_data = mut_data@data %>% dplyr::filter(Variant_Type %in% c("DEL", "INS"))

SNP_matrix_input = SNP_data %>% dplyr::select(Hugo_Symbol, Start_Position, Tumor_Sample_Barcode)
str_sub(SNP_matrix_input$Tumor_Sample_Barcode, 16, -1) = ""

INDEL_matrix_input = INDEL_data %>% 
  dplyr::select(Hugo_Symbol, Start_Position, End_Position, Tumor_Sample_Barcode)
str_sub(INDEL_matrix_input$Tumor_Sample_Barcode, 16, -1) = ""

# Set the type of features
mut_features = "Gene" # can also be "Gene+Position"

if (mut_features == "Gene") {
  #SNPs
  SNP_matrix_input = SNP_matrix_input %>% dplyr::select(Hugo_Symbol,
                                                        Tumor_Sample_Barcode) %>%
    distinct()
  SNP_matrix <- dcast(SNP_matrix_input, 
                      Hugo_Symbol ~ Tumor_Sample_Barcode, fun.aggregate = NULL, 
                      value.var = "Hugo_Symbol")
  
  # INDELs
  INDEL_matrix_input = INDEL_matrix_input %>% dplyr::select(Hugo_Symbol,
                                                            Tumor_Sample_Barcode) %>%
    distinct()
  INDEL_matrix <- dcast(INDEL_matrix_input, 
                        Hugo_Symbol ~ Tumor_Sample_Barcode, fun.aggregate = NULL, 
                        value.var = "Hugo_Symbol")
  
} else if (mut_features == "Gene+Position") {
  # SNPs
  SNP_matrix_input$feature = paste0(SNP_matrix_input$Hugo_Symbol, "_", 
                                    SNP_matrix_input$Start_Position)
  SNP_matrix <- dcast(SNP_matrix_input, 
                      feature ~ Tumor_Sample_Barcode, fun.aggregate = NULL, 
                      value.var = "feature")
  
  # INDELs
  INDEL_matrix_input$feature = paste0(INDEL_matrix_input$Hugo_Symbol, "_", 
                                      INDEL_matrix_input$Start_Position, "_",
                                      INDEL_matrix_input$End_Position)
  INDEL_matrix <- dcast(INDEL_matrix_input, 
                        feature ~ Tumor_Sample_Barcode, fun.aggregate = NULL, 
                        value.var = "feature")
}

# Set rownames
SNP_matrix = as.matrix(SNP_matrix)
rownames(SNP_matrix) = SNP_matrix[, 1]
SNP_matrix = SNP_matrix[, 2:ncol(SNP_matrix)]

INDEL_matrix = as.matrix(INDEL_matrix)
rownames(INDEL_matrix) = INDEL_matrix[, 1]
INDEL_matrix = INDEL_matrix[, 2:ncol(INDEL_matrix)]

# Replace NAs with 0 and gene names with 1
SNP_matrix[is.na(SNP_matrix)] <- 0
SNP_matrix[SNP_matrix != 0] <- 1

INDEL_matrix[is.na(INDEL_matrix)] <- 0
INDEL_matrix[INDEL_matrix != 0] <- 1

# Numeric matrices
class(SNP_matrix) = "numeric"
class(INDEL_matrix) = "numeric"

# Pick one matrix for MOVICS analysis
data_object[["SNPs"]] = as.data.frame(SNP_matrix)

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
feature_column = c("X", "X", "X", "X", FALSE)

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

# Keep female samples that have measurements in all modalities
male_samples = paste0(clinical_data$bcr_patient_barcode[clinical_data$gender == "MALE"],
                      "-01")
if (alg_feature_pref == "rows") {
  # Sample names are in the columns
  overlap = setdiff(Reduce(intersect, lapply(input, colnames)),
                    male_samples)
  
  # Filter inputs
  input = lapply(input, function(x) {
    x = x[, overlap]
  })
} else {
  # Sample names are in the rows
  overlap = setdiff(Reduce(intersect, lapply(input, rownames)),
                    male_samples)
  
  # Filter inputs
  input = lapply(input, function(x) {
    x = x[overlap, ]
  })
}

# Run algorithm #####
library(MOVICS)

# identify optimal clustering number (may take a while)
optk = getClustNum(data = input,
                   is.binary = c(F,F,F,F,T),
                   try.N.clust = 2:10,
                   center = FALSE,
                   scale = FALSE, # default: FALSE
                   fig.path = "Results/MOVICS_baseline",
                   fig.name = paste0("optimal_k_plot_", data_source,
                                     "_", data_types))

# Perform multi-omic clustering with 9 available methods using default parameters
# iClusterBayes will be run on the cluster
moic.res.list = getMOIC(data = input,
                        methodslist = list("SNF", "CIMLR", "PINSPlus", "NEMO", 
                                           "COCA", # "MoCluster",
                                           "LRAcluster", "ConsensusClustering", 
                                           "IntNMF"),
                                           # , "iClusterBayes"),
                        N.clust = optk$N.clust,
                        type = c("gaussian", "gaussian", "gaussian", "gaussian",
                                 "binomial"))

# Save results to local file
save(moic.res.list, file = paste0(home, "/Results/MOVICS_baseline/", 
                                  algorithm, "_", data_source, "_",
                                  data_types, "_eval_on_", evaluation_source,
                                  "_moic.res.list.rda"))

# iClusterBayes with lower burnin and draw parameter values
iClusterBayes.res = getMOIC(data        = input,
                            N.clust     = 3,
                            methodslist = "iClusterBayes",
                            type        = c("gaussian",
                                            "gaussian",
                                            "gaussian",
                                            "gaussian",
                                            "binomial"),
                            n.burnin    = 1800,
                            n.draw      = 1200,
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
                            N.clust     = 3,
                            methodslist = "MoCluster",
                            type        = c("gaussian",
                                            "gaussian",
                                            "gaussian",
                                            "gaussian",
                                            "binomial"),
                            center      = FALSE,
                            scale       = FALSE)
moic.res.list = append(moic.res.list, 
                       list("MoCluster" = moCluster.res))
rm(moCluster.res); gc()

# Save final moic.res.list object
save(moic.res.list, file = paste0(home, "/Results/MOVICS_baseline/", 
                                  algorithm, "_", data_source, "_",
                                  data_types, "_eval_on_", evaluation_source,
                                  "_moic.res.list.rda"))

# Inspect output for similarities/differences across clusterings #####
library(dplyr)

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
getSilhouette(sil      = consensus$sil,
              fig.path = "Results/MOVICS_baseline",
              fig.name = "Silhouette",
              height   = 5.5,
              width    = 5)

# Downstream comparisons #####

# Create data frame for comparisons
cc_res = consensus$clust.res
colnames(cc_res) = c("Sample.ID", "Consensus Subtype")
cc_res$`Consensus Subtype` = factor(cc_res$`Consensus Subtype`,
                                    levels = c(1, 2, 3),
                                    labels = c("CS1", "CS2", "CS3"))
subset_of_interest = clinical_data[, c(1, 7:12, 16, 22, 24, 26, 38, 40, 41,
                                       51:54, 67:69, 107)] %>%
  dplyr::rename(Sample.ID = bcr_patient_barcode) %>%
  dplyr::mutate(Sample.ID = paste0(Sample.ID, "-01")) %>%
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

# Statistical comparisons
clin_comp = compClinvar(moic.res = consensus,
                        var2comp = var2comp,
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
                        includeNA = FALSE,
                        doWord = TRUE,
                        tab.name = "Summary_of_clinical_variables",
                        res.path = paste0(home, "/Results/MOVICS_baseline/"))

# Color annotation #####
# set color for each omics data
# if no color list specified all subheatmaps will be unified to green and red color pattern
mRNA.col   <- c("#00ff00", "#000000", "#ff0000")
CNV.col <- c("#6699CC", "white", "#FF3C38")
mut.col   <- c("#EFE5AF", "white", "#780A43")
methylation.col    <- c("#57087C", "#000000", "#FF3C38")
miRNA.col <- c("#D3ACEF", "#000000", "#FA076B")
col.list   <- list(CNV.col, methylation.col, mRNA.col, miRNA.col, mut.col)

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
# Save environment
save.image(paste0(home, "/Results/MOVICS_baseline/", 
                  algorithm, "_", data_source, "_",
                  data_types, "_eval_on_", evaluation_source,
                  "_env.RData"))
