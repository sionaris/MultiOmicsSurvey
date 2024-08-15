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


