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

# Generate hyperparameters report text
generate_hyperparams_text = function (algorithm, hyperparameters) {
  hyperparameters_text = paste(readLines(paste0("Resources/algorithm_descriptions/", algorithm,
                         "_hyperparams_text.Rmd")),
        collapse = "\n")
  #hyperparameters_text = glue::glue(hyperparameters_text)
  return(hyperparameters_text)
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
