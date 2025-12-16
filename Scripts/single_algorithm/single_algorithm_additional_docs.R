# single_algorithm_additional_docs.R
# ==============================================================================
# Roxygen-style documentation for custom functions defined in the single_algorithm/
# scripts. This file serves as a centralized reference and does NOT modify the
# original source files. Each function is documented with its source file(s).
# ==============================================================================

################################################################################
# MAD-BASED WEIGHTING FUNCTIONS
################################################################################

#' Compute MAD-Based Feature Weights
#'
#' @description
#' Calculates Median Absolute Deviation (MAD) based weights for features in
#' a data matrix. Used to assign higher importance to features with greater
#' variability across samples.
#'
#' @param X Numeric matrix with features in rows and samples in columns.
#'
#' @return Named numeric vector of normalized weights (sum = 1), with names
#'   corresponding to row names of the input matrix. Features with zero or

#'   non-finite MAD receive zero weight; if all features are constant, uniform
#'   weights are assigned.
#'
#' @details
#' The function computes MAD for each feature (row), handles non-finite values
#' by setting them to zero, and normalizes so that all weights sum to 1.
#'
#' @seealso Used for weighted distance calculations in similarity network
#'   fusion algorithms.
#'
#' @source Found in:
#'   - Scripts/single_algorithm/ab-SNF.R (line ~108)
#'   - Scripts/single_algorithm/wMKL.R (line ~160)
mad_based_weights <- function(X) {
  v <- apply(X, 2, mad, na.rm = TRUE) # MAD of each feature (features in columns)
  v[!is.finite(v)] <- 0
  if (sum(v) == 0) {
    # all features are constant, fallback: uniform weighting or return zeros
    w <- rep(0, length(v))
  } else {
    w <- v / sum(v) # normalize so sum of all weights is 1
  }
  return(w)
}


#' Compute Weighted Euclidean Distance Matrix
#'
#' @description
#' Calculates a weighted Euclidean distance matrix between samples, where each
#' feature contribution is scaled by the square root of its weight.
#'
#' @param X Numeric matrix with features in rows and samples in columns.
#' @param w Named numeric vector of feature weights (length = nrow(X)).
#'
#' @return A square distance matrix of dimension ncol(X) x ncol(X).
#'
#' @details
#' For features with non-zero weights, the distance calculation uses
#' sqrt(w) * (x_i - x_j), then aggregates via sum of squares.
#'
#' @source Found in:
#'   - Scripts/single_algorithm/ab-SNF.R (line ~120)
weighted_euclidean_dist <- function(X, w) {
  Xw <- sweep(X, 2, sqrt(w), `*`) # multiply each feature column by sqrt(w_k)
  dmat <- abSNF::dist2(as.matrix(Xw), as.matrix(Xw))
  return(dmat)
}


#' Compute Weighted Hamming Fraction Distance
#'
#' @description
#' Calculates a weighted Hamming distance matrix for binary data, assigning
#' different weights to driver genes vs. non-driver genes.
#'
#' @param binary_mat Binary (0/1) matrix with features in rows and samples in columns.
#' @param features Character vector of feature names (rownames of binary_mat).
#' @param drivers Character vector of driver gene identifiers.
#' @param driver_weight Numeric weight for driver gene mismatches.
#' @param nondriver_weight Numeric weight for non-driver gene mismatches.
#'
#' @return A square distance matrix of dimension ncol(binary_mat) x ncol(binary_mat).
#'
#' @details
#' Computes weighted Hamming distance as the fraction of disagreeing positions,
#' with driver genes weighted more heavily than non-drivers.
#'
#' @source Found in:
#'   - Scripts/single_algorithm/ab-SNF.R (line ~133)
weighted_hamming_fraction <- function(binary_mat, features, drivers,
                                      driver_weight, nondriver_weight) {
  w <- ifelse(features %in% drivers, driver_weight, nondriver_weight)
  
  # Compute for each pair (i, j):
  #   d(i, j) = [ sum_{k}  w_k * |SNP_mat[i,k] - SNP_mat[j,k]| ] / sum_{k} w_k
  # This yields the *fraction* of mismatches, weighted by w_k.
  
  total_w <- sum(w)
  N <- nrow(binary_mat)
  dmat <- matrix(0, nrow = N, ncol = N)
  
  for (i in seq_len(N)) {
    for (j in seq_len(N)) {
      # Weighted sum of absolute differences across features
      mismatches_ij <- sum(w * abs(binary_mat[i, ] - binary_mat[j, ]))
      # Divide by total weight to get fraction of weighted mismatches
      dmat[i, j] <- mismatches_ij / total_w
    }
  }
  
  return(dmat)
}


#' Compute Weighted Mutation Weights
#'
#' @description
#' Assigns weights to features based on whether they are known driver genes,
#' normalizing so the total weight sums to 1.
#'
#' @param features Character vector of all feature names.
#' @param drivers Character vector of known driver gene identifiers.
#' @param driver_weight Numeric weight to assign to driver genes.
#' @param nondriver_weight Numeric weight to assign to non-driver genes.
#'
#' @return Numeric vector of normalized weights (length = length(features), sum = 1).
#'
#' @source Found in:
#'   - Scripts/single_algorithm/wMKL.R (line ~176)
weighted_mutations <- function(features, drivers, driver_weight, nondriver_weight) {
  w <- rep(NA, length(features))
  w[which(features %in% drivers)] = driver_weight
  w[which(!features %in% drivers)] = nondriver_weight
  w <- w / sum(w) # normalize so sum of all weights is 1
  return(w)
}


################################################################################
# MATRIX SIMILARITY FUNCTIONS
################################################################################

#' Compute Pairwise Matrix Similarity Metrics
#'
#' @description
#' Calculates both Frobenius norm distance and Pearson correlation between all
#' pairs of matrices in a list. Used to assess similarity between affinity
#' matrices computed under different hyperparameter settings.
#'
#' @param matrices Named list of numeric matrices of identical dimensions.
#'
#' @return List with two elements:
#'   \describe{
#'     \item{Frobenius}{Square matrix of pairwise Frobenius norm distances}
#'     \item{Pearson}{Square matrix of pairwise Pearson correlations}
#'   }
#'   Row and column names correspond to the names of the input matrices.
#'
#' @details
#' Diagonal elements are zero. Relies on helper functions `frobenius_norm()`
#' and `pearson_correlation()` for individual pairwise computations.
#'
#' @source Found in:
#'   - Scripts/single_algorithm/ab-SNF.R (line ~254)
#'   - Scripts/single_algorithm/ANF.R (line ~173)
#'   - Scripts/single_algorithm/RWR-F.R
#'   - Scripts/single_algorithm/RWR-NF.R
#'   - Scripts/single_algorithm/SNF.R
#'   - Scripts/single_algorithm/Spectrum.R
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


################################################################################
# SPECTRAL CLUSTERING FUNCTIONS
################################################################################

#' Spectral Clustering with Eigenvector Extraction
#'
#' @description
#' Performs spectral clustering on an affinity matrix and returns both cluster
#' assignments and the eigenvector embedding. This is a modified version of
#' `SNFtool::spectralClustering` that exposes the intermediate eigen-space.
#'
#' @param affinity Square symmetric affinity matrix (n x n).
#' @param K Integer, the number of clusters.
#' @param type Integer (1, 2, or 3) specifying the Laplacian normalization:
#'   - 1: Unnormalized Laplacian (L = D - A)
#'   - 2: Random walk Laplacian (D^{-1} L)
#'   - 3: Symmetric normalized Laplacian (D^{-1/2} L D^{-1/2}), default
#'
#' @return Data frame containing:
#'   \describe{
#'     \item{eig1...eigK}{Numeric columns for each of the K eigenvectors}
#'     \item{Cluster}{Integer cluster assignment (1 to K)}
#'     \item{Sample.ID}{Character sample identifiers from affinity column names}
#'   }
#'
#' @details
#' Computes the graph Laplacian, extracts the K smallest eigenvectors (excluding
#' the trivial eigenvector), row-normalizes if type = 3, then applies k-means
#' via `SNFtool:::.discretisation` to obtain discrete cluster assignments.
#'
#' @source Found in:
#'   - Scripts/single_algorithm/ab-SNF.R (line ~758)
#'   - Scripts/single_algorithm/SNF.R (line ~693)
#'   - Scripts/single_algorithm/NEMO.R (line ~575)
spectralClustering_eig <- function(affinity, K, type = 3) {
  library(cluster)
  d = rowSums(affinity)
  d[d == 0] = .Machine$double.eps
  D = diag(d)
  L = D - affinity
  if (type == 1) {
    NL = L
  }
  else if (type == 2) {
    Di = diag(1/d)
    NL = Di %*% L
  }
  else if (type == 3) {
    Di = diag(1/sqrt(d))
    NL = Di %*% L %*% Di
  }
  eig = eigen(NL)
  res = sort(abs(eig$values), index.return = TRUE)
  U = eig$vectors[, res$ix[1:K]]
  normalize <- function(x) x/sqrt(sum(x^2))
  if (type == 3) {
    U = t(apply(U, 1, normalize))
  }
  eigDiscrete = SNFtool:::.discretisation(U)
  eigDiscrete = eigDiscrete$discrete
  labels = apply(eigDiscrete, 1, which.max)
  U = as.data.frame(cbind(U, labels))
  U$Sample.ID = colnames(affinity)
  colnames(U)[(ncol(U)-1):ncol(U)] = c("Cluster", "Sample.ID")
  return(U)
}


#' ANF-Style Spectral Clustering with Eigenvector Extraction
#'
#' @description
#' Adapted spectral clustering for ANF that uses a random-walk normalized
#' Laplacian and returns eigenvectors along with cluster assignments.
#'
#' @param affinity Square symmetric affinity matrix (n x n).
#' @param type Character, normalization type. Currently only "rw" (random walk)
#'   is implemented.
#' @param optk Integer, the number of clusters.
#'
#' @return Data frame containing:
#'   \describe{
#'     \item{eigenvector columns}{Numeric columns for the optk eigenvectors}
#'     \item{Cluster}{Integer cluster assignment}
#'     \item{Sample.ID}{Character sample identifiers}
#'   }
#'
#' @details
#' Uses `ANF::pod` for discretization instead of SNFtool's internal method.
#'
#' @source Found in:
#'   - Scripts/single_algorithm/ANF.R (line ~315)
spectral_clustering_eig <- function(affinity, type = "rw", optk) {
  # Adapted source code of ANF::spectral_clustering to return eigenvectors as well
  library(cluster)
  n = nrow(affinity)
  d = rowSums(affinity)
  d[d == 0] = .Machine$double.eps
  L = diag(d) - affinity
  NL = diag(1/d) %*% L
  eig = eigen(NL)
  Y = Re(eig$vectors[, -1:-(n - optk)])
  labels = apply(pod(Y), 1, which.max)
  U = as.data.frame(cbind(Y, labels))
  U$Sample.ID = colnames(affinity)
  colnames(U)[(ncol(U)-1):ncol(U)] = c("Cluster", "Sample.ID")
  return(U)
}


#' Extract Eigenspace from Affinity Matrix
#'
#' @description
#' Extracts the spectral embedding (eigenvectors) and cluster labels from an
#' affinity matrix. Similar to `spectralClustering_eig` but uses `wMKL`
#' discretization internally.
#'
#' @param affinity Square symmetric affinity matrix (n x n).
#' @param K Integer, the number of clusters.
#' @param type Integer (1, 2, or 3) specifying Laplacian normalization.
#'
#' @return List with elements:
#'   \describe{
#'     \item{labels}{Integer vector of cluster assignments}
#'     \item{U}{Matrix of eigenvectors (n x K)}
#'     \item{eigDiscrete}{Discrete cluster indicator matrix}
#'   }
#'
#' @source Found in:
#'   - Scripts/single_algorithm/wMKL.R (line ~340)
extract_eigenspace <- function(affinity, K, type = 3) {
  d <- rowSums(affinity)
  d[d == 0] <- .Machine$double.eps
  D <- diag(d)
  L <- D - affinity
  if (type == 1) {
    NL <- L
  }
  else if (type == 2) {
    Di <- diag(1/d)
    NL <- Di %*% L
  }
  else if (type == 3) {
    Di <- diag(1/sqrt(d))
    NL <- Di %*% L %*% Di
  }
  eig <- eigen(NL)
  res <- sort(abs(eig$values), index.return = TRUE)
  U <- eig$vectors[, res$ix[1:K]]
  normalize <- function(x) x/sqrt(sum(x^2))
  if (type == 3) {
    U <- t(apply(U, 1, normalize))
  }
  eigDiscrete <- wMKL:::.discretisation(U)
  eigDiscrete <- eigDiscrete$discrete
  labels <- apply(eigDiscrete, 1, which.max)
  return(list(labels = labels, U = U, eigDiscrete = eigDiscrete))
}


################################################################################
# STATISTICAL ANALYSIS FUNCTIONS
################################################################################

#' Bias-Corrected Cramer's V Test
#'
#' @description
#' Computes the bias-corrected Cramer's V (or Phi coefficient for 2x2 tables)
#' as a measure of association between two categorical variables.
#'
#' @param x Contingency table (matrix or table object).
#' @param string Character description of the comparison for output labeling.
#' @param digits Integer, number of decimal places to round to (default: 3).
#'
#' @return List with two elements:
#'   \describe{
#'     \item{text}{Formatted string describing the result}
#'     \item{value}{Rounded numeric Cramer's V value}
#'   }
#'
#' @details
#' Uses `rcompanion::cramerV` with bias correction enabled. This provides a
#' more accurate effect size estimate especially for smaller samples.
#'
#' @source Found in (identical definition in all files):
#'   - Scripts/single_algorithm/ab-SNF.R (line ~2142)
#'   - Scripts/single_algorithm/ANF.R (line ~1688)
#'   - Scripts/single_algorithm/CIMLR.R (line ~1574)
#'   - Scripts/single_algorithm/COCA.R (line ~1192)
#'   - Scripts/single_algorithm/iClusterBayes.R (line ~1282)
#'   - Scripts/single_algorithm/KLIC.R (line ~1445)
#'   - Scripts/single_algorithm/LRAcluster.R (line ~1407)
#'   - Scripts/single_algorithm/MDICC.R (line ~1438)
#'   - Scripts/single_algorithm/MFA.R (line ~1860)
#'   - Scripts/single_algorithm/MOFA.R (line ~1909)
#'   - Scripts/single_algorithm/MONET.R (line ~1641)
#'   - Scripts/single_algorithm/MSNE.R (line ~1209)
#'   - Scripts/single_algorithm/NEMO.R (line ~1908)
#'   - Scripts/single_algorithm/RWR-F.R (line ~1971)
#'   - Scripts/single_algorithm/RWR-NF.R (line ~1970)
#'   - Scripts/single_algorithm/SNF.R (line ~2078)
#'   - Scripts/single_algorithm/Spectrum.R (line ~1392)
#'   - Scripts/single_algorithm/wMKL.R (line ~1494)
unbiased.cv.test <- function(x, string, digits = 3) {
  CV = rcompanion::cramerV(x, bias.correct = TRUE)
  return(list(text = paste0("Bias-corrected Cramer's V / Phi for ", 
                            string, ": ", round(as.numeric(CV), digits)),
              value = round(as.numeric(CV), digits)))
}


################################################################################
# GRAPH/NETWORK ANALYSIS FUNCTIONS
################################################################################

#' Determine Edge Modality Support
#'
#' @description
#' Identifies which data modalities primarily support a given edge in a fused
#' affinity network, based on relative edge weights across modality-specific
#' networks.
#'
#' @param edge_weights Numeric vector of edge weights from each modality network.
#'
#' @return Integer vector of indices indicating which modalities support the edge:
#'   - Single index if one modality dominates (>10% higher than others)
#'   - Multiple indices if several modalities contribute similarly
#'   - All indices if all modalities contribute within 10% of the maximum
#'
#' @details
#' Uses a 10% threshold to distinguish dominant modalities from those with
#' similar contributions. Useful for visualizing which data types drive
#' specific sample relationships in the fused network.
#'
#' @source Found in:
#'   - Scripts/single_algorithm/ab-SNF.R (line ~2407)
#'   - Scripts/single_algorithm/ANF.R (line ~1952)
#'   - Scripts/single_algorithm/MDICC.R (line ~1705)
#'   - Scripts/single_algorithm/RWR-F.R (line ~2230)
#'   - Scripts/single_algorithm/RWR-NF.R (line ~2235)
#'   - Scripts/single_algorithm/SNF.R (line ~2342)
determine_edge_support <- function(edge_weights) {
  sorted_weights <- sort(edge_weights, decreasing = TRUE)
  # If the highest weight is more than 10% greater than all others, it is supported by a single modality
  if (sorted_weights[1] > sorted_weights[2] * 1.1) {
    return(which(edge_weights == sorted_weights[1]))
  }
  # If the difference between the two highest weights is less than 10%, it is supported by those two modalities
  else if (sorted_weights[1] <= sorted_weights[2] * 1.1 && sorted_weights[2] > sorted_weights[3] * 1.1) {
    return(which(edge_weights >= sorted_weights[2]))
  }
  # If the difference between all weights is less than 10%, it is supported by all modalities
  else {
    return(which(edge_weights >= sorted_weights[1] * 0.9))
  }
}


################################################################################
# MONET-SPECIFIC FUNCTIONS
################################################################################

#' Apply Offset to Correlation Matrix
#'
#' @description
#' Transforms a correlation matrix by subtracting the global mean and an
#' additional offset, then setting diagonal to zero. This is the core
#' preprocessing step for MONET's network construction.
#'
#' @param cor.matrix Square correlation matrix.
#' @param offset Numeric offset value to subtract (e.g., 0.2 as in MONET paper).
#'
#' @return Transformed correlation matrix with diagonal set to zero.
#'
#' @details
#' The formula is: weighted = cor.matrix - mean(cor.matrix) - offset
#' This centers the distribution and shifts it so that only stronger-than-average
#' correlations remain positive.
#'
#' @source Found in:
#'   - Scripts/single_algorithm/MONET.R (line ~106)
offset_correlation <- function(cor.matrix, offset) {
  avg_val = mean(mean(cor.matrix))
  weighted = cor.matrix - avg_val - offset
  diag(weighted) = 0
  return(weighted)
}


#' Compute Mean Absolute Edge Weight
#'
#' @description
#' Calculates the mean of absolute edge weights in a NetworkX graph object.
#' Used to establish baseline connectivity levels for MONET modules.
#'
#' @param nx_graph A reticulate-wrapped NetworkX graph object with weighted edges.
#'
#' @return Numeric scalar, the mean absolute weight across all edges.
#'
#' @details
#' Iterates over all edges in the graph, extracts absolute weights, and
#' computes their arithmetic mean.
#'
#' @source Found in:
#'   - Scripts/single_algorithm/MONET.R (line ~450)
edge_absolute_mean <- function(nx_graph) {
  w <- vapply(iterate(nx_graph$edges(data = TRUE)),
              function(e) abs(as.numeric(e[[3]]$weight)),
              numeric(1))
  mean(w)  # just return the scalar mean
}


#' Get Module-Specific Mean Absolute Edge Weight
#'
#' @description
#' Computes the mean absolute edge weight within a MONET module for a specific
#' omic network. Used to quantify how strongly each omic supports a given module.
#'
#' @param mod MONET module object with patient names accessible via
#'   `get_patients_names_as_list()`.
#' @param omic_name Character name of the omic (e.g., "RNAseq", "CNV").
#'
#' @return Numeric scalar, the mean absolute weight of edges between module
#'   members in the specified omic network. Returns NA if the module has fewer
#'   than 2 members.
#'
#' @details
#' Examines all pairwise combinations of samples in the module and retrieves
#' edge weights from the corresponding omic graph, computing the mean of
#' absolute values.
#'
#' @source Found in:
#'   - Scripts/single_algorithm/MONET.R (line ~464)
get_mod_edge_absmean <- function(mod, omic_name) {
  g   <- MONET$glob_var$omics[[omic_name]]$graph
  pts <- unlist(mod$get_patients_names_as_list())
  if (length(pts) < 2) return(NA)
  
  combs <- t(combn(pts, 2))
  w <- apply(combs, 1, function(x)
    if (g$has_edge(x[1], x[2]))
      abs(as.numeric(g$get_edge_data(x[1], x[2])$weight))
    else 0)
  mean(w)
}


################################################################################
# COCA-SPECIFIC FUNCTIONS
################################################################################
  
#' Hierarchical Clustering Gap Statistic Wrapper
#'
#' @description
#' Wrapper function for use with `cluster::clusGap` that returns precomputed
#' hierarchical clustering results for a given k.
#'
#' @param x Data matrix (not directly used; gap statistic framework requirement).
#' @param k Integer, the number of clusters to return.
#'
#' @return Clustering result object from `list_of_k[[paste0("k = ", k)]]`,
#'   which must be predefined in the calling environment.
#'
#' @details
#' This function serves as a `FUNcluster` argument to `clusGap`, allowing
#' precomputed clusterings (e.g., from hierarchical cutting at various k)
#' to be evaluated via the gap statistic framework.
#'
#' @source Found in:
#'   - Scripts/single_algorithm/COCA.R (line ~136)
hclust_gap <- function(x, k) {
  return(list_of_k[[paste0("k = ", k)]])
}


################################################################################
# PAMOGK PREPROCESSING FUNCTIONS
################################################################################

#' Clean Row Names by Removing Version Suffixes
#'
#' @description
#' Removes version number suffixes (e.g., ".1", ".2") from row names of a data
#' frame. Commonly needed when dealing with gene/transcript identifiers.
#'
#' @param df Data frame with row names containing version suffixes.
#'
#' @return Data frame with cleaned row names.
#'
#' @source Found in:
#'   - Scripts/single_algorithm/PAMOGK_prep.R (line ~15)
clean_rownames <- function(df) {
  rownames(df) <- sub("\\.\\d+$", "", rownames(df))
  return(df)
}


#' Remove Rows with NA Row Names
#'
#' @description
#' Filters out rows that have NA as their row name, with reporting of how many
#' rows are removed.
#'
#' @param df Data frame potentially containing rows with NA row names.
#'
#' @return Data frame with only rows having non-NA row names.
#'
#' @details
#' Prints a message indicating the number of rows removed if any have NA names.
#'
#' @source Found in:
#'   - Scripts/single_algorithm/PAMOGK_prep.R (line ~30)
remove_na_rownames <- function(df) {
  # Identify rows with non-NA row names
  valid_rows <- !is.na(rownames(df))
  
  # Report the number of rows being removed
  removed_rows <- sum(!valid_rows)
  if (removed_rows > 0) {
    message(sprintf("Removing %d rows with NA row names.", removed_rows))
  }
  
  # Subset the data frame to keep only valid rows
  cleaned_df <- df[valid_rows, , drop = FALSE]
  
  return(cleaned_df)
}


#' Get Dimensions or NULL Status
#'
#' @description
#' Returns a formatted string describing the dimensions of an object, or "NULL"
#' if the object is NULL.
#'
#' @param x Object to check dimensions of.
#'
#' @return Character string: either "NULL" or "nrow rows x ncol columns".
#'
#' @source Found in:
#'   - Scripts/single_algorithm/PAMOGK_prep.R (line ~56)
get_dimensions <- function(x) {
  if (is.null(x)) {
    return("NULL")
  } else {
    dims <- dim(x)
    return(paste(dims[1], "rows x", dims[2], "columns"))
  }
}


#' Fix UniProt Row Names from List Format
#'
#' @description
#' Converts row names that were stored as R vector strings (e.g., 
#' 'c("A8K052","P04217")') back into semicolon-separated strings.
#'
#' @param df Data frame with problematic row name format.
#'
#' @return Data frame with row names converted to semicolon-separated format.
#'
#' @details
#' Parses row names that look like `c("ID1","ID2",...)` using `eval(parse())`
#' and collapses them with semicolons.
#'
#' @source Found in:
#'   - Scripts/single_algorithm/PAMOGK_prep.R (line ~88)
fix_uniprot_row_names <- function(df) {
  old_rownames <- rownames(df)
  new_rownames <- sapply(old_rownames, function(x) {
    # If it literally starts with c("..."), parse it:
    # Safeguard: if the row name doesn't match that pattern, skip or return it as-is.
    if (length(strsplit(x, split = ",")[[1]]) == 1) {
      paste(x)
    } else {
      # Otherwise parse:
      val <- eval(parse(text = x))  # yields a character vector c("A8K052","P04217",...)
      paste(val, collapse = ";")    # turn c("A8K052","P04217") into "A8K052;P04217"
    }
  })
  
  rownames(df) <- new_rownames
  df
}


#' Filter to Valid PAMOGK UniProt IDs
#'
#' @description
#' Filters data frames to retain only rows with UniProt IDs that exist in
#' PAMOGK's pathway/interaction database, handling multi-ID row names.
#'
#' @param up_input List of data frames with UniProt IDs as row names.
#' @param PAMOGK_UP_IDs Character vector of valid UniProt IDs from PAMOGK.
#'
#' @return List of filtered data frames, with duplicates resolved by:
#'   - For SNPs: keeping the row with highest sum (most mutations)
#'   - For continuous data: keeping the row with highest MAD
#'
#' @details
#' For row names containing multiple semicolon-separated IDs, only the first
#' valid ID (matching PAMOGK_UP_IDs) is retained.
#'
#' @source Found in:
#'   - Scripts/single_algorithm/PAMOGK_prep.R (line ~142)
filter_uniprot_ids <- function(up_input, PAMOGK_UP_IDs) {
  filtered_up_input <- lapply(names(up_input), function(datatype) {
    df <- up_input[[datatype]]
    filtered_rows <- sapply(rownames(df), function(row) {
      ids <- strsplit(row, ";")[[1]]
      valid_ids <- intersect(ids, PAMOGK_UP_IDs)
      if (length(valid_ids) > 0) {
        return(valid_ids[1])  # Keep only the first valid ID
      } else {
        return(NA)
      }
    })
    
    valid_rows <- !is.na(filtered_rows)
    df_filtered <- df[valid_rows, ]
    rownames(df_filtered) <- filtered_rows[valid_rows]
    
    # Handle duplicates
    if (datatype == "SNPs") {
      # For SNPs, keep the row with the highest sum
      df_filtered <- df_filtered[order(rowSums(df_filtered), decreasing = TRUE), ]
      df_filtered <- df_filtered[!duplicated(rownames(df_filtered)), ]
    } else {
      # For continuous data, keep the row with the highest MAD
      row_mads <- matrixStats::rowMads(as.matrix(df_filtered))
      df_filtered <- df_filtered[order(row_mads, decreasing = TRUE), ]
      df_filtered <- df_filtered[!duplicated(rownames(df_filtered)), ]
    }
    
    return(df_filtered)
  })
  
  names(filtered_up_input) <- names(up_input)
  return(filtered_up_input)
}


################################################################################
# iCLUSTERBAYES TUNING FUNCTIONS
################################################################################

#' Prepare Global Data for iClusterBayes Visualization
#'
#' @description
#' Prepares a tidy data frame of Z.ar (latent variable acceptance rates) across
#' all K values for density plotting.
#'
#' @param data_list List of numeric vectors (Z.ar values for each K).
#' @param K_values Integer vector of K indices corresponding to data_list.
#' @param var_name Character label for the variable (e.g., "Z.ar").
#'
#' @return Data frame with columns: value, K (factor), variable.
#'
#' @source Found in:
#'   - Scripts/single_algorithm/iCB_inspect_HPC_tune_results.R (line ~28)
prepare_global_data <- function(data_list, K_values, var_name) {
  df <- do.call(rbind, lapply(seq_along(data_list), function(i) {
    data.frame(value = data_list[[i]], K = K_values[i] + 1)
  }))
  df$K <- as.factor(df$K)
  df$variable <- var_name
  return(df)
}


#' Prepare Modality-Specific Data for iClusterBayes Visualization
#'
#' @description
#' Prepares a tidy data frame of Beta/Gamma acceptance rates broken down by
#' modality and K value for visualization.
#'
#' @param data_list List of lists, where each inner list contains vectors per
#'   modality.
#' @param K_values Integer vector of K indices corresponding to data_list.
#' @param var_name Character label for the variable (e.g., "beta.ar").
#'
#' @return Data frame with columns: value, modality, K, variable.
#'
#' @details
#' Assumes `modalities` variable exists in the calling environment.
#'
#' @source Found in:
#'   - Scripts/single_algorithm/iCB_inspect_HPC_tune_results.R (line ~38)
prepare_modality_data <- function(data_list, K_values, var_name) {
  do.call(rbind, lapply(seq_along(data_list), function(i) {
    data.frame(
      value    = unlist(data_list[[i]]),
      modality = rep(modalities, sapply(data_list[[i]], length)),
      K        = K_values[i] + 1,
      variable = var_name
    )
  }))
}


#' Get Summary Statistics Text Output
#'
#' @description
#' Generates formatted text output summarizing Z.ar, beta.ar, and gamma.ar
#' statistics across K values for iClusterBayes tuning results.
#'
#' @param fit_list List of fit results, indexed by K.
#' @param K_values Integer vector of K values to summarize.
#'
#' @return Character vector of captured output (one element per line).
#'
#' @source Found in:
#'   - Scripts/single_algorithm/iCB_inspect_HPC_tune_results.R (line ~50)
get_summary_statistics_text <- function(fit_list, K_values) {
  summary_text <- capture.output({
    for (K in K_values) {
      cat("\n--- Summary for K =", K, "---\n")
      
      if (!is.null(fit_list[[K]]$Z.ar)) {
        cat("Z.ar:\n")
        print(summary(fit_list[[K]]$Z.ar))
      }
      
      if (!is.null(fit_list[[K]]$beta.ar)) {
        cat("\nbeta.ar:\n")
        for (m in seq_along(fit_list[[K]]$beta.ar)) {
          cat("  Modality", m, ":\n")
          print(summary(fit_list[[K]]$beta.ar[[m]]))
        }
      }
      
      if (!is.null(fit_list[[K]]$gamma.ar)) {
        cat("\ngamma.ar:\n")
        for (m in seq_along(fit_list[[K]]$gamma.ar)) {
          cat("  Modality", m, ":\n")
          print(summary(fit_list[[K]]$gamma.ar[[m]]))
        }
      }
    }
  })
  return(summary_text)
}


#' Create Global Density Plot
#'
#' @description
#' Creates a ggplot2 density plot of Z.ar values colored by K.
#'
#' @param data Data frame from `prepare_global_data()`.
#' @param title Character plot title.
#'
#' @return ggplot2 object.
#'
#' @source Found in:
#'   - Scripts/single_algorithm/iCB_inspect_HPC_tune_results.R (line ~81)
create_global_plot <- function(data, title) {
  ggplot(data, aes(x = value, color = K, fill = K)) +
    geom_density(alpha = 0.3, linewidth = 0.2) +
    labs(title = title, x = "Value", y = "Density") +
    theme_minimal(base_size = 7) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom"
    )
}


#' Create Modality-Specific Density Plot
#'
#' @description
#' Creates a ggplot2 density plot of Beta/Gamma acceptance rates colored by
#' modality, for a specific K value.
#'
#' @param data Data frame from `prepare_modality_data()`.
#' @param K Integer, the K value to plot.
#' @param title Character plot title prefix.
#'
#' @return ggplot2 object.
#'
#' @details
#' Assumes `modality_colours` variable exists in the calling environment.
#'
#' @source Found in:
#'   - Scripts/single_algorithm/iCB_inspect_HPC_tune_results.R (line ~93)
create_modality_plot <- function(data, K, title) {
  ggplot(data[data$K == K+1, ], aes(x = value, color = modality, fill = modality)) +
    geom_density(alpha = 0.3, linewidth = 0.2) +
    labs(title = paste(title, "- K =", K+1), x = "Value", y = "Density") +
    scale_color_manual(values = modality_colours) +
    scale_fill_manual(values = modality_colours) +
    theme_minimal(base_size = 7) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom"
    )
}


#' Compute Robust Metrics and Score for iClusterBayes Tuning
#'
#' @description
#' Computes median acceptance rates and a composite score based on distance
#' from the optimal acceptance rate (0.234) for each K value.
#'
#' @param tune_results List containing `fit` element with results per K.
#' @param K_values Integer vector of K values to evaluate.
#'
#' @return List with:
#'   \describe{
#'     \item{per_K_metrics}{List of metrics for each K (median values and score)}
#'     \item{final_score}{Sum of scores across all K (lower is better)}
#'   }
#'
#' @details
#' The target acceptance rate of 0.234 comes from Roberts et al.'s work on
#' optimal scaling for MCMC algorithms. The score sums |median - 0.234| for
#' Z.ar, beta.ar, and gamma.ar.
#'
#' @references
#' Roberts, G.O. et al. (1997). Weak convergence and optimal scaling of random
#' walk Metropolis algorithms. Ann. Appl. Probab., 7(1), 110-120.
#'
#' @source Found in:
#'   - Scripts/single_algorithm/iCB_inspect_HPC_tune_results.R (line ~110)
compute_robust_metrics_and_score <- function(tune_results, K_values) {
  metrics_list <- list()
  final_score  <- 0
  
  for (K in K_values) {
    z_ar     <- tune_results$fit[[K]]$Z.ar
    beta_ar  <- tune_results$fit[[K]]$beta.ar
    gamma_ar <- tune_results$fit[[K]]$gamma.ar
    
    if (is.null(z_ar)) next
    
    median_z   <- median(z_ar, na.rm = TRUE)
    all_beta   <- unlist(beta_ar)
    median_beta <- median(all_beta, na.rm = TRUE)
    all_gamma  <- unlist(gamma_ar)
    median_gamma <- median(all_gamma, na.rm = TRUE)
    
    # Scoring: sum of distances from 0.234 for each median
    # why 0.234?: https://www.maths.lancs.ac.uk/~sherlocc/Publications/rwm.final.pdf
    # Lower is "better"
    score_k <- (abs(median_z - 0.234) +
                  abs(median_beta - 0.234) +
                  abs(median_gamma - 0.234))
    
    metrics_list[[K]] <- list(
      K               = K,
      median_z_ar     = median_z,
      median_beta_ar  = median_beta,
      median_gamma_ar = median_gamma,
      score_k         = score_k
    )
    final_score <- final_score + score_k
  }
  
  return(list(
    per_K_metrics = metrics_list,
    final_score   = final_score
  ))
}


#' Penalize Extreme Acceptance Rates
#'
#' @description
#' Returns a penalty if an acceptance rate falls outside acceptable bounds.
#' Used to discourage hyperparameter combinations that lead to poor MCMC mixing.
#'
#' @param acc_rate Numeric acceptance rate (typically between 0 and 1).
#' @param low_thres Numeric lower threshold (default: 0.1).
#' @param high_thres Numeric upper threshold (default: 0.9).
#' @param penalty Numeric penalty to apply if out of bounds (default: 0.05).
#'
#' @return Numeric: 0 if within bounds, `penalty` otherwise.
#'
#' @source Found in:
#'   - Scripts/single_algorithm/iCB_inspect_HPC_tune_results.R (line ~424)
penalize_acceptance <- function(acc_rate, low_thres = 0.1, high_thres = 0.9, 
                                 penalty = 0.05) {
  if (acc_rate < low_thres || acc_rate > high_thres) {
    return(penalty)
  } else {
    return(0)
  }
}


################################################################################
# HPC FOLDER FUNCTIONS: LRAcluster_HPC/full_LRAcluster_source.R
# These functions implement the LRAcluster algorithm for multi-omics integration
################################################################################

#' Check if Element is a Matrix
#'
#' @description
#' Simple utility to check if an element is NOT a matrix. Used in input validation
#' for LRAcluster data processing.
#'
#' @param x Any R object to check.
#'
#' @return Logical: TRUE if x is NOT a matrix, FALSE if x is a matrix.
#'
#' @source Found in:
#'   - Scripts/single_algorithm/LRAcluster_HPC/full_LRAcluster_source.R (line 1)
check.matrix.element <- function(x) {
  if (!is.matrix(x)) {
    return(T)
  } else {
    return(F)
  }
}


#' Get Number of Columns of Matrix Element
#'
#' @description
#' Wrapper around ncol() for use with sapply on list of matrices.
#'
#' @param x A matrix.
#'
#' @return Integer number of columns.
#'
#' @source Found in:
#'   - Scripts/single_algorithm/LRAcluster_HPC/full_LRAcluster_source.R (line 13)
ncol.element <- function(x) {
  return(ncol(x))
}


#' Get Number of Rows of Matrix Element
#'
#' @description
#' Wrapper around nrow() for use with sapply on list of matrices.
#'
#' @param x A matrix.
#'
#' @return Integer number of rows.
#'
#' @source Found in:
#'   - Scripts/single_algorithm/LRAcluster_HPC/full_LRAcluster_source.R (line 18)
nrow.element <- function(x) {
  return(nrow(x))
}


#' Check Data Matrix Based on Type
#'
#' @description
#' Dispatches to type-specific checking functions for binary, gaussian, or poisson
#' data matrices in LRAcluster.
#'
#' @param mat Data matrix to check.
#' @param type Character string: "binary", "gaussian", or "poisson".
#' @param name Character name for error messages.
#'
#' @return The validated/cleaned data matrix.
#'
#' @source Found in:
#'   - Scripts/single_algorithm/LRAcluster_HPC/full_LRAcluster_source.R (line 23)
check <- function(mat, type, name) {
  if (type == "binary") {
    return(check.binary(mat, name))
  } else if (type == "gaussian") {
    return(check.gaussian(mat, name))
  } else if (type == "poisson") {
    return(check.poisson(mat, name))
  } else {
    e <- paste("unknown type ", type, sep = "")
    stop(e)
  }
}


#' LRAcluster Main Function
#'
#' @description
#' Low-Rank Approximation clustering for integrating multiple data types. 
#' Performs dimension reduction while accounting for different data distributions
#' (binary, gaussian, poisson).
#'
#' @param data List of data matrices, each with features in rows and samples in columns.
#' @param types Character vector specifying the data type for each matrix 
#'   ("binary", "gaussian", or "poisson").
#' @param dimension Integer number of dimensions for the low-rank approximation (default: 2).
#' @param names Character vector of names for each data type (default: 1:length(data)).
#'
#' @return List containing:
#'   \item{coordinate}{Matrix of sample coordinates in the reduced space}
#'   \item{potential}{Ratio of explained vs unexplained likelihood}
#'
#' @details
#' The algorithm iteratively updates the low-rank approximation using type-specific
#' update rules until convergence. Handles missing values (NA) appropriately.
#'
#' @source Found in:
#'   - Scripts/single_algorithm/LRAcluster_HPC/full_LRAcluster_source.R (line 44)
LRAcluster <- function(data, types, dimension = 2, names = as.character(1:length(data))) {
  eps <- 0.0
  if (!is.list(data)) {
    stop("the input data must be a list!")
  }
  c <- sapply(data, check.matrix.element)
  if (sum(c) > 0) {
    stop("each element of input list must be a matrix!")
  }
  c <- sapply(data, ncol.element)
  if (length(levels(factor(c))) > 1) {
    stop("each element of input list must have the same column number!")
  }
  if (length(data) != length(types)) {
    stop("data and types must be the same length!")
  }
  nSample <- c[1]
  loglmin <- 0
  loglmax <- 0
  loglu <- 0.0
  nData <- length(data)
  for (i in 1:nData) {
    data[[i]] <- check(data[[i]], types[[i]], names[[i]])
  }
  nGeneArr <- sapply(data, nrow.element)
  nGene <- sum(nGeneArr)
  indexData <- list()
  k <- 1
  for (i in 1:nData) {
    indexData[[i]] <- (k):(k + nGeneArr[i] - 1)
    k <- k + nGeneArr[i]
  }
  base <- matrix(0, nGene, nSample)
  now <- matrix(0, nGene, nSample)
  update <- matrix(0, nGene, nSample)
  thr <- array(0, nData)
  for (i in 1:nData) {
    if (types[[i]] == "binary") {
      base[indexData[[i]], ] <- base.binary(data[[i]])
      loglmin <- loglmin + LLmin.binary(data[[i]], base[indexData[[i]], ])
      loglmax <- loglmax + LLmax.binary(data[[i]])
    } else if (types[[i]] == "gaussian") {
      base[indexData[[i]], ] <- base.gaussian(data[[i]])
      loglmin <- loglmin + LLmin.gaussian(data[[i]], base[indexData[[i]], ])
      loglmax <- loglmax + LLmax.gaussian(data[[i]])
    } else if (types[[i]] == "poisson") {
      base[indexData[[i]], ] <- base.poisson(data[[i]])
      loglmin <- loglmin + LLmin.poisson(data[[i]], base[indexData[[i]], ])
      loglmax <- loglmax + LLmax.poisson(data[[i]])
    }
  }
  for (i in 1:nData) {
    if (types[[i]] == "binary") {
      update[indexData[[i]], ] <- update.binary(data[[i]], base[indexData[[i]], ], now[indexData[[i]], ], exp(eps))
    } else if (types[[i]] == "gaussian") {
      update[indexData[[i]], ] <- update.gaussian(data[[i]], base[indexData[[i]], ], now[indexData[[i]], ], exp(eps))
    } else if (types[[i]] == "poisson") {
      update[indexData[[i]], ] <- update.poisson(data[[i]], base[indexData[[i]], ], now[indexData[[i]], ], exp(eps))
    }
  }
  update <- nuclear_approximation(update, dimension)
  nIter <- 0
  thres <- array(Inf, 3)
  epsN <- array(Inf, 2)
  while (T) {
    for (i in 1:nData) {
      if (types[[i]] == "binary") {
        thr[i] <- stop.binary(data[[i]], base[indexData[[i]], ], now[indexData[[i]], ], update[indexData[[i]], ])
      } else if (types[[i]] == "gaussian") {
        thr[i] <- stop.gaussian(data[[i]], base[indexData[[i]], ], now[indexData[[i]], ], update[indexData[[i]], ])
      } else if (types[[i]] == "poisson") {
        thr[i] <- stop.poisson(data[[i]], base[indexData[[i]], ], now[indexData[[i]], ], update[indexData[[i]], ])
      }
    }
    nIter <- nIter + 1
    thres[1] <- thres[2]
    thres[2] <- thres[3]
    thres[3] <- sum(thr)
    epsN[1] <- epsN[2]
    epsN[2] <- eps
    if (nIter > 5) {
      if (runif(1) < thres[1] * thres[3] / (thres[2] * thres[2] + thres[1] * thres[3])) {
        eps <- epsN[1] + 0.05 * runif(1) - 0.025
      } else {
        eps <- epsN[2] + 0.05 * runif(1) - 0.025
      }
      if (eps < -0.7) {
        eps <- 0
        epsN <- c(0, 0)
      }
      if (eps > 1.4) {
        eps <- 0
        epsN <- c(0, 0)
      }
    }
    if (sum(thr) < nData * 0.2) {
      break
    }
    now <- update
    for (i in 1:nData) {
      if (types[[i]] == "binary") {
        update[indexData[[i]], ] <- update.binary(data[[i]], base[indexData[[i]], ], now[indexData[[i]], ], exp(eps))
      } else if (types[[i]] == "gaussian") {
        update[indexData[[i]], ] <- update.gaussian(data[[i]], base[indexData[[i]], ], now[indexData[[i]], ], exp(eps))
      } else if (types[[i]] == "poisson") {
        update[indexData[[i]], ] <- update.poisson(data[[i]], base[indexData[[i]], ], now[indexData[[i]], ], exp(eps))
      }
    }
    update <- nuclear_approximation(update, dimension)
  }
  for (i in 1:nData) {
    if (types[[i]] == "binary") {
      loglu <- loglu + LL.binary(data[[i]], base[indexData[[i]], ], update[indexData[[i]], ])
    } else if (types[[i]] == "gaussian") {
      loglu <- loglu + LL.gaussian(data[[i]], base[indexData[[i]], ], update[indexData[[i]], ])
    } else if (types[[i]] == "poisson") {
      loglu <- loglu + LL.poisson(data[[i]], base[indexData[[i]], ], update[indexData[[i]], ])
    }
  }
  sv <- svd(update, nu = 0, nv = dimension)
  coordinate <- diag(c(sv$d[1:dimension], 0))[1:dimension, 1:dimension] %*% t(sv$v)
  colnames(coordinate) <- colnames(data[[1]])
  rownames(coordinate) <- paste("PC ", as.character(1:dimension), sep = "")
  ratio <- (loglu - loglmin) / (loglmax - loglmin)
  return(list("coordinate" = coordinate, "potential" = ratio))
}


#' Check Gaussian Data Row for Validity
#'
#' @description
#' Checks if a row of gaussian data contains at least one non-NA value.
#'
#' @param arr Numeric vector (one row of data).
#'
#' @return Logical: TRUE if valid, FALSE if all NA.
#'
#' @source Found in:
#'   - Scripts/single_algorithm/LRAcluster_HPC/full_LRAcluster_source.R (line 219)
check.gaussian.row <- function(arr) {
  if (sum(!is.na(arr)) == 0) {
    return(F)
  } else {
    return(T)
  }
}


#' Check and Clean Gaussian Data Matrix
#'
#' @description
#' Validates a gaussian data matrix, removing rows that are entirely NA.
#'
#' @param mat Numeric data matrix.
#' @param name Character name for warning messages.
#'
#' @return Cleaned matrix with invalid rows removed.
#'
#' @source Found in:
#'   - Scripts/single_algorithm/LRAcluster_HPC/full_LRAcluster_source.R (line 230)
check.gaussian <- function(mat, name) {
  index <- array(T, nrow(mat))
  for (i in 1:nrow(mat)) {
    if (sum(is.na(mat[i, ]) == ncol(mat))) {
      war <- paste("Warning: ", name, "'s ", as.character(i), " line is all NA. Delete this line", sep = "")
      warning(war)
      index[i] <- F
    }
  }
  mat_c <- mat[index, ]
  rownames(mat_c) <- rownames(mat)[index]
  colnames(mat_c) <- colnames(mat)
  return(mat_c)
}


#' Compute Base Value for Gaussian Row
#'
#' @description
#' Computes the mean of non-NA values in a row for gaussian data baseline.
#'
#' @param arr Numeric vector.
#'
#' @return Numeric mean of non-NA values.
#'
#' @source Found in:
#'   - Scripts/single_algorithm/LRAcluster_HPC/full_LRAcluster_source.R (line 248)
base.gaussian.row <- function(arr) {
  idx <- !is.na(arr)
  return(mean(arr[idx]))
}


#' Compute Baseline Matrix for Gaussian Data
#'
#' @description
#' Computes row-wise means as the baseline for gaussian data in LRAcluster.
#'
#' @param mat Numeric data matrix.
#'
#' @return Matrix of same dimensions with row means repeated across columns.
#'
#' @source Found in:
#'   - Scripts/single_algorithm/LRAcluster_HPC/full_LRAcluster_source.R (line 254)
base.gaussian <- function(mat) {
  mat_b <- matrix(0, nrow(mat), ncol(mat))
  ar_b <- apply(mat, 1, base.gaussian.row)
  mat_b[1:nrow(mat_b), ] <- ar_b
  return(mat_b)
}


#' Update Step for Gaussian Data
#'
#' @description
#' Performs the gradient update step for gaussian data type in LRAcluster.
#'
#' @param mat Original data matrix.
#' @param mat_b Baseline matrix.
#' @param mat_now Current approximation matrix.
#' @param eps Step size parameter.
#'
#' @return Updated matrix.
#'
#' @source Found in:
#'   - Scripts/single_algorithm/LRAcluster_HPC/full_LRAcluster_source.R (line 262)
update.gaussian <- function(mat, mat_b, mat_now, eps) {
  mat_p <- mat_b + mat_now
  mat_u <- matrix(0, nrow(mat), ncol(mat))
  index <- !is.na(mat)
  mat_u[index] <- mat_now[index] + eps * epsilon.gaussian * (mat[index] - mat_p[index])
  index <- is.na(mat)
  mat_u[index] <- mat_now[index]
  return(mat_u)
}


#' Stopping Criterion for Gaussian Data
#'
#' @description
#' Computes the log-likelihood improvement for convergence checking.
#'
#' @param mat Original data matrix.
#' @param mat_b Baseline matrix.
#' @param mat_now Current approximation.
#' @param mat_u Updated approximation.
#'
#' @return Numeric log-likelihood difference.
#'
#' @source Found in:
#'   - Scripts/single_algorithm/LRAcluster_HPC/full_LRAcluster_source.R (line 273)
stop.gaussian <- function(mat, mat_b, mat_now, mat_u) {
  index <- !is.na(mat)
  mn <- mat_b + mat_now
  mu <- mat_b + mat_u
  ren <- mat[index] - mn[index]
  reu <- mat[index] - mu[index]
  lgn <- -0.5 * sum(ren * ren)
  lgu <- -0.5 * sum(reu * reu)
  return(lgu - lgn)
}


#' Log-Likelihood for Gaussian Data
#'
#' @description
#' Computes the log-likelihood of the current approximation for gaussian data.
#'
#' @param mat Original data matrix.
#' @param mat_b Baseline matrix.
#' @param mat_u Current approximation matrix.
#'
#' @return Numeric log-likelihood value.
#'
#' @source Found in:
#'   - Scripts/single_algorithm/LRAcluster_HPC/full_LRAcluster_source.R (line 285)
LL.gaussian <- function(mat, mat_b, mat_u) {
  index <- !is.na(mat)
  mu <- mat_b + mat_u
  reu <- mat[index] - mu[index]
  lgu <- -0.5 * sum(reu * reu)
  return(lgu)
}


#' Maximum Log-Likelihood for Gaussian Data
#'
#' @description
#' Returns the maximum achievable log-likelihood for gaussian data (saturated model).
#'
#' @param mat Original data matrix.
#'
#' @return Numeric: 0 (perfect fit has zero residual).
#'
#' @source Found in:
#'   - Scripts/single_algorithm/LRAcluster_HPC/full_LRAcluster_source.R (line 294)
LLmax.gaussian <- function(mat) {
  return(0.0)
}


#' Minimum Log-Likelihood for Gaussian Data
#'
#' @description
#' Computes the minimum (null model) log-likelihood for gaussian data.
#'
#' @param mat Original data matrix.
#' @param mat_b Baseline matrix.
#'
#' @return Numeric log-likelihood of the null model.
#'
#' @source Found in:
#'   - Scripts/single_algorithm/LRAcluster_HPC/full_LRAcluster_source.R (line 299)
LLmin.gaussian <- function(mat, mat_b) {
  index <- !is.na(mat)
  reu <- mat[index] - mat_b[index]
  lgu <- -0.5 * sum(reu * reu)
  return(lgu)
}


#' Gaussian Type Base Function
#'
#' @description
#' Standalone function to fit gaussian data using LRAcluster methodology.
#'
#' @param data Numeric data matrix.
#' @param dimension Target dimension (default: 2).
#' @param name Character name for messages.
#'
#' @return Low-rank approximation matrix.
#'
#' @source Found in:
#'   - Scripts/single_algorithm/LRAcluster_HPC/full_LRAcluster_source.R (line 307)
gaussian_base <- function(data, dimension = 2, name = "test") {
  data <- check.gaussian(data, name)
  data_b <- base.gaussian(data)
  data_now <- matrix(0, nrow(data), ncol(data))
  data_u <- update.gaussian(data, data_b, data_now)
  data_u <- nuclear_approximation(data_u, dimension)
  while (T) {
    thr <- stop.gaussian(data, data_b, data_now, data_u)
    message(thr)
    if (thr < 0.2) {
      break
    }
    data_now <- data_u
    data_u <- update.gaussian(data, data_b, data_now)
    data_u <- nuclear_approximation(data_u, dimension)
  }
  return(data_now)
}


#' Check Poisson Data Row for Validity
#'
#' @description
#' Validates a row for poisson data (non-NA and non-negative values).
#'
#' @param arr Numeric vector.
#'
#' @return Logical: TRUE if valid, FALSE otherwise.
#'
#' @source Found in:
#'   - Scripts/single_algorithm/LRAcluster_HPC/full_LRAcluster_source.R (line 333)
check.poisson.row <- function(arr) {
  if (sum(!is.na(arr)) == 0) {
    return(F)
  } else {
    idx <- !is.na(arr)
    if (sum(arr[idx] < 0) > 0) {
      return(F)
    } else {
      return(T)
    }
  }
}


#' Check and Clean Poisson Data Matrix
#'
#' @description
#' Validates poisson data matrix, adds 1 to all counts (for log stability),
#' and removes invalid rows.
#'
#' @param mat Count data matrix.
#' @param name Character name for messages.
#'
#' @return Cleaned matrix with +1 added to counts.
#'
#' @source Found in:
#'   - Scripts/single_algorithm/LRAcluster_HPC/full_LRAcluster_source.R (line 353)
check.poisson <- function(mat, name) {
  w <- paste(name, " is poisson type. Add 1 to all counts", sep = "")
  message(w)
  index <- apply(mat, 1, check.poisson.row)
  n <- sum(!index)
  if (n > 0) {
    w <- paste("Warning: ", name, " have ", as.character(n), " invalid lines!", sep = "")
    warning(w)
  }
  mat_c <- mat[index, ] + 1
  rownames(mat_c) <- rownames(mat)[index]
  colnames(mat_c) <- colnames(mat)
  return(mat_c)
}


#' Compute Base Value for Poisson Row
#'
#' @description
#' Computes log-mean of non-NA values for poisson baseline.
#'
#' @param arr Numeric vector.
#'
#' @return Numeric mean of log values.
#'
#' @source Found in:
#'   - Scripts/single_algorithm/LRAcluster_HPC/full_LRAcluster_source.R (line 370)
base.poisson.row <- function(arr) {
  idx <- !is.na(arr)
  m <- sum(log(arr[idx]))
  n <- sum(idx)
  return(m / n)
}


#' Compute Baseline Matrix for Poisson Data
#'
#' @description
#' Computes row-wise log-means as baseline for poisson data.
#'
#' @param mat Count data matrix.
#'
#' @return Baseline matrix with log-means.
#'
#' @source Found in:
#'   - Scripts/single_algorithm/LRAcluster_HPC/full_LRAcluster_source.R (line 378)
base.poisson <- function(mat) {
  mat_b <- matrix(0, nrow(mat), ncol(mat))
  ar_b <- apply(mat, 1, base.poisson.row)
  mat_b[1:nrow(mat_b), ] <- ar_b
  return(mat_b)
}


#' Update Step for Poisson Data
#'
#' @description
#' Performs gradient update for poisson data in LRAcluster.
#'
#' @param mat Original count matrix.
#' @param mat_b Baseline matrix.
#' @param mat_now Current approximation.
#' @param eps Step size parameter.
#'
#' @return Updated matrix.
#'
#' @source Found in:
#'   - Scripts/single_algorithm/LRAcluster_HPC/full_LRAcluster_source.R (line 386)
update.poisson <- function(mat, mat_b, mat_now, eps) {
  mat_p <- mat_b + mat_now
  mat_u <- matrix(0, nrow(mat), ncol(mat))
  index <- !is.na(mat)
  mat_u[index] <- mat_now[index] + eps * epsilon.poisson * (log(mat[index]) - mat_p[index])
  index <- is.na(mat)
  mat_u[index] <- mat_now[index]
  return(mat_u)
}


#' Stopping Criterion for Poisson Data
#'
#' @description
#' Computes log-likelihood improvement for poisson data convergence check.
#'
#' @param mat Original count matrix.
#' @param mat_b Baseline matrix.
#' @param mat_now Current approximation.
#' @param mat_u Updated approximation.
#'
#' @return Numeric log-likelihood difference.
#'
#' @source Found in:
#'   - Scripts/single_algorithm/LRAcluster_HPC/full_LRAcluster_source.R (line 397)
stop.poisson <- function(mat, mat_b, mat_now, mat_u) {
  index <- !is.na(mat)
  mn <- mat_b + mat_now
  mu <- mat_b + mat_u
  lgn <- sum(mat[index] * mn[index] - exp(mn[index]))
  lgu <- sum(mat[index] * mu[index] - exp(mu[index]))
  return(lgu - lgn)
}


#' Log-Likelihood for Poisson Data
#'
#' @description
#' Computes poisson log-likelihood for current approximation.
#'
#' @param mat Original count matrix.
#' @param mat_b Baseline matrix.
#' @param mat_u Current approximation.
#'
#' @return Numeric log-likelihood.
#'
#' @source Found in:
#'   - Scripts/single_algorithm/LRAcluster_HPC/full_LRAcluster_source.R (line 407)
LL.poisson <- function(mat, mat_b, mat_u) {
  index <- !is.na(mat)
  mu <- mat_b + mat_u
  lgu <- sum(mat[index] * mu[index] - exp(mu[index]))
  return(lgu)
}


#' Maximum Log-Likelihood for Poisson Data
#'
#' @description
#' Computes maximum (saturated model) log-likelihood for poisson data.
#'
#' @param mat Original count matrix.
#'
#' @return Numeric maximum log-likelihood.
#'
#' @source Found in:
#'   - Scripts/single_algorithm/LRAcluster_HPC/full_LRAcluster_source.R (line 415)
LLmax.poisson <- function(mat) {
  index <- !is.na(mat)
  lgu <- sum(mat[index] * log(mat[index]) - mat[index])
  return(lgu)
}


#' Minimum Log-Likelihood for Poisson Data
#'
#' @description
#' Computes minimum (null model) log-likelihood for poisson data.
#'
#' @param mat Original count matrix.
#' @param mat_b Baseline matrix.
#'
#' @return Numeric minimum log-likelihood.
#'
#' @source Found in:
#'   - Scripts/single_algorithm/LRAcluster_HPC/full_LRAcluster_source.R (line 422)
LLmin.poisson <- function(mat, mat_b) {
  index <- !is.na(mat)
  lgu <- sum(mat[index] * mat_b[index] - exp(mat_b[index]))
  return(lgu)
}


#' Poisson Type Base Function
#'
#' @description
#' Standalone function to fit poisson data using LRAcluster methodology.
#'
#' @param data Count data matrix.
#' @param dimension Target dimension (default: 2).
#' @param name Character name for messages.
#'
#' @return Low-rank approximation matrix.
#'
#' @source Found in:
#'   - Scripts/single_algorithm/LRAcluster_HPC/full_LRAcluster_source.R (line 429)
poisson_type_base <- function(data, dimension = 2, name = "test") {
  data <- check.poisson(data, name)
  data_b <- base.poisson(data)
  data_now <- matrix(0, nrow(data), ncol(data))
  data_u <- update.poisson(data, data_b, data_now)
  data_u <- nuclear_approximation(data_u, dimension)
  while (T) {
    thr <- stop.poisson(data, data_b, data_now, data_u)
    message(thr)
    if (thr < 0.2) {
      break
    }
    data_now <- data_u
    data_u <- update.poisson(data, data_b, data_now)
    data_u <- nuclear_approximation(data_u, dimension)
  }
  return(data_now)
}


#' Check Binary Data Row for Validity
#'
#' @description
#' Validates a binary data row (must have variation - not all 0 or all 1).
#'
#' @param arr Binary vector.
#'
#' @return Logical: TRUE if valid, FALSE if constant or all NA.
#'
#' @source Found in:
#'   - Scripts/single_algorithm/LRAcluster_HPC/full_LRAcluster_source.R (line 454)
check.binary.row <- function(arr) {
  if (sum(!is.na(arr)) == 0) {
    return(F)
  } else {
    idx <- !is.na(arr)
    if (sum(arr[idx]) == 0 || sum(arr[idx]) == sum(idx)) {
      return(F)
    } else {
      return(T)
    }
  }
}


#' Check and Clean Binary Data Matrix
#'
#' @description
#' Validates binary data matrix and removes rows without variation.
#'
#' @param mat Binary data matrix.
#' @param name Character name for messages.
#'
#' @return Cleaned matrix with constant rows removed.
#'
#' @source Found in:
#'   - Scripts/single_algorithm/LRAcluster_HPC/full_LRAcluster_source.R (line 473)
check.binary <- function(mat, name) {
  index <- apply(mat, 1, check.binary.row)
  n <- sum(!index)
  if (n > 0) {
    w <- paste("Warning: ", name, " have ", as.character(n), " invalid lines!", sep = "")
    warning(w)
  }
  mat_c <- mat[index, ]
  rownames(mat_c) <- rownames(mat)[index]
  colnames(mat_c) <- colnames(mat)
  return(mat_c)
}


#' Compute Base Value for Binary Row
#'
#' @description
#' Computes log-odds baseline for a binary row.
#'
#' @param arr Binary vector.
#'
#' @return Numeric log-odds value.
#'
#' @source Found in:
#'   - Scripts/single_algorithm/LRAcluster_HPC/full_LRAcluster_source.R (line 488)
base.binary.row <- function(arr) {
  idx <- !is.na(arr)
  n <- sum(idx)
  m <- sum(arr[idx])
  return(log(m / (n - m)))
}


#' Compute Baseline Matrix for Binary Data
#'
#' @description
#' Computes row-wise log-odds as baseline for binary data.
#'
#' @param mat Binary data matrix.
#'
#' @return Baseline matrix with log-odds.
#'
#' @source Found in:
#'   - Scripts/single_algorithm/LRAcluster_HPC/full_LRAcluster_source.R (line 496)
base.binary <- function(mat) {
  mat_b <- matrix(0, nrow(mat), ncol(mat))
  ar_b <- apply(mat, 1, base.binary.row)
  mat_b[1:nrow(mat_b), ] <- ar_b
  return(mat_b)
}


#' Update Step for Binary Data
#'
#' @description
#' Performs logistic regression gradient update for binary data.
#'
#' @param mat Original binary matrix.
#' @param mat_b Baseline matrix.
#' @param mat_now Current approximation.
#' @param eps Step size parameter.
#'
#' @return Updated matrix.
#'
#' @source Found in:
#'   - Scripts/single_algorithm/LRAcluster_HPC/full_LRAcluster_source.R (line 504)
update.binary <- function(mat, mat_b, mat_now, eps) {
  mat_p <- mat_b + mat_now
  mat_u <- matrix(0, nrow(mat), ncol(mat))
  idx1 <- !is.na(mat) & mat == 1
  idx0 <- !is.na(mat) & mat == 0
  index <- is.na(mat)
  arr <- exp(mat_p)
  mat_u[index] <- mat_now[index]
  mat_u[idx1] <- mat_now[idx1] + eps * epsilon.binary / (1.0 + arr[idx1])
  mat_u[idx0] <- mat_now[idx0] - eps * epsilon.binary * arr[idx0] / (1.0 + arr[idx0])
  return(mat_u)
}


#' Stopping Criterion for Binary Data
#'
#' @description
#' Computes log-likelihood improvement for binary data convergence check.
#'
#' @param mat Original binary matrix.
#' @param mat_b Baseline matrix.
#' @param mat_now Current approximation.
#' @param mat_u Updated approximation.
#'
#' @return Numeric log-likelihood difference.
#'
#' @source Found in:
#'   - Scripts/single_algorithm/LRAcluster_HPC/full_LRAcluster_source.R (line 518)
stop.binary <- function(mat, mat_b, mat_now, mat_u) {
  index <- !is.na(mat)
  mn <- mat_b + mat_now
  mu <- mat_b + mat_u
  arn <- exp(mn)
  aru <- exp(mu)
  idx1 <- !is.na(mat) & mat == 1
  idx0 <- !is.na(mat) & mat == 0
  lgn <- sum(log(arn[idx1] / (1 + arn[idx1]))) + sum(log(1 / (1 + arn[idx0])))
  lgu <- sum(log(aru[idx1] / (1 + aru[idx1]))) + sum(log(1 / (1 + aru[idx0])))
  return(lgu - lgn)
}


#' Log-Likelihood for Binary Data
#'
#' @description
#' Computes binary (logistic) log-likelihood for current approximation.
#'
#' @param mat Original binary matrix.
#' @param mat_b Baseline matrix.
#' @param mat_u Current approximation.
#'
#' @return Numeric log-likelihood.
#'
#' @source Found in:
#'   - Scripts/single_algorithm/LRAcluster_HPC/full_LRAcluster_source.R (line 532)
LL.binary <- function(mat, mat_b, mat_u) {
  index <- !is.na(mat)
  mu <- mat_b + mat_u
  aru <- exp(mu)
  idx1 <- !is.na(mat) & mat == 1
  idx0 <- !is.na(mat) & mat == 0
  lgu <- sum(log(aru[idx1] / (1 + aru[idx1]))) + sum(log(1 / (1 + aru[idx0])))
  return(lgu)
}


#' Maximum Log-Likelihood for Binary Data
#'
#' @description
#' Returns maximum (saturated model) log-likelihood for binary data.
#'
#' @param mat Original binary matrix.
#'
#' @return Numeric: 0 (perfect prediction has zero loss).
#'
#' @source Found in:
#'   - Scripts/single_algorithm/LRAcluster_HPC/full_LRAcluster_source.R (line 543)
LLmax.binary <- function(mat) {
  return(0)
}


#' Minimum Log-Likelihood for Binary Data
#'
#' @description
#' Computes minimum (null model) log-likelihood for binary data.
#'
#' @param mat Original binary matrix.
#' @param mat_b Baseline matrix.
#'
#' @return Numeric minimum log-likelihood.
#'
#' @source Found in:
#'   - Scripts/single_algorithm/LRAcluster_HPC/full_LRAcluster_source.R (line 548)
LLmin.binary <- function(mat, mat_b) {
  index <- !is.na(mat)
  aru <- exp(mat_b)
  idx1 <- !is.na(mat) & mat == 1
  idx0 <- !is.na(mat) & mat == 0
  lgu <- sum(log(aru[idx1] / (1 + aru[idx1]))) + sum(log(1 / (1 + aru[idx0])))
  return(lgu)
}


#' Binary Type Base Function
#'
#' @description
#' Standalone function to fit binary data using LRAcluster methodology.
#'
#' @param data Binary data matrix.
#' @param dimension Target dimension (default: 2).
#' @param name Character name for messages.
#'
#' @return Low-rank approximation matrix.
#'
#' @source Found in:
#'   - Scripts/single_algorithm/LRAcluster_HPC/full_LRAcluster_source.R (line 558)
binary_type_base <- function(data, dimension = 2, name = "test") {
  data <- check.binary(data, name)
  data_b <- base.binary(data)
  data_now <- matrix(0, nrow(data), ncol(data))
  data_u <- update.binary(data, data_b, data_now)
  data_u <- nuclear_approximation(data_u, dimension)
  while (T) {
    thr <- stop.binary(data, data_b, data_now, data_u)
    message(thr)
    if (thr < 0.2) {
      break
    }
    data_now <- data_u
    data_u <- update.binary(data, data_b, data_now)
    data_u <- nuclear_approximation(data_u, dimension)
  }
  return(data_now)
}


#' Nuclear Norm Approximation
#'
#' @description
#' Computes the nuclear norm (low-rank) approximation of a matrix using SVD.
#' Projects the matrix onto the space of matrices with rank at most 'dimension'.
#'
#' @param mat Input matrix to approximate.
#' @param dimension Target rank for the approximation.
#'
#' @return Low-rank approximation of the input matrix.
#'
#' @details
#' Uses soft-thresholding on singular values: subtracts the (dimension+1)th
#' singular value from all singular values and truncates at zero.
#'
#' @source Found in:
#'   - Scripts/single_algorithm/LRAcluster_HPC/full_LRAcluster_source.R (line 580)
nuclear_approximation <- function(mat, dimension) {
  svd <- svd(mat, nu = 0, nv = 0)
  if (dimension < length(svd$d)) {
    lambda <- svd$d[dimension + 1]
    svd <- svd(mat, nu = dimension, nv = dimension)
    indexh <- svd$d > lambda
    indexm <- svd$d < lambda
    dia <- array(svd$d, length(svd$d))
    dia[indexh] <- dia[indexh] - lambda
    dia[indexm] <- 0
    mat_low <- svd$u %*% diag(c(dia[1:dimension], 0))[1:dimension, 1:dimension] %*% t(svd$v)
  } else {
    mat_low <- mat
  }
  return(mat_low)
}


################################################################################
# HPC FOLDER FUNCTIONS: RWR-F_HPC/RWR-F_source.R
# Random Walk with Restart based Fusion for multi-omics integration
################################################################################

#' Random Walk with Restart Fusion (Neighbor-based)
#'
#' @description
#' Implements Random Walk with Restart (RWR) based similarity network fusion
#' with neighbor-based transitions. Integrates multiple similarity networks
#' into a single fused similarity matrix.
#'
#' @param sim_list List of similarity matrices (one per omics type). All matrices
#'   must have the same dimensions and sample ordering.
#' @param iteration_max Maximum number of iterations for convergence (default: 1000).
#' @param gama Restart probability parameter (default: 0.7). Higher values favor
#'   returning to the starting node.
#' @param neighbor_num Number of nearest neighbors to consider (default: 10).
#' @param alpha Weight for self-network vs cross-network transitions (default: 0.9).
#' @param beta Self-transition probability within neighbor set (default: 0.9).
#'
#' @return A fused similarity matrix combining information from all input networks.
#'   Matrix is normalized to [0,1] range and made symmetric.
#'
#' @details
#' Uses parallel processing via foreach. The algorithm:
#' 1. Constructs neighbor-based adjacency matrices
#' 2. Builds a transition probability matrix across all networks
#' 3. Performs random walk with restart until convergence
#' 4. Aggregates results across all networks
#'
#' @note Requires the foreach and doParallel packages for parallel execution.
#'
#' @source Found in:
#'   - Scripts/single_algorithm/RWR-F_HPC/RWR-F_source.R (line 1)
RWR_fusion_neighbor <- function(sim_list, iteration_max = 1000, gama = 0.7,
                                 neighbor_num = 10, alpha = 0.9, beta = 0.9) {
  len_net <- length(sim_list)
  len_node <- nrow(sim_list[[1]])
  
  # Parameters
  lamda <- 1 / len_net
  
  # Output matrix
  result_matrix <- matrix(data = NA, nrow = len_net * len_node, ncol = len_net * len_node)
  
  # Adjacency matrix construction
  adj_list <- NULL
  for (i in 1:len_net) {
    temp <- sim_list[[i]]
    temp[which(temp != 0)] <- 0
    for (j in 1:len_node) {
      num <- setdiff(which(rank(-sim_list[[i]][j, ]) <= neighbor_num + 1), j)
      temp[j, num] <- (1 - beta) / neighbor_num
      temp[j, j] <- beta
      rm(num)
    }
    adj_list[[i]] <- temp
    rm(temp)
  }
  
  # Random walk (parallel)
  start_time_all <- Sys.time()
  result_matrix <- foreach(i = 1:nrow(result_matrix), .combine = rbind) %dopar% {
    index <- ceiling(i / len_node)
    
    # Set initial value
    p0 <- rep(0, len_net * len_node)
    p0[i] <- alpha
    temp <- i %% len_node
    for (j in 1:len_net) {
      if (j != index) {
        if (temp == 0) {
          p0[(1 + (j - 1) * len_node):(j * len_node)] <- (1 - alpha) * (adj_list[[index]][len_node, ]) / (len_net - 1)
        } else {
          p0[(1 + (j - 1) * len_node):(j * len_node)] <- (1 - alpha) * (adj_list[[index]][temp, ]) / (len_net - 1)
        }
      }
    }
    
    # Build transition probability matrix
    W_list <- NULL
    for (m in 1:len_net) {
      for (n in 1:len_net) {
        if (m == n) {
          temp <- NULL
          for (j in 1:len_node) {
            temp <- rbind(temp, lamda * sim_list[[m]][j, ] / sum(sim_list[[m]][j, ]))
          }
          W_list <- c(W_list, list(temp))
          rm(temp)
        } else {
          temp <- NULL
          for (j in 1:len_node) {
            temp <- rbind(temp, lamda * adj_list[[m]][j, ] / sum(adj_list[[m]][j, ]))
          }
          W_list <- c(W_list, list(temp))
          rm(temp)
        }
      }
    }
    W <- NULL
    k <- 1
    for (m in 1:len_net) {
      temp <- NULL
      for (n in 1:len_net) {
        temp <- cbind(temp, W_list[[k]])
        k <- k + 1
      }
      W <- rbind(W, temp)
    }
    
    # Iterate
    temp <- p0
    for (iteration_num in 1:iteration_max) {
      pt0 <- temp
      pt1 <- (1 - gama) * t(W) %*% pt0 + gama * p0
      temp <- pt1
      end_clock <- sum(abs(pt1 - pt0))
      if (end_clock <= 1e-10) {
        break
      }
    }
    rm(temp)
    rm(end_clock)
    pt1 <- pt1 / sum(pt1)
    result_matrix[i, ] <- pt1
    t(pt1)
  }
  end_time_all <- Sys.time()
  print(end_time_all - start_time_all)
  
  rownames(result_matrix) <- rep(rownames(sim_list[[1]]), len_net)
  colnames(result_matrix) <- rownames(result_matrix)
  
  RWR_similarity <- result_matrix[1:len_node, 1:len_node]
  RWR_similarity[which(RWR_similarity != 0)] <- 0
  for (m in 1:len_net) {
    for (n in 1:len_net) {
      RWR_similarity <- result_matrix[(1 + (m - 1) * len_node):(m * len_node),
                                       (1 + (m - 1) * len_node):(m * len_node)] + RWR_similarity
    }
  }
  RWR_similarity <- RWR_similarity / (len_net)
  RWR_similarity <- (RWR_similarity + t(RWR_similarity)) / 2
  RWR_similarity <- RWR_similarity / max(RWR_similarity)
  
  return(RWR_similarity)
}


#' Random Walk with Restart Fusion (Standard)
#'
#' @description
#' Implements standard Random Walk with Restart based similarity network fusion.
#' A simpler variant without neighbor-based transitions.
#'
#' @param sim_list List of similarity matrices (one per omics type).
#' @param iteration_max Maximum iterations for convergence (default: 1000).
#' @param gama Restart probability parameter (default: 0.7).
#'
#' @return A fused similarity matrix. Normalized to [0,1] and symmetric.
#'
#' @details
#' Simpler than RWR_fusion_neighbor - uses identity matrices for cross-network
#' transitions instead of neighbor-based adjacency matrices.
#'
#' @note Requires foreach and doParallel packages for parallel execution.
#'
#' @source Found in:
#'   - Scripts/single_algorithm/RWR-F_HPC/RWR-F_source.R (line 176)
RWR_fusion <- function(sim_list, iteration_max = 1000, gama = 0.7) {
  len_net <- length(sim_list)
  len_node <- nrow(sim_list[[1]])
  
  # Parameters
  alpha_list <- rep(1 / len_net, len_net)
  lamda <- 1 / len_net
  
  # Output matrix
  result_matrix <- matrix(data = NA, nrow = len_net * len_node, ncol = len_net * len_node)
  
  # Random walk (parallel)
  start_time_all <- Sys.time()
  result_matrix <- foreach(i = 1:nrow(result_matrix), .combine = rbind) %dopar% {
    start_time <- Sys.time()
    index <- ceiling(i / len_node)
    
    # Set initial value
    p0 <- rep(0, len_net * len_node)
    temp <- i %% len_node
    if (temp == 0) {
      temp <- len_node
    }
    for (j in 1:len_net) {
      p0[temp + (j - 1) * len_node] <- alpha_list[j]
    }
    
    # Build transition probability matrix
    W_list <- NULL
    for (m in 1:len_net) {
      for (n in 1:len_net) {
        if (m == n) {
          temp <- NULL
          for (j in 1:len_node) {
            temp <- rbind(temp, lamda * sim_list[[m]][j, ] / sum(sim_list[[m]][j, ]))
          }
          W_list <- c(W_list, list(temp))
          rm(temp)
        } else {
          W_list <- c(W_list, list(diag(x = lamda, nrow = len_node, ncol = len_node)))
        }
      }
    }
    W <- NULL
    k <- 1
    for (m in 1:len_net) {
      temp <- NULL
      for (n in 1:len_net) {
        temp <- cbind(temp, W_list[[k]])
        k <- k + 1
      }
      W <- rbind(W, temp)
    }
    
    # Iterate
    temp <- p0
    for (iteration_num in 1:iteration_max) {
      pt0 <- temp
      pt1 <- (1 - gama) * t(W) %*% pt0 + gama * p0
      temp <- pt1
      end_clock <- sum(abs(pt1 - pt0))
      if (end_clock <= 1e-10) {
        break
      }
    }
    rm(temp)
    rm(end_clock)
    pt1 <- pt1 / sum(pt1)
    result_matrix[i, ] <- pt1
    end_time <- Sys.time()
    t(pt1)
  }
  end_time_all <- Sys.time()
  print(end_time_all - start_time_all)
  
  rownames(result_matrix) <- rep(rownames(sim_list[[1]]), len_net)
  colnames(result_matrix) <- rownames(result_matrix)
  
  RWR_similarity <- result_matrix[1:len_node, 1:len_node]
  RWR_similarity[which(RWR_similarity != 0)] <- 0
  for (m in 1:len_net) {
    for (n in 1:len_net) {
      RWR_similarity <- result_matrix[(1 + (m - 1) * len_node):(m * len_node),
                                       (1 + (m - 1) * len_node):(m * len_node)] + RWR_similarity
    }
  }
  RWR_similarity <- RWR_similarity / (len_net)
  RWR_similarity <- (RWR_similarity + t(RWR_similarity)) / 2
  RWR_similarity <- RWR_similarity / max(RWR_similarity)
  
  return(RWR_similarity)
}


################################################################################
# SCRIPTS/CONSENSUS FUNCTIONS
# Functions for consensus clustering analysis and post-processing
################################################################################

#' Calculate ARI Between Two Pathway Sets
#'
#' @description
#' Computes Adjusted Rand Index (ARI) between two sets of pathways by converting
#' them to binary membership vectors over a common universe of pathways.
#'
#' @param set1 Character vector of pathway IDs in the first set.
#' @param set2 Character vector of pathway IDs in the second set.
#' @param universe Character vector containing all possible pathway IDs.
#'
#' @return Numeric ARI value between -1 and 1.
#'
#' @details
#' Converts each set to a binary vector indicating membership in the universe,
#' then applies mclust::adjustedRandIndex() to compare the two vectors.
#'
#' @source Found in:
#'   - Scripts/Consensus/Post_consensus.R (line 169)
#'   - Scripts/Consensus/Post_consensus_more_than_2_clusters.R (line 157)
calc_ari <- function(set1, set2, universe) {
  v1 <- as.integer(universe %in% set1)
  v2 <- as.integer(universe %in% set2)
  adjustedRandIndex(v1, v2)
}


################################################################################
# SCRIPTS/MOVICS FUNCTIONS
# Functions for MOVICS baseline analysis and output inspection
################################################################################

#' Check Missing Values with Row Indices
#'
#' @description
#' Examines each matrix in a list for missing values, returning both the count
#' of NA values and the row indices where they occur.
#'
#' @param input List of matrices to check for missing values.
#'
#' @return Named list (one element per input matrix), each containing:
#'   \item{na_count}{Total number of NA values in the matrix}
#'   \item{missing_indices}{Integer vector of row indices containing NA values}
#'
#' @details
#' Useful for data quality control before running multi-omics integration
#' algorithms that may not handle missing values gracefully.
#'
#' @source Found in:
#'   - Scripts/MOVICS/MOVICS_baseline.R (line 227)
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


#' Preprocess Text for Word Frequency Analysis
#'
#' @description
#' Converts a vector of text strings into a term-document matrix and extracts
#' term frequencies for wordcloud generation.
#'
#' @param textVector Character vector of text strings to process.
#' @param removeStopwords Logical; if TRUE, removes English stopwords (default: FALSE).
#'
#' @return Data frame with columns:
#'   \item{term}{Character terms extracted from text}
#'   \item{freq}{Integer frequency count for each term}
#'
#' @details
#' Uses the tm package to create a corpus, optionally remove stopwords, and
#' build a term-document matrix. Frequencies are sorted in descending order.
#'
#' @source Found in:
#'   - Scripts/MOVICS/MOVICS_consensus_output_inspection.R (line 580)
#'   - Scripts/method_comparisons/Comparisons.R (line 505)
preprocessText <- function(textVector, removeStopwords = FALSE) {
  # Create a text corpus
  corp <- Corpus(VectorSource(textVector))
  
  # Optionally remove English stopwords
  if (removeStopwords) {
    corp <- tm_map(corp, removeWords, stopwords("english"))
  }
  
  # Create a term-document matrix
  tdm <- TermDocumentMatrix(corp)
  
  # Convert the matrix to a data frame of terms and their frequencies
  m <- as.matrix(tdm)
  termFrequency <- sort(rowSums(m), decreasing = TRUE)
  dfFrequency <- data.frame(term = names(termFrequency), freq = termFrequency)
  
  return(dfFrequency)
}


#' Calculate Overlap Coefficient Between Two Sets
#'
#' @description
#' Computes the overlap coefficient (Szymkiewicz-Simpson coefficient) between
#' two sets, defined as the intersection size divided by the smaller set size.
#'
#' @param set1 First set (vector).
#' @param set2 Second set (vector).
#'
#' @return Numeric value between 0 and 1. Returns 0 if either set is empty.
#'
#' @details
#' The overlap coefficient is less sensitive to set size differences than
#' the Jaccard index, making it useful for comparing pathway sets of
#' different sizes.
#'
#' @source Found in:
#'   - Scripts/MOVICS/MOVICS_consensus_output_inspection.R (line 644)
#'   - Scripts/method_comparisons/Comparisons.R (line 383)
overlap_coefficient <- function(set1, set2) {
  # Handle potential empty sets:
  if (length(set1) == 0 || length(set2) == 0) return(0)
  length(intersect(set1, set2)) / min(length(set1), length(set2))
}


#' Calculate Similarity Matrix Using Overlap Coefficient
#'
#' @description
#' Computes a pairwise similarity matrix between algorithms based on their
#' pathway sets using the overlap coefficient.
#'
#' @param algorithms Character vector of algorithm names.
#' @param pathway_column Unquoted column name containing pathway lists.
#' @param pathway_list Data frame with Algorithm column and pathway list column.
#'
#' @return Square numeric matrix with overlap coefficients between all pairs
#'   of algorithms.
#'
#' @details
#' Uses dplyr programming with {{ }} for column selection. The matrix is
#' symmetric with 1s on the diagonal.
#'
#' @source Found in:
#'   - Scripts/MOVICS/MOVICS_consensus_output_inspection.R (line 649)
calculate_similarities <- function(algorithms, pathway_column, pathway_list) {
  # Create a matrix to store the results
  similarity_matrix <- matrix(0, nrow = length(algorithms), ncol = length(algorithms),
                              dimnames = list(algorithms, algorithms))
  
  for (i in 1:length(algorithms)) {
    for (j in i:length(algorithms)) {
      # Correctly filter and extract the pathway list for each algorithm
      path1 <- pathway_list %>% 
        filter(Algorithm == algorithms[i]) %>% 
        pull({{pathway_column}})
      path2 <- pathway_list %>% 
        filter(Algorithm == algorithms[j]) %>% 
        pull({{pathway_column}})
      
      # Calculate the overlap coefficient
      similarity_score <- overlap_coefficient(path1[[1]], path2[[1]])
      
      # Assign the computed score to both symmetric positions in the matrix
      similarity_matrix[i, j] <- similarity_score
      if (i != j) {
        similarity_matrix[j, i] <- similarity_score
      }
    }
  }
  
  return(similarity_matrix)
}


################################################################################
# SCRIPTS/METHOD_COMPARISONS/COMPARISONS.R FUNCTIONS
# Functions for method comparison analysis
################################################################################

#' Calculate Pathway Similarity Matrix
#'
#' @description
#' Computes a pairwise similarity matrix between algorithms based on their
#' aggregated pathway sets using the overlap coefficient.
#'
#' @param agg_list Named list of pathway vectors, one per algorithm.
#'
#' @return Square numeric matrix with overlap coefficients between all pairs
#'   of algorithms.
#'
#' @details
#' Simpler interface than calculate_similarities() - takes a named list directly
#' rather than a data frame with pathway columns.
#'
#' @source Found in:
#'   - Scripts/method_comparisons/Comparisons.R (line 389)
calculate_pathway_similarity_matrix <- function(agg_list) {
  algos <- names(agg_list)
  sim_matrix <- matrix(0, nrow = length(algos), ncol = length(algos),
                       dimnames = list(algos, algos))
  for (i in seq_along(algos)) {
    for (j in i:length(algos)) {
      sim_val <- overlap_coefficient(agg_list[[algos[i]]], agg_list[[algos[j]]])
      sim_matrix[i, j] <- sim_val
      sim_matrix[j, i] <- sim_val
    }
  }
  return(sim_matrix)
}


#' Calculate ARI Between Pathway Sets
#'
#' @description
#' Computes Adjusted Rand Index between two aggregated pathway sets by treating
#' pathway membership as a binary clustering assignment.
#'
#' @param agg_path1 Character vector of pathway IDs for first algorithm.
#' @param agg_path2 Character vector of pathway IDs for second algorithm.
#' @param all_pathways Character vector of all possible pathway IDs (universe).
#'
#' @return Numeric ARI value between -1 and 1.
#'
#' @details
#' Similar to calc_ari() but with different parameter naming. Creates binary
#' membership vectors and computes ARI using mclust::adjustedRandIndex().
#'
#' @source Found in:
#'   - Scripts/method_comparisons/Comparisons.R (line 544)
calculate_ari_pathway <- function(agg_path1, agg_path2, all_pathways) {
  # Create binary membership vectors for the union of all pathways
  vec1 <- as.integer(all_pathways %in% agg_path1)
  vec2 <- as.integer(all_pathways %in% agg_path2)
  
  # Load mclust and compute the ARI between these two binary vectors
  library(mclust)
  adjustedRandIndex(vec1, vec2)
}


#' Custom Cell Function for Combined ARI Heatmap
#'
#' @description
#' Cell drawing function for ComplexHeatmap that handles diagonal cells,
#' NA values, and dual-color schemes for cluster vs pathway ARI values.
#'
#' @param j Column index.
#' @param i Row index.
#' @param x X position for cell center.
#' @param y Y position for cell center.
#' @param width Cell width.
#' @param height Cell height.
#' @param fill Fill color (from color mapping).
#'
#' @return NULL (draws directly to graphics device).
#'
#' @details
#' Requires external variables: M (data matrix), type_mat (cluster/pathway indicator),
#' col_fun_cluster, col_fun_pathway (color mapping functions).
#' - Diagonal cells: white, no text
#' - NA cells: white
#' - .Machine$double.eps cells: white with "NA" text
#' - Other cells: colored by type with value text
#'
#' @note This is an inline function that depends on variables from the enclosing scope.
#'
#' @source Found in:
#'   - Scripts/method_comparisons/Comparisons.R (line 677)
cell_fun <- function(j, i, x, y, width, height, fill) {
  # For diagonal cells, draw white tiles and no text
  if (i == j) {
    grid.rect(x = x, y = y, width = width, height = height,
              gp = gpar(fill = "white", col = NA))
    return()
  }
  
  val <- M[i, j]
  if (is.na(val)) {
    grid.rect(x = x, y = y, width = width, height = height,
              gp = gpar(fill = "white", col = NA))
  } else if (val == .Machine$double.eps) {
    grid.rect(x = x, y = y, width = width, height = height,
              gp = gpar(fill = "white", col = NA))
    grid.text("NA", x, y, gp = gpar(col = "black", fontsize = 6))
  } else {
    cell_color <- if (type_mat[i, j] == "cluster") {
      col_fun_cluster(val)
    } else {
      col_fun_pathway(val)
    }
    grid.rect(x = x, y = y, width = width, height = height,
              gp = gpar(fill = cell_color, col = NA))
    grid.text(sprintf("%.2f", val), x, y, 
              gp = gpar(col = "black", fontsize = 6))
  }
}


#' PCA from ARI Matrix
#'
#' @description
#' Performs kernel PCA on an ARI similarity matrix and creates a 2D visualization
#' of algorithm relationships with category-based coloring.
#'
#' @param sim_matrix Square numeric similarity/ARI matrix.
#' @param clust_res Data frame with algorithm metadata including Category column.
#' @param cluster_colors Named vector of colors for each category.
#' @param output_path Character path for saving output (currently unused).
#' @param pointsize Numeric size for plot points (default: 2).
#' @param geom_label_size Numeric size for text labels (default: 2).
#'
#' @return A ggplot object showing the kernel PCA projection.
#'
#' @details
#' Modified version of pca_from_sim_matrix() from custom_functions.R.
#' Uses kernel PCA (eigendecomposition of centered kernel matrix) rather than
#' standard PCA. Points are colored by algorithm category and labeled with
#' algorithm names using ggrepel for non-overlapping labels.
#'
#' @source Found in:
#'   - Scripts/method_comparisons/Comparisons.R (line 1157)
pca_from_ari_matrix <- function(sim_matrix = NULL, clust_res = NULL,
                                cluster_colors = NULL, output_path = NULL,
                                pointsize = NULL,
                                geom_label_size = NULL) {
  
  # Kernel PCA of the clusters (modified from kernel_pca() from Spectrum)
  km <- sim_matrix
  m <- nrow(km)
  kc <- t(t(km - colSums(km)/m) - rowSums(km)/m) + sum(km)/m^2
  res <- eigen(kc/m, symmetric = TRUE)
  features <- m
  ret <- suppressWarnings(t(t(res$vectors[, 1:features]) / sqrt(res$values[1:features])))
  scores <- data.frame(ret)
  rownames(scores) <- colnames(km)
  colnames(scores)[1:2] <- c("PC1", "PC2")
  scores$algorithm <- rownames(scores)
  
  # Filter the scores dataset for annotation by joining with clust_res
  plot_df <- scores %>% inner_join(clust_res, by = "algorithm")
  rownames(plot_df) <- plot_df$algorithm
  
  if(is.null(pointsize)) { pointsize <- 2 }
  if(is.null(geom_label_size)) { geom_label_size <- 2 }
  
  # Actual plot: points are now colored by Category.
  kernelPCA <- ggplot(plot_df, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = Category), size = pointsize) +
    ggrepel::geom_text_repel(aes(label = algorithm),
                    size = pointsize,
                    max.overlaps = Inf, min.segment.length = 0,
                    segment.size = 0.1) +
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


################################################################################
# TOP-LEVEL SCRIPTS FUNCTIONS
# Functions from Scripts/Download_TCGA_data.R and Scripts/RNAseq_clean_counts.R
################################################################################

#' Check for Duplicate Row Names in Data Object
#'
#' @description
#' Examines each matrix in a list (data object) to identify any duplicated
#' row names. Returns a list of duplicates found in each matrix.
#'
#' @param data_object List of matrices/data frames to check.
#'
#' @return Named list where each element contains the duplicated row names
#'   found in that matrix. Empty list if no duplicates found.
#'
#' @details
#' Iterates through each element of the input list and checks for duplicated
#' row names using base R's duplicated() function. Only reports matrices
#' where duplicates were found.
#'
#' @source Found in:
#'   - Scripts/Download_TCGA_data.R (line 497)
#'   - Scripts/RNAseq_clean_counts.R (line 14)
check_duplicate_rownames <- function(data_object) {
  duplicate_info <- list()  # Initialize an empty list to store results
  
  for (i in seq_along(data_object)) {
    df <- data_object[[i]]
    if (!is.null(rownames(df))) {  # Check if the matrix has row names
      duplicates <- duplicated(rownames(df))  # Check for duplicated row names
      if (any(duplicates)) {
        duplicate_info[[paste("Matrix", i)]] <- rownames(df)[duplicates]
      }
    }
  }
  
  return(duplicate_info)
}


#' Filter to Most Variant Rows for Duplicate Row Names
#'
#' @description
#' Handles duplicate row names in multi-omics data by keeping the row with
#' highest variance. In case of ties, uses averaging for RNAseq/Methylation
#' and median for CNV data.
#'
#' @param data_object Named list of matrices where names indicate data type
#'   (e.g., "RNAseq", "CNV", "Methylation").
#'
#' @return Named list of filtered matrices with unique row names.
#'
#' @details
#' For each unique row name:
#' - If only one row exists: keep it
#' - If multiple rows: keep the one with maximum variance
#' - If variance ties exist:
#'   - RNAseq/Methylation: average the tied rows
#'   - CNV: take the median of tied rows
#'
#' @note Requires the matrixStats package for rowVars() and colMedians().
#'
#' @source Found in:
#'   - Scripts/Download_TCGA_data.R (line 528)
#'   - Scripts/RNAseq_clean_counts.R (line 43)
filter_most_variant_rows <- function(data_object) {
  # Initialize a list to store the processed matrices
  filtered_data_object <- list()
  
  for (matrix_name in names(data_object)) {
    mat <- data_object[[matrix_name]]
    
    # Find unique rownames and their indices
    unique_rownames <- unique(rownames(mat))
    filtered_mat <- do.call(rbind, lapply(unique_rownames, function(rname) {
      # Get the indices of rows with this rowname
      idx <- which(rownames(mat) == rname)
      
      if (length(idx) == 1) {
        # If there's only one row with this name, keep it
        return(mat[idx, , drop = FALSE])
      } else {
        # Calculate the variance for each row
        row_vars <- rowVars(as.matrix(mat[idx, ]))
        
        # Identify the rows with the maximum variance
        max_var_rows <- idx[row_vars == max(row_vars)]
        
        if (length(max_var_rows) == 1) {
          # If only one row has the max variance, keep it
          return(mat[max_var_rows, , drop = FALSE])
        } else {
          # Handle ties differently based on matrix type
          if (matrix_name %in% c("RNAseq", "Methylation")) {
            # Average the rows for RNAseq and Methylation
            return(colMeans(mat[max_var_rows, , drop = FALSE]))
          } else if (matrix_name == "CNV") {
            # Use median for CNV
            return(colMedians(as.matrix(mat[max_var_rows, , drop = FALSE])))
          }
        }
      }
    }))
    
    # Assign rownames to the filtered matrix
    rownames(filtered_mat) <- unique_rownames
    # Store the filtered matrix back in the list
    filtered_data_object[[matrix_name]] <- filtered_mat
  }
  
  return(filtered_data_object)
}


# ==============================================================================
# END OF DOCUMENTATION FILE
# ==============================================================================
