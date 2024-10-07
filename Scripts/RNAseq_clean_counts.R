# Load data
RNAfull = readRDS("Resources/TCGA/RNA_full.rds")
data_object = readRDS("Resources/TCGA/norm_data_object.rds")

# Extract counts, give names, filter for zero rows
RNAcounts = RNAfull@assays@data@listData[["unstranded"]]
rownames(RNAcounts) = RNAfull@rowRanges@elementMetadata@listData[["gene_name"]]
colnames(RNAcounts) = colnames(RNAfull)
which(is.na(RNAcounts)) # No NA
not_all_zeros = which(rowSums(RNAcounts) != 0) # indices of rows that have non-zero sum of counts
RNAcounts = RNAcounts[not_all_zeros, ]; rm(not_all_zeros); gc()

# Check duplicate rownames
check_duplicate_rownames <- function(rna_object) {
  duplicate_info <- list()  # Initialize an empty list to store results
  
  for (i in seq_along(rna_object)) {
    df <- rna_object[[i]]
    if (!is.null(rownames(df))) {  # Check if the matrix has row names
      duplicates <- duplicated(rownames(df))  # Check for duplicated row names
      if (any(duplicates)) {
        duplicate_info[[paste("Matrix", i)]] <- rownames(df)[duplicates]
      }
    }
  }
  
  return(duplicate_info)
}

# Filter for the colnames of interest
rna_cols = colnames(RNAcounts)
rna_cols = substr(rna_cols, 0, 16)
colnames(RNAcounts) = rna_cols
RNAcounts = RNAcounts[, colnames(data_object$RNAseq)]

rna_object = list(RNAseq = RNAcounts)
res = check_duplicate_rownames(rna_object)

# Load necessary libraries
library(matrixStats)

# Function to handle duplicates and filtering
filter_most_variant_rows <- function(rna_object) {
  # Initialize a list to store the processed matrices
  filtered_data_object <- list()
  
  for (matrix_name in names(rna_object)) {
    mat <- rna_object[[matrix_name]]
    
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

filtered_data_object <- filter_most_variant_rows(rna_object)

# Check for missing values
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

res = check_missing_values_with_indices(rna_object)

# Check if there are rows with only zeros
# Assuming rna_object is a list of matrices
zero_rows_info <- lapply(rna_object, function(mat) {
  if (!is.matrix(mat)) {
    stop("All elements of rna_object should be matrices.")
  }
  
  # Find row indices where all values are zero
  zero_indices <- which(rowSums(mat == 0) == ncol(mat))
  
  # Return the indices of rows that contain only zeros
  return(zero_indices)
})

# Assign names to the list for clarity
names(zero_rows_info) <- names(rna_object)

# Remove the identified rows
for (name in names(zero_rows_info)) {
  zero_indices <- zero_rows_info[[name]]
  
  if (length(zero_indices) > 0) {
    # Remove the rows with only zeros from the corresponding matrix
    rna_object[[name]] <- rna_object[[name]][-zero_indices, , drop = FALSE]
    message(paste("Removed", length(zero_indices), "rows from", name))
  } else {
    message(paste("No rows with only zeros in", name))
  }
}

# Check consistency with data object
nrow(rna_object$RNAseq)
nrow(data_object$RNAseq)

# Keep the same genes as in the normalized case
rna_object$RNAseq = rna_object$RNAseq[rownames(data_object$RNAseq), ]

# Export counts object
saveRDS(rna_object$RNAseq, "Resources/TCGA/counts_for_DESeq2.rds")
