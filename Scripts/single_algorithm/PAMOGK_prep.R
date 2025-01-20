input = readRDS("Resources/TCGA/mm_input.rds")
up_input = input

library(org.Hs.eg.db)

# Base R for renaming miRNAs
mirnas = rownames(input[["miRNA"]])
substs = gsub("^hsa-mir-", "MIR", mirnas)
substs = gsub("^hsa-let-", "MIRLET", substs)
substs <- gsub("-", "", substs)
substs <- gsub("([a-z])", "\\U\\1", substs, perl=TRUE)
rownames(input[["miRNA"]]) = substs

# Function to clean row names
clean_rownames <- function(df) {
  rownames(df) <- sub("\\.\\d+$", "", rownames(df))
  return(df)
}

input <- lapply(input, clean_rownames)

for (i in 1:length(input)) {
  gene_symbols = rownames(input[[i]])
  uniprot_ids = mapIds(org.Hs.eg.db, keys=gene_symbols,
                       column="UNIPROT", keytype="SYMBOL", multiVals="first")
  rownames(up_input[[i]]) = uniprot_ids
}

# Function to remove rows with NA row names
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

# Apply the function to each element in 'up_input'
up_input <- lapply(up_input, remove_na_rownames)

# Removing 267 rows with NA row names. - SNPs
# Removing 37847 rows with NA row names. - RNAseq
# Removing 22484 rows with NA row names. - CNV
# Removing 1568 rows with NA row names. - miRNA (no proteins in miRNA)
# Removing 22022 rows with NA row names. - Methylation

# Function to get dimensions or indicate NULL
get_dimensions <- function(x) {
  if (is.null(x)) {
    return("NULL")
  } else {
    dims <- dim(x)
    return(paste(dims[1], "rows x", dims[2], "columns"))
  }
}

# Apply the function to each element in the list
dimensions <- sapply(up_input, get_dimensions)

# Convert to a data frame for better readability
dimensions_df <- data.frame(
  Element = names(dimensions),
  Dimensions = dimensions,
  stringsAsFactors = FALSE
)

# Print the data frame
print(dimensions_df)

# Element               Dimensions
# SNPs               SNPs 12659 rows x 625 columns
# RNAseq           RNAseq 18870 rows x 625 columns
# CNV                 CNV 13374 rows x 625 columns
# miRNA             miRNA     0 rows x 625 columns
# Methylation Methylation 14739 rows x 625 columns