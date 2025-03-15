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
                       column="UNIPROT", keytype="SYMBOL", multiVals="list")
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

# Remove miRNAs
up_input = up_input[c("RNAseq", "CNV", "SNPs", "Methylation")]

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

up_input = lapply(up_input, fix_uniprot_row_names)

# # Due to large running times keep the top 10% of features ranked by MAD
# # PAMOGK however works with pathway/interaction network information from UniProt
# # meaning we might be dropping a significant amount of biologically relevant information here
# 
# # for continuous datasets
# up_input[c("RNAseq", "CNV", "Methylation")] = lapply(up_input[c("RNAseq", "CNV", "Methylation")], function(x) {
#   rowMAD <- apply(x, 1, mad, na.rm = TRUE)
#   ordered_rows <- order(rowMAD, decreasing = TRUE)
#   n_top <- ceiling(0.1 * nrow(x))
#   x[ordered_rows[1:n_top], ]
# })
# 
# # Keep cancer drivers and the top 10% most mutated genes of the rest for SNPs
# COSMIC_BC_drivers = read.csv("Resources/COSMIC_CGC_Breast_somatic.csv")$Gene.Symbol %>%
#   as.character()
# 
# snp_mat <- up_input[["SNPs"]]
# mutation_counts <- rowSums(snp_mat, na.rm = TRUE)
# non_drivers <- setdiff(rownames(snp_mat), COSMIC_BC_drivers)
# ordered_non_drivers <- non_drivers[order(mutation_counts[non_drivers], decreasing = TRUE)]
# n_top <- ceiling(0.1 * length(non_drivers))
# top_non_drivers <- ordered_non_drivers[1:n_top]
# genes_to_keep <- unique(c(COSMIC_BC_drivers, top_non_drivers))
# up_input[["SNPs"]] <- snp_mat[intersect(rownames(snp_mat), genes_to_keep), ]


################################################################################
# Find intersect of TCGA features with the json files in PAMOGK
PAMOGK_UP_IDs = as.character(data.table::fread("Python/PAMOGK/unique_uniprot_ids.txt")$V1)
# TCGA_UP_IDs = unique(unlist(lapply(up_input, rownames)))
# intersection = intersect(PAMOGK_UP_IDs, TCGA_UP_IDs)

library(matrixStats)

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
      row_mads <- rowMads(as.matrix(df_filtered))
      df_filtered <- df_filtered[order(row_mads, decreasing = TRUE), ]
      df_filtered <- df_filtered[!duplicated(rownames(df_filtered)), ]
    }
    
    return(df_filtered)
  })
  
  names(filtered_up_input) <- names(up_input)
  return(filtered_up_input)
}

# Apply the filter
up_input_filt <- filter_uniprot_ids(up_input, PAMOGK_UP_IDs)

# Export
for (modality in names(up_input_filt)) {
  write.csv(as.data.frame(up_input_filt[[modality]]), paste0("Python/PAMOGK/", modality, ".csv"))
}
