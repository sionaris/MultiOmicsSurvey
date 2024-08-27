# Download BRCA mRNA, miRNA, CNV, Methylation and Mutation data
source("Scripts/automated_scripts/custom_functions.R")
library(TCGAbiolinks)
library(dplyr)

# Exportable object with TCGA data initiation
data_object = list()

# Transcription data #####
query.exp.hg38 <- GDCquery(
  project = "TCGA-BRCA", 
  data.category = "Transcriptome Profiling", 
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts",
  access = "open",
  sample.type = "Primary Tumor"
)
GDCdownload(query.exp.hg38)
expdat <- GDCprepare(
  query = query.exp.hg38#,
  #save = TRUE, 
  #save.filename = "exp.rda"
)
saveRDS(expdat, "Resources/TCGA/RNA_full.rds")
rm(expdat, query.exp.hg38); gc()

# CNV data #####
CNV_query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Copy Number Variation",
  data.type = "Gene Level Copy Number",
  access = "open",
  sample.type = "Primary Tumor" # ABSOLUTE LiftOver data
)
GDCdownload(CNV_query)
CNV_data <- GDCprepare(CNV_query)

saveRDS(CNV_data, "Resources/TCGA/CNV_full.rds")
rm(CNV_data, CNV_query); gc()

# miRNA data #####
query.mirna <- GDCquery(
  project = "TCGA-BRCA", 
  experimental.strategy = "miRNA-Seq",
  data.category = "Transcriptome Profiling", 
  data.type = "miRNA Expression Quantification",
  access = "open"
)
GDCdownload(query.mirna)
mirna.dat <- GDCprepare(
  query = query.mirna
)

saveRDS(mirna.dat, "Resources/TCGA/miRNA_full.rds")
rm(mirna.dat, query.mirna); gc()

# Methylation data #####
query_met.hg38 <- GDCquery(
  project= "TCGA-BRCA", 
  data.category = "DNA Methylation", 
  data.type = "Methylation Beta Value",
  platform = "Illumina Human Methylation 450", 
  access = "open"
)
GDCdownload(query_met.hg38)
data.met.hg38 <- GDCprepare(query_met.hg38)

saveRDS(data.met.hg38, "Resources/TCGA/Methyl_full.rds")
rm(data.met.hg38, query_met.hg38); gc()

# Download categorical mutation data #####
library(maftools)
query_mut = GDCquery(
  project = "TCGA-BRCA",  # Replace with the appropriate TCGA project ID
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)

# Execute the query
GDCdownload(query_mut)

# Prepare the clinical data
mut_data = GDCprepare(query_mut)
mut_data = mut_data %>% maftools::read.maf()

# The @data object is filtered for silent mutations. We use it to create a 
# binary matrix

# Split for SNP and indels
library(stringr)
library(data.table)
SNP_data = mut_data@data %>% dplyr::filter(Variant_Type == "SNP")
INDEL_data = mut_data@data %>% dplyr::filter(Variant_Type %in% c("DEL", "INS"))

SNP_matrix_input = SNP_data %>% dplyr::select(Hugo_Symbol, Start_Position, Tumor_Sample_Barcode)
str_sub(SNP_matrix_input$Tumor_Sample_Barcode, 17, -1) = ""

INDEL_matrix_input = INDEL_data %>% 
  dplyr::select(Hugo_Symbol, Start_Position, End_Position, Tumor_Sample_Barcode)
str_sub(INDEL_matrix_input$Tumor_Sample_Barcode, 17, -1) = ""

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
rm(list=setdiff(ls(), "data_object")); gc()
source("Scripts/automated_scripts/custom_functions.R")

tcga_samples = str_sub(colnames(data_object[["SNPs"]]), 1, 12)

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

# Male samples (remove)
male_samples = clinical_data$bcr_patient_barcode[clinical_data$gender == "MALE"]
female_indices = which(!tcga_samples %in% male_samples)
data_object$SNPs = data_object$SNPs[, female_indices] # Only keep female samples from now on

# Load the rest of the objects and filter for common samples #####
RNAseq = readRDS("Resources/TCGA/RNA_full.rds")
miRNA = readRDS("Resources/TCGA/miRNA_full.rds")
count_cols = which(grepl("count", colnames(miRNA)))
miRNA = miRNA[, c(1, count_cols)]

# RNA and miRNA samples and SNPs
RNAsamples = RNAseq@colData@rownames
stringr::str_sub(RNAsamples, 17, 28) = ""
mirnasamples = gsub("read_count_", "", colnames(miRNA)[2:ncol(miRNA)])
str_sub(mirnasamples, 17, 28) = ""
overlap = intersect(mirnasamples, intersect(RNAsamples, colnames(data_object$SNPs)))

# Filter RNA
colnames(RNAseq) = RNAsamples
RNAseq = RNAseq[, overlap]
exprs = RNAseq@assays@data@listData[["unstranded"]]
rownames(exprs) = RNAseq@rowRanges@elementMetadata@listData[["gene_name"]]
colnames(exprs) = colnames(RNAseq)
rm(RNAseq); gc()

which(is.na(exprs)) # No NA
not_all_zeros = which(rowSums(exprs) != 0) # indices of rows that have non-zero sum of counts
exprs = exprs[not_all_zeros, ]; rm(not_all_zeros); gc()

# Make miRNA neat
colnames(miRNA) = c("miR_ID", mirnasamples)
miRNA = miRNA[, c("miR_ID", overlap)]; gc()
mirs = miRNA$miR_ID
miRNA = as.matrix(miRNA[, 2:ncol(miRNA)])
rownames(miRNA) = mirs
rm(mirnasamples, mirs, RNAsamples, count_cols); gc()

which(is.na(miRNA)) # No NA
not_all_zeros = which(rowSums(miRNA) != 0) # indices of rows that have non-zero sum of counts
miRNA = miRNA[not_all_zeros, ]; rm(not_all_zeros); gc()

# Load CNV data
cnv = readRDS("Resources/TCGA/CNV_full.rds")
cnv_dat = cnv@assays@data@listData[["copy_number"]]
colnames(cnv_dat) = cnv@colData@listData[["sample"]]
rownames(cnv_dat) = cnv@rowRanges@elementMetadata@listData[["gene_name"]]

# Filter for samples with all modalities
overlap = intersect(colnames(cnv_dat), overlap)
cnv_dat = cnv_dat[, overlap]
rows_with_na = apply(cnv_dat, 1, function(row) any(is.na(row)))
no_of_nas = apply(cnv_dat[rows_with_na, ], 1, function(row) sum(is.na(row)))

table(no_of_nas)

# All rows with >10 missing values are removed
cnv_dat = cnv_dat[apply(cnv_dat, 1, function(row) sum(is.na(row)) <= 10), ]

# The rest of the rows could be imputed with "2" which is the normal value
# or be removed as well. Here we choose the strict way

cnv_dat = na.omit(cnv_dat)
rm(cnv, no_of_nas, rows_with_na)
gc()

# Methylation data
Methyl = readRDS("Resources/TCGA/Methyl_full.rds")
methyl_dat = Methyl@assays@data@listData[[1]]
colnames(methyl_dat) = Methyl@colData@listData[["sample"]]
rownames(methyl_dat) = Methyl@rowRanges@elementMetadata@listData[["gene"]]

# Filter for samples with all modalities
overlap = intersect(colnames(methyl_dat), overlap)
methyl_dat = methyl_dat[, overlap]
rows_with_na = apply(methyl_dat, 1, function(row) any(is.na(row))) # 485577 rows with missing values
no_of_nas = apply(methyl_dat[rows_with_na, ], 1, function(row) sum(is.na(row)))

table(no_of_nas)

# Remove all rows with NAs in more than 25% of the samples
# Calculate the threshold for 25% missing values
threshold <- round(ncol(methyl_dat) * 0.25)

# Identify rows with less than 25% missing values
rows_to_keep <- apply(methyl_dat, 1, function(row) sum(is.na(row)) < threshold)

# Subset the matrix to keep only the rows with less than 25% missing values
methyl_dat <- methyl_dat[rows_to_keep, ]

# Clean up
rm(threshold, rows_to_keep, no_of_nas, rows_with_na, Methyl); gc()

# Impute missing values using the impute package
library(impute)
RNGversion("4.2.2")
methyl_dat = impute.knn(methyl_dat, k = 10, maxp = 1500,
                        rng.seed = 123, rowmax = 0.25, colmax = 0.8)$data

# Store everything in data_object
data_object$RNAseq = exprs; rm(exprs); gc()
data_object$CNV = cnv_dat; rm(cnv_dat); gc()
data_object$miRNA = miRNA; rm(miRNA); gc()
data_object$Methylation = methyl_dat; rm(methyl_dat); gc()

# Convert everything to matrices and filter for overlap samples
data_object = lapply(data_object, as.matrix)
data_object = lapply(data_object, function(matrix) matrix[, overlap, drop = FALSE])

# Total: 625 samples

# Preprocessing #####

# Methylation
# We will convert to M-values using the sesame package (which is used by TCGA). Here is why:
# https://life-epigenetics-methylprep.readthedocs-hosted.com/en/latest/docs/introduction/introduction.html#:~:text=The%20M%2Dvalue%20does%20not,because%20they%20are%20more%20intuitive.

library(sesame)
data_object$Methylation = BetaValueToMValue(data_object$Methylation)

# RNAseq and miRNAseq will be converted to normalized counts using DESeq2
# Why? Here:
# 1. https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2016.00164/full
# 2. https://hbctraining.github.io/DGE_workshop_salmon/lessons/02_DGE_count_normalization.html

library(DESeq2)

req_norm = c("RNAseq", "miRNA")
size_factors = list()

for (i in req_norm) {
  count_matrix = data_object[[i]]
  
  # mock sample conditions
  sample_conditions <- data.frame(
    row.names = colnames(count_matrix),
    condition = factor(rep("condition", ncol(count_matrix)))
  )
  
  # Create DESeqDataSet object
  dds <- DESeqDataSetFromMatrix(countData = count_matrix, 
                                colData = sample_conditions, 
                                design = ~ 1) 
  
  dds <- estimateSizeFactors(dds)
  
  # Save size factors for later
  size_factors[[i]] = dds$sizeFactor
  
  # Get the normalized counts
  normalized_counts <- counts(dds, normalized = TRUE)
  rownames(normalized_counts) = rownames(data_object[[i]])
  colnames(normalized_counts) = colnames(data_object[[i]])
  
  data_object[[i]] = normalized_counts
}

rm(count_matrix, dds, normalized_counts, i, req_norm, sample_conditions); gc()

# Check densities prior to standardization #####
library(ggplot2)
library(tidyr)

# Density curve for RNA and miRNA normalized counts
create_density_curve(data_object$RNAseq)+
  scale_x_continuous(limits = c(-1, 30), breaks = seq(0, 30, 5))
ggsave(filename = "TCGA_RNAseq_normalized_counts_density_curves_xlim30.tiff",
       path = "Results/TCGA_exploratory",
       width = 1920, height = 1080, device = 'tiff', units = "px",
       dpi = 150, compression = "lzw")
dev.off()

create_density_curve(data_object$miRNA)+
  scale_x_continuous(limits = c(-1, 30), breaks = seq(0, 30, 5))
ggsave(filename = "TCGA_miRNA_normalized_counts_density_curves_xlim30.tiff",
       path = "Results/TCGA_exploratory",
       width = 1920, height = 1080, device = 'tiff', units = "px",
       dpi = 150, compression = "lzw")
dev.off()

# Color density plots
create_density_plot_color(data_object$RNAseq)+
  scale_x_continuous(limits = c(-1, 30), breaks = seq(0, 30, 5))
ggsave(filename = "TCGA_RNAseq_normalized_counts_color_densities_xlim30.tiff",
       path = "Results/TCGA_exploratory",
       width = 1920, height = 1080, device = 'tiff', units = "px",
       dpi = 150, compression = "lzw")
dev.off()

create_density_plot_color(data_object$miRNA)+
  scale_x_continuous(limits = c(-1, 30), breaks = seq(0, 30, 5))
ggsave(filename = "TCGA_miRNA_normalized_counts_color_densities_xlim30.tiff",
       path = "Results/TCGA_exploratory",
       width = 1920, height = 1080, device = 'tiff', units = "px",
       dpi = 150, compression = "lzw")
dev.off()

# To avoid the extreme effect of outliers on the standardization we first
# log-transform with an added pseudocount determined as Lun et al. suggest here:
# https://www.biorxiv.org/content/10.1101/404962v1.full

# Clean up clinical data
clinical_data = clinical_data[clinical_data$bcr_patient_barcode %in%
                                substring(overlap, 0, 12), ] %>%
  dplyr::rename(Patient.ID = bcr_patient_barcode) %>%
  group_by(Patient.ID) %>%
  arrange(Patient.ID, rowSums(is.na(across(-Patient.ID)))) %>%  # Arrange by Sample.ID and NA count
  distinct(Patient.ID, .keep_all = TRUE) %>%  # Keep the first occurrence in case of ties
  ungroup()


# Add a Sample.ID column
samples_df = as.data.frame(list(Sample.ID = overlap))
samples_df$Patient.ID = substring(samples_df$Sample.ID, 0, 12)

clinical_data = clinical_data %>%
  inner_join(samples_df, by = "Patient.ID") %>%
  dplyr::select(Sample.ID, Patient.ID, everything())

openxlsx::write.xlsx(clinical_data, "Resources/TCGA/clinical_data.xlsx",
                     overwrite = TRUE)
rm(query, samples_df, male_samples, female_indices, tcga_samples); gc()

# We need a grouping for the breast cancer samples. We are going to go for 
# ER+, HER2+, TNBC
her2_samples = clinical_data$Sample.ID[clinical_data$breast_carcinoma_estrogen_receptor_status == "Negative" &
                                       clinical_data$lab_proc_her2_neu_immunohistochemistry_receptor_status %in%
                                         c("Equivocal", "Positive")]
tnbc_samples = clinical_data$Sample.ID[clinical_data$breast_carcinoma_estrogen_receptor_status == "Negative" &
                                         clinical_data$breast_carcinoma_progesterone_receptor_status != "Positive" &
                                         clinical_data$lab_proc_her2_neu_immunohistochemistry_receptor_status == "Negative"]
er_samples = setdiff(clinical_data$Sample.ID, c(her2_samples, tnbc_samples))

# We can therefore split in ER+ and non-ER+ for the optimal pseudocount estimation

# miRNA
avg_er_size_factor_miRNA = mean(size_factors$miRNA[er_samples])
avg_else_size_factor_miRNA = mean(size_factors$miRNA[c(her2_samples, tnbc_samples)])
sf_miRNA = c(avg_er_size_factor_miRNA, avg_else_size_factor_miRNA)

# RNAseq
avg_er_size_factor_RNAseq = mean(size_factors$RNAseq[er_samples])
avg_else_size_factor_RNAseq = mean(size_factors$RNAseq[c(her2_samples, tnbc_samples)])
sf_RNA = c(avg_er_size_factor_RNAseq, avg_else_size_factor_RNAseq)

# Lun et al. suggest pseudocount = max{1, r|1/smin - 1/smax|}, where r = 1 (suggestion)
pseudocount_miRNA = max(c(1, abs(1/min(sf_miRNA) - 1/max(sf_miRNA)))) # 1
pseudocount_RNA = max(c(1, abs(1/min(sf_RNA) - 1/max(sf_RNA)))) # 1

# We proceed with the log2(norm.counts + 1) transformation
data_object$RNAseq = log2(data_object$RNAseq + pseudocount_RNA)
data_object$miRNA = log2(data_object$miRNA + pseudocount_miRNA)

# Clean up
rm(avg_else_size_factor_miRNA, avg_else_size_factor_RNAseq,
   avg_er_size_factor_miRNA, avg_er_size_factor_RNAseq,
   sf_miRNA, sf_RNA, pseudocount_miRNA, pseudocount_RNA,
   er_samples, her2_samples, tnbc_samples); gc()

# Which modalities have duplicated rownames?
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

res = check_duplicate_rownames(data_object)

# Remove rows with NA rownames from Methylation
noname = which(is.na(rownames(data_object$Methylation)))
data_object$Methylation = data_object$Methylation[-noname,]

# There are duplicate rownames in RNAseq, CNV and Methylation
# We keep the most variant rows
# In the case of ties, we use averaging for RNAseq, Methylation
# Keep the median in CNV

# Load necessary libraries
library(matrixStats)

# Function to handle duplicates and filtering
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

filtered_data_object <- filter_most_variant_rows(data_object[c("RNAseq", "CNV", "Methylation")])
data_object$RNAseq = filtered_data_object$RNAseq
data_object$CNV = filtered_data_object$CNV
data_object$Methylation = filtered_data_object$Methylation
rm(filtered_data_object, noname, res, size_factors); gc()

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

res = check_missing_values_with_indices(data_object)

# Methylation has 72 rows with missing values. We remove those
data_object$Methylation = data_object$Methylation[-res$Methylation$missing_indices, ]
rm(res); gc()

# All data types apart from SNPs will be standardized feature-wise
req_standardize = c("RNAseq", "Methylation", "CNV", "miRNA")

# Check if there are rows with only zeros
# Assuming data_object is a list of matrices
zero_rows_info <- lapply(data_object, function(mat) {
  if (!is.matrix(mat)) {
    stop("All elements of data_object should be matrices.")
  }
  
  # Find row indices where all values are zero
  zero_indices <- which(rowSums(mat == 0) == ncol(mat))
  
  # Return the indices of rows that contain only zeros
  return(zero_indices)
})

# Assign names to the list for clarity
names(zero_rows_info) <- names(data_object)

# Remove the identified rows
for (name in names(zero_rows_info)) {
  zero_indices <- zero_rows_info[[name]]
  
  if (length(zero_indices) > 0) {
    # Remove the rows with only zeros from the corresponding matrix
    data_object[[name]] <- data_object[[name]][-zero_indices, , drop = FALSE]
    message(paste("Removed", length(zero_indices), "rows from", name))
  } else {
    message(paste("No rows with only zeros in", name))
  }
}

for (matrix_name in names(data_object)) {
  if (matrix_name %in% req_standardize) {
    # Standardize the rows for matrices that need to be standardized
    data_object[[matrix_name]] <- standardize_rows(data_object[[matrix_name]])
  }
}

# Check the data_object for any missing values, one last time
res = check_missing_values_with_indices(data_object)

# 2 rows with missing values in Methylation. We remove them
data_object$Methylation = data_object$Methylation[-res$Methylation$missing_indices, ]
rm(res); gc()

# Export object
saveRDS(data_object, "Resources/TCGA/norm_data_object.rds"); gc()
rm(zero_rows_info, matrix_name, name, zero_indices); gc()
save.image("Resources/TCGA/tcga_env.RData")

# Export session info as .txt
writeLines(capture.output(sessionInfo()), "Resources/TCGA/TCGA_preprocess_sessionInfo.txt")
