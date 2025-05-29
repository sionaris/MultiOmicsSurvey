library(MOVICS)
library(TCGAbiolinks)
library(dplyr)
library(DESeq2)
library(matrixStats)
library(crayon)
library(survminer)
library(survival)
library(openxlsx)

# Ensure reproducibility
RNGversion("4.2.2")
set.seed(123)

# Load custom helper functions
source("Scripts/automated_scripts/custom_functions.R")
source("Scripts/automated_scripts/modified_MOVICS_functions.R")

# Set home directory
home = getwd()

# This RDS object was produced using the Scripts/MOVICS/MOVICS_baseline.R script
input = readRDS("Resources/TCGA/mm_input.rds")
train_samples = colnames(input$SNPs) # These are not barcodes, but we will work around it

# Query the TCGA-BRCA project for all samples with expression data in RNAseq
query.exp.hg38 <- GDCquery(
  project = "TCGA-BRCA", 
  data.category = "Transcriptome Profiling", 
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts",
  access = "open",
  sample.type = "Primary Tumor"
)

# Filter query for non-overlapping samples
query_filt = query.exp.hg38
query_filt[[1]][[1]]$surv_samples <- sapply(query_filt[[1]][[1]]$cases, function(x) {
  parts <- strsplit(x, "-")[[1]]
  paste(head(parts, 4), collapse = "-")
})
query_filt[[1]][[1]] = query_filt[[1]][[1]] %>%
  dplyr::filter(!surv_samples %in% train_samples)

# Download data
# system("subst x: \"C:/Users/user/path/to/dir\"") - to avoid the character 260 limit
GDCdownload(query_filt, directory = paste0("x://", "Surv"))

# Prepare data
expdat <- GDCprepare(
  query = query_filt,
  directory = paste0("x://", "Surv")#,
  #save = TRUE, 
  #save.filename = "exp.rda"
)
rm(query.exp.hg38, query_filt); gc()

# Create a query to retrieve clinical data for the specified samples #####
query <- GDCquery(
  project = "TCGA-BRCA",  # Replace with the appropriate TCGA project ID
  data.category = "Clinical",
  data.type = "Clinical Supplement",
  barcode = expdat@colData@listData[["patient"]]
)

# Execute the query
GDCdownload(query, directory = paste0("x://", "Surv"))

# Prepare the clinical data
holdout_clinical_data = GDCprepare_clinic(query, clinical.info = "patient",
                                  directory = paste0("x://", "Surv"))

# Add a Sample.ID column
samples_df = as.data.frame(list(Sample.ID = expdat@colData@listData[["sample"]],
                                Long.barcode = expdat@colData@listData[["barcode"]]))
samples_df$Patient.ID = expdat@colData@listData[["patient"]]

holdout_clinical_data = holdout_clinical_data %>%
  dplyr::rename(Patient.ID = bcr_patient_barcode) %>%
  inner_join(samples_df, by = "Patient.ID", relationship = "many-to-many") %>%
  dplyr::select(Sample.ID, Patient.ID, Long.barcode, everything())

# Clean up clinical data
holdout_clinical_data = holdout_clinical_data %>%
  dplyr::filter(gender == "FEMALE") %>%
  group_by(Patient.ID) %>%
  arrange(Patient.ID, rowSums(is.na(across(-Patient.ID)))) %>%  # Arrange by Sample.ID and NA count
  distinct(Patient.ID, .keep_all = TRUE) %>%  # Keep the first occurrence in case of ties
  ungroup()

# Filter expdat for these samples
expdat = expdat[, holdout_clinical_data$Long.barcode]
colnames(expdat) = holdout_clinical_data$Sample.ID

# Export
saveRDS(expdat, "Resources/TCGA/Surv_RNA_full.rds")
rm(query, samples_df); gc()

# Bring system back to normal
# system("subst x: /D")

# Extract expression matrix
transcr_pre = expdat@assays@data@listData[["unstranded"]]
rownames(transcr_pre) = expdat@rowRanges@elementMetadata@listData[["gene_id"]] # temporary Ensembl IDs
colnames(transcr_pre) = colnames(expdat)

# Normalize counts using DESeq2 #####

# mock sample conditions
sample_conditions <- data.frame(
  row.names = colnames(transcr_pre),
  condition = factor(rep("condition", ncol(transcr_pre)))
)

# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = transcr_pre, 
                              colData = sample_conditions, 
                              design = ~ 1) 

dds <- estimateSizeFactors(dds)

# Save size factors for later
size_factors = dds$sizeFactor

# Get the normalized counts
normalized_counts <- counts(dds, normalized = TRUE)
rownames(normalized_counts) = rownames(transcr_pre)
colnames(normalized_counts) = colnames(transcr_pre)

rm(dds, sample_conditions); gc()

# To avoid the extreme effect of outliers on the standardization we first
# log-transform with an added pseudocount determined as Lun et al. suggest here:
# https://www.biorxiv.org/content/10.1101/404962v1.full

# We need a grouping for the breast cancer samples. We are going to go for 
# ER+, HER2+, TNBC
her2_samples = holdout_clinical_data$Sample.ID[holdout_clinical_data$breast_carcinoma_estrogen_receptor_status == "Negative" &
                                         holdout_clinical_data$lab_proc_her2_neu_immunohistochemistry_receptor_status %in%
                                         c("Equivocal", "Positive")]
tnbc_samples = holdout_clinical_data$Sample.ID[holdout_clinical_data$breast_carcinoma_estrogen_receptor_status == "Negative" &
                                         holdout_clinical_data$breast_carcinoma_progesterone_receptor_status != "Positive" &
                                         holdout_clinical_data$lab_proc_her2_neu_immunohistochemistry_receptor_status == "Negative"]
er_samples = setdiff(holdout_clinical_data$Sample.ID, c(her2_samples, tnbc_samples))

# We can therefore split in ER+ and non-ER+ for the optimal pseudocount estimation

# RNAseq
avg_er_size_factor_RNAseq = mean(size_factors[er_samples])
avg_else_size_factor_RNAseq = mean(size_factors[c(her2_samples, tnbc_samples)])
sf_RNA = c(avg_er_size_factor_RNAseq, avg_else_size_factor_RNAseq)

# Lun et al. suggest pseudocount = max{1, r|1/smin - 1/smax|}, where r = 1 (suggestion)
pseudocount_RNA = max(c(1, abs(1/min(sf_RNA) - 1/max(sf_RNA)))) # 1

# We proceed with the log2(norm.counts + opt.pseudocount) transformation
log2_normalized_counts = log2(normalized_counts + pseudocount_RNA)

# Clean up
rm(avg_else_size_factor_RNAseq, avg_er_size_factor_RNAseq,
   sf_RNA, pseudocount_RNA,
   er_samples, her2_samples, tnbc_samples); gc()

# Give gene names to the rows and filter for most variant rows in case of duplicates #####
rownames(log2_normalized_counts) = expdat@rowRanges@elementMetadata@listData[["gene_name"]]

# There are a few duplicate rownames (1233) in the transcr object. We keep the most 
# variant row in each case
length(which(duplicated(rownames(log2_normalized_counts))))

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

log2_filt = filter_most_variant_rows(list(RNAseq = log2_normalized_counts))
log2_filt = as.matrix(log2_filt[[1]])

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

res = check_missing_values_with_indices(list(RNAseq = log2_filt)) # no missing values

# Check if there are rows with only zeros
# Assuming data_object is a list of matrices
zero_rows_info <- lapply(list(RNAseq = log2_filt), function(mat) {
  if (!is.matrix(mat)) {
    stop("All elements of data_object should be matrices.")
  }
  
  # Find row indices where all values are zero
  zero_indices <- which(rowSums(mat == 0) == ncol(mat))
  
  # Return the indices of rows that contain only zeros
  return(zero_indices)
})

names(zero_rows_info) <- "RNAseq" # 2722 zero only rows

# Remove the identified rows
if (length(zero_rows_info$RNAseq) > 0) {
  # Remove the rows with only zeros from the corresponding matrix
  log2_filt <- log2_filt[-zero_rows_info$RNAseq, , drop = FALSE]
  message(paste("Removed", length(zero_rows_info$RNAseq), "rows from RNAseq"))
} else {
  message(paste("No rows with only zeros in RNAseq"))
}

z_transcr <- standardize_rows(log2_filt)

# Check the data_object for any missing values, one last time
res = check_missing_values_with_indices(list(RNAseq = z_transcr)) # No missing values

# Export transcr
saveRDS(z_transcr, "Resources/TCGA/Surv_standardized_transcr.rds"); gc()


# Export all input required by the NTP/PAM mapping algorithms and run the 
# label mapping in an HPC environment
algorithms = c("ab-SNF", "ANF", "CIMLR", "COCA", "iClusterBayes", "KLIC",
               "LRAcluster", "MDICC", "MFA",
               "NEMO", "RWR-F", "RWR-NF", "SNF", "Spectrum", "wMKL",
               "MONET", "MSNE", "MOFA")
env_objects = list()

for (algorithm in algorithms) {
  
  # .RData object location
  env_file = paste0("Results/single_algorithm/", 
                    algorithm, "/", algorithm,
                    "_TCGA_RNAseq-CNV-Methylation-miRNA-SNPs_eval_on_transNEO_env.RData")
  
  # check if file exists
  if (file.exists(env_file)) {
    
    # load environment
    load(env_file)
    rm(list=setdiff(ls(), c("z_transcr", "home", "plot_object", "algorithm", "algorithms",
                            "plotdata", "train_samples", "holdout_clinical_data", "expdat",
                            "log2_filt", "env_objects")))
    
    env_objects[[algorithm]] = list(plotdata = plotdata,
                                    plot_object = plot_object)
    
    rm(list=setdiff(ls(), c("home", "algorithm", "algorithms", "z_transcr", "train_samples",
                            "holdout_clinical_data", "expdat",
                            "log2_filt", "env_objects")))
    cat(blue("Done with", algorithm, "\n"))
  } else {
    cat(blue("Skipped", algorithm, "\n"))
  }
  invisible(gc())
}

# Export object for usage in an HPC environment
saveRDS(env_objects, "Results/Survival_evaluations/env_objects.rds")

# Now filter env_objects only for RNA data (only ones we use for NTP/PAM)
env_objects <- lapply(env_objects, function(x) {
    x$plotdata <- x$plotdata["RNAseq"]
    return(x)
})

# Clean up
rm(input, log2_normalized_counts, normalized_counts, res, transcr_pre, 
   zero_rows_info, size_factors); gc()

# Label mapping loop ###
label_maps = list()
for (algorithm in algorithms) {
  
  # Load custom helper functions
  source("Scripts/automated_scripts/custom_functions.R")
  source("Scripts/automated_scripts/modified_MOVICS_functions.R")
  
  # get as many templates as possible
  dgea.marker.up_1000 <- runMarker_single_algorithm_no_export(algorithm_name = algorithm,
                                                              moic.res = env_objects[[algorithm]]$plot_object,
                                                              n.marker = 1000,
                                                              dea.method    = "limma", # name of DEA method
                                                              prefix        = "dgea_", # MUST be the same of argument in runDEA()
                                                              dat.path      = paste0(home, "/Results/single_algorithm/", algorithm), # path of DEA files
                                                              p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                                                              p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                                                              norm.expr = env_objects[[algorithm]]$plotdata$RNAseq,
                                                              dirct         = "up" # direction of dysregulation in expression
  )
  
  # If the directory of the algorithm does not exist, create it
  newdir = paste0(home, "/Results/Survival_evaluations/", algorithm)
  if (!dir.exists(newdir)) {
    dir.create(newdir)
  }
  
  # If NTP worked for the algorithm in the initial runs, use NTP;
  # otherwise use PAM
  if (TRUE %in% grepl("ntp_expr_up", 
                      list.files(paste0("Results/single_algorithm/",
                                        algorithm)))) {
    # NTP
    TCGA_ntp_expr_up = runNTP(
      expr = z_transcr[intersect(dgea.marker.up_1000$templates$probe,
                                 rownames(z_transcr)), ],
      templates = dgea.marker.up_1000$templates,
      scaleFlag = TRUE,
      centerFlag = TRUE,
      nPerm = 10000,
      seed = 123,
      distance = "cosine", # default
      doPlot = TRUE,
      height = 8,
      width = 12,
      fig.path = paste0(home, "/Results/Survival_evaluations/", algorithm),
      fig.name = paste0(algorithm, "_ntp_expr_up_heatmap_TCGA"))
    
    df = as.data.frame(TCGA_ntp_expr_up$clust.res) %>%
      dplyr::select(samID, clust)
    colnames(df) = c("Sample.ID", algorithm)
    label_maps[[algorithm]] = df %>%
      inner_join(holdout_clinical_data, by = "Sample.ID")
  } else {
    # PAM
    TCGA_pam = runPAM_single_algorithm(algorithm_name = algorithm,
                                       train.expr = env_objects[[algorithm]]$plotdata$RNAseq,
                                       moic.res   = env_objects[[algorithm]]$plot_object,
                                       test.expr  = z_transcr)
    df = as.data.frame(TCGA_pam$clust.res) %>%
      dplyr::select(samID, everything())
    colnames(df) = c("Sample.ID", algorithm)
    label_maps[[algorithm]] = df %>%
      inner_join(holdout_clinical_data, by = "Sample.ID")
  }
  
  rm(list=setdiff(ls(), c("home", "algorithm", "algorithms", "z_transcr", "train_samples",
                          "holdout_clinical_data", "expdat",
                          "log2_filt", "label_maps", "env_objects", "t_start")))
  cat(blue("Done with", algorithm, "\n"))
  invisible(gc())
}

# Clean cache
save.image(paste0(home, "/Results/Survival_evaluations/Surv.RData"))
gc()

# Export
write.xlsx(holdout_clinical_data, "Resources/TCGA/Surv_clinical_data.xlsx",
                     overwrite = TRUE)

# Survival analysis #####
# Load custom helper functions
source("Scripts/automated_scripts/custom_functions.R")
source("Scripts/automated_scripts/modified_MOVICS_functions.R")

# Coloring schemes
cluster_colors = c("#2EC4B6", "#E71D36", "#FF9F1C", "#BDD5EA", 
                   "#FFA5AB", "#011627", "#023E8A", "#9D4EDD")

# Set up survival analysis
surv = list()

# Import survival data from cBioBortal
cBioPortal = read.xlsx("Resources/cBioPortal_surv.xlsx")

# If we use XENA USC data

# xena <- read.table("Resources/Xena_USC_survival_BRCA.txt",
#                    header = TRUE, sep = "\t") %>%
#   rename(Patient.ID = X_PATIENT, OS_DAYS = OS.time, OS_STATUS = OS) %>%
#   filter(Redaction != "Redacted") %>%
#   arrange(Patient.ID, desc(OS_DAYS)) %>%
#   distinct(Patient.ID, .keep_all = TRUE) %>%
#   mutate(vital_status = ifelse(OS_STATUS == 0, "Alive", "Dead"),
#          OS_MONTHS = OS_DAYS/30)

survival_data = cBioPortal

# Set correct_survival to NULL, "BH", or "Bonferroni" before the loop
correct_survival <- "BH"  # Or NULL or "Bonferroni"

dashfile <- file.path(home, "Results", "Survival_evaluations", "survival_checks.txt")
unlink(dashfile)

for (algorithm in algorithms) {
  
  # Get the corresponding label-mapping data on the holdout TCGA data
  surv_df = label_maps[[algorithm]] %>%
    mutate(Source = "Holdout") %>%
    #dplyr::select(Sample.ID, !!sym(algorithm), Patient.ID) %>%
    dplyr::select(Sample.ID, !!sym(algorithm), Patient.ID, everything()) %>%
    dplyr::select(-vital_status, -days_to_death, -days_to_last_known_alive) %>%
    inner_join(survival_data, by = "Patient.ID") %>%
    dplyr::rename(LN_status = primary_lymph_node_presentation_assessment,
                  ER_status = breast_carcinoma_estrogen_receptor_status,
                  HER2_status = lab_proc_her2_neu_immunohistochemistry_receptor_status,
                  # stage = stage_event_pathologic_stage,
                  dist_metastasis = distant_metastasis_present_ind2) %>%
    mutate(age = days_to_birth/365)
  
  surv_df = surv_df[!is.na(surv_df[, algorithm]), ]
  surv_df$OS_MONTHS = as.numeric(surv_df$OS_MONTHS)
  surv_df$OS_DAYS = as.numeric(surv_df$OS_DAYS)
  
  # Create a format expected by TCGAanalyze_survival()
  surv_df$days_to_death = ifelse(surv_df$vital_status == "Dead", surv_df$OS_DAYS,
                                 surv_df$days_to_last_followup)
  
  # Convert grouping variable to factor
  surv_df[, algorithm] = as.factor(surv_df[, algorithm])
  
  # Capture warnings, Run TCGAanalyze_survival
  surv[[algorithm]] <- TCGAanalyze_survival_custom3(
        data         = surv_df,
        clusterCol   = algorithm,
        adjustVars   = c("age"),
        main         = bquote( bold(.(algorithm) ~ "on TCGA-BRCA holdout set") ),
        title.size   = 18,
        xlab         = expression(bold("Time since diagnosis (days)")),
        ylab         = expression(bold("Survival probability")),
        color        = cluster_colors[1:length(unique(surv_df[, algorithm]))],
        legend       = expression(bold("Legend")),
        save.filename= paste0(home, "/Results/Survival_evaluations/", algorithm, "/",
                              algorithm, "_survival_plot.pdf"),
        save.width   = 10,
        save.height  = 10*seq(1, 1.25, length.out = 9)[length(unique(surv_df[, algorithm]))-1],
        save.dpi     = 700,
        ph_threshold = 0.05, # Set Proportional Hazard threshold
        vif_cutoff = 5, # Set VIF threshold
        summary.filename = dashfile
  )
  
  # Skip if results for TCGAanalyze_survival_custom3 where not calculated
  if (is.null(surv[[algorithm]]$pvalue)) {
    # Exporting
    write.xlsx(surv_df,
               paste0(home, "/Results/Survival_evaluations/", algorithm, "/",
                      algorithm, "_surv_data.xlsx"),
               overwrite = TRUE)
    
    # Save in R
    surv[[algorithm]][["df"]] = surv_df
    next
  }
  
  # Perform p-value correction if requested
  if (!is.null(correct_survival) && !is.na(surv[[algorithm]]$pvalue)) {
    if (correct_survival == "BH") {
      surv[[algorithm]]$pvalue_adjusted <- p.adjust(surv[[algorithm]]$pvalue,
                                                    method = "BH",
                                                    n = length(algorithms))
    } else if (correct_survival == "Bonferroni") {
      surv[[algorithm]]$pvalue_adjusted <- p.adjust(surv[[algorithm]]$pvalue,
                                                    method = "bonferroni",
                                                    n = length(algorithms))
    }
  }
  
  # Exporting
  write.xlsx(surv_df,
             paste0(home, "/Results/Survival_evaluations/", algorithm, "/",
                    algorithm, "_surv_data.xlsx"),
             overwrite = TRUE)
  
  # Save in R
  surv[[algorithm]][["df"]] = surv_df
}

rm(algorithm, surv_df, i, newdir); gc()

# Clinical comparisons #####
# The following parameters were extracted using the 
# Scripts/Survival analysis/Extract_pdf_parameters_for_clincomp.py script
pdf_parameters = read.table("Resources/extracted_parameters.txt", header = TRUE,
                            sep = "\t", quote = "")
# Convert the string representations to R objects:
pdf_parameters$pdf_level_col_width <- sapply(pdf_parameters$pdf_level_col_width, 
                                             function(x) eval(parse(text = x)), 
                                             simplify = FALSE)

pdf_parameters$pdf_count_col_width <- sapply(pdf_parameters$pdf_count_col_width,
                                             function(x) eval(parse(text = x)))

pdf_parameters$pdf_pval_col_width <- sapply(pdf_parameters$pdf_pval_col_width,
                                            function(x) eval(parse(text = x)))

pdf_parameters$pdf_test_col_width <- sapply(pdf_parameters$pdf_test_col_width,
                                            function(x) eval(parse(text = x)))

# Import styles and helpful objects
scheme = readRDS("Resources/scheme.rds")
cols_of_interest = colnames(scheme$var2comp)[2:ncol(scheme$var2comp)]

clincomps = list()
# Loop for clinical comparisons
for (algorithm in algorithms) {
  var2comp = holdout_clinical_data[, c("Sample.ID", cols_of_interest)] %>%
    inner_join(surv[[algorithm]]$df[, c("Sample.ID", algorithm)], by = "Sample.ID") %>%
    tibble::column_to_rownames(var = "Sample.ID")
  
  cols_to_numeric <- c("days_to_birth", "days_to_death",
                       "days_to_last_known_alive", 
                       "days_to_last_followup",
                       "age_at_initial_pathologic_diagnosis",
                       "er_level_cell_percentage_category",
                       "progesterone_receptor_level_cell_percent_category",
                       "number_of_lymphnodes_positive_by_ihc",
                       "number_of_lymphnodes_positive_by_he")
  
  # Convert the specified columns to numeric
  var2comp[, cols_to_numeric] <- lapply(var2comp[, cols_to_numeric, drop = FALSE],
                                        function(x) as.numeric(x))
  
  # Remove unknown levels for statistical tests
  var2comp_nonas = var2comp
  for (i in 1:ncol(var2comp)) {
    nas = which(var2comp[, i] == "Unknown")
    var2comp_nonas[nas, i] = NA
    empties = which(var2comp[, i] == "")
    var2comp_nonas[empties, i] = NA
  }
  rm(nas, empties); gc()
  
  moic.res = list()
  moic.res$clust.res = label_maps[[algorithm]][, c("Sample.ID", algorithm)]
  colnames(moic.res$clust.res) = c("samID", "clust")
  rownames(moic.res$clust.res) = moic.res$clust.res$samID
  
  # The function below expects clusters to not have the algorithm label, but just
  # numeric labels
  moic.res$clust.res$clust = gsub(algorithm, "", moic.res$clust.res$clust)
  
  # Statistical comparisons
  clin_comp = compClinvar_single_algorithm(algorithm_name = algorithm,
                                           moic.res = moic.res,
                                           var2comp = var2comp_nonas,
                                           strata = algorithm,
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
                                           res.path = paste0(home, "/Results/Survival_evaluations/", algorithm),
                                           output_pdf = TRUE,
                                           pdf_level_col_width = pdf_parameters$pdf_level_col_width[pdf_parameters$Algorithm == algorithm][[1]],
                                           pdf_count_col_width = pdf_parameters$pdf_count_col_width[pdf_parameters$Algorithm == algorithm],
                                           pdf_pval_col_width = pdf_parameters$pdf_pval_col_width[pdf_parameters$Algorithm == algorithm],
                                           pdf_test_col_width = pdf_parameters$pdf_test_col_width[pdf_parameters$Algorithm == algorithm],
                                           pdf_tab_font_size = pdf_parameters$pdf_tab_font_size[pdf_parameters$Algorithm == algorithm])
  
  clin_ordinal_comp = compClinvar_ordinal_single_algorithm(algorithm_name = algorithm,
                                                           moic.res = moic.res,
                                                           var2comp = var2comp_nonas[, c("number_of_lymphnodes_positive_by_ihc",
                                                                           "number_of_lymphnodes_positive_by_he",
                                                                           algorithm)],
                                                           strata = algorithm,
                                                           ordinalVars = c("number_of_lymphnodes_positive_by_ihc",
                                                                           "number_of_lymphnodes_positive_by_he"),
                                                           includeNA = FALSE,
                                                           tab.name = "Summary of ordinal clinical variables",
                                                           res.path = paste0(home, "/Results/Survival_evaluations/", algorithm),
                                                           output_pdf = TRUE,
                                                           pdf_template_loc = paste0(home, "/Scripts/automated_scripts/clincomp_template.Rmd"),
                                                           pdf_level_col_width = pdf_parameters$pdf_level_col_width[pdf_parameters$Algorithm == algorithm][[1]],
                                                           pdf_count_col_width = pdf_parameters$pdf_count_col_width[pdf_parameters$Algorithm == algorithm],
                                                           pdf_pval_col_width = pdf_parameters$pdf_pval_col_width[pdf_parameters$Algorithm == algorithm],
                                                           pdf_test_col_width = pdf_parameters$pdf_test_col_width[pdf_parameters$Algorithm == algorithm],
                                                           pdf_tab_font_size = pdf_parameters$pdf_tab_font_size[pdf_parameters$Algorithm == algorithm])
  
  # Save
  clincomps[[algorithm]] = list(var2comp = var2comp,
                                var2comp_nonas = var2comp_nonas,
                                clincomp = clin_comp,
                                ordinal = clin_ordinal_comp)
}

rm(algorithm, moic.res); gc()

# annCol = scheme$annCol
# annColors = scheme$annColors

# Bar charts with variables of interest #####



# Save environment
save.image(paste0(home, "/Results/Survival_evaluations/Surv.RData"))

# Export session info
writeLines(capture.output(sessionInfo()), 
           paste0(home, "/Results/Survival_evaluations/sessionInfo.txt"))