# Exports stratified sample subsets (per fraction and replicate) as a single TSV.GZ.
library(data.table)
library(openxlsx)
library(dplyr)

# This RDS object was produced using the Scripts/MOVICS/MOVICS_baseline.R script
input <- readRDS("Resources/TCGA/mm_input.rds")

# Import clinical data
clinical_data <- openxlsx::read.xlsx("Resources/TCGA/clinical_data.xlsx") %>%
  dplyr::rename(ER_status = breast_carcinoma_estrogen_receptor_status)
unknown = which(clinical_data$ER_status == "")
clinical_data$ER_status[unknown] = "Unknown"

# Reproducibility
RNGversion("4.2.2")
set.seed(123)

out_dir <- "Resources/Performance/Sample_perturbations"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_dir, "sample_subsets.tsv.gz")

# Map ER status to sample IDs and keep only samples that exist in the input
sid <- clinical_data[["Sample.ID"]]
er  <- clinical_data[["ER_status"]]
names(er) <- sid

# Use RNAseq columns as the "universe" of valid sample IDs
ref_mod <- if (!is.null(input[["RNAseq"]])) "RNAseq" else names(input)[which(!vapply(input, is.null, logical(1)))[1]]
valid_ids <- colnames(input[[ref_mod]])

er <- er[names(er) %in% valid_ids]
er <- er[!is.na(er)]
er <- droplevels(factor(er))

fractions <- c(10, 20, 50, 70, 90)
B <- 10

strat_sample_ids <- function(frac, B, er_factor_named) {
  ids_all <- names(er_factor_named)
  levs <- levels(er_factor_named)
  out <- vector("list", B)
  
  for (b in seq_len(B)) {
    picked <- character(0)
    for (lv in levs) {
      ids_lv <- ids_all[er_factor_named == lv]
      n_lv <- max(1L, round(length(ids_lv) * frac / 100))
      n_lv <- min(n_lv, length(ids_lv))
      picked <- c(picked, sample(ids_lv, size = n_lv, replace = FALSE))
    }
    out[[b]] <- picked
  }
  out
}

all_rows <- list()
k <- 1L
for (frac in fractions) {
  ids_list <- strat_sample_ids(frac, B, er)
  for (b in seq_len(B)) {
    all_rows[[k]] <- data.table(
      fraction = frac,
      replicate = b,
      sample_id = ids_list[[b]]
    )
    k <- k + 1L
  }
}

sample_subsets <- rbindlist(all_rows, use.names = TRUE)

# Write compressed TSV (portable across R/Python)
fwrite(
  sample_subsets, out_file,
  sep = "\t", quote = FALSE, na = "NA",
  compress = "gzip"
)
