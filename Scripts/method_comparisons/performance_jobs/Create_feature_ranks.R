library(data.table)

# This RDS object was produced using the Scripts/MOVICS/MOVICS_baseline.R script
input = readRDS("Resources/TCGA/mm_input.rds")

# Ensure reproducibility
RNGversion("4.2.2")
set.seed(123)

# Create feature subsets
out_dir <- "Resources/Performance/Feature_perturbations"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_dir, "feature_rankings.tsv.gz")

rank_mad_ids <- function(mat) {
  stopifnot(!is.null(rownames(mat)))
  m <- apply(mat, 1, function(x) mad(x, na.rm = TRUE))
  keep <- which(!is.na(m))
  ord <- keep[order(m[keep], decreasing = TRUE)]
  rownames(mat)[ord]
}

rows <- list()

# Continuous modalities by MAD
cont_modalities <- c("RNAseq", "CNV", "Methylation", "miRNA")
for (nm in cont_modalities) {
  if (!is.null(input[[nm]])) {
    ids <- rank_mad_ids(input[[nm]])
    rows[[nm]] <- data.table(modality = nm, rank = seq_along(ids), feature_id = ids)
  }
}

# Binary modality (SNPs): Rank by variance p(1-p)
if (!is.null(input[["SNPs"]])) {
  snp <- input[["SNPs"]]
  
  # prevalence of non-zero / mutated
  p <- rowMeans(snp != 0, na.rm = TRUE)
  
  # variance for Bernoulli variable
  v <- p * (1 - p)
  
  # deterministic ordering on ties
  ord <- order(v, p, decreasing = TRUE, rownames(snp))
  
  ids <- rownames(snp)[ord]
  rows[["SNPs"]] <- data.table(
    modality = "SNPs",
    rank = seq_along(ids),
    feature_id = ids
  )
}

feat_rank <- rbindlist(rows, use.names = TRUE, fill = TRUE)

# Write compressed TSV (portable across R/Python)
fwrite(
  feat_rank, out_file,
  sep = "\t", quote = FALSE, na = "NA",
  compress = "gzip"
)
