library(TCGAbiolinks)
library(dplyr)
library(survminer)
library(survival)
library(openxlsx)

# Ensure reproducibility
RNGversion("4.2.2")
set.seed(123)

# Load custom helper functions
source("Scripts/automated_scripts/custom_functions.R")
source("Scripts/automated_scripts/modified_MOVICS_functions.R")

# List of algorithms
algorithms = c("ab-SNF", "ANF", "CIMLR", "COCA", "iClusterBayes", "KLIC",
               "LRAcluster", "MDICC", "MFA",
               "NEMO", "RWR-F", "RWR-NF", "SNF", "Spectrum", "wMKL",
               "MONET", "MSNE", "MOFA")

# Set home directory
home = getwd()

# Import coloring schemes etc
scheme = readRDS("Resources/scheme.rds")
cluster_colors = c("#2EC4B6", "#E71D36", "#FF9F1C", "#BDD5EA", 
                   "#FFA5AB", "#011627", "#023E8A", "#9D4EDD")

# Set up survival analysis
train_clinical_data = read.xlsx("Resources/TCGA/clinical_data.xlsx")
surv = list()

# Import survival data from cBioBortal
cBioPortal = read.xlsx("Resources/cBioPortal_surv.xlsx")
for (algorithm in algorithms) {
  alg_clusterings_train = read.xlsx(paste0("Results/single_algorithm/",
                                           algorithm, "/", algorithm,
                                           "_TCGA_RNAseq-CNV-Methylation-miRNA-SNPs_eval_on_transNEO_clusterings.xlsx"))
  
  surv_df = train_clinical_data %>%
                    dplyr::select(Sample.ID, Patient.ID, vital_status,
                                  days_to_death, days_to_last_followup) %>%
                    inner_join(alg_clusterings_train, by = "Sample.ID") %>%
                    mutate(Source = "Train") %>%
                    dplyr::select(Sample.ID, !!sym(algorithm) := Cluster, Patient.ID, 
                                  Source) %>%
    inner_join(cBioPortal, by = "Patient.ID")
  
  surv_df = surv_df[!is.na(surv_df[, algorithm]), ]
  surv_df$OS_MONTHS = as.numeric(surv_df$OS_MONTHS)
  surv_df$OS_DAYS = as.numeric(surv_df$OS_DAYS)
  
  # Create a format expected by TCGAanalyze_survival()
  surv_df$days_to_death = ifelse(surv_df$vital_status == "Dead", surv_df$OS_DAYS,
                                 NA)
  surv_df$days_to_last_follow_up = surv_df$OS_DAYS
  
  # Convert grouping variable to factor
  surv_df[, algorithm] = as.factor(surv_df[, algorithm])
  
  # Run TCGAanalyze_survival
  surv[[algorithm]] = TCGAanalyze_survival_custom(
    data = surv_df,
    clusterCol = algorithm,
    main = paste(algorithm, "on TCGA-BRCA training set: survival analysis"),
    ylab = expression(bold("Survival probability")),
    xlab = expression(bold("Time since diagnosis (days)")),
    filename = paste0(home, "/Results/Survival_evaluations/", algorithm, "/",
                      algorithm,
                      "_survival_plot_training_data.pdf"),
    legend = expression(bold("Legend")),
    risk.table.height = 0.2*seq(1, 1.25, length.out = 9)[length(unique(surv_df[, algorithm]))-1],
    color = cluster_colors[1:length(unique(surv_df[, algorithm]))],
    main_fontsize = 18,
    height = 10*seq(1, 1.25, length.out = 9)[length(unique(surv_df[, algorithm]))-1],
    width = 10,
    dpi = 700
  )
  
  # Exporting
  write.xlsx(surv_df,
             paste0(home, "/Results/Survival_evaluations/", algorithm, "/",
                    algorithm, "_surv_TRAINING_data.xlsx"),
             overwrite = TRUE)
  
  # Save in R
  surv[[algorithm]][["df"]] = surv_df
}

rm(algorithm, surv_df); gc()

# Clinical comparisons #####

# Bar charts with variables of interest #####


# Save environment
save.image(paste0(home, "/Results/Survival_evaluations/Surv_training.RData"))

# Export session info
writeLines(capture.output(sessionInfo()), 
           paste0(home, "/Results/Survival_evaluations/sessionInfo_training.txt"))