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

dashfile <- file.path(home, "Results", "Survival_evaluations", "training_survival_checks.txt")
unlink(dashfile)

for (algorithm in algorithms) {
  
  alg_clusterings_train = read.xlsx(paste0("Results/single_algorithm/",
                                           algorithm, "/", algorithm,
                                           "_TCGA_RNAseq-CNV-Methylation-miRNA-SNPs_eval_on_transNEO_clusterings.xlsx"))
  
  # Get the corresponding label-mapping data on the holdout TCGA data
  surv_df = train_clinical_data %>%
    dplyr::select(Sample.ID, Patient.ID, days_to_last_followup, everything()) %>%
    inner_join(alg_clusterings_train, by = "Sample.ID") %>%
    mutate(Source = "Train") %>%
    dplyr::select(Sample.ID, !!sym(algorithm) := Cluster, Patient.ID,
                  everything()) %>%
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
        main         = bquote( bold(.(algorithm) ~ "on TCGA-BRCA training set") ),
        title.size   = 18,
        xlab         = expression(bold("Time since diagnosis (days)")),
        ylab         = expression(bold("Survival probability")),
        color        = cluster_colors[1:length(unique(surv_df[, algorithm]))],
        legend       = expression(bold("Legend")),
        save.filename= paste0(home, "/Results/Survival_evaluations/", algorithm, "/",
                              algorithm,
                              "_survival_plot_training_data.pdf"),
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
                      algorithm, "_surv_TRAINING_data.xlsx"),
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
                    algorithm, "_surv_TRAINING_data.xlsx"),
             overwrite = TRUE)
  
  # Save in R
  surv[[algorithm]][["df"]] = surv_df
}

rm(algorithm, surv_df); gc()

# Save environment
save.image(paste0(home, "/Results/Survival_evaluations/Surv_training.RData"))

# Export session info
writeLines(capture.output(sessionInfo()), 
           paste0(home, "/Results/Survival_evaluations/sessionInfo_training.txt"))