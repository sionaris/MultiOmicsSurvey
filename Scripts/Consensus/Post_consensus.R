# Libraries
library(openxlsx)
library(dplyr)
library(ggplot2)

# Import the full environment for consensus but only keep plot_object and plotdata
load("Results/Consensus/CC/CC_TCGA_RNAseq-CNV-Methylation-miRNA-SNPs_eval_on_transNEO_env.RData")
rm(list=setdiff(ls(), c("plotdata", "plot_object"))); gc()

# Ensure reproducibility
RNGversion("4.2.2")
set.seed(123)

# Load custom helper functions
source("Scripts/automated_scripts/custom_functions.R")
source("Scripts/automated_scripts/modified_MOVICS_functions.R")

# Preamble
home = getwd()
data_source = "TCGA" # e.g. TCGA, TCGA-transNEO, transNEO-PARTNER
data_types = "RNAseq-CNV-Methylation-miRNA-SNPs" # e.g. RNAseq, RNAseq-CNV-miRNA
evaluation_source = "transNEO" # e.g. PARTNER, transNEO-PARTNER 

# Import clusterings
R_algorithms = c("ab-SNF", "ANF", "CIMLR", "COCA", "iClusterBayes", "KLIC",
                 "LRAcluster", "MDICC", "MFA", # "mixKernel", #"MOFA", 
                 "NEMO", "RWR-F", "RWR-NF", "SNF", "Spectrum", "wMKL")
Python_algorithms = c("MONET", "MSNE", "MOFA") #, "PAMOGK",)
algorithms = c(R_algorithms, Python_algorithms)
algorithm_languages = c(rep("R", length(R_algorithms)),
                        rep("Python", length(Python_algorithms)))
names(algorithm_languages) = algorithms
algorithm_languages["MOFA"] = "R & Python"
algorithm_languages["MDICC"] = "R & Python"
algorithm_languages["MixKernel"] = "R & Python"

clusterings = list()

# R methods
for (R_algorithm in R_algorithms) {
  res_dir = paste0(home, "/Results/single_algorithm/", R_algorithm, "/")
  res_dir_files = list.files(path = res_dir, pattern = ".*_clusterings\\.xlsx$", 
                             full.names = TRUE)
  if (length(res_dir_files == 1)) {
    clusterings[[R_algorithm]] = read.xlsx(res_dir_files[1])
  } else {
    clusterings[[R_algorithm]] = NA
  }
}

# Python methods
for (Python_algorithm in Python_algorithms) {
  res_dir = paste0(home, "/Results/single_algorithm/", Python_algorithm, "/")
  res_dir_files = list.files(path = res_dir, pattern = ".*_clusterings\\.xlsx$", 
                             full.names = TRUE)
  if (length(res_dir_files == 1)) {
    clusterings[[Python_algorithm]] = read.xlsx(res_dir_files[1])
  } else {
    clusterings[[Python_algorithm]] = NA
  }
}
names(clusterings) = algorithms

# Set up method categories
similarity_network_methods = c("ab-SNF", "ANF", "MDICC", "MSNE", "NEMO", "RWR-F", "RWR-NF", "SNF", "Spectrum")
multiple_kernel_learning = c("CIMLR", "KLIC", "wMKL") # "mixKernel"
matrix_factorization_latent_variables = c("MFA", "MOFA", "LRAcluster")
graph_methods = c("MONET") #, "PAMOGK")
bayesian = c("iClusterBayes")
# cca_methods = c("RGCCA", "SGCCA")
# low_rank_methods = c("LRAcluster") #, moCluster, PINSPlus
cc_ensemble = c("COCA")

# Primary annotation
primary_annotation_rag = c(rep("Similarity Network", length(similarity_network_methods)),
                           rep("Multiple Kernel Learning", length(multiple_kernel_learning)),
                           rep("Matrix Factorization/Latent Variables", length(matrix_factorization_latent_variables)),
                           rep("Graph-based Methods", length(graph_methods)),
                           rep("Bayesian", length(bayesian)),
                           rep("Consensus/Ensemble Clustering", length(cc_ensemble)))
# rep("Canonical Correlation", length(cca_methods)),
# rep("Low-rank Projection", length(low_rank_methods)),
# rep("Miscellaneous", length(misc)))
names(primary_annotation_rag) = c(similarity_network_methods, multiple_kernel_learning,
                                  matrix_factorization_latent_variables, graph_methods, bayesian,
                                  cc_ensemble)
# cca_methods, low_rank_methods, 
# misc)

clusterings = clusterings[names(primary_annotation_rag)]
algorithm_languages = algorithm_languages[names(primary_annotation_rag)]

# # Secondary annotation
# secondary_annotation_rag = c(rep("Graph-based methods", 8),
#                              rep("Low-rank Projection", 7),
#                              rep("Miscellaneous", 7))
# names(secondary_annotation_rag) = c(similarity_network_methods, "Spectrum",
#                                     cca_methods, bayesian, matrix_factorization_latent_variables,
#                                     multiple_kernel_learning, graph_methods,
#                                     low_rank_methods, "COCA")

# Convert individual cluster labels from "algorithm#" to just #
generic_clusterings = clusterings
for (i in 1:length(generic_clusterings)) {
  if (!is.null(ncol(generic_clusterings[[i]]))) { # temporary error control
    generic_clusterings[[i]]$Cluster <- gsub("[^0-9]", "", generic_clusterings[[i]]$Cluster)
  }
}

rm(res_dir, res_dir_files, R_algorithm, Python_algorithm); gc()


# ARI agreement #####
library(mclust)

# Initialize a matrix to store ARI values
ari_vector <- setNames(rep(NA, length(algorithms)), algorithms)

# Import the consensus clustering
consensus = read.xlsx("Results/Consensus/CC/CC_TCGA_RNAseq-CNV-Methylation-miRNA-SNPs_eval_on_transNEO_clusterings.xlsx") %>%
  dplyr::mutate(Cluster = gsub("CC", "", Cluster))

# Populate the ari_vector
for (algorithm in algorithms) {
  ari_vector[[algorithm]] = adjustedRandIndex(generic_clusterings[[algorithm]]$Cluster,
                          consensus$Cluster)
}

# Print ARI vector
print(ari_vector)

# Create data frame
ARI_df = as.data.frame(list(algorithm = names(ari_vector),
                            ARI = ari_vector)) %>%
  inner_join(as.data.frame(list(algorithm = names(primary_annotation_rag),
                                Category = primary_annotation_rag)),
             by = "algorithm") %>%
  arrange(Category, algorithm) %>%
  mutate(algorithm = factor(algorithm, levels = unique(algorithm)),
         annot_x = -0.035)

# Bar chart for ARI
library(ggnewscale)
library(rcartocolor)
category_colors <- c(
  "Similarity Network" = carto_pal("Bold", n = 12)[1],
  "Multiple Kernel Learning" = carto_pal("Bold", n = 12)[2],
  "Matrix Factorization/Latent Variables" = carto_pal("Antique", n = 12)[5],
  "Graph-based Methods" = carto_pal("Bold", n = 12)[4],
  "Bayesian" = carto_pal("Bold", n = 12)[11],
  "Consensus/Ensemble Clustering" = carto_pal("Bold", n = 12)[9]
  # "Canonical Correlation" = carto_pal("Bold", n = 12)[9],
  # "Low-rank Projection" = carto_pal("Bold", n = 12)[10]
)

p <- ggplot(ARI_df, aes(y = algorithm)) +
  # Category annotation tile with legend
  geom_tile(aes(x = annot_x, fill = Category), width = 0.015, 
            height = 0.8, color = "grey50", linewidth = 0.2) +
  scale_fill_manual(values = category_colors, name = "Category") +
  new_scale_fill() +
  # ARI bars with colorbar legend
  geom_col(aes(x = ARI, fill = ARI), color = NA) +
  scale_fill_carto_c(palette = "Teal", name = "ARI") +
  coord_cartesian(xlim = c(-0.05, 1)) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
  labs(
    title = "ARI between consensus and individual methods",
    x = "ARI",
    y = NULL
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 7),
    axis.title.x = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )

# Save plot
ggsave(file.path(home, "Results/Consensus/Post/ARI_to_conensus_bar_chart.png"), p, dpi = 700, width = 8, height = 6)
ggsave(file.path(home, "Results/Consensus/Post/ARI_to_conensus_bar_chart.pdf"), p, dpi = 700, width = 8, height = 6)

# NTP #####
dgea.marker.up_1000 <- runMarker_single_algorithm_no_export(algorithm_name = "CC",
                                                              moic.res = plot_object,
                                                              n.marker = 1000,
                                                              dea.method    = "limma", # name of DEA method
                                                              prefix        = "dgea_", # MUST be the same of argument in runDEA()
                                                              dat.path      = paste0(home, "/Results/Consensus/CC"), # path of DEA files
                                                              p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                                                              p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                                                              norm.expr = plotdata$RNAseq,
                                                              dirct         = "up" # direction of dysregulation in expression
)

# If the directory of the algorithm does not exist, create it
newdir = paste0(home, "/Results/Consensus/Post")
if (!dir.exists(newdir)) {
  dir.create(newdir)
}

# NTP
library(MOVICS)
z_transcr = readRDS("Resources/TCGA/Surv_standardized_transcr.rds")
holdout_clinical_data = read.xlsx("Resources/TCGA/Surv_clinical_data.xlsx")

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
  fig.path = paste0(home, "/Results/Consensus/Post"),
  fig.name = paste0("CC_ntp_expr_up_heatmap_TCGA"))

holdout_df = as.data.frame(TCGA_ntp_expr_up$clust.res) %>%
  dplyr::select(samID, clust)
colnames(holdout_df) = c("Sample.ID", "CC")
holdout_df = holdout_df %>%
  inner_join(holdout_clinical_data, by = "Sample.ID")

# Consensus and Survival #####
# Coloring schemes
cluster_colors = c("#2EC4B6", "#E71D36")

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

# Survival in the holdout ###

# Get the corresponding label-mapping data on the holdout TCGA data
surv_df = holdout_df %>%
  mutate(Source = "Holdout") %>%
  #dplyr::select(Sample.ID, CC, Patient.ID) %>%
  dplyr::select(Sample.ID, CC, Patient.ID, everything()) %>%
  dplyr::select(-vital_status) %>%
  inner_join(survival_data, by = "Patient.ID") %>%
  dplyr::rename(LN_status = primary_lymph_node_presentation_assessment,
                ER_status = breast_carcinoma_estrogen_receptor_status,
                HER2_status = lab_proc_her2_neu_immunohistochemistry_receptor_status,
                # stage = stage_event_pathologic_stage,
                dist_metastasis = distant_metastasis_present_ind2) %>%
  mutate(age = days_to_birth/365)

surv_df = surv_df[!is.na(surv_df[, "CC"]), ]
surv_df$OS_MONTHS = as.numeric(surv_df$OS_MONTHS)
surv_df$OS_DAYS = as.numeric(surv_df$OS_DAYS)

# Create a format expected by TCGAanalyze_survival()
surv_df$days_to_death = ifelse(surv_df$vital_status == "Dead", surv_df$OS_DAYS,
                               NA)
surv_df$days_to_last_follow_up = surv_df$OS_DAYS

# Convert grouping variable to factor
surv_df[, "CC"] = as.factor(surv_df[, "CC"])

# Capture warnings, Run TCGAanalyze_survival
library(survival)
withCallingHandlers({
  surv <- tryCatch({
    TCGAanalyze_survival_custom3(
      data         = surv_df,
      clusterCol   = "CC",
      adjustVars   = c("age", "histological_type", "LN_status", "menopause_status",
                       "ER_status", "HER2_status", "dist_metastasis"),
      main         = paste("CC on TCGA-BRCA holdout set"),
      title.size   = 18,
      xlab         = expression(bold("Time since diagnosis (days)")),
      ylab         = expression(bold("Survival probability")),
      color        = cluster_colors[1:length(unique(surv_df[, "CC"]))],
      legend       = expression(bold("Legend")),
      save.filename= paste0(home, "/Results/Consensus/Post/CC_survival_plot.pdf"),
      save.width   = 10,
      save.height  = 10*seq(1, 1.25, length.out = 9)[length(unique(surv_df[, "CC"]))-1],
      save.dpi     = 700,
      ph_threshold = 0.05, # Set Proportional Hazard threshold
      vif_cutoff = 5 # Set VIF threshold
    )
  }, error = function(e) {
    # If an error occurs, return a list with error information
    list(error = e$message, pvalue = NA) # important to add pvalue = NA so code does not crash
  })
}, warning = function(w) {
  # Capture warnings here
  holdout_warnings <- conditionMessage(w)
  invokeRestart("muffleWarning")
})

# Store convergence info from warnings
if (is.list(surv) && !is.null(surv$convergence)) {
  surv[["convergence"]] <- list(
    singular = any(grepl("computationally singular", holdout_warnings)),
    loglik_failed = any(grepl("Loglik converged before", holdout_warnings))
  )
} else {
  surv[["convergence"]] <- list(
    singular = NA,
    loglik_failed = NA
  )
}

# Exporting
write.xlsx(surv_df,
           paste0(home, "/Results/Consensus/Post/CC_surv_data.xlsx"),
           overwrite = TRUE)

# Save in R
surv[["df"]] = surv_df

# Survival in training set ###
train_clinical_data = read.xlsx("Resources/TCGA/clinical_data.xlsx")

# Get the corresponding label-mapping data on the training TCGA data
surv_training_df = holdout_df %>%
  mutate(Source = "Holdout") %>%
  #dplyr::select(Sample.ID, CC, Patient.ID) %>%
  dplyr::select(Sample.ID, CC, Patient.ID, everything()) %>%
  dplyr::select(-vital_status) %>%
  inner_join(survival_data, by = "Patient.ID") %>%
  dplyr::rename(LN_status = primary_lymph_node_presentation_assessment,
                ER_status = breast_carcinoma_estrogen_receptor_status,
                HER2_status = lab_proc_her2_neu_immunohistochemistry_receptor_status,
                # stage = stage_event_pathologic_stage,
                dist_metastasis = distant_metastasis_present_ind2) %>%
  mutate(age = days_to_birth/365)

surv_training_df = surv_training_df[!is.na(surv_training_df[, "CC"]), ]
surv_training_df$OS_MONTHS = as.numeric(surv_training_df$OS_MONTHS)
surv_training_df$OS_DAYS = as.numeric(surv_training_df$OS_DAYS)

# Create a format expected by TCGAanalyze_survival()
surv_training_df$days_to_death = ifelse(surv_training_df$vital_status == "Dead", surv_training_df$OS_DAYS,
                                        NA)
surv_training_df$days_to_last_follow_up = surv_training_df$OS_DAYS

# Convert grouping variable to factor
surv_training_df[, "CC"] = as.factor(surv_training_df[, "CC"])

# Capture warnings, Run TCGAanalyze_survival
withCallingHandlers({
  surv_training <- tryCatch({
    TCGAanalyze_survival_custom3(
      data         = surv_training_df,
      clusterCol   = "CC",
      adjustVars   = c("age", "histological_type", "LN_status", "menopause_status",
                       "ER_status", "HER2_status", "dist_metastasis"),
      main         = paste("CC on TCGA-BRCA holdout set"),
      title.size   = 18,
      xlab         = expression(bold("Time since diagnosis (days)")),
      ylab         = expression(bold("survival probability")),
      color        = cluster_colors[1:length(unique(surv_training_df[, "CC"]))],
      legend       = expression(bold("Legend")),
      save.filename= paste0(home, "/Results/Consensus/Post/CC_training_survival_plot.pdf"),
      save.width   = 10,
      save.height  = 10*seq(1, 1.25, length.out = 9)[length(unique(surv_training_df[, "CC"]))-1],
      save.dpi     = 700,
      ph_threshold = 0.05, # Set Proportional Hazard threshold
      vif_cutoff = 5 # Set VIF threshold
    )
  }, error = function(e) {
    # If an error occurs, return a list with error information
    list(error = e$message, pvalue = NA) # important to add pvalue = NA so code does not crash
  })
}, warning = function(w) {
  # Capture warnings here
  training_warnings <- conditionMessage(w)
  invokeRestart("muffleWarning")
})

# Store convergence info from warnings
if (is.list(surv_training) && !is.null(surv_training$convergence)) {
  surv_training[["convergence"]] <- list(
    singular = any(grepl("computationally singular", training_warnings)),
    loglik_failed = any(grepl("Loglik converged before", training_warnings))
  )
} else {
  surv_training[["convergence"]] <- list(
    singular = NA,
    loglik_failed = NA
  )
}

# Exporting
write.xlsx(surv_training_df,
           paste0(home, "/Results/Consensus/Post/CC_surv_training_data.xlsx"),
           overwrite = TRUE)

# Save in R
surv_training[["df"]] = surv_training_df

# Save environment
save.image(paste0(home, "/Results/Consensus/Post/CC_Post.RData"))

# Export session info
writeLines(capture.output(sessionInfo()), 
           paste0(home, "/Results/Consensus/Post/CC_Post_sessionInfo_training.txt"))