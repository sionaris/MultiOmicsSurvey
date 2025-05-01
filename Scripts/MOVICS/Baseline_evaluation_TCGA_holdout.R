library(MOVICS)
library(openxlsx)
library(dplyr)

# Import environment objects
home = getwd()
load(paste0(home, "/Results/MOVICS_baseline/MOVICS_TCGA_RNAseq-CNV-Methylation-miRNA-SNPs_eval_on_transNEO_env.RData"))
z_transcr = readRDS(paste0(home, "/Resources/TCGA/Surv_standardized_transcr.rds"))
TCGA_holdout = read.xlsx(paste0(home, "/Resources/TCGA/Surv_clinical_data.xlsx"))

# Up-regulated expression features
RNGversion("4.2.2")
timestamp()
TCGA_holdout_ntp_expr_up = runNTP(
  expr = z_transcr,
  templates = dgea.marker.up_1000$templates,
  scaleFlag = TRUE,
  centerFlag = TRUE,
  nPerm = 10000,
  seed = 123,
  distance = "cosine", # default
  doPlot = TRUE,
  height = 8,
  width = 12,
  fig.path = paste0(home, "/Results/MOVICS_baseline"),
  fig.name = "ntp_expr_up_heatmap_TCGA_holdout")
timestamp() # 12.5 min

# down-regulated
RNGversion("4.2.2")
timestamp()
TCGA_holdout_ntp_expr_down = runNTP(
  expr = z_transcr,
  templates = dgea.marker.down_1000$templates,
  scaleFlag = TRUE,
  centerFlag = TRUE,
  nPerm = 10000,
  seed = 123,
  distance = "cosine", # default
  doPlot = TRUE,
  height = 8,
  width = 12,
  fig.path = paste0(home, "/Results/MOVICS_baseline"),
  fig.name = "ntp_expr_down_heatmap_TCGA_holdout")
timestamp() # 12.5 min

# Check concordance
expr_conc = as.data.frame(TCGA_holdout_ntp_expr_down$clust.res) %>%
  dplyr::rename(clust_down = clust) %>%
  inner_join(as.data.frame(TCGA_holdout_ntp_expr_up$clust.res) %>%
               dplyr::rename(clust_up = clust),
             by = "samID")

# This is counter-intuitive but due to opposite directions of deregulation this
# is how it works (perhaps this was expected)
expr_conc$agreement = ifelse(expr_conc$clust_down!=expr_conc$clust_up, "Yes", "No")
paste("Agremeent of NTP subtypes with respect to expression data from the external cohort is: ",
      length(which(expr_conc$agreement == "Yes"))/nrow(expr_conc)*100, "% (", nrow(expr_conc),
      " samples).")

# Import styles and helpful objects
scheme = readRDS("Resources/scheme.rds")
cols_of_interest = colnames(scheme$var2comp)[2:ncol(scheme$var2comp)]
TCGA_holdout_var2comp = TCGA_holdout[, c("Sample.ID", cols_of_interest)] %>%
  inner_join(expr_conc %>% dplyr::select(Sample.ID = samID, `Consensus Subtype` = clust_up),
             by = "Sample.ID") %>%
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
TCGA_holdout_var2comp[, cols_to_numeric] <- lapply(TCGA_holdout_var2comp[, cols_to_numeric, drop = FALSE],
                                      function(x) as.numeric(x))

# Remove unknown levels for statistical tests
TCGA_holdout_var2comp_nonas = TCGA_holdout_var2comp
for (i in 1:ncol(TCGA_holdout_var2comp)) {
  nas = which(TCGA_holdout_var2comp[, i] == "Unknown")
  TCGA_holdout_var2comp_nonas[nas, i] = NA
  empties = which(TCGA_holdout_var2comp[, i] == "")
  TCGA_holdout_var2comp_nonas[empties, i] = NA
}
rm(nas, empties); gc()

# Statistical comparisons
clin_comp = compClinvar_single_algorithm(algorithm_name = "CS",
                                         moic.res = TCGA_holdout_ntp_expr_up,
                                         var2comp = TCGA_holdout_var2comp_nonas %>%
                                           mutate(`Consensus Subtype` = paste0("CS", `Consensus Subtype`)),,
                                         strata = "Consensus Subtype",
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
                                         tab.name = "TCGA_holdout_Summary_of_clinical_variables",
                                         res.path = paste0(home, "/Results/MOVICS_baseline"),
                                         output_pdf = TRUE,
                                         pdf_level_col_width = c("7em", "10em"),
                                         pdf_count_col_width = "10em",
                                         pdf_pval_col_width = "3em",
                                         pdf_test_col_width = "8em",
                                         pdf_tab_font_size = 9)

clin_ordinal_comp = compClinvar_ordinal_single_algorithm(algorithm_name = "CS",
                                                         moic.res = TCGA_holdout_ntp_expr_up,
                                                         var2comp = TCGA_holdout_var2comp_nonas %>%
                                                           select(number_of_lymphnodes_positive_by_ihc, 
                                                                  number_of_lymphnodes_positive_by_he, `Consensus Subtype`),
                                                         strata = "Consensus Subtype",
                                                         ordinalVars = c("number_of_lymphnodes_positive_by_ihc",
                                                                         "number_of_lymphnodes_positive_by_he"),
                                                         includeNA = FALSE,
                                                         tab.name = "Summary of ordinal clinical variables - TCGA holdout",
                                                         res.path = paste0(home, "/Results/MOVICS_baseline"),
                                                         output_pdf = TRUE,
                                                         pdf_template_loc = paste0(home, "/Scripts/automated_scripts/clincomp_template.Rmd"),
                                                         pdf_level_col_width = c("7em", "10em"),
                                                         pdf_count_col_width = "10em",
                                                         pdf_pval_col_width = "3em",
                                                         pdf_test_col_width = "8em",
                                                         pdf_tab_font_size = 9)

# Run PAM
RNGversion("4.2.2")
set.seed(123)
TCGA_holdout_pam = runPAM(train.expr = plotdata$RNAseq,
                      moic.res   = consensus,
                      test.expr  = as.matrix(z_transcr))

# NTP TCGA_holdout vs PAM TCGA_holdout
runKappa(subt1 = as.numeric(TCGA_holdout_ntp_expr_up$clust.res$clust),
         subt2 = as.numeric(TCGA_holdout_pam$clust.res$clust),
         subt1.lab = "TCGA_holdout NTP",
         subt2.lab = "TCGA_holdout PAM",
         height = 8,
         width = 8,
         fig.path = paste0(home, "/Results/MOVICS_baseline"),
         fig.name = "kappa_NTP_vs_PAM_TCGA_holdout")

save.image(paste0(home, "/Results/MOVICS_baseline/TCGA_holdout_evaluation.RData"))