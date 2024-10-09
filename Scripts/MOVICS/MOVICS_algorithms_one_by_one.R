# Libraries and imports #####
library(MOVICS)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(ComplexHeatmap)
library(igraph)
library(ggraph)
library(tidygraph)
library(openxlsx)
library(ggpubr)
library(plotly)
library(data.table)
library(colorspace)

load("Results/MOVICS_baseline/MO_comparisons/MO_comparisons_MOVICS_TCGA_RNAseq-CNV-Methylation-miRNA-SNPs_eval_on_transNEO_env.RData")
source("Scripts/automated_scripts/custom_functions.R")
source("Scripts/automated_scripts/modified_MOVICS_functions.R")

# Ensure reproducibility
RNGversion("4.2.2")
set.seed(123)

# Chi-square tests between clusterings and clinical variables #####
# Bias-corrected Cramer's V calculation using package rcompanion:
unbiased.cv.test = function(x, string, digits = 3) {
  CV = rcompanion::cramerV(x, bias.correct = TRUE)
  return(list(text = paste0("Bias-corrected Cramer's V / Phi for ", 
                            string, ": ", round(as.numeric(CV), digits)),
              value = round(as.numeric(CV), digits)))
}

chisq_outputs = list()
algorithms = colnames(clust_annot_pheno)[c(2:11)]

# Modify the clust_annot_pheno object for plotting
clust_annot_pheno = annCol %>% mutate(samID = rownames(.)) %>%
  inner_join(clust_annot, by = "samID")
rownames(clust_annot_pheno) = clust_annot_pheno$samID
for(algorithm in algorithms) {
  clust_annot_pheno[[algorithm]] = paste0(algorithm, clust_annot_pheno[[algorithm]])
}

clust_annot_pheno_nonas = clust_annot_pheno
for(i in 1:ncol(clust_annot_pheno_nonas)) {
  clust_annot_pheno_nonas[, i] = as.character(clust_annot_pheno_nonas[, i])
  nas = which(clust_annot_pheno_nonas[, i] == "Unknown")
  clust_annot_pheno_nonas[nas, i] = NA
  clust_annot_pheno_nonas[, i] = factor(clust_annot_pheno_nonas[, i])
}
rm(nas); gc()

# Variables of interest
voi = setdiff(colnames(clust_annot_pheno_nonas), c(algorithms, "samID"))

for (a in 1:length(algorithms)) {
  output = as.data.frame(matrix(NA, nrow = 0, ncol = 4))
  for (v in 1:length(voi)){
    keepers = which(!is.na(clust_annot_pheno_nonas[, voi[v]]))
    test = suppressWarnings(chisq.test(table(clust_annot_pheno_nonas[keepers, algorithms[a]], 
                                             clust_annot_pheno_nonas[keepers, voi[v]])))
    chifit_p = test$p.value
    chifit_xsq = test$statistic
    chifit_cv = suppressWarnings(unbiased.cv.test(table(clust_annot_pheno_nonas[keepers, algorithms[a]], 
                                                        clust_annot_pheno_nonas[keepers, voi[v]]),
                                                  string = voi[i],
                                                  digits = 3)$value)
    comparison = paste0(voi[v], " vs ", algorithms[a], " cluster")
    output = rbind(output, c(comparison, chifit_p, chifit_xsq, chifit_cv))
    rm(test, comparison, chifit_p, chifit_xsq, chifit_cv, keepers)
  }
  colnames(output) = c("Comparison", "p-value", "Statistic", "Cramer's V")
  chisq_outputs[[a]] = output
  rm(output)
}

rm(a, v); gc()

chisq_wb = createWorkbook()
for (i in 1:length(chisq_outputs)) {
  addWorksheet(chisq_wb, algorithms[i])
  writeData(chisq_wb, algorithms[i], chisq_outputs[[i]])
}
saveWorkbook(chisq_wb, file = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/chisq_tables.xlsx"),
             overwrite = TRUE); rm(chisq_wb)

names(chisq_outputs) = algorithms

# Create output directories in the MO_comparisons subdir
for (algorithm in algorithms) {
  if (!dir.exists(paste0(home, "/Results/MOVICS_baseline/MO_comparisons/", 
                         algorithm, "_extra"))) {
    dir.create(paste0(home, "/Results/MOVICS_baseline/MO_comparisons/", 
                      algorithm, "_extra"))
  }
}

# SNF #####
aff_CNV = normalize_affinity_matrix(
  SNFtool::affinityMatrix(
    SNFtool::dist2(t(input$CNV),
                   t(input$CNV)
    ),
    K = 30, sigma = 0.5)
)
colnames(aff_CNV) = rownames(aff_CNV) = colnames(input$CNV)

aff_rna = normalize_affinity_matrix(
  SNFtool::affinityMatrix(
    SNFtool::dist2(t(input$RNAseq),
                   t(input$RNAseq)
    ),
    K = 30, sigma = 0.5)
)
colnames(aff_rna) = rownames(aff_rna) = colnames(input$RNAseq)

aff_mut = normalize_affinity_matrix(
  SNFtool::affinityMatrix(
    as.matrix(dist(as.matrix(t(input$SNPs)),
                   as.matrix(t(input$SNPs)),
                   method = "binary")),
    K = 30, sigma = 0.5)
)
colnames(aff_mut) = rownames(aff_mut) = colnames(input$SNPs)

aff_Methyl = normalize_affinity_matrix(
  SNFtool::affinityMatrix(
    SNFtool::dist2(t(input$Methylation),
                   t(input$Methylation)
    ),
    K = 30, sigma = 0.5)
)
colnames(aff_Methyl) = rownames(aff_Methyl) = colnames(input$Methylation)

aff_miRNA = normalize_affinity_matrix(
  SNFtool::affinityMatrix(
    SNFtool::dist2(t(input$miRNA),
                   t(input$miRNA)
    ),
    K = 30, sigma = 0.5)
)
colnames(aff_miRNA) = rownames(aff_miRNA) = colnames(input$miRNA)

aff_final = moic.res.list[["SNF"]][["fit"]]
colnames(aff_final) = rownames(aff_final) = moic.res.list[["SNF"]][["clust.res"]]$samID

# Heatmaps ###
colors_heatmap = rev(colorRampPalette(viridisLite::magma(10))(255))
cluster_colors_heatmap = c("#2EC4B6", "#E71D36")
afh_colnames = voi

# CNV
create_MO_heatmap(matrix = aff_CNV, algorithm = "SNF", 
                  need.diag.zero = TRUE, 
                  clust_annot_pheno = clust_annot_pheno %>%
                    mutate(SNF = paste0("MOVICS_", SNF)) %>%
                    select(samID, all_of(afh_colnames), SNF),
                  afh_colnames = afh_colnames, 
                  colors = colors_heatmap,
                  annColors = annColors,
                  heatmap_title = "CNV first affinity heatmap",
                  cluster_colors = cluster_colors_heatmap,
                  legend_title = "Normalized affinity",
                  cluster_cols_flag = FALSE,
                  cluster_rows_flag = FALSE,
                  splits_flag = TRUE,
                  output_file_name = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/SNF_extra/aff_CNV_heatmap.png"))

# RNAseq
create_MO_heatmap(matrix = aff_rna, algorithm = "SNF", 
                  need.diag.zero = TRUE, 
                  clust_annot_pheno = clust_annot_pheno %>%
                    mutate(SNF = paste0("MOVICS_", SNF)) %>%
                    select(samID, all_of(afh_colnames), SNF),
                  afh_colnames = afh_colnames, 
                  colors = colors_heatmap,
                  annColors = annColors,
                  heatmap_title = "RNAseq first affinity heatmap",
                  cluster_colors = cluster_colors_heatmap,
                  legend_title = "Normalized affinity",
                  cluster_cols_flag = FALSE,
                  cluster_rows_flag = FALSE,
                  splits_flag = TRUE,
                  output_file_name = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/SNF_extra/aff_RNAseq_heatmap.png"))

# miRNA
create_MO_heatmap(matrix = aff_miRNA, algorithm = "SNF", 
                  need.diag.zero = TRUE, 
                  clust_annot_pheno = clust_annot_pheno %>%
                    mutate(SNF = paste0("MOVICS_", SNF)) %>%
                    select(samID, all_of(afh_colnames), SNF),
                  afh_colnames = afh_colnames, 
                  colors = colors_heatmap,
                  annColors = annColors,
                  heatmap_title = "miRNA first affinity heatmap",
                  cluster_colors = cluster_colors_heatmap,
                  legend_title = "Normalized affinity",
                  cluster_cols_flag = FALSE,
                  cluster_rows_flag = FALSE,
                  splits_flag = TRUE,
                  output_file_name = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/SNF_extra/aff_miRNA_heatmap.png"))

# Methylation
create_MO_heatmap(matrix = aff_Methyl, algorithm = "SNF", 
                  need.diag.zero = TRUE, 
                  clust_annot_pheno = clust_annot_pheno %>%
                    mutate(SNF = paste0("MOVICS_", SNF)) %>%
                    select(samID, all_of(afh_colnames), SNF),
                  afh_colnames = afh_colnames, 
                  colors = colors_heatmap,
                  annColors = annColors,
                  heatmap_title = "Methylation first affinity heatmap",
                  cluster_colors = cluster_colors_heatmap,
                  legend_title = "Normalized affinity",
                  cluster_cols_flag = FALSE,
                  cluster_rows_flag = FALSE,
                  splits_flag = TRUE,
                  output_file_name = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/SNF_extra/aff_Methylation_heatmap.png"))

# SNPs
create_MO_heatmap(matrix = aff_mut, algorithm = "SNF", 
                  need.diag.zero = TRUE, 
                  clust_annot_pheno = clust_annot_pheno %>%
                    mutate(SNF = paste0("MOVICS_", SNF)) %>%
                    select(samID, all_of(afh_colnames), SNF),
                  afh_colnames = afh_colnames, 
                  colors = colors_heatmap,
                  annColors = annColors,
                  heatmap_title = "SNPs first affinity heatmap",
                  cluster_colors = cluster_colors_heatmap,
                  legend_title = "Normalized affinity",
                  cluster_cols_flag = FALSE,
                  cluster_rows_flag = FALSE,
                  splits_flag = TRUE,
                  output_file_name = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/SNF_extra/aff_SNPs_heatmap.png"))

# Final affinity matrix
create_MO_heatmap(matrix = aff_final, algorithm = "SNF", 
                  need.diag.zero = TRUE, 
                  clust_annot_pheno = clust_annot_pheno %>%
                    mutate(SNF = paste0("MOVICS_", SNF)) %>%
                    select(samID, all_of(afh_colnames), SNF),
                  afh_colnames = afh_colnames, 
                  colors = colors_heatmap,
                  annColors = annColors,
                  heatmap_title = "Final affinity heatmap",
                  cluster_colors = cluster_colors_heatmap,
                  legend_title = "Normalized affinity",
                  cluster_cols_flag = FALSE,
                  cluster_rows_flag = FALSE,
                  splits_flag = TRUE,
                  output_file_name = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/SNF_extra/aff_final_affinity_heatmap.png"))

# Final affinity matrix with clustered rows and columns
create_MO_heatmap(matrix = aff_final, algorithm = "SNF", 
                  need.diag.zero = TRUE, 
                  clust_annot_pheno = clust_annot_pheno %>%
                    mutate(SNF = paste0("MOVICS_", SNF)) %>%
                    select(samID, all_of(afh_colnames), SNF),
                  afh_colnames = afh_colnames, 
                  colors = colors_heatmap,
                  annColors = annColors,
                  heatmap_title = "Final affinity heatmap",
                  cluster_colors = cluster_colors_heatmap,
                  legend_title = "Normalized affinity",
                  cluster_cols_flag = TRUE,
                  cluster_rows_flag = TRUE,
                  splits_flag = FALSE,
                  output_file_name = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/SNF_extra/hclust_aff_final_affinity_heatmap.png"))


# PCA plots ###
snf_clust_res = clust_annot_pheno %>% dplyr::select(samID, SNF) %>%
  mutate(SNF = paste0("MOVICS_", SNF))

# CNV
pca_from_sim_matrix(sim_matrix = aff_CNV, algorithm = "SNF", clust_res = snf_clust_res,
                    cluster_colors = c("#2EC4B6", "#E71D36"), 
                    output_path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/SNF_extra"), 
                    title_add = "CNV")

# RNA
pca_from_sim_matrix(sim_matrix = aff_rna, algorithm = "SNF", clust_res = snf_clust_res,
                    cluster_colors = c("#2EC4B6", "#E71D36"), 
                    output_path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/SNF_extra"), 
                    title_add = "RNA-seq")

# SNPs
pca_from_sim_matrix(sim_matrix = aff_mut, algorithm = "SNF", clust_res = snf_clust_res,
                    cluster_colors = c("#2EC4B6", "#E71D36"),
                    output_path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/SNF_extra"), 
                    title_add = "SNPs")

# Methylation
pca_from_sim_matrix(sim_matrix = aff_Methyl, algorithm = "SNF", clust_res = snf_clust_res,
                    cluster_colors = c("#2EC4B6", "#E71D36"), 
                    output_path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/SNF_extra"), 
                    title_add = "Methylation")

# miRNA
pca_from_sim_matrix(sim_matrix = aff_miRNA, algorithm = "SNF", clust_res = snf_clust_res,
                    cluster_colors = c("#2EC4B6", "#E71D36"), 
                    output_path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/SNF_extra"), 
                    title_add = "miRNA")

# Final matrix
pca_from_sim_matrix(sim_matrix = aff_final, algorithm = "SNF", clust_res = snf_clust_res,
                    cluster_colors = c("#2EC4B6", "#E71D36"), 
                    output_path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/SNF_extra"), 
                    title_add = "Final Fusion")

# Setup for barcharts ###
# Stage
scale_fill_stage = scale_fill_manual(values = c(`Stage I` = "#00C9FF", 
                                                `Stage II` = "#099CF5", 
                                                `Stage III` = "#097BF5", 
                                                `Stage IV` = "#0B5684", 
                                                `Unknown` = "grey40"))

# Lymph node status
scale_fill_lymph_node_status = scale_fill_manual(values = c(No = "grey75", 
                                                            Yes = "#4A0558", 
                                                            Unknown = "grey40"))

# ER status
scale_fill_ER_status = scale_fill_manual(values = c(Negative = "#C11D9C", 
                                                    Positive = "#0F1682", 
                                                    Unknown = "grey40"))

# PR status
scale_fill_PR_status = scale_fill_manual(values = c(Indeterminate = "aliceblue", 
                                                    Positive = "dodgerblue4", 
                                                    Negative = "#F0C6C3", 
                                                    Unknown = "grey40"))

# HER2 status
scale_fill_HER2_status = scale_fill_manual(values = c(Negative = "#0B9EF8", 
                                                      Positive = "#560DA7", 
                                                      Indeterminate = "mistyrose1", 
                                                      Equivocal = "hotpink4", 
                                                      Unknown = "grey40"))

# Vital status
scale_fill_vital_status = scale_fill_manual(values = c(Alive = "lightpink1", 
                                                       Dead = "black", 
                                                       Unknown = "grey40"))

# Ethnicity
scale_fill_ethnicity = scale_fill_manual(values = c(`Hispanic or latino` = "#E58606", 
                                                    `Not hispanic or latino` = "#24796C", 
                                                    Unknown = "grey40"))

# Race
scale_fill_race = scale_fill_manual(values = c(`American indian or alaska native` = "#E73F74", 
                                               Asian = "#3969AC", 
                                               `Black or african american` = "#666666", 
                                               White = "beige", 
                                               Unknown = "grey40"))

# Metastasis
scale_fill_metastasis = scale_fill_manual(values = c(Yes = "deeppink4", 
                                                     No = "cadetblue2", 
                                                     Unknown = "grey40"))

# Histology
scale_fill_histology = scale_fill_manual(values = c(`Infiltrating Carcinoma NOS` = "#88CCEE", 
                                                    `Infiltrating Ductal Carcinoma` = "#CC6677", 
                                                    `Infiltrating Lobular Carcinoma` = "#DDCC77", 
                                                    `Medullary Carcinoma` = "#117733", 
                                                    `Metaplastic Carcinoma` = "#332288", 
                                                    Mixed = "#AA4499", 
                                                    `Mucinous Carcinoma` = "#44AA99", 
                                                    Other = "#999933", 
                                                    Unknown = "grey40"))

# Menopausal status
scale_fill_menopausal_status = scale_fill_manual(values = c(Indeterminate = "mistyrose2", 
                                                            `Pre-menopausal` = "#FAA476", 
                                                            Perimenopausal = "#DC3977", 
                                                            `Post-menopausal` = "#7C1D6F", 
                                                            Unknown = "grey40"))

# Combine all scales into a list
barchart_scales = list(scale_fill_stage, scale_fill_lymph_node_status, scale_fill_ER_status, 
                       scale_fill_PR_status, scale_fill_HER2_status, scale_fill_vital_status, 
                       scale_fill_ethnicity, scale_fill_race, scale_fill_metastasis, 
                       scale_fill_histology, scale_fill_menopausal_status)

# Name the scales accordingly
names(barchart_scales) = c("Stage", "Lymph node status", "ER status", "PR status", "HER2 status", 
                           "Vital status", "Ethnicity", "Race", "Metastasis", "Histology", 
                           "Menopausal status")

# Bar charts with clinical variables of interest ###
SNF_barcharts = list()
plotdata_bar = clust_annot_pheno_nonas %>%
  dplyr::mutate(SNF = paste0("MOVICS_", SNF))
plotdata_bar$SNF = factor(plotdata_bar$SNF)

for (i in 1:length(voi)) {
  chifit = chisq_outputs[["SNF"]]
  loc = which(grepl(voi[i], chifit$Comparison))
  chifit = chifit[loc, ]
  SNF_barcharts[[i]] = create_annot_barchart(plotdata = plotdata_bar, fill = voi[i],
                                             chifit = chifit,
                                             na.action = "na.omit",
                                             algorithm = "SNF",
                                             barchart_ylim = 650,
                                             text_y = 630, rect_ymin = 530,
                                             rect_ymax = 650, x_annot = 1.5,
                                             v_gap = 35, rect_xmin = 1,
                                             rect_xmax = 2, 
                                             annot_text_size = 2.25,
                                             legend.text.size = 5,
                                             x.axis.text.size = 5) +
    barchart_scales[[voi[i]]]
  print(SNF_barcharts[[i]])
  ggsave(filename = paste0("SNF_", voi[i], "_barchart.png"),
         path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/SNF_extra"), 
         width = 2320, height = 2320, device = 'png', units = "px",
         dpi = 700)
  dev.off()
}
names(SNF_barcharts) = voi
rm(loc, chifit)

# Multiplot (PNG) - bar charts
ggarrange(SNF_barcharts[[1]], SNF_barcharts[[2]], SNF_barcharts[[3]],
          SNF_barcharts[[4]], SNF_barcharts[[5]], SNF_barcharts[[6]],
          SNF_barcharts[[7]], SNF_barcharts[[8]], SNF_barcharts[[9]],
          SNF_barcharts[[10]], SNF_barcharts[[11]],
          ncol = 3, nrow = 4, labels = c("A", "B", "C", "D", "E", "F", "G", "H",
                                         "I", "J", "K"),
          font.label = list(size = 8, face = "bold", color ="black"))
ggsave(filename = "Multiplot_SNF_barcharts.png",
       path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/SNF_extra"), 
       width = 7000, height = 8000, device = 'png', units = "px",
       dpi = 700)
dev.off()

# Just significant ones now
SNF_barcharts_sig = list()
plotdata_bar_sig = clust_annot_pheno_nonas %>% dplyr::select(SNF, Race, Histology, 
                                                       `ER status`, `PR status`, Stage) %>%
  dplyr::mutate(SNF = paste0("MOVICS_", SNF))
plotdata_bar_sig$SNF = factor(plotdata_bar_sig$SNF)
voi_sig = setdiff(colnames(plotdata_bar_sig), "SNF")
for (i in 1:length(voi_sig)) {
  chifit = chisq_outputs[["SNF"]]
  loc = which(grepl(voi_sig[i], chifit$Comparison))
  chifit = chifit[loc, ]
  SNF_barcharts_sig[[i]] = create_annot_barchart(plotdata = plotdata_bar_sig, fill = voi_sig[i],
                                             chifit = chifit,
                                             na.action = "na.omit",
                                             algorithm = "SNF",
                                             barchart_ylim = 650,
                                             text_y = 630, rect_ymin = 530,
                                             rect_ymax = 650, x_annot = 1.5,
                                             v_gap = 35, rect_xmin = 1,
                                             rect_xmax = 2, 
                                             annot_text_size = 2.25,
                                             legend.text.size = 5,
                                             x.axis.text.size = 5) +
    barchart_scales[[voi_sig[i]]]
  print(SNF_barcharts_sig[[i]])
  ggsave(filename = paste0("sig_SNF_", voi_sig[i], "_barchart.png"),
         path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/SNF_extra"), 
         width = 2320, height = 2320, device = 'png', units = "px",
         dpi = 700)
  dev.off()
}
names(SNF_barcharts_sig) = voi_sig
rm(loc, chifit)

# Multiplot (PNG) - bar charts
ggarrange(SNF_barcharts_sig[[1]], SNF_barcharts_sig[[2]], SNF_barcharts_sig[[3]],
          SNF_barcharts_sig[[4]], SNF_barcharts_sig[[5]], 
          ncol = 2, nrow = 3, labels = c("A", "B", "C", "D", "E"),
          font.label = list(size = 8, face = "bold", color ="black"))
ggsave(filename = "sig_Multiplot_SNF_barcharts.png",
       path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/SNF_extra"), 
       width = 5500, height = 7000, device = 'png', units = "px",
       dpi = 700)
dev.off()

# Sunburst plot ###
library(plotly)
Pheno_sunburst_SNF = clust_annot_pheno
Pheno_sunburst_SNF$`ER status` = gsub("Unknown", "Unkn ER status", Pheno_sunburst_SNF$`ER status`)
Pheno_sunburst_SNF$Stage = gsub("Unknown", "Unkn stage", 
                                Pheno_sunburst_SNF$Stage)
Pheno_sunburst_SNF = Pheno_sunburst_SNF %>%
  dplyr::select(SNF, `ER status`, Stage) %>%
  group_by(SNF, `ER status`, Stage) %>%
  summarise(Counts = n()) %>%
  as.data.frame()

sunburst_coloring_SNF = data.frame(stringsAsFactors = FALSE,
                                   colors = tolower(gplots::col2hex(c("#2EC4B6", "#E71D36", 
                                                                      "#C11D9C", "#0F1682",  "grey40",
                                                                      "#00C9FF", "#099CF5", "#097BF5", 
                                                                      "#0B5684", "grey40"))),
                                   labels = c("SNF1", "SNF2",
                                              "Negative", "Positive", "Unkn ER status",
                                              "Stage I", "Stage II", "Stage III",
                                              "Stage IV", "Unkn stage"))

sunburstDF_SNF = as.sunburstDF(Pheno_sunburst_SNF, value_column = "Counts", add_root = FALSE) %>%
  inner_join(sunburst_coloring_SNF, by = "labels")

pie_SNF = plot_ly() %>%
  add_trace(ids = sunburstDF_SNF$ids, labels= sunburstDF_SNF$labels, 
            parents = sunburstDF_SNF$parents, 
            values= sunburstDF_SNF$values, type='sunburst', branchvalues = 'total',
            insidetextorientation='radial', maxdepth = 5,
            marker = list(colors = sunburstDF_SNF$colors)) %>%
  layout(
    grid = list(columns =1, rows = 1),
    margin = list(l = 0, r = 0, b = 0, t = 0)
  )
pie_SNF
rm(Pheno_sunburst_SNF, sunburstDF_SNF, sunburst_coloring_SNF, pie_SNF); gc()

# Draw graphs from affinity matrices
aff_CNV_S = calculate_S(aff_CNV)
aff_rna_S = calculate_S(aff_rna)
aff_miRNA_S = calculate_S(aff_miRNA)
aff_Methyl_S = calculate_S(aff_Methyl)
aff_mut_S = calculate_S(aff_mut)
aff_final_S = calculate_S(aff_final)

list_aff_S = list(aff_CNV_S, aff_rna_S, aff_miRNA_S, aff_Methyl_S, 
                  aff_mut_S, aff_final_S)
names(list_aff_S) = c(paste0("CNV Original Affinity 30-NN Graph"),
                      paste0("RNAseq Original Affinity 30-NN Graph"),
                      paste0("miRNA Original Affinity 30-NN Graph"),
                      paste0("Methylation Original Affinity 30-NN Graph"),
                      paste0("SNPs Original Affinity 30-NN Graph"),
                      paste0("Final Fused Affinity 30-NN Graph"))

for (i in 1:length(list_aff_S)) {
  
  # Prepare the graph object
  g <- graph_from_adjacency_matrix(list_aff_S[[i]], 
                                   mode = "directed", weighted = TRUE, diag = FALSE)
  g <- delete_edges(g, E(g)[weight == 0])
  E(g)$width <- sqrt(E(g)$weight) * 5  # Example transformation for visibility
  nodes_data <- data.frame(name = V(g)$name) %>%
    inner_join(clust_annot_pheno %>% dplyr::select(samID, SNF), by = c("name" = "samID"))
  
  # Set SNF as a factor for coloring
  nodes_data$SNF <- as.factor(paste0("MOVICS_", nodes_data$SNF))
  V(g)$SNF <- nodes_data$SNF
  
  # Set color based on SNF
  V(g)$color <- ifelse(V(g)$SNF == "MOVICS_SNF1", "#2EC4B6", "#E71D36")
  
  png(paste0(home, "/Results/MOVICS_baseline/MO_comparisons/SNF_extra/",
             names(list_aff_S)[i], ".png"),
      width = 6000, height = 6000, res = 700)
  
  par(mar = c(2, 2, 2, 5))  # Adjust right margin to accommodate legend
  
  # Plot the graph with a layout that spreads nodes well
  plot(g, vertex.color = V(g)$color,
       edge.width = E(g)$width,
       vertex.size = 4, 
       vertex.label = NA, 
       edge.color = "gray85",
       layout = layout_with_fr(g),  # Use Fruchterman-Reingold layout
       main = "",
       edge.arrow.size=0.05)
  
  # Add title with reduced size using title() function
  title(main = names(list_aff_S)[i], cex.main = 1.7)
  
  # Add a legend to the right of the plot
  legend("bottomright", 
         title="Node Color Legend",    
         legend=c("MOVICS_SNF1", "MOVICS_SNF2"), 
         fill=c("#2EC4B6", "#E71D36"),  
         cex=0.7,      
         box.lwd=1)  
  
  dev.off() 
}
rm(g, nodes_data)

# See concordance with the final consensus
clust_annot_pheno2 = clust_annot_pheno %>%
  inner_join(as.data.frame(consensus$clust.res))

# Create a list of cross-tabulations between the CS and individual algorithm labels
CS_comp_list = vector("list", 10)
names(CS_comp_list) = algorithms
for(i in 1:length(CS_comp_list)) {
  CS_comp_list[[i]] <- list(
    table = NULL,
    NMI = NULL,
    ARI = NULL
  )
}

# Compare MOVICS SNF to MOVICS consensus
CS_comp_list[["SNF"]]$table = table(clust_annot_pheno2$SNF, 
                              paste0("CS", clust_annot_pheno2$clust))

CS_comp_list[["SNF"]]$ARI = calculate_ari_index(cluster_df1 = snf_clust_res %>%
                                                  dplyr::rename(Cluster = SNF) %>%
                                                  mutate(Cluster = gsub("MOVICS_SNF", "", Cluster)),
                                                cluster_df2 = as.data.frame(consensus$clust.res) %>%
                                                  dplyr::rename(Cluster = clust),
                                                sample_col = "samID",
                                                clust_col = "Cluster",
                                                suffixes = c("_MOVICS_SNF", "_CS"))

CS_comp_list[["SNF"]]$NMI = calculate_nmi_index(cluster_df1 = snf_clust_res %>%
                                                  dplyr::rename(Cluster = SNF) %>%
                                                  mutate(Cluster = gsub("MOVICS_SNF", "", Cluster)),
                                                cluster_df2 = as.data.frame(consensus$clust.res) %>%
                                                  dplyr::rename(Cluster = clust),
                                                sample_col = "samID",
                                                clust_col = "Cluster",
                                                suffixes = c("_MOVICS_SNF", "_CS"))

# Print all comparison data
print(CS_comp_list$SNF)

# CIMLR #####

# PCA from original matrices ###
CIMLR_clust_res = clust_annot_pheno %>% dplyr::select(samID, CIMLR) %>%
  mutate(CIMLR = paste0("MOVICS_", CIMLR))

# RNA
pca_from_original_matrix(mydata = plotdata$RNAseq, 
                         algorithm = "CIMLR", 
                         clust_res = CIMLR_clust_res,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/CIMLR_extra"),
                         title_add = "RNA-seq")

# miRNA
pca_from_original_matrix(mydata = plotdata$miRNA, 
                         algorithm = "CIMLR", 
                         clust_res = CIMLR_clust_res,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/CIMLR_extra"),
                         title_add = "miRNA")

# CNV
pca_from_original_matrix(mydata = plotdata$CNV, 
                         algorithm = "CIMLR", 
                         clust_res = CIMLR_clust_res,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/CIMLR_extra"),
                         title_add = "CNV")

# Use multidimensional scaling for SNPs
# Features must be in rows
mds_from_original_matrix(matrix = plotdata$SNPs, dist_method = "binary",
                         algorithm = "CIMLR", 
                         clust_res = CIMLR_clust_res,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/CIMLR_extra"),
                         title_add = "SNPs")

# Methylation
pca_from_original_matrix(mydata = plotdata$Methylation, 
                         algorithm = "CIMLR", 
                         clust_res = CIMLR_clust_res,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/CIMLR_extra"),
                         title_add = "Methylation")

# Draw a heatmap of the final S matrix ###
cimlr_matrix = moic.res.list$CIMLR$fit$S
dimnames(cimlr_matrix) = dimnames(aff_final)
create_MO_heatmap(matrix = cimlr_matrix, algorithm = "CIMLR", 
                  need.diag.zero = FALSE, # already zero
                  clust_annot_pheno = clust_annot_pheno %>%
                    mutate(CIMLR = paste0("MOVICS_", CIMLR)) %>%
                    select(samID, all_of(afh_colnames), CIMLR),
                  afh_colnames = afh_colnames, 
                  colors = colors_heatmap,
                  annColors = annColors,
                  heatmap_title = "Final CIMLR similarity heatmap",
                  cluster_colors = cluster_colors_heatmap,
                  cluster_rows_flag = FALSE,
                  cluster_cols_flag = FALSE,
                  splits_flag = TRUE,
                  legend_title = "Final kernel similarity",
                  output_file_name = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/CIMLR_extra/CIMLR_final_S_matrix_heatmap.png"))

# Bar charts with clinical variables of interest ###
CIMLR_barcharts = list()
plotdata_bar = clust_annot_pheno_nonas %>%
  dplyr::mutate(CIMLR = paste0("MOVICS_", CIMLR))
plotdata_bar$CIMLR = factor(plotdata_bar$CIMLR)
for (i in 1:length(voi)) {
  chifit = chisq_outputs[["CIMLR"]]
  loc = which(grepl(voi[i], chifit$Comparison))
  chifit = chifit[loc, ]
  CIMLR_barcharts[[i]] = create_annot_barchart(plotdata = plotdata_bar, fill = voi[i],
                                               chifit = chifit,
                                               na.action = "na.omit",
                                               algorithm = "CIMLR",
                                               barchart_ylim = 650,
                                               text_y = 630, rect_ymin = 530,
                                               rect_ymax = 650, x_annot = 1.5,
                                               v_gap = 35, rect_xmin = 1,
                                               rect_xmax = 2, 
                                               annot_text_size = 2.25,
                                               legend.text.size = 5,
                                               x.axis.text.size = 5) +
    barchart_scales[[voi[i]]]
  print(CIMLR_barcharts[[i]])
  ggsave(filename = paste0("CIMLR_", voi[i], "_barchart.png"),
         path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/CIMLR_extra"), 
         width = 2320, height = 2320, device = 'png', units = "px",
         dpi = 700)
  dev.off()
}
names(CIMLR_barcharts) = voi
rm(loc, chifit)

# Multiplot (PNG) - bar charts
ggarrange(CIMLR_barcharts[[1]], CIMLR_barcharts[[2]], CIMLR_barcharts[[3]],
          CIMLR_barcharts[[4]], CIMLR_barcharts[[5]], CIMLR_barcharts[[6]],
          CIMLR_barcharts[[7]], CIMLR_barcharts[[8]], CIMLR_barcharts[[9]],
          CIMLR_barcharts[[10]], CIMLR_barcharts[[11]],
          ncol = 3, nrow = 4, labels = c("A", "B", "C", "D", "E", "F", "G", "H",
                                         "I", "J", "K"),
          font.label = list(size = 8, face = "bold", color ="black"))
ggsave(filename = "Multiplot_CIMLR_barcharts.png",
       path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/CIMLR_extra"), 
       width = 7000, height = 8000, device = 'png', units = "px",
       dpi = 700)
dev.off()

# Just significant ones now
CIMLR_barcharts_sig = list()
plotdata_bar_sig = clust_annot_pheno_nonas %>% dplyr::select(CIMLR, Race, `Menopausal status`, 
                                                       `ER status`, `PR status`, `HER2 status`) %>%
  dplyr::mutate(CIMLR = paste0("MOVICS_", CIMLR))
plotdata_bar_sig$CIMLR = factor(plotdata_bar_sig$CIMLR)
voi_sig = setdiff(colnames(plotdata_bar_sig), "CIMLR")
for (i in 1:length(voi_sig)) {
  chifit = chisq_outputs[["CIMLR"]]
  loc = which(grepl(voi_sig[i], chifit$Comparison))
  chifit = chifit[loc, ]
  CIMLR_barcharts_sig[[i]] = create_annot_barchart(plotdata = plotdata_bar_sig, fill = voi_sig[i],
                                                   chifit = chifit,
                                                   na.action = "na.omit",
                                                   algorithm = "CIMLR",
                                                   barchart_ylim = 650,
                                                   text_y = 630, rect_ymin = 530,
                                                   rect_ymax = 650, x_annot = 1.5,
                                                   v_gap = 35, rect_xmin = 1,
                                                   rect_xmax = 2, 
                                                   annot_text_size = 2.25,
                                                   legend.text.size = 5,
                                                   x.axis.text.size = 5) +
    barchart_scales[[voi_sig[i]]]
  print(CIMLR_barcharts_sig[[i]])
  ggsave(filename = paste0("sig_CIMLR_", voi_sig[i], "_barchart.png"),
         path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/CIMLR_extra"), 
         width = 2320, height = 2320, device = 'png', units = "px",
         dpi = 700)
  dev.off()
}
names(CIMLR_barcharts_sig) = voi_sig
rm(loc, chifit)

# Multiplot (PNG) - bar charts
ggarrange(CIMLR_barcharts_sig[[1]], CIMLR_barcharts_sig[[2]], CIMLR_barcharts_sig[[3]],
          CIMLR_barcharts_sig[[4]], CIMLR_barcharts_sig[[5]], 
          ncol = 2, nrow = 3, labels = c("A", "B", "C", "D", "E"),
          font.label = list(size = 8, face = "bold", color ="black"))
ggsave(filename = "sig_Multiplot_CIMLR_barcharts.png",
       path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/CIMLR_extra"), 
       width = 5500, height = 7000, device = 'png', units = "px",
       dpi = 700)
dev.off()

# Sunburst plot ###
Pheno_sunburst_CIMLR = clust_annot_pheno
Pheno_sunburst_CIMLR$`ER status` = gsub("Unknown", "Unkn ER status", Pheno_sunburst_CIMLR$`ER status`)
Pheno_sunburst_CIMLR$`ER status` = gsub("Positive", "ER+", Pheno_sunburst_CIMLR$`ER status`)
Pheno_sunburst_CIMLR$`ER status` = gsub("Negative", "ER-", Pheno_sunburst_CIMLR$`ER status`)
Pheno_sunburst_CIMLR$`HER2 status` = gsub("Unknown", "Unkn HER2 status", 
                                          Pheno_sunburst_CIMLR$`HER2 status`)
Pheno_sunburst_CIMLR$`HER2 status` = gsub("Positive", "HER2+", Pheno_sunburst_CIMLR$`HER2 status`)
Pheno_sunburst_CIMLR$`HER2 status` = gsub("Negative", "HER2-", Pheno_sunburst_CIMLR$`HER2 status`)
Pheno_sunburst_CIMLR$`Menopausal status` = gsub("Unknown", "Unkn Meno status", 
                                                Pheno_sunburst_CIMLR$`Menopausal status`)
Pheno_sunburst_CIMLR$`Menopausal status` = gsub("Indeterminate", "Indeterminate Meno", 
                                                Pheno_sunburst_CIMLR$`Menopausal status`)
Pheno_sunburst_CIMLR = Pheno_sunburst_CIMLR %>%
  dplyr::select(CIMLR, `ER status`, `HER2 status`, `Menopausal status`) %>%
  group_by(CIMLR, `ER status`, `HER2 status`, `Menopausal status`) %>%
  summarise(Counts = n()) %>%
  as.data.frame()

sunburst_coloring_CIMLR = data.frame(stringsAsFactors = FALSE,
                                     colors = tolower(gplots::col2hex(c("#2EC4B6", "#E71D36", 
                                                                        "#C11D9C", "#0F1682",  "grey40",
                                                                        "#0B9EF8", "#560DA7", "mistyrose1", 
                                                                        "hotpink4", "grey40",
                                                                        "mistyrose2", "#FAA476",
                                                                        "#DC3977", "#7C1D6F", "grey40"))),
                                     labels = c("CIMLR1", "CIMLR2",
                                                "ER-", "ER+", "Unkn ER status",
                                                "HER2-", "HER2+", "Indeterminate",
                                                "Equivocal", "Unkn HER2 status",
                                                "Indeterminate Meno", "Pre-menopausal", "Perimenopausal",
                                                "Post-menopausal", "Unkn Meno status"))

sunburstDF_CIMLR = as.sunburstDF(Pheno_sunburst_CIMLR, value_column = "Counts", add_root = FALSE) %>%
  inner_join(sunburst_coloring_CIMLR, by = "labels")

pie_CIMLR = plot_ly() %>%
  add_trace(ids = sunburstDF_CIMLR$ids, labels= sunburstDF_CIMLR$labels, 
            parents = sunburstDF_CIMLR$parents, 
            values= sunburstDF_CIMLR$values, type='sunburst', branchvalues = 'total',
            insidetextorientation='radial', maxdepth = 5,
            marker = list(colors = sunburstDF_CIMLR$colors)) %>%
  layout(
    grid = list(columns =1, rows = 1),
    margin = list(l = 0, r = 0, b = 0, t = 0)
  )
pie_CIMLR
rm(Pheno_sunburst_CIMLR, sunburstDF_CIMLR, sunburst_coloring_CIMLR, pie_CIMLR); gc()

# See concordance with the final consensus
CS_comp_list[["CIMLR"]]$table = table(clust_annot_pheno2$CIMLR, 
                                      paste0("CS", clust_annot_pheno2$clust))

CS_comp_list[["CIMLR"]]$ARI = calculate_ari_index(cluster_df1 = CIMLR_clust_res %>%
                                                    dplyr::rename(Cluster = CIMLR) %>%
                                                    mutate(Cluster = gsub("MOVICS_CIMLR", "", Cluster)),
                                                  cluster_df2 = as.data.frame(consensus$clust.res) %>%
                                                    dplyr::rename(Cluster = clust),
                                                  sample_col = "samID",
                                                  clust_col = "Cluster",
                                                  suffixes = c("_MOVICS_CIMLR", "_CS"))

CS_comp_list[["CIMLR"]]$NMI = calculate_nmi_index(cluster_df1 = CIMLR_clust_res %>%
                                                    dplyr::rename(Cluster = CIMLR) %>%
                                                    mutate(Cluster = gsub("MOVICS_CIMLR", "", Cluster)),
                                                  cluster_df2 = as.data.frame(consensus$clust.res) %>%
                                                    dplyr::rename(Cluster = clust),
                                                  sample_col = "samID",
                                                  clust_col = "Cluster",
                                                  suffixes = c("_MOVICS_CIMLR", "_CS"))

# Print all comparison data
print(CS_comp_list$CIMLR)

# PINSPlus #####
# PCA from original matrices ###
PINSPlus_clust_res = clust_annot_pheno %>% dplyr::select(samID, PINSPlus) %>%
  mutate(PINSPlus = paste0("MOVICS_", PINSPlus))

# RNA
pca_from_original_matrix(mydata = plotdata$RNAseq, 
                         algorithm = "PINSPlus", 
                         clust_res = PINSPlus_clust_res,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/PINSPlus_extra"),
                         title_add = "RNA-seq")

# miRNA
pca_from_original_matrix(mydata = plotdata$miRNA, 
                         algorithm = "PINSPlus", 
                         clust_res = PINSPlus_clust_res,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/PINSPlus_extra"),
                         title_add = "miRNA")

# CNV
pca_from_original_matrix(mydata = plotdata$CNV, 
                         algorithm = "PINSPlus", 
                         clust_res = PINSPlus_clust_res,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/PINSPlus_extra"),
                         title_add = "CNV")

# Use multidimensional scaling for SNPs
# Features must be in rows
mds_from_original_matrix(matrix = plotdata$SNPs, dist_method = "binary",
                         algorithm = "PINSPlus", 
                         clust_res = PINSPlus_clust_res,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/PINSPlus_extra"),
                         title_add = "SNPs")

# Methylation
pca_from_original_matrix(mydata = plotdata$Methylation, 
                         algorithm = "PINSPlus", 
                         clust_res = PINSPlus_clust_res,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/PINSPlus_extra"),
                         title_add = "Methylation")

# Bar charts with clinical variables of interest ###
PINSPlus_barcharts = list()
plotdata_bar = clust_annot_pheno_nonas %>%
  dplyr::mutate(PINSPlus = paste0("MOVICS_", PINSPlus))
plotdata_bar$PINSPlus = factor(plotdata_bar$PINSPlus)
for (i in 1:length(voi)) {
  chifit = chisq_outputs[["PINSPlus"]]
  loc = which(grepl(voi[i], chifit$Comparison))
  chifit = chifit[loc, ]
  PINSPlus_barcharts[[i]] = create_annot_barchart(plotdata = plotdata_bar, fill = voi[i],
                                                  chifit = chifit,
                                                  na.action = "na.omit",
                                                  algorithm = "PINSPlus",
                                                  barchart_ylim = 650,
                                                  text_y = 630, rect_ymin = 530,
                                                  rect_ymax = 650, x_annot = 1.5,
                                                  v_gap = 35, rect_xmin = 1,
                                                  rect_xmax = 2, 
                                                  annot_text_size = 2.25,
                                                  legend.text.size = 5,
                                                  x.axis.text.size = 5) +
    barchart_scales[[voi[i]]]
  print(PINSPlus_barcharts[[i]])
  ggsave(filename = paste0("PINSPlus_", voi[i], "_barchart.png"),
         path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/PINSPlus_extra"), 
         width = 2320, height = 2320, device = 'png', units = "px",
         dpi = 700)
  dev.off()
}
names(PINSPlus_barcharts) = voi
rm(loc, chifit)

# Multiplot (PNG) - bar charts
ggarrange(PINSPlus_barcharts[[1]], PINSPlus_barcharts[[2]], PINSPlus_barcharts[[3]],
          PINSPlus_barcharts[[4]], PINSPlus_barcharts[[5]], PINSPlus_barcharts[[6]],
          PINSPlus_barcharts[[7]], PINSPlus_barcharts[[8]], PINSPlus_barcharts[[9]],
          PINSPlus_barcharts[[10]], PINSPlus_barcharts[[11]],
          ncol = 3, nrow = 4, labels = c("A", "B", "C", "D", "E", "F", "G", "H",
                                         "I", "J", "K"),
          font.label = list(size = 8, face = "bold", color ="black"))
ggsave(filename = "Multiplot_PINSPlus_barcharts.png",
       path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/PINSPlus_extra"), 
       width = 7000, height = 8000, device = 'png', units = "px",
       dpi = 700)
dev.off()

# Just significant ones now
PINSPlus_barcharts_sig = list()
plotdata_bar_sig = clust_annot_pheno_nonas %>% dplyr::select(PINSPlus, `Lymph node status`,
                                                             `ER status`, `PR status`, `HER2 status`,
                                                       Histology, Stage) %>%
  dplyr::mutate(PINSPlus = paste0("MOVICS_", PINSPlus))
plotdata_bar_sig$PINSPlus = factor(plotdata_bar_sig$PINSPlus)
voi_sig = setdiff(colnames(plotdata_bar_sig), "PINSPlus")
for (i in 1:length(voi_sig)) {
  chifit = chisq_outputs[["PINSPlus"]]
  loc = which(grepl(voi_sig[i], chifit$Comparison))
  chifit = chifit[loc, ]
  PINSPlus_barcharts_sig[[i]] = create_annot_barchart(plotdata = plotdata_bar_sig, fill = voi_sig[i],
                                                      chifit = chifit,
                                                      na.action = "na.omit",
                                                      algorithm = "PINSPlus",
                                                      barchart_ylim = 650,
                                                      text_y = 630, rect_ymin = 530,
                                                      rect_ymax = 650, x_annot = 1.5,
                                                      v_gap = 35, rect_xmin = 1,
                                                      rect_xmax = 2, 
                                                      annot_text_size = 2.25,
                                                      legend.text.size = 5,
                                                      x.axis.text.size = 5) +
    barchart_scales[[voi_sig[i]]]
  print(PINSPlus_barcharts_sig[[i]])
  ggsave(filename = paste0("sig_PINSPlus_", voi_sig[i], "_barchart.png"),
         path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/PINSPlus_extra"), 
         width = 2320, height = 2320, device = 'png', units = "px",
         dpi = 700)
  dev.off()
}
names(PINSPlus_barcharts_sig) = voi_sig
rm(loc, chifit)

# Multiplot (PNG) - bar charts
ggarrange(PINSPlus_barcharts_sig[[1]], PINSPlus_barcharts_sig[[2]], PINSPlus_barcharts_sig[[3]],
          PINSPlus_barcharts_sig[[4]], PINSPlus_barcharts_sig[[5]], PINSPlus_barcharts_sig[[6]],
          ncol = 2, nrow = 3, labels = c("A", "B", "C", "D", "E", "F"),
          font.label = list(size = 8, face = "bold", color ="black"))
ggsave(filename = "sig_Multiplot_PINSPlus_barcharts.png",
       path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/PINSPlus_extra"), 
       width = 5500, height = 8000, device = 'png', units = "px",
       dpi = 700)
dev.off()

# Sunburst plot ###
Pheno_sunburst_PINSPlus = clust_annot_pheno
Pheno_sunburst_PINSPlus$`ER status` = gsub("Unknown", "Unkn ER status", Pheno_sunburst_PINSPlus$`ER status`)
Pheno_sunburst_PINSPlus$`ER status` = gsub("Positive", "ER+", Pheno_sunburst_PINSPlus$`ER status`)
Pheno_sunburst_PINSPlus$`ER status` = gsub("Negative", "ER-", Pheno_sunburst_PINSPlus$`ER status`)
Pheno_sunburst_PINSPlus$`HER2 status` = gsub("Unknown", "Unkn HER2 status", 
                                             Pheno_sunburst_PINSPlus$`HER2 status`)
Pheno_sunburst_PINSPlus$`HER2 status` = gsub("Positive", "HER2+", Pheno_sunburst_PINSPlus$`HER2 status`)
Pheno_sunburst_PINSPlus$`HER2 status` = gsub("Negative", "HER2-", Pheno_sunburst_PINSPlus$`HER2 status`)
Pheno_sunburst_PINSPlus = Pheno_sunburst_PINSPlus %>%
  dplyr::select(PINSPlus, `ER status`, `HER2 status`) %>%
  group_by(PINSPlus, `ER status`, `HER2 status`) %>%
  summarise(Counts = n()) %>%
  as.data.frame()

sunburst_coloring_PINSPlus = data.frame(stringsAsFactors = FALSE,
                                        colors = tolower(gplots::col2hex(c("#2EC4B6", "#E71D36", 
                                                                           "#C11D9C", "#0F1682",  "grey40",
                                                                           "#0B9EF8", "#560DA7", "mistyrose1", 
                                                                           "hotpink4", "grey40"))),
                                        labels = c("PINSPlus1", "PINSPlus2",
                                                   "ER-", "ER+", "Unkn ER status",
                                                   "HER2-", "HER2+", "Indeterminate",
                                                   "Equivocal", "Unkn HER2 status"))

sunburstDF_PINSPlus = as.sunburstDF(Pheno_sunburst_PINSPlus, value_column = "Counts", add_root = FALSE) %>%
  inner_join(sunburst_coloring_PINSPlus, by = "labels")

pie_PINSPlus = plot_ly() %>%
  add_trace(ids = sunburstDF_PINSPlus$ids, labels= sunburstDF_PINSPlus$labels, 
            parents = sunburstDF_PINSPlus$parents, 
            values= sunburstDF_PINSPlus$values, type='sunburst', branchvalues = 'total',
            insidetextorientation='radial', maxdepth = 5,
            marker = list(colors = sunburstDF_PINSPlus$colors)) %>%
  layout(
    grid = list(columns =1, rows = 1),
    margin = list(l = 0, r = 0, b = 0, t = 0)
  )
pie_PINSPlus
rm(Pheno_sunburst_PINSPlus, sunburstDF_PINSPlus, sunburst_coloring_PINSPlus, pie_PINSPlus); gc()

# Compare MOVICS PINSPlus to MOVICS consensus
CS_comp_list[["PINSPlus"]]$table = table(clust_annot_pheno2$PINSPlus, 
                                         paste0("CS", clust_annot_pheno2$clust))

CS_comp_list[["PINSPlus"]]$ARI = calculate_ari_index(cluster_df1 = PINSPlus_clust_res %>%
                                                       dplyr::rename(Cluster = PINSPlus) %>%
                                                       mutate(Cluster = gsub("MOVICS_PINSPlus", "", Cluster)),
                                                     cluster_df2 = as.data.frame(consensus$clust.res) %>%
                                                       dplyr::rename(Cluster = clust),
                                                     sample_col = "samID",
                                                     clust_col = "Cluster",
                                                     suffixes = c("_MOVICS_PINSPlus", "_CS"))

CS_comp_list[["PINSPlus"]]$NMI = calculate_nmi_index(cluster_df1 = PINSPlus_clust_res %>%
                                                       dplyr::rename(Cluster = PINSPlus) %>%
                                                       mutate(Cluster = gsub("MOVICS_PINSPlus", "", Cluster)),
                                                     cluster_df2 = as.data.frame(consensus$clust.res) %>%
                                                       dplyr::rename(Cluster = clust),
                                                     sample_col = "samID",
                                                     clust_col = "Cluster",
                                                     suffixes = c("_MOVICS_PINSPlus", "_CS"))

# Print all comparison data
print(CS_comp_list$PINSPlus)

# NEMO #####
# PCA from original matrices ###
NEMO_clust_res = clust_annot_pheno %>% dplyr::select(samID, NEMO) %>%
  mutate(NEMO = paste0("MOVICS_", NEMO))

# RNA
pca_from_original_matrix(mydata = plotdata$RNAseq, 
                         algorithm = "NEMO", 
                         clust_res = NEMO_clust_res,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/NEMO_extra"),
                         title_add = "RNA-seq")

# miRNA
pca_from_original_matrix(mydata = plotdata$miRNA, 
                         algorithm = "NEMO", 
                         clust_res = NEMO_clust_res,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/NEMO_extra"),
                         title_add = "miRNA")

# CNV
pca_from_original_matrix(mydata = plotdata$CNV, 
                         algorithm = "NEMO", 
                         clust_res = NEMO_clust_res,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/NEMO_extra"),
                         title_add = "CNV")

# Use multidimensional scaling for SNPs
# Features must be in rows
mds_from_original_matrix(matrix = plotdata$SNPs, dist_method = "binary",
                         algorithm = "NEMO", 
                         clust_res = NEMO_clust_res,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/NEMO_extra"),
                         title_add = "SNPs")

# Methylation
pca_from_original_matrix(mydata = plotdata$Methylation, 
                         algorithm = "NEMO", 
                         clust_res = NEMO_clust_res,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/NEMO_extra"),
                         title_add = "Methylation")

# Bar charts with clinical variables of interest ###
NEMO_barcharts = list()
plotdata_bar = clust_annot_pheno_nonas %>%
  dplyr::mutate(NEMO = paste0("MOVICS_", NEMO))
plotdata_bar$NEMO = factor(plotdata_bar$NEMO)
for (i in 1:length(voi)) {
  chifit = chisq_outputs[["NEMO"]]
  loc = which(grepl(voi[i], chifit$Comparison))
  chifit = chifit[loc, ]
  NEMO_barcharts[[i]] = create_annot_barchart(plotdata = plotdata_bar, fill = voi[i],
                                              chifit = chifit,
                                              na.action = "na.omit",
                                              algorithm = "NEMO",
                                              barchart_ylim = 650,
                                              text_y = 630, rect_ymin = 530,
                                              rect_ymax = 650, x_annot = 1.5,
                                              v_gap = 35, rect_xmin = 1,
                                              rect_xmax = 2, 
                                              annot_text_size = 2.25,
                                              legend.text.size = 5,
                                              x.axis.text.size = 5) +
    barchart_scales[[voi[i]]]
  print(NEMO_barcharts[[i]])
  ggsave(filename = paste0("NEMO_", voi[i], "_barchart.png"),
         path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/NEMO_extra"), 
         width = 2320, height = 2320, device = 'png', units = "px",
         dpi = 700)
  dev.off()
}
names(NEMO_barcharts) = voi
rm(loc, chifit)

# Multiplot (PNG) - bar charts
ggarrange(NEMO_barcharts[[1]], NEMO_barcharts[[2]], NEMO_barcharts[[3]],
          NEMO_barcharts[[4]], NEMO_barcharts[[5]], NEMO_barcharts[[6]],
          NEMO_barcharts[[7]], NEMO_barcharts[[8]], NEMO_barcharts[[9]],
          NEMO_barcharts[[10]], NEMO_barcharts[[11]],
          ncol = 3, nrow = 4, labels = c("A", "B", "C", "D", "E", "F", "G", "H",
                                         "I", "J", "K"),
          font.label = list(size = 8, face = "bold", color ="black"))
ggsave(filename = "Multiplot_NEMO_barcharts.png",
       path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/NEMO_extra"), 
       width = 7000, height = 8000, device = 'png', units = "px",
       dpi = 700)
dev.off()

# Just significant ones now
NEMO_barcharts_sig = list()
plotdata_bar_sig = clust_annot_pheno_nonas %>% dplyr::select(NEMO, Race, Histology, 
                                                       `ER status`, `PR status`, `Stage`) %>%
  dplyr::mutate(NEMO = paste0("MOVICS_", NEMO))
plotdata_bar_sig$NEMO = factor(plotdata_bar_sig$NEMO)
voi_sig = setdiff(colnames(plotdata_bar_sig), "NEMO")
for (i in 1:length(voi_sig)) {
  chifit = chisq_outputs[["NEMO"]]
  loc = which(grepl(voi_sig[i], chifit$Comparison))
  chifit = chifit[loc, ]
  NEMO_barcharts_sig[[i]] = create_annot_barchart(plotdata = plotdata_bar_sig, fill = voi_sig[i],
                                                  chifit = chifit,
                                                  na.action = "na.omit",
                                                  algorithm = "NEMO",
                                                  barchart_ylim = 650,
                                                  text_y = 630, rect_ymin = 530,
                                                  rect_ymax = 650, x_annot = 1.5,
                                                  v_gap = 35, rect_xmin = 1,
                                                  rect_xmax = 2, 
                                                  annot_text_size = 2.25,
                                                  legend.text.size = 5,
                                                  x.axis.text.size = 5) +
    barchart_scales[[voi_sig[i]]]
  print(NEMO_barcharts_sig[[i]])
  ggsave(filename = paste0("sig_NEMO_", voi_sig[i], "_barchart.png"),
         path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/NEMO_extra"), 
         width = 2320, height = 2320, device = 'png', units = "px",
         dpi = 700)
  dev.off()
}
names(NEMO_barcharts_sig) = voi_sig
rm(loc, chifit)

# Multiplot (PNG) - bar charts
ggarrange(NEMO_barcharts_sig[[1]], NEMO_barcharts_sig[[2]], NEMO_barcharts_sig[[3]],
          NEMO_barcharts_sig[[4]], NEMO_barcharts_sig[[5]], 
          ncol = 2, nrow = 3, labels = c("A", "B", "C", "D", "E"),
          font.label = list(size = 8, face = "bold", color ="black"))
ggsave(filename = "sig_Multiplot_NEMO_barcharts.png",
       path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/NEMO_extra"), 
       width = 5500, height = 7000, device = 'png', units = "px",
       dpi = 700)
dev.off()

# Sunburst plot ###
Pheno_sunburst_NEMO = clust_annot_pheno
Pheno_sunburst_NEMO$`ER status` = gsub("Unknown", "Unkn ER status", Pheno_sunburst_NEMO$`ER status`)
Pheno_sunburst_NEMO$`ER status` = gsub("Positive", "ER+", Pheno_sunburst_NEMO$`ER status`)
Pheno_sunburst_NEMO$`ER status` = gsub("Negative", "ER-", Pheno_sunburst_NEMO$`ER status`)
Pheno_sunburst_NEMO$`HER2 status` = gsub("Unknown", "Unkn HER2 status", 
                                         Pheno_sunburst_NEMO$`HER2 status`)

Pheno_sunburst_NEMO = Pheno_sunburst_NEMO %>%
  dplyr::select(NEMO, `ER status`, Stage) %>%
  group_by(NEMO, `ER status`, Stage) %>%
  summarise(Counts = n()) %>%
  as.data.frame()

sunburst_coloring_NEMO = data.frame(stringsAsFactors = FALSE,
                                    colors = tolower(gplots::col2hex(c("#2EC4B6", "#E71D36", 
                                                                       "#C11D9C", "#0F1682",  "grey40",
                                                                       "#00C9FF", "#099CF5", "#097BF5", 
                                                                       "#0B5684", "grey40"))),
                                    labels = c("NEMO1", "NEMO2",
                                               "ER-", "ER+", "Unkn ER status",
                                               "Stage I", "Stage II", "Stage III",
                                               "Stage IV", "Unkn stage"))

sunburstDF_NEMO = as.sunburstDF(Pheno_sunburst_NEMO, value_column = "Counts", add_root = FALSE) %>%
  inner_join(sunburst_coloring_NEMO, by = "labels")

pie_NEMO = plot_ly() %>%
  add_trace(ids = sunburstDF_NEMO$ids, labels= sunburstDF_NEMO$labels, 
            parents = sunburstDF_NEMO$parents, 
            values= sunburstDF_NEMO$values, type='sunburst', branchvalues = 'total',
            insidetextorientation='radial', maxdepth = 5,
            marker = list(colors = sunburstDF_NEMO$colors)) %>%
  layout(
    grid = list(columns =1, rows = 1),
    margin = list(l = 0, r = 0, b = 0, t = 0)
  )
pie_NEMO
rm(Pheno_sunburst_NEMO, sunburstDF_NEMO, sunburst_coloring_NEMO, pie_NEMO); gc()

# Compare MOVICS NEMO to MOVICS consensus
CS_comp_list[["NEMO"]]$table = table(clust_annot_pheno2$NEMO, 
                                     paste0("CS", clust_annot_pheno2$clust))

CS_comp_list[["NEMO"]]$ARI = calculate_ari_index(cluster_df1 = NEMO_clust_res %>%
                                                   dplyr::rename(Cluster = NEMO) %>%
                                                   mutate(Cluster = gsub("MOVICS_NEMO", "", Cluster)),
                                                 cluster_df2 = as.data.frame(consensus$clust.res) %>%
                                                   dplyr::rename(Cluster = clust),
                                                 sample_col = "samID",
                                                 clust_col = "Cluster",
                                                 suffixes = c("_MOVICS_NEMO", "_CS"))

CS_comp_list[["NEMO"]]$NMI = calculate_nmi_index(cluster_df1 = NEMO_clust_res %>%
                                                   dplyr::rename(Cluster = NEMO) %>%
                                                   mutate(Cluster = gsub("MOVICS_NEMO", "", Cluster)),
                                                 cluster_df2 = as.data.frame(consensus$clust.res) %>%
                                                   dplyr::rename(Cluster = clust),
                                                 sample_col = "samID",
                                                 clust_col = "Cluster",
                                                 suffixes = c("_MOVICS_NEMO", "_CS"))

# Print all comparison data
print(CS_comp_list$NEMO)

# COCA #####
# Get the final Jaccard matrix
coca_jaccard = as.matrix(as.dist(vegan::vegdist(as.matrix(moic.res.list$COCA$fit$moc), method = "jaccard")))

# Final dissimilarity matrix
create_MO_heatmap(matrix = coca_jaccard, algorithm = "COCA", 
                  need.diag.zero = TRUE, 
                  clust_annot_pheno %>%
                    mutate(COCA = paste0("MOVICS_", COCA)) %>%
                    select(samID, all_of(afh_colnames), COCA),
                  afh_colnames = afh_colnames, 
                  colors = colors_heatmap,
                  annColors = annColors,
                  cluster_colors = cluster_colors_heatmap,
                  cluster_rows_flag = FALSE,
                  cluster_cols_flag = FALSE,
                  splits_flag = TRUE,
                  heatmap_title = "COCA final dissimilarity matrix heatmap",
                  legend_title = "Dissimilarity",
                  output_file_name = paste0(home, 
                                            "/Results/MOVICS_baseline/MO_comparisons/COCA_extra/COCA_dissimilarity_heatmap.png"))

# PCA from original matrices ###
COCA_clust_res = clust_annot_pheno %>% dplyr::select(samID, COCA) %>%
  mutate(COCA = paste0("MOVICS_", COCA))

# RNA
pca_from_original_matrix(mydata = plotdata$RNAseq, 
                         algorithm = "COCA", 
                         clust_res = COCA_clust_res,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/COCA_extra"),
                         title_add = "RNA-seq")

# miRNA
pca_from_original_matrix(mydata = plotdata$miRNA, 
                         algorithm = "COCA", 
                         clust_res = COCA_clust_res,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/COCA_extra"),
                         title_add = "miRNA")

# CNV
pca_from_original_matrix(mydata = plotdata$CNV, 
                         algorithm = "COCA", 
                         clust_res = COCA_clust_res,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/COCA_extra"),
                         title_add = "CNV")

# Use multidimensional scaling for SNPs
# Features must be in rows
mds_from_original_matrix(matrix = plotdata$SNPs, dist_method = "binary",
                         algorithm = "COCA", 
                         clust_res = COCA_clust_res,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/COCA_extra"),
                         title_add = "SNPs")

# Methylation
pca_from_original_matrix(mydata = plotdata$Methylation, 
                         algorithm = "COCA", 
                         clust_res = COCA_clust_res,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/COCA_extra"),
                         title_add = "Methylation")

# Bar charts with clinical variables of interest ###
COCA_barcharts = list()
plotdata_bar = clust_annot_pheno_nonas %>%
  dplyr::mutate(COCA = paste0("MOVICS_", COCA))
plotdata_bar$COCA = factor(plotdata_bar$COCA)
for (i in 1:length(voi)) {
  chifit = chisq_outputs[["COCA"]]
  loc = which(grepl(voi[i], chifit$Comparison))
  chifit = chifit[loc, ]
  COCA_barcharts[[i]] = create_annot_barchart(plotdata = plotdata_bar, fill = voi[i],
                                              chifit = chifit,
                                              na.action = "na.omit",
                                              algorithm = "COCA",
                                              barchart_ylim = 700,
                                              text_y = 680, rect_ymin = 580,
                                              rect_ymax = 700, x_annot = 1.5,
                                              v_gap = 35, rect_xmin = 1,
                                              rect_xmax = 2, 
                                              annot_text_size = 2.25,
                                              legend.text.size = 5,
                                              x.axis.text.size = 5) +
    barchart_scales[[voi[i]]]
  print(COCA_barcharts[[i]])
  ggsave(filename = paste0("COCA_", voi[i], "_barchart.png"),
         path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/COCA_extra"), 
         width = 2320, height = 2320, device = 'png', units = "px",
         dpi = 700)
  dev.off()
}
names(COCA_barcharts) = voi
rm(loc, chifit)

# Multiplot (PNG) - bar charts
ggarrange(COCA_barcharts[[1]], COCA_barcharts[[2]], COCA_barcharts[[3]],
          COCA_barcharts[[4]], COCA_barcharts[[5]], COCA_barcharts[[6]],
          COCA_barcharts[[7]], COCA_barcharts[[8]], COCA_barcharts[[9]],
          COCA_barcharts[[10]], COCA_barcharts[[11]],
          ncol = 3, nrow = 4, labels = c("A", "B", "C", "D", "E", "F", "G", "H",
                                         "I", "J", "K"),
          font.label = list(size = 8, face = "bold", color ="black"))
ggsave(filename = "Multiplot_COCA_barcharts.png",
       path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/COCA_extra"), 
       width = 7000, height = 8000, device = 'png', units = "px",
       dpi = 700)
dev.off()

# Only significantly associated variable is HER2

# Sunburst plot ###
Pheno_sunburst_COCA = clust_annot_pheno
Pheno_sunburst_COCA$`HER2 status` = gsub("Unknown", "Unkn HER2 status", 
                                             Pheno_sunburst_COCA$`HER2 status`)
Pheno_sunburst_COCA$`HER2 status` = gsub("Positive", "HER2+", Pheno_sunburst_COCA$`HER2 status`)
Pheno_sunburst_COCA$`HER2 status` = gsub("Negative", "HER2-", Pheno_sunburst_COCA$`HER2 status`)
Pheno_sunburst_COCA = Pheno_sunburst_COCA %>%
  dplyr::select(COCA, `HER2 status`) %>%
  group_by(COCA, `HER2 status`) %>%
  summarise(Counts = n()) %>%
  as.data.frame()

sunburst_coloring_COCA = data.frame(stringsAsFactors = FALSE,
                                    colors = tolower(gplots::col2hex(c("#2EC4B6", "#E71D36", 
                                                                       "#0B9EF8", "#560DA7", "mistyrose1", 
                                                                       "hotpink4", "grey40"))),
                                    labels = c("COCA1", "COCA2",
                                               "HER2-", "HER2+", "Indeterminate",
                                               "Equivocal", "Unkn HER2 status"))

sunburstDF_COCA = as.sunburstDF(Pheno_sunburst_COCA, value_column = "Counts", add_root = FALSE) %>%
  inner_join(sunburst_coloring_COCA, by = "labels")

pie_COCA = plot_ly() %>%
  add_trace(ids = sunburstDF_COCA$ids, labels= sunburstDF_COCA$labels, 
            parents = sunburstDF_COCA$parents, 
            values= sunburstDF_COCA$values, type='sunburst', branchvalues = 'total',
            insidetextorientation='radial', maxdepth = 5,
            marker = list(colors = sunburstDF_COCA$colors)) %>%
  layout(
    grid = list(columns =1, rows = 1),
    margin = list(l = 0, r = 0, b = 0, t = 0)
  )
pie_COCA
rm(Pheno_sunburst_COCA, sunburstDF_COCA, sunburst_coloring_COCA, pie_COCA); gc()

# Compare MOVICS COCA to MOVICS consensus
CS_comp_list[["COCA"]]$table = table(clust_annot_pheno2$COCA, 
                                     paste0("CS", clust_annot_pheno2$clust))

CS_comp_list[["COCA"]]$ARI = calculate_ari_index(cluster_df1 = COCA_clust_res %>%
                                                   dplyr::rename(Cluster = COCA) %>%
                                                   mutate(Cluster = gsub("MOVICS_COCA", "", Cluster)),
                                                 cluster_df2 = as.data.frame(consensus$clust.res) %>%
                                                   dplyr::rename(Cluster = clust),
                                                 sample_col = "samID",
                                                 clust_col = "Cluster",
                                                 suffixes = c("_MOVICS_COCA", "_CS"))

CS_comp_list[["COCA"]]$NMI = calculate_nmi_index(cluster_df1 = COCA_clust_res %>%
                                                   dplyr::rename(Cluster = COCA) %>%
                                                   mutate(Cluster = gsub("MOVICS_COCA", "", Cluster)),
                                                 cluster_df2 = as.data.frame(consensus$clust.res) %>%
                                                   dplyr::rename(Cluster = clust),
                                                 sample_col = "samID",
                                                 clust_col = "Cluster",
                                                 suffixes = c("_MOVICS_COCA", "_CS"))

# Print all comparison data
print(CS_comp_list$COCA)

# MoCluster #####
mocluster_mat_2d = moic.res.list[["MoCluster"]][["fit"]]@fac.scr
dist_mocluster_2d = as.matrix(dist(mocluster_mat_2d, method = "euclidean"))

# Final 2D CPCA Euclidean distance heatmap
create_MO_heatmap(matrix = dist_mocluster_2d, algorithm = "MoCluster", 
                  need.diag.zero = TRUE, 
                  clust_annot_pheno = clust_annot_pheno %>%
                    mutate(COCA = paste0("MOVICS_", MoCluster)) %>%
                    select(samID, all_of(afh_colnames), MoCluster),
                  afh_colnames = afh_colnames, 
                  colors = colors_heatmap,
                  annColors = annColors,
                  cluster_colors = cluster_colors_heatmap,
                  cluster_rows_flag = FALSE,
                  cluster_cols_flag = FALSE,
                  splits_flag = TRUE,
                  heatmap_title = "MoCluster 2D CPCA Euclidean distance heatmap",
                  legend_title = "Euclidean distance",
                  output_file_name = paste0(home, 
                                            "/Results/MOVICS_baseline/MO_comparisons/MoCluster_extra/MoCluster_2D_CPCA_Euclidean_distance_heatmap.png"))

# PCA from original matrices ###
MoCluster_clust_res = clust_annot_pheno %>% dplyr::select(samID, MoCluster) %>%
  mutate(MoCluster = paste0("MOVICS_", MoCluster))

# RNA
pca_from_original_matrix(mydata = plotdata$RNAseq, 
                         algorithm = "MoCluster", 
                         clust_res = MoCluster_clust_res,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/MoCluster_extra"),
                         title_add = "RNA-seq")

# miRNA
pca_from_original_matrix(mydata = plotdata$miRNA, 
                         algorithm = "MoCluster", 
                         clust_res = MoCluster_clust_res,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/MoCluster_extra"),
                         title_add = "miRNA")

# CNV
pca_from_original_matrix(mydata = plotdata$CNV, 
                         algorithm = "MoCluster", 
                         clust_res = MoCluster_clust_res,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/MoCluster_extra"),
                         title_add = "CNV")

# Use multidimensional scaling for SNPs
# Features must be in rows
mds_from_original_matrix(matrix = plotdata$SNPs, dist_method = "binary",
                         algorithm = "MoCluster", 
                         clust_res = MoCluster_clust_res,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/MoCluster_extra"),
                         title_add = "SNPs")

# Methylation
pca_from_original_matrix(mydata = plotdata$Methylation, 
                         algorithm = "MoCluster", 
                         clust_res = MoCluster_clust_res,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/MoCluster_extra"),
                         title_add = "Methylation")

# Bar charts with clinical variables of interest ###
MoCluster_barcharts = list()
plotdata_bar = clust_annot_pheno_nonas %>%
  dplyr::mutate(MoCluster = paste0("MOVICS_", MoCluster))
plotdata_bar$MoCluster = factor(plotdata_bar$MoCluster)
for (i in 1:length(voi)) {
  chifit = chisq_outputs[["MoCluster"]]
  loc = which(grepl(voi[i], chifit$Comparison))
  chifit = chifit[loc, ]
  MoCluster_barcharts[[i]] = create_annot_barchart(plotdata = plotdata_bar, fill = voi[i],
                                                   chifit = chifit,
                                                   na.action = "na.omit",
                                                   algorithm = "MoCluster",
                                                   barchart_ylim = 650,
                                                   text_y = 630, rect_ymin = 530,
                                                   rect_ymax = 650, x_annot = 1.5,
                                                   v_gap = 35, rect_xmin = 1,
                                                   rect_xmax = 2, 
                                                   annot_text_size = 2.25,
                                                   legend.text.size = 5,
                                                   x.axis.text.size = 5) +
    barchart_scales[[voi[i]]]
  print(MoCluster_barcharts[[i]])
  ggsave(filename = paste0("MoCluster_", voi[i], "_barchart.png"),
         path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/MoCluster_extra"), 
         width = 2320, height = 2320, device = 'png', units = "px",
         dpi = 700)
  dev.off()
}
names(MoCluster_barcharts) = voi
rm(loc, chifit)

# Multiplot (PNG) - bar charts
ggarrange(MoCluster_barcharts[[1]], MoCluster_barcharts[[2]], MoCluster_barcharts[[3]],
          MoCluster_barcharts[[4]], MoCluster_barcharts[[5]], MoCluster_barcharts[[6]],
          MoCluster_barcharts[[7]], MoCluster_barcharts[[8]], MoCluster_barcharts[[9]],
          MoCluster_barcharts[[10]], MoCluster_barcharts[[11]],
          ncol = 3, nrow = 4, labels = c("A", "B", "C", "D", "E", "F", "G", "H",
                                         "I", "J", "K"),
          font.label = list(size = 8, face = "bold", color ="black"))
ggsave(filename = "Multiplot_MoCluster_barcharts.png",
       path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/MoCluster_extra"), 
       width = 7000, height = 8000, device = 'png', units = "px",
       dpi = 700)
dev.off()

# Just significant ones now
MoCluster_barcharts_sig = list()
plotdata_bar_sig = clust_annot_pheno_nonas %>% dplyr::select(MoCluster, Race, Histology,
                                                       `ER status`, `PR status`) %>%
  dplyr::mutate(MoCluster = paste0("MOVICS_", MoCluster))
plotdata_bar_sig$MoCluster = factor(plotdata_bar_sig$MoCluster)
voi_sig = setdiff(colnames(plotdata_bar_sig), "MoCluster")
for (i in 1:length(voi_sig)) {
  chifit = chisq_outputs[["MoCluster"]]
  loc = which(grepl(voi_sig[i], chifit$Comparison))
  chifit = chifit[loc, ]
  MoCluster_barcharts_sig[[i]] = create_annot_barchart(plotdata = plotdata_bar_sig, fill = voi_sig[i],
                                                       chifit = chifit,
                                                       na.action = "na.omit",
                                                       algorithm = "MoCluster",
                                                       barchart_ylim = 650,
                                                       text_y = 630, rect_ymin = 530,
                                                       rect_ymax = 650, x_annot = 1.5,
                                                       v_gap = 35, rect_xmin = 1,
                                                       rect_xmax = 2, 
                                                       annot_text_size = 2.25,
                                                       legend.text.size = 5,
                                                       x.axis.text.size = 5) +
    barchart_scales[[voi_sig[i]]]
  print(MoCluster_barcharts_sig[[i]])
  ggsave(filename = paste0("sig_MoCluster_", voi_sig[i], "_barchart.png"),
         path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/MoCluster_extra"), 
         width = 2320, height = 2320, device = 'png', units = "px",
         dpi = 700)
  dev.off()
}
names(MoCluster_barcharts_sig) = voi_sig
rm(loc, chifit)

# Multiplot (PNG) - bar charts
ggarrange(MoCluster_barcharts_sig[[1]], MoCluster_barcharts_sig[[2]], MoCluster_barcharts_sig[[3]],
          MoCluster_barcharts_sig[[4]], MoCluster_barcharts_sig[[5]], MoCluster_barcharts_sig[[6]],
          ncol = 2, nrow = 3, labels = c("A", "B", "C", "D", "E", "F"),
          font.label = list(size = 8, face = "bold", color ="black"))
ggsave(filename = "sig_Multiplot_MoCluster_barcharts.png",
       path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/MoCluster_extra"), 
       width = 5500, height = 7000, device = 'png', units = "px",
       dpi = 700)
dev.off()

# Sunburst plot ###
Pheno_sunburst_MoCluster = clust_annot_pheno
Pheno_sunburst_MoCluster$`ER status` = gsub("Unknown", "Unkn ER status", Pheno_sunburst_MoCluster$`ER status`)
Pheno_sunburst_MoCluster$`ER status` = gsub("Positive", "ER+", Pheno_sunburst_MoCluster$`ER status`)
Pheno_sunburst_MoCluster$`ER status` = gsub("Negative", "ER-", Pheno_sunburst_MoCluster$`ER status`)

Pheno_sunburst_MoCluster = Pheno_sunburst_MoCluster %>%
  dplyr::select(MoCluster, `ER status`) %>%
  group_by(MoCluster, `ER status`) %>%
  summarise(Counts = n()) %>%
  as.data.frame()

sunburst_coloring_MoCluster = data.frame(stringsAsFactors = FALSE,
                                         colors = tolower(gplots::col2hex(c("#2EC4B6", "#E71D36", 
                                                                            "#C11D9C", "#0F1682",  "grey40"))),
                                         labels = c("MoCluster1", "MoCluster2",
                                                    "ER-", "ER+", "Unkn ER status"))

sunburstDF_MoCluster = as.sunburstDF(Pheno_sunburst_MoCluster, value_column = "Counts", add_root = FALSE) %>%
  inner_join(sunburst_coloring_MoCluster, by = "labels")

pie_MoCluster = plot_ly() %>%
  add_trace(ids = sunburstDF_MoCluster$ids, labels= sunburstDF_MoCluster$labels, 
            parents = sunburstDF_MoCluster$parents, 
            values= sunburstDF_MoCluster$values, type='sunburst', branchvalues = 'total',
            insidetextorientation='radial', maxdepth = 5,
            marker = list(colors = sunburstDF_MoCluster$colors)) %>%
  layout(
    grid = list(columns =1, rows = 1),
    margin = list(l = 0, r = 0, b = 0, t = 0)
  )
pie_MoCluster
rm(Pheno_sunburst_MoCluster, sunburstDF_MoCluster, sunburst_coloring_MoCluster, pie_MoCluster); gc()

# Compare MOVICS MoCluster to MOVICS consensus
CS_comp_list[["MoCluster"]]$table = table(clust_annot_pheno2$MoCluster, 
                                          paste0("CS", clust_annot_pheno2$clust))

CS_comp_list[["MoCluster"]]$ARI = calculate_ari_index(cluster_df1 = MoCluster_clust_res %>%
                                                        dplyr::rename(Cluster = MoCluster) %>%
                                                        mutate(Cluster = gsub("MOVICS_MoCluster", "", Cluster)),
                                                      cluster_df2 = as.data.frame(consensus$clust.res) %>%
                                                        dplyr::rename(Cluster = clust),
                                                      sample_col = "samID",
                                                      clust_col = "Cluster",
                                                      suffixes = c("_MOVICS_MoCluster", "_CS"))

CS_comp_list[["MoCluster"]]$NMI = calculate_nmi_index(cluster_df1 = MoCluster_clust_res %>%
                                                        dplyr::rename(Cluster = MoCluster) %>%
                                                        mutate(Cluster = gsub("MOVICS_MoCluster", "", Cluster)),
                                                      cluster_df2 = as.data.frame(consensus$clust.res) %>%
                                                        dplyr::rename(Cluster = clust),
                                                      sample_col = "samID",
                                                      clust_col = "Cluster",
                                                      suffixes = c("_MOVICS_MoCluster", "_CS"))

# Print all comparison data
print(CS_comp_list$MoCluster)

# LRAcluster #####
LRA_ld_coordinates = t(moic.res.list[["LRAcluster"]][["fit"]][["coordinate"]])
dist_LRA_2d = as.matrix(dist(LRA_ld_coordinates, method = "euclidean"))

# Final 2D Euclidean distance heatmap
create_MO_heatmap(matrix = dist_LRA_2d, algorithm = "LRAcluster", 
                  need.diag.zero = TRUE,
                  clust_annot_pheno = clust_annot_pheno %>%
                    mutate(LRAcluster = paste0("MOVICS_", LRAcluster)) %>%
                    select(samID, all_of(afh_colnames), LRAcluster),
                  afh_colnames = afh_colnames, 
                  colors = colors_heatmap,
                  annColors = annColors,
                  cluster_colors = cluster_colors_heatmap,
                  cluster_rows_flag = FALSE,
                  cluster_cols_flag = FALSE,
                  splits_flag = TRUE,
                  heatmap_title = "LRAcluster 2D Euclidean distance heatmap",
                  legend_title = "Euclidean distance",
                  output_file_name = paste0(home, 
                                            "/Results/MOVICS_baseline/MO_comparisons/LRAcluster_extra/LRAcluster_2D_Euclidean_distance_heatmap.png"))

# PCA from original matrices ###
LRAcluster_clust_res = clust_annot_pheno %>% dplyr::select(samID, LRAcluster) %>%
  mutate(LRAcluster = paste0("MOVICS_", LRAcluster))

# RNA
pca_from_original_matrix(mydata = plotdata$RNAseq, 
                         algorithm = "LRAcluster", 
                         clust_res = LRAcluster_clust_res,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/LRAcluster_extra"),
                         title_add = "RNA-seq")

# miRNA
pca_from_original_matrix(mydata = plotdata$miRNA, 
                         algorithm = "LRAcluster", 
                         clust_res = LRAcluster_clust_res,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/LRAcluster_extra"),
                         title_add = "miRNA")

# CNV
pca_from_original_matrix(mydata = plotdata$CNV, 
                         algorithm = "LRAcluster", 
                         clust_res = LRAcluster_clust_res,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/LRAcluster_extra"),
                         title_add = "CNV")

# Use multidimensional scaling for SNPs
# Features must be in rows
mds_from_original_matrix(matrix = plotdata$SNPs, dist_method = "binary",
                         algorithm = "LRAcluster", 
                         clust_res = LRAcluster_clust_res,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/LRAcluster_extra"),
                         title_add = "SNPs")

# Methylation
pca_from_original_matrix(mydata = plotdata$Methylation, 
                         algorithm = "LRAcluster", 
                         clust_res = LRAcluster_clust_res,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/LRAcluster_extra"),
                         title_add = "Methylation")

# Bar charts with clinical variables of interest ###
LRAcluster_barcharts = list()
plotdata_bar = clust_annot_pheno_nonas %>%
  dplyr::mutate(LRAcluster = paste0("MOVICS_", LRAcluster))
plotdata_bar$LRAcluster = factor(plotdata_bar$LRAcluster)
for (i in 1:length(voi)) {
  chifit = chisq_outputs[["LRAcluster"]]
  loc = which(grepl(voi[i], chifit$Comparison))
  chifit = chifit[loc, ]
  LRAcluster_barcharts[[i]] = create_annot_barchart(plotdata = plotdata_bar, fill = voi[i],
                                                    chifit = chifit,
                                                    na.action = "na.omit",
                                                    algorithm = "LRAcluster",
                                                    barchart_ylim = 650,
                                                    text_y = 630, rect_ymin = 530,
                                                    rect_ymax = 650, x_annot = 1.5,
                                                    v_gap = 35, rect_xmin = 1,
                                                    rect_xmax = 2, 
                                                    annot_text_size = 2.25,
                                                    legend.text.size = 5,
                                                    x.axis.text.size = 5) +
    barchart_scales[[voi[i]]]
  print(LRAcluster_barcharts[[i]])
  ggsave(filename = paste0("LRAcluster_", voi[i], "_barchart.png"),
         path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/LRAcluster_extra"), 
         width = 2320, height = 2320, device = 'png', units = "px",
         dpi = 700)
  dev.off()
}
names(LRAcluster_barcharts) = voi
rm(loc, chifit)

# Multiplot (PNG) - bar charts
ggarrange(LRAcluster_barcharts[[1]], LRAcluster_barcharts[[2]], LRAcluster_barcharts[[3]],
          LRAcluster_barcharts[[4]], LRAcluster_barcharts[[5]], LRAcluster_barcharts[[6]],
          LRAcluster_barcharts[[7]], LRAcluster_barcharts[[8]], LRAcluster_barcharts[[9]],
          LRAcluster_barcharts[[10]], LRAcluster_barcharts[[11]],
          ncol = 3, nrow = 4, labels = c("A", "B", "C", "D", "E", "F", "G", "H",
                                         "I", "J", "K"),
          font.label = list(size = 8, face = "bold", color ="black"))
ggsave(filename = "Multiplot_LRAcluster_barcharts.png",
       path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/LRAcluster_extra"), 
       width = 7000, height = 8000, device = 'png', units = "px",
       dpi = 700)
dev.off()

# Just significant ones now
LRAcluster_barcharts_sig = list()
plotdata_bar_sig = clust_annot_pheno_nonas %>% dplyr::select(LRAcluster, `ER status`, `PR status`, 
                                                       `HER2 status`, Histology, Stage) %>%
  dplyr::mutate(LRAcluster = paste0("MOVICS_", LRAcluster))
plotdata_bar_sig$LRAcluster = factor(plotdata_bar_sig$LRAcluster)
voi_sig = setdiff(colnames(plotdata_bar_sig), "LRAcluster")
for (i in 1:length(voi_sig)) {
  chifit = chisq_outputs[["LRAcluster"]]
  loc = which(grepl(voi_sig[i], chifit$Comparison))
  chifit = chifit[loc, ]
  LRAcluster_barcharts_sig[[i]] = create_annot_barchart(plotdata = plotdata_bar_sig, fill = voi_sig[i],
                                                        chifit = chifit,
                                                        na.action = "na.omit",
                                                        algorithm = "LRAcluster",
                                                        barchart_ylim = 650,
                                                        text_y = 630, rect_ymin = 530,
                                                        rect_ymax = 650, x_annot = 1.5,
                                                        v_gap = 35, rect_xmin = 1,
                                                        rect_xmax = 2, 
                                                        annot_text_size = 2.25,
                                                        legend.text.size = 5,
                                                        x.axis.text.size = 5) +
    barchart_scales[[voi_sig[i]]]
  print(LRAcluster_barcharts_sig[[i]])
  ggsave(filename = paste0("sig_LRAcluster_", voi_sig[i], "_barchart.png"),
         path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/LRAcluster_extra"), 
         width = 2320, height = 2320, device = 'png', units = "px",
         dpi = 700)
  dev.off()
}
names(LRAcluster_barcharts_sig) = voi_sig
rm(loc, chifit)

# Multiplot (PNG) - bar charts
ggarrange(LRAcluster_barcharts_sig[[1]], LRAcluster_barcharts_sig[[2]], LRAcluster_barcharts_sig[[3]],
          LRAcluster_barcharts_sig[[4]], LRAcluster_barcharts_sig[[5]],
          ncol = 2, nrow = 3, labels = c("A", "B", "C", "D", "E"),
          font.label = list(size = 8, face = "bold", color ="black"))
ggsave(filename = "sig_Multiplot_LRAcluster_barcharts.png",
       path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/LRAcluster_extra"), 
       width = 5500, height = 8000, device = 'png', units = "px",
       dpi = 700)
dev.off()

# Sunburst plot ###
Pheno_sunburst_LRAcluster = clust_annot_pheno
Pheno_sunburst_LRAcluster$`ER status` = gsub("Unknown", "Unkn ER status", Pheno_sunburst_LRAcluster$`ER status`)
Pheno_sunburst_LRAcluster$`ER status` = gsub("Positive", "ER+", Pheno_sunburst_LRAcluster$`ER status`)
Pheno_sunburst_LRAcluster$`ER status` = gsub("Negative", "ER-", Pheno_sunburst_LRAcluster$`ER status`)
Pheno_sunburst_LRAcluster$`HER2 status` = gsub("Unknown", "Unkn HER2 status", 
                                               Pheno_sunburst_LRAcluster$`HER2 status`)
Pheno_sunburst_LRAcluster$`HER2 status` = gsub("Positive", "HER2+", Pheno_sunburst_LRAcluster$`HER2 status`)
Pheno_sunburst_LRAcluster$`HER2 status` = gsub("Negative", "HER2-", Pheno_sunburst_LRAcluster$`HER2 status`)
Pheno_sunburst_LRAcluster = Pheno_sunburst_LRAcluster %>%
  dplyr::select(LRAcluster, `ER status`, `HER2 status`) %>%
  group_by(LRAcluster, `ER status`, `HER2 status`) %>%
  summarise(Counts = n()) %>%
  as.data.frame()

sunburst_coloring_LRAcluster = data.frame(stringsAsFactors = FALSE,
                                          colors = tolower(gplots::col2hex(c("#2EC4B6", "#E71D36", 
                                                                             "#C11D9C", "#0F1682",  "grey40",
                                                                             "#0B9EF8", "#560DA7", "mistyrose1", 
                                                                             "hotpink4", "grey40"))),
                                          labels = c("LRAcluster1", "LRAcluster2",
                                                     "ER-", "ER+", "Unkn ER status",
                                                     "HER2-", "HER2+", "Indeterminate",
                                                     "Equivocal", "Unkn HER2 status"))

sunburstDF_LRAcluster = as.sunburstDF(Pheno_sunburst_LRAcluster, value_column = "Counts", add_root = FALSE) %>%
  inner_join(sunburst_coloring_LRAcluster, by = "labels")

pie_LRAcluster = plot_ly() %>%
  add_trace(ids = sunburstDF_LRAcluster$ids, labels= sunburstDF_LRAcluster$labels, 
            parents = sunburstDF_LRAcluster$parents, 
            values= sunburstDF_LRAcluster$values, type='sunburst', branchvalues = 'total',
            insidetextorientation='radial', maxdepth = 5,
            marker = list(colors = sunburstDF_LRAcluster$colors)) %>%
  layout(
    grid = list(columns =1, rows = 1),
    margin = list(l = 0, r = 0, b = 0, t = 0)
  )
pie_LRAcluster
rm(Pheno_sunburst_LRAcluster, sunburstDF_LRAcluster, sunburst_coloring_LRAcluster, pie_LRAcluster); gc()

# Compare MOVICS LRAcluster to MOVICS consensus
CS_comp_list[["LRAcluster"]]$table = table(clust_annot_pheno2$LRAcluster, 
                                           paste0("CS", clust_annot_pheno2$clust))

CS_comp_list[["LRAcluster"]]$ARI = calculate_ari_index(cluster_df1 = LRAcluster_clust_res %>%
                                                         dplyr::rename(Cluster = LRAcluster) %>%
                                                         mutate(Cluster = gsub("MOVICS_LRAcluster", "", Cluster)),
                                                       cluster_df2 = as.data.frame(consensus$clust.res) %>%
                                                         dplyr::rename(Cluster = clust),
                                                       sample_col = "samID",
                                                       clust_col = "Cluster",
                                                       suffixes = c("_MOVICS_LRAcluster", "_CS"))

CS_comp_list[["LRAcluster"]]$NMI = calculate_nmi_index(cluster_df1 = LRAcluster_clust_res %>%
                                                         dplyr::rename(Cluster = LRAcluster) %>%
                                                         mutate(Cluster = gsub("MOVICS_LRAcluster", "", Cluster)),
                                                       cluster_df2 = as.data.frame(consensus$clust.res) %>%
                                                         dplyr::rename(Cluster = clust),
                                                       sample_col = "samID",
                                                       clust_col = "Cluster",
                                                       suffixes = c("_MOVICS_LRAcluster", "_CS"))

# Print all comparison data
print(CS_comp_list$LRAcluster)

# Consensus Clustering #####
final_cc_matrix = 1 - moic.res.list[["ConsensusClustering"]][["fit"]][[2]][["consensusMatrix"]]
dimnames(final_cc_matrix) = list(names(moic.res.list[["ConsensusClustering"]][["fit"]][[2]][["consensusClass"]]),
                                 names(moic.res.list[["ConsensusClustering"]][["fit"]][[2]][["consensusClass"]]))

# Final 2D Euclidean distance heatmap
create_MO_heatmap(matrix = final_cc_matrix, algorithm = "ConsensusClustering", 
                  need.diag.zero = TRUE, 
                  clust_annot_pheno = clust_annot_pheno %>%
                    mutate(ConsensusClustering = paste0("MOVICS_", ConsensusClustering)) %>%
                    select(samID, all_of(afh_colnames), ConsensusClustering),
                  afh_colnames = afh_colnames, 
                  colors = colors_heatmap,
                  annColors = annColors,
                  cluster_colors = cluster_colors_heatmap,
                  cluster_rows_flag = FALSE,
                  cluster_cols_flag = FALSE,
                  splits_flag = TRUE,
                  heatmap_title = "Consensus Clustering final connectivity heatmap",
                  legend_title = "1 - final connectivity",
                  output_file_name = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/CC_extra/CC_final_connectivity_heatmap.png"))

# PCA from original matrices ###
ConsensusClustering_clust_res = clust_annot_pheno %>% dplyr::select(samID, ConsensusClustering) %>%
  mutate(ConsensusClustering = paste0("MOVICS_", ConsensusClustering))

# RNA
pca_from_original_matrix(mydata = plotdata$RNAseq, 
                         algorithm = "ConsensusClustering", 
                         clust_res = ConsensusClustering_clust_res,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/CC_extra"),
                         title_add = "RNA-seq")

# miRNA
pca_from_original_matrix(mydata = plotdata$miRNA, 
                         algorithm = "ConsensusClustering", 
                         clust_res = ConsensusClustering_clust_res,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/CC_extra"),
                         title_add = "miRNA")

# CNV
pca_from_original_matrix(mydata = plotdata$CNV, 
                         algorithm = "ConsensusClustering", 
                         clust_res = ConsensusClustering_clust_res,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/CC_extra"),
                         title_add = "CNV")

# Use multidimensional scaling for SNPs
# Features must be in rows
mds_from_original_matrix(matrix = plotdata$SNPs, dist_method = "binary",
                         algorithm = "ConsensusClustering", 
                         clust_res = ConsensusClustering_clust_res,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/CC_extra"),
                         title_add = "SNPs")

# Methylation
pca_from_original_matrix(mydata = plotdata$Methylation, 
                         algorithm = "ConsensusClustering", 
                         clust_res = ConsensusClustering_clust_res,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/CC_extra"),
                         title_add = "Methylation")

# Bar charts with clinical variables of interest ###
ConsensusClustering_barcharts = list()
plotdata_bar = clust_annot_pheno_nonas %>%
  dplyr::mutate(ConsensusClustering = gsub("ConsensusClustering", "CC", ConsensusClustering)) %>%
  dplyr::mutate(ConsensusClustering = paste0("MOVICS_", ConsensusClustering))
plotdata_bar$ConsensusClustering = factor(plotdata_bar$ConsensusClustering)
for (i in 1:length(voi)) {
  chifit = chisq_outputs[["ConsensusClustering"]]
  loc = which(grepl(voi[i], chifit$Comparison))
  chifit = chifit[loc, ]
  ConsensusClustering_barcharts[[i]] = create_annot_barchart(plotdata = plotdata_bar, fill = voi[i],
                                                             chifit = chifit,
                                                             na.action = "na.omit",
                                                             algorithm = "ConsensusClustering",
                                                             barchart_ylim = 650,
                                                             text_y = 630, rect_ymin = 530,
                                                             rect_ymax = 650, x_annot = 1.5,
                                                             v_gap = 35, rect_xmin = 1,
                                                             rect_xmax = 2, 
                                                             annot_text_size = 2.25,
                                                             legend.text.size = 5,
                                                             x.axis.text.size = 5) +
    barchart_scales[[voi[i]]]
  print(ConsensusClustering_barcharts[[i]])
  ggsave(filename = paste0("ConsensusClustering_", voi[i], "_barchart.png"),
         path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/CC_extra"), 
         width = 2320, height = 2320, device = 'png', units = "px",
         dpi = 700)
  dev.off()
}
names(ConsensusClustering_barcharts) = voi
rm(loc, chifit)

# Multiplot (PNG) - bar charts
ggarrange(ConsensusClustering_barcharts[[1]], ConsensusClustering_barcharts[[2]], ConsensusClustering_barcharts[[3]],
          ConsensusClustering_barcharts[[4]], ConsensusClustering_barcharts[[5]], ConsensusClustering_barcharts[[6]],
          ConsensusClustering_barcharts[[7]], ConsensusClustering_barcharts[[8]], ConsensusClustering_barcharts[[9]],
          ConsensusClustering_barcharts[[10]], ConsensusClustering_barcharts[[11]],
          ncol = 3, nrow = 4, labels = c("A", "B", "C", "D", "E", "F", "G", "H",
                                         "I", "J", "K"),
          font.label = list(size = 8, face = "bold", color ="black"))
ggsave(filename = "Multiplot_ConsensusClustering_barcharts.png",
       path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/CC_extra"), 
       width = 7000, height = 8000, device = 'png', units = "px",
       dpi = 700)
dev.off()

# Just significant ones now
ConsensusClustering_barcharts_sig = list()
plotdata_bar_sig = clust_annot_pheno_nonas %>% dplyr::select(ConsensusClustering, `ER status`, 
                                                       `PR status`, `HER2 status`, Histology) %>%
  dplyr::mutate(ConsensusClustering = gsub("ConsensusClustering", "CC", ConsensusClustering)) %>%
  dplyr::mutate(ConsensusClustering = paste0("MOVICS_", ConsensusClustering))
plotdata_bar_sig$ConsensusClustering = factor(plotdata_bar_sig$ConsensusClustering)
voi_sig = setdiff(colnames(plotdata_bar_sig), "ConsensusClustering")
for (i in 1:length(voi_sig)) {
  chifit = chisq_outputs[["ConsensusClustering"]]
  loc = which(grepl(voi_sig[i], chifit$Comparison))
  chifit = chifit[loc, ]
  ConsensusClustering_barcharts_sig[[i]] = create_annot_barchart(plotdata = plotdata_bar_sig, fill = voi_sig[i],
                                                                 chifit = chifit,
                                                                 na.action = "na.omit",
                                                                 algorithm = "ConsensusClustering",
                                                                 barchart_ylim = 650,
                                                                 text_y = 630, rect_ymin = 530,
                                                                 rect_ymax = 650, x_annot = 1.5,
                                                                 v_gap = 35, rect_xmin = 1,
                                                                 rect_xmax = 2, 
                                                                 annot_text_size = 2.25,
                                                                 legend.text.size = 5,
                                                                 x.axis.text.size = 5) +
    barchart_scales[[voi_sig[i]]]
  print(ConsensusClustering_barcharts_sig[[i]])
  ggsave(filename = paste0("sig_ConsensusClustering_", voi_sig[i], "_barchart.png"),
         path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/CC_extra"), 
         width = 2320, height = 2320, device = 'png', units = "px",
         dpi = 700)
  dev.off()
}
names(ConsensusClustering_barcharts_sig) = voi_sig
rm(loc, chifit)

# Multiplot (PNG) - bar charts
ggarrange(ConsensusClustering_barcharts_sig[[1]], ConsensusClustering_barcharts_sig[[2]], ConsensusClustering_barcharts_sig[[3]],
          ConsensusClustering_barcharts_sig[[4]],  
          ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"),
          font.label = list(size = 8, face = "bold", color ="black"))
ggsave(filename = "sig_Multiplot_ConsensusClustering_barcharts.png",
       path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/CC_extra"), 
       width = 5500, height = 5500, device = 'png', units = "px",
       dpi = 700)
dev.off()

# Sunburst plot ###
Pheno_sunburst_ConsensusClustering = clust_annot_pheno
Pheno_sunburst_ConsensusClustering$`ER status` = gsub("Unknown", "Unkn ER status", Pheno_sunburst_ConsensusClustering$`ER status`)
Pheno_sunburst_ConsensusClustering$`ER status` = gsub("Positive", "ER+", Pheno_sunburst_ConsensusClustering$`ER status`)
Pheno_sunburst_ConsensusClustering$`ER status` = gsub("Negative", "ER-", Pheno_sunburst_ConsensusClustering$`ER status`)
Pheno_sunburst_ConsensusClustering$`HER2 status` = gsub("Unknown", "Unkn HER2 status", 
                                                        Pheno_sunburst_ConsensusClustering$`HER2 status`)
Pheno_sunburst_ConsensusClustering$`HER2 status` = gsub("Positive", "HER2+", Pheno_sunburst_ConsensusClustering$`HER2 status`)
Pheno_sunburst_ConsensusClustering$`HER2 status` = gsub("Negative", "HER2-", Pheno_sunburst_ConsensusClustering$`HER2 status`)
Pheno_sunburst_ConsensusClustering = Pheno_sunburst_ConsensusClustering %>%
  dplyr::mutate(ConsensusClustering = gsub("ConsensusClustering", "CC", ConsensusClustering)) %>%
  dplyr::select(ConsensusClustering, `ER status`, `HER2 status`) %>%
  group_by(ConsensusClustering, `ER status`, `HER2 status`) %>%
  summarise(Counts = n()) %>%
  as.data.frame()

sunburst_coloring_ConsensusClustering = data.frame(stringsAsFactors = FALSE,
                                                   colors = tolower(gplots::col2hex(c("#2EC4B6", "#E71D36", 
                                                                                      "#C11D9C", "#0F1682",  "grey40",
                                                                                      "#0B9EF8", "#560DA7", "mistyrose1", 
                                                                                      "hotpink4", "grey40"))),
                                                   labels = c("CC1", "CC2",
                                                              "ER-", "ER+", "Unkn ER status",
                                                              "HER2-", "HER2+", "Indeterminate",
                                                              "Equivocal", "Unkn HER2 status"))

sunburstDF_ConsensusClustering = as.sunburstDF(Pheno_sunburst_ConsensusClustering, value_column = "Counts", add_root = FALSE) %>%
  inner_join(sunburst_coloring_ConsensusClustering, by = "labels")

pie_ConsensusClustering = plot_ly() %>%
  add_trace(ids = sunburstDF_ConsensusClustering$ids, labels= sunburstDF_ConsensusClustering$labels, 
            parents = sunburstDF_ConsensusClustering$parents, 
            values= sunburstDF_ConsensusClustering$values, type='sunburst', branchvalues = 'total',
            insidetextorientation='radial', maxdepth = 5,
            marker = list(colors = sunburstDF_ConsensusClustering$colors)) %>%
  layout(
    grid = list(columns =1, rows = 1),
    margin = list(l = 0, r = 0, b = 0, t = 0)
  )
pie_ConsensusClustering
rm(Pheno_sunburst_ConsensusClustering, sunburstDF_ConsensusClustering, 
   sunburst_coloring_ConsensusClustering, pie_ConsensusClustering); gc()

# Compare MOVICS ConsensusClustering to MOVICS consensus
CS_comp_list[["ConsensusClustering"]]$table = table(clust_annot_pheno2$ConsensusClustering, 
                                                    paste0("CS", clust_annot_pheno2$clust))

CS_comp_list[["ConsensusClustering"]]$ARI = calculate_ari_index(cluster_df1 = ConsensusClustering_clust_res %>%
                                                                  dplyr::rename(Cluster = ConsensusClustering) %>%
                                                                  mutate(Cluster = gsub("MOVICS_ConsensusClustering", "", Cluster)),
                                                                cluster_df2 = as.data.frame(consensus$clust.res) %>%
                                                                  dplyr::rename(Cluster = clust),
                                                                sample_col = "samID",
                                                                clust_col = "Cluster",
                                                                suffixes = c("_MOVICS_ConsensusClustering", "_CS"))

CS_comp_list[["ConsensusClustering"]]$NMI = calculate_nmi_index(cluster_df1 = ConsensusClustering_clust_res %>%
                                                                  dplyr::rename(Cluster = ConsensusClustering) %>%
                                                                  mutate(Cluster = gsub("MOVICS_ConsensusClustering", "", Cluster)),
                                                                cluster_df2 = as.data.frame(consensus$clust.res) %>%
                                                                  dplyr::rename(Cluster = clust),
                                                                sample_col = "samID",
                                                                clust_col = "Cluster",
                                                                suffixes = c("_MOVICS_ConsensusClustering", "_CS"))

# Print all comparison data
print(CS_comp_list$ConsensusClustering)

# IntNMF #####
# PCA from original matrices ###
IntNMF_clust_res = clust_annot_pheno %>% dplyr::select(samID, IntNMF) %>%
  mutate(IntNMF = paste0("MOVICS_", IntNMF))

# RNA
pca_from_original_matrix(mydata = plotdata$RNAseq, 
                         algorithm = "IntNMF", 
                         clust_res = IntNMF_clust_res,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/IntNMF_extra"),
                         title_add = "RNA-seq")

# miRNA
pca_from_original_matrix(mydata = plotdata$miRNA, 
                         algorithm = "IntNMF", 
                         clust_res = IntNMF_clust_res,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/IntNMF_extra"),
                         title_add = "miRNA")

# CNV
pca_from_original_matrix(mydata = plotdata$CNV, 
                         algorithm = "IntNMF", 
                         clust_res = IntNMF_clust_res,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/IntNMF_extra"),
                         title_add = "CNV")

# Use multidimensional scaling for SNPs
# Features must be in rows
mds_from_original_matrix(matrix = plotdata$SNPs, dist_method = "binary",
                         algorithm = "IntNMF", 
                         clust_res = IntNMF_clust_res,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/IntNMF_extra"),
                         title_add = "SNPs")

# Methylation
pca_from_original_matrix(mydata = plotdata$Methylation, 
                         algorithm = "IntNMF", 
                         clust_res = IntNMF_clust_res,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/IntNMF_extra"),
                         title_add = "Methylation")

# Bar charts with clinical variables of interest ###
IntNMF_barcharts = list()
plotdata_bar = clust_annot_pheno_nonas %>%
  dplyr::mutate(IntNMF = paste0("MOVICS_", IntNMF))
plotdata_bar$IntNMF = factor(plotdata_bar$IntNMF)
for (i in 1:length(voi)) {
  chifit = chisq_outputs[["IntNMF"]]
  loc = which(grepl(voi[i], chifit$Comparison))
  chifit = chifit[loc, ]
  IntNMF_barcharts[[i]] = create_annot_barchart(plotdata = plotdata_bar, fill = voi[i],
                                                chifit = chifit,
                                                na.action = "na.omit",
                                                algorithm = "IntNMF",
                                                barchart_ylim = 650,
                                                text_y = 630, rect_ymin = 530,
                                                rect_ymax = 650, x_annot = 1.5,
                                                v_gap = 35, rect_xmin = 1,
                                                rect_xmax = 2, 
                                                annot_text_size = 2.25,
                                                legend.text.size = 5,
                                                x.axis.text.size = 5) +
    barchart_scales[[voi[i]]]
  print(IntNMF_barcharts[[i]])
  ggsave(filename = paste0("IntNMF_", voi[i], "_barchart.png"),
         path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/IntNMF_extra"), 
         width = 2320, height = 2320, device = 'png', units = "px",
         dpi = 700)
  dev.off()
}
names(IntNMF_barcharts) = voi
rm(loc, chifit)

# Multiplot (PNG) - bar charts
ggarrange(IntNMF_barcharts[[1]], IntNMF_barcharts[[2]], IntNMF_barcharts[[3]],
          IntNMF_barcharts[[4]], IntNMF_barcharts[[5]], IntNMF_barcharts[[6]],
          IntNMF_barcharts[[7]], IntNMF_barcharts[[8]], IntNMF_barcharts[[9]],
          IntNMF_barcharts[[10]], IntNMF_barcharts[[11]],
          ncol = 3, nrow = 4, labels = c("A", "B", "C", "D", "E", "F", "G", "H",
                                         "I", "J", "K"),
          font.label = list(size = 8, face = "bold", color ="black"))
ggsave(filename = "Multiplot_IntNMF_barcharts.png",
       path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/IntNMF_extra"), 
       width = 7000, height = 8000, device = 'png', units = "px",
       dpi = 700)
dev.off()

# Just significant ones now
IntNMF_barcharts_sig = list()
plotdata_bar_sig = clust_annot_pheno_nonas %>% dplyr::select(IntNMF, Histology, Stage, `Lymph node status`,
                                                       `ER status`, `PR status`, `HER2 status`) %>%
  dplyr::mutate(IntNMF = paste0("MOVICS_", IntNMF))
plotdata_bar_sig$IntNMF = factor(plotdata_bar_sig$IntNMF)
voi_sig = setdiff(colnames(plotdata_bar_sig), "IntNMF")
for (i in 1:length(voi_sig)) {
  chifit = chisq_outputs[["IntNMF"]]
  loc = which(grepl(voi_sig[i], chifit$Comparison))
  chifit = chifit[loc, ]
  IntNMF_barcharts_sig[[i]] = create_annot_barchart(plotdata = plotdata_bar_sig, fill = voi_sig[i],
                                                    chifit = chifit,
                                                    na.action = "na.omit",
                                                    algorithm = "IntNMF",
                                                    barchart_ylim = 650,
                                                    text_y = 630, rect_ymin = 530,
                                                    rect_ymax = 650, x_annot = 1.5,
                                                    v_gap = 35, rect_xmin = 1,
                                                    rect_xmax = 2, 
                                                    annot_text_size = 2.25,
                                                    legend.text.size = 5,
                                                    x.axis.text.size = 5) +
    barchart_scales[[voi_sig[i]]]
  print(IntNMF_barcharts_sig[[i]])
  ggsave(filename = paste0("sig_IntNMF_", voi_sig[i], "_barchart.png"),
         path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/IntNMF_extra"), 
         width = 2320, height = 2320, device = 'png', units = "px",
         dpi = 700)
  dev.off()
}
names(IntNMF_barcharts_sig) = voi_sig
rm(loc, chifit)

# Multiplot (PNG) - bar charts
ggarrange(IntNMF_barcharts_sig[[1]], IntNMF_barcharts_sig[[2]], IntNMF_barcharts_sig[[3]],
          IntNMF_barcharts_sig[[4]], IntNMF_barcharts_sig[[5]], IntNMF_barcharts_sig[[6]],
          ncol = 2, nrow = 3, labels = c("A", "B", "C", "D", "E", "F"),
          font.label = list(size = 8, face = "bold", color ="black"))
ggsave(filename = "sig_Multiplot_IntNMF_barcharts.png",
       path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/IntNMF_extra"), 
       width = 5500, height = 7000, device = 'png', units = "px",
       dpi = 700)
dev.off()

# Sunburst plot ###
Pheno_sunburst_IntNMF = clust_annot_pheno
Pheno_sunburst_IntNMF$`ER status` = gsub("Unknown", "Unkn ER status", Pheno_sunburst_IntNMF$`ER status`)
Pheno_sunburst_IntNMF$`ER status` = gsub("Positive", "ER+", Pheno_sunburst_IntNMF$`ER status`)
Pheno_sunburst_IntNMF$`ER status` = gsub("Negative", "ER-", Pheno_sunburst_IntNMF$`ER status`)
Pheno_sunburst_IntNMF$`HER2 status` = gsub("Unknown", "Unkn HER2 status", 
                                           Pheno_sunburst_IntNMF$`HER2 status`)
Pheno_sunburst_IntNMF$`HER2 status` = gsub("Positive", "HER2+", Pheno_sunburst_IntNMF$`HER2 status`)
Pheno_sunburst_IntNMF$`HER2 status` = gsub("Negative", "HER2-", Pheno_sunburst_IntNMF$`HER2 status`)
Pheno_sunburst_IntNMF = Pheno_sunburst_IntNMF %>%
  dplyr::select(IntNMF, `ER status`, `HER2 status`) %>%
  group_by(IntNMF, `ER status`, `HER2 status`) %>%
  summarise(Counts = n()) %>%
  as.data.frame()

sunburst_coloring_IntNMF = data.frame(stringsAsFactors = FALSE,
                                      colors = tolower(gplots::col2hex(c("#2EC4B6", "#E71D36", 
                                                                         "#C11D9C", "#0F1682",  "grey40",
                                                                         "#0B9EF8", "#560DA7", "mistyrose1", 
                                                                         "hotpink4", "grey40"))),
                                      labels = c("IntNMF1", "IntNMF2",
                                                 "ER-", "ER+", "Unkn ER status",
                                                 "HER2-", "HER2+", "Indeterminate",
                                                 "Equivocal", "Unkn HER2 status"))

sunburstDF_IntNMF = as.sunburstDF(Pheno_sunburst_IntNMF, value_column = "Counts", add_root = FALSE) %>%
  inner_join(sunburst_coloring_IntNMF, by = "labels")

pie_IntNMF = plot_ly() %>%
  add_trace(ids = sunburstDF_IntNMF$ids, labels= sunburstDF_IntNMF$labels, 
            parents = sunburstDF_IntNMF$parents, 
            values= sunburstDF_IntNMF$values, type='sunburst', branchvalues = 'total',
            insidetextorientation='radial', maxdepth = 5,
            marker = list(colors = sunburstDF_IntNMF$colors)) %>%
  layout(
    grid = list(columns =1, rows = 1),
    margin = list(l = 0, r = 0, b = 0, t = 0)
  )
pie_IntNMF
rm(Pheno_sunburst_IntNMF, sunburstDF_IntNMF, sunburst_coloring_IntNMF, pie_IntNMF); gc()

# Compare MOVICS IntNMF to MOVICS consensus
CS_comp_list[["IntNMF"]]$table = table(clust_annot_pheno2$IntNMF, 
                                       paste0("CS", clust_annot_pheno2$clust))

CS_comp_list[["IntNMF"]]$ARI = calculate_ari_index(cluster_df1 = IntNMF_clust_res %>%
                                                     dplyr::rename(Cluster = IntNMF) %>%
                                                     mutate(Cluster = gsub("MOVICS_IntNMF", "", Cluster)),
                                                   cluster_df2 = as.data.frame(consensus$clust.res) %>%
                                                     dplyr::rename(Cluster = clust),
                                                   sample_col = "samID",
                                                   clust_col = "Cluster",
                                                   suffixes = c("_MOVICS_IntNMF", "_CS"))

CS_comp_list[["IntNMF"]]$NMI = calculate_nmi_index(cluster_df1 = IntNMF_clust_res %>%
                                                     dplyr::rename(Cluster = IntNMF) %>%
                                                     mutate(Cluster = gsub("MOVICS_IntNMF", "", Cluster)),
                                                   cluster_df2 = as.data.frame(consensus$clust.res) %>%
                                                     dplyr::rename(Cluster = clust),
                                                   sample_col = "samID",
                                                   clust_col = "Cluster",
                                                   suffixes = c("_MOVICS_IntNMF", "_CS"))

# Print all comparison data
print(CS_comp_list$IntNMF)

# iClusterBayes #####
iCB_feature_ranks = as.data.frame(moic.res.list[["iClusterBayes"]][["feat.res"]])
write.xlsx(iCB_feature_ranks,
           paste0(home, "/Results/MOVICS_baseline/MO_comparisons/iClusterBayes_extra/iCB_feature_ranks.xlsx"),
           overwrite = TRUE)

# PCA from original matrices ###
iClusterBayes_clust_res = clust_annot_pheno %>% dplyr::select(samID, iClusterBayes) %>%
  mutate(iClusterBayes = paste0("MOVICS_", iClusterBayes))

# RNA
pca_from_original_matrix(mydata = plotdata$RNAseq, 
                         algorithm = "iClusterBayes", 
                         clust_res = iClusterBayes_clust_res,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/iClusterBayes_extra"),
                         title_add = "RNA-seq")

# miRNA
pca_from_original_matrix(mydata = plotdata$miRNA, 
                         algorithm = "iClusterBayes", 
                         clust_res = iClusterBayes_clust_res,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/iClusterBayes_extra"),
                         title_add = "miRNA")

# CNV
pca_from_original_matrix(mydata = plotdata$CNV, 
                         algorithm = "iClusterBayes", 
                         clust_res = iClusterBayes_clust_res,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/iClusterBayes_extra"),
                         title_add = "CNV")

# Use multidimensional scaling for SNPs
# Features must be in rows
mds_from_original_matrix(matrix = plotdata$SNPs, dist_method = "binary",
                         algorithm = "iClusterBayes", 
                         clust_res = iClusterBayes_clust_res,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/iClusterBayes_extra"),
                         title_add = "SNPs")

# Methylation
pca_from_original_matrix(mydata = plotdata$Methylation, 
                         algorithm = "iClusterBayes", 
                         clust_res = iClusterBayes_clust_res,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/iClusterBayes_extra"),
                         title_add = "Methylation")

# Bar charts with clinical variables of interest ###
iClusterBayes_barcharts = list()
plotdata_bar = clust_annot_pheno_nonas %>%
  dplyr::mutate(iClusterBayes = paste0("MOVICS_", iClusterBayes))
plotdata_bar$iClusterBayes = factor(plotdata_bar$iClusterBayes)
for (i in 1:length(voi)) {
  chifit = chisq_outputs[["iClusterBayes"]]
  loc = which(grepl(voi[i], chifit$Comparison))
  chifit = chifit[loc, ]
  iClusterBayes_barcharts[[i]] = create_annot_barchart(plotdata = plotdata_bar, fill = voi[i],
                                                       chifit = chifit,
                                                       na.action = "na.omit",
                                                       algorithm = "iClusterBayes",
                                                       barchart_ylim = 650,
                                                       text_y = 630, rect_ymin = 530,
                                                       rect_ymax = 650, x_annot = 1.5,
                                                       v_gap = 35, rect_xmin = 1,
                                                       rect_xmax = 2, 
                                                       annot_text_size = 2.25,
                                                       legend.text.size = 5,
                                                       x.axis.text.size = 5) +
    barchart_scales[[voi[i]]]
  print(iClusterBayes_barcharts[[i]])
  ggsave(filename = paste0("iClusterBayes_", voi[i], "_barchart.png"),
         path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/iClusterBayes_extra"), 
         width = 2320, height = 2320, device = 'png', units = "px",
         dpi = 700)
  dev.off()
}
names(iClusterBayes_barcharts) = voi
rm(loc, chifit)

# Multiplot (PNG) - bar charts
ggarrange(iClusterBayes_barcharts[[1]], iClusterBayes_barcharts[[2]], iClusterBayes_barcharts[[3]],
          iClusterBayes_barcharts[[4]], iClusterBayes_barcharts[[5]], iClusterBayes_barcharts[[6]],
          iClusterBayes_barcharts[[7]], iClusterBayes_barcharts[[8]], iClusterBayes_barcharts[[9]],
          iClusterBayes_barcharts[[10]], iClusterBayes_barcharts[[11]],
          ncol = 3, nrow = 4, labels = c("A", "B", "C", "D", "E", "F", "G", "H",
                                         "I", "J", "K"),
          font.label = list(size = 8, face = "bold", color ="black"))
ggsave(filename = "Multiplot_iClusterBayes_barcharts.png",
       path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/iClusterBayes_extra"), 
       width = 7000, height = 8000, device = 'png', units = "px",
       dpi = 700)
dev.off()

# Just the significant ones
iClusterBayes_barcharts_sig = list()
plotdata_bar_sig = clust_annot_pheno_nonas %>% dplyr::select(iClusterBayes, Histology, Stage,
                                                       `ER status`, `PR status`, `HER2 status`) %>%
  dplyr::mutate(iClusterBayes = paste0("MOVICS_", iClusterBayes))
plotdata_bar_sig$iClusterBayes = factor(plotdata_bar_sig$iClusterBayes)
voi_sig = setdiff(colnames(plotdata_bar_sig), "iClusterBayes")
for (i in 1:length(voi_sig)) {
  chifit = chisq_outputs[["iClusterBayes"]]
  loc = which(grepl(voi_sig[i], chifit$Comparison))
  chifit = chifit[loc, ]
  iClusterBayes_barcharts_sig[[i]] = create_annot_barchart(plotdata = plotdata_bar_sig, fill = voi_sig[i],
                                                           chifit = chifit,
                                                           na.action = "na.omit",
                                                           algorithm = "iClusterBayes",
                                                           barchart_ylim = 650,
                                                           text_y = 630, rect_ymin = 530,
                                                           rect_ymax = 650, x_annot = 1.5,
                                                           v_gap = 35, rect_xmin = 1,
                                                           rect_xmax = 2, 
                                                           annot_text_size = 2.25,
                                                           legend.text.size = 5,
                                                           x.axis.text.size = 5) +
    barchart_scales[[voi_sig[i]]]
  print(iClusterBayes_barcharts_sig[[i]])
  ggsave(filename = paste0("sig_iClusterBayes_", voi_sig[i], "_barchart.png"),
         path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/iClusterBayes_extra"), 
         width = 2320, height = 2320, device = 'png', units = "px",
         dpi = 700)
  dev.off()
}
names(iClusterBayes_barcharts_sig) = voi_sig
rm(loc, chifit)

# Multiplot (PNG) - bar charts
ggarrange(iClusterBayes_barcharts_sig[[1]], iClusterBayes_barcharts_sig[[2]], iClusterBayes_barcharts_sig[[3]],
          iClusterBayes_barcharts_sig[[4]], iClusterBayes_barcharts_sig[[5]], 
          ncol = 2, nrow = 3, labels = c("A", "B", "C", "D", "E"),
          font.label = list(size = 8, face = "bold", color ="black"))
ggsave(filename = "sig_Multiplot_iClusterBayes_barcharts.png",
       path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/iClusterBayes_extra"), 
       width = 5500, height = 7000, device = 'png', units = "px",
       dpi = 700)
dev.off()

# Sunburst plot ###
Pheno_sunburst_iClusterBayes = clust_annot_pheno
Pheno_sunburst_iClusterBayes$`ER status` = gsub("Unknown", "Unkn ER status", Pheno_sunburst_iClusterBayes$`ER status`)
Pheno_sunburst_iClusterBayes$`ER status` = gsub("Positive", "ER+", Pheno_sunburst_iClusterBayes$`ER status`)
Pheno_sunburst_iClusterBayes$`ER status` = gsub("Negative", "ER-", Pheno_sunburst_iClusterBayes$`ER status`)
Pheno_sunburst_iClusterBayes$`HER2 status` = gsub("Unknown", "Unkn HER2 status", 
                                                  Pheno_sunburst_iClusterBayes$`HER2 status`)
Pheno_sunburst_iClusterBayes$`HER2 status` = gsub("Positive", "HER2+", Pheno_sunburst_iClusterBayes$`HER2 status`)
Pheno_sunburst_iClusterBayes$`HER2 status` = gsub("Negative", "HER2-", Pheno_sunburst_iClusterBayes$`HER2 status`)
Pheno_sunburst_iClusterBayes = Pheno_sunburst_iClusterBayes %>%
  dplyr::select(iClusterBayes, `ER status`, `HER2 status`) %>%
  group_by(iClusterBayes, `ER status`, `HER2 status`) %>%
  summarise(Counts = n()) %>%
  as.data.frame()

sunburst_coloring_iClusterBayes = data.frame(stringsAsFactors = FALSE,
                                             colors = tolower(gplots::col2hex(c("#2EC4B6", "#E71D36", 
                                                                                "#C11D9C", "#0F1682",  "grey40",
                                                                                "#0B9EF8", "#560DA7", "mistyrose1", 
                                                                                "hotpink4", "grey40"))),
                                             labels = c("iClusterBayes1", "iClusterBayes2",
                                                        "ER-", "ER+", "Unkn ER status",
                                                        "HER2-", "HER2+", "Indeterminate",
                                                        "Equivocal", "Unkn HER2 status"))

sunburstDF_iClusterBayes = as.sunburstDF(Pheno_sunburst_iClusterBayes, value_column = "Counts", add_root = FALSE) %>%
  inner_join(sunburst_coloring_iClusterBayes, by = "labels")

pie_iClusterBayes = plot_ly() %>%
  add_trace(ids = sunburstDF_iClusterBayes$ids, labels= sunburstDF_iClusterBayes$labels, 
            parents = sunburstDF_iClusterBayes$parents, 
            values= sunburstDF_iClusterBayes$values, type='sunburst', branchvalues = 'total',
            insidetextorientation='radial', maxdepth = 5,
            marker = list(colors = sunburstDF_iClusterBayes$colors)) %>%
  layout(
    grid = list(columns =1, rows = 1),
    margin = list(l = 0, r = 0, b = 0, t = 0)
  )
pie_iClusterBayes
rm(Pheno_sunburst_iClusterBayes, sunburstDF_iClusterBayes, sunburst_coloring_iClusterBayes, pie_iClusterBayes); gc()

# Compare MOVICS iClusterBayes to MOVICS consensus
CS_comp_list[["iClusterBayes"]]$table = table(clust_annot_pheno2$iClusterBayes, 
                                              paste0("CS", clust_annot_pheno2$clust))

CS_comp_list[["iClusterBayes"]]$ARI = calculate_ari_index(cluster_df1 = iClusterBayes_clust_res %>%
                                                            dplyr::rename(Cluster = iClusterBayes) %>%
                                                            mutate(Cluster = gsub("MOVICS_iClusterBayes", "", Cluster)),
                                                          cluster_df2 = as.data.frame(consensus$clust.res) %>%
                                                            dplyr::rename(Cluster = clust),
                                                          sample_col = "samID",
                                                          clust_col = "Cluster",
                                                          suffixes = c("_MOVICS_iClusterBayes", "_CS"))

CS_comp_list[["iClusterBayes"]]$NMI = calculate_nmi_index(cluster_df1 = iClusterBayes_clust_res %>%
                                                            dplyr::rename(Cluster = iClusterBayes) %>%
                                                            mutate(Cluster = gsub("MOVICS_iClusterBayes", "", Cluster)),
                                                          cluster_df2 = as.data.frame(consensus$clust.res) %>%
                                                            dplyr::rename(Cluster = clust),
                                                          sample_col = "samID",
                                                          clust_col = "Cluster",
                                                          suffixes = c("_MOVICS_iClusterBayes", "_CS"))

# Print all comparison data
print(CS_comp_list$iClusterBayes)
