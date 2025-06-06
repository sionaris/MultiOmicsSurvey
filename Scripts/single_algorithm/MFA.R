# Import data from gitignored "Resources/BRCA complete/" folder #####

# This RDS object was produced using the Scripts/MOVICS/MOVICS_baseline.R script
input = readRDS("Resources/TCGA/mm_input.rds")

# Setup environment variables for markdown #####

# Ensure reproducibility
RNGversion("4.2.2")
set.seed(123)

# Load custom helper functions
source("Scripts/automated_scripts/custom_functions.R")
source("Scripts/automated_scripts/modified_MOVICS_functions.R")

# Preamble
home = getwd()
algorithm = "MFA"
alg_feature_pref = "cols" # Where does the algorithm expect the features to be
citation = fetch_citation(algorithm = algorithm)
data_source = "TCGA" # e.g. TCGA, TCGA-transNEO, transNEO-PARTNER
data_types = "RNAseq-CNV-Methylation-miRNA-SNPs" # e.g. RNAseq, RNAseq-CNV-miRNA
evaluation_source = "transNEO" # e.g. PARTNER, transNEO-PARTNER 
title = paste0("Results from ", algorithm)
subtitle = paste0("<b>Train</b>: ", data_source, " ", data_types, 
                  " | <b>Evaluation</b>: ", evaluation_source)
in_a_nutshell = fetch_in_a_nutshell(algorithm = algorithm)
ground_truth_labels = openxlsx::read.xlsx("Results/MOVICS_baseline/MOVICS_TCGA_RNAseq-CNV-Methylation-miRNA-SNPs_eval_on_transNEO_clusterings.xlsx")
ground_truth_k = 2 # optk from MOVICS
optk_boolean = "FALSE" # either TRUE or FALSE. Answers whether the algorithm suggests an optimal k
optk_text = ifelse(optk_boolean == TRUE,
                   "<u>suggests</u> an estimate of the optimal number of multi-omic clusters $k$",
                   "<u>does not suggest</u> an optimal number of multi-omic clusters $k$")

# Detailed description of the algorithm
description = paste(readLines(paste0("Resources/algorithm_descriptions/", algorithm,
                                     "_description.Rmd")),
                    collapse = "\n") # File path to .Rmd file within Resources/algorithm_descriptions

# Create algorithm directory if it doesn't exist
if (!dir.exists(paste0(home, "/Results/single_algorithm/", 
                       algorithm))) {
  dir.create(paste0(home, "/Results/single_algorithm/", 
                    algorithm))
}

# Preprocessing flags and code #####
library(stringr)
library(dplyr)

# All preprocessing for this input has already been performed using the 
# Scripts/MOVICS/MOVICS_baseline.R script

# Extract the names of the modalities that will be used
modalities = unlist(strsplit(data_types, "-"))

# Replace with TRUE wherever features are in rows
features_in_rows = rep(TRUE, length(modalities))

# Run algorithm #####
library(FactoMineR)

# Check what kind of arrangement the algorithm requires as input 
# (i.e. features in rows or columns)?
if (alg_feature_pref == "rows") {
  rogue_indices = which(features_in_rows == FALSE)
  if (length(rogue_indices >= 1)) {
    for (index in rogue_indices) {
      cols = colnames(input[[index]])
      rows = rownames(input[[index]])
      input[[index]] = t(input[[index]])
      rownames(input[[index]]) = cols # transpose names
      colnames(input[[index]]) = rows # transpose names
      rm(rows, cols)
    }
  }
} else if (alg_feature_pref == "cols") {
  rogue_indices = which(features_in_rows == TRUE)
  if (length(rogue_indices >= 1)) {
    for (index in rogue_indices) {
      cols = colnames(input[[index]])
      rows = rownames(input[[index]])
      input[[index]] = t(input[[index]])
      rownames(input[[index]]) = cols # transpose names
      colnames(input[[index]]) = rows # transpose names
      rm(rows, cols)
    }
  }
}
rm(rogue_indices, index); gc()

# Import clinical data for the TCGA samples of interest
clinical_data = openxlsx::read.xlsx("Resources/TCGA/clinical_data.xlsx")

# Keep the top 33% features for each continuous dataset based on MAD
mfa_input = input

# for continuous datasets
mfa_input[c("RNAseq", "CNV", "Methylation", "miRNA")] = lapply(mfa_input[c("RNAseq", "CNV", "Methylation", "miRNA")], function(x) {
  colMAD <- apply(x, 2, mad, na.rm = TRUE)
  ordered_cols <- order(colMAD, decreasing = TRUE)
  n_top <- ceiling(0.33 * ncol(x))
  x[, ordered_cols[1:n_top]]
})

# Keep cancer drivers and the top 33%  most mutated genes of the rest for SNPs
COSMIC_BC_drivers = read.csv("Resources/COSMIC_CGC_Breast_somatic.csv")$Gene.Symbol %>%
  as.character()

snp_mat <- input[["SNPs"]]
mutation_counts <- colSums(snp_mat, na.rm = TRUE)
non_drivers <- setdiff(colnames(snp_mat), COSMIC_BC_drivers)
ordered_non_drivers <- non_drivers[order(mutation_counts[non_drivers], decreasing = TRUE)]
n_top <- ceiling(0.33 * length(non_drivers))
top_non_drivers <- ordered_non_drivers[1:n_top]
genes_to_keep <- unique(c(COSMIC_BC_drivers, top_non_drivers))
mfa_input[["SNPs"]] <- snp_mat[, intersect(colnames(snp_mat), genes_to_keep)]

# Export input for HPC
for (i in 1:length(mfa_input)) {
  colnames(mfa_input[[i]]) = paste0(names(mfa_input)[i], "_", colnames(mfa_input[[i]]))
}
snps_no = as.numeric(length(intersect(colnames(snp_mat), genes_to_keep)))

mfa_input = do.call(cbind, mfa_input)
snp_cols = sum(grepl("SNPs_", colnames(mfa_input)))
rna_cols = sum(grepl("RNAseq_", colnames(mfa_input)))
cnv_cols = sum(grepl("CNV_", colnames(mfa_input)))
mirna_cols = sum(grepl("miRNA_", colnames(mfa_input)))
methyl_cols = sum(grepl("Methylation_", colnames(mfa_input)))

mfa_input = as.data.frame(mfa_input)
mfa_input[, 1:snps_no] <- lapply(mfa_input[, 1:snps_no], as.factor)

saveRDS(mfa_input, "Resources/MFA_input.rds")
rm(i, snps_no, snp_mat, mutation_counts, non_drivers, ordered_non_drivers,
   n_top, top_non_drivers, genes_to_keep); gc()

# HPC: We run the job with scripts in the Scripts/single_algorithm/MFA_HPC/ directory
# Hyperparameter setup
# The only hyperparameter here is the number of factors
NCOMP = 200

# The results are in the Resources/HPC output/MFA_HPC/ directory (gitignored, 5GB)
library(factoextra)
mfa200 = readRDS(paste0("Resources/HPC output/MFA_HPC/MFA_ncp_", NCOMP, "_results.rds"))

# Create subdirectory for supplementary plots
if (!dir.exists(paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement"))) {
  dir.create(paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement"))
}

# Inspect eigenvalues and variance explained
eig.val = get_eigenvalue(mfa200$MFA)
scree = fviz_screeplot(mfa200$MFA, addlabels = TRUE, ylim = c(0, 7.5),
               ggtheme = theme_classic(), ncp = NCOMP) +
  theme(plot.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 0.1))

for (i in seq_along(scree$layers)) {
  if (inherits(scree$layers[[i]]$geom, "GeomText")) {
    scree$layers[[i]]$aes_params$size <- 0.1
    scree$layers[[i]]$aes_params$vjust <- -2
    scree$layers[[i]]$aes_params$hjust <- 0.5
  }
}

ggsave(plot = scree,
  filename = paste0("screeplot_MCA_ann_df_ncp", NCOMP, ".png"),
  path = paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement"),
  dpi = 700, height = 3400, 
  width = 5400, units = "px", device = "png"
)
ggsave(plot = scree,
       filename = paste0("screeplot_MCA_ann_df_ncp", NCOMP, ".pdf"),
       path = paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement"),
       dpi = 700, height = 3400, 
       width = 5400, units = "px", device = "pdf",
)
gc()

# Just the top 20 ncp:
top_ncp = 20
scree_top_ncp = fviz_screeplot(mfa200$MFA, addlabels = TRUE, ylim = c(0, 7.5),
                       ggtheme = theme_classic(), ncp = top_ncp) +
  theme(plot.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"))

for (i in seq_along(scree_top_ncp$layers)) {
  if (inherits(scree_top_ncp$layers[[i]]$geom, "GeomText")) {
    scree_top_ncp$layers[[i]]$aes_params$size <- 2.5
    scree_top_ncp$layers[[i]]$aes_params$vjust <- -2
    scree_top_ncp$layers[[i]]$aes_params$hjust <- 0.5
  }
}

ggsave(plot = scree_top_ncp,
       filename = paste0("screeplot_MCA_ann_df_ncp", top_ncp, ".png"),
       path = paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement"),
       dpi = 700, height = 3400, 
       width = 5400, units = "px", device = "png"
)
ggsave(plot = scree_top_ncp,
       filename = paste0("screeplot_MCA_ann_df_ncp", top_ncp, ".pdf"),
       path = paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement"),
       dpi = 700, height = 3400, 
       width = 5400, units = "px", device = "pdf",
)
gc()

# Based on the visual inspection of the  plot, there is a knee/elbow after 
# dimension 11. So we therefore pick 11 dimensions.
# Additionally, after dimension 11, we see that the %var.explained by each
# new dimension is less than 1%

# inflection::uik() identifies dimension 37 as the knee point, but by then
# %var.explained is 0.37% and in terms of cumulative variance explained, though,
# dim 11 is at ~27% while dim 37 is at ~42%
visual_optdim = 11

# Additionally, inflection identifies dim 134 as the knee/elbow point for 
# cumulative variance explained

# We will use all three options with M3C clustering and determine on the best result
# by examining stability statistics along with parsimony
library(inflection)
uik(1:nrow(eig.val), eig.val[, 1]) # eigenvalues: 37
uik(1:nrow(eig.val), eig.val[, 2]) # %var.explained: 37
uik(1:nrow(eig.val), eig.val[, 3]) # %cumul.var.explained: 134

uiks = list(eigenvalue = uik(1:nrow(eig.val), eig.val[, 1]), # eigenvalues: 37)
            perc.var.expl. = uik(1:nrow(eig.val), eig.val[, 2]), # %var.explained: 37
            perc.cumul.var.expl. = uik(1:nrow(eig.val), eig.val[, 3]) # %cumul.var.explained: 134
)
check_dims = sort(unlist(c(visual_optdim, c(unique(uiks)))))

# M3C k-means clustering ###
library(M3C)

# Import resources
scheme = readRDS("Resources/scheme.rds")
annCol = scheme$annCol
annColors = scheme$annColors
cluster_colors = scheme$clust.colors
col.list = scheme$col.list
var2comp = scheme$var2comp
rm(scheme); gc()

# Here we create a class column for ER status
m3c_des = annCol
m3c_des$class = m3c_des$`ER status`
m3c_des$ID = rownames(m3c_des)

# Run a loop of M3C for all results. In the end pick the one that achieves the
# best combo of entropy, RCSI, p-value

M3C_clusterings = list()
for (n.dim in check_dims) {
  m3c_input = t(mfa200[["MFA"]][["ind"]][["coord"]][, 1:n.dim]) %>% as.data.frame()
  rownames(m3c_input) = paste0("MFAdim", 1:nrow(m3c_input))
  # colnames(m3c_input) = colnames(input[["SNPs"]])
  
  RNGversion("4.2.2")
  M3C_clusterings[[paste0("n.dim = ", n.dim)]] = M3C(m3c_input, des = m3c_des,
                                                       iters = 100, repsref = 250, 
                                                       repsreal = 250, seed = 123, fsize = 18, lthick = 2, dotsize = 1.25,
                                                       clusteralg = "km", maxK = 10)
  cat("Done with n.dim =", n.dim, ".", "\n")
}

# console logs:
# n.dim = 11, optimal K: 3
# n.dim = 37, optimal K: 8
# n.dim = 134, optimal K: 8

# Inspection of clustering results #####
scores_df = as.data.frame(cbind(list(n.dim = unlist(lapply(check_dims, 
                                                           function (x) rep(paste0("n.dim = ", x), 9)))), 
                                rbind(M3C_clusterings[[paste0("n.dim = ", check_dims[1])]][["scores"]],
                                      M3C_clusterings[[paste0("n.dim = ", check_dims[2])]][["scores"]],
                                      M3C_clusterings[[paste0("n.dim = ", check_dims[3])]][["scores"]])))

# RCSI plot
rcsi = ggplot(scores_df, aes(x = K, y = RCSI, group = n.dim, color = factor(n.dim)))+
  geom_line(linewidth = 0.3*2)+
  # geom_errorbar(aes(ymin = RCSI - RCSI_SE,
  #                   ymax = RCSI + RCSI_SE,
  #                   color = "grey"), width = 0.05, linewidth = 0.1*2)+
  geom_point(size = 1)+
  scale_x_continuous(limits = c(1.9, 10.1), breaks = seq(2, 10, 1))+
  scale_y_continuous(limits = c(min(scores_df$RCSI - scores_df$RCSI_SE) - 0.15, 
                                max(scores_df$RCSI + scores_df$RCSI_SE) + 0.15), 
                     breaks = c(-rev(seq(0, abs(round(min(scores_df$RCSI - scores_df$RCSI_SE), 1)), 0.5)), 
                                seq(0, round(max(scores_df$RCSI + scores_df$RCSI_SE), 1), 0.5)))+
  scale_colour_manual(values=rcartocolor::carto_pal(n = 9, "Safe"), name="n.dim") +
  theme(plot.title = element_text(size = 5*2, face = "bold"),
        axis.title.x = element_text(size = 4*2, face = "bold"),
        axis.title.y = element_text(size = 4*2, face = "bold"),
        axis.ticks = element_line(linewidth = 0.15*2),
        axis.text.x = element_text(size = 4*2),
        axis.text.y = element_text(size = 4*2),
        legend.position = "right",
        legend.text = element_text(size = 5),
        legend.title = element_text(face = "bold", size = 6.5, hjust = 0.5),
        legend.key.spacing = unit(2, "mm"),
        panel.background = element_rect(fill = "white", 
                                        colour = "white"),
        panel.grid = element_blank(),
        axis.line = element_line(linewidth = 0.2*2))+
  labs(title = "RCSI vs. number of clusters K for different n.dim values",
       x = "K", y = "RCSI")
rcsi
ggsave(filename = "RCSI.pdf",
       path = "Results/single_algorithm/MFA", 
       width = 140, height = 100, device = 'pdf', units = "mm",
       dpi = 350)
dev.off()

# Statistical significance of clusters
library(ggnewscale)
statsig_clust = ggplot(scores_df, aes(x = K, y = P_SCORE, color = factor(n.dim)))+
  geom_point(size = 1.5*2, alpha = 0.6)+
  scale_color_manual(values=c(rcartocolor::carto_pal(n = 9, "Safe")),
                     name="n.dim") +
  new_scale("color") +
  geom_hline(linetype = "dashed", linewidth = 0.2*2,
             aes(yintercept = -log10(0.05), color = "grey40"))+
  scale_color_manual(values="grey40", labels = expression(-log[10]("0.05")),
                     name="Statistical \nsignificance") +
  scale_x_continuous(limits = c(1.9, 10.1), breaks = seq(2, 10, 1))+
  scale_y_continuous(limits = c(min(scores_df$P_SCORE) - 0.15, 
                                max(scores_df$P_SCORE) + 0.15), 
                     breaks = seq(0, max(scores_df$P_SCORE) + 0.1, 1))+
  theme_bw()+
  theme(plot.title = element_text(size = 5*2, face = "bold"),
        axis.title.x = element_text(size = 4*2, face = "bold"),
        axis.title.y = element_text(size = 4*2, face = "bold"),
        axis.ticks = element_line(linewidth = 0.15*2),
        axis.text.x = element_text(size = 4*2),
        axis.text.y = element_text(size = 4*2),
        legend.position = "right",
        legend.text = element_text(size = 5),
        legend.title = element_text(face = "bold", size = 6.5, hjust = 0.5),
        legend.key.spacing = unit(2, "mm"),
        panel.background = element_rect(fill = "white", 
                                        colour = "white"),
        panel.grid = element_blank(),
        axis.line = element_line(linewidth = 0.2*2))+
  labs(title = "Statistical significance of different values of K",
       y = bquote(bold(-log[10]("p"))))
statsig_clust
ggsave(filename = "Stat_Sig.pdf",
       path = "Results/single_algorithm/MFA", 
       width = 140, height = 100, device = 'pdf', units = "mm",
       dpi = 350)
dev.off()

# Entropy
entropy = ggplot(scores_df, aes(x = K, y = ENTROPY_REAL, alpha = 0.85, 
                                group = n.dim, color = factor(n.dim)))+
  geom_line(linewidth = 0.6)+
  geom_point(size = 1) +
  scale_x_continuous(limits = c(1.9, 10.1), breaks = seq(2, 10, 1))+
  scale_y_continuous(limits = c(0, 100000), breaks = seq(0, 100000, 10000),
                     labels = scales::label_comma()) +
  scale_colour_manual(values=rcartocolor::carto_pal(n = 9, "Safe"), name="n.dim") +
  theme_bw()+
  theme(panel.border = element_rect(linewidth = 0.2*2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 5*2, face = "bold", vjust = 0.5, hjust = 0.5),
        legend.title = element_text(face = "bold", size = 4*2, hjust = 0.5),
        legend.text = element_text(size = 3*2),
        legend.key.size = unit(0.2*2, "cm"),
        legend.margin = ggplot2::margin(0, 0, 0, 0, unit = "mm"),
        legend.spacing.y = unit(0.5*2, units = "mm"),
        axis.title.x = element_text(size = 4*2, face = "bold"),
        axis.title.y = element_text(size = 4*2, face = "bold"),
        axis.ticks = element_line(linewidth = 0.15*2),
        axis.text.x = element_text(size = 4*2),
        axis.text.y = element_text(size = 4*2))+
  labs(y = "Entropy",
       x = "K",
       title = "Entropy in Real Data")+
  guides(alpha = "none")
entropy
ggsave(filename = "Entropy.pdf",
       path = "Results/single_algorithm/MFA", 
       width = 140, height = 100, device = 'pdf', units = "mm",
       dpi = 350)
dev.off()

# Determine the best clustering based on p-value filtering and then examining the RCSI
# Entropy is biased for higher K
best_clusterings = scores_df[scores_df$NORM_P < 0.05, ] %>%
  dplyr::arrange(desc(RCSI))

# According to these criteria the best clustering is:
print(best_clusterings[1, ]) # n.dim = 11, K = 3, lowest entropy, p-value, highest RCSI

# Determine optn.dim
optn.dim = as.numeric(str_sub(best_clusterings$n.dim[1], -2, -1))
conclusion = paste0("Best clustering based on statistical significance and RCSI is for $n.dim = ", 
                     optn.dim, "$. We therefore proceed with $n.dim = ",
                     optn.dim, "$.")

optk = as.numeric(best_clusterings[1, ]$K)

# Exploratory plot of the results (just the first two dimensions)
assignments = as.data.frame(M3C_clusterings[[best_clusterings$n.dim[1]]]$realdataresults[[optk]]$assignments) %>%
  tibble::rownames_to_column(var = "Sample.ID")
colnames(assignments)[2] = "Cluster"
cluster_DF = as.data.frame(mfa200[["MFA"]][["ind"]][["coord"]][, 1:optn.dim]) %>%
  tibble::rownames_to_column(var = "Sample.ID") %>%
  inner_join(assignments, by = "Sample.ID")
colnames(cluster_DF)[grep("Dim", colnames(cluster_DF))] = paste0("MFA_", 
                                                                 colnames(cluster_DF)[grep("Dim", colnames(cluster_DF))])
cluster_DF$Cluster = paste0("MFA", cluster_DF$Cluster)
cluster_scatter = ggplot(data = cluster_DF, aes(x = MFA_Dim.1, y = MFA_Dim.2, group = Cluster,
                                                color = factor(Cluster))) +
  geom_point(size = 0.5) +
  scale_color_manual(values=c(rcartocolor::carto_pal(n = 12, "Safe")[1:3]),
                     name="Cluster") +
  theme_bw()+
  theme(plot.title = element_text(size = 5*2, face = "bold"),
        axis.title.x = element_text(size = 4*2, face = "bold"),
        axis.title.y = element_text(size = 4*2, face = "bold"),
        axis.ticks = element_line(linewidth = 0.15*2),
        axis.text.x = element_text(size = 4*2),
        axis.text.y = element_text(size = 4*2),
        legend.position = "right",
        legend.text = element_text(size = 8),
        legend.title = element_text(face = "bold", size = 6.5, hjust = 0.5),
        legend.key.spacing = unit(2, "mm"),
        panel.background = element_rect(fill = "white", 
                                        colour = "white"),
        panel.grid = element_blank(),
        axis.line = element_line(linewidth = 0.2*2))+
  labs(title = paste0("Cluster scatter plot for the first two dimensions and K = ",
                      best_clusterings$K[1], " (total dims = ",
                      optn.dim, ")"),
       x = "MFA dim1", y = "MFA dim2")
cluster_scatter
ggsave(filename = "Cluster_scatter_2D.pdf",
       path = "Results/single_algorithm/MFA", 
       width = 140, height = 100, device = 'pdf', units = "mm",
       dpi = 350)
dev.off()

# Explore MFA results / Interpretability #####
# Variables in the 11 dimensions ###

varplots = list()
for (dim1 in 1:optn.dim) {
  for (dim2 in 1:optn.dim) {
    varplots[[paste0("Dim", dim1, "_Dim", dim2)]] = fviz_mfa_var(mfa200$MFA, "group",
                                                                 axes = c(dim1, dim2),
                                                                 repel = TRUE,
                                                                 col.var = "black",
                                                                 title = paste0("MFA dims ", dim1, " vs. ", dim2))
    varplots[[paste0("Dim", dim1, "_Dim", dim2)]] = varplots[[paste0("Dim", dim1, "_Dim", dim2)]] +
      theme_classic() +
      theme(plot.title = element_text(size = 5*2, face = "bold"),
            axis.title.x = element_text(size = 4*2, face = "bold"),
            axis.title.y = element_text(size = 4*2, face = "bold"),
            axis.ticks = element_line(linewidth = 0.15*2),
            axis.text.x = element_text(size = 4*2),
            axis.text.y = element_text(size = 4*2),
            axis.line = element_line(linewidth = 0.2*2))
  }
}

# Arrange in a 11x11 grid
library(ggpubr)
FIG_varplot = ggarrange(plotlist = varplots,
                        ncol = optn.dim, nrow = optn.dim,
                        common.legend = TRUE,
                        align = "hv", labels = NULL)

ggsave(plot = FIG_varplot,
       filename = paste0("MFA_varplots_", optn.dim, "_dims.png"),
       path = paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement"),
       dpi = 400, height = 10400, 
       width = 10400, units = "px", device = "png"
)

rm(dim1, dim2); gc()

# Modality contribution bar charts ###
barplots = list()
for (dim in 1:optn.dim) {
    barplots[[paste0("Dim", dim)]] = fviz_contrib(mfa200$MFA, "group", axes = dim)
    barplots[[paste0("Dim", dim)]] = barplots[[paste0("Dim", dim)]] +
      theme_classic() +
      theme(plot.title = element_text(size = 5*2, face = "bold"),
            axis.title.x = element_text(size = 4*2, face = "bold"),
            axis.title.y = element_text(size = 4*2, face = "bold"),
            axis.ticks = element_line(linewidth = 0.15*2),
            axis.text.x = element_text(size = 4*2),
            axis.text.y = element_text(size = 4*2),
            axis.line = element_line(linewidth = 0.2*2))
}

FIG_barplot = ggarrange(plotlist = barplots,
                        ncol = 3, nrow = 4,
                        common.legend = TRUE,
                        align = "hv", labels = NULL)

ggsave(plot = FIG_barplot,
       filename = paste0("MFA_contrib_barplots_", optn.dim, "_dims.png"),
       path = paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement"),
       dpi = 700, height = 10400, 
       width = 10400, units = "px", device = "png"
)

rm(dim); gc()

# Feature contribution plots ###
featplots_10000 = list()
for (dim in 1:optn.dim) {
  featplots_10000[[paste0("Dim", dim)]] = fviz_contrib(mfa200$MFA, "quanti.var", axes = dim, top = 10000,
                                                 palette = "jco")
  featplots_10000[[paste0("Dim", dim)]] = featplots_10000[[paste0("Dim", dim)]] +
    theme_classic() +
    theme(plot.title = element_text(size = 5*2, face = "bold"),
          axis.title.x = element_text(size = 4*2, face = "bold"),
          axis.title.y = element_text(size = 4*2, face = "bold"),
          axis.ticks = element_line(linewidth = 0.15*2),
          axis.text.x = element_text(size = 1, angle = 45),
          axis.text.y = element_text(size = 4*2),
          axis.line = element_line(linewidth = 0.2*2))
}

FIG_featplot_top10000 = ggarrange(plotlist = featplots_10000,
                        ncol = 3, nrow = 4,
                        common.legend = TRUE,
                        align = "hv", labels = NULL)

ggsave(plot = FIG_featplot_top10000,
       filename = paste0("MFA_contrib_featplots_", optn.dim, "_dims_top10000.png"),
       path = paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement"),
       dpi = 700, height = 10400, 
       width = 13400, units = "px", device = "png"
)

# top 100
featplots_100 = list()
for (dim in 1:optn.dim) {
  featplots_100[[paste0("Dim", dim)]] = fviz_contrib(mfa200$MFA, "quanti.var", axes = dim, top = 100,
                                                     palette = "jco")
  featplots_100[[paste0("Dim", dim)]] = featplots_100[[paste0("Dim", dim)]] +
    theme_classic() +
    theme(plot.title = element_text(size = 5*2, face = "bold"),
          axis.title.x = element_text(size = 4*2, face = "bold"),
          axis.title.y = element_text(size = 4*2, face = "bold"),
          axis.ticks = element_line(linewidth = 0.15*2),
          axis.text.x = element_text(size = 1, angle = 45),
          axis.text.y = element_text(size = 4*2),
          axis.line = element_line(linewidth = 0.2*2))
}

FIG_featplot_top100 = ggarrange(plotlist = featplots_100,
                                ncol = 3, nrow = 4,
                                common.legend = TRUE,
                                align = "hv", labels = NULL)

ggsave(plot = FIG_featplot_top100,
       filename = paste0("MFA_contrib_featplots_", optn.dim, "_dims_top100.png"),
       path = paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement"),
       dpi = 700, height = 10400, 
       width = 13400, units = "px", device = "png"
)

rm(dim); gc()

# Plots for ER, HER2 and Stage ###
ER_palette = c(Negative = "#C11D9C", 
               Positive = "#0F1682", 
               Unknown = "grey40")

HER2_palette = c(Negative = "#0B9EF8", 
                 Positive = "#560DA7", 
                 Indeterminate = "mistyrose1", 
                 Equivocal = "hotpink4", 
                 Unknown = "grey40")

stage_palette = c(`Stage I` = "#7cc6ad", 
                  `Stage II` = "#0a6da5", 
                  `Stage III` = "#9a9afc", 
                  `Stage IV` = "#5d032d", 
                  `Unknown` = "grey40")

# ER status
ER_plots = list()
for (dim1 in 1:optn.dim) {
  for (dim2 in 1:optn.dim) {
    
    if (dim1 != dim2) {
      ER_plots[[paste0("Dim", dim1, "_Dim", dim2)]] = fviz_mfa_ind(mfa200$MFA, label = "none",
                                                                   habillage = "ER_status",
                                                                   palette = ER_palette,
                                                                   axes = c(dim1, dim2),
                                                                   addEllipses = TRUE, ellipse.type = "confidence", 
                                                                   repel = TRUE # Avoid text overlapping
      ) 
      ER_plots[[paste0("Dim", dim1, "_Dim", dim2)]] = ER_plots[[paste0("Dim", dim1, "_Dim", dim2)]] +
        theme_classic() +
        ggtitle(paste0("ER status: dims ", dim1, "-", dim2)) +
        theme(plot.title = element_text(size = 5*2, face = "bold"),
              axis.title.x = element_text(size = 4*2, face = "bold"),
              axis.title.y = element_text(size = 4*2, face = "bold"),
              axis.ticks = element_line(linewidth = 0.15*2),
              axis.text.x = element_text(size = 4*2),
              axis.text.y = element_text(size = 4*2),
              axis.line = element_line(linewidth = 0.2*2))
    } else {
      scores <- get_mfa_ind(mfa200$MFA)$coord[, dim1]
      
      dens_df <- data.frame(
        Dim       = scores,
        ER_status = mfa200$MFA$call$X$ER_status
      )
      
      ER_plots[[paste0("Dim", dim1, "_Dim", dim1)]] <-
        ggplot(dens_df, aes(x = Dim, fill = ER_status)) +
        geom_density(alpha = 0.4, adjust = 1, color = "grey50") +
        scale_fill_manual(values = ER_palette) +
        theme_classic() +
        ggtitle(sprintf("Dim %d ER densities", dim1)) +
        xlab(sprintf("Dimension %d coordinate", dim1)) +
        ylab("Density") +
        theme(plot.title = element_text(size = 5*2, face = "bold"),
              axis.title.x = element_text(size = 4*2, face = "bold"),
              axis.title.y = element_text(size = 4*2, face = "bold"),
              axis.ticks = element_line(linewidth = 0.15*2),
              axis.text.x = element_text(size = 4*2),
              axis.text.y = element_text(size = 4*2))
    }
  }
}

FIG_ER_plot = ggarrange(plotlist = ER_plots,
                        ncol = optn.dim, nrow = optn.dim,
                        common.legend = TRUE, legend = "bottom",
                        align = "hv", labels = NULL)

ggsave(plot = FIG_ER_plot,
       filename = paste0("MFA_ER_plots_", optn.dim, "_dims.png"),
       path = paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement"),
       dpi = 400, height = 10800, 
       width = 10400, units = "px", device = "png"
)

rm(dim1, dim2); gc()

# HER2 status
mfa200$MFA$call$X$HER2_status = annCol[rownames(mfa200$MFA$call$X), "HER2 status"]
HER2_plots = list()
for (dim1 in 1:optn.dim) {
  for (dim2 in 1:optn.dim) {
    
    if (dim1 != dim2) {
      HER2_plots[[paste0("Dim", dim1, "_Dim", dim2)]] = fviz_mfa_ind(mfa200$MFA, label = "none",
                                                                     habillage = "HER2_status",
                                                                     palette = HER2_palette,
                                                                     axes = c(dim1, dim2),
                                                                     addEllipses = TRUE, ellipse.type = "confidence", 
                                                                     repel = TRUE # Avoid text overlapping
      ) 
      HER2_plots[[paste0("Dim", dim1, "_Dim", dim2)]] = HER2_plots[[paste0("Dim", dim1, "_Dim", dim2)]] +
        theme_classic() +
        ggtitle(paste0("HER2 status: dims ", dim1, "-", dim2)) +
        theme(plot.title = element_text(size = 5*2, face = "bold"),
              axis.title.x = element_text(size = 4*2, face = "bold"),
              axis.title.y = element_text(size = 4*2, face = "bold"),
              axis.ticks = element_line(linewidth = 0.15*2),
              axis.text.x = element_text(size = 4*2),
              axis.text.y = element_text(size = 4*2),
              axis.line = element_line(linewidth = 0.2*2))
    } else {
      scores <- get_mfa_ind(mfa200$MFA)$coord[, dim1]
      
      dens_df <- data.frame(
        Dim       = scores,
        HER2_status = mfa200$MFA$call$X$HER2_status
      )
      
      HER2_plots[[paste0("Dim", dim1, "_Dim", dim1)]] <-
        ggplot(dens_df, aes(x = Dim, fill = HER2_status)) +
        geom_density(alpha = 0.4, adjust = 1, color = "grey50") +
        scale_fill_manual(values = HER2_palette) +
        theme_classic() +
        ggtitle(sprintf("Dim %d HER2 densities", dim1)) +
        xlab(sprintf("Dimension %d coordinate", dim1)) +
        ylab("Density") +
        theme(plot.title = element_text(size = 5*2, face = "bold"),
              axis.title.x = element_text(size = 4*2, face = "bold"),
              axis.title.y = element_text(size = 4*2, face = "bold"),
              axis.ticks = element_line(linewidth = 0.15*2),
              axis.text.x = element_text(size = 4*2),
              axis.text.y = element_text(size = 4*2))
    }
  }
}

FIG_HER2_plot = ggarrange(plotlist = HER2_plots,
                          ncol = optn.dim, nrow = optn.dim,
                          common.legend = TRUE, legend = "bottom",
                          align = "hv", labels = NULL)

ggsave(plot = FIG_HER2_plot,
       filename = paste0("MFA_HER2_plots_", optn.dim, "_dims.png"),
       path = paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement"),
       dpi = 400, height = 10800, 
       width = 10400, units = "px", device = "png"
)

rm(dim1, dim2); gc()

# Stage
mfa200$MFA$call$X$Stage = annCol[rownames(mfa200$MFA$call$X), "Stage"]
Stage_plots = list()
for (dim1 in 1:optn.dim) {
  for (dim2 in 1:optn.dim) {
    
    if (dim1 != dim2) {
      Stage_plots[[paste0("Dim", dim1, "_Dim", dim2)]] = fviz_mfa_ind(mfa200$MFA, label = "none",
                                                                      habillage = "Stage",
                                                                      palette = stage_palette,
                                                                      axes = c(dim1, dim2),
                                                                      addEllipses = TRUE, ellipse.type = "confidence", 
                                                                      repel = TRUE # Avoid text overlapping
      ) 
      Stage_plots[[paste0("Dim", dim1, "_Dim", dim2)]] = Stage_plots[[paste0("Dim", dim1, "_Dim", dim2)]] +
        theme_classic() +
        ggtitle(paste0("Stage: dims ", dim1, "-", dim2)) +
        theme(plot.title = element_text(size = 5*2, face = "bold"),
              axis.title.x = element_text(size = 4*2, face = "bold"),
              axis.title.y = element_text(size = 4*2, face = "bold"),
              axis.ticks = element_line(linewidth = 0.15*2),
              axis.text.x = element_text(size = 4*2),
              axis.text.y = element_text(size = 4*2),
              axis.line = element_line(linewidth = 0.2*2))
    } else {
      scores <- get_mfa_ind(mfa200$MFA)$coord[, dim1]
      
      dens_df <- data.frame(
        Dim       = scores,
        Stage = mfa200$MFA$call$X$Stage
      )
      
      Stage_plots[[paste0("Dim", dim1, "_Dim", dim1)]] <-
        ggplot(dens_df, aes(x = Dim, fill = Stage)) +
        geom_density(alpha = 0.4, adjust = 1, color = "grey50") +
        scale_fill_manual(values = stage_palette) +
        theme_classic() +
        ggtitle(sprintf("Dim %d Stage densities", dim1)) +
        xlab(sprintf("Dimension %d coordinate", dim1)) +
        ylab("Density") +
        theme(plot.title = element_text(size = 5*2, face = "bold"),
              axis.title.x = element_text(size = 4*2, face = "bold"),
              axis.title.y = element_text(size = 4*2, face = "bold"),
              axis.ticks = element_line(linewidth = 0.15*2),
              axis.text.x = element_text(size = 4*2),
              axis.text.y = element_text(size = 4*2))
    }
  }
}

FIG_Stage_plot = ggarrange(plotlist = Stage_plots,
                           ncol = optn.dim, nrow = optn.dim,
                           common.legend = TRUE, legend = "bottom",
                           align = "hv", labels = NULL)

ggsave(plot = FIG_Stage_plot,
       filename = paste0("MFA_Stage_plots_", optn.dim, "_dims.png"),
       path = paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement"),
       dpi = 400, height = 10800, 
       width = 10400, units = "px", device = "png"
)

rm(dim1, dim2); gc()

# Main results #####
# Examine cluster similarity to MOVICS by measuring NMI and ARI indices #####
# (Jaccard may be misleading)
MFA_clusters = cluster_DF %>% dplyr::select(Sample.ID, Cluster)
MFA_clusters$Sample.ID = gsub("\\.", "-", MFA_clusters$Sample.ID)
rownames(MFA_clusters) = MFA_clusters$Sample.ID

# Calculate ARI and NMI
library(mclust)
library(clue)

ARI_to_MOVICS = calculate_ari_index(cluster_df1 = ground_truth_labels,
                                    cluster_df2 = MFA_clusters,
                                    sample_col = "Sample.ID",
                                    clust_col = "Cluster",
                                    suffixes = c("_MOVICS_CS",
                                                 paste0("_", algorithm)))

NMI_to_MOVICS = calculate_nmi_index(cluster_df1 = ground_truth_labels,
                                    cluster_df2 = MFA_clusters,
                                    sample_col = "Sample.ID",
                                    clust_col = "Cluster",
                                    suffixes = c("_MOVICS_CS",
                                                 paste0("_", algorithm)))

# Very low statistics when compared to the MOVICS. Results differ

# Import coloring scheme
scheme = readRDS("Resources/scheme.rds")
annCol = scheme$annCol
annColors = scheme$annColors
cluster_colors = c("#2EC4B6", "#E71D36", "#FF9F1C")
col.list = scheme$col.list
var2comp = scheme$var2comp %>%
  dplyr::select(-`Consensus Subtype`) %>%
  mutate(Sample.ID = rownames(.)) %>%
  inner_join(MFA_clusters, by = "Sample.ID") %>%
  tibble::column_to_rownames(var = "Sample.ID") %>%
  mutate(MFA = paste0(algorithm, Cluster)) %>%
  dplyr::select(MFA, everything()) %>%
  dplyr::select(-Cluster)
rm(scheme); gc()

# Silhouette
library(MOVICS)
library(cluster)
silhouette = silhouette(as.integer(gsub("MFA", "", cluster_DF$Cluster)),
                        dist = Rfast::Dist(cluster_DF[, grep("MFA_Dim", colnames(cluster_DF))],
                                           method = "euclidean"))

getSilhouette_ggplot(sil      = silhouette,
                     fig.path = paste0(home, "/Results/single_algorithm/", algorithm),
                     fig.name = "Silhouette",
                     height   = 5.5,
                     width    = 5.5,
                     axis_label_size = 12,
                     axis_label_font = "bold",
                     text_size = 1.5,
                     title_size = 16,
                     algorithm = algorithm,
                     save_plot = TRUE)
dev.off()

# MOVICS-like analysis #####
library(ComplexHeatmap)

# Heatmap prep
plotdata <- lapply(lapply(input, as.matrix), 
                   function(mat) mat[, colSums(mat != 0) > 0])
plotdata <- lapply(plotdata, t)
plotdata = getStdiz(
  data = plotdata,
  halfwidth = c(NA, 3, 3, 3, 3), # No halfwidth for SNPs
  centerFlag = c(F, F, F, F, F),
  scaleFlag = c(F, F, F, F, F)
)

MFA_clusters$Cluster = gsub(algorithm, "", MFA_clusters$Cluster)
plot_object = list(clust.res = MFA_clusters %>%
                     dplyr::rename(samID = Sample.ID, clust = Cluster))

# Export consensus clustering object
clust = as.data.frame(plot_object$clust.res)
colnames(clust) = c("Sample.ID", "Cluster")
clust$Cluster = paste0(algorithm, clust$Cluster)
openxlsx::write.xlsx(clust, paste0(home, "/Results/single_algorithm/", algorithm, "/", 
                                   algorithm, "_", data_source, "_",
                                   data_types, "_eval_on_", evaluation_source,
                                   "_clusterings.xlsx"))

# Order features
feature_orders = readRDS("Resources/TCGA/mm_feature_orders.rds")
for (i in 1:length(plotdata)) {
  plotdata[[i]] = plotdata[[i]][feature_orders[[names(plotdata)[i]]], , drop = FALSE]
}

getMoHeatmap_single_algorithm2(algorithm_name = algorithm,
                               data          = plotdata,
                               row.title     = names(plotdata),
                               is.binary     = c(T,F,F,F,F), 
                               legend.name   = c("SNPs",
                                                 "Standardized RNAseq norm. counts",
                                                 "Standardized CNV",
                                                 "Standardized miRNA norm. counts",
                                                 "Standardized Methylation M-values"
                               ),
                               cluster_rows = rep(F, length(plotdata)),
                               cluster_cols = rep(F, length(plotdata)),
                               show.col.dend = rep(F, length(plotdata)),
                               show.colnames = FALSE,
                               show.row.dend = rep(F, length(plotdata)),
                               show.rownames = rep(F, length(plotdata)),
                               clust.res     = plot_object$clust.res, # consensusMOIC-like results
                               # clust.dist.row = c("manhattan", rep("euclidean", 4)),
                               # clust.method.row = rep("ward.D", length(plotdata)),
                               annRow        = NULL, # no selected features
                               color         = col.list,
                               annCol        = annCol, # annotation for samples
                               annColors     = annColors, # annotation color
                               width         = 20, # width of each subheatmap
                               height        = 10, # height of each subheatmap
                               fig.path      = paste0(home, "/Results/single_algorithm/", algorithm),
                               fig.name      = paste0("default_", algorithm, "_Comprehensive_heatmap"))
dev.off()
gc()

# Clinical variables ###
# Remove unknown levels for statistical tests
var2comp_nonas = var2comp
for (i in 1:ncol(var2comp)) {
  nas = which(var2comp[, i] == "Unknown")
  var2comp_nonas[nas, i] = NA
  empties = which(var2comp[, i] == "")
  var2comp_nonas[empties, i] = NA
}
rm(nas, empties); gc()

# Statistical comparisons
clin_comp = compClinvar_single_algorithm(algorithm_name = algorithm,
                                         moic.res = plot_object,
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
                                         res.path = paste0(home, "/Results/single_algorithm/", algorithm, "/"),
                                         output_pdf = TRUE,
                                         pdf_level_col_width = c("7em", "10em"),
                                         pdf_count_col_width = "10em",
                                         pdf_pval_col_width = "3em",
                                         pdf_test_col_width = "8em",
                                         pdf_tab_font_size = 9)

clin_ordinal_comp = compClinvar_ordinal_single_algorithm(algorithm_name = algorithm,
                                                         moic.res = plot_object,
                                                         var2comp = var2comp_nonas %>%
                                                           dplyr::select(number_of_lymphnodes_positive_by_ihc,
                                                                         number_of_lymphnodes_positive_by_he,
                                                                         MFA),
                                                         strata = algorithm,
                                                         ordinalVars = c("number_of_lymphnodes_positive_by_ihc",
                                                                         "number_of_lymphnodes_positive_by_he"),
                                                         includeNA = FALSE,
                                                         tab.name = "Summary of ordinal clinical variables",
                                                         res.path = paste0(home, "/Results/single_algorithm/", algorithm, "/"),
                                                         output_pdf = TRUE,
                                                         pdf_template_loc = paste0(home, "/Scripts/automated_scripts/clincomp_template.Rmd"),
                                                         pdf_level_col_width = c("7em", "10em"),
                                                         pdf_count_col_width = "10em",
                                                         pdf_pval_col_width = "3em",
                                                         pdf_test_col_width = "8em",
                                                         pdf_tab_font_size = 9)

# Oncoprint ###
oncoprint <- compMut_single_algorithm(algorithm_name = algorithm,
                                      moic.res  = plot_object,
                                      mut.matrix   = plotdata$SNPs, # binary somatic mutation matrix
                                      doWord       = TRUE, # generate table in .docx format
                                      doPlot       = TRUE, # draw OncoPrint
                                      freq.cutoff  = 0.05, # keep those genes that mutated in at least 5% of samples
                                      p.adj.cutoff = 0.05, # keep those genes with adjusted p value < 0.05 to draw OncoPrint
                                      innerclust   = TRUE, # perform clustering within each subtype
                                      annCol       = annCol, # same annotation for heatmap
                                      annColors    = annColors, # same annotation color for heatmap
                                      width        = 12, 
                                      height       = 6,
                                      fig.name     = paste0(algorithm, "_", data_source, "_",
                                                            data_types, "_eval_on_", evaluation_source,
                                                            "_oncoprint"),
                                      tab.name     = "Independent test between subtype and mutation",
                                      fig.path     = paste0(home, "/Results/single_algorithm/", algorithm),
                                      res.path     = paste0(home, "/Results/single_algorithm/", algorithm))

# Drug sensitivity comparison ###
drug_sensitivity <- compDrugsen_single_algorithm(algorithm_name = algorithm,
                                                 moic.res    = plot_object,
                                                 norm.expr   = plotdata$RNAseq,
                                                 drugs       = c("Cisplatin", "Paclitaxel", "Lapatinib",
                                                                 "Doxorubicin", "5-Fluorouracil",
                                                                 "Sorafenib"), # a vector of names of drug in GDSC
                                                 tissueType  = "breast", # choose specific tissue type to construct ridge regression model
                                                 test.method = "nonparametric", # statistical testing method
                                                 prefix      = "Violin_plot_of_IC50",
                                                 seed = 123,
                                                 fig.path = paste0(home, "/Results/single_algorithm/", algorithm))

# Agreement with other subtypes ###
subtype_agreement <- compAgree_single_algorithm(algorithm_name = algorithm,
                                                moic.res  = plot_object,
                                                subt2comp = annCol[, c("ER status", "PR status",
                                                                       "HER2 status", "Metastasis", "Stage")],
                                                doPlot    = TRUE,
                                                box.width = 0.2,
                                                fig.name  = "Classification_agreement",
                                                fig.path  = paste0(home, "/Results/single_algorithm/", algorithm),
                                                width     = 12)
dev.off()

# DGEA ###
dgea = runDEA_mod(dea.method = "limma", # we use normalized data as input
                  expr = plotdata$RNAseq,
                  moic.res = plot_object,
                  prefix = "dgea_",
                  sort.p = TRUE,
                  overwt = TRUE,
                  verbose = TRUE,
                  res.path = paste0(home, "/Results/single_algorithm/", algorithm),
                  algorithm = algorithm)

# # Identify unique subtype biomarkers
# # 1. Up-regulated markers
dgea.marker.up <- runMarker_single_algorithm(algorithm_name = algorithm,
                                             moic.res = plot_object,
                                             dea.method    = "limma", # name of DEA method
                                             prefix        = "dgea_", # MUST be the same of argument in runDEA()
                                             dat.path      = paste0(home, "/Results/single_algorithm/", algorithm), # path of DEA files
                                             res.path      = paste0(home, "/Results/single_algorithm/", algorithm), # path to save marker files
                                             p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                                             p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                                             dirct         = "up", # direction of dysregulation in expression
                                             n.marker      = 100, # number of biomarkers for each subtype
                                             doplot        = TRUE, # generate diagonal heatmap
                                             norm.expr     = plotdata$RNAseq, # use normalized expression as heatmap input
                                             annCol        = annCol, # sample annotation in heatmap
                                             annColors     = annColors, # colors for sample annotation
                                             show_rownames = TRUE, # show no rownames (biomarker name)
                                             centerFlag = F,
                                             scaleFlag = F,
                                             halfwidth = 3,
                                             fig.name      = "upregulated_biomarkers_heatmap",
                                             fig.path = paste0(home, "/Results/single_algorithm/", algorithm),
                                             width = 14,
                                             height = 12,
                                             fontsize_row = 3,
                                             name = "normalized RNA-seq")
dev.off()

# # 2. Down-regulated markers
dgea.marker.down <- runMarker_single_algorithm(algorithm_name = algorithm,
                                               moic.res = plot_object,
                                               dea.method    = "limma", # name of DEA method
                                               prefix        = "dgea_", # MUST be the same of argument in runDEA()
                                               dat.path      = paste0(home, "/Results/single_algorithm/", algorithm), # path of DEA files
                                               res.path      = paste0(home, "/Results/single_algorithm/", algorithm), # path to save marker files
                                               p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                                               p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                                               dirct         = "down", # direction of dysregulation in expression
                                               n.marker      = 100, # number of biomarkers for each subtype
                                               doplot        = TRUE, # generate diagonal heatmap
                                               norm.expr     = plotdata$RNAseq, # use normalized expression as heatmap input
                                               annCol        = annCol, # sample annotation in heatmap
                                               annColors     = annColors, # colors for sample annotation
                                               show_rownames = TRUE, # show no rownames (biomarker name)
                                               centerFlag = F,
                                               scaleFlag = F,
                                               halfwidth = 3,
                                               fig.name      = "downregulated_biomarkers_heatmap",
                                               fig.path = paste0(home, "/Results/single_algorithm/", algorithm),
                                               width = 14,
                                               height = 12,
                                               fontsize_row = 3,
                                               name = "normalized RNA-seq")
dev.off()

# DMEA ###
dmea = runDEA_mod(dea.method = "limma", # we use normalized data as input
                  expr = plotdata$Methylation,
                  moic.res = plot_object,
                  prefix = "dmea_",
                  sort.p = TRUE,
                  overwt = TRUE,
                  verbose = TRUE,
                  res.path = paste0(home, "/Results/single_algorithm/", algorithm),
                  algorithm = algorithm)

# # Identify unique subtype biomarkers
# # 1. Up-regulated markers
methyl.marker.up <- runMarker_single_algorithm(algorithm_name = algorithm,
                                               moic.res = plot_object,
                                               dea.method    = "limma", # name of DEA method
                                               prefix        = "dmea_", # MUST be the same of argument in runDEA()
                                               dat.path      = paste0(home, "/Results/single_algorithm/", algorithm), # path of DEA files
                                               res.path      = paste0(home, "/Results/single_algorithm/", algorithm), # path to save marker files
                                               p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                                               p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                                               dirct         = "up", # direction of dysregulation in expression
                                               n.marker      = 100, # number of biomarkers for each subtype
                                               doplot        = TRUE, # generate diagonal heatmap
                                               norm.expr     = plotdata$Methylation, # use normalized expression as heatmap input
                                               annCol        = annCol, # sample annotation in heatmap
                                               annColors     = annColors, # colors for sample annotation
                                               show_rownames = TRUE, # show no rownames (biomarker name)
                                               centerFlag = F,
                                               scaleFlag = F,
                                               halfwidth = 3,
                                               fig.name      = "hypermethylated_biomarkers_heatmap",
                                               fig.path = paste0(home, "/Results/single_algorithm/", algorithm),
                                               width = 14,
                                               height = 12,
                                               fontsize_row = 0, # 3 default
                                               name = "normalized Methylation")
dev.off()

# # 2. Down-regulated markers
methyl.marker.down <- runMarker_single_algorithm(algorithm_name = algorithm,
                                                 moic.res = plot_object,
                                                 dea.method    = "limma", # name of DEA method
                                                 prefix        = "dmea_", # MUST be the same of argument in runDEA()
                                                 dat.path      = paste0(home, "/Results/single_algorithm/", algorithm), # path of DEA files
                                                 res.path      = paste0(home, "/Results/single_algorithm/", algorithm), # path to save marker files
                                                 p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                                                 p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                                                 dirct         = "down", # direction of dysregulation in expression
                                                 n.marker      = 100, # number of biomarkers for each subtype
                                                 doplot        = TRUE, # generate diagonal heatmap
                                                 norm.expr     = plotdata$Methylation, # use normalized expression as heatmap input
                                                 annCol        = annCol, # sample annotation in heatmap
                                                 annColors     = annColors, # colors for sample annotation
                                                 show_rownames = TRUE, # show no rownames (biomarker name)
                                                 centerFlag = F,
                                                 scaleFlag = F,
                                                 halfwidth = 3,
                                                 fig.name      = "hypomethylated_biomarkers_heatmap",
                                                 fig.path = paste0(home, "/Results/single_algorithm/", algorithm),
                                                 width = 14,
                                                 height = 12,
                                                 fontsize_row = 0, # 3 default
                                                 name = "normalized Methylation")
dev.off()

# DmiREA ###
dmiRea = runDEA_mod(dea.method = "limma", # we use normalized data as input
                    expr = plotdata$miRNA,
                    moic.res = plot_object,
                    prefix = "dmiRea_",
                    sort.p = TRUE,
                    overwt = TRUE,
                    verbose = TRUE,
                    res.path = paste0(home, "/Results/single_algorithm/", algorithm),
                    algorithm = algorithm)

# # Identify unique subtype biomarkers
# # 1. Up-regulated markers
miRNA.marker.up <- runMarker_single_algorithm(algorithm_name = algorithm,
                                              moic.res = plot_object,
                                              dea.method    = "limma", # name of DEA method
                                              prefix        = "dmiRea_", # MUST be the same of argument in runDEA()
                                              dat.path      = paste0(home, "/Results/single_algorithm/", algorithm), # path of DEA files
                                              res.path      = paste0(home, "/Results/single_algorithm/", algorithm), # path to save marker files
                                              p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                                              p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                                              dirct         = "up", # direction of dysregulation in expression
                                              n.marker      = 100, # number of biomarkers for each subtype
                                              doplot        = TRUE, # generate diagonal heatmap
                                              norm.expr     = plotdata$miRNA, # use normalized expression as heatmap input
                                              annCol        = annCol, # sample annotation in heatmap
                                              annColors     = annColors, # colors for sample annotation
                                              show_rownames = TRUE, # show no rownames (biomarker name)
                                              centerFlag = F,
                                              scaleFlag = F,
                                              halfwidth = 3,
                                              fig.name      = "upregulated_miRNA_biomarkers_heatmap",
                                              fig.path = paste0(home, "/Results/single_algorithm/", algorithm),
                                              width = 14,
                                              height = 12,
                                              fontsize_row = 0, # 3 default
                                              name = "normalized miRNA")
dev.off()

# # 2. Down-regulated markers
miRNA.marker.down <- runMarker_single_algorithm(algorithm_name = algorithm,
                                                moic.res = plot_object,
                                                dea.method    = "limma", # name of DEA method
                                                prefix        = "dmiRea_", # MUST be the same of argument in runDEA()
                                                dat.path      = paste0(home, "/Results/single_algorithm/", algorithm), # path of DEA files
                                                res.path      = paste0(home, "/Results/single_algorithm/", algorithm), # path to save marker files
                                                p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                                                p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                                                dirct         = "down", # direction of dysregulation in expression
                                                n.marker      = 100, # number of biomarkers for each subtype
                                                doplot        = TRUE, # generate diagonal heatmap
                                                norm.expr     = plotdata$miRNA, # use normalized expression as heatmap input
                                                annCol        = annCol, # sample annotation in heatmap
                                                annColors     = annColors, # colors for sample annotation
                                                show_rownames = TRUE, # show no rownames (biomarker name)
                                                centerFlag = F,
                                                scaleFlag = F,
                                                halfwidth = 3,
                                                fig.name      = "downregulated_miRNA_biomarkers_heatmap",
                                                fig.path = paste0(home, "/Results/single_algorithm/", algorithm),
                                                width = 14,
                                                height = 12,
                                                fontsize_row = 0, # 3 default
                                                name = "normalized miRNA")
dev.off()

# GSEA ###
# Load MSigDb file
MSIGDB.FILE <- paste0(home, "/Resources/Pathways/GO-BP_c5.go.bp.v2024.1.Hs.symbols.gmt")

# GSEA up-regulated
RNGversion("4.2.2")
set.seed(123)
gsea.up <- runGSEA_mod_4.4_single_algorithm(algorithm_name = algorithm,
                                            moic.res     = plot_object,
                                            dea.method   = "limma", # name of DEA method
                                            prefix       = "dgea_", # MUST be the same of argument in runDEA()
                                            dat.path      = paste0(home, "/Results/single_algorithm/", algorithm), # path of DEA files
                                            res.path      = paste0(home, "/Results/single_algorithm/", algorithm), # path to save marker files
                                            msigdb.path  = MSIGDB.FILE, # MUST be the ABSOLUTE path of msigdb file
                                            norm.expr    = plotdata$RNAseq, # use normalized expression to calculate enrichment score
                                            dirct        = "up", # direction of dysregulation in pathway
                                            n.path       = 20,
                                            p.cutoff     = 0.05, # p cutoff to identify significant pathways
                                            p.adj.cutoff = 0.05, # padj cutoff to identify significant pathways
                                            gsva.method  = "gsva", # method to calculate single sample enrichment score
                                            name         = "GSVA scores", # name for colorbar
                                            norm.method  = "mean", # normalization method to calculate subtype-specific enrichment score
                                            fig.name     = "upregulated_pathway_heatmap",
                                            nPerm = 10000,
                                            minGSSize = 10,
                                            maxGSSize = 500,
                                            fig.path = paste0(home, "/Results/single_algorithm/", algorithm),
                                            width = 14, height = 12)

# GSEA down-regulated
RNGversion("4.2.2")
set.seed(123)
gsea.down <- runGSEA_mod_4.4_single_algorithm(algorithm_name = algorithm,
                                              moic.res     = plot_object,
                                              dea.method   = "limma", # name of DEA method
                                              prefix       = "dgea_", # MUST be the same of argument in runDEA()
                                              dat.path      = paste0(home, "/Results/single_algorithm/", algorithm), # path of DEA files
                                              res.path      = paste0(home, "/Results/single_algorithm/", algorithm), # path to save marker files
                                              msigdb.path  = MSIGDB.FILE, # MUST be the ABSOLUTE path of msigdb file
                                              norm.expr    = plotdata$RNAseq, # use normalized expression to calculate enrichment score
                                              dirct        = "down", # direction of dysregulation in pathway
                                              n.path       = 20,
                                              p.cutoff     = 0.05, # p cutoff to identify significant pathways
                                              p.adj.cutoff = 0.05, # padj cutoff to identify significant pathways
                                              gsva.method  = "gsva", # method to calculate single sample enrichment score
                                              name         = "GSVA scores", # name for colorbar
                                              norm.method  = "mean", # normalization method to calculate subtype-specific enrichment score
                                              fig.name     = "downregulated_pathway_heatmap",
                                              nPerm = 10000,
                                              minGSSize = 10,
                                              maxGSSize = 500,
                                              fig.path = paste0(home, "/Results/single_algorithm/", algorithm),
                                              width = 14, height = 12)

# Gene set variation analysis #####
# locate ABSOLUTE path of gene set file
GSET.FILE <- paste0(home, "/Resources/Pathways/gene_sets_of_interest.gmt")

RNGversion("4.2.2")
set.seed(123)
gsva.res = runGSVA_mod_4.4_single_algorithm(algorithm_name = algorithm,
                                            moic.res      = plot_object,
                                            norm.expr     = plotdata$RNAseq,
                                            gset.gmt.path = GSET.FILE, # ABSOLUTE path of gene set file
                                            gsva.method   = "gsva", # method to calculate single sample enrichment score
                                            annCol        = annCol,
                                            annColors     = annColors,
                                            fig.path      = paste0(home, "/Results/single_algorithm/", algorithm),
                                            fig.name      = "gene_sets_of_interest_heatmap",
                                            centerFlag    = F,
                                            scaleFlag     = F,
                                            distance      = 'euclidean',
                                            linkage       = 'average',
                                            show_rownames = TRUE,
                                            show_colnames = FALSE,
                                            height        = 8,
                                            width         = 12,
                                            name          = "GSVA scores")
dev.off()

# Hierarchical clustering of pathways
library(pathfindR)
library(fastcluster)

# Get unique pathways for each subtype
GSEAfiles_up <- sort(dir(paste0(home, "/Results/single_algorithm/", algorithm), 
                         pattern = "unique_upexpr_pathway.txt$"))
GSEAfiles_down <- sort(dir(paste0(home, "/Results/single_algorithm/", algorithm), 
                           pattern = "unique_downexpr_pathway.txt$"))

unique_upexpr_pathways = list()
for (i in 1:length(gsea.up$gsea.list)) {
  unique_upexpr_pathways[[i]] = data.table::fread(paste0(paste0(home, "/Results/single_algorithm/", algorithm), 
                                                         "/", GSEAfiles_up[i]),
                                                  header = TRUE, sep = "\t")
}

unique_downexpr_pathways = list()
for (i in 1:length(gsea.down$gsea.list)) {
  unique_downexpr_pathways[[i]] = data.table::fread(paste0(paste0(home, "/Results/single_algorithm/", algorithm), 
                                                           "/", GSEAfiles_down[i]),
                                                    header = TRUE, sep = "\t")
}

names(unique_downexpr_pathways) = names(unique_upexpr_pathways) = names(gsea.up$gsea.list)

# Filter GSEA input
gsea.up_unique = gsea.up
for (i in 1:length(unique_upexpr_pathways)) {
  unq = unique_upexpr_pathways[[i]]$V1
  gsea.up_unique$gsea.list[[i]]@result = gsea.up_unique$gsea.list[[i]]@result[gsea.up_unique$gsea.list[[i]]@result$ID %in%
                                                                                unq, ]
}

gsea.down_unique = gsea.down
for (i in 1:length(unique_downexpr_pathways)) {
  unq = unique_downexpr_pathways[[i]]$V1
  gsea.down_unique$gsea.list[[i]]@result = gsea.down_unique$gsea.list[[i]]@result[gsea.down_unique$gsea.list[[i]]@result$ID %in%
                                                                                    unq, ]
}

rm(unq); gc()

hclust_input_up = prepare_gsea_output_for_hclust(gsea_output = gsea.up_unique, 
                                                 dgea_output_name_style = "dgea_", 
                                                 dea.method = "limma", 
                                                 mo.method = "",
                                                 dat.path = paste0(home, "/Results/single_algorithm/", algorithm), 
                                                 dgea_padj_cutoff = 0.05,
                                                 logfc_cutoff = 0,
                                                 pathway_padj_cutoff = 0.05)

hclust_input_down = prepare_gsea_output_for_hclust(gsea_output = gsea.down_unique, 
                                                   dgea_output_name_style = "dgea_", 
                                                   dea.method = "limma", 
                                                   mo.method = "",
                                                   dat.path = paste0(home, "/Results/single_algorithm/", algorithm), 
                                                   dgea_padj_cutoff = 0.05,
                                                   logfc_cutoff = 0,
                                                   pathway_padj_cutoff = 0.05)

hclust_input = c(hclust_input_up, hclust_input_down)
names(hclust_input) = c(paste0(rep("up_", length(hclust_input_up)), 
                               names(hclust_input_up)),
                        paste0(rep("down_", length(hclust_input_down)), 
                               names(hclust_input_down)))

rm(hclust_input_down, hclust_input_up); gc()

# Load doParallel if not already loaded
library(parallel)
library(foreach)
library(doParallel)

# Set up the number of cores to use: minimum of length(hclust_input) or 5
cl <- makeCluster(min(length(hclust_input), 5))
registerDoParallel(cl)

# Use foreach with parallel processing
timestamp()
hclust_output <- foreach(i = 1:length(hclust_input), .packages = c("pathfindR", "fastcluster")) %dopar% {
  RNGversion("4.2.2")
  set.seed(123)
  source("Scripts/automated_scripts/fast_pathfindR_hclust.R")
  
  # Perform clustering, handle errors
  result <- cluster_enriched_terms_fast(hclust_input[[i]], method = "hierarchical", plot_clusters_graph = FALSE,
                                        use_description = FALSE, use_active_snw_genes = FALSE)
  if (is.character(result) && result == "hclust impossible") {
    return("hclust impossible")
  } else {
    return(result)
  }
}

timestamp() # ~2.5 mins
stopCluster(cl)
gc()
names(hclust_output) <- names(hclust_input)

# Are there any null sets?
which(hclust_output == "hclust impossible")

# Export
library(openxlsx)
wb = createWorkbook()
for (j in 1:length(hclust_output)) {
  addWorksheet(wb, names(hclust_output)[j])
  writeData(wb, names(hclust_output)[j], hclust_output[[j]])
}
saveWorkbook(wb, file = paste0(home, "/Results/single_algorithm/", algorithm, "/", 
                               algorithm, "_representative_pathways.xlsx"),
             overwrite = TRUE); rm(wb)

# Plot pathway heatmaps
hclust_pathway_plots_up = plot_pathway_heatmaps(gsea.lists = hclust_output[grepl("up", names(hclust_output))], 
                                                norm.expr = plotdata$RNAseq, 
                                                present_clusters = c("MFA1", "MFA2"),
                                                representative = TRUE, moic.res = plot_object,
                                                subtype_prefix = algorithm, n.path = 20, msigdb.path = MSIGDB.FILE,
                                                norm.method = "mean", dirct = "up",
                                                fig.name = "upregulated_pathway_heatmap",
                                                name = "GSVA scores",
                                                fig.path = paste0(home, "/Results/single_algorithm/", algorithm), 
                                                width = 15, height = 12, gsva.method = "gsva")

hclust_pathway_plots_down = plot_pathway_heatmaps(gsea.lists = hclust_output[grepl("down", names(hclust_output))], 
                                                  norm.expr = plotdata$RNAseq, 
                                                  present_clusters = c("MFA1", "MFA2"),
                                                  representative = TRUE, moic.res = plot_object,
                                                  subtype_prefix = algorithm, n.path = 20, msigdb.path = MSIGDB.FILE,
                                                  norm.method = "mean", dirct = "down",
                                                  fig.name = "downregulated_pathway_heatmap",
                                                  name = "GSVA scores",
                                                  fig.path = paste0(home, "/Results/single_algorithm/", algorithm), 
                                                  width = 15, height = 12, gsva.method = "gsva")

# Fraction Genome Altered ###
fga_df = readRDS("Resources/TCGA/fga_df.rds"); gc()

fga.MFA <- compFGA_optimized(moic.res     = plot_object,
                             segment      = fga_df,
                             iscopynumber = TRUE, 
                             test.method  = "nonparametric", # statistical testing method (Wilcoxon with asymptotic approximation. Consider Kruskall Wallis?)
                             fig.path     = paste0(home, "/Results/single_algorithm/", algorithm),
                             fig.name     = paste0("FGA_barplot_", algorithm),
                             prefix = algorithm,
                             width = 16,
                             ga_column = "ga", # genome altered column
                             clust.col = cluster_colors,
                             title = paste0(algorithm, " FGA plot: simple criteria"))

fga.MFA.COSMIC <- compFGA_optimized(moic.res     = plot_object,
                                    segment      = fga_df,
                                    iscopynumber = TRUE, 
                                    test.method  = "nonparametric", # statistical testing method (Wilcoxon with asymptotic approximation. Consider Kruskall Wallis?)
                                    fig.path     = paste0(home, "/Results/single_algorithm/", algorithm),
                                    fig.name     = paste0("COSMIC_criteria_FGA_barplot_", algorithm),
                                    prefix = algorithm,
                                    width = 16,
                                    ga_column = "COSMIC_ga", # genome altered column
                                    clust.col = cluster_colors,
                                    title = paste0(algorithm, " FGA plot: COSMIC criteria"))

rm(fga_df); gc()

# Evaluation #####
# Run Nearest Template Prediction in transNEO cohort ###
# Load transNEO data
transNEO_mm_inputs = readRDS("Resources/transNEO/transNEO_multimodal_inputs.rds")
transcr = readRDS("Resources/transNEO/log2.norm.counts.plus1_transNEO.rds")

# get as many templates as possible
dgea.marker.up_1000 <- runMarker_single_algorithm_no_export(algorithm_name = algorithm,
                                                            moic.res = plot_object,
                                                            n.marker = 1000,
                                                            dea.method    = "limma", # name of DEA method
                                                            prefix        = "dgea_", # MUST be the same of argument in runDEA()
                                                            dat.path      = paste0(home, "/Results/single_algorithm/", algorithm), # path of DEA files
                                                            p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                                                            p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                                                            norm.expr = plotdata$RNAseq,
                                                            dirct         = "up" # direction of dysregulation in expression
)

# 2. Down-regulated markers
dgea.marker.down_1000 <- runMarker_single_algorithm_no_export(algorithm_name = algorithm,
                                                              moic.res = plot_object,
                                                              n.marker = 1000,
                                                              dea.method    = "limma", # name of DEA method
                                                              prefix        = "dgea_", # MUST be the same of argument in runDEA()
                                                              dat.path      = paste0(home, "/Results/single_algorithm/", algorithm), # path of DEA files
                                                              p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                                                              p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                                                              norm.expr = plotdata$RNAseq,
                                                              dirct         = "down" # direction of dysregulation in expression
)

# Up-regulated expression features
RNGversion("4.2.2")
timestamp()
transNEO_ntp_expr_up = runNTP(
  expr = transcr,
  templates = dgea.marker.up_1000$templates,
  scaleFlag = TRUE,
  centerFlag = TRUE,
  nPerm = 10000,
  seed = 123,
  distance = "cosine", # default
  doPlot = TRUE,
  height = 8,
  width = 12,
  fig.path = paste0(home, "/Results/single_algorithm/", algorithm),
  fig.name = "ntp_expr_up_heatmap_transNEO")
timestamp() # 12 min

# down-regulated
RNGversion("4.2.2")
timestamp()
transNEO_ntp_expr_down = runNTP(
  expr = transcr,
  templates = dgea.marker.down_1000$templates,
  scaleFlag = TRUE,
  centerFlag = TRUE,
  nPerm = 10000,
  seed = 123,
  distance = "cosine", # default
  doPlot = TRUE,
  height = 8,
  width = 12,
  fig.path = paste0(home, "/Results/single_algorithm/", algorithm),
  fig.name = "ntp_expr_down_heatmap_transNEO")
timestamp() # 12 min

# Check concordance
expr_conc = as.data.frame(transNEO_ntp_expr_down$clust.res) %>%
  dplyr::rename(clust_down = clust) %>%
  inner_join(as.data.frame(transNEO_ntp_expr_up$clust.res) %>%
               dplyr::rename(clust_up = clust),
             by = "samID")

# This is counter-intuitive but due to opposite directions of deregulation this
# is how it works (perhaps this was expected)
expr_conc$agreement = ifelse(expr_conc$clust_down!=expr_conc$clust_up, "Yes", "No")
paste("Agremeent of NTP subtypes with respect to expression data from the external cohort is: ",
      length(which(expr_conc$agreement == "Yes"))/nrow(expr_conc)*100, "% (", nrow(expr_conc),
      " samples).")

# Compare clinical variables of interest across clusters
transNEO_var2comp = transNEO_mm_inputs$`Full pheno` %>%
  dplyr::select(LN.status.at.diagnosis, ER.status, HER2.status,
                Grade.pre.NAT, pCR.RD, Age, T.stage, PAM50, iC10,
                NAT.regimen, Chemo.cycles,
                aHER2.cycles, RCB.score, STAT1.gsva,
                GGI.gsva, ESC.gsva, TMB, HRD.sum, Donor.ID) %>%
  inner_join(expr_conc %>% dplyr::select(Donor.ID = samID, MFA = clust_up),
             by = "Donor.ID")
rownames(transNEO_var2comp) = transNEO_var2comp$Donor.ID
transNEO_var2comp = transNEO_var2comp %>% dplyr::select(-Donor.ID)

# Convert to factors
transNEO_var2comp$LN.status.at.diagnosis = factor(transNEO_var2comp$LN.status.at.diagnosis,
                                                  levels = c("NEG", "POS"),
                                                  labels = c("Negative", "Positive"))
transNEO_var2comp$ER.status = factor(transNEO_var2comp$ER.status,
                                     levels = c("NEG", "POS"),
                                     labels = c("Negative", "Positive"))
transNEO_var2comp$HER2.status = factor(transNEO_var2comp$HER2.status,
                                       levels = c("NEG", "POS"),
                                       labels = c("Negative", "Positive"))
transNEO_var2comp$Grade.pre.NAT = factor(transNEO_var2comp$Grade.pre.NAT,
                                         levels = c(1, 2, 3, 4),
                                         labels = c("Grade 1", "Grade 2", "Grade 3", "Grade 4"))
transNEO_var2comp$pCR.RD = factor(transNEO_var2comp$pCR.RD,
                                  levels = c("pCR", "RD"),
                                  labels = c("pCR", "Residual Disease"))
transNEO_var2comp$PAM50 = factor(transNEO_var2comp$PAM50,
                                 levels = c("Basal", "Her2", "LumB", "LumA", "Normal", "Unk"),
                                 labels = c("Basal-like", "HER2+", "Luminal B", "Luminal A",
                                            "Normal-like", "Unknown"))
transNEO_var2comp$iC10 = factor(transNEO_var2comp$iC10,
                                levels = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
                                labels = paste("iC", seq(1, 10, 1), sep = ""))

# Remove unknown levels for statistical tests
transNEO_var2comp_nonas = transNEO_var2comp
for (i in 1:ncol(transNEO_var2comp)) {
  nas = which(transNEO_var2comp[, i] == "Unknown")
  transNEO_var2comp_nonas[nas, i] = NA
  empties = which(transNEO_var2comp[, i] == "")
  transNEO_var2comp_nonas[empties, i] = NA
}
rm(nas, empties); gc()

transNEO_clincomp = compClinvar_single_algorithm(algorithm_name = algorithm,
                                                 moic.res = transNEO_ntp_expr_up,
                                                 var2comp = transNEO_var2comp_nonas,
                                                 strata = algorithm,
                                                 factorVars = c("ER.status", "HER2.status",
                                                                "NAT.regimen", 
                                                                "pCR.RD", "LN.status.at.diagnosis"),
                                                 nonnormalVars = c("Age",
                                                                   "RCB.score", "STAT1.gsva", "GGI.gsva",
                                                                   "ESC.gsva", "TMB", "HRD.sum",
                                                                   "Grade.pre.NAT", "Chemo.cycle", "aHER2.cycles"),
                                                 includeNA = FALSE,
                                                 doWord = TRUE,
                                                 tab.name = "transNEO_Summary_of_clinical_variables",
                                                 res.path = paste0(home, "/Results/single_algorithm/", algorithm, "/"),
                                                 output_pdf = TRUE,
                                                 pdf_level_col_width = c("7em", "10em"),
                                                 pdf_count_col_width = "10em",
                                                 pdf_pval_col_width = "3em",
                                                 pdf_test_col_width = "8em",
                                                 pdf_tab_font_size = 9)

transNEO_ntp_expr_up_ord = transNEO_ntp_expr_up
transNEO_ntp_expr_up_ord$clust.res$clust = gsub(algorithm, "", transNEO_ntp_expr_up_ord$clust.res$clust)
transNEO_ordinal_clincomp = compClinvar_ordinal_single_algorithm(algorithm_name = algorithm,
                                                                 moic.res = transNEO_ntp_expr_up_ord,
                                                                 var2comp = transNEO_var2comp_nonas %>%
                                                                   dplyr::select(Grade.pre.NAT, 
                                                                                 Chemo.cycles, 
                                                                                 aHER2.cycles,
                                                                                 MFA),
                                                                 strata = algorithm,
                                                                 ordinalVars = c("Grade.pre.NAT",
                                                                                 "Chemo.cycles",
                                                                                 "aHER2.cycles"),
                                                                 includeNA = FALSE,
                                                                 tab.name = "transNEO Summary of ordinal clinical variables",
                                                                 res.path = paste0(home, "/Results/single_algorithm/", algorithm, "/"),
                                                                 output_pdf = TRUE,
                                                                 pdf_template_loc = paste0(home, "/Scripts/automated_scripts/clincomp_template.Rmd"),
                                                                 pdf_level_col_width = c("7em", "10em"),
                                                                 pdf_count_col_width = "10em",
                                                                 pdf_pval_col_width = "3em",
                                                                 pdf_test_col_width = "8em",
                                                                 pdf_tab_font_size = 9)

# Run PAM ###
RNGversion("4.2.2")
set.seed(123)
transNEO_pam = runPAM_single_algorithm(algorithm_name = algorithm,
                                       train.expr = plotdata$RNAseq,
                                       moic.res   = plot_object,
                                       test.expr  = transcr)

# Check consistency across methods

# Get predictions for TCGA (discovery cohort)
RNGversion("4.2.2")
set.seed(123)
TCGA.ntp.pred = runNTP(expr = plotdata$RNAseq[, plot_object$clust.res$samID],
                       templates = dgea.marker.up_1000$templates, distance = "cosine",
                       doPlot = F, nPerm = 10000)

TCGA.pam.pred = runPAM_single_algorithm(algorithm_name = algorithm,
                                        train.expr = plotdata$RNAseq[, plot_object$clust.res$samID],
                                        moic.res = plot_object,
                                        test.expr = plotdata$RNAseq[, plot_object$clust.res$samID])

# consensus TCGA vs NTP TCGA
runKappa_single_algorithm(algorithm_name = algorithm,
                          subt1 = plot_object$clust.res$clust,
                          subt2 = gsub(algorithm, "", TCGA.ntp.pred$clust.res$clust),
                          subt1.lab = algorithm,
                          subt2.lab = "NTP TCGA",
                          height = 8,
                          width = 8,
                          fig.path = paste0(home, "/Results/single_algorithm/", algorithm),
                          fig.name = paste0("kappa_", algorithm, "_vs_NTP_TCGA"))

# consensus TCGA vs PAM TCGA
runKappa_single_algorithm(algorithm_name = algorithm,
                          subt1 = plot_object$clust.res$clust,
                          subt2 = gsub(algorithm, "", TCGA.pam.pred$clust.res$clust),
                          subt1.lab = algorithm,
                          subt2.lab = "PAM TCGA",
                          height = 8,
                          width = 8,
                          fig.path = paste0(home, "/Results/single_algorithm/", algorithm),
                          fig.name = paste0("kappa_", algorithm, "_vs_PAM_TCGA"))

# NTP transNEO vs PAM transNEO
runKappa_single_algorithm(algorithm_name = algorithm,
                          subt1 = as.numeric(gsub(algorithm, "",
                                                  transNEO_ntp_expr_up$clust.res$clust)),
                          subt2 = as.numeric(transNEO_pam$clust.res$clust),
                          subt1.lab = "transNEO NTP",
                          subt2.lab = "transNEO PAM",
                          height = 8,
                          width = 8,
                          fig.path = paste0(home, "/Results/single_algorithm/", algorithm),
                          fig.name = "kappa_NTP_vs_PAM_transNEO")

# Supplementary results #####
# Setup for heatmaps
colors_heatmap = rev(colorRampPalette(viridisLite::magma(10))(255))
cluster_colors_heatmap = cluster_colors
clust_annot_pheno = annCol %>% mutate(Sample.ID = rownames(.)) %>%
  inner_join(clust, by = "Sample.ID") %>%
  dplyr::rename(MFA = Cluster, samID = "Sample.ID")
rownames(clust_annot_pheno) = clust_annot_pheno$samID
afh_colnames = colnames(annCol)

# Same data frame. Different columns. Just for easiness
MFA_clust_res = MFA_clusters %>% dplyr::rename(samID = Sample.ID, 
                                               MFA = Cluster) %>%
  dplyr::mutate(MFA = gsub(algorithm, "", MFA))

# PCA ###
# RNA
pca_from_original_matrix(mydata = plotdata$RNAseq, 
                         algorithm = algorithm, 
                         clust_res = MFA_clust_res,
                         cluster_colors = cluster_colors_heatmap, 
                         output_path = paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement"),
                         title_add = "RNAseq")

# miRNA
pca_from_original_matrix(mydata = plotdata$miRNA, 
                         algorithm = algorithm, 
                         clust_res = MFA_clust_res,
                         cluster_colors = cluster_colors_heatmap, 
                         output_path = paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement"),
                         title_add = "miRNA")

# CNV
pca_from_original_matrix(mydata = plotdata$CNV, 
                         algorithm = algorithm, 
                         clust_res = MFA_clust_res,
                         cluster_colors = cluster_colors_heatmap, 
                         output_path = paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement"),
                         title_add = "CNV")

# Use multidimensional scaling for SNPs
# Features must be in rows
mds_from_original_matrix(matrix = plotdata$SNPs, dist_method = "binary",
                         algorithm = algorithm, 
                         clust_res = MFA_clust_res,
                         cluster_colors = cluster_colors_heatmap, 
                         output_path = paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement"),
                         title_add = "SNPs")

# Methylation
pca_from_original_matrix(mydata = plotdata$Methylation, 
                         algorithm = algorithm, 
                         clust_res = MFA_clust_res,
                         cluster_colors = cluster_colors_heatmap,
                         output_path = paste0(home, "/Results/single_algorithm/", algorithm, "/Supplement"),
                         title_add = "Methylation")

# Draw a heatmap of the final S matrix ###
write.xlsx(cluster_DF, paste0(home, "/Results/single_algorithm/", algorithm, "/",
                              algorithm, "_embeddings.xlsx"))
dist_embeddings = as.matrix(Rfast::Dist(cluster_DF[, grep("MFA_Dim", colnames(cluster_DF))],
                                        method = "euclidean"))
dimnames(dist_embeddings) = list(cluster_DF$Sample.ID, cluster_DF$Sample.ID)
create_MO_heatmap(matrix = dist_embeddings, algorithm = algorithm, 
                  need.diag.zero = FALSE, # already zero
                  clust_annot_pheno = clust_annot_pheno ,
                  afh_colnames = afh_colnames, 
                  colors = colors_heatmap,
                  annColors = annColors,
                  heatmap_title = "MFA embeddings distance matrix",
                  cluster_colors = cluster_colors_heatmap,
                  cluster_rows_flag = FALSE,
                  cluster_cols_flag = FALSE,
                  splits_flag = TRUE,
                  legend_title = "Distance",
                  output_file_name = paste0(home, "/Results/single_algorithm/", algorithm, 
                                            "/Supplement/MFA_embeddings_distance_matrix_heatmap.png"))

# Setup for barcharts ###
# Lymph node status
scale_fill_lymph_node_status = scale_fill_manual(values = c(No = "grey75", 
                                                            Yes = "#4A0558", 
                                                            Unknown = "grey40"))

# PR status
scale_fill_PR_status = scale_fill_manual(values = c(Indeterminate = "aliceblue", 
                                                    Positive = "dodgerblue4", 
                                                    Negative = "#F0C6C3", 
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

# ER status
scale_fill_ER_status = scale_fill_manual(values = c(Negative = "#C11D9C", 
                                                    Positive = "#0F1682", 
                                                    Unknown = "grey40"))

# HER2 status
scale_fill_HER2_status = scale_fill_manual(values = c(Negative = "#0B9EF8", 
                                                      Positive = "#560DA7", 
                                                      Indeterminate = "mistyrose1", 
                                                      Equivocal = "hotpink4", 
                                                      Unknown = "grey40"))

# Stage
scale_fill_stage = scale_fill_manual(values = c(`Stage I` = "#00C9FF", 
                                                `Stage II` = "#099CF5", 
                                                `Stage III` = "#097BF5", 
                                                `Stage IV` = "#0B5684", 
                                                `Unknown` = "grey40"))

# Combine all scales into a list
barchart_scales = list(scale_fill_stage, scale_fill_lymph_node_status, scale_fill_ER_status, 
                       scale_fill_PR_status, scale_fill_HER2_status, scale_fill_vital_status, 
                       scale_fill_ethnicity, scale_fill_race, scale_fill_metastasis, 
                       scale_fill_histology, scale_fill_menopausal_status)

# Name the scales accordingly
names(barchart_scales) = c("Stage", "Lymph node status", "ER status", "PR status", "HER2 status", 
                           "Vital status", "Ethnicity", "Race", "Metastasis", "Histology", 
                           "Menopausal status")

# Chi-square tests ###
# Bias-corrected Cramer's V calculation using package rcompanion:
unbiased.cv.test = function(x, string, digits = 3) {
  CV = rcompanion::cramerV(x, bias.correct = TRUE)
  return(list(text = paste0("Bias-corrected Cramer's V / Phi for ", 
                            string, ": ", round(as.numeric(CV), digits)),
              value = round(as.numeric(CV), digits)))
}

clust_annot_pheno_nonas = clust_annot_pheno
for(i in 1:ncol(clust_annot_pheno_nonas)) {
  clust_annot_pheno_nonas[, i] = as.character(clust_annot_pheno_nonas[, i])
  nas = which(clust_annot_pheno_nonas[, i] == "Unknown")
  clust_annot_pheno_nonas[nas, i] = NA
  clust_annot_pheno_nonas[, i] = factor(clust_annot_pheno_nonas[, i])
}
rm(nas); gc()

voi = colnames(clust_annot_pheno_nonas)[1:11]
output = as.data.frame(matrix(NA, nrow = 0, ncol = 4))
for (v in 1:length(voi)){
  keepers = which(!is.na(clust_annot_pheno_nonas[, voi[v]]))
  test = suppressWarnings(chisq.test(table(clust_annot_pheno_nonas[keepers, algorithm], 
                                           clust_annot_pheno_nonas[keepers, voi[v]])))
  chifit_p = test$p.value
  chifit_xsq = test$statistic
  chifit_cv = suppressWarnings(unbiased.cv.test(table(clust_annot_pheno_nonas[keepers, algorithm], 
                                                      clust_annot_pheno_nonas[keepers, voi[v]]),
                                                string = voi[v],
                                                digits = 3)$value)
  comparison = paste0(voi[v], " vs ", algorithm, " cluster")
  output = rbind(output, c(comparison, chifit_p, chifit_xsq, chifit_cv))
  rm(test, comparison, chifit_p, chifit_xsq, chifit_cv, keepers)
}
colnames(output) = c("Comparison", "p-value", "Statistic", "Cramer's V")

rm(v); gc()
openxlsx::write.xlsx(output, 
                     paste0(home, 
                            "/Results/single_algorithm/", algorithm, "/Supplement/Chisq_tests.xlsx"),
                     overwrite = TRUE)

# Bar chart generation
MFA_barcharts = list()
plotdata_bar = clust_annot_pheno_nonas
plotdata_bar[[algorithm]] = factor(plotdata_bar[[algorithm]])
for (i in 1:length(voi)) {
  chifit = output
  loc = which(grepl(voi[i], chifit$Comparison))
  chifit = chifit[loc, ]
  MFA_barcharts[[i]] = create_annot_barchart(plotdata = plotdata_bar, fill = voi[i],
                                              chifit = chifit,
                                              na.action = "na.omit",
                                              algorithm = algorithm,
                                              barchart_ylim = 450,
                                              text_y = 400, rect_ymin = 320,
                                              rect_ymax = 420, x_annot = 2,
                                              v_gap = 25, rect_xmin = 1.5,
                                              rect_xmax = 2.5, 
                                              annot_text_size = 2.25,
                                              legend.text.size = 5,
                                              x.axis.text.size = 5) +
    barchart_scales[[voi[i]]]
  print(MFA_barcharts[[i]])
  ggsave(filename = paste0(algorithm, "_", voi[i], "_barchart.png"),
         path = paste0(home, 
                       "/Results/single_algorithm/", algorithm, "/Supplement"), 
         width = 2620, height = 2320, device = 'png', units = "px",
         dpi = 700)
  dev.off()
}
names(MFA_barcharts) = voi
rm(loc, chifit)

# Multiplot (PNG) - bar charts
library(ggpubr)
ggarrange(MFA_barcharts[[1]], MFA_barcharts[[2]], MFA_barcharts[[3]],
          MFA_barcharts[[4]], MFA_barcharts[[5]], MFA_barcharts[[6]],
          MFA_barcharts[[7]], MFA_barcharts[[8]], MFA_barcharts[[9]],
          MFA_barcharts[[10]], MFA_barcharts[[11]],
          ncol = 3, nrow = 4, labels = c("A", "B", "C", "D", "E", "F", "G", "H",
                                         "I", "J", "K"),
          font.label = list(size = 8, face = "bold", color ="black"))
ggsave(filename = paste0("Multiplot_", algorithm, "_barcharts.png"),
       path = paste0(home, 
                     "/Results/single_algorithm/", algorithm, "/Supplement"), 
       width = 9000, height = 8000, device = 'png', units = "px",
       dpi = 700)
dev.off()

# Just significant ones now
MFA_barcharts_sig = list()
plotdata_bar_sig = clust_annot_pheno_nonas %>% dplyr::select(MFA, Race, Histology, 
                                                             `ER status`, `PR status`, 
                                                             `HER2 status`, `Menopausal status`)
plotdata_bar_sig$MFA = factor(plotdata_bar_sig$MFA)
voi_sig = setdiff(colnames(plotdata_bar_sig), algorithm)
for (i in 1:length(voi_sig)) {
  chifit = output
  loc = which(grepl(voi_sig[i], chifit$Comparison))
  chifit = chifit[loc, ]
  MFA_barcharts_sig[[i]] = create_annot_barchart(plotdata = plotdata_bar_sig, fill = voi_sig[i],
                                                  chifit = chifit,
                                                  na.action = "na.omit",
                                                  algorithm = algorithm,
                                                  barchart_ylim = 450,
                                                  text_y = 400, rect_ymin = 320,
                                                  rect_ymax = 420, x_annot = 2,
                                                  v_gap = 25, rect_xmin = 1.5,
                                                  rect_xmax = 2.5, 
                                                  annot_text_size = 2.25,
                                                  legend.text.size = 5,
                                                  x.axis.text.size = 5) +
    barchart_scales[[voi_sig[i]]]
  print(MFA_barcharts_sig[[i]])
  ggsave(filename = paste0("sig_", algorithm, "_", voi_sig[i], "_barchart.png"),
         path = paste0(home, 
                       "/Results/single_algorithm/", algorithm, "/Supplement"), 
         width = 2620, height = 2320, device = 'png', units = "px",
         dpi = 700)
  dev.off()
}
names(MFA_barcharts_sig) = voi_sig
rm(loc, chifit)

# Multiplot (PNG) - bar charts
ggarrange(MFA_barcharts_sig[[1]], MFA_barcharts_sig[[2]], MFA_barcharts_sig[[3]],
          MFA_barcharts_sig[[4]], MFA_barcharts_sig[[5]], MFA_barcharts_sig[[6]],
          ncol = 2, nrow = 3, labels = c("A", "B", "C", "D", "E", "F"),
          font.label = list(size = 8, face = "bold", color ="black"))
ggsave(filename = paste0("sig_Multiplot_", algorithm, "_barcharts.png"),
       path = paste0(home, 
                     "/Results/single_algorithm/", algorithm, "/Supplement"), 
       width = 5500, height = 2320*3, device = 'png', units = "px",
       dpi = 700)
dev.off()

# Sunburst plot ###
library(plotly)
Pheno_sunburst_MFA = clust_annot_pheno
Pheno_sunburst_MFA$`ER status` = gsub("Unknown", "Unkn ER status", Pheno_sunburst_MFA$`ER status`)
Pheno_sunburst_MFA$`ER status` = gsub("Positive", "ER+", Pheno_sunburst_MFA$`ER status`)
Pheno_sunburst_MFA$`ER status` = gsub("Negative", "ER-", Pheno_sunburst_MFA$`ER status`)
Pheno_sunburst_MFA$`HER2 status` = gsub("Unknown", "Unkn HER2 status", 
                                         Pheno_sunburst_MFA$`HER2 status`)
Pheno_sunburst_MFA$`HER2 status` = gsub("Positive", "HER2+", Pheno_sunburst_MFA$`HER2 status`)
Pheno_sunburst_MFA$`HER2 status` = gsub("Negative", "HER2-", Pheno_sunburst_MFA$`HER2 status`)
Pheno_sunburst_MFA = Pheno_sunburst_MFA %>%
  dplyr::select(MFA, `ER status`, `HER2 status`) %>%
  group_by(MFA, `ER status`, `HER2 status`) %>%
  summarise(Counts = n()) %>%
  as.data.frame()

sunburst_coloring_MFA = data.frame(stringsAsFactors = FALSE,
                                    colors = tolower(gplots::col2hex(c("#2EC4B6", "#E71D36", "#FF9F1C",
                                                                       "#C11D9C", "#0F1682",  "grey40",
                                                                       "#0B9EF8", "#560DA7", "mistyrose1", 
                                                                       "hotpink4", "grey40"))),
                                    labels = c("MFA1", "MFA2",
                                               "MFA3", 
                                               "ER-", "ER+", "Unkn ER status",
                                               "HER2-", "HER2+", "Indeterminate",
                                               "Equivocal", "Unkn HER2 status"))

sunburstDF_MFA = as.sunburstDF(Pheno_sunburst_MFA, value_column = "Counts", add_root = FALSE) %>%
  inner_join(sunburst_coloring_MFA, by = "labels")

pie_MFA = plot_ly() %>%
  add_trace(ids = sunburstDF_MFA$ids, labels= sunburstDF_MFA$labels, 
            parents = sunburstDF_MFA$parents, 
            values= sunburstDF_MFA$values, type='sunburst', branchvalues = 'total',
            insidetextorientation='radial', maxdepth = 5,
            marker = list(colors = sunburstDF_MFA$colors)) %>%
  layout(
    grid = list(columns =1, rows = 1),
    margin = list(l = 0, r = 0, b = 0, t = 0)
  )
pie_MFA
rm(Pheno_sunburst_MFA, sunburstDF_MFA, sunburst_coloring_MFA, pie_MFA); gc()

# Wrap up #####
hyperparameters = list(ncp = 200,
                       conclusion = conclusion,
                       optk = optk,
                       optn.dim = optn.dim,
                       uiks = uiks
)

# Put all parameters in a list
params = list(algorithm = algorithm, data_source = data_source, data_types = data_types,
              evaluation_source = evaluation_source, title = title, subtitle = subtitle,
              description = description, in_a_nutshell = in_a_nutshell, optk_text = optk_text,
              citation = citation, NMI_to_MOVICS = NMI_to_MOVICS, ARI_to_MOVICS = ARI_to_MOVICS,
              hyperparameters = hyperparameters, ground_truth_k = ground_truth_k,
              transNEO_var2comp = transNEO_var2comp,
              sessionInfo = sessionInfo(), home = home)

# Render the R Markdown document with the parameters
rmarkdown::render(paste0(getwd(), "/Results/single_algorithm/", algorithm,
                         "/", algorithm, "_report.Rmd"), 
                  params = params, 
                  output_file = paste0(home, "/Results/single_algorithm/", 
                                       algorithm, "/", algorithm, "_report_",
                                       data_source, "_",
                                       data_types, "_eval_on_", evaluation_source,
                                       ".html"))

# Export session info as .txt
writeLines(capture.output(sessionInfo()), paste0("sessionInfo/",
                                                 algorithm, "_", data_source, "_",
                                                 data_types, "_eval_on_", evaluation_source,
                                                 "_sessionInfo.txt"))

# Save environment
save.image(paste0(home, "/Results/single_algorithm/", 
                  algorithm, "/", algorithm, "_", data_source, "_",
                  data_types, "_eval_on_", evaluation_source,
                  "_env.RData"))
