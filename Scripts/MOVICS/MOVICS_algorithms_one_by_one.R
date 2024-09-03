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

# Chi-square tests between clusterings and clinical variables #####
chisq_outputs = list()
algorithms = colnames(clust_annot_pheno)[c(2:11)]

for (a in 1:length(algorithms)) {
  output = as.data.frame(matrix(NA, nrow = 0, ncol = 4))
  for (v in 1:length(voi)){
    test = suppressWarnings(chisq.test(table(clust_annot_pheno[, algorithms[a]], 
                                             clust_annot_pheno[, voi[v]])))
    chifit_p = test$p.value
    chifit_xsq = test$statistic
    chifit_cv = suppressWarnings(unbiased.cv.test(table(clust_annot_pheno[, algorithms[a]], 
                                                        clust_annot_pheno[, voi[v]]),
                                                  string = voi[i],
                                                  digits = 3)$value)
    comparison = paste0(voi[v], " vs ", algorithms[a], " cluster")
    output = rbind(output, c(comparison, chifit_p, chifit_xsq, chifit_cv))
    rm(test, comparison, chifit_p, chifit_xsq, chifit_cv)
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
saveWorkbook(chisq_wb, file = "new_code/output/MOVICS/MO_comparisons/chisq_tables.xlsx",
             overwrite = TRUE); rm(chisq_wb)

names(chisq_outputs) = algorithms

# SNF #####
# Normalize affinity matrix code:
normalize_affinity_matrix <- function(W) {
  N <- nrow(W)  # Determine the size of the matrix
  P <- matrix(0, nrow = N, ncol = N)  # Initialize P with zeros
  
  # Calculate row sums of W, excluding the diagonal elements
  row_sums <- rowSums(W) - diag(W)
  
  # Fill P matrix based on the conditions
  for (i in 1:N) {
    for (j in 1:N) {
      if (i != j) {
        P[i, j] = W[i, j] / (2 * row_sums[i])
      } else {
        P[i, j] = 0.5
      }
    }
  }
  
  return(P)  # Return the normalized matrix
}

aff_lymph = normalize_affinity_matrix(
  SNFtool::affinityMatrix(
    SNFtool::dist2(t(MOVICS_inputs$`Digital Pathology`),
                   t(MOVICS_inputs$`Digital Pathology`)
    ),
    K = 30, sigma = 0.5)
)
colnames(aff_lymph) = rownames(aff_lymph) = colnames(MOVICS_inputs$`Digital Pathology`)

aff_rna = normalize_affinity_matrix(
  SNFtool::affinityMatrix(
    SNFtool::dist2(t(MOVICS_inputs$RNA),
                   t(MOVICS_inputs$RNA)
    ),
    K = 30, sigma = 0.5)
)
colnames(aff_rna) = rownames(aff_rna) = colnames(MOVICS_inputs$RNA)

aff_mut = normalize_affinity_matrix(
  SNFtool::affinityMatrix(
    SNFtool::dist2(t(MOVICS_inputs$`Mutational Signatures`),
                   t(MOVICS_inputs$`Mutational Signatures`)
    ),
    K = 30, sigma = 0.5)
)
colnames(aff_mut) = rownames(aff_mut) = colnames(MOVICS_inputs$`Mutational Signatures`)

aff_immune = normalize_affinity_matrix(
  SNFtool::affinityMatrix(
    SNFtool::dist2(t(MOVICS_inputs$Immunophenoscore),
                   t(MOVICS_inputs$Immunophenoscore)
    ),
    K = 30, sigma = 0.5)
)
colnames(aff_immune) = rownames(aff_immune) = colnames(MOVICS_inputs$Immunophenoscore)

aff_final = moic.res.list[["SNF"]][["fit"]]

# Heatmaps ###
# Create heatmap for lymph data
create_MO_heatmap(matrix = aff_lymph, algorithm = "SNF", 
                  need.diag.zero = TRUE, 
                  clust_annot_pheno = clust_annot_pheno,
                  afh_colnames = afh_colnames_movics, 
                  colors = colors_heatmap,
                  heatmap_title = "DigPath first affinity heatmap",
                  cluster_colors = cluster_colors_heatmap,
                  legend_title = "Normalized affinity",
                  output_file_name = "new_code/output/MOVICS/MO_comparisons/SNF_extra/aff_lymph_heatmap.png")

# RNA
create_MO_heatmap(matrix = aff_rna, algorithm = "SNF", 
                  need.diag.zero = TRUE, 
                  clust_annot_pheno = clust_annot_pheno,
                  afh_colnames = afh_colnames_movics, 
                  colors = colors_heatmap,
                  heatmap_title = "RNA first affinity heatmap",
                  cluster_colors = cluster_colors_heatmap,
                  legend_title = "Normalized affinity",
                  output_file_name = "new_code/output/MOVICS/MO_comparisons/SNF_extra/aff_rna_heatmap.png")

# Mutational signatures
create_MO_heatmap(matrix = aff_mut, algorithm = "SNF", 
                  need.diag.zero = TRUE, 
                  clust_annot_pheno = clust_annot_pheno,
                  afh_colnames = afh_colnames_movics, 
                  colors = colors_heatmap,
                  heatmap_title = "Mutational signatures first affinity heatmap",
                  cluster_colors = cluster_colors_heatmap,
                  legend_title = "Normalized affinity",
                  output_file_name = "new_code/output/MOVICS/MO_comparisons/SNF_extra/aff_mut_heatmap.png")

# Immunophenoscore
create_MO_heatmap(matrix = aff_immune, algorithm = "SNF", 
                  need.diag.zero = TRUE, 
                  clust_annot_pheno = clust_annot_pheno,
                  afh_colnames = afh_colnames_movics, 
                  colors = colors_heatmap,
                  heatmap_title = "Immunophenoscore first affinity heatmap",
                  cluster_colors = cluster_colors_heatmap,
                  legend_title = "Normalized affinity",
                  output_file_name = "new_code/output/MOVICS/MO_comparisons/SNF_extra/aff_immmune_heatmap.png")

# Final SNF fused matrix
create_MO_heatmap(matrix = aff_final, algorithm = "SNF", 
                  need.diag.zero = TRUE, 
                  clust_annot_pheno = clust_annot_pheno,
                  afh_colnames = afh_colnames_movics, 
                  colors = colors_heatmap,
                  heatmap_title = "Final fused affinity heatmap",
                  cluster_colors = cluster_colors_heatmap,
                  legend_title = "Normalized affinity",
                  output_file_name = "new_code/output/MOVICS/MO_comparisons/SNF_extra/aff_final_heatmap.png")

# PCA plots ###
snf_clust_res = clust_annot_pheno %>% dplyr::select(samID, SNF)

# Lymph
pca_from_sim_matrix(sim_matrix = aff_lymph, algorithm = "SNF", clust_res = snf_clust_res,
                    cluster_colors = c("#2EC4B6", "#E71D36"), 
                    output_path = "new_code/output/MOVICS/MO_comparisons/SNF_extra", 
                    title_add = "lymph data")

# RNA
pca_from_sim_matrix(sim_matrix = aff_rna, algorithm = "SNF", clust_res = snf_clust_res,
                    cluster_colors = c("#2EC4B6", "#E71D36"), 
                    output_path = "new_code/output/MOVICS/MO_comparisons/SNF_extra", 
                    title_add = "RNA data")

# Mutational signatures
pca_from_sim_matrix(sim_matrix = aff_mut, algorithm = "SNF", clust_res = snf_clust_res,
                    cluster_colors = c("#2EC4B6", "#E71D36"), 
                    output_path = "new_code/output/MOVICS/MO_comparisons/SNF_extra", 
                    title_add = "Mutational Signatures")

# Immunophenoscore
pca_from_sim_matrix(sim_matrix = aff_immune, algorithm = "SNF", clust_res = snf_clust_res,
                    cluster_colors = c("#2EC4B6", "#E71D36"), 
                    output_path = "new_code/output/MOVICS/MO_comparisons/SNF_extra", 
                    title_add = "Immunophenoscore")

# Final matrix
pca_from_sim_matrix(sim_matrix = aff_final, algorithm = "SNF", clust_res = snf_clust_res,
                    cluster_colors = c("#2EC4B6", "#E71D36"), 
                    output_path = "new_code/output/MOVICS/MO_comparisons/SNF_extra", 
                    title_add = "Final Fusion")

# Bar charts with clinical variables of interest ###
SNF_barcharts = list()
cols_to_factor <- c(2:15, 18:21)
plotdata = clust_annot_pheno
plotdata[cols_to_factor] <- lapply(plotdata[cols_to_factor], as.factor)
for (i in 1:length(voi)) {
  chifit = chisq_outputs[["SNF"]]
  loc = which(grepl(voi[i], chifit$Comparison))
  chifit = chifit[loc, ]
  SNF_barcharts[[i]] = create_annot_barchart(plotdata = plotdata, fill = voi[i],
                                             chifit = chifit,
                                             algorithm = "SNF",
                                             text_y = 137, rect_ymin = 112,
                                             rect_ymax = 145) +
    barchart_scales[[voi[i]]]
  print(SNF_barcharts[[i]])
  ggsave(filename = paste0("SNF_", voi[i], "_barchart.png"),
         path = "new_code/output/MOVICS/MO_comparisons/SNF_extra", 
         width = 1920, height = 1620, device = 'png', units = "px",
         dpi = 700)
  dev.off()
}
names(SNF_barcharts) = voi
rm(loc, chifit)

# Multiplot (PNG) - bar charts
ggarrange(SNF_barcharts[[1]], SNF_barcharts[[2]], SNF_barcharts[[3]],
          SNF_barcharts[[4]], SNF_barcharts[[5]], SNF_barcharts[[6]],
          SNF_barcharts[[7]], SNF_barcharts[[8]],
          ncol = 2, nrow = 4, labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
          font.label = list(size = 8, face = "bold", color ="black"))
ggsave(filename = "Multiplot_SNF_barcharts.png",
       path = "new_code/output/MOVICS/MO_comparisons/SNF_extra", 
       width = 4612, height = 6000, device = 'png', units = "px",
       dpi = 700)
dev.off()

# Sunburst plot ###
Pheno_sunburst_SNF = clust_annot_pheno %>%
  dplyr::select(SNF, pCR.RD, PAM50, T.stage) %>%
  group_by(SNF, pCR.RD, PAM50, T.stage) %>%
  summarise(Counts = n()) %>%
  as.data.frame()
Pheno_sunburst_SNF$SNF = paste0("SNF", Pheno_sunburst_SNF$SNF)

sunburst_coloring_SNF = data.frame(stringsAsFactors = FALSE,
                                   colors = tolower(gplots::col2hex(c("#2EC4B6", "#E71D36", 
                                                                      "deeppink4", "dodgerblue4",
                                                                      "red4", "violet", "darkblue", "skyblue", "lightgreen","grey",
                                                                      "#00C9FF", "#099CF5", "#097BF5", "#0B5684"))),
                                   labels = c("SNF1", "SNF2", "RD", "pCR",
                                              "Basal", "Her2", "LumA", "LumB", "Normal", "Unk",
                                              "T1", "T2", "T3", "T4"))

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
# We are using the original S matrices because they contain fewer edges
calculate_S <- function(W) {
  n <- nrow(W)  # Assuming W is a square matrix
  S <- matrix(0, n, n)  # Initialize S as a zero matrix of the same size as W
  
  # Preserve row and column names
  rownames(S) <- rownames(W)
  colnames(S) <- colnames(W)
  
  for (i in 1:n) {
    # Find the indices of the 30 highest values in W[i, ]
    nn_indices <- order(W[i, ], decreasing = TRUE)[1:30]
    
    # Set S[i, nn_indices] to W[i, nn_indices]
    S[i, nn_indices] <- W[i, nn_indices]
  }
  
  return(S)
}

aff_lymph_S = calculate_S(aff_lymph)
aff_rna_S = calculate_S(aff_rna)
aff_mut_S = calculate_S(aff_mut)
aff_immune_S = calculate_S(aff_immune)
aff_final_S = calculate_S(aff_final)

list_aff_S = list(aff_lymph_S, aff_rna_S, aff_mut_S, aff_immune_S, aff_final_S)
names(list_aff_S) = c("Lymph Original Affinity 30-NN Graph",
                      "RNA Original Affinity 30-NN Graph",
                      "Mutational signatures Original Affinity 30-NN Graph",
                      "Immunophenoscore Original Affinity 30-NN Graph",
                      "Final Fused Affinity 30-NN Graph")

for (i in 1:length(list_aff_S)) {
  
  # Prepare the graph object
  g <- graph_from_adjacency_matrix(list_aff_S[[i]], 
                                   mode = "undirected", weighted = TRUE, diag = FALSE)
  g <- delete_edges(g, E(g)[weight == 0])
  E(g)$width <- sqrt(E(g)$weight) * 5  # Example transformation for visibility
  nodes_data <- data.frame(name = V(g)$name) %>%
    inner_join(clust_annot_pheno %>% dplyr::select(samID, SNF), by = c("name" = "samID"))
  
  # Set SNF as a factor for coloring
  nodes_data$SNF <- as.factor(nodes_data$SNF)
  V(g)$SNF <- nodes_data$SNF
  
  # Set color based on SNF
  V(g)$color <- ifelse(V(g)$SNF == 1, "#2EC4B6", "#E71D36")
  
  png(paste0("new_code/output/MOVICS/MO_comparisons/SNF_extra/",
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
       main = "")
  
  # Add title with reduced size using title() function
  title(main = names(list_aff_S)[i], cex.main = 1.7)
  
  # Add a legend to the right of the plot
  legend("bottomright", 
         title="Node Color Legend",    
         legend=c("SNF1", "SNF2"), 
         fill=c("#2EC4B6", "#E71D36"),  
         cex=0.7,      
         box.lwd=1)  
  
  dev.off() 
}
rm(g, nodes_data)

# See concordance with the final consensus
clust_annot_pheno2 = clust_annot_pheno %>%
  inner_join(as.data.frame(consensus$clust.res))
table(paste0("SNF", clust_annot_pheno2$SNF), paste0("CS", clust_annot_pheno2$clust))

# CIMLR #####

# PCA from original matrices ###

# RNA
pca_from_original_matrix(mydata = MOVICS_inputs$RNA, 
                         algorithm = "CIMLR", 
                         clust_res = clust_annot_pheno,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = "new_code/output/MOVICS/MO_comparisons/CIMLR_extra",
                         title_add = "RNA")

# Digital Pathology
pca_from_original_matrix(mydata = MOVICS_inputs$`Digital Pathology`, 
                         algorithm = "CIMLR", 
                         clust_res = clust_annot_pheno,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = "new_code/output/MOVICS/MO_comparisons/CIMLR_extra",
                         title_add = "Digital Pathology")

# Immune
pca_from_original_matrix(mydata = MOVICS_inputs$Immunophenoscore, 
                         algorithm = "CIMLR", 
                         clust_res = clust_annot_pheno,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = "new_code/output/MOVICS/MO_comparisons/CIMLR_extra",
                         title_add = "Immunophenoscore")

# Mutational signatures
pca_from_original_matrix(mydata = MOVICS_inputs$`Mutational Signatures`, 
                         algorithm = "CIMLR", 
                         clust_res = clust_annot_pheno,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = "new_code/output/MOVICS/MO_comparisons/CIMLR_extra",
                         title_add = "Mutational Signatures")

# Draw a heatmap of the final S matrix ###
cimlr_matrix = moic.res.list$CIMLR$fit$S
dimnames(cimlr_matrix) = dimnames(aff_final)
create_MO_heatmap(matrix = cimlr_matrix, algorithm = "CIMLR", 
                  need.diag.zero = FALSE, # already zero
                  clust_annot_pheno = clust_annot_pheno,
                  afh_colnames = afh_colnames_movics, 
                  colors = colors_heatmap,
                  heatmap_title = "Final CIMLR similarity heatmap",
                  cluster_colors = cluster_colors_heatmap,
                  legend_title = "Final kernel similarity",
                  output_file_name = "new_code/output/MOVICS/MO_comparisons/CIMLR_extra/CIMLR_final_S_matrix_heatmap.png")

# Bar charts with clinical variables of interest ###
CIMLR_barcharts = list()
cols_to_factor <- c(2:15, 18:21)
plotdata = clust_annot_pheno
plotdata[cols_to_factor] <- lapply(plotdata[cols_to_factor], as.factor)
for (i in 1:length(voi)) {
  chifit = chisq_outputs[["CIMLR"]]
  loc = which(grepl(voi[i], chifit$Comparison))
  chifit = chifit[loc, ]
  CIMLR_barcharts[[i]] = create_annot_barchart(plotdata = plotdata, fill = voi[i],
                                               chifit = chifit,
                                               algorithm = "CIMLR",
                                               text_y = 137, rect_ymin = 112,
                                               rect_ymax = 145) +
    barchart_scales[[voi[i]]]
  print(CIMLR_barcharts[[i]])
  ggsave(filename = paste0("CIMLR_", voi[i], "_barchart.png"),
         path = "new_code/output/MOVICS/MO_comparisons/CIMLR_extra", 
         width = 1920, height = 1620, device = 'png', units = "px",
         dpi = 700)
  dev.off()
}
names(CIMLR_barcharts) = voi
rm(loc, chifit)

# Multiplot (PNG) - bar charts
ggarrange(CIMLR_barcharts[[1]], CIMLR_barcharts[[2]], CIMLR_barcharts[[3]],
          CIMLR_barcharts[[4]], CIMLR_barcharts[[5]], CIMLR_barcharts[[6]],
          CIMLR_barcharts[[7]], CIMLR_barcharts[[8]],
          ncol = 2, nrow = 4, labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
          font.label = list(size = 8, face = "bold", color ="black"))
ggsave(filename = "Multiplot_CIMLR_barcharts.png",
       path = "new_code/output/MOVICS/MO_comparisons/CIMLR_extra", 
       width = 4612, height = 6000, device = 'png', units = "px",
       dpi = 700)
dev.off()

# Sunburst plot ###
Pheno_sunburst_CIMLR = clust_annot_pheno %>%
  dplyr::select(CIMLR, pCR.RD, PAM50, T.stage) %>%
  group_by(CIMLR, pCR.RD, PAM50, T.stage) %>%
  summarise(Counts = n()) %>%
  as.data.frame()
Pheno_sunburst_CIMLR$CIMLR = paste0("CIMLR", Pheno_sunburst_CIMLR$CIMLR)

sunburst_coloring_CIMLR = data.frame(stringsAsFactors = FALSE,
                                     colors = tolower(gplots::col2hex(c("#2EC4B6", "#E71D36", 
                                                                        "deeppink4", "dodgerblue4",
                                                                        "red4", "violet", "darkblue", "skyblue", "lightgreen","grey",
                                                                        "#00C9FF", "#099CF5", "#097BF5", "#0B5684"))),
                                     labels = c("CIMLR1", "CIMLR2", "RD", "pCR",
                                                "Basal", "Her2", "LumA", "LumB", "Normal", "Unk",
                                                "T1", "T2", "T3", "T4"))

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
table(paste0("CIMLR", clust_annot_pheno2$CIMLR), paste0("CS", clust_annot_pheno2$clust))

# PINSPlus #####
# PCA from original matrices ###
# RNA
pca_from_original_matrix(mydata = MOVICS_inputs$RNA, 
                         algorithm = "PINSPlus", 
                         clust_res = clust_annot_pheno,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = "new_code/output/MOVICS/MO_comparisons/PINSPlus_extra",
                         title_add = "RNA")

# Digital Pathology
pca_from_original_matrix(mydata = MOVICS_inputs$`Digital Pathology`, 
                         algorithm = "PINSPlus", 
                         clust_res = clust_annot_pheno,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = "new_code/output/MOVICS/MO_comparisons/PINSPlus_extra",
                         title_add = "Digital Pathology")

# Immune
pca_from_original_matrix(mydata = MOVICS_inputs$Immunophenoscore, 
                         algorithm = "PINSPlus", 
                         clust_res = clust_annot_pheno,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = "new_code/output/MOVICS/MO_comparisons/PINSPlus_extra",
                         title_add = "Immunophenoscore")

# Mutational signatures
pca_from_original_matrix(mydata = MOVICS_inputs$`Mutational Signatures`, 
                         algorithm = "PINSPlus", 
                         clust_res = clust_annot_pheno,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = "new_code/output/MOVICS/MO_comparisons/PINSPlus_extra",
                         title_add = "Mutational Signatures")

# Bar charts with clinical variables of interest ###
PINSPlus_barcharts = list()
cols_to_factor <- c(2:15, 18:21)
plotdata = clust_annot_pheno
plotdata[cols_to_factor] <- lapply(plotdata[cols_to_factor], as.factor)
for (i in 1:length(voi)) {
  chifit = chisq_outputs[["PINSPlus"]]
  loc = which(grepl(voi[i], chifit$Comparison))
  chifit = chifit[loc, ]
  PINSPlus_barcharts[[i]] = create_annot_barchart(plotdata = plotdata, fill = voi[i],
                                                  chifit = chifit,
                                                  algorithm = "PINSPlus",
                                                  text_y = 137, rect_ymin = 112,
                                                  rect_ymax = 145) +
    barchart_scales[[voi[i]]]
  print(PINSPlus_barcharts[[i]])
  ggsave(filename = paste0("PINSPlus_", voi[i], "_barchart.png"),
         path = "new_code/output/MOVICS/MO_comparisons/PINSPlus_extra", 
         width = 1920, height = 1620, device = 'png', units = "px",
         dpi = 700)
  dev.off()
}
names(PINSPlus_barcharts) = voi
rm(loc, chifit)

# Multiplot (PNG) - bar charts
ggarrange(PINSPlus_barcharts[[1]], PINSPlus_barcharts[[2]], PINSPlus_barcharts[[3]],
          PINSPlus_barcharts[[4]], PINSPlus_barcharts[[5]], PINSPlus_barcharts[[6]],
          PINSPlus_barcharts[[7]], PINSPlus_barcharts[[8]],
          ncol = 2, nrow = 4, labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
          font.label = list(size = 8, face = "bold", color ="black"))
ggsave(filename = "Multiplot_PINSPlus_barcharts.png",
       path = "new_code/output/MOVICS/MO_comparisons/PINSPlus_extra", 
       width = 4612, height = 6000, device = 'png', units = "px",
       dpi = 700)
dev.off()

# Sunburst plot ###
Pheno_sunburst_PINSPlus = clust_annot_pheno %>%
  dplyr::select(PINSPlus, pCR.RD, PAM50, T.stage) %>%
  group_by(PINSPlus, pCR.RD, PAM50, T.stage) %>%
  summarise(Counts = n()) %>%
  as.data.frame()
Pheno_sunburst_PINSPlus$PINSPlus = paste0("PINSPlus", Pheno_sunburst_PINSPlus$PINSPlus)

sunburst_coloring_PINSPlus = data.frame(stringsAsFactors = FALSE,
                                        colors = tolower(gplots::col2hex(c("#2EC4B6", "#E71D36", 
                                                                           "deeppink4", "dodgerblue4",
                                                                           "red4", "violet", "darkblue", "skyblue", "lightgreen","grey",
                                                                           "#00C9FF", "#099CF5", "#097BF5", "#0B5684"))),
                                        labels = c("PINSPlus1", "PINSPlus2", "RD", "pCR",
                                                   "Basal", "Her2", "LumA", "LumB", "Normal", "Unk",
                                                   "T1", "T2", "T3", "T4"))

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


# See concordance with the final consensus
table(paste0("PINSPlus", clust_annot_pheno2$PINSPlus), paste0("CS", clust_annot_pheno2$clust))

# NEMO #####
# PCA from original matrices ###
# RNA
pca_from_original_matrix(mydata = MOVICS_inputs$RNA, 
                         algorithm = "NEMO", 
                         clust_res = clust_annot_pheno,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = "new_code/output/MOVICS/MO_comparisons/NEMO_extra",
                         title_add = "RNA")

# Digital Pathology
pca_from_original_matrix(mydata = MOVICS_inputs$`Digital Pathology`, 
                         algorithm = "NEMO", 
                         clust_res = clust_annot_pheno,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = "new_code/output/MOVICS/MO_comparisons/NEMO_extra",
                         title_add = "Digital Pathology")

# Immune
pca_from_original_matrix(mydata = MOVICS_inputs$Immunophenoscore, 
                         algorithm = "NEMO", 
                         clust_res = clust_annot_pheno,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = "new_code/output/MOVICS/MO_comparisons/NEMO_extra",
                         title_add = "Immunophenoscore")

# Mutational signatures
pca_from_original_matrix(mydata = MOVICS_inputs$`Mutational Signatures`, 
                         algorithm = "NEMO", 
                         clust_res = clust_annot_pheno,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = "new_code/output/MOVICS/MO_comparisons/NEMO_extra",
                         title_add = "Mutational Signatures")

# Bar charts with clinical variables of interest ###
NEMO_barcharts = list()
cols_to_factor <- c(2:15, 18:21)
plotdata = clust_annot_pheno
plotdata[cols_to_factor] <- lapply(plotdata[cols_to_factor], as.factor)
for (i in 1:length(voi)) {
  chifit = chisq_outputs[["NEMO"]]
  loc = which(grepl(voi[i], chifit$Comparison))
  chifit = chifit[loc, ]
  NEMO_barcharts[[i]] = create_annot_barchart(plotdata = plotdata, fill = voi[i],
                                              chifit = chifit,
                                              algorithm = "NEMO",
                                              text_y = 137, rect_ymin = 112,
                                              rect_ymax = 145) +
    barchart_scales[[voi[i]]]
  print(NEMO_barcharts[[i]])
  ggsave(filename = paste0("NEMO_", voi[i], "_barchart.png"),
         path = "new_code/output/MOVICS/MO_comparisons/NEMO_extra", 
         width = 1920, height = 1620, device = 'png', units = "px",
         dpi = 700)
  dev.off()
}
names(NEMO_barcharts) = voi
rm(loc, chifit)

# Multiplot (PNG) - bar charts
ggarrange(NEMO_barcharts[[1]], NEMO_barcharts[[2]], NEMO_barcharts[[3]],
          NEMO_barcharts[[4]], NEMO_barcharts[[5]], NEMO_barcharts[[6]],
          NEMO_barcharts[[7]], NEMO_barcharts[[8]],
          ncol = 2, nrow = 4, labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
          font.label = list(size = 8, face = "bold", color ="black"))
ggsave(filename = "Multiplot_NEMO_barcharts.png",
       path = "new_code/output/MOVICS/MO_comparisons/NEMO_extra", 
       width = 4612, height = 6000, device = 'png', units = "px",
       dpi = 700)
dev.off()

# Sunburst plot ###
Pheno_sunburst_NEMO = clust_annot_pheno %>%
  dplyr::select(NEMO, pCR.RD, PAM50, T.stage) %>%
  group_by(NEMO, pCR.RD, PAM50, T.stage) %>%
  summarise(Counts = n()) %>%
  as.data.frame()
Pheno_sunburst_NEMO$NEMO = paste0("NEMO", Pheno_sunburst_NEMO$NEMO)

sunburst_coloring_NEMO = data.frame(stringsAsFactors = FALSE,
                                    colors = tolower(gplots::col2hex(c("#2EC4B6", "#E71D36", 
                                                                       "deeppink4", "dodgerblue4",
                                                                       "red4", "violet", "darkblue", "skyblue", "lightgreen","grey",
                                                                       "#00C9FF", "#099CF5", "#097BF5", "#0B5684"))),
                                    labels = c("NEMO1", "NEMO2", "RD", "pCR",
                                               "Basal", "Her2", "LumA", "LumB", "Normal", "Unk",
                                               "T1", "T2", "T3", "T4"))

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

# See concordance with the final consensus
table(paste0("NEMO", clust_annot_pheno2$NEMO), paste0("CS", clust_annot_pheno2$clust))

# COCA #####
# Get the final Jaccard matrix
coca_jaccard = as.matrix(as.dist(vegan::vegdist(as.matrix(moic.res.list$COCA$fit$moc), method = "jaccard")))

# Final dissimilarity matrix
create_MO_heatmap(matrix = coca_jaccard, algorithm = "COCA", 
                  need.diag.zero = TRUE, 
                  clust_annot_pheno = clust_annot_pheno,
                  afh_colnames = afh_colnames_movics, 
                  colors = colors_heatmap,
                  heatmap_title = "COCA final dissimilarity matrix heatmap",
                  cluster_colors = cluster_colors_heatmap,
                  legend_title = "Dissimilarity",
                  output_file_name = "new_code/output/MOVICS/MO_comparisons/COCA_extra/COCA_dissimilarity_heatmap.png")

# PCA from original matrices ###
# RNA
pca_from_original_matrix(mydata = MOVICS_inputs$RNA, 
                         algorithm = "COCA", 
                         clust_res = clust_annot_pheno,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = "new_code/output/MOVICS/MO_comparisons/COCA_extra",
                         title_add = "RNA")

# Digital Pathology
pca_from_original_matrix(mydata = MOVICS_inputs$`Digital Pathology`, 
                         algorithm = "COCA", 
                         clust_res = clust_annot_pheno,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = "new_code/output/MOVICS/MO_comparisons/COCA_extra",
                         title_add = "Digital Pathology")

# Immune
pca_from_original_matrix(mydata = MOVICS_inputs$Immunophenoscore, 
                         algorithm = "COCA", 
                         clust_res = clust_annot_pheno,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = "new_code/output/MOVICS/MO_comparisons/COCA_extra",
                         title_add = "Immunophenoscore")

# Mutational signatures
pca_from_original_matrix(mydata = MOVICS_inputs$`Mutational Signatures`, 
                         algorithm = "COCA", 
                         clust_res = clust_annot_pheno,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = "new_code/output/MOVICS/MO_comparisons/COCA_extra",
                         title_add = "Mutational Signatures")

# Bar charts with clinical variables of interest ###
COCA_barcharts = list()
cols_to_factor <- c(2:15, 18:21)
plotdata = clust_annot_pheno
plotdata[cols_to_factor] <- lapply(plotdata[cols_to_factor], as.factor)
for (i in 1:length(voi)) {
  chifit = chisq_outputs[["COCA"]]
  loc = which(grepl(voi[i], chifit$Comparison))
  chifit = chifit[loc, ]
  COCA_barcharts[[i]] = create_annot_barchart(plotdata = plotdata, fill = voi[i],
                                              chifit = chifit,
                                              algorithm = "COCA",
                                              text_y = 137, rect_ymin = 112,
                                              rect_ymax = 145) +
    barchart_scales[[voi[i]]]
  print(COCA_barcharts[[i]])
  ggsave(filename = paste0("COCA_", voi[i], "_barchart.png"),
         path = "new_code/output/MOVICS/MO_comparisons/COCA_extra", 
         width = 1920, height = 1620, device = 'png', units = "px",
         dpi = 700)
  dev.off()
}
names(COCA_barcharts) = voi
rm(loc, chifit)

# Multiplot (PNG) - bar charts
ggarrange(COCA_barcharts[[1]], COCA_barcharts[[2]], COCA_barcharts[[3]],
          COCA_barcharts[[4]], COCA_barcharts[[5]], COCA_barcharts[[6]],
          COCA_barcharts[[7]], COCA_barcharts[[8]],
          ncol = 2, nrow = 4, labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
          font.label = list(size = 8, face = "bold", color ="black"))
ggsave(filename = "Multiplot_COCA_barcharts.png",
       path = "new_code/output/MOVICS/MO_comparisons/COCA_extra", 
       width = 4612, height = 6000, device = 'png', units = "px",
       dpi = 700)
dev.off()

# Sunburst plot ###
Pheno_sunburst_COCA = clust_annot_pheno %>%
  dplyr::select(COCA, pCR.RD, PAM50, T.stage) %>%
  group_by(COCA, pCR.RD, PAM50, T.stage) %>%
  summarise(Counts = n()) %>%
  as.data.frame()
Pheno_sunburst_COCA$COCA = paste0("COCA", Pheno_sunburst_COCA$COCA)

sunburst_coloring_COCA = data.frame(stringsAsFactors = FALSE,
                                    colors = tolower(gplots::col2hex(c("#2EC4B6", "#E71D36", 
                                                                       "deeppink4", "dodgerblue4",
                                                                       "red4", "violet", "darkblue", "skyblue", "lightgreen","grey",
                                                                       "#00C9FF", "#099CF5", "#097BF5", "#0B5684"))),
                                    labels = c("COCA1", "COCA2", "RD", "pCR",
                                               "Basal", "Her2", "LumA", "LumB", "Normal", "Unk",
                                               "T1", "T2", "T3", "T4"))

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

# See concordance with the final consensus
table(paste0("COCA", clust_annot_pheno2$COCA), paste0("CS", clust_annot_pheno2$clust))

# MoCluster #####
mocluster_mat_2d = moic.res.list[["MoCluster"]][["fit"]]@fac.scr
dist_mocluster_2d = as.matrix(dist(mocluster_mat_2d, method = "euclidean"))

# Final 2D CPCA Euclidean distance heatmap
create_MO_heatmap(matrix = dist_mocluster_2d, algorithm = "MoCluster", 
                  need.diag.zero = TRUE, 
                  clust_annot_pheno = clust_annot_pheno,
                  afh_colnames = afh_colnames_movics, 
                  colors = colors_heatmap,
                  heatmap_title = "MoCluster 2D CPCA Euclidean distance heatmap",
                  cluster_colors = cluster_colors_heatmap,
                  legend_title = "Euclidean distance",
                  output_file_name = "new_code/output/MOVICS/MO_comparisons/MoCluster_extra/MoCluster_2D_CPCA_Euclidean_distance_heatmap.png")

# PCA from original matrices ###
# RNA
pca_from_original_matrix(mydata = MOVICS_inputs$RNA, 
                         algorithm = "MoCluster", 
                         clust_res = clust_annot_pheno,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = "new_code/output/MOVICS/MO_comparisons/MoCluster_extra",
                         title_add = "RNA")

# Digital Pathology
pca_from_original_matrix(mydata = MOVICS_inputs$`Digital Pathology`, 
                         algorithm = "MoCluster", 
                         clust_res = clust_annot_pheno,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = "new_code/output/MOVICS/MO_comparisons/MoCluster_extra",
                         title_add = "Digital Pathology")

# Immune
pca_from_original_matrix(mydata = MOVICS_inputs$Immunophenoscore, 
                         algorithm = "MoCluster", 
                         clust_res = clust_annot_pheno,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = "new_code/output/MOVICS/MO_comparisons/MoCluster_extra",
                         title_add = "Immunophenoscore")

# Mutational signatures
pca_from_original_matrix(mydata = MOVICS_inputs$`Mutational Signatures`, 
                         algorithm = "MoCluster", 
                         clust_res = clust_annot_pheno,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = "new_code/output/MOVICS/MO_comparisons/MoCluster_extra",
                         title_add = "Mutational Signatures")

# Bar charts with clinical variables of interest ###
MoCluster_barcharts = list()
cols_to_factor <- c(2:15, 18:21)
plotdata = clust_annot_pheno
plotdata[cols_to_factor] <- lapply(plotdata[cols_to_factor], as.factor)
for (i in 1:length(voi)) {
  chifit = chisq_outputs[["MoCluster"]]
  loc = which(grepl(voi[i], chifit$Comparison))
  chifit = chifit[loc, ]
  MoCluster_barcharts[[i]] = create_annot_barchart(plotdata = plotdata, fill = voi[i],
                                                   chifit = chifit,
                                                   algorithm = "MoCluster",
                                                   text_y = 137, rect_ymin = 112,
                                                   rect_ymax = 145) +
    barchart_scales[[voi[i]]]
  print(MoCluster_barcharts[[i]])
  ggsave(filename = paste0("MoCluster_", voi[i], "_barchart.png"),
         path = "new_code/output/MOVICS/MO_comparisons/MoCluster_extra", 
         width = 1920, height = 1620, device = 'png', units = "px",
         dpi = 700)
  dev.off()
}
names(MoCluster_barcharts) = voi
rm(loc, chifit)

# Multiplot (PNG) - bar charts
ggarrange(MoCluster_barcharts[[1]], MoCluster_barcharts[[2]], MoCluster_barcharts[[3]],
          MoCluster_barcharts[[4]], MoCluster_barcharts[[5]], MoCluster_barcharts[[6]],
          MoCluster_barcharts[[7]], MoCluster_barcharts[[8]],
          ncol = 2, nrow = 4, labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
          font.label = list(size = 8, face = "bold", color ="black"))
ggsave(filename = "Multiplot_MoCluster_barcharts.png",
       path = "new_code/output/MOVICS/MO_comparisons/MoCluster_extra", 
       width = 4612, height = 6000, device = 'png', units = "px",
       dpi = 700)
dev.off()

# Sunburst plot ###
Pheno_sunburst_MoCluster = clust_annot_pheno %>%
  dplyr::select(MoCluster, pCR.RD, PAM50, T.stage) %>%
  group_by(MoCluster, pCR.RD, PAM50, T.stage) %>%
  summarise(Counts = n()) %>%
  as.data.frame()
Pheno_sunburst_MoCluster$MoCluster = paste0("MoCluster", Pheno_sunburst_MoCluster$MoCluster)

sunburst_coloring_MoCluster = data.frame(stringsAsFactors = FALSE,
                                         colors = tolower(gplots::col2hex(c("#2EC4B6", "#E71D36", 
                                                                            "deeppink4", "dodgerblue4",
                                                                            "red4", "violet", "darkblue", "skyblue", "lightgreen","grey",
                                                                            "#00C9FF", "#099CF5", "#097BF5", "#0B5684"))),
                                         labels = c("MoCluster1", "MoCluster2", "RD", "pCR",
                                                    "Basal", "Her2", "LumA", "LumB", "Normal", "Unk",
                                                    "T1", "T2", "T3", "T4"))

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

# See concordance with the final consensus
table(paste0("MoCluster", clust_annot_pheno2$MoCluster), paste0("CS", clust_annot_pheno2$clust))

# LRAcluster #####
LRA_ld_coordinates = t(moic.res.list[["LRAcluster"]][["fit"]][["coordinate"]])
dist_LRA_2d = as.matrix(dist(LRA_ld_coordinates, method = "euclidean"))

# Final 2D Euclidean distance heatmap
create_MO_heatmap(matrix = dist_LRA_2d, algorithm = "LRAcluster", 
                  need.diag.zero = TRUE, 
                  clust_annot_pheno = clust_annot_pheno,
                  afh_colnames = afh_colnames_movics, 
                  colors = colors_heatmap,
                  heatmap_title = "LRAcluster 2D Euclidean distance heatmap",
                  cluster_colors = cluster_colors_heatmap,
                  legend_title = "Euclidean distance",
                  output_file_name = "new_code/output/MOVICS/MO_comparisons/LRAcluster_extra/LRAcluster_2D_Euclidean_distance_heatmap.png")

# PCA from original matrices ###
# RNA
pca_from_original_matrix(mydata = MOVICS_inputs$RNA, 
                         algorithm = "LRAcluster", 
                         clust_res = clust_annot_pheno,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = "new_code/output/MOVICS/MO_comparisons/LRAcluster_extra",
                         title_add = "RNA")

# Digital Pathology
pca_from_original_matrix(mydata = MOVICS_inputs$`Digital Pathology`, 
                         algorithm = "LRAcluster", 
                         clust_res = clust_annot_pheno,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = "new_code/output/MOVICS/MO_comparisons/LRAcluster_extra",
                         title_add = "Digital Pathology")

# Immune
pca_from_original_matrix(mydata = MOVICS_inputs$Immunophenoscore, 
                         algorithm = "LRAcluster", 
                         clust_res = clust_annot_pheno,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = "new_code/output/MOVICS/MO_comparisons/LRAcluster_extra",
                         title_add = "Immunophenoscore")

# Mutational signatures
pca_from_original_matrix(mydata = MOVICS_inputs$`Mutational Signatures`, 
                         algorithm = "LRAcluster", 
                         clust_res = clust_annot_pheno,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = "new_code/output/MOVICS/MO_comparisons/LRAcluster_extra",
                         title_add = "Mutational Signatures")

# Bar charts with clinical variables of interest ###
LRAcluster_barcharts = list()
cols_to_factor <- c(2:15, 18:21)
plotdata = clust_annot_pheno
plotdata[cols_to_factor] <- lapply(plotdata[cols_to_factor], as.factor)
for (i in 1:length(voi)) {
  chifit = chisq_outputs[["LRAcluster"]]
  loc = which(grepl(voi[i], chifit$Comparison))
  chifit = chifit[loc, ]
  LRAcluster_barcharts[[i]] = create_annot_barchart(plotdata = plotdata, fill = voi[i],
                                                    chifit = chifit,
                                                    algorithm = "LRAcluster",
                                                    text_y = 137, rect_ymin = 112,
                                                    rect_ymax = 145) +
    barchart_scales[[voi[i]]]
  print(LRAcluster_barcharts[[i]])
  ggsave(filename = paste0("LRAcluster_", voi[i], "_barchart.png"),
         path = "new_code/output/MOVICS/MO_comparisons/LRAcluster_extra", 
         width = 1920, height = 1620, device = 'png', units = "px",
         dpi = 700)
  dev.off()
}
names(LRAcluster_barcharts) = voi
rm(loc, chifit)

# Multiplot (PNG) - bar charts
ggarrange(LRAcluster_barcharts[[1]], LRAcluster_barcharts[[2]], LRAcluster_barcharts[[3]],
          LRAcluster_barcharts[[4]], LRAcluster_barcharts[[5]], LRAcluster_barcharts[[6]],
          LRAcluster_barcharts[[7]], LRAcluster_barcharts[[8]],
          ncol = 2, nrow = 4, labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
          font.label = list(size = 8, face = "bold", color ="black"))
ggsave(filename = "Multiplot_LRAcluster_barcharts.png",
       path = "new_code/output/MOVICS/MO_comparisons/LRAcluster_extra", 
       width = 4612, height = 6000, device = 'png', units = "px",
       dpi = 700)
dev.off()

# Sunburst plot ###
Pheno_sunburst_LRAcluster = clust_annot_pheno %>%
  dplyr::select(LRAcluster, pCR.RD, PAM50, T.stage) %>%
  group_by(LRAcluster, pCR.RD, PAM50, T.stage) %>%
  summarise(Counts = n()) %>%
  as.data.frame()
Pheno_sunburst_LRAcluster$LRAcluster = paste0("LRAcluster", Pheno_sunburst_LRAcluster$LRAcluster)

sunburst_coloring_LRAcluster = data.frame(stringsAsFactors = FALSE,
                                          colors = tolower(gplots::col2hex(c("#2EC4B6", "#E71D36", 
                                                                             "deeppink4", "dodgerblue4",
                                                                             "red4", "violet", "darkblue", "skyblue", "lightgreen","grey",
                                                                             "#00C9FF", "#099CF5", "#097BF5", "#0B5684"))),
                                          labels = c("LRAcluster1", "LRAcluster2", "RD", "pCR",
                                                     "Basal", "Her2", "LumA", "LumB", "Normal", "Unk",
                                                     "T1", "T2", "T3", "T4"))

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

# See concordance with the final consensus
table(paste0("LRAcluster", clust_annot_pheno2$LRAcluster), paste0("CS", clust_annot_pheno2$clust))

# Consensus Clustering #####
final_cc_matrix = 1 - moic.res.list[["ConsensusClustering"]][["fit"]][[2]][["consensusMatrix"]]
dimnames(final_cc_matrix) = list(names(moic.res.list[["ConsensusClustering"]][["fit"]][[2]][["consensusClass"]]),
                                 names(moic.res.list[["ConsensusClustering"]][["fit"]][[2]][["consensusClass"]]))

# Final 2D Euclidean distance heatmap
create_MO_heatmap(matrix = final_cc_matrix, algorithm = "ConsensusClustering", 
                  need.diag.zero = TRUE, 
                  clust_annot_pheno = clust_annot_pheno,
                  afh_colnames = afh_colnames_movics, 
                  colors = colors_heatmap,
                  heatmap_title = "Consensus Clustering final connectivity heatmap",
                  cluster_colors = cluster_colors_heatmap,
                  legend_title = "1 - final connectivity",
                  output_file_name = "new_code/output/MOVICS/MO_comparisons/CC_extra/CC_final_connectivity_heatmap.png")

# PCA from original matrices ###
# RNA
pca_from_original_matrix(mydata = MOVICS_inputs$RNA, 
                         algorithm = "ConsensusClustering", 
                         clust_res = clust_annot_pheno,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = "new_code/output/MOVICS/MO_comparisons/CC_extra",
                         title_add = "RNA")

# Digital Pathology
pca_from_original_matrix(mydata = MOVICS_inputs$`Digital Pathology`, 
                         algorithm = "ConsensusClustering", 
                         clust_res = clust_annot_pheno,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = "new_code/output/MOVICS/MO_comparisons/CC_extra",
                         title_add = "Digital Pathology")

# Immune
pca_from_original_matrix(mydata = MOVICS_inputs$Immunophenoscore, 
                         algorithm = "ConsensusClustering", 
                         clust_res = clust_annot_pheno,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = "new_code/output/MOVICS/MO_comparisons/CC_extra",
                         title_add = "Immunophenoscore")

# Mutational signatures
pca_from_original_matrix(mydata = MOVICS_inputs$`Mutational Signatures`, 
                         algorithm = "ConsensusClustering", 
                         clust_res = clust_annot_pheno,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = "new_code/output/MOVICS/MO_comparisons/CC_extra",
                         title_add = "Mutational Signatures")

# Bar charts with clinical variables of interest ###
CC_barcharts = list()
cols_to_factor <- c(2:15, 18:21)
plotdata = clust_annot_pheno
plotdata[cols_to_factor] <- lapply(plotdata[cols_to_factor], as.factor)
for (i in 1:length(voi)) {
  chifit = chisq_outputs[["ConsensusClustering"]]
  loc = which(grepl(voi[i], chifit$Comparison))
  chifit = chifit[loc, ]
  CC_barcharts[[i]] = create_annot_barchart(plotdata = plotdata, fill = voi[i],
                                            chifit = chifit,
                                            algorithm = "ConsensusClustering",
                                            text_y = 137, rect_ymin = 112,
                                            rect_ymax = 145) +
    barchart_scales[[voi[i]]]
  print(CC_barcharts[[i]])
  ggsave(filename = paste0("CC_", voi[i], "_barchart.png"),
         path = "new_code/output/MOVICS/MO_comparisons/CC_extra", 
         width = 1920, height = 1620, device = 'png', units = "px",
         dpi = 700)
  dev.off()
}
names(CC_barcharts) = voi
rm(loc, chifit)

# Multiplot (PNG) - bar charts
ggarrange(CC_barcharts[[1]], CC_barcharts[[2]], CC_barcharts[[3]],
          CC_barcharts[[4]], CC_barcharts[[5]], CC_barcharts[[6]],
          CC_barcharts[[7]], CC_barcharts[[8]],
          ncol = 2, nrow = 4, labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
          font.label = list(size = 8, face = "bold", color ="black"))
ggsave(filename = "Multiplot_CC_barcharts.png",
       path = "new_code/output/MOVICS/MO_comparisons/CC_extra", 
       width = 4612, height = 6000, device = 'png', units = "px",
       dpi = 700)
dev.off()

# Sunburst plot ###
Pheno_sunburst_CC = clust_annot_pheno %>%
  dplyr::select(ConsensusClustering, pCR.RD, PAM50, T.stage) %>%
  group_by(ConsensusClustering, pCR.RD, PAM50, T.stage) %>%
  summarise(Counts = n()) %>%
  as.data.frame()
Pheno_sunburst_CC$ConsensusClustering = paste0("CC", Pheno_sunburst_CC$ConsensusClustering)

sunburst_coloring_CC = data.frame(stringsAsFactors = FALSE,
                                  colors = tolower(gplots::col2hex(c("#2EC4B6", "#E71D36", 
                                                                     "deeppink4", "dodgerblue4",
                                                                     "red4", "violet", "darkblue", "skyblue", "lightgreen","grey",
                                                                     "#00C9FF", "#099CF5", "#097BF5", "#0B5684"))),
                                  labels = c("CC1", "CC2", "RD", "pCR",
                                             "Basal", "Her2", "LumA", "LumB", "Normal", "Unk",
                                             "T1", "T2", "T3", "T4"))

sunburstDF_CC = as.sunburstDF(Pheno_sunburst_CC, value_column = "Counts", add_root = FALSE) %>%
  inner_join(sunburst_coloring_CC, by = "labels")

pie_CC = plot_ly() %>%
  add_trace(ids = sunburstDF_CC$ids, labels= sunburstDF_CC$labels, 
            parents = sunburstDF_CC$parents, 
            values= sunburstDF_CC$values, type='sunburst', branchvalues = 'total',
            insidetextorientation='radial', maxdepth = 5,
            marker = list(colors = sunburstDF_CC$colors)) %>%
  layout(
    grid = list(columns =1, rows = 1),
    margin = list(l = 0, r = 0, b = 0, t = 0)
  )
pie_CC
rm(Pheno_sunburst_CC, sunburstDF_CC, sunburst_coloring_CC, pie_CC); gc()

# See concordance with the final consensus
table(paste0("CC", clust_annot_pheno2$ConsensusClustering), paste0("CS", clust_annot_pheno2$clust))

# IntNMF #####
# PCA from original matrices ###
# RNA
pca_from_original_matrix(mydata = MOVICS_inputs$RNA, 
                         algorithm = "IntNMF", 
                         clust_res = clust_annot_pheno,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = "new_code/output/MOVICS/MO_comparisons/IntNMF_extra",
                         title_add = "RNA")

# Digital Pathology
pca_from_original_matrix(mydata = MOVICS_inputs$`Digital Pathology`, 
                         algorithm = "IntNMF", 
                         clust_res = clust_annot_pheno,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = "new_code/output/MOVICS/MO_comparisons/IntNMF_extra",
                         title_add = "Digital Pathology")

# Immune
pca_from_original_matrix(mydata = MOVICS_inputs$Immunophenoscore, 
                         algorithm = "IntNMF", 
                         clust_res = clust_annot_pheno,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = "new_code/output/MOVICS/MO_comparisons/IntNMF_extra",
                         title_add = "Immunophenoscore")

# Mutational signatures
pca_from_original_matrix(mydata = MOVICS_inputs$`Mutational Signatures`, 
                         algorithm = "IntNMF", 
                         clust_res = clust_annot_pheno,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = "new_code/output/MOVICS/MO_comparisons/IntNMF_extra",
                         title_add = "Mutational Signatures")

# Bar charts with clinical variables of interest ###
IntNMF_barcharts = list()
cols_to_factor <- c(2:15, 18:21)
plotdata = clust_annot_pheno
plotdata[cols_to_factor] <- lapply(plotdata[cols_to_factor], as.factor)
for (i in 1:length(voi)) {
  chifit = chisq_outputs[["IntNMF"]]
  loc = which(grepl(voi[i], chifit$Comparison))
  chifit = chifit[loc, ]
  IntNMF_barcharts[[i]] = create_annot_barchart(plotdata = plotdata, fill = voi[i],
                                                chifit = chifit,
                                                algorithm = "IntNMF",
                                                text_y = 137, rect_ymin = 112,
                                                rect_ymax = 145) +
    barchart_scales[[voi[i]]]
  print(IntNMF_barcharts[[i]])
  ggsave(filename = paste0("IntNMF_", voi[i], "_barchart.png"),
         path = "new_code/output/MOVICS/MO_comparisons/IntNMF_extra", 
         width = 1920, height = 1620, device = 'png', units = "px",
         dpi = 700)
  dev.off()
}
names(IntNMF_barcharts) = voi
rm(loc, chifit)

# Multiplot (PNG) - bar charts
ggarrange(IntNMF_barcharts[[1]], IntNMF_barcharts[[2]], IntNMF_barcharts[[3]],
          IntNMF_barcharts[[4]], IntNMF_barcharts[[5]], IntNMF_barcharts[[6]],
          IntNMF_barcharts[[7]], IntNMF_barcharts[[8]],
          ncol = 2, nrow = 4, labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
          font.label = list(size = 8, face = "bold", color ="black"))
ggsave(filename = "Multiplot_IntNMF_barcharts.png",
       path = "new_code/output/MOVICS/MO_comparisons/IntNMF_extra", 
       width = 4612, height = 6000, device = 'png', units = "px",
       dpi = 700)
dev.off()

# Sunburst plot ###
Pheno_sunburst_IntNMF = clust_annot_pheno %>%
  dplyr::select(IntNMF, pCR.RD, PAM50, T.stage) %>%
  group_by(IntNMF, pCR.RD, PAM50, T.stage) %>%
  summarise(Counts = n()) %>%
  as.data.frame()
Pheno_sunburst_IntNMF$IntNMF = paste0("IntNMF", Pheno_sunburst_IntNMF$IntNMF)

sunburst_coloring_IntNMF = data.frame(stringsAsFactors = FALSE,
                                      colors = tolower(gplots::col2hex(c("#2EC4B6", "#E71D36", 
                                                                         "deeppink4", "dodgerblue4",
                                                                         "red4", "violet", "darkblue", "skyblue", "lightgreen","grey",
                                                                         "#00C9FF", "#099CF5", "#097BF5", "#0B5684"))),
                                      labels = c("IntNMF1", "IntNMF2", "RD", "pCR",
                                                 "Basal", "Her2", "LumA", "LumB", "Normal", "Unk",
                                                 "T1", "T2", "T3", "T4"))

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

# See concordance with the final consensus
table(paste0("IntNMF", clust_annot_pheno2$IntNMF), paste0("CS", clust_annot_pheno2$clust))

# iClusterBayes #####
iCB_feature_ranks = as.data.frame(moic.res.list[["iClusterBayes"]][["feat.res"]])
write.xlsx(iCB_feature_ranks,
           "new_code/output/MOVICS/MO_comparisons/iClusterBayes_extra/iCB_feature_ranks.xlsx",
           overwrite = TRUE)

# PCA from original matrices ###
# RNA
pca_from_original_matrix(mydata = MOVICS_inputs$RNA, 
                         algorithm = "iClusterBayes", 
                         clust_res = clust_annot_pheno,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = "new_code/output/MOVICS/MO_comparisons/iClusterBayes_extra",
                         title_add = "RNA")

# Digital Pathology
pca_from_original_matrix(mydata = MOVICS_inputs$`Digital Pathology`, 
                         algorithm = "iClusterBayes", 
                         clust_res = clust_annot_pheno,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = "new_code/output/MOVICS/MO_comparisons/iClusterBayes_extra",
                         title_add = "Digital Pathology")

# Immune
pca_from_original_matrix(mydata = MOVICS_inputs$Immunophenoscore, 
                         algorithm = "iClusterBayes", 
                         clust_res = clust_annot_pheno,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = "new_code/output/MOVICS/MO_comparisons/iClusterBayes_extra",
                         title_add = "Immunophenoscore")

# Mutational signatures
pca_from_original_matrix(mydata = MOVICS_inputs$`Mutational Signatures`, 
                         algorithm = "iClusterBayes", 
                         clust_res = clust_annot_pheno,
                         cluster_colors = c("#2EC4B6", "#E71D36"), 
                         output_path = "new_code/output/MOVICS/MO_comparisons/iClusterBayes_extra",
                         title_add = "Mutational Signatures")

# Bar charts with clinical variables of interest ###
iClusterBayes_barcharts = list()
cols_to_factor <- c(2:15, 18:21)
plotdata = clust_annot_pheno
plotdata[cols_to_factor] <- lapply(plotdata[cols_to_factor], as.factor)
for (i in 1:length(voi)) {
  chifit = chisq_outputs[["iClusterBayes"]]
  loc = which(grepl(voi[i], chifit$Comparison))
  chifit = chifit[loc, ]
  iClusterBayes_barcharts[[i]] = create_annot_barchart(plotdata = plotdata, fill = voi[i],
                                                       chifit = chifit,
                                                       algorithm = "iClusterBayes",
                                                       text_y = 137, rect_ymin = 112,
                                                       rect_ymax = 145) +
    barchart_scales[[voi[i]]]
  print(iClusterBayes_barcharts[[i]])
  ggsave(filename = paste0("iClusterBayes_", voi[i], "_barchart.png"),
         path = "new_code/output/MOVICS/MO_comparisons/iClusterBayes_extra", 
         width = 1920, height = 1620, device = 'png', units = "px",
         dpi = 700)
  dev.off()
}
names(iClusterBayes_barcharts) = voi
rm(loc, chifit)

# Multiplot (PNG) - bar charts
ggarrange(iClusterBayes_barcharts[[1]], iClusterBayes_barcharts[[2]], iClusterBayes_barcharts[[3]],
          iClusterBayes_barcharts[[4]], iClusterBayes_barcharts[[5]], iClusterBayes_barcharts[[6]],
          iClusterBayes_barcharts[[7]], iClusterBayes_barcharts[[8]],
          ncol = 2, nrow = 4, labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
          font.label = list(size = 8, face = "bold", color ="black"))
ggsave(filename = "Multiplot_iClusterBayes_barcharts.png",
       path = "new_code/output/MOVICS/MO_comparisons/iClusterBayes_extra", 
       width = 4612, height = 6000, device = 'png', units = "px",
       dpi = 700)
dev.off()

# Sunburst plot ###
Pheno_sunburst_iClusterBayes = clust_annot_pheno %>%
  dplyr::select(iClusterBayes, pCR.RD, PAM50, T.stage) %>%
  group_by(iClusterBayes, pCR.RD, PAM50, T.stage) %>%
  summarise(Counts = n()) %>%
  as.data.frame()
Pheno_sunburst_iClusterBayes$iClusterBayes = paste0("iClusterBayes", Pheno_sunburst_iClusterBayes$iClusterBayes)

sunburst_coloring_iClusterBayes = data.frame(stringsAsFactors = FALSE,
                                             colors = tolower(gplots::col2hex(c("#2EC4B6", "#E71D36", 
                                                                                "deeppink4", "dodgerblue4",
                                                                                "red4", "violet", "darkblue", "skyblue", "lightgreen","grey",
                                                                                "#00C9FF", "#099CF5", "#097BF5", "#0B5684"))),
                                             labels = c("iClusterBayes1", "iClusterBayes2", "RD", "pCR",
                                                        "Basal", "Her2", "LumA", "LumB", "Normal", "Unk",
                                                        "T1", "T2", "T3", "T4"))

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

# See concordance with the final consensus
table(paste0("iClusterBayes", clust_annot_pheno2$iClusterBayes), paste0("CS", clust_annot_pheno2$clust))
