# Load environment from consensus MOVICS analysis
load("Results/MOVICS_baseline/MOVICS_TCGA_RNAseq-CNV-Methylation-miRNA-SNPs_eval_on_transNEO_env.RData")
source("Scripts/automated_scripts/modified_MOVICS_functions.R")
source("Scripts/automated_scripts/custom_functions.R")

# Libraries #####
library(MOVICS)
library(dplyr)
library(ComplexHeatmap)
library(ggplot2)
library(tidyverse)
library(wordcloud)
library(tm)

# Examine the differences in clinical variable distribution across clusterings #####
clust_annot = as.data.frame(moic.res.list$SNF$clust.res) %>%
  dplyr::rename(SNF = clust) %>%
  inner_join(as.data.frame(moic.res.list$CIMLR$clust.res) %>%
               dplyr::rename(CIMLR = clust), by = "samID") %>%
  inner_join(as.data.frame(moic.res.list$PINSPlus$clust.res) %>%
               dplyr::rename(PINSPlus = clust), by = "samID") %>%
  inner_join(as.data.frame(moic.res.list$NEMO$clust.res) %>%
               dplyr::rename(NEMO = clust), by = "samID") %>%
  inner_join(as.data.frame(moic.res.list$COCA$clust.res) %>%
               dplyr::rename(COCA = clust), by = "samID") %>%
  inner_join(as.data.frame(moic.res.list$MoCluster$clust.res) %>%
               dplyr::rename(MoCluster = clust), by = "samID") %>%
  inner_join(as.data.frame(moic.res.list$LRAcluster$clust.res) %>%
               dplyr::rename(LRAcluster = clust), by = "samID") %>%
  inner_join(as.data.frame(moic.res.list$ConsensusClustering$clust.res) %>%
               dplyr::rename(ConsensusClustering = clust), by = "samID") %>%
  inner_join(as.data.frame(moic.res.list$IntNMF$clust.res) %>%
               dplyr::rename(IntNMF = clust), by = "samID") %>%
  inner_join(as.data.frame(moic.res.list$iClusterBayes$clust.res) %>%
               dplyr::rename(iClusterBayes = clust), by = "samID")

# Convert to long format
clust_annot_long <- clust_annot %>%
  pivot_longer(cols = -samID, names_to = "Algorithm", values_to = "Cluster")
counts_by_alg = ggplot(clust_annot_long, aes(x = Algorithm, fill = as.factor(Cluster), 
                                             group = interaction(Algorithm, Cluster))) +
  geom_bar(position = position_dodge(width = 0.9), stat = "count") +
  scale_fill_manual(values = c("#2EC4B6", "#E71D36"), 
                    labels = c("Cluster 1", "Cluster 2"),
                    name = "Cluster") +
  labs(title = "Cluster Counts by Algorithm",
       x = "Clustering Algorithm",
       y = "Count") +
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.4),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        plot.title = element_text(face = "bold"))
counts_by_alg
ggsave(filename = "counts_by_alg.png",
       path = "Results/MOVICS_baseline/MO_comparisons/", 
       width = 192, height = 150, device = 'png', units = "mm",
       dpi = 700)
dev.off()

# Join pheno data with clusterings
clust_annot_pheno = clust_annot %>%
  inner_join(var2comp %>%
               mutate(samID = rownames(.)), by = "samID")

# ER status:
# Convert to long format and adjust factor levels
clust_annot_pheno_long_ER <- clust_annot_pheno %>%
  pivot_longer(cols = c("SNF", "CIMLR", "PINSPlus", "NEMO", "COCA", "MoCluster", "LRAcluster", 
                        "ConsensusClustering", "IntNMF", "iClusterBayes"),
               names_to = "Algorithm", values_to = "Cluster") %>%
  select(Algorithm, Cluster, samID, ER.status = breast_carcinoma_estrogen_receptor_status) %>%
  mutate(Cluster_ER = interaction(Cluster, ER.status, sep = " - "),
         # Ordering the levels as desired
         Cluster_ER = factor(Cluster_ER, labels = c("1 - ER+", "1 - ER-", "2 - ER+", "2 - ER-",
                                                    "1 - Unknown", "2 - Unknown"),
                             levels = c("1 - Positive", "1 - Negative", "2 - Positive", "2 - Negative",
                                        "1 - ", "2 - ")))

# Create the plot with the desired order of bars
comp_ER_clust_alg = ggplot(clust_annot_pheno_long_ER, aes(x = Algorithm, fill = Cluster_ER)) +
  geom_bar(position = position_dodge(width = 0.9), stat = "count") +
  scale_fill_manual(values = c("1 - ER+" = "#0F1682", "1 - ER-" = "#C11D9C", "1 - Unknown" = "grey40", 
                               "2 - ER+" = "#0F1682", "2 - ER-" = "#C11D9C", "2 - Unknown" = "grey40"),
                    name = "Cluster - ER status",
                    labels = c("1 - ER+", "1 - ER-", "1 - Unknown", 
                               "2 - ER+", "2 - ER-", "2 - Unknown")) +
  labs(title = "Cluster and ER Status by Algorithm",
       x = "Clustering Algorithm",
       y = "Count") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.4),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        plot.title = element_text(face = "bold"))

print(comp_ER_clust_alg)
ggsave(filename = "comp_ER_clust_alg.png",
       path = "Results/MOVICS_baseline/MO_comparisons/", 
       width = 192, height = 150, device = 'png', units = "mm",
       dpi = 700)


# Examine the differences in oncoprints #####

# Create oncoprints for each algorithm
# comp_oncoprint_full = comp_oncoprint_coding = comp_oncoprint_driver_coding = list()
# comp_oncoprint_full = comp_oncoprint_coding = list()
# 
# for (i in 1:length(moic.res.list)) {
#   # mutational frequency comparison
#   comp_oncoprint_full[[i]] <- compMut(moic.res  = moic.res.list[[i]],
#                             mut.matrix   = bmm_full, # binary somatic mutation matrix
#                             doWord       = TRUE, # generate table in .docx format
#                             doPlot       = TRUE, # draw OncoPrint
#                             freq.cutoff  = 0.05, # keep those genes that mutated in at least 5% of samples
#                             p.adj.cutoff = 0.05, # keep those genes with adjusted p value < 0.05 to draw OncoPrint
#                             innerclust   = TRUE, # perform clustering within each subtype
#                             annCol       = annCol, # same annotation for heatmap
#                             annColors    = annColors, # same annotation color for heatmap
#                             width        = 12, 
#                             height       = 6,
#                             fig.name     = paste0("oncoprint_full_", names(moic.res.list)[i]),
#                             tab.name     = paste0("Independent test between ", names(moic.res.list)[i],
#                                             " subtype and mutation"),
#                             fig.path     = "new_code/output/MOVICS/MO_comparisons/oncoprints",
#                             res.path     = "new_code/output/MOVICS/MO_comparisons/oncoprints")
#   
#   comp_oncoprint_coding[[i]] <- compMut(moic.res  = moic.res.list[[i]],
#                               mut.matrix   = bmm_coding, # binary somatic mutation matrix
#                               doWord       = TRUE, # generate table in .docx format
#                               doPlot       = TRUE, # draw OncoPrint
#                               freq.cutoff  = 0.05, # keep those genes that mutated in at least 5% of samples
#                               p.adj.cutoff = 0.05, # keep those genes with adjusted p value < 0.05 to draw OncoPrint
#                               innerclust   = TRUE, # perform clustering within each subtype
#                               annCol       = annCol, # same annotation for heatmap
#                               annColors    = annColors, # same annotation color for heatmap
#                               width        = 12, 
#                               height       = 6,
#                               fig.name     = paste0("oncoprint_coding_", names(moic.res.list)[i]),
#                               tab.name     = paste0("Independent test between ", names(moic.res.list)[i],
#                                                     " subtype and coding mutation"),
#                               fig.path     = "new_code/output/MOVICS/MO_comparisons/oncoprints",
#                               res.path     = "new_code/output/MOVICS/MO_comparisons/oncoprints")
#   
#   # driv = moic.res.list[[i]]
#   # driv$fit = driv$fit[colnames(bmm_driver_coding), colnames(bmm_driver_coding)]
#   # driv$clust.res = driv$clust.res[colnames(bmm_driver_coding), ]
#   # comp_oncoprint_driver_coding[[i]] <- compMut(moic.res  = driv,
#   #                                   mut.matrix   = bmm_driver_coding, # binary somatic mutation matrix
#   #                                   doWord       = TRUE, # generate table in .docx format
#   #                                   doPlot       = TRUE, # draw OncoPrint
#   #                                   freq.cutoff  = 0.05, # keep those genes that mutated in at least 5% of samples
#   #                                   p.adj.cutoff = 0.05, # keep those genes with adjusted p value < 0.05 to draw OncoPrint
#   #                                   innerclust   = TRUE, # perform clustering within each subtype
#   #                                   annCol       = annCol, # same annotation for heatmap
#   #                                   annColors    = annColors, # same annotation color for heatmap
#   #                                   width        = 12, 
#   #                                   height       = 6,
#   #                                   fig.name     = paste0("oncoprint_coding_driver_", names(moic.res.list)[i]),
#   #                                   tab.name     = paste0("Independent test between ", names(moic.res.list)[i],
#   #                                                         " subtype and coding driver mutation"),
#   #                                   fig.path     = "new_code/output/MOVICS/MO_comparisons/oncoprints",
#   #                                   res.path     = "new_code/output/MOVICS/MO_comparisons/oncoprints")
# }
# 
# names(comp_oncoprint_full) = names(comp_oncoprint_coding) = names(moic.res.list)

# Examine the differences in terms of gene expression and pathways across clusterings #####

# DGEA
comp_dgea = comp_dgea.marker.up = comp_dgea.marker.down = list()
# No significant genes for COCA (just 1 down-regulated, for i = 6)
for (i in 1:length(moic.res.list)) {
  comp_dgea[[i]] = runDEA(dea.method = "limma",
                          expr = input$RNAseq,
                          moic.res = moic.res.list[[i]],
                          prefix = "dgea_",
                          sort.p = TRUE,
                          overwt = TRUE,
                          verbose = TRUE,
                          res.path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/DGEA"))
  
  # Identify unique subtype biomarkers
  # 1. Up-regulated markers
  comp_dgea.marker.up[[i]] <- runMarker_single_algorithm(algorithm_name = names(moic.res.list)[i],
                                                         moic.res = moic.res.list[[i]],
                                                         dea.method    = "limma", # name of DEA method
                                                         prefix        = "dgea_", # MUST be the same of argument in runDEA()
                                                         dat.path      = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/DGEA"), # path of DEA files
                                                         res.path      = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/DGEA"), # path to save marker files
                                                         p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                                                         p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                                                         dirct         = "up", # direction of dysregulation in expression
                                                         n.marker      = 100, # number of biomarkers for each subtype
                                                         doplot        = TRUE, # generate diagonal heatmap
                                                         norm.expr     = input$RNAseq, # use normalized expression as heatmap input
                                                         annCol        = annCol, # sample annotation in heatmap
                                                         annColors     = annColors, # colors for sample annotation
                                                         show_rownames = TRUE, # show no rownames (biomarker name)
                                                         centerFlag = F,
                                                         scaleFlag = F,
                                                         halfwidth = 3,
                                                         fig.name      = paste0(names(moic.res.list)[i], 
                                                                                "_upregulated_biomarkers_heatmap"),
                                                         fig.path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/DGEA"),
                                                         width = 14,
                                                         height = 12,
                                                         fontsize_row = 3,
                                                         name = "normalised RNA-seq")
  
  # 2. Down-regulated markers
  comp_dgea.marker.down[[i]] <- runMarker_single_algorithm(algorithm_name = names(moic.res.list)[i],
                                                           moic.res = moic.res.list[[i]],
                                                           dea.method    = "limma", # name of DEA method
                                                           prefix        = "dgea_", # MUST be the same of argument in runDEA()
                                                           dat.path      = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/DGEA"), # path of DEA files
                                                           res.path      = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/DGEA"), # path to save marker files
                                                           p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                                                           p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                                                           dirct         = "down", # direction of dysregulation in expression
                                                           n.marker      = 100, # number of biomarkers for each subtype
                                                           doplot        = TRUE, # generate diagonal heatmap
                                                           norm.expr     = input$RNAseq, # use normalized expression as heatmap input
                                                           annCol        = annCol, # sample annotation in heatmap
                                                           annColors     = annColors, # colors for sample annotation
                                                           show_rownames = TRUE, # show no rownames (biomarker name)
                                                           centerFlag = F,
                                                           scaleFlag = F,
                                                           halfwidth = 3,
                                                           fig.name      = paste0(names(moic.res.list)[i], 
                                                                                  "_downregulated_biomarkers_heatmap"),
                                                           fig.path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/DGEA"),
                                                           width = 14,
                                                           height = 12,
                                                           fontsize_row = 3,
                                                           name = "normalised RNA-seq")
}

# Pathways
MSIGDB.FILE <- system.file("extdata", "c5.bp.v7.1.symbols.xls", package = "MOVICS", mustWork = TRUE)

comp_gsea.up = comp_gsea.down = list()
# Exclude COCA
for (i in c(1:4, 6:10)) {
  # GSEA up-regulated
  RNGversion("4.2.2")
  set.seed(123)
  comp_gsea.up[[i]] <- runGSEA_mod_4.4_single_algorithm(algorithm_name = names(moic.res.list)[i],
                                                        moic.res     = moic.res.list[[i]],
                                                        dea.method   = "limma", # name of DEA method
                                                        prefix       = "dgea_", # MUST be the same of argument in runDEA()
                                                        dat.path      = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/DGEA"), # path of DEA files
                                                        res.path      = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/GSEA"), # path to save marker files
                                                        msigdb.path  = MSIGDB.FILE, # MUST be the ABSOLUTE path of msigdb file
                                                        norm.expr    = input$RNAseq, # use normalized expression to calculate enrichment score
                                                        dirct        = "up", # direction of dysregulation in pathway
                                                        n.path       = 10,
                                                        p.cutoff     = 0.05, # p cutoff to identify significant pathways
                                                        p.adj.cutoff = 0.1, # padj cutoff to identify significant pathways
                                                        gsva.method  = "gsva", # method to calculate single sample enrichment score
                                                        name         = "GSVA scores", # name for colorbar
                                                        norm.method  = "mean", # normalization method to calculate subtype-specific enrichment score
                                                        fig.name     = paste0(names(moic.res.list)[i], 
                                                                              "_upregulated_pathway_heatmap"),
                                                        nPerm = 10000,
                                                        minGSSize = 10,
                                                        maxGSSize = 500,
                                                        fig.path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/GSEA"),
                                                        width = 14, height = 12)
  
  RNGversion("4.2.2")
  set.seed(123)
  comp_gsea.down[[i]] <- runGSEA_mod_4.4_single_algorithm(algorithm_name = names(moic.res.list)[i],
                                                          moic.res     = moic.res.list[[i]],
                                                          dea.method   = "limma", # name of DEA method
                                                          prefix       = "dgea_", # MUST be the same of argument in runDEA()
                                                          dat.path      = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/DGEA"), # path of DEA files
                                                          res.path      = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/GSEA"), # path to save marker files
                                                          msigdb.path  = MSIGDB.FILE, # MUST be the ABSOLUTE path of msigdb file
                                                          norm.expr    = input$RNAseq, # use normalized expression to calculate enrichment score
                                                          dirct        = "down", # direction of dysregulation in pathway
                                                          n.path       = 10,
                                                          p.cutoff     = 0.05, # p cutoff to identify significant pathways
                                                          p.adj.cutoff = 0.1, # padj cutoff to identify significant pathways
                                                          gsva.method  = "gsva", # method to calculate single sample enrichment score
                                                          name         = "GSVA scores", # name for colorbar
                                                          norm.method  = "mean", # normalization method to calculate subtype-specific enrichment score
                                                          fig.name     = paste0(names(moic.res.list)[i], 
                                                                                "_downregulated_pathway_heatmap"),
                                                          nPerm = 10000,
                                                          minGSSize = 10,
                                                          maxGSSize = 500,
                                                          fig.path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/GSEA"),
                                                          width = 14, height = 12)
}

names(comp_gsea.up) = names(comp_gsea.down) = names(moic.res.list)

# Exclude COCA slot
comp_gsea.up = comp_gsea.up[setdiff(names(comp_gsea.up), "COCA")]
comp_gsea.down = comp_gsea.down[setdiff(names(comp_gsea.down), "COCA")]

# Jaccard heatmap for pathways (see if the clusterings uncover similar biology)

# Up-regulated patwhays in CS1
cs1_upreg_path = as.data.frame(matrix(data = NA, nrow = 0, ncol = 3))
colnames(cs1_upreg_path) = c("ID", "NES", "Algorithm")
for (i in 1:length(comp_gsea.up)) {
  cs1_upreg_path_alg = as.data.frame(comp_gsea.up[[i]][["gsea.list"]][[paste0(names(comp_gsea.up)[i], 
                                                                              "_dgea__limma_test_result.CS1")]]@result) %>%
    dplyr::filter(qvalue < 0.05) %>%
    dplyr::filter(NES > 0) %>%
    dplyr::select(ID, NES) %>%
    dplyr::mutate(Algorithm = names(comp_gsea.up)[i])
  cs1_upreg_path = rbind(cs1_upreg_path, cs1_upreg_path_alg)
}
rm(cs1_upreg_path_alg)
table(cs1_upreg_path$Algorithm)

# Number of distinct up-regulated pathways
length(unique(cs1_upreg_path$ID)) # 2939

# Down-regulated pathways in CS1
cs1_downreg_path = as.data.frame(matrix(data = NA, nrow = 0, ncol = 3))
colnames(cs1_downreg_path) = c("ID", "NES", "Algorithm")
for (i in 1:length(comp_gsea.up)) {
  cs1_downreg_path_alg = as.data.frame(comp_gsea.down[[i]][["gsea.list"]][[paste0(names(comp_gsea.down)[i], 
                                                                                "_dgea__limma_test_result.CS1")]]@result) %>%
    dplyr::filter(qvalue < 0.05) %>%
    dplyr::filter(NES < 0) %>%
    dplyr::select(ID, NES) %>%
    dplyr::mutate(Algorithm = names(comp_gsea.up)[i])
  cs1_downreg_path = rbind(cs1_downreg_path, cs1_downreg_path_alg)
}
rm(cs1_downreg_path_alg)
table(cs1_downreg_path$Algorithm)

# Number of distinct down-regulated pathways
length(unique(cs1_downreg_path$ID)) # 2702

# Bar plot of top 10 most frequent pathways
top10_upreg_pathways <- cs1_upreg_path %>%
  group_by(ID) %>%
  summarize(Frequency = n(),                  
            Avg_NES = mean(NES), .groups = 'drop') %>%  
  arrange(desc(Frequency), desc(Avg_NES)) %>%  
  slice_head(n = 10)

top10_downreg_pathways <- cs1_downreg_path %>%
  group_by(ID) %>%
  summarize(Frequency = n(),                  
            Avg_NES = mean(NES), .groups = 'drop') %>%  
  arrange(desc(Frequency), desc(-Avg_NES)) %>%  
  slice_head(n = 10)

# Create a bar plot for the counts of the top 10 pathways (up-regulated)
barplot_top10_upreg_path_cs1 <- ggplot(top10_upreg_pathways, aes(x = reorder(ID, -Frequency), y = Frequency, fill = Avg_NES)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = ID), position = position_stack(vjust = 0.5),
            angle = 90, color = "black", size = 4) +
  labs(title = "Top 10 Pathways by Frequency", x = "Pathway", y = "Frequency") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.4),
        axis.text.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14)) +
  scale_fill_gradient(low = "yellow", high = "red",
                      name = "Average NES",
                      limits = range(top10_upreg_pathways$Avg_NES),
                      breaks = scales::pretty_breaks(n = 5),
                      labels = scales::number_format(accuracy = 0.1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
print(barplot_top10_upreg_path_cs1)
ggsave(filename = "top10_upregulated_pathways_CS1.png",
       path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/GSEA"), 
       width = 450, height = 390, device = 'png', units = "mm",
       dpi = 700)

# Create a bar plot for the counts of the top 10 pathways (down-regulated)
barplot_top10_downreg_path_cs1 <- ggplot(top10_downreg_pathways,
                                         aes(x = reorder(ID, -Frequency),
                                             y = Frequency, fill = Avg_NES)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = ID), position = position_stack(vjust = 0.5),
            angle = 90, color = "black", size = 4) +
  labs(title = "Top 10 Pathways by Frequency", x = "Pathway", y = "Frequency") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.4),
        axis.text.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14)) +
  scale_fill_gradient(low = "dodgerblue4", high = "skyblue",
                      name = "Average NES",
                      limits = range(top10_downreg_pathways$Avg_NES),
                      breaks = scales::pretty_breaks(n = 5),
                      labels = scales::number_format(accuracy = 0.1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
print(barplot_top10_downreg_path_cs1)
ggsave(filename = "top10_downregulated_pathways_CS1.png",
       path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/GSEA"), 
       width = 450, height = 390, device = 'png', units = "mm",
       dpi = 700)

# Create a wordcloud plot for each data frame

# Function to create wordcloud input:
preprocessText <- function(textVector, removeStopwords = FALSE) {
  # Create a text corpus
  corp <- Corpus(VectorSource(textVector))
  
  # Clean the corpus by removing punctuation, numbers, and excessive whitespace
  # corp <- tm_map(corp, removePunctuation)
  # corp <- tm_map(corp, removeNumbers)
  # corp <- tm_map(corp, stripWhitespace)
  
  # Optionally remove English stopwords
  if (removeStopwords) {
    corp <- tm_map(corp, removeWords, stopwords("english"))
  }
  
  # Create a term-document matrix
  tdm <- TermDocumentMatrix(corp)
  
  # Convert the matrix to a data frame of terms and their frequencies
  m <- as.matrix(tdm)
  termFrequency <- sort(rowSums(m), decreasing = TRUE)
  dfFrequency <- data.frame(term = names(termFrequency), freq = termFrequency)
  
  return(dfFrequency)
}

upreg.freq <- preprocessText(cs1_upreg_path$ID, removeStopwords = FALSE)
downreg.freq <- preprocessText(cs1_downreg_path$ID, removeStopwords = FALSE)

# Up-regulated pathways CS1 wordcloud
png(filename = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/GSEA/CS1_upreg_wordcloud.png"),
    width = 5500, height = 5500, res = 700)
wordcloud(words = upreg.freq$term, freq = upreg.freq$freq, min.freq = 1,
          max.words = 200, random.order = FALSE, rot.per = 0.35, 
          colors = brewer.pal(8, "Dark2"))
dev.off()

# Down-regulated pathways CS1 wordcloud
png(filename = paste0(home, "/Results/MOVICS_baseline/MO_comparisons/GSEA/CS1_downreg_wordcloud.png"),
    width = 5500, height = 5500, res = 700)
wordcloud(words = downreg.freq$term, freq = downreg.freq$freq, min.freq = 1,
          max.words = 200, random.order = FALSE, rot.per = 0.35, 
          colors = brewer.pal(8, "Dark2"))
dev.off()

# Examine similarity indices in terms of pathways across clusterings
# Group by Algorithm and collect all IDs into a list
algorithm_upreg_pathways <- cs1_upreg_path %>%
  group_by(Algorithm) %>%
  summarise(Pathways_upreg = list(ID), .groups = 'drop')

algorithm_downreg_pathways <- cs1_downreg_path %>%
  group_by(Algorithm) %>%
  summarise(Pathways_downreg = list(ID), .groups = 'drop')

algorithm_all_pathways <- rbind(cs1_upreg_path, cs1_downreg_path) %>%
  group_by(Algorithm) %>%
  summarise(Pathways_all = list(ID), .groups = 'drop')
algorithms <- unique(cs1_upreg_path$Algorithm)

# The Jaccard index is very sensitive to large differences in the sizes
# of the compared sets. Here we use the overlap coefficient which divides
# the intersection by the size of the smaller set

# Function to calculate Overlap Coefficient
overlap_coefficient <- function(set1, set2) {
  length(intersect(set1, set2)) / min(length(set1), length(set2))
}

# Function to calculate similarities using the Overlap Coefficient
calculate_similarities <- function(algorithms, pathway_column, pathway_list) {
  # Create a matrix to store the results
  similarity_matrix <- matrix(0, nrow = length(algorithms), ncol = length(algorithms),
                              dimnames = list(algorithms, algorithms))
  
  for (i in 1:length(algorithms)) {
    for (j in i:length(algorithms)) {
      # Correctly filter and extract the pathway list for each algorithm
      path1 <- pathway_list %>% 
        filter(Algorithm == algorithms[i]) %>% 
        pull({{pathway_column}})
      path2 <- pathway_list %>% 
        filter(Algorithm == algorithms[j]) %>% 
        pull({{pathway_column}})
      
      # Calculate the overlap coefficient
      similarity_score <- overlap_coefficient(path1[[1]], path2[[1]])
      
      # Assign the computed score to both symmetric positions in the matrix
      similarity_matrix[i, j] <- similarity_score
      if (i != j) {
        similarity_matrix[j, i] <- similarity_score
      }
    }
  }
  
  return(similarity_matrix)
}

# Assuming the pathway dataframes are already prepared as described
upreg_similarity_matrix <- calculate_similarities(algorithms, "Pathways_upreg", algorithm_upreg_pathways)
downreg_similarity_matrix <- calculate_similarities(algorithms, "Pathways_downreg", algorithm_downreg_pathways)
all_similarity_matrix <- calculate_similarities(algorithms, "Pathways_all", algorithm_all_pathways)

# Create heatmaps
sim_matrices = list(upreg_similarity_matrix, downreg_similarity_matrix,
                    all_similarity_matrix)
names(sim_matrices) = c("Up-regulated sets in CS1", 
                        "Down-regulated sets in CS1", 
                        "All sets in CS1")
overlap_heatmap = list()

for (k in 1:length(sim_matrices)) {
  mat = as.matrix(sim_matrices[[k]])
  class(mat) = "numeric"
  
  # Define the color palette using viridis
  color_palette <- viridis::viridis(100)
  
  # Draw the heatmap
  overlap_heatmap[[k]] = Heatmap(mat, 
                                 name = "Overlap coefficient", 
                                 column_title = paste0(names(sim_matrices)[k], " overlap heatmap"), 
                                 column_title_gp = gpar(fontsize = 8, fontface = "bold"),
                                 col = color_palette, 
                                 cluster_rows = FALSE, 
                                 cluster_columns = FALSE, 
                                 show_row_names = TRUE, 
                                 show_column_names = TRUE,
                                 row_names_gp = grid::gpar(fontsize = 6, fontface = "bold"), 
                                 column_names_gp = grid::gpar(fontsize = 6, fontface = "bold"),
                                 cell_fun = function(j, i, x, y, width, height, fill) {
                                   grid::grid.text(sprintf("%.2f", sim_matrices[[k]][i, j]), x, y, 
                                                   gp = grid::gpar(col = "black", fontsize = 6))
                                 },
                                 heatmap_legend_param = list(
                                   title = "Overlap coefficient",
                                   title_gp = grid::gpar(fontsize = 6, fontface = "bold"), 
                                   labels_gp = grid::gpar(fontsize = 6),
                                   legend_height = unit(4, "cm"),
                                   grid_width = unit(0.25, "cm"),
                                   title_position = "leftcenter-rot"
                                 ))
  
  png(paste0(home, "/Results/MOVICS_baseline/MO_comparisons/", names(sim_matrices)[k],
             " overlap heatmap.png"),
      width = 4300, height = 4300, res = 700)
  draw(overlap_heatmap[[k]])
  dev.off()
}

names(overlap_heatmap) = names(sim_matrices)

# Create a comprehensive plot that will both display similarities in:
# a) clusterings and b) identified pathways
library(reshape2)

# Average of NMI and ARI for cluster similarity
average_cluster_similarity <- (nmi_matrix + ari_matrix) / 2
overlap_matrix <- sim_matrices[["All sets in CS1"]]

# Create a long format data frame from the matrices
nmi_df <- melt(average_cluster_similarity)
colnames(nmi_df) <- c("Algorithm1", "Algorithm2", "Cluster similarity")
overlap_df <- melt(overlap_matrix)
colnames(overlap_df) <- c("Algorithm1", "Algorithm2", "Pathway overlap")

# Merge the two data frames
comparison_df <- merge(nmi_df, overlap_df, by = c("Algorithm1", "Algorithm2"))

# Plot using ggplot2
alg_biol_comp = ggplot(comparison_df, aes(x = Algorithm1, y = Algorithm2)) +
  geom_tile(aes(fill = `Cluster similarity`), color = "white") +
  geom_point(aes(size = `Pathway overlap`), shape = 16) + 
  scale_fill_gradientn(colors = rev(colorRampPalette(viridisLite::magma(10))(255))[3:225]) +
  scale_size_continuous(range = c(1, 10)) + 
  ggtitle("Algorithm comparisons: clusterings and underlying biology") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
        panel.background = element_blank(),
        legend.title = element_text(face = "bold", size = 6.5),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.background = element_rect(fill = "white", linetype = "solid"),
        legend.text = element_text(size = 5.5)) +
  labs(x = "", y = "", 
       fill = "Cluster Similarity (NMI + ARI)", shape = "Pathway overlap")
alg_biol_comp
ggsave(filename = "biology_and_clust_comparisons_across_algorithms.pdf",
       path = paste0(home, "/Results/MOVICS_baseline/MO_comparisons"), 
       width = 5300, height = 4300, device = 'pdf', units = "px",
       dpi = 700)
dev.off()
