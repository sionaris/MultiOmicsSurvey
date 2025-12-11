# MultiOmicsSurvey Function Documentation

This document provides comprehensive documentation for all custom functions used in the MultiOmicsSurvey project. Functions are organized by source file and presented in alphabetical order within each section.

---

## Table of Contents

- [custom_functions.R](#custom_functionsr)
- [fast_pathfindR_hclust.R](#fast_pathfindr_hclustr)
- [modified_MOVICS_functions.R](#modified_movics_functionsr)
- [Spectrum_functions.R](#spectrum_functionsr)
- [Alphabetical Index](#alphabetical-index)

---

## custom_functions.R

Functions for multi-omics data analysis, clustering evaluation, visualization, and data processing.

### as.sunburstDF

```r
as.sunburstDF(DF, value_column = NULL, add_root = FALSE)
```

**Description:** Converts a hierarchical data frame into a sunburst-compatible format for visualization.

| Parameter | Type | Description |
|-----------|------|-------------|
| `DF` | data.frame | Hierarchical data frame with category columns |
| `value_column` | character | Column name containing values for sizing |
| `add_root` | logical | Whether to add a root node (default: FALSE) |

**Returns:** A data frame formatted for sunburst charts with `ids`, `labels`, `parents`, and `values` columns.

---

### calculate_ari_index

```r
calculate_ari_index(cluster_df1, cluster_df2, 
                    cluster_col1 = "Cluster", cluster_col2 = "Cluster",
                    sample_col1 = "samID", sample_col2 = "samID")
```

**Description:** Calculates the Adjusted Rand Index (ARI) between two clustering results.

| Parameter | Type | Description |
|-----------|------|-------------|
| `cluster_df1` | data.frame | First clustering result |
| `cluster_df2` | data.frame | Second clustering result |
| `cluster_col1` | character | Column name for clusters in df1 (default: "Cluster") |
| `cluster_col2` | character | Column name for clusters in df2 (default: "Cluster") |
| `sample_col1` | character | Column name for sample IDs in df1 (default: "samID") |
| `sample_col2` | character | Column name for sample IDs in df2 (default: "samID") |

**Returns:** Numeric value representing the Adjusted Rand Index (-1 to 1).

---

### calculate_nmi_index

```r
calculate_nmi_index(cluster_df1, cluster_df2,
                    cluster_col1 = "Cluster", cluster_col2 = "Cluster",
                    sample_col1 = "samID", sample_col2 = "samID")
```

**Description:** Calculates the Normalized Mutual Information (NMI) between two clustering results.

**Returns:** Numeric value representing NMI (0 to 1).

---

### calculate_S

```r
calculate_S(W)
```

**Description:** Calculates the S matrix from an affinity matrix W using the SNF transition probability formula.

| Parameter | Type | Description |
|-----------|------|-------------|
| `W` | matrix | Affinity/similarity matrix |

**Returns:** Normalized transition probability matrix S.

---

### CIMLR_mod

```r
CIMLR_mod(X, c, no.dim = NA, k = 10, cores.ratio = 1, binary_flags = NULL,
          binary_distance = "binary", max_iter = 100, use_rsvd = FALSE)
```

**Description:** Modified CIMLR (Cancer Integration via Multi-kernel Learning) algorithm supporting mixed data types.

| Parameter | Type | Description |
|-----------|------|-------------|
| `X` | list | List of data matrices (features × samples) |
| `c` | integer | Number of clusters |
| `no.dim` | integer | Number of dimensions for embedding (default: NA) |
| `k` | integer | Number of neighbors (default: 10) |
| `cores.ratio` | numeric | Ratio of cores for parallel processing (default: 1) |
| `binary_flags` | character | Vector indicating binary data types |
| `binary_distance` | character | Distance metric for binary data (default: "binary") |

**Returns:** List containing cluster assignments, similarity matrix, and convergence info.

---

### CIMLR_Estimate_Number_of_Clusters_weight_mod

```r
CIMLR_Estimate_Number_of_Clusters_weight_mod(all_data, NUMC = 2:5, 
                                              cores.ratio = 0, weight,
                                              binary_flags = NULL, binary_distance = NULL)
```

**Description:** Estimates optimal number of clusters for weighted CIMLR using eigengap and rotation cost.

**Returns:** List with optimal K estimates from eigengap and rotation methods.

---

### CIMLR.weight_mod

```r
CIMLR.weight_mod(X, c, no.dim = NA, k = 10, cores.ratio = 0, weight,
                 binary_flags = NULL, binary_distance = NULL, max_iter = 100)
```

**Description:** Weighted CIMLR implementation with feature importance weights.

**Returns:** Clustering result with weighted kernel integration.

---

### coca_cc_mod

```r
coca_cc_mod(data = NULL, K = 2, B = 100, pItem = 0.8, clMethod = "hclust",
            innerLinkage = "complete", maxK = NULL)
```

**Description:** Modified COCA (Cluster of Clusters Analysis) consensus clustering.

| Parameter | Type | Description |
|-----------|------|-------------|
| `data` | list | List of clustering matrices |
| `K` | integer | Number of final clusters (default: 2) |
| `B` | integer | Number of bootstrap iterations (default: 100) |
| `pItem` | numeric | Proportion of items to sample (default: 0.8) |
| `clMethod` | character | Clustering method (default: "hclust") |

**Returns:** Consensus clustering result with cluster assignments.

---

### compute_matrix_similarity

```r
compute_matrix_similarity(matrices)
```

**Description:** Computes pairwise similarity between multiple matrices using Frobenius norm and Pearson correlation.

| Parameter | Type | Description |
|-----------|------|-------------|
| `matrices` | list | Named list of matrices to compare |

**Returns:** List with Frobenius and Pearson similarity matrices.

---

### compute_silhouette

```r
compute_silhouette(cluster_df, similarity_matrix, normalize_matrix = FALSE)
```

**Description:** Computes silhouette scores from a similarity matrix and cluster assignments.

| Parameter | Type | Description |
|-----------|------|-------------|
| `cluster_df` | data.frame | Data frame with sample IDs and cluster assignments |
| `similarity_matrix` | matrix | Sample similarity matrix |
| `normalize_matrix` | logical | Whether to normalize the matrix (default: FALSE) |

**Returns:** Silhouette object with per-sample scores.

---

### concordanceNetworkNMI_ANF

```r
concordanceNetworkNMI_ANF(Wall, C, type)
```

**Description:** Calculates concordance between networks and clustering using NMI for ANF method.

**Returns:** NMI-based concordance scores.

---

### create_annot_barchart

```r
create_annot_barchart(plotdata = NULL, fill = NULL, title = NULL,
                      xlab = NULL, ylab = NULL, legend.position = "right")
```

**Description:** Creates annotated bar charts for visualizing cluster distributions.

**Returns:** ggplot2 object.

---

### create_density_curve

```r
create_density_curve(matrix)
```

**Description:** Creates a density curve plot for matrix values.

**Returns:** ggplot2 density plot object.

---

### create_density_plot_color

```r
create_density_plot_color(matrix)
```

**Description:** Creates a colored density plot for matrix value distributions.

**Returns:** ggplot2 density plot with color gradient.

---

### create_MO_heatmap

```r
create_MO_heatmap(matrix = NULL, algorithm = NULL, clust_res = NULL,
                  data_name = NULL, fig.path = getwd(), ...)
```

**Description:** Creates multi-omics heatmaps with cluster annotations.

**Returns:** ComplexHeatmap object saved to file.

---

### dist2.cimlr_mod2

```r
dist2.cimlr_mod2(x, c = NA, name = NULL, method = "sqeuclidean")
```

**Description:** Modified squared Euclidean distance calculation for CIMLR.

**Returns:** Distance matrix.

---

### dist2.cimlr.weight_mod

```r
dist2.cimlr.weight_mod(x, weight, method = "sqeuclidean")
```

**Description:** Weighted distance calculation for CIMLR.

**Returns:** Weighted distance matrix.

---

### divide_by_quantile_normalize

```r
divide_by_quantile_normalize(matrix, norm_quant = 0.05)
```

**Description:** Normalizes matrix by dividing by a quantile value.

**Returns:** Normalized matrix.

---

### estimateNumberOfClustersGivenGraph_mod

```r
estimateNumberOfClustersGivenGraph_mod(W, NUMC = 2:5)
```

**Description:** Estimates optimal cluster number using eigengap and rotation cost heuristics.

**Returns:** List with optimal K estimates.

---

### fetch_citation

```r
fetch_citation(algorithm)
```

**Description:** Retrieves citation information for a given algorithm.

**Returns:** Character string with citation.

---

### fetch_in_a_nutshell

```r
fetch_in_a_nutshell(algorithm)
```

**Description:** Retrieves brief description for a given algorithm.

**Returns:** Character string with algorithm summary.

---

### frobenius_norm

```r
frobenius_norm(mat1, mat2)
```

**Description:** Calculates the Frobenius norm distance between two matrices.

**Returns:** Numeric Frobenius norm value.

---

### getMoHeatmap_prenorm_quart

```r
getMoHeatmap_prenorm_quart(data = NULL, is.binary = c(FALSE, FALSE, FALSE, FALSE, FALSE),
                           clust.res = NULL, ...)
```

**Description:** Generates multi-omics heatmap with quartile-based pre-normalization.

**Returns:** ComplexHeatmap object.

---

### linear_quantile_normalize

```r
linear_quantile_normalize(matrix, norm_quant = 0.05)
```

**Description:** Linear quantile normalization scaling matrix to [0,1] based on quantiles.

**Returns:** Normalized matrix.

---

### mds_from_original_matrix

```r
mds_from_original_matrix(matrix = NULL, dist_method = NULL, algorithm = NULL,
                         clust_res = NULL, fig.path = getwd(), ...)
```

**Description:** Performs MDS from original data matrix for visualization.

**Returns:** MDS plot saved to file.

---

### MOVICS_jaccard_index

```r
MOVICS_jaccard_index(clust1, clust2)
```

**Description:** Calculates Jaccard index between two cluster assignments.

**Returns:** Numeric Jaccard index (0 to 1).

---

### multiple.kernel.cimlr_mod

```r
multiple.kernel.cimlr_mod(x, cores.ratio = 1, name = NULL, method = "sqeuclidean")
```

**Description:** Generates multiple kernels for CIMLR integration.

**Returns:** List of kernel matrices.

---

### multiple.kernel.cimlr.weight_mod

```r
multiple.kernel.cimlr.weight_mod(x, cores.ratio = 0, weight, method = "sqeuclidean")
```

**Description:** Weighted multiple kernel generation for CIMLR.

**Returns:** Weighted kernel list.

---

### nemo.affinity.graph_mod

```r
nemo.affinity.graph_mod(raw.data, k = NA, sigma = 0.5,
                        binary_flags = rep("No", length(raw.data)),
                        binary_distance = NULL)
```

**Description:** Modified NEMO affinity graph construction supporting mixed data types.

**Returns:** Affinity matrix.

---

### normalize_affinity_matrix

```r
normalize_affinity_matrix(W)
```

**Description:** Normalizes an affinity matrix to have row sums of 1.

**Returns:** Row-normalized affinity matrix.

---

### pca_from_original_matrix

```r
pca_from_original_matrix(mydata = NULL, algorithm = NULL, clust_res = NULL,
                         data_name = NULL, fig.path = getwd(), ...)
```

**Description:** Performs PCA on original data and creates visualization with cluster colors.

**Returns:** PCA plot saved to file.

---

### pca_from_sim_matrix

```r
pca_from_sim_matrix(sim_matrix = NULL, algorithm = NULL, clust_res = NULL,
                    data_name = NULL, fig.path = getwd(), ...)
```

**Description:** Performs PCA-based visualization from a similarity matrix.

**Returns:** PCA plot from similarity matrix.

---

### pearson_correlation

```r
pearson_correlation(mat1, mat2)
```

**Description:** Calculates Pearson correlation between vectorized matrices.

**Returns:** Pearson correlation coefficient.

---

### plot_MOFA_cumul_var_exp

```r
plot_MOFA_cumul_var_exp(object, factor_range = NULL, ...)
```

**Description:** Plots cumulative variance explained by MOFA factors.

**Returns:** ggplot2 cumulative variance plot.

---

### plot_MOFA_var_exp

```r
plot_MOFA_var_exp(object, factor_range = NULL, ...)
```

**Description:** Plots variance explained per view and factor for MOFA models.

**Returns:** ggplot2 variance explained heatmap.

---

### rankFeaturesByNMI_parallely

```r
rankFeaturesByNMI_parallely(data, W, ncores = detectCores() - 1, 
                            binary = FALSE, nn = NULL, sigma)
```

**Description:** Ranks features by NMI contribution using parallel processing.

**Returns:** Data frame with feature NMI scores.

---

### rankFeaturesByNMI_parallely_ANF

```r
rankFeaturesByNMI_parallely_ANF(data, W, ncores = detectCores() - 1,
                                 binary = FALSE, ...)
```

**Description:** ANF-specific feature ranking by NMI.

**Returns:** Feature ranking data frame.

---

### reshape_Pearson_for_tests

```r
reshape_Pearson_for_tests(similarities)
```

**Description:** Reshapes Pearson correlation matrix for statistical testing.

**Returns:** Reshaped data frame for testing.

---

### standardize_rows

```r
standardize_rows(mat)
```

**Description:** Z-score standardizes each row of a matrix.

**Returns:** Row-standardized matrix.

---

### TCGAanalyze_survival_custom3

```r
TCGAanalyze_survival_custom3(clinical_data, cluster_col, time_col, 
                              event_col, algorithm = NULL, ...)
```

**Description:** Custom survival analysis for TCGA data with multi-group comparisons.

**Returns:** Survival analysis results with Kaplan-Meier plots.

---

## fast_pathfindR_hclust.R

Optimized functions for pathway enrichment term clustering.

### cluster_enriched_terms_fast

```r
cluster_enriched_terms_fast(enrichment_res, method = "hierarchical",
                            plot_clusters_graph = TRUE, ...)
```

**Description:** Main function for clustering enriched pathway terms using optimized algorithms.

| Parameter | Type | Description |
|-----------|------|-------------|
| `enrichment_res` | data.frame | Enrichment results from pathfindR |
| `method` | character | Clustering method: "hierarchical" or "fuzzy" |
| `plot_clusters_graph` | logical | Whether to plot cluster graph |

**Returns:** List with clustered terms and visualization.

---

### cluster_graph_vis_fast

```r
cluster_graph_vis_fast(clu_obj, kappa_mat, enrichment_res, 
                       kappa_threshold = 0.35, ...)
```

**Description:** Creates fast graph visualization of term clusters.

**Returns:** igraph visualization object.

---

### create_kappa_matrix_fast

```r
create_kappa_matrix_fast(enrichment_res, use_description = FALSE, 
                         use_active_snw_genes = FALSE)
```

**Description:** Creates kappa statistic matrix for pathway term similarity using optimized algorithm.

**Returns:** Kappa matrix for term similarity.

---

### fuzzy_term_clustering_fast

```r
fuzzy_term_clustering_fast(kappa_mat, enrichment_res, 
                           kappa_threshold = 0.35, ...)
```

**Description:** Performs fuzzy clustering on pathway terms.

**Returns:** Fuzzy clustering assignments.

---

### hierarchical_term_clustering_fast

```r
hierarchical_term_clustering_fast(kappa_mat, enrichment_res, 
                                  num_clusters = NULL, ...)
```

**Description:** Performs hierarchical clustering on pathway terms.

**Returns:** Hierarchical clustering result.

---

## modified_MOVICS_functions.R

Modified MOVICS package functions for multi-omics analysis.

### compAgree2

```r
compAgree2(moic.res = NULL, subt2comp = NULL, doPlot = TRUE, ...)
```

**Description:** Compares agreement between two clustering results with visualization.

**Returns:** Agreement statistics and optional Sankey plot.

---

### compAgree_single_algorithm

```r
compAgree_single_algorithm(algorithm_name = "CS", moic.res = NULL, 
                           subt2comp = NULL, ...)
```

**Description:** Single-algorithm version of clustering agreement comparison.

**Returns:** Agreement metrics.

---

### compClinvar2

```r
compClinvar2(moic.res = NULL, var2comp = NULL, strata = NULL, 
             factorVars = NULL, ...)
```

**Description:** Compares clinical variables across subtypes.

**Returns:** Clinical comparison table with p-values.

---

### compClinvar_ordinal_single_algorithm

```r
compClinvar_ordinal_single_algorithm(algorithm_name = "CS", moic.res = NULL,
                                     var2comp = NULL, ordinal_vars = NULL, ...)
```

**Description:** Compares ordinal clinical variables using appropriate tests.

**Returns:** Ordinal comparison results.

---

### compClinvar_single_algorithm

```r
compClinvar_single_algorithm(algorithm_name = "CS", moic.res = NULL,
                             var2comp = NULL, ...)
```

**Description:** Single-algorithm clinical variable comparison.

**Returns:** Comparison table.

---

### compDrugsen_single_algorithm

```r
compDrugsen_single_algorithm(algorithm_name = "CS", moic.res = NULL,
                             norm.expr = NULL, drugs = c("Cisplatin", "Paclitaxel"), ...)
```

**Description:** Compares predicted drug sensitivity across subtypes.

**Returns:** Drug sensitivity comparison.

---

### compFGA_mod

```r
compFGA_mod(moic.res = NULL, segment = NULL, iscopynumber = FALSE, 
            ga_column = NULL, ...)
```

**Description:** Compares Fraction of Genome Altered (FGA) across subtypes.

**Returns:** FGA comparison with visualization.

---

### compFGA_optimized

```r
compFGA_optimized(moic.res = NULL, segment = NULL, iscopynumber = FALSE,
                  ga_column = NULL, ...)
```

**Description:** Optimized FGA comparison with improved performance.

**Returns:** FGA comparison results.

---

### compMut_single_algorithm

```r
compMut_single_algorithm(algorithm_name = "CS", moic.res = NULL, 
                         mut.matrix = NULL, freq.cutoff = 0.05, ...)
```

**Description:** Compares mutation frequencies across subtypes.

**Returns:** Mutation comparison with oncoprint.

---

### compMut_single_algorithm_sc

```r
compMut_single_algorithm_sc(algorithm_name = "CS", moic.res = NULL,
                            mut.matrix = NULL, ...)
```

**Description:** Single-cell adapted mutation comparison.

**Returns:** Mutation comparison for single-cell data.

---

### compSurv_ext

```r
compSurv_ext(moic.res = NULL, surv.info = NULL, convt.time = "d", ...)
```

**Description:** Extended survival comparison with multiple tests.

**Returns:** Comprehensive survival analysis.

---

### getMoHeatmap_single_algorithm

```r
getMoHeatmap_single_algorithm(algorithm_name = "CS", data = NULL,
                               is.binary = c(FALSE), ...)
```

**Description:** Generates multi-omics heatmap for single algorithm result.

**Returns:** ComplexHeatmap object.

---

### getMoHeatmap_single_algorithm2

```r
getMoHeatmap_single_algorithm2(algorithm_name = "CS", data = NULL, ...)
```

**Description:** Alternative multi-omics heatmap generation.

**Returns:** Heatmap visualization.

---

### getSilhouette_ggplot

```r
getSilhouette_ggplot(sil = NULL, clust.col = NULL, fig.path = getwd(),
                     fig.name = "silhouette_ggplot", algorithm = "", ...)
```

**Description:** Creates ggplot2-based silhouette plot.

**Returns:** ggplot2 silhouette visualization.

---

### ntp_mod

```r
ntp_mod(emat, templates, nPerm = 1000, distance = "cosine", ...)
```

**Description:** Modified Nearest Template Prediction algorithm.

**Returns:** NTP predictions with p-values.

---

### plot_pathway_heatmaps

```r
plot_pathway_heatmaps(gsea.lists, norm.expr = NULL, algorithm_name = NULL, ...)
```

**Description:** Creates pathway-level heatmaps from GSEA results.

**Returns:** Pathway heatmap visualizations.

---

### prepare_gsea_output_for_hclust

```r
prepare_gsea_output_for_hclust(gsea_output, dgea_output_name_style = "", ...)
```

**Description:** Prepares GSEA output for hierarchical clustering analysis.

**Returns:** Formatted data for pathway clustering.

---

### runDEA_mod

```r
runDEA_mod(dea.method = c("deseq2", "edger", "limma"), expr = NULL,
           moic.res = NULL, ...)
```

**Description:** Modified differential expression analysis wrapper.

**Returns:** DEA results for all pairwise comparisons.

---

### runGSEA_mod

```r
runGSEA_mod(moic.res = NULL, dea.method = c("deseq2", "edger", "limma"), ...)
```

**Description:** Modified GSEA analysis for MOVICS results.

**Returns:** GSEA enrichment results.

---

### runGSEA_mod_4.4

```r
runGSEA_mod_4.4(moic.res = NULL, dea.method = c("deseq2", "edger", "limma"), ...)
```

**Description:** GSEA modified for R 4.4 compatibility.

**Returns:** R 4.4 compatible GSEA results.

---

### runGSEA_mod_4.4_single_algorithm

```r
runGSEA_mod_4.4_single_algorithm(algorithm_name = "CS", moic.res = NULL, ...)
```

**Description:** Single-algorithm GSEA for R 4.4.

**Returns:** GSEA results with algorithm prefix.

---

### runGSVA_mod_4.4

```r
runGSVA_mod_4.4(moic.res = NULL, norm.expr = NULL, gset.gmt.path = NULL, ...)
```

**Description:** Modified GSVA analysis for R 4.4.

**Returns:** GSVA enrichment scores.

---

### runGSVA_mod_4.4_single_algorithm

```r
runGSVA_mod_4.4_single_algorithm(algorithm_name = "CS", moic.res = NULL, ...)
```

**Description:** Single-algorithm GSVA for R 4.4.

**Returns:** GSVA results with algorithm prefix.

---

### runKappa_single_algorithm

```r
runKappa_single_algorithm(algorithm_name = "CS", moic.res = NULL, ...)
```

**Description:** Calculates Cohen's Kappa for clustering agreement.

**Returns:** Kappa statistics.

---

### runMarker_mod_4.4

```r
runMarker_mod_4.4(moic.res = NULL, dea.method = c("deseq2", "edger", "limma"), ...)
```

**Description:** Marker gene identification for R 4.4.

**Returns:** Marker genes per subtype.

---

### runMarker_single_algorithm

```r
runMarker_single_algorithm(algorithm_name = "CS", moic.res = NULL,
                           dea.method = c("deseq2", "edger", "limma"), ...)
```

**Description:** Single-algorithm marker identification.

**Returns:** Marker genes with algorithm prefix.

---

### runMarker_single_algorithm_no_export

```r
runMarker_single_algorithm_no_export(algorithm_name = "CS", moic.res = NULL, ...)
```

**Description:** Marker identification without file export.

**Returns:** In-memory marker results.

---

### runNTP_mod

```r
runNTP_mod(expr = NULL, templates = NULL, scaleFlag = TRUE, 
           centerFlag = TRUE, ...)
```

**Description:** Modified Nearest Template Prediction.

**Returns:** NTP classification results.

---

### runPAM_single_algorithm

```r
runPAM_single_algorithm(algorithm_name = "CS", moic.res = NULL, 
                        norm.expr = NULL, ...)
```

**Description:** PAM (Prediction Analysis for Microarrays) for subtype prediction.

**Returns:** PAM classifier and predictions.

---

### subHeatmap_mod

```r
subHeatmap_mod(emat, res, templates, keepN = TRUE, labRow = NULL, ...)
```

**Description:** Modified subtype heatmap visualization.

**Returns:** Subtype heatmap.

---

### twoclassdeseq2_mod

```r
twoclassdeseq2_mod(moic.res = NULL, countsTable = NULL, prefix = NULL, ...)
```

**Description:** Two-class DESeq2 differential expression.

**Returns:** DESeq2 results.

---

### twoclassedger_mod

```r
twoclassedger_mod(moic.res = NULL, countsTable = NULL, prefix = NULL, ...)
```

**Description:** Two-class edgeR differential expression.

**Returns:** edgeR results.

---

### twoclasslimma_mod

```r
twoclasslimma_mod(moic.res = NULL, norm.expr = NULL, prefix = NULL, ...)
```

**Description:** Two-class limma differential expression.

**Returns:** limma results.

---

## Spectrum_functions.R

Modified Spectrum spectral clustering functions.

### CNN_kernel_mod

```r
CNN_kernel_mod(mat, NN = 3, NN2 = 7, distance = "euclidean")
```

**Description:** Modified Continuous Nearest Neighbor kernel construction.

| Parameter | Type | Description |
|-----------|------|-------------|
| `mat` | matrix | Data matrix |
| `NN` | integer | Number of nearest neighbors (default: 3) |
| `NN2` | integer | Extended neighbors (default: 7) |
| `distance` | character | Distance metric |

**Returns:** CNN kernel matrix.

---

### kernfinder_local_mod

```r
kernfinder_local_mod(data, maxk = 10, fontsize = 12, silent = FALSE, ...)
```

**Description:** Local kernel parameter finder for Spectrum.

**Returns:** Optimal local kernel parameters.

---

### kernfinder_mine_mod

```r
kernfinder_mine_mod(data, maxk = 10, fontsize = 12, silent = FALSE, ...)
```

**Description:** Mining-based kernel parameter finder.

**Returns:** Optimal kernel parameters.

---

### rbfkernel_b_mod

```r
rbfkernel_b_mod(mat, K = 3, sigma = 1, distance = "euclidean")
```

**Description:** Modified RBF kernel construction.

| Parameter | Type | Description |
|-----------|------|-------------|
| `mat` | matrix | Data matrix |
| `K` | integer | Number of neighbors (default: 3) |
| `sigma` | numeric | Kernel width (default: 1) |
| `distance` | character | Distance metric |

**Returns:** RBF kernel matrix.

---

### Spectrum_bin_and_par

```r
Spectrum_bin_and_par(data, method = "CNN", maxk = 10, ...)
```

**Description:** Main Spectrum clustering with binary data support and parallel processing.

| Parameter | Type | Description |
|-----------|------|-------------|
| `data` | list | List of data matrices |
| `method` | character | Kernel method ("CNN", "RBF") |
| `maxk` | integer | Maximum clusters to test |

**Returns:** Spectrum clustering results with cluster assignments and metrics.

---

## Alphabetical Index

| Function | Source File |
|----------|-------------|
| as.sunburstDF | custom_functions.R |
| calculate_ari_index | custom_functions.R |
| calculate_nmi_index | custom_functions.R |
| calculate_S | custom_functions.R |
| CIMLR_mod | custom_functions.R |
| CIMLR_Estimate_Number_of_Clusters_weight_mod | custom_functions.R |
| CIMLR.weight_mod | custom_functions.R |
| cluster_enriched_terms_fast | fast_pathfindR_hclust.R |
| cluster_graph_vis_fast | fast_pathfindR_hclust.R |
| CNN_kernel_mod | Spectrum_functions.R |
| coca_cc_mod | custom_functions.R |
| compAgree2 | modified_MOVICS_functions.R |
| compAgree_single_algorithm | modified_MOVICS_functions.R |
| compClinvar2 | modified_MOVICS_functions.R |
| compClinvar_ordinal_single_algorithm | modified_MOVICS_functions.R |
| compClinvar_single_algorithm | modified_MOVICS_functions.R |
| compDrugsen_single_algorithm | modified_MOVICS_functions.R |
| compFGA_mod | modified_MOVICS_functions.R |
| compFGA_optimized | modified_MOVICS_functions.R |
| compMut_single_algorithm | modified_MOVICS_functions.R |
| compMut_single_algorithm_sc | modified_MOVICS_functions.R |
| compSurv_ext | modified_MOVICS_functions.R |
| compute_matrix_similarity | custom_functions.R |
| compute_silhouette | custom_functions.R |
| concordanceNetworkNMI_ANF | custom_functions.R |
| create_annot_barchart | custom_functions.R |
| create_density_curve | custom_functions.R |
| create_density_plot_color | custom_functions.R |
| create_kappa_matrix_fast | fast_pathfindR_hclust.R |
| create_MO_heatmap | custom_functions.R |
| dist2.cimlr_mod2 | custom_functions.R |
| dist2.cimlr.weight_mod | custom_functions.R |
| divide_by_quantile_normalize | custom_functions.R |
| estimateNumberOfClustersGivenGraph_mod | custom_functions.R |
| fetch_citation | custom_functions.R |
| fetch_in_a_nutshell | custom_functions.R |
| frobenius_norm | custom_functions.R |
| fuzzy_term_clustering_fast | fast_pathfindR_hclust.R |
| getMoHeatmap_prenorm_quart | custom_functions.R |
| getMoHeatmap_single_algorithm | modified_MOVICS_functions.R |
| getMoHeatmap_single_algorithm2 | modified_MOVICS_functions.R |
| getSilhouette_ggplot | modified_MOVICS_functions.R |
| hierarchical_term_clustering_fast | fast_pathfindR_hclust.R |
| kernfinder_local_mod | Spectrum_functions.R |
| kernfinder_mine_mod | Spectrum_functions.R |
| linear_quantile_normalize | custom_functions.R |
| mds_from_original_matrix | custom_functions.R |
| MOVICS_jaccard_index | custom_functions.R |
| multiple.kernel.cimlr_mod | custom_functions.R |
| multiple.kernel.cimlr.weight_mod | custom_functions.R |
| nemo.affinity.graph_mod | custom_functions.R |
| normalize_affinity_matrix | custom_functions.R |
| ntp_mod | modified_MOVICS_functions.R |
| pca_from_original_matrix | custom_functions.R |
| pca_from_sim_matrix | custom_functions.R |
| pearson_correlation | custom_functions.R |
| plot_MOFA_cumul_var_exp | custom_functions.R |
| plot_MOFA_var_exp | custom_functions.R |
| plot_pathway_heatmaps | modified_MOVICS_functions.R |
| prepare_gsea_output_for_hclust | modified_MOVICS_functions.R |
| rankFeaturesByNMI_parallely | custom_functions.R |
| rankFeaturesByNMI_parallely_ANF | custom_functions.R |
| rbfkernel_b_mod | Spectrum_functions.R |
| reshape_Pearson_for_tests | custom_functions.R |
| runDEA_mod | modified_MOVICS_functions.R |
| runGSEA_mod | modified_MOVICS_functions.R |
| runGSEA_mod_4.4 | modified_MOVICS_functions.R |
| runGSEA_mod_4.4_single_algorithm | modified_MOVICS_functions.R |
| runGSVA_mod_4.4 | modified_MOVICS_functions.R |
| runGSVA_mod_4.4_single_algorithm | modified_MOVICS_functions.R |
| runKappa_single_algorithm | modified_MOVICS_functions.R |
| runMarker_mod_4.4 | modified_MOVICS_functions.R |
| runMarker_single_algorithm | modified_MOVICS_functions.R |
| runMarker_single_algorithm_no_export | modified_MOVICS_functions.R |
| runNTP_mod | modified_MOVICS_functions.R |
| runPAM_single_algorithm | modified_MOVICS_functions.R |
| Spectrum_bin_and_par | Spectrum_functions.R |
| standardize_rows | custom_functions.R |
| subHeatmap_mod | modified_MOVICS_functions.R |
| TCGAanalyze_survival_custom3 | custom_functions.R |
| twoclassdeseq2_mod | modified_MOVICS_functions.R |
| twoclassedger_mod | modified_MOVICS_functions.R |
| twoclasslimma_mod | modified_MOVICS_functions.R |

---

*Documentation auto-generated for MultiOmicsSurvey project.*
