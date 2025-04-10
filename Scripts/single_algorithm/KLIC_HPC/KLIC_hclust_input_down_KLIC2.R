library(pathfindR)
library(fastcluster)
source("fast_pathfindR_hclust.R")

RNGversion("4.2.2")
set.seed(123)

input = readRDS("hclust_input_down_KLIC2.rds")

t1 = Sys.time()
result = cluster_enriched_terms_fast(input, method = "hierarchical", plot_clusters_graph = FALSE,
                                     use_description = FALSE, use_active_snw_genes = FALSE)
dt = Sys.time() - t1

export = list(clustering = result, dt = dt)

# Exports
saveRDS(export, "hclust_output_down_KLIC2.rds")
writeLines(capture.output(sessionInfo()), "sessionInfo_hclust_input_down_KLIC2.txt")
