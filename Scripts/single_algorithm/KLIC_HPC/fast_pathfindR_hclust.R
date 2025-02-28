# This script makes hierarchical clustering faster when using pathfindR

library(pathfindR)
library(fastcluster)

create_kappa_matrix_fast <- function(enrichment_res, use_description = FALSE, use_active_snw_genes = FALSE) {
  ### Argument checks
  if (!is.logical(use_description)) {
    stop("`use_description` should be TRUE or FALSE")
  }
  
  if (!is.logical(use_active_snw_genes)) {
    stop("`use_active_snw_genes` should be TRUE or FALSE")
  }
  
  if (!is.data.frame(enrichment_res)) {
    stop("`enrichment_res` should be a data frame of enrichment results")
  }
  if (nrow(enrichment_res) < 2) {
    stop("`enrichment_res` should contain at least 2 rows")
  }
  
  nec_cols <- c("Down_regulated", "Up_regulated")
  if (use_description) {
    nec_cols <- c("Term_Description", nec_cols)
  } else {
    nec_cols <- c("ID", nec_cols)
  }
  if (use_active_snw_genes) {
    nec_cols <- c(nec_cols, "non_Signif_Snw_Genes")
  }
  
  if (!all(nec_cols %in% colnames(enrichment_res))) {
    stop("`enrichment_res` should contain all of ", paste(dQuote(nec_cols), collapse = ", "))
  }
  
  ### Initial steps Column to use for gene set names
  chosen_id <- ifelse(use_description, which(colnames(enrichment_res) == "Term_Description"),
                      which(colnames(enrichment_res) == "ID"))
  
  # list of genes
  down_idx <- which(colnames(enrichment_res) == "Down_regulated")
  up_idx <- which(colnames(enrichment_res) == "Up_regulated")
  
  genes_lists <- apply(enrichment_res, 1, function(x) {
    base::toupper(c(unlist(strsplit(as.character(x[up_idx]), ", ")), unlist(strsplit(as.character(x[down_idx]),
                                                                                     ", "))))
  })
  
  if (use_active_snw_genes) {
    active_idx <- which(colnames(enrichment_res) == "non_Signif_Snw_Genes")
    
    genes_lists <- mapply(function(x, y) {
      c(x, unlist(strsplit(as.character(y), ", ")))
    }, genes_lists, enrichment_res[, active_idx])
  }
  
  # Exclude zero-length gene sets
  excluded_idx <- which(vapply(genes_lists, length, 1) == 0)
  if (length(excluded_idx) != 0) {
    genes_lists <- genes_lists[-excluded_idx]
    enrichment_res <- enrichment_res[-excluded_idx, ]
  }
  
  ### Create Kappa Matrix
  all_genes <- unique(unlist(genes_lists, use.names = FALSE))
  N <- nrow(enrichment_res)
  term_names <- enrichment_res[, chosen_id]
  
  kappa_mat <- matrix(0, nrow = N, ncol = N, dimnames = list(term_names, term_names))
  diag(kappa_mat) <- 1
  
  total <- length(all_genes)
  for (i in 1:(N - 1)) {
    for (j in (i + 1):N) {
      genes_i <- genes_lists[[i]]
      genes_j <- genes_lists[[j]]
      
      both <- length(intersect(genes_i, genes_j))
      term_i <- length(base::setdiff(genes_i, genes_j))
      term_j <- length(base::setdiff(genes_j, genes_i))
      no_terms <- total - sum(both, term_i, term_j)
      
      observed <- (both + no_terms)/total
      chance <- (both + term_i) * (both + term_j)
      chance <- chance + (term_j + no_terms) * (term_i + no_terms)
      chance <- chance/total^2
      kappa_mat[j, i] <- kappa_mat[i, j] <- (observed - chance)/(1 - chance)
    }
  }
  
  return(kappa_mat)
}

hierarchical_term_clustering_fast <- function(kappa_mat, enrichment_res, num_clusters = NULL,
                                         use_description = FALSE, clu_method = "average", plot_hmap = FALSE, plot_dend = TRUE) {
  ### Set ID/Name index
  chosen_id <- ifelse(use_description, which(colnames(enrichment_res) == "Term_Description"),
                      which(colnames(enrichment_res) == "ID"))
  
  ### Argument checks
  if (!isSymmetric.matrix(kappa_mat)) {
    stop("`kappa_mat` should be a symmetric matrix")
  }
  
  if (!all(colnames(kappa_mat) %in% enrichment_res[, chosen_id])) {
    stop("All terms in `kappa_mat` should be present in `enrichment_res`")
  }
  
  if (!is.logical(plot_hmap)) {
    stop("`plot_hmap` should be TRUE or FALSE")
  }
  
  if (!is.logical(plot_dend)) {
    stop("`plot_dend` should be TRUE or FALSE")
  }
  
  ### Add excluded (zero-length) genes
  kappa_mat2 <- kappa_mat
  cond <- !enrichment_res[, chosen_id] %in% rownames(kappa_mat2)
  outliers <- enrichment_res[cond, chosen_id]
  outliers_mat <- matrix(-1, nrow = nrow(kappa_mat2), ncol = length(outliers),
                         dimnames = list(rownames(kappa_mat2), outliers))
  kappa_mat2 <- cbind(kappa_mat2, outliers_mat)
  outliers_mat <- matrix(-1, nrow = length(outliers), ncol = ncol(kappa_mat2),
                         dimnames = list(outliers, colnames(kappa_mat2)))
  kappa_mat2 <- rbind(kappa_mat2, outliers_mat)
  
  ### Perform hierarchical clustering
  clu <- fastcluster::hclust(stats::as.dist(1 - kappa_mat2), method = clu_method)
  
  if (plot_hmap) {
    stats::heatmap(kappa_mat2, distfun = function(x) stats::as.dist(1 - x), hclustfun = function(x) stats::hclust(x,
                                                                                                                  method = clu_method))
  }
  
  ### Choose optimal k (if not specified)
  if (is.null(num_clusters)) {
    kmax <- max(nrow(kappa_mat2)%/%2, 2)
    
    # sequence of k (number of clusters) to try
    if (kmax <= 20) {
      kseq <- 2:kmax
    } else if (kmax <= 100) {
      kseq <- c(2:19, seq(20, kmax%/%10 * 10, 10))
    } else {
      kseq <- c(2:19, seq(20, 99, 10), seq(100, kmax%/%50 * 50, 50))
    }
    
    # calculate average silhouette width per k in sequence
    avg_sils <- c()
    for (k in kseq) {
      avg_sils <- c(avg_sils, fpc::cluster.stats(stats::as.dist(1 - kappa_mat2),
                                                 stats::cutree(clu, k = k), silhouette = TRUE)$avg.silwidth)
    }
    
    k_opt <- kseq[which.max(avg_sils)]
    
    message(paste("The maximum average silhouette width was", round(max(avg_sils),
                                                                    2), "for k =", k_opt, "\n\n"))
  } else {
    k_opt <- num_clusters
  }
  
  
  if (plot_dend) {
    graphics::plot(clu)
    stats::rect.hclust(clu, k = k_opt)
  }
  
  clusters <- stats::cutree(clu, k = k_opt)
  
  return(clusters)
}

fuzzy_term_clustering_fast <- function(kappa_mat, enrichment_res, kappa_threshold = 0.35,
                                  use_description = FALSE) {
  ### Set ID/Name index
  chosen_id <- ifelse(use_description, which(colnames(enrichment_res) == "Term_Description"),
                      which(colnames(enrichment_res) == "ID"))
  
  ### Argument checks
  if (!isSymmetric.matrix(kappa_mat)) {
    stop("`kappa_mat` should be a symmetric matrix")
  }
  
  if (!all(colnames(kappa_mat) %in% enrichment_res[, chosen_id])) {
    stop("All terms in `kappa_mat` should be present in `enrichment_res`")
  }
  
  if (!is.numeric(kappa_threshold)) {
    stop("`kappa_threshold` should be numeric")
  }
  
  if (kappa_threshold > 1) {
    stop("`kappa_threshold` should be at most 1 as kappa statistic is always <= 1")
  }
  
  ### Find Qualified Seeds
  qualified_seeds <- list()
  j <- 1
  for (i in base::seq_len(nrow(kappa_mat))) {
    current_term <- rownames(kappa_mat)[i]
    current_term_kappa <- kappa_mat[i, ]
    
    
    init_membership_cond <- current_term_kappa >= kappa_threshold
    if (sum(init_membership_cond) > 3) {
      related_terms <- names(current_term_kappa)[init_membership_cond]
      terms <- c(current_term, related_terms)
      related_kappa <- kappa_mat[rownames(kappa_mat) %in% terms, colnames(kappa_mat) %in%
                                   terms]
      diag(related_kappa) <- 0
      tight_relationship_cond <- sum(related_kappa >= kappa_threshold)/(nrow(related_kappa)^2) >=
        0.5
      
      if (tight_relationship_cond) {
        qualified_seeds[[j]] <- related_terms
        names(qualified_seeds)[j] <- current_term
        j <- j + 1
      }
    }
  }
  
  ### Fuzzy Clustering
  clusters <- unique(qualified_seeds)
  i <- 1
  j <- i + 1
  while (i < length(clusters)) {
    common_terms <- intersect(clusters[[i]], clusters[[j]])
    all_terms <- union(clusters[[i]], clusters[[j]])
    
    if (length(common_terms)/length(all_terms) > 0.5 & i != j) {
      clusters[[i]] <- all_terms
      clusters[[j]] <- NULL
      i <- 1
      j <- i + 1
    } else if (j < length(clusters)) {
      j <- j + 1
    } else {
      i <- i + 1
      j <- 1
    }
  }
  
  ### Find Outliers
  cond <- !enrichment_res[, chosen_id] %in% c(names(clusters), unlist(clusters))
  outliers <- enrichment_res[cond, chosen_id]
  for (outlier in outliers) {
    clusters[[outlier]] <- outlier
  }
  ### Return Cluster Matrix
  names(clusters) <- base::seq_len(length(clusters))
  
  cluster_mat <- matrix(FALSE, nrow = nrow(enrichment_res), ncol = length(clusters),
                        dimnames = list(enrichment_res[, chosen_id], names(clusters)))
  for (clu in names(clusters)) {
    clu_terms <- clusters[[clu]]
    cluster_mat[clu_terms, clu] <- TRUE
  }
  
  return(cluster_mat)
}

cluster_graph_vis_fast <- function(clu_obj, kappa_mat, enrichment_res, kappa_threshold = 0.35,
                              use_description = FALSE, vertex.label.cex = 0.7, vertex.size.scaling = 2.5) {
  ### Set ID/Name index
  chosen_id <- ifelse(use_description, which(colnames(enrichment_res) == "Term_Description"),
                      which(colnames(enrichment_res) == "ID"))
  
  ### For coloring nodes
  all_cols <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33",
                "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3",
                "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", "#8DD3C7", "#FFFFB3", "#BEBADA",
                "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD",
                "#CCEBC5", "#FFED6F", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
                "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928")
  
  if (is.matrix(clu_obj)) {
    ### Argument checks
    if (!all(rownames(clu_obj) %in% colnames(kappa_mat))) {
      stop("Not all terms in `clu_obj` present in `kappa_mat`!")
    }
    
    ### Prep data Remove weak links
    kappa_mat2 <- kappa_mat
    diag(kappa_mat2) <- 0
    kappa_mat2 <- ifelse(kappa_mat2 < kappa_threshold, 0, kappa_mat2)
    
    # Add missing terms
    missing <- rownames(clu_obj)[!rownames(clu_obj) %in% colnames(kappa_mat2)]
    missing_mat <- matrix(0, nrow = nrow(kappa_mat2), ncol = length(missing),
                          dimnames = list(rownames(kappa_mat2), missing))
    kappa_mat2 <- cbind(kappa_mat2, missing_mat)
    missing <- rownames(clu_obj)[!rownames(clu_obj) %in% rownames(kappa_mat2)]
    missing_mat <- matrix(0, nrow = length(missing), ncol = ncol(kappa_mat2),
                          dimnames = list(missing, colnames(kappa_mat2)))
    kappa_mat2 <- rbind(kappa_mat2, missing_mat)
    
    ### Create Graph, Set Color, Size and Percentages
    values <- apply(clu_obj, 1, function(x) which(x))
    percs <- list()
    for (i in base::seq_len(length(values))) {
      percs[[i]] <- rep(1/length(values[[i]]), length(values[[i]]))
    }
    
    g <- igraph::graph_from_adjacency_matrix(kappa_mat2, weighted = TRUE)
    
    if (length(all_cols) < max(as.integer(colnames(clu_obj)))) {
      num_extra <- max(as.integer(colnames(clu_obj))) - length(all_cols)
      extra_colors <- grDevices::rainbow(num_extra)
      all_cols <- c(all_cols, extra_colors)
    }
    
    # Node shapes are either circle (single cluster) or pie (multiple
    # clusters)
    igraph::V(g)$shape <- ifelse(vapply(percs, length, 1) > 1, "pie", "circle")
    
    # Node colors are cluster memberships
    cols <- lapply(values, function(x) all_cols[x])
    igraph::V(g)$color <- vapply(cols, function(x) x[1], "")
    
    # Node sizes are -log(lowest_p)
    p_idx <- match(names(igraph::V(g)), enrichment_res[, chosen_id])
    transformed_p <- -log10(enrichment_res$lowest_p[p_idx])
    igraph::V(g)$size <- transformed_p * vertex.size.scaling
    
    ### Plot Graph
    igraph::plot.igraph(g, vertex.pie = percs, vertex.pie.color = cols, layout = igraph::layout_nicely(g),
                        edge.curved = FALSE, vertex.label.dist = 0, vertex.label.color = "black",
                        asp = 1, vertex.label.cex = vertex.label.cex, edge.width = igraph::E(g)$weight,
                        edge.arrow.mode = 0)
  } else if (is.integer(clu_obj)) {
    ### Argument checks
    if (!all(names(clu_obj) %in% colnames(kappa_mat))) {
      stop("Not all terms in `clu_obj` present in `kappa_mat`!")
    }
    
    ### Prep data Remove weak links
    kappa_mat2 <- kappa_mat
    diag(kappa_mat2) <- 0
    kappa_mat2 <- ifelse(kappa_mat2 > kappa_threshold, kappa_mat2, 0)
    
    # Add missing terms
    missing <- names(clu_obj)[!names(clu_obj) %in% colnames(kappa_mat2)]
    missing_mat <- matrix(0, nrow = nrow(kappa_mat2), ncol = length(missing),
                          dimnames = list(rownames(kappa_mat2), missing))
    kappa_mat2 <- cbind(kappa_mat2, missing_mat)
    missing <- names(clu_obj)[!names(clu_obj) %in% rownames(kappa_mat2)]
    missing_mat <- matrix(0, nrow = length(missing), ncol = ncol(kappa_mat2),
                          dimnames = list(missing, colnames(kappa_mat2)))
    kappa_mat2 <- rbind(kappa_mat2, missing_mat)
    
    ### Create Graph, Set Colors and Sizes
    g <- igraph::graph_from_adjacency_matrix(kappa_mat2, weighted = TRUE)
    
    igraph::V(g)$Clu <- clu_obj[match(igraph::V(g)$name, names(clu_obj))]
    
    if (length(all_cols) < max(as.integer(igraph::V(g)$Clu))) {
      num_extra <- max(clu_obj) - length(all_cols)
      extra_colors <- grDevices::rainbow(num_extra)
      all_cols <- c(all_cols, extra_colors)
    }
    
    # Node colors are cluster memberships
    igraph::V(g)$color <- all_cols[as.integer(igraph::V(g)$Clu)]
    
    # Node sizes are -log(lowest_p)
    p_idx <- match(names(igraph::V(g)), enrichment_res[, chosen_id])
    transformed_p <- -log10(enrichment_res$lowest_p[p_idx])
    igraph::V(g)$size <- transformed_p * vertex.size.scaling
    
    ### Plot graph
    igraph::plot.igraph(g, layout = igraph::layout_nicely(g), edge.curved = FALSE,
                        vertex.label.dist = 0, vertex.label.color = "black", asp = 0, vertex.label.cex = vertex.label.cex,
                        edge.width = igraph::E(g)$weight, edge.arrow.mode = 0)
  } else {
    stop("Invalid class for `clu_obj`!")
  }
}

cluster_enriched_terms_fast <- function(enrichment_res, method = "hierarchical", plot_clusters_graph = TRUE,
                                        use_description = FALSE, use_active_snw_genes = FALSE, ...) {
  ### Argument Checks
  if (!method %in% c("hierarchical", "fuzzy")) {
    stop("The clustering `method` must either be \"hierarchical\" or \"fuzzy\"")
  }
  
  if (!is.logical(plot_clusters_graph)) {
    stop("`plot_clusters_graph` must be logical!")
  }
  
  ### Set ID/Name index for logging
  chosen_id <- if ("Term_Description" %in% colnames(enrichment_res) && use_description) {
    which(colnames(enrichment_res) == "Term_Description")
  } else if ("ID" %in% colnames(enrichment_res)) {
    which(colnames(enrichment_res) == "ID")
  } else {
    return("hclust impossible")  # Return "hclust impossible" if neither column is found
  }
  
  if (is.null(chosen_id) || length(chosen_id) == 0) {
    return("hclust impossible")  # Return "hclust impossible" if chosen_id is invalid
  }
  
  ### Handle cases with fewer than 2 rows
  if (nrow(enrichment_res) < 2) {
    return("hclust impossible")  # Return "hclust impossible" if there are fewer than 2 rows
  }
  
  ### Create Kappa Matrix
  kappa_mat <- tryCatch({
    create_kappa_matrix_fast(enrichment_res = enrichment_res, use_description = use_description,
                             use_active_snw_genes = use_active_snw_genes)
  }, error = function(e) {
    return("hclust impossible")  # Return "hclust impossible" if kappa matrix creation fails
  })
  
  kappa_mat[is.na(kappa_mat)] <- 0
  
  ### Cluster Terms
  clu_obj <- NULL
  if (method == "hierarchical") {
    clu_obj <- tryCatch({
      hierarchical_term_clustering_fast(kappa_mat = kappa_mat, enrichment_res = enrichment_res, use_description = use_description, ...)
    }, error = function(e) {
      return("hclust impossible")  # Return "hclust impossible" if hierarchical clustering fails
    })
  } else if (method == "fuzzy") {
    clu_obj <- tryCatch({
      fuzzy_term_clustering_fast(kappa_mat = kappa_mat, enrichment_res = enrichment_res, use_description = use_description, ...)
    }, error = function(e) {
      return("hclust impossible")  # Return "hclust impossible" if fuzzy clustering fails
    })
  }
  
  if (is.character(clu_obj) && clu_obj == "hclust impossible") {
    return("hclust impossible")  # If clustering failed, return "hclust impossible"
  }
  
  ### Graph Visualization of Clusters
  if (plot_clusters_graph) {
    tryCatch({
      cluster_graph_vis_fast(clu_obj = clu_obj, kappa_mat = kappa_mat, enrichment_res = enrichment_res, use_description = use_description, ...)
    }, error = function(e) {
      message("Error in plotting cluster graph: ", e$message)
    })
  }
  
  ### Returned Data Frame with Cluster Information
  clustered_df <- enrichment_res
  
  ### Assign Clusters and Representatives for hierarchical method
  if (method == "hierarchical") {
    clu_idx <- match(clustered_df[, chosen_id], names(clu_obj))
    
    if (any(is.na(clu_idx))) {
      warning("Some terms in clustered_df could not be matched to clu_obj. They will be excluded.")
      clustered_df <- clustered_df[!is.na(clu_idx), ]
      clu_idx <- clu_idx[!is.na(clu_idx)]
    }
    
    clustered_df$Cluster <- clu_obj[clu_idx]
    clustered_df <- clustered_df[order(clustered_df$Cluster, clustered_df$lowest_p, decreasing = FALSE), ]
    
    tmp <- tapply(clustered_df[, chosen_id], clustered_df$Cluster, function(x) x[1])
    stat_cond <- clustered_df[, chosen_id] %in% tmp
    clustered_df$Status <- ifelse(stat_cond, "Representative", "Member")
  }
  
  return(list(clustered_df = clustered_df))  # Removed `skipped` slot
}
