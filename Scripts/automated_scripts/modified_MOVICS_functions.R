# runMarker #####

# Just changing the CS labels when focusing on a single algorithm
runMarker_single_algorithm = function (algorithm_name = "CS", moic.res = NULL, dea.method = c("deseq2", "edger", 
                                          "limma"), prefix = NULL, dat.path = getwd(), res.path = getwd(), 
          p.cutoff = 0.05, p.adj.cutoff = 0.05, dirct = "up", n.marker = 200, 
          doplot = TRUE, norm.expr = NULL, annCol = NULL, annColors = NULL, 
          clust.col = c("#2EC4B6", "#E71D36", "#FF9F1C", "#BDD5EA", 
                        "#FFA5AB", "#011627", "#023E8A", "#9D4EDD"), halfwidth = 3, 
          centerFlag = TRUE, scaleFlag = TRUE, show_rownames = FALSE, 
          show_colnames = FALSE, color = c("#5bc0eb", "black", "#ECE700"), 
          fig.path = getwd(), fig.name = NULL, width = 8, height = 8, 
          ...) 
{
  n.moic <- length(unique(moic.res$clust.res$clust))
  mo.method <- moic.res$mo.method
  DEpattern <- paste(mo.method, "_", ifelse(is.null(prefix), 
                                            "", paste0(prefix, "_")), dea.method, ".*._vs_Others.txt$", 
                     sep = "")
  DEfiles <- dir(dat.path, pattern = DEpattern)
  if (length(DEfiles) == 0) {
    stop("no DEfiles!")
  }
  if (length(DEfiles) != n.moic) {
    stop("not all the multi-omics clusters have DEfile!")
  }
  if (!is.element(dirct, c("up", "down"))) {
    stop("dirct type error! Allowed value contains c('up', 'down').")
  }
  if (dirct == "up") {
    outlabel <- "unique_upexpr_marker.txt"
  }
  if (dirct == "down") {
    outlabel <- "unique_downexpr_marker.txt"
  }
  genelist <- c()
  for (filek in DEfiles) {
    DEres <- read.table(file.path(dat.path, filek), header = TRUE, 
                        row.names = NULL, sep = "\t", quote = "", stringsAsFactors = FALSE)
    DEres <- DEres[!duplicated(DEres[, 1]), ]
    DEres <- DEres[!is.na(DEres[, 1]), ]
    rownames(DEres) <- DEres[, 1]
    DEres <- DEres[, -1]
    if (dirct == "up") {
      genelist <- c(genelist, rownames(DEres[!is.na(DEres$padj) & 
                                               DEres$pvalue < p.cutoff & DEres$padj < p.adj.cutoff & 
                                               !is.na(DEres$log2fc) & DEres$log2fc > 0, ]))
    }
    if (dirct == "down") {
      genelist <- c(genelist, rownames(DEres[!is.na(DEres$padj) & 
                                               DEres$pvalue < p.cutoff & DEres$padj < p.adj.cutoff & 
                                               !is.na(DEres$log2fc) & DEres$log2fc < 0, ]))
    }
  }
  unqlist <- setdiff(genelist, genelist[duplicated(genelist)])
  marker <- list()
  for (filek in DEfiles) {
    DEres <- read.table(file.path(dat.path, filek), header = TRUE, 
                        row.names = NULL, sep = "\t", quote = "", stringsAsFactors = FALSE)
    DEres <- DEres[!duplicated(DEres[, 1]), ]
    DEres <- DEres[!is.na(DEres[, 1]), ]
    rownames(DEres) <- DEres[, 1]
    DEres <- DEres[, -1]
    if (dirct == "up") {
      outk <- intersect(unqlist, rownames(DEres[!is.na(DEres$padj) & 
                                                  DEres$pvalue < p.cutoff & DEres$padj < p.adj.cutoff & 
                                                  !is.na(DEres$log2fc) & DEres$log2fc > 0, ]))
      outk <- DEres[outk, ]
      outk <- outk[order(outk$log2fc, decreasing = TRUE), 
      ]
      if (nrow(outk) > n.marker) {
        marker[[filek]] <- outk[1:n.marker, ]
      }
      else {
        marker[[filek]] <- outk
      }
      marker$dirct <- "up"
    }
    if (dirct == "down") {
      outk <- intersect(unqlist, rownames(DEres[!is.na(DEres$padj) & 
                                                  DEres$pvalue < p.cutoff & DEres$padj < p.adj.cutoff & 
                                                  !is.na(DEres$log2fc) & DEres$log2fc < 0, ]))
      outk <- DEres[outk, ]
      outk <- outk[order(outk$log2fc, decreasing = FALSE), 
      ]
      if (nrow(outk) > n.marker) {
        marker[[filek]] <- outk[1:n.marker, ]
      }
      else {
        marker[[filek]] <- outk
      }
      marker$dirct <- "down"
    }
    write.table(outk, file = file.path(res.path, paste(gsub("_vs_Others.txt", 
                                                            "", filek, fixed = TRUE), outlabel, sep = "_")), 
                row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE)
  }
  templates <- NULL
  for (filek in DEfiles) {
    tmp <- data.frame(probe = rownames(marker[[filek]]), 
                      class = sub("_vs_Others.txt", "", sub(".*.result.", 
                                                            "", filek)), dirct = marker$dirct, stringsAsFactors = FALSE)
    templates <- rbind.data.frame(templates, tmp, stringsAsFactors = FALSE)
  }
  write.table(templates, file = file.path(res.path, paste0(mo.method, 
                                                           "_", dea.method, "_", dirct, "regulated_marker_templates.txt")), 
              row.names = FALSE, sep = "\t", quote = FALSE)
  if (doplot) {
    if (is.null(norm.expr)) {
      stop("please provide a matrix or data.frame of normalized expression data with rows for genes and columns for samples; FPKM or TPM without log2 transformation is recommended.")
    }
    comsam <- intersect(moic.res$clust.res$samID, colnames(norm.expr))
    if (length(comsam) == nrow(moic.res$clust.res)) {
      message("--all samples matched.")
    }
    else {
      message(paste0("--", (nrow(moic.res$clust.res) - 
                              length(comsam)), " samples mismatched from current subtypes."))
    }
    moic.res$clust.res <- moic.res$clust.res[comsam, , drop = FALSE]
    norm.expr <- norm.expr[, comsam]
    if (is.null(fig.name)) {
      outFig <- paste0("markerheatmap_using_", dirct, 
                       "regulated_genes.pdf")
    }
    else {
      outFig <- paste0(fig.name, "_using_", dirct, "regulated_genes.pdf")
    }
    sam.order <- moic.res$clust.res[order(moic.res$clust.res$clust, 
                                          decreasing = FALSE), "samID"]
    colvec <- clust.col[1:n.moic]
    names(colvec) <- paste0(algorithm_name, 1:n.moic)
    if (!is.null(annCol) & !is.null(annColors)) {
      annCol <- annCol[sam.order, , drop = FALSE]
      annCol$Subtype <- paste0(algorithm_name, moic.res$clust.res[sam.order, 
                                                        "clust"])
      annColors[["Subtype"]] <- colvec
    }
    else {
      annCol <- data.frame(Subtype = paste0(algorithm_name, moic.res$clust.res[sam.order, 
                                                                     "clust"]), row.names = sam.order, stringsAsFactors = FALSE)
      annColors <- list(Subtype = colvec)
    }
    if (max(norm.expr) < 25 | (max(norm.expr) >= 25 & min(norm.expr) < 
                               0)) {
      message("--expression profile seems to have been standardised (z-score or log transformation), no more action will be performed.")
      gset <- norm.expr
    }
    if (max(norm.expr) >= 25 & min(norm.expr) >= 0) {
      message("--log2 transformation done for expression data.")
      gset <- log2(norm.expr + 1)
    }
    standarize.fun <- function(indata = NULL, halfwidth = NULL, 
                               centerFlag = TRUE, scaleFlag = TRUE) {
      outdata = t(scale(t(indata), center = centerFlag, 
                        scale = scaleFlag))
      if (!is.null(halfwidth)) {
        outdata[outdata > halfwidth] = halfwidth
        outdata[outdata < (-halfwidth)] = -halfwidth
      }
      return(outdata)
    }
    plotdata <- standarize.fun(gset[intersect(templates$probe, 
                                              rownames(gset)), sam.order], halfwidth = halfwidth, 
                               centerFlag = centerFlag, scaleFlag = scaleFlag)
    if (!is.null(annCol) & !is.null(annColors)) {
      for (i in names(annColors)) {
        if (is.function(annColors[[i]])) {
          annColors[[i]] <- annColors[[i]](pretty(range(annCol[, 
                                                               i]), n = 64))
        }
      }
    }
    hm <- ComplexHeatmap::pheatmap(mat = plotdata, border_color = NA, 
                                   cluster_cols = FALSE, cluster_rows = FALSE, annotation_col = annCol, 
                                   annotation_colors = annColors, legend_breaks = pretty(c(-halfwidth, 
                                                                                           halfwidth)), legend_labels = pretty(c(-halfwidth, 
                                                                                                                                 halfwidth)), show_rownames = show_rownames, 
                                   show_colnames = show_colnames, treeheight_col = 0, 
                                   treeheight_row = 0, color = (grDevices::colorRampPalette(color))(64), 
                                   ...)
    pdf(file.path(fig.path, outFig), width = width, height = height)
    draw(hm, annotation_legend_side = "left", heatmap_legend_side = "left")
    invisible(dev.off())
    draw(hm, annotation_legend_side = "left", heatmap_legend_side = "left")
    return(list(unqlist = unqlist, templates = templates, 
                dirct = dirct, heatmap = hm))
  }
  else {
    return(list(unqlist = unqlist, templates = templates, 
                dirct = dirct))
  }
}


# runGSEA #####

# Modifications:
# add a name argument for colorbars
runGSEA_mod <- function (moic.res = NULL, dea.method = c("deseq2", "edger", 
                                          "limma"), norm.expr = NULL, prefix = NULL, dat.path = getwd(), 
          res.path = getwd(), dirct = "up", n.path = 10, msigdb.path = NULL, 
          nPerm = 1000, minGSSize = 10, maxGSSize = 500, p.cutoff = 0.05, 
          p.adj.cutoff = 0.05, gsva.method = "gsva", norm.method = "mean", 
          clust.col = c("#2EC4B6", "#E71D36", "#FF9F1C", "#BDD5EA", 
                        "#FFA5AB", "#011627", "#023E8A", "#9D4EDD"), color = NULL, 
          fig.name = NULL, fig.path = getwd(), width = 15, height = 10, name = NULL) 
{
  comsam <- intersect(moic.res$clust.res$samID, colnames(norm.expr))
  if (length(comsam) == nrow(moic.res$clust.res)) {
    message("--all samples matched.")
  }
  else {
    message(paste0("--", (nrow(moic.res$clust.res) - length(comsam)), 
                   " samples mismatched from current subtypes."))
  }
  moic.res$clust.res <- moic.res$clust.res[comsam, , drop = FALSE]
  norm.expr <- norm.expr[, comsam]
  n.moic <- length(unique(moic.res$clust.res$clust))
  mo.method <- moic.res$mo.method
  DEpattern <- paste(mo.method, "_", ifelse(is.null(prefix), 
                                            "", paste0(prefix, "_")), dea.method, ".*._vs_Others.txt$", 
                     sep = "")
  DEfiles <- dir(dat.path, pattern = DEpattern)
  if (length(DEfiles) == 0) {
    stop("no DEfiles!")
  }
  if (length(DEfiles) != n.moic) {
    stop("not all multi-omics clusters have DEfile!")
  }
  if (!is.element(dirct, c("up", "down"))) {
    stop("dirct type error! Allowed value contains c('up', 'down').")
  }
  if (!is.element(gsva.method, c("gsva", "ssgsea", "zscore", 
                                 "plage"))) {
    stop("GSVA method error! Allowed value contains c('gsva', 'ssgsea', 'zscore', 'plage').")
  }
  if (dirct == "up") {
    outlabel <- "unique_upexpr_pathway.txt"
  }
  if (dirct == "down") {
    outlabel <- "unique_downexpr_pathway.txt"
  }
  gsea.list <- list()
  gseaidList <- c()
  for (filek in DEfiles) {
    DEres <- read.table(file.path(dat.path, filek), header = TRUE, 
                        row.names = NULL, sep = "\t", quote = "", stringsAsFactors = FALSE)
    DEres <- DEres[!duplicated(DEres[, 1]), ]
    DEres <- DEres[!is.na(DEres[, 1]), ]
    rownames(DEres) <- DEres[, 1]
    DEres <- DEres[, -1]
    geneList <- DEres$log2fc
    names(geneList) <- rownames(DEres)
    geneList <- sort(geneList, decreasing = TRUE)
    msigdb <- try(clusterProfiler::read.gmt(msigdb.path), 
                  silent = TRUE)
    if (class(msigdb) == "try-error") {
      stop("please provide correct ABSOLUTE PATH for MSigDB file.")
    }
    moic.lab <- paste0("CS", gsub("_vs_Others.txt", "", 
                                  sub(".*CS", "", filek)))
    gsea.list[[moic.lab]] <- suppressWarnings(clusterProfiler::GSEA(geneList = geneList, 
                                                                    TERM2GENE = msigdb, nPerm = nPerm, minGSSize = minGSSize, 
                                                                    maxGSSize = maxGSSize, seed = TRUE, verbose = FALSE, 
                                                                    pvalueCutoff = 1))
    gsea.dat <- as.data.frame(gsea.list[[moic.lab]])
    write.table(gsea.dat[, setdiff(colnames(gsea.dat), "ID")], 
                file = file.path(res.path, paste(gsub("_vs_Others.txt", 
                                                      "", filek, fixed = TRUE), "gsea_all_results.txt", 
                                                 sep = "_")), sep = "\t", row.names = FALSE, 
                col.names = TRUE, quote = FALSE)
    if (dirct == "up") {
      gseaidList <- c(gseaidList, rownames(gsea.dat[which(gsea.dat$NES > 
                                                            0 & gsea.dat$pvalue < p.cutoff & gsea.dat$p.adjust < 
                                                            p.adj.cutoff), ]))
    }
    if (dirct == "down") {
      gseaidList <- c(gseaidList, rownames(gsea.dat[which(gsea.dat$NES < 
                                                            0 & gsea.dat$pvalue < p.cutoff & gsea.dat$p.adjust < 
                                                            p.adj.cutoff), ]))
    }
    unqlist <- setdiff(gseaidList, gseaidList[duplicated(gseaidList)])
  }
  message("GSEA done...")
  GSEApattern <- paste(mo.method, "_", ifelse(is.null(prefix), 
                                              "", paste0(prefix, "_")), dea.method, ".*._gsea_all_results.txt$", 
                       sep = "")
  GSEAfiles <- dir(res.path, pattern = GSEApattern)
  pathway <- pathcore <- list()
  pathnum <- c()
  for (filek in GSEAfiles) {
    GSEAres <- read.table(file.path(res.path, filek), header = TRUE, 
                          row.names = 1, sep = "\t", quote = "", stringsAsFactors = FALSE)
    if (dirct == "up") {
      outk <- intersect(unqlist, rownames(GSEAres[which(GSEAres$NES > 
                                                          0 & GSEAres$pvalue < p.cutoff & GSEAres$p.adjust < 
                                                          p.adj.cutoff), ]))
      outk <- GSEAres[outk, ]
      outk <- outk[order(outk$NES, decreasing = TRUE), 
      ]
      if (nrow(outk) > n.path) {
        pathway[[filek]] <- outk[1:n.path, ]
      }
      else {
        pathway[[filek]] <- outk
      }
      pathnum <- c(pathnum, nrow(pathway[[filek]]))
      pathway$dirct <- "up"
      for (i in rownames(pathway[[filek]])) {
        pathcore[[i]] <- msigdb[which(msigdb[, 1] %in% 
                                        i), "gene"]
      }
    }
    if (dirct == "down") {
      outk <- intersect(unqlist, rownames(GSEAres[which(GSEAres$NES < 
                                                          0 & GSEAres$pvalue < p.cutoff & GSEAres$p.adjust < 
                                                          p.adj.cutoff), ]))
      outk <- GSEAres[outk, ]
      outk <- outk[order(outk$NES, decreasing = FALSE), 
      ]
      if (nrow(outk) > n.path) {
        pathway[[filek]] <- outk[1:n.path, ]
      }
      else {
        pathway[[filek]] <- outk
      }
      pathnum <- c(pathnum, nrow(pathway[[filek]]))
      pathway$dirct <- "down"
      for (i in rownames(pathway[[filek]])) {
        pathcore[[i]] <- msigdb[which(msigdb[, 1] %in% 
                                        i), "gene"]
      }
    }
    write.table(outk, file = file.path(res.path, paste(gsub("_gsea_all_results.txt", 
                                                            "", filek, fixed = TRUE), outlabel, sep = "_")), 
                row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE)
  }
  standarize.fun <- function(indata = NULL, halfwidth = NULL, 
                             centerFlag = TRUE, scaleFlag = TRUE) {
    outdata = t(scale(t(indata), center = centerFlag, scale = scaleFlag))
    if (!is.null(halfwidth)) {
      outdata[outdata > halfwidth] = halfwidth
      outdata[outdata < (-halfwidth)] = -halfwidth
    }
    return(outdata)
  }
  rowmean <- function(x) {
    return(apply(x, 1, mean))
  }
  rowmedian <- function(x) {
    return(apply(x, 1, median))
  }
  if (max(norm.expr) < 25 | (max(norm.expr) >= 25 & min(norm.expr) < 
                             0)) {
    message("--expression profile seems to have been standardised (z-score or log transformation), no more action will be performed.")
    gset <- norm.expr
  }
  if (max(norm.expr) >= 25 & min(norm.expr) >= 0) {
    message("--log2 transformation done for expression data.")
    gset <- log2(norm.expr + 1)
  }
  sam.order <- moic.res$clust.res[order(moic.res$clust.res$clust, 
                                        decreasing = FALSE), "samID"]
  colvec <- clust.col[1:n.moic]
  names(colvec) <- paste0("CS", 1:n.moic)
  annCol <- data.frame(Subtype = paste0("CS", moic.res$clust.res[sam.order, 
                                                                 "clust"]), row.names = sam.order, stringsAsFactors = FALSE)
  annColors <- list(Subtype = colvec)
  es <- GSVA::gsva(expr = as.matrix(gset[, rownames(annCol), 
                                         drop = FALSE]), gset.idx.list = pathcore, method = gsva.method, 
                   parallel.sz = 1)
  es.backup <- es
  es <- standarize.fun(es, halfwidth = 1, centerFlag = TRUE, 
                       scaleFlag = TRUE)
  message(gsva.method, " done...")
  esm <- data.frame(row.names = rownames(es))
  if (norm.method == "mean") {
    for (i in paste0("CS", 1:n.moic)) {
      esm <- cbind.data.frame(esm, data.frame(rowmean(es[, 
                                                         rownames(annCol[which(annCol$Subtype == i), 
                                                                         , drop = FALSE])])))
    }
  }
  if (norm.method == "median") {
    for (i in paste0("CS", 1:n.moic)) {
      esm <- cbind.data.frame(esm, data.frame(rowmedian(es[, 
                                                           rownames(annCol[which(annCol$Subtype == i), 
                                                                           , drop = FALSE])])))
    }
  }
  colnames(esm) <- paste0("CS", 1:n.moic)
  annRow <- data.frame(Subtype = rep(paste0("CS", 1:n.moic), 
                                     pathnum), row.names = rownames(esm), stringsAsFactors = FALSE)
  if (is.null(color)) {
    mapcolor <- (grDevices::colorRampPalette(c("#0000FF", 
                                               "#8080FF", "#FFFFFF", "#FF8080", "#FF0000")))(64)
  }
  else {
    mapcolor <- (grDevices::colorRampPalette(color))(64)
  }
  hm <- ComplexHeatmap::pheatmap(mat = as.matrix(esm), cluster_rows = FALSE, name = name,
                                 cluster_cols = FALSE, show_rownames = TRUE, show_colnames = TRUE, 
                                 annotation_row = annRow, annotation_colors = annColors, 
                                 annotation_names_row = FALSE, legend = TRUE, color = mapcolor, 
                                 border_color = "black", legend_breaks = c(-1, -0.5, 
                                                                           0, 0.5, 1), legend_labels = c(-1, -0.5, 0, 0.5, 
                                                                                                         1), cellwidth = 15, cellheight = 10)
  ComplexHeatmap::draw(hm, annotation_legend_side = "left", heatmap_legend_side = "left")
  if (is.null(fig.name)) {
    outFig <- paste0("gseaheatmap_using_", dirct, "regulated_pathways.pdf")
  }
  else {
    outFig <- paste0(fig.name, "_using_", dirct, "regulated_pathways.pdf")
  }
  pdf(file.path(fig.path, outFig), width = width, height = height)
  ComplexHeatmap::draw(hm, annotation_legend_side = "left", heatmap_legend_side = "left")
  invisible(dev.off())
  message("heatmap done...")
  return(list(gsea.list = gsea.list, raw.es = es.backup, scaled.es = es, 
              grouped.es = esm, heatmap = hm))
}

# runGSEA for newer GSVA versions
runGSEA_mod_4.4 <- function (moic.res = NULL, dea.method = c("deseq2", "edger", 
                                                             "limma"), norm.expr = NULL, prefix = NULL, dat.path = getwd(), 
                             res.path = getwd(), dirct = "up", n.path = 10, msigdb.path = NULL, 
                             nPerm = 1000, minGSSize = 10, maxGSSize = 500, p.cutoff = 0.05, 
                             p.adj.cutoff = 0.05, gsva.method = "gsva", norm.method = "mean", 
                             clust.col = c("#2EC4B6", "#E71D36", "#FF9F1C", "#BDD5EA", 
                                           "#FFA5AB", "#011627", "#023E8A", "#9D4EDD"), color = NULL, 
                             fig.name = NULL, fig.path = getwd(), width = 15, height = 10, name = NULL) 
{
  comsam <- intersect(moic.res$clust.res$samID, colnames(norm.expr))
  if (length(comsam) == nrow(moic.res$clust.res)) {
    message("--all samples matched.")
  }
  else {
    message(paste0("--", (nrow(moic.res$clust.res) - length(comsam)), 
                   " samples mismatched from current subtypes."))
  }
  moic.res$clust.res <- moic.res$clust.res[comsam, , drop = FALSE]
  norm.expr <- norm.expr[, comsam]
  n.moic <- length(unique(moic.res$clust.res$clust))
  mo.method <- moic.res$mo.method
  DEpattern <- paste(mo.method, "_", ifelse(is.null(prefix), 
                                            "", paste0(prefix, "_")), dea.method, ".*._vs_Others.txt$", 
                     sep = "")
  DEfiles <- dir(dat.path, pattern = DEpattern)
  if (length(DEfiles) == 0) {
    stop("no DEfiles!")
  }
  if (length(DEfiles) != n.moic) {
    stop("not all multi-omics clusters have DEfile!")
  }
  if (!is.element(dirct, c("up", "down"))) {
    stop("dirct type error! Allowed value contains c('up', 'down').")
  }
  if (!is.element(gsva.method, c("gsva", "ssgsea", "zscore", 
                                 "plage"))) {
    stop("GSVA method error! Allowed value contains c('gsva', 'ssgsea', 'zscore', 'plage').")
  }
  if (dirct == "up") {
    outlabel <- "unique_upexpr_pathway.txt"
  }
  if (dirct == "down") {
    outlabel <- "unique_downexpr_pathway.txt"
  }
  gsea.list <- list()
  gseaidList <- c()
  for (filek in DEfiles) {
    DEres <- read.table(file.path(dat.path, filek), header = TRUE, 
                        row.names = NULL, sep = "\t", quote = "", stringsAsFactors = FALSE)
    DEres <- DEres[!duplicated(DEres[, 1]), ]
    DEres <- DEres[!is.na(DEres[, 1]), ]
    rownames(DEres) <- DEres[, 1]
    DEres <- DEres[, -1]
    geneList <- DEres$log2fc
    names(geneList) <- rownames(DEres)
    geneList <- sort(geneList, decreasing = TRUE)
    msigdb <- try(clusterProfiler::read.gmt(msigdb.path), 
                  silent = TRUE)
    if (class(msigdb) == "try-error") {
      stop("please provide correct ABSOLUTE PATH for MSigDB file.")
    }
    moic.lab <- paste0("CS", gsub("_vs_Others.txt", "", 
                                  sub(".*CS", "", filek)))
    gsea.list[[moic.lab]] <- suppressWarnings(clusterProfiler::GSEA(geneList = geneList, 
                                                                    TERM2GENE = msigdb, nPerm = nPerm, minGSSize = minGSSize, 
                                                                    maxGSSize = maxGSSize, seed = TRUE, verbose = FALSE, 
                                                                    pvalueCutoff = 1))
    gsea.dat <- as.data.frame(gsea.list[[moic.lab]])
    write.table(gsea.dat[, setdiff(colnames(gsea.dat), "ID")], 
                file = file.path(res.path, paste(gsub("_vs_Others.txt", 
                                                      "", filek, fixed = TRUE), "gsea_all_results.txt", 
                                                 sep = "_")), sep = "\t", row.names = FALSE, 
                col.names = TRUE, quote = FALSE)
    if (dirct == "up") {
      gseaidList <- c(gseaidList, rownames(gsea.dat[which(gsea.dat$NES > 
                                                            0 & gsea.dat$pvalue < p.cutoff & gsea.dat$p.adjust < 
                                                            p.adj.cutoff), ]))
    }
    if (dirct == "down") {
      gseaidList <- c(gseaidList, rownames(gsea.dat[which(gsea.dat$NES < 
                                                            0 & gsea.dat$pvalue < p.cutoff & gsea.dat$p.adjust < 
                                                            p.adj.cutoff), ]))
    }
    unqlist <- setdiff(gseaidList, gseaidList[duplicated(gseaidList)])
  }
  message("GSEA done...")
  GSEApattern <- paste(mo.method, "_", ifelse(is.null(prefix), 
                                              "", paste0(prefix, "_")), dea.method, ".*._gsea_all_results.txt$", 
                       sep = "")
  GSEAfiles <- dir(res.path, pattern = GSEApattern)
  pathway <- pathcore <- list()
  pathnum <- c()
  for (filek in GSEAfiles) {
    GSEAres <- read.table(file.path(res.path, filek), header = TRUE, 
                          row.names = 1, sep = "\t", quote = "", stringsAsFactors = FALSE)
    if (dirct == "up") {
      outk <- intersect(unqlist, rownames(GSEAres[which(GSEAres$NES > 
                                                          0 & GSEAres$pvalue < p.cutoff & GSEAres$p.adjust < 
                                                          p.adj.cutoff), ]))
      outk <- GSEAres[outk, ]
      outk <- outk[order(outk$NES, decreasing = TRUE), 
      ]
      if (nrow(outk) > n.path) {
        pathway[[filek]] <- outk[1:n.path, ]
      }
      else {
        pathway[[filek]] <- outk
      }
      pathnum <- c(pathnum, nrow(pathway[[filek]]))
      pathway$dirct <- "up"
      for (i in rownames(pathway[[filek]])) {
        pathcore[[i]] <- msigdb[which(msigdb[, 1] %in% 
                                        i), "gene"]
      }
    }
    if (dirct == "down") {
      outk <- intersect(unqlist, rownames(GSEAres[which(GSEAres$NES < 
                                                          0 & GSEAres$pvalue < p.cutoff & GSEAres$p.adjust < 
                                                          p.adj.cutoff), ]))
      outk <- GSEAres[outk, ]
      outk <- outk[order(outk$NES, decreasing = FALSE), 
      ]
      if (nrow(outk) > n.path) {
        pathway[[filek]] <- outk[1:n.path, ]
      }
      else {
        pathway[[filek]] <- outk
      }
      pathnum <- c(pathnum, nrow(pathway[[filek]]))
      pathway$dirct <- "down"
      for (i in rownames(pathway[[filek]])) {
        pathcore[[i]] <- msigdb[which(msigdb[, 1] %in% 
                                        i), "gene"]
      }
    }
    write.table(outk, file = file.path(res.path, paste(gsub("_gsea_all_results.txt", 
                                                            "", filek, fixed = TRUE), outlabel, sep = "_")), 
                row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE)
  }
  standarize.fun <- function(indata = NULL, halfwidth = NULL, 
                             centerFlag = TRUE, scaleFlag = TRUE) {
    outdata = t(scale(t(indata), center = centerFlag, scale = scaleFlag))
    if (!is.null(halfwidth)) {
      outdata[outdata > halfwidth] = halfwidth
      outdata[outdata < (-halfwidth)] = -halfwidth
    }
    return(outdata)
  }
  rowmean <- function(x) {
    return(apply(x, 1, mean))
  }
  rowmedian <- function(x) {
    return(apply(x, 1, median))
  }
  if (max(norm.expr) < 25 | (max(norm.expr) >= 25 & min(norm.expr) < 
                             0)) {
    message("--expression profile seems to have been standardised (z-score or log transformation), no more action will be performed.")
    gset <- norm.expr
  }
  if (max(norm.expr) >= 25 & min(norm.expr) >= 0) {
    message("--log2 transformation done for expression data.")
    gset <- log2(norm.expr + 1)
  }
  sam.order <- moic.res$clust.res[order(moic.res$clust.res$clust, 
                                        decreasing = FALSE), "samID"]
  colvec <- clust.col[1:n.moic]
  names(colvec) <- paste0("CS", 1:n.moic)
  annCol <- data.frame(Subtype = paste0("CS", moic.res$clust.res[sam.order, 
                                                                 "clust"]), row.names = sam.order, stringsAsFactors = FALSE)
  annColors <- list(Subtype = colvec)
  es <- GSVA::gsva(param = GSVA::gsvaParam(exprData = as.matrix(gset[, rownames(annCol), 
                                                                     drop = FALSE]),
                                           geneSets = pathcore
                                           ))
  es.backup <- es
  es <- standarize.fun(es, halfwidth = 1, centerFlag = TRUE, 
                       scaleFlag = TRUE)
  message(gsva.method, " done...")
  esm <- data.frame(row.names = rownames(es))
  if (norm.method == "mean") {
    for (i in paste0("CS", 1:n.moic)) {
      esm <- cbind.data.frame(esm, data.frame(rowmean(es[, 
                                                         rownames(annCol[which(annCol$Subtype == i), 
                                                                         , drop = FALSE])])))
    }
  }
  if (norm.method == "median") {
    for (i in paste0("CS", 1:n.moic)) {
      esm <- cbind.data.frame(esm, data.frame(rowmedian(es[, 
                                                           rownames(annCol[which(annCol$Subtype == i), 
                                                                           , drop = FALSE])])))
    }
  }
  colnames(esm) <- paste0("CS", 1:n.moic)
  annRow <- data.frame(Subtype = rep(paste0("CS", 1:n.moic), 
                                     pathnum), row.names = rownames(esm), stringsAsFactors = FALSE)
  if (is.null(color)) {
    mapcolor <- (grDevices::colorRampPalette(c("#0000FF", 
                                               "#8080FF", "#FFFFFF", "#FF8080", "#FF0000")))(64)
  }
  else {
    mapcolor <- (grDevices::colorRampPalette(color))(64)
  }
  hm <- ComplexHeatmap::pheatmap(mat = as.matrix(esm), cluster_rows = FALSE, name = name,
                                 cluster_cols = FALSE, show_rownames = TRUE, show_colnames = TRUE, 
                                 annotation_row = annRow, annotation_colors = annColors, 
                                 annotation_names_row = FALSE, legend = TRUE, color = mapcolor, 
                                 border_color = "black", legend_breaks = c(-1, -0.5, 
                                                                           0, 0.5, 1), legend_labels = c(-1, -0.5, 0, 0.5, 
                                                                                                         1), cellwidth = 15, cellheight = 10)
  ComplexHeatmap::draw(hm, annotation_legend_side = "left", heatmap_legend_side = "left")
  if (is.null(fig.name)) {
    outFig <- paste0("gseaheatmap_using_", dirct, "regulated_pathways.pdf")
  }
  else {
    outFig <- paste0(fig.name, "_using_", dirct, "regulated_pathways.pdf")
  }
  pdf(file.path(fig.path, outFig), width = width, height = height)
  ComplexHeatmap::draw(hm, annotation_legend_side = "left", heatmap_legend_side = "left")
  invisible(dev.off())
  message("heatmap done...")
  return(list(gsea.list = gsea.list, raw.es = es.backup, scaled.es = es, 
              grouped.es = esm, heatmap = hm))
}

# runGSEA for newer GSVA versions
runGSEA_mod_4.4_single_algorithm <- function (algorithm_name = "CS", moic.res = NULL, dea.method = c("deseq2", "edger", 
                                                             "limma"), norm.expr = NULL, prefix = NULL, dat.path = getwd(), 
                             res.path = getwd(), dirct = "up", n.path = 10, msigdb.path = NULL, 
                             nPerm = 1000, minGSSize = 10, maxGSSize = 500, p.cutoff = 0.05, 
                             p.adj.cutoff = 0.05, gsva.method = "gsva", norm.method = "mean", 
                             clust.col = c("#2EC4B6", "#E71D36", "#FF9F1C", "#BDD5EA", 
                                           "#FFA5AB", "#011627", "#023E8A", "#9D4EDD"), color = NULL, 
                             fig.name = NULL, fig.path = getwd(), width = 15, height = 10, name = NULL) 
{
  comsam <- intersect(moic.res$clust.res$samID, colnames(norm.expr))
  if (length(comsam) == nrow(moic.res$clust.res)) {
    message("--all samples matched.")
  }
  else {
    message(paste0("--", (nrow(moic.res$clust.res) - length(comsam)), 
                   " samples mismatched from current subtypes."))
  }
  moic.res$clust.res <- moic.res$clust.res[comsam, , drop = FALSE]
  norm.expr <- norm.expr[, comsam]
  n.moic <- length(unique(moic.res$clust.res$clust))
  mo.method <- moic.res$mo.method
  DEpattern <- paste(mo.method, "_", ifelse(is.null(prefix), 
                                            "", paste0(prefix, "_")), dea.method, ".*._vs_Others.txt$", 
                     sep = "")
  DEfiles <- dir(dat.path, pattern = DEpattern)
  if (length(DEfiles) == 0) {
    stop("no DEfiles!")
  }
  if (length(DEfiles) != n.moic) {
    stop("not all multi-omics clusters have DEfile!")
  }
  if (!is.element(dirct, c("up", "down"))) {
    stop("dirct type error! Allowed value contains c('up', 'down').")
  }
  if (!is.element(gsva.method, c("gsva", "ssgsea", "zscore", 
                                 "plage"))) {
    stop("GSVA method error! Allowed value contains c('gsva', 'ssgsea', 'zscore', 'plage').")
  }
  if (dirct == "up") {
    outlabel <- "unique_upexpr_pathway.txt"
  }
  if (dirct == "down") {
    outlabel <- "unique_downexpr_pathway.txt"
  }
  gsea.list <- list()
  gseaidList <- c()
  for (filek in DEfiles) {
    DEres <- read.table(file.path(dat.path, filek), header = TRUE, 
                        row.names = NULL, sep = "\t", quote = "", stringsAsFactors = FALSE)
    DEres <- DEres[!duplicated(DEres[, 1]), ]
    DEres <- DEres[!is.na(DEres[, 1]), ]
    rownames(DEres) <- DEres[, 1]
    DEres <- DEres[, -1]
    geneList <- DEres$log2fc
    names(geneList) <- rownames(DEres)
    geneList <- sort(geneList, decreasing = TRUE)
    msigdb <- try(clusterProfiler::read.gmt(msigdb.path), 
                  silent = TRUE)
    if (class(msigdb) == "try-error") {
      stop("please provide correct ABSOLUTE PATH for MSigDB file.")
    }
    moic.lab <- paste0(algorithm_name, gsub("_vs_Others.txt", "", 
                                  sub(paste0(".*", algorithm_name), "", filek)))
    gsea.list[[moic.lab]] <- suppressWarnings(clusterProfiler::GSEA(geneList = geneList, 
                                                                    TERM2GENE = msigdb, nPerm = nPerm, minGSSize = minGSSize, 
                                                                    maxGSSize = maxGSSize, seed = TRUE, verbose = FALSE, 
                                                                    pvalueCutoff = 1))
    gsea.dat <- as.data.frame(gsea.list[[moic.lab]])
    write.table(gsea.dat[, setdiff(colnames(gsea.dat), "ID")], 
                file = file.path(res.path, paste(gsub("_vs_Others.txt", 
                                                      "", filek, fixed = TRUE), "gsea_all_results.txt", 
                                                 sep = "_")), sep = "\t", row.names = FALSE, 
                col.names = TRUE, quote = FALSE)
    if (dirct == "up") {
      gseaidList <- c(gseaidList, rownames(gsea.dat[which(gsea.dat$NES > 
                                                            0 & gsea.dat$pvalue < p.cutoff & gsea.dat$p.adjust < 
                                                            p.adj.cutoff), ]))
    }
    if (dirct == "down") {
      gseaidList <- c(gseaidList, rownames(gsea.dat[which(gsea.dat$NES < 
                                                            0 & gsea.dat$pvalue < p.cutoff & gsea.dat$p.adjust < 
                                                            p.adj.cutoff), ]))
    }
    unqlist <- setdiff(gseaidList, gseaidList[duplicated(gseaidList)])
  }
  message("GSEA done...")
  GSEApattern <- paste(mo.method, "_", ifelse(is.null(prefix), 
                                              "", paste0(prefix, "_")), dea.method, ".*._gsea_all_results.txt$", 
                       sep = "")
  GSEAfiles <- dir(res.path, pattern = GSEApattern)
  pathway <- pathcore <- list()
  pathnum <- c()
  for (filek in GSEAfiles) {
    GSEAres <- read.table(file.path(res.path, filek), header = TRUE, 
                          row.names = 1, sep = "\t", quote = "", stringsAsFactors = FALSE)
    if (dirct == "up") {
      outk <- intersect(unqlist, rownames(GSEAres[which(GSEAres$NES > 
                                                          0 & GSEAres$pvalue < p.cutoff & GSEAres$p.adjust < 
                                                          p.adj.cutoff), ]))
      outk <- GSEAres[outk, ]
      outk <- outk[order(outk$NES, decreasing = TRUE), 
      ]
      if (nrow(outk) > n.path) {
        pathway[[filek]] <- outk[1:n.path, ]
      }
      else {
        pathway[[filek]] <- outk
      }
      pathnum <- c(pathnum, nrow(pathway[[filek]]))
      pathway$dirct <- "up"
      for (i in rownames(pathway[[filek]])) {
        pathcore[[i]] <- msigdb[which(msigdb[, 1] %in% 
                                        i), "gene"]
      }
    }
    if (dirct == "down") {
      outk <- intersect(unqlist, rownames(GSEAres[which(GSEAres$NES < 
                                                          0 & GSEAres$pvalue < p.cutoff & GSEAres$p.adjust < 
                                                          p.adj.cutoff), ]))
      outk <- GSEAres[outk, ]
      outk <- outk[order(outk$NES, decreasing = FALSE), 
      ]
      if (nrow(outk) > n.path) {
        pathway[[filek]] <- outk[1:n.path, ]
      }
      else {
        pathway[[filek]] <- outk
      }
      pathnum <- c(pathnum, nrow(pathway[[filek]]))
      pathway$dirct <- "down"
      for (i in rownames(pathway[[filek]])) {
        pathcore[[i]] <- msigdb[which(msigdb[, 1] %in% 
                                        i), "gene"]
      }
    }
    write.table(outk, file = file.path(res.path, paste(gsub("_gsea_all_results.txt", 
                                                            "", filek, fixed = TRUE), outlabel, sep = "_")), 
                row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE)
  }
  standarize.fun <- function(indata = NULL, halfwidth = NULL, 
                             centerFlag = TRUE, scaleFlag = TRUE) {
    outdata = t(scale(t(indata), center = centerFlag, scale = scaleFlag))
    if (!is.null(halfwidth)) {
      outdata[outdata > halfwidth] = halfwidth
      outdata[outdata < (-halfwidth)] = -halfwidth
    }
    return(outdata)
  }
  rowmean <- function(x) {
    return(apply(x, 1, mean))
  }
  rowmedian <- function(x) {
    return(apply(x, 1, median))
  }
  if (max(norm.expr) < 25 | (max(norm.expr) >= 25 & min(norm.expr) < 
                             0)) {
    message("--expression profile seems to have been standardised (z-score or log transformation), no more action will be performed.")
    gset <- norm.expr
  }
  if (max(norm.expr) >= 25 & min(norm.expr) >= 0) {
    message("--log2 transformation done for expression data.")
    gset <- log2(norm.expr + 1)
  }
  sam.order <- moic.res$clust.res[order(moic.res$clust.res$clust, 
                                        decreasing = FALSE), "samID"]
  colvec <- clust.col[1:n.moic]
  names(colvec) <- paste0(algorithm_name, 1:n.moic)
  annCol <- data.frame(Subtype = paste0(algorithm_name, moic.res$clust.res[sam.order, 
                                                                 "clust"]), row.names = sam.order, stringsAsFactors = FALSE)
  annColors <- list(Subtype = colvec)
  es <- GSVA::gsva(param = GSVA::gsvaParam(exprData = as.matrix(gset[, rownames(annCol), 
                                                                     drop = FALSE]),
                                           geneSets = pathcore
  ))
  es.backup <- es
  es <- standarize.fun(es, halfwidth = 1, centerFlag = TRUE, 
                       scaleFlag = TRUE)
  message(gsva.method, " done...")
  esm <- data.frame(row.names = rownames(es))
  if (norm.method == "mean") {
    for (i in paste0(algorithm_name, 1:n.moic)) {
      esm <- cbind.data.frame(esm, data.frame(rowmean(es[, 
                                                         rownames(annCol[which(annCol$Subtype == i), 
                                                                         , drop = FALSE])])))
    }
  }
  if (norm.method == "median") {
    for (i in paste0(algorithm_name, 1:n.moic)) {
      esm <- cbind.data.frame(esm, data.frame(rowmedian(es[, 
                                                           rownames(annCol[which(annCol$Subtype == i), 
                                                                           , drop = FALSE])])))
    }
  }
  colnames(esm) <- paste0(algorithm_name, 1:n.moic)
  annRow <- data.frame(Subtype = rep(paste0(algorithm_name, 1:n.moic), 
                                     pathnum), row.names = rownames(esm), stringsAsFactors = FALSE)
  if (is.null(color)) {
    mapcolor <- (grDevices::colorRampPalette(c("#0000FF", 
                                               "#8080FF", "#FFFFFF", "#FF8080", "#FF0000")))(64)
  }
  else {
    mapcolor <- (grDevices::colorRampPalette(color))(64)
  }
  hm <- ComplexHeatmap::pheatmap(mat = as.matrix(esm), cluster_rows = FALSE, name = name,
                                 cluster_cols = FALSE, show_rownames = TRUE, show_colnames = TRUE, 
                                 annotation_row = annRow, annotation_colors = annColors, 
                                 annotation_names_row = FALSE, legend = TRUE, color = mapcolor, 
                                 border_color = "black", legend_breaks = c(-1, -0.5, 
                                                                           0, 0.5, 1), legend_labels = c(-1, -0.5, 0, 0.5, 
                                                                                                         1), cellwidth = 15, cellheight = 10)
  ComplexHeatmap::draw(hm, annotation_legend_side = "left", heatmap_legend_side = "left")
  if (is.null(fig.name)) {
    outFig <- paste0("gseaheatmap_using_", dirct, "regulated_pathways.pdf")
  }
  else {
    outFig <- paste0(fig.name, "_using_", dirct, "regulated_pathways.pdf")
  }
  pdf(file.path(fig.path, outFig), width = width, height = height)
  ComplexHeatmap::draw(hm, annotation_legend_side = "left", heatmap_legend_side = "left")
  invisible(dev.off())
  message("heatmap done...")
  return(list(gsea.list = gsea.list, raw.es = es.backup, scaled.es = es, 
              grouped.es = esm, heatmap = hm))
}


# Redefine compAgree for subtype comparisons across different classifications #####

# Modification:
# Added three library calls in the beginning of the function definition
compAgree2 = function (moic.res = NULL, subt2comp = NULL, doPlot = TRUE, 
                       clust.col = c("#2EC4B6", "#E71D36", "#FF9F1C", "#BDD5EA", 
                                     "#FFA5AB", "#011627", "#023E8A", "#9D4EDD"), box.width = 0.1, 
                       fig.name = NULL, fig.path = getwd(), width = 6, height = 5) 
{
  library(ggplot2)
  library(ggalluvial)
  library(cowplot)
  dat <- moic.res$clust.res
  colnames(dat)[2] <- "Subtype"
  dat$Subtype <- paste0("CS", dat$Subtype)
  comsam <- intersect(dat$samID, rownames(subt2comp))
  if (length(comsam) == nrow(dat)) {
    message("--all samples matched.")
  }
  else {
    message(paste0("--", (nrow(dat) - length(comsam)), " samples mismatched from current subtypes."))
  }
  dat <- cbind.data.frame(Subtype = dat[comsam, "Subtype", 
                                        drop = FALSE], subt2comp[comsam, , drop = FALSE])
  dat <- as.data.frame(na.omit(dat))
  if (nrow(dat) != nrow(moic.res$clust.res)) {
    message("--removed NA values in subt2comp.")
  }
  var <- colnames(dat)
  n.var <- length(var)
  if (n.var > 6) {
    stop("please indicate less than 6 subtypes (including current subtypes) that need to compare.")
  }
  outTab <- NULL
  c1 <- as.vector(as.numeric(factor(dat[, 1])))
  for (i in 2:ncol(dat)) {
    c2 <- as.vector(as.numeric(factor(dat[, i])))
    RI <- flexclust::comPart(c1, c2, type = c("RI"))
    AMI <- aricode::AMI(c1, c2)
    JI <- flexclust::comPart(c1, c2, type = c("J"))
    FM <- flexclust::comPart(c1, c2, type = c("FM"))
    outTab <- rbind.data.frame(outTab, data.frame(current.subtype = colnames(dat)[1], 
                                                  other.subtype = colnames(dat)[i], RI = as.numeric(RI), 
                                                  AMI = as.numeric(AMI), JI = as.numeric(JI), FM = as.numeric(FM), 
                                                  stringsAsFactors = FALSE), stringsAsFactors = FALSE)
  }
  assign("StatStratum", ggalluvial::StatStratum, envir = globalenv())
  if (doPlot) {
    if (is.null(fig.name)) {
      outFig <- "Agreement between current subtype and other classifications.pdf"
    }
    else {
      outFig <- paste0(fig.name, ".pdf")
    }
    agreement <- reshape2::melt(outTab[, 2:ncol(outTab)], 
                                id.vars = "other.subtype", variable.name = "Method")
    b <- ggplot(data = agreement, aes(x = Method, y = value, 
                                      fill = other.subtype)) + geom_bar(stat = "identity", 
                                                                        position = position_dodge()) + scale_fill_brewer(palette = "Set1") + 
      ggplot2::labs(x = "", y = "Scalar") + scale_y_continuous(limits = c(0, 
                                                                          1), expand = c(0, 0)) + theme_bw() + theme(legend.position = "top", 
                                                                                                                     legend.title = element_blank(), panel.grid = element_blank(), 
                                                                                                                     axis.ticks = element_blank(), axis.text.x = element_text(color = "black", 
                                                                                                                                                                              size = 12, face = "bold", vjust = -5), axis.text.y = element_text(color = "black", 
                                                                                                                                                                                                                                                size = 12, face = "bold"), axis.title.x = element_text(color = "black", 
                                                                                                                                                                                                                                                                                                       size = 12, face = "bold"), axis.title.y = element_text(color = "black", 
                                                                                                                                                                                                                                                                                                                                                              size = 12, face = "bold")) + ggtitle("")
    col = clust.col[1:length(unique(dat$Subtype))]
    var1 <- var[1]
    if (n.var == 2) {
      var2 <- var[2]
      subdf <- dat[, 1:n.var]
      colnames(subdf) <- c("Subtype", paste0("Subtype", 
                                             1:(n.var - 1)))
      subdf <- subdf %>% group_by(Subtype, Subtype1) %>% 
        tally(name = "Freq") %>% as.data.frame()
      p <- ggplot(subdf, aes(y = Freq, axis1 = Subtype, 
                             axis2 = Subtype1)) + scale_fill_manual(values = col) + 
        geom_flow(stat = "alluvium", width = 1/8, aes(fill = Subtype)) + 
        geom_stratum(width = 1/8, reverse = TRUE) + 
        geom_text(stat = "stratum", aes(label = after_stat(stratum)), 
                  reverse = TRUE) + scale_x_continuous(breaks = 1:n.var, 
                                                       labels = c(var1, var2)) + theme_bw() + theme(legend.position = "top", 
                                                                                                    legend.title = element_blank(), panel.grid = element_blank(), 
                                                                                                    panel.border = element_blank(), axis.title.x = element_blank(), 
                                                                                                    axis.title.y = element_blank(), axis.text.y = element_blank(), 
                                                                                                    axis.ticks = element_blank(), axis.text.x = element_text(size = 12, 
                                                                                                                                                             face = "bold", color = "black")) + ggtitle("")
    }
    if (n.var == 3) {
      var2 <- var[2]
      var3 <- var[3]
      subdf <- dat[, 1:n.var]
      colnames(subdf) <- c("Subtype", paste0("Subtype", 
                                             1:(n.var - 1)))
      subdf <- subdf %>% group_by(Subtype, Subtype1, Subtype2) %>% 
        tally(name = "Freq") %>% as.data.frame()
      p <- ggplot(subdf, aes(y = Freq, axis1 = Subtype, 
                             axis2 = Subtype1, axis3 = Subtype2)) + scale_fill_manual(values = col) + 
        geom_flow(stat = "alluvium", width = 1/8, aes(fill = Subtype)) + 
        geom_stratum(width = box.width, reverse = TRUE) + 
        geom_text(stat = "stratum", aes(label = after_stat(stratum)), 
                  reverse = TRUE) + scale_x_continuous(breaks = 1:n.var, 
                                                       labels = c(var1, var2, var3)) + theme_bw() + 
        theme(legend.position = "top", legend.title = element_blank(), 
              panel.grid = element_blank(), panel.border = element_blank(), 
              axis.title.x = element_blank(), axis.title.y = element_blank(), 
              axis.text.y = element_blank(), axis.ticks = element_blank(), 
              axis.text.x = element_text(size = 12, face = "bold", 
                                         color = "black")) + ggtitle("")
    }
    if (n.var == 4) {
      var2 <- var[2]
      var3 <- var[3]
      var4 <- var[4]
      subdf <- dat[, 1:n.var]
      colnames(subdf) <- c("Subtype", paste0("Subtype", 
                                             1:(n.var - 1)))
      subdf <- subdf %>% group_by(Subtype, Subtype1, Subtype2, 
                                  Subtype3) %>% tally(name = "Freq") %>% as.data.frame()
      p <- ggplot(subdf, aes(y = Freq, axis1 = Subtype, 
                             axis2 = Subtype1, axis3 = Subtype2, axis4 = Subtype3)) + 
        scale_fill_manual(values = col) + geom_flow(stat = "alluvium", 
                                                    width = 1/8, aes(fill = Subtype)) + geom_stratum(width = 1/8, 
                                                                                                     reverse = TRUE) + geom_text(stat = "stratum", 
                                                                                                                                 aes(label = after_stat(stratum)), reverse = TRUE) + 
        scale_x_continuous(breaks = 1:n.var, labels = c(var1, 
                                                        var2, var3, var4)) + theme_bw() + theme(legend.position = "top", 
                                                                                                legend.title = element_blank(), panel.grid = element_blank(), 
                                                                                                panel.border = element_blank(), axis.title.x = element_blank(), 
                                                                                                axis.title.y = element_blank(), axis.text.y = element_blank(), 
                                                                                                axis.ticks = element_blank(), axis.text.x = element_text(size = 12, 
                                                                                                                                                         face = "bold", color = "black")) + ggtitle("")
    }
    if (n.var == 5) {
      var2 <- var[2]
      var3 <- var[3]
      var4 <- var[4]
      var5 <- var[5]
      subdf <- dat[, 1:n.var]
      colnames(subdf) <- c("Subtype", paste0("Subtype", 
                                             1:(n.var - 1)))
      subdf <- subdf %>% group_by(Subtype, Subtype1, Subtype2, 
                                  Subtype3, Subtype4) %>% tally(name = "Freq") %>% 
        as.data.frame()
      p <- ggplot(subdf, aes(y = Freq, axis1 = Subtype, 
                             axis2 = Subtype1, axis3 = Subtype2, axis4 = Subtype3, 
                             axis5 = Subtype4)) + scale_fill_manual(values = col) + 
        geom_flow(stat = "alluvium", width = 1/8, aes(fill = Subtype)) + 
        geom_stratum(width = 1/8, reverse = TRUE) + 
        geom_text(stat = "stratum", aes(label = after_stat(stratum)), 
                  reverse = TRUE) + scale_x_continuous(breaks = 1:n.var, 
                                                       labels = c(var1, var2, var3, var4, var5)) + 
        theme_bw() + theme(legend.position = "top", 
                           legend.title = element_blank(), panel.grid = element_blank(), 
                           panel.border = element_blank(), axis.title.x = element_blank(), 
                           axis.title.y = element_blank(), axis.text.y = element_blank(), 
                           axis.ticks = element_blank(), axis.text.x = element_text(size = 12, 
                                                                                    face = "bold", color = "black")) + ggtitle("")
    }
    if (n.var == 6) {
      var2 <- var[2]
      var3 <- var[3]
      var4 <- var[4]
      var5 <- var[5]
      var6 <- var[6]
      subdf <- dat[, 1:n.var]
      colnames(subdf) <- c("Subtype", paste0("Subtype", 
                                             1:(n.var - 1)))
      subdf <- subdf %>% group_by(Subtype, Subtype1, Subtype2, 
                                  Subtype3, Subtype4, Subtype5) %>% tally(name = "Freq") %>% 
        as.data.frame()
      p <- ggplot(subdf, aes(y = Freq, axis1 = Subtype, 
                             axis2 = Subtype1, axis3 = Subtype2, axis4 = Subtype3, 
                             axis5 = Subtype4, axis6 = Subtype5)) + scale_fill_manual(values = col) + 
        geom_flow(stat = "alluvium", width = 1/8, aes(fill = Subtype)) + 
        geom_stratum(width = 1/8, reverse = TRUE) + 
        geom_text(stat = "stratum", aes(label = after_stat(stratum)), 
                  reverse = TRUE) + scale_x_continuous(breaks = 1:n.var, 
                                                       labels = c(var1, var2, var3, var4, var5, var6)) + 
        theme_bw() + theme(legend.position = "top", 
                           legend.title = element_blank(), panel.grid = element_blank(), 
                           panel.border = element_blank(), axis.title.x = element_blank(), 
                           axis.title.y = element_blank(), axis.text.y = element_blank(), 
                           axis.ticks = element_blank(), axis.text.x = element_text(size = 12, 
                                                                                    face = "bold", color = "black")) + ggtitle("")
    }
    agree <- list(b, p)
    bp <- plot_grid(plotlist = agree, ncol = 2)
    ggsave(file.path(fig.path, outFig), width = width, height = height)
    print(bp)
  }
  return(outTab)
}

# compClinvar modification #####
# Changed the warn variable to character(0) due to unnecessary error stops
compClinvar2 = function (moic.res = NULL, var2comp = NULL, strata = NULL, factorVars = NULL, 
          nonnormalVars = NULL, exactVars = NULL, includeNA = FALSE, 
          doWord = TRUE, tab.name = NULL, res.path = getwd(), ...) 
{
  dat <- moic.res$clust.res
  colnames(dat)[which(colnames(dat) == "clust")] <- "Subtype"
  dat$Subtype <- paste0("CS", dat$Subtype)
  com_sam <- intersect(dat$samID, rownames(var2comp))
  if (length(com_sam) == nrow(dat)) {
    message("--all samples matched.")
  }
  else {
    message(paste0("--", (nrow(dat) - length(com_sam)), 
                   " samples mismatched from current subtypes."))
  }
  dat <- cbind.data.frame(Subtype = dat[com_sam, "Subtype", 
                                        drop = FALSE], var2comp[com_sam, , drop = FALSE])
  if (is.null(strata)) {
    strata <- "Subtype"
  }
  if (!is.element(strata, colnames(dat))) {
    stop("fail to find this strata in var2comp. Consider using NULL by default.")
  }
  warn <- character(0)  # Initialize `warn` as an empty character vector
  
  tryCatch(stabl <- jstable::CreateTableOne2(vars = setdiff(colnames(dat), 
                                                            strata), strata = strata, data = dat, factorVars = factorVars, 
                                             nonnormal = nonnormalVars, exact = exactVars, includeNA = includeNA, 
                                             showAllLevels = TRUE, ...), warning = function(w) {
                                               warn <<- append(warn, conditionMessage(w))
                                             })
  
  if (length(warn) > 0 && grepl("NA", warn, fixed = TRUE)) {
    set.seed(19991018)
    stabl <- jstable::CreateTableOne2(vars = setdiff(colnames(dat), 
                                                     strata), strata = strata, data = dat, factorVars = factorVars, 
                                      nonnormal = nonnormalVars, exact = exactVars, includeNA = includeNA, 
                                      showAllLevels = TRUE, argsExact = list(simulate.p.value = T), 
                                      ...)
  }
  else {
    stabl <- jstable::CreateTableOne2(vars = setdiff(colnames(dat), 
                                                     strata), strata = strata, data = dat, factorVars = factorVars, 
                                      nonnormal = nonnormalVars, exact = exactVars, includeNA = includeNA, 
                                      showAllLevels = TRUE, ...)
  }
  comtable <- as.data.frame(stabl)
  comtable <- cbind.data.frame(var = rownames(stabl), comtable)
  rownames(comtable) <- NULL
  colnames(comtable)[1] <- " "
  comtable[is.na(comtable)] <- ""
  comtable <- comtable[, setdiff(colnames(comtable), "sig")]
  if (is.null(tab.name)) {
    outFile <- "summarization of clinical variables stratified by current subtype.txt"
  }
  else {
    outFile <- paste0(tab.name, ".txt")
  }
  write.table(stabl, file.path(res.path, outFile), sep = "\t", 
              quote = FALSE)
  if (doWord) {
    table_subtitle <- colnames(comtable)
    title_name <- paste0("Table *. ", gsub(".txt", "", outFile, 
                                           fixed = TRUE))
    mynote <- "Note: ..."
    my_doc <- officer::read_docx()
    my_doc %>% officer::body_add_par(value = title_name, 
                                     style = "table title") %>% officer::body_add_table(value = comtable, 
                                                                                        style = "table_template") %>% officer::body_add_par(value = mynote) %>% 
      print(target = file.path(res.path, paste0("TABLE ", 
                                                gsub(".txt", "", outFile, fixed = TRUE), ".docx")))
  }
  return(list(compTab = comtable))
}
