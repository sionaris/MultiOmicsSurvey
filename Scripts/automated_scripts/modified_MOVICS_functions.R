# runDEA: custom output name modification #####
runDEA_mod = function (dea.method = c("deseq2", "edger", "limma"), expr = NULL, 
                       moic.res = NULL, prefix = NULL, overwt = FALSE, sort.p = TRUE, 
                       verbose = TRUE, res.path = getwd(), algorithm = "CS") 
{
  if (!is.element(dea.method, c("deseq2", "edger", "limma"))) {
    stop("unsupported algorithm: dea.method should be one of 'deseq2', 'edger', or 'limma'.")
  }
  method <- dea.method[1]
  comsam <- intersect(moic.res$clust.res$samID, colnames(expr))
  if (length(comsam) == nrow(moic.res$clust.res)) {
    message("--all samples matched.")
  }
  else {
    message(paste0("--", (nrow(moic.res$clust.res) - length(comsam)), 
                   " samples mismatched from current subtypes."))
  }
  moic.res$clust.res <- moic.res$clust.res[comsam, , drop = FALSE]
  expr <- expr[, comsam]
  rundea_mod <- switch(method, deseq2 = twoclassdeseq2_mod, edger = twoclassedger_mod, 
                   limma = twoclasslimma_mod)
  if (method %in% c("deseq2", "edger")) {
    message(paste0("--you choose ", method, " and please make sure an RNA-Seq count data was provided."))
    rundea_mod(moic.res = moic.res, countsTable = expr, prefix = prefix, 
           overwt = overwt, sort.p = sort.p, verbose = verbose, 
           res.path = res.path, algorithm = algorithm)
  }
  else {
    message(paste0("--you choose ", method, " and please make sure a microarray profile or a normalized expression data [FPKM or TPM without log2 transformation is recommended] was provided."))
    rundea_mod(moic.res = moic.res, norm.expr = expr, prefix = prefix, 
           overwt = overwt, sort.p = sort.p, verbose = verbose, 
           res.path = res.path, algorithm = algorithm)
  }
}

twoclassdeseq2_mod = function (moic.res = NULL, countsTable = NULL, prefix = NULL, 
                               overwt = FALSE, sort.p = TRUE, verbose = TRUE, res.path = getwd(),
                               algorithm = "CS") 
{
  createList <- function(moic.res = NULL) {
    mo.method <- moic.res$mo.method
    moic.res <- moic.res$clust.res
    tumorsam <- moic.res$samID
    sampleList = list()
    treatsamList = list()
    treatnameList <- c()
    ctrlnameList <- c()
    n.moic <- length(unique(moic.res$clust))
    for (i in 1:n.moic) {
      sampleList[[i]] <- tumorsam
      treatsamList[[i]] = intersect(tumorsam, moic.res[which(moic.res$clust == 
                                                               i), "samID"])
      treatnameList[i] <- paste0(algorithm, i)
      ctrlnameList[i] <- "Others"
    }
    return(list(sampleList, treatsamList, treatnameList, 
                ctrlnameList, mo.method))
  }
  complist <- createList(moic.res = moic.res)
  sampleList <- complist[[1]]
  treatsamList <- complist[[2]]
  treatnameList <- complist[[3]]
  ctrlnameList <- complist[[4]]
  mo.method <- complist[[5]]
  allsamples <- colnames(countsTable)
  options(warn = 1)
  for (k in 1:length(sampleList)) {
    samples <- sampleList[[k]]
    treatsam <- treatsamList[[k]]
    treatname <- treatnameList[k]
    ctrlname <- ctrlnameList[k]
    compname <- paste(treatname, "_vs_", ctrlname, sep = "")
    tmp <- rep("others", times = length(allsamples))
    names(tmp) <- allsamples
    tmp[samples] <- "control"
    tmp[treatsam] <- "treatment"
    if (!is.null(prefix)) {
      outfile <- file.path(res.path, paste(mo.method, 
                                           "_", prefix, "_deseq2_test_result.", compname, 
                                           ".txt", sep = ""))
    }
    else {
      outfile <- file.path(res.path, paste(mo.method, 
                                           "_deseq2_test_result.", compname, ".txt", sep = ""))
    }
    if (file.exists(outfile) & (overwt == FALSE)) {
      cat(paste0("deseq2 of ", compname, " exists and skipped...\n"))
      next
    }
    saminfo <- data.frame(Type = as.factor(tmp[samples]), 
                          SampleID = samples, stringsAsFactors = FALSE)
    cts <- countsTable[, samples]
    coldata <- saminfo[samples, ]
    dds <- quiet(DESeq2::DESeqDataSetFromMatrix(countData = cts, 
                                                colData = coldata, design = as.formula("~ Type")))
    dds$Type <- relevel(dds$Type, ref = "control")
    dds <- quiet(DESeq2::DESeq(dds))
    res <- DESeq2::results(dds, contrast = c("Type", "treatment", 
                                             "control"))
    if (sort.p) {
      resData <- as.data.frame(res[order(res$padj), ])
    }
    else {
      resData <- as.data.frame(res)
    }
    resData$id <- rownames(resData)
    resData <- resData[, c("id", "baseMean", "log2FoldChange", 
                           "lfcSE", "stat", "pvalue", "padj")]
    colnames(resData) <- c("id", "baseMean", "log2fc", "lfcSE", 
                           "stat", "pvalue", "padj")
    resData$fc <- 2^resData$log2fc
    if (verbose) {
      resData <- resData[, c("id", "fc", "log2fc", "pvalue", 
                             "padj")]
    }
    else {
      resData <- resData[, c("id", "fc", "log2fc", "baseMean", 
                             "lfcSE", "stat", "pvalue", "padj")]
    }
    write.table(resData, file = outfile, row.names = FALSE, 
                sep = "\t", quote = FALSE)
    cat(paste0("deseq2 of ", compname, " done...\n"))
  }
  options(warn = 0)
}

twoclassedger_mod = function (moic.res = NULL, countsTable = NULL, prefix = NULL, 
                              overwt = FALSE, sort.p = TRUE, verbose = TRUE, res.path = getwd(),
                              algorithm = "CS") 
{
  createList <- function(moic.res = NULL) {
    mo.method <- moic.res$mo.method
    moic.res <- moic.res$clust.res
    tumorsam <- moic.res$samID
    sampleList <- list()
    treatsamList <- list()
    treatnameList <- c()
    ctrlnameList <- c()
    n.moic <- length(unique(moic.res$clust))
    for (i in 1:n.moic) {
      sampleList[[i]] <- tumorsam
      treatsamList[[i]] <- intersect(tumorsam, moic.res[which(moic.res$clust == 
                                                                i), "samID"])
      treatnameList[i] <- paste0(algorithm, i)
      ctrlnameList[i] <- "Others"
    }
    return(list(sampleList, treatsamList, treatnameList, 
                ctrlnameList, mo.method))
  }
  complist <- createList(moic.res = moic.res)
  sampleList <- complist[[1]]
  treatsamList <- complist[[2]]
  treatnameList <- complist[[3]]
  ctrlnameList <- complist[[4]]
  mo.method <- complist[[5]]
  allsamples <- colnames(countsTable)
  options(warn = 1)
  for (k in 1:length(sampleList)) {
    samples <- sampleList[[k]]
    treatsam <- treatsamList[[k]]
    treatname <- treatnameList[k]
    ctrlname <- ctrlnameList[k]
    compname <- paste(treatname, "_vs_", ctrlname, sep = "")
    tmp <- rep("others", times = length(allsamples))
    names(tmp) <- allsamples
    tmp[samples] <- "control"
    tmp[treatsam] <- "treatment"
    if (!is.null(prefix)) {
      outfile <- file.path(res.path, paste(mo.method, 
                                           "_", prefix, "_edger_test_result.", compname, 
                                           ".txt", sep = ""))
    }
    else {
      outfile <- file.path(res.path, paste(mo.method, 
                                           "_edger_test_result.", compname, ".txt", sep = ""))
    }
    if (file.exists(outfile) & (overwt == FALSE)) {
      cat(paste0("edger of ", compname, " exists and skipped...\n"))
      next
    }
    saminfo <- data.frame(Type = tmp[samples], SampleID = samples, 
                          stringsAsFactors = FALSE)
    group = factor(saminfo$Type, levels = c("control", "treatment"))
    design <- model.matrix(~group)
    rownames(design) <- samples
    y <- edgeR::DGEList(counts = countsTable[, samples], 
                        group = saminfo$Type)
    y <- edgeR::calcNormFactors(y)
    y <- edgeR::estimateDisp(y, design, robust = TRUE)
    fit <- edgeR::glmFit(y, design)
    lrt <- edgeR::glmLRT(fit)
    ordered_tags <- edgeR::topTags(lrt, n = 1e+05)
    allDiff <- ordered_tags$table
    allDiff <- allDiff[is.na(allDiff$FDR) == FALSE, ]
    diff <- allDiff
    diff$id <- rownames(diff)
    resData <- diff[, c("id", "logFC", "logCPM", "LR", "PValue", 
                        "FDR")]
    colnames(resData) <- c("id", "log2fc", "logCPM", "LR", 
                           "pvalue", "padj")
    resData$fc <- 2^resData$log2fc
    if (sort.p) {
      resData <- as.data.frame(resData[order(resData$padj), 
      ])
    }
    else {
      resData <- as.data.frame(resData)
    }
    if (verbose) {
      resData <- resData[, c("id", "fc", "log2fc", "pvalue", 
                             "padj")]
    }
    else {
      resData <- resData[, c("id", "fc", "log2fc", "logCPM", 
                             "LR", "pvalue", "padj")]
    }
    write.table(resData, file = outfile, row.names = FALSE, 
                sep = "\t", quote = FALSE)
    cat(paste0("edger of ", compname, " done...\n"))
  }
  options(warn = 0)
}

twoclasslimma_mod = function (moic.res = NULL, norm.expr = NULL, prefix = NULL, 
                              overwt = FALSE, sort.p = TRUE, verbose = TRUE, res.path = getwd(),
                              algorithm = "CS") 
{
  createList <- function(moic.res = NULL) {
    mo.method <- moic.res$mo.method
    moic.res <- moic.res$clust.res
    tumorsam <- moic.res$samID
    sampleList <- list()
    treatsamList <- list()
    treatnameList <- c()
    ctrlnameList <- c()
    n.moic <- length(unique(moic.res$clust))
    for (i in 1:n.moic) {
      sampleList[[i]] <- tumorsam
      treatsamList[[i]] <- intersect(tumorsam, moic.res[which(moic.res$clust == 
                                                                i), "samID"])
      treatnameList[i] <- paste0(algorithm, i)
      ctrlnameList[i] <- "Others"
    }
    return(list(sampleList, treatsamList, treatnameList, 
                ctrlnameList, mo.method))
  }
  complist <- createList(moic.res = moic.res)
  sampleList <- complist[[1]]
  treatsamList <- complist[[2]]
  treatnameList <- complist[[3]]
  ctrlnameList <- complist[[4]]
  mo.method <- complist[[5]]
  allsamples <- colnames(norm.expr)
  if (max(norm.expr) < 25 | (max(norm.expr) >= 25 & min(norm.expr) < 
                             0)) {
    message("--expression profile seems to have been standardised (z-score or log transformation), no more action will be performed.")
    gset <- norm.expr
  }
  if (max(norm.expr) >= 25 & min(norm.expr) >= 0) {
    message("--log2 transformation done for expression data.")
    gset <- log2(norm.expr + 1)
  }
  options(warn = 1)
  for (k in 1:length(sampleList)) {
    samples <- sampleList[[k]]
    treatsam <- treatsamList[[k]]
    treatname <- treatnameList[k]
    ctrlname <- ctrlnameList[k]
    compname <- paste(treatname, "_vs_", ctrlname, sep = "")
    tmp <- rep("others", times = length(allsamples))
    names(tmp) <- allsamples
    tmp[samples] <- "control"
    tmp[treatsam] <- "treatment"
    if (!is.null(prefix)) {
      outfile <- file.path(res.path, paste(mo.method, 
                                           "_", prefix, "_limma_test_result.", compname, 
                                           ".txt", sep = ""))
    }
    else {
      outfile <- file.path(res.path, paste(mo.method, 
                                           "_limma_test_result.", compname, ".txt", sep = ""))
    }
    if (file.exists(outfile) & (overwt == FALSE)) {
      cat(paste0("limma of ", compname, " exists and skipped...\n"))
      next
    }
    pd <- data.frame(Samples = names(tmp), Group = as.character(tmp), 
                     stringsAsFactors = FALSE)
    design <- model.matrix(~-1 + factor(pd$Group, levels = c("treatment", 
                                                             "control")))
    colnames(design) <- c("treatment", "control")
    fit <- limma::lmFit(gset, design = design)
    contrastsMatrix <- limma::makeContrasts(treatment - 
                                              control, levels = c("treatment", "control"))
    fit2 <- limma::contrasts.fit(fit, contrasts = contrastsMatrix)
    fit2 <- limma::eBayes(fit2, 0.01)
    resData <- limma::topTable(fit2, adjust = "fdr", sort.by = "B", 
                               number = 1e+05)
    resData <- as.data.frame(subset(resData, select = c("logFC", 
                                                        "t", "B", "P.Value", "adj.P.Val")))
    resData$id <- rownames(resData)
    colnames(resData) <- c("log2fc", "t", "B", "pvalue", 
                           "padj", "id")
    resData$fc <- 2^resData$log2fc
    if (sort.p) {
      resData <- resData[order(resData$padj), ]
    }
    else {
      resData <- as.data.frame(resData)
    }
    if (verbose) {
      resData <- resData[, c("id", "fc", "log2fc", "pvalue", 
                             "padj")]
    }
    else {
      resData <- resData[, c("id", "fc", "log2fc", "t", 
                             "B", "pvalue", "padj")]
    }
    write.table(resData, file = outfile, row.names = FALSE, 
                sep = "\t", quote = FALSE)
    cat(paste0("limma of ", compname, " done...\n"))
  }
  options(warn = 0)
}

# runMarker single algorithm #####

# Just changing the CS labels when focusing on a single algorithm
runMarker_single_algorithm = function (algorithm_name = "CS", moic.res = NULL, dea.method = c("deseq2", "edger", 
                                          "limma"), prefix = NULL, dat.path = getwd(), res.path = getwd(), 
          p.cutoff = 0.05, p.adj.cutoff = 0.05, dirct = "up", n.marker = 200, 
          doplot = TRUE, norm.expr = NULL, annCol = NULL, annColors = NULL, 
          clust.col = c("#2EC4B6", "#E71D36", "#FF9F1C", "#BDD5EA", 
                        "#FFA5AB", "#011627", "#023E8A", "#9D4EDD", "#f09c6c", "#09f3b3"), halfwidth = 3, 
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
    
    # Check if the DEres is empty or has only one gene
    if (nrow(DEres) < 2) {
      stop(paste("Skipping file", filek, "because it has less than two rows of data."))
    }
    
    # Proceed with the existing code
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
    
    # Check if the DEres is empty or has only one gene
    if (nrow(DEres) < 2) {
      stop(paste("Skipping file", filek, "because it has less than two rows of data."))
    }
    
    # Proceed with the existing code
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
    tmp <- NULL
    if (!is.null(marker[[filek]]) && nrow(marker[[filek]]) > 0) {
      tmp <- data.frame(probe = rownames(marker[[filek]]), 
                        class = sub("_vs_Others.txt", "", sub(".*.result.", "", filek)), 
                        dirct = marker$dirct, stringsAsFactors = FALSE)
    }
    if (!is.null(tmp)) {
      templates <- rbind.data.frame(templates, tmp, stringsAsFactors = FALSE)
    }
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

runMarker_single_algorithm_no_export = function (algorithm_name = "CS", moic.res = NULL, dea.method = c("deseq2", "edger", 
                                                                                                        "limma"), prefix = NULL, dat.path = getwd(), 
                                                 p.cutoff = 0.05, p.adj.cutoff = 0.05, dirct = "up", n.marker = 200, 
                                                 doplot = TRUE, norm.expr = NULL, annCol = NULL, annColors = NULL, 
                                                 clust.col = c("#2EC4B6", "#E71D36", "#FF9F1C", "#BDD5EA", 
                                                               "#FFA5AB", "#011627", "#023E8A", "#9D4EDD", "#f09c6c", "#09f3b3"), halfwidth = 3, 
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
  
  genelist <- c()
  for (filek in DEfiles) {
    DEres <- read.table(file.path(dat.path, filek), header = TRUE, 
                        row.names = NULL, sep = "\t", quote = "", stringsAsFactors = FALSE)
    
    # Check if the DEres is empty or has only one gene
    if (nrow(DEres) < 2) {
      stop(paste("Skipping file", filek, "because it has less than two rows of data."))
    }
    
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
    
    if (nrow(DEres) < 2) {
      stop(paste("Skipping file", filek, "because it has less than two rows of data."))
    }
    
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
  }
  
  templates <- NULL
  for (filek in DEfiles) {
    tmp <- NULL
    if (!is.null(marker[[filek]]) && nrow(marker[[filek]]) > 0) {
      tmp <- data.frame(probe = rownames(marker[[filek]]), 
                        class = sub("_vs_Others.txt", "", sub(".*.result.", "", filek)), 
                        dirct = marker$dirct, stringsAsFactors = FALSE)
    }
    if (!is.null(tmp)) {
      templates <- rbind.data.frame(templates, tmp, stringsAsFactors = FALSE)
    }
  }
  
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
    
    return(list(unqlist = unqlist, templates = templates, 
                dirct = dirct, plotdata = plotdata, annCol = annCol, annColors = annColors))
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
                        "#FFA5AB", "#011627", "#023E8A", "#9D4EDD", "#f09c6c", "#09f3b3"), color = NULL, 
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

# runGSEA for newer GSVA versions #####
runGSEA_mod_4.4 <- function (moic.res = NULL, dea.method = c("deseq2", "edger", 
                                                             "limma"), norm.expr = NULL, prefix = NULL, dat.path = getwd(), 
                             res.path = getwd(), dirct = "up", n.path = 10, msigdb.path = NULL, 
                             nPerm = 1000, minGSSize = 10, maxGSSize = 500, p.cutoff = 0.05, 
                             p.adj.cutoff = 0.05, gsva.method = "gsva", norm.method = "mean", 
                             clust.col = c("#2EC4B6", "#E71D36", "#FF9F1C", "#BDD5EA", 
                                           "#FFA5AB", "#011627", "#023E8A", "#9D4EDD", "#f09c6c", "#09f3b3"), color = NULL, 
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

# runGSEA for newer GSVA versions #####
runGSEA_mod_4.4_single_algorithm <- function (algorithm_name = "CS", moic.res = NULL, dea.method = c("deseq2", "edger", 
                                                             "limma"), norm.expr = NULL, prefix = NULL, dat.path = getwd(), 
                             res.path = getwd(), dirct = "up", n.path = 10, msigdb.path = NULL, 
                             nPerm = 1000, minGSSize = 10, maxGSSize = 500, p.cutoff = 0.05, 
                             p.adj.cutoff = 0.05, gsva.method = "gsva", norm.method = "mean", 
                             clust.col = c("#2EC4B6", "#E71D36", "#FF9F1C", "#BDD5EA", 
                                           "#FFA5AB", "#011627", "#023E8A", "#9D4EDD", "#f09c6c", "#09f3b3"), color = NULL, 
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
  if (length(pathcore) > 0) {
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
              grouped.es = esm, heatmap = hm)) } else {
                # Stop the function and print a message that no deregulated pathways were foundmessage("No deregulated pathways were found in the GSEA results.")
                return(NULL)
              }
}


# Redefine compAgree for subtype comparisons across different classifications #####

# Modification:
# Added three library calls in the beginning of the function definition
compAgree2 = function (moic.res = NULL, subt2comp = NULL, doPlot = TRUE, 
                       clust.col = c("#2EC4B6", "#E71D36", "#FF9F1C", "#BDD5EA", 
                                     "#FFA5AB", "#011627", "#023E8A", "#9D4EDD", "#f09c6c", "#09f3b3"), box.width = 0.1, 
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

# compClinvar single algorithm #####
compClinvar_single_algorithm_wordOnly = function (algorithm_name = "CS",
                                         moic.res = NULL, var2comp = NULL, strata = NULL, factorVars = NULL, 
                         nonnormalVars = NULL, exactVars = NULL, includeNA = FALSE, 
                         doWord = TRUE, tab.name = NULL, res.path = getwd(), ...) 
{
  dat <- moic.res$clust.res
  colnames(dat)[which(colnames(dat) == "clust")] <- "Subtype"
  dat$Subtype <- paste0(algorithm_name, dat$Subtype)
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

compClinvar_single_algorithm <- function(algorithm_name = "CS",
                                         moic.res = NULL,
                                         var2comp = NULL,
                                         strata = NULL,
                                         factorVars = NULL,
                                         nonnormalVars = NULL,
                                         exactVars = NULL,
                                         includeNA = FALSE,
                                         doWord = TRUE,
                                         tab.name = NULL,
                                         res.path = getwd(),
                                         output_pdf = FALSE,  
                                         pdf_level_col_width = c("15em", "15em"),  
                                         pdf_count_col_width = "10em",  
                                         pdf_pval_col_width = "8em",  
                                         pdf_test_col_width = "8em",  
                                         pdf_tab_font_size = 8,
                                         omit_warnings = FALSE
) {
  library(knitr)
  library(kableExtra)
  library(rmarkdown)
  library(officer)  
  library(stringr)
  
  # Prepare data for the table
  dat <- moic.res$clust.res
  colnames(dat)[which(colnames(dat) == "clust")] <- "Subtype"
  dat$Subtype <- paste0(algorithm_name, dat$Subtype)
  com_sam <- intersect(dat$samID, rownames(var2comp))
  
  if (length(com_sam) == nrow(dat)) {
    message("--all samples matched.")
  } else {
    message(paste0(
      "--",
      (nrow(dat) - length(com_sam)),
      " samples mismatched from current subtypes."
    ))
  }
  
  dat <- cbind.data.frame(Subtype = dat[com_sam, "Subtype", drop = FALSE], var2comp[com_sam, , drop = FALSE])
  
  if (is.null(strata)) {
    strata <- "Subtype"
  }
  if (!is.element(strata, colnames(dat))) {
    stop("fail to find this strata in var2comp. Consider using NULL by default.")
  }
  
  # Run the table creation with error handling for warnings
  warn <- character(0)
  tryCatch({
    stabl <- jstable::CreateTableOne2(
      vars = setdiff(colnames(dat), strata),
      strata = strata,
      data = dat,
      factorVars = factorVars,
      nonnormal = nonnormalVars,
      exact = exactVars,
      includeNA = includeNA,
      showAllLevels = TRUE
    )
  }, warning = function(w) {
    warn <<- append(warn, conditionMessage(w))
  })
  
  if(!omit_warnings) {
    if (length(warn) > 0 && grepl("NA", warn, fixed = TRUE)) {
      set.seed(19991018)
      stabl <- jstable::CreateTableOne2(
        vars = setdiff(colnames(dat), strata),
        strata = strata,
        data = dat,
        factorVars = factorVars,
        nonnormal = nonnormalVars,
        exact = exactVars,
        includeNA = includeNA,
        showAllLevels = TRUE,
        argsExact = list(simulate.p.value = TRUE)
      )
    }
  } else {
    set.seed(19991018)
    stabl <- jstable::CreateTableOne2(
      vars = setdiff(colnames(dat), strata),
      strata = strata,
      data = dat,
      factorVars = factorVars,
      nonnormal = nonnormalVars,
      exact = exactVars,
      includeNA = includeNA,
      showAllLevels = TRUE,
      argsExact = list(simulate.p.value = TRUE)
    )
  }
  
  # Prepare the table
  comtable <- as.data.frame(stabl)
  comtable <- cbind.data.frame(var = rownames(stabl), comtable)
  rownames(comtable) <- NULL
  colnames(comtable)[1] <- " "
  comtable[is.na(comtable)] <- ""
  comtable <- comtable[, setdiff(colnames(comtable), "sig")]
  comtable$` ` <- gsub("_", " ", comtable$` `)
  
  if (is.null(tab.name)) {
    outFile <- "summarization_of_clinical_variables_stratified_by_current_subtype"
  } else {
    outFile <- tab.name
  }
  
  # Create the PDF output using rmarkdown
  if (output_pdf) {
    # Create a temporary Rmarkdown file
    rmd_file <- tempfile(fileext = ".Rmd")
    rmd_content <- paste0(
      "---\n",
      "title: '", outFile, "'\n",
      "output:\n  pdf_document:\n    toc: true\n    number_sections: true\n    latex_engine: xelatex\n",
      "header-includes:\n  - \\usepackage{fontspec}\n  - \\setmainfont{Arial}\n",
      "---\n\n",
      "```{r, echo=FALSE, message=FALSE, warning=FALSE}\n",
      "library(knitr)\n",
      "library(kableExtra)\n\n",
      
      "# Display table with custom column width for specific columns\n",
      "kable(comtable, align = c('l', 'c', rep('c', ncol(comtable)-2)), caption = 'Summarization of clinical variables stratified by subtype') %>%\n",
      "  kable_styling(latex_options = c('striped', 'scale_down', 'repeat_header', 'hold_position'), full_width = FALSE, font_size = ", pdf_tab_font_size, ") %>%\n",
      
      # Set width for the first two columns
      "  column_spec(1, width = '", pdf_level_col_width[1], "') %>%\n",
      "  column_spec(2, width = '", pdf_level_col_width[2], "') %>%\n",
      
      # Apply pdf_count_col_width to all columns except the last two
      "  column_spec(3:", ncol(comtable)-2, ", width = '", pdf_count_col_width, "') %>%\n",
      
      # Set width for the second-to-last column (p-value column)
      "  column_spec(", ncol(comtable)-1, ", width = '", pdf_pval_col_width, "') %>%\n",
      
      # Set width for the last column (test column)
      "  column_spec(", ncol(comtable), ", width = '", pdf_test_col_width, "')\n",
      "```\n"
    )
    
    # Write the Rmarkdown content to a file
    writeLines(rmd_content, rmd_file)
    
    # Render the PDF from Rmarkdown
    pdf_file <- file.path(res.path, paste0(outFile, ".pdf"))
    rmarkdown::render(rmd_file, output_file = pdf_file)
    
    message("PDF created at: ", pdf_file)
  }
  
  # Create a Word document if doWord is TRUE
  if (doWord) {
    my_doc <- officer::read_docx() %>%
      officer::body_add_par(value = paste0("Table: ", outFile), style = "heading 1") %>%
      officer::body_add_table(value = comtable, style = "table_template") %>%
      officer::body_add_par(value = "Note: Clinical variables stratified by current subtype.")
    
    word_file <- file.path(res.path, paste0(outFile, ".docx"))
    print(my_doc, target = word_file)
    message("Word document created at: ", word_file)
  }
  
  return(list(compTab = comtable))
}

compClinvar_ordinal_single_algorithm <- function(algorithm_name = "CS",
                                                 moic.res = NULL,
                                                 var2comp = NULL,
                                                 strata = NULL,
                                                 ordinalVars = NULL,
                                                 includeNA = FALSE,
                                                 tab.name = NULL,
                                                 res.path = getwd(),
                                                 output_pdf = TRUE,
                                                 pdf_template_loc = getwd(),  # New argument for the template location
                                                 pdf_level_col_width = c("10em", "12em"),  
                                                 pdf_count_col_width = "8em",  
                                                 pdf_pval_col_width = "5em",  
                                                 pdf_test_col_width = "5em",  
                                                 pdf_tab_font_size = 8 
) {
  library(knitr)
  library(kableExtra)
  library(rmarkdown)
  library(clinfun)
  
  # Prepare data for the table
  dat <- moic.res$clust.res
  colnames(dat)[which(colnames(dat) == "clust")] <- "Subtype"
  dat$Subtype <- paste0(algorithm_name, dat$Subtype)
  com_sam <- intersect(dat$samID, rownames(var2comp))
  
  indices <- which(colnames(var2comp) %in% ordinalVars)
  ordinalVars <- gsub("_", " ", ordinalVars)
  colnames(var2comp)[indices] <- ordinalVars
  
  if (length(com_sam) == nrow(dat)) {
    message("--all samples matched.")
  } else {
    message(paste0("--", (nrow(dat) - length(com_sam)), " samples mismatched from current subtypes."))
  }
  
  dat <- cbind.data.frame(Subtype = dat[com_sam, "Subtype", drop = FALSE], var2comp[com_sam, , drop = FALSE])
  
  if (is.null(strata)) {
    strata <- "Subtype"
  }
  if (!is.element(strata, colnames(dat))) {
    stop("fail to find this strata in var2comp. Consider using NULL by default.")
  }
  
  # Prepare unique subtypes and explicitly assign them to columns
  unique_subtypes <- sort(unique(dat$Subtype))
  
  # Convert Subtype column to an ordered factor
  dat$Subtype <- factor(dat$Subtype, levels = unique_subtypes, ordered = TRUE)
  
  # Create a results table with proper column names for subtypes
  results_table <- data.frame(Variable = character(),
                              Levels = character(),
                              stringsAsFactors = FALSE)
  
  # Add columns for each unique subtype (SNF1, SNF2, etc.)
  for (subtype in unique_subtypes) {
    results_table[[subtype]] <- character(0)  # Add empty columns for each subtype
  }
  
  # Add columns for p-value and test result
  results_table$p <- character(0)  # Use character type to handle both numbers and empty strings
  results_table$test <- character(0)
  
  # Perform JT test for each ordinal variable
  for (var in ordinalVars) {
    # Filter the dataset to exclude rows where the current variable is NA
    filtered_dat <- dat[!is.na(dat[[var]]), ]
    
    # Now calculate total counts after filtering NA for each subtype
    total_counts <- table(filtered_dat$Subtype)
    
    levels_var <- sort(unique(filtered_dat[[var]]))
    test_result <- jonckheere.test(x = as.numeric(filtered_dat[[var]]), g = filtered_dat$Subtype, nperm = 10000)
    
    # Add rows to the table for each level of the ordinal variable
    first_row <- TRUE
    for (level in levels_var) {
      row <- list(Variable = ifelse(first_row, paste0(var, " (%)"), ""),  # Add "(%)" to the variable name
                  Levels = level)
      
      # Count occurrences for each group and calculate percentages based on **filtered data**
      for (group in unique_subtypes) {
        count <- sum(filtered_dat$Subtype == group & filtered_dat[[var]] == level, na.rm = TRUE)
        total_in_group <- total_counts[group]  # Total count for this subtype in the filtered data
        percentage <- if (total_in_group > 0) (count / total_in_group) * 100 else 0  # Calculate percentage safely
        row[[group]] <- paste0(count, " (", sprintf("%.1f", percentage), ")")  # Format: count (percentage)
      }
      
      # Add p-value and test only in the first row for this variable
      row$p <- ifelse(first_row, round(test_result$p.value, 4), "")  # Use empty string instead of NA
      row$test <- ifelse(first_row, "JT", "")
      
      results_table <- rbind(results_table, as.data.frame(row, stringsAsFactors = FALSE))
      first_row <- FALSE  # After the first row, the variable name is not repeated
    }
  }
  
  # Reorder the columns to be: Variable, Levels, SNF1, SNF2, ..., p, test
  colnames(results_table) = c("Variable", "Levels", unique_subtypes, "p", "test")
  results_table <- results_table[, c("Variable", "Levels", unique_subtypes, "p", "test")]
  
  # PDF generation
  if (output_pdf) {
    # Check if the template file exists in the provided location
    rmd_template <- file.path(pdf_template_loc)
    if (!file.exists(rmd_template)) {
      stop("The Rmarkdown template was not found at the specified location: ", pdf_template_loc)
    }
    
    # Pass parameters dynamically to the Rmarkdown render function
    pdf_file <- file.path(res.path, paste0(tab.name, ".pdf"))
    rmarkdown::render(rmd_template, 
                      output_file = pdf_file,
                      params = list(
                        tab_name = tab.name,
                        results_table = results_table,
                        n_subtypes = length(unique_subtypes),
                        pdf_level_col_width = pdf_level_col_width[1],
                        pdf_count_col_width = pdf_count_col_width,
                        pdf_pval_col_width = pdf_pval_col_width,
                        pdf_test_col_width = pdf_test_col_width,
                        pdf_tab_font_size = pdf_tab_font_size
                      ))
    
    message("PDF created at: ", pdf_file)
  }
  
  return(list(compTab = results_table))
}

# compMut single algorithm #####
compMut_single_algorithm  = function (algorithm_name = "CS", moic.res = NULL, mut.matrix = NULL, freq.cutoff = 0.05, 
                                     test.method = "fisher", p.adj.method = "BH", doWord = TRUE, 
                                     doPlot = TRUE, innerclust = TRUE, res.path = getwd(), tab.name = NULL, 
                                     fig.path = getwd(), fig.name = NULL, annCol = NULL, annColors = NULL, 
                                     mut.col = "#21498D", bg.col = "#dcddde", p.cutoff = 0.05, 
                                     p.adj.cutoff = 0.05, clust.col = c("#2EC4B6", "#E71D36", 
                                                                        "#FF9F1C", "#BDD5EA", "#FFA5AB", "#011627", "#023E8A", 
                                                                        "#9D4EDD", "#f09c6c", "#09f3b3"), width = 8, height = 4,
                                     simulate.p.value = FALSE, return.binary = FALSE) 
{
  library(MOVICS)
  library(grid)
  library(ComplexHeatmap)
  if (!is.element(test.method, c("fisher", "chisq"))) {
    stop("test.method for independency can be one of fisher or chisq.\n")
  }
  if (!is.element(p.adj.method, c("holm", "hochberg", "hommel", 
                                  "bonferroni", "BH", "BY", "fdr"))) {
    stop("p.adj.method can be one of holm, hochberg, hommel, bonferroni, BH, BY, fdr.\n")
  }
  create.anntrack <- function(samples = NULL, subtype = NULL, 
                              typesToPlot = NULL) {
    if (is.null(samples)) {
      stop("samples can not be NULL!")
    }
    if (length(samples) != length(subtype)) {
      stop("samples and subtype do not have equal length!")
    }
    if (is.null(typesToPlot)) {
      typesToPlot <- levels(factor(unique(subtype)))
    }
    names(subtype) <- samples
    return(list(subtype = subtype, typesToPlot = typesToPlot))
  }
  createMutSubtype <- function(indata = NULL, samples = NULL, 
                               genename = NULL) {
    if (!is.element(genename, rownames(indata))) {
      stop(paste(genename, "not found in indata!", sep = " "))
    }
    comsam <- intersect(colnames(indata), samples)
    out <- rep("Not Available", times = length(samples))
    names(out) <- samples
    out[comsam] <- as.character(indata[genename, comsam])
    out[is.na(out)] <- "Not Available"
    out[out == "0"] <- "Normal"
    out[out == "1"] <- "Mutated"
    return(out)
  }
  clust.res <- moic.res$clust.res
  clust.res <- clust.res[order(clust.res$clust, decreasing = FALSE), 
                         , drop = FALSE]
  comsam <- intersect(clust.res$samID, colnames(mut.matrix))
  clust.res <- clust.res[comsam, , drop = FALSE]
  mut.matrix <- mut.matrix[, comsam]
  n.moic <- length(unique(clust.res$clust))
  if (length(comsam) == nrow(moic.res$clust.res)) {
    message("--all samples matched.")
  }
  else {
    message("--", paste0((nrow(moic.res$clust.res) - length(comsam)), 
                         " samples mismatched from current subtypes."))
  }
  ans <- rep(paste0(algorithm_name, 1:n.moic), as.numeric(table(clust.res$clust)))
  names(ans) <- clust.res$samID
  genelist <- rownames(mut.matrix[rowSums(mut.matrix) > freq.cutoff * 
                                    nrow(clust.res), ])
  binarymut <- as.data.frame(matrix(0, nrow = nrow(clust.res), 
                                    ncol = length(genelist)))
  rownames(binarymut) <- names(ans)
  colnames(binarymut) <- genelist
  for (k in 1:length(genelist)) {
    res <- create.anntrack(samples = names(ans), subtype = createMutSubtype(mut.matrix, 
                                                                            names(ans), genelist[k]))
    binarymut[, genelist[k]] <- res$subtype
  }
  out <- matrix(0, nrow = length(genelist), ncol = n.moic + 
                  1)
  colnames(out) <- c(levels(factor(ans)), "pvalue")
  rownames(out) <- paste(genelist, "Mutated", sep = "_")
  for (k in 1:length(genelist)) {
    genek <- genelist[k]
    x <- ans
    y <- binarymut[names(x), genek, drop = TRUE]
    y <- as.character(y)
    names(y) <- names(x)
    tmp <- setdiff(names(x), names(y)[y == "Not Available"])
    x <- x[tmp]
    y <- y[tmp]
    res <- table(y, x)
    if (!all(colnames(res) == colnames(out)[1:(ncol(out) - 
                                               1)])) {
      stop(paste("colnames mismatch for ", k, sep = ""))
    }
    pct <- paste0("(", format(round(res["Mutated", ]/as.numeric(table(clust.res$clust)) * 
                                      100, 1), digits = 3), "%)")
    freqpct <- paste(res["Mutated", ], pct, sep = " ")
    out[k, 1:(ncol(out) - 1)] <- freqpct
    if (test.method == "fisher") {
      out[k, "pvalue"] <- as.numeric(fisher.test(x, y, 
                                                 workspace = 2e+09,
                                                 simulate.p.value = simulate.p.value)$p.value)
    }
    else {
      out[k, "pvalue"] <- as.numeric(chisq.test(x, y)$p.value)
    }
  }
  out <- as.data.frame(out)
  out$pvalue <- formatC(as.numeric(as.character(out$pvalue)), 
                        format = "e", digits = 2)
  out$padj <- formatC(p.adjust(as.numeric(out$pvalue), method = p.adj.method), 
                      format = "e", digits = 2)
  if (is.null(tab.name)) {
    outFile <- "Independent test between subtype and mutation.txt"
  }
  else {
    outFile <- paste0(tab.name, ".txt")
  }
  tmb <- rowSums(mut.matrix[genelist, ])
  pct <- paste0("(", format(round(tmb/ncol(mut.matrix) * 100, 
                                  1), digits = 1), "%)")
  freqpct <- paste(tmb, pct, sep = " ")
  out <- cbind.data.frame(data.frame(`Gene (Mutated)` = gsub("_Mutated", 
                                                             "", rownames(out)), TMB = freqpct, check.names = FALSE), 
                          out)
  rownames(out) <- NULL
  write.table(out, file.path(res.path, outFile), row.names = FALSE, 
              col.names = TRUE, sep = "\t", quote = FALSE)
  if (doWord) {
    comtable <- out
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
  if (doPlot) {
    if (is.null(fig.name)) {
      outFig <- paste0("oncoprint for mutations with frequency over than ", 
                       freq.cutoff * 100, " pct.pdf")
    }
    else {
      outFig <- paste0(fig.name, ".pdf")
    }
    sam.order <- moic.res$clust.res[order(moic.res$clust.res$clust, 
                                          decreasing = FALSE), "samID"]
    colvec <- clust.col[1:length(unique(moic.res$clust.res$clust))]
    names(colvec) <- paste0(algorithm_name, sort(unique(moic.res$clust.res$clust)))
    if (!is.null(annCol) & !is.null(annColors)) {
      annCol <- annCol[sam.order, , drop = FALSE]
      annCol$Subtype <- paste0(algorithm_name, moic.res$clust.res[sam.order, 
                                                        "clust"])
      annColors[["Subtype"]] <- colvec
    }
    else {
      annCol <- data.frame(Subtype = paste0(algorithm_name, moic.res$clust.res[sam.order, 
                                                                     "clust"]), row.names = sam.order)
      annColors <- list(Subtype = colvec)
    }
    sig.mut <- as.character(out[which(as.numeric(out$pvalue) < 
                                        p.cutoff & as.numeric(out$padj) < p.adj.cutoff), 
                                "Gene (Mutated)"])
    if (length(sig.mut) == 0) {
      
    } else {
      onco_dat <- t(binarymut[rownames(annCol), sig.mut, drop = FALSE])
      onco_dat[onco_dat == "Normal"] <- ""
      onco_dat <- as.data.frame(onco_dat)
      alter_fun = list(background = function(x, y, w, h) {
        grid::grid.rect(x, y, w - unit(0.5, "mm"), h - unit(0.5, 
                                                            "mm"), gp = gpar(fill = bg.col, col = NA))
      }, Mutated = function(x, y, w, h) {
        grid::grid.rect(x, y, w - unit(0.5, "mm"), h - unit(0.5, 
                                                            "mm"), gp = gpar(fill = mut.col, col = NA))
      })
      col = c(Mutated = mut.col)
      if (innerclust) {
        sam.reorder <- c()
        for (i in 1:n.moic) {
          sam <- moic.res$clust.res[which(moic.res$clust.res$clust == 
                                            i), "samID"]
          tmp <- MOVICS:::quiet(ComplexHeatmap::oncoPrint(onco_dat[, 
                                                                   sam], get_type = function(x) x, alter_fun = alter_fun, 
                                                          col = col, remove_empty_columns = FALSE, show_pct = FALSE, 
                                                          bottom_annotation = NULL, top_annotation = NULL, 
                                                          show_heatmap_legend = FALSE))
          sam.reorder <- c(sam.reorder, sam[tmp@column_order])
        }
        my_annotation = ComplexHeatmap::HeatmapAnnotation(df = annCol[sam.reorder, 
                                                                      , drop = FALSE], col = annColors)
        p <- MOVICS:::quiet(ComplexHeatmap::oncoPrint(onco_dat[, 
                                                               sam.reorder], get_type = function(x) x, alter_fun = alter_fun, 
                                                      col = col, remove_empty_columns = FALSE, column_order = sam.reorder, 
                                                      show_pct = TRUE, bottom_annotation = my_annotation, 
                                                      top_annotation = NULL, show_heatmap_legend = FALSE))
      }
      else {
        my_annotation = ComplexHeatmap::HeatmapAnnotation(df = annCol, 
                                                          col = annColors)
        p <- MOVICS:::quiet(ComplexHeatmap::oncoPrint(onco_dat, get_type = function(x) x, 
                                                      alter_fun = alter_fun, col = col, remove_empty_columns = FALSE, 
                                                      column_order = colnames(onco_dat), show_pct = TRUE, 
                                                      bottom_annotation = my_annotation, top_annotation = NULL, 
                                                      show_heatmap_legend = FALSE))
      }
      pdf(file.path(fig.path, outFig), width = width, height = height)
      draw(p)
      invisible(dev.off())
      draw(p)
    } 
    }
  if (return.binary) {
    binarymut$Sample.ID <- rownames(binarymut)   # keeps the IDs in a column
    binarymut <- binarymut[, c("Sample.ID", setdiff(names(binarymut), "Sample.ID"))]
    return(list(summary = out, sample_binary = binarymut))
  } else {
    return(out)
  }
}

# Special case for ER
compMut_single_algorithm_sc  <- function(
    algorithm_name = "CS",
    moic.res      = NULL,
    mut.matrix    = NULL,
    freq.cutoff   = 0.05,
    test.method   = "fisher",
    p.adj.method  = "BH",
    doWord        = TRUE,
    doPlot        = TRUE,
    innerclust    = TRUE,
    res.path      = getwd(),
    tab.name      = NULL,
    fig.path      = getwd(),
    fig.name      = NULL,
    annCol        = NULL,
    annColors     = NULL,
    mut.col       = "#21498D",
    bg.col        = "#dcddde",
    p.cutoff      = 0.05,
    p.adj.cutoff  = 0.05,
    clust.col     = c("#2EC4B6", "#E71D36", "#FF9F1C", "#BDD5EA", "#FFA5AB",
                      "#011627", "#023E8A", "#9D4EDD", "#f09c6c", "#09f3b3"),
    width         = 8,
    height        = 4,
    simulate.p.value = FALSE,
    return.binary = FALSE
) {
  library(MOVICS)
  library(grid)
  library(ComplexHeatmap)
  
  if (!is.element(test.method, c("fisher", "chisq"))) {
    stop("test.method for independency can be one of fisher or chisq.\n")
  }
  if (!is.element(p.adj.method, c("holm", "hochberg", "hommel",
                                  "bonferroni", "BH", "BY", "fdr"))) {
    stop("p.adj.method can be one of holm, hochberg, hommel, bonferroni, BH, BY, fdr.\n")
  }
  
  create.anntrack <- function(samples = NULL, subtype = NULL, typesToPlot = NULL) {
    if (is.null(samples)) stop("samples can not be NULL!")
    if (length(samples) != length(subtype)) {
      stop("samples and subtype do not have equal length!")
    }
    if (is.null(typesToPlot)) {
      typesToPlot <- levels(factor(unique(subtype)))
    }
    names(subtype) <- samples
    return(list(subtype = subtype, typesToPlot = typesToPlot))
  }
  
  createMutSubtype <- function(indata = NULL, samples = NULL, genename = NULL) {
    if (!is.element(genename, rownames(indata))) {
      stop(paste(genename, "not found in indata!", sep = " "))
    }
    comsam <- intersect(colnames(indata), samples)
    out <- rep("Not Available", times = length(samples))
    names(out) <- samples
    out[comsam] <- as.character(indata[genename, comsam])
    out[is.na(out)] <- "Not Available"
    out[out == "0"] <- "Normal"
    out[out == "1"] <- "Mutated"
    return(out)
  }
  
  ## --- Prepare clustering / mutation data -----------------------------------
  
  clust.res <- moic.res$clust.res
  clust.res <- clust.res[order(clust.res$clust, decreasing = FALSE), , drop = FALSE]
  
  comsam <- intersect(clust.res$samID, colnames(mut.matrix))
  clust.res <- clust.res[comsam, , drop = FALSE]
  mut.matrix <- mut.matrix[, comsam, drop = FALSE]
  
  if (length(comsam) == nrow(moic.res$clust.res)) {
    message("--all samples matched.")
  } else {
    message("--", paste0(nrow(moic.res$clust.res) - length(comsam),
                         " samples mismatched from current subtypes."))
  }
  
  ## robust subtype labelling (no rep(), no assumptions about raw clust values)
  clust_fac      <- factor(clust.res$clust)              # drop unused, sort
  n.moic         <- nlevels(clust_fac)
  subtype_names  <- paste0(algorithm_name, seq_len(n.moic))    # ER1, ER2, ...
  names(subtype_names) <- levels(clust_fac)                     # map raw -> ER1/ER2
  
  ans <- subtype_names[as.character(clust_fac)]          # one per sample
  names(ans) <- clust.res$samID
  
  ## --- Filter genes by frequency --------------------------------------------
  
  genelist <- rownames(mut.matrix[rowSums(mut.matrix) > freq.cutoff * nrow(clust.res), , drop = FALSE])
  
  if (length(genelist) == 0) {
    warning("No genes pass freq.cutoff; nothing to test.")
    if (return.binary) return(list(summary = NULL, sample_binary = NULL))
    return(NULL)
  }
  
  ## --- Build mutation status matrix (Normal / Mutated / Not Available) ------
  
  binarymut <- as.data.frame(matrix(0, nrow = nrow(clust.res),
                                    ncol = length(genelist)))
  rownames(binarymut) <- names(ans)  # sample IDs
  colnames(binarymut) <- genelist
  
  for (k in seq_along(genelist)) {
    res_k <- create.anntrack(
      samples = names(ans),
      subtype = createMutSubtype(mut.matrix, names(ans), genelist[k])
    )
    binarymut[, genelist[k]] <- res_k$subtype
  }
  
  ## --- Test association + compute per-subtype mutation frequencies ----------
  
  out <- matrix(0, nrow = length(genelist), ncol = n.moic + 1)
  colnames(out) <- c(subtype_names, "pvalue")
  rownames(out) <- paste(genelist, "Mutated", sep = "_")
  
  for (k in seq_along(genelist)) {
    genek <- genelist[k]
    
    x <- ans                                        # subtype label per sample
    y <- binarymut[names(x), genek, drop = TRUE]   # mutation status per sample
    
    y <- as.character(y)
    names(y) <- names(x)
    
    ## drop NA / "Not Available"
    tmp <- setdiff(names(x), names(y)[y == "Not Available"])
    x <- x[tmp]
    y <- y[tmp]
    
    res <- table(y, x)   # rows: Normal/Mutated; cols: ER1, ER2, ...
    
    if (!all(colnames(res) == colnames(out)[1:(ncol(out) - 1)])) {
      stop(paste("colnames mismatch for ", genek, " (index ", k, ")", sep = ""))
    }
    
    ## make sure "Mutated" row exists
    if (!("Mutated" %in% rownames(res))) {
      mut_counts <- rep(0, ncol(res))
      names(mut_counts) <- colnames(res)
    } else {
      mut_counts <- res["Mutated", ]
    }
    
    denom <- colSums(res)   # total samples per subtype for this gene
    pct   <- paste0(
      "(",
      format(round(mut_counts / denom * 100, 1), digits = 3),
      "%)"
    )
    freqpct <- paste(mut_counts, pct, sep = " ")
    
    out[k, 1:(ncol(out) - 1)] <- freqpct
    
    if (test.method == "fisher") {
      out[k, "pvalue"] <- as.numeric(fisher.test(x, y,
                                                 workspace = 2e+09,
                                                 simulate.p.value = simulate.p.value)$p.value)
    } else {
      out[k, "pvalue"] <- as.numeric(chisq.test(x, y)$p.value)
    }
  }
  
  ## --- Format p-values, add TMB columns ------------------------------------
  
  out <- as.data.frame(out, stringsAsFactors = FALSE)
  out$pvalue <- formatC(as.numeric(as.character(out$pvalue)),
                        format = "e", digits = 2)
  out$padj <- formatC(p.adjust(as.numeric(out$pvalue), method = p.adj.method),
                      format = "e", digits = 2)
  
  if (is.null(tab.name)) {
    outFile <- "Independent test between subtype and mutation.txt"
  } else {
    outFile <- paste0(tab.name, ".txt")
  }
  
  tmb <- rowSums(mut.matrix[genelist, , drop = FALSE])
  tmb_pct <- paste0(
    "(",
    format(round(tmb / ncol(mut.matrix) * 100, 1), digits = 1),
    "%)"
  )
  tmb_freqpct <- paste(tmb, tmb_pct, sep = " ")
  
  out <- cbind.data.frame(
    data.frame(
      `Gene (Mutated)` = gsub("_Mutated", "", rownames(out)),
      TMB              = tmb_freqpct,
      check.names      = FALSE
    ),
    out,
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL
  
  write.table(out, file.path(res.path, outFile),
              row.names = FALSE, col.names = TRUE,
              sep = "\t", quote = FALSE)
  
  ## --- Word export ----------------------------------------------------------
  
  if (doWord) {
    comtable   <- out
    title_name <- paste0("Table *. ", gsub(".txt", "", outFile, fixed = TRUE))
    mynote     <- "Note: ..."
    
    my_doc <- officer::read_docx()
    my_doc %>%
      officer::body_add_par(value = title_name, style = "table title") %>%
      officer::body_add_table(value = comtable, style = "table_template") %>%
      officer::body_add_par(value = mynote) %>%
      print(target = file.path(
        res.path,
        paste0("TABLE ", gsub(".txt", "", outFile, fixed = TRUE), ".docx")
      ))
  }
  
  ## --- OncoPrint ------------------------------------------------------------
  
  if (doPlot) {
    if (is.null(fig.name)) {
      outFig <- paste0("oncoprint for mutations with frequency over than ",
                       freq.cutoff * 100, " pct.pdf")
    } else {
      outFig <- paste0(fig.name, ".pdf")
    }
    
    sam.order <- moic.res$clust.res[order(moic.res$clust.res$clust,
                                          decreasing = FALSE), "samID"]
    
    colvec <- clust.col[1:length(unique(moic.res$clust.res$clust))]
    names(colvec) <- paste0(algorithm_name,
                            sort(unique(moic.res$clust.res$clust)))
    
    ## ensure annotation data frame has Subtype and is ordered
    annCol <- if (!is.null(annCol)) annCol[sam.order, , drop = FALSE] else
      data.frame(row.names = sam.order)
    
    annCol$Subtype <- paste0(algorithm_name, moic.res$clust.res[sam.order, "clust"])
    
    ## build a clean color mapping only for columns present in annCol
    annColors_use <- list()
    for (nm in colnames(annCol)) {
      if (nm == "Subtype") {
        sub_levels <- sort(unique(annCol$Subtype))
        cols <- clust.col[seq_along(sub_levels)]
        names(cols) <- sub_levels             # VERY IMPORTANT: names = levels
        annColors_use[[nm]] <- cols
      } else if (!is.null(annColors) && nm %in% names(annColors)) {
        annColors_use[[nm]] <- annColors[[nm]]
      }
    }
    
    ## sanity check: all entries must be named vectors
    bad <- vapply(annColors_use, function(v) is.null(names(v)), logical(1))
    if (any(bad)) {
      stop("These entries in annColors_use are not named vectors: ",
           paste(names(annColors_use)[bad], collapse = ", "))
    }
    
    sig.mut <- as.character(out[
      which(as.numeric(out$pvalue) < p.cutoff &
              as.numeric(out$padj)   < p.adj.cutoff),
      "Gene (Mutated)"
    ])
    
    if (length(sig.mut) == 0) {
      message("No mutations pass p < ", p.cutoff,
              " and padj < ", p.adj.cutoff,
              " — skipping OncoPrint.")
    } else {
      onco_dat <- t(binarymut[rownames(annCol), sig.mut, drop = FALSE])
      onco_dat[onco_dat == "Normal"] <- ""
      onco_dat <- as.data.frame(onco_dat)
      
      alter_fun <- list(
        background = function(x, y, w, h) {
          grid::grid.rect(x, y, w - unit(0.5, "mm"), h - unit(0.5, "mm"),
                          gp = gpar(fill = bg.col, col = NA))
        },
        Mutated = function(x, y, w, h) {
          grid::grid.rect(x, y, w - unit(0.5, "mm"), h - unit(0.5, "mm"),
                          gp = gpar(fill = mut.col, col = NA))
        }
      )
      col <- c(Mutated = mut.col)
      
      if (innerclust) {
        sam.reorder <- c()
        for (i in seq_len(n.moic)) {
          sam <- moic.res$clust.res[which(moic.res$clust.res$clust == i), "samID"]
          tmp_onco <- MOVICS:::quiet(
            ComplexHeatmap::oncoPrint(
              onco_dat[, sam, drop = FALSE],
              get_type = function(x) x,
              alter_fun = alter_fun,
              col = col,
              remove_empty_columns = FALSE,
              show_pct = FALSE,
              bottom_annotation = NULL,
              top_annotation = NULL,
              show_heatmap_legend = FALSE
            )
          )
          sam.reorder <- c(sam.reorder, sam[tmp_onco@column_order])
        }
        my_annotation <- ComplexHeatmap::HeatmapAnnotation(
          df  = annCol[sam.reorder, , drop = FALSE],
          col = annColors_use
        )
        p <- MOVICS:::quiet(
          ComplexHeatmap::oncoPrint(
            onco_dat[, sam.reorder, drop = FALSE],
            get_type = function(x) x,
            alter_fun = alter_fun,
            col = col,
            remove_empty_columns = FALSE,
            column_order = sam.reorder,
            show_pct = TRUE,
            bottom_annotation = my_annotation,
            top_annotation = NULL,
            show_heatmap_legend = FALSE
          )
        )
      } else {
        my_annotation <- ComplexHeatmap::HeatmapAnnotation(
          df  = annCol,
          col = annColors
        )
        p <- MOVICS:::quiet(
          ComplexHeatmap::oncoPrint(
            onco_dat,
            get_type = function(x) x,
            alter_fun = alter_fun,
            col = col,
            remove_empty_columns = FALSE,
            column_order = colnames(onco_dat),
            show_pct = TRUE,
            bottom_annotation = my_annotation,
            top_annotation = NULL,
            show_heatmap_legend = FALSE
          )
        )
      }
      
      pdf(file.path(fig.path, outFig), width = width, height = height)
      draw(p)
      invisible(dev.off())
      draw(p)
    }
  }
  
  if (return.binary) {
    binarymut$Sample.ID <- rownames(binarymut)
    binarymut <- binarymut[, c("Sample.ID", setdiff(names(binarymut), "Sample.ID"))]
    return(list(summary = out, sample_binary = binarymut))
  } else {
    return(out)
  }
}


# compDrugsen single algorithm #####
compDrugsen_single_algorithm = function (algorithm_name = "CS", moic.res = NULL, norm.expr = NULL, drugs = c("Cisplatin", 
                                                                                      "Paclitaxel"), tissueType = "all", test.method = "nonparametric", 
                                         clust.col = c("#2EC4B6", "#E71D36", "#FF9F1C", "#BDD5EA", 
                                                       "#FFA5AB", "#011627", "#023E8A", "#9D4EDD", "#f09c6c", "#09f3b3"), prefix = NULL, 
                                         seed = 123456, fig.path = getwd(), width = 5, height = 5, notch = TRUE) 
{
  library(MOVICS)
  library(ggplot2)
  if (!is.element(test.method, c("nonparametric", "parametric"))) {
    stop("test.method can be one of nonparametric or parametric.")
  }
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
  sam.order <- moic.res$clust.res[order(moic.res$clust.res$clust, 
                                        decreasing = FALSE), "samID"]
  colvec <- clust.col[1:length(unique(moic.res$clust.res$clust))]
  names(colvec) <- paste0(algorithm_name, sort(unique(moic.res$clust.res$clust)))
  annCol <- data.frame(Subtype = paste0(algorithm_name, moic.res$clust.res[sam.order, 
                                                                 "clust"]), samID = sam.order, row.names = sam.order, 
                       stringsAsFactors = FALSE)
  if (max(norm.expr) < 25 | (max(norm.expr) >= 25 & min(norm.expr) < 
                             0)) {
    message("--expression profile seems to have veen standardised (z-score or log transformation), no more action will be performed.")
    gset <- norm.expr
  }
  if (max(norm.expr) >= 25 & min(norm.expr) >= 0) {
    message("--log2 transformation done for expression data.")
    gset <- log2(norm.expr + 1)
  }
  predictedPtype <- predictedBoxdat <- list()
  for (drug in drugs) {
    set.seed(seed)
    predictedPtype[[drug]] <- MOVICS:::quiet(MOVICS:::pRRopheticPredict(testMatrix = as.matrix(gset[, 
                                                                                  rownames(annCol)]), drug = drug, tissueType = tissueType, 
                                                      dataset = "cgp2016", minNumSamples = 5, selection = 1))
    if (!all(names(predictedPtype[[drug]]) == rownames(annCol))) {
      stop("name mismatched!\n")
    }
    predictedBoxdat[[drug]] <- data.frame(Est.IC50 = predictedPtype[[drug]], 
                                          Subtype = as.character(annCol$Subtype), row.names = names(predictedPtype[[drug]]), 
                                          stringsAsFactors = FALSE)
    message(drug, " done...")
    if (n.moic == 2 & test.method == "nonparametric") {
      statistic = "wilcox.test"
      ic50.test <- wilcox.test(predictedBoxdat[[drug]]$Est.IC50 ~ 
                                 predictedBoxdat[[drug]]$Subtype)$p.value
      cat(paste0("Wilcoxon rank sum test p value = ", 
                 formatC(ic50.test, format = "e", digits = 2), 
                 " for ", drug))
    }
    if (n.moic == 2 & test.method == "parametric") {
      statistic = "t.test"
      ic50.test <- t.test(predictedBoxdat[[drug]]$Est.IC50 ~ 
                            predictedBoxdat[[drug]]$Subtype)$p.value
      cat(paste0("Student's t test p value = ", formatC(ic50.test, 
                                                        format = "e", digits = 2), " for ", drug))
    }
    if (n.moic > 2 & test.method == "nonparametric") {
      statistic = "kruskal.test"
      ic50.test <- kruskal.test(predictedBoxdat[[drug]]$Est.IC50 ~ 
                                  predictedBoxdat[[drug]]$Subtype)$p.value
      pairwise.ic50.test <- pairwise.wilcox.test(predictedBoxdat[[drug]]$Est.IC50, 
                                                 predictedBoxdat[[drug]]$Subtype, p.adjust.method = "BH")
      cat(paste0(drug, ": Kruskal-Wallis rank sum test p value = ", 
                 formatC(ic50.test, format = "e", digits = 2), 
                 "\npost-hoc pairwise wilcoxon rank sum test with Benjamini-Hochberg adjustment presents below:\n"))
      print(formatC(pairwise.ic50.test$p.value, format = "e", 
                    digits = 2))
    }
    if (n.moic > 2 & test.method == "parametric") {
      statistic = "anova"
      ic50.test <- summary(aov(predictedBoxdat[[drug]]$Est.IC50 ~ 
                                 predictedBoxdat[[drug]]$Subtype))[[1]][["Pr(>F)"]][1]
      pairwise.ic50.test <- pairwise.t.test(predictedBoxdat[[drug]]$Est.IC50, 
                                            predictedBoxdat[[drug]]$Subtype, p.adjust.method = "BH")
      cat(paste0(drug, ": One-way anova test p value = ", 
                 formatC(ic50.test, format = "e", digits = 2), 
                 "\npost-hoc pairwise Student's t test with Benjamini-Hochberg adjustment presents below:\n"))
      print(formatC(pairwise.ic50.test$p.value, format = "e", 
                    digits = 2))
    }
    p <- ggplot(data = predictedBoxdat[[drug]], aes(x = Subtype, 
                                                    y = Est.IC50, fill = Subtype)) + scale_fill_manual(values = colvec) + 
      geom_violin(alpha = 0.4, position = position_dodge(width = 0.75), 
                  size = 0.8, color = "black") + geom_boxplot(notch = notch, 
                                                              outlier.size = -1, color = "black", lwd = 0.8, alpha = 0.7) + 
      geom_point(shape = 21, size = 2, position = position_jitterdodge(), 
                 color = "black", alpha = 1) + theme_classic() + 
      ylab(bquote("Estimated IC"[50] ~ "of" ~ .(drug))) + 
      xlab("") + theme(axis.text.x = element_text(angle = 45, 
                                                  hjust = 1, size = 12), axis.ticks = element_line(linewidth = 0.2, 
                                                                                                   color = "black"), axis.ticks.length = unit(0.2, 
                                                                                                                                              "cm"), legend.position = "none", axis.title = element_text(size = 15), 
                       axis.text = element_text(size = 10)) + ggpubr::stat_compare_means(method = statistic, 
                                                                                 hjust = ifelse(n.moic%%2 == 0, 0.5, 0), label.x = ifelse(n.moic%%2 == 
                                                                                                                                            0, n.moic/2 + 0.5, n.moic/2), label.y = min(predictedBoxdat[[drug]]$Est.IC50))
    if (is.null(prefix)) {
      outFig <- paste0("boxviolin of estimated ic50 for ", 
                       drug, ".pdf")
    }
    else {
      outFig <- paste0(prefix, " for ", drug, ".pdf")
    }
    ggsave(file.path(fig.path, outFig), width = width, height = height)
    print(p)
  }
  return(predictedBoxdat)
}

# compAgree single algorithm #####
compAgree_single_algorithm = function (algorithm_name = "CS",
                                       moic.res = NULL, subt2comp = NULL, doPlot = TRUE, 
                       clust.col = c("#2EC4B6", "#E71D36", "#FF9F1C", "#BDD5EA", 
                                     "#FFA5AB", "#011627", "#023E8A", "#9D4EDD", "#f09c6c", "#09f3b3"), box.width = 0.1, 
                       fig.name = NULL, fig.path = getwd(), width = 6, height = 5) 
{
  library(ggplot2)
  library(ggalluvial)
  library(cowplot)
  dat <- moic.res$clust.res
  colnames(dat)[which(colnames(dat) == "clust")] <- "Subtype"
  dat$Subtype <- paste0(algorithm_name, dat$Subtype)
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
# runGSVA_mod_4.4 #####
runGSVA_mod_4.4 <- function (moic.res = NULL, norm.expr = NULL, gset.gmt.path = NULL, 
          gsva.method = "gsva", centerFlag = TRUE, scaleFlag = TRUE, 
          halfwidth = 1, annCol = NULL, annColors = NULL, clust.col = c("#2EC4B6", 
                                                                        "#E71D36", "#FF9F1C", "#BDD5EA", "#FFA5AB", "#011627", 
                                                                        "#023E8A", "#9D4EDD", "#f09c6c", "#09f3b3"), distance = "euclidean", linkage = "ward.D", 
          show_rownames = TRUE, show_colnames = FALSE, color = c("#366A9B", 
                                                                 "#4E98DE", "#DDDDDD", "#FBCFA7", "#F79C4A"), fig.path = getwd(), 
          fig.name = NULL, width = 8, height = 8, ...) 
{
  standarize.fun <- function(indata = NULL, halfwidth = NULL, 
                             centerFlag = TRUE, scaleFlag = TRUE) {
    outdata = t(scale(t(indata), center = centerFlag, scale = scaleFlag))
    if (!is.null(halfwidth)) {
      outdata[outdata > halfwidth] = halfwidth
      outdata[outdata < (-halfwidth)] = -halfwidth
    }
    return(outdata)
  }
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
  gset <- try(clusterProfiler::read.gmt(gset.gmt.path), silent = TRUE)
  if (class(gset) == "try-error") {
    stop("please provide correct ABSOLUTE PATH for gene sets of interest.")
  }
  term <- unique(gset[, 1])
  gset.list <- list()
  for (i in term) {
    gset.list[[i]] <- gset[which(gset[, 1] == i), 2]
  }
  if (max(norm.expr) < 25 | (max(norm.expr) >= 25 & min(norm.expr) < 
                             0)) {
    message("--expression profile seems to have been standardised (z-score or log transformation), no more action will be performed.")
  }
  if (max(norm.expr) >= 25 & min(norm.expr) >= 0) {
    message("--log2 transformation done for expression data.")
    norm.expr <- log2(norm.expr + 1)
  }
  es <- GSVA::gsva(param = GSVA::gsvaParam(exprData = as.matrix(norm.expr),
                                     geneSets = gset.list
  ))
  es.backup <- es
  es <- standarize.fun(es, halfwidth = halfwidth, centerFlag = centerFlag, 
                       scaleFlag = scaleFlag)
  message(gsva.method, " done...")
  if (is.null(fig.name)) {
    outFig <- paste0("enrichment_heatmap_using_", gsva.method, 
                     ".pdf")
  }
  else {
    outFig <- paste0(fig.name, "_", gsva.method, ".pdf")
  }
  sam.order <- moic.res$clust.res[order(moic.res$clust.res$clust, 
                                        decreasing = FALSE), "samID"]
  colvec <- clust.col[1:n.moic]
  names(colvec) <- paste0("CS", 1:n.moic)
  if (!is.null(annCol) & !is.null(annColors)) {
    annCol <- annCol[sam.order, , drop = FALSE]
    annCol$Subtype <- paste0("CS", moic.res$clust.res[sam.order, 
                                                      "clust"])
    annColors[["Subtype"]] <- colvec
  }
  else {
    annCol <- data.frame(Subtype = paste0("CS", moic.res$clust.res[sam.order, 
                                                                   "clust"]), row.names = sam.order, stringsAsFactors = FALSE)
    annColors <- list(Subtype = colvec)
  }
  if (!is.null(annCol) & !is.null(annColors)) {
    for (i in names(annColors)) {
      if (is.function(annColors[[i]])) {
        annColors[[i]] <- annColors[[i]](pretty(range(annCol[, 
                                                             i]), n = 64))
      }
    }
  }
  ht_opt$message = FALSE
  if (is.null(distance) | is.null(linkage)) {
    hcg <- FALSE
  }
  else {
    hcg <- fastcluster::hclust(ClassDiscovery::distanceMatrix(t(as.matrix(es[, 
                                                                sam.order])), distance), linkage)
  }
  hm <- ComplexHeatmap::pheatmap(mat = es[, sam.order], border_color = NA, 
                                 cluster_cols = FALSE, cluster_rows = hcg, annotation_col = annCol, 
                                 annotation_colors = annColors, show_rownames = show_rownames, 
                                 show_colnames = show_colnames, color = (grDevices::colorRampPalette(color))(64), 
                                 ...)
  pdf(file.path(fig.path, outFig), width = width, height = height)
  draw(hm)
  invisible(dev.off())
  draw(hm)
  return(list(gset.list = gset.list, raw.es = es.backup, scaled.es = es))
}

# runGSVA_mod_4.4 single algorithm #####
runGSVA_mod_4.4_single_algorithm <- function (algorithm_name = "CS",
                             moic.res = NULL, norm.expr = NULL, gset.gmt.path = NULL, 
                             gsva.method = "gsva", centerFlag = TRUE, scaleFlag = TRUE, 
                             halfwidth = 1, annCol = NULL, annColors = NULL, clust.col = c("#2EC4B6", 
                                                                                           "#E71D36", "#FF9F1C", "#BDD5EA", "#FFA5AB", "#011627", 
                                                                                           "#023E8A", "#9D4EDD", "#f09c6c", "#09f3b3"), distance = "euclidean", linkage = "ward.D", 
                             show_rownames = TRUE, show_colnames = FALSE, color = c("#366A9B", 
                                                                                    "#4E98DE", "#DDDDDD", "#FBCFA7", "#F79C4A"), fig.path = getwd(), 
                             fig.name = NULL, width = 8, height = 8, ...) 
{
  standarize.fun <- function(indata = NULL, halfwidth = NULL, 
                             centerFlag = TRUE, scaleFlag = TRUE) {
    outdata = t(scale(t(indata), center = centerFlag, scale = scaleFlag))
    if (!is.null(halfwidth)) {
      outdata[outdata > halfwidth] = halfwidth
      outdata[outdata < (-halfwidth)] = -halfwidth
    }
    return(outdata)
  }
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
  gset <- try(clusterProfiler::read.gmt(gset.gmt.path), silent = TRUE)
  if (class(gset) == "try-error") {
    stop("please provide correct ABSOLUTE PATH for gene sets of interest.")
  }
  term <- unique(gset[, 1])
  gset.list <- list()
  for (i in term) {
    gset.list[[i]] <- gset[which(gset[, 1] == i), 2]
  }
  if (max(norm.expr) < 25 | (max(norm.expr) >= 25 & min(norm.expr) < 
                             0)) {
    message("--expression profile seems to have been standardised (z-score or log transformation), no more action will be performed.")
  }
  if (max(norm.expr) >= 25 & min(norm.expr) >= 0) {
    message("--log2 transformation done for expression data.")
    norm.expr <- log2(norm.expr + 1)
  }
  es <- GSVA::gsva(param = GSVA::gsvaParam(exprData = as.matrix(norm.expr),
                                           geneSets = gset.list
  ))
  es.backup <- es
  es <- standarize.fun(es, halfwidth = halfwidth, centerFlag = centerFlag, 
                       scaleFlag = scaleFlag)
  message(gsva.method, " done...")
  if (is.null(fig.name)) {
    outFig <- paste0("enrichment_heatmap_using_", gsva.method, 
                     ".pdf")
  }
  else {
    outFig <- paste0(fig.name, "_", gsva.method, ".pdf")
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
  if (!is.null(annCol) & !is.null(annColors)) {
    for (i in names(annColors)) {
      if (is.function(annColors[[i]])) {
        annColors[[i]] <- annColors[[i]](pretty(range(annCol[, 
                                                             i]), n = 64))
      }
    }
  }
  ht_opt$message = FALSE
  if (is.null(distance) | is.null(linkage)) {
    hcg <- FALSE
  }
  else {
    hcg <- fastcluster::hclust(ClassDiscovery::distanceMatrix(t(as.matrix(es[, 
                                                                sam.order])), distance), linkage)
  }
  hm <- ComplexHeatmap::pheatmap(mat = es[, sam.order], border_color = NA, 
                                 cluster_cols = FALSE, cluster_rows = hcg, annotation_col = annCol, 
                                 annotation_colors = annColors, show_rownames = show_rownames, 
                                 show_colnames = show_colnames, color = (grDevices::colorRampPalette(color))(64), 
                                 ...)
  pdf(file.path(fig.path, outFig), width = width, height = height)
  draw(hm)
  invisible(dev.off())
  draw(hm)
  return(list(gset.list = gset.list, raw.es = es.backup, scaled.es = es))
}

# runMarker_mod_4.4 #####
runMarker_mod_4.4 <- function (moic.res = NULL, dea.method = c("deseq2", "edger", "limma"), 
                       prefix = NULL, dat.path = getwd(), res.path = getwd(), 
                       p.cutoff = 0.05, p.adj.cutoff = 0.05, dirct = "up", 
                       n.marker = 200, doplot = TRUE, norm.expr = NULL, 
                       annCol = NULL, annColors = NULL, clust.col = c("#2EC4B6", 
                                                                      "#E71D36", "#FF9F1C", "#BDD5EA", "#FFA5AB", "#011627", 
                                                                      "#023E8A", "#9D4EDD", "#f09c6c", "#09f3b3"), halfwidth = 3, centerFlag = TRUE, 
                       scaleFlag = TRUE, show_rownames = FALSE, show_colnames = FALSE, 
                       color = c("#5bc0eb", "black", "#ECE700"), fig.path = getwd(), 
                       fig.name = NULL, width = 8, height = 8, ...) {
  n.moic <- length(unique(moic.res$clust.res$clust))
  mo.method <- moic.res$mo.method
  DEpattern <- paste(mo.method, "_", ifelse(is.null(prefix), "", paste0(prefix, "_")), dea.method, ".*._vs_Others.txt$", sep = "")
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
    DEres <- read.table(file.path(dat.path, filek), header = TRUE, row.names = NULL, sep = "\t", quote = "", stringsAsFactors = FALSE)
    DEres <- DEres[!duplicated(DEres[, 1]), ]
    DEres <- DEres[!is.na(DEres[, 1]), ]
    rownames(DEres) <- DEres[, 1]
    DEres <- DEres[, -1]
    if (dirct == "up") {
      genelist <- c(genelist, rownames(DEres[!is.na(DEres$padj) & DEres$pvalue < p.cutoff & DEres$padj < p.adj.cutoff & !is.na(DEres$log2fc) & DEres$log2fc > 0, ]))
    }
    if (dirct == "down") {
      genelist <- c(genelist, rownames(DEres[!is.na(DEres$padj) & DEres$pvalue < p.cutoff & DEres$padj < p.adj.cutoff & !is.na(DEres$log2fc) & DEres$log2fc < 0, ]))
    }
  }
  unqlist <- setdiff(genelist, genelist[duplicated(genelist)])
  marker <- list()
  for (filek in DEfiles) {
    DEres <- read.table(file.path(dat.path, filek), header = TRUE, row.names = NULL, sep = "\t", quote = "", stringsAsFactors = FALSE)
    DEres <- DEres[!duplicated(DEres[, 1]), ]
    DEres <- DEres[!is.na(DEres[, 1]), ]
    rownames(DEres) <- DEres[, 1]
    DEres <- DEres[, -1]
    if (dirct == "up") {
      outk <- intersect(unqlist, rownames(DEres[!is.na(DEres$padj) & DEres$pvalue < p.cutoff & DEres$padj < p.adj.cutoff & !is.na(DEres$log2fc) & DEres$log2fc > 0, ]))
      outk <- DEres[outk, ]
      outk <- outk[order(outk$log2fc, decreasing = TRUE), ]
      if (nrow(outk) > n.marker) {
        marker[[filek]] <- outk[1:n.marker, ]
      } else {
        marker[[filek]] <- outk
      }
      marker$dirct <- "up"
    }
    if (dirct == "down") {
      outk <- intersect(unqlist, rownames(DEres[!is.na(DEres$padj) & DEres$pvalue < p.cutoff & DEres$padj < p.adj.cutoff & !is.na(DEres$log2fc) & DEres$log2fc < 0, ]))
      outk <- DEres[outk, ]
      outk <- outk[order(outk$log2fc, decreasing = FALSE), ]
      if (nrow(outk) > n.marker) {
        marker[[filek]] <- outk[1:n.marker, ]
      } else {
        marker[[filek]] <- outk
      }
      marker$dirct <- "down"
    }
    write.table(outk, file = file.path(res.path, paste(gsub("_vs_Others.txt", "", filek, fixed = TRUE), outlabel, sep = "_")), row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE)
  }
  templates <- NULL
  for (filek in DEfiles) {
    if (nrow(marker[[filek]]) > 0) { # Check if the marker list is not empty
      tmp <- data.frame(probe = rownames(marker[[filek]]), class = sub("_vs_Others.txt", "", sub(".*.result.", "", filek)), dirct = marker$dirct, stringsAsFactors = FALSE)
      templates <- rbind.data.frame(templates, tmp, stringsAsFactors = FALSE)
    } else {
      message(paste("No significant genes found in file:", filek))
    }
  }
  write.table(templates, file = file.path(res.path, paste0(mo.method, "_", dea.method, "_", dirct, "regulated_marker_templates.txt")), row.names = FALSE, sep = "\t", quote = FALSE)
  if (doplot) {
    if (is.null(norm.expr)) {
      stop("please provide a matrix or data.frame of normalized expression data with rows for genes and columns for samples; FPKM or TPM without log2 transformation is recommended.")
    }
    comsam <- intersect(moic.res$clust.res$samID, colnames(norm.expr))
    if (length(comsam) == nrow(moic.res$clust.res)) {
      message("--all samples matched.")
    } else {
      message(paste0("--", (nrow(moic.res$clust.res) - length(comsam)), " samples mismatched from current subtypes."))
    }
    moic.res$clust.res <- moic.res$clust.res[comsam, , drop = FALSE]
    norm.expr <- norm.expr[, comsam]
    if (is.null(fig.name)) {
      outFig <- paste0("markerheatmap_using_", dirct, "regulated_genes.pdf")
    } else {
      outFig <- paste0(fig.name, "_using_", dirct, "regulated_genes.pdf")
    }
    sam.order <- moic.res$clust.res[order(moic.res$clust.res$clust, decreasing = FALSE), "samID"]
    colvec <- clust.col[1:n.moic]
    names(colvec) <- paste0("CS", 1:n.moic)
    if (!is.null(annCol) & !is.null(annColors)) {
      annCol <- annCol[sam.order, , drop = FALSE]
      annCol$Subtype <- paste0("CS", moic.res$clust.res[sam.order, "clust"])
      annColors[["Subtype"]] <- colvec
    } else {
      annCol <- data.frame(Subtype = paste0("CS", moic.res$clust.res[sam.order, "clust"]), row.names = sam.order, stringsAsFactors = FALSE)
      annColors <- list(Subtype = colvec)
    }
    if (max(norm.expr) < 25 | (max(norm.expr) >= 25 & min(norm.expr) < 0)) {
      message("--expression profile seems to have been standardised (z-score or log transformation), no more action will be performed.")
      gset <- norm.expr
    }
    if (max(norm.expr) >= 25 & min(norm.expr) >= 0) {
      message("--log2 transformation done for expression data.")
      gset <- log2(norm.expr + 1)
    }
    standarize.fun <- function(indata = NULL, halfwidth = NULL, centerFlag = TRUE, scaleFlag = TRUE) {
      outdata = t(scale(t(indata), center = centerFlag, scale = scaleFlag))
      if (!is.null(halfwidth)) {
        outdata[outdata > halfwidth] = halfwidth
        outdata[outdata < (-halfwidth)] = -halfwidth
      }
      return(outdata)
    }
    plotdata <- standarize.fun(gset[intersect(templates$probe, rownames(gset)), sam.order], halfwidth = halfwidth, centerFlag = centerFlag, scaleFlag = scaleFlag)
    if (!is.null(annCol) & !is.null(annColors)) {
      for (i in names(annColors)) {
        if (is.function(annColors[[i]])) {
          annColors[[i]] <- annColors[[i]](pretty(range(annCol[, i]), n = 64))
        }
      }
    }
    hm <- ComplexHeatmap::pheatmap(mat = plotdata, border_color = NA, cluster_cols = FALSE, cluster_rows = FALSE, annotation_col = annCol, annotation_colors = annColors, legend_breaks = pretty(c(-halfwidth, halfwidth)), legend_labels = pretty(c(-halfwidth, halfwidth)), show_rownames = show_rownames, show_colnames = show_colnames, treeheight_col = 0, treeheight_row = 0, color = (grDevices::colorRampPalette(color))(64), ...)
    pdf(file.path(fig.path, outFig), width = width, height = height)
    draw(hm, annotation_legend_side = "left", heatmap_legend_side = "left")
    invisible(dev.off())
    draw(hm, annotation_legend_side = "left", heatmap_legend_side = "left")
    return(list(unqlist = unqlist, templates = templates, dirct = dirct, heatmap = hm))
  } else {
    return(list(unqlist = unqlist, templates = templates, dirct = dirct))
  }
}

# Modified survival function (for extend = TRUE) #####
compSurv_ext <- function (moic.res = NULL, surv.info = NULL, convt.time = "d", 
                          surv.cut = NULL, xyrs.est = NULL, clust.col = c("#2EC4B6", 
                                                                          "#E71D36", "#FF9F1C", "#BDD5EA", "#FFA5AB", "#011627", 
                                                                          "#023E8A", "#9D4EDD", "#f09c6c", "#09f3b3"), p.adjust.method = "BH", surv.median.line = "none", 
                          fig.name = NULL, fig.path = getwd()) 
{
  if (!all(is.element(c("futime", "fustat"), colnames(surv.info)))) {
    stop("fail to find variables of futime and fustat.")
  }
  if (!all(is.element(convt.time, c("d", "m", "y")))) {
    stop("unsupported time conversion. Allowed values contain c('d', 'm', 'y').")
  }
  comsam <- intersect(rownames(surv.info), rownames(moic.res$clust.res))
  mosurv.res <- cbind.data.frame(surv.info[comsam, c("futime", 
                                                     "fustat")], moic.res$clust.res[comsam, "clust", drop = FALSE])
  message(paste0("--a total of ", length(comsam), " samples are identified."))
  if (sum(c(is.na(mosurv.res$futime), is.na(mosurv.res$fustat))) > 
      0) {
    message("--removed missing values.")
    mosurv.res <- as.data.frame(na.omit(mosurv.res))
    message(paste0("--leaving ", nrow(mosurv.res), " observations."))
  }
  if (max(mosurv.res$futime) < 365) {
    warning("it seems the 'futime' might not in [day] unit, please make sure you provide the correct survival information.")
  }
  mosurv.res$Subtype <- paste0("CS", mosurv.res$clust)
  mosurv.res <- mosurv.res[order(mosurv.res$Subtype), ]
  if (!is.null(xyrs.est)) {
    if (max(xyrs.est) * 365 < max(mosurv.res$futime)) {
      xyrs <- summary(survfit(Surv(futime, fustat) ~ Subtype, 
                              data = mosurv.res), times = xyrs.est * 365,
                      extend = TRUE)
    }
    else {
      stop("the maximal survival time is less than the time point you want to estimate!")
    }
  }
  else {
    xyrs <- "[Not Available]: argument of xyrs.est was not specified."
  }
  mosurv.res$futime <- switch(convt.time, d = mosurv.res$futime, 
                              m = round(mosurv.res$futime/30.5, 4), y = round(mosurv.res$futime/365, 
                                                                              4))
  date.lab <- switch(convt.time, d = "Days", m = "Months", 
                     y = "Years")
  if (date.lab == "Days" & max(mosurv.res$futime) > 3650) {
    brk = 365
  }
  if (date.lab == "Days" & max(mosurv.res$futime) <= 3650) {
    brk = floor(max(mosurv.res$futime)/10)
  }
  if (date.lab == "Months" & max(mosurv.res$futime) > 120) {
    brk = 12
  }
  if (date.lab == "Months" & max(mosurv.res$futime) <= 120) {
    brk = floor(max(mosurv.res$futime)/10)
  }
  if (date.lab == "Years" & max(mosurv.res$futime) > 10) {
    brk = 1
  }
  if (date.lab == "Years" & max(mosurv.res$futime) <= 10) {
    brk = 1
  }
  if (is.null(surv.cut)) {
    xlim = c(0, max(mosurv.res$futime))
  }
  else {
    message(paste0("--cut survival curve up to ", surv.cut, 
                   " ", tolower(date.lab)))
    xlim = c(0, surv.cut)
  }
  n.moic <- length(unique(mosurv.res$Subtype))
  fitd <- survdiff(Surv(futime, fustat) ~ Subtype, data = mosurv.res, 
                   na.action = na.exclude)
  p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
  fit <- survfit(Surv(futime, fustat) ~ Subtype, data = mosurv.res, 
                 type = "kaplan-meier", error = "greenwood", conf.type = "plain", 
                 na.action = na.exclude, extend = TRUE)
  names(fit$strata) <- gsub("Subtype=", "", names(fit$strata))
  p <- suppressWarnings(ggsurvplot(fit = fit, conf.int = FALSE, 
                                   risk.table = TRUE, risk.table.col = "strata", palette = clust.col[1:n.moic], 
                                   data = mosurv.res, size = 1, xlim = xlim, break.time.by = brk, 
                                   legend.title = "", surv.median.line = surv.median.line, 
                                   xlab = paste0("Time (", date.lab, ")"), ylab = "Survival probability (%)", 
                                   risk.table.y.text = FALSE))
  p$plot <- suppressWarnings(p$plot + scale_y_continuous(breaks = seq(0, 
                                                           1, 0.25), labels = seq(0, 100, 25)))
  if (n.moic > 2) {
    p.lab <- paste0("Overall P", ifelse(p.val < 0.001, " < 0.001", 
                                        paste0(" = ", round(p.val, 3))))
    p$plot <- p$plot + annotate("text", x = 0, y = 0.55, 
                                hjust = 0, fontface = 4, label = p.lab)
    ps <- pairwise_survdiff(Surv(futime, fustat) ~ Subtype, 
                            data = mosurv.res, p.adjust.method = p.adjust.method)
    addTab <- as.data.frame(as.matrix(ifelse(round(ps$p.value, 
                                                   3) < 0.001, "<0.001", round(ps$p.value, 3))))
    addTab[is.na(addTab)] <- "-"
    df <- tibble(x = 0, y = 0, tb = list(addTab))
    p$plot <- p$plot + geom_table(data = df, aes(x = x, 
                                                 y = y, label = tb), table.rownames = TRUE)
  }
  else {
    p.lab <- paste0("P", ifelse(p.val < 0.001, " < 0.001", 
                                paste0(" = ", round(p.val, 3))))
    p$plot <- p$plot + annotate("text", x = 0, y = 0.55, 
                                hjust = 0, fontface = 4, label = p.lab)
  }
  if (!is.null(fig.name)) {
    outFile <- file.path(fig.path, paste0(fig.name, ".pdf"))
  }
  else {
    outFile <- file.path(fig.path, paste0("km_curve_", moic.res$mo.method, 
                                          ".pdf"))
  }
  pdf.options(reset = TRUE, onefile = FALSE)
  pdf(outFile, width = 6, height = 7)
  print(p)
  dev.off()
  print(p)
  if (n.moic > 2) {
    return(list(fitd = fitd, fit = fit, xyrs.est = xyrs, 
                overall.p = p.val, pairwise.p = ps))
  }
  else {
    return(list(fitd = fitd, fit = fit, xyrs.est = xyrs, 
                overall.p = p.val))
  }
}

# runPAM single algorithm #####
runPAM_single_algorithm = function (algorithm_name = "CS",
                                    train.expr = NULL, moic.res = NULL, test.expr = NULL, 
                                    gene.subset = NULL) 
{
  comsam <- intersect(moic.res$clust.res$samID, colnames(train.expr))
  if (length(comsam) == nrow(moic.res$clust.res)) {
    message("--all samples matched.")
  }
  else {
    message(paste0("--", (nrow(moic.res$clust.res) - length(comsam)), 
                   " samples mismatched from current subtypes."))
  }
  moic.res$clust.res <- moic.res$clust.res[comsam, , drop = FALSE]
  train.expr <- train.expr[, comsam]
  if (is.null(gene.subset)) {
    comgene <- intersect(rownames(train.expr), rownames(test.expr))
    message(paste0("--a total of ", length(comgene), " genes shared and used."))
  }
  else {
    comgene <- intersect(intersect(rownames(train.expr), 
                                   rownames(test.expr)), gene.subset)
    message(paste0("--a total of ", length(comgene), " genes shared in the gene subset and used."))
  }
  train.expr <- train.expr[comgene, ]
  test.expr <- test.expr[comgene, ]
  train.subt <- paste0(algorithm_name, moic.res$clust.res$clust)
  if (max(train.expr) < 25 | (max(train.expr) >= 25 & min(train.expr) < 
                              0)) {
    message("--training expression profile seems to have been standardised (z-score or log transformation), no more action will be performed.")
    train.expr <- train.expr
  }
  if (max(train.expr) >= 25 & min(train.expr) >= 0) {
    message("--log2 transformation done for training expression data.")
    train.expr <- log2(train.expr + 1)
  }
  train.expr <- as.data.frame(t(scale(t(train.expr))))
  if (max(test.expr) < 25 | (max(test.expr) >= 25 & min(test.expr) < 
                             0)) {
    message("--testing expression profile seems to have been standardised (z-score or log transformation), no more action will be performed.")
    test.expr <- test.expr
  }
  if (max(test.expr) >= 25 & min(test.expr) >= 0) {
    message("--log2 transformation done for testing expression data.")
    test.expr <- log2(test.expr + 1)
  }
  test.expr <- as.data.frame(t(scale(t(test.expr))))
  mylist <- list(x = as.matrix(train.expr), y = as.vector(train.subt))
  pamr.classifier <- MOVICS:::quiet(pamr::pamr.train(mylist))
  pamr.pred.test <- clusterRepro::IGP.clusterRepro(Centroids = pamr.classifier$centroids, 
                                     Data = as.matrix(test.expr))
  IGP <- pamr.pred.test$IGP
  names(IGP) <- paste0(algorithm_name, 1:length(IGP))
  ex.moic.res <- data.frame(samID = names(pamr.pred.test$Class), 
                            clust = as.character(pamr.pred.test$Class), row.names = names(pamr.pred.test$Class), 
                            stringsAsFactors = FALSE)
  return(list(IGP = IGP, clust.res = ex.moic.res, mo.method = "PAM"))
}

# runKappa single algorithm #####
runKappa_single_algorithm = function (algorithm_name = "CS",
                                      subt1 = NULL, subt2 = NULL, subt1.lab = NULL, subt2.lab = NULL, 
                                      fig.path = getwd(), fig.name = "constheatmap", width = 5, 
                                      height = 5) 
{
  subt1 <- as.vector(as.character(subt1))
  subt2 <- as.vector(as.character(subt2))
  if (length(subt1) != length(subt2)) {
    stop("subtypes identified from different cohorts.")
  }
  if (!identical(sort(unique(subt1)), sort(unique(subt2)))) {
    stop("subtypes fail to matched from two appraisements.")
  }
  if (is.null(subt1.lab) | is.null(subt2.lab)) {
    stop("label for subtype1 and subtype2 must be both indicated.")
  }
  comb.subt <- data.frame(subt1 = paste0(algorithm_name, subt1), subt2 = paste0(algorithm_name, 
                                                                      subt2), stringsAsFactors = F)
  tab_classify <- as.data.frame.array(table(comb.subt$subt1, 
                                            comb.subt$subt2))
  x <- table(comb.subt$subt1, comb.subt$subt2)
  nr <- nrow(x)
  nc <- ncol(x)
  N <- sum(x)
  Po <- sum(diag(x))/N
  Pe <- sum(rowSums(x) * colSums(x)/N)/N
  kappa <- (Po - Pe)/(1 - Pe)
  seK0 <- sqrt(Pe/(N * (1 - Pe)))
  p.v <- 1 - pnorm(kappa/seK0)
  p.lab <- ifelse(p.v < 0.001, "P < 0.001", paste0("P = ", 
                                                   format(round(p.v, 3), nsmall = 3)))
  blue <- "#204F8D"
  lblue <- "#498EB9"
  dwhite <- "#B6D1E8"
  white <- "#E6EAF7"
  par(bty = "n", mgp = c(2, 0.5, 0), mar = c(4.1, 4.1, 4.1, 
                                             2.1), tcl = -0.25, font.main = 3)
  par(xpd = NA)
  plot(c(0, ncol(tab_classify)), c(0, nrow(tab_classify)), 
       col = "white", xlab = "", xaxt = "n", ylab = "", yaxt = "n")
  title(paste0("Consistency between ", subt1.lab, " and ", 
               subt2.lab, "\nKappa = ", format(round(kappa, 3), nsmall = 3), 
               "\n", p.lab), adj = 0, line = 0)
  axis(2, at = 0.5:(nrow(tab_classify) - 0.5), labels = FALSE)
  text(y = 0.5:(nrow(tab_classify) - 0.5), par("usr")[1], 
       labels = rownames(tab_classify)[nrow(tab_classify):1], 
       srt = 0, pos = 2, xpd = TRUE)
  mtext(paste0("Subtypes derived from ", subt1.lab), side = 2, 
        line = 3)
  axis(1, at = 0.5:(ncol(tab_classify) - 0.5), labels = FALSE)
  text(x = 0.5:(ncol(tab_classify) - 0.5), par("usr")[1] - 
         0.2, labels = colnames(tab_classify), srt = 45, pos = 1, 
       xpd = TRUE)
  mtext(paste0("Subtypes derived from ", subt2.lab), side = 1, 
        line = 3)
  input_matrix <- as.matrix(tab_classify)
  mat.max = max(input_matrix)
  unq.value <- unique(sort(as.vector(input_matrix)))
  rbPal <- colorRampPalette(c(white, dwhite, lblue, blue))
  col.vec <- rbPal(max(unq.value) + 1)
  col.mat <- matrix(NA, byrow = T, ncol = ncol(input_matrix), 
                    nrow = nrow(input_matrix))
  for (i in 1:nrow(col.mat)) {
    for (j in 1:ncol(col.mat)) {
      col.mat[i, j] <- col.vec[input_matrix[i, j] + 1]
    }
  }
  x_size <- ncol(input_matrix)
  y_size <- nrow(input_matrix)
  my_xleft = rep(c(0:(x_size - 1)), each = x_size)
  my_xright = my_xleft + 1
  my_ybottom = rep(c((y_size - 1):0), y_size)
  my_ytop = my_ybottom + 1
  rect(xleft = my_xleft, ybottom = my_ybottom, xright = my_xright, 
       ytop = my_ytop, col = col.mat, border = F)
  text(my_xleft + 0.5, my_ybottom + 0.5, input_matrix, cex = 1.3)
  outFig <- paste0(fig.name, ".pdf")
  invisible(dev.copy2pdf(file = file.path(fig.path, outFig), 
                         width = width, height = height))
}

# getMoHeatmap_single_algorithm #####
getMoHeatmap_single_algorithm = function (algorithm_name = "CS", data = NULL, is.binary = c(FALSE, FALSE, FALSE, FALSE, 
                                                        FALSE, FALSE), row.title = c("Data1", "Data2", "Data3", 
                                                                                     "Data4", "Data5", "Data6"), legend.name = c("Data1", "Data2", 
                                                                                                                                 "Data3", "Data4", "Data5", "Data6"), clust.res = NULL, clust.dend = NULL, 
                             show.col.dend = TRUE, show.colnames = FALSE, show.row.dend = c(TRUE, 
                                                                                            TRUE, TRUE, TRUE, TRUE, TRUE), show.rownames = c(FALSE, 
                                                                                                                                             FALSE, FALSE, FALSE, FALSE, FALSE), clust.dist.row = c("pearson", 
                                                                                                                                                                                                    "pearson", "pearson", "pearson", "pearson", "pearson"), 
                             clust.method.row = c("ward.D", "ward.D", "ward.D", "ward.D", 
                                                  "ward.D", "ward.D"), clust.col = c("#2EC4B6", "#E71D36", 
                                                                                     "#FF9F1C", "#BDD5EA", "#FFA5AB", "#011627", "#023E8A", 
                                                                                     "#9D4EDD", "#f09c6c", "#09f3b3"), color = rep(list(c("#00FF00", "#000000", 
                                                                                                                    "#FF0000")), length(data)), annCol = NULL, annColors = NULL, 
                             annRow = NULL, width = 6, height = 4, fig.path = getwd(), 
                             fig.name = "moheatmap") 
{
  ht_opt$message = FALSE
  defaultW <- getOption("warn")
  options(warn = -1)
  if (is.null(names(data))) {
    names(data) <- sprintf("dat%s", 1:length(data))
  }
  n_dat <- length(data)
  if (n_dat > 6) {
    stop("current verision of MOVICS can support up to 6 datasets.")
  }
  if (n_dat < 2) {
    stop("current verision of MOVICS needs at least 2 omics data.")
  }
  colvec <- clust.col[1:length(unique(clust.res$clust))]
  names(colvec) <- paste0(algorithm_name, sort(unique(clust.res$clust)))
  if (!is.null(annCol) & !is.null(annColors)) {
    annCol <- annCol[colnames(data[[1]]), , drop = FALSE]
    annCol$Subtype <- paste0(algorithm_name, clust.res[colnames(data[[1]]), "clust"])
    annColors[["Subtype"]] <- colvec
    if (is.null(clust.dend)) {
      clust.res <- clust.res[order(clust.res$clust), ]
      annCol <- annCol[clust.res$samID, , drop = FALSE]
    }
    ha <- ComplexHeatmap::HeatmapAnnotation(df = annCol, 
                                            col = annColors, border = FALSE)
  }
  else {
    annCol <- data.frame(Subtype = paste0(algorithm_name, clust.res[colnames(data[[1]]), 
                                                                    "clust"]), row.names = colnames(data[[1]]), stringsAsFactors = FALSE)
    annColors <- list(Subtype = colvec)
    if (is.null(clust.dend)) {
      clust.res <- clust.res[order(clust.res$clust), ]
      annCol <- annCol[clust.res$samID, , drop = FALSE]
    }
    ha <- ComplexHeatmap::HeatmapAnnotation(df = annCol, 
                                            col = annColors, border = FALSE)
  }
  if (!is.null(annRow)) {
    if (!is.list(annRow)) {
      stop("argument of annRow should be a list!")
    }
  }
  ht <- list()
  for (i in 1:n_dat) {
    hcg <- fastcluster::hclust(ClassDiscovery::distanceMatrix(as.matrix(t(data[[i]])), 
                                                 clust.dist.row[i]), clust.method.row[i])
    if (is.null(annRow[[i]][1])) {
      rowlab <- ""
      rowlab.index <- 0
    }
    else if (is.na(annRow[[i]][1])) {
      rowlab <- ""
      rowlab.index <- 0
    }
    else {
      rowlab <- intersect(rownames(data[[i]]), annRow[[i]])
      rowlab.index <- match(rowlab, rownames(data[[i]]))
    }
    if (is.null(clust.dend)) {
      data <- lapply(data, function(x) x[, clust.res$samID])
      if (!is.binary[i]) {
        ht[[i]] <- ComplexHeatmap::Heatmap(matrix = as.matrix(data[[i]]), 
                                           row_title = row.title[i], name = legend.name[i], 
                                           cluster_columns = FALSE, cluster_rows = hcg, 
                                           show_column_dend = FALSE, show_column_names = show.colnames, 
                                           show_row_dend = show.row.dend[i], show_row_names = show.rownames[i], 
                                           col = (grDevices::colorRampPalette(color[[i]]))(64), 
                                           top_annotation = switch((i == 1) + 1, NULL, 
                                                                   ha), width = grid::unit(width, "cm"), height = grid::unit(height, 
                                                                                                                             "cm"), heatmap_legend_param = list(at = pretty(range(data[[i]])), 
                                                                                                                                                                labels = pretty(range(data[[i]]))), right_annotation = ComplexHeatmap::rowAnnotation(link = anno_mark(at = rowlab.index, 
                                                                                                                                                                                                                                                                      labels = rowlab, which = "row", lines_gp = grid::gpar(fontsize = 5), 
                                                                                                                                                                                                                                                                      link_width = grid::unit(3, "mm"), padding = grid::unit(0.8, 
                                                                                                                                                                                                                                                                                                                             "mm"), labels_gp = grid::gpar(fontsize = 7))))
      }
      else {
        col_fun = circlize::colorRamp2(c(0, 1), color[[i]])
        ht[[i]] <- ComplexHeatmap::Heatmap(matrix = as.matrix(data[[i]]), 
                                           row_title = row.title[i], name = legend.name[i], 
                                           cluster_columns = FALSE, cluster_rows = hcg, 
                                           show_column_dend = FALSE, show_column_names = show.colnames, 
                                           show_row_dend = show.row.dend[i], show_row_names = show.rownames[i], 
                                           col = color[[i]], top_annotation = switch((i == 
                                                                                        1) + 1, NULL, ha), width = grid::unit(width, 
                                                                                                                              "cm"), height = grid::unit(height, "cm"), 
                                           heatmap_legend_param = list(at = c(0, 1), 
                                                                       legend_gp = grid::gpar(fill = col_fun(c(0, 
                                                                                                               1))), labels = c("0", "1")), right_annotation = ComplexHeatmap::rowAnnotation(link = anno_mark(at = rowlab.index, 
                                                                                                                                                                                                              labels = rowlab, which = "row", lines_gp = grid::gpar(fontsize = 5), 
                                                                                                                                                                                                              link_width = grid::unit(3, "mm"), padding = grid::unit(0.8, 
                                                                                                                                                                                                                                                                     "mm"), labels_gp = grid::gpar(fontsize = 7))))
      }
    }
    else {
      if (!is.binary[i]) {
        ht[[i]] <- ComplexHeatmap::Heatmap(matrix = as.matrix(data[[i]]), 
                                           row_title = row.title[i], name = legend.name[i], 
                                           cluster_columns = clust.dend, cluster_rows = hcg, 
                                           show_column_dend = show.col.dend, show_column_names = show.colnames, 
                                           show_row_dend = show.row.dend[i], show_row_names = show.rownames[i], 
                                           col = (grDevices::colorRampPalette(color[[i]]))(64), 
                                           top_annotation = switch((i == 1) + 1, NULL, 
                                                                   ha), width = grid::unit(width, "cm"), height = grid::unit(height, 
                                                                                                                             "cm"), heatmap_legend_param = list(at = pretty(range(data[[i]])), 
                                                                                                                                                                labels = pretty(range(data[[i]]))), right_annotation = ComplexHeatmap::rowAnnotation(link = anno_mark(at = rowlab.index, 
                                                                                                                                                                                                                                                                      labels = rowlab, which = "row", lines_gp = grid::gpar(fontsize = 5), 
                                                                                                                                                                                                                                                                      link_width = grid::unit(3, "mm"), padding = grid::unit(0.8, 
                                                                                                                                                                                                                                                                                                                             "mm"), labels_gp = grid::gpar(fontsize = 7))))
      }
      else {
        col_fun = circlize::colorRamp2(c(0, 1), color[[i]])
        ht[[i]] <- ComplexHeatmap::Heatmap(matrix = as.matrix(data[[i]]), 
                                           row_title = row.title[i], name = legend.name[i], 
                                           cluster_columns = clust.dend, cluster_rows = hcg, 
                                           show_column_dend = show.col.dend, show_column_names = show.colnames, 
                                           show_row_dend = show.row.dend[i], show_row_names = show.rownames[i], 
                                           col = color[[i]], top_annotation = switch((i == 
                                                                                        1) + 1, NULL, ha), width = grid::unit(width, 
                                                                                                                              "cm"), height = grid::unit(height, "cm"), 
                                           heatmap_legend_param = list(at = c(0, 1), 
                                                                       legend_gp = grid::gpar(fill = col_fun(c(0, 
                                                                                                               1))), labels = c("0", "1")), right_annotation = ComplexHeatmap::rowAnnotation(link = anno_mark(at = rowlab.index, 
                                                                                                                                                                                                              labels = rowlab, which = "row", lines_gp = grid::gpar(fontsize = 5), 
                                                                                                                                                                                                              link_width = grid::unit(3, "mm"), padding = grid::unit(0.8, 
                                                                                                                                                                                                                                                                     "mm"), labels_gp = grid::gpar(fontsize = 7))))
      }
    }
  }
  if (n_dat == 1) {
    ht_list <- ht[[1]]
  }
  if (n_dat == 2) {
    ht_list <- ht[[1]] %v% ht[[2]]
  }
  if (n_dat == 3) {
    ht_list <- ht[[1]] %v% ht[[2]] %v% ht[[3]]
  }
  if (n_dat == 4) {
    ht_list <- ht[[1]] %v% ht[[2]] %v% ht[[3]] %v% ht[[4]]
  }
  if (n_dat == 5) {
    ht_list <- ht[[1]] %v% ht[[2]] %v% ht[[3]] %v% ht[[4]] %v% 
      ht[[5]]
  }
  if (n_dat == 6) {
    ht_list <- ht[[1]] %v% ht[[2]] %v% ht[[3]] %v% ht[[4]] %v% 
      ht[[5]] %v% ht[[6]]
  }
  outFile <- file.path(fig.path, paste0(fig.name, ".pdf"))
  if (is.null(annCol)) {
    pdf(outFile, width = width, height = height * n_dat/2)
  }
  else {
    pdf(outFile, width = width, height = height * n_dat/1.5)
  }
  draw(ht_list, merge_legend = TRUE, heatmap_legend_side = "right")
  invisible(dev.off())
  draw(ht_list, merge_legend = TRUE, heatmap_legend_side = "right")
  options(warn = defaultW)
}

getMoHeatmap_single_algorithm2 = function(
    algorithm_name = "CS",
    data = NULL,
    is.binary = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
    row.title = c("Data1", "Data2", "Data3", "Data4", "Data5", "Data6"),
    legend.name = c("Data1", "Data2", "Data3", "Data4", "Data5", "Data6"),
    # Preexisting clustering (data frame with columns 'samID' and 'clust')
    clust.res = NULL,
    # Back-compat (unused) 
    clust.dend = NULL,
    # Control row & column clustering
    cluster_rows = rep(FALSE, 6),
    cluster_cols = rep(FALSE, 6),
    # Show/hide the dendrogram once computed
    show.row.dend = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
    show.col.dend = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
    show.rownames = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
    show.colnames = FALSE,
    # Dist/method for row clustering
    clust.dist.row = c("pearson","pearson","pearson","pearson","pearson","pearson"),
    clust.method.row = c("ward.D","ward.D","ward.D","ward.D","ward.D","ward.D"),
    # Dist/method for column clustering
    clust.dist.col = c("pearson","pearson","pearson","pearson","pearson","pearson"),
    clust.method.col = c("ward.D","ward.D","ward.D","ward.D","ward.D","ward.D"),
    # Colors for the subtypes in annotation
    clust.col = c("#2EC4B6","#E71D36","#FF9F1C","#BDD5EA","#FFA5AB",
                  "#011627","#023E8A","#9D4EDD","#f09c6c","#09f3b3"),
    # Color palettes for each dataset
    color = rep(list(c("#00FF00", "#000000", "#FF0000")), 6),
    annCol = NULL,
    annColors = NULL,
    annRow = NULL,
    width = 6,
    height = 4,
    fig.path = getwd(),
    fig.name = "moheatmap"
) {
  # Turn off messages from ComplexHeatmap
  ht_opt$message = FALSE
  
  # Temporarily suppress warnings
  defaultW <- getOption("warn")
  options(warn = -1)
  
  # Ensure data has names
  if (is.null(names(data))) {
    names(data) <- sprintf("dat%s", seq_along(data))
  }
  n_dat <- length(data)
  if (n_dat > 6) {
    stop("Current version can support up to 6 datasets.")
  }
  if (n_dat < 2) {
    stop("Need at least 2 omics data.")
  }
  
  #---------------------------------------------
  # 1) Possibly reorder columns by 'clust.res'
  #---------------------------------------------
  # We do a single reorder if ANY dataset has cluster_cols[i] = FALSE
  # and we have a non-null clust.res. This ensures consistent sample
  # order across all datasets in that scenario.
  
  reorder_needed <- any(cluster_cols == FALSE) && !is.null(clust.res)
  
  if (reorder_needed) {
    # sort clust.res by cluster
    # (assuming 'clust.res' has columns: samID, clust)
    clust.res <- clust.res[order(clust.res$clust), , drop = FALSE]
    # reorder all data sets to match that sample order
    data <- lapply(data, function(mat) {
      # keep only samples that are in clust.res$samID (and in the same order)
      mat[, clust.res$samID, drop = FALSE]
    })
  }
  
  #---------------------------------------------
  # 2) Build color mapping for the clusters
  #---------------------------------------------
  # If user gave us a clust.res, define colvec accordingly
  if (!is.null(clust.res)) {
    unique_clusters <- sort(unique(clust.res$clust))
    colvec <- clust.col[seq_along(unique_clusters)]
    names(colvec) <- paste0(algorithm_name, unique_clusters)
  } else {
    # fallback if no clust.res
    colvec <- clust.col
  }
  
  #---------------------------------------------
  # 3) Validate annRow
  #---------------------------------------------
  if (!is.null(annRow) && !is.list(annRow)) {
    stop("Argument 'annRow' should be a list if provided!")
  }
  
  #---------------------------------------------
  # 4) Helper to compute distances via Rfast or cor
  #---------------------------------------------
  compute_distance <- function(mat, dist_method, do_column = FALSE) {
    if (dist_method %in% c("pearson","spearman")) {
      # correlation-based
      if (do_column) {
        # columns => cor(mat)
        cmat <- cor(mat, method = dist_method, use="pairwise.complete.obs")
      } else {
        # rows => cor(t(mat))
        cmat <- cor(t(mat), method = dist_method, use="pairwise.complete.obs")
      }
      dist_mat <- as.dist(1 - cmat)
    } else {
      if (!requireNamespace("Rfast", quietly = TRUE)) {
        stop("Package 'Rfast' not installed. Please install it or change dist method.")
      }
      if (do_column) {
        dist_full <- Rfast::Dist(t(mat), method=dist_method)
      } else {
        dist_full <- Rfast::Dist(mat, method=dist_method)
      }
      dist_mat <- as.dist(dist_full)
    }
    dist_mat
  }
  
  #---------------------------------------------
  # 5) Build Heatmaps for Each Data Set
  #---------------------------------------------
  ht_list <- vector("list", n_dat)
  
  for (i in seq_len(n_dat)) {
    
    #--- Column Clustering or None
    if (cluster_cols[i]) {
      # do a new column dendrogram
      col_dist_mat <- compute_distance(data[[i]], clust.dist.col[i], do_column=TRUE)
      col_hclust   <- fastcluster::hclust(col_dist_mat, method = clust.method.col[i])
    } else {
      # no new column clustering
      col_hclust <- FALSE
      # The columns are already reordered above if reorder_needed=TRUE
    }
    
    #--- Row Clustering or None
    if (cluster_rows[i]) {
      row_dist_mat <- compute_distance(data[[i]], clust.dist.row[i], do_column=FALSE)
      row_hclust   <- fastcluster::hclust(row_dist_mat, method = clust.method.row[i])
    } else {
      row_hclust <- FALSE
    }
    
    #--- Figure out row labels for anno_mark
    if (!is.null(annRow) && !is.null(annRow[[i]]) && !is.na(annRow[[i]][1])) {
      rowlab       <- intersect(rownames(data[[i]]), annRow[[i]])
      rowlab.index <- match(rowlab, rownames(data[[i]]))
    } else {
      rowlab       <- character(0)
      rowlab.index <- integer(0)
    }
    
    #--- top_annotation only on the first heatmap
    #---------------------------------------------
    # Build Sample Annotation AFTER reorder
    #---------------------------------------------
    # Now annCol matches the new column order in data[[1]] if we did reorder
    if (!is.null(annCol) && !is.null(annColors)) {
      # subset annCol to match the new colnames of data[[1]]
      annCol <- annCol[colnames(data[[1]]), , drop=FALSE]
      if (!is.null(clust.res)) {
        # Add Subtype info
        annCol$Subtype <- paste0(algorithm_name, clust.res[colnames(data[[1]]), "clust"])
        annColors[["Subtype"]] <- colvec
      }
      ha <- ComplexHeatmap::HeatmapAnnotation(df = annCol,
                                              col = annColors,
                                              border = FALSE)
    } else {
      # minimal annotation
      if (!is.null(clust.res)) {
        # build minimal df with cluster info
        minimalAnn <- data.frame(
          Subtype = paste0(algorithm_name, clust.res[colnames(data[[1]]), "clust"]),
          row.names = colnames(data[[1]]),
          stringsAsFactors = FALSE
        )
        annColorsLocal <- list(Subtype = colvec)
        ha <- ComplexHeatmap::HeatmapAnnotation(df = minimalAnn,
                                                col = annColorsLocal,
                                                border = FALSE)
      } else {
        ha <- NULL
      }
    }
    top_annotation <- if (i == 1) ha else NULL
    
    #--- Build the Heatmap object
    if (!is.binary[i]) {
      # continuous data
      rng <- range(data[[i]], na.rm=TRUE)
      cont_colors <- grDevices::colorRampPalette(color[[i]])(64)
      ht_list[[i]] <- ComplexHeatmap::Heatmap(
        matrix            = as.matrix(data[[i]]),
        row_title         = row.title[i],
        name              = legend.name[i],
        cluster_rows      = row_hclust,
        cluster_columns   = col_hclust,
        show_row_dend     = show.row.dend[i],
        show_column_dend  = show.col.dend[i],
        show_row_names    = show.rownames[i],
        show_column_names = show.colnames,
        col = cont_colors,
        top_annotation = top_annotation,
        width  = grid::unit(width, "cm"),
        height = grid::unit(height, "cm"),
        heatmap_legend_param = list(
          at     = pretty(rng),
          labels = pretty(rng)
        ),
        right_annotation = ComplexHeatmap::rowAnnotation(
          link = ComplexHeatmap::anno_mark(
            at = rowlab.index,
            labels = rowlab,
            which = "row",
            lines_gp = grid::gpar(fontsize = 5),
            link_width = grid::unit(3, "mm"),
            padding = grid::unit(0.8, "mm"),
            labels_gp = grid::gpar(fontsize = 7)
          )
        )
      )
    } else {
      # binary data
      col_fun_bin <- circlize::colorRamp2(c(0, 1), color[[i]])
      ht_list[[i]] <- ComplexHeatmap::Heatmap(
        matrix            = as.matrix(data[[i]]),
        row_title         = row.title[i],
        name              = legend.name[i],
        cluster_rows      = row_hclust,
        cluster_columns   = col_hclust,
        show_row_dend     = show.row.dend[i],
        show_column_dend  = show.col.dend[i],
        show_row_names    = show.rownames[i],
        show_column_names = show.colnames,
        col = color[[i]],   # or col_fun_bin
        top_annotation = top_annotation,
        width  = grid::unit(width, "cm"),
        height = grid::unit(height, "cm"),
        heatmap_legend_param = list(
          at = c(0, 1),
          legend_gp = grid::gpar(fill = col_fun_bin(c(0,1))),
          labels = c("0", "1")
        ),
        right_annotation = ComplexHeatmap::rowAnnotation(
          link = ComplexHeatmap::anno_mark(
            at = rowlab.index,
            labels = rowlab,
            which = "row",
            lines_gp = grid::gpar(fontsize = 5),
            link_width = grid::unit(3, "mm"),
            padding = grid::unit(0.8, "mm"),
            labels_gp = grid::gpar(fontsize = 7)
          )
        )
      )
    }
  }
  
  #---------------------------------------------
  # 6) Combine & Output
  #---------------------------------------------
  ht_combined <- ht_list[[1]]
  if (n_dat > 1) {
    for (k in 2:n_dat) {
      ht_combined <- ht_combined %v% ht_list[[k]]
    }
  }
  
  outFile <- file.path(fig.path, paste0(fig.name, ".pdf"))
  # PDF height depends on presence of annotation
  if (is.null(annCol)) {
    pdf(outFile, width = width, height = height * n_dat / 2)
  } else {
    pdf(outFile, width = width, height = height * n_dat / 1.5)
  }
  
  draw(ht_combined, merge_legend = TRUE, heatmap_legend_side = "right")
  invisible(dev.off())
  
  # Also draw on current device
  draw(ht_combined, merge_legend = TRUE, heatmap_legend_side = "right")
  
  # Restore warning level
  options(warn = defaultW)
}

# NTP modification for subscript out of bounds error #####
runNTP_mod = function (expr = NULL, templates = NULL, scaleFlag = TRUE, centerFlag = TRUE, 
                       nPerm = 1000, distance = "cosine", seed = 123456, verbose = TRUE, 
                       doPlot = FALSE, fig.path = getwd(), fig.name = "ntpheatmap", 
                       width = 5, height = 5, algorithm_name = "CS") 
{
  library(MOVICS)
  if (!is.element(distance, c("cosine", "pearson", "spearman", 
                              "kendall"))) {
    stop("the argument of distance should be one of cosine, pearson, spearman, or kendall.")
  }
  com_feat <- intersect(rownames(expr), templates$probe)
  message(paste0("--original template has ", nrow(templates), 
                 " biomarkers and ", length(com_feat), " are matched in external expression profile."))
  expr <- expr[com_feat, , drop = FALSE]
  rownames(expr) = com_feat
  templates <- templates[which(templates$probe %in% com_feat), 
                         , drop = FALSE]
  if (is.element(0, as.numeric(table(templates$class)))) {
    stop("at least one class has no probes/genes matched in template file!")
  }
  emat <- t(scale(t(expr), scale = scaleFlag, center = centerFlag))
  if (doPlot) {
    outFig <- paste0(fig.name, ".pdf")
    ntp.res <- ntp_mod(emat = emat, templates = templates, doPlot = doPlot, 
                   nPerm = nPerm, distance = distance, nCores = 1, 
                   seed = seed, verbose = verbose)
    invisible(dev.copy2pdf(file = file.path(fig.path, outFig), 
                           width = width, height = height))
  } else {
    ntp.res <- ntp_mod(emat = emat, templates = templates, doPlot = doPlot, 
                   nPerm = nPerm, distance = distance, nCores = 1, 
                   seed = seed, verbose = verbose)
  }
  ntp.res[, setdiff(colnames(ntp.res), "prediction")] <- round(ntp.res[, 
                                                                       setdiff(colnames(ntp.res), "prediction")], 4)
  ex.moic.res <- data.frame(samID = rownames(ntp.res), clust = gsub(algorithm_name, 
                                                                    "", ntp.res$prediction), row.names = rownames(ntp.res), 
                            stringsAsFactors = FALSE)
  return(list(ntp.res = ntp.res, clust.res = ex.moic.res, 
              mo.method = "NTP"))
}

ntp_mod = function (emat, templates, nPerm = 1000, distance = "cosine", 
                    nCores = 1, seed = NULL, verbose = getOption("verbose"), 
                    doPlot = FALSE) 
{
  library(CMScaller)
  if (class(emat)[1] == "ExpressionSet") 
    emat <- suppressPackageStartupMessages(Biobase::exprs(emat))
  if (is.data.frame(emat) | is.vector(emat)) 
    emat <- as.matrix(emat)
  if (is.null(rownames(emat))) 
    stop("missing emat rownames, check input")
  if (is.null(templates$class) | is.null(templates$probe)) {
    stop("missing columns in templates, check input")
  }
  if (is.character(templates$class)) 
    templates$class <- factor(templates$class)
  if (is.factor(templates$probe) | is.integer(templates$probe)) {
    warning("templates$probe coerced to character", call. = FALSE)
    templates$probe <- as.character(templates$probe)
  }
  if (is.character(distance) & isTRUE(verbose)) {
    message(paste0(distance, " correlation distance"))
  }
  keepP <- stats::complete.cases(emat)
  if (sum(!keepP) > 0) {
    if (isTRUE(verbose)) 
      message(paste0(sum(!keepP), "/", length(keepP), 
                     " features NA, discarded"))
    emat <- emat[keepP, , drop = FALSE]
  }
  keepT <- templates$probe %in% rownames(emat)
  if (sum(!keepT) > 0) {
    if (isTRUE(verbose)) 
      message(paste0(sum(!keepT), "/", length(keepT), 
                     " templates features not in emat, discarded"))
    templates <- templates[keepT, ]
  }
  if (min(table(templates$class)) < 2) {
    message("<2 matched features/class")
    stop("check templates$probe is matchable against rownames(emat)", 
         call. = FALSE)
  }
  if (min(table(templates$class)) < 5) {
    warning("<5 matched features/class - unstable predictions", 
            call. = FALSE)
  }
  N <- ncol(emat)
  K <- nlevels(templates$class)
  S <- nrow(templates)
  P <- nrow(emat)
  # cat("N = ", N, ", K = ", K, ", S = ", S, ", P = ", P, "\n", sep = "")
  class.names <- levels(templates$class)
  templates$class <- as.numeric(templates$class)
  emat.mean <- round(mean(emat), 2)
  if (abs(emat.mean) > 1) {
    isnorm <- " <- check feature centering!"
    emat.sd <- round(stats::sd(emat), 2)
    warning(paste0("emat mean=", emat.mean, "; sd=", emat.sd, 
                   isnorm), call. = FALSE)
  }
  feat.class <- paste(range(table(templates$class)), collapse = "-")
  if (isTRUE(verbose)) 
    message(paste0(N, " samples; ", K, " classes; ", feat.class, 
                   " features/class"))
  mm <- match(templates$probe, rownames(emat), nomatch = 0)
  mm <- mm[mm > 0]  # Remove unmatched probes
  # cat("The length of mm is ", length(mm), "\n")
  if (!all(rownames(emat)[mm] == templates$probe)) {
    stop("error matching probes, check rownames(emat) and templates$probe")
  }
  pReplace <- length(templates$probe) > length(unique(templates$probe))
  tmat <- matrix(rep(templates$class, K), ncol = K)
  for (k in seq_len(K)) tmat[, k] <- as.numeric(tmat[, k] == 
                                                  k)
  if (K == 2) 
    tmat[tmat == 0] <- -1
  if (!distance %in% c("cosine", "pearson", "spearman", "kendall")) 
    stop("invalid distance method")
  if (distance == "cosine") {
    simFun <- function(x, y) corCosine(x, y)
  } else {
    simFun <- function(x, y) {
      stats::cor(x, y, method = distance)
    }
  }
  # cat("The dimensions of tmat are ", paste(dim(tmat), sep = " x "), "\n")
  ntpFUN <- function(n) {
    n.sim <- as.vector(simFun(emat[mm, n, drop = FALSE], 
                              tmat))
    n.sim.perm.max <- apply(simFun(matrix(emat[, n][sample.int(P, 
                                                               S * nPerm, replace = TRUE)], ncol = nPerm), tmat), 
                            1, max)
    n.ntp <- which.max(n.sim)
    n.sim.ranks <- rank(-c(n.sim[n.ntp], (n.sim.perm.max)))
    n.pval <- n.sim.ranks[1]/length(n.sim.ranks)
    return(c(n.ntp, CMScaller:::simToDist(n.sim), n.pval))
  }
  # ntpFUN <- function(n, iter) {
  #   tryCatch({
  #     # Ensure n is within valid column index range
  #     if (n > ncol(emat)) {
  #       stop(paste("Column index out of bounds at iteration:", iter))
  #     }
  #     
  #     # Debugging: Print dimensions of emat
  #     message(paste("Iteration:", iter))
  #     cat("The dimensions of the subset are: ", print(dim(emat[mm, n, drop = FALSE])), "\n")  # Print dimensions of the subsetting
  #     
  #     # Adjust simFun to work with conformable matrices
  #     # Now we are only selecting the nth sample column
  #     n.sim <- as.vector(simFun(emat[mm, n, drop = FALSE], tmat))  # Correct dimension match
  #     
  #     # Adjust the permutation step to ensure valid matrix dimensions
  #     n.sim.perm.max <- apply(simFun(matrix(emat[, n][sample.int(P, min(S, P) * nPerm, replace = TRUE)], ncol = nPerm), tmat), 1, max)
  #     
  #     # Compute the maximum similarity and rank
  #     n.ntp <- which.max(n.sim)
  #     n.sim.ranks <- rank(-c(n.sim[n.ntp], (n.sim.perm.max)))
  #     n.pval <- n.sim.ranks[1] / length(n.sim.ranks)
  #     
  #     return(c(n.ntp, CMScaller:::simToDist(n.sim), n.pval))
  #     
  #   }, error = function(e) {
  #     message(paste("Error encountered during iteration:", iter))
  #     message("Error message:", conditionMessage(e))
  #     return(NULL)  # Return NULL or handle the error gracefully
  #   })
  # }
  if (!is.null(seed)) {
    set.seed(seed)
    nCores <- 1
  }
  existParallel <- CMScaller:::packageExists("parallel")
  existSnow <- CMScaller:::packageExists("snow")
  if ((existParallel | existParallel) & nCores != 1) {
    funVal <- vector(mode = "numeric", length = 2 + K)
    if (.Platform$OS.type == "windows") {
      nSlaves <- ifelse(nCores == 0, as.numeric(Sys.getenv("NUMBER_OF_PROCESSORS")), 
                        nCores)
      if (nSlaves == 0) 
        nSlaves <- parallel::detectCores()
      if (isTRUE(verbose)) 
        message(paste0("parallel; ", nSlaves, " sockets pkg:snow", 
                       "; ", nPerm, " permutation(s)..."))
      nParts <- split(seq_len(N), cut(seq_len(N), nSlaves, 
                                      labels = FALSE))
      cl <- snow::makeCluster(nSlaves, type = "SOCK")
      res <- snow::parLapply(cl, nParts, function(n) vapply(n, 
                                                            ntpFUN, funVal))
      snow::stopCluster(cl)
      res <- data.frame(t(do.call(cbind, res)))
    }
    else {
      nCores <- ifelse(nCores == 0, parallel::detectCores(), 
                       nCores)
      if (isTRUE(verbose)) 
        message(paste0("parallel; ", nCores, " cores pkg:parallel", 
                       "; ", nPerm, " permutation(s)..."))
      options(mc.cores = nCores)
      nParts <- split(seq_len(N), cut(seq_len(N), nCores, 
                                      labels = FALSE))
      res <- parallel::mclapply(nParts, function(n) vapply(n, 
                                                           ntpFUN, funVal))
      res <- data.frame(t(do.call(cbind, res)))
    }} else {
      if (isTRUE(verbose)) 
        message(paste0("serial processing; ", nPerm, " permutation(s)..."))
      res <- lapply(seq_len(N), ntpFUN)
      res <- data.frame(do.call(rbind, res))
    }
  # } else {
  #   if (isTRUE(verbose)) 
  #     message(paste0("serial processing; ", nPerm, " permutation(s)..."))
  #   # Perform matching and ensure templates$probe matches rownames(emat)
  #   templates$probe <- tolower(trimws(templates$probe))
  #   rownames(emat) <- tolower(trimws(rownames(emat)))
  #   # Wrap the entire lapply call inside tryCatch for additional error handling
  #   res <- tryCatch({
  #     lapply(seq_len(N), function(i) {
  #       message(paste("Running iteration:", i))
  #       ntpFUN(i, iter = i)
  #     })
  #   }, error = function(e) {
  #     message("Error in the lapply call:")
  #     message(conditionMessage(e))
  #     return(NULL)
  #   })
  #   res <- data.frame(do.call(rbind, res))
  # }
  colnames(res) <- c("prediction", paste0("d.", class.names), 
                     "p.value")
  res$prediction <- factor(class.names[res$prediction], levels = class.names)
  rownames(res) <- colnames(emat)
  res$p.value[res$p.value < 1/nPerm] <- 1/nPerm
  res$FDR <- stats::p.adjust(res$p.value, "fdr")
  if (isTRUE(doPlot)) {
    subHeatmap_mod(emat = emat, res = res, templates = templates)
  }
  if (isTRUE(verbose)) {
    message("predicted samples/class (FDR<0.05)")
    print(table(suppressMessages(CMScaller::subSetNA(res, FDR = 0.05))$prediction, 
                useNA = "always"))
  }
  return(res)
}

subHeatmap_mod <- function (emat, res, templates, keepN = TRUE, labRow = NULL, 
                            classCol = getOption("subClassCol"), 
                            heatCol = colorRampPalette(c("dodgerblue4", "white", "deeppink4"))(100),
                            # heatCol = rev(colorRampPalette(viridisLite::magma(10))(255)),
                            N = ncol(emat), ...) 
{
  library(MOVICS)
  library(CMScaller)
  
  # Check if a package is installed, if not, install it
  if (!requireNamespace("dlfUtils", quietly = TRUE)) {
    
    # Check if the 'remotes' package is installed
    if (!requireNamespace("remotes", quietly = TRUE)) {
      message("Installing 'remotes' package...")
      install.packages("remotes")
    }
    
    # Install 'dlfUtils' from GitHub
    message("Installing 'dlfUtils' from GitHub...")
    remotes::install_github("daynefiler/dlfUtils")
    
    # Load the installed 'dlfUtils' package
    library(dlfUtils)
  } else {
    # If 'dlfUtils' is already installed, load it
    library(dlfUtils)
  }
  
  # Convert rownames and templates$probe to lowercase and trim whitespace
  message("Standardizing rownames of emat and templates...")
  rownames(emat) <- tolower(trimws(rownames(emat)))
  templates$probe <- tolower(trimws(templates$probe))
  
  # Check if emat/res match
  message("Checking if emat and res dimensions match...")
  if (!all(colnames(emat) == rownames(res))) {
    stop("Error: emat and res do not match in dimensions or names")
  }
  
  # Initial subsetting of emat and res based on keepN
  emat <- emat[, keepN]
  res <- res[keepN, ]
  
  class <- res$prediction
  K <- nlevels(class)
  
  # Check for duplicates in templates
  message("Checking for duplicated rows in templates...")
  if (length(which(duplicated(templates$probe))) > 1) {
    message(paste("Number of duplicated rows to be removed:", sum(duplicated(templates$probe))))
  }
  templates <- templates[!duplicated(templates$probe), ]
  rownames(templates) <- templates$probe
  
  
  # Intersect templates and emat rownames
  P <- intersect(rownames(emat), rownames(templates))
  
  # Debug: Check if intersection is working as expected
  if (length(P) == 0) {
    stop("Error: No intersecting rownames between emat and templates")
  }
  
  templates <- templates[P, ]
  rownames(templates) = P
  
  # Ensure emat has the same rownames as ordered templates
  message("Subsetting and adjusting emat using rownames of templates...")
  
  # Debugging: Find mismatched rownames
  unmatched_in_templates <- setdiff(rownames(templates), rownames(emat))
  unmatched_in_emat <- setdiff(rownames(emat), rownames(templates))
  
  if (length(unmatched_in_templates) > 0) {
    message("Unmatched rownames in templates that are not in emat:")
    print(unmatched_in_templates)
  }
  
  if (length(unmatched_in_emat) > 0) {
    message("Unmatched rownames in emat that are not in templates:")
    print(unmatched_in_emat)
  }
  
  if (length(unmatched_in_templates) > 0 || length(unmatched_in_emat) > 0) {
    stop("Error: some rownames in templates do not match rownames in emat")
  }
  
  emat <- ematAdjust(emat[rownames(templates), , drop = FALSE])
  
  # Ensure heatCol is valid
  if (is.null(heatCol)) {
    heatCol <- subData[["hmCol"]]
    message("Using default heatCol from subData")
  }
  
  # Dynamically calculate breaks to have one more than the number of colors
  pMax <- 3
  breaks <- seq(-pMax, pMax, length.out = length(heatCol) + 1)  # Ensure correct number of breaks
  
  # Ensure the data is within the pMax range
  emat[emat > pMax] <- pMax
  emat[emat < -pMax] <- -pMax
  
  # Plot heatmap and class predictions
  message("Preparing heatmap plot...")
  xx <- seq(0, 1, length.out = ncol(emat) + 1)
  yy <- seq(0, 1, length.out = nrow(emat) + 1)
  graphics::image(x = xx, y = yy, z = t(emat), yaxt = "n", 
                  xaxt = "n", useRaster = TRUE, col = heatCol, breaks = breaks, 
                  xlab = "class predictions", ylab = "template features", ...)
  
  # Plot class rectangles
  xb <- cumsum(c(0, (sapply(split(res$prediction, res$prediction), length) / length(res$prediction))))
  xl <- xb[-(K + 1)]
  xr <- xb[-1]
  yy <- dlfUtils::line2user(line = 1:0, side = 1)
  yy <- seq(yy[1], yy[2], length = 3)
  
  graphics::rect(xleft = xl, xright = xr, ybottom = yy[1], 
                 ytop = yy[2], col = classCol, xpd = TRUE, border = FALSE, 
                 lwd = 0.75)
  
  # Plot class labels
  xxl <- xl + (xr - xl) / 2
  graphics::text(xxl, yy[1], pos = 1, levels(class), adj = 0, 
                 xpd = TRUE, cex = 0.75)
  
  # Plot p-values if available
  if (min(res$p.value) < 1) {
    py <- abs(res$p.value[N])
    px <- (seq_along(N) / length(N)) - 1 / length(N) / 2
    py <- (py * (yy[2] - yy[1])) + yy[1]
    
    graphics::rect(xleft = xx[-length(xx)], xright = xx[-1], 
                   ybottom = yy[1], ytop = py, border = FALSE, xpd = NA, 
                   col = "white", lwd = 2)
    
    graphics::rect(xleft = xl, xright = xr, ybottom = yy[1], 
                   ytop = yy[2], col = NA, xpd = TRUE, lwd = 0.75)
    
    graphics::text(0, yy[1], pos = 2, expression(italic(p) - 
                                                   value), xpd = TRUE, cex = 0.75)
    graphics::text(x = 1, y = c(yy[1], yy[2]), c(0, 1), 
                   pos = c(4), xpd = TRUE, cex = 0.75)
    graphics::text(4, 2, paste(bquote(italic(p) - value)))
  }
  
  # Plot the class bars on the side of the heatmap
  xx <- dlfUtils::line2user(line = 1:0, side = 2)
  xx <- seq(xx[1], xx[2], length = 3)
  yy <- c(0, cumsum(sapply(split(templates$probe, templates$class), 
                           length) / length(templates$probe)))
  yb <- yy[-(K + 1)]
  yt <- yy[-1]
  
  graphics::rect(xleft = xx[1], xright = xx[2], ybottom = yb, 
                 ytop = yt, col = classCol, xpd = NA, lwd = 0.75)
  
  # Add row labels if labRow is not NULL
  if (!is.null(labRow)) {
    message("Adding row labels...")
    textfun <- function(..., cex.text) graphics::text(..., 
                                                      xpd = TRUE, cex = 0.75)
    id <- templates[, labRow]
    yy <- seq(0 + (1 / nrow(templates) / 2), 1 - (1 / nrow(templates) / 2), 
              length.out = nrow(templates))
    textfun(1, yy, id, pos = 4, adj = 1)
  }
  
  message("Heatmap plotting complete.")
}

# Modified gain definition for FGA #####
compFGA_mod = function (moic.res = NULL, segment = NULL, iscopynumber = FALSE, ga_column = NULL,
                        cnathreshold = 0.2, test.method = "nonparametric", barcolor = c("#008B8A", 
                                                                                        "#F2042C", "#21498D"), clust.col = c("#2EC4B6", "#E71D36", 
                                                                                                                             "#FF9F1C", "#BDD5EA", "#FFA5AB", "#011627", "#023E8A", 
                                                                                                                             "#9D4EDD", "#f09c6c", "#09f3b3"), 
                        fig.path = getwd(), fig.name = NULL, width = 8, 
                        height = 4, prefix = "", title = NULL) 
{
  library(patchwork)
  
  if (!all(is.element(c("sample", "chrom", "start", "end", "value"), colnames(segment)))) {
    stop("segment data must have the following columns: sample, chrom, start, end, value.")
  }
  
  if (iscopynumber) {
    segment$value <- log2(segment$value / 2)
  }
  
  comsam <- intersect(moic.res$clust.res$samID, unique(segment$sample))
  
  if (length(comsam) == nrow(moic.res$clust.res)) {
    message("--all samples matched.")
  } else {
    message(paste0("--", (nrow(moic.res$clust.res) - length(comsam)), 
                   " samples mismatched from current subtypes."))
  }
  
  if (!is.element(test.method, c("nonparametric", "parametric"))) {
    stop("test.method can be one of nonparametric or parametric.")
  }
  
  clust.res <- moic.res$clust.res[comsam, , drop = FALSE]
  segment <- segment[which(segment$sample %in% comsam), ]
  n.moic <- length(unique(clust.res$clust))
  segment$bases <- segment$end - segment$start
  
  display.progress <- function(index, totalN, breakN = 20) {
    if (index %% ceiling(totalN / breakN) == 0) {
      cat(paste(round(index * 100 / totalN), "% ", sep = ""))
    }
  }
  
  std <- function(x, na.rm = TRUE) {
    if (na.rm) {
      x <- as.numeric(na.omit(x))
      sd(x) / sqrt(length(x))
    } else {
      sd(x) / sqrt(length(x))
    }
  }
  
  outTab <- data.frame()
  
  for (i in 1:length(unique(segment$sample))) {
    display.progress(index = i, totalN = length(unique(segment$sample)))
    tmp <- segment[segment$sample == names(table(segment$sample))[i], ]
    
    # If ga_column is not NULL, use it to classify gains/losses/normal
    if (!is.null(ga_column) && ga_column %in% colnames(segment)) {
      tmp$classification <- tmp[[ga_column]]
      
      # Use the classification in the specified column
      FGA <- sum(tmp[tmp$classification %in% c("gain", "loss"), "bases"]) / sum(tmp[, "bases"])
      FGG <- sum(tmp[tmp$classification == "gain", "bases"]) / sum(tmp[, "bases"])
      FGL <- sum(tmp[tmp$classification == "loss", "bases"]) / sum(tmp[, "bases"])
      
    } else {
      # Proceed with the original logic if ga_column is NULL
      if (length(tmp[abs(tmp$value) > cnathreshold, "bases"][6]) == 0) {
        FGA = 0
      } else {
        FGA = sum(tmp[abs(tmp$value) > cnathreshold, "bases"]) / sum(tmp[, "bases"])
      }
      
      if (length(tmp[tmp$value > cnathreshold, "bases"][6]) == 0) {
        FGG = 0
      } else {
        FGG = sum(tmp[tmp$value > cnathreshold, "bases"]) / sum(tmp[, "bases"])
      }
      
      if (length(tmp[tmp$value < (-cnathreshold), "bases"][6]) == 0) {
        FGL = 0
      } else {
        FGL = sum(tmp[tmp$value < (-cnathreshold), "bases"]) / sum(tmp[, "bases"])
      }
    }
    
    # Store results for each sample
    tmp_out <- data.frame(samID = names(table(segment$sample))[i], 
                          FGA = FGA, FGG = FGG, FGL = FGL, stringsAsFactors = FALSE)
    outTab <- rbind.data.frame(outTab, tmp_out, stringsAsFactors = FALSE)
  }
  
  outTab$Subtype <- paste0(prefix, clust.res[outTab$samID, "clust"])
  
  # Summarize FGA, FGG, and FGL by Subtype
  summaryFGA <- outTab %>% group_by(Subtype) %>% dplyr::summarize(mean = mean(FGA, na.rm = TRUE), se = std(FGA, na.rm = TRUE))
  summaryFGG <- outTab %>% group_by(Subtype) %>% dplyr::summarize(mean = mean(FGG, na.rm = TRUE), se = std(FGG, na.rm = TRUE))
  summaryFGL <- outTab %>% group_by(Subtype) %>% dplyr::summarize(mean = mean(FGL, na.rm = TRUE), se = std(FGL, na.rm = TRUE))
  
  summaryFGGL <- data.frame(rbind.data.frame(summaryFGG, summaryFGL), 
                            class = rep(c("FGG", "FGL"), c(nrow(summaryFGG), nrow(summaryFGL))), 
                            stringsAsFactors = FALSE)
  
  # Statistical tests based on number of subtypes
  if (n.moic == 2 & test.method == "nonparametric") {
    statistic <- "wilcox.test"
    FGA.test <- wilcox.test(outTab$FGA ~ outTab$Subtype)$p.value
    FGG.test <- wilcox.test(outTab$FGG ~ outTab$Subtype)$p.value
    FGL.test <- wilcox.test(outTab$FGL ~ outTab$Subtype)$p.value
  } else if (n.moic == 2 & test.method == "parametric") {
    statistic <- "t.test"
    FGA.test <- t.test(outTab$FGA ~ outTab$Subtype)$p.value
    FGG.test <- t.test(outTab$FGG ~ outTab$Subtype)$p.value
    FGL.test <- t.test(outTab$FGL ~ outTab$Subtype)$p.value
  } else if (n.moic > 2 & test.method == "nonparametric") {
    statistic <- "kruskal.test"
    FGA.test <- kruskal.test(outTab$FGA ~ outTab$Subtype)$p.value
    FGG.test <- kruskal.test(outTab$FGG ~ outTab$Subtype)$p.value
    FGL.test <- kruskal.test(outTab$FGL ~ outTab$Subtype)$p.value
    pairwise.FGA.test <- pairwise.wilcox.test(outTab$FGA, outTab$Subtype, p.adjust.method = "BH")
    pairwise.FGG.test <- pairwise.wilcox.test(outTab$FGG, outTab$Subtype, p.adjust.method = "BH")
    pairwise.FGL.test <- pairwise.wilcox.test(outTab$FGL, outTab$Subtype, p.adjust.method = "BH")
  } else if (n.moic > 2 & test.method == "parametric") {
    statistic <- "anova"
    FGA.test <- summary(aov(outTab$FGA ~ outTab$Subtype))[[1]][["Pr(>F)"]][1]
    FGG.test <- summary(aov(outTab$FGG ~ outTab$Subtype))[[1]][["Pr(>F)"]][1]
    FGL.test <- summary(aov(outTab$FGL ~ outTab$Subtype))[[1]][["Pr(>F)"]][1]
    pairwise.FGA.test <- pairwise.t.test(outTab$FGA, outTab$Subtype, p.adjust.method = "BH")
    pairwise.FGG.test <- pairwise.t.test(outTab$FGG, outTab$Subtype, p.adjust.method = "BH")
    pairwise.FGL.test <- pairwise.t.test(outTab$FGL, outTab$Subtype, p.adjust.method = "BH")
  }
  FGA.col <- barcolor[1]
  FGG.col <- barcolor[2]
  FGL.col <- barcolor[3]
  p1 <- ggplot(summaryFGA, aes(x = Subtype, y = mean, fill = rep("0", nrow(summaryFGA)))) + 
    geom_bar(stat = "identity") + 
    geom_errorbar(aes(ymax = mean + se, ymin = mean - se), position = position_dodge(0.9), width = 0.15) + 
    geom_text(aes(label = cut(FGA.test, c(0, 0.001, 0.01, 0.05, 0.1, 1), labels = c("****", "***", "**", "*", "."))), 
              x = n.moic / 2 + 0.5, 
              y = as.numeric(summaryFGA[which.max(summaryFGA$mean), "mean"]), 
              size = 8, angle = 90, fontface = "bold") + 
    scale_x_discrete(name = "", position = "top") + 
    theme_bw() + 
    theme(axis.line.y = element_line(linewidth = 0.8), 
          axis.ticks.y = element_line(linewidth = 0.2), 
          axis.text.y = element_blank(), 
          axis.title.x = element_text(vjust = -0.3, size = 12), 
          axis.text.x = element_text(size = 10, color = "black"), 
          plot.margin = unit(c(0.3, -1.7, 0.3, 0.3), "lines"), 
          legend.title = element_blank()) + 
    coord_flip() + 
    scale_fill_manual(values = FGA.col, breaks = c("0"), labels = c("Copy number-altered genome")) + 
    scale_y_reverse(expand = c(0.01, 0), name = "FGA (Fraction of Genome Altered)", position = "left")
  
  p2 <- ggplot(summaryFGGL, aes(x = Subtype, y = ifelse(class == "FGG", mean, -mean), fill = class)) + 
    geom_bar(stat = "identity") + 
    geom_errorbar(data = summaryFGGL[summaryFGGL$class == "FGG", ], 
                  aes(ymax = mean + se, ymin = mean - se), 
                  position = position_dodge(0.9), width = 0.15) + 
    geom_errorbar(data = summaryFGGL[summaryFGGL$class == "FGL", ], 
                  aes(ymax = -mean - se, ymin = -mean + se), 
                  position = position_dodge(0.9), width = 0.15) + 
    geom_text(aes(label = cut(FGL.test, c(0, 0.001, 0.01, 0.05, 0.1, 1), labels = c("****", "***", "**", "*", "."))), 
              x = n.moic / 2 + 0.5, 
              y = -as.numeric(summaryFGL[which.max(summaryFGL$mean), "mean"]), 
              size = 8, angle = 90, fontface = "bold") + 
    geom_text(aes(label = cut(FGG.test, c(0, 0.001, 0.01, 0.05, 0.1, 1), labels = c("****", "***", "**", "*", "."))), 
              x = n.moic / 2 + 0.5, 
              y = as.numeric(summaryFGG[which.max(summaryFGG$mean), "mean"]), 
              size = 8, angle = 90, fontface = "bold") + 
    scale_x_discrete(name = "") + 
    theme_bw() + 
    theme(axis.line.y = element_line(linewidth = 0.8), 
          axis.ticks.y = element_line(linewidth = 0.2), 
          axis.text.y = element_blank(), 
          axis.title.x = element_text(vjust = -0.3, size = 12), 
          axis.text.x = element_text(size = 10, color = "black"), 
          plot.margin = unit(c(0.3, 0.3, 0.3, -1), "lines"), 
          legend.title = element_blank()) + 
    coord_flip() + 
    scale_fill_manual(values = c(FGL.col, FGG.col), 
                      breaks = c("FGL", "FGG"), 
                      labels = c("Copy number-lost genome", "Copy number-gained genome")) + 
    scale_y_continuous(expand = c(0.01, 0), name = "FGL or FGG (Fraction of Genome Lost or Gained)")
  
  pp <- ggplot() + geom_label(data = summaryFGGL, aes(label = Subtype, 
                                                      x = Subtype, fill = Subtype), y = 0.5, color = "white", 
                              size = 0.9 * 11/.pt, hjust = 0.4, vjust = 0.5) + scale_fill_manual(values = clust.col) + 
    theme_minimal() + theme(axis.line.y = element_blank(), 
                            axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
                            axis.title.y = element_blank(), axis.title.x = element_blank(), 
                            plot.margin = unit(c(0.3, 0, 0.3, 0), "lines")) + guides(fill = "none") + 
    coord_flip() + scale_y_reverse()
  pal <- p1 + pp + p2 + plot_layout(widths = c(7, 1, 7), guides = "collect") & 
    theme(legend.position = "top")
  if (!is.null(title)) {
    pal = pal +
      theme(plot.title = element_text(face = "bold", size = 16, hjust = -0.5, vjust = 0))+
      labs(title = title)
  } else {
    pal = pal +
      theme(plot.title = element_text(face = "bold", size = 16, hjust = -0.5, vjust = 0))+
      labs(title = "FGA barplot")
  }
  if (is.null(fig.name)) {
    outFig <- "barplot of FGA.pdf"
  }
  else {
    outFig <- paste0(fig.name, ".pdf")
  }
  ggsave(file.path(fig.path, outFig), width = width, height = height)
  print(pal)
  if (n.moic > 2) {
    return(list(summary = outTab, FGA.p.value = FGA.test, 
                pairwise.FGA.test = pairwise.FGA.test, FGG.p.value = FGG.test, 
                pairwise.FGG.test = pairwise.FGG.test, FGL.p.value = FGL.test, 
                pairwise.FGL.test = pairwise.FGL.test, test.method = statistic))
  }
  else {
    return(list(summary = outTab, FGA.p.value = FGA.test, 
                FGG.p.value = FGG.test, FGL.p.value = FGL.test, 
                test.method = statistic))
  }
}

compFGA_optimized <- function(moic.res = NULL, segment = NULL, iscopynumber = FALSE, ga_column = NULL,
                              cnathreshold = 0.2, test.method = "nonparametric", barcolor = c("#008B8A", 
                                                                                              "#F2042C", "#21498D"), 
                              clust.col = c("#2EC4B6", "#E71D36", "#FF9F1C", "#BDD5EA", "#FFA5AB", 
                                            "#011627", "#023E8A", "#9D4EDD", "#f09c6c", "#09f3b3"), 
                              fig.path = getwd(), fig.name = NULL, width = 8, height = 4, 
                              prefix = "", title = NULL) {
  
  library(patchwork)
  library(dplyr)
  
  # Check required columns
  required_cols <- c("sample", "chrom", "start", "end", "value")
  if (!all(required_cols %in% colnames(segment))) {
    stop("segment data must have the following columns: sample, chrom, start, end, value.")
  }
  
  n.moic <- length(unique(moic.res$clust.res$clust))
  
  # Convert copy number to log2 ratio if needed
  if (iscopynumber) {
    segment$value <- log2(segment$value / 2)
  }
  
  # Match samples and subset data
  comsam <- intersect(moic.res$clust.res$samID, unique(segment$sample))
  clust.res <- moic.res$clust.res[comsam, , drop = FALSE]
  segment <- segment[segment$sample %in% comsam, ]
  
  # Add 'bases' column for segment length
  segment <- segment %>%
    mutate(bases = end - start)
  
  # Classification based on the given column or value thresholds
  if (!is.null(ga_column) && ga_column %in% colnames(segment)) {
    segment <- segment %>%
      mutate(classification = .data[[ga_column]])
    
    fga_summary <- segment %>%
      group_by(sample) %>%
      summarize(FGA = sum(bases[classification %in% c("gain", "loss")]) / sum(bases),
                FGG = sum(bases[classification == "gain"]) / sum(bases),
                FGL = sum(bases[classification == "loss"]) / sum(bases))
  } else {
    fga_summary <- segment %>%
      group_by(sample) %>%
      summarize(FGA = sum(bases[abs(value) > cnathreshold]) / sum(bases),
                FGG = sum(bases[value > cnathreshold]) / sum(bases),
                FGL = sum(bases[value < -cnathreshold]) / sum(bases))
  }
  
  # Error bars
  std <- function(x, na.rm = TRUE) {
    if (na.rm) {
      x <- as.numeric(na.omit(x))
      sd(x) / sqrt(length(x))
    } else {
      sd(x) / sqrt(length(x))
    }
  }
  
  # Add Subtype information
  outTab <- fga_summary %>%
    mutate(Subtype = paste0(prefix, clust.res[sample, "clust"]))
  
  # Summarize by Subtype
  summaryFGA <- outTab %>%
    group_by(Subtype) %>%
    summarize(mean = mean(FGA, na.rm = TRUE), se = std(FGA, na.rm = TRUE))
  
  summaryFGG <- outTab %>%
    group_by(Subtype) %>%
    summarize(mean = mean(FGG, na.rm = TRUE), se = std(FGG, na.rm = TRUE))
  
  summaryFGL <- outTab %>%
    group_by(Subtype) %>%
    summarize(mean = mean(FGL, na.rm = TRUE), se = std(FGL, na.rm = TRUE))
  
  summaryFGGL <- data.frame(rbind.data.frame(summaryFGG, summaryFGL), 
                            class = rep(c("FGG", "FGL"), c(nrow(summaryFGG), nrow(summaryFGL))), 
                            stringsAsFactors = FALSE)
  
  # Statistical tests
  if (n.moic == 2 & test.method == "nonparametric") {
    statistic <- "wilcox.test"
    FGA.test <- wilcox.test(outTab$FGA ~ outTab$Subtype)$p.value
    FGG.test <- wilcox.test(outTab$FGG ~ outTab$Subtype)$p.value
    FGL.test <- wilcox.test(outTab$FGL ~ outTab$Subtype)$p.value
  } else if (n.moic == 2 & test.method == "parametric") {
    statistic <- "t.test"
    FGA.test <- t.test(outTab$FGA ~ outTab$Subtype)$p.value
    FGG.test <- t.test(outTab$FGG ~ outTab$Subtype)$p.value
    FGL.test <- t.test(outTab$FGL ~ outTab$Subtype)$p.value
  } else if (n.moic > 2 & test.method == "nonparametric") {
    statistic <- "kruskal.test"
    FGA.test <- kruskal.test(outTab$FGA ~ outTab$Subtype)$p.value
    FGG.test <- kruskal.test(outTab$FGG ~ outTab$Subtype)$p.value
    FGL.test <- kruskal.test(outTab$FGL ~ outTab$Subtype)$p.value
  } else if (n.moic > 2 & test.method == "parametric") {
    statistic <- "anova"
    FGA.test <- summary(aov(outTab$FGA ~ outTab$Subtype))[[1]][["Pr(>F)"]][1]
    FGG.test <- summary(aov(outTab$FGG ~ outTab$Subtype))[[1]][["Pr(>F)"]][1]
    FGL.test <- summary(aov(outTab$FGL ~ outTab$Subtype))[[1]][["Pr(>F)"]][1]
  }
  
  # Plotting
  FGA.col <- barcolor[1]
  FGG.col <- barcolor[2]
  FGL.col <- barcolor[3]
  
  p1 <- ggplot(summaryFGA, aes(x = Subtype, y = mean, fill = rep("0", nrow(summaryFGA)))) + 
    geom_bar(stat = "identity") + 
    geom_errorbar(aes(ymax = mean + se, ymin = mean - se), position = position_dodge(0.9), width = 0.15) + 
    geom_text(aes(x = n.moic / 2 + 0.5, y = max(mean, na.rm = TRUE), 
                  label = cut(FGA.test, c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                              labels = c("****", "***", "**", "*", "."))), 
              size = 8, angle = 90, fontface = "bold") + 
    scale_x_discrete(name = "", position = "top") + 
    theme_bw() + 
    theme(axis.line.y = element_line(linewidth = 0.8), 
          axis.ticks.y = element_line(linewidth = 0.2), 
          axis.text.y = element_blank(), 
          axis.title.x = element_text(vjust = -0.3, size = 12), 
          axis.text.x = element_text(size = 10, color = "black"), 
          plot.margin = unit(c(0.3, -1.7, 0.3, 0.3), "lines"), 
          legend.title = element_blank()) + 
    coord_flip() + 
    scale_fill_manual(values = FGA.col, breaks = c("0"), labels = c("Copy number-altered genome")) + 
    scale_y_reverse(expand = c(0.01, 0), name = "FGA (Fraction of Genome Altered)", position = "left")
  
  p2 <- ggplot(summaryFGGL, aes(x = Subtype, y = ifelse(class == "FGG", mean, -mean), fill = class)) + 
    geom_bar(stat = "identity") + 
    geom_errorbar(data = summaryFGGL[summaryFGGL$class == "FGG", ], 
                  aes(ymax = mean + se, ymin = mean - se), position = position_dodge(0.9), width = 0.15) + 
    geom_errorbar(data = summaryFGGL[summaryFGGL$class == "FGL", ], 
                  aes(ymax = -mean - se, ymin = -mean + se), position = position_dodge(0.9), width = 0.15) + 
    geom_text(aes(x = n.moic / 2 + 0.5, y = max(mean, na.rm = TRUE), 
                  label = cut(FGG.test, c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                              labels = c("****", "***", "**", "*", "."))), 
              size = 8, angle = 90, fontface = "bold") + 
    geom_text(aes(x = n.moic / 2 + 0.5, y = -max(mean, na.rm = TRUE), 
                  label = cut(FGL.test, c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                              labels = c("****", "***", "**", "*", "."))), 
              size = 8, angle = 90, fontface = "bold") + 
    scale_x_discrete(name = "") + 
    theme_bw() + 
    theme(axis.line.y = element_line(linewidth = 0.8), 
          axis.ticks.y = element_line(linewidth = 0.2), 
          axis.text.y = element_blank(), 
          axis.title.x = element_text(vjust = -0.3, size = 12), 
          axis.text.x = element_text(size = 10, color = "black"), 
          plot.margin = unit(c(0.3, 0.3, 0.3, -1), "lines"), 
          legend.title = element_blank()) + 
    coord_flip() + 
    scale_fill_manual(values = c(FGL.col, FGG.col), 
                      breaks = c("FGL", "FGG"), 
                      labels = c("Copy number-lost genome", "Copy number-gained genome")) + 
    scale_y_continuous(expand = c(0.01, 0), name = "FGL or FGG (Fraction of Genome Lost or Gained)")
  
  # Combine plots
  pp <- ggplot() + geom_label(data = summaryFGGL, aes(label = Subtype, 
                                                      x = Subtype, fill = Subtype), y = 0.5, color = "white", 
                              size = 0.9 * 11/.pt, hjust = 0.4, vjust = 0.5) + 
    scale_fill_manual(values = clust.col) + 
    theme_minimal() + 
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), 
          axis.text.y = element_blank(), axis.title.y = element_blank(), 
          axis.title.x = element_blank(), plot.margin = unit(c(0.3, 0, 0.3, 0), "lines")) + 
    guides(fill = "none") + coord_flip() + scale_y_reverse()
  
  pal <- p1 + pp + p2 + plot_layout(widths = c(7, 1, 7), guides = "collect") & 
    theme(legend.position = "top")
  
  if (!is.null(title)) {
    pal <- pal + theme(plot.title = element_text(face = "bold", size = 16, hjust = -0.5, vjust = 0)) + 
      labs(title = title)
  } else {
    pal <- pal + theme(plot.title = element_text(face = "bold", size = 16, hjust = -0.5, vjust = 0)) + 
      labs(title = "FGA barplot")
  }
  
  if (is.null(fig.name)) {
    outFig <- "barplot of FGA.pdf"
  } else {
    outFig <- paste0(fig.name, ".pdf")
  }
  
  ggsave(file.path(fig.path, outFig), width = width, height = height)
  print(pal)
  
  if (n.moic > 2) {
    return(list(summary = outTab, FGA.p.value = FGA.test, FGG.p.value = FGG.test, FGL.p.value = FGL.test))
  } else {
    return(list(summary = outTab, FGA.p.value = FGA.test, FGG.p.value = FGG.test, FGL.p.value = FGL.test))
  }
}

# Prepare GSEA output for hierarchical clustering #####
prepare_gsea_output_for_hclust = function (gsea_output, dgea_output_name_style = "", 
                                           dea.method = "", mo.method = "",
                                           dat.path = "", dgea_padj_cutoff = 0.05,
                                           pathway_padj_cutoff = 0.05,
                                           logfc_cutoff = 0) {
  
  # Load corresponding DGEA files
  n.moic = length(gsea_output$gsea.list)
  DEpattern <- paste(mo.method, "_", ifelse(is.null(dgea_output_name_style), 
                                            "", paste0(dgea_output_name_style, "_")), dea.method, ".*._vs_Others.txt$", 
                     sep = "")
  DEfiles <- dir(dat.path, pattern = DEpattern)
  
  # Input validity checks
  if (length(DEfiles) == 0) {
    stop("no DEfiles!")
  }
  if (length(DEfiles) != n.moic) {
    stop("not all multi-omics clusters have DEfile!")
  }
  
  # DGEA input list
  dgea_sets = list()
  for (i in 1:length(DEfiles)) {
    dgea_sets[[i]] = data.table::fread(paste0(dat.path, "/", DEfiles[i]),
                                       header = TRUE, sep = "\t")
    # dgea_sets[[i]][, 2:ncol(dgea_sets[[i]])] = lapply(dgea_sets[[i]][, 2:ncol(dgea_sets[[i]])], as.numeric)
    # dgea_sets[[i]][, 1] = as.character(dgea_sets[[i]][, 1])
  }
  names(dgea_sets) = names(gsea_output$gsea.list)
  
  # Prepare output list
  hclust_input_dfs = list()
  
  for (i in 1:n.moic) {
    
    subtype_id = names(dgea_sets)[i]
    
    # Significant pathway results for a particular subtype
    result_object = gsea_output$gsea.list[[subtype_id]]@result %>%
      dplyr::filter(p.adjust < pathway_padj_cutoff)
    
    # signficantly deregulated genes
    sig_up = dgea_sets[[subtype_id]]$id[dgea_sets[[subtype_id]]$log2fc > abs(logfc_cutoff) & 
                                          dgea_sets[[subtype_id]]$padj < dgea_padj_cutoff]
    sig_down = dgea_sets[[subtype_id]]$id[dgea_sets[[subtype_id]]$log2fc < -abs(logfc_cutoff) &
                                            dgea_sets[[subtype_id]]$padj < dgea_padj_cutoff]
    
    # Prepare output data frame
    hclust_input_dfs[[subtype_id]] = data.frame(
      ID = result_object$ID,
      Down_regulated = sapply(result_object$ID, function(pathway) {
        
        # Filter genes in the pathway with logFC < 0 and p-values passing cutoff
        genes <- unlist(strsplit(gsea_output$gsea.list[[subtype_id]]@geneSets[[pathway]], "/"))
        down_genes <- intersect(genes, sig_down)
        paste(down_genes, collapse = ", ")
      }),
      Up_regulated = sapply(result_object$ID, function(pathway) {
        # Filter genes in the pathway with logFC > 0 and p-values passing cutoff
        genes <- unlist(strsplit(gsea_output$gsea.list[[subtype_id]]@geneSets[[pathway]], "/"))
        up_genes <- intersect(genes, sig_up)
        paste(up_genes, collapse = ", ")
      }),
      NES = result_object$NES,
      lowest_p = result_object$p.adjust
    )
  }
  
  return(hclust_input_dfs)
}

plot_pathway_heatmaps = function(gsea.lists, norm.expr = NULL, 
                                 representative = TRUE, moic.res = NULL,
                                 clust.col = c("#2EC4B6", "#E71D36", "#FF9F1C", "#BDD5EA", 
                                               "#FFA5AB", "#011627", "#023E8A", "#9D4EDD", "#f09c6c", "#09f3b3"),
                                 present_clusters = NULL, # New argument for cluster labels
                                 subtype_prefix = "CS", n.path = 10, msigdb.path = NULL,
                                 norm.method = "mean", dirct = NULL, color = NULL,
                                 fig.name = NULL, fig.path = getwd(), width = 15, height = 10, name = NULL,
                                 gsva.method = "gsva") {
  
  if (!(is.logical(representative) && length(representative) == 1)) {
    stop("The 'representative' argument must either be TRUE or FALSE")
  }
  
  if (length(present_clusters) != length(gsea.lists)) {
    gsea.lists = gsea.lists[paste0(dirct, "_", present_clusters)]
  }
  
  if (representative) {
    fig.name = paste0("representative_", fig.name)

    gsea.lists = lapply(gsea.lists, function (x) {
      x = x[["clustered_df"]] %>%
        dplyr::filter(Status == "Representative")
    })
    
    if (dirct == "up") {
      gsea.lists = lapply(gsea.lists, function (x) {
        x = x %>%
          dplyr::arrange(lowest_p, desc(NES))
      })
    } else {
      gsea.lists = lapply(gsea.lists, function (x) {
        x = x %>%
          dplyr::arrange(lowest_p, NES)
      })
    }
  }
  
  msigdb <- try(clusterProfiler::read.gmt(msigdb.path), 
                silent = TRUE)
  if (class(msigdb) == "try-error") {
    stop("please provide correct ABSOLUTE PATH for MSigDB file.")
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
  
  # If present_clusters is NULL, use all clusters
  if (is.null(present_clusters)) {
    present_clusters <- paste0(subtype_prefix, 1:length(clust.col))
  }
  
  # Extract cluster numbers from present_clusters
  cluster_numbers <- as.numeric(gsub("[^0-9]", "", present_clusters))
  
  # Filter colors based on the extracted cluster numbers
  colvec <- clust.col[cluster_numbers]
  names(colvec) <- present_clusters
  
  # Filter moic.res and gsea.lists to include only present clusters
  moic.res$clust.res <- moic.res$clust.res[moic.res$clust.res$clust %in% cluster_numbers, ]
  gsea.lists <- gsea.lists[sapply(gsea.lists, function(x) any(as.numeric(gsub("[^0-9]", "", x$Cluster)) %in% cluster_numbers))]
  
  n.moic <- length(cluster_numbers)
  
  pathway <- pathcore <- list()
  pathnum <- c()
  for (filek in 1:length(gsea.lists)) {
    if (nrow(gsea.lists[[filek]]) > n.path) {
      pathway[[filek]] <- gsea.lists[[filek]][1:n.path, ]
    }
    else {
      pathway[[filek]] <- gsea.lists[[filek]]
    }
    
    pathnum <- c(pathnum, nrow(pathway[[filek]]))
    pathway[[filek]]$dirct <- dirct
    for (i in pathway[[filek]]$ID) {
      pathcore[[i]] <- msigdb[which(msigdb[, 1] %in% 
                                      i), "gene"]
    }
  }
  
  sam.order <- moic.res$clust.res[order(moic.res$clust.res$clust, 
                                        decreasing = FALSE), "samID"]
  annCol <- data.frame(Subtype = paste0(subtype_prefix, moic.res$clust.res[sam.order, 
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
    for (i in present_clusters) {
      esm <- cbind.data.frame(esm, data.frame(rowmean(es[, 
                                                         rownames(annCol[which(annCol$Subtype == i), 
                                                                         , drop = FALSE])])))
    }
  }
  if (norm.method == "median") {
    for (i in present_clusters) {
      esm <- cbind.data.frame(esm, data.frame(rowmedian(es[, 
                                                           rownames(annCol[which(annCol$Subtype == i), 
                                                                           , drop = FALSE])])))
    }
  }
  colnames(esm) <- present_clusters
  annRow <- data.frame(Subtype = rep(present_clusters, 
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
  return(list(gsea.list = gsea.lists, raw.es = es.backup, scaled.es = es, 
              grouped.es = esm, heatmap = hm))
}

# Silhouette
getSilhouette_mod = function (sil = NULL, 
                              clust.col = c("#2EC4B6", "#E71D36", "#FF9F1C", 
                                            "#BDD5EA", "#FFA5AB", "#011627", 
                                            "#023E8A", "#9D4EDD", "#f09c6c", 
                                            "#09f3b3"), 
                              fig.path = getwd(), 
                              fig.name = "silhouette", 
                              width = 5.5, 
                              height = 5,
                              axis_label_size = 1, 
                              axis_label_font = 1, 
                              annotation_size = 1.2,
                              title_text_size = 1
) 
{
  N.clust <- length(unique(sil[, 1]))
  colvec <- clust.col[1:N.clust]
  outFig <- paste0(fig.name, ".pdf")
  
  # Set graphical parameters, including text sizes and fonts
  par(bty = "o", mgp = c(2.5, 0.33, 0), 
      mar = c(5.1, 2.1, 3.1, 2.1) + 0.1, 
      las = 1, 
      tcl = -0.25, 
      cex.axis = axis_label_size,    
      cex.lab = axis_label_size,     
      font.lab = axis_label_font,    
      cex.main = title_text_size
  )
  
  # Plot the silhouette with the specified colors
  plot(sil, border = NA, col = colvec, cex = annotation_size)
  
  # Save the plot to a PDF
  dev.copy2pdf(file = file.path(fig.path, outFig), width = width, height = height)
}

getSilhouette_ggplot = function(sil = NULL, 
                                clust.col = c("#2EC4B6", "#E71D36", "#FF9F1C", 
                                              "#BDD5EA", "#FFA5AB", "#011627", 
                                              "#023E8A", "#9D4EDD", "#f09c6c", 
                                              "#09f3b3"), 
                                fig.path = getwd(), 
                                fig.name = "silhouette_ggplot", 
                                width = 7,
                                height = 5,
                                axis_label_size = 12,
                                axis_label_font = "plain",
                                text_size = 3.5,
                                title_size = 16,
                                algorithm = "",
                                save_plot = TRUE
) 
{
  library(ggplot2)
  
  # Prepare the data frame from silhouette object
  sil_df <- as.data.frame(sil)
  sil_df$samID <- rownames(sil_df)
  
  # Calculate the average silhouette width for each cluster
  cluster_stats <- aggregate(sil_df$sil_width, by = list(sil_df$cluster), FUN = mean)
  colnames(cluster_stats) <- c("cluster", "Average_Silhouette_Width")
  
  # Add cluster stats to the sil_df for annotation
  sil_df <- merge(sil_df, cluster_stats, by.x = "cluster", by.y = "cluster")
  
  # Create a new column for prefixed cluster labels using the 'algorithm' argument
  sil_df$prefixed_cluster <- paste0(algorithm, sil_df$cluster)
  
  # Assign colors to the prefixed clusters
  unique_prefixed_clusters <- unique(sil_df$prefixed_cluster)
  cluster_colors <- setNames(clust.col[1:length(unique_prefixed_clusters)], unique_prefixed_clusters)
  
  # Create a separate data frame just for the average silhouette widths
  avg_sil_text_df <- unique(sil_df[, c("cluster", "prefixed_cluster", "Average_Silhouette_Width")])
  
  # Create the ggplot2 silhouette plot
  p <- ggplot(sil_df, aes(x = sil_width, y = reorder(samID, sil_width), fill = prefixed_cluster)) +
    geom_bar(stat = "identity", size = 0.25, show.legend = FALSE) +  # Remove legend
    scale_fill_manual(values = cluster_colors) +  # Use colors for prefixed clusters
    facet_wrap(~prefixed_cluster, scales = "free_y", ncol = 1) +  # Use the prefixed cluster label in facets
    labs(title = "Silhouette Plot", x = "Silhouette Width", y = "") +  # X-axis title
    theme_minimal() +  # White background, minimal theme
    theme(
      axis.text.x = element_text(size = axis_label_size),
      axis.text.y = element_blank(),  # Remove sample names from y-axis
      axis.ticks.y = element_blank(), # Remove y-axis ticks
      axis.title.x = element_text(size = axis_label_size, face = axis_label_font),  # Bold axis label
      axis.line.x = element_line(linewidth = 0.35),
      axis.ticks.x = element_line(linewidth = 0.15),
      plot.title = element_text(size = title_size, hjust = 0.5, face = "bold"),
      strip.text = element_text(size = 14),  # Prefixed cluster labels in facets
      panel.grid.major = element_blank(),  # Remove gridlines
      panel.grid.minor = element_blank(),  # Remove minor gridlines
      panel.background = element_rect(fill = "white", color = NA)  # Ensure white background
    ) +
    # Print average silhouette width once per facet
    geom_text(data = avg_sil_text_df, aes(x = 0, y = 0,  # Adjust vertical position
                                          label = paste("Avg. sil =", round(Average_Silhouette_Width, 2))),
              inherit.aes = FALSE, hjust = 2, vjust = -3, size = axis_label_size / 4, color = "black", fontface = "bold")
  
  # Display plot in the RStudio viewer
  print(p)
  
  # Optionally save the plot as a PDF
  if (save_plot) {
    ggsave(filename = file.path(fig.path, paste0(fig.name, ".pdf")), plot = p, width = width, height = height)
  }
  
  return(p)
}
