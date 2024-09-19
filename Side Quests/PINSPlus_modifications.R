SubtypingOmicsData_mod = function (dataList, kMin = 2, kMax = 5, k = NULL, agreementCutoff = 0.5, 
                                   ncore = 1, verbose = T, sampledSetSize = 2000, knn.k = NULL,
                                   binary_flags = rep("No", length(dataList)), binary_distance = NULL,
                                   ...)
{
  now = Sys.time()
  mlog <- if (!verbose) 
    function(...) {
    }
  else function(...) {
    message(...)
    flush.console()
  }
  dataListComplete <- dataList
  
  # Check if there are binary_modalities
  if ("Yes" %in% binary_flags) {
    binary_modalities_indices = which(binary_flags == "Yes")
    binary_modalities = names(dataList)[binary_modalities_indices]
    continuous_modalities_indices <- which(binary_flags == "No")
    continuous_modalities = names(dataList)[continuous_modalities_indices]
  }
  
  commonSamples <- Reduce(f = intersect, x = lapply(dataList, 
                                                    rownames))
  dataList <- lapply(dataList, function(d) d[commonSamples, 
  ])
  notCommonData <- lapply(dataListComplete, function(d) {
    rn <- rownames(d)[(!rownames(d) %in% commonSamples)]
    d <- matrix(d[rn, ], ncol = ncol(d))
    rownames(d) <- rn
    d
  })
  dataListTrain <- NULL
  dataListTest <- NULL
  seed = round(rnorm(1) * 10^6)
  dataList <- lapply(dataList, as.data.frame)
  if (nrow(dataList[[1]]) > sampledSetSize) {
    n_samples <- nrow(dataList[[1]])
    ind <- sample.int(n_samples, size = sampledSetSize)
    dataListTrain <- lapply(dataList, function(x) x[ind, 
    ])
    dataListTest <- lapply(dataList, function(x) x[-ind, 
                                                   , drop = F])
    dataList <- dataListTrain
  }
  runPerturbationClustering_mod <- function(dataList, kMin, kMax, 
                                            stage = 1, forceSplit = FALSE, k = NULL) {
    dataTypeResult <- lapply(seq_along(dataList), function(idx) {
      data <- as.matrix(dataList[[idx]])
      binary_check <- ifelse(binary_flags[idx] == "Yes", TRUE, FALSE)
      
      # Perform clustering with binary or continuous data based on the flag
      set.seed(seed)
      data <- data[rowSums(is.na(data)) == 0, ]
      PerturbationClustering_mod(data, kMin, kMax, ncore = ncore, 
                                 binary = binary_check, verbose = verbose, ...)
    })
    allSamples <- unique(unlist(lapply(dataList, rownames)))
    origList <- lapply(dataTypeResult, function(r) r$origS[[r$k]])
    origMerged <- do.call(what = rbind, args = lapply(origList, 
                                                      function(o) {
                                                        o <- as.data.frame(o)
                                                        as.numeric(as.matrix(t(as.data.frame(t(o[allSamples, 
                                                        ]))[allSamples, ])))
                                                      }))
    orig = matrix(colMeans(origMerged, na.rm = T), nrow = length(allSamples))
    rownames(orig) <- colnames(orig) <- allSamples
    orig <- impute::impute.knn(orig)$data
    orig[is.na(orig)] <- 0
    PW = matrix(as.numeric(colSums(origMerged == 0, na.rm = T) == 
                             0), nrow = length(allSamples))
    rownames(PW) <- colnames(PW) <- allSamples
    agreement = (sum(orig == 0) + sum(orig == 1) - nrow(orig))/(nrow(orig)^2 - 
                                                                  nrow(orig))
    pertList <- lapply(dataTypeResult, function(r) r$pertS[[r$k]])
    pertMerged <- do.call(what = rbind, args = lapply(pertList, 
                                                      function(p) {
                                                        p <- as.data.frame(p)
                                                        as.numeric(as.matrix(t(as.data.frame(t(p[allSamples, 
                                                        ]))[allSamples, ])))
                                                      }))
    pert = matrix(colMeans(pertMerged, na.rm = T), nrow = length(allSamples))
    rownames(pert) <- colnames(pert) <- allSamples
    pert <- impute::impute.knn(pert)$data
    pert[is.na(pert)] <- 0
    groups <- NULL
    mlog("STAGE : ", stage, "\t Agreement : ", agreement)
    if (agreement >= agreementCutoff | forceSplit) {
      hcW <- hclust(dist(PW))
      maxK = min(kMax * 2, dim(unique(PW, MARGIN = 2))[2] - 
                   (stage - 1))
      maxHeight = FindMaxHeight(hcW, maxK = min(2 * maxK, 
                                                10))
      groups <- cutree(hcW, maxHeight)
      if (!is.null(k) && max(groups) > k) {
        groups <- cutree(hcW, k)
      }
    }
    list(dataTypeResult = dataTypeResult, orig = orig, pert = pert, 
         PW = PW, groups = groups, agreement = agreement)
  }
  pResult <- runPerturbationClustering_mod(dataList, kMin, kMax, 
                                           k = k)
  groups <- pResult$groups
  groups2 <- NULL
  if (!is.null(groups)) {
    groups2 <- groups
    if (is.null(k)) {
      for (g in sort(unique(groups))) {
        miniGroup <- names(groups[groups == g])
        if (length(miniGroup) > 30) {
          groupsM <- runPerturbationClustering_mod(dataList = lapply(dataList, 
                                                                     function(d) d[miniGroup, ]), kMin = kMin, 
                                                   kMax = min(kMax, 5), stage = 2)$groups
          if (!is.null(groupsM)) 
            groups2[miniGroup] <- paste(g, groupsM, 
                                        sep = "-")
        }
      }
    }
    else {
      agreements <- rep(1, length(groups))
      names(agreements) <- names(groups)
      tbl <- sort(table(groups2), decreasing = T)
      minGroupSize <- 30
      while (length(unique(groups2)) < k) {
        if (all(tbl <= 30)) {
          minGroupSize <- 10
        }
        if (all(tbl <= 10)) 
          (break)()
        for (g in names(tbl)) {
          miniGroup <- names(groups2[groups2 == g])
          if (length(miniGroup) > minGroupSize) {
            splitRes <- runPerturbationClustering_mod(dataList = lapply(dataList, 
                                                                        function(d) d[miniGroup, ]), kMin = kMin, 
                                                      kMax = min(kMax, 5), stage = 2, forceSplit = T)
            groupsM <- splitRes$groups
            if (!is.null(groupsM)) {
              groups2[miniGroup] <- paste(g, groupsM, 
                                          sep = "-")
              agreements[miniGroup] <- splitRes$agreement
            }
          }
        }
        tbl <- sort(table(groups2), decreasing = T)
      }
      agreements.unique = unique(agreements)
      for (aggr in sort(unique(agreements))) {
        if (length(unique(groups2)) == k) 
          (break)()
        merge.group <- agreements == aggr
        k.smallGroup <- length(unique(groups2[merge.group]))
        k.need <- k - (length(unique(groups2)) - k.smallGroup + 
                         1) + 1
        groups2[merge.group] <- unlist(lapply(strsplit(groups2[merge.group], 
                                                       "-"), function(g) {
                                                         paste0(g[1:(length(g) - 1)], collapse = "-")
                                                       }))
        if (k.need > 1) {
          splitRes <- runPerturbationClustering_mod(dataList = lapply(dataList, 
                                                                      function(d) d[miniGroup, ]), kMin = kMin, 
                                                    kMax = min(kMax, 5), stage = 2, forceSplit = T, 
                                                    k = k.need)
          groupsM <- splitRes$groups
          groups2[merge.group] <- paste(groups2[merge.group], 
                                        groupsM, sep = "-")
        }
        agreements[merge.group] <- 1
      }
    }
  }
  else {
    set.seed(seed)
    orig <- pResult$orig
    dataTypeResult <- pResult$dataTypeResult
    clusteringAlgorithm = GetClusteringAlgorithm(...)$fun
    groupings <- lapply(dataTypeResult, function(r) clusteringAlgorithm(data = r$origS[[r$k]], 
                                                                        k = r$k))
    pGroups <- ClusterUsingPAM(orig = orig, kMax = kMax * 
                                 2, groupings = groupings)
    hGroups <- ClusterUsingHierarchical(orig = orig, kMax = kMax * 
                                          2, groupings = groupings)
    pAgree = pGroups$agree
    hAgree = hGroups$agree
    groups <- (if (pAgree > hAgree) 
      pGroups
      else if (hAgree > pAgree) 
        hGroups
      else {
        pAgree = ClusterUsingPAM(orig = pResult$pert, kMax = kMax, 
                                 groupings = groupings)$agree
        hAgree = ClusterUsingHierarchical(orig = pResult$pert, 
                                          kMax = kMax, groupings = groupings)$agree
        if (hAgree - pAgree >= 0.001) 
          hGroups
        else pGroups
      })$cluster
    names(groups) <- rownames(orig)
    groups2 <- groups
    if (is.null(k)) {
      mlog("Check if can proceed to stage II")
      normalizedEntropy = entropy::entropy(table(groups))/log(length(unique(groups)), 
                                                              exp(1))
      if (normalizedEntropy < 0.5) {
        for (g in sort(unique(groups))) {
          miniGroup <- names(groups[groups == g])
          if (length(miniGroup) > 30) {
            groupsM <- runPerturbationClustering_mod(dataList = lapply(dataList, 
                                                                       function(d) d[miniGroup, ]), kMin = kMin, 
                                                     kMax = min(kMax, 5), stage = 2, forceSplit = T, 
                                                     k = NULL)$groups
            if (!is.null(groupsM)) 
              groups2[miniGroup] <- paste(g, groupsM, 
                                          sep = "-")
          }
        }
      }
    }
    else {
      if (length(unique(groups2)) > k) {
        pGroups <- ClusterUsingPAM(orig = orig, kMax = kMax * 
                                     2, groupings = groupings, k)
        hGroups <- ClusterUsingHierarchical(orig = orig, 
                                            kMax = kMax * 2, groupings = groupings, k)
        pAgree = pGroups$agree
        hAgree = hGroups$agree
        groups <- (if (pAgree > hAgree) 
          pGroups
          else if (hAgree > pAgree) 
            hGroups
          else {
            pAgree = ClusterUsingPAM(orig = pResult$pert, 
                                     kMax = kMax, groupings = groupings, k)$agree
            hAgree = ClusterUsingHierarchical(orig = pResult$pert, 
                                              kMax = kMax, groupings = groupings, k)$agree
            if (hAgree - pAgree >= 0.001) 
              hGroups
            else pGroups
          })$cluster
        names(groups) <- rownames(orig)
        groups2 <- groups
      }
      else if (length(unique(groups2)) < k) {
        normalizedEntropy = entropy::entropy(table(groups))/log(length(unique(groups)), 
                                                                exp(1))
        agreements <- rep(1, length(groups))
        names(agreements) <- names(groups)
        if (normalizedEntropy < 0.5) {
          for (g in sort(unique(groups))) {
            miniGroup <- names(groups[groups == g])
            if (length(miniGroup) > 30) {
              splitRes <- runPerturbationClustering_mod(dataList = lapply(dataList, 
                                                                          function(d) d[miniGroup, ]), kMin = kMin, 
                                                        kMax = min(kMax, 5), stage = 2, forceSplit = T)
              groupsM <- splitRes$groups
              if (!is.null(groupsM)) {
                groups2[miniGroup] <- paste(g, groupsM, 
                                            sep = "-")
                agreements[miniGroup] <- splitRes$agreement
              }
            }
          }
        }
        if (length(unique(groups2)) < k) {
          tbl <- sort(table(groups2), decreasing = T)
          minGroupSize <- 30
          while (length(unique(groups2)) < k) {
            if (all(tbl <= 30)) {
              minGroupSize <- 10
            }
            if (all(tbl <= 10)) 
              (break)()
            for (g in names(tbl)) {
              miniGroup <- names(groups2[groups2 == 
                                           g])
              if (length(miniGroup) > minGroupSize) {
                splitRes <- runPerturbationClustering_mod(dataList = lapply(dataList, 
                                                                            function(d) d[miniGroup, ]), kMin = kMin, 
                                                          kMax = min(kMax, 5), stage = 2, forceSplit = T)
                groupsM <- splitRes$groups
                if (!is.null(groupsM)) {
                  groups2[miniGroup] <- paste(g, groupsM, 
                                              sep = "-")
                  agreements[miniGroup] <- splitRes$agreement
                }
              }
            }
            tbl <- sort(table(groups2), decreasing = T)
          }
        }
        agreements.unique = unique(agreements)
        for (aggr in sort(unique(agreements))) {
          if (length(unique(groups2)) == k) 
            (break)()
          merge.group <- agreements == aggr
          k.smallGroup <- length(unique(groups2[merge.group]))
          k.need <- k - (length(unique(groups2)) - k.smallGroup + 
                           1) + 1
          groups2[merge.group] <- unlist(lapply(strsplit(groups2[merge.group], 
                                                         "-"), function(g) {
                                                           paste0(g[1:(length(g) - 1)], collapse = "-")
                                                         }))
          if (k.need > 1) {
            splitRes <- runPerturbationClustering_mod(dataList = lapply(dataList, 
                                                                        function(d) d[merge.group, ]), kMin = kMin, 
                                                      kMax = min(kMax, 5), stage = 2, forceSplit = T, 
                                                      k = k.need)
            groupsM <- splitRes$groups
            groups2[merge.group] <- paste(groups2[merge.group], 
                                          groupsM, sep = "-")
          }
          agreements[merge.group] <- 1
        }
      }
    }
  }
  {
    train_y <- groups
    train_y2 <- groups2
    if (!is.null(dataListTest)) {
      set.seed(seed)
      RcppParallel::setThreadOptions(ncore)
      test_prob <- matrix(0, nrow = n_samples - sampledSetSize, 
                          ncol = length(unique(groups)))
      if (!is.null(train_y2)) {
        test_prob2 <- matrix(0, nrow = n_samples - sampledSetSize, 
                             ncol = length(unique(groups2)))
      }
      for (i in 1:length(dataListTrain)) {
        train <- dataListTrain[[i]]
        test <- dataListTest[[i]]
        
        if (binary_flags[i] == TRUE) {
          if (ncol(train) * nrow(train) > 2e+07) {
            pca <- list()
            pca$x <- cmdscale(d = stats::dist(train, method = binary_distance),
                              k = min(nrow(train), 20))
            colnames(pca$x) <- paste0("PC", 1:ncol(pca$x))
          } else {
            pca <- list()
            pca$x <- cmdscale(d = stats::dist(train, method = binary_distance),
                              k = min(nrow(train), 200))
            colnames(pca$x) <- paste0("PC", 1:ncol(pca$x))
          }
          
          # Project test data into the MDS space
          test_distances <- as.matrix(stats::dist(rbind(train, test), method = binary_distance))
          test_train_distances <- test_distances[(nrow(train) + 1):nrow(test_distances), 1:nrow(train)]
          
          # Use Gower interpolation to project test data
          test_coords <- mds_project(train_coords = pca$x, 
                                     train_distances = as.matrix(stats::dist(train, method = binary_distance)),
                                     test_distances = test_train_distances)
          
          # Use `test_coords` for further analysis
          test <- test_coords
        } else {
          # Handle PCA as usual for continuous data
          if (ncol(train) * nrow(train) > 2e+07) {
            pca <- rpca.para(train, min(nrow(train), 20), scale = F, p = 10)
          } else {
            pca <- prcomp(train, rank. = min(nrow(train), 200))
          }
          test <- predict.rpca.para(pca, test)
        }
        train <- pca$x
        test_prob <- test_prob + classifierProb(train, 
                                                groups, test, knn.k)
        if (!is.null(train_y2)) {
          test_prob2 <- test_prob2 + classifierProb(train, 
                                                    groups2, test, knn.k)
        }
      }
      test_y <- apply(test_prob, 1, which.max)
      groups <- rep(0, n_samples)
      groups[ind] <- train_y
      groups[-ind] <- test_y
      if (!is.null(train_y2)) {
        test_y2 <- apply(test_prob2, 1, which.max)
        groups2 <- rep(0, n_samples)
        groups2[ind] <- train_y2
        groups2[-ind] <- test_y2
      }
    }
  }
  {
    train_y <- groups
    train_y2 <- groups2
    allSamples <- unique(unlist(lapply(dataListComplete, 
                                       rownames)))
    n_samples <- length(allSamples)
    notCommonSamples <- allSamples[!(allSamples %in% commonSamples)]
    if (length(allSamples) > length(commonSamples)) {
      set.seed(seed)
      RcppParallel::setThreadOptions(ncore)
      if (!is.null(dataListTrain)) {
        dataListTrain <- lapply(1:length(dataList), 
                                function(i) {
                                  as.matrix(rbind(dataList[[i]], dataListTrain[[i]]))
                                })
      }
      else {
        dataListTrain <- dataList
      }
      dataListTest <- notCommonData
      test_prob <- matrix(0, nrow = length(notCommonSamples), 
                          ncol = length(unique(groups)))
      rownames(test_prob) <- notCommonSamples
      colnames(test_prob) <- unique(groups)
      if (!is.null(train_y2)) {
        test_prob2 <- matrix(0, nrow = length(notCommonSamples), 
                             ncol = length(unique(groups2)))
        rownames(test_prob2) <- notCommonSamples
        colnames(test_prob2) <- unique(groups2)
      }
      for (i in 1:length(dataListTrain)) {
        train <- dataListTrain[[i]]
        test <- dataListTest[[i]]
        if (nrow(test) == 0) 
          (next)()
        if (binary_flags[i] == "Yes") {
          # Handling binary data using MDS
          pca <- list()
          pca$x <- cmdscale(d = stats::dist(train, method = binary_distance), 
                            k = min(nrow(train), 20))
          colnames(pca$x) <- paste0("PC", 1:ncol(pca$x))
          
          # Project test data into the MDS space
          test_distances <- as.matrix(stats::dist(rbind(train, test), method = binary_distance))
          test_train_distances <- test_distances[(nrow(train) + 1):nrow(test_distances), 1:nrow(train)]
          
          # Use Gower interpolation to project test data
          test_coords <- mds_project(train_coords = pca$x, 
                                     train_distances = as.matrix(stats::dist(train, method = binary_distance)),
                                     test_distances = test_train_distances)
          
          test <- test_coords  # Use projected test coordinates
          
        } else {
          # Handling continuous data using PCA
          pca <- rpca.para(train, min(nrow(train), 20), scale = F)
          test <- predict.rpca.para(pca, test)  # Project test data for continuous case
        }
        
        train <- pca$x  # Set train to transformed PCA or MDS coordinates
        test_prob_tmp <- classifierProb(train, groups, 
                                        test, knn.k)
        test_prob[rownames(test_prob_tmp), ] <- test_prob[rownames(test_prob_tmp), 
        ] + test_prob_tmp
        if (!is.null(train_y2)) {
          test_prob_tmp2 <- classifierProb(train, groups2, 
                                           test, knn.k)
          test_prob2[rownames(test_prob_tmp2), ] <- test_prob2[rownames(test_prob_tmp2), 
          ] + test_prob_tmp2
        }
      }
      test_y <- colnames(test_prob)[apply(test_prob, 1, 
                                          which.max)]
      groups <- rep(0, n_samples)
      names(groups) <- allSamples
      groups[commonSamples] <- train_y
      groups[notCommonSamples] <- test_y
      if (!is.null(train_y2)) {
        test_y2 <- colnames(test_prob2)[apply(test_prob2, 
                                              1, which.max)]
        groups2 <- rep(0, n_samples)
        names(groups2) <- allSamples
        groups2[commonSamples] <- train_y2
        groups2[notCommonSamples] <- test_y2
      }
    }
  }
  timediff = Sys.time() - now
  mlog("Done in ", timediff, " ", units(timediff), ".\n")
  list(cluster1 = groups, cluster2 = groups2, dataTypeResult = pResult$dataTypeResult)
}

PerturbationClustering_mod = function (data, kMin = 2, kMax = 5, k = NULL, verbose = T, ncore = 1, 
                                       clusteringMethod = "kmeans", clusteringFunction = NULL, 
                                       clusteringOptions = NULL, perturbMethod = "noise", perturbFunction = NULL, 
                                       perturbOptions = NULL, PCAFunction = NULL, iterMin = 20, 
                                       iterMax = 200, madMin = 0.001, msdMin = 1e-06, sampledSetSize = 2000, 
                                       knn.k = NULL, binary = FALSE) 
{
  RcppParallel::setThreadOptions(ncore)
  if (nrow(data) <= sampledSetSize) {
    now = Sys.time()
    log <- if (!verbose) 
      function(...) {
      }
    else function(...) {
      message(...)
      flush.console()
    }
    clusteringAlgorithm = GetClusteringAlgorithm(clusteringMethod = clusteringMethod, 
                                                 clusteringFunction = clusteringFunction, clusteringOptions = clusteringOptions)
    perturbationAlgorithm = GetPerturbationAlgorithm(data = data, 
                                                     perturbMethod = perturbMethod, perturbFunction = perturbFunction, 
                                                     perturbOptions = perturbOptions)
    log("Clustering method: ", clusteringAlgorithm$name)
    log("Perturbation method: ", perturbationAlgorithm$name)
    seed = round(rnorm(1) * 10^6)
    if (is.null(PCAFunction)) {
      if (binary) {
        if (ncol(train) * nrow(train) > 2e+07) {
          pca = list()
          pca$x = cmdscale(d = stats::dist(train, method = binary_distance),
                           k = min(nrow(train), 20))
          colnames(pca$x) = paste0("PC", 1:ncol(pca$x))
        }
        else {
          pca = list()
          pca$x = cmdscale(d = stats::dist(train, method = binary_distance),
                           k = min(nrow(train), 200))
          colnames(pca$x) = paste0("PC", 1:ncol(pca$x))
        }
      } else {
        if (ncol(train) * nrow(train) > 2e+07) {
          pca <- rpca.para(train, min(nrow(train), 20), 
                           scale = F, p = 10)
        }
        else {
          pca <- prcomp(train, rank. = min(nrow(train), 
                                           200))
        }
      }
    }
    else {
      pca <- list(x = PCAFunction(train))
    }
    if (!is.null(k) || kMin == kMax) {
      if (is.null(k)) {
        k = kMin
      }
      set.seed(seed)
      cluster <- kmeans(pca$x, centers = k, iter.max = 1000, 
                        nstart = 1000)$cluster
      return(list(k = k, cluster = cluster, origS = NULL, 
                  pertS = NULL, Discrepancy = NULL, pca = pca))
    }
    log("Building original connectivity matrices")
    set.seed(seed)
    origPartition <- GetOriginalSimilarity(data = pca$x, 
                                           clusRange = kMin:kMax, clusteringAlgorithm = clusteringAlgorithm$fun, 
                                           showProgress = verbose, ncore)
    origS <- origPartition$origS
    listAUC <- list()
    for (k in kMin:kMax) {
      listAUC[k] <- list(c())
    }
    stoppingCriteriaHandler <- DiffDevianceStoppingHandler(kMin = kMin, 
                                                           kMax = kMax, origS = origS, iterMin = iterMin, madMin = madMin, 
                                                           msdMin = msdMin, onExcute = function(k, AUCs) {
                                                             listAUC[[k]] <<- AUCs
                                                           })
    log("\nBuilding perturbed connectivity matrices")
    set.seed(seed)
    pertS <- GetPerturbedSimilarity(data = pca$x, clusRange = kMin:kMax, 
                                    iterMax = iterMax, iterMin = iterMin, origS = origS, 
                                    clusteringAlgorithm = clusteringAlgorithm$fun, perturbedFunction = perturbationAlgorithm$fun, 
                                    stoppingCriteriaHandler = stoppingCriteriaHandler, 
                                    showProgress = verbose, ncore = ncore)
    for (k in kMin:kMax) {
      AUCs <- listAUC[[k]]
      pert <- matrix(0, nrow(data), nrow(data))
      for (i in 1:length(AUCs)) {
        pert <- pert + pertS[[k]][[i]] * length(which(AUCs == 
                                                        AUCs[i]))
      }
      pertS[[k]] <- pert/sum(as.matrix(table(AUCs))[, 
                                                    1]^2)
    }
    Discrepancy <- CalcPerturbedDiscrepancy(origS, pertS, 
                                            clusRange = kMin:kMax)
    Discrepancy$AUC = round(Discrepancy$AUC, digits = 4)
    clus <- min(which(Discrepancy$AUC == max(Discrepancy$AUC[kMin:kMax])))
    timediff = Sys.time() - now
    log("Done in ", timediff, " ", units(timediff), ".\n")
    list(k = clus, cluster = origPartition$groupings[[clus]], 
         origS = origS, pertS = pertS, Discrepancy = Discrepancy, 
         pca = pca)
  }
  else {
    log <- if (!verbose) 
      function(...) {
      }
    else function(...) {
      message(...)
      flush.console()
    }
    log("\nUsing knn...")
    now = Sys.time()
    names <- rownames(data)
    set.seed(1)
    ind <- sample.int(nrow(data), sampledSetSize)
    train <- data[ind, ]
    test <- data[-ind, , drop = F]
    clusteringAlgorithm = GetClusteringAlgorithm(clusteringMethod = clusteringMethod, 
                                                 clusteringFunction = clusteringFunction, clusteringOptions = clusteringOptions)
    perturbationAlgorithm = GetPerturbationAlgorithm(data = data, 
                                                     perturbMethod = perturbMethod, perturbFunction = perturbFunction, 
                                                     perturbOptions = perturbOptions)
    log("Clustering method: ", clusteringAlgorithm$name)
    log("Perturbation method: ", perturbationAlgorithm$name)
    seed = round(rnorm(1) * 10^6)
    if (is.null(PCAFunction)) {
      if (binary) {
        # Handling binary data using MDS
        if (ncol(train) * nrow(train) > 2e+07) {
          pca <- list()
          pca$x <- cmdscale(d = stats::dist(train, method = binary_distance), 
                            k = min(nrow(train), 20))
          colnames(pca$x) <- paste0("PC", 1:ncol(pca$x))
        } else {
          pca <- list()
          pca$x <- cmdscale(d = stats::dist(train, method = binary_distance), 
                            k = min(nrow(train), 200))
          colnames(pca$x) <- paste0("PC", 1:ncol(pca$x))
        }
        
        # Project test data into the MDS space
        test_distances <- as.matrix(stats::dist(rbind(train, test), method = binary_distance))
        test_train_distances <- test_distances[(nrow(train) + 1):nrow(test_distances), 1:nrow(train)]
        
        # Use Gower interpolation to project test data
        test_coords <- mds_project(train_coords = pca$x, 
                                   train_distances = as.matrix(stats::dist(train, method = binary_distance)),
                                   test_distances = test_train_distances)
        
        test <- test_coords  # Use projected test coordinates
      } else {
        # Handling continuous data using PCA
        if (ncol(train) * nrow(train) > 2e+07) {
          pca <- rpca.para(train, min(nrow(train), 20), scale = F, p = 10)
        } else {
          pca <- prcomp(train, rank. = min(nrow(train), 200))
        }
        
        test <- predict.rpca.para(pca, test)  # Project test data for continuous case
      }
    } else {
      # If a custom PCA function is provided
      pca <- list(x = PCAFunction(train))
    }
    
    # Set train to transformed PCA or MDS coordinates
    train <- pca$x
    if (!is.null(k) || kMin == kMax) {
      if (is.null(k)) {
        k = kMin
      }
      set.seed(seed)
      cluster.train <- kmeans(pca$x, centers = k, iter.max = 1000, 
                              nstart = 1000)$cluster
      clus <- k
      origS <- NULL
      pertS <- NULL
      Discrepancy <- NULL
    }
    else {
      log("Building original connectivity matrices")
      set.seed(seed)
      origPartition <- GetOriginalSimilarity(data = train, 
                                             clusRange = kMin:kMax, clusteringAlgorithm = clusteringAlgorithm$fun, 
                                             showProgress = verbose, ncore)
      origS <- origPartition$origS
      listAUC <- list()
      for (k in kMin:kMax) {
        listAUC[k] <- list(c())
      }
      stoppingCriteriaHandler <- DiffDevianceStoppingHandler(kMin = kMin, 
                                                             kMax = kMax, origS = origS, iterMin = iterMin, 
                                                             madMin = madMin, msdMin = msdMin, onExcute = function(k, 
                                                                                                                   AUCs) {
                                                               listAUC[[k]] <<- AUCs
                                                             })
      log("Building perturbed connectivity matrices")
      set.seed(seed)
      pertS <- GetPerturbedSimilarity(data = train, clusRange = kMin:kMax, 
                                      iterMax = iterMax, iterMin = iterMin, origS = origS, 
                                      clusteringAlgorithm = clusteringAlgorithm$fun, 
                                      perturbedFunction = perturbationAlgorithm$fun, 
                                      stoppingCriteriaHandler = stoppingCriteriaHandler, 
                                      showProgress = verbose, ncore = ncore)
      for (k in kMin:kMax) {
        AUCs <- listAUC[[k]]
        pert <- matrix(0, nrow(train), nrow(train))
        for (i in 1:length(AUCs)) {
          pert <- pert + pertS[[k]][[i]] * length(which(AUCs == 
                                                          AUCs[i]))
        }
        pertS[[k]] <- pert/sum(as.matrix(table(AUCs))[, 
                                                      1]^2)
      }
      Discrepancy <- CalcPerturbedDiscrepancy(origS, pertS, 
                                              clusRange = kMin:kMax)
      Discrepancy$AUC = round(Discrepancy$AUC, digits = 4)
      clus <- min(which(Discrepancy$AUC == max(Discrepancy$AUC[kMin:kMax])))
      cluster.train <- origPartition$groupings[[clus]]
    }
    cluster.test <- as.integer(classify(train, cluster.train, 
                                        test, knn.k))
    cluster <- rep(0, nrow(data))
    cluster[ind] <- cluster.train
    cluster[-ind] <- cluster.test
    timediff = Sys.time() - now
    log("Done in ", timediff, " ", units(timediff), ".\n")
    list(k = clus, cluster = cluster, origS = origS, pertS = pertS, 
         Discrepancy = Discrepancy, pca = pca)
  }
}

logisticPCA_optimized <- function(x, k = 2, m = 4, quiet = TRUE, partial_decomp = FALSE,
                                  max_iters = 1000, conv_criteria = 1e-05, random_start = FALSE,
                                  start_U = NULL, start_mu = NULL, main_effects = TRUE, validation = NULL,
                                  M = NULL, use_irlba = NULL) {
  
  library(Matrix)  # Sparse matrices support
  library(RcppParallel)  # Parallelization for matrix operations
  
  # Utility functions (Assumes inv.logit.mat and log_like_Bernoulli are pre-defined)
  inv.logit.mat <- function(x) {
    return(1 / (1 + exp(-x)))
  }
  
  log_like_Bernoulli <- function(q, theta) {
    return(sum(q * log(inv.logit.mat(theta)) + (1 - q) * log(1 - inv.logit.mat(theta))))
  }
  
  # Compatibility for deprecated arguments
  if (!missing(M)) {
    m = M
    warning("M is deprecated. Use m instead. Using m = ", m)
  }
  if (!missing(use_irlba)) {
    partial_decomp = use_irlba
    warning("use_irlba is deprecated. Use partial_decomp instead.")
  }
  
  # Convert data to sparse matrix
  q <- as.matrix(2 * x - 1)
  q <- as(q, "dgCMatrix")
  
  missing_mat <- is.na(q)
  q[is.na(q)] <- 0  # Handling missing values by replacing NAs with 0
  n <- nrow(q)
  d <- ncol(q)
  
  if (k >= d & partial_decomp) {
    message("k >= dimension. Setting partial_decomp = FALSE")
    partial_decomp <- FALSE
    k <- d
  }
  
  if (main_effects) {
    mu <- if (!is.null(start_mu)) start_mu else colMeans(m * q)
  } else {
    mu <- rep(0, d)
  }
  
  if (!is.null(start_U)) {
    U <- sweep(start_U, 2, sqrt(colSums(start_U^2)), "/")
  } else if (random_start) {
    U <- matrix(rnorm(d * k), d, k)
    U <- qr.Q(qr(U))
  } else {
    # SVD initialization (Use rARPACK for partial decomposition when applicable)
    if (partial_decomp && requireNamespace("rARPACK", quietly = TRUE)) {
      udv <- rARPACK::svds(scale(q, center = main_effects, scale = FALSE), k = k)
    } else {
      udv <- svd(scale(q, center = main_effects, scale = FALSE))
    }
    U <- matrix(udv$v[, 1:k], d, k)
  }
  
  qTq <- crossprod(q)  # Sparse crossproduct
  loss_trace <- numeric(max_iters + 1)
  eta <- m * q + missing_mat * outer(rep(1, n), mu)
  theta <- outer(rep(1, n), mu) + scale(eta, center = mu, scale = FALSE) %*% tcrossprod(U)
  loglike <- log_like_Bernoulli(q = q, theta = theta)
  loss_trace[1] <- (-loglike) / sum(q != 0)
  
  ptm <- proc.time()
  
  for (i in 1:max_iters) {
    last_U <- U
    last_m <- m
    last_mu <- mu
    
    Z <- as.matrix(theta + 4 * q * (1 - inv.logit.mat(q * theta)))
    
    if (main_effects) {
      mu <- as.numeric(colMeans(Z - eta %*% tcrossprod(U)))
    }
    
    eta <- m * q + missing_mat * outer(rep(1, n), mu)
    
    mat_temp <- crossprod(scale(eta, center = mu, scale = FALSE), Z)
    mat_temp <- mat_temp + t(mat_temp) - crossprod(eta) + n * outer(mu, mu)
    
    # Use partial decomposition if applicable
    repeat {
      if (partial_decomp) {
        eig <- rARPACK::eigs_sym(mat_temp, k = min(k + 2, d))
      }
      if (!partial_decomp || any(eig$values[1:k] < 0)) {
        eig <- eigen(mat_temp, symmetric = TRUE)
      }
      U <- matrix(eig$vectors[, 1:k], d, k)
      theta <- outer(rep(1, n), mu) + scale(eta, center = mu, scale = FALSE) %*% tcrossprod(U)
      this_loglike <- log_like_Bernoulli(q = q, theta = theta)
      
      if (!partial_decomp | this_loglike >= loglike) {
        loglike <- this_loglike
        break
      } else {
        partial_decomp <- FALSE
        warning("rARPACK::eigs_sym was too inaccurate in iteration ", i, ". Switched to base::eigen")
      }
    }
    
    loss_trace[i + 1] <- (-loglike) / sum(q != 0)
    
    if (i > 4 && abs(loss_trace[i] - loss_trace[i + 1]) < conv_criteria) {
      break
    }
  }
  
  # Construct output object
  object <- list(mu = mu, U = U, PCs = scale(eta, center = mu, scale = FALSE) %*% U, 
                 m = m, M = m, iters = i, loss_trace = loss_trace[1:(i + 1)], 
                 prop_deviance_expl = 1 - loglike / sum(q != 0))
  
  class(object) <- "lpca"
  return(object)
}
