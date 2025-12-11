#' Modified Spectrum Clustering with Binary Data and Parallelization Support
#'
#' @description Modified version of Spectrum clustering algorithm with support for
#'   binary/categorical data via alternative distance metrics and parallel processing.
#'   Implements spectral clustering with adaptive kernel parameter selection.
#'
#' @param data A list of data matrices (one per view/modality) with samples as columns
#'   and features as rows, or a single matrix.
#' @param method Integer; method for selecting optimal K. 1 = eigenvalue gap,
#'   2 = eigenvector multimodality, 3 = fixed K. Default is 1.
#' @param silent Logical; if TRUE, suppress progress messages. Default is FALSE.
#' @param showres Logical; if TRUE, display diagnostic plots. Default is TRUE.
#' @param diffusion Logical; if TRUE, apply diffusion process on similarity matrix.
#'   Default is TRUE.
#' @param kerneltype Character; kernel type, either "density" (CNN kernel) or "stsc"
#'   (self-tuning spectral clustering). Default is "density".
#' @param maxk Integer; maximum number of clusters to consider. Default is 10.
#' @param NN Integer; number of nearest neighbors for kernel construction. Default is 3.
#' @param NN2 Integer; extended neighborhood size for CNN kernel. Default is 7.
#' @param showpca Logical; if TRUE, show PCA plot of results. Default is FALSE.
#' @param frac Numeric; fraction parameter for method 2. Default is 2.
#' @param thresh Numeric; threshold parameter for method 2. Default is 7.
#' @param fontsize Numeric; font size for plots. Default is 18.
#' @param dotsize Numeric; point size for plots. Default is 3.
#' @param tunekernel Logical; if TRUE, automatically tune kernel parameters. Default is FALSE.
#' @param clusteralg Character; clustering algorithm, "GMM" or "km". Default is "GMM".
#' @param FASP Logical; if TRUE, use Fast Approximate Spectral Clustering. Default is FALSE.
#' @param FASPk Integer; number of centroids for FASP. Required if FASP is TRUE.
#' @param fixk Integer; fixed K value for method 3. Required if method is 3.
#' @param krangemax Integer; maximum K for range clustering. Default is 10.
#' @param runrange Logical; if TRUE, cluster over range of K values. Default is FALSE.
#' @param diffusion_iters Integer; number of diffusion iterations. Default is 4.
#' @param KNNs_p Integer; KNN parameter for diffusion. Default is 10.
#' @param missing Logical; if TRUE, impute missing data. Default is FALSE.
#' @param distances Character; distance metric(s) for each view. Can be a single value
#'   or vector. Supports: "euclidean", "manhattan", "cosine", "binary", etc.
#' @param cores Integer; number of CPU cores for parallel processing. Default is 1.
#'
#' @return List containing cluster assignments, eigenvector analysis, optimal K,
#'   similarity matrix, and eigensystem decomposition.
#'
#' @export
Spectrum_bin_and_par <- function (
    data,
    method = 1,
    silent = FALSE,
    showres = TRUE,
    diffusion = TRUE,
    kerneltype = c("density", "stsc"),
    maxk = 10,
    NN = 3,
    NN2 = 7,
    showpca = FALSE,
    frac = 2,
    thresh = 7,
    fontsize = 18,
    dotsize = 3,
    tunekernel = FALSE,
    clusteralg = "GMM",
    FASP = FALSE,
    FASPk = NULL,
    fixk = NULL,
    krangemax = 10,
    runrange = FALSE,
    diffusion_iters = 4,
    KNNs_p = 10,
    missing = FALSE,
    distances = "euclidean",  # Added distances argument
    cores = 1                 # Added cores argument
)
{
  # Load necessary libraries
  required_packages <- c("foreach", "doParallel", "diptest", "Rfast", "ggplot2")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("The", pkg, "package is required but not installed. Please install it before proceeding."))
    }
  }
  
  library(foreach)
  library(doParallel)
  
  kerneltype <- match.arg(kerneltype)
  
  # Handle the 'distances' argument
  if (length(distances) == 1) {
    distances <- rep(distances, length(data))
  } else if (length(distances) != length(data)) {
    stop("Error: 'distances' must be either a single string or a vector with length equal to the number of data views.")
  }
  
  # Validate that all distance types are supported
  supported_distances <- c(
    "euclidean",
    "manhattan",
    "minimum",
    "maximum",
    "minkowski",
    "bhattacharyya",
    "hellinger",
    "kullback_leibler",
    "jensen_shannon",
    "haversine",
    "canberra",
    "chi_square",
    "soergel",
    "sorensen",
    "cosine",
    "wave_hedges",
    "motyka",
    "harmonic_mean",
    "jeffries_matusita",
    "gower",
    "kulczynski",
    "itakura_saito"
  )
  if (!all(distances %in% supported_distances)) {
    stop(paste("Error: Unsupported distance type detected. Supported distances are:",
               paste(supported_distances, collapse = ", ")))
  }
  
  # Validate 'cores' argument
  if (!is.numeric(cores) || length(cores) != 1 || cores < 1 || cores != floor(cores)) {
    stop("Error: 'cores' must be a positive integer.")
  }
  
  available_cores <- parallel::detectCores()
  if (cores > available_cores) {
    warning(paste("Requested number of cores (", cores, ") exceeds available cores (", available_cores, "). Using ", available_cores, " cores instead.", sep = ""))
    cores <- available_cores
  }
  
  # Convert data to list if not already
  if (!inherits(data, "list")) {
    datalist <- list(data)
    distances <- distances[1]  # Ensure distances aligns with datalist
  }
  else {
    datalist <- data
  }
  
  # Check FASP constraints
  if (length(datalist) > 1 & FASP == TRUE) {
    stop("Error: FASP method works for only a single view")
  }
  if (is.null(FASPk) == TRUE & FASP == TRUE) {
    stop("Error: FASP method requires a number of centroids to compute")
  }
  if (runrange == TRUE & method == 3) {
    stop("Error: cannot run a range of K whilst method=3")
  }
  if (is.null(fixk) == TRUE & method == 3) {
    stop("Error: need to set the value of K using the fixk parameter for method 3")
  }
  
  # Informational messages
  if (silent == FALSE) {
    message("***Spectrum***")
    message(paste("Detected views:", length(datalist)))
    message(paste("Method:", method))
    message(paste("Kernel type:", kerneltype))
    message(paste("Using", cores, "core(s) for parallel processing"))
  }
  
  # FASP data compression
  if (FASP) {
    message("Running with FASP data compression")
    cs <- kmeans(t(datalist[[1]]), centers = FASPk)
    csx <- cs$centers
    cas <- cs$cluster
    datalist[[1]] <- data.frame(t(csx))
  }
  
  # Set up parallel backend for kernel computation
  cl_kernel <- makeCluster(cores)
  registerDoParallel(cl_kernel)
  
  # List of helper functions to export
  helper_functions <- c("CNN_kernel_mod", "kernfinder_mine_mod", "kernfinder_local_mod")
  #,
  # "rbfkernel_b_mod", "EM_finder", "findk",
  # "plot_egap", "plot_multigap", "pca",
  # "harmonise_ids", "mean_imputation")
  
  # Parallelized kernel computation using foreach
  kernellist <- foreach(platform = seq_along(datalist),
                        .packages = c("Spectrum", "Rfast", "ggplot2", "diptest"),
                        .export = helper_functions) %dopar% {
                          if (silent == FALSE) {
                            message(paste("Calculating similarity matrix for modality", platform))
                          }
                          
                          # Retrieve the distance for the current modality
                          current_distance <- distances[platform]
                          
                          # Initialize kerneli
                          kerneli <- NULL
                          
                          # Modify kernel computation based on the specified distance
                          if (kerneltype == "stsc") {
                            if (method == 2 && tunekernel) {
                              NN_current <- kernfinder_local_mod(
                                datalist[[platform]],
                                maxk = maxk,
                                silent = silent,
                                fontsize = fontsize,
                                dotsize = dotsize,
                                showres = showres,
                                distance = current_distance  # Passing the specified distance
                              )
                              NN <- NN_current
                            }
                            
                            # Compute the RBF kernel based on the specified distance
                            kerneli <- rbfkernel_b_mod(
                              datalist[[platform]],
                              K = NN,
                              sigma = 1,
                              distance = current_distance  # Passing the specified distance
                            )
                          }
                          else if (kerneltype == "density") {
                            if (method == 2 && tunekernel) {
                              NN_current <- kernfinder_mine_mod(
                                datalist[[platform]],
                                maxk = maxk,
                                silent = silent,
                                showres = showres,
                                fontsize = fontsize,
                                dotsize = dotsize,
                                distance = current_distance  # Passing the specified distance
                              )
                              NN <- NN_current
                            }
                            
                            # Compute the CNN kernel based on the specified distance
                            kerneli <- CNN_kernel_mod(
                              datalist[[platform]],
                              NN = NN,
                              NN2 = NN2,
                              distance = current_distance  # Passing the specified distance
                            )
                          }
                          
                          if (silent == FALSE) {
                            message("Done.")
                          }
                          
                          return(kerneli)
                        }
  
  # Stop the kernel computation cluster
  stopCluster(cl_kernel)
  registerDoSEQ()
  
  # Handle missing data if required
  if (missing) {
    message("Imputing missing data...")
    kernellist <- harmonise_ids(kernellist)
    kernellist <- mean_imputation(kernellist)
    message("Done.")
  }
  
  # Combine similarity matrices
  if (silent == FALSE) {
    message("Combining similarity matrices and creating kNN graph...")
  }
  A <- Reduce("+", kernellist)
  
  # Diffusion process
  if (diffusion == TRUE) {
    for (col in seq_len(ncol(A))) {
      KNNs <- head(rev(sort(A[, col])), (KNNs_p + 1))
      tokeep <- names(KNNs)
      A[!(rownames(A) %in% tokeep), col] <- 0
    }
    A <- A / rowSums(A)
    if (silent == FALSE) {
      message("Done.")
    }
    if (silent == FALSE) {
      message("Diffusing on tensor product graph...")
    }
    Qt <- A
    im <- matrix(0, ncol = ncol(A), nrow = ncol(A))
    diag(im) <- 1
    for (t in seq_len(diffusion_iters)) {
      Qt <- A %*% Qt %*% t(A) + im
    }
    A2 <- t(Qt)
    if (silent == FALSE) {
      message("Done.")
    }
  }
  else if (diffusion == FALSE) {
    A2 <- A / length(datalist)
  }
  
  # Calculating graph Laplacian (L)
  if (silent == FALSE) {
    message("Calculating graph Laplacian (L)...")
  }
  dv <- 1 / sqrt(rowSums(A2))
  l <- dv * A2 %*% diag(dv)
  
  # Eigen decomposition and selecting optimal K
  if (method == 1) {
    if (silent == FALSE) {
      message("Performing eigen decomposition of L...")
    }
    decomp <- eigen(l)
    if (silent == FALSE) {
      message("Done.")
      message("Examining eigenvalues to select K...")
    }
    evals <- as.numeric(decomp$values)
    diffs <- diff(evals)
    diffs <- diffs[-1]
    if (maxk - 1 < length(diffs)) {
      diffs_subset <- diffs[1:(maxk - 1)]
    } else {
      diffs_subset <- diffs
    }
    optk <- which.max(abs(diffs_subset)) + 1
    if (silent == FALSE) {
      message(paste("Optimal K:", optk))
    }
    nn <- maxk + 1
    d <- data.frame(K = seq_len(nn), evals = evals[1:nn])
    if (showres == TRUE) {
      plot_egap(d, maxk = maxk, dotsize = dotsize, fontsize = fontsize)
    }
  }
  else if (method == 2) {
    if (silent == FALSE) {
      message("Performing eigen decomposition of L...")
    }
    decomp <- eigen(l)
    if (silent == FALSE) {
      message("Done.")
      message("Examining eigenvector distributions to select K...")
    }
    xi <- decomp$vectors[, 1:(maxk + 1)]
    res <- Spectrum:::EM_finder(xi, silent = silent)
    d <- data.frame(K = seq_len(maxk + 1), Z = res[1:(maxk + 1), 2])
    if (showres == TRUE) {
      Spectrum:::plot_multigap(d, maxk = maxk, dotsize = dotsize, fontsize = fontsize)
    }
    optk <- Spectrum:::findk(res, maxk = maxk, frac = frac, thresh = thresh)
    if (silent == FALSE) {
      message(paste("Optimal K:", optk))
    }
  }
  else if (method == 3) {
    decomp <- eigen(l)
    optk <- fixk
  }
  
  results <- list()
  
  # Parallelization setup for clustering over a range of K if runrange is TRUE
  if (runrange) {
    if (silent == FALSE) {
      message("Clustering over a range of K values...")
    }
    
    # Set up parallel backend for clustering
    cl_clustering <- makeCluster(cores)
    registerDoParallel(cl_clustering)
    
    # Parallelized clustering over range of K
    results <- foreach(tk = 2:krangemax,
                       .packages = c("Spectrum"),  # Replace "ClusterR" with "Spectrum"
                       .export = c(helper_functions, "cas", "d", "findk")) %dopar% {
                         xi <- decomp$vectors[, 1:tk]
                         yi <- xi / sqrt(rowSums(xi^2))
                         yi[!is.finite(yi)] <- 0
                         
                         if (clusteralg == "GMM") {
                           # Assuming Spectrum has GMM and predict_GMM functions
                           gmm <- ClusterR:::GMM(
                             yi,
                             tk,
                             verbose = FALSE,
                             seed_mode = "random_spread"
                           )
                           pr <- ClusterR::predict_GMM(
                             yi,
                             gmm$centroids,
                             gmm$covariance_matrices,
                             gmm$weights
                           )
                           names(pr)[3] <- "cluster"
                           if (0 %in% pr$cluster) {
                             pr$cluster <- pr$cluster + 1
                           }
                           if (silent == FALSE) {
                             message(paste("Clustered for K =", tk))
                           }
                         }
                         else if (clusteralg == "km") {
                           pr <- kmeans(yi, tk)
                           if (silent == FALSE) {
                             message(paste("Clustered for K =", tk))
                           }
                         }
                         
                         # Assemble the results as in the original function
                         if (method != 3) {
                           if (FASP) {
                             casn <- cas
                             casn <- casn[seq_along(casn)] <- pr$cluster[as.numeric(casn[seq_along(casn)])]
                             names(casn) <- names(cas)
                             list(
                               allsample_assignments = casn,
                               centroid_assignments = pr$cluster,
                               eigenvector_analysis = d,
                               K = tk,
                               similarity_matrix = A2,
                               eigensystem = decomp
                             )
                           }
                           else {
                             list(
                               assignments = pr$cluster,
                               eigenvector_analysis = d,
                               K = tk,
                               similarity_matrix = A2,
                               eigensystem = decomp
                             )
                           }
                         }
                       }
    
    # Stop the clustering cluster
    stopCluster(cl_clustering)
    registerDoSEQ()
  }
  else {
    xi <- decomp$vectors[, 1:optk]
    yi <- xi / sqrt(rowSums(xi^2))
    yi[!is.finite(yi)] <- 0
    if (clusteralg == "GMM") {
      if (silent == FALSE) {
        message("Performing GMM clustering...")
      }
      
      gmm <-ClusterR::GMM(
        yi,
        optk,
        verbose = FALSE,
        seed_mode = "random_spread"
      )
      pr <- ClusterR::predict_GMM(
        yi,
        gmm$centroids,
        gmm$covariance_matrices,
        gmm$weights
      )
      names(pr)[3] <- "cluster"
      if (0 %in% pr$cluster) {
        pr$cluster <- pr$cluster + 1
      }
      if (silent == FALSE) {
        message("Done.")
      }
    }
    else if (clusteralg == "km") {
      if (silent == FALSE) {
        message("Performing k-means clustering...")
      }
      pr <- kmeans(yi, optk)
      if (silent == FALSE) {
        message("Done.")
      }
    }
    if (length(datalist) == 1 && showres == TRUE) {
      if (showpca == TRUE) {
        pca(
          datalist[[1]],
          labels = as.factor(pr$cluster),
          axistextsize = fontsize,
          legendtextsize = fontsize,
          dotsize = dotsize
        )
      }
    }
    if (method != 3) {
      if (FASP) {
        casn <- cas
        casn <- casn[seq_along(casn)] <- pr$cluster[as.numeric(casn[seq_along(casn)])]
        names(casn) <- names(cas)
        results <- list(
          allsample_assignments = casn,
          centroid_assignments = pr$cluster,
          eigenvector_analysis = d,
          K = optk,
          similarity_matrix = A2,
          eigensystem = decomp
        )
      }
      else {
        results <- list(
          assignments = pr$cluster,
          eigenvector_analysis = d,
          K = optk,
          similarity_matrix = A2,
          eigensystem = decomp
        )
      }
    }
    else if (method == 3) {
      if (FASP) {
        casn <- cas
        casn <- casn[seq_along(casn)] <- pr$cluster[as.numeric(casn[seq_along(casn)])]
        names(casn) <- names(cas)
        results <- list(
          allsample_assignments = casn,
          centroid_assignments = pr$cluster,
          K = optk,
          similarity_matrix = A2,
          eigensystem = decomp
        )
      }
      else {
        results <- list(
          assignments = pr$cluster,
          K = optk,
          similarity_matrix = A2,
          eigensystem = decomp
        )
      }
    }
  }
  
  if (silent == FALSE) {
    message("Finished.")
  }
  return(results)
}


#' Common Nearest Neighbors Kernel with Modified Distance Support
#'
#' @description Computes a similarity matrix using the Common Nearest Neighbors (CNN)
#'   kernel. Extends the original Spectrum CNN kernel to support multiple distance
#'   metrics including those suitable for binary/categorical data.
#'
#' @param mat Data matrix with samples as columns and features as rows.
#' @param NN Integer; number of nearest neighbors for local sigma estimation. Default is 3.
#' @param NN2 Integer; extended neighborhood size for CNN similarity computation. Default is 7.
#' @param distance Character; distance metric to use. Supported values include
#'   "euclidean", "manhattan", "cosine", "binary", "canberra", etc.
#'
#' @return A symmetric similarity matrix with dimensions n x n where n is the number
#'   of samples. Values range from 0 to 1.
#'
#' @export
CNN_kernel_mod <- function(mat, NN = 3, NN2 = 7, distance = "euclidean")
{
  # Validate distance parameter
  supported_distances <- c(
    "euclidean",
    "manhattan",
    "minimum",
    "maximum",
    "minkowski",
    "bhattacharyya",
    "hellinger",
    "kullback_leibler",
    "jensen_shannon",
    "haversine",
    "canberra",
    "chi_square",
    "soergel",
    "sorensen",
    "cosine",
    "wave_hedges",
    "motyka",
    "harmonic_mean",
    "jeffries_matusita",
    "gower",
    "kulczynski",
    "itakura_saito"
  )
  if (!(distance %in% supported_distances)) {
    stop(paste("Unsupported distance type:", distance, ". Supported distances are:",
               paste(supported_distances, collapse = ", ")))
  }
  
  n <- ncol(mat)
  nbs <- list()
  
  # Compute distance matrix with specified distance metric
  dm <- Rfast::Dist(t(mat), method = distance)
  dimnames(dm) <- list(colnames(mat), colnames(mat))
  
  kn <- c()
  for (i in seq_len(n)) {
    sortedvec <- sort.int(dm[i, ], index.return = FALSE)
    kn <- c(kn, sortedvec[NN + 1])
    nbs[[i]] <- names(sortedvec[2:(NN2 + 1)])
    names(nbs)[i] <- names(sortedvec)[1]
  }
  
  sigmamatrix <- kn %o% kn
  out <- matrix(nrow = n, ncol = n)
  upper <- -dm^2
  
  for (i in 2:n) {
    for (j in 1:(i - 1)) {
      cnns <- length(intersect(nbs[[i]], nbs[[j]]))
      upperval <- upper[i, j]
      localsigma <- sigmamatrix[i, j]
      out[i, j] <- exp(upperval / (localsigma * (cnns + 1)))
    }
  }
  
  out <- pmax(out, t(out), na.rm = TRUE)
  diag(out) <- 1
  colnames(out) <- colnames(mat)
  row.names(out) <- colnames(mat)
  
  return(out)
}

#' Kernel Parameter Finder Using Eigenvector Multimodality (Modified)
#'
#' @description Finds the optimal nearest neighbor (NN) parameter for the CNN kernel
#'   by examining eigenvector distributions using dip tests for multimodality.
#'   Modified to support alternative distance metrics.
#'
#' @param data Data matrix with samples as columns and features as rows.
#' @param maxk Integer; maximum number of clusters to consider. Default is 10.
#' @param fontsize Numeric; font size for diagnostic plots. Default is 12.
#' @param silent Logical; if TRUE, suppress progress messages. Default is FALSE.
#' @param showres Logical; if TRUE, display diagnostic plot. Default is TRUE.
#' @param dotsize Numeric; point size for plots. Default is 2.
#' @param distance Character; distance metric to use. Default is "euclidean".
#'
#' @return Integer; optimal NN parameter value.
#'
#' @export
kernfinder_mine_mod <- function(data, maxk = 10, fontsize = 12, silent = FALSE,
                                showres = TRUE, dotsize = 2, distance = "euclidean")
{
  # Validate distance parameter
  supported_distances <- c(
    "euclidean",
    "manhattan",
    "minimum",
    "maximum",
    "minkowski",
    "bhattacharyya",
    "hellinger",
    "kullback_leibler",
    "jensen_shannon",
    "haversine",
    "canberra",
    "chi_square",
    "soergel",
    "sorensen",
    "cosine",
    "wave_hedges",
    "motyka",
    "harmonic_mean",
    "jeffries_matusita",
    "gower",
    "kulczynski",
    "itakura_saito"
  )
  if (!(distance %in% supported_distances)) {
    stop(paste("Unsupported distance type:", distance, ". Supported distances are:",
               paste(supported_distances, collapse = ", ")))
  }
  
  if (silent == FALSE) {
    message("Finding optimal NN kernel parameter by examining eigenvector distributions")
  }
  
  rr <- c()
  for (param in seq(1, 10)) {
    if (silent == FALSE) {
      message(paste("Tuning kernel NN parameter:", param))
    }
    
    # Pass 'distance' to CNN_kernel
    kern <- CNN_kernel_mod(data, NN = param, NN2 = 7, distance = distance)
    kern[which(!is.finite(kern))] <- 0
    
    dv <- 1 / sqrt(rowSums(kern))
    l <- dv * kern %*% diag(dv)
    xi <- eigen(l)$vectors
    
    res <- matrix(nrow = ncol(xi), ncol = 2)
    for (ii in seq(1, ncol(xi))) {
      r <- diptest::dip.test(xi[, ii], simulate.p.value = FALSE, B = 2000)
      res[ii, 1] <- r$p.value
      res[ii, 2] <- r$statistic
    }
    
    diffs <- diff(res[, 2])
    diffs <- diffs[-1]
    tophit <- diffs[1:(maxk + 1)][which.min(diffs[1:(maxk + 1)])]
    rr <- c(rr, tophit)
  }
  
  optimalparam <- which.min(rr)
  if (silent == FALSE) {
    message(paste("Optimal NN:", optimalparam))
  }
  
  d <- data.frame(x = seq(1, 10), y = rr)
  py <- ggplot2::ggplot(data = d, aes(x = x, y = y)) +
    ggplot2::geom_point(colour = "black", size = dotsize) +
    ggplot2::theme_bw() +
    ggplot2::geom_line() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = fontsize, colour = "black"),
      axis.text.x = ggplot2::element_text(size = fontsize, colour = "black"),
      axis.title.x = ggplot2::element_text(size = fontsize),
      axis.title.y = ggplot2::element_text(size = fontsize),
      legend.text = ggplot2::element_text(size = fontsize),
      legend.title = ggplot2::element_text(size = fontsize),
      plot.title = ggplot2::element_text(size = fontsize, colour = "black", hjust = 0.5),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::ylab("D") +
    ggplot2::xlab("NN") +
    ggplot2::scale_x_continuous(limits = c(1, 10), breaks = seq(1, 10, by = 1))
  
  if (showres == TRUE) {
    print(py)
  }
  
  return(optimalparam)
}

#' Kernel Parameter Finder for Local Scaling (Modified)
#'
#' @description Finds the optimal nearest neighbor (NN) parameter for the RBF kernel
#'   with local scaling by examining eigenvector distributions using dip tests.
#'   Modified to support alternative distance metrics for STSC kernel type.
#'
#' @param data Data matrix with samples as columns and features as rows.
#' @param maxk Integer; maximum number of clusters to consider. Default is 10.
#' @param fontsize Numeric; font size for diagnostic plots. Default is 12.
#' @param silent Logical; if TRUE, suppress progress messages. Default is FALSE.
#' @param showres Logical; if TRUE, display diagnostic plot. Default is TRUE.
#' @param dotsize Numeric; point size for plots. Default is 2.
#' @param distance Character; distance metric to use. Default is "euclidean".
#'
#' @return Integer; optimal NN parameter value.
#'
#' @export
kernfinder_local_mod <- function(data, maxk = 10, fontsize = 12, silent = FALSE,
                                 showres = TRUE, dotsize = 2, distance = "euclidean")
{
  # Validate distance parameter
  supported_distances <- c(
    "euclidean",
    "manhattan",
    "minimum",
    "maximum",
    "minkowski",
    "bhattacharyya",
    "hellinger",
    "kullback_leibler",
    "jensen_shannon",
    "haversine",
    "canberra",
    "chi_square",
    "soergel",
    "sorensen",
    "cosine",
    "wave_hedges",
    "motyka",
    "harmonic_mean",
    "jeffries_matusita",
    "gower",
    "kulczynski",
    "itakura_saito"
  )
  if (!(distance %in% supported_distances)) {
    stop(paste("Unsupported distance type:", distance, ". Supported distances are:",
               paste(supported_distances, collapse = ", ")))
  }
  
  if (silent == FALSE) {
    message("Finding optimal kernel NN parameter by examining eigenvectors")
  }
  
  rr <- c()
  for (param in seq(1, 10)) {
    if (silent == FALSE) {
      message(paste("Tuning NN parameter:", param))
    }
    
    # Pass 'distance' to rbfkernel_b
    kern <- rbfkernel_b_mod(data, K = param, sigma = 1, distance = distance)
    
    # Handle non-finite values
    kern[!is.finite(kern)] <- 0
    
    # Compute graph Laplacian
    dv <- 1 / sqrt(rowSums(kern))
    l <- dv * kern %*% diag(dv)
    
    # Eigen decomposition
    xi <- eigen(l)$vectors
    
    # Perform Dip Test on each eigenvector
    res <- matrix(nrow = ncol(xi), ncol = 2)
    for (ii in seq_len(ncol(xi))) {
      r <- diptest::dip.test(xi[, ii], simulate.p.value = FALSE, B = 2000)
      res[ii, 1] <- r$p.value
      res[ii, 2] <- r$statistic
    }
    
    # Calculate differences in dip statistics
    diffs <- diff(res[, 2])
    diffs <- diffs[-1]
    
    # Identify the parameter with the minimum dip statistic difference
    tophit <- diffs[1:(maxk + 1)][which.min(diffs[1:(maxk + 1)])]
    rr <- c(rr, tophit)
  }
  
  # Determine the optimal NN parameter
  optimalparam <- which.min(rr)
  
  if (silent == FALSE) {
    message(paste("Optimal NN:", optimalparam))
  }
  
  # Plot the dip statistic differences
  d <- data.frame(x = seq(1, 10), y = rr)
  py <- ggplot2::ggplot(data = d, aes(x = x, y = y)) +
    ggplot2::geom_point(colour = "black", size = dotsize) +
    ggplot2::theme_bw() +
    ggplot2::geom_line() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = fontsize, colour = "black"),
      axis.text.x = ggplot2::element_text(size = fontsize, colour = "black"),
      axis.title.x = ggplot2::element_text(size = fontsize),
      axis.title.y = ggplot2::element_text(size = fontsize),
      legend.text = ggplot2::element_text(size = fontsize),
      legend.title = ggplot2::element_text(size = fontsize),
      plot.title = ggplot2::element_text(size = fontsize, colour = "black", hjust = 0.5),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::ylab("D") +
    ggplot2::xlab("NN") +
    ggplot2::scale_x_continuous(limits = c(1, 10), breaks = seq(1, 10, by = 1))
  
  if (showres == TRUE) {
    print(py)
  }
  
  return(optimalparam)
}

#' RBF Kernel with Local Scaling (Modified)
#'
#' @description Computes a Radial Basis Function (RBF) kernel with local scaling
#'   using K nearest neighbors. Modified to support alternative distance metrics
#'   for use with binary or categorical data.
#'
#' @param mat Data matrix with samples as columns and features as rows.
#' @param K Integer; number of nearest neighbors for local sigma estimation. Default is 3.
#' @param sigma Numeric; scaling parameter for the kernel. Default is 1.
#' @param distance Character; distance metric to use. Supported values include
#'   "euclidean", "manhattan", "cosine", "binary", "canberra", etc.
#'
#' @return A symmetric similarity matrix with dimensions n x n where n is the number
#'   of samples. Values range from 0 to 1.
#'
#' @export
rbfkernel_b_mod <- function(mat, K = 3, sigma = 1, distance = "euclidean")
{
  # Validate distance parameter
  supported_distances <- c(
    "euclidean",
    "manhattan",
    "minimum",
    "maximum",
    "minkowski",
    "bhattacharyya",
    "hellinger",
    "kullback_leibler",
    "jensen_shannon",
    "haversine",
    "canberra",
    "chi_square",
    "soergel",
    "sorensen",
    "cosine",
    "wave_hedges",
    "motyka",
    "harmonic_mean",
    "jeffries_matusita",
    "gower",
    "kulczynski",
    "itakura_saito"
  )
  
  if (!(distance %in% supported_distances)) {
    stop(paste("Unsupported distance type:", distance, ". Supported distances are:",
               paste(supported_distances, collapse = ", ")))
  }
  
  n <- ncol(mat)
  NN <- K
  nbs <- list()
  
  # Compute distance matrix with specified distance metric
  dm <- Rfast::Dist(t(mat), method = distance)
  dimnames(dm) <- list(colnames(mat), colnames(mat))
  
  kn <- c()
  for (i in seq_len(n)) {
    sortedvec <- as.numeric(sort.int(dm[i, ]))
    sortedvec <- sortedvec[sortedvec != 0] # Exclude zero distances
    if (length(sortedvec) < NN) {
      stop(paste("Not enough neighbors for sample", colnames(mat)[i],
                 "with NN =", NN))
    }
    kn <- c(kn, sortedvec[NN])
  }
  
  sigmamatrix <- kn %o% kn
  upper <- -dm^2
  out <- matrix(nrow = n, ncol = n)
  
  for (i in 2:n) {
    for (j in 1:(i - 1)) {
      lowerval <- sigmamatrix[i, j]
      upperval <- upper[i, j]
      out[i, j] <- exp(upperval / (lowerval * sigma))
    }
  }
  
  out <- pmax(out, t(out), na.rm = TRUE)
  diag(out) <- 1
  colnames(out) <- colnames(mat)
  row.names(out) <- colnames(mat)
  
  return(out)
}
