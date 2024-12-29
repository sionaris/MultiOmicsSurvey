# Selection of optimum number of clusters (k)
nmf.opt.k_parallel <- function(dat = dat, 
                      n.runs = 30, 
                      n.fold = 5, 
                      k.range = 2:8, 
                      result = TRUE, 
                      make.plot = TRUE, 
                      progress = TRUE,
                      st.count = 10, 
                      maxiter = 100, 
                      wt = if (is.list(dat)) rep(1, length(dat)) else 1,
                      num_cores = 1) {  # Added num_cores as an argument
  
  library(IntNMF)
  library(parallel)
  library(doParallel)
  library(foreach)
  library(MASS)
  library(mclust)
  library(doRNG)
  
  # ------------------------------
  # Input Validation
  # ------------------------------
  
  if (!is.list(dat) & length(wt) > 1) {
    stop("Weight is not applicable for single data")
  }
  
  if (!is.list(dat)) {
    dat <- list(dat)
  }
  
  n.dat <- length(dat)
  
  if (n.dat != length(wt)) {
    stop("Number of weights must match number of data")
  }
  
  # Check for non-negative data
  for (i in 1:n.dat) {
    if (any(dat[[i]] < 0)) {
      stop(paste("All values must be positive. There are -ve entries in dat", i))
    }
  }
  
  set.seed(12345)  # For reproducibility
  n.sample <- nrow(dat[[1]])
  
  # Initialize CPI matrix
  CPI <- matrix(NA, nrow = length(k.range), ncol = n.runs)
  dimnames(CPI) <- list(paste("k", k.range, sep = ""), paste("run", 1:n.runs, sep = ""))
  
  # ------------------------------
  # Parallel Backend Setup
  # ------------------------------
  
  if (!is.numeric(num_cores) || num_cores < 1) {
    stop("num_cores must be a positive integer")
  }
  
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  # Ensure proper cleanup in case of error
  on.exit(stopCluster(cl))
  
  # Export necessary variables and load required packages on workers
  clusterExport(cl, varlist = c("dat", "n.dat", "k.range", "n.fold", "st.count", "maxiter", "wt", "n.sample"), envir = environment())
  clusterEvalQ(cl, {
    library(IntNMF)
    library(MASS)
    library(mclust)
  })
  
  # ------------------------------
  # Parallel Computation
  # ------------------------------
  
  # To ensure reproducibility across parallel runs, set different seeds for each run
  # Use a reproducible parallel RNG
  # Install and use the doRNG package for reproducible foreach loops
  if (!requireNamespace("doRNG", quietly = TRUE)) {
    install.packages("doRNG", repos = "https://cran.r-project.org")
  }
  library(doRNG)
  
  CPI_results <- foreach(run = 1:n.runs, 
                         .combine = 'cbind', 
                         .packages = c("IntNMF", "MASS", "mclust"),
                         .options.RNG = 12345) %dorng% {  # Using doRNG for reproducibility
                           
                           run_CPI <- numeric(length(k.range))
                           names(run_CPI) <- paste("k", k.range, sep = "")
                           
                           for (ki in seq_along(k.range)) {
                             k <- k.range[ki]
                             R.ind <- numeric(n.fold)
                             
                             random.sample <- sample(n.sample)
                             
                             for (j in 1:n.fold) {
                               # Define test and train indices
                               test_start <- floor((j - 1) * n.sample / n.fold) + 1
                               test_end <- floor(j * n.sample / n.fold)
                               test.sample <- random.sample[test_start:test_end]
                               train.sample <- setdiff(random.sample, test.sample)
                               
                               # Split data
                               d.train <- lapply(dat, function(x) x[train.sample, , drop = FALSE])
                               d.test <- lapply(dat, function(x) x[test.sample, , drop = FALSE])
                               
                               # Fit NMF on training data
                               fit.train <- nmf.mnnals(dat = d.train, 
                                                       k = k, 
                                                       maxiter = maxiter, 
                                                       st.count = st.count, 
                                                       n.ini = 1,
                                                       ini.nndsvd = FALSE, 
                                                       seed = FALSE, 
                                                       wt = wt)
                               
                               # Compute XHt and HHt
                               XHt <- Reduce(`+`, lapply(1:n.dat, function(m) {
                                 sqrt(wt[m]) * d.test[[m]] %*% t(fit.train$H[[m]])
                               }))
                               
                               HHt <- Reduce(`+`, lapply(1:n.dat, function(m) {
                                 sqrt(wt[m]) * fit.train$H[[m]] %*% t(fit.train$H[[m]])
                               }))
                               
                               # Predict W
                               W.predict <- XHt %*% MASS::ginv(HHt)
                               W.predict <- W.predict + abs(min(W.predict))
                               predicted.cluster.mem <- apply(W.predict, 1, which.max)
                               
                               # Fit NMF on test data to compute cluster membership
                               fit.test <- nmf.mnnals(dat = d.test, 
                                                      k = k, 
                                                      maxiter = maxiter, 
                                                      st.count = st.count, 
                                                      n.ini = 1,
                                                      ini.nndsvd = FALSE, 
                                                      seed = FALSE, 
                                                      wt = wt)
                               computed.cluster.mem <- fit.test$clusters
                               
                               # Compute Adjusted Rand Index
                               R.ind[j] <- mclust::adjustedRandIndex(predicted.cluster.mem, computed.cluster.mem)
                             }
                             
                             # Store mean R.ind for this k
                             run_CPI[ki] <- mean(R.ind)
                           }
                           
                           run_CPI
                         }
  
  # ------------------------------
  # Assign results to CPI
  # ------------------------------
  
  CPI <- CPI_results
  
  # ------------------------------
  # Plotting
  # ------------------------------
  
  if (make.plot) {
    # Open a new plotting device
    # For non-interactive environments, use a specific device like png or pdf
    if (interactive()) {
      if (.Platform$OS.type == "windows") {
        windows(width = 4, height = 5)
      } else if (Sys.info()["sysname"] == "Darwin") {
        quartz(width = 4, height = 5)
      } else {
        x11(width = 4, height = 5)
      }
    } else {
      # If not interactive, save plot to a file
      png(filename = "nmf_opt_k_plot.png", width = 800, height = 600)
    }
    
    plot(k.range, CPI[,1], 
         ylim = c(min(CPI, na.rm = TRUE), max(CPI, na.rm = TRUE)), 
         pch = 20, 
         main = "", 
         xlab = "k", 
         ylab = "CPI")
    
    for (m in 2:n.runs) {
      points(k.range, CPI[,m], pch = 20)
    }
    
    lines(k.range, apply(CPI, 1, mean), col = "red", lwd = 2)
    mtext("Optimum k", outer = TRUE, cex = 1, line = -2)
    
    if (!interactive()) {
      dev.off()
    }
  }
  
  # ------------------------------
  # Return Results
  # ------------------------------
  
  if (result) return(CPI)
}
