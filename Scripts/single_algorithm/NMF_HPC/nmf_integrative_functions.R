nmf.opt.k.integrative <- function(dat, 
                                  is.binary,
                                  n.runs = 30,
                                  n.fold = 5,
                                  k.range = 2:8,
                                  result = TRUE,
                                  make.plot = TRUE,
                                  progress = TRUE,
                                  maxiter = 100,
                                  lr = 1e-3,
                                  tol = 1e-6,
                                  allowParallel = FALSE,
                                  n.cores = NULL,
                                  seed = 12345) {
  # Required libraries
  library(mclust)    # for adjustedRandIndex
  library(foreach)
  library(doParallel)
  library(MASS)      # for ginv, if needed
  
  # Basic checks
  if (!is.list(dat)) stop("Input 'dat' must be a list of matrices.")
  M <- length(dat)
  if (length(is.binary) != M) 
    stop("Length of 'is.binary' must match the number of modalities in 'dat'.")
  
  # Ensure non-negativity
  for (i in seq_len(M)) {
    if (min(dat[[i]]) < 0) {
      dat[[i]] <- pmax(dat[[i]] + abs(min(dat[[i]])), 0) + .Machine$double.eps
    }
  }
  
  # Weights for binary vs continuous
  M_b <- sum(is.binary)
  M_c <- M - M_b
  if (M_b > 0 && M_c > 0) {
    total_binary_weight <- M_b / M
    total_cont_weight   <- M_c / M
    wt <- numeric(M)
    wt[is.binary]  <- total_binary_weight / M_b
    wt[!is.binary] <- total_cont_weight   / M_c
  } else {
    wt <- rep(1/M, M)  # All binary or all continuous
  }
  
  sigmoid <- function(x) 1 / (1 + exp(-x))
  
  #--------------------------------------#
  # 1) Training Fit Function
  #    Now returns (W, H, clusters, iteration_count)
  #--------------------------------------#
  nmf.integrative.fit <- function(dat.list, k, maxiter, lr, tol, wt, is.binary, seed) {
    n <- nrow(dat.list[[1]])
    set.seed(seed)
    W <- matrix(runif(n * k, min = 0, max = 1), n, k)
    
    H.list <- vector("list", M)
    for (m in seq_len(M)) {
      H.list[[m]] <- matrix(runif(k * ncol(dat.list[[m]]), min = 0, max = 1), k, ncol(dat.list[[m]]))
    }
    
    prev_W <- W
    iteration_count <- 0
    for (iter in seq_len(maxiter)) {
      iteration_count <- iter
      
      # Update H
      for (m in seq_len(M)) {
        X_m <- dat.list[[m]]
        pred_m <- W %*% H.list[[m]]
        
        if (is.binary[m]) {
          # Logistic update
          Sigm <- sigmoid(pred_m)
          Grad_H <- wt[m] * t(W) %*% (Sigm - X_m)
          H.list[[m]] <- pmax(H.list[[m]] - lr * Grad_H, 0)
        } else {
          # Continuous update
          Grad_H <- wt[m] * t(W) %*% ((W %*% H.list[[m]]) - X_m)
          H.list[[m]] <- pmax(H.list[[m]] - lr * Grad_H, 0)
        }
      }
      
      # Update W
      W_grad <- matrix(0, nrow = n, ncol = k)
      for (m in seq_len(M)) {
        X_m <- dat.list[[m]]
        pred_m <- W %*% H.list[[m]]
        if (is.binary[m]) {
          Sigm <- sigmoid(pred_m)
          W_grad <- W_grad + wt[m] * ((Sigm - X_m) %*% t(H.list[[m]]))
        } else {
          W_grad <- W_grad + wt[m] * (((W %*% H.list[[m]]) - X_m) %*% t(H.list[[m]]))
        }
      }
      W <- pmax(W - lr * W_grad, 0)
      
      rel_change <- sum(abs(W - prev_W)) / (sum(abs(prev_W)) + .Machine$double.eps)
      if (rel_change < tol) break
      prev_W <- W
    }
    
    clusters <- apply(W, 1, which.max)
    list(W = W, H = H.list, clusters = clusters, iteration_count = iteration_count)
  }
  
  #--------------------------------------#
  # 2) Test W Estimation
  #    Now returns both W_test and iteration_count
  #--------------------------------------#
  nmf.test.W <- function(d.test, H.list, wt, is.binary, lr, tol, maxiter, seed) {
    n.test <- nrow(d.test[[1]])
    k <- nrow(H.list[[1]])
    set.seed(seed + 1)
    W_test <- matrix(runif(n.test * k, min=0, max=1), n.test, k)
    
    prev_W <- W_test
    iteration_count <- 0
    for (iter in seq_len(maxiter)) {
      iteration_count <- iter
      W_grad <- matrix(0, n.test, k)
      for (m in seq_len(M)) {
        X_m <- d.test[[m]]
        pred_m <- W_test %*% H.list[[m]]
        if (is.binary[m]) {
          Sigm <- sigmoid(pred_m)
          W_grad <- W_grad + wt[m] * ((Sigm - X_m) %*% t(H.list[[m]]))
        } else {
          # continuous
          W_grad <- W_grad + wt[m] * (((W_test %*% H.list[[m]]) - X_m) %*% t(H.list[[m]]))
        }
      }
      W_test <- pmax(W_test - lr * W_grad, 0)
      
      rel_change <- sum(abs(W_test - prev_W)) / (sum(abs(prev_W)) + .Machine$double.eps)
      if (rel_change < tol) break
      prev_W <- W_test
    }
    list(W_test = W_test, iteration_count = iteration_count)
  }
  
  #--------------------------------------#
  # Optionally Setup Parallelization
  #--------------------------------------#
  if (allowParallel) {
    if (is.null(n.cores)) n.cores <- parallel::detectCores() - 1
    cl <- parallel::makeCluster(n.cores)
    doParallel::registerDoParallel(cl)
  }
  
  #--------------------------------------#
  # Prepare to collect results
  # We'll do a custom .combine to store CPI, plus iteration counts
  #--------------------------------------#
  # Each run returns:
  #   $cpi            (length(k.range)) numeric
  #   $train_iter     (length(k.range)) numeric
  #   $test_iter      (length(k.range)) numeric
  #
  # We'll combine them into three matrices cpi, train_iter, test_iter, each size
  # (length(k.range), n.runs).
  
  combineRuns <- function(x, y) {
    list(
      cpi        = cbind(x$cpi,        y$cpi),
      train_iter = cbind(x$train_iter, y$train_iter),
      test_iter  = cbind(x$test_iter,  y$test_iter)
    )
  }
  
  initValue <- list(
    cpi        = matrix(, nrow=length(k.range), ncol=0),
    train_iter = matrix(, nrow=length(k.range), ncol=0),
    test_iter  = matrix(, nrow=length(k.range), ncol=0)
  )
  
  start_time <- Sys.time()
  
  res_list <- foreach(i = 1:n.runs, 
                      .combine = combineRuns, 
                      .init = initValue,
                      .packages = c("mclust")) %dopar% 
    {
      run_cpi        <- numeric(length(k.range))
      run_train_iter <- numeric(length(k.range))
      run_test_iter  <- numeric(length(k.range))
      
      for (ki in seq_along(k.range)) {
        k <- k.range[ki]
        R.ind <- c()
        sum_train_iter <- 0
        sum_test_iter  <- 0
        
        random.sample <- sample(seq_len(nrow(dat[[1]])), nrow(dat[[1]]))
        fold_size <- floor(nrow(dat[[1]]) / n.fold)
        
        for (j in 1:n.fold) {
          test_idx <- ((j-1)*fold_size+1):min(j*fold_size, nrow(dat[[1]]))
          test.sample  <- random.sample[test_idx]
          train.sample <- setdiff(random.sample, test.sample)
          
          d.train <- lapply(dat, function(x) x[train.sample, , drop = FALSE])
          d.test  <- lapply(dat, function(x) x[test.sample,  , drop = FALSE])
          
          # Train model
          fit.train <- nmf.integrative.fit(d.train, k=k, maxiter=maxiter,
                                           lr=lr, tol=tol, wt=wt,
                                           is.binary=is.binary, seed=seed)
          sum_train_iter <- sum_train_iter + fit.train$iteration_count
          
          # Compute test W
          testW_res <- nmf.test.W(d.test, fit.train$H, wt, is.binary, lr, tol, maxiter, seed)
          predicted.cluster.mem <- apply(testW_res$W_test, 1, which.max)
          sum_test_iter <- sum_test_iter + testW_res$iteration_count
          
          # Fit model on test set for ground truth clusters
          fit.test <- nmf.integrative.fit(d.test, k=k, maxiter=maxiter,
                                          lr=lr, tol=tol, wt=wt,
                                          is.binary=is.binary, seed=seed+10)
          computed.cluster.mem <- fit.test$clusters
          
          R.ind <- c(R.ind, mclust::adjustedRandIndex(predicted.cluster.mem, computed.cluster.mem))
          
          if (progress) {
            done <- ((i-1)*length(k.range)*n.fold + (ki-1)*n.fold + j) / (n.runs*length(k.range)*n.fold)
            pct  <- round(done * 100)
            if (pct %in% seq(5, 100, by=5)) {
              message(pct, "% complete")
              flush.console()
            }
          }
        }
        # Store mean ARI across folds in run_cpi for this k
        run_cpi[ki] <- mean(R.ind)
        
        # Average iteration counts across folds
        run_train_iter[ki] <- sum_train_iter / n.fold
        run_test_iter[ki]  <- sum_test_iter  / n.fold
      }
      
      list(
        cpi        = matrix(run_cpi,        nrow=length(k.range), ncol=1),
        train_iter = matrix(run_train_iter, nrow=length(k.range), ncol=1),
        test_iter  = matrix(run_test_iter,  nrow=length(k.range), ncol=1)
      )
    }
  
  total_time <- Sys.time() - start_time
  
  if (allowParallel) {
    parallel::stopCluster(cl)
  }
  
  # res_list is now a list with $cpi, $train_iter, $test_iter, each dimension
  # (length(k.range), n.runs)
  final_cpi        <- res_list$cpi
  final_train_iter <- res_list$train_iter
  final_test_iter  <- res_list$test_iter
  
  dimnames(final_cpi)        <- list(paste0("k", k.range), paste0("run", 1:n.runs))
  dimnames(final_train_iter) <- list(paste0("k", k.range), paste0("run", 1:n.runs))
  dimnames(final_test_iter)  <- list(paste0("k", k.range), paste0("run", 1:n.runs))
  
  if (make.plot) {
    # Simple plot of CPI
    dev.new(width = 4, height = 5)
    plot(k.range, final_cpi[, 1], 
         ylim = c(min(final_cpi, na.rm=TRUE), max(final_cpi, na.rm=TRUE)), 
         pch = 20, main = "", xlab = "k", ylab = "CPI")
    for (m in 2:n.runs) {
      points(k.range, final_cpi[, m], pch = 20)
    }
    lines(k.range, apply(final_cpi, 1, mean, na.rm=TRUE), col = "red", lwd = 2)
    mtext("Optimum k", outer=TRUE, cex=1, line=-2)
  }
  
  # Construct final result object
  results_out <- list(
    cpi            = final_cpi,
    train_iter     = final_train_iter,
    test_iter      = final_test_iter,
    hyperparams    = list(
      n.runs        = n.runs,
      n.fold        = n.fold,
      k.range       = k.range,
      maxiter       = maxiter,
      lr            = lr,
      tol           = tol,
      allowParallel = allowParallel,
      n.cores       = n.cores,
      seed          = seed
    ),
    total_runtime  = total_time
  )
  
  if (result) {
    return(results_out)
  } else {
    return(invisible(results_out))
  }
}