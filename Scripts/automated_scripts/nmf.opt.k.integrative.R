nmf.opt.k.integrative <- function(dat, is.binary, n.runs = 30, n.fold = 5, k.range = 2:8, 
                                  result = TRUE, make.plot = TRUE, progress = TRUE, 
                                  maxiter = 100, lr = 1e-3, tol = 1e-6, 
                                  allowParallel = FALSE, n.cores = NULL, seed = 12345) {
  # Required libraries:
  library(mclust) # for adjustedRandIndex
  library(foreach)
  library(doParallel)
  library(MASS) # for ginv
  
  if (!is.list(dat)) stop("Input 'dat' must be a list of matrices.")
  M <- length(dat)
  if (length(is.binary) != M) stop("Length of 'is.binary' must match the number of modalities in 'dat'.")
  
  for (i in seq_len(M)) {
    if (min(dat[[i]]) < 0) {
      dat[[i]] <- pmax(dat[[i]] + abs(min(dat[[i]])), 0) + .Machine$double.eps
    }
  }
  
  # Calculate weights based on proportions
  M_b <- sum(is.binary)
  M_c <- M - M_b
  if (M_b > 0 && M_c > 0) {
    total_binary_weight <- M_b / M
    total_cont_weight <- M_c / M
    wt <- numeric(M)
    wt[is.binary] <- total_binary_weight / M_b
    wt[!is.binary] <- total_cont_weight / M_c
  } else {
    # All binary or all continuous
    wt <- rep(1/M, M)
  }
  
  sigmoid <- function(x) 1 / (1 + exp(-x))
  
  # Training fit function
  nmf.integrative.fit <- function(dat.list, k, maxiter, lr, tol, wt, is.binary, seed) {
    n <- nrow(dat.list[[1]])
    set.seed(seed)
    W <- matrix(runif(n * k, min = 0, max = 1), n, k)
    H.list <- vector("list", M)
    for (m in seq_len(M)) {
      H.list[[m]] <- matrix(runif(k * ncol(dat.list[[m]]), min = 0, max = 1), k, ncol(dat.list[[m]]))
    }
    
    prev_W <- W
    for (iter in seq_len(maxiter)) {
      # Update H
      for (m in seq_len(M)) {
        X_m <- dat.list[[m]]
        pred_m <- W %*% H.list[[m]]
        if (is.binary[m]) {
          Sigm <- sigmoid(pred_m)
          Grad_H <- wt[m] * t(W) %*% (Sigm - X_m) 
          H.list[[m]] <- pmax(H.list[[m]] - lr * Grad_H, 0)
        } else {
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
    list(W = W, H = H.list, clusters = clusters)
  }
  
  # Test W estimation given H from training:
  # We fix H and solve for W by gradient descent, mixing binary and continuous losses.
  nmf.test.W <- function(d.test, H.list, wt, is.binary, lr, tol, maxiter, seed) {
    n.test <- nrow(d.test[[1]])
    k <- nrow(H.list[[1]])
    set.seed(seed+1)
    W_test <- matrix(runif(n.test * k, min=0, max=1), n.test, k)
    prev_W <- W_test
    
    for (iter in seq_len(maxiter)) {
      W_grad <- matrix(0, n.test, k)
      for (m in seq_len(M)) {
        X_m <- d.test[[m]]
        pred_m <- W_test %*% H.list[[m]]
        if (is.binary[m]) {
          Sigm <- 1/(1+exp(-pred_m))
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
    W_test
  }
  
  if (allowParallel) {
    if (is.null(n.cores)) n.cores <- parallel::detectCores() - 1
    cl <- parallel::makeCluster(n.cores)
    doParallel::registerDoParallel(cl)
  }
  
  set.seed(seed)
  n.sample <- nrow(dat[[1]])
  CPI <- matrix(NA, length(k.range), n.runs)
  dimnames(CPI) <- list(paste("k", k.range, sep = ""), paste("run", 1:n.runs, sep = ""))
  
  res_list <- foreach::foreach(i = 1:n.runs, .combine = 'cbind', .packages = c("mclust")) %dopar% {
    run_results <- numeric(length(k.range))
    for (ki in seq_along(k.range)) {
      k <- k.range[ki]
      R.ind <- NULL
      random.sample <- sample(seq(n.sample), n.sample)
      fold_size <- floor(n.sample/n.fold)
      
      for (j in 1:n.fold) {
        test_idx <- ((j-1)*fold_size+1):min(j*fold_size, n.sample)
        test.sample <- random.sample[test_idx]
        train.sample <- setdiff(random.sample, test.sample)
        
        d.train <- lapply(dat, function(x) x[train.sample, , drop = FALSE])
        d.test <- lapply(dat, function(x) x[test.sample, , drop = FALSE])
        
        # Train model on training set
        fit.train <- nmf.integrative.fit(d.train, k = k, maxiter = maxiter, lr = lr, tol = tol, wt = wt, is.binary = is.binary, seed = seed)
        
        # Compute test W by fixing H and doing iterative updates
        W.predict <- nmf.test.W(d.test, fit.train$H, wt, is.binary, lr, tol, maxiter, seed)
        predicted.cluster.mem <- apply(W.predict, 1, which.max)
        
        # Fit model on test set for ground truth clusters
        fit.test <- nmf.integrative.fit(d.test, k = k, maxiter = maxiter, lr = lr, tol = tol, wt = wt, is.binary = is.binary, seed = seed+10)
        computed.cluster.mem <- fit.test$clusters
        
        R.ind <- c(R.ind, mclust::adjustedRandIndex(predicted.cluster.mem, computed.cluster.mem))
        
        if (progress) {
          done <- ((i-1)*length(k.range)*n.fold + (ki-1)*n.fold + j) / (n.runs*length(k.range)*n.fold)
          pct <- round(done*100)
          if (pct %in% seq(5,100,by=5)) {
            message(pct, "% complete")
            flush.console()
          }
        }
      }
      run_results[ki] <- mean(R.ind)
    }
    run_results
  }
  
  if (allowParallel) {
    parallel::stopCluster(cl)
  }
  
  for (i in 1:ncol(res_list)) {
    CPI[, i] <- res_list[, i]
  }
  
  if (make.plot) {
    dev.new(width = 4, height = 5)
    plot(k.range, CPI[, 1], ylim = c(min(CPI, na.rm=TRUE), max(CPI, na.rm=TRUE)), pch = 20, main = "", xlab = "k", ylab = "CPI")
    for (m in 2:n.runs) points(k.range, CPI[, m], pch = 20)
    lines(k.range, apply(CPI, 1, mean, na.rm=TRUE), col = "red", lwd = 2)
    mtext("Optimum k", outer = TRUE, cex = 1, line = -2)
  }
  
  if (result) return(CPI)
}
