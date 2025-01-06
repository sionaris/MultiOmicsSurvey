# Adapted from MONET's official repository:
# https://github.com/Shamir-Lab/MONET/blob/master/R_code/monet_exp.R
run.clusternomics <- function(omics.list, num.clusters=NULL, 
                              num.clusters.per.omic=NULL, dataDistributions = NULL,
                              ncores = NULL) {
  library(clusternomics)
  if (is.null(num.clusters)) {
    stop("Provide a valid non-NULL/NA num.clusters value!")
  }
  if (is.null(num.clusters.per.omic)) {
    stop("Provide a valid non-NULL/NA num.clusters.per.omic value! It must be a list of two vectors.")
  }
  if(is.null(ncores)) {
    stop("ncores must be a numeric value.")
  }
  if (is.null(dataDistributions)) {
    stop("Provide a valid non-NULL/NA value for dataDistributions.")
  } else if (!is.null(dataDistributions) && length(dataDistributions) != length(omics.list)) {
    stop("length(omics.list) == length(dataDistributions) is essential.")
  }
  all.param.options = expand.grid(1:length(num.clusters), 1:length(num.clusters.per.omic))
  all.rets = mclapply(1:nrow(all.param.options), function(i) {
    set.seed(123 + i)
    cur.num.clusters = num.clusters[all.param.options[i, 1]]
    cur.num.clusters.per.omic = num.clusters.per.omic[[all.param.options[i, 2]]]
    start = Sys.time()
    if (length(cur.num.clusters.per.omic) == 1) {
      num.clusters.per.omic = rep(num.clusters.per.omic, length(omics.list))
    }
    cluster.counts = list(global=cur.num.clusters, context=cur.num.clusters.per.omic)
    
    # Hack due to a false assertion in clusternomics package
    if (length(omics.list) > 2) {
      cluster.counts = c(cluster.counts, rep('UNUSED', length(omics.list) - 2))
    }
    # parameters are as used in the clusternomics publication.
    results = clusternomics::contextCluster(omics.list.trans, cluster.counts, maxIter=1e4, 
                             burnin=5e3, lag=3, dataDistributions=dataDistributions, verbose=T)
    cur.clustering = results$samples[[length(results$samples)]]$Global
    dic = results$DIC
    time.taken.per.param = as.numeric(Sys.time() - start, units='secs')
    return(list(clustering=cur.clustering, dic=dic, timing=time.taken.per.param, clusternomics.ret=results))
  }, mc.cores=ncores)
  
  best.sol.index = which.min(sapply(all.rets, function(x) x$dic))
  wall.timing = max(sapply(all.rets, function(x) x$timing)) + time.taken.normalization
  return(list(clustering=all.rets[[best.sol.index]]$clustering, timing=wall.timing, all.rets=all.rets))
}