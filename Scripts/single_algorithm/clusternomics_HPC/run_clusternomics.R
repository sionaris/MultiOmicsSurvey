library(clusternomics)
source("modified_clusternomics.R")

# load input
input = readRDS("clusternomics_input.rds")

# Run the modified algorithm
clusternomics = run.clusternoics(omics.list = input,
                                 num.clusters = 2:10,
                                 num.clusters.per.omic = 5, # MONET uses  list(rep(3, length(omics.list)), rep(5, length(omics.list))) as default
                                 dataDistributions = c('binary', 'diagNormal', 
                                                       'diagNormal', 'diagNormal', 
                                                       'diagNormal'),
                                 ncores = 30)

# Write out output
saveRDS("clusternomics_output.rds")
writeLines(capture.output(sessionInfo()), "clusternomics.sessionInfo.txt")
