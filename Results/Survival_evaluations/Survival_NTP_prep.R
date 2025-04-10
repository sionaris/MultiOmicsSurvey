library(MOVICS)
library(TCGAbiolinks)

# Ensure reproducibility
RNGversion("4.2.2")
set.seed(123)

# Load custom helper functions
source("Scripts/automated_scripts/custom_functions.R")
source("Scripts/automated_scripts/modified_MOVICS_functions.R")

# NTP loop #####
for (algorithm in length(algorithms)) {
  # load environment
  
  
  # get as many templates as possible
  dgea.marker.up_1000 <- runMarker_single_algorithm_no_export(algorithm_name = algorithm,
                                                              moic.res = plot_object,
                                                              n.marker = 1000,
                                                              dea.method    = "limma", # name of DEA method
                                                              prefix        = "dgea_", # MUST be the same of argument in runDEA()
                                                              dat.path      = paste0(home, "/Results/single_algorithm/", algorithm), # path of DEA files
                                                              p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                                                              p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                                                              norm.expr = plotdata$RNAseq,
                                                              dirct         = "up" # direction of dysregulation in expression
  )
  
  # NTP
  transNEO_ntp_expr_up = runNTP(
    expr = transcr,
    templates = dgea.marker.up_1000$templates,
    scaleFlag = TRUE,
    centerFlag = TRUE,
    nPerm = 10000,
    seed = 123,
    distance = "cosine", # default
    doPlot = TRUE,
    height = 8,
    width = 12,
    fig.path = paste0(home, "/Results/Survival_evaluations"),
    fig.name = paste0(algorithm, "_ntp_expr_up_heatmap_TCGA"))
}
