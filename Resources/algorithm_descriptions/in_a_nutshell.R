# A collection of the high-level descriptions of algorithms
desc_list = list()

desc_list[["MOVICS"]] = "MOVICS is used to perform multi-omics integrative clustering and visualization for cancer subtyping research. It provides a unified interface for 10 state-of-the-art multi-omics clustering algorithms, and standardizes the output for each algorithm so as to form a pipeline for downstream analyses. Ten algorithms are CIMLR, iClusterBayes, MoCluster, COCA, ConsensusClustering, IntNMF, LRAcluster, NEMO, PINSPlus, and SNF, where the former three methods can also perform the process of feature selection. For cancer subtyping studies, MOVICS also forms a pipeline for most commonly used downstream analyses for further subtype characterization and creates editable publication-quality illustrations."

# SNF
desc_list[["SNF"]] = "SNF approaches the multi-omic problem by constructing networks of samples (e.g., patients) for each available data type and then efficiently fusing these into one network that represents the full spectrum of underlying data."

# CIMLR
desc_list[["CIMLR"]] = ""

# PINSPlus
desc_list[["PINSPlus"]] = ""

# NEMO
desc_list[["NEMO"]] = "NEMO processes data by taking matrices from multiple omics types and calculating similarity matrices for each type using a radial basis function kernel. These matrices consider the Euclidean distance between sample profiles and adjust for density differences using a normalization factor (identical to the SNF approach until here). NEMO then creates a relative similarity matrix for each omic to adjust similarities based on local neighborhood data, making the comparison between different omics more consistent. The algorithm averages these relative similarity matrices to form an **Average Relative Similarity matrix**. This Average Relative Similarity matrix, which treats similarities as transition probabilities (akin to a random walk on a graph), is used for spectral clustering. The number of clusters is determined using a variant of the eigengap method, optimized to enhance the prognostic value by possibly suggesting a higher number of clusters."

# IntNMF
desc_list[["IntNMF"]] = "Integrative Non-negative Matrix Factorization (IntNMF) is an NMF approach to multi-omics where one common basis is used for all modalities, but separate sets of coefficients for each modality's features. Through iterative updates, a final set of basis vectors (each of which represents a sample's association with a distinct cluster) is produced, and each sample is allocated a label corresponding to the index of the basis vector with the highest value for that sample. The first $k$ basis vectors are examined, where $k$ is specified by the user or is determined based on the Cluster Prediction Index (CPI) of the examined $k$'s."

# iClusterBayes
desc_list[["iClusterBayes"]] = "iClusterBayes uses a few latent variables to capture the inherent structure of multiple omics datasets to achieve joint dimension reduction. As a result, the tumor samples can be clustered in the latent variable space and relevant omics features that drive the sample clustering are identified through Bayesian variable selection."