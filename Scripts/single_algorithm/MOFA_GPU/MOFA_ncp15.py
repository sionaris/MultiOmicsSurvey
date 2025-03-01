
import glob
import os
import pandas as pd
import numpy as np
from mofapy2.run.entry_point import entry_point
import cupy as cp

print(cp.cuda.runtime.runtimeGetVersion())
print(cp.show_config())

# ------------------------------------------------------------------------------
# 1. Import CSV files, each file is treated as one 'view' in MOFA
# ------------------------------------------------------------------------------
input_folder = "Omics_input"
csv_files = glob.glob(os.path.join(input_folder, "*.csv"))

# Read each CSV: first column is feature names and subsequent columns are samples
data_matrices = []
view_names = []
for file_path in csv_files:
    view_name = os.path.splitext(os.path.basename(file_path))[0] # give appropriate names to views
    view_names.append(view_name)

    # Read the CSV file; first column = row labels (features), rest of columns = sample measurements
    df = pd.read_csv(file_path, header=0, index_col=0)

    # Transpose so that rows = samples, columns = features
    mat = df.T.values  # shape: [number_of_samples, number_of_features]
    data_matrices.append(mat)

# Wrap in a nested list to signal that we have M views and 1 group: data_mat[m][g]
# For one group (g=0), each "view" is data_mat[m][0]
data_mat_nested = [[m] for m in data_matrices]

# ------------------------------------------------------------------------------
# 2. Initialize the MOFA entry point and set data options
# ------------------------------------------------------------------------------
ent = entry_point()
ent.set_data_options(
    scale_views=False,
    scale_groups=False,
    center_groups=False
)

# ------------------------------------------------------------------------------
# 3. Add data to the model
# ------------------------------------------------------------------------------
# Provide one likelihood for each view
likelihoods = ["gaussian","gaussian","bernoulli","gaussian","gaussian"]

ent.set_data_matrix(
    data_mat_nested,
    likelihoods=likelihoods,
    views_names=view_names
)

# ------------------------------------------------------------------------------
# 4. Set model options
# ------------------------------------------------------------------------------
ent.set_model_options(
    factors=15,
    spikeslab_weights=False,
    spikeslab_factors=False,
    ard_weights=True,
    ard_factors=False
)

# ------------------------------------------------------------------------------
# 5. Set training options
# ------------------------------------------------------------------------------
ent.set_train_options(
    convergence_mode="slow",
    dropR2=-1,
    iter=20000,
    gpu_mode=True,
    freqELBO=5,
    startELBO=1,
    seed=123
)

# ------------------------------------------------------------------------------
# 6. Build and train the model
# ------------------------------------------------------------------------------
ent.build()
ent.run()

# ------------------------------------------------------------------------------
# 7. Save the trained model output
# ------------------------------------------------------------------------------
ent.save(outfile="MOFA_output_factors15_iter20000_R2_-1.hdf5")

print("MOFA training complete. Results saved to \"MOFA_output_factors15_iter20000_R2_-1.hdf5\"")
