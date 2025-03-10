from scipy import stats
import pandas as pd
import numpy as np

def norm2distance(df):
    temp = (df ** 2).sum(axis=1, keepdims=True).repeat(df.shape[0], axis=1)
    distance = temp + temp.T - 2 * df @ df.T
    return distance

def similarity(df, k=20, sig=0.5, ktop=True, high=True, input_type="omic_matrices"):
    """
    similarity() can operate in two modes:
      1) input_type="omic_matrices" (default): 
         'df' is a (samples x features) data frame; Eucl. distance is computed.
      2) input_type="symmetric_relationships": 
         'df' is already a (samples x samples) symmetrical distance matrix. 
         We'll skip norm2distance() and proceed with 
         top-k filtering, Gaussian-like weighting, etc.

    df: pd.DataFrame
    k: top k neighbors
    sig: scaling factor
    ktop: whether to keep only k neighbors
    high: whether to do high_order(s)
    input_type: "omic_matrices" or "symmetric_relationships"
    """
    assert isinstance(df, pd.DataFrame)
    samples = df.index

    if input_type == "omic_matrices":
        # df is shape (n_samples, n_features)
        df = df.values
        distance = norm2distance(df)
        # Make sure it's symmetric
        distance = (distance + distance.T) / 2

    elif input_type == "symmetric_relationships":
        # df is shape (n_samples, n_samples) (already a distance or relationship matrix)
        # We assume it's already symmetrical, but let's symmetrize just to be safe:
        distance = df.values
        distance = (distance + distance.T) / 2

    else:
        raise ValueError("input_type must be either 'omic_matrices' or 'symmetric_relationships'.")

    # Next lines are the local-scaling / Gaussian weighting logic
    knear = distance.argsort()[:, 1:k+1]
    x = np.repeat(np.arange(distance.shape[0]).reshape(-1,1), k, axis=1)
    delta0 = distance[x, knear].mean(axis=1)
    delta = delta0.reshape(-1,1).repeat(distance.shape[0], axis=1)
    delta = sig * (delta + delta.T + distance) / 3.0

    # Convert distances to similarities via Gaussian pdf
    s = stats.norm(0, delta).pdf(distance)

    # Keep only top-k neighbors if ktop=True
    if ktop:
        ag = s.argsort()[:, :-(k + 1)]
        x = np.repeat(np.arange(ag.shape[0]).reshape(-1, 1), ag.shape[1], axis=1)
        s[x, ag] = 0
    s = (s + s.T) / 2

    # Now row-normalize (assuming you want each row to sum to 1)
    s = normalize(s)

    # Possibly do high-order transitions
    if high:
        s = high_order(s)

    s = pd.DataFrame(s, index=samples, columns=samples)
    return s

def normalize(s):
    ispd=isinstance(s,pd.DataFrame)
    if(ispd):
        index=s.index
        col=s.columns
        s=s.values
    s=s/s.sum(axis=1).reshape(-1,1)
    if(ispd):
        s=pd.DataFrame(s,index=index,columns=col)
    return (s+s.T)/2

def high_order(s,n=5):
    for i in range(n):
        s=s@s
    s=(s+s.T)/2
    return normalize(s)

def partition_num(num, workers):
    if num % workers == 0:
        return [num//workers]*workers
    else:
        return [num//workers]*workers + [num % workers]
