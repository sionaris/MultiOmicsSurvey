#!/usr/bin/env bash
module load rhel8/default-amp
module load cuda/11.4

source ~/miniconda/etc/profile.d/conda.sh

conda create -n mofa_env_gpu -y python=3.10 pip
conda activate mofa_env_gpu

# CuPy matched for CUDA 11.x
pip install cupy-cuda11x

# MOFA + GPU mem querying
pip install mofapy2==0.7.2 nvidia-ml-py3 pandas numpy

