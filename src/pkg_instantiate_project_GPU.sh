#!/bin/bash

#PBS -N instantiate_GPU
#PBS -P y99
#PBS -l mem=47GB
#PBS -q gpuvolta
#PBS -l walltime=01:00:00
#PBS -l ncpus=12
#PBS -l ngpus=1
#PBS -l storage=scratch/y99+gdata/y99
#PBS -l jobfs=4GB
#PBS -o scratch_output/PBS/
#PBS -e scratch_output/PBS/
#PBS -l wd

repo_root=/home/561/bp3051/Projects/TMIP/ACCESS-TMIP
echo "REPO_ROOT=$repo_root"

echo "Instantiating packages on compute node on GPU"
module load cuda/12.9.0
export JULIA_CUDA_USE_COMPAT=false
julia --project -e 'using Pkg; Pkg.instantiate()'
echo "Done instantiating packages on compute node on GPU"
