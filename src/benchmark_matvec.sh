#!/bin/bash

#PBS -P y99
#PBS -N benchmark_ncpusplaceholder_sparseBLASplaceholder
#PBS -l ncpus=ncpusplaceholder
#PBS -q normal
#PBS -l mem=190GB
#PBS -l jobfs=4GB
#PBS -l walltime=00:10:00
#PBS -l storage=scratch/gh0+scratch/y99+gdata/xp65
#PBS -l wd
#PBS -o output/PBS/
#PBS -j oe

# export OPENBLAS_NUM_THREADS=1 # or set the number of OpenBLAS threads?
export JULIA_NUM_THREADS=$PBS_NCPUS   # Set the number of Julia threads?
export MKL_NUM_THREADS=$PBS_NCPUS   # or set the number of MKL threads?

export sparseBLAS=sparseBLASplaceholder
# export sparseBLAS=Julia  # Julia default
# export sparseBLAS=MKL    # Use MKLSparse
# export sparseBLAS=CSR    # Use ThreadedSparseCSR (row-major format)

echo "Going into ACCESS-TMIP"
cd ~/Projects/TMIP/ACCESS-TMIP

echo "Running script"
julia src/benchmark_matvec.jl &> output/benchmark9_matvec.$PBS_JOBID.$sparseBLAS-BLAS.$PBS_NCPUS-CPUs.out
