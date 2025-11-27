#!/bin/bash

#PBS -P y99
#PBS -N age_OM2-1_GMRES_ILU
#PBS -l ncpus=24
#PBS -q normal
#PBS -l mem=24GB
#PBS -l jobfs=4GB
#PBS -l walltime=03:00:00
#PBS -l storage=scratch/gh0+scratch/y99+gdata/xp65
#PBS -l wd
#PBS -o output/PBS/
#PBS -j oe

# export JULIA_NUM_THREADS=$PBS_NCPUS    # Set the number of Julia threads?
# export OPENBLAS_NUM_THREADS=1 # or set the number of OpenBLAS threads?
export MKL_NUM_THREADS=$PBS_NCPUS    # or set the number of MKL threads?

echo "Going into ACCESS-TMIP"
cd ~/Projects/TMIP/ACCESS-TMIP

echo "Running script"
julia src/compute_age_ACCESS-OM2-1_GMRES_ILU.jl &> output/compute_age_ACCESS-OM2-1_GMRES_ILU.$PBS_JOBID.out
