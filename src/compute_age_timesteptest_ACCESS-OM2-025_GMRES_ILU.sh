#!/bin/bash

#PBS -P y99
#PBS -N age_step_test_GMRES_ILU
#PBS -l ncpus=48
#PBS -q normal
#PBS -l mem=190GB
#PBS -l jobfs=4GB
#PBS -l walltime=02:00:00
#PBS -l storage=scratch/gh0+scratch/y99+gdata/xp65
#PBS -l wd
#PBS -o output/PBS/
#PBS -j oe

# export OPENBLAS_NUM_THREADS=1 # or set the number of OpenBLAS threads?
export JULIA_NUM_THREADS=$PBS_NCPUS    # Set the number of Julia threads?
export MKL_NUM_THREADS=$PBS_NCPUS    # or set the number of MKL threads?

export ILUTAU=1e-8 # drop tolerance for ILU preconditioner

echo "Going into ACCESS-TMIP"
cd ~/Projects/TMIP/ACCESS-TMIP

echo "Running script"
julia src/compute_age_timesteptest_ACCESS-OM2-025_GMRES_ILU.jl &> output/compute_age_timesteptest_ACCESS-OM2-025_GMRES_ILU.$PBS_JOBID.$ILUTAU.out
