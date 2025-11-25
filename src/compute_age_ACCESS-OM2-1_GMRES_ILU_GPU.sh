#!/bin/bash

#PBS -P y99
#PBS -N GPU_test
#PBS -q gpuvolta
#PBS -l ncpus=24
#PBS -l ngpus=2
#PBS -l mem=192GB
#PBS -l jobfs=4GB
#PBS -l walltime=01:00:00
#PBS -l storage=scratch/y99+gdata/xp65
#PBS -l wd
#PBS -o output/PBS/
#PBS -j oe

# module load cuda/12.9.0 # Maybe I don;t need that (downloaded by Julia CUDA package?)

echo "Going into ACCESS-TMIP"
cd ~/Projects/TMIP/ACCESS-TMIP

echo "Running script"
julia src/compute_age_ACCESS-OM2-1_GMRES_ILU_GPU.jl &> output/compute_age_ACCESS-OM2-1_GMRES_ILU_GPU.$PBS_JOBID.out
