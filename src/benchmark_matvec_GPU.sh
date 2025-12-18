#!/bin/bash

#PBS -P y99
#PBS -N benchmark_GPU_matvec
#PBS -q gpuvolta
#PBS -l ncpus=12
#PBS -l ngpus=1
#PBS -l mem=96GB
#PBS -l jobfs=4GB
#PBS -l walltime=00:10:00
#PBS -l storage=scratch/gh0+scratch/y99+gdata/xp65
#PBS -l wd
#PBS -o output/PBS/
#PBS -j oe

echo "Going into ACCESS-TMIP"
cd ~/Projects/TMIP/ACCESS-TMIP

echo "Running script"
julia src/benchmark_matvec_GPU.jl &> output/benchmark9_matvec_GPU.$PBS_JOBID.out
