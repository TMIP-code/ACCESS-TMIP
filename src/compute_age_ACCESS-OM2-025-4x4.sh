#!/bin/bash

#PBS -P y99
#PBS -N age_OM2-025-4x4
#PBS -l ncpus=24
#PBS -q normal
#PBS -l mem=95GB
#PBS -l jobfs=4GB
#PBS -l walltime=01:00:00
#PBS -l storage=scratch/gh0+scratch/y99+scratch/p66
#PBS -l wd
#PBS -o output/PBS/
#PBS -j oe



echo "Going into ACCESS-TMIP"
cd ~/Projects/TMIP/ACCESS-TMIP

echo "Running script"
julia src/compute_age_ACCESS-OM2-025-4x4.jl &> output/$PBS_JOBID.compute_age_ACCESS-OM2-025-4x4.out
