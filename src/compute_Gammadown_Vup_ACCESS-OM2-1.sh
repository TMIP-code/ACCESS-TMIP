#!/bin/bash

#PBS -P y99
#PBS -N Gdown_Vup_OM2-1
#PBS -l ncpus=48
#PBS -q normal
#PBS -l mem=190GB
#PBS -l jobfs=4GB
#PBS -l walltime=02:00:00
#PBS -l storage=scratch/gh0+scratch/y99+gdata/xp65
#PBS -l wd
#PBS -o output/PBS/
#PBS -j oe


echo "Going into ACCESS-TMIP"
cd ~/Projects/TMIP/ACCESS-TMIP

echo "Running script"
julia src/compute_Gammadown_Vup_ACCESS-OM2-1.jl &> output/compute_Gammadown_Vup_ACCESS-OM2-1.$PBS_JOBID.out
