#!/bin/bash

#PBS -P y99
#PBS -N Gup_Vdown_OM2-025
#PBS -l ncpus=48
#PBS -q megamem
#PBS -l mem=2900GB
#PBS -l jobfs=4GB
#PBS -l walltime=48:00:00
#PBS -l storage=scratch/gh0+scratch/y99+gdata/xp65
#PBS -l wd
#PBS -o output/PBS/
#PBS -j oe


echo "Going into ACCESS-TMIP"
cd ~/Projects/TMIP/ACCESS-TMIP

echo "Running script"
julia src/compute_Gammaup_Vdown_ACCESS-OM2-025_periodic.jl &> output/compute_Gammaup_Vdown_ACCESS-OM2-025_periodic.$PBS_JOBID.out
