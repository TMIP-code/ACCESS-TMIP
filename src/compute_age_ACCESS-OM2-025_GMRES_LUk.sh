#!/bin/bash

#PBS -P y99
#PBS -N age_OM2-025_GMRES_LUk
#PBS -l ncpus=24
#PBS -q hugemem
#PBS -l mem=735GB
#PBS -l jobfs=4GB
#PBS -l walltime=24:00:00
#PBS -l storage=scratch/gh0+scratch/y99+gdata/xp65
#PBS -l wd
#PBS -o output/PBS/
#PBS -j oe

echo "Going into ACCESS-TMIP"
cd ~/Projects/TMIP/ACCESS-TMIP

echo "Running script"
julia src/compute_age_ACCESS-OM2-025_GMRES_LUk.jl &> output/compute_age_ACCESS-OM2-025_GMRES_LUk.$PBS_JOBID.out
