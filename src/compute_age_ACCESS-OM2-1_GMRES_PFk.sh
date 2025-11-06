#!/bin/bash

#PBS -P y99
#PBS -N age_OM2-1_GMRES_PFk
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
julia src/compute_age_ACCESS-OM2-1_GMRES_PFk.jl &> output/compute_age_ACCESS-OM2-1_GMRES_PFk.$PBS_JOBID.out
