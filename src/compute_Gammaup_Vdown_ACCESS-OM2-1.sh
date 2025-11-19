#!/bin/bash

#PBS -P y99
#PBS -N Gup_Vdown_OM2-1
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
# julia src/compute_Gammaup_Vdown_ACCESS-OM2-1.jl &> output/compute_Gammaup_Vdown_ACCESS-OM2-1.$PBS_JOBID.out
julia src/compute_Gammaup_Vdown_ACCESS-OM2-1_nosplit.jl &> output/compute_Gammaup_Vdown_ACCESS-OM2-1_nosplit.$PBS_JOBID.out
