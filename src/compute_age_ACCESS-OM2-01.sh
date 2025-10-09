#!/bin/bash

#PBS -P y99
#PBS -N age_OM2-01
#PBS -l ncpus=48
#PBS -q megamem
#PBS -l mem=2900GB
#PBS -l jobfs=4GB
#PBS -l walltime=48:00:00
#PBS -l storage=scratch/gh0+scratch/y99+scratch/p66
#PBS -l wd
#PBS -o output/PBS/
#PBS -j oe



echo "Going into ACCESS-TMIP"
cd ~/Projects/TMIP/ACCESS-TMIP

echo "Running script"
julia src/compute_age_ACCESS-OM2-01.jl &> output/$PBS_JOBID.compute_age_ACCESS-OM2-01.out
