#!/bin/bash

#PBS -P xv83
#PBS -N ACCESS_Gamma
#PBS -l ncpus=28
#PBS -l mem=180GB
#PBS -l jobfs=4GB
#PBS -l walltime=12:00:00
#PBS -l storage=scratch/xv83
#PBS -l wd
#PBS -o output/PBS/
#PBS -j oe

echo "Going into ACCESS-TMIP"
cd ~/Projects/TMIP/ACCESS-TMIP

echo "Running script"
# julia src/compute_Gammas.jl
julia src/compute_and_plot_Gammas_varied_diffusivities.jl &> output/$PBS_JOBID.compute_and_plot_Gammas_varied_diffusivities.out
