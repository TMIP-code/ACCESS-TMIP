#!/bin/bash

#PBS -P xv83
#PBS -N Gamma
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

# model=ACCESS1-0
# model=ACCESS1-3
# model=ACCESS-CM2
# model=ACCESS-ESM1-5
# model=CESM2
# model=CESM2-FV2
# model=CESM2-WACCM-FV2
model=CMCC-CM2-HR4
# model=CMCC-CM2-SR5
# model=CMCC-ESM2
# model=FGOALS-f3-L
# model=FGOALS-g3
# model=MPI-ESM-1-2-HAM
# model=MPI-ESM1-2-HR
# model=MPI-ESM1-2-LR
# model=NorCPM1
# model=NorESM2-LM
# model=NorESM2-MM

echo "Running script"
julia src/compute_Gammas.jl $model &> output/$PBS_JOBID.$model.compute_Gammas.out
# julia src/compute_and_plot_Gammas_varied_diffusivities.jl &> output/$PBS_JOBID.compute_and_plot_Gammas_varied_diffusivities.out
