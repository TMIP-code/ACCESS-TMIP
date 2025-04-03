#!/bin/bash

#PBS -P xv83
#PBS -N plot_TTD
#PBS -l ncpus=6
#PBS -l mem=32GB
#PBS -l jobfs=4GB
#PBS -l walltime=1:00:00
#PBS -l storage=scratch/gh0+scratch/xv83
#PBS -l wd
#PBS -o output/PBS/
#PBS -j oe

echo "Going into ACCESS-TMIP"
cd ~/Projects/TMIP/ACCESS-TMIP

echo "Running script"
julia src/plot_timeseries_map_TTD_seqeff.jl &> output/plot_timeseries_map_TTD_seqeff.$PBS_JOBID.out
