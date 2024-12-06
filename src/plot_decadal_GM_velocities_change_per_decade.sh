#!/bin/bash

#PBS -P xv83
#PBS -N plot_decadal_GM_archive
#PBS -l ncpus=24
#PBS -l mem=90GB
#PBS -l jobfs=4GB
#PBS -l walltime=24:00:00
#PBS -l storage=scratch/gh0+gdata/xv83+scratch/xv83
#PBS -l wd
#PBS -o output/PBS/
#PBS -j oe

echo "Going into ACCESS-TMIP"
cd ~/Projects/TMIP/ACCESS-TMIP

echo "Running script"
julia src/plot_decadal_GM_velocities_change_per_decade.jl &> output/$PBS_JOBID.plot_decadal_GM_velocities_change_per_decade.out
