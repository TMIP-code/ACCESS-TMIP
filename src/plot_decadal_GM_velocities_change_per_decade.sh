#!/bin/bash

#PBS -P xv83
#PBS -N plot_decadal_GM_archive
#PBS -l ncpus=28
#PBS -l mem=180GB
#PBS -l jobfs=4GB
#PBS -l walltime=12:00:00
#PBS -l storage=gdata/xv83
#PBS -l wd
#PBS -o output/PBS/
#PBS -j oe

echo "Going into ACCESS-TMIP"
cd ~/Projects/TMIP/ACCESS-TMIP

echo "Running script"
julia src/plot_decadal_GM_velocities_change_per_decade.jl &> output/$PBS_JOBID.plot_decadal_GM_velocities_change_per_decade.out
