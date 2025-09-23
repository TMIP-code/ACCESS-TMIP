#!/bin/bash

#PBS -P y99
#PBS -N yrlyTM_OM2-025-CMIP6
#PBS -l ncpus=28
#PBS -l mem=120GB
#PBS -l jobfs=4GB
#PBS -l walltime=3:00:00
#PBS -l storage=scratch/gh0+gdata/xv83+scratch/xv83+scratch/y99+gdata/cj50+gdata/ik11
#PBS -l wd
#PBS -o output/PBS/
#PBS -j oe

# qsub -I -P y99 -N yearlyTM_OM2-1 -l ncpus=28 -l mem=120GB -l jobfs=4GB -l walltime=3:00:00 -l storage=scratch/gh0+gdata/xv83+scratch/xv83+scratch/y99+gdata/cj50+gdata/ik11 -l wd -o output/PBS/ -j oe

echo "Going into ACCESS-TMIP"
cd ~/Projects/TMIP/ACCESS-TMIP

julia src/build_yearlyTM_ACCESS-OM2-025_CMIP6.jl &> output/build_yearlyTM_ACCESS-OM2-025_CMIP6.$PBS_JOBID.out


