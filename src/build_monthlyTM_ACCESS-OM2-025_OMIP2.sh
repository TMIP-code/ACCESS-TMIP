#!/bin/bash

#PBS -P y99
#PBS -N monthlyTM_OM2-025
#PBS -l ncpus=24
#PBS -l mem=95GB
#PBS -l jobfs=4GB
#PBS -l walltime=1:00:00
#PBS -l storage=scratch/gh0+gdata/xp65+gdata/xv83+scratch/xv83+scratch/y99+gdata/cj50+gdata/ik11
#PBS -l wd
#PBS -o output/PBS/
#PBS -j oe

echo "Going into ACCESS-TMIP"
cd ~/Projects/TMIP/ACCESS-TMIP

julia src/build_monthlyTM_ACCESS-OM2-025_OMIP2.jl &> output/build_monthlyTM_ACCESS-OM2-025_OMIP2.$PBS_JOBID.out


