#!/bin/bash

#PBS -P xv83
#PBS -N time_step_injection
#PBS -l ncpus=48
#PBS -q hugemem
#PBS -l mem=360GB
#PBS -l jobfs=4GB
#PBS -l walltime=05:00:00
#PBS -l storage=scratch/gh0+scratch/xv83
#PBS -l wd
#PBS -o output/PBS/
#PBS -j oe

experiment=historical
member=r6i1p1f1
time_window=Jan1990-Dec1999

echo "Going into ACCESS-TMIP"
cd ~/Projects/TMIP/ACCESS-TMIP

echo "Running script"
julia src/timestep_Carbon_injection.jl $experiment $member $time_window $lonlat &> output/timestep_Carbon_injection.$PBS_JOBID.out
