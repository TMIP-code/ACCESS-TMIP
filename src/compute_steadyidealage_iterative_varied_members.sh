#!/bin/bash

#PBS -P xv83
#PBS -N age_member_placeholder
#PBS -l ncpus=24
#PBS -l mem=100GB
#PBS -l jobfs=4GB
#PBS -l walltime=0:30:00
#PBS -l storage=scratch/gh0+scratch/xv83
#PBS -l wd
#PBS -o output/PBS/
#PBS -j oe

echo "Going into ACCESS-TMIP"
cd ~/Projects/TMIP/ACCESS-TMIP

experiment=historical
# time_window=Jan1990-Dec1999
time_window=Jan1850-Dec1859
# experiment=ssp370
# time_window=Jan2030-Dec2039
# time_window=Jan2090-Dec2099

member=member_placeholder

# average=monthlymatrices
average=flow

echo "Computing age for $experiment $member $time_window (mean $average)"
julia src/compute_steadyidealage_iterative_varied_members.jl $experiment $member $time_window $average &> output/$PBS_JOBID.compute_steadyidealage_iterative_varied_members.$experiment.$member.$time_window.mean$average.out


