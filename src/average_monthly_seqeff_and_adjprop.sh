#!/bin/bash

#PBS -P xv83
#PBS -N avg_member_placeholder_time_window_placeholder
#PBS -l ncpus=16
#PBS -q express
#PBS -l mem=63GB
#PBS -l jobfs=4GB
#PBS -l walltime=0:10:00
#PBS -l storage=scratch/gh0+scratch/xv83
#PBS -l wd
#PBS -o output/PBS/
#PBS -j oe

experiment=ssp370
time_window=time_window_placeholder
member=member_placeholder

echo "Going into ACCESS-TMIP"
cd ~/Projects/TMIP/ACCESS-TMIP

echo "Running script"
julia src/average_monthly_seqeff_and_adjprop.jl $experiment $member $time_window &> output/average_monthly_seqeff_and_adjprop.$experiment.$member.$time_window.$PBS_JOBID.out
