#!/bin/bash

#PBS -P xv83
#PBS -N cycloage
#PBS -l ncpus=48
#PBS -q hugemem
#PBS -l mem=360GB
#PBS -l jobfs=4GB
#PBS -l walltime=1:00:00
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

for member in r{20..20}i1p1f1; do
    echo "Computing age for $experiment $member $time_window"
    julia src/compute_cycloidealage_iterative.jl $experiment $member $time_window &> output/$PBS_JOBID.compute_cycloidealage_iterative..$experiment.$member.$time_window.out
done


