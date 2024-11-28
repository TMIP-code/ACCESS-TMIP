#!/bin/bash

#PBS -P xv83
#PBS -N cyclomonthmatrices
#PBS -l ncpus=48
#PBS -l mem=180GB
#PBS -l jobfs=4GB
#PBS -l walltime=1:00:00
#PBS -l storage=gdata/gh0+gdata/xv83+scratch/xv83
#PBS -l wd
#PBS -o output/PBS/
#PBS -j oe

echo "Going into ACCESS-TMIP"
cd ~/Projects/TMIP/ACCESS-TMIP

experiment=historical
time_window=Jan1990-Dec1999

for member in r{26..40}i1p1f1; do
    echo "building $experiment $member $time_window"
    julia src/build_cyclomonthmatrices.jl $experiment $member $time_window &> output/$PBS_JOBID.build_cyclomonthmatrices.$experiment.$member.$time_window.out
done


