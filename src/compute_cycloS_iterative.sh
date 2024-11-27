#!/bin/bash

#PBS -P xv83
#PBS -N cycloS
#PBS -l ncpus=48
#PBS -q hugemem
#PBS -l mem=360GB
#PBS -l jobfs=4GB
#PBS -l walltime=2:00:00
#PBS -l storage=gdata/gh0+scratch/xv83
#PBS -l wd
#PBS -o output/PBS/
#PBS -j oe

echo "Going into ACCESS-TMIP"
cd ~/Projects/TMIP/ACCESS-TMIP

echo "Running script"
julia src/compute_cycloS_iterative.jl &> output/$PBS_JOBID.$model.compute_cycloS_iterative.out
