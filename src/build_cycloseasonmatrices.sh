#!/bin/bash

#PBS -P xv83
#PBS -N cycloseasonmatrices
#PBS -l ncpus=48
#PBS -l mem=180GB
#PBS -l jobfs=4GB
#PBS -l walltime=12:00:00
#PBS -l storage=scratch/gh0+gdata/xv83+scratch/xv83
#PBS -l wd
#PBS -o output/PBS/
#PBS -j oe

echo "Going into ACCESS-TMIP"
cd ~/Projects/TMIP/ACCESS-TMIP

echo "Running script"
julia src/build_cycloseasonmatrices.jl &> output/$PBS_JOBID.build_cycloseasonmatrices.out
