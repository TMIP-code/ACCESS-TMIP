#!/bin/bash

#PBS -P xv83
#PBS -N ACCESS-ESM1-5_diag
#PBS -l ncpus=28
#PBS -l mem=180GB
#PBS -l jobfs=4GB
#PBS -l walltime=24:00:00
#PBS -l storage=gdata/xv83
#PBS -l wd
#PBS -o output/PBS/
#PBS -j oe

echo "Going into ACCESS-TMIP"
cd ~/Projects/TMIP/ACCESS-TMIP

echo "Running script"
julia src/compute_ideal_ages.jl
