#!/bin/bash

#PBS -P xv83
#PBS -N Gamma_monthlytransport
#PBS -l ncpus=28
#PBS -l mem=180GB
#PBS -l jobfs=4GB
#PBS -l walltime=12:00:00
#PBS -l storage=scratch/p66+gdata/p73
#PBS -l storage=gdata/fs38+gdata/xv83+scratch/xv83+gdata/dk92+gdata/hh5+gdata/xp65+gdata/p73+scratch/p66
#PBS -l wd
#PBS -o output/PBS/
#PBS -j oe

echo "Loading conda/analysis3-24.04 module"
module use /g/data/hh5/public/modules
module load conda/analysis3-24.04
conda activate conda/analysis3-24.04
conda info

echo "Loading python3/3.12.1"
module load python3/3.12.1

echo "Going into ACCESS-TMIP"
cd ~/Projects/TMIP/ACCESS-TMIP

echo "Running script"
julia src/compute_Gammas_frommonthlymatrices.jl &> output/$PBS_JOBID.compute_Gammas_frommonthlymatrices.out
