#!/bin/bash

#PBS -P xv83
#PBS -N Gamma_maxMLD
#PBS -l ncpus=28
#PBS -l mem=180GB
#PBS -l jobfs=4GB
#PBS -l walltime=12:00:00
#PBS -l storage=scratch/p66+gdata/p73
#PBS -l wd
#PBS -o output/PBS/
#PBS -j oe

echo "Going into ACCESS-TMIP"
cd ~/Projects/TMIP/ACCESS-TMIP

echo "Running script"
julia src/compute_Gammas_GMcomparison_maxMLD.jl &> output/$PBS_JOBID.compute_Gammas_GMcomparison_maxMLD.out
