#!/bin/bash

#PBS -P y99
#PBS -N age_OM2-1
#PBS -l ncpus=48
#PBS -q normal
#PBS -l mem=190GB
#PBS -l jobfs=4GB
#PBS -l walltime=02:00:00
#PBS -l storage=scratch/gh0+scratch/y99+gdata/xp65
#PBS -l wd
#PBS -o output/PBS/
#PBS -j oe

echo "Going into ACCESS-TMIP"
cd ~/Projects/TMIP/ACCESS-TMIP

export MKL_PARDISO_OOC_PATH = /scratch/y99/TMIP/data/ACCESS-OM2-1/1deg_jra55_iaf_omip2_cycle6/Jan1960-Dec1979/
export MKL_PARDISO_OOC_KEEP_FILE = 1

echo "Running script"
julia src/compute_age_ACCESS-OM2-1.jl &> output/compute_age_ACCESS-OM2-1.$PBS_JOBID.out
