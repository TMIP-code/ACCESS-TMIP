#!/bin/bash

#PBS -P y99
#PBS -N yrlyTM_OM2-01-qian
#PBS -q hugemem
#PBS -l ncpus=16
#PBS -l mem=490GB
#PBS -l jobfs=4GB
#PBS -l walltime=3:00:00
#PBS -l storage=scratch/gh0+gdata/xv83+scratch/xv83+scratch/p66+scratch/y99+gdata/cj50+gdata/ik11
#PBS -l wd
#PBS -o output/PBS/
#PBS -j oe

# qsub -I -P y99 -N yrlyTM_OM2-01-qian -q hugemem -l ncpus=24 -l mem=735GB -l jobfs=4GB -l walltime=3:00:00 -l storage=scratch/gh0+gdata/xv83+scratch/xv83+scratch/p66+scratch/y99+gdata/cj50+gdata/ik11 -l wd -o output/PBS/ -j oe


echo "Going into ACCESS-TMIP"
cd ~/Projects/TMIP/ACCESS-TMIP

julia src/build_yearlyTM_ACCESS-OM2-01_qian.jl &> output/build_yearlyTM_ACCESS-OM2-01_qian.$PBS_JOBID.out


