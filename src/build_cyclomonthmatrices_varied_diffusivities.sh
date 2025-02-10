#!/bin/bash

#PBS -P xv83
#PBS -N cyclomonthmatrices
#PBS -l ncpus=28
#PBS -l mem=120GB
#PBS -l jobfs=4GB
#PBS -l walltime=3:00:00
#PBS -l storage=scratch/gh0+gdata/xv83+scratch/xv83
#PBS -l wd
#PBS -o output/PBS/
#PBS -j oe

echo "Going into ACCESS-TMIP"
cd ~/Projects/TMIP/ACCESS-TMIP

experiment=historical
time_window=Jan1850-Dec1859
# time_window=Jan1990-Dec1999
# experiment=ssp370
# time_window=Jan2030-Dec2039
# time_window=Jan2090-Dec2099

# for member in r{1..40}i1p1f1; do
for member in r{20..20}i1p1f1; do
    echo "building $experiment $member $time_window"
    julia src/build_cyclomonthmatrices_varied_diffusivities.jl $experiment $member $time_window &> output/build_cyclomonthmatrices_varied_diffusivities.$experiment.$member.$time_window.$PBS_JOBID.out
done


