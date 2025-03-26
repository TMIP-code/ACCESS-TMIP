#!/bin/bash

#PBS -P xv83
#PBS -N age_kVdeep_placeholder_kVML_placeholder_kH_placeholder
#PBS -q express
#PBS -l ncpus=24
#PBS -l mem=100GB
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

kVdeep=kVdeep_placeholder
kVML=kVML_placeholder
kH=kH_placeholder

# average=monthlymatrices
average=flow

# for member in r{20..20}i1p1f1; do
for member in AA; do
    echo "Computing age for $experiment $member $time_window (mean $average) kVdeep=$kVdeep kVML=$kVML kH=$kH"
    julia src/compute_steadyidealage_iterative_varied_diffusivities.jl $experiment $member $time_window $average $kVdeep $kVML $kH &> output/$PBS_JOBID.compute_steadyidealage_iterative_varied_members.$experiment.$member.$time_window.mean$average.$kVdeep.$kVML.$kH.out
done


