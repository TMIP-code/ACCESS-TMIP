#!/bin/bash

#PBS -P xv83
#PBS -N adjoint_member_placeholder_finalmonth_placeholder
#PBS -l ncpus=48
#PBS -q hugemem
#PBS -l mem=500GB
#PBS -l jobfs=4GB
#PBS -l walltime=15:00:00
#PBS -l storage=scratch/gh0+scratch/xv83
#PBS -l wd
#PBS -o output/PBS/
#PBS -j oe

# experiment=historical
# time_window=Jan1850-Dec1859
# time_window=Jan1990-Dec1999
experiment=ssp370
# time_window=Jan2030-Dec2039
time_window=Jan2090-Dec2099
member=member_placeholder
finalmonth=finalmonth_placeholder
WRITEDATA=true

echo "Going into ACCESS-TMIP"
cd ~/Projects/TMIP/ACCESS-TMIP

echo "Running script"
julia src/timestep_adjoint_propagator.jl $experiment $member $time_window $finalmonth $WRITEDATA &> output/timestep_adjoint_propagator.$experiment.$member.$time_window.$finalmonth.$PBS_JOBID.out
