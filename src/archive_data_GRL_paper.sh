#!/bin/bash
#PBS -P xv83
#PBS -l ncpus=1
#PBS -l mem=30GB
#PBS -q copyq
#PBS -l walltime=03:00:00
#PBS -l wd
#PBS -l storage=scratch/xv83+gdata/xv83
#PBS -o output/PBS/
#PBS -j oe


SOURCEDIR='/scratch/xv83/TMIP/data/ACCESS-ESM1-5/'
DESTDIR='/g/data/xv83/TMIP/data/ACCESS-ESM1-5/'

echo "Archiving monthly climatologies"
rsync -a -m --include='ssp370/r*i1p1f1/Jan20*-Dec20*/cyclomonth/t*_trans_s*.nc' --include='*/' --exclude='*' ${SOURCEDIR} ${DESTDIR}

echo "Archiving monthly matrices"
rsync -a -m --include='ssp370/r*i1p1f1/Jan20*-Dec20*/cyclomonth/cyclo_matrix_centered_kVdeep3e-05_kH300_kVML1e+00_*.jld2' --include='*/' --exclude='*' ${SOURCEDIR} ${DESTDIR}

echo "Archiving mean reemergence time"
rsync -a -m --include='ssp370/r*i1p1f1/Jan20*-Dec20*/cyclomonth/reemergence_time_centered_kVdeep3e-05_kH300_kVML1e+00.nc' --include='*/' --exclude='*' ${SOURCEDIR} ${DESTDIR}

echo "Archiving mean TTD"
rsync -a -m --include='ssp370/r*i1p1f1/Jan20*-Dec20*/calgtilde_centered_kVdeep3e-05_kH300_kVML1e+00.nc' --include='*/' --exclude='*' ${SOURCEDIR} ${DESTDIR}

echo "Archiving mean sequestration efficiency"
rsync -a -m --include='ssp370/r*i1p1f1/Jan20*-Dec20*/calE_centered_kVdeep3e-05_kH300_kVML1e+00.nc' --include='*/' --exclude='*' ${SOURCEDIR} ${DESTDIR}

echo "Done archiving!"