# command used to get the job resources
# qsub -I -P y99 -l mem=47GB -l storage=scratch/gh0+scratch/xv83+scratch/y99+gdata/cj50 -l walltime=01:00:00 -l ncpus=12

module load netcdf
module load nco

inputdir=/g/data/cj50/access-om2/raw-output/access-om2-01/01deg_jra55v13_ryf9091
outputdir=/scratch/y99/bp3051/hiresTM

# just take the mean to start with
# do it for both tx_trans and ty_trans

# ocean.nc is yearly data
# ocean_month.nc is monthly data
# 1/4 year (= 3 months, but averaged into a single time step) per ocean.nc file at 0.1°
# 3 months per ocean_month.nc file at 0.1°

ncdump -h ${inputdir}/output000/ocean/ocean.nc | grep 't[xy]_trans.*(.*)'
# shows that 0.1° yearly output (ocean.nc) includes
# - tx_trans
# - ty_trans
# - ty_trans_submeso # WHY?
# - tx_trans_rho
# - ty_trans_rho
# - ty_trans_nrho_submeso # WHY?

ncdump -h ${inputdir}/output000/ocean/ocean_month.nc | grep 't[xy]_trans.*(.*)'
# shows that 0.1° monthly output (ocean_month.nc) includes only
# - tx_trans_int_z
# - ty_trans_int_z
# So I will probably need to diff these along z if I use monthly data # CHECK

# So let's start taking a decadal mean of tx_trans and ty_trans
# Since
for v in tx_trans ty_trans; do
    time ncra -v ${v} -O ${inputdir}/output0[0-3]?/ocean/ocean.nc ${outputdir}/${v}_nco_test.nc
done

ncdump -h /scratch/y99/bp3051/hiresTM/tx_trans_nco_test.nc
ls -lh /scratch/y99/bp3051/hiresTM/

# exit

# /g/data/cj50/access-om2/raw-output/access-om2-01/01deg_jra55v13_ryf9091/output000/ocean/ocean.nc
# /g/data/cj50/access-om2/raw-output/access-om2-01/01deg_jra55v13_ryf9091/output000/ocean/ocean.nc