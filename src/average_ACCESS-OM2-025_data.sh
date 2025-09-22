# command used to get the job resources
# qsub -I -P y99 -l mem=47GB -l storage=scratch/gh0+scratch/xv83+scratch/y99+gdata/cj50 -l walltime=01:00:00 -l ncpus=12

module load netcdf
module load nco

inputdir=/g/data/cj50/access-om2/raw-output/access-om2-025/025deg_jra55v13_ryf9091_gmredi6
outputdir=/scratch/y99/bp3051/hiresTM/access-om2-025
mkdir -p ${outputdir}

ls ${inputdir}/output0[0-3]?/ocean/ocean.nc

# ocean.nc is yearly data
# ocean_month.nc is monthly data
# 2 years per file at 0.25° (2 yearly steps = 24 monthly steps)

ncdump -h ${inputdir}/output000/ocean/ocean.nc | grep 't[xy]_trans.*(.*)'
# shows that 0.25° yearly output (ocean.nc) includes
# - tx_trans
# - ty_trans
# - ty_trans_submeso # WHY?
# - tx_trans_rho
# - ty_trans_rho
# - ty_trans_rho_gm # WHY?
# - ty_trans_nrho_submeso # WHY?

ncdump -h ${inputdir}/output000/ocean/ocean_month.nc | grep 't[xy]_trans.*(.*)'
# shows that 0.25° monthly output (ocean_month.nc) includes only
# - tx_trans_int_z
# - ty_trans_int_z
# So I will probably need to diff these along z if I use monthly data # CHECK

# just take the mean to start with
# do it for both tx_trans and ty_trans and parameterized terms
for v in tx_trans ty_trans tx_trans_gm ty_trans_gm tx_trans_submeso ty_trans_submeso; do
    time ncra -v ${v} -O ${inputdir}/output0[0-3]?/ocean/ocean.nc ${outputdir}/${v}_nco_test.nc
done

ncdump -h /scratch/y99/bp3051/hiresTM/access-om2-025/tx_trans_nco_test.nc
ls -lh /scratch/y99/bp3051/hiresTM/access-om2-025/

exit
