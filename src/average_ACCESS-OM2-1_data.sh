# command used to get the job resources
# qsub -I -P y99 -l mem=47GB -l storage=scratch/gh0+scratch/xv83+scratch/y99+gdata/cj50+gdata/ik11 -l walltime=01:00:00 -l ncpus=12

# I have built matrices for 1° simulations before,
# from ACCESS1.0, ACCESS1.3, ACCESSACCESS-ESM1.5, and ACCESS-CM2
# but not from ACCESS-OM2 at 1°, let alone finer resolutions 0.25° and 0.1°.
# So let's try that here.

# For both 1° and 0.25°, I need the transport output from parameterized
# transport terms. In CMIP it should be all in umo and vmo, but
# in ACCESS output it's split into tx_trans, tx_trans_gm,
# and tx_trans_submeso (and same for ty_trans). The parameterized terms
# are also summed vertically, so they need to be diffed vertically
# to get back the corresponding transport on the cell face.
# No need for gm terms for the 0.1° simulations, but need submeso terms? # CHECK
#
# I think i need the following variables:
# - tx_trans
# - ty_trans
# - tx_trans_submeso
# - ty_trans_submeso
# - area_t
# - dzt
# and these for 1° and 0.25°
# - mld
# - tx_trans_gm
# - ty_trans_gm
#
# But unfortunately it looks like only these runs have it:
# - 01deg_jra55v13_ryf9091
# - 01deg_jra55v13_ryf9091_easterlies_down10
# - 01deg_jra55v13_ryf9091_easterlies_up10
# - 01deg_jra55v13_ryf9091_easterlies_up10_meridional
# - 01deg_jra55v13_ryf9091_easterlies_up10_zonal
# - 01deg_jra55v13_ryf9091_weddell_down2
# - 01deg_jra55v13_ryf9091_weddell_up1
# - 1deg_jra55_ryf9091_gadi/output031

# So as a test, let's try 1deg_jra55_ryf9091_gadi

module load netcdf
module load nco



# inputdir=/g/data/cj50/access-om2/raw-output/access-om2/1deg_jra55_ryf9091_gadi
inputdir=/g/data/ik11/outputs/access-om2/1deg_jra55_ryf9091_gadi
outputdir=/scratch/y99/bp3051/hiresTM/access-om2-1/1deg_jra55_ryf9091_gadi
mkdir -p ${outputdir}

ls ${inputdir}/output0[0-3]?/ocean/ocean.nc

# ocean.nc is yearly data
# ocean_month.nc is monthly data
# 10 years per file at 1° (10 yearly steps = 120 monthly steps)

ncdump -h ${inputdir}/output000/ocean/ocean.nc | grep 't[xy]_trans.*(.*)'
# shows that 1° yearly output (ocean.nc) includes
# - tx_trans
# - ty_trans
# - tx_trans_rho
# - ty_trans_rho
# - ty_trans_rho_gm # WHY?

ncdump -h ${inputdir}/output000/ocean/ocean_month.nc | grep 't[xy]_trans.*(.*)'
# shows that 1° monthly output (ocean_month.nc) includes only
# - tx_trans_int_z
# - ty_trans_int_z
# So I will probably need to diff these along z if I use monthly data # CHECK

# just take the mean to start with
# do it for both tx_trans and ty_trans and parameterized terms
for v in tx_trans ty_trans tx_trans_submeso ty_trans_submeso area_t dzt mld tx_trans_gm ty_trans_gm; do
    echo ${v}
    ncra -v ${v} -O  ${inputdir}/output031/ocean/ocean.nc ${outputdir}/${v}.nc
    echo ${v}_grid
    ncks -v ${v} -O ${inputdir}/output031/ocean/ocean_grid.nc ${outputdir}/${v}_grid.nc
    echo ${v}_month
    ncra -v ${v} -O ${inputdir}/output031/ocean/ocean_month.nc ${outputdir}/${v}_month.nc
done

ls -lh ${outputdir}/

# ncdump -h /scratch/y99/bp3051/hiresTM/access-om2-1/tx_trans_nco_test.nc
# ls -lh /scratch/y99/bp3051/hiresTM/access-om2-1/

# exit
