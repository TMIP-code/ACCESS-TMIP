# for member in r{1..40}i1p1f1; do
for member in r{1..1}i1p1f1; do
    # for srcname in Marlo; do
    for finalmonth in {1..12}; do
        # change member_placeholders in timestep_adjoint_propagator.sh and create file and remove it
        sed "s/member_placeholder/$member/g" src/timestep_adjoint_propagator.sh > src/tmp.sh
        sed "s/finalmonth_placeholder/$finalmonth/g" src/tmp.sh > src/tmp2.sh
        qsub src/tmp2.sh
        rm src/tmp.sh
        rm src/tmp2.sh
    done
done
