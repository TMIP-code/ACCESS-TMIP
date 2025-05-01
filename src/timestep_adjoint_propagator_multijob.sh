# for member in r{1..40}i1p1f1; do
for member in r{40..40}i1p1f1; do
    for finalmonth in 3 4 5 6; do
        # change member_placeholders in timestep_adjoint_propagator.sh and create file and remove it
        sed "s/member_placeholder/$member/g" src/timestep_adjoint_propagator.sh > src/tmp.sh
        sed -i "s/finalmonth_placeholder/$finalmonth/g" src/tmp.sh
        qsub src/tmp.sh
        rm src/tmp.sh
    done
done

