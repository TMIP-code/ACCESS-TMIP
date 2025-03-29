for member in r{1..40}i1p1f1; do
    # change member_placeholders in timestep_adjoint_propagator.sh and create file and remove it
    sed "s/member_placeholder/$member/g" src/timestep_adjoint_propagator_yearly.sh > src/tmp.sh
    qsub src/tmp.sh
    rm src/tmp.sh
done

