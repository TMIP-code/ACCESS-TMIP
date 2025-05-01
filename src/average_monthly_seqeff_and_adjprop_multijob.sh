for member in r{14..14}i1p1f1; do
    for time_window in Jan2030-Dec2039 Jan2090-Dec2099; do
        # change member_placeholders in timestep_adjoint_propagator.sh and create file and remove it
        sed "s/member_placeholder/$member/g" src/average_monthly_seqeff_and_adjprop.sh > src/tmp.sh
        sed -i "s/time_window_placeholder/$time_window/g" src/tmp.sh
        qsub src/tmp.sh
        rm src/tmp.sh
    done
done

