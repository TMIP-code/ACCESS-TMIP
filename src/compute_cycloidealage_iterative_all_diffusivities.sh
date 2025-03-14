for kVdeep in 3e-8 1e-7 3e-7 1e-6 3e-6; do
    for kVML in 0 1e-4 3e-4 1e-3 3e-3 1e-2; do
        for kH in 0 1 3 10 30 100; do
            sed "s/kVdeep_placeholder/$kVdeep/g" src/compute_cycloidealage_iterative_varied_diffusivities.sh > src/tmp.sh
            sed "s/kVML_placeholder/$kVML/g" src/tmp.sh > src/tmp2.sh
            sed "s/kH_placeholder/$kH/g" src/tmp2.sh > src/tmp3.sh
            qsub src/tmp3.sh
            rm src/tmp.sh
            rm src/tmp2.sh
            rm src/tmp3.sh
        done
    done
done