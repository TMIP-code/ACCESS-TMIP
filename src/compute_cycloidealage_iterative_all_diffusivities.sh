for kVdeep in 3e-8 1e-7 3e-7 1e-6 3e-6; do
    for kVML in 0 1e-4 3e-4 1e-3 3e-3 1e-2; do
        for kH in 0 1 3 10 30 100; do
            sed "s/kVdeep_placeholder/$kVdeep/g" src/compute_cycloidealage_iterative_varied_diffusivities.sh > tmp.sh
            sed -i "s/kVML_placeholder/$kVML/g" tmp.sh
            sed -i "s/kH_placeholder/$kH/g" tmp.sh
            qsub tmp.sh
            rm tmp.sh
        done
    done
done