for kVdeep in 2e-5 3e-5 4e-5; do
    for kVML in 0.001 0.01 0.1 1 10; do
        for kH in 50 100 200 300 400 500 600; do
            sed "s/kVdeep_placeholder/$kVdeep/g" src/compute_cycloidealage_iterative_varied_diffusivities.sh > tmp.sh
            sed -i "s/kVML_placeholder/$kVML/g" tmp.sh
            sed -i "s/kH_placeholder/$kH/g" tmp.sh
            qsub tmp.sh
            rm tmp.sh
        done
    done
done