# for kVdeep in 1e-7 3e-7 1e-6 3e-6 1e-5 3e-5 1e-4 3e-4; do
for kVML in 0.001 0.01 0.1 1; do
    # for kH in 5 15; do
    # sed "s/kVdeep_placeholder/$kVdeep/g" src/compute_cycloidealage_iterative_varied_diffusivities.sh > src/tmp.sh
    # sed "s/kH_placeholder/$kH/g" src/tmp.sh > src/tmp2.sh
    sed "s/kVML_placeholder/$kVML/g" src/compute_cycloidealage_iterative_varied_diffusivities.sh > src/tmp2.sh
    qsub src/tmp2.sh
    # rm src/tmp.sh
    rm src/tmp2.sh
    # done
done