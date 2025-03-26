# for kVdeep in 1e-8 1e-7 1e-6 1e-5 1e-4; do
#     for kVML in 1e-8 1e-6 1e-4 1e-2 1e-0; do
#         for kH in 1 10 100 1000; do
# for kVdeep in 1e-5 2e-5 3e-5 4e-5 7e-5 1e-4; do
#     for kVML in 1e-8 1e-6 1e-4 1e-2 1e-0; do
#         for kH in 50 100 200 500 1000 2000; do
for kVdeep in 2e-5 3e-5 4e-5; do
    for kVML in 1e-8 1e-6 1e-4 1e-2 1e-0; do
        for kH in 50 100 200 500 1000 2000; do
            sed "s/kVdeep_placeholder/$kVdeep/g" src/compute_steadyidealage_iterative_varied_diffusivities.sh > tmp.sh
            sed -i "s/kVML_placeholder/$kVML/g" tmp.sh
            sed -i "s/kH_placeholder/$kH/g" tmp.sh
            qsub tmp.sh
            rm tmp.sh
        done
    done
done