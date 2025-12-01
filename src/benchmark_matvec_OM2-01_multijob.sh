for NCPUS in 1 2 3 4 6 8 12 16 20 24 36 48 ; do
# for NCPUS in 1 2 3 4 6 8 12 16 20 24 36 48 ; do
# for NCPUS in 1 2 ; do
    for sparseBLAS in Julia MKL CSR-tmul CSR-bmul; do
    # for sparseBLAS in Julia MKL; do
        # change member_placeholders in timestep_adjoint_propagator.sh and create file and remove it
        sed "s/ncpusplaceholder/$NCPUS/g" src/benchmark_matvec_OM2-01.sh > src/tmp.sh
        sed -i "s/sparseBLASplaceholder/$sparseBLAS/g" src/tmp.sh
        qsub src/tmp.sh
        rm src/tmp.sh
    done
done

