for member in r{2..40}i1p1f1; do
    sed "s/member_placeholder/$member/g" src/compute_steadyidealage_iterative_varied_members.sh > src/tmp.sh
    qsub src/tmp.sh
    rm src/tmp.sh
done