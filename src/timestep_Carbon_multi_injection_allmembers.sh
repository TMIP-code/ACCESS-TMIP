for member in r{28..40}i1p1f1; do
# for member in r{1..1}i1p1f1; do
    # for srcname in Marlo; do
    for srcname in Karratha Portland Marlo; do
        # change member_placeholders in timestep_Carbon_multi_injection.sh and create file and remove it
        sed "s/member_placeholder/$member/g" src/timestep_Carbon_multi_injection.sh > src/tmp.sh
        sed "s/srcname_placeholder/$srcname/g" src/tmp.sh > src/tmp2.sh
        qsub src/tmp2.sh
        rm src/tmp.sh
        rm src/tmp2.sh
    done
done