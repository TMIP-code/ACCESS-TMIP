#!/bin/bash

#PBS -P y99
#PBS -N age_OM2-025_PETSc
#PBS -l ncpus=48
#PBS -q normal
#PBS -l mem=190GB
#PBS -l jobfs=4GB
#PBS -l walltime=24:00:00
#PBS -l storage=scratch/gh0+scratch/y99+scratch/p66
#PBS -l wd
#PBS -o output/PBS/
#PBS -j oe



echo "Going into ACCESS-TMIP"
cd ~/Projects/TMIP/ACCESS-TMIP

echo "Running script"
julia src/compute_age_ACCESS-OM2-025_PETSc.jl &> output/$PBS_JOBID.compute_age_ACCESS-OM2-025_PETSc.out


everything = [
    "01deg_jra55_ryf_Control"
    "01deg_jra55_ryf_ENFull"
    "01deg_jra55_ryf_LNFull"
    "01deg_jra55v13_ryf9091"
    "01deg_jra55v13_ryf9091_easterlies_down10"
    "01deg_jra55v13_ryf9091_easterlies_up10"
    "01deg_jra55v13_ryf9091_easterlies_up10_meridional"
    "01deg_jra55v13_ryf9091_easterlies_up10_zonal"
    "01deg_jra55v13_ryf9091_qian_wthmp"
    "01deg_jra55v13_ryf9091_qian_wthp"
    "01deg_jra55v13_ryf9091_weddell_down2"
    "01deg_jra55v13_ryf9091_weddell_up1"
    "01deg_jra55v140_iaf"
    "01deg_jra55v140_iaf_cycle2"
    "01deg_jra55v140_iaf_cycle3"
    "01deg_jra55v140_iaf_cycle4"
    "01deg_jra55v140_iaf_cycle4_jra55v150_extension"
    "01deg_jra55v150_iaf_cycle1"
    "025deg_jra55_iaf_omip2_cycle1"
    "025deg_jra55_iaf_omip2_cycle2"
    "025deg_jra55_iaf_omip2_cycle3"
    "025deg_jra55_iaf_omip2_cycle4"
    "025deg_jra55_iaf_omip2_cycle5"
    "025deg_jra55_iaf_omip2_cycle6"
    "1deg_jra55_iaf_omip2_cycle1"
    "1deg_jra55_iaf_omip2_cycle2"
    "1deg_jra55_iaf_omip2_cycle3"
    "1deg_jra55_iaf_omip2_cycle4"
    "1deg_jra55_iaf_omip2_cycle5"
    "1deg_jra55_iaf_omip2_cycle6"
    "1deg_jra55_iaf_omip2spunup_cycle1"
    "1deg_jra55_iaf_omip2spunup_cycle2"
    "1deg_jra55_iaf_omip2spunup_cycle3"
    "1deg_jra55_iaf_omip2spunup_cycle34"
    "1deg_jra55_iaf_omip2spunup_cycle35"
    "1deg_jra55_iaf_omip2spunup_cycle36"
    "1deg_jra55_iaf_omip2spunup_cycle37"
    "1deg_jra55_iaf_omip2spunup_cycle38"
    "1deg_jra55_iaf_omip2spunup_cycle39"
    "1deg_jra55_iaf_omip2spunup_cycle4"
    "1deg_jra55_iaf_omip2spunup_cycle5"
    "1deg_jra55_iaf_omip2spunup_cycle6"
    "1deg_jra55_ryf9091_gadi"
]


- `01deg_jra55v13_ryf9091`
- `01deg_jra55v13_ryf9091_easterlies_down10`
- `01deg_jra55v13_ryf9091_easterlies_up10`
- `01deg_jra55v13_ryf9091_easterlies_up10_meridional`
- `01deg_jra55v13_ryf9091_easterlies_up10_zonal`
- `01deg_jra55v13_ryf9091_weddell_down2`
- `01deg_jra55v13_ryf9091_weddell_up1`
- `025deg_jra55_iaf_omip2_cycle1`
- `025deg_jra55_iaf_omip2_cycle2`
- `025deg_jra55_iaf_omip2_cycle3`
- `025deg_jra55_iaf_omip2_cycle4`
- `025deg_jra55_iaf_omip2_cycle5`
- `025deg_jra55_iaf_omip2_cycle6`
- `1deg_jra55_iaf_omip2_cycle1`
- `1deg_jra55_iaf_omip2_cycle2`
- `1deg_jra55_iaf_omip2_cycle3`
- `1deg_jra55_iaf_omip2_cycle4`
- `1deg_jra55_iaf_omip2_cycle5`
- `1deg_jra55_iaf_omip2_cycle6`
- `1deg_jra55_iaf_omip2spunup_cycle1`
- `1deg_jra55_iaf_omip2spunup_cycle2`
- `1deg_jra55_iaf_omip2spunup_cycle3`
- `1deg_jra55_iaf_omip2spunup_cycle34`
- `1deg_jra55_iaf_omip2spunup_cycle35`
- `1deg_jra55_iaf_omip2spunup_cycle36`
- `1deg_jra55_iaf_omip2spunup_cycle37`
- `1deg_jra55_iaf_omip2spunup_cycle38`
- `1deg_jra55_iaf_omip2spunup_cycle39`
- `1deg_jra55_iaf_omip2spunup_cycle4`
- `1deg_jra55_iaf_omip2spunup_cycle5`
- `1deg_jra55_iaf_omip2spunup_cycle6`
- `1deg_jra55_ryf9091_gadi`


no GM
21         025deg_jra55_iaf_omip2_cycle1
22         025deg_jra55_iaf_omip2_cycle2
23         025deg_jra55_iaf_omip2_cycle3
24         025deg_jra55_iaf_omip2_cycle4
25         025deg_jra55_iaf_omip2_cycle5
26         025deg_jra55_iaf_omip2_cycle6
32           1deg_jra55_iaf_omip2_cycle1
33           1deg_jra55_iaf_omip2_cycle2
34           1deg_jra55_iaf_omip2_cycle3
35           1deg_jra55_iaf_omip2_cycle4
36           1deg_jra55_iaf_omip2_cycle5
37           1deg_jra55_iaf_omip2_cycle6
38     1deg_jra55_iaf_omip2spunup_cycle1
49     1deg_jra55_iaf_omip2spunup_cycle2
60     1deg_jra55_iaf_omip2spunup_cycle3
65    1deg_jra55_iaf_omip2spunup_cycle34
66    1deg_jra55_iaf_omip2spunup_cycle35
67    1deg_jra55_iaf_omip2spunup_cycle36
68    1deg_jra55_iaf_omip2spunup_cycle37
69    1deg_jra55_iaf_omip2spunup_cycle38
70    1deg_jra55_iaf_omip2spunup_cycle39
71     1deg_jra55_iaf_omip2spunup_cycle4
78     1deg_jra55_iaf_omip2spunup_cycle5
79     1deg_jra55_iaf_omip2spunup_cycle6



no GM no submeso

21         025deg_jra55_iaf_omip2_cycle1
22         025deg_jra55_iaf_omip2_cycle2
23         025deg_jra55_iaf_omip2_cycle3
24         025deg_jra55_iaf_omip2_cycle4
25         025deg_jra55_iaf_omip2_cycle5
26         025deg_jra55_iaf_omip2_cycle6
32           1deg_jra55_iaf_omip2_cycle1
33           1deg_jra55_iaf_omip2_cycle2
34           1deg_jra55_iaf_omip2_cycle3
35           1deg_jra55_iaf_omip2_cycle4
36           1deg_jra55_iaf_omip2_cycle5
37           1deg_jra55_iaf_omip2_cycle6
38     1deg_jra55_iaf_omip2spunup_cycle1
49     1deg_jra55_iaf_omip2spunup_cycle2
60     1deg_jra55_iaf_omip2spunup_cycle3
65    1deg_jra55_iaf_omip2spunup_cycle34
66    1deg_jra55_iaf_omip2spunup_cycle35
67    1deg_jra55_iaf_omip2spunup_cycle36
68    1deg_jra55_iaf_omip2spunup_cycle37
69    1deg_jra55_iaf_omip2spunup_cycle38
70    1deg_jra55_iaf_omip2spunup_cycle39
71     1deg_jra55_iaf_omip2spunup_cycle4
78     1deg_jra55_iaf_omip2spunup_cycle5
79     1deg_jra55_iaf_omip2spunup_cycle6
0                               01deg_jra55_ryf_Control
1                                01deg_jra55_ryf_ENFull
2                                01deg_jra55_ryf_LNFull
3                                01deg_jra55v13_ryf9091
4              01deg_jra55v13_ryf9091_easterlies_down10
5                01deg_jra55v13_ryf9091_easterlies_up10
6     01deg_jra55v13_ryf9091_easterlies_up10_meridional
7          01deg_jra55v13_ryf9091_easterlies_up10_zonal
8                     01deg_jra55v13_ryf9091_qian_wthmp
9                      01deg_jra55v13_ryf9091_qian_wthp
10                 01deg_jra55v13_ryf9091_weddell_down2
11                   01deg_jra55v13_ryf9091_weddell_up1
12                                  01deg_jra55v140_iaf
13                           01deg_jra55v140_iaf_cycle2
14                           01deg_jra55v140_iaf_cycle3
15                           01deg_jra55v140_iaf_cycle4
16       01deg_jra55v140_iaf_cycle4_jra55v150_extension
17                           01deg_jra55v150_iaf_cycle1
27                            025deg_jra55_ryf9091_gadi
29                                        1deg_era5_iaf
31                        1deg_jra55_iaf_era5comparison
83                              1deg_jra55_ryf9091_gadi