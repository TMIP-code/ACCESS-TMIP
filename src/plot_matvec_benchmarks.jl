using GLMakie
using DataFrames
using MakieExtra

df = DataFrame()
[
    push!(df, row) for row in
        [
            # (transpose = false, NCPUs = 1, BLAS = "Julia", meantime = 4.662299106363636e8, mediantime = 4.56166099e8, mintime = 4.55542125e8, maxtime = 4.95921563e8, stdtime = 1.65606501166579e7)
            # (transpose = true, NCPUs = 1, BLAS = "Julia", meantime = 4.382445339166667e8, mediantime = 4.290010905e8, mintime = 4.28296656e8, maxtime = 4.65980272e8, stdtime = 1.6341199280377736e7)
            # (transpose = false, NCPUs = 1, BLAS = "MKL", meantime = 6.013181693333334e8, mediantime = 5.89961543e8, mintime = 5.88482177e8, maxtime = 6.2797813e8, stdtime = 1.81569644223452e7)
            # (transpose = true, NCPUs = 1, BLAS = "MKL", meantime = 7.336929502857143e8, mediantime = 7.122526e8, mintime = 7.10131867e8, maxtime = 7.64447654e8, stdtime = 2.788763778033029e7)
            # (transpose = false, NCPUs = 2, BLAS = "Julia", meantime = 4.519010435833333e8, mediantime = 4.42098666e8, mintime = 4.41733024e8, maxtime = 4.81810373e8, stdtime = 1.7717541346496236e7)
            # (transpose = true, NCPUs = 2, BLAS = "Julia", meantime = 4.2975022775e8, mediantime = 4.23078797e8, mintime = 4.22742384e8, maxtime = 4.63263268e8, stdtime = 1.5249999590234056e7)
            # (transpose = false, NCPUs = 2, BLAS = "MKL", meantime = 4.360637025e8, mediantime = 4.3391359e8, mintime = 4.19969825e8, maxtime = 4.55956046e8, stdtime = 1.5719693465315193e7)
            # (transpose = true, NCPUs = 2, BLAS = "MKL", meantime = 3.6039191621428573e8, mediantime = 3.570793145e8, mintime = 3.55594178e8, maxtime = 3.75740453e8, stdtime = 6.912871292603298e6)
            # (transpose = false, NCPUs = 3, BLAS = "Julia", meantime = 4.646093421818182e8, mediantime = 4.56531746e8, mintime = 4.56146862e8, maxtime = 4.95530269e8, stdtime = 1.551938251921635e7)
            # (transpose = true, NCPUs = 3, BLAS = "Julia", meantime = 4.351914885e8, mediantime = 4.286962945e8, mintime = 4.28438762e8, maxtime = 4.67393303e8, stdtime = 1.365574520751625e7)
            # (transpose = false, NCPUs = 3, BLAS = "MKL", meantime = 4.175496208333333e8, mediantime = 4.12719332e8, mintime = 4.10695773e8, maxtime = 4.30614317e8, stdtime = 8.702372882689063e6)
            # (transpose = true, NCPUs = 3, BLAS = "MKL", meantime = 2.5177187457142857e8, mediantime = 2.39833517e8, mintime = 2.39258666e8, maxtime = 2.93626763e8, stdtime = 2.10352000665863e7)
            # (transpose = false, NCPUs = 4, BLAS = "Julia", meantime = 4.497061046666667e8, mediantime = 4.464430025e8, mintime = 4.45991151e8, maxtime = 4.7517058e8, stdtime = 8.545211506795364e6)
            # (transpose = true, NCPUs = 4, BLAS = "Julia", meantime = 4.3389692925e8, mediantime = 4.24121931e8, mintime = 4.23982728e8, maxtime = 4.6355172e8, stdtime = 1.7501447826466262e7)
            # (transpose = false, NCPUs = 4, BLAS = "MKL", meantime = 3.9262944138461536e8, mediantime = 3.82069897e8, mintime = 3.79016795e8, maxtime = 4.11211244e8, stdtime = 1.4346107156566316e7)
            # (transpose = true, NCPUs = 4, BLAS = "MKL", meantime = 1.8652908214814815e8, mediantime = 1.80681527e8, mintime = 1.80252515e8, maxtime = 2.33217621e8, stdtime = 1.3003649225391341e7)
            # (transpose = false, NCPUs = 6, BLAS = "Julia", meantime = 4.607998128181818e8, mediantime = 4.49364234e8, mintime = 4.48951491e8, maxtime = 4.91564664e8, stdtime = 1.780873712000416e7)
            # (transpose = true, NCPUs = 6, BLAS = "Julia", meantime = 4.421662735e8, mediantime = 4.38996131e8, mintime = 4.27338717e8, maxtime = 4.77700362e8, stdtime = 1.6183889257554242e7)
            # (transpose = false, NCPUs = 6, BLAS = "MKL", meantime = 3.893798229230769e8, mediantime = 3.59108014e8, mintime = 3.58922598e8, maxtime = 5.6550377e8, stdtime = 7.396704739305362e7)
            # (transpose = true, NCPUs = 6, BLAS = "MKL", meantime = 1.2993595838461539e8, mediantime = 1.20178876e8, mintime = 1.19797196e8, maxtime = 3.36842085e8, stdtime = 3.746025588050141e7)
            # (transpose = false, NCPUs = 8, BLAS = "Julia", meantime = 4.486334979166667e8, mediantime = 4.45288893e8, mintime = 4.45003086e8, maxtime = 4.82943841e8, stdtime = 1.0821020649974244e7)
            # (transpose = true, NCPUs = 8, BLAS = "Julia", meantime = 4.36907548e8, mediantime = 4.269023675e8, mintime = 4.26547847e8, maxtime = 4.67518358e8, stdtime = 1.824483044735164e7)
            # (transpose = false, NCPUs = 8, BLAS = "MKL", meantime = 7.593063531428572e8, mediantime = 6.59573499e8, mintime = 6.01654061e8, maxtime = 1.141431649e9, stdtime = 2.0265096231295466e8)
            # (transpose = true, NCPUs = 8, BLAS = "MKL", meantime = 2.1067869054166666e8, mediantime = 1.81774866e8, mintime = 1.78098619e8, maxtime = 5.25841194e8, stdtime = 7.507262968650267e7)
            # (transpose = false, NCPUs = 12, BLAS = "Julia", meantime = 4.6273400345454544e8, mediantime = 4.56189498e8, mintime = 4.55820918e8, maxtime = 4.93954559e8, stdtime = 1.411300164251971e7)
            # (transpose = true, NCPUs = 12, BLAS = "Julia", meantime = 4.558107316363636e8, mediantime = 4.52085216e8, mintime = 4.35661153e8, maxtime = 4.76736971e8, stdtime = 1.7975141515259467e7)
            # (transpose = false, NCPUs = 12, BLAS = "MKL", meantime = 4.224576985e8, mediantime = 4.194326665e8, mintime = 4.19017637e8, maxtime = 4.39155899e8, stdtime = 6.8465953122291835e6)
            # (transpose = true, NCPUs = 12, BLAS = "MKL", meantime = 1.478351086470588e8, mediantime = 1.286892225e8, mintime = 1.24138784e8, maxtime = 4.7547798e8, stdtime = 7.133122254958108e7)
            # (transpose = false, NCPUs = 16, BLAS = "Julia", meantime = 4.6000529445454544e8, mediantime = 4.59697238e8, mintime = 4.59161303e8, maxtime = 4.6214509e8, stdtime = 1.0435163653225917e6)
            # (transpose = true, NCPUs = 16, BLAS = "Julia", meantime = 4.469761928333333e8, mediantime = 4.36953207e8, mintime = 4.36762817e8, maxtime = 4.77925892e8, stdtime = 1.6369546538163424e7)
            # (transpose = false, NCPUs = 16, BLAS = "MKL", meantime = 4.292986104166667e8, mediantime = 4.04913057e8, mintime = 4.01455751e8, maxtime = 6.42039078e8, stdtime = 6.788290640751022e7)
            # (transpose = true, NCPUs = 16, BLAS = "MKL", meantime = 1.4853937970588234e8, mediantime = 1.0230638e8, mintime = 9.0928573e7, maxtime = 3.55159278e8, stdtime = 8.693982964312722e7)
            # (transpose = false, NCPUs = 20, BLAS = "Julia", meantime = 4.573211471818182e8, mediantime = 4.52035373e8, mintime = 4.51447991e8, maxtime = 4.88610971e8, stdtime = 1.1953679495507982e7)
            # (transpose = true, NCPUs = 20, BLAS = "Julia", meantime = 4.473714365e8, mediantime = 4.402625075e8, mintime = 4.28183657e8, maxtime = 4.71916532e8, stdtime = 1.917841557881981e7)
            # (transpose = false, NCPUs = 20, BLAS = "MKL", meantime = 5.407032969e8, mediantime = 4.85515017e8, mintime = 4.31955919e8, maxtime = 7.68051841e8, stdtime = 1.31720719035713e8)
            # (transpose = true, NCPUs = 20, BLAS = "MKL", meantime = 1.1821755760465117e8, mediantime = 8.6905782e7, mintime = 6.2529388e7, maxtime = 2.70172965e8, stdtime = 5.93575581419265e7)
            # (transpose = false, NCPUs = 24, BLAS = "Julia", meantime = 4.45926959e8, mediantime = 4.41228727e8, mintime = 4.40965224e8, maxtime = 4.78314826e8, stdtime = 1.0995497286179168e7)
            # (transpose = true, NCPUs = 24, BLAS = "Julia", meantime = 4.214945328333333e8, mediantime = 4.21170112e8, mintime = 4.20844407e8, maxtime = 4.23900531e8, stdtime = 877187.4122648875)
            # (transpose = false, NCPUs = 24, BLAS = "MKL", meantime = 4.6575263454545456e8, mediantime = 4.64615954e8, mintime = 4.64340919e8, maxtime = 4.7732185e8, stdtime = 3.844497948986951e6)
            # (transpose = true, NCPUs = 24, BLAS = "MKL", meantime = 1.078366069361702e8, mediantime = 1.0107288e8, mintime = 9.0110459e7, maxtime = 2.87826077e8, stdtime = 3.5226037921233065e7)
            # (transpose = false, NCPUs = 36, BLAS = "Julia", meantime = 4.531115790833333e8, mediantime = 4.525883615e8, mintime = 4.52128507e8, maxtime = 4.56282466e8, stdtime = 1.266477446152298e6)
            # (transpose = true, NCPUs = 36, BLAS = "Julia", meantime = 4.50372326e8, mediantime = 4.506963725e8, mintime = 4.29180808e8, maxtime = 4.71497003e8, stdtime = 2.056734779081425e7)
            # (transpose = false, NCPUs = 36, BLAS = "MKL", meantime = 4.5813660790909094e8, mediantime = 4.58375475e8, mintime = 4.44048091e8, maxtime = 4.69507888e8, stdtime = 6.672245444456813e6)
            # (transpose = true, NCPUs = 36, BLAS = "MKL", meantime = 1.0926559076595744e8, mediantime = 9.6175866e7, mintime = 9.0141027e7, maxtime = 2.8552446e8, stdtime = 4.366129822598922e7)
            # (transpose = false, NCPUs = 48, BLAS = "Julia", meantime = 4.392508005e8, mediantime = 4.391397825e8, mintime = 4.38751625e8, maxtime = 4.40350589e8, stdtime = 426710.2623721713)
            # (transpose = true, NCPUs = 48, BLAS = "Julia", meantime = 4.295984261666667e8, mediantime = 4.184454935e8, mintime = 4.17783303e8, maxtime = 4.56113981e8, stdtime = 1.7315670960982036e7)
            # (transpose = false, NCPUs = 48, BLAS = "MKL", meantime = 7.603238917142857e8, mediantime = 7.89538035e8, mintime = 4.48626537e8, maxtime = 1.145136288e9, stdtime = 3.1043062912782615e8)
            # (transpose = true, NCPUs = 48, BLAS = "MKL", meantime = 8.517274671186441e7, mediantime = 4.4519008e7, mintime = 2.8534268e7, maxtime = 3.15918892e8, stdtime = 8.45964133899557e7)
            # (NCPUs = 1, BLAS = "Julia", meantime = 4.62278333e8, mediantime = 4.55552711e8, mintime = 4.5467243e8, maxtime = 4.92527927e8, stdtime = 1.4817770561039506e7)
            # (NCPUs = 1, BLAS = "MKL", meantime = 6.04560891e8, mediantime = 6.03554815e8, mintime = 5.89523154e8, maxtime = 6.29752089e8, stdtime = 1.6405261745905e7)
            # (NCPUs = 1, BLAS = "CSR-tmul", meantime = 4.028326890769231e8, mediantime = 4.02279553e8, mintime = 4.01374638e8, maxtime = 4.05886627e8, stdtime = 1.5097978904339846e6)
            # (NCPUs = 1, BLAS = "CSR-bmul", meantime = 4.219883654166667e8, mediantime = 4.21441053e8, mintime = 4.20897836e8, maxtime = 4.25269667e8, stdtime = 1.4590923944203658e6)
            # (NCPUs = 2, BLAS = "Julia", meantime = 4.500307463333333e8, mediantime = 4.433706255e8, mintime = 4.43070197e8, maxtime = 4.84170604e8, stdtime = 1.5227198780427491e7)
            # (NCPUs = 2, BLAS = "MKL", meantime = 4.2699047975e8, mediantime = 4.152132185e8, mintime = 4.1378855e8, maxtime = 4.49414778e8, stdtime = 1.5787783019653404e7)
            # (NCPUs = 2, BLAS = "CSR-tmul", meantime = 2.0390314932e8, mediantime = 2.02295413e8, mintime = 2.01529769e8, maxtime = 2.22375889e8, stdtime = 4.212425236350262e6)
            # (NCPUs = 2, BLAS = "CSR-bmul", meantime = 2.175719264347826e8, mediantime = 2.18423457e8, mintime = 2.10732132e8, maxtime = 2.20955385e8, stdtime = 2.5461932044544606e6)
            (NCPUs = 1, BLAS = "Julia", meantime = 4.6583576145454544e8, mediantime = 4.5908871e8, mintime = 4.51912021e8, maxtime = 5.00148546e8, stdtime = 1.7573298204491504e7)
            (NCPUs = 1, BLAS = "MKL", meantime = 6.169911312222222e8, mediantime = 6.00700313e8, mintime = 5.99170842e8, maxtime = 6.45242046e8, stdtime = 2.0795401942313038e7)
            (NCPUs = 1, BLAS = "CSR-tmul", meantime = 4.14499527e8, mediantime = 4.13977536e8, mintime = 4.13233859e8, maxtime = 4.16624129e8, stdtime = 1.2726475089008347e6)
            (NCPUs = 1, BLAS = "CSR-bmul", meantime = 3.8553272553846157e8, mediantime = 3.85125397e8, mintime = 3.84322126e8, maxtime = 3.87642279e8, stdtime = 1.1686182321752624e6)
            (NCPUs = 2, BLAS = "Julia", meantime = 4.500202596666667e8, mediantime = 4.40334494e8, mintime = 4.3973049e8, maxtime = 4.79781954e8, stdtime = 1.7710026872141983e7)
            (NCPUs = 2, BLAS = "MKL", meantime = 4.308758598333333e8, mediantime = 4.18918103e8, mintime = 4.1828469e8, maxtime = 4.54099345e8, stdtime = 1.6024961103006408e7)
            (NCPUs = 2, BLAS = "CSR-tmul", meantime = 2.426639337142857e8, mediantime = 2.12501578e8, mintime = 2.0813021e8, maxtime = 4.20639346e8, stdtime = 7.347194228101255e7)
            (NCPUs = 2, BLAS = "CSR-bmul", meantime = 2.1671208520833334e8, mediantime = 2.166827695e8, mintime = 2.1556805e8, maxtime = 2.1800694e8, stdtime = 586265.892535711)
            (NCPUs = 3, BLAS = "Julia", meantime = 4.7233825372727275e8, mediantime = 4.65597959e8, mintime = 4.53985926e8, maxtime = 4.94518392e8, stdtime = 1.8574126068093356e7)
            (NCPUs = 3, BLAS = "MKL", meantime = 3.9327069084615386e8, mediantime = 3.82272547e8, mintime = 3.8155852e8, maxtime = 4.15311982e8, stdtime = 1.4679500213839415e7)
            (NCPUs = 3, BLAS = "CSR-tmul", meantime = 1.4583319894285715e8, mediantime = 1.45757912e8, mintime = 1.44322793e8, maxtime = 1.48322237e8, stdtime = 772642.6577086302)
            (NCPUs = 3, BLAS = "CSR-bmul", meantime = 1.8345215539285713e8, mediantime = 1.83419241e8, mintime = 1.83149848e8, maxtime = 1.83924514e8, stdtime = 184876.03058392828)
            (NCPUs = 4, BLAS = "Julia", meantime = 4.7430965809090906e8, mediantime = 4.67614759e8, mintime = 4.67173969e8, maxtime = 5.08673157e8, stdtime = 1.4932092276927615e7)
            (NCPUs = 4, BLAS = "MKL", meantime = 3.9834018153846157e8, mediantime = 3.92201046e8, mintime = 3.91971752e8, maxtime = 4.19208978e8, stdtime = 1.117040703242977e7)
            (NCPUs = 4, BLAS = "CSR-tmul", meantime = 1.2410687053658536e8, mediantime = 1.18221731e8, mintime = 1.17211893e8, maxtime = 1.72103891e8, stdtime = 1.515312541711425e7)
            (NCPUs = 4, BLAS = "CSR-bmul", meantime = 1.2973867707692307e8, mediantime = 1.16114065e8, mintime = 1.13201356e8, maxtime = 3.5416585e8, stdtime = 4.315727956011334e7)
            (NCPUs = 6, BLAS = "Julia", meantime = 4.5514322209090906e8, mediantime = 4.48654395e8, mintime = 4.48262756e8, maxtime = 4.87832046e8, stdtime = 1.44755971369501e7)
            (NCPUs = 6, BLAS = "MKL", meantime = 6.713506965e8, mediantime = 6.41598363e8, mintime = 5.92214387e8, maxtime = 8.22732235e8, stdtime = 9.256930717392786e7)
            (NCPUs = 6, BLAS = "CSR-tmul", meantime = 1.121697738478261e8, mediantime = 9.5276458e7, mintime = 8.8200223e7, maxtime = 3.23382881e8, stdtime = 4.871659409631042e7)
            (NCPUs = 6, BLAS = "CSR-bmul", meantime = 1.006501365e8, mediantime = 1.005192455e8, mintime = 1.00126112e8, maxtime = 1.02291686e8, stdtime = 463179.178196765)
            (NCPUs = 8, BLAS = "Julia", meantime = 4.6431017372727275e8, mediantime = 4.5496795e8, mintime = 4.4622815e8, maxtime = 4.86232469e8, stdtime = 1.854058042518827e7)
            (NCPUs = 8, BLAS = "MKL", meantime = 7.883905407142857e8, mediantime = 6.6586202e8, mintime = 5.87080366e8, maxtime = 1.321166058e9, stdtime = 2.617112311549547e8)
            (NCPUs = 8, BLAS = "CSR-tmul", meantime = 9.259070253703703e7, mediantime = 9.14083675e7, mintime = 9.0039023e7, maxtime = 1.24657237e8, stdtime = 5.966292222889591e6)
            (NCPUs = 8, BLAS = "CSR-bmul", meantime = 9.334718896296297e7, mediantime = 9.3115348e7, mintime = 9.2820029e7, maxtime = 9.5096868e7, stdtime = 596516.4276446401)
            (NCPUs = 12, BLAS = "Julia", meantime = 4.55957178e8, mediantime = 4.43053414e8, mintime = 4.4135856e8, maxtime = 4.80643858e8, stdtime = 1.7691859800031427e7)
            (NCPUs = 12, BLAS = "MKL", meantime = 5.362515351e8, mediantime = 5.32125366e8, mintime = 4.22008107e8, maxtime = 6.52514142e8, stdtime = 8.288033353422347e7)
            (NCPUs = 12, BLAS = "CSR-tmul", meantime = 9.99984858e7, mediantime = 8.5255936e7, mintime = 6.4741637e7, maxtime = 2.0423732e8, stdtime = 3.544362587387413e7)
            (NCPUs = 12, BLAS = "CSR-bmul", meantime = 9.15022222e7, mediantime = 9.1266564e7, mintime = 9.1060545e7, maxtime = 9.3728432e7, stdtime = 705418.8990191468)
            (NCPUs = 16, BLAS = "Julia", meantime = 4.6250807172727275e8, mediantime = 4.49493409e8, mintime = 4.48217917e8, maxtime = 4.87019784e8, stdtime = 1.897821567554414e7)
            (NCPUs = 16, BLAS = "MKL", meantime = 4.708429616363636e8, mediantime = 4.65929516e8, mintime = 4.63237569e8, maxtime = 4.84596507e8, stdtime = 8.407438770930959e6)
            (NCPUs = 16, BLAS = "CSR-tmul", meantime = 6.78760578108108e7, mediantime = 6.0549516e7, mintime = 5.0710359e7, maxtime = 1.63814213e8, stdtime = 2.4369239572309755e7)
            (NCPUs = 16, BLAS = "CSR-bmul", meantime = 6.1019996951219514e7, mediantime = 5.64277625e7, mintime = 5.6314897e7, maxtime = 1.52231926e8, stdtime = 1.6408329284242325e7)
            (NCPUs = 20, BLAS = "Julia", meantime = 4.519027624166667e8, mediantime = 4.51948113e8, mintime = 4.51542306e8, maxtime = 4.52137754e8, stdtime = 187214.65051668184)
            (NCPUs = 20, BLAS = "MKL", meantime = 5.10571902e8, mediantime = 4.65378662e8, mintime = 4.26856512e8, maxtime = 7.50507751e8, stdtime = 1.1480584481393579e8)
            (NCPUs = 20, BLAS = "CSR-tmul", meantime = 6.540467172727273e7, mediantime = 5.9663719e7, mintime = 5.1639302e7, maxtime = 1.50550288e8, stdtime = 1.829230241743437e7)
            (NCPUs = 20, BLAS = "CSR-bmul", meantime = 6.897503306849316e7, mediantime = 5.9899562e7, mintime = 5.2272845e7, maxtime = 2.04177594e8, stdtime = 3.0893855803945493e7)
            (NCPUs = 24, BLAS = "Julia", meantime = 4.55101717e8, mediantime = 4.52037076e8, mintime = 4.51759255e8, maxtime = 4.8377893e8, stdtime = 9.527181465267438e6)
            (NCPUs = 24, BLAS = "MKL", meantime = 5.023315307e8, mediantime = 4.441035685e8, mintime = 4.20634365e8, maxtime = 6.85390554e8, stdtime = 1.0818446505240606e8)
            (NCPUs = 24, BLAS = "CSR-tmul", meantime = 5.7985023302325584e7, mediantime = 5.2344588e7, mintime = 5.2045106e7, maxtime = 1.07448573e8, stdtime = 8.806190353342475e6)
            (NCPUs = 24, BLAS = "CSR-bmul", meantime = 9.177558947272727e7, mediantime = 8.925485e7, mintime = 8.3190799e7, maxtime = 1.62312992e8, stdtime = 1.4380460664546315e7)
            (NCPUs = 36, BLAS = "Julia", meantime = 5.638422828888888e8, mediantime = 4.87026616e8, mintime = 4.45008496e8, maxtime = 7.3329078e8, stdtime = 1.232456099203042e8)
            (NCPUs = 36, BLAS = "MKL", meantime = 5.434447704e8, mediantime = 5.19041098e8, mintime = 5.14858685e8, maxtime = 7.5843983e8, stdtime = 7.569725862001485e7)
            (NCPUs = 36, BLAS = "CSR-tmul", meantime = 7.362073208695652e7, mediantime = 6.2785466e7, mintime = 5.8108043e7, maxtime = 2.20269305e8, stdtime = 3.588675063651513e7)
            (NCPUs = 36, BLAS = "CSR-bmul", meantime = 7.441944882089552e7, mediantime = 6.5105618e7, mintime = 4.9333828e7, maxtime = 2.42889888e8, stdtime = 3.824778127422448e7)
            (NCPUs = 48, BLAS = "Julia", meantime = 4.553303191666667e8, mediantime = 4.4938369e8, mintime = 4.49042279e8, maxtime = 4.86874418e8, stdtime = 1.3618035837851485e7)
            (NCPUs = 48, BLAS = "MKL", meantime = 5.037321015e8, mediantime = 4.90758936e8, mintime = 4.48678888e8, maxtime = 5.97816192e8, stdtime = 5.030351385792964e7)
            (NCPUs = 48, BLAS = "CSR-tmul", meantime = 5.804969129069767e7, mediantime = 5.4469401e7, mintime = 4.5689475e7, maxtime = 1.70802124e8, stdtime = 2.0123236253490612e7)
            (NCPUs = 48, BLAS = "CSR-bmul", meantime = 9.761288211320755e7, mediantime = 7.8427757e7, mintime = 3.4658509e7, maxtime = 2.79920545e8, stdtime = 7.1115256439898e7)
        ]
]

using GLMakie

iJulia = findall(df.BLAS .== "Julia")
iMKL = findall(df.BLAS .== "MKL")

# fig = Figure(size = (600, 1200))
fig = Figure(size = (600, 600))

# for (itranspose, transpose) in enumerate((false, true))
itranspose = 1
    ax = Axis(
        fig[itranspose, 1];
        xlabel = "Number of CPUs",
        ylabel = "Benchmark Time (s)",
        # xticks = unique(df.NCPUs),
        # yticks = [0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000],
        yticks = BaseMulTicks([1, 2, 3, 4, 5, 7]),
        yminorticks = BaseMulTicks(2:9),
        xticks = BaseMulTicks([1, 2, 3, 4, 5, 7]),
        xminorticks = BaseMulTicks(2:9),
        # limits = (nothing, nothing, 0, nothing),
        xscale = log10,
        yscale = log10,
        # title = transpose ? "Transposed Matrix-Vector Multiplication" : "Matrix-Vector Multiplication",
        xminorticksvisible = true,
        yminorticksvisible = true,
    )

    # df2 = filter(row -> row.transpose == transpose, df)
    df2 = df

    for (iblas, blas) in enumerate(unique(df2.BLAS))
        df3 = filter(row -> row.BLAS == blas, df2)
        lines!(df3.NCPUs, 1.0e-9 * df3.mediantime, label = blas, color = Cycled(iblas))
        lines!(df3.NCPUs, 1.0e-9 * df3.meantime, color = Cycled(iblas), linestyle = :dash)
        band!(df3.NCPUs, 1.0e-9 * df3.mintime, 1.0e-9 * df3.maxtime, color = Cycled(iblas), alpha = 0.3)
    end

    axislegend(ax; position = :rt)

# end

# linkyaxes!([ax for ax in fig.content if ax isa Axis]...)

fig
save("/Users/benoitpasquier/tmp/MKLSparse_Benchmark.png", fig)
