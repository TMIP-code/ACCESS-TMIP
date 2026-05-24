# `compute_age_ACCESS-OM2-01_NxN` runs — log survey

Survey of all PBS submissions of the plain `NxN` variants of
`compute_age_ACCESS-OM2-01` (i.e. the base script and its `2x2` / `5x5`
coarsening counterparts). The hemispheric-mask variants such as `2x235S`,
`2x240S`, …, `5x535S` are intentionally excluded.

Scripts:

- `1x1` (no coarsening): [src/compute_age_ACCESS-OM2-01.jl](src/compute_age_ACCESS-OM2-01.jl) / [src/compute_age_ACCESS-OM2-01.sh](src/compute_age_ACCESS-OM2-01.sh)
- `2x2`: [src/compute_age_ACCESS-OM2-01_2x2.jl](src/compute_age_ACCESS-OM2-01_2x2.jl) / [src/compute_age_ACCESS-OM2-01_2x2.sh](src/compute_age_ACCESS-OM2-01_2x2.sh)
- `5x5`: [src/compute_age_ACCESS-OM2-01_5x5.jl](src/compute_age_ACCESS-OM2-01_5x5.jl) / [src/compute_age_ACCESS-OM2-01_5x5.sh](src/compute_age_ACCESS-OM2-01_5x5.sh)

## Summary by coarsening factor

| NxN | Attempts | Converged? | Best wall time | Best mem used | Failure mode of failed runs |
|-----|---------:|------------|----------------|---------------|-----------------------------|
| 1x1 | 6 | No | — | — | 5 early script/data errors; 1 reached the Pardiso solve and failed with `InexactError: trunc(Int32, 2147483653)` — MKL Pardiso 32-bit index overflow on a 3.5 × 10⁸ × 3.5 × 10⁸ matrix with 2.4 × 10⁹ nnz. **Not OOM** (used 169 GB / 2.83 TB requested). |
| 2x2 | 1 | **Yes** (exit 0) | **09:16:31** | **2.34 TB** (of 2.83 TB requested) | n/a |
| 5x5 | 3 | No | — | — | All three reached `LUMP and SPRAY: Matrix size reduction: 96.0%` and then failed with `UndefVarError` (`v`, `volcello_ds`, `experiment`) — script bugs in the post-LUMP code paths, **not OOM and not solver convergence failures**. Peak mem reached 454 GB / 735 GB requested. |

Caveats:

- For the failed `5x5` runs the linear solve may or may not have actually
  completed before the `UndefVarError` — the logs show the LUMP/SPRAY reduction
  message but do not clearly bracket the Pardiso solve. Treat "did 5x5 converge?"
  as unknown from these logs.
- The single `2x2` success used 2.34 TB out of 2.83 TB requested, i.e. ~83 % of
  the memory ceiling.

## Per-job details

Columns: NxN, PBS job id (→ PBS resource log), julia stdout/stderr log,
walltime used, memory used, exit status, and a one-line outcome.

| NxN | Job (PBS log) | Julia log | Walltime | Mem used | Exit | Outcome |
|-----|---------------|-----------|---------:|---------:|-----:|---------|
| 1x1 | [151967081](output/PBS/151967081.gadi-pbs.OU) | [151967081 julia](output/151967081.gadi-pbs.compute_age_ACCESS-OM2-01.out) | 00:02:18 | 2.63 GB | 1 | `NetCDF error code 2: No such file or directory` at [compute_age_ACCESS-OM2-01.jl:50](src/compute_age_ACCESS-OM2-01.jl#L50) |
| 1x1 | [151967713](output/PBS/151967713.gadi-pbs.OU) | [151967713 julia](output/151967713.gadi-pbs.compute_age_ACCESS-OM2-01.out) | 00:00:48 | 11.87 GB | 1 | `UndefVarError: 'model' not defined` at [compute_age_ACCESS-OM2-01.jl:59](src/compute_age_ACCESS-OM2-01.jl#L59) |
| 1x1 | [151972267](output/PBS/151972267.gadi-pbs.OU) | [151972267 julia](output/151972267.gadi-pbs.compute_age_ACCESS-OM2-01.out) | 00:00:52 | 8.18 GB | 1 | `KeyError: key "_FillValue" not found` in `makegridmetrics` at [compute_age_ACCESS-OM2-01.jl:72](src/compute_age_ACCESS-OM2-01.jl#L72) |
| 1x1 | [151994252](output/PBS/151994252.gadi-pbs.OU) | [151994252 julia](output/151994252.gadi-pbs.compute_age_ACCESS-OM2-01.out) | 00:01:25 | 15.33 GB | 1 | Same `KeyError: "_FillValue"` at [compute_age_ACCESS-OM2-01.jl:72](src/compute_age_ACCESS-OM2-01.jl#L72) |
| 1x1 | [151994733](output/PBS/151994733.gadi-pbs.OU) | [151994733 julia](output/151994733.gadi-pbs.compute_age_ACCESS-OM2-01.out) | 00:00:51 | 7.44 GB | 1 | Same `KeyError: "_FillValue"` at [compute_age_ACCESS-OM2-01.jl:72](src/compute_age_ACCESS-OM2-01.jl#L72) |
| 1x1 | [151994957](output/PBS/151994957.gadi-pbs.OU) | [151994957 julia](output/151994957.gadi-pbs.compute_age_ACCESS-OM2-01.out) | 00:13:31 | 169.32 GB | 1 | Reached Pardiso solve, failed with `InexactError: trunc(Int32, 2147483653)` (32-bit index overflow) at [compute_age_ACCESS-OM2-01.jl:99](src/compute_age_ACCESS-OM2-01.jl#L99) |
| 2x2 | [152034390](output/PBS/152034390.gadi-pbs.OU) | [152034390 julia](output/152034390.gadi-pbs.compute_age_ACCESS-OM2-01_2x2.out) | **09:16:31** | **2.34 TB** | **0** | **Converged.** Saved `steady_age_kVdeep3e-05_kH300_kVML1e+00_2x2.nc`. |
| 5x5 | [151997332](output/PBS/151997332.gadi-pbs.OU) | [151997332 julia](output/151997332.gadi-pbs.compute_age_ACCESS-OM2-01_5x5.out) | 00:14:38 | 174.73 GB | 1 | `UndefVarError: 'v' not defined` at [compute_age_ACCESS-OM2-01_5x5.jl:97](src/compute_age_ACCESS-OM2-01_5x5.jl#L97) (failed inside/around the `lump_and_spray` call) |
| 5x5 | [152018949](output/PBS/152018949.gadi-pbs.OU) | [152018949 julia](output/152018949.gadi-pbs.compute_age_ACCESS-OM2-01_5x5.out) | 00:57:27 | 441.17 GB | 1 | Got past LUMP/SPRAY (96.0 % size reduction). `UndefVarError: 'volcello_ds' not defined` at [compute_age_ACCESS-OM2-01_5x5.jl:111](src/compute_age_ACCESS-OM2-01_5x5.jl#L111) |
| 5x5 | [152027539](output/PBS/152027539.gadi-pbs.OU) | [152027539 julia](output/152027539.gadi-pbs.compute_age_ACCESS-OM2-01_5x5.out) | 00:45:20 | 454.20 GB | 1 | Got past LUMP/SPRAY (96.0 % size reduction). `UndefVarError: 'experiment' not defined` at [compute_age_ACCESS-OM2-01_5x5.jl:112](src/compute_age_ACCESS-OM2-01_5x5.jl#L112) |

### Resource requests

For context (taken from the PBS resource-usage block at the foot of each PBS log):

- `1x1`: 48 NCPUs, mem request 2.83 TB (`megamem`), walltime request 48:00:00.
- `2x2`: 48 NCPUs, mem request 2.83 TB (`megamem`), walltime request 48:00:00.
- `5x5`: 24 NCPUs, mem request 735.0 GB, walltime request 48:00:00.

Job timestamps (PBS "Resource Usage on …"):

- 1x1: 2025-10-09 14:05 → 22:53 (six consecutive attempts the same day)
- 5x5: 2025-10-10 00:28 → 11:53 (three attempts)
- 2x2: 2025-10-11 05:01 (single successful run)
