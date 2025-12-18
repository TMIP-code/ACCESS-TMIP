using Pkg
Pkg.activate(".")
Pkg.instantiate()

using SparseArrays
using LinearAlgebra
using Statistics
using Format
using FileIO
using BenchmarkTools
using CUDA
CUDA.set_runtime_version!(v"12.9.1")
@show CUDA.versioninfo()
using CUDA.CUSPARSE, CUDA.CUSOLVER

model = "ACCESS-OM2-025"
experiment = "025deg_jra55_iaf_omip2_cycle6"
time_window = "Jan1960-Dec1979"
inputdir = "/scratch/y99/TMIP/data/$model/$experiment/$time_window"

# preferred diffusivities
κVdeep = 3.0e-5 # m^2/s
κVML = 1.0      # m^2/s
κH = 300.0 / 4  # m^2/s (grid-scaling by sqrt(area))

κVdeep_str = "kVdeep" * format(κVdeep, conversion = "e")
κVML_str = "kVML" * format(κVML, conversion = "e")
κH_str = "kH" * format(κH, conversion = "d")

TMfile = joinpath(inputdir, "yearly_matrix_$(κVdeep_str)_$(κH_str)_$(κVML_str).jld2")
T = load(TMfile, "T")
N = size(T, 1)

A_gpu = CuSparseMatrixCSR(T)
u_gpu = CuVector(ones(N))
du_gpu = CuVector(ones(N))

@time "mul!(du, A, u) (1st)" mul!(du_gpu, A_gpu, u_gpu);
@time "mul!(du, A, u)" mul!(du_gpu, A_gpu, u_gpu);
@time "du .= A * u (1st)" du_gpu .= A_gpu * u_gpu;
@time "du .= A * u" du_gpu .= A_gpu * u_gpu;

println("Done.")
