using Pkg
Pkg.activate(".")
Pkg.instantiate()

using SparseArrays
using LinearAlgebra
using Statistics
using Format
using FileIO
using BenchmarkTools
if ENV["sparseBLAS"] == "MKL"
    using MKLSparse
elseif contains(ENV["sparseBLAS"], "CSR")
    using SparseMatricesCSR
    using ThreadedSparseCSR
    if contains(ENV["sparseBLAS"], "tmul")
        ThreadedSparseCSR.multithread_matmul(BaseThreads())
    elseif contains(ENV["sparseBLAS"], "bmul")
        ThreadedSparseCSR.multithread_matmul(PolyesterThreads())
    end
end

println("queue: ", ENV["PBS_O_QUEUE"])
println("VMEM: ", parse(Int, ENV["PBS_VMEM"]) ÷ 10^9, "GB")
println("NCPUS: ", parse(Int, ENV["PBS_NCPUS"]))
println("SparseBLAS: ", ENV["sparseBLAS"])

model = "ACCESS-OM2-01"
experiment = "01deg_jra55v140_iaf_cycle4"
time_window = "Jan1960-Dec1979"
inputdir = "/scratch/y99/TMIP/data/$model/$experiment/$time_window"

# preferred diffusivities
κVdeep = 3.0e-5 # m^2/s
κVML = 1.0      # m^2/s
κH = 300.0 / 10 # m^2/s (grid-scaling by sqrt(area))

κVdeep_str = "kVdeep" * format(κVdeep, conversion = "e")
κVML_str = "kVML" * format(κVML, conversion = "e")
κH_str = "kH" * format(κH, conversion = "d")

TMfile = joinpath(inputdir, "yearly_matrix_$(κVdeep_str)_$(κH_str)_$(κVML_str).jld2")
T = load(TMfile, "T")
if contains(ENV["sparseBLAS"], "CSR")
    i, j, v = findnz(T)
    T = sparsecsr(i, j, v, size(T)...)
end

v = ones(size(T, 2))
dv = similar(v)
b1 = @benchmark mul!($dv, $T, $v);

row1 = (;
    NCPUs = parse(Int, ENV["PBS_NCPUS"]),
    BLAS = ENV["sparseBLAS"],
    meantime = mean(b1).time,
    mediantime = median(b1).time,
    mintime = minimum(b1).time,
    maxtime = maximum(b1).time,
    stdtime = std(b1).time,
)

@show row1

# T2 = transpose(T)
# b2 = @benchmark mul!($dv, $T2, $v);

# row2 = (;
#     NCPUs = parse(Int, ENV["PBS_NCPUS"]),
#     BLAS = ENV["sparseBLAS"],
#     meantime = mean(b2).time,
#     mediantime = median(b2).time,
#     mintime = minimum(b2).time,
#     maxtime = maximum(b2).time,
#     stdtime = std(b2).time,
# )

# @show row2

println("Done.")
