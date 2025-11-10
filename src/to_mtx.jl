using Pkg
Pkg.instantiate()
using FileIO
using MatrixMarket
using Format
using SparseArrays

κVdeep = 3.0e-5 # m^2/s
κVML = 1.0      # m^2/s
κH = 300.0      # m^2/s (grid-scaling by sqrt(area))

inputdirs = [
    "/scratch/y99/TMIP/data/ACCESS-OM2-1/1deg_jra55_iaf_omip2_cycle6/Jan1960-Dec1979"     # 1 degree
    "/scratch/y99/TMIP/data/ACCESS-OM2-025/025deg_jra55_iaf_omip2_cycle6/Jan1960-Dec1979" # 0.25 degree
]

κH_gridscalings = [
    1      # 1 degree
    4      # 0.25 degree
]

# Remove surface layer for matrix market save,
# so it is smaller and then the RHS is just ones.
Nsrfs = [
    69809      # 1 degree
    970921     # 0.25 degree
]

for (inputdir, κH_gridscaling, Nsrf) in zip(inputdirs, κH_gridscalings, Nsrfs)
    κVdeep_str = "kVdeep" * format(κVdeep, conversion = "e")
    κVML_str = "kVML" * format(κVML, conversion = "e")
    κH_str = "kH" * format(κH / κH_gridscaling, conversion = "d")

    TMfile = joinpath(inputdir, "yearly_matrix_$(κVdeep_str)_$(κH_str)_$(κVML_str).jld2")
    @info "Loading matrix from $TMfile"
    T = load(TMfile, "T")
    @info "Matrix size: $(size(T)), nnz = $(nnz(T))"

    T = T[(Nsrf + 1):end, (Nsrf + 1):end]
    @info """Removing surface layer:
    New matrix size: $(size(T)), nnz = $(nnz(T))"""

    MMfile = joinpath(inputdir, "yearly_matrix_$(κVdeep_str)_$(κH_str)_$(κVML_str).mtx")
    @info "Saving matrix to Matrix Market file:\n  $(MMfile)"
    MatrixMarket.mmwrite(MMfile, T)
end
