# qsub -I -P xv83 -l mem=64GB -l storage=scratch/gh0+scratch/xv83 -l walltime=01:00:00 -l ncpus=48

using Pkg
Pkg.activate(".")
Pkg.instantiate()


ENV["JULIA_CONDAPKG_BACKEND"] = "Null"
using OceanTransportMatrixBuilder
using NetCDF
using YAXArrays
using DataFrames
using DimensionalData
using SparseArrays
using LinearAlgebra
using Unitful
using Unitful: s, yr, d
using Statistics
using Format
using Dates
using FileIO
using ProgressMeter
try
    using CairoMakie
catch
    using CairoMakie
end
using GeoMakie
using OceanBasins
using NaNStatistics



# include("plotting_functions.jl")


# script options
@show model = "ACCESS-ESM1-5"
if isempty(ARGS)
    member = "r1i1p1f1"
    # experiment = "historical"
    # time_window = "Jan1850-Dec1859"
    # time_window = "Jan1990-Dec1999"
    experiment = "ssp370"
    # time_window = "Jan2030-Dec2039"
    time_window = "Jan2090-Dec2099"
    WRITEDATA = "true"
else
    experiment, member, time_window, finalmonth, WRITEDATA = ARGS
end
WRITEDATA = parse(Bool, WRITEDATA)
@show experiment
@show member
@show time_window


# Load areacello and volcello for grid geometry
fixedvarsinputdir = "/scratch/xv83/TMIP/data/$model"
volcello_ds = open_dataset(joinpath(fixedvarsinputdir, "volcello.nc"))
areacello_ds = open_dataset(joinpath(fixedvarsinputdir, "areacello.nc"))

members = ["r$(i)i1p1f1" for i in 1:3]

for member in members

    # Gadi directory for input files
    lumpby = "month"
    inputdir = "/scratch/xv83/TMIP/data/$model/$experiment/$member/$(time_window)"
    cycloinputdir = joinpath(inputdir, "cyclo$lumpby")



    # Load fixed variables in memory
    areacello = readcubedata(areacello_ds.areacello)
    volcello = readcubedata(volcello_ds.volcello)
    lon = readcubedata(volcello_ds.lon)
    lat = readcubedata(volcello_ds.lat)
    lev = volcello_ds.lev
    # Identify the vertices keys (vary across CMIPs / models)
    volcello_keys = propertynames(volcello_ds)
    lon_vertices_key = volcello_keys[findfirst(x -> occursin("lon", x) & occursin("vert", x), string.(volcello_keys))]
    lat_vertices_key = volcello_keys[findfirst(x -> occursin("lat", x) & occursin("vert", x), string.(volcello_keys))]
    lon_vertices = readcubedata(getproperty(volcello_ds, lon_vertices_key))
    lat_vertices = readcubedata(getproperty(volcello_ds, lat_vertices_key))

    # Make makegridmetrics
    gridmetrics = makegridmetrics(; areacello, volcello, lon, lat, lev, lon_vertices, lat_vertices)
    (; lon_vertices, lat_vertices, v3D) = gridmetrics

    # Make indices
    indices = makeindices(v3D)
    (; N, wet3D) = indices

    finalmonths = 1:12


    @time "building mean_days_in_months" mean_days_in_months = map(finalmonths) do m
        inputfile = joinpath(cycloinputdir, "cyclo_matrix_$m.jld2")
        load(inputfile, "mean_days_in_month")
    end
    # So the δt that multiplies M̃ₜ is δ(t..t+1)
    # which is 0.5 of the mean days in finalmonths k and k+1
    δts = map(finalmonths) do m
        ustrip(s, (mean_days_in_months[mod1(m + 1, 12)] + mean_days_in_months[m]) / 2 * d)
    end
    # Weights to be used to combine the ℊ̃ from each final month
    ωs = δts / sum(δts)

    # # quick check that it makes sense
    Δt = sum(δts)

    # Take the average of ℰ over the 12 possible final finalmonths
    # (weighted by the duration of that month)
    ℰ = sum(map(finalmonths) do finalmonth
        # Load from NetCDF file
        finalmonthstr = format(finalmonth, width = 2, zeropadding = true)
        # FIXME this below will need to be uncommented and the line below removed
        # because I messed up the file name to save to (overwrote ℊ̃ with ℰ)
        outputfile = joinpath(inputdir, "seqeff_$(finalmonthstr).nc")
        @info "Loading adjoint propagrator as netCDF file:\n  $(outputfile)"
        ds = open_dataset(outputfile)
        ωs[finalmonth] * readcubedata(ds.seqeff)
    end)

    # save data
    cube3D = rebuild(volcello_ds["volcello"];
        data = ℰ.data,
        dims = dims(ℰ),
        metadata = Dict(
            "origin" => "cyclo-stationary sequestration efficiency, calE, averaged over all final months",
            "model" => model,
            "experiment" => experiment,
            "member" => member,
            "time window" => time_window,
            "units" => "",
            "Ti unit" => "yr",
        )
    )

    arrays = Dict(:calE => cube3D, :lat => volcello_ds.lat, :lon => volcello_ds.lon)
    ds = Dataset(; volcello_ds.properties, arrays...)

    # Save to netCDF file
    outputfile = joinpath(inputdir, "calE.nc")
    @info "Saving mean sequestration efficiency as netCDF file:\n  $(outputfile)"
    savedataset(ds, path = outputfile, driver = :netcdf, overwrite = true)


    # Do the same for 𝒢̃
    𝒢̃ = sum(map(finalmonths) do finalmonth
        # Load from NetCDF file
        finalmonthstr = format(finalmonth, width = 2, zeropadding = true)
        # FIXME this below will need to be uncommented and the line below removed
        # because I messed up the file name to save to (overwrote ℊ̃ with ℰ)
        outputfile = joinpath(inputdir, "calgtilde_$(finalmonthstr).nc")
        @info "Loading adjoint propagrator as netCDF file:\n  $(outputfile)"
        ds = open_dataset(outputfile)
        ωs[finalmonth] * readcubedata(ds.calgtilde)
    end)

    # save data
    cube3D = rebuild(volcello_ds["volcello"];
        data = 𝒢̃.data,
        dims = dims(𝒢̃),
        metadata = Dict(
            "origin" => "cyclo-stationary adjoint propagator, calgtilde, averaged over all final months",
            "model" => model,
            "experiment" => experiment,
            "member" => member,
            "time window" => time_window,
            "units" => "",
            "Ti unit" => "yr",
        )
    )

    arrays = Dict(:calgtilde => cube3D, :lat => volcello_ds.lat, :lon => volcello_ds.lon)
    ds = Dataset(; volcello_ds.properties, arrays...)

    # Save to netCDF file
    outputfile = joinpath(inputdir, "calgtilde.nc")
    @info "Saving mean adjoint propagator as netCDF file:\n  $(outputfile)"
    savedataset(ds, path = outputfile, driver = :netcdf, overwrite = true)


end

