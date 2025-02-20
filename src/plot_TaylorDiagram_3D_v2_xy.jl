# # qsub -I -P xv83 -l mem=64GB -l storage=scratch/gh0+scratch/xv83 -l walltime=02:00:00 -l ncpus=12

# using Pkg
# Pkg.activate(".")
# Pkg.instantiate()

# using OceanTransportMatrixBuilder
# using NetCDF
# using YAXArrays
# using DataFrames
# using DimensionalData
# # using SparseArrays
# # using LinearAlgebra
# using Unitful
# using Unitful: s, yr
# try
#     using CairoMakie
# catch
#     using CairoMakie
# end
# using GeoMakie
# using Interpolations
# using OceanBasins
# using Statistics
# using NaNStatistics
# using StatsBase
# using FileIO
# using Contour
# using GeometryBasics
# using GeometryOps
# using LibGEOS
# # using LaTeXStrings
# using Format

# model = "ACCESS-ESM1-5"

# time_window = "Jan1850-Dec1859"
# experiment = "historical"
# member = "r20i1p1f1"

# # Gadi directory for input files
# # inputdirfun(member) = "/scratch/xv83/TMIP/data/$model/$experiment/all members/$(time_window)"
# inputdir = "/scratch/xv83/TMIP/data/$model/$experiment/$(member)/$(time_window)/cyclomonth"
# outputdir = inputdir

# # Load areacello and volcello for grid geometry
# fixedvarsinputdir = "/scratch/xv83/TMIP/data/$model"
# volcello_ds = open_dataset(joinpath(fixedvarsinputdir, "volcello.nc"))
# areacello_ds = open_dataset(joinpath(fixedvarsinputdir, "areacello.nc"))

# # Load fixed variables in memory
# areacello = readcubedata(areacello_ds.areacello)
# volcello = readcubedata(volcello_ds.volcello)
# lon = readcubedata(volcello_ds.lon)
# lat = readcubedata(volcello_ds.lat)
# lev = volcello_ds.lev
# # Identify the vertices keys (vary across CMIPs / models)
# volcello_keys = propertynames(volcello_ds)
# lon_vertices_key = volcello_keys[findfirst(x -> occursin("lon", x) & occursin("vert", x), string.(volcello_keys))]
# lat_vertices_key = volcello_keys[findfirst(x -> occursin("lat", x) & occursin("vert", x), string.(volcello_keys))]
# lon_vertices = readcubedata(getproperty(volcello_ds, lon_vertices_key))
# lat_vertices = readcubedata(getproperty(volcello_ds, lat_vertices_key))
# # Make makegridmetrics
# gridmetrics = makegridmetrics(; areacello, volcello, lon, lat, lev, lon_vertices, lat_vertices)
# (; lon_vertices, lat_vertices, lon, lat, zt, v3D, thkcello, Z3D) = gridmetrics
# lev = zt
# # Make indices
# indices = makeindices(gridmetrics.v3D)
# (; wet3D, N) = indices




# # # Matrix ages for varied diffusivities
# # @info "Loading age computed from matrices with different diffusivities"
# # κVdeeps = [3e-8, 1e-7, 3e-7, 1e-6, 3e-6] # m^2/s
# # κVMLs = [0, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2] # m^2/s
# # κHs = [0, 1, 3, 10, 30, 100] # m^2/s

# # κs = Iterators.product(κVdeeps, κVMLs, κHs)
# # model_data = map(κs) do κ
# #     κVdeep, κVML, κH = κ
# #     κVdeep_str = "kVdeep" * format(κVdeep, conversion="e")
# #     κVML_str = "kVML" * format(κVML, conversion="e")
# #     κH_str = "kH" * format(κH, conversion="d")
# #     f = "ideal_mean_age_$(κVdeep_str)_$(κH_str)_$(κVML_str).nc"
# #     age_ds = open_dataset(joinpath(inputdir, f))
# #     age3D = (age_ds.age[Ti=1].data + age_ds.age[Ti=12].data) / 2
# #     age3D[wet3D]
# # end

# isvalidagefile(f) = startswith(f, "ideal_mean_age") && contains(f, "kVdeep") && contains(f, "kVML") && contains(f, "kH")
# function parse_κs(f)
#     κVdeep_str, κH_str, κVML_str = split(f, "_")[4:6]
#     κVdeep = parse(Float64, κVdeep_str[7:11])
#     κVML = parse(Float64, κVML_str[5:9])
#     κH = parse(Float64, κH_str[3:end])
#     return κVdeep, κVML, κH
# end
# files = [f for f in readdir(inputdir) if isvalidagefile(f)]
# model_data = map(files) do f
#     age_ds = open_dataset(joinpath(inputdir, f))
#     age3D = (age_ds.age[Ti=1].data + age_ds.age[Ti=12].data) / 2
#     age3D[wet3D]
# end
# κs = parse_κs.(files)
# κVdeeps = map(x -> x[1], κs)
# κVMLs = map(x -> x[2], κs)
# κHs = map(x -> x[3], κs)


# @info "Loading Anderson Acceleration age"
# # Anderson Accelerated age
# AAfile = "/scratch/xv83/bp3051/access-esm/archive/andersonacceleration_test-n10-5415f621/age_output/ocean_age.res_0035.nc"
# obs_ds = open_dataset(AAfile)
# obs_data = obs_ds.age_global[Time=1].data[wet3D]
# # iters = 0:35
# # obs_data2 = map(iters) do iter
# #     AAfile = joinpath(AAdir, "ocean_age.res_$(format(iter; width = 4, zeropadding = true)).nc")
# #     obs_ds = open_dataset(AAfile)
# #     obs_data = obs_ds.age_global[Time=1].data[wet3D]
# # end

# @info "Loading age with same constants as Matt"
# model_data2 = begin
#     f = "/scratch/xv83/TMIP/data/ACCESS-ESM1-5/historical/r20i1p1f1/Jan1850-Dec1859/cyclomonth/ideal_mean_age_kVdeep3e-05_kH500.nc"
#     age_ds = open_dataset(joinpath(inputdir, f))
#     age3D = (age_ds.age[Ti=1].data + age_ds.age[Ti=12].data) / 2
#     age3D[wet3D]
# end


# # Taylor diagram function that returns all the required values
# # notation taken from the original Taylor paper
# # TODO check that the identity holds when using weights!
# function taylordiagramvalues(f, r, args...)

#     # STDs and means
#     σf = std(f, args...; corrected = false)
#     σr = std(r, args...; corrected = false)
#     f̄ = mean(f, args...)
#     r̄ = mean(r, args...)

#     # Correlation coefficient
#     R = cor([f r], args...)[2]

#     # Root Mean Square Difference
#     E = sqrt(mean((f .- r) .^ 2, args...))

#     # Bias
#     Ē = f̄ - r̄

#     # Centered Root Mean Square Difference
#     E′ = sqrt(mean(((f .- f̄) - (r .- r̄)) .^ 2, args...))

#     # Full Mean Square Difference
#     E² = E′^2 + Ē^2

#     # Normalized values (maybe that needs to be a kwarg)
#     Ê′ = E′ / σr
#     σ̂f = σf / σr
#     σ̂r = 1.0

#     return (; σr, σf, R, E, Ē, E′, E², Ê′, σ̂f, σ̂r, f̄, r̄)
# end

# # Calculate the Taylor diagram values
# w = weights(v3D[wet3D])
# TDvals = [taylordiagramvalues(data, obs_data, w) for data in model_data]
# TDvals2 = taylordiagramvalues(model_data2, obs_data, w)
# # TDvalsobs = [taylordiagramvalues(data, obs_data, w) for data in obs_data2]
σr = TDvals[1].σr
σmax = 1.25TDvals[1].σr
σmax2 = (1.25 + 0.1) * TDvals[1].σr

σfs = [vals.σf for vals in TDvals]
Rs = [vals.R for vals in TDvals]
E′s = [vals.E′ for vals in TDvals]
Ēs = [vals.Ē for vals in TDvals]
Es = [vals.E for vals in TDvals]

# taken from MakieExtra (not using Pkg because it has outdated Makie dep)
using Makie.IntervalSets
Makie.project(s, r::HyperRectangle) = HyperRectangle(Makie.project(s, r.origin), Makie.project(s, r.origin + r.widths) - Makie.project(s, r.origin))
corner(r::HyperRectangle{2}, which::NTuple{2,Integer}) = Makie.Point(extrema(r)[_which_to_ix(which[1])][1], extrema(r)[_which_to_ix(which[2])][2])
_which_to_ix(which::Integer) = which == -1 ? 1 : which == 1 ? 2 : error("which must be -1 or 1, got $which")
fullproject(ax, p) = Makie.project(Makie.get_scene(ax), Makie.apply_transform(Makie.transform_func(ax), p)) + viewport(ax)[].origin
Base.:(⊆)(a::HyperRectangle, b::HyperRectangle) = all(map(⊆, intervals(a), intervals(b)))
intervals(r::HyperRectangle) = Interval.(r.origin, r.origin + r.widths)
function zoom_lines!(ax1, ax2; strokewidth=1.5, strokecolor=:black, color=(:black, 0), rectattrs=(;), lineattrs=(;))
    pscene = parent(parent(Makie.parent_scene(ax1)))
    @assert parent(parent(Makie.parent_scene(ax2))) === pscene
    obs = lift(ax1.finallimits, ax2.finallimits, ax1.scene.viewport, ax2.scene.viewport, ax1.scene.camera.projectionview, ax2.scene.camera.projectionview, Makie.transform_func(ax1), Makie.transform_func(ax2)) do _...
        lims = [ax1.finallimits[], ax2.finallimits[]]
        axs = lims[1] ⊆ lims[2] ? (ax1, ax2) :
              lims[2] ⊆ lims[1] ? (ax2, ax1) :
              nothing
        slines = if isnothing(axs)
            nothing
        else
            r1 = fullproject(axs[1], axs[1].finallimits[])
            r2 = fullproject(axs[2], axs[1].finallimits[])
            # cornsets = [
            #     ((corner(r1, (1,1)), corner(r2, (-1,1))), (corner(r1, (1,-1)), corner(r2, (-1,-1)))),
            #     ((corner(r1, (1,-1)), corner(r2, (1,1))), (corner(r1, (-1,-1)), corner(r2, (-1,1)))),
            #     ((corner(r1, (-1,-1)), corner(r2, (1,-1))), (corner(r1, (-1,1)), corner(r2, (1,1)))),
            #     ((corner(r1, (-1,1)), corner(r2, (-1,-1))), (corner(r1, (1,1)), corner(r2, (1,-1)))),
            # ]
            # argmin(cornsets) do ((a1, a2), (b1, b2))
            #     min(norm(a1-a2), norm(b1-b2))
            # end
            # BP: below is my zoom lines that does not work in general
            # cornsets2 = [
            #     ((corner(r1, (1,1)), corner(r2, (1,1))), (corner(r1, (-1,-1)), corner(r2, (-1,-1)))),
            #     ((corner(r1, (1,-1)), corner(r2, (1,-1))), (corner(r1, (-1,1)), corner(r2, (-1,1)))),
            # ]
            # argmin(cornsets2) do ((a1, a2), (b1, b2))
            #     max(norm(a1-a2), norm(b1-b2))
            # end
            [
                corner(r1, (1,1)), corner(r2, (1,1)),
                corner(r1, (-1,-1)), corner(r2, (-1,-1)),
                corner(r1, (1,-1)), corner(r2, (1,-1)),
                corner(r1, (-1,1)), corner(r2, (-1,1)),
            ]
            # BP: below is my zoom lines that does not work in general
            # (corner(r1, (1,-1)), corner(r2, (1,-1))), (corner(r1, (-1,1)), corner(r2, (-1,1)))
            # (corner(r1, (1,1)), corner(r2, (1,1))), (corner(r1, (-1,-1)), corner(r2, (-1,-1)))
        end
        (
            rect1=ax2.finallimits[],
            rect2=ax1.finallimits[],
            # slines=isnothing(slines) ? Point2{Float32}[] : Point2{Float32}[slines[1]..., slines[2]...],
            slines=isnothing(slines) ? Point2{Float32}[] : Point2{Float32}[slines...],
        )
    end

    rectattrs = (; strokewidth, strokecolor, color, xautolimits=false, yautolimits=false, rectattrs...)
    p1 = poly!(ax1, (@lift $obs.rect1); rectattrs...)
    p2 = poly!(ax2, (@lift $obs.rect2); rectattrs...)
    translate!(p1, 0, 0, 200)
    translate!(p2, 0, 0, 200)
    plt = linesegments!(pscene, (@lift $obs.slines); color=strokecolor, linewidth=strokewidth, linestyle=:dot, lineattrs...)
    translate!(plt, 0, 0, 200)
    return nothing
end

# Functions to transform data to and from Cartesian
#   x = r cos(θ)
#   y = r sin(θ)
#   r = √(x² + y²)
#   θ = acos(x / r) = asin(x / r)
# To get correlation (R) and STD (σ), we have
#   σ = r
#   R = cos(θ)
# so
#   σ = √(x² + y²)
#   R = cos(θ) = x / √(x² + y²) = x / σ
# and
#   x = R * σ
#   y = √(σ² - x²)
function xy_from_Rσ(R, σ)
    x = R * σ
    Point2(x, sqrt(σ^2 - x^2))
end
function Rσ_from_xy(x, y)
    σ = √(x^2 + y^2)
    Point2(x / σ, σ)
end

# Do the actual plotting now
# First, construct the figure and a polar axis on the first quadrant
fig = Figure(size=(1000, 600))

# Corrticks for Taylor diagram
# corrticks = [-1; -0.99; -0.95; -0.9:0.1:-0.7; -0.6:0.2:0.6; 0.7:0.1:0.9; 0.95; 0.99; 1.0]
# corrticks = [0:0.2:0.6; 0.7:0.1:0.9; 0.95; 0.99; 1.0]
corrticks = [0.2:0.2:0.6; 0.7:0.1:0.9; 0.95; 0.99]
# corrticks = [0.7:0.1:0.9; 0.95; 0.99; 1.0]
function myformat(corrtick; printat=0)
    if corrtick == printat
    # if corrtick == 0.7
        return rich(rich("R", font=:italic), " = $corrtick")
    else
        isinteger(corrtick) && (corrtick = Int(corrtick))
        str = string(corrtick)
        return replace(str, "-" => "−")
    end
end

σref = "σ"
rticks = 0:σr/4:σmax

Rlimits = (0, 1)
thetalimits = acos.(Rlimits[[2,1]]) # zoom in
rlimits = (0, σmax)

thetaticks = (acos.(corrticks), myformat.(corrticks; printat=0))
xlow, ylow = xy_from_Rσ(Rlimits[2], rlimits[1])
xhigh, yhigh = xy_from_Rσ(Rlimits[1], rlimits[2])
axa = Axis(fig[1, 1];
    spinewidth = 1,
    aspect = DataAspect(),
    limits = (0, σmax2, 0, σmax2),
    xgridvisible = false,
    ygridvisible = false,
    topspinevisible = false,
    rightspinevisible = false,
    xtrimspine = true,
    ytrimspine = true,
    xlabel = "STD",
    ylabel = "STD (years)",
    xticks = (rticks, ["0", rich("0.25", σref), rich("0.5", σref), rich("0.75", σref), rich("1", σref), rich("1.25", σref)]),
    yticks = round.(rticks),
)





# Plot polar axis grid
for (θ, text) in zip(thetaticks...)
    lines!(axa, σmax * cos(θ) * [0, 1], σmax * sin(θ) * [0, 1]; color = (:black, 0.2), linewidth = 1)
    lines!(axa, σmax * cos(θ) * [1, 1.01], σmax * sin(θ) * [1, 1.01]; color = :black, linewidth = 1)
    # ticks
    text!(axa, 1.02σmax * cos(θ), 1.02σmax * sin(θ); text, color = :black, rotation = θ, align = (:left, :center))
end

text!(axa, σmax2 * cos(π/4), σmax2 * sin(π/4); text = "correlation", color = :black, rotation = -π/4, align = (:center, :bottom))
# θs = acos(Rlimits[2]):0.01:acos(Rlimits[1])
θs = 0:0.01:π/2
for r in rticks
    lines!(axa, r * cos.(θs), r * sin.(θs); color = (:black, 0.2), linewidth = 1)
end
lines!(axa, σmax * cos.(θs), σmax * sin.(θs); color = :black, linewidth = 1)


# Create isolines of root centered mean squared difference (E′) and labels
levels = (0.25:0.25:4) .* σr
rgrid = (0:0.01:1) .* σmax
θgrid = (0:0.01:π)
E′fun(σf, σr, R) = sqrt(σf^2 + σr^2 - 2 * σf * σr * R)
# E′grid = [sqrt(r^2 + σr^2 - 2 * σr * r * cos(θ)) for θ in θgrid, r in rgrid]
E′grid = [E′fun(r, σr, cos(θ)) for θ in θgrid, r in rgrid]

# labelformatter(E′s) = map(E′ -> rich("$(E′/σr)", rich(" σ", subscript("ref"))), E′s)
# labelformatter(E′s) = map(E′ -> "$(format(round(10E′/σr)/10, stripzeros = true)) $σref", E′s)
function labelformatterfun(E′; printat=0.5)
    if isapprox(E′, printat * σr; rtol = 1e-4)
        # "E′ = 0.5$σref"
        "RMSD = $printat$σref"
    elseif (mod(E′/σr, 1) ≈ 0) || (mod(E′/σr, 1) ≈ 1)
        format(round(E′/σr), stripzeros = true) * "$σref"
    elseif (mod(E′/σr, 0.5) ≈ 0) || (mod(E′/σr, 0.5) ≈ 1)
        format(round(10E′/σr)/10, stripzeros = true) * "$σref"
    else
        format(round(100E′/σr)/100, stripzeros = true) * "$σref"
    end
end
labelformatter1(E′s) = map(labelformatterfun, E′s)
contour!(axa, rgrid' .* cos.(θgrid), rgrid' .* sin.(θgrid), E′grid;
    levels,
    linestyle = :dot,
    labels = true,
    linewidth = 1,
    labelformatter = labelformatter1,
    # color = cgrad(:Archambault, categorical = true)[3],
    color = :black,
    # color = :black,
)


# skill score isolines
R₀ = 0.9940801590730742 # maximum correlation obtainable from ensemble
S(σf, σr, R) = 4 * (1 + R) / ((σf/σr + σr/σf)^2 * (1 + R₀))
# S(σf, σr, R) = 4 * (1 + R)^4 / ((σf/σr + σr/σf)^2 * (1 + R₀)^4)
Sgrid = [S(r, σr, cos(θ)) for θ in θgrid, r in rgrid]
Slevels = [0:0.1:0.9; 0.95; 0.99]
# labelformatter2(v) = map(x -> (x ≈ 0.99) ? " score = $x" : "$x", v)
contour!(axa, rgrid' .* cos.(θgrid), rgrid' .* sin.(θgrid), Sgrid;
    levels = Slevels,
    labels = true,
    linewidth = 1,
    linestyle = :dash,
    # labelformatter = labelformatter2,
    # color = cgrad(:Archambault, categorical = true)[4]
    color = :black
)
# Slevels2 = [0:0.1:0.5; 0.6:0.05:0.85; 0.9:0.01:0.99]
Slevels2 = [0.7:0.05:0.85; 0.9:0.01:0.99]
colorscale = ReversibleScale(x -> 1 - 2acos(x) / π, x -> cos(π/2 * (1 - x)), limits = (0.001, 0.999))
ctrf = contourf!(axa, rgrid' .* cos.(θgrid), rgrid' .* sin.(θgrid), Sgrid;
    levels = Slevels2,
    colormap = :nuuk,
    # colormap = cgrad(cgrad(:watermelon, categorical = true)[5:7], rev = true),
    extendhigh = :auto,
    extendlow = :auto,
    colorscale,
)
translate!(ctrf, 0, 0, -100)


# Now, plot the actual data




Ps = xy_from_Rσ.(Rs, σfs)


# Plot reference (AA age)
xAA, yAA = xy_from_Rσ(1, σr)
# Add lines for AA Interation
# σfsobs = [vals.σf for vals in TDvalsobs]
# Rsobs = [vals.R for vals in TDvalsobs]
# xobs, yobs = collect.(zip(xy_from_Rσ.(Rsobs, σfsobs)...) |> collect)
# scatterlines!(axa, xobs, yobs; color = :black, linewidth = 1, markersize = 3)
# text!(axa, xobs[1], yobs[1]; text = "AA sequence", align = (:left, :bottom), offset = (0, 3), fontsize)
offset = 20
txtline = [offset, 0]
fontsize = 12
lines!(axa, xAA .+ txtline, yAA .+ txtline; linewidth = 1, color = :black)
scatter!(axa, xAA, yAA;
    color = :black,
)
text!(axa, xAA + offset, yAA + offset; text = "AA age", align = (:left, :bottom), offset = (0, 0), fontsize)






# Plot Age with Matt's constants
xMatt, yMatt = xy_from_Rσ(TDvals2.R, TDvals2.σf)
# No need for a cross if there is already a dot
# scatter!(axa, xMatt, yMatt;
#     color = :black,
#     marker = :xcross,
# )
lines!(axa, xMatt .- txtline, yMatt .- txtline; linewidth = 1, color = :black)
text!(axa, xMatt - offset, yMatt - offset; text = "C19", align = (:right, :top), offset = (0, 0), fontsize)
# text!(axa, xMatt - offset, yMatt - offset; text = "diffusivity constants of\nChamberlain et al. (2019)", align = (:right, :top), offset = (0, 0), fontsize)

ikeep = σfs .< σmax
markersize = 5
scatter!(axa, Ps[ikeep];
    color = Ēs,
    colorrange = (-200, 200),
    colormap = cgrad(:tableau_red_blue, rev = true),
    # marker = :cross,
    markersize,
    strokewidth = 2,
    strokecolor = :black,
)
sc = scatter!(axa, Ps[ikeep];
    color = Ēs,
    colorrange = (-200, 200),
    colormap = cgrad(:tableau_red_blue, rev = true),
    highclip = cgrad(:tableau_red_blue, rev = true)[end],
    lowclip =  cgrad(:tableau_red_blue, rev = true)[1],
    # marker = :cross,
    markersize,
)
# scatterlines!(axa, cube_lines;
#     color = (:black, 0.1),
#     markersize = 1,
#     markercolor = :black,
#     linewidth = 1,
# )

# Picking the right diffusivity
# I'm going to base this off the skill score mostly, then look at the bias.
df = DataFrame(
    skillscore = S.(σfs, σr, Rs)[:],
    bias = Ēs[:],
    correlation = Rs[:],
    RMSD = Es[:],
    CRMSD = E′s[:],
    ΔSTD = σfs[:] .- σr,
    κVdeep = κVdeeps,
    κVML = κVMLs,
    κH = κHs,
)
r̄ = TDvals[1].r̄
α = 0.066
sum((df.skillscore .> 0.96) .& (-α * r̄ .< df.bias .< α * r̄) .& (-α * σr .< df.ΔSTD .< α * σr))
iopt = only(findall((df.skillscore .> 0.96) .& (-α * r̄ .< df.bias .< α * r̄) .& (-α * σr .< df.ΔSTD .< α * σr)))
# iopt = only(findall((df.skillscore .> 0.96) .& (-1 .< df.bias .< 1)))
# iopt = only(findall((df.skillscore .> 0.96) .& (-15 .< df.bias .< 15)))
maxbias = 0.05r̄
dfsel = (df.skillscore .> 0.96) .& (-maxbias .< df.bias .< maxbias)
ioptsub = argmax(df.skillscore[dfsel])
iopt = findall(dfsel)[ioptsub]
@show df[iopt, :]
xopt, yopt = xy_from_Rσ(Rs[iopt], σfs[iopt])
scatter!(axa, [xopt], [yopt];
    color = Ēs[iopt],
    colorrange = (-200, 200),
    colormap = cgrad(:tableau_red_blue, rev = true),
    marker = :star5,
    markersize = 2markersize,
    strokewidth = 2,
    strokecolor = :black,
)
txtline = [2offset, 0]
lines!(axa, xopt .+ txtline, yopt .- 0 .* txtline; linewidth = 1, color = :black)
text!(axa, xopt + 2offset, yopt; text = "preferred", align = (:left, :center), offset = (3, 0), fontsize)
scatter!(axa, [xopt], [yopt];
    color = Ēs[iopt],
    colorrange = (-200, 200),
    colormap = cgrad(:tableau_red_blue, rev = true),
    marker = :star5,
    markersize = 2markersize,
)





labeloptions = (
    font = :bold,
    align = (:left, :top),
    # align = (:center, :center),
    offset = (5, -2),
    space = :relative,
    fontsize = 24,
)
txt1 = text!(axa, 0, 1; text = "a", labeloptions..., strokecolor = :white, strokewidth = 3)
txt2 = text!(axa, 0, 1; text = "a", labeloptions...)
translate!(txt1, 0, 0, 100)
translate!(txt2, 0, 0, 100)







Rlimits = (0.92, 0.96)
thetalimits = acos.(Rlimits[[2,1]]) # zoom in
rlimits = (0.82σr, 1.02σr)
corrtickszoom = Rlimits[1]:0.01:Rlimits[2]

# plot zoomed in
xlow, ylow = xy_from_Rσ(Rlimits[2], rlimits[1])
xhigh, yhigh = xy_from_Rσ(Rlimits[1], rlimits[2])
limits = (xlow, xhigh, ylow, yhigh)
axb = Axis(fig[1, 2];
    spinewidth = 1,
    aspect = DataAspect(),
    limits,
    xticksvisible = false,
    yticksvisible = false,
    xticklabelsvisible = false,
    yticklabelsvisible = false,
    xgridvisible = false,
    ygridvisible = false,
)

# Plot polar axis grid
for θ in thetaticks[1]
    ablines!(axb, 0, sin(θ); color = (:black, 0.2), linewidth = 1)
end
# θs = acos(Rlimits[2]):0.01:acos(Rlimits[1])
θs = 0:0.01:π/2
for r in rticks
    lines!(axb, r * cos.(θs), r * sin.(θs); color = (:black, 0.2), linewidth = 1)
end


xgrid = range(xlow, xhigh, length = 200)
ygrid = range(ylow, yhigh, length = 200)
E′gridxy = [E′fun(σf, σr, R) for (R, σf) in Rσ_from_xy.(xgrid, ygrid')]

labelformatter1bis(E′s) = map(x -> labelformatterfun(x; printat=0.25), E′s)
contour!(axb, xgrid, ygrid, E′gridxy;
    levels,
    linestyle = :dot,
    labels = true,
    linewidth = 1,
    labelformatter = labelformatter1bis,
    color = :black,
)

Sgridxy = [S(σf, σr, R) for (R, σf) in Rσ_from_xy.(xgrid, ygrid')]
contour!(axb, xgrid, ygrid, Sgridxy;
    levels = Slevels,
    labels = true,
    linewidth = 1,
    labelformatter = v -> map(x -> "score = " * format(x, conversion = "f", stripzeros = true), v),
    linestyle = :dash,
    color = :black
)
ctrf = contourf!(axb, xgrid, ygrid, Sgridxy;
    levels = Slevels2,
    colormap = :nuuk,
    extendhigh = :auto,
    extendlow = :auto,
    colorscale,
)
translate!(ctrf, 0, 0, -100)

markersize = 10
scatter!(axb, Ps[ikeep];
    color = Ēs,
    colorrange = (-200, 200),
    colormap = cgrad(:tableau_red_blue, rev = true),
    # marker = :cross,
    markersize,
    strokewidth = 2,
    strokecolor = :black,
)
sc = scatter!(axb, Ps[ikeep];
    color = Ēs,
    colorrange = (-200, 200),
    colormap = cgrad(:tableau_red_blue, rev = true),
    highclip = cgrad(:tableau_red_blue, rev = true)[end],
    lowclip =  cgrad(:tableau_red_blue, rev = true)[1],
    # marker = :cross,
    markersize,
)

scatter!(axb, [xopt], [yopt];
    color = Ēs[iopt],
    colorrange = (-200, 200),
    colormap = cgrad(:tableau_red_blue, rev = true),
    marker = :star5,
    markersize = 2markersize,
    strokewidth = 2,
    strokecolor = :black,
)
offset = 10
txtline = [offset, 0]
lines!(axb, xopt .+ txtline, yopt .- 0 .* txtline; linewidth = 1, color = :black)
text!(axb, xopt + offset, yopt; text = "preferred", align = (:left, :center), offset = (3, 0), fontsize)
sc2 = scatter!(axb, [xopt], [yopt];
    color = Ēs[iopt],
    colorrange = (-200, 200),
    colormap = cgrad(:tableau_red_blue, rev = true),
    marker = :star5,
    markersize = 2markersize,
)
translate!(sc2, 0, 0, 100)




cb = Colorbar(fig[2, 1], ctrf;
    vertical = false,
    flipaxis = false,
    tellheight = true,
    ticks = corrticks,
    tickformat = v -> map(x -> format(x, conversion = "f", stripzeros = true), v),
    # ticks = (Slevels2, map(x -> x ∈ corrticks ? string(x) : "", Slevels2)),
    label = rich("skill score"),
)
cb.width = Relative(0.6)

cb = Colorbar(fig[2, 2], sc;
    vertical = false,
    flipaxis = false,
    tellheight = true,
    # ticks = (, cbarticklabelformat.(levels)),
    label = rich("bias (years)"),
)
cb.width = Relative(0.6)


labeloptions = (
    font = :bold,
    align = (:left, :top),
    offset = (5, -2),
    space = :relative,
    fontsize = 24,
)
txt1 = text!(axb, 0, 1; text = "b", labeloptions..., strokecolor = :white, strokewidth = 3)
txt2 = text!(axb, 0, 1; text = "b", labeloptions...)
translate!(txt1, 0, 0, 100)
translate!(txt2, 0, 0, 100)

rectattrs = (strokecolor = (:black, 0.2), linestyle = :solid)
lineattrs = (color = (:black, 0.2), linestyle = :solid)
zoom_lines!(axa, axb; rectattrs, lineattrs)

outputfile = joinpath(outputdir, "Taylor_diagram_3D_v3_xy.png")
@info "Saving image file:\n  $(outputfile)"
save(outputfile, fig)
outputfile = joinpath(outputdir, "Taylor_diagram_3D_v3_xy.pdf")
@info "Saving image file:\n  $(outputfile)"
save(outputfile, fig)




# fig = Figure(size=(1200, 1200))

# options = (
#     xlabel = "κVML",
#     xticks = (1:length(κVMLs), format.(κVMLs, conversion = "g")),
#     ylabel = "κH",
#     yticks = (1:length(κHs), format.(κHs, conversion = "e")),
# )

# colorscale = ReversibleScale(x -> 1 - 2acos(x) / π, x -> cos(π/2 * (1 - x)), limits = (0.001, 0.999))



# for irow in 1:5
#     # Correlation
#     ax = Axis(fig[irow, 1]; options...)
#     kwargs = (
#         colorrange = (0.5, 1),
#         colormap = :viridis,
#         colorscale,
#     )
#     hm1 = heatmap!(ax, Rs[irow, :, :]; kwargs...)
#     (irow ≠ 5) && hidexdecorations!(ax)

#     # STD - STDref
#     ax = Axis(fig[irow, 2]; options...)
#     kwargs = (
#         colorrange = (-σr/2, σr/2),
#         colormap = cgrad(:RdBu, rev = true),
#     )
#     hm2 = heatmap!(ax, σfs[irow, :, :] .- σr; kwargs...)
#     hideydecorations!(ax)
#     (irow ≠ 5) && hidexdecorations!(ax)


#     # RMS
#     ax = Axis(fig[irow, 3]; options...)
#     kwargs = (
#         colorrange = (0, σr/2),
#         colormap = :plasma,
#     )
#     hm3 = heatmap!(ax, E′s[irow, :, :]; kwargs...)
#     hideydecorations!(ax)
#     (irow ≠ 5) && hidexdecorations!(ax)

#     # bias
#     ax = Axis(fig[irow, 4]; options...)
#     kwargs = (
#         colorrange = (-200, 200),
#         colormap = cgrad(:tableau_red_blue, rev = true),
#     )
#     hm4 = heatmap!(ax, Ēs[irow, :, :]; kwargs...)
#     hideydecorations!(ax)
#     (irow ≠ 5) && hidexdecorations!(ax)

#     # score
#     ax = Axis(fig[irow, 5]; options...)
#     kwargs = (
#         colorrange = (0.9, 0.97),
#         colormap = :cividis,
#         colorscale,
#     )
#     hm5 = heatmap!(ax, S.(σfs, σr, Rs)[irow, :, :]; kwargs...)
#     hideydecorations!(ax)
#     (irow ≠ 5) && hidexdecorations!(ax)

#     Label(fig[irow, 0]; text = "κVDeep = " * format.(κVdeeps[irow], conversion = "e"), rotation = π/2, tellheight = false)

#     if irow == 5
#         cbopt = (
#             vertical = false,
#             flipaxis = false,
#             tellheight = true,
#             width = Relative(0.85),
#         )
#         cb = Colorbar(fig[6, 1], hm1; cbopt..., label = rich("correlation"), ticks = [0, 0.5, 0.8, 0.95, 1])
#         cb = Colorbar(fig[6, 2], hm2; cbopt..., label = rich("STD - STDref"))
#         cb = Colorbar(fig[6, 3], hm3; cbopt..., label = rich("RMS"))
#         cb = Colorbar(fig[6, 4], hm4; cbopt..., label = rich("bias"))
#         cb = Colorbar(fig[6, 5], hm5; cbopt..., label = rich("skill score"), ticks = [0.9, 0.95, 0.97])
#     end

# end



# Label(fig[0,1]; text = "correlation", tellwidth = false)
# Label(fig[0,2]; text = "STD - STDref", tellwidth = false)
# Label(fig[0,3]; text = "RMS", tellwidth = false)
# Label(fig[0,4]; text = "bias", tellwidth = false)
# Label(fig[0,5]; text = "skill score", tellwidth = false)

# outputfile = joinpath(outputdir, "Taylor_diagram_3D_heatmaps_v2.png")
# @info "Saving image file:\n  $(outputfile)"
# save(outputfile, fig)






# @show sort(df, :skillscore, rev = true)
# @show sort(df, :RMSD)



