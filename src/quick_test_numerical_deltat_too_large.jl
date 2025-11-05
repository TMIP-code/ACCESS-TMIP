# using Pkg
# Pkg.activate(; temp=true)
using GLMakie
using LinearAlgebra
using SparseArrays
using Unitful
using Unitful: m, km, s, minute, hr, d, yr
using Statistics
using Format

κVML = 0.1m^2 / s
κVdeep = 1.0e-5m^2 / s # TODO rerun with this value (1e-5)
a = (100km)^2
h = 10m
vol = a * h

TVML = ustrip(s^-1, κVML * a / (h * vol))
TVdeep = ustrip(s^-1, κVdeep * a / (h * vol))
N = 50
z = range(start = h / 2, length = N, step = h)
MLD(t) = 100.1m - 33.3m * cos(π / 4 + 2π * t / ustrip(s, 1yr))
iMLD(t) = searchsortedfirst(z, MLD(t))
offdiag(iMLD) = [TVML * (i < iMLD) + TVdeep * (iMLD ≤ i) for i in 1:(N - 1)]
ondiag(iMLD) = [0; offdiag(iMLD)] + [offdiag(iMLD); 0]
# ΩMLD(iMLD) = [i ≤ iMLD for i in 1:N]
T(iMLD) = Tridiagonal(
    -offdiag(iMLD),
    ondiag(iMLD),
    -offdiag(iMLD)
)
v = fill(vol, N)
V = Diagonal(v)
V⁻¹ = Diagonal(1 ./ v)
# T = [Tcoeff -Tcoeff; -Tcoeff Tcoeff]

τ₀ = ustrip(s, 1yr)


Ω = [i == 1 for i in 1:N]

M(iMLD) = T(iMLD) + Diagonal(Ω) / τ₀

M̃(iMLD) = V⁻¹ * M(iMLD)' * V


function Ã₊(iMLD, δt)
    # println("memoizing")
    return lu(I + δt * M̃(iMLD))
end
# Ã₋(iMLD, δt) = I - δt * M̃(iMLD)
# Ã₊2(iMLD, δt) = I + δt/2 * M̃(iMLD)
# Ã₋2(iMLD, δt) = I - δt/2 * M̃(iMLD)


g0 = Ω / τ₀

# Tsim = 0.5yr
# uδts = Any[10s, 1hr, 1d, 1yr/12 |> d]
# Tsim = 10yr
# uδts = Any[1hr, 1d, 1yr/12 |> d, 1yr]
# Tsim = 10yr
# uδts = Any[1minute, 1hr, 1d, 30d]

Tsim = 50yr

# uδts = Any[1minute, 1hr, 1d]
resolutions = [12, 4]
# resolutions = [4, 2, 1]
resolutions = [48, 12, 4, 2, 1]
resolutions = [12]
uδts = Any[1yr / res for res in resolutions]

fig = Figure(size = (1000, 400 * length(uδts)))

for (iresolution, resolution) in enumerate(resolutions)
    uδt = uδts[iresolution]
    δt = ustrip(s, uδt)

    ts = 0yr:uδt:(Tsim + 1yr)
    tsav = ts[1:resolution:end]


    maxisim = length(ts) - 1
    @show maxisim, maxisim ÷ resolution
    ℊ̃HF = Array{Float64}(undef, length(g0), length(ts), resolution)
    ℊ̃ = Array{Float64}(undef, length(g0), maxisim ÷ resolution + 1, resolution)
    ℊ̃tmp = Vector{Float64}(undef, length(g0))
    # CNg = Array{Float64}(undef, length(g0), length(ts))
    # CNtmp = Vector{Float64}(undef, length(g0))
    # EFg = Array{Float64}(undef, length(g0), maxisim + 1)

    for itf in 1:resolution
        tf = ustrip(s, itf * uδt)
        ℊ̃[:, 1, itf] .= g0
        ℊ̃HF[:, 1, itf] .= g0
        ℊ̃tmp[:] .= g0

        isav = 1

        for (i, tu) in enumerate(ts[1:(end - 1)])
            t = ustrip(s, tu)
            ldiv!(Ã₊(iMLD(tf - t), δt), ℊ̃tmp)
            if mod(i, resolution) == 0
                ℊ̃[:, isav + 1, itf] .= ℊ̃tmp
                isav += 1
            end
            ℊ̃HF[:, i + 1, itf] .= ℊ̃tmp
        end
    end

    # Stich together the ℊ̃ from each final month
    meanℊ̃ = dropdims(mean(ℊ̃, dims = 3), dims = 3)
    meanℊ̃HF = dropdims(mean(ℊ̃HF, dims = 3), dims = 3)

    # limits = (0, ustrip(yr, min(5yr, Tsim)), nothing, nothing)
    limits = (0, ustrip(yr, Tsim), 0, nothing)
    δtstr = "δt = 1/$(resolution) yr"
    ax = Axis(fig[iresolution, 1]; xlabel = "time (years)", ylabel = "meanℊ̃", limits)
    # ax2 = Axis(fig[iresolution, 1]; xlabel = "time (years)", ylabel = "ℊ̃ (δt = $uδt)", limits, yaxisposition = :right, yticklabelcolor = :red, )
    # izs = 2:length(z)
    izs = length(z)
    dodge = false
    colors = collect(cgrad(:viridis, length(z), categorical = true, rev = true))
    for iz in izs
        # lines!(ax, ustrip.(yr, ts), ℊ̃[iz, :])
        # scatterlines!(ax, ustrip.(yr, ts), ℊ̃[iz, :], markersize = 10 * (iz .≤ iMLD.(ustrip.(s, ts))))
        color = colors[iz]
        # scatterlines!(ax, ustrip.(yr, ts), meanℊ̃HF[iz, :] ./ maximum(meanℊ̃HF[iz, :]) .+ dodge * 0.02iz; color, markersize = 0 * (iz .≤ iMLD.(ustrip.(s, ts))))
        # scatterlines!(ax, ustrip.(yr, tsav), meanℊ̃[iz, :] ./ maximum(meanℊ̃[iz, :]) .+ dodge * 0.02iz; color, markersize = 0 * (iz .≤ iMLD.(ustrip.(s, tsav))), linestyle = :dash)
        scl1 = scatterlines!(ax, ustrip.(yr, ts), meanℊ̃HF[iz, :] .+ dodge * 0.02iz; color = :black, markersize = 10, label = "monthly values")
        # scatterlines!(ax, ustrip.(yr, tsav), meanℊ̃[iz, :] .+ dodge * 0.02iz; color, markersize = 0 * (iz .≤ iMLD.(ustrip.(s, tsav))), linestyle = :dash)
        scl2 = scatterlines!(ax, ustrip.(yr, tsav), meanℊ̃[iz, :] .+ dodge * 0.02iz; color = :red, markersize = 10, linestyle = :dash, label = "yearly saves")
        # scatterlines!(ax, ustrip.(yr, ts), ℊ̃[iz, :]; markersize = 7 * (iz .≤ iMLD.(ustrip.(s, ts))), color)
        # scatter!(ax2, ustrip.(yr, ts), iMLD.(ustrip.(s, ts)), color = :red)
        axislegend(ax)
    end
    hidexdecorations!(ax, label = iresolution < length(uδts), ticklabels = iresolution < length(uδts), ticks = iresolution < length(uδts), grid = false)
    # hidexdecorations!(ax2, label = iresolution < length(uδts), ticklabels = iresolution < length(uδts), ticks = iresolution < length(uδts), grid = false)
    limits = (0, ustrip(yr, Tsim), 0, nothing)
    ax = Axis(fig[iresolution, 2]; xlabel = "time (years)", ylabel = "ℰ = 1 − ∫ℊ̃dt", limits)
    for iz in izs
        color = colors[iz]
        scatterlines!(ax, ustrip.(yr, ts), 1 .- cumsum(meanℊ̃HF[iz, :] * δt) .+ dodge * 0.02iz; color = :black, markersize = 10)
        scatterlines!(ax, ustrip.(yr, tsav), 1 .- cumsum(meanℊ̃[iz, :] * ustrip(s, 1yr)) .+ dodge * 0.02iz; color = :red, markersize = 10, linestyle = :dash)
    end
    hidexdecorations!(ax, label = iresolution < length(uδts), ticklabels = iresolution < length(uδts), ticks = iresolution < length(uδts), grid = false)
    # xlims!(ax, 0, ustrip(yr, 5yr))
    # for (iy, y) in enumerate(eachslice(ℊ̃, dims = 1)[2:end])
    # end
    # for (iy, y) in enumerate(eachslice(CNg, dims = 1)[2:end])
    #     lines!(ax, ustrip.(yr, 0yr:uδt:Tsim), y, alpha = 0.5)
    # end
    # lines!(ax, ustrip.(yr, 0yr:uδt:Tsim), cumsum(CNg[end, :] * δt |> Vector); color = Cycled(2), label = "CN", alpha = 0.5)
    # lines!(ax, ustrip.(yr, 0yr:uδt:Tsim), cumsum(EFg[end, :] * δt |> Vector); color = Cycled(3), label = "EF", linewidth = 2)
    # ylims!(ax, -0.02, 1.02)
    Label(fig[iresolution, 0], δtstr; tellheight = false, fontsize = 14, rotation = π / 2)
end
tfstr = rich(rich("t", font = :italic), subscript("f"))
Label(fig[0, 1], rich("ℊ̃ (averaged over the ", tfstr, ")"); tellwidth = false, fontsize = 14)
Label(fig[0, 2], "ℰ = 1 − ∫ℊ̃dt"; tellwidth = false, fontsize = 14)

# for t in ustrip.(s, range(0yr, stop = Tsim, step = 1yr/12))
#     @show ustrip(yr, t * s) * 12
#     @show MLD(t)
#     @show iMLD(t)
# end
ones(N)'T(ustrip(s, 0.1yr))
TVML
T(4)
Ã₊(4, ustrip(s, uδts[1]))
Ã₊⁻¹ = inv(Tridiagonal(I + ustrip(s, uδts[1]) * M̃(4)))
Ã₊⁻¹[abs.(Ã₊⁻¹) .< 1.0e-4] .= 0
sparse(Ã₊⁻¹)
Ã₊⁻¹ * ones(N)
M(4)
fig
