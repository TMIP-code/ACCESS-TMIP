# using Pkg
# Pkg.activate(; temp=true)
using GLMakie
using LinearAlgebra
using SparseArrays
using Unitful
using Unitful: m, km, s, minute, hr, d, yr
using Memoize

let
    κVML = 0.1m^2/s
    κVdeep = 3e-5m^2/s
    a = (100km)^2
    h = 10m
    vol = a * h

    TVML = ustrip(s^-1, κVML * a / (h * vol))
    TVdeep = ustrip(s^-1, κVdeep * a / (h * vol))
    N = 20
    z = range(start = h/2, length = N, step = h)
    MLD(t) = 40.1m - 33m * sin(2π * t / ustrip(s, 1yr))
    iMLD(t) = searchsortedfirst(z, MLD(t))
    offdiag(iMLD) = [TVML * (i < iMLD) + TVdeep * (iMLD ≤ i) for i in 1:N-1]
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

    τ₀ = ustrip(s, 1s)


    Ω = [i == 1 for i in 1:N]

    M(iMLD) = T(iMLD) + Diagonal(Ω) / τ₀

    M̃(iMLD) = V⁻¹ * M(iMLD)' * V



    function Ã₊(iMLD, δt)
        # println("memoizing")
        lu(I + δt * M̃(iMLD))
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
    uδts = Any[1yr//12, 1yr//4]

    fig = Figure(size=(1000, 200 * length(uδts)))

    for (iδt, uδt) in enumerate(uδts)
        @show δt = ustrip(s, uδt)

        ts = 0yr:uδt:Tsim
        tskip = uδt < 1d ? 1d ÷ uδt : 1
        # tskip = uδt < 1hr ? 1hr ÷ uδt : 1
        # tskip = uδt < 1minute ? 1minute ÷ uδt : 1
        ts_saves = ts[1:tskip:end]

        maxisim = length(ts) - 1
        EBg = Array{Float64}(undef, length(g0), length(ts_saves))
        ∫ℊ̃ = Array{Float64}(undef, length(g0), length(ts_saves))
        EBtmp = Vector{Float64}(undef, length(g0))
        ∫ℊ̃tmp = Vector{Float64}(undef, length(g0))
        # CNg = Array{Float64}(undef, length(g0), length(ts_saves))
        # CNtmp = Vector{Float64}(undef, length(g0))
        # EFg = Array{Float64}(undef, length(g0), maxisim + 1)

        EBtmp[:] .= g0
        ∫ℊ̃tmp[:] .= 0
        # CNtmp[:] .= g0
        EBg[:, 1] .= g0
        ∫ℊ̃[:, 1] .= 0
        # CNg[:, 1] .= g0
        # EFg[:] .= g0

        # FIXME: loop over ts instead and write every tskip with a counter
        # for i in eachindex(ts_saves)[1:end-1]
        #     for iskip in 1:tskip
        #         t = ustrip(s, ts_saves[i] + iskip * uδt)
        #         EBtmp .= Ã₊(t, δt) \ EBtmp
        #         # CNtmp .= Ã₊2(t, δt) \ (Ã₋2(t, δt) * CNtmp)
        #         # EFg[:, i + 1] .= Ã₋(t, δt) * EFg[:, i]
        #     end
        #     EBg[:, i + 1] .= EBtmp
        #     # CNg[:, i + 1] .= CNtmp
        # end
        isav = 1
        for (i, tu) in enumerate(ts[1:end-1])
            t = ustrip(s, tu)
            ∫ℊ̃tmp .+= EBtmp * δt
            ldiv!(Ã₊(iMLD(t), δt), EBtmp)
            if mod(i, tskip) == 0
                isav += 1
                EBg[:, isav] .= EBtmp
                ∫ℊ̃[:, isav] .= ∫ℊ̃tmp
            end
        end

        # limits = (0, ustrip(yr, min(5yr, Tsim)), nothing, nothing)
        limits = (0, ustrip(yr, Tsim), nothing, nothing)
        ax = Axis(fig[iδt, 1]; xlabel = "time (years)", ylabel = "ℊ̃ (δt = $uδt)", limits, yscale = Makie.pseudolog10)
        # ax2 = Axis(fig[iδt, 1]; xlabel = "time (years)", ylabel = "ℊ̃ (δt = $uδt)", limits, yaxisposition = :right, yticklabelcolor = :red, )
        izs = 2:length(z)
        dodge = false
        colors = collect(cgrad(:viridis, length(z), categorical = true, rev = true))
        for iz in izs
            # lines!(ax, ustrip.(yr, ts), EBg[iz, :])
            # scatterlines!(ax, ustrip.(yr, ts), EBg[iz, :], markersize = 10 * (iz .≤ iMLD.(ustrip.(s, ts))))
            color = colors[iz]
            scatterlines!(ax, ustrip.(yr, ts_saves), EBg[iz, :] ./ maximum(EBg[iz, :]) .+ dodge * 0.02iz; color, markersize = 0 * (iz .≤ iMLD.(ustrip.(s, ts_saves))))
            # scatterlines!(ax, ustrip.(yr, ts_saves), EBg[iz, :]; markersize = 7 * (iz .≤ iMLD.(ustrip.(s, ts_saves))), color)
            # scatter!(ax2, ustrip.(yr, ts), iMLD.(ustrip.(s, ts)), color = :red)
        end
        hidexdecorations!(ax, label = iδt < length(uδts), ticklabels = iδt < length(uδts), ticks = iδt < length(uδts), grid = false)
        # hidexdecorations!(ax2, label = iδt < length(uδts), ticklabels = iδt < length(uδts), ticks = iδt < length(uδts), grid = false)
        ax = Axis(fig[iδt, 2]; xlabel = "time (years)", ylabel = "1 − ℰ = ∫ℊ̃dt (δt = $uδt)", limits)
        for iz in izs
            color = colors[iz]
            scatterlines!(ax, ustrip.(yr, ts_saves), ∫ℊ̃[iz, :] .+ dodge * 0.02iz; color, markersize = 0 * (iz .≤ iMLD.(ustrip.(s, ts_saves))))
        end
        hidexdecorations!(ax, label = iδt < length(uδts), ticklabels = iδt < length(uδts), ticks = iδt < length(uδts), grid = false)
        # xlims!(ax, 0, ustrip(yr, 5yr))
        # for (iy, y) in enumerate(eachslice(EBg, dims = 1)[2:end])
        # end
        # for (iy, y) in enumerate(eachslice(CNg, dims = 1)[2:end])
        #     lines!(ax, ustrip.(yr, 0yr:uδt:Tsim), y, alpha = 0.5)
        # end
        # lines!(ax, ustrip.(yr, 0yr:uδt:Tsim), cumsum(CNg[end, :] * δt |> Vector); color = Cycled(2), label = "CN", alpha = 0.5)
        # lines!(ax, ustrip.(yr, 0yr:uδt:Tsim), cumsum(EFg[end, :] * δt |> Vector); color = Cycled(3), label = "EF", linewidth = 2)
        # ylims!(ax, -0.02, 1.02)
    end

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
    Ã₊⁻¹[abs.(Ã₊⁻¹) .< 1e-4] .= 0
    sparse(Ã₊⁻¹)
    Ã₊⁻¹ * ones(N)
    M(4)
    fig
end


