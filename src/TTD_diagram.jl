# # using Pkg
# # Pkg.activate(; temp=true)
# # using GLMakie
# using CairoMakie
# using LinearAlgebra
# using SparseArrays
# using Unitful
# using Unitful: m, km, s, minute, hr, d, yr
# using Statistics
# using Format

# # κVML = 0.1m^2/s
# κVML = 0.001m^2/s
# # κVdeep = 1e-5m^2/s
# κVdeep = 1e-4m^2/s
# a = (100km)^2
# h = 10m
# vol = a * h

# TVML = ustrip(s^-1, κVML * a / (h * vol))
# TVdeep = ustrip(s^-1, κVdeep * a / (h * vol))

# N = 20

# z = range(start = h/2, length = N, step = h)
# MLD(t) = 100.1m - 33.3m * cos(π/4 + 2π * t / ustrip(s, 1yr))
# iMLD(t) = searchsortedfirst(z, MLD(t))
# offdiag(iMLD) = [TVML * (i < iMLD) + TVdeep * (iMLD ≤ i) for i in 1:N-1]
# ondiag(iMLD) = [0; offdiag(iMLD)] + [offdiag(iMLD); 0]
# # ΩMLD(iMLD) = [i ≤ iMLD for i in 1:N]
# T(iMLD) = Tridiagonal(
#     -offdiag(iMLD),
#     ondiag(iMLD),
#     -offdiag(iMLD)
# )
# v = fill(vol, N)
# V = Diagonal(v)
# V⁻¹ = Diagonal(1 ./ v)
# # T = [Tcoeff -Tcoeff; -Tcoeff Tcoeff]

# τ₀ = ustrip(s, 1yr)


# Ω = [i == 1 for i in 1:N]

# M(iMLD) = T(iMLD) + Diagonal(Ω) / τ₀

# M̃(iMLD) = V⁻¹ * M(iMLD)' * V

# function Ã₊(iMLD, δt)
#     # println("memoizing")
#     lu(I + δt * M̃(iMLD))
# end
# # Ã₋(iMLD, δt) = I - δt * M̃(iMLD)
# # Ã₊2(iMLD, δt) = I + δt/2 * M̃(iMLD)
# # Ã₋2(iMLD, δt) = I - δt/2 * M̃(iMLD)

# g0 = Ω / τ₀


# Tsim = 13yr

# resolution = 48

# uδt = 1yr/resolution
# δt = ustrip(s, uδt)

# ts = 0yr:uδt:Tsim + 1yr

# maxisim = length(ts) - 1
# 𝒢̃ = Array{Float64}(undef, length(g0), length(ts), resolution)
# 𝒢̃tmp = Vector{Float64}(undef, length(g0))

# for itf in 1:resolution
#     tf = ustrip(s, itf * uδt)
#     𝒢̃[:, 1, itf] .= g0
#     𝒢̃tmp[:] .= g0


#     for (i, tu) in enumerate(ts[1:end-1])
#         t = ustrip(s, tu)
#         ldiv!(Ã₊(iMLD(tf - t), δt), 𝒢̃tmp)
#         𝒢̃[:, i + 1, itf] .= 𝒢̃tmp
#     end
# end

# iz = length(z) - 1

# Stich together the 𝒢̃ to make a map
Nmap𝒢̃ = length(ts)
# map𝒢̃ = zeros(Nmap𝒢̃, Nmap𝒢̃)
map𝒢̃ = fill(NaN, Nmap𝒢̃, Nmap𝒢̃)
for i in axes(map𝒢̃, 1)
    for f in axes(map𝒢̃, 2)
        if i ≤ f
            map𝒢̃[i, f] = 𝒢̃[iz, f - i + 1, mod1(f, resolution)]
        end
    end
end


fig = Figure(size = (400, 450))
𝓉 = rich("t", font = :italic)
𝓉i = rich(𝓉, subscript("i", offset = (0.1, 0)))
𝓉f = rich(𝓉, subscript("f", offset = (0.1, 0)))

xticks = 0:1:ustrip(yr, Tsim + 1yr)
yticks = 0:1:ustrip(yr, Tsim + 1yr)

Δt = ustrip(yr, Tsim / 2 + 1yr)
limits = (0, Δt + 1 + 1, ustrip(yr, Tsim) - Δt - 1, ustrip(yr, Tsim + 1yr))

ax = Axis(
    fig[1, 1];
    xlabel = rich("injection time, ", 𝓉i, " (yr)"),
    ylabel = rich("reemergence time, ", 𝓉f, " (yr)"),
    # limits = (nothing, nothing, nothing, nothing),
    limits,
    xticks,
    yticks,
    # xticklabelsvisible = false,
    # xticksvisible = false,
    # yticklabelsvisible = false,
    # yticksvisible = false,
    aspect = DataAspect(),
    backgroundcolor = :lightgray,
)

ctrf = contourf!(
    ax, ustrip.(yr, ts), ustrip.(yr, ts), map𝒢̃;
    colormap = cgrad(:devon, rev = true),
)
# ctrf = heatmap!(ax, ustrip.(yr, ts), ustrip.(yr, ts), map𝒢̃)

# move the contourf under the grid
translate!(ctrf, 0, 0, -100)

# Draw the window of interest
# lines!(ax, collect(limits)[[1, 2, 2, 1, 1]], collect(limits)[[3, 3, 4, 4, 2]], color = :blue)

# draw ti = tf line
# ablines!(ax, 0, 1, color = :black)
brokenabline = [t0 - 10, t0 - 1.4, NaN, t0 - 0.6, t0 + 1.6, NaN, t0 + 2.4, t0 + 10]
lines!(ax, brokenabline, brokenabline, color = :black)
# lines!(ax, [, t0 - 0.7], [t0 - 1.3, t0 - 0.7], color = :white)
text!(ax, t0 - 1, t0 - 1, text = rich(𝓉i, " = ", 𝓉f), rotation = π / 4, align = (:center, :center), color = :black)
# lines!(ax, [t0 + 1.7, t0 + 2.3], [t0 + 1.7, t0 + 2.3], color = :white)
text!(ax, t0 + 2, t0 + 2, text = "τ = 0", rotation = π / 4, align = (:center, :center), color = :black)

# colors of polygons and labels
colors = cgrad(:Egypt, categorical = true)[[3, 1]]

# text for tf = ti + τ
text!(ax, t0 + 0.5 - 0.2, t0 + 0.5 + τ + 0.2, text = rich(𝓉f, " = ", 𝓉i, " + τ"), rotation = π / 4, align = (:center, :center), color = colors[1])
text!(ax, t0 + 0.5 - 0.2 - τ, t0 + 0.5 + 0.2, text = rich(𝓉i, " = ", 𝓉f, " − τ"), rotation = π / 4, align = (:center, :center), color = colors[2])

# Draw
τ = Δt - 2.1
t0 = Δt - 0.7
# scatter!(ax, [t0], [t0], color=:red)
# Add lines to delineate the matching patches
for τ2 in 0:floor(τ)
    lines!(ax, t0 .+ [1 - min(1, τ - τ2), 1], t0 + 1 .+ [τ2, τ2], color = colors[1])
    lines!(ax, t0 .- [τ2, τ2], t0 .+ [0, min(1, τ - τ2)], color = colors[2])
    # lines!(ax, t0 .+ [0, 1], t0 .+ [-τ2, -τ2], color = :black, linestyle = :dash)
end
redpoly = (t0 .+ [0, 1, 1, 0, 0], t0 .+ [0, 1, 1 + τ, τ, 0])
poly!(ax, redpoly...; color = (colors[1], 0.1))
lines!(ax, redpoly..., color = colors[1])
bluepoly = (t0 .+ [0, 1, 1 - τ, -τ, 0], t0 .+ [0, 1, 1, 0, 0])
poly!(ax, bluepoly...; color = (colors[2], 0.1))
lines!(ax, bluepoly...; color = colors[2])

# Add line for τ
# bracket!(ax, t0 - τ, t0, t0, t0, text = "τ", offset = 2, color = colors[2], textcolor = colors[2], orientation = :down)
offset = 0.3
# bracket!(ax, t0 + 1 + offset, t0 + 1, t0 + 1 + offset, t0 + 1 + τ, text = rich("τ = ", 𝓉f, " − ", 𝓉i), offset = 2, color = colors[1], textcolor = colors[1], orientation = :down)
# lines!(ax, [t0 + 1, t0 + 1 + offset, NaN, t0 + 1, t0 + 1 + offset], [t0 + 1, t0 + 1, NaN, t0 + 1 + τ, t0 + 1 + τ], color = colors[1], linestyle = :dot)
bracket!(ax, t0, t0 - offset, t0 + 1, t0 - offset, text = "1 yr", offset = 2, color = :black, textcolor = :black, orientation = :down)
lines!(ax, [t0, t0, NaN, t0 + 1, t0 + 1], [t0 - offset, t0, NaN, t0 - offset, t0 + 1], color = :black, linestyle = :dot)

# Equality
𝒢̃str = rich("𝒢", superscript("~", offset = (-0.5, 0.2)))
𝒢̃funstr = rich(𝒢̃str, "(", 𝓉i, ",", 𝓉f, ")")
# intstr(sub, sup) = rich("∫", subsup(sub, sup))
intstr(sub, sup, offsub, offsup) = rich("∫", subscript(sub, offset = (-offsub, -1)), superscript(sup, offset = (-offsup, 1)))
text = rich(intstr("0", "1", 0.6, 0.5), "d", 𝓉i, "  ", intstr(𝓉i, rich(𝓉i, " + τ"), 0.6, 1.2), "d", 𝓉f, " ", 𝒢̃funstr)
text!(ax, t0 - 0.2, t0 + 4; text, align = (:right, :center), color = colors[1], fontsize = 20)
text = rich(intstr("0", "1", 0.6, 0.5), "d", 𝓉f, "  ", intstr(rich(𝓉f, " − τ"), 𝓉f, 1.2, 1.2), " d", 𝓉i, " ", 𝒢̃funstr)
text!(ax, t0 - 2.5, t0 + 1.2; text, align = (:center, :bottom), color = colors[2], fontsize = 20)

# Colorbar
cb = Colorbar(
    fig[2, 1], ctrf;
    label = rich("Cyclostationary ", 𝒢̃str, "(", 𝓉i, ",", 𝓉f, ")"),
    # label = L"$\tilde{\mathcal{G}}(t_\mathrm{i}, t_\mathrm{f})$",
    width = Relative(3 / 4),
    vertical = false,
    flipaxis = false,
    ticklabelsvisible = false,
    ticksvisible = false,
)

colgap!(fig.layout, 10)
rowgap!(fig.layout, 10)
# Label(fig[0, 1]; text = "")


save("TTD_diagram.pdf", fig)
fig
