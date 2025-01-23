using GLMakie
vertices = [
    0.0 0.0;
    1.0 0.0;
    1.0 1.0;
    0.0 1.0;
]

fig = Figure()


colors = [:red, :green, :blue, :orange]

ax = Axis(fig[1, 1])
faces = [
    1 2 3;
    3 4 1;
]
mesh!(ax, vertices, faces, color = colors, shading = NoShading)

ax = Axis(fig[2, 1])
faces = [
    4 1 2;
    2 3 4;
]
mesh!(ax, vertices, faces, color = colors, shading = NoShading)

fig