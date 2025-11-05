using Pkg
Pkg.activate(".")
Pkg.instantiate()

using DataFrames
using PrettyTables

# Directory for data
if gethostname() == "benoits-MacBook-Pro.local"
    TMIPDIR = "/Users/benoitpasquier/Data/TMIP"
else # on Gadi.
    TMIPDIR = "/scratch/xv83/TMIP"
end
DATADIR = joinpath(TMIPDIR, "data")

# Go through all the directories to check which models are available
MASSTRANSPORT_VARS = ["umo", "vmo"]
VELOCITY_VARS = ["uo", "vo"]
MLD_VARS = ["mlotst"]
GRID_VARS = ["areacello", "volcello"]
TEMPSALT_VARS = ["thetao", "so"]
hasrequiredvariable(dirpath, var) = isfile(joinpath(dirpath, "$var.nc"))
hasallrequiredvariables(dirpath, vars) = all(hasrequiredvariable(dirpath, v) for v in vars)

df = DataFrame(
    model = String[],
    experiment = String[],
    member = String[],
    time_window = String[],
    ϕ = Union{String, Missing}[],
    u = Union{String, Missing}[],
    mld = Union{String, Missing}[],
    grid = Union{String, Missing}[],
    TS = Union{String, Missing}[]
)

# files to ignore
ignore(f) = f ∈ [".DS_Store"]

models = [m for m in readdir(DATADIR) if !ignore(m)]

for model in models
    modeldir = joinpath(DATADIR, model)
    experiments = [e for e in readdir(modeldir) if !ignore(e)]
    for experiment in experiments
        experimentdir = joinpath(modeldir, experiment)
        members = [m for m in readdir(experimentdir) if !ignore(m)]
        for member in members
            memberdir = joinpath(experimentdir, member)
            time_windows = [t for t in readdir(memberdir) if !ignore(t)]
            for time_window in time_windows
                inputdir = joinpath(memberdir, time_window)
                row = (
                    model = model,
                    experiment = experiment,
                    member = member,
                    time_window = time_window,
                    ϕ = hasallrequiredvariables(inputdir, MASSTRANSPORT_VARS) ? "ϕᵢ + ϕⱼ" : missing,
                    u = hasallrequiredvariables(inputdir, VELOCITY_VARS) ? "uᵢ + uⱼ" : missing,
                    mld = hasallrequiredvariables(inputdir, MLD_VARS) ? "MLD" : missing,
                    grid = hasallrequiredvariables(inputdir, GRID_VARS) ? "grid" : missing,
                    TS = hasallrequiredvariables(inputdir, TEMPSALT_VARS) ? "T + S" : missing,
                )
                # Don't print if required data is missing
                ismissing(row.grid) || ismissing(row.mld) || (ismissing(row.ϕ) && ismissing(row.u)) && continue
                push!(df, row)
            end
        end
    end
end

hnotOK = Highlighter(
    (data, i, j) -> true,
    bold = false,
    foreground = :light_gray
)
hOK = Highlighter(
    (data, i, j) -> !(ismissing(data[i, 8]) || ismissing(data[i, 7]) || (ismissing(data[i, 5]) && ismissing(data[i, 6]))),
    bold = true,
    foreground = :green
)
hOK2 = Highlighter(
    (data, i, j) -> !(ismissing(data[i, 8]) || ismissing(data[i, 7]) || ismissing(data[i, 5]) || ismissing(data[i, 6])),
    bold = true,
    foreground = :blue
)
hOK3 = Highlighter(
    (data, i, j) -> !(ismissing(data[i, 8]) || ismissing(data[i, 7]) || ismissing(data[i, 5]) || ismissing(data[i, 6]) || ismissing(data[i, 9])),
    bold = true,
    foreground = :magenta
)

pretty_table(df, highlighters = (hOK3, hOK2, hOK, hnotOK), formatters = ft_nomissing, crop = :none)
