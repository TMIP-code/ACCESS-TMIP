using Pkg
Pkg.activate(".")
Pkg.instantiate()

using DataFrames

# Directory for data
if gethostname() == "benoits-MacBook-Pro.local"
    TMIPDIR = "/Users/benoitpasquier/Data/TMIP"
else # on Gadi.
    TMIPDIR = "/scratch/xv83/TMIP"
end
DATADIR = joinpath(TMIPDIR, "data")

# Go through all the directories to check which models are available
REQUIREDVARIABLES = ["umo", "vmo", "mlotst", "volcello", "areacello"]
hasrequireddata(dirpath, var) = isfile(joinpath(dirpath, "$var.nc"))
hasrequireddata(dirpath) = all(var -> hasrequireddata(dirpath, var), REQUIREDVARIABLES)

df = DataFrame(model=String[], experiment=String[], member=String[], time_window=String[])

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
                if hasrequireddata(inputdir)
                    push!(df, (model, experiment, member, time_window))
                end
            end
        end
    end
end

@show df

