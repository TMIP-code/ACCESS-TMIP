# This script just defines some simple functions to fetch the GM terms for different models
# At this stage it's a bespoke script tailored to the few ACCESS model members that have GM data
# stored on Gadi

# function ϕGM(model, experiment, member, time_window)

#     if (model == "ACCESS-ESM1-5") && (experiment == "historical") && (member == "r5i1p1f1")
#         return ϕGMfromACCESSESM15(member, time_window)
#     elseif model == "ACCESS-CM2"
#         return ϕGMfromACCESSCM2(experiment, member, time_window)
#     elseif model == "ACCESS1-3"
#         return ϕGMfromACCESS13(experiment, member, time_window)
#     elseif model == "ACCESS1-0"
#         return ϕGMfromACCESS10(experiment, member, time_window)
#     else
#         error("$model $experiment $member not supported")
#     end

# end

function CMIP6member2CSIROmember(member)
    # Convert to the member number to be used in the file path by CSIRO
    # I.e., r1i1p1f1 -> HI-05 (so add 4 to the r number and format it as a 2 digit number)
    member_r = parse(Int, split(member, "r")[2][1])
    CSIRO_member = "HI-$(format(member_r + 4, width=2, zeropadding=true))"
    return CSIRO_member
end

function parse_time_window(time_window)
    # parse the time_window string to get the years
    yearstart = parse(Int, split(time_window, "-")[1][4:end])
    yearend = parse(Int, split(time_window, "-")[2][4:end])
    return yearstart, yearend
end

function timeaverage(files, varname)
    FILLVALUE = ncgetatt(first(files), varname, "_FillValue")
    varmonthly = replace(ncread(first(files), varname) |> Array .|> Float64, FILLVALUE => 0.0)
    weights = diff(ncread(first(files), "time_bounds"), dims=1)
    cumweights = sum(weights)
    varmean = sum(weights) * dropdims(nanmean(varmonthly, reshape(weights, (1, 1, 1, length(weights))), dims=4), dims=4)
    for filepath in files[2:end]
        weights = diff(ncread(filepath, "time_bounds"), dims=1)
        cumweights += sum(weights)
        varmonthly .= replace(ncread(filepath, varname) |> Array .|> Float64, FILLVALUE => 0.0)
        varmean .+= sum(weights) * dropdims(nanmean(varmonthly, reshape(weights, (1, 1, 1, length(weights))), dims=4), dims=4)
    end
    varmean ./= cumweights
end


function ϕfromACCESSESM15(gm_or_submeso, member, time_window)

    @assert gm_or_submeso ∈ ("gm", "submeso")

    # Convert to the member number to be used in the file path by CSIRO
    CSIRO_member = CMIP6member2CSIROmember(member)

    # parse the time_window string to get the years
    yearstart, yearend = parse_time_window(time_window)

    # Load the GM terms for ACCESS-ESM1-5
    files = ["/g/data/p73/archive/CMIP6/ACCESS-ESM1-5/$CSIRO_member/history/ocn/ocean_month.nc-$(year)1231" for year in yearstart:yearend]
    @info "Loading GM velocities from ACCESS-ESM1-5 data"
    @info "  u GM"
    ψᵢmean = timeaverage(files, "tx_trans_$gm_or_submeso")
    @info "  v GM"
    ψⱼmean = timeaverage(files, "ty_trans_$gm_or_submeso")

    # These are transport diagnostics at the bottom of the face cells, so need to diff to get mass transport
    @info "  Taking vertical diff"
    (nx, ny, _) = size(ψᵢmean)
    ϕᵢmean = diff([fill(0.0, nx, ny, 1);;; ψᵢmean], dims=3)
    ϕⱼmean = diff([fill(0.0, nx, ny, 1);;; ψⱼmean], dims=3)

    return ϕᵢmean, ϕⱼmean

end





function ϕfromACCESSESM15_Tilo(gm_or_submeso, member, time_window)

    @assert gm_or_submeso ∈ ("gm", "submeso")

    # Convert to the member number to be used in the file path by CSIRO
    CSIRO_member = CMIP6member2CSIROmember(member)

    # parse the time_window string to get the years
    yearstart, yearend = parse_time_window(time_window)

    # Load the GM terms for ACCESS-ESM1-5
    files = ["/g/data/p73/archive/CMIP6/ACCESS-ESM1-5/$CSIRO_member/history/ocn/ocean_month.nc-$(year)1231" for year in yearstart:yearend]
    @info "Loading GM velocities from ACCESS-ESM1-5 data"
    @info "  u GM"
    ψᵢmean = timeaverage(files, "tx_trans_$gm_or_submeso")
    @info "  v GM"
    ψⱼmean = timeaverage(files, "ty_trans_$gm_or_submeso")

    # These are transport diagnostics at the bottom of the face cells, so need to diff to get mass transport
    @info "  Taking vertical diff"
    (nx, ny, _) = size(ψᵢmean)
    ϕᵢmean = diff([fill(0.0, nx, ny, 1);;; ψᵢmean], dims=3)
    ϕⱼmean = diff([fill(0.0, nx, ny, 1);;; ψⱼmean], dims=3)

    return ϕᵢmean, ϕⱼmean

end
