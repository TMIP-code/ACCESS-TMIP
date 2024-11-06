using Format

model = "ACCESS-ESM1-5"

members = 1:40

# decades = 1850:10:2010
decades = [1990]
time_window = "Jan1990-Dec1999"
# time_window = "Jan1071-Dec1100" # <- last 30 years of ACCESS-ESM1-5 piControl
# time_window = "Jan1420-Dec1449" # <- last 30 years of ACCESS-CM2 piControl

experiment = "historical"
# experiment = "piControl"



scratchdir = "/scratch/xv83/TMIP/data"
@assert isdir(scratchdir)

datadir = "/g/data/xv83/TMIP/data"
@assert isdir(datadir)

# Check scracth dir for decadally averaged CMIP6 output of transport
modeldir = joinpath(scratchdir, model)
for member in members
    CMIP6_member = "r$(member)i1p1f1"
    memberdir = joinpath(modeldir, CMIP6_member)
    tx_trans_gm = [isfile(joinpath(memberdir, "month_tx_trans_gm_$(decade)s.nc")) for decade in decades]
    ty_trans_gm = [isfile(joinpath(memberdir, "month_ty_trans_gm_$(decade)s.nc")) for decade in decades]
    tx_trans_submeso = [isfile(joinpath(memberdir, "month_tx_trans_submeso_$(decade)s.nc")) for decade in decades]
    ty_trans_submeso = [isfile(joinpath(memberdir, "month_ty_trans_submeso_$(decade)s.nc")) for decade in decades]
    membersdone = tx_trans_gm .& ty_trans_gm .& tx_trans_submeso .& ty_trans_submeso
    numTODO = sum(.!membersdone)
    decadesTODO = decades[.!membersdone]
    numdecadeTODOmax = 10
    postfix = reduce((x,y) -> x * " " * y, string.(decadesTODO[1:min(end,numdecadeTODOmax)]), init="")
    postfix = all(membersdone) ? "" : all(.!membersdone) ? " missing all!" : " missing:" * postfix * (numTODO > numdecadeTODOmax ? "..." : "")
    prefix = all(membersdone) ? "DONE" : "TODO"
    println("$prefix $CMIP6_member = r$(member)i1f1p1$postfix")
end





# Gadi directory for input files
inputdirfun(member) = "/scratch/xv83/TMIP/data/$model/$experiment/$member/$(time_window)"

# find all members for which the inputdir contains umo.nc, vmo.nc, mlotst.nc, volcello.nc, and areacello.nc
requiredvariables = ["umo", "vmo", "mlotst", "volcello", "areacello"]
hasrequireddata(member, variable_name) = isfile(joinpath(inputdirfun(member), "$variable_name.nc"))
hasrequireddata(member) = all(variable_name -> hasrequireddata(member, variable_name), requiredvariables)
members = readdir("/scratch/xv83/TMIP/data/$model/$experiment")

# sort members by r, i, p[, f]

member_regex = CMIP_version == "CMIP6" ? r"r(\d+)i(\d+)p(\d+)f(\d+)" : r"r(\d+)i(\d+)p(\d+)"
parse_member(member) = parse.(Int, match(member_regex, member).captures)
members = sort(members, by = x -> parse_member(x))
dataavailability = DataFrame(
    :member => members,
    :has_it_all => hasrequireddata.(members),
    [Symbol(var) => [hasrequireddata(member, var) for member in members] for var in requiredvariables]...,
)
show(dataavailability; allrows = true)

all(.!dataavailability.has_it_all) && @warn "Nothing to do: missing something for all members"


# Check gdata for decadally averaged GM outputs
modeldir = joinpath(datadir, model)
for member in members
    CSIRO_member = "HI-$(format(member + 4, width=2, zeropadding=true))"
    memberdir = joinpath(modeldir, CSIRO_member)
    tx_trans_gm = [isfile(joinpath(memberdir, "month_tx_trans_gm_$(decade)s.nc")) for decade in decades]
    ty_trans_gm = [isfile(joinpath(memberdir, "month_ty_trans_gm_$(decade)s.nc")) for decade in decades]
    tx_trans_submeso = [isfile(joinpath(memberdir, "month_tx_trans_submeso_$(decade)s.nc")) for decade in decades]
    ty_trans_submeso = [isfile(joinpath(memberdir, "month_ty_trans_submeso_$(decade)s.nc")) for decade in decades]
    membersdone = tx_trans_gm .& ty_trans_gm .& tx_trans_submeso .& ty_trans_submeso
    numTODO = sum(.!membersdone)
    decadesTODO = decades[.!membersdone]
    numdecadeTODOmax = 10
    postfix = reduce((x,y) -> x * " " * y, string.(decadesTODO[1:min(end,numdecadeTODOmax)]), init="")
    postfix = all(membersdone) ? "" : all(.!membersdone) ? " missing all!" : " missing:" * postfix * (numTODO > numdecadeTODOmax ? "..." : "")
    prefix = all(membersdone) ? "DONE" : "TODO"
    println("$prefix $CSIRO_member = r$(member)i1f1p1$postfix")
end