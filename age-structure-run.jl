using CSV,DataFrames
using Distributed
addprocs(8)
begin
    using ProgressMeter
end
#################
# testing purpose
profile = "example_1_b"
overwrite = false

#################
if length(ARGS) == 1
    profile = ARGS[1]
    overwrite = false
elseif length(ARGS) == 2
    profile = ARGS[1]
    if ARGS[2] == "overwrite"
        overwrite = true
    else
        overwrite = false
    end
else
    error("usage: age-structure-run.jl profile [overwrite]")
end
println("loading profile $(profile)...")
if isfile("profile/$(profile).jl")
    include("profile/$(profile).jl")
else
    error("profile/$(profile).jl not found")
end
