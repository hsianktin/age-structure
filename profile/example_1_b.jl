βs = [i/10 for i in 15:5:25]
μs = [i/10 for i in 5:5:15]
ks = [1]
type = "semi-x′"
label = "example_1_b"
cmds = []
for β in βs, μ in μs, k in ks
    cmd = `julia age-structure-base.jl $(β) $(μ) $(k) $(type)`
    push!(cmds,cmd)
end
if isfile("fig/$(label).csv") && !overwrite
    println("data exists, skipping simulation...")
    df = CSV.read("fig/$(label).csv",DataFrame)
else
    println("starting simulation...")
    @showprogress pmap(run,cmds)
    println("simulation completed, saving data...")
    df = DataFrame(
        β_0 = Float64[],
        μ_0 = Float64[],
        k_0 = Float64[],
        N = Float64[],
        T = Float64[],
        # ∂ᵤN = Float64[],
        # ∂ᵤlnN = Float64[],
    )
    for f in readdir("./data/",join=true)
        temp_df = CSV.read(f,DataFrame)
        for i in 1:length(temp_df[:,1])
            push!(df,temp_df[i,:])
        end
    end
    sort!(df,[:β_0,:μ_0])
    CSV.write("fig/$(label).csv",df)
    for f in readdir("./data/",join=true)
        rm(f)
    end
end
println("creating TeX figures...")
include("../age-structure-plot-bm.jl")
println("done.")