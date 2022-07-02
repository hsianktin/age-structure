βs = [i/10 for i in 5:1:7]
μs = [i/10 for i in 2:2:10]
ks = [1]
type = "semi-x′"
label = "test_1_b"
cmds = []
for β in βs, μ in μs, k in ks
    cmd = `julia age-structure-base.jl $(β) $(μ) $(k) $(type)`
    push!(cmds,cmd)
end
if isfile("fig/$(label).csv") && !overwrite
    println("data exists, skipping calculation...")
    df = CSV.read("fig/$(label).csv",DataFrame)
else
    println("starting calculation...")
    @showprogress pmap(run,cmds)
    println("calculation completed, saving data...")
    df = DataFrame(
        β_0 = Float64[],
        μ_0 = Float64[],
        k_0 = Float64[],
        N = Float64[],
        Λ = Float64[],
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