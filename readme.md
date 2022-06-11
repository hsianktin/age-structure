# Environment setup
- Install [The Julia Programming Language](https://julialang.org/), make sure `julia` is added to PATH.
- In the julialang REPL, run 
```julia
using Pkg;Pkg.add(["CSV","DataFrames","PGFPlotsX","ProgressMeter"]
```
- To compile the tex into svg plots, need full LaTeX installation in PATH.

# Instructions
The code is an implementation of sampling the positive steady state with constant $\beta$, $\mu$, and two-compartment asymmetric interaction kernel $k$. It tries to vary different $\beta$ and $\mu$, and determine if overcompensation exists.
- In the shell environment (Powershell/Cmd/Bash), run
```shell
julia age-structure-run.jl "Profile Name" [overwrite]
```
- It will read the profile with the name `Profile Name.jl` under the directory `profile/`.
- In the profile, one specifies the necessary arguments needed for simulation, and provide them in the array of commands `cmds`. It also provides a string `label` and the method for storing data.
	- A sample profile `general.jl` is given as follows.
	- To run it, use `julia age-structure-run.jl general`.
```julia
βs = [i/10 for i in 1:20]
μs = [i/10 for i in 1:0.5:10]
ks = [1]
label = "general"
cmds = []
for β in βs, μ in μs, k in ks
    cmd = `julia age-structure-base.jl $(β) $(μ) $(k)`
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
        ∂ᵤN = Float64[],
        ∂ᵤlnN = Float64[],
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
```
-  `age-structure-utils.jl` contains the general implementation of functions.
- `age-structure-base.jl` receives the commands and run the simulation.
	- It also handles the parameters passed by commandline arguments.
	- To use different types of dynamics, please change the way $\beta$, $\mu$, and $k$ are defined through parameters. 

# Multithreading
By default, in `age-structure-run.jl`, it uses 8 processors by `addprocs(8)`. Change `8` to another positive integer according to your needs and your machine.