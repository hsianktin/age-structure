# simple finite element method for solving
#    PDE: ∂ₜn = - (μ + ∫ₓkndx′)n - ∂ₓn 
using CSV,DataFrames

β_0 = 2
μ_0 = 0.4
k_0 = 2
type = "semi-x′"
if length(ARGS) == 3
    β_0 = parse(Float64, ARGS[1])
    μ_0 = parse(Float64, ARGS[2])
    k_0 = parse(Float64, ARGS[3])
    type = "semi-constant" # default type of kernel
elseif length(ARGS) == 4
    β_0 = parse(Float64, ARGS[1])
    μ_0 = parse(Float64, ARGS[2])
    k_0 = parse(Float64, ARGS[3])
    type = ARGS[4]
else
    error!("usage: age-structure-base.jl β_0 μ_0 k_0")
end

# step size for age (x)
dx = 0.02

# truncation of age. In the usual node, n(x) represents density.
#     However, in the last node, its value n(x=A) = ∫ₓndx' for x > A.
A = 10 # this selection is special, due to the specific form

# birth rate, a function.
β(x) = β_0 # constant birth rate

# death rate μ
μ(x) = μ_0 # constant death rate

# the interaction k
K = A # cut off age for interaction, not used in fact
function k(x′,x, type="semi-constant")
    if type == "semi-constant"
        # semi-constant specification
        a = 2 # critical age of interaction 
        # k(x',x) = k_0 if x' > a > x, 0 otherwise
        if 0 ≤ x < a < x′ ≤ K 
            return k_0
        else
            return 0
        end
    elseif type == "semi-x′"
        # semi-constant specification
        a = 2 # critical age of interaction 
        # k(x',x) = k_0 if x' > a > x, 0 otherwise
        if 0 ≤ x < a < x′ ≤ K 
            return k_0 * x′
        else
            return 0
        end
    end
end
## import the necessary functions
include("age-structure-utils.jl")

## Total population N from distribution n
N(n) = ∑(n[1:end-1])*dx+n[end]

## function to obtain equilibrium density via time-forward method
function n₊(β, μ, k)
    Nₜ = Array{Float64,1}()
    Tₜ = Array{Float64,1}()
    t = 0
    # set initial condition, a proper initial condition helps
    # to converge faster
    n = ones(1+convert(Int,A/dx))
    for i in 2:convert(Int,A/dx)+1
        n[i] = n[i-1] * exp(-μ(x(i-1))*dx)
    end
    # run
    # count = 0
    while t < 200
        push!(Nₜ,N(n))
        push!(Tₜ,t)
        # the hard upper limit t = 1000 is necessary because the actual solution
        # might be oscillating... In practice it is often the case
        # n,t = ∂ₜndt(n,t,β,μ,k)
        n,t = euler(n,t,β,μ,k)
        # count += 1
        # if count % 10 == 0
        #     print("t = $(t), ∫βndx/dx = $(n[1])\r")
        # end
    end
    # CSV.write("traces/$(β_0)_$(μ_0)_$(k_0).csv", DataFrame(N=Nₜ,T=Tₜ))
    return (n,t)
end
# obtain the current value
(n,t) = n₊(β, μ, k);
# obtain the perturbed value
# δμ(x) = μ_0 + 0.01
# (δn,δt) = n₊(β, δμ, k)

# # inspect the equilibrium density
# using PGFPlotsX
# @pgf axis = SemiLogYAxis(
#     {
#         xlabel = "age \$x\$",
#         ylabel = "density \$n(x)\$",
#     },
# )

# push!(axis,    Plot(Coordinates([x(i) for i in 1:i(A)-1],n[1:i(A)-1])))
# push!(axis,LegendEntry("numerical equilibrium density"))
# plot_entry = Plot(Coordinates([x(i) for i in 1:i(A)-1],n̄(n)))
# legend_entry = LegendEntry("estimated equilibrium density")
# push!(axis,plot_entry)
# push!(axis,legend_entry)
# axis
# μ̃(n,3)

# save data
df = DataFrame(
    β_0 = Float64[],
    μ_0 = Float64[],
    k_0 = Float64[],
    N = Float64[],
    T = Float64[],
    # ∂ᵤN = Float64[],
    # ∂ᵤlnN = Float64[],
)
push!(df,[
        β_0,
        μ_0,
        k_0,
        ∑(n[1:end-1])*dx+n[end],
        t,
        # (∑(δn[1:end-1])*dx+δn[end]-∑(n[1:end-1])*dx-n[end])/0.01,
        # (∑(δn[1:end-1])*dx+δn[end]-∑(n[1:end-1])*dx-n[end])/(0.01*∑(n[1:end-1])*dx+n[end])
    ])
if !isdir("data")
    mkdir("data")
end

CSV.write("data/$(rand()).csv",df)
