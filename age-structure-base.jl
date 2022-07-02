# simple finite element method for solving
#    PDE: ∂ₜn = - (μ + ∫ₓkndx′)n - ∂ₓn 
using CSV,DataFrames
using LinearAlgebra
using Roots
β_0 = 2.5
μ_0 = 0.4
k_0 = 1
type = "semi-(x′-x)"
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
    println(ARGS)
else
    error("usage: age-structure-base.jl β_0 μ_0 k_0")
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
K = A - dx # cut off age for interaction, not used in fact
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
    elseif type == "semi-(x′-x)"
        # semi-constant specification
        a = 2 # critical age of interaction 
        # k(x',x) = k_0 if x' > a > x, 0 otherwise
        if 0 ≤ x < a < x′ ≤ K 
            return k_0 * (x′-x)
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

function n⁻(β, μ, k, nᵩ)
    n = [nᵩ for i in 1:i(A)] # density n[i]=n(x=i*dx)
    # for the last node, μ(A) * n[i(A)] = n[i(A)-1], with n[i(A)] = nᵩ
    n[i(A) - 1] = μ(A) * nᵩ
    for j in i(A)-1:-1:2 # backwards iteration
        # (μ(x(j)) + ∫ₓkndx′(k,n,x(j))) * n[j] = (n[j-1] - n[j])/dx
        n[j-1] = n[j] + ((μ(x(j)) + ∫ₓkndx′(k,n,x(j))) * n[j]) * dx
    end
    # our loss function is ∂ₜlog(n(0)) 
    ℓ = (∫βndx(β,n)/dx # birth 
            - (μ(x(1)) + ∫ₓkndx′(k,n,x(1))) * n[1] # death
            - n[1]/dx)/n[1] # age transition 
    return (n,ℓ)
end

function n⁻(β, μ, k)
    ℓ(nᵩ) = n⁻(β, μ, k, nᵩ)[2]
    lower_bound = -10.0
    upper_bound = 2.0
    it_counter = 0
    while isnan(ℓ(exp(lower_bound))) && it_counter < 1e2
        it_counter += 1
        lower_bound += 0.1
    end
    while ℓ(exp(lower_bound)) < 0 && it_counter < 1e2
        lower_bound += 0.1
        lower_bound -= 0.1
    end
    while isnan(ℓ(exp(upper_bound))) && it_counter < 1e2
        it_counter += 1
        upper_bound -= 0.1
    end
    while ℓ(exp(upper_bound)) > 0 && it_counter < 1e3
        it_counter += 1
        upper_bound += 0.1
    end
    if it_counter == 1e2
        println("iteration limit reached")
        return [0.0 for i in 1:i(A)]
    end
    if isnan(ℓ(exp(lower_bound)))
        error("lower bound is nan")
    end
    if isnan(ℓ(exp(upper_bound)))
        error("upper bound is nan")
    end
    println("ℓ_ = $(ℓ(exp(lower_bound)))")
    println("ℓ^ = $(ℓ(exp(upper_bound)))")
    nᵩ = exp(find_zero(ln_a -> ℓ(exp(ln_a)), (lower_bound, upper_bound), Bisection() ))
    return n⁻(β, μ, k, nᵩ)[1]
end

function λ₀(β, μ, k, n)
    J = zeros(i(A),i(A))
    e = [[convert(Float64,j==k) for j in 1:i(A)] for k in 1:i(A)]
    for j in 1:i(A)
        J[:,j] .= δ∂ₜn(β, μ, k,n,e[j])
    end
    return maximum(real.(eigvals(J)))
end

# obtain the current value
@time n = n⁻(β, μ, k);
@time λ = λ₀(β, μ, k, n);
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
    Λ = Float64[],
    # ∂ᵤN = Float64[],
    # ∂ᵤlnN = Float64[],
)
push!(df,[
        β_0,
        μ_0,
        k_0,
        N(n),
        λ,
        # (∑(δn[1:end-1])*dx+δn[end]-∑(n[1:end-1])*dx-n[end])/0.01,
        # (∑(δn[1:end-1])*dx+δn[end]-∑(n[1:end-1])*dx-n[end])/(0.01*∑(n[1:end-1])*dx+n[end])
    ])
if !isdir("data")
    mkdir("data")
end

CSV.write("data/$(rand()).csv",df)
