## utils
i(x) = round(Int,x/dx) + 1 # x to index
x(i) = (i-1)*dx # index to x
∑(x) = sum(x)
# calculate the interaction experienced at age x=a
# under current density profile n and interaction k
# here k is a function of x′ and x
function ∫ₓkndx′(k::Function,n::Array,a::Number)
    ∫ₓkndx = 0
    for j in 1:i(A)-1
        ∫ₓkndx += k(x(j),a,type) * n[j] * dx
    end
    ∫ₓkndx += k(A,a) * n[i(A)]
    return ∫ₓkndx
end

# birth function
# here β is a function β(x)
function ∫βndx(β::Function,n::Array)
    ∫βnda = 0
    for j in 1:i(A)-1
        ∫βnda += β(x(j)) * n[j] * dx
    end
    ∫βnda += β(A) * n[i(A)]
    return ∫βnda
end

# obtain ∂ₜn
# here β, μ, k are functions
# n is the current density profile
function ∂ₜn(β::Function,μ::Function,k::Function,n::Array)
    V = zeros(1+convert(Int,A/dx))
    V[1] = (∫βndx(β,n)/dx # birth 
            - (μ(x(1)) + ∫ₓkndx′(k,n,x(1))) * n[1] # death
             - n[1]/dx) # age transition 
    for j in 2:i(A)-1
        V[j] = (- (μ(x(j)) + ∫ₓkndx′(k,n,x(j))) * n[j] # death
         + (n[j-1] - n[j])/dx) # age transition
    end
    V[i(A)] = (- μ(A) * n[i(A)] # death 
                + (n[i(A)-1])) # age transition
                # note that n[i(A)-1] = n(A-dx) represents the density
                # however, n[i(A)] = ∫ₓndx' for x > A represents the total number
                # actual number of people at n[i(A)-1] is n(A-dx)dx.
                # dx cancels the rate 1/dx
    return V
end

## main method
function euler(n::Array,t::Number,β::Function,μ::Function,k::Function)
	# adaptive version of euler method
	dndt = ∂ₜn(β, μ, k, n)
    # the step size should not change 1% of the current n(x) ∀x
	norm_dt = maximum(abs.(dndt./n))
	ε = 0.01/norm_dt
    nₜ₊₁ = n + ε * dndt
    return (nₜ₊₁, t + ε)
    # push!(N,∑(n)*dx)
    # push!(T,t)
end

function ∂ₜndt(n::Array,t::Number,β::Function,μ::Function,k::Function)
	# euler method with fixed step size
	dndt = ∂ₜn(β, μ, k, n)
	ε = 1/500;
    nₜ₊₁ = n + ε * dndt
    return (nₜ₊₁, t + ε)
    # push!(N,∑(n)*dx)
    # push!(T,t)
end


## quality control
N(n) = ∑(n[1:i(A)-1])*dx + n[i(A)]
μ̃(n,x) = μ(x) + ∫ₓkndx′(k,n,x)
function n̄(n)
    n̄ = ones(i(A)-1)
    n̄[1] = n[1]
    for i in 2:i(A)-1
        n̄[i] = n̄[i-1] * exp(-μ̃(n,x(i-1))*dx)
    end
    return n̄
end
function loss(n)
    return ∑(abs.(n̄(n)-n[1:i(A)-1]) ./ n[1:i(A)-1])/(i(A)-1)
end


## deprecated
function euler!(β,μ,k,ε)
    global n
    global t
    n += ε * ∂ₜn(β, μ, k, n)
    t += ε
    # push!(N,∑(n)*dx)
    # push!(T,t)
end

###### deprecated ################################################
### specific function
### assumption:
###   1. k(x′,x) = const. for x' > a > x, 0 otherwise.
###   2. β, μ are constants.
# function n₋(β_0, μ_0, k_0, a)
#     dx = 0.1 
#     A = 40 
#     β(x) = β_0
#     K = 42
#     μ(x) = μ_0  
#     a = 10
#     function k(x′,x)
#         # interaction
#         if x′ > a && x' < K && x < a
#             return k_0
#         else
#             return 0
#         end
#     end
#     n = ones(1+convert(Int,A/dx))
#     t = 0
#     # set initial condition
#     for i in 2:convert(Int,A/dx)+1
#         n[i] = n[i-1] * exp(-μ(x(i-1))*dx)
#     end
#     # run
#     while maximum(∂ₜn(β,μ,k,n)./n) > 1e-2 && t < 100
#         euler!(β,μ,k)
#     end
#     return n
# end

# function ∂ᵤN(n)
#     μ̃ = μ_0 + ∫ₓkndx′(k,n,0)
#     nₐ = n[i(a)]
#     ∂ᵤnₐ = (μ_0/k_0) * ((-exp(-μ̃*a)/μ_0^2)/(a*exp(-μ̃*a)*(1/μ_0 - 1/μ̃) + (1-exp(-μ̃*a))/μ̃^2) - 1 + k_0*nₐ/μ_0^2)
#     ∂ᵤμ̃ = 1 - k_0*nₐ/μ_0^2 + (k_0/μ_0) * ∂ᵤnₐ
#     return (exp(μ̃*a)/β_0) * ∂ᵤnₐ + (a*nₐ*exp(μ̃*a)/β_0) * ∂ᵤμ̃
# end

###################################################################