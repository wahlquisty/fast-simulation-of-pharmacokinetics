# ForwardDiff with closure for fast PK simulation as well as update using x(k+1) = Ad*xk + Bd*uk
# Date: 221017

using BenchmarkTools, ControlSystemsBase, ForwardDiff

include("simulation.jl") # fcts to simulate state and compute output
include("pkmodels.jl") # calculations of λ, R from parameter vector θ
include("../expm.jl")
# Since we cannot differentiate through VectorizationBase.vexp, we cannot use SLEEFPirates.exp. Instead, use expm.julia. Adds approximately 1 us in time per patient.
# Diagonal of discrete time system matrix
@inline @fastmath function getΦdiag(λ, h)
    return @. exp(λ * h) 
end

# Simulate with closure
function simulate(θ, infusionrate, bolusdose, h, youts, order)
    function simulate_inner(θ::Array{T}) where {T<:Real}
        model = PK(θ, order)
        j = 1 # counter to keep track of next free spot in y
        x = @SVector zeros(eltype(u), 3) # initial state
        y = zeros(T,length(youts))
        for i in eachindex(u, hs, v)
            if i in youts # if we want to compute output
                x, yi = @inbounds updatestateoutput(x, hs[i], model.V1inv, model.λ, model.λinv, model.R, u[i], v[i]) # update state and compute output
                y[j] = yi
                j += 1
            else
                x = @inbounds updatestate(x, hs[i], model.λ, model.λinv, u[i], v[i]) # update state
            end
        end
        return y
    end
    return simulate_inner
end

# Simulation in discrete time, x(k+1) = Ad*x(k) + Bd*u(k)
function simAdBd(θ, infusionrate, bolusdose, hs, youts)
    function simAdBd_inner(θ::Array{T}) where {T<:Real}
        k10, k12, k13, k21, k31, V1 = θ
        A = [-(k10 + k12 + k13) k12 k13
            k21 -k21 0.0
            k31 0.0 -k31]

        B = [1 / V1; 0; 0]
        C = [1 0 0]
        D = 0

        PK = ss(A, B, C, D)

        nt = length(hs)
        Cp_s = Array{T}(zeros(nt))
        x = [0, 0, 0] # Initial state vector

        for i = 1:nt-1
            if !iszero(bolusdose[i])
                x = x + PK.B * bolusdose[i] # bolus
                Cp_s[i] = x[1]
            end
            if hs[i] == 0.0
                Cp_s[i+1] = Cp_s[i]
                continue
            end
            PK_d = c2d(PK, hs[i])
            x = PK_d.A * x + PK_d.B * infusionrate[i] # infusion
            Cp_s[i+1] = x[1]
        end
        return Cp_s[youts]
    end
    return simAdBd_inner
end

# Example
θ = rand(Float32, 6); # model parameters, k10, k12, k13, k21, k31, V1
nt = 100 # length of time series
u = rand(Float32, nt); # infusionrates
v = rand(Float32, nt); # bolusdoses
hs = sort(rand(Float32, nt)); # sampling times
youts = unique(sort(Int32.(floor.(rand(10,) * 10) .+ 1))) # observation indices


# jacobian for new simulation
simulate_inner = simulate(θ, u, v, hs, youts,3)
simulate_inner(θ)
ForwardDiff.jacobian(simulate_inner, θ) # works! (without SLEEFPirates)

# jacobian for the other simulation
simAdBd_inner = simAdBd(θ, u, v, hs, youts)
simAdBd_inner(θ)
ForwardDiff.jacobian(simAdBd_inner, θ) # works!

# btiming
@btime ForwardDiff.jacobian(simulate_inner, $θ) # 757 allocations, 30 us
@btime ForwardDiff.jacobian(simAdBd_inner, $θ) # 6400 allocations, 1.2 ms
# comparison: old simulation is x slower, x more allocations
