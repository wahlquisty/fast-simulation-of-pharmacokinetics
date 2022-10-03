# Simulation of three compartment mamillary model for full state output (full x)
# Random simulation data.
# Date: 220929

# computes all states (full R)
using Pkg
Pkg.activate("..")

using ControlSystems, StaticArrays, LinearAlgebra
using BenchmarkTools

include("pksimulation_fullstate.jl") # fcts to simulate state and compute output
include("pkmodels_fullstate.jl") # calculations of λ, R from parameter vector θ

# Example
θ = rand(Float32, 6); # model parameters, k10, k12, k13, k21, k31, V1
nt = 100 # length of time series
u = 10 * rand(Float32, nt); # infusionrates
v = 10 * rand(Float32, nt); # bolusdoses
hs = sort(rand(Float32, nt)); # sampling times
youts = unique(sort(Int32.(floor.(rand(10,) * nt) .+ 1))) # observation indices

# Simulation
y = zeros(Float32, length(youts))
PKsim!(y, θ, u, v, hs, youts) # simulation, not interested in states
y, states = PKsim!(y, θ, u, v, hs, youts, xout=true) # simulation, gets all states


# btiming
@btime PKsim!($y, $θ, $u, $v, $hs, $youts, xout=false) # 1.4 us, 0 allocations
@btime PKsim!($y, $θ, $u, $v, $hs, $youts, xout=true) # 1.6 us, 0 allocations