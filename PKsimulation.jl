# Date: 220629
# Here, we assume that we are only looking for the concentration of the central compartment volume

using StaticArrays, SLEEFPirates

####

# update state for bolus
@inline function bolus(x, v)
    return x .+ v
end

# update state with input step
@inline function step(x, λ, Φdiag, u)
    Γdiag = @. 1 / λ * (Φdiag - 1)
    return @. Φdiag * x + Γdiag * u
end

# output computation
function gety(x, V1, R)
    return (1/V1*R'*x)[1]
end

# Diagonal of discrete time system matrix
@inline @fastmath function getΦdiag(λ, h)
    Φdiag = @. SLEEFPirates.exp(λ * h) # cannot differentiate through VectorizationBase.vexp
    # return @. exp(λ * h) 
    return Φdiag
end

# Initiate/update state
@inline @fastmath function updatestate(x, h, λ, u=0.0, v=0.0)
    Φdiag = getΦdiag(λ, h) # compute Φ
    x = bolus(x, v) # update state for bolus
    x = step(x, λ, Φdiag, u) # infusion affect next sample
    return x
end

# Update state and compute output
@inline @fastmath function updatestateoutput(x, h, V1, λ, R, u=0.0, v=0.0)
    Φdiag = getΦdiag(λ, h) # compute Φ
    x = bolus(x, v) # update state for bolus
    y = gety(x, V1, R) # compute output
    x = step(x, λ, Φdiag, u) # infusion affect next sample
    return x, y
end

# Fast PK simulation for 3 compartment model. Modifies first argument, output y (prealloated)
function PKsim!(y, θ, u, v, hs, youts)
    λ, R = update(θ) # Setting up simulator
    j = 1 # counter to keep track of next free spot in y
    x = @SVector [0.0f0, 0.0f0, 0.0f0] # initial state
    for i in eachindex(u)
        if i in youts # if we want to compute output
            x, yi = updatestateoutput(x, hs[i], θ[6], λ, R, u[i], v[i]) # update state and compute output
            y[j] = yi
            j += 1
        else
            x = updatestate(x, hs[i], λ, u[i], v[i]) # update state
        end
    end
    return y
end