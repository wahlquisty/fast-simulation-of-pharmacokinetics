# Date: 220926
# Full state vector simulation (returns all three states of x)

using StaticArrays, LinearAlgebra, SLEEFPirates

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
    Φdiag = @. SLEEFPirates.exp(λ * h) # cannot differentiate through VectorizationBase.vexp. Remove 1 us in compilation time
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

# Initiate/update state
@inline @fastmath function updateandgetstate(x, h, V1, λ, R, u=0.0, v=0.0)
    Φdiag = getΦdiag(λ, h) # compute Φ
    x = bolus(x, v) # update state for bolus
    xtrue = getstates(x, V1, R)
    x = step(x, λ, Φdiag, u) # infusion affect next sample
    return x, xtrue
end

# Update state and compute output
@inline @fastmath function updatestateoutput(x, h, V1, λ, R, u=0.0, v=0.0; xout=false)
    Φdiag = getΦdiag(λ, h) # compute Φ
    x = bolus(x, v) # update state for bolus
    y = gety(x, V1, R) # compute output
    x = step(x, λ, Φdiag, u) # infusion affect next sample
    return x, y
end

# Update state and compute output
@inline @fastmath function updateandgetstateoutput(x, h, V1, λ, R, u=0.0, v=0.0; xout=false)
    Φdiag = getΦdiag(λ, h) # compute Φ
    x = bolus(x, v) # update state for bolus
    xtrue = getstates(x, V1, R)
    y = gety(x, V1, R) # compute output
    x = step(x, λ, Φdiag, u) # infusion affect next sample
    return x, xtrue, y
end

function getstates(x, V1, R)
    return 1 / V1 * R' * x
end


function PKsim!(Cp, θ, u, v, h, youts; xout=false)
    λ, R = update(θ) # Setting up simulator
    nyouts = length(youts)
    nt = length(u)

    x = @SVector [0.0f0, 0.0f0, 0.0f0] # initial state

    if xout
        xouts = zeros(Float32, nyouts, 3)
        j = 1 # next empty index in Cp
        for i = 1:nt
            if i in youts # if we want to compute output
                x, xtrue, yi = updateandgetstateoutput(x, h[i], θ[6], λ, R, u[i], v[i]) # update state and compute output
                Cp[j] = yi
                xouts[j, :] = xtrue
                j += 1 # update next empty index in Cp
            else
                x = updatestate(x, h[i], λ, u[i], v[i]) # update state
            end
        end
        return Cp, xouts
    else
        j = 1
        for i = 1:nt
            if i in youts # if we want to compute output
                x, yi = updatestateoutput(x, h[i], θ[6], λ, R, u[i], v[i]) # update state and compute output
                Cp[j] = yi
                j += 1
            else
                x = updatestate(x, h[i], λ, u[i], v[i]) # update state
            end
        end
        return Cp
    end
end
