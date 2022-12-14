# Mamillary three-compartment model functions
# Date: 220926

# Computes R for all states (note, specific for 3 compartment mammillary model)

using StaticArrays, LinearAlgebra

# Initiate/update parameters
function update(θ)
    λ = getλ(θ)
    R = getR(θ, λ)
    return λ, R
end

# Compute λ for 3 compartment mammillary model
@inline @fastmath function getλ(θ)
    k10, k12, k13, k21, k31, _ = θ
    b1 = k10 + k12 + k13 + k21 + k31
    b2 = k21 * (k10 + k13 + k31) + k31 * (k10 + k12)
    b3 = k10 * k21 * k31

    # Wengert list used to compute λ.
    a1 = b1 / 3
    a2 = a1^2
    a3 = a1 * a2
    a4 = b2 / 3
    a5 = a4 - a2
    a6 = (b1 * a4 - b3) / 2
    a7 = 2(a6 + sqrt(complex(a5^3 + (a3 - a6)^2)) - a3)^(1 / 3.0f0)
    a8 = -real(a7)
    a9 = imag(a7)
    a10 = a9 * sqrt(3) / 2
    a11 = a1 - a8 / 2

    return @SVector [-a1 - a8, -a10 - a11, a10 - a11] # The eigenvalues of the continuous-time system matrix
end

# Compute R for all states (note, specific for 3 compartment mammillary model)
@inline function getR(θ, λ)
    _, _, _, k21, k31, _ = θ
    l1, l2, l3 = λ
    # b = @SVector [1, k21 + k31, k21 * k31]     # Quite often we would only be interested in the first column
    a1 = l2 - l1
    a2 = l3 - l1
    a3 = l3 - l2
    d1 = a2 * a1
    d2 = -a3 * a1
    d3 = a2 * a3

    # Qinv = @SMatrix [[l1^2 l1 1] / d1; [l2^2 l2 1] / d2; [l3^2 l3 1] / d3]
    Qinv = @SMatrix [l1^2/d1 l1/d1 1/d1; l2^2/d2 l2/d2 1/d2; l3^2/d3 l3/d3 1/d3]
    b = @SMatrix [1 0 0; k21+k31 k21 k31; k21*k31 k21*k31 k21*k31] # Each column in b corresponds to the b of one output.
    return Qinv * b
end