# Simulation of population Eleveld PK model from paper "Pharmacokinetic-pharmacodynamic model for propofol for broad application in anaesthesia and sedation"
# Simulates using computation of the exponential matrix in discrete time at every time instance. Three-compartment mammillary model.
# x(k+1) = A_d*x(k) + B_d u(k), y(k) = C*x(k).
# Date: 220926

# TODO: run full benchmark with julia 1.8

using Pkg
Pkg.activate("..")

using ControlSystems, CSV, DataFrames, StaticArrays
using BenchmarkTools

include("../expm.jl") # Exponential matrix computation 

# Load datafiles
modeldf = CSV.read("../inputs/eleveld_modelparams.csv", DataFrame) # Model parameters V1,V2,... for PK model
inputdf = CSV.read("../inputs/eleveld_infusiondata.csv", DataFrame) # Input vector for simulation, with time and rates/ doses for infusions and boluses
idxdf = CSV.read("../inputs/eleveld_youts.csv", DataFrame) # Indices of time points for final concentration vector with time points (relates to indices in eleveld_infusiondata.csv)


include("../getdata.jl") # Functions to get input data from files (model parameters and input data)

# Simulation using c2d in ControlSystems.jl and updating state using x(k+1) = A_d*x(k) + B_d u(k), where A_d requires computation of the exponential matrix. Simulates blood plasma concentration of the central compartment.
function simulate(θ, u, v, h, youts)
    k10, k12, k13, k21, k31, V1 = θ
    A = [-(k10 + k12 + k13) k12 k13
        k21 -k21 0.0
        k31 0.0 -k31]

    B = [1 / V1; 0; 0]
    C = [1 0 0]
    D = 0

    PK = ss(A, B, C, D)

    nt = length(u)
    y = zeros(nt,)
    x = [0, 0, 0] # Initial state vector

    for i = 1:nt-1
        if !iszero(v[i])
            x = x + PK.B * v[i] # bolus
            y[i] = x[1]
        end
        if h[i] == 0.0
            y[i+1] = y[i]
            continue
        end
        PK_d = c2d(PK, h[i])
        x = PK_d.A * x + PK_d.B * u[i] # infusion
        y[i+1] = x[1] # infusion affects next sample
    end
    y = y[youts] # observed outputs
    return y
end


# Simulating one patient
studynbr = 29
id = 842 # long time series
study_df, first, _ = getstudydata(studynbr)
θ, infusionrate, bolusdose, time, h, youts = getpatientdata(id, study_df)
y = simulate(θ, infusionrate, bolusdose, h, youts)


# Benchmark one patient
@btime simulate($θ, $infusionrate, $bolusdose, $h, $youts) # 618 us, 3713 allocations


## Benchmarking simulation
# Simulation times of each patient are added together
nstudy = 30 # nbr of studies
benchtime = 0.0 # ns
nallocs = 0
for studynbr = 1:nstudy
    studynbr = 1
    @show studynbr
    study_df, firstid, lastid = getstudydata(studynbr) # get dataframe for this study
    for id = firstid:lastid
        @show id
        if id in [893, 897] # no measurements exists for these patients
            continue
        end
        θ, infusionrate, bolusdose, time, h, youts = getpatientdata(id, study_df) # get patient data
        bm = @benchmark simulate($θ, $infusionrate, $bolusdose, $h, $youts) # add samples = 100 (or similar?)
        benchtime += median(bm.times) # median simulation time
        nallocs += bm.allocs
    end
end

benchtime # 5.481491855e8 ns = 548 ns (julia 1.7)
nallocs # 2.027384e6
