# Fast simulation of population Eleveld PK model from paper "Pharmacokinetic-pharmacodynamic model for propofol for broad application in anaesthesia and sedation"
# Date: 220926

# Simulates the three-compartment mammillary model. Does not compute all states, only first state as output (simpler R). Simulates blood plasma concentration of the central compartment.

# TODO: Compute and compare median or minimum time of simulation? @btime returns minimum time, with @benchmark we can get any measure.
# TODO: How many iterations for @benchmark?
# TODO: run full benchmark with julia 1.8

using Pkg
Pkg.activate(".")

using CSV, DataFrames, StaticArrays
using BenchmarkTools

include("PKsimulation.jl") # functions to simulate state and compute output
include("PKmodels.jl") # calculations of λ, R from parameter vector θ

# Load datafiles
modeldf = CSV.read("inputs/eleveld_modelparams.csv", DataFrame) # Model parameters V1,V2,... for PK model
inputdf = CSV.read("inputs/eleveld_infusiondata.csv", DataFrame) # Input vector for simulation, with time and rates/ doses for infusions and boluses
idxdf = CSV.read("inputs/eleveld_youts.csv", DataFrame) # Indices of time points for final concentration vector with time points (relates to indices in eleveld_infusiondata.csv)

include("getdata.jl") # Functions to get input data from files (model parameters and input data)

# Simulate one patient
studynbr = 29
id = 842 # long time series
study_df, first, _ = getstudydata(studynbr)
θ, infusionrate, bolusdose, time, h, youts = getpatientdata(id, study_df)
y = zeros(Float32, length(youts)) # create empty output vector
PKsim!(y, θ, infusionrate, bolusdose, h, youts)

# Benchmark one patient
@btime PKsim!($y, $θ, $infusionrate, $bolusdose, $h, $youts) # 2.6 us, 0 allocations


## Simulate all patients
nstudy = 30 # nbr of studies
for studynbr = 1:nstudy
    study_df, firstid, lastid = getstudydata(studynbr) # get dataframe for this study
    for id = firstid:lastid
        if id in (893, 897) # no measurements exists for these patients
            continue
        end
        θ, infusionrate, bolusdose, time, h, youts = getpatientdata(id, study_df) # get patient data
        y = zeros(Float32, length(youts)) # create empty output vector
        PKsim!(y, θ, infusionrate, bolusdose, h, youts)
    end
end


## Benchmarking for fast simulation
# Simulation times of each patient are added together
# Total number of allocations: 0.
function runsim()
    nstudy = 30 # nbr of studies
    benchtime = 0.0 # ns
    for studynbr = 1:nstudy
        # @show studynbr
        study_df, firstid, lastid = getstudydata(studynbr) # get dataframe for this study
        for id = firstid:lastid
            @show id
            if id in (893, 897) # no measurements exists for these patients
                continue
            end
            θ, infusionrate, bolusdose, time, h, youts = getpatientdata(id, study_df) # get patient data
            yset = Set(youts) 
            y = zeros(Float32, length(youts)) # create empty output vector
            bm = @benchmark PKsim!($y, $θ, $infusionrate, $bolusdose, $h, $yset) samples=5 evals=5 gctrial=false
            benchtime += median(bm.times) # median simulation time
            # bm = @elapsed PKsim!(y, θ, infusionrate, bolusdose, h, yset)
            # benchtime += bm# median simulation time
        end
    end
    benchtime
end

benchtime = runsim() # 2.14804786e6 nanoseconds = 2.148 ms (for julia 1.7)


