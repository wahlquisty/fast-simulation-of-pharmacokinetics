# Fast simulation of population Eleveld PK model from paper "Pharmacokinetic-pharmacodynamic model for propofol for broad application in anaesthesia and sedation"

# Simulates the three-compartment mammillary model. Does not compute all states, only first state as output (simpler R). Simulates blood plasma concentration of the central compartment.

using Pkg
cd(@__DIR__)
Pkg.activate(".")
Pkg.add(url="https://github.com/wahlquisty/FastPKSim.jl") # FastPKSim.jl from github

using FastPKSim
using CSV, DataFrames, StaticArrays
using BenchmarkTools

# # Load datafiles
modeldf = CSV.read("inputs/eleveld_modelparams.csv", DataFrame) # Model parameters V1,V2,... for PK model
inputdf = CSV.read("inputs/eleveld_infusiondata.csv", DataFrame) # Input vector for simulation, with time and rates/ doses for infusions and boluses
idxdf = CSV.read("inputs/eleveld_youts.csv", DataFrame) # Indices of time points for final concentration vector with time points (relates to indices in eleveld_infusiondata.csv)

include("getdata.jl") # Functions to get input data from files (model parameters and input data)

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
        pk3sim!(y, θ, infusionrate, bolusdose, h, youts)
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
            # @show id
            if id in (893, 897) # no measurements exists for these patients
                continue
            end
            θ, infusionrate, bolusdose, time, h, youts = getpatientdata(id, study_df) # get patient data
            yset = Set(youts)
            y = zeros(Float32, length(youts)) # create empty output vector
            bm = @benchmark pk3sim!($y, $θ, $infusionrate, $bolusdose, $h, $yset) samples = 10 evals = 10 gctrial = false
            benchtime += median(bm.times) # median simulation time
        end
    end
    benchtime
end

benchtime = runsim() # 1.5833 e6 nanoseconds = 1.5833 ms (for julia 1.8.2)


