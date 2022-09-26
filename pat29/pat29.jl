# generate data for patient 29 to paper on pksim
# date 220906


using Pkg
Pkg.activate(".")

using CSV, DataFrames, StaticArrays
using BenchmarkTools

include("../PKsimulation.jl") # functions to simulate state and compute output
include("../PKmodels.jl") # calculations of λ, R from parameter vector θ

# Load datafiles
modeldf = CSV.read("inputs/eleveld_modelparams.csv", DataFrame) # Model parameters V1,V2,... for PK model
inputdf = CSV.read("inputs/eleveld_infusiondata.csv", DataFrame) # Input vector for simulation, with time and rates/ doses for infusions and boluses
idxdf = CSV.read("inputs/eleveld_youts.csv", DataFrame) # Indices of time points for final concentration vector with time points (relates to indices in eleveld_infusiondata.csv)

include("../getdata.jl") # Functions to get input data from files (model parameters and input data)

# Simulate one patient
studynbr = 29
id = 842 # long time series
study_df, first, _ = getstudydata(studynbr)
θ, infusionrate, bolusdose, time, h, youts = getpatientdata(id, study_df)
y = zeros(Float32, length(youts)) # create empty output vector
PKsim!(y, θ, infusionrate, bolusdose, h, youts)
yobs = y
y = zeros(Float32, length(infusionrate)) # create empty output vector
PKsim!(y, θ, infusionrate, bolusdose, h, 1:length(infusionrate)) # full simulation, all steps
yfull = y

observationdata = [time[youts] ./ 60 yold]
fulldata = [time ./ 60 yfull]
inputdata = [time ./ 60 infusionrate]

using DelimitedFiles

writedlm("csv/observationdata_pat29.csv", observationdata, ',')
writedlm("csv/fulldata_pat29.csv", fulldata, ',')
writedlm("csv/inputdata_pat29.csv", inputdata, ',')