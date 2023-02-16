using Pkg
Pkg.activate(".")

using CSV, DataFrames, StaticArrays
using BenchmarkTools

include("pksimulation.jl") # functions to simulate state and compute output
include("pkmodels.jl") # calculations of λ, R from parameter vector θ

# Load datafiles
modeldf = CSV.read("../inputs/eleveld_modelparams.csv", DataFrame) # Model parameters V1,V2,... for PK model
inputdata = CSV.read("pat29.csv", DataFrame, header=0) # Model parameters V1,V2,... for PK model

# inputdf = CSV.read("../inputs/eleveld_infusiondata.csv", DataFrame) # Input vector for simulation, with time and rates/ doses for infusions and boluses
# idxdf = CSV.read("inputs/eleveld_youts.csv", DataFrame) # Indices of time points for final concentration vector with time points (relates to indices in eleveld_infusiondata.csv)

# include("../getdata.jl") # Functions to get input data from files (model parameters and input data)

# Simulate one patient
studynbr = 29
id = 842 # long time series
# study_df, first, _ = getstudydata(studynbr)
# θ, infusionrate, bolusdose, time, h, youts = getpatientdata(id, study_df)
study_df = modeldf[in([studynbr]).(modeldf.StudyNbr), :] # get all rows with nbr = studynbr
subject = study_df[in([id]).(study_df.ID), :] # get all rows with nbr = id

t = inputdata[:,2]
infusionrate = inputdata[:, 3]
bolusdose = inputdata[:, 4]
ismeasurements = inputdata[:,5]

V1 = subject[!, :V1][1]
V2 = subject[!, :V2][1]
V3 = subject[!, :V3][1]
CL = subject[!, :CL][1]
Q2 = subject[!, :Q2][1]
Q3 = subject[!, :Q3][1]

# drug transfer rate constants
CL /= 60 # [l/min]
Q2 /= 60
Q3 /= 60
k10 = CL / V1 # [1/s]
k12 = Q2 / V1
k13 = Q3 / V1
k21 = Q2 / V2
k31 = Q3 / V3

θ = @SVector [k10, k12, k13, k21, k31, V1] # PK model parameter vector
h = diff(t) # durations between inputs
h = [h; h[end]] # assume last h is same as second last.
# Indices of measured concentrations
youts = findall(x->x!=0.0,ismeasurements)
# t = t[youts]


y = zeros(Float32, length(youts)) # create empty output vector
pk3sim!(y, θ, infusionrate, bolusdose, h, youts)
yobs = y
y = zeros(Float32, length(infusionrate)) # create empty output vector
pk3sim!(y, θ, infusionrate, bolusdose, h, 1:length(infusionrate)) # full simulation, all steps
yfull = y

observationdata = [t[youts] ./ 60 yobs]
fulldata = [t ./ 60 yfull]
infdata = [t ./ 60 infusionrate]

using Plots
plot(t./60, infusionrate)
plot(t ./ 60, yfull)
plot!(t[youts] ./ 60, yobs)
