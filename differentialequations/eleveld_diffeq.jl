# First attempt of simulating the Eleveld model using DifferentialEquations.jl
# Use callbacks for change of infusion rates and bolus doses

# FIXME: Does not generate correct output! Compare with output from eleveld_pksim.jl! For the example patient, only the first samples are correct.

using Pkg
Pkg.activate("..")

using ControlSystemsBase, CSV, DataFrames, LinearAlgebra
using BenchmarkTools
using DifferentialEquations
using StaticArrays

# Load datafiles
modeldf = CSV.read("../inputs/eleveld_modelparams.csv", DataFrame) # Model parameters V1,V2,... for PK model
inputdf = CSV.read("../inputs/eleveld_infusiondata.csv", DataFrame) # Input vector for simulation, with time and amounts for infusions and boluses
idxdf = CSV.read("../inputs/eleveld_youts.csv", DataFrame) # Indices of time points for final concentration vector with time points (relates to inputdf)

include("../getdata.jl") # Functions to get input data from files (model parameters and input data)

# In this case, we are not interested in the indicies of the measured outputs, but rather the specific time points. Therefore, this function is modified
function getpatientdata(id, study_df)
    subject = study_df[in([id]).(study_df.ID), :] # get all rows with nbr = id

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

    θ = [k10, k12, k13, k21, k31, V1] # PK model parameter vector

    # Input data
    inputs = inputdf[in([id]).(inputdf.ID), :]
    youtstime = idxdf[in([id]).(idxdf.ID), :]

    time = Float32.(inputs[!, :Time]) #[s]
    infusionrate = Float32.(inputs[!, :InfusionRate]) # unit?
    bolusdose = Float32.(inputs[!, :Bolus]) # unit?
    h = diff(time)
    h = [h; h[end]]
    # Indices of measured concentrations
    youts = Int32.(youtstime[!, :Idx])
    #time = time[youts]

    θdiffeq = θ[1:6]
    tol = 1e-6
    tspan = [time[1], time[end] + tol]
    youtstime = time[youts] .+ tol

    return Float32.(θ), infusionrate, bolusdose, time, h, youts, youtstime, tspan
end


# PK model
function pk3!(dx, x, p, t)
    k10, k12, k13, k21, k31, V1, w = p
    dx[1] = -(k10 + k12 + k13) * x[1] + k12 * x[2] + k13 * x[3] + (1 / V1) * w
    dx[2] = k21 * x[1] - k21 * x[2]
    dx[3] = k31 * x[1] - k31 * x[3]
end

function sim_diffeq(θ, uinf, ubol, utime, youtstime, tspan)
    ibolstart = 1 # if no bolus at start
    iinfstart = 1
    if ubol[1] != 0.0 # bolus in first time point t=0
        x0 = [popfirst!(ubol) / θ[6], 0.0, 0.0]
        ibolstart = 2
    else
        x0 = [0.0, 0.0, 0.0]
    end
    if uinf[1] != 0.0
        w = popfirst!(uinf)
        iinfstart = 2
    else
        w = 0.0
    end
    # end
    p = [θ; w]
    prob = ODEProblem(pk3!, x0, tspan, p)

    # do these need to be inside function?
    affectbolus!(integrator) = integrator.u[1] += popfirst!(ubol) / θ[6] # can probably do this in a better way
    cbbolus = PresetTimeCallback(utime[ibolstart:end-1], affectbolus!, save_positions=(false, false))

    affectinfusion!(integrator) = integrator.p[7] = popfirst!(uinf)
    cbinfusion = PresetTimeCallback(utime[iinfstart:end], affectinfusion!, save_positions=(false, false))

    cbs = CallbackSet(cbbolus, cbinfusion)

    # Solve
    sol = solve(prob, callback=cbs, saveat=youtstime)
    # sol = solve(prob, callback=cbbolus, saveat=0.01)
    # return sol
    y = reduce(hcat, sol.u)[1, :]
    return y
end


# Simulate one patient
# studynbr = 29
# id = 842 # long time series
# study_df, first, _ = getstudydata(studynbr)
# θ, infusionrate, bolusdose, time, h, youts, youtstime, tspan = getpatientdata(id, study_df)
# ydiffeq = sim_diffeq(θ, copy(infusionrate), copy(bolusdose), time, youtstime, tspan)

# Output does not match result from eleveld_pksim.jl at all time instances! why?

# @btime sim_diffeq($θ, $copy(infusionrate), $copy(bolusdose), $time, $youtstime, $tspan)


## Simulate all patients
nstudy = 30 # nbr of studies
for studynbr = 1:nstudy
    study_df, firstid, lastid = getstudydata(studynbr) # get dataframe for this study
    for id = firstid:lastid
        if id in [893, 897] # no measurements exists for these patients
            continue
        end
        θ, infusionrate, bolusdose, time, h, youts, youtstime, tspan = getpatientdata(id, study_df)
        sim_diffeq(θ, infusionrate, bolusdose, time, youtstime, tspan)
    end
end

## Benchmarking simulation
# Simulation times of each patient are added together
function runsim()
    nstudy = 30 # nbr of studies
    benchtime = 0.0 # ns
    nallocs = 0
    for studynbr = 1:nstudy
        # @show studynbr
        study_df, firstid, lastid = getstudydata(studynbr) # get dataframe for this study
        for id = firstid:lastid
            @show id
            if id in [893, 897] # no measurements exists for these patients
                continue
            end
            θ, infusionrate, bolusdose, time, h, youts, youtstime, tspan = getpatientdata(id, study_df)
            # bm = @benchmark sim_diffeq($θ, $deepcopy(infusionrate), $deepcopy(bolusdose), $time, $youtstime, $tspan) samples = 10 evals = 10 gctrial = false # does not work, why?
            nit = 10 # nbr of iterations per patient
            bit = zeros(nit)
            for it = 1:nit # nbr of trials
                bit[it] = @elapsed sim_diffeq(θ, copy(infusionrate), copy(bolusdose), time, youtstime, tspan)
            end
            benchtime += median(bit)
            nallocs += @allocated sim_diffeq(θ, infusionrate, bolusdose, time, youtstime, tspan)
            # θ, infusionrate, bolusdose, time, h, youts = getpatientdata(id, study_df) # get patient data
            # bm = @benchmark simulate($θ, $infusionrate, $bolusdose, $h, $youts) 
            # benchtime += median(bm.times) # median simulation time
            # nallocs += bm.allocs
        end
    end
    benchtime, nallocs
end

benchtime, nallocs = runsim()
# 0.23 s = 230.3 ms, 27869280 allocations
