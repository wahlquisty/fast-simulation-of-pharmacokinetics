# Simulating the Eleveld model using DifferentialEquations.jl
# Use callbacks for change of infusion rates and bolus doses

using Pkg
cd(@__DIR__)
Pkg.activate("..")

using CSV, DataFrames, DifferentialEquations, LinearAlgebra, StaticArrays
using BenchmarkTools

# Load datafiles
modeldf = CSV.read("../inputs/eleveld_modelparams.csv", DataFrame) # Model parameters V1,V2,... for PK model
inputdf = CSV.read("../inputs/eleveld_infusiondata.csv", DataFrame) # Input vector for simulation, with time and amounts for infusions and boluses
idxdf = CSV.read("../inputs/eleveld_youts.csv", DataFrame) # Indices of time points for final concentration vector with time points (relates to inputdf)

# Modified functions from eleveld_pksim.jl
function getstudydata(studynbr)
    study_df = modeldf[in([studynbr]).(modeldf.StudyNbr), :] # get all rows with nbr = studynbr
    idnbrs = study_df[!, :ID] # get array with all id-numbers in current study
    last = idnbrs[length(idnbrs)] # last id-number of current study
    first = idnbrs[1] # first id-number of current study
    return study_df, first, last
end

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

    affectbolus!(integrator) = integrator.u[1] += popfirst!(ubol) / θ[6] # can probably do this in a better way
    cbbolus = PresetTimeCallback(utime[ibolstart:end-1], affectbolus!, save_positions=(false, false))

    affectinfusion!(integrator) = integrator.p[7] = popfirst!(uinf)
    cbinfusion = PresetTimeCallback(utime[iinfstart:end], affectinfusion!, save_positions=(false, false))

    cbs = CallbackSet(cbbolus, cbinfusion)

    # Solve
    sol = solve(prob, callback=cbs, saveat=youtstime)
    y = reduce(hcat, sol.u)[1, :]
    return y
end


## Simulate all patients
nstudy = 30 # nbr of studies
@time for studynbr = 1:nstudy
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
            θ, infusionrate, bolusdose, time, _, _, youtstime, tspan = getpatientdata(id, study_df)
            ni = 10 # nbr of iterations per patient
            ti = zeros(ni)
            for it = 1:ni # nbr of trials
                ti[it] = @elapsed sim_diffeq(θ, copy(infusionrate), copy(bolusdose), time, youtstime, tspan)
            end
            benchtime += median(ti)
            nallocs += @allocated sim_diffeq(θ, infusionrate, bolusdose, time, youtstime, tspan)
        end
    end
    benchtime, nallocs
end

benchtime, nallocs = runsim()
# 0.2304 s = 230.4 ms, 27869280 allocations
