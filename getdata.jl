
# Functions to get data from data files put in global variables modeldf, inputdf, idxdf
# modeldf: Model parameters V1,V2,... for PK model
# inputdf: Input vector for simulation, with time and rates/ doses for infusions and boluses
# idxdf: Indices of time points for final concentration vector with time points (relates to indices in inputdf)

# Functions to extract input data and model parameters from files
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

    θ = @SVector [k10, k12, k13, k21, k31, V1] # PK model parameter vector

    # Input data
    inputs = inputdf[in([id]).(inputdf.ID), :]
    youtstime = idxdf[in([id]).(idxdf.ID), :]

    time = Float32.(inputs[!, :Time]) #[s]
    infusionrate = Float32.(inputs[!, :InfusionRate]) # unit?
    bolusdose = Float32.(inputs[!, :Bolus]) # unit?
    h = diff(time) # durations between inputs
    h = [h; h[end]] # assume last h is same as second last.
    # Indices of measured concentrations
    youts = Int32.(youtstime[!, :Idx])
    time = time[youts]
    return Float32.(θ), infusionrate, bolusdose, time, h, youts
end