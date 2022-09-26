% Eleveld simulation
% Date: 220912
% first test

clear all
close all
warning('off','all')

modelparams = csvread("../inputs/eleveld_modelparams.csv",1); % order: patient ID, studynbr, V1, V2, V3, CL, Q2, Q3
inputs = csvread("../inputs/eleveld_infusiondata.csv",1); % order: ID,Time,InfusionRate,Bolus
youts_idxs = csvread("../inputs/eleveld_youts.csv",1); % order: ID,Idx,Time
%meas = csvread("../inputs/eleveld_measurements.csv",1); % order: ID,Time,Cp,PRED

% struct with patient data
S = struct();


% create struct for data
for studynbr = 1:30
study_ind = find(modelparams(:,2)==studynbr);
modelparams_study = modelparams(study_ind,:);

v1 = modelparams_study(:,3);
v2 = modelparams_study(:,4);
v3 = modelparams_study(:,5);
cl = modelparams_study(:,6);
q2 = modelparams_study(:,7);
q3 = modelparams_study(:,8);
cl = cl/60; % [l/min]
q2 = q2/60;
q3 = q3/60;
k10 = cl ./ v1;                     % [1/s]
k12 = q2 ./ v1;                     % [1/s]
k13 = q3 ./ v1;                     % [1/s]
k21 = q2 ./ v2;                     % [1/s]
k31 = q3 ./ v3;                     % [1/s]

k = 1; % local id
for ID = modelparams_study(:,1)' % global id
    % model parameters
    S(studynbr).Patient(k).V1 = v1(k);
    S(studynbr).Patient(k).V2 = v2(k);
    S(studynbr).Patient(k).V3 = v3(k);
    S(studynbr).Patient(k).CL = cl(k);
    S(studynbr).Patient(k).Q2 = q2(k);
    S(studynbr).Patient(k).Q3 = q3(k);
    S(studynbr).Patient(k).k10 = k10(k);
    S(studynbr).Patient(k).k12 = k12(k);
    S(studynbr).Patient(k).k13 = k13(k);
    S(studynbr).Patient(k).k21 = k21(k);
    S(studynbr).Patient(k).k31 = k31(k);
    
    % inputs
    ind_k = find(inputs(:,1)==ID);
    input_k = inputs(ind_k,:);
    S(studynbr).Patient(k).InputTime = input_k(:,2);
    S(studynbr).Patient(k).Infusion = input_k(:,3);
    S(studynbr).Patient(k).Bolus = input_k(:,4);
    
    % observations
    ind_meas_k = find(youts_idxs(:,1)==ID);
    youts_k = youts_idxs(ind_meas_k,:);
    S(studynbr).Patient(k).youtsTime = youts_k(:,3);
    S(studynbr).Patient(k).youtsIdxs = youts_k(:,2);
    

    % precompute PK model
    A = [-(k10(k) + k12(k) + k13(k)) k12(k) k13(k)
    k21(k) -k21(k) 0.0
    k31(k) 0.0 -k31(k)];

    B = [1 / v1(k); 0; 0];
    C = [1 0 0];
    D = 0;

    PK = ss(A, B, C, D);
    S(studynbr).Patient(k).PK = PK;
    
    k = k+1;
end

end

%% simulate using lsim
k = 1;
studynbr = 29;

PK = S(studynbr).Patient(k).PK;
uinf = S(studynbr).Patient(k).Infusion;
ubol = S(studynbr).Patient(k).Bolus;
tu = S(studynbr).Patient(k).InputTime;
youts = S(studynbr).Patient(k).youtsIdxs;

f = @() simulate_lsim(PK,uinf,ubol,tu,youts);
t1 = timeit(f); % 0.2069

%y = simulate_lsim(PK,uinf,ubol,tu,youts);

%% simulate all patients in the eleveld data set

totaltime = 0.0;
for studynbr = 1:30
    lastid = length(S(studynbr));
    for k = 1:lastid
        PK = S(studynbr).Patient(k).PK;
        uinf = S(studynbr).Patient(k).Infusion;
        ubol = S(studynbr).Patient(k).Bolus;
        tu = S(studynbr).Patient(k).InputTime;
        youts = S(studynbr).Patient(k).youtsIdxs;

        f = @() simulate_lsim(PK,uinf,ubol,tu,youts);
        t1 = timeit(f); 
        totaltime = totaltime + t1;
    end
end

totaltime
