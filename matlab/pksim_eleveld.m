% Eleveld simulation using lsim
% Date: 221128

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

    v1 = modelparams_study(:,3); % [l]
    v2 = modelparams_study(:,4);
    v3 = modelparams_study(:,5);
    cl = modelparams_study(:,6); % [l/s]
    q2 = modelparams_study(:,7);
    q3 = modelparams_study(:,8);
    cl = cl/60; % [l/min]
    q2 = q2/60;
    q3 = q3/60;
    k10 = cl ./ v1;                     % [1/min]
    k12 = q2 ./ v1;                     % [1/min]
    k13 = q3 ./ v1;                     % [1/min]
    k21 = q2 ./ v2;                     % [1/min]
    k31 = q3 ./ v3;                     % [1/min]

    k = 1; % local id
    for ID = modelparams_study(:,1)' % global id
        % model parameters
        S(studynbr).Patient(k).V1 = v1(k); % [l]
        S(studynbr).Patient(k).V2 = v2(k);
        S(studynbr).Patient(k).V3 = v3(k);
        S(studynbr).Patient(k).CL = cl(k); % [l/min]
        S(studynbr).Patient(k).Q2 = q2(k);
        S(studynbr).Patient(k).Q3 = q3(k);
        S(studynbr).Patient(k).k10 = k10(k); % [1/min]
        S(studynbr).Patient(k).k12 = k12(k);
        S(studynbr).Patient(k).k13 = k13(k);
        S(studynbr).Patient(k).k21 = k21(k);
        S(studynbr).Patient(k).k31 = k31(k);


        % inputs
        ind_k = find(inputs(:,1)==ID);
        input_k = inputs(ind_k,:);
        S(studynbr).Patient(k).InputTime = input_k(:,2); % [s]
        S(studynbr).Patient(k).Infusion = input_k(:,3); % [ml/min]
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

        B = [1 / (v1(k)); 0; 0]; % use [ml] as unit
        C = [1 0 0];
        D = 0;

        PK = ss(A, B, C, D);
        S(studynbr).Patient(k).PK = PK;

        k = k+1;
    end

end

%% simulate one patient using lsim
k = 1;
studynbr = 29;

PK = S(studynbr).Patient(k).PK;
uinf = S(studynbr).Patient(k).Infusion;
ubol = S(studynbr).Patient(k).Bolus;
tu = S(studynbr).Patient(k).InputTime;
youts = S(studynbr).Patient(k).youtsIdxs;

f = @() simulate_lsim(PK,uinf,ubol,tu,youts);
t1 = timeit(f); % 0.2137 s

y = simulate_lsim(PK,uinf,ubol,tu,youts);
plot(S(studynbr).Patient(k).youtsTime./60,y,'*')
title('Smiulated output of patient 842 in study 29')
xlabel('Time (min)')
ylabel('Blood plasma concentration (mg/ml)')

%% Simulate all patients in the Eleveld data set

neval = 10; % Number of evaluations per patient
totaltime = 0.0;
telapsed = zeros(neval,1);
for studynbr = 1:30
    studynbr
    lastid = length(S(studynbr).Patient);
    for k = 1:lastid
        k
        PK = S(studynbr).Patient(k).PK;
        uinf = S(studynbr).Patient(k).Infusion;
        ubol = S(studynbr).Patient(k).Bolus;
        tu = S(studynbr).Patient(k).InputTime;
        youts = S(studynbr).Patient(k).youtsIdxs;
        
        for i = 1:neval
            tstart = tic;
            y = simulate_lsim(PK,uinf,ubol,tu,youts);
            telapsed(i) = toc(tstart);
        end
        t1 = median(telapsed); % Compute median computation time
        totaltime = totaltime + t1;
    end
end

totaltime % 64.5848 seconds
