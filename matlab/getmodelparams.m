function [k10,k12,k13,k21,k31,V1] = getmodelparams(k,modelparams_study)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here



% Create PK model
% drug transfer rate constants
CL = CL/60; % [l/min]
Q2 = Q2/60;
Q3 = Q3/60;
k10 = CL / V1;                     % [1/s]
k12 = Q2 / V1;                     % [1/s]
k13 = Q3 / V1;                     % [1/s]
k21 = Q2 / V2;                     % [1/s]
k31 = Q3 / V3;                     % [1/s]

end

