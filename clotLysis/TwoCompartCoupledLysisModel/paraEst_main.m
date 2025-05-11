%% ---------------------------------------------------------------------- %
% Parameter Estimation Main Script.                                       %
% ----------------------------------------------------------------------- %

clc; clear all; format long;

global initValueSave fileID

fileID = fopen('paraEstLog.txt','w');

[k1f,k1r,k2,k3f,k3r,k4,k5f,k5r,k6,kMG,kPAI] = kinetic_para();
initValue = [k1f,k1r,k2,k3f,k3r,k4,k5f,k5r,k6,kMG,kPAI];
initValueSave = initValue;
options = [];
% [x,fval,exitflag,output] = fminsearch('paraEst_expData_func',initValue,options);
lb = initValue*1e-3;
ub = initValue*1e3;
[x,fval,exitflag,output] = fmincon('paraEst_expData_func',initValue,[],[],[],[],lb,ub,[],options);

fclose(fileID);