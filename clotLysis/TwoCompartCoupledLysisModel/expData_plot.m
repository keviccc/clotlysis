%% Experimental data plot for PKPD
% units: time in [min], tPA conc. in [micro mol/L] and others in [%]
% variables: N_patient, dose_0X (X: number), duration_0X, conc <--str of proteins names
% optional var.: tPA_patient_0X, FBG_patient_0X, AP_patient_0X and PLG_patient_0X
clc; clear all; format long; close all;
% --- Verstraete et al. Journal of Pharmacology and Experimental Therapeutics 235 (1985)
% tPA, FBG and AP, 3 patients and 1 healthy (only tPA)
load('Tanswell1989.mat');
for ii=1:N_patient
    figureName_temp = ['Group ', num2str(ii)];
    for jj=1:length(conc)
        figureName = [figureName_temp,'-',conc{jj}];
        nameOfdata2plot = [conc{jj},'_patient_0',num2str(ii)];
        data2plot = eval(nameOfdata2plot);
        figure('Name',figureName)
        plot(data2plot(:,1),data2plot(:,2),'--o');
        xlabel('Time [min]');
        if strcmp(conc{jj},'tPA')
            ylabel('Conc. [\muM]');
        else; ylabel('Conc. [%]');
        end
    end
end