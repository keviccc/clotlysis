%% ---------------------------------------------------------------------- %
% Comparison of compartmental model and experimental data                 %
%  - Plot simulation results and experimental data from the literature    %
%  - Estimation of kinetics parameters                                    %
% programmed by BG, 20/12/2018                                            %
% ----------------------------------------------------------------------- %

clc; clear all; format long; close all;
global kel_PLG kel_PLS kel_FBG kel_AP kel_MG kel_PAI k12 k21 k10
global KM_PLG k_PLG_cat k_AP_f k_AP_r k_AP_cat KM_FBG k_FBG_cat kMG kPAI
global C_tPA0 C_PLG0 C_PLS0 C_FBG0 C_AP0 C_MG0 C_PAI0 

% fileNames = {'Verstraete1985','Tienfenbrunn1986','Collen1986',...
%     'Noe1987','Tebbe1989','Tanswell1989'}; verstraete, noe excluded
fileNames = {'Tanswell1989'};
N_dataSet = length(fileNames);
% ----- Model parameters --------------------------------------------------
[k12,k21,k10,kel_PLG,kel_PLS,kel_FBG,kel_AP,kel_MG,kel_PAI] = elimination_para();
[KM_PLG,k_PLG_cat,k_AP_f,k_AP_r,k_AP_cat,KM_FBG,k_FBG_cat,kMG,kPAI] = kinetic_para();
[C_tPA0,C_PLG0,C_PLS0,C_FBG0,C_AP0,C_MG0,C_PAI0] = initConc_values();
% -------------------------------------------------------------------------

for nn = 1:N_dataSet
    % Load data set containing variables...check MAT-files for their names
    clear *patient* dose* duration*
    load([fileNames{nn},'.mat']);
    fprintf('Data set from %s\n',fileNames{nn});
    % Loop over the patient groups 
    for ii = 1:N_patient % N_patient: No of patient groups in each data set
        fprintf('--- Patient group %d out of %d\n',ii,N_patient);
        varName_dose = ['dose_0',num2str(ii)];
        varName_infTime = ['duration_0',num2str(ii)];
        % --- INPUT Variables ---------------------------------------------
        % tPAdose, total tPA dose in [mg]
        % N_stage, Number of infusion stage
        % t_infusion, infusion times [s]
        % t_delay, time delay between infusion modes [s]
        % perc_infusion, percentage of each infusion of total amount [-]
        % -----------------------------------------------------------------
        % tPA dose
        tPAdose_each = eval(varName_dose);
        if length(tPAdose_each)>1
            tPAdose = sum(tPAdose_each);
            perc_infusion = tPAdose_each/tPAdose;
        else
            tPAdose = tPAdose_each;
            perc_infusion = 1;
        end
        % infusion & delay times
        t_infusion = 60*eval(varName_infTime);
        N_stage = length(t_infusion);
        if N_stage>1
            t_delay = zeros(1,N_stage-1);
        else
            t_delay = [];
        end
        
        % Concentration data
        tPAexpData = [];	PLGexpData = [];
        PLSexpData = [];	FBGexpData = [];
        APexpData = [];     MGexpData = [];     PAIexpData = [];  
        maxTime = 0; % [min]
        for rr=1:length(conc) % conc: type of data, what concentration?
            concName = conc{rr};
            switch concName
                case 'tPA'; indRes = 1; 
                    tPAexpData = eval([concName,'_patient_0',num2str(ii)]);
                    maxTime = max([maxTime,tPAexpData(end,1)]);
                case 'PLG'; indRes = 3; 
                    PLGexpData = eval([concName,'_patient_0',num2str(ii)]);
                    maxTime = max([maxTime,PLGexpData(end,1)]);
                case 'FBG'; indRes = 5; 
                    FBGexpData = eval([concName,'_patient_0',num2str(ii)]);
                    maxTime = max([maxTime,FBGexpData(end,1)]);
                case 'AP'; indRes = 6; 
                    APexpData = eval([concName,'_patient_0',num2str(ii)]);
                    maxTime = max([maxTime,APexpData(end,1)]);
                otherwise; error('Something wrong with the name of data set!!! Data set, %s',fileNames{nn});
            end
        end
        
        % Dosage regimen
        [tPAdose_array, timePoints, t_ramping] = dosageRegimen_func(tPAdose,N_stage,t_infusion,t_delay,perc_infusion);

        %% Solving the compartmental model
        dt = 1; % [s]
        tfinal = round(maxTime*60); % in [s]
        tspan = 0:dt:max([(sum(timePoints)),tfinal]);       % simulation duration [s]
        Cinit = [C_tPA0, 0, C_PLG0, C_PLS0, C_FBG0, C_AP0, C_MG0, C_PAI0, 0];
        
        % ODE solver
        options = odeset('RelTol',1e-8,'AbsTol',1e-12);
        [tsolve,Csolve] = ode15s(@compartmentalODE_NewFunc,tspan,Cinit,options,tPAdose_array,timePoints,t_ramping);
        allres(nn,ii).tsolve = tsolve;
        allres(nn,ii).Csolve = Csolve;
        
        % ---- Plot results -----------------------------------------------
        % Data to plot infusion profile
        allNames = {'tPA','tPA pp','PLG','PLS','FBG','AP','MG','PAI'};
        IR_tPA = zeros(1,length(tspan));
        for jj=1:length(tspan)
            IR_tPA(jj) = infusionRate_func_general(tspan(jj), tPAdose_array, timePoints, t_ramping);
        end
        figureName = [fileNames{nn}, '-Group', num2str(ii)];
        figureSize = [-30 5 30 15];
        figure('Name',figureName,'Units','centimeters','Position',figureSize);
        subplot(241) % Infusion rate
        plot(tspan/60,1e3*IR_tPA*60,'LineWidth',2); xlabel('Time [min]'); ylabel('Infusion rate [mg/min]');
        
        for ppp=[1,3:(length(Cinit)-1)]
            subN = (1)*(ppp==1) + ppp; 
            subplot(2,4,subN) 
            % Simulation data plotting
            if ppp==1 % tPA only, central &peripheral tPA
                plot(tsolve/60,Csolve(:,1),'LineWidth',2); hold on;
%                 plot(tsolve/60,Csolve(:,2));
                xlabel('Time [min]'); ylabel([allNames{ppp},' Conc. [\muM]']);
%                 legend('Central','Peripheral');
            else
                plot(tsolve/60,Csolve(:,ppp),'LineWidth',2); hold on;
                xlabel('Time [min]'); ylabel([allNames{ppp},' Conc. [\muM]']);
            end
            % Experimental data plotting
            if ppp==1 && ~isempty(tPAexpData)
                plot(tPAexpData(:,1),tPAexpData(:,2),'o','LineWidth',2);
                legend('Model','Exp. data')
            elseif ppp==3 && ~isempty(PLGexpData)
                plot(PLGexpData(:,1),PLGexpData(:,2)*C_PLG0/100,'o','LineWidth',2);
                legend('Model','Exp. data')
            elseif ppp==5 && ~isempty(FBGexpData)
                plot(FBGexpData(:,1),FBGexpData(:,2)*C_FBG0/100,'o','LineWidth',2);
                legend('Model','Exp. data')
            elseif ppp==6 && ~isempty(APexpData) 
                plot(APexpData(:,1),APexpData(:,2)*C_AP0/100,'o','LineWidth',2);
                legend('Model','Exp. data')
            end
            hold off;
            pause(0.5)
        end
    end
end