function obj = paraEst_expData_func(x)
% ----------------------------------------------------------------------- %
% Comparison of compartmental model and experimental data                 %
%  - Plot simulation results and experimental data from the literature    %
%  - Estimation of kinetics parameters                                    %
% programmed by BG, 20/12/2018                                            %
% ----------------------------------------------------------------------- %

global kel_PLG kel_PLS kel_FBG kel_AP kel_MG kel_PAI k12 k21 k10
global k1f k1r k2 k3f k3r k4 k5f k5r k6 kMG kPAI
global C_tPA0 C_PLG0 C_PLS0 C_FBG0 C_AP0 C_MG0 C_PAI0 
global initValueSave fileID

fileNames = {'Verstraete1985','Tienfenbrunn1986','Collen1986',...
    'Noe1987','Tebbe1989','Tanswell1989'};
% fileNames = {'Verstraete1985'};
N_dataSet = length(fileNames);
% ----- Model parameters --------------------------------------------------
[k12,k21,k10,kel_PLG,kel_PLS,kel_FBG,kel_AP,kel_MG,kel_PAI] = elimination_para();
% [k1f,k1r,k2,k3f,k3r,k4,k5f,k5r,k6,kMG,kPAI] = kinetic_para();
[C_tPA0,C_PLG0,C_PLS0,C_FBG0,C_AP0,C_MG0,C_PAI0] = initConc_values();
% -------------------------------------------------------------------------
k1f = x(1);
k1r = x(2);
k2 = x(3);
k3f = x(4);
k3r = x(5);
k4 = x(6);
k5f = x(7);
k5r = x(8);
k6 = x(9);
kMG = x(10);
kPAI = x(11);
fprintf(fileID, 'Discrepancy between initial and current values\n');
for qq=1:length(x)
    fprintf(fileID, ' - Parameter %d: %2.2f%%\n ',qq,(x(qq)-initValueSave(qq))/initValueSave(qq)*100);
end
%% Calculation
err_C_all = []; % errors to build up
errInd = 0; allDataLength = 0;
for nn = 1:N_dataSet
    % Load data set containing variables...check MAT-files for their names
    clear *patient* dose* duration*
    load([fileNames{nn},'.mat']);
    fprintf(fileID, '===== Data set from %s ======\n',fileNames{nn});
    % Loop over the patient groups 
    for ii = 1:N_patient % N_patient: No of patient groups in each data set
        fprintf(fileID, '--- Patient group %d out of %d ---\n',ii,N_patient);
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
        FBGexpData = [];    APexpData = []; 
        maxTime = 0; % [min]
        for rr=1:length(conc) % conc: type of data, what concentration?
            concName = conc{rr};
            fprintf(fileID, '* Extracting Data of %s\n',concName);
            switch concName
                case 'tPA'; indRes = 1;
                    tPAexpData = eval([concName,'_patient_0',num2str(ii)]);
                    maxTime = max([maxTime,tPAexpData(end,1)]);
                    % -----------------------------------------------------
                    allDataLength = allDataLength + size(tPAexpData,1);
                    if tPAexpData(1,1)==0; allDataLength = allDataLength-1; end
                    % -----------------------------------------------------
                case 'PLG'; indRes = 3;
                    PLGexpData = eval([concName,'_patient_0',num2str(ii)]);
                    maxTime = max([maxTime,PLGexpData(end,1)]);
                    % -----------------------------------------------------
                    allDataLength = allDataLength + size(PLGexpData,1);
                    if PLGexpData(1,1)==0; allDataLength = allDataLength-1; end
                    % -----------------------------------------------------
                case 'FBG'; indRes = 5;
                    FBGexpData = eval([concName,'_patient_0',num2str(ii)]);
                    maxTime = max([maxTime,FBGexpData(end,1)]);
                    % -----------------------------------------------------
                    allDataLength = allDataLength + size(FBGexpData,1);
                    if FBGexpData(1,1)==0; allDataLength = allDataLength-1; end
                    % -----------------------------------------------------
                case 'AP'; indRes = 6; 
                    APexpData = eval([concName,'_patient_0',num2str(ii)]);
                    maxTime = max([maxTime,APexpData(end,1)]);
                    % -----------------------------------------------------
                    allDataLength = allDataLength + size(APexpData,1);
                    if APexpData(1,1)==0; allDataLength = allDataLength-1; end
                    % -----------------------------------------------------
                otherwise; error('Something wrong with the name of data set!!! Data set, %s',fileNames{nn});
            end
        end
        
        % Dosage regimen
        [tPAdose_array, timePoints, t_ramping] = dosageRegimen_func(tPAdose,N_stage,t_infusion,t_delay,perc_infusion);

        %% Solving the compartmental model
        dt = 1; % [s]
        tfinal = round(maxTime*60); % in [s]
        tspan = 0:dt:max([(sum(timePoints)),tfinal]);       % simulation duration [s]
        Cinit = [C_tPA0, C_tPA0, C_PLG0, C_PLS0, C_FBG0, C_AP0, C_MG0, C_PAI0, 0, 0, 0];
        
        % ODE solver
        options = odeset('RelTol',1e-8,'AbsTol',1e-10);
        [tsolve,Csolve] = ode15s(@compartmentalODE_NewFunc,tspan,Cinit,options,tPAdose_array,timePoints,t_ramping);
        allres(nn,ii).tsolve = tsolve;
        allres(nn,ii).Csolve = Csolve;
        
        % Calculate errors
        C_tPA_solve = [tsolve, Csolve(:,1)];
        C_PLG_solve = [tsolve, Csolve(:,3)];
        C_FBG_solve = [tsolve, Csolve(:,5)];
        C_AP_solve = [tsolve, Csolve(:,6)];
        
        namesConc = {'tPA','PLG','FBG','AP'};
        for sss = 1:length(namesConc)
            expdataName = [namesConc{sss},'expData'];
            simdataName = ['C_',namesConc{sss},'_solve'];
            expdata = eval(expdataName);
            simdata = eval(simdataName);
            
            if ~isempty(expdata)
                fprintf(fileID, '* Calculating Errors for Data of %s\n',namesConc{sss});
                err_C = calculateError_btw_simNexp(expdata,simdata);
                if ~isempty(err_C); errInd = errInd+1; end
                err_C_all = [err_C_all, err_C];
            end
        end
        
        % ---- Plot results -----------------------------------------------
        % Data to plot infusion profile
        allNames = {'tPA','tPA pp','PLG','PLS','FBG','AP','MG','PAI'};
        IR_tPA = zeros(1,length(tspan));
        for jj=1:length(tspan)
            IR_tPA(jj) = infusionRate_func_general(tspan(jj), tPAdose_array, timePoints, t_ramping);
        end
        figureName = [fileNames{nn}, '-Group', num2str(ii)];
        figureSize = [-30 5 30 15];
%         figure('Name',figureName,'Units','centimeters','Position',figureSize);
%         subplot(241) % Infusion rate
%         plot(tspan/60,1e3*IR_tPA*60); xlabel('Time [min]'); ylabel('Infusion rate [mg/min]');
%         
%         for ppp=[1,3:(length(Cinit)-3)]
%             subN = (1)*(ppp==1) + ppp; 
%             subplot(2,4,subN) 
%             % Simulation data plotting
%             if ppp==1 % tPA only, central &peripheral tPA
%                 plot(tsolve/60,Csolve(:,1)); hold on;
%                 plot(tsolve/60,Csolve(:,2));
%                 xlabel('Time [min]'); ylabel([allNames{ppp},' Conc. [\muM]']);
%                 legend('Central','Peripheral');
%             else
%                 plot(tsolve/60,Csolve(:,ppp)); hold on;
%                 xlabel('Time [min]'); ylabel([allNames{ppp},' Conc. [\muM]']);
%             end
%             % Experimental data plotting
%             if ppp==1 && ~isempty(tPAexpData)
%                 plot(tPAexpData(:,1),tPAexpData(:,2),'o');
%                 legend('Model','Exp. data')
%             elseif ppp==3 && ~isempty(PLGexpData)
%                 plot(PLGexpData(:,1),PLGexpData(:,2)*C_PLG0/100,'o');
%                 legend('Model','Exp. data')
%             elseif ppp==5 && ~isempty(FBGexpData)
%                 plot(FBGexpData(:,1),FBGexpData(:,2)*C_FBG0/100,'o');
%                 legend('Model','Exp. data')
%             elseif ppp==6 && ~isempty(APexpData) 
%                 plot(APexpData(:,1),APexpData(:,2)*C_AP0/100,'o');
%                 legend('Model','Exp. data')
%             end
%             hold off;
%             pause(0.5)
%         end
        % ---- End Plot ---------------------------------------------------
    end
end
% Check data
errSize1 = length(err_C_all);
errSize2 = errInd;
dataSize = allDataLength;
fprintf(fileID, 'Error matrix size = %d or %d and Total data size = %d\n',errSize1,errSize2,dataSize);

%% Objective function
obj = err_C_all*err_C_all';
fprintf(fileID, '==============================================\n');
fprintf(fileID, 'Objective function value = %3.3e\n',obj);
fprintf(fileID, '==============================================\n');


