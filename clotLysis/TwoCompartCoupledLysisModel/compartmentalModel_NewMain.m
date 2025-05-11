%% ---------------------------------------------------------------------- %
% Two Compartmental model                                                 %
% programmed by B Gu, 18/08/2018, modified 20/12/2018                     %
% ----------------------------------------------------------------------- %

clc; clear all; format long; %close all;
global kel_PLG kel_PLS kel_FBG kel_AP kel_MG kel_PAI k12 k21 k10
% global kcat_tPA KM_tPA kcat_PLS KM_PLS kMG kPAI kcat_AP KM_AP
global C_tPA0 C_PLG0 C_PLS0 C_FBG0 C_AP0 C_MG0 C_PAI0 
% global k1f k1r k2 k3f k3r k4 k5f k5r k6 kMG kPAI
global KM_PLG k_PLG_cat k_AP_f k_AP_r k_AP_cat KM_FBG k_FBG_cat kMG kPAI

%% Properties Parameters
[k12,k21,k10,kel_PLG,kel_PLS,kel_FBG,kel_AP,kel_MG,kel_PAI] = elimination_para();

%% Kinetic parameters
% [kcat_tPA,KM_tPA,kcat_PLS,KM_PLS,kcat_AP,KM_AP,kMG,kPAI] = kinetic_para();
% [k1f,k1r,k2,k3f,k3r,k4,k5f,k5r,k6,kMG,kPAI] = kinetic_para();
[KM_PLG,k_PLG_cat,k_AP_f,k_AP_r,k_AP_cat,KM_FBG,k_FBG_cat,kMG,kPAI] = kinetic_para();

%% initial concentrations in plasma
[C_tPA0,C_PLG0,C_PLS0,C_FBG0,C_AP0,C_MG0,C_PAI0] = initConc_values();

%% Treatment specific parameters
% === INPUT ===============================================================
tPAdose = 0.9*80;                     % total tPA dose in [mg]
N_stage = 2;
t_infusion = [1,60]*60;    % infusion duration [s]
t_delay = [0]*60;           % time delay between infusion modes [s]
perc_infusion = [0.1,0.9];    % percentage of each infusion of total amount [-]
% =========================================================================
[tPAdose_array, timePoints, t_ramping] = dosageRegimen_func(tPAdose,N_stage,t_infusion,t_delay,perc_infusion);

%% Calculate infusion rate
dt = 0.01;%0.075;
% tspan = 0:dt:(sum(timePoints)+60*100);       % simulation duration [s]
% tspan = 0:dt:3600+60+30*60;       % simulation duration [s]
tspan = 0:dt:200*60
for ii=1:length(tspan)
    IR_tPA(ii) = infusionRate_func_general(tspan(ii), tPAdose_array, timePoints, t_ramping);
end

figure(4)
plot(tspan/60,1e3*IR_tPA*60,'LineWidth',2); grid on; hold on;
ax = gca; ax.FontSize = 18;
xlabel('time [min]'); ylabel('Infusion rate [mg/min]');
% dd
%% Solve ODE for the compartmental model
% --- initial conc. [central tPA, pheripheral tPA, PLG, PLS, FBG, AP, MG, PAI] [micro M] - 9 factors
Cintit = [C_tPA0, 0, C_PLG0, C_PLS0, C_FBG0, C_AP0, C_MG0, C_PAI0, 0];
[tsolve,Csolve] = ode23s(@compartmentalODE_NewFunc,tspan,Cintit,[],tPAdose_array,timePoints,t_ramping);
filename = ['compartSol_low','.csv'];
% dlmwrite(filename,[tsolve,Csolve(:,[1,3:4])],'precision','%4.2f')
% 
% dd
%% Plot Sim results
subplotOFF = 2;
if subplotOFF==1
%     figure('Name','tPA');
    figure(2);
%     yyaxis left
    plot(tsolve/60,Csolve(:,1),'LineWidth',2); grid on; hold on; 
    ax = gca; ax.FontSize = 18;
    ylabel('C_{tPA,sys} [{\mu}M]'); xlabel('Time [min]'); 
    legend('1.5:5:3.5','1.5:3.5:5','1:6:3','1:3:6');
    
%     figure(21);
% yyaxis right
%     plot(tsolve/60,Csolve(:,2),'LineWidth',2); grid on; hold on;
%     ylabel('Peripheral C_{tPA} [{\mu}M]'); xlabel('Time [min]'); 
%     ax = gca; ax.FontSize = 18;
%     legend('Low','Medium','High');
    
%     figure('Name','PLG');
    figure(3);
    plot(tsolve/60,Csolve(:,3),'LineWidth',2); grid on; hold on;
    xlabel('Time [min]'); ylabel('C_{PLG,sys} [{\mu}M]');
    ax = gca; ax.FontSize = 18;
%     legend('0:10','1:9','2:8','3:7','5:5');
    legend('1.5:5:3.5','1.5:3.5:5','1:6:3','1:3:6');
%     figure('Name','PLS');
    figure(4);
    plot(tsolve/60,Csolve(:,4),'LineWidth',2); grid on; hold on;
    xlabel('Time [min]'); ylabel('C_{PLS,sys} [{\mu}M]');
    ax = gca; ax.FontSize = 18;
%     legend('0:10','1:9','2:8','3:7','5:5');
legend('1.5:5:3.5','1.5:3.5:5','1:6:3','1:3:6');
%     figure('Name','FBG');
    figure(5);
    plot(tsolve/60,Csolve(:,5),'LineWidth',2); grid on; hold on;
    xlabel('Time [min]'); ylabel('C_{FBG,sys} [{\mu}M]');
    ax = gca; ax.FontSize = 18;
%     legend('0:10','1:9','2:8','3:7','5:5');
legend('1.5:5:3.5','1.5:3.5:5','1:6:3','1:3:6');
%     figure('Name','AP');
    figure(6);
    plot(tsolve/60,Csolve(:,6),'LineWidth',2); grid on; hold on;
    xlabel('Time [min]'); ylabel('C_{AP,sys} [{\mu}M]');
    ax = gca; ax.FontSize = 18;
%     legend('Low','Medium','High');
    legend('1.5:5:3.5','1.5:3.5:5','1:6:3','1:3:6');
%     figure('Name','MG');
    figure(7);
    plot(tsolve/60,Csolve(:,7),'LineWidth',2); grid on; hold on;
    xlabel('Time [min]'); ylabel('C_{MG,sys} [{\mu}M]');
    ax = gca; ax.FontSize = 18;
%     legend('Low','Medium','High');
    legend('1.5:5:3.5','1.5:3.5:5','1:6:3','1:3:6');
%     figure('Name','PAI'); 
    figure(8);
    plot(tsolve/60,Csolve(:,8),'LineWidth',2); grid on; ylim([0,C_PAI0]); hold on;
    xlabel('Time [min]'); ylabel('C_{PAI,sys} [{\mu}M]');
    ax = gca; ax.FontSize = 18;
%     legend('Low','Medium','High');
legend('1.5:5:3.5','1.5:3.5:5','1:6:3','1:3:6');
elseif subplotOFF==2
    figure('Name','Results','Units','centimeter','Position',[1,1,30,10]);
    subplot(241);
    plot(tsolve/60,Csolve(:,1)); 
    xlabel('Time [min]'); ylabel('C_{tPA} [{\mu}M]');
    subplot(242);
    plot(tsolve/60,Csolve(:,2)); 
    xlabel('Time [min]'); ylabel('C_{tPA,peri} [{\mu}M]');
    subplot(243);
    plot(tsolve/60,Csolve(:,3));
    xlabel('Time [min]'); ylabel('C_{PLG} [{\mu}M]');
    subplot(244);
    plot(tsolve/60,Csolve(:,4));
    xlabel('Time [min]'); ylabel('C_{PLS} [{\mu}M]');
    subplot(245);
    plot(tsolve/60,Csolve(:,5));
    xlabel('Time [min]'); ylabel('C_{FBG} [{\mu}M]');
    subplot(246);
    plot(tsolve/60,Csolve(:,6));
    xlabel('Time [min]'); ylabel('C_{AP} [{\mu}M]');
    subplot(247);
    plot(tsolve/60,Csolve(:,7));
    xlabel('Time [min]'); ylabel('C_{MG} [{\mu}M]');
    subplot(248);
    plot(tsolve/60,Csolve(:,8)); ylim([0,C_PAI0])
    xlabel('Time [min]'); ylabel('C_{PAI} [{\mu}M]');
end



