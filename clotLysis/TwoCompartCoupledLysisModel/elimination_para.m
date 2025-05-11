function [k12,k21,k10,kel_PLG,kel_PLS,kel_FBG,kel_AP,kel_MG,kel_PAI] = elimination_para()
global Vc
%% Parameter specification
% Half life for each species & 2-compartment model parameters for tPA

%% tPA: taken from the literature
% --- For drug, central and peripheral compartment
% =========================================================================
alp = 9.45;     % elimination constant of the alpha phase [1/hr]
bet = 1.039;    % elimination constant of the beta phase [1/hr]
Vc = 0.057;     % central compartmental volume, [L/kg]
V_ss = 0.11;    % volume of distribution at steady state [L/kg]
% =========================================================================
t_hf_alp = log(2)/alp;
t_hf_bet = log(2)/bet;
% fprintf('alpha phase half life of tPA, %2.2f [min]\n',t_hf_alp*60);
% fprintf('beta phase half life of tPA, %2.2f [min]\n',t_hf_bet*60);
% k21 = 1.34;     % [1/hr]
% k10 = alp*bet/k21  % [1/hr]
% k12 = alp + bet - k10 - k21 % [1/hr]
% dd
XX = V_ss/Vc - 1;   % [-]
k21_temp1 = ((alp + bet) + sqrt((alp+bet)^2 - 4*(XX+1)*alp*bet))/(2*(XX+1));
k21_temp2 = ((alp + bet) - sqrt((alp+bet)^2 - 4*(XX+1)*alp*bet))/(2*(XX+1));
k21_hr = k21_temp2;        % in [hr]
k12_hr = XX*k21_hr;           % in [hr]
k10_hr = alp*bet/k21_hr;       % in [hr]
fprintf('Calculated half life of tPA, %2.2f [min]\n',log(2)/k10_hr*60)
k21 = k21_hr/3600; % in [s]
k12 = k12_hr/3600; % in [s]
k10 = k10_hr/3600; % in [s]

%% Half life of others
% =========================================================================
t_PLG = 2.2;    % half life of PLG [days]
t_PLS = 0.1;    % half life of PLS [s]
t_FBG = 4.14;   % half life of fibrinogen [days]
t_AP = 2.64;    % half life of AP [days]
t_MG = 3;       % half life of MG [mins]
t_PAI = 1.5;    % half life of PAI [hr]
% =========================================================================
% --- unit conversion of half life parameters and decay constants
t_HL = [t_PLG*24*3600, t_PLS, t_FBG*24*3600, t_AP*24*3600, t_MG*60, t_PAI*3600]; % s
k_el = log(2)./t_HL;     % decay constants [1/s]
kel_PLG = k_el(1);
kel_PLS = k_el(2);
kel_FBG = k_el(3);
kel_AP = k_el(4);
kel_MG = k_el(5);
kel_PAI = k_el(6);
