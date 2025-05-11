function dy = compartmentalODE_func(t,y,tPAdose_tot,perc_bolus,t0,tb,td,tc,pb,pc)
global Mw_tPA V_plasma
global kel_PLG kel_PLS kel_FBG kel_AP kel_MG kel_PAI k12 k21 k10
global kcat_tPA KM_tPA kcat_PLS KM_PLS kAP kMG kPAI kcat_AP KM_AP
global C_tPA0 C_PLG0 C_PLS0 C_FBG0 C_AP0 C_MG0 C_PAI0
% C0: initial conc. array in the order of tPA, PLG, PLS, FBG, AP, MG and PAI[micron M]
% k_el: decay constants, array [1/s]
% kinpara: structure variable to contain kinetic parameters [variable unit]
% specpara: structure variable to contain treatment specific parameters
% y: concentrations to be solved in micro M
C_tPA = y(1);
C_tPA_peri = y(2);
C_PLG = y(3);
C_PLS = y(4);
C_FBG = y(5);
C_AP = y(6);
C_MG = y(7);
C_PAI = y(8);

% generation terms, [micro M /s]
C0 = [C_tPA0, C_PLG0,C_PLS0,C_FBG0, C_AP0, C_MG0, C_PAI0];
S_init = [k10, kel_PLG, kel_PLS, kel_FBG, kel_AP, kel_MG, kel_PAI].*C0;
% S_init = 0.5*[k10, kel_PLG, kel_PLS, kel_FBG, kel_AP, kel_MG, kel_PAI].*(C0([1,3:end])-y([1,3:end]));


% IR: infusion rate [g/s]
IR_tPA = infusionRate_func_modified(perc_bolus,t,tPAdose_tot,t0,tb,td,tc,pb,pc);

fprintf('time = %4.2f\n',t)

%% ODEs
% --- Michaelis Metnen kinetics
PLG2PLSbyTPA = (kcat_tPA*C_tPA*C_PLG/(KM_tPA + C_PLG))*(C_tPA>C_tPA0);
FBG2FDPbyPLS = (kcat_PLS*C_PLS*C_FBG/(KM_PLS + C_FBG))*(C_PLS>C_PLS0);
% FBG2FDPbyPLS = 8*C_PLS*C_FBG;
% --- tPA inhibition by PAI
tPA_inhibited = kPAI*C_tPA*C_PAI*(C_tPA>C_tPA0);
% --- PLS inactivation by AP and MG
PLS_inactByAP = (kcat_AP*C_PLS*C_AP/(KM_AP + C_PLS))*(C_PLS>C_PLS0);
% PLS_inactByAP = kAP*C_PLS*C_AP*(C_PLS>C_PLS0)
PLS_inactByMG = kMG*C_PLS*C_MG*(C_PLS>C_PLS0);

% --- Generation term
S = S_init;

% --- component balance
dC_tPA = 1e6*(IR_tPA/Mw_tPA)/V_plasma - k10*C_tPA - k12*C_tPA + k21*C_tPA_peri - tPA_inhibited + S(1);
dC_tPA_peri = k12*C_tPA - k21*C_tPA_peri;
dC_PLG = -kel_PLG*C_PLG - PLG2PLSbyTPA + S(2);
dC_PLS = -kel_PLS*C_PLS + PLG2PLSbyTPA - PLS_inactByAP - PLS_inactByMG;
dC_FBG = -kel_FBG*C_FBG - FBG2FDPbyPLS + S(4);
dC_AP = -kel_AP*C_AP - PLS_inactByAP + S(5);
dC_MG = -kel_MG*C_MG - PLS_inactByMG + S(6);
dC_PAI = -kel_PAI*C_PAI - tPA_inhibited + S(7);
  
%% Return
dy = [dC_tPA; dC_tPA_peri; dC_PLG; dC_PLS; dC_FBG; dC_AP; dC_MG; dC_PAI];





