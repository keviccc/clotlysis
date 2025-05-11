function dy = compartmentalODE_NewFunc(t,y,tPAdose_array, timePoints, t_ramping)
global Mw_tPA V_plasma
global kel_PLG kel_PLS kel_FBG kel_AP kel_MG kel_PAI k12 k21 k10
% global kcat_tPA KM_tPA kcat_PLS KM_PLS kAP kMG kPAI kcat_AP KM_AP
% global k1f k1r k2 k3f k3r k4 k5f k5r k6 kMG kPAI
global KM_PLG k_PLG_cat k_AP_f k_AP_r k_AP_cat KM_FBG k_FBG_cat kMG kPAI
global C_tPA0 C_PLG0 C_PLS0 C_FBG0 C_AP0 C_MG0 C_PAI0
% C0: initial conc. array in the order of tPA, PLG, PLS, FBG, AP, MG and PAI[micron M]
% k_el: decay constants, array [1/s]
% kinpara: structure variable to contain kinetic parameters [variable unit]
% specpara: structure variable to contain treatment specific parameters
% y: concentrations to be solved in micro M

% Nonphysical values
% if sum(y<0)>0
%     y(y<0) = 0;
% end
C_tPA = y(1);
C_tPA_peri = y(2);
C_PLG = y(3);
C_PLS = y(4);
C_FBG = y(5);
C_AP = y(6);
C_MG = y(7);
C_PAI = y(8);
% C_tPA_PLG = y(9);
C_PLS_AP = y(9);
% C_PLS_FBG = y(11);

% generation terms, [micro M /s]
C0 = [C_tPA0, C_PLG0,C_PLS0,C_FBG0, C_AP0, C_MG0, C_PAI0];
S_init = [k10, kel_PLG, kel_PLS, kel_FBG, kel_AP, kel_MG, kel_PAI].*C0;
% S_init = 0.5*[k10, kel_PLG, kel_PLS, kel_FBG, kel_AP, kel_MG, kel_PAI].*(C0([1,3:end])-y([1,3:end]));


% IR: infusion rate [g/s]
IR_tPA = infusionRate_func_general(t, tPAdose_array, timePoints, t_ramping);
% timeMin = t/60;
% if mod(t,60)==0
%     fprintf('time = %4.2f [min]\n',timeMin)
% end
%% ODEs
% --- Generation term
S = S_init;

% --- Reaction terms: Written for forward reactions
% 1: tPA + PLG <--> tPA*PLG --> PLS + tPA
R1 = k_PLG_cat*C_tPA.*C_PLG./( KM_PLG + C_PLG );
% 2: PLS + FBG <--> PLS*FBG --> PLS + FDP
R2 = k_FBG_cat*C_PLS.*C_FBG./( KM_FBG + C_FBG );
% 3: PLS + AP <--> PLS*AP
R3 = k_AP_f*C_PLS.*C_AP - k_AP_r*C_PLS_AP;
% 4: PLS*AP --> inactive
R4 = k_AP_cat*C_PLS_AP;
% 5: PLS + MG --> inactive
R5 = kMG*C_PLS*C_MG;
% 6: tPA + PAI --> inactive
R6 = kPAI*C_tPA*C_PAI;

% For each component
R_tPA = -R6;
R_PLG = -R1;
R_PLS = R1 - R3 - R5;
R_FBG = -R2;
R_AP = -R3;
R_MG = -R5;
R_PAI = -R6;
R_PLS_AP = R3 - R4;

% Component balance
dC_tPA = 1e6*(IR_tPA/Mw_tPA)/V_plasma - k10*C_tPA - k12*C_tPA + k21*C_tPA_peri + R_tPA  + S(1);
dC_tPA_peri = k12*C_tPA - k21*C_tPA_peri + 0;
dC_PLG = -kel_PLG*C_PLG + R_PLG  + S(2);
dC_PLS = -kel_PLS*C_PLS + R_PLS;
dC_FBG = -kel_FBG*C_FBG + R_FBG + S(4);
dC_AP = -kel_AP*C_AP + R_AP + S(5);
dC_MG = -kel_MG*C_MG + R_MG + S(6);
dC_PAI = -kel_PAI*C_PAI + R_PAI + S(7);
dC_PLS_AP = R_PLS_AP;

%% Return
dy = [dC_tPA; dC_tPA_peri; dC_PLG; dC_PLS; dC_FBG; dC_AP; dC_MG; dC_PAI; dC_PLS_AP];





