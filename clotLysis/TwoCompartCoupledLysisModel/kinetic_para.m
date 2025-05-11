function[KM_PLG,k_PLG_cat,k_AP_f,k_AP_r,k_AP_cat,KM_FBG,k_FBG_cat,kMG,kPAI] = kinetic_para()

% [k1f,k1r,k2,k3f,k3r,k4,k5f,k5r,k6,kMG,kPAI] = kinetic_para()

%% Kinetic parameters
% Reaction 1. PLG --> PLS by tPA, Michaelis Menten
k1f = 10;               % 1/(micro M*s]
k1r = 280;              % 1/s
k2 = 0.3;               % 1/s
% -------------------------------------
k_PLG_cat = k2;
KM_PLG = (k1r + k2)/k1f;  % micro M
% -------------------------------------

% Reaction 2. FBG-->FDP by PLS, Michaelis-Menten
k5f = 10 ;              % 1/(micro M*s]
k5r = 300;              % 1/s
k6 = 25*10;             % 1/s *****************************************
% -------------------------------------
k_FBG_cat = k6;
KM_FBG = (k5r + k6)/k5f;% micro M
% -------------------------------------

% Reaction 3. PLS inactivation by AP, second order Michaelis-Menten
k_AP_f = 10;
k_AP_r = 0.0021;
k_AP_cat = 0.004;

% Reaction 4. PLS inactivation by MG, second order
kMG = 0.35;         % [1/(microM s)]

% Reaction 5. tPA inhibiation by PAI, second order
kPAI = 37;          % [1/(micro M s)]



% %% Kinetic parameters
% % Reaction 1. PLG --> PLS by tPA, Michaelis Menten
% % kcat_tPA = 0.3;     % [1/s]
% % KM_tPA = 28.03;     % [micro M]
% % kcat_tPA = 0.1;     % [1/s]
% % KM_tPA = 65;     % [micro M]
% k1f = 10;
% k1r = 280;
% k2 = 0.3;
% % Reaction 2. FBG-->FDP by PLS, Michaelis-Menten
% % kcat_PLS = 20*25;      % [1/s]
% % KM_PLS = 32.5;      % [micro M]
% k5f = 10 ;
% k5r = 300;
% k6 = 25;
% % Reaction 3. PLS inactivation by AP, second order Michaelis-Menten
% % kAP = 10;           % [1/(microM s)]
% % kcat_AP = 2e-3;
% % KM_AP = (2e-3 + 2.2e-3)/22;
% k3f = 10;
% k3r = 0.0021;
% k4 = 0.004;
% 
% % Reaction 4. PLS inactivation by MG, second order
% kMG = 0.35;         % [1/(microM s)]
% 
% % Reaction 5. tPA inhibiation by PAI, second order
% kPAI = 37;          % [1/(micro M s)]

