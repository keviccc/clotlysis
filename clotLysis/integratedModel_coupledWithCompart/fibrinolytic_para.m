function [ka_tPA,kd_tPA,ka_PLG,kd_PLG,ka_PLS,kd_PLS,k_MM,KM,k_deg,gamma] = fibrinolytic_para()

%% Kinetic parameters for fibrinolysis
% adsoprtion and desorption of fibrinolytic proteins
ka_tPA = 0.01;      % adsorption constant for tPA, [1/(microM*s)]
kd_tPA = 0.0058;    % desorption constant for tPA, [1/s]
ka_PLG = 0.1;       % adroption constant for PLG, [1/(microM*s)]
kd_PLG = 3.8;       % desorption constant for PLG, [1/s]
ka_PLS = 0.1;       % adroption constant for PLS, [1/(microM*s)]
kd_PLS = 0.05;      % desorption consntat for PLS, [1/s]
% kf_PLS = kd_PLS;    % desoption of PLS from lysed binding site, [1/s]

% lysis related parameters
k_MM = 0.3;         % Michaelis reaction constant, [1/s]
KM = 0.19;          % Michaelis reaction constant, [microM]
k_deg = 2.178;      % Lysis coefficient, [1/s]
gamma = 10;         % no. of cuts needed for 1 fibrin unit
