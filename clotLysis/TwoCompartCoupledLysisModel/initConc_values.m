function [C_tPA0,C_PLG0,C_PLS0,C_FBG0,C_AP0,C_MG0,C_PAI0] = initConc_values()
%% Initial concentrations
C_tPA0 = 0.05e-3;   % tPA [micro M]
C_PLG0 = 2.2;       % PLG [micro M]
C_PLS0 = 0;         % PLS [micro M]
C_FBG0 = 8;        % FBG [micro M]
C_AP0 = 1;          % AP [micro M]
C_MG0 = 3;          % MG [micro M]
C_PAI0 = 22.5/(43e3);      % PAI, between 5/43e3 and 40/43e3 [micro M]
