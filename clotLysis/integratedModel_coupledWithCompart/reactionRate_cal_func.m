function SourceTerm = reactionRate_cal_func(C0,n0,epsilon_clot)
global ka_tPA kd_tPA ka_PLG kd_PLG ka_PLS kd_PLS k_MM KM k_deg gamma
global KM_PLG k_PLG_cat k_AP_f k_AP_r k_AP_cat KM_FBG k_FBG_cat kMG kPAI
global plasmaReactOn

% C0: matrix [nx by 8]
% selectComp: name of component in String
%% Convect the matrix to vectors (nx by 1)
% free phase
C_tPA_vec = C0(:,1);
C_PLG_vec = C0(:,2);
C_PLS_vec = C0(:,3);
C_FBG_vec = C0(:,4);
C_AP_vec = C0(:,5);
C_MG_vec = C0(:,6);
C_PAI_vec = C0(:,7);
C_PLS_AP_vec = C0(:,8);

% bound phase
n_tPA_vec = n0(:,1);
n_PLG_vec = n0(:,2);
n_PLS_vec = n0(:,3);
n_PLSlysed_vec = n0(:,4);
n_tot_vec = n0(:,5);
n_free_vec = n_tot_vec - (n_tPA_vec+n_PLG_vec+n_PLS_vec);

%% Reactions in the clot
% --- fibrinolytic reactions within the clot
% Each reaction is of a size of [nx,1]
r1 = ka_tPA.*C_tPA_vec.*n_free_vec - kd_tPA*n_tPA_vec;
r2 = ka_PLG*C_PLG_vec.*n_free_vec - kd_PLG*n_PLG_vec;
r3 = ka_PLS*C_PLS_vec.*n_free_vec - kd_PLS*n_PLS_vec;
r4 = k_MM*n_tPA_vec.*n_PLG_vec./( KM*(1-epsilon_clot) + n_PLG_vec );
r5 = k_deg*gamma*n_PLS_vec;
r6 = kd_PLS.*n_PLSlysed_vec;
if sum(isnan(r4))~=0
    r4(isnan(r4)) = 0;
end
% --- Bound phase species reaction rates ----------------------------------
R_tPA_S = r1;
R_PLG_S = r2 - r4;
R_PLS_S = r3 + r4 - r5;
R_PLS_Slysed = r5 - r6;
R_tot = -r5;
% -------------------------------------------------------------------------
% --- Free phase species reaction rates -----------------------------------
% tPA, PLG and PLS react with the fibre.
% AP, FBG, MG, PAI and intermediate complex do not react with the fibre.
R_tPA_wtClot = -r1./epsilon_clot;
R_PLG_wtClot = -r2./epsilon_clot;
R_PLS_wtClot = (-r3 - r5 + r6)./epsilon_clot;

% -------------------------------------------------------------------------
%% Reactions in the plasma
% --- fibrinolytic reactions in PLASMA!!!!!!!
% --- Reaction terms: Written for forward reactions
% 1: tPA + PLG <--> tPA*PLG --> PLS + tPA
R1 = k_PLG_cat*C_tPA_vec.*C_PLG_vec./( KM_PLG + C_PLG_vec );
% 2: PLS + FBG <--> PLS*FBG --> PLS + FDP
R2 = k_FBG_cat*C_PLS_vec.*C_FBG_vec./( KM_FBG + C_FBG_vec );
% 3: PLS + AP <--> PLS*AP
R3 = k_AP_f*C_PLS_vec.*C_AP_vec - k_AP_r*C_PLS_AP_vec;
% 4: PLS*AP --> inactive
R4 = k_AP_cat*C_PLS_AP_vec;
% 5: PLS + MG --> inactive
R5 = kMG*C_PLS_vec.*C_MG_vec;
% 6: tPA + PAI --> inactive
R6 = kPAI*C_tPA_vec.*C_PAI_vec;


% --- Free phase species reaction rates -----------------------------------
% Only free phase species are involved in the reactions, no bound phase.
% tPA, PLG and PLS reactions occur in both the clot and plasma.
R_tPA_plasma = -R6;
R_PLG_plasma = -R1;
R_PLS_plasma = R1 - R3 - R5;
R_FBG_plasma = -R2;
R_AP_plasma = -R3;
R_MG_plasma = -R5;
R_PAI_plasma = -R6;
R_PLS_AP_plasma = R3 - R4;


%% reaction terms for each C_xx and n_xx
% --- Bound phase species
% Defined above.
% --- Free phase species
R_tPA = R_tPA_wtClot;
R_PLG = R_PLG_wtClot;
R_PLS = R_PLS_wtClot;
R_FBG = zeros(size(C0,1),1);
R_AP = zeros(size(C0,1),1);
R_MG = zeros(size(C0,1),1);
R_PAI = zeros(size(C0,1),1);
R_PLS_AP = zeros(size(C0,1),1);
if plasmaReactOn
    R_tPA = R_tPA + R_tPA_plasma.*(n_tot_vec==0);
    R_PLG = R_PLG + R_PLG_plasma.*(n_tot_vec==0);
    R_PLS = R_PLS + R_PLS_plasma.*(n_tot_vec==0);
    R_FBG = R_FBG + R_FBG_plasma.*(n_tot_vec==0);
    R_AP = R_AP + R_AP_plasma.*(n_tot_vec==0);
    R_MG = R_MG + R_MG_plasma.*(n_tot_vec==0);
    R_PAI = R_PAI + R_PAI_plasma.*(n_tot_vec==0);
    R_PLS_AP = R_PLS_AP + R_PLS_AP_plasma.*(n_tot_vec==0);    
end
% pause
%% return
SourceTerm = [R_tPA, R_PLG, R_PLS, R_FBG, R_AP, R_MG, R_PAI, R_PLS_AP,...
    R_tPA_S, R_PLG_S, R_PLS_S, R_PLS_Slysed, R_tot];

end



























