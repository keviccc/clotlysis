function ret = solveArteries_func(initCond,inletCond)
% Temporal discretisation, dt and nt same for ICA, MCA and ACA
% inletCond: concentrations at inlet, [1x (8 free + 5 bound) comp]
% Spatial discretisation: even distribution

global dt nt dx nx R_clot_0 phi_0 couplingON

%% Initialise matrices of concentrations
nx_tot = sum(nx);
ret = intialise_func(initCond,nx);
Cinit_mat = [ret.C_tPA_vec, ret.C_PLG_vec, ret.C_PLS_vec, ret.C_FBG_vec...
    ret.C_AP_vec, ret.C_MG_vec, ret.C_PAI_vec, ret.C_PLS_AP_vec]; % initial values, size of [(nx_tot)x8]
ninit_mat = [ret.n_tPA_vec, ret.n_PLG_vec, ret.n_PLS_vec, ret.n_PLSlysed_vec,...
    ret.n_tot_vec];             % initial values, size of [(nx_tot)x5] only for MAC

% Create vectors for clot property variables
% --- preallocation
clotInit.R_clot_vec = zeros(nx_tot,1);   % resistance
clotInit.volFrac_vec = zeros(nx_tot,1);  % fibrin volume fraction

% --- assigning values
indClot = (nx(1)+1):(nx(1)+nx(2));
clotInit.R_clot_vec(indClot) = R_clot_0;
clotInit.volFrac_vec(indClot) = phi_0;

%% Return
if couplingON
    ret = convectionDiffusionReaction_mat_timeVaryingInlet_func(nx_tot,nt,dx,dt,Cinit_mat,ninit_mat,inletCond,clotInit);
else
    ret = convectionDiffusionReaction_mat_func(nx_tot,nt,dx,dt,Cinit_mat,ninit_mat,inletCond,clotInit);
end

