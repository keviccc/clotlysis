function ret = intialise_func(initCond,Nseg)
% This function is to assign intial values of variables to each segment 
% that are used in the convection-diffusion-reaction equation.
% It can be used for both clot-free arteries and clotted arteries.
global n_0 N_free

%% Initial concentrations
% --- free phase concentrations -------------------------------------------
if length(initCond)==N_free
    C_tPA_vec = initCond(1)*ones(1,sum(Nseg));        % initial free tPA conc. [microM]
    C_PLG_vec = initCond(2)*ones(1,sum(Nseg));        % initial free PLG conc. [microM]
    C_PLS_vec = initCond(3)*ones(1,sum(Nseg));        % initial free PLS conc. [microM]
    C_FBG_vec = initCond(4)*ones(1,sum(Nseg));        % initial fibrinogen conc. [microM]
    C_AP_vec = initCond(5)*ones(1,sum(Nseg));         % initial free AP conc. [microM]
    C_MG_vec = initCond(6)*ones(1,sum(Nseg));         % initial free MG conc. [microM]
    C_PAI_vec = initCond(7)*ones(1,sum(Nseg));        % initial free PAI conc. [microM]
    C_PLS_AP_vec = initCond(8)*ones(1,sum(Nseg));     % initial free PLS-AP complex conc. [microM]
else
    error('Error: Size of initCond NOT correct!!!');
end
% --- Bound phase concentrations ------------------------------------------
n_tPA_vec = zeros(1,sum(Nseg));
n_PLG_vec = zeros(1,sum(Nseg));
n_PLS_vec = zeros(1,sum(Nseg));
n_PLSlysed_vec = zeros(1,sum(Nseg));
n_tot_vec =  zeros(1,sum(Nseg));
n_tot_vec((Nseg(1)+1):(Nseg(1)+Nseg(2))) = n_0;

%% Return
ret.C_tPA_vec = C_tPA_vec';
ret.C_PLG_vec = C_PLG_vec';
ret.C_PLS_vec = C_PLS_vec';
ret.C_FBG_vec = C_FBG_vec';
ret.C_AP_vec = C_AP_vec';
ret.C_MG_vec = C_MG_vec';
ret.C_PAI_vec = C_PAI_vec';
ret.C_PLS_AP_vec = C_PLS_AP_vec';
ret.n_tPA_vec = n_tPA_vec';
ret.n_PLG_vec = n_PLG_vec';
ret.n_PLS_vec = n_PLS_vec';
ret.n_PLSlysed_vec = n_PLSlysed_vec';
ret.n_tot_vec = n_tot_vec';
