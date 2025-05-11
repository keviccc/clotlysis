function ret = modelODE_solve_func(Nseg,initCond)
%% Initial conditions 
% --- free phase concentration
C_tPA = initCond(1);        % initial free tPA conc. [microM]
C_PLG = initCond(2);        % initial free PLG conc. [microM] 
C_PLS = initCond(3);        % initial free PLS conc. [microM]
C_AP = initCond(4);         % initial free AP conc. [microM]
C_FBG = initCond(5);        % initial fibrinogen conc. [microM]
% --- bound phase concentration
n_tPA = initCond(6);        % initial bound tPA conc. [microM]
n_PLG = initCond(7);        % initial bound PLG conc. [microM] 
n_PLS = initCond(8);        % initial bound PLS conc. [microM]
n_PLSlysed = initCond(9);   % initial bound PLS conc. on the lysed binding site [microM]
n_tot = initCond(10);       % total conc. of active binding site [microM]
% --- convert them to vectors
C_tPA_vec = C_tPA*ones(1,Nseg);
C_PLG_vec = C_PLG*ones(1,Nseg);
C_PLS_vec = C_PLS*ones(1,Nseg);
C_AP_vec = C_AP*ones(1,Nseg);
C_FBG_vec = C_FBG*ones(1,Nseg);
% --- bound phase concentration
% These values can be comprised of zeros when they are for clot-free
% arteries.
if n_tot==0
    n_tPA_vec = n_tPA*ones(1,Nseg);
    n_PLG_vec = n_PLG*ones(1,Nseg);
    n_PLS_vec = n_PLS*ones(1,Nseg);
    n_PLSlysed_vec = n_PLSlysed*ones(1,Nseg);
    n_tot_vec = n_tot*ones(1,Nseg);
else
    
    
end
