function ret =  calculateClotProperties_func(ntot)
global R_f0 EL_crit Q_ICA D_ACA D_MCA n_0 phi_0
% ntot: a vector with a size of [nx by 1]
ExtLysis = 1 - ntot/n_0;
if sum(isnan(ExtLysis))~=0
    ExtLysis(isnan(ExtLysis)) = 1;
end

% Update clot properties
volFrac = phi_0*(1-ExtLysis);  % updated fibrin fraction in clot
permb = permeability_func(R_f0,volFrac);
R_clot_cal = 1./permb;  % in [1/m2]
ind_clot_new = find(ExtLysis<EL_crit);      % clot indices
ind_trimmedClot = find(ExtLysis>=EL_crit);  % to be trimmed out as EL>ELcrit
if isempty(ind_clot_new)
     fprintf('No clot left!! \n');
end

% Get the apparent values after applying EL_crit
n_tot_app = zeros(length(ntot),1);
n_tot_app(ind_clot_new) = ntot(ind_clot_new);
R_clot_app = zeros(length(R_clot_cal),1);
R_clot_app(ind_clot_new) =  R_clot_cal(ind_clot_new);
volFrac_app = zeros(length(volFrac),1);
volFrac_app(ind_clot_new) = volFrac(ind_clot_new);

% Average values
R_clot_ave = mean(R_clot_app(ind_clot_new));
ExtLysis_ave =  mean(ExtLysis(ind_clot_new));
volFrac_ave = mean(volFrac(ind_clot_new));
epsilon_ave = mean(1-volFrac(ind_clot_new));

% Average clot resistance & Update flowrate
if isempty(ind_clot_new)
    % clot dissolved and MCA opened
    R_clot_ave = 0;
    Q_MCA = Q_ICA*(D_MCA^2/(D_MCA^2 + D_ACA^2));
else
    Q_MCA =  DarcyLaw_func(R_clot_ave);
end

%% Return
ret.Q_MCA = Q_MCA;
ret.indClot = ind_clot_new;
ret.ind_trimmedClot = ind_trimmedClot;
ret.NindClot = length(ind_clot_new);
% --- raw values
ret.ExtLysis = ExtLysis;
ret.volFrac = volFrac;
ret.R_clot = R_clot_cal;
ret.n_tot_app = n_tot_app;
% --- after applying EL_crit
ret.volFrac_app = volFrac_app;
ret.R_clot_app = R_clot_app;
% --- average
ret.ExtLysis_ave = ExtLysis_ave;
ret.volFrac_ave = volFrac_ave;
ret.epsilon_ave = epsilon_ave;
ret.R_clot_ave = R_clot_ave;




