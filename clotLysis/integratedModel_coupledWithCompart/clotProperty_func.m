function [n_0,phi_0,R_clot_0] = clotProperty_func(N_AV,R_f0,rho_0,rho_fibre,dr,dth,L_M,epsilon_0)

% Total number of layers per cross-section
N_layer = R_f0/dr;

% Calculating the number of binding sites per cross-section
it = 1:N_layer;
BS_CS = sum(pi./asin(dth./(2*it*dr)));

Mw = 340;   % kDa = kg/mole

% Fibre length density
LtVt = (rho_0/(rho_fibre*1000))*1E21/(pi*R_f0*R_f0) ;

% Density of cross-sections
CS_V = 2*LtVt/(1E12)*L_M;


%% Return
% Calculate clot properties
n_0 = 2*CS_V * BS_CS*1E6*1E15*(1-epsilon_0)/N_AV;
phi_0 = 1 - epsilon_0;

permb = permeability_func(R_f0,phi_0);
% permb = (1e-18*R_f0*R_f0*4)/(70*(phi_0^1.5)*( 1+52*phi_0^1.5 ));
R_clot_0 = 1/permb;  % in [1/m2]

