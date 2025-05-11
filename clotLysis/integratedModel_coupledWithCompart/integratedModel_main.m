%% --------------------------------------------------------------------- %%
% Integrated model for                                                    %
% - 1D blood flow by accounting for only artery length and clot position  %
% - protein transport using convection-diffusion-reaction                 %
% - clot dissolution as fibrinolysis kinetics model                       %
% - clot property estimation such as resistance and porosity              %
% - compartmental model to obtain dynamic changes at inlet                %
% based on the previous model, programmed by B Gu, 17/01/2018             %
% ----------------------------------------------------------------------- %

clc; clear all; format long;

%% Fluid and protiens parameters
propPara.density = 1030;             % density of fluid [kg/m3]
propPara.viscosity = 1;              % dynamic viscosity [Pa s]
propPara.kinVis = propPara.viscosity/propPara.density; % kinematic viscosity [m2/s]
propPara.Dcoeff = 5e-11;             % dispersion coeff of proteins [m2/s]
propPara.Q_ICA = 4.31e-6;            % flowrate at ICA [m3/s]
%% Arterial geometry
propPara.D_ICA = 3.35e-3;            % diameter of ICA [m]
propPara.D_ACA = 2.20e-3;            % diameter of ACA [m]
propPara.D_MCA = 3.00e-3;            % diameter of MCA [m]
propPara.L_ICA = 5*propPara.D_ICA;   % length of ICA [m]
propPara.L_ACA = 5*propPara.D_ACA;   % length of ACA [m]
propPara.L_MCA = 5*propPara.D_MCA;   % length of MCA [m]

%% Information on clot
clot.L_clot = 3e-3;              % length of clot [m]
clot.clotPos = 0.1*propPara.L_MCA;        % clot position, distance from the MCA start [m]
clot.dpdx_clot = 300;              % pressure drop per length in the clot [Pa/m]
% **** revise the units ***
clot.R_f0 = 85;
clot.rho_0 = 3;
clot.rho_fibre = 0.28; 
clot.dr = 6;
clot.dth = 6;
clot.L_M = 0.04444;
clot.N_AV = 6.02E23;
% *************************

%% Kinetic parameters for fibrinolysis
kinPara.ka_tPA = 0.01;      % adsorption constant for tPA, [1/(microM*s)]
kinPara.kd_tPA = 0.0058;    % desorption constant for tPA, [1/s]
kinPara.ka_PLG = 0.1;       % adroption constant for PLG, [1/(microM*s)]
kinPara.kd_PLG = 3.5;       % desorption constant for PLG, [1/s]
kinPara.ka_PLS = 0.1;       % adroption constant for PLS, [1/(microM*s)]
kinPara.kd_PLS = 0.05;      % desorption consntat for PLS, [1/s]
kinPara.k2 = 0.3;           % Michaelis reaction constant, [1/s]
kinPara.KM = 0.16;          % Michaelis reaction constant, [microM]
kinPara.k_AP = 10;          % reaction rate constant of AP inhibition, [1/(microM*s)]
kinPara.k_cat = 2.178;      % Lysis coefficient, [1/s]
kinPara.kf_PLS = kinPara.kd_PLS;    % desoption of PLS from lysed binding site, [1/s]
kinPara.gamma = 10;         % no. of cuts needed for 1 fibrin unit

%% Discretisation for FDM
n_ICA = 100;             % No. of discretised segments in ICA
n_ACA = 50;             % No. of discretised segments in ACA
n_MCA1 = 20;            % No. of discretised segments in MCA region 1 (clot free region proximal to bifurcation)
n_MCA2 = 50;            % No. of discretised segments in MCA region 2 (clot)
n_MCA3 = 50;            % No. of discretised segments in MCA region 3 (clot free region distal to bifurcation)
dx_ICA = propPara.L_ICA/n_ICA;	% ICA segment size [m]
dx_ACA = propPara.L_ACA/n_ACA;  % ACA segment size [m] 
dx_MCA1 = clot.clotPos/n_MCA1;  % MCA segment size 1 [m]
dx_MCA2 = clot.L_clot/n_MCA2;   % MCA segment size 2 [m]
dx_MCA3 = (propPara.L_MCA-clot.L_clot-clot.clotPos)/n_MCA2;    % MCA segment size 3 [m]
% --- vectors of size of each segment
dx_ICA_vec = dx_ICA*ones(n_ICA,1);
dx_ACA_vec = dx_ACA*ones(n_ACA,1);
dx_MCA_vec = [dx_MCA1*ones(n_MCA1,1); dx_MCA2*ones(n_MCA2,1); dx_MCA3*ones(n_MCA3,1)];

%% Variables to be solved
% --- free phase concentration
% Initial concentrations
C_tPA_t0 = 0;      % free tPA conc. [microM]
C_PLG_t0 = 2.2;    % free PLG conc. [microM] 
C_PLS_t0 = 0;      % free PLS conc. [microM]
C_AP_t0 = 1;       % free AP conc. [microM]
C_FBG_t0 = 8;      % fibrinogen conc. [microM]
% Inlet concentrations
C_tPA_x0 = 40e-3;       % Inlet free tPA conc. [microM]
C_PLG_x0 = 2.2;    % Inlet free PLG conc. [microM] 
C_PLS_x0 = 0;      % Inlet free PLS conc. [microM]
C_AP_x0 = 1;       % Inlet free AP conc. [microM]
C_FBG_x0 = 8;      % Inlet fibrinogen conc. [microM]


initCond = [C_tPA_t0,C_PLG_t0,C_PLS_t0,C_AP_t0,C_FBG_t0];
% Nseg = [n_MCA1,n_MCA2,n_MCA3];
whichArt = 'ICA';
initVal = intialise_func(clot,initCond,n_ICA,whichArt);
% Store initial clot properties
clotInit.epsilon_0 = initVal.epsilon_0;
clotInit.n_tot0 = initVal.n_tot0;
clotInit.phi_0 = initVal.phi_0;
clotInit.R_clot0 = initVal.R_clot0;


inletCond = [C_tPA_x0,C_PLG_x0,C_PLS_x0,C_AP_x0,C_FBG_x0];
y0_ode = [initVal.C_tPA_vec, initVal.C_PLG_vec, initVal.C_PLS_vec,...
    initVal.C_AP_vec, initVal.C_FBG_vec, initVal.n_tPA_vec,...
    initVal.n_PLG_vec, initVal.n_PLS_vec, initVal.n_PLSlysed_vec,...
    initVal.n_tot_vec]'; % initial values, size of [(Nseg*10)x1]
dt = 0.00001;
tspan = 0:dt:1;
% dy = convectionDiffusionReaction_ode_func(0,y0_ode,whichArt,propPara,kinPara,clot,n_ICA,dx_ICA_vec,clotInit,inletCond)
% convectionDiffusionReaction_ode_func(t,y,whichArt,propPara,KP,clot,Nseg,dx_vec,clotInit,inletCond)
% [tsolv,Csolv] = ode15s(@convectionDiffusionReaction_ode_func,tspan,y0_ode,[],whichArt,propPara,kinPara,clot,n_ICA,dx_ICA_vec,clotInit,inletCond);

% y_mat = zeros(length(y0_ode),length(tspan));

Nseg = n_ICA;
y_pre = y0_ode;
Csolv = zeros(length(y0_ode),length(tspan));
for tt = 1:length(tspan)
    runtime = tspan(tt);
    dy = convectionDiffusionReaction_ode_func(runtime,y_pre,whichArt,propPara,kinPara,clot,n_ICA,dx_ICA_vec,clotInit,inletCond);
    disp(max(dy))
    y_new = y_pre + dy*dt;
    Csolv(:,tt) = y_new;
    y_pre = y_new;
end
tsolv = tspan;
Csolv = Csolv';

xvec = [0; cumsum(dx_ICA_vec)];
C_tPA = [C_tPA_x0*ones(length(tsolv),1), Csolv(:,1:n_ICA)];
C_PLG = [C_PLG_x0*ones(length(tsolv),1), Csolv(:,(1:n_ICA)+n_ICA)];
C_PLS = [C_PLS_x0*ones(length(tsolv),1), Csolv(:,(1:n_ICA)+2*n_ICA)];
C_AP = [C_AP_x0*ones(length(tsolv),1), Csolv(:,(1:n_ICA)+3*n_ICA)];
C_FBG = [C_FBG_x0*ones(length(tsolv),1), Csolv(:,(1:n_ICA)+4*n_ICA)];
figure(1)
mesh(xvec,tsolv,C_tPA); ylabel('Time [s]'); xlabel('Length [m]'); zlabel('C_{tPA} [\muM]');
figure(2)
subplot(151)
mesh(xvec,tsolv,C_PLG); ylabel('Time [s]'); xlabel('Length [m]'); zlabel('C_{PLG} [\muM]');
subplot(152)
mesh(xvec,tsolv,C_PLS); ylabel('Time [s]'); xlabel('Length [m]'); zlabel('C_{PLS} [\muM]');
subplot(153)
mesh(xvec,tsolv,C_AP); ylabel('Time [s]'); xlabel('Length [m]'); zlabel('C_{AP} [\muM]');
subplot(154)
mesh(xvec,tsolv,C_FBG); ylabel('Time [s]'); xlabel('Length [m]'); zlabel('C_{FBG} [\muM]');













