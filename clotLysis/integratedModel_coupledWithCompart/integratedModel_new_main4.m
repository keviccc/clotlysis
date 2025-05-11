%% --------------------------------------------------------------------- %%
% Integrated model for                                                    %
% - 1D blood flow by accounting for only artery length and clot position  %
% - protein transport using convection-diffusion-reaction                 %
% - clot dissolution as fibrinolysis kinetics model                       %
% - clot property estimation such as resistance and porosity              %
% - compartmental model to obtain dynamic changes at inlet                %
% based on the previous model, programmed by B Gu, 17/01/2018             %
%                               modified 31/12/2018                       %
% ----------------------------------------------------------------------- %

clc; clear all; format long; close all;
global density viscosity Dcoeff Q_ICA D_ICA D_ACA D_MCA
global R_f0 epsilon_0 n_0 phi_0 R_clot_0 dpdx_clot EL_crit
global dt dx nt nx
global ka_tPA kd_tPA ka_PLG kd_PLG ka_PLS kd_PLS k_MM KM k_deg gamma
global KM_PLG k_PLG_cat k_AP_f k_AP_r k_AP_cat KM_FBG k_FBG_cat kMG kPAI
global finishTimeInd plasmaReactOn N_free N_bound
tic;
%% Fluid and protiens parameters
density = 1060;             % density of fluid [kg/m3]
kinVis = 3.5e-6;            % kinematic viscosity [m2/s]
viscosity = kinVis*density;              % dynamic viscosity [Pa s]
Dcoeff = 5e-11;             % dispersion coeff of proteins [m2/s]
Q_ICA = 4.31e-6;            % flowrate at ICA [m3/s]

%% Arterial geometry
D_ICA = 3.35e-3;            % diameter of ICA [m]
D_ACA = 2.20e-3;            % diameter of ACA [m]
D_MCA = 3.00e-3;            % diameter of MCA [m]
L_MCA = 15e-3;              % length of MCA [m]
%% Information on clot
L_clot = 5e-3;              % length of clot [m]
clotPos = 0.0e-3;%0.1*propPara.L_MCA;        % clot position, distance from the MCA start [m]
dpdx_clot_mmHg_cm = 60;                    % in mmHg/cm
dpdx_clot = dpdx_clot_mmHg_cm*133.32*1e2;  % in Pa/m

% *************************
R_f0 = 300;             % fibre radius [nm]
epsilon_0 = 0.95;      % initial voidage
N_AV = 6.02E23;        % avogaro's number
rho_0 = 3;
rho_fibre = 0.28;       % 0.28 0.21
dr = 10;
dth = 10;
L_M = 0.04444;
% *************************
[n_0,phi_0,R_clot_0] = clotProperty_func(N_AV,R_f0,rho_0,rho_fibre,dr,dth,L_M,epsilon_0);

%% Kinetic parameters for fibrinolysis
[ka_tPA,kd_tPA,ka_PLG,kd_PLG,ka_PLS,kd_PLS,k_MM,KM,k_deg,gamma] = fibrinolytic_para();
[KM_PLG,k_PLG_cat,k_AP_f,k_AP_r,k_AP_cat,KM_FBG,k_FBG_cat,kMG,kPAI] = kinetic_para();
EL_crit = 0.95;
plasmaReactOn = 1; % 1: plasma reactions to be inluded
%% Discretisation for FDM
% -------------------------------------------------------------------------
dx = 0.01e-3;           % discreised segment size [m]
% -------------------------------------------------------------------------
if clotPos>0
    if dx>L_clot || dx>clotPos
        warning('Selected dxx > Length of clot or distance between clot front and birfucation!!');
        fprintf('The value of dxx is changed to "clotPos" %2.2e...\n',clotPos); pause
        fprintf('Wish to continue with the new dxx value? If yes, hit any key, otherwise, ctrl+C.\n');
    end
end
n_MCA1 = round(clotPos/dx);    % No. of discretised segments in MCA region 1 (clot free region proximal to bifurcation)
n_MCA2 = round(L_clot/dx);     % No. of discretised segments in MCA region 2 (clot)
nx = [n_MCA1, n_MCA2];

%% Temporal discretisation: even distribution
% -------------------------------------------------------------------------
dt = 0.05;
% -------------------------------------------------------------------------
tf = 60*180;
nt = tf/dt;

%% Variables to be solved
N_free = 8;     % Number of free phase species
N_bound = 5;    % Number of bound phase species
% --- free phase concentration
% Initial concentrations
C_tPA_t0 = 0.00005;     % free tPA conc. [microM]
C_PLG_t0 = 2.2;         % free PLG conc. [microM] 
C_PLS_t0 = 0;           % free PLS conc. [microM]
C_FBG_t0 = 8;           % fibrinogen conc. [microM]
C_AP_t0 = 1;            % free AP conc. [microM]
C_MG_t0 = 3;            % free MG conc. [microM]
C_PAI_t0 = 22.5/(43e3); % free PAI conc. [microM]
C_PLS_AP_t0 = 0;        % free PLS-AP complex conc. [microM]
% Inlet concentrations
C_tPA_x0 = 25e-3;       % Inlet free tPA conc. [microM]
C_PLG_x0 = 2.2;         % Inlet free PLG conc. [microM] 
C_PLS_x0 = 0;           % Inlet free PLS conc. [microM]
C_FBG_x0 = 8;           % Inlet fibrinogen conc. [microM]
C_AP_x0 = 1;            % Inlet free AP conc. [microM]
C_MG_x0 = 3;            % Inlet MG conc. [microM]
C_PAI_x0 = 22.5/(43e3); % Inlet PAI conc. [microM]
C_PLS_AP_x0 = 0;        % Inlet PLS-AP complex conc. [microM]

% vectors
initCond = [C_tPA_t0,C_PLG_t0,C_PLS_t0,C_FBG_t0,C_AP_t0,C_MG_t0,...
    C_PAI_t0,C_PLS_AP_t0]; % 8x1
inletCond = [C_tPA_x0,C_PLG_x0,C_PLS_x0,C_FBG_x0,C_AP_x0,C_MG_x0,...
    C_PAI_x0,C_PLS_AP_x0]; % 8x1

ret = solveArteries_func(initCond,inletCond); toc;
if isempty(finishTimeInd)
    fprintf('Clot NOT completely dissolved within Time = %5.2f [min]...\n',tf/60);
else
    fprintf('Clot completly dissolved at Time = %5.2f [min]...\n',(finishTimeInd-1)*dt/60);
end


%% Plot
if isempty(finishTimeInd)
    tvec = 0:dt:dt*nt;
else
    tvec = 0:dt:dt*(finishTimeInd-1);
end

xvec = 0:dx(1):dx(1)*sum(nx);
tselect = 1:100:length(tvec);
% figure('Name','Free Phase','Units','centimeter','Position',[-35,15,45,10]);
% subplot(2,4,1)
% mesh(tvec(tselect)/60,xvec,ret.C_tPA(:,tselect)); xlabel('Time [min]'); ylabel('Length [m]'); title('C_{tPA} [\muM]');
% view([90,-90]); colorbar;
% subplot(2,4,2)
% mesh(tvec(tselect)/60,xvec,ret.C_PLG(:,tselect)); xlabel('Time [min]'); ylabel('Length [m]');  title('C_{PLG} [\muM]');
% view([90,-90]); colorbar;
% subplot(2,4,3)
% mesh(tvec(tselect)/60,xvec,ret.C_PLS(:,tselect)); xlabel('Time [min]'); ylabel('Length [m]');  title('C_{PLS} [\muM]');
% view([90,-90]); colorbar;
% subplot(2,4,4)
% mesh(tvec(tselect)/60,xvec,ret.C_FBG(:,tselect)); xlabel('Time [min]'); ylabel('Length [m]');  title('C_{FBG} [\muM]');
% view([90,-90]); colorbar;
% subplot(2,4,5)
% mesh(tvec(tselect)/60,xvec,ret.C_AP(:,tselect)); xlabel('Time [min]'); ylabel('Length [m]');  title('C_{AP} [\muM]');
% view([90,-90]); colorbar;
% subplot(2,4,6)
% mesh(tvec(tselect)/60,xvec,ret.C_MG(:,tselect)); xlabel('Time [min]'); ylabel('Length [m]');  title('C_{MG} [\muM]');
% view([90,-90]); colorbar;
% subplot(2,4,7)
% mesh(tvec(tselect)/60,xvec,ret.C_PAI(:,tselect)); xlabel('Time [min]'); ylabel('Length [m]');  title('C_{PAI} [\muM]');
% view([90,-90]); colorbar;
% subplot(2,4,8)
% mesh(tvec(tselect)/60,xvec,ret.C_PLS_AP(:,tselect)); xlabel('Time [min]'); ylabel('Length [m]');  title('C_{PLS-AP} [\muM]');
% view([90,-90]); colorbar;
% 
% figure('Name','Bound Phase','Units','centimeter','Position',[-25,1,25,10]);
% subplot(231)
% mesh(tvec(tselect)/60,xvec,ret.n_tot(:,tselect)); xlabel('Time [min]'); ylabel('Length [m]'); title('n_{tot} [\muM]');
% view([90,-90]); colorbar;
% subplot(232)
% mesh(tvec(tselect)/60,xvec,ret.n_tPA(:,tselect)); xlabel('Time [min]'); ylabel('Length [m]'); title('n_{tPA} [\muM]');
% view([90,-90]); colorbar;
% subplot(233)
% mesh(tvec(tselect)/60,xvec,ret.n_PLG(:,tselect)); xlabel('Time [min]'); ylabel('Length [m]'); title('n_{PLG} [\muM]');
% view([90,-90]); colorbar;
% subplot(234)
% mesh(tvec(tselect)/60,xvec,ret.n_PLS(:,tselect)); xlabel('Time [min]'); ylabel('Length [m]'); title('n_{PLS} [\muM]');
% view([90,-90]); colorbar;
% subplot(235)
% mesh(tvec(tselect)/60,xvec,ret.n_PLSlysed(:,tselect)); xlabel('Time [min]'); ylabel('Length [m]'); title('n_{lysed} [\muM]');
% view([90,-90]); colorbar;
% 
% figure('Name','Extent of Lysis');
% EL = 1 - ret.n_tot(:,tselect)/n_0;
% mesh(tvec(tselect)/60,xvec,EL*100); xlabel('Time [min]'); ylabel('Length [m]'); title('E_{L} [%]');
% view([90,-90]); colorbar;
% 
% figure('Name','Resistance');
% mesh(tvec(tselect)/60,xvec(2:end),ret.Rclot_vec(:,tselect)); xlabel('Time [min]'); ylabel('Length [m]'); title('R_{clot} [m^{-2}]');
% view([90,-90]); colorbar;
% 
% % figure('Name','Resistance diff');
% % mesh(tvec(tselect)/60,xvec(2:end),ret.Rclot_vec(:,tselect)-ret.Rclot_app_vec(:,tselect)); xlabel('Time [min]'); ylabel('Length [m]'); title('R_{clot} [m^{-2}]');
% % view([90,-90]); colorbar;
% 
% 
% 
