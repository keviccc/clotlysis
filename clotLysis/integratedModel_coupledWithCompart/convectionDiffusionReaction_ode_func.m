function dy = convectionDiffusionReaction_ode_func(t,y,whichArt,propPara,KP,clot,Nseg,dx_vec,clotInit,inletCond)
dy = zeros(length(y),1);
inletCond = [inletCond, zeros(1,5)]; % for 10 variables
NsegTot = sum(Nseg);
%% Preprocess
% Divide the intial values to several vectors
% Each vector is a size of [1,NsegTot]
ind = 1:NsegTot;
var_names = {'C_tPA','C_PLG','C_PLS','C_AP','C_FBG','n_tPA','n_PLG','n_PLS','n_PLSlysed','n_tot'};
for ii = 1:length(var_names)
   var_each = [var_names{ii},'_vec']; 
   eval([var_each,' = y(ind);']);
   ind = ind + NsegTot;
end
% for inlet values
for kk = 1:length(var_names)
   var_each = [var_names{kk},'_x0']; 
   eval([var_each,' = inletCond(kk);']);
end
n_free_vec = n_tot_vec - (n_tPA_vec + n_PLG_vec + n_PLS_vec);
% assign dC_tPA, ... and so on ...
for jj = 1:length(var_names)
   var_each = ['d',var_names{jj},'_vec']; 
   eval([var_each,' = zeros(NsegTot,1);']);
end

%% Clot properties
ExtLysis = 1 - n_tot_vec./clotInit.n_tot0;
if isnan(ExtLysis)
    ExtLysis = 1;
end
phi_clot = clotInit.phi_0*(1 - ExtLysis);
epsilon_clot = 1 - phi_clot;
R_clot = ((16*(phi_clot.^1.5)).*(1 + 56*phi_clot.^3))*1e18/(0.286*(clot.R_f0*clot.R_f0));

if epsilon_clot==1; 
    epsilon_clot = 0; 
end

%% Flowrates and velocity
if strcmp(whichArt,'ICA')
    Ac = pi*propPara.D_ICA*propPara.D_ICA/4;
    Q = propPara.Q_ICA; 
else
    Ac_MCA = pi*propPara.D_MCA*propPara.D_MCA/4;
    Ac_ACA = pi*propPara.D_ACA*propPara.D_ACA/4;
    
    
end



% Q_MCA = (pi*propPara.D_MCA*propPara.D_MCA/4)/(R_clot*propPara.viscosity)*clot.dpdx_clot; % flowrate through the porous media (clot) calcuated by Darcy's law
% Q_ACA = propPara.Q_ICA - Q_MCA;    % flowrate at ACA [m3/s]
% switch whichArt
%     case 'ICA'; Q = propPara.Q_ICA; Ac = pi*propPara.D_ICA*propPara.D_ICA/4;
%     case 'MCA'; Q = Q_MCA;          Ac = pi*propPara.D_MCA*propPara.D_MCA/4;
%     case 'ACA'; Q = Q_ACA;          Ac = pi*propPara.D_ACA*propPara.D_ACA/4;
%     otherwise; error('Error: Wrong artery name!!!')
% end
U = Q/Ac;
%% Reactions
% Each reaction is of a size of [1,NsegTot]
r1 = KP.ka_tPA.*C_tPA_vec.*n_free_vec - KP.kd_tPA*n_tPA_vec;
r2 = KP.ka_PLG*C_PLG_vec.*n_free_vec - KP.kd_PLG*n_PLG_vec;
r3 = KP.ka_PLS*C_PLS_vec.*n_free_vec - KP.kd_PLS*n_PLS_vec;
r4 = KP.k2*n_tPA_vec.*n_PLG_vec./( KP.KM*(1-epsilon_clot) + n_PLG_vec );
r5 = KP.k_AP*C_PLS_vec.*C_AP_vec;
r6 = KP.k_cat*KP.gamma*n_PLS_vec;
r7 = KP.kf_PLS.*n_PLSlysed_vec; 
% reaction terms for each C_xx and n_xx
R_tPA = -r1;
R_PLG = -r2;
R_PLS = -r3 - r5 + r7;
R_AP = -r5;
R_FBG = zeros(NsegTot,1);
R_tPA_S = r1;
R_PLG_S = r2 - r4;
R_PLS_S = r3 + r4 - r6;
R_PLS_Slysed = r6 - r7;
R_tot = -r6;

%% Component balance in the free phase
t
dC_tPA = convDiffReact_ODE_func(C_tPA_x0,C_tPA_vec,dx_vec,U,propPara.Dcoeff,R_tPA);
% pause
dC_PLG = convDiffReact_ODE_func(C_PLG_x0,C_PLG_vec,dx_vec,U,propPara.Dcoeff,R_PLG);
% pause
dC_PLS = convDiffReact_ODE_func(C_PLS_x0,C_PLS_vec,dx_vec,U,propPara.Dcoeff,R_PLS); 
% pause
dC_AP = convDiffReact_ODE_func(C_AP_x0,C_AP_vec,dx_vec,U,propPara.Dcoeff,R_AP); 
% pause
dC_FBG = convDiffReact_ODE_func(C_FBG_x0,C_FBG_vec,dx_vec,U,propPara.Dcoeff,R_FBG); 
% pause
%% Component balance in the clot
dn_tPA = R_tPA_S;
dn_PLG = R_PLG_S;
dn_PLS = R_PLS_S;
dn_PLSonLysed = R_PLS_Slysed;
dn_tot = R_tot;

%% Return
dy = [dC_tPA; dC_PLG; dC_PLS; dC_AP; dC_FBG;
    dn_tPA; dn_PLG; dn_PLS; dn_PLSonLysed; dn_tot];


end

function dC_vec = convDiffReact_ODE_func(C_x0,C_vec,dx_vec,U,D,R)
dC_vec = zeros(length(C_vec),1);
dCdx = zeros(length(C_vec),1);
d2Cdx2 = zeros(length(C_vec),1);
% for the inlet segment, ii=1
dCdx(1) = ( -C_x0 + C_vec(2) ) / (2*dx_vec(1));
d2Cdx2(1) = ( C_x0 - 2*C_vec(1) + C_vec(2) ) / (dx_vec(1)*dx_vec(1));
dC_vec(1) = -U*dCdx(1) + D*d2Cdx2(1) + R(1);
% for the rest, ii=2:end-1
indrange = 2:(length(dx_vec)-1);
dCdx(indrange) = ( -C_vec(indrange-1) + C_vec(indrange+1) ) ./ (2*dx_vec(indrange)); 
d2Cdx2(indrange) = ( C_vec(indrange-1) - 2*C_vec(indrange) + C_vec(indrange+1) ) ./ (dx_vec(indrange).^2);
dC_vec(indrange) = -U*dCdx(indrange) + D*d2Cdx2(indrange) + R(indrange);
% for the last segment, ii=end
% dCdx(end) = ( 3*C_vec(end) - 4*C_vec(end-1) + C_vec(end-2) ) / (2*dx_vec(end));
dCdx(end) = 0;
d2Cdx2(end) = ( -2*C_vec(end) + 5*C_vec(end-1) - 4*C_vec(end-2) + C_vec(end-3) ) / (dx_vec(end)*dx_vec(end));
dC_vec(end) = -U*dCdx(end) + D*d2Cdx2(end) + R(end);

end



