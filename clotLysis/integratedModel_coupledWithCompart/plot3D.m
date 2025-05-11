close all;

%% Free phase
figure('Name','CtPA');
mesh(tvec(1:finishTimeInd)/60,xvec(2:end)*1e3,ret.C_tPA(2:end,1:finishTimeInd)); xlabel('Time [min]'); ylabel('Length [mm]'); 
title('C_{tPA} [\muM]');
view([90,-90]); colorbar; xlim([0,tvec(finishTimeInd)/60]); ylim([0,xvec(end)*1e3]);
ax = gca; ax.FontSize = 18;

figure('Name','CPLG');
mesh(tvec(1:finishTimeInd)/60,xvec(2:end)*1e3,ret.C_PLG(2:end,1:finishTimeInd)); xlabel('Time [min]'); ylabel('Length [mm]'); 
title('C_{PLG} [\muM]');
view([90,-90]); colorbar; xlim([0,tvec(finishTimeInd)/60]); ylim([0,xvec(end)*1e3]);
ax = gca; ax.FontSize = 18;

figure('Name','CPLS');
mesh(tvec(1:finishTimeInd)/60,xvec(2:end)*1e3,ret.C_PLS(2:end,1:finishTimeInd)); xlabel('Time [min]'); ylabel('Length [mm]'); 
title('C_{PLS} [\muM]');
view([90,-90]); colorbar; xlim([0,tvec(finishTimeInd)/60]); ylim([0,xvec(end)*1e3]);
ax = gca; ax.FontSize = 18;

figure('Name','CFBG');
mesh(tvec(1:finishTimeInd)/60,xvec(2:end)*1e3,ret.C_FBG(2:end,1:finishTimeInd)); xlabel('Time [min]'); ylabel('Length [mm]'); 
title('C_{FBG} [\muM]');
view([90,-90]); colorbar; xlim([0,tvec(finishTimeInd)/60]); ylim([0,xvec(end)*1e3]);
ax = gca; ax.FontSize = 18;


%% Bound phase
figure('Name','ntot');
mesh(tvec(1:finishTimeInd)/60,xvec(2:end)*1e3,ret.n_tot(2:end,1:finishTimeInd)); xlabel('Time [min]'); ylabel('Length [mm]'); 
title('n_{tot} [\muM]');
view([90,-90]); colorbar; xlim([0,tvec(finishTimeInd)/60]); ylim([0,xvec(end)*1e3]);
ax = gca; ax.FontSize = 18;

figure('Name','ntPA');
mesh(tvec(1:finishTimeInd)/60,xvec(2:end)*1e3,ret.n_tPA(2:end,1:finishTimeInd)); xlabel('Time [min]'); ylabel('Length [mm]'); 
title('n_{tPA} [\muM]');
view([90,-90]); colorbar; xlim([0,tvec(finishTimeInd)/60]); ylim([0,xvec(end)*1e3]);
ax = gca; ax.FontSize = 18;

figure('Name','nPLG');
mesh(tvec(1:finishTimeInd)/60,xvec(2:end)*1e3,ret.n_PLG(2:end,1:finishTimeInd)); xlabel('Time [min]'); ylabel('Length [mm]'); 
title('n_{PLG} [\muM]');
view([90,-90]); colorbar; xlim([0,tvec(finishTimeInd)/60]); ylim([0,xvec(end)*1e3]);
ax = gca; ax.FontSize = 18;

figure('Name','nPLS');
mesh(tvec(1:finishTimeInd)/60,xvec(2:end)*1e3,ret.n_PLS(2:end,1:finishTimeInd)); xlabel('Time [min]'); ylabel('Length [mm]'); 
title('n_{PLS} [\muM]');
view([90,-90]); colorbar; xlim([0,tvec(finishTimeInd)/60]); ylim([0,xvec(end)*1e3]);
ax = gca; ax.FontSize = 18;

% figure('Name','Extent of Lysis');
% EL = 1 - ret.n_tot(2:end,1:finishTimeInd)/n_0;
% mesh(tvec(1:finishTimeInd)/60,xvec(2:end)*1e3,EL*100); xlabel('Time [min]'); ylabel('Length [mm]'); 
% title('E_{L} [%]');
% view([90,-90]); colorbar; xlim([0,tvec(finishTimeInd)/60]); ylim([0,xvec(end)*1e3]);
% ax = gca; ax.FontSize = 18;
% 
% figure('Name','Resistance');
% mesh(tvec(1:finishTimeInd)/60,xvec(2:end)*1e3,ret.Rclot_vec(:,1:finishTimeInd)); xlabel('Time [min]'); ylabel('Length [mm]'); 
% title('Clot Resistance [m^{-2}]');
% view([90,-90]); colorbar; xlim([0,tvec(finishTimeInd)/60]); ylim([0,xvec(end)*1e3]);
% ax = gca; ax.FontSize = 18;