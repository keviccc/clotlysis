function ret = convectionDiffusionReaction_mat_timeVaryingInlet_func(nx,nt,dx,dt,Cinit_mat,ninit_mat,C_x0,clotInit)
global R_clot_0 D_MCA epsilon_0 Dcoeff finishTimeInd
% clotInit: structural variable containing clot properties
% Cinit_mat: initial values = nx by No of free components (10), 
% ninit_mat: initial bound values = nx by No of bound components (5), 
% C_x0: inlet values, [1 x 8]
% Convection-diffusion Free phase only
% : C_tPA, C_PLG, C_PLS, C_FBG, C_AP, C_MG, C_PAI, C_PLS_AP
% Bound phase: n_tPA, n_PLG, n_PLS, n_PLSlysed, n_tot
compType = {'C_tPA', 'C_PLG', 'C_PLS', 'C_FBG', 'C_AP', 'C_MG', 'C_PAI', 'C_PLS_AP'};
boundType = {'n_tPA', 'n_PLG', 'n_PLS', 'n_PLSlysed', 'n_tot'};
%% Preallocation
% --- initial conditions
C_mat_tot = zeros(nx,nt+1,length(compType)); % 3D matrix, nx by nt+1 by 5
for comp = 1:length(compType)
    C_mat_tot(:,1,comp) = Cinit_mat(:,comp);
end
n_mat_tot = zeros(nx,nt+1,length(boundType));       % 3D matrix, nx by nt+1 by 5
for bound = 1:length(boundType)
    n_mat_tot(:,1,bound) = ninit_mat(:,bound);
end

% --- clot properties
Rclot_vec = [clotInit.R_clot_vec, zeros(nx,nt)];        % raw values
Rclot_app_vec = [clotInit.R_clot_vec, zeros(nx,nt)];    % after EL_crit
volFrac_vec = [clotInit.volFrac_vec, zeros(nx,nt)];     % raw values
volFrac_app_vec = [clotInit.volFrac_vec, zeros(nx,nt)]; % after EL_crit
epsilon_vec = 1-volFrac_vec;            % raw values
epsilon_app_vec = 1-volFrac_app_vec;    % after EL_crit
indClot = find(clotInit.R_clot_vec>0);              % where the clot is, indices
NindClot = [length(indClot), zeros(1,nt)];          % number of discretised segments with clot

%% Initial velocity at t=0, (tt=1)
Q_MCA0 = DarcyLaw_func(R_clot_0);
U = Q_MCA0/(pi*D_MCA*D_MCA/4)/epsilon_0;
fprintf('Initial Q_MCA = %3.3e [m3/s], U_MCA = %3.3e [m/s]\n',Q_MCA0,U);
fprintf('Courant number U*dt/dx = %3.2e\n',U*dt/dx);
% fprintf('Press any key to continue...\n');
pause(1)
%% Integrate over time
% initial clot size
C_mat = Cinit_mat;  % initial conc. for all components
n_mat = ninit_mat;  % initial bound phase conc.
for tt = 2:(nt+1)
    if mod((tt-1),(1/dt))==0
        runtime = (tt-1)/(1/dt); 
        clc; fprintf('Time %d sec: U = %2.3e, Rclot = %2.3e, k = %2.3e\n',runtime,U,newClot.R_clot_ave,1/newClot.R_clot_ave); 
    end
    % Source term calculated using (t-dt) information
    SourceTerm = reactionRate_cal_func(C_mat,n_mat,epsilon_app_vec(:,tt-1));
    % ---------------------------------------------------------------------
    % Free phaseFor each component (a total of 10)
    for ss = 1:length(compType)
%         fprintf(' %s ',compType{ss}\n); 
        C_pre = C_mat_tot(:,tt-1,ss);     % at t=t-dt (nx by 1)
        R_pre = SourceTerm(:,ss);       % reaction source term
        
        % -----------------------------------------------------------------
        % Change depending on the time
        C_inlet = C_x0(tt-1,ss);             % inlet conc. (scalar)
        % -----------------------------------------------------------------
        
        % Constuct the matrix on the left hand side
        Amat = construct_Amat_back_func(nx,dt,dx,U,Dcoeff);
        Bvec = construct_Bvec_back_func(dt,dx,nx,C_pre,C_inlet,Dcoeff,U,R_pre);
        C_t_new_temp = Amat\Bvec;
        C_n = (4*C_t_new_temp(end) - C_t_new_temp(end-1))/3;
        C_t_new = [C_t_new_temp; C_n];
        negC_ind = find(C_t_new<0);
        if ~isempty(negC_ind)
%            fprintf('Negative %s conc. appears at time step = %d\n',compType{ss},tt);
%            fprintf('Negative values changed to 0...');
%            pause
%            fprintf('Index: %d\n',negC_ind);
           C_t_new(negC_ind) = 0;
        end

        % Update the matrix for the next time step
        C_mat_tot(:,tt,ss) = C_t_new; 
        
        % Updated the matrix to be used in the next time step
        C_mat(:,ss) = C_t_new;
    end
    % ---------------------------------------------------------------------
    % Bound phase
    % Modify bound phase after applying E_L,crit
    n_tot_ind = 5;
    n_tot_new = n_mat(:,end) + dt*SourceTerm(:, (length(compType) + n_tot_ind));
    % Calculate new clot properties
    newClot =  calculateClotProperties_func(n_tot_new);
    indClot = newClot.indClot;
    indClotRemoved = newClot.ind_trimmedClot;
    % Save modified n_tot
    n_mat_tot(:,tt,n_tot_ind) = newClot.n_tot_app;
    n_mat(:,n_tot_ind) = newClot.n_tot_app;
    % for other bound species
    for bb = 1:(length(boundType)-1)
       n_t_new = n_mat(:,bb) + dt*SourceTerm(:, (bb + length(compType)));
       n_t_new(indClotRemoved) = 0;
       n_mat_tot(:,tt,bb) = n_t_new;
       n_mat(:,bb) = n_t_new;
    end
    
    % ---------------------------------------------------------------------
    if isempty(indClot) % when clot completely dissolved
        finishTimeInd = tt; 
        break
    end
    % ---------------------------------------------------------------------
    % Update flowrate
    Q_MCA_updated = newClot.Q_MCA;
    if isnan(newClot.epsilon_ave); newClot.epsilon_ave = 1; end
    U = (Q_MCA_updated/(pi*D_MCA*D_MCA/4))/newClot.epsilon_ave ;
    if isnan(U)
        fprintf('NaN in U, epsClotAve = %e\n',newClot.epsilon_ave);
        pause
    end
    % ---------------------------------------------------------------------
   
    % Store results
    % --- resistance
    Rclot_app_vec(:,tt) = newClot.R_clot_app;
    Rclot_vec(:,tt) = newClot.R_clot;
    % --- volume fraction of fibres
    volFrac_vec(:,tt) = newClot.volFrac;
    volFrac_app_vec(:,tt) = newClot.volFrac_app;
    % --- porosity
    epsilon_vec(:,tt) = 1-volFrac_vec(:,tt);
    epsilon_app_vec(:,tt) = 1-volFrac_app_vec(:,tt);
    % --- clot location
    NindClot(tt) = newClot.NindClot;
    clotInfo(tt-1) = newClot;
end


%% Return
ret.C_tPA = [C_x0(1)*ones(1,nt+1); C_mat_tot(:,:,1)];
ret.C_PLG = [C_x0(2)*ones(1,nt+1); C_mat_tot(:,:,2)];
ret.C_PLS = [C_x0(3)*ones(1,nt+1); C_mat_tot(:,:,3)];
ret.C_FBG = [C_x0(4)*ones(1,nt+1); C_mat_tot(:,:,4)];
ret.C_AP = [C_x0(5)*ones(1,nt+1); C_mat_tot(:,:,5)];
ret.C_MG = [C_x0(6)*ones(1,nt+1); C_mat_tot(:,:,6)];
ret.C_PAI = [C_x0(7)*ones(1,nt+1); C_mat_tot(:,:,7)];
ret.C_PLS_AP = [C_x0(8)*ones(1,nt+1); C_mat_tot(:,:,8)];
ret.n_tPA = [zeros(1,nt+1); n_mat_tot(:,:,1)];
ret.n_PLG = [zeros(1,nt+1); n_mat_tot(:,:,2)];
ret.n_PLS = [zeros(1,nt+1); n_mat_tot(:,:,3)];
ret.n_PLSlysed = [zeros(1,nt+1); n_mat_tot(:,:,4)];
ret.n_tot = [zeros(1,nt+1); n_mat_tot(:,:,5)];
ret.clotInfo = clotInfo;
ret.Rclot_vec = [Rclot_vec];
ret.Rclot_app_vec = Rclot_app_vec;
ret.epsilon_vec = epsilon_vec;
ret.NindClot = NindClot;
end





