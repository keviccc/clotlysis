function C_t_new = matrix_cal_func(R_pre,C_inlet,C_pre,Amat,dt,dx,nx,D,U)
% Find the correct ghost value to make dCdx(end) = 0
C_nxp1_guess = C_pre(end);  % initial guess
dCdx_end = 1;                   % to enter the while-loop
tol = 1e-8;                    % tolerance
iter = 0;                       % iteration No.
while(abs(dCdx_end)>tol)
    iter = iter+1;
    if iter==300; tol = tol*100; fprintf('\nIteration No. 300\n'); end
    if iter==800; tol = tol*100; fprintf('\nIteration No. 800\n'); end
    if iter>800; disp([dCdx_end,C_nxp1_guess]); pause; end
    % construct the vector on the right hand side
    Bvec_base = construct_Bvec_func(dt,dx,nx,C_pre,C_inlet,D,U,C_nxp1_guess);
    
    Bvec = Bvec_base + R_pre;
    
    % Matrix inversion: C_vec (now) = inverse(A)*B
    C_t_new = Amat\Bvec;
    
    % Check any sub zero values and correct them.
    ind_zero = find(C_t_new<0);
    if ~isempty(ind_zero)
%         fprintf('C below zero: index %d\n', ind_zero);
%         warning('Warning! C found below zero!');
%         pause
    end
    C_t_new(ind_zero) = 0;
    
    % check dCdx(end)=0 (outlet boundary condition)
    dCdx_end = (3*C_t_new(nx)-4*C_t_new(nx-1)+C_t_new(nx-2))/(2*dx);
%     dCdx_end = (3*C_nxp1_guess-4*C_t_new(nx)+C_t_new(nx-1))/(2*dx)
    C_nxp1_guess
    C_t_new(nx)
    C_t_new(nx-1)
    pause
%     dCdx_end = (C_nxp1_guess-C_t_new(nx-1))/(2*dx);
    % update C(nx+1)
    C_nxp1_guess = C_nxp1_guess + 0.1*dCdx_end*(2*dx);
    if C_nxp1_guess<0; C_nxp1_guess = 1e-16; end
%             disp([dCdx_end,C_nxp1_guess])
%             pause
end
if abs(C_nxp1_guess-C_t_new(end-1))>1e-12
fprintf('error %e\n',C_nxp1_guess-C_t_new(end-1)); pause
end

end