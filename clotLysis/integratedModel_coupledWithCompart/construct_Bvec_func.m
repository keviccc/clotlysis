function Bvec = construct_Bvec_func(dt,dx,nx,C_t0,C_x0,D,U,C_nxp1_guess)
Bvec = zeros(nx,1);
Bvec(1) = C_t0(1)/dt + U*C_x0/(2*dx) + D*C_x0/(dx*dx);
for jj = 2:(nx-1)
    Bvec(jj) =  C_t0(jj)/dt;
end
Bvec(nx) = C_t0(nx)/dt + C_nxp1_guess*(D/(dx*dx) - U/(2*dx));
% Bvec(nx) = C_t0(nx)/dt + C_nxp1_guess*(D/(dx*dx) - U/(2*dx));
% Bvec_final = Bvec + R_pre;
end