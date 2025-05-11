function Bvec = construct_Bvec_back_func(dt,dx,nx,C_t0,C_x0,D,U,R)
Bvec = zeros(nx-1,1);
Bvec(1) = R(1) + C_t0(1)/dt + 3*U*C_x0/(2*dx) + D*C_x0/(dx*dx);
Bvec(2) = R(2) + C_t0(2)/dt - U*C_x0/(2*dx);
for jj = 3:(nx-1)
    Bvec(jj) =  R(jj) + C_t0(jj)/dt;
end
