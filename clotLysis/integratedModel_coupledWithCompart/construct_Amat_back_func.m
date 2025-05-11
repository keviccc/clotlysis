function Amat = construct_Amat_back_func(nx,dt,dx,U,D)
%% Construct the matrix on the left side
Amat = zeros(nx-1);
M11 = 1/dt + 4*U/(2*dx) + 2*D/(dx*dx);
M12 = -U/(2*dx) - D/(dx*dx);
M21 = -4*U/(2*dx) - D/(dx*dx);
M22 = 1/dt + 3*U/(2*dx) + 2*D/(dx*dx);
M23 = -D/(dx*dx);
M31 = U/(2*dx);
repUnit = [M31, M21, M22, M23];
Amat(1,1:2) = [M11, M12];
Amat(2,1:3) = [M21, M22, M23];
indupdate = 1:4;
for ii = 3:(nx-2)
    Amat(ii,indupdate) = repUnit;
    indupdate = indupdate+1;
end
Amat(nx-1,end-2:end) = [M31, M21+D/(3*dx*dx), M22-4*D/(3*dx*dx)];
