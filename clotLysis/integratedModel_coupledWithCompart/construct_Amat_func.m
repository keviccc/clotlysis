function Amat = construct_Amat_func(nx,dt,dx,U,D)
%% Construct the matrix on the left side
Amat = zeros(nx);
const1 = 1/dt + 2*D/(dx*dx);
const2 = U/(2*dx);
const3 = D/(dx*dx);
const4 = 2*D/(dx*dx);
Amat(1,1:2) = [const1, const2 - const3];
indupdate = 1:3;
for ii = 2:(nx-1)
    Amat(ii,indupdate) = [-const2 - const3, const1, const2 - const3];
    indupdate = indupdate+1;
end
Amat(nx,end-1:end) = [-const2-const3, const1];
% Amat(nx,end-1:end) = [-const4, const1];
end