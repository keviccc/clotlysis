function Q = DarcyLaw_func(R)
global dpdx_clot D_MCA viscosity
% D: diameter [m]
% R: resistance [m^-2]
% vis: viscosity [Pa s]
% dpdx: dP/dx [Pa/m]
% flowrate through the porous media (clot) calcuated by Darcy's law

Q = (pi*D_MCA*D_MCA/4)/(R*viscosity)*dpdx_clot; 
