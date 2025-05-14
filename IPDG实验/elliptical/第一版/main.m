% main.m 
clc
clear

left = 0;
right = pi;
n_base = 20;
mp = 3;
nt = 3;

% 
Gauss_coefficient = generate_Gauss_lobatto(mp);
% 
inv_mass = generate_1D_invmass(mp);
% 

error_order = zeros(nt-1,2);
error = zeros(nt,2);
for k = 1:nt

    ng = n_base*2^(k-1);

    % partion
    [h,mid_points]=generate_1D_uniform_grid(left,right,ng);
    
    % Constant Matrices
    generate_1D_constant_Matrices_IPDG;

    % caluate
    [uh,error2,errors]=solve_1D_elliptical(Gauss_coefficient,inv_mass,h,mid_points,mp,ng);
    
    error(k,1) = error2;
    error(k,2) = errors;

end