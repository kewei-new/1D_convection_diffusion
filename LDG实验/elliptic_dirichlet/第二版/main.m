% main.m 

left = 0;
right = pi;
n_base = 20;
mp = 2;
nt = 4;

% 
Gauss_coefficient = generate_Gauss_lobatto(mp);
% 
inv_mass = generate_1D_invmass(mp);
% 

error1_last = 1;
error2_last = 1;
error_order = zeros(nt,2);

for k = 1:nt
    ng = n_base*2^(k-1);
    % 
    [h,mid_points]=generate_1D_uniform_grid(left,right,ng);
    %
    [H,R,L] = generate_1D_LDG_matrix(Gauss_coefficient,mp,inv_mass,h);
    % 
    [M1,M2,Mn] = assemble_LDG_matrix(inv_mass,ng,mp,mid_points,h,H,R,L);
    
    % 
    [phih,qh]=solve_1D_LDG(Gauss_coefficient,inv_mass,M1,M2,Mn,h,mid_points,mp,ng);
    % [phi,E]=test(inv_mass,ng,mid_points,h,mp);
    
    error1 = evaluate_error_absolute('exact_fun',0,phih,h,mid_points,ng,mp);
    error2 = evaluate_error_absolute('exact_fun_x',0,qh,h,mid_points,ng,mp);
    error_order(k,:) = [log(error1_last/error1)/log(2),log(error2_last/error2)/log(2)];
    error1_last = error1;
    error2_last = error2;

end
% 
E = -qh;
subplot(2,1,1)
plot_DG(phih,ng,mp,mid_points);
subplot(2,1,2)
plot_DG(E,ng,mp,mid_points);
error_order