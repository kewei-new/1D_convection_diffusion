% main.m 

left = 0;
right = pi;
ng = 40;
mp = 2;

% 
Gauss_coefficient = generate_Gauss_lobatto(mp);

% 
[h,mid_points]=generate_1D_uniform_grid(left,right,ng);

% 
inv_mass = generate_1D_invmass(mp);

% 
[H,R,L] = generate_1D_LDG_matrix(Gauss_coefficient,mp,inv_mass,h);
[LDG_matrix,f_values,Q_matrix] = assemble_LDG_matrix(inv_mass,ng,mp,mid_points,h,H,R,L);

% 
[phih,qh]=solve_1D_LDG(Gauss_coefficient,LDG_matrix,Q_matrix,f_values,h,mid_points,mp,ng);
% [phi,E]=test(inv_mass,ng,mid_points,h,mp);

error1 = evaluate_error_absolute('exact_fun',0,phih,h,mid_points,ng,mp);
error2 = evaluate_error_absolute('exact_fun_x',0,qh,h,mid_points,ng,mp);
% 
subplot(2,1,1)
plot_DG(phih,ng,mp,mid_points);
subplot(2,1,2)
plot_DG(E,ng,mp,mid_points);