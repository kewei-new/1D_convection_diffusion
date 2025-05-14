function uh = solve_1D_elliptical(Gauss_coefficient,inv_mass,h,mid_points,mp,ng)

% assemble matrices
Matrices = assemble_1D_Matrices_IPDG();

% generate vector
Vector = generate_1D_Vector();

% to solve it
uh = Matrices\Vector;

% pot

