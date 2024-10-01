function solve_poisson_1D_FEM()
% solve \phi_{xx} = e/{\epsilon}*(n - nd)
% we know n(x,0), and \phi satisfy the dirichlet boundary condition,so we
% could get \phi0, and after in the io(1,2,3,4)th step in Rk, we need this
% proceding to get the right hand side f(x+h/2,t+h/2);

% To reduce the number of calculations which we will use many many times, we
% use the same segment, the matrix had been computed

matrix_R = generate_1D_matrix_FEM_r()
