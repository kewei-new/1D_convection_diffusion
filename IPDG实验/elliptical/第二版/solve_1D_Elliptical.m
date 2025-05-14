function [error2,errors,uh,solution] = solve_1D_Elliptical(mp)
% solve -uxx = sin(x)

% initial
global lambda weight
[lambda,weight] = get_gauss(mp);
get_basis;

% assemble matirces
Matrices = assemble_1D_Matrices_IPDG;

% generate vector
Vector = generate_1D_vector;

% boundary treatment
[Matrices,Vector] = boundary_treatment(Matrices,Vector);

% calculate algebra system
% 设置参数
tol = 1e-6;     % 容差
maxit = size(Matrices,1);   % 最大迭代次数
M = diag(diag(Matrices)); % 简单预处理矩阵（对角线缩放）

% 调用 GMRES
[uh, flag] = gmres(Matrices, Vector, [], tol, maxit, M);

% 检查收敛性
% if flag == 0
%     disp('GMRES 收敛');
% else
%     disp('GMRES 未收敛，尝试调整预处理或参数。');
% end

% size reshape
global dimPk nx
uh = reshape(uh,[dimPk,nx]);

% plot
[solution,X_points] = plot_1D_DG(uh);

% error evaluate
error2 = evaluate_L2_norm_error(solution,X_points);
errors = evaluate_absolute_norm_error(solution,X_points);

% fprintf('%d points errors is %4f \n',nx,errors)


