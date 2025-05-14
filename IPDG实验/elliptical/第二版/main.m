% main.m
% 1D Elliptical Problem With Dirichlet Condition
clear
clc

xa = 0;
xb = 2*pi;

n_base = 20;

global alpha beta dimPk

dimPk = 3;
mp = dimPk-1;

if mp == 1
    alpha = 1;
elseif mp == 2
    alpha = 8;
elseif mp == 3
    alpha = 20;
end
beta = 1;

% 加密次数
nt = 2;

% main
error = zeros(nt,2);
error_order = zeros(nt-1,2);

global nx hx Xm
for k = 1:nt
    
    nx = n_base*2^(k-1);
    %% 剖分
    [hx,Xm] = generate_1D_uniform_grid(xa,xb);
       
    %% 计算
    [error2,errors,uh,solution] = solve_1D_Elliptical(mp);
    error(k,1) = error2;
    error(k,2) = errors;

end


for k = 1:nt-1

    error_order(k,1) = log(error(k,1)/error(k+1,1))/log(2);
    error_order(k,2) = log(error(k,2)/error(k+1,2))/log(2);

end