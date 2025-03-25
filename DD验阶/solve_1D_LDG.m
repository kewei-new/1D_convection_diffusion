function [phih,qh]=solve_1D_LDG(Gauss_coefficient,LDG_matrix,Q_matrix,f_values,uh,h,mid_points,mp,ng,num,t)
% phi_xx = xe/epsilon * (n - nd) = f , let q=phi_x
% 以下为推导过程，periodic边界和dirichlet边界的flux不一样
% (q_x,v) = (f,v)
% (q,p) = (phi_x,p)
% -(q,v_x) + (q,v)e = (f,v)
% -(phi,p_x) + (phi,p)e = (q,p)

% 除了两个端点处其余的和periodic一样采用交替，如果Q用+，phi用-则左端点正常
% 右端点采用flux(Q_N) =Q(N-) - (mp/h)*phi(N-) + (mp/h)*g_b,(其中g_b是phi在右边界的值)
% 同时U采用虚拟节点 phi(0-) = g_a(g_a表示phi在左边界的值)
% 由u1h = xe*nh可得
nh = uh;

%% evaluate (f_fun,v) value
% num==0表示正在初始化，初值右端n-nd=0
if num > 0
    f_values = generate_1D_LDG_vector(f_values,Gauss_coefficient,ng,mp,h,mid_points,nh,t);
end
%% 需要组装矩阵联立求方程组 这里少乘了质量矩阵的逆矩阵，需要完善
% 组装矩阵可以放到外面，只需要在里面把边界端元处理的边界值加上


%% boundary treatment
% 也放在外面了


%% 计算phih
f_values = reshape(f_values,(mp+1)*ng,1);
phih = LDG_matrix\f_values;

%% 计算qh
M2 = Q_matrix(:,1:ng*(mp+1));
b2 = Q_matrix(:,end);
qh = M2*phih + b2;

phih = reshape(phih,mp+1,ng);
qh = reshape(qh,mp+1,ng);
