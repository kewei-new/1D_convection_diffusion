function result = evaluate_local_FEM(x,phi,basis_type,basis_derivative,mp,P,T,Tb)
% Ask for phi or phi_x.
% phi:Store the finite element nodal coefficients calculated by solve_poisson_1D_FEM.m

% nargin
if nargin == 4
    % 当x是输入参考单元时
    mp = length(phi)-1;
    result = 0;
    for i = 1:mp+1
        result = result + phi(i)*basis_reference_FEM(x,basis_derivative,i,basis_type);
    end

else
    % 当x是输入局部单元，此时就需要转换到相应的有限单元上去计算
    % 这应该放在外面，最终和参考单元输入应该是一样的，但是因为还没写计算误差和画图的函数现在这里放一下
    % index:Determine within which finite element x belongs, if x is at the breakpoint, take the left interval
    index = find(((P-x)>=0)==0,1,'last');
    local_points = P(:,T(:,index));
    h_local = max(local_points)-min(local_points);
    mid_point = mean(local_points);
    % x_ref:Scale x from finite element to [-0.5,0.5]
    x_ref = (x-mid_point)/h_local;
    
    result = 0;
    for i = 1:mp+1
    
        result = result + phi(Tb(i,index))*basis_reference_FEM(x_ref,basis_derivative,i,basis_type)*(h_local)^(-basis_derivative);
    
    end

end



